#include "netclient.h"
#include "animals.h"
#include "combat.h"
#include "playersession.h"
#include "hostserver.h"
#include "net.h"
#include "world.h"
#include "inventory.h"
#include "voxel.h"
#include "plots.h"
#include "buildtool.h"
#include "player.h"
#include "scene.h"
#include "simclock.h"
#include "worldtypes.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <unordered_map>

#include <raylib.h>
#include <cstdio>
#include <vector>
#include <set>

extern const char* kGameVersion; // "<branch> <commit>" (#501)

namespace {
    // A player's state from the server - its own message, or decoded from a snapshot (#769): the same handling either way.
    void OnPlayerState(NetClient& c, const MsgPlayerState& s) {
        if (s.session == c.welcome.session) {
            if (!c.haveOwn || s.tick > c.own.tick) { c.own = s; c.haveOwn = true; } // (unreliable: may arrive out of order)
            ++c.states;
        } else {
            auto& q = c.others[s.session];
            if (q.empty() || (s.tick > q.back().tick && s.ackInput > q.back().ackInput)) q.push_back(s); // (newer, and moved on: a state per step taken)
            while (q.size() > 32) q.pop_front();
        }
    }
    template <class M> void SendMsg(const NetClient& c, const M& m) {
        NetWriter w; NetWrite(w, m);
        Net_Send(c.conn, w.bytes.data(), (uint32_t)w.bytes.size(), true);
    }
    std::string Commit() { const std::string v = kGameVersion; const size_t sp = v.rfind(' '); return sp == std::string::npos ? v : v.substr(sp + 1); }
    void End(NetClient& c, NetClient::State s) {
        if (c.conn) Net_Close(c.conn);
        c.conn = 0;
        c.state = s;
    }
    void OnMessage(NetClient& c, const NetMessage& m) {
        c.lastHeard = GetTime();
        const NetMsg kind = NetPeek(m.data);
        if (c.syncHold && (kind == NetMsg::SyncIndex || kind == NetMsg::SyncData || kind == NetMsg::SyncDone)) { c.syncHeld.push_back(m); return; } // (#776)
        NetReader r(m.data.data() + (m.data.empty() ? 0 : 1), m.data.empty() ? 0 : m.data.size() - 1);
        switch (kind) {
            case NetMsg::Welcome: {
                MsgWelcome w;
                if (c.state != NetClient::State::AwaitingWelcome || !NetRead(r, w)) return;
                c.welcome = w; c.serverTick = w.serverTick; c.state = NetClient::State::Joined;
                c.syncing = true; c.syncStartedAt = GetTime(); c.syncDoneAt = -1.0; // (#776) the world's edits come next
                if (gNetRole == NetRole::Client) SetWorldTime(w.worldTime); // #339: the server's clock (spin, orbits, sun, weather derive from it)
                return;
            }
            case NetMsg::Reject: {
                MsgReject j;
                if (!NetRead(r, j)) return;
                c.reject = j;
                End(c, NetClient::State::Rejected);
                return;
            }
            case NetMsg::Bye: { MsgBye b; if (NetRead(r, b)) { c.bye = b.reason; End(c, NetClient::State::Closed); } return; }
            case NetMsg::Pong: {
                MsgPong p;
                if (!NetRead(r, p)) return;
                c.rttMs = (GetTime() - p.clientTime) * 1000.0; c.serverTick = p.serverTick; ++c.pongs;
                { // #339: keep following the server's world clock (slewed: the sky never jumps); a host's own clients only measure it
                    // Only pongs near the best recent round trip: a resent or stalled one makes the trip lopsided (NTP-style).
                    c.rttFloor = c.rttFloor < 0.0 ? c.rttMs : std::min(c.rttMs, c.rttFloor + 2.0);
                    if (c.rttMs > c.rttFloor * 1.25 + 5.0) return;
                    const double est = p.worldTime + 0.0005 * c.rttMs;
                    c.clockError = est - WorldTime(); if (c.pongs > 1) c.worstClockError = std::max(c.worstClockError, std::abs(c.clockError));
                    if (gNetRole == NetRole::Client) WorldTime_Follow(est);
                }
                return;
            }
            case NetMsg::Shower: { MsgShower s; if (NetRead(r, s)) c.shower = s; return; }
            case NetMsg::Chat: { MsgChat m; if (c.state == NetClient::State::Joined && NetRead(r, m)) { c.chat.push_back(m); if (c.chat.size() > 64) c.chat.pop_front(); } return; } // (#326)
            case NetMsg::Roster: {
                MsgRoster ro;
                if (c.state != NetClient::State::Joined || !NetRead(r, ro)) return;
                c.roster.clear();
                for (const auto& e : ro.players) c.roster[e.session] = e.name;
                // Gone: forget their states and loadout (a reconnect starts them fresh).
                for (auto it = c.others.begin(); it != c.others.end();) it = c.roster.count(it->first) ? std::next(it) : c.others.erase(it);
                for (auto it = c.loadouts.begin(); it != c.loadouts.end();) it = c.roster.count(it->first) ? std::next(it) : c.loadouts.erase(it);
                ++c.rosters;
                return;
            }
            case NetMsg::Loadout: {
                MsgLoadout l;
                if (c.state != NetClient::State::Joined || !NetRead(r, l) || l.session == 0) return;
                if (l.session != c.welcome.session) c.loadouts[l.session] = l;
                return;
            }
            case NetMsg::TerrainEdit: { // (#329) in the server's order; a host in this process already changed the world
                MsgTerrainEdit t;
                if (c.state != NetClient::State::Joined || !NetRead(r, t)) return;
                if (c.syncing) { c.syncQueue.push_back(m); ++c.syncQueued; return; } // (#776) after the join snapshot
                ++c.editsSeen;
                c.editOrder = c.editOrder * 1099511628211ull ^ (((uint64_t)t.session << 32) | t.seq); // the order, to compare clients
                if (t.session == c.welcome.session) ++c.editsOwn;
                if (!(t.flags & 3) && t.strength > 0.0f) ++c.placesSeen;
                if (Host_Listening()) return;
                Planetoid* p = Planetoid::findNet(t.body);
                if (!p) return;
                TerrainEdit e{};
                e.center = t.center.toFloat(); e.radius = t.radius; e.strength = t.strength; e.paint = t.paint;
                e.carve = (t.flags & 1) != 0; e.fill = (t.flags & 2) != 0; e.box = (t.flags & 4) != 0; e.softOnly = (t.flags & 8) != 0;
                if (e.box) {
                    e.half = { t.hx, t.hy, t.hz };
                    e.ax = { t.ax[0], t.ax[1], t.ax[2] }; e.ay = { t.ay[0], t.ay[1], t.ay[2] };
                    e.az = Vector3CrossProduct(e.ax, e.ay);
                }
                ApplyTerrainEdit(*p, e);
                return;
            }
            case NetMsg::EditRefused: {
                MsgEditRefused f;
                if (c.state != NetClient::State::Joined || !NetRead(r, f)) return;
                ++c.editsRefused; c.lastRefusal = f.why; c.refusedWhy |= 1u << (unsigned)f.why;
                return;
            }
            case NetMsg::InventoryDelta: { // the server's inventory is the truth: ours follows it exactly (hundredths)
                MsgInventoryDelta d;
                if (c.state != NetClient::State::Joined || !NetRead(r, d)) return;
                ++c.inventoryDeltas;
                if (c.inventory) { if (d.amount > 0.0) c.inventory->addAmount(d.item, d.amount); else c.inventory->removeAmount(d.item, -d.amount); }
                return;
            }
            case NetMsg::PlotEdit: { // (#330) in the server's order; a host in this process already changed the world
                MsgPlotEdit e;
                if (c.state != NetClient::State::Joined || !NetRead(r, e)) return;
                if (c.syncing) { c.syncQueue.push_back(m); ++c.syncQueued; return; } // (#776) after the join snapshot
                ++c.plotEdits;
                c.plotOrder = c.plotOrder * 1099511628211ull ^ (((uint64_t)e.session << 40) ^ ((uint64_t)e.seq << 8) ^ (uint64_t)e.op);
                if (e.session == c.welcome.session) { ++c.plotEditsOwn; if (e.op == NetPartOp::NewPlot) c.lastNewPlot = e.plot; }
                if (c.keepPlotLog) c.plotLog.push_back(e);
                if (Host_Listening()) return;
                Planetoid* p = Planetoid::findNet(e.body);
                if (!p) return;
                if (e.op == NetPartOp::NewPlot) {
                    if (!p->plots) p->plots = std::make_unique<PlotSet>();
                    p->plots->createAt(e.plot, Vector3{ e.anchor[0], e.anchor[1], e.anchor[2] });
                    return;
                }
                Plot* plot = p->plots ? p->plots->byId(e.plot) : nullptr;
                if (!plot) return;
                if (e.op == NetPartOp::Place) {
                    Part part{}; part.type = (PartType)e.part.type; part.origin = { e.part.x, e.part.y, e.part.z }; part.orient = e.part.orient; part.material = e.part.material;
                    BuildTool_ApplyPart(*p, *plot, part);
                } else plot->remove((int)e.partIndex);
                return;
            }
            case NetMsg::SyncIndex: { // (#776) a page of the join snapshot's entries; all in: say which we need
                MsgSyncIndex x;
                if (c.state != NetClient::State::Joined || !NetRead(r, x) || !c.syncing) return;
                if (x.snapshot != c.syncSnapshot || c.syncPages == 0) { c.syncSnapshot = x.snapshot; c.syncPages = x.pages; c.syncPagesGot = 0; c.syncIndex.assign(x.total, NetSyncEntry{}); }
                const size_t base = (size_t)x.page * kSyncIndexPage;
                for (size_t i = 0; i < x.entries.size() && base + i < c.syncIndex.size(); ++i) c.syncIndex[base + i] = x.entries[i];
                if (++c.syncPagesGot < c.syncPages) return;
                // Which do we lack? An entry we already hold with the same hash is skipped (a cached client save).
                const bool install = !Host_Listening(); // (a host in this process: same world; fetched to check, not installed)
                MsgSyncWant want; want.snapshot = c.syncSnapshot; want.total = (uint32_t)c.syncIndex.size();
                want.bits.assign((c.syncIndex.size() + 7) / 8, 0);
                c.syncWanted.assign(c.syncIndex.size(), 0);
                std::map<uint64_t, std::set<uint64_t>> listed; // regions in the snapshot, by body
                for (size_t i = 0; i < c.syncIndex.size(); ++i) {
                    const NetSyncEntry& e = c.syncIndex[i];
                    Planetoid* p = Planetoid::findNet(e.body);
                    uint64_t mine = 0;
                    if (e.kind == NetSyncKind::Region) {
                        listed[e.body].insert(e.key);
                        if (p && p->edits && !c.syncAssumeEmpty) { const std::vector<uint8_t> b = Region_Encode(*p->edits, e.key); if (b.size() > 8) mine = Region_Hash(b); }
                    } else if (p && p->plots && !c.syncAssumeEmpty) mine = Region_Hash(Plots_Blob(*p->plots));
                    if (mine == e.hash) continue;
                    want.bits[i / 8] |= (uint8_t)(1u << (i % 8));
                    c.syncWanted[i] = 1; ++c.syncEntriesWanted; c.syncBytesWanted += e.size;
                }
                // Regions we hold that the server has no edits in: back to the generated ground.
                if (install && !c.syncAssumeEmpty)
                    for (Planetoid* p : Planetoid::All()) {
                        if (!p->edits || !p->netId) continue;
                        for (TerrainDeltas::Key k : Regions_List(*p->edits))
                            if (!listed[p->netId].count(k) && Region_Install(*p, k, nullptr, 0)) ++c.syncCleared;
                    }
                SendMsg(c, want);
                return;
            }
            case NetMsg::SyncData: {
                MsgSyncData d;
                if (c.state != NetClient::State::Joined || !NetRead(r, d) || d.snapshot != c.syncSnapshot || d.entry >= c.syncIndex.size() || !c.syncWanted[d.entry]) return;
                const NetSyncEntry& e = c.syncIndex[d.entry]; // (e.size: as stored; d.size: as sent, compressed)
                std::vector<uint8_t>& buf = c.syncBuf[d.entry];
                if (c.syncGot[d.entry] == 0 && buf.empty()) buf.assign(d.size, 0);
                if (buf.size() != d.size) return; // (pieces of one entry all carry its sent size)
                std::copy(d.bytes.begin(), d.bytes.end(), buf.begin() + d.offset);
                c.syncBytes += d.bytes.size();
                if ((c.syncGot[d.entry] += (uint32_t)d.bytes.size()) < d.size) return;
                // Complete: check it, install it.
                c.syncWanted[d.entry] = 0; ++c.syncEntriesDone;
                std::vector<uint8_t> plain; // (sent DEFLATE-compressed; hashed as stored)
                if (!buf.empty()) { int n = 0; unsigned char* u = DecompressData(buf.data(), (int)buf.size(), &n); if (u) { plain.assign(u, u + n); MemFree(u); } }
                if (Region_Hash(plain) != e.hash) ++c.syncHashBad;
                else if (!Host_Listening()) {
                    if (Planetoid* p = Planetoid::findNet(e.body)) {
                        const bool ok = e.kind == NetSyncKind::Region ? Region_Install(*p, e.key, plain.data(), plain.size()) : Plots_Install(*p, plain);
                        if (ok) ++c.syncInstalled; else ++c.syncHashBad;
                    }
                }
                c.syncBuf.erase(d.entry); c.syncGot.erase(d.entry);
                return;
            }
            case NetMsg::SyncDone: { // everything asked for is here: the live edits held since the snapshot apply now, in order
                MsgSyncDone d;
                if (c.state != NetClient::State::Joined || !NetRead(r, d) || d.snapshot != c.syncSnapshot) return;
                c.syncing = false; c.syncDoneAt = GetTime();
                std::vector<NetMessage> held;
                held.swap(c.syncQueue);
                for (const NetMessage& q : held) OnMessage(c, q);
                return;
            }
            case NetMsg::PartRefused: {
                MsgPartRefused f;
                if (c.state != NetClient::State::Joined || !NetRead(r, f)) return;
                ++c.partsRefused; c.partRefusedWhy |= 1u << (unsigned)f.why;
                return;
            }
            case NetMsg::PlayerState: {
                MsgPlayerState s;
                if (c.state != NetClient::State::Joined || !NetRead(r, s)) return;
                OnPlayerState(c, s);
                return;
            }
            case NetMsg::Animals: { // (#334) the server's animals near us; only a client process keeps copies (a host here owns the real ones)
                MsgAnimals am;
                if (c.state != NetClient::State::Joined || !NetRead(r, am)) return;
                ++c.animalMsgs;
                if (!gAnimalsReplicated) return;
                static std::vector<AnimalNet> in; in.clear();
                for (const MsgAnimals::Entry& e : am.entries) in.push_back({ e.netId, e.species, e.body, e.pos.toFloat(), e.heading, e.state, e.health, e.age, e.flags });
                Animals_ApplyNet(in, am.removed);
                return;
            }
            case NetMsg::Damage: { // (#334 slice 2) a hit on our player (an animal's bite): the damage here (shield, health), the shove the server's
                MsgDamage dm;
                if (c.state != NetClient::State::Joined || !NetRead(r, dm)) return;
                ++c.damageTaken;
                c.lastDamage = dm.amount;
                if (gNetRole != NetRole::Client) return; // (a host's in-process test client: the host's own player isn't the one hit)
                DamageEvent e; e.amount = dm.amount; e.kind = (DamageKind)dm.kind; e.hasSource = dm.hasSource; e.knockback = 0.0f; e.sourceId = dm.sourceId;
                if (Planetoid* b = Planetoid::findNet(dm.body)) e.source = b->toWorldD(dm.local); else e.hasSource = false;
                const DamageResult res = Combat_Damage(LocalSession(), e);
                c.lastDamageBlocked = res.blocked;
                return;
            }
            case NetMsg::PlayerSnapshot: { // (#769) the states of one snapshot, each against one we acknowledged (or full)
                if (c.state != NetClient::State::Joined) return;
                MsgPlayerSnapshot snap;
                const bool ok = NetRead(r, snap, [&](uint32_t session, uint32_t seq, PlayerStateQ& out) {
                    auto h = c.snapHistory.find(session);
                    if (h == c.snapHistory.end()) return false;
                    for (const auto& [s, q] : h->second) if (s == seq) { out = q; return true; }
                    return false;
                });
                if (!ok) { ++c.snapUndecodable; return; } // (a base we no longer have: not acked, so the server keeps to older bases or resends full)
                // Acknowledge it (the newest sequence and a mask of the 32 before), keep its states as future bases, apply them.
                if (snap.seq > c.snapNewest) {
                    const uint32_t d = snap.seq - c.snapNewest;
                    c.snapMask = d >= 33 ? 0u : (d == 32 ? 0u : (c.snapMask << d)) | (c.snapNewest && d <= 32 ? 1u << (d - 1) : 0u);
                    c.snapNewest = snap.seq;
                } else if (snap.seq < c.snapNewest && c.snapNewest - snap.seq <= 32) c.snapMask |= 1u << (c.snapNewest - snap.seq - 1);
                for (const MsgPlayerSnapshot::Entry& e : snap.entries) {
                    auto& h = c.snapHistory[e.s.session];
                    h.push_back({ snap.seq, e.s });
                    while (h.size() > 32) h.pop_front();
                    OnPlayerState(c, NetDequantize(e.s, snap.tick));
                }
                ++c.snapshots;
                return;
            }
            default: return;
        }
    }
}

const char* NetClient_StateName(NetClient::State s) {
    switch (s) {
        case NetClient::State::Idle: return "idle"; case NetClient::State::Connecting: return "connecting";
        case NetClient::State::AwaitingWelcome: return "awaiting welcome"; case NetClient::State::Joined: return "joined";
        case NetClient::State::Rejected: return "rejected"; case NetClient::State::Closed: return "closed";
    }
    return "?";
}

bool NetClient_Connect(NetClient& c, const char* address, const std::string& name, const std::string& password) {
    if (!Host_Listening()) { std::string err; if (!Net_Init(err)) { c.state = NetClient::State::Closed; return false; } }
    c.name = name; c.password = password;
    c.conn = Net_Connect(address);
    if (!c.conn) { c.state = NetClient::State::Closed; return false; }
    if (Host_Listening()) Host_SetLocalClient(c.conn); // a host in this process receives for us
    c.state = NetClient::State::Connecting;
    c.started = c.lastHeard = GetTime();
    return true;
}

void NetClient_Update(NetClient& c) {
    if (c.state == NetClient::State::Idle || c.state == NetClient::State::Rejected || c.state == NetClient::State::Closed) return;
    if (!c.syncHold && !c.syncHeld.empty()) { // (#776) the world is ready now: the sync's messages held so far, in order
        std::vector<NetMessage> held;
        held.swap(c.syncHeld);
        for (const NetMessage& m : held) OnMessage(c, m);
    }
    // Each other player's playback clock, in that player's own steps (ackInput: the server steps a player only by its
    // inputs, so by the server's tick it moves in bursts as they arrive; by its own steps evenly, at 60 a second). It runs
    // at real time and eases toward kInterpDelay behind the newest state heard (a 1 s time constant, per second: per
    // frame it was frame-rate bound); 2 s off, it jumps (half a second snapped after one long frame: the ease catches up at most ~1.5x).
    for (auto it = c.viewTicks.begin(); it != c.viewTicks.end();) it = c.others.count(it->first) ? std::next(it) : c.viewTicks.erase(it);
    for (const auto& [session, q] : c.others) {
        if (q.empty() || q.back().ackInput == 0) continue;
        const double target = q.back().ackInput - kInterpDelay * kSimHz;
        double& v = c.viewTicks[session];
        v += GetFrameTime() * kSimHz;
        if (v <= 0.0 || fabs(v - target) > 120.0) v = target;
        else v += (target - v) * (1.0 - exp(-GetFrameTime() / 1.0));
    }
    std::vector<NetMessage> in;
    if (Host_Listening()) in = Host_TakeLocalClientMessages(c.conn);
    else { // several clients in one process (#348: a bot swarm): one poll, each message kept for its own client
        static std::unordered_map<uint32_t, std::vector<NetMessage>> inbox;
        Net_Poll();
        std::vector<NetMessage> got; Net_Receive(got, 1024);
        for (NetMessage& m : got) inbox[m.conn].push_back(std::move(m));
        auto it = inbox.find(c.conn); if (it != inbox.end()) { in.swap(it->second); inbox.erase(it); }
    }
    for (const NetMessage& m : in) {
        if (m.conn != c.conn) continue;
        OnMessage(c, m);
        if (c.state == NetClient::State::Rejected || c.state == NetClient::State::Closed) return;
    }
    const double now = GetTime();
    if (c.state == NetClient::State::Connecting) {
        if (Net_IsConnected(c.conn)) {
            SendMsg(c, MsgHello{ c.protocol, c.build.empty() ? Commit() : c.build, c.name, c.password });
            c.state = NetClient::State::AwaitingWelcome;
            c.lastHeard = now;
        } else if (now - c.started > kClientConnectTimeout) { c.bye = NetBye::TimedOut; End(c, NetClient::State::Closed); }
        return;
    }
    if (c.state == NetClient::State::Joined && now - c.lastPing >= 1.0) { SendMsg(c, MsgPing{ ++c.pingSeq, now }); c.lastPing = now; }
    if (now - c.lastHeard > kClientSilentTimeout || !Net_IsConnected(c.conn)) { c.bye = NetBye::TimedOut; End(c, NetClient::State::Closed); }
}

void NetClient_Disconnect(NetClient& c) {
    if (!c.conn) return;
    SendMsg(c, MsgBye{ NetBye::Quit });
    c.bye = NetBye::Quit;
    End(c, NetClient::State::Closed);
}

// ---- Movement (#327) ----

MsgPlayerState NetClient_ToMsg(const Player::NetState& n) {
    MsgPlayerState s;
    s.body = n.body; s.local = n.local;
    s.vx = n.vel.x; s.vy = n.vel.y; s.vz = n.vel.z;
    s.jumpBuffer = n.jumpBuffer; s.coyote = n.coyote; s.jumpLock = n.jumpLock; s.jetpackFuel = n.jetpackFuel;
    s.flags = (n.grounded ? 1 : 0) | (n.hovering ? 16 : 0) | ((n.posture & 3) << 5);
    return s;
}
Player::NetState NetClient_FromMsg(const MsgPlayerState& s) {
    Player::NetState n;
    n.body = s.body; n.local = s.local;
    n.vel = { s.vx, s.vy, s.vz };
    n.grounded = (s.flags & 1) != 0; n.hovering = (s.flags & 16) != 0; n.posture = (uint8_t)((s.flags >> 5) & 3);
    n.jumpBuffer = s.jumpBuffer; n.coyote = s.coyote; n.jumpLock = s.jumpLock; n.jetpackFuel = s.jetpackFuel;
    return n;
}

void NetClient_Step(NetClient& c, Scene& world, Player& self, const PlayerInput& in) {
    if (c.state != NetClient::State::Joined) return;
    // Reconcile first: where were we, by our own reckoning, when the server applied that input?
    if (c.haveOwn) {
        c.haveOwn = false;
        Planetoid* body = Planetoid::findNet(c.own.body);
        while (!c.history.empty() && c.history.front().tick < c.own.ackInput) c.history.pop_front();
        if (body && !c.history.empty() && c.history.front().tick == c.own.ackInput) {
            const DVec3 d = c.history.front().local - c.own.local;
            c.lastError = d.length();
            c.worstError = std::max(c.worstError, c.lastError);
            c.history.pop_front(); // acknowledged
            if (c.lastError > kReconcileError) {
                ++c.corrections;
                // Back to the server's state, then replay what it hasn't applied yet, each step with the view it was
                // sent with (the view now is the newest; replaying with it turned every replayed step).
                const Vector3 viewLocal = body->dirToLocal(self.forward);
                self.applyNetState(world, NetClient_FromMsg(c.own));
                for (NetPredicted& h : c.history) {
                    self.forward = body->dirToWorld({ h.fx, h.fy, h.fz });
                    RenderOriginScope frame(world, self.position);
                    self.update(world, h.in, SimDt());
                    h.local = body->toLocalD(self.position);
                }
                self.forward = body->dirToWorld(viewLocal);
            }
        }
    }
    // Then this step: send it, and predict it.
    Planetoid* body = self.nearest;
    NetPredicted step;
    step.tick = ++c.tick; step.in = in;
    if (body) { const Vector3 f = body->dirToLocal(self.forward); step.fx = f.x; step.fy = f.y; step.fz = f.z; }
    {
        RenderOriginScope frame(world, self.position);
        self.update(world, in, SimDt());
    }
    if (body) step.local = body->toLocalD(self.position);
    c.history.push_back(step);
    while (c.history.size() > 600) c.history.pop_front(); // 10 s unacknowledged: something is wrong; don't grow forever
    // The newest kInputRedundancy steps, unreliable (a lost packet is covered by the next).
    MsgInputs m;
    m.action = (uint8_t)std::clamp(c.action, 0, 3); m.actionPhase = std::clamp(c.actionPhase, 0.0f, 1.0f);
    m.ackSnap = c.snapNewest; m.ackMask = c.snapMask; // (#769)
    for (auto it = c.history.rbegin(); it != c.history.rend() && (int)m.steps.size() < kInputRedundancy; ++it) {
        if (!m.steps.empty() && it->tick != m.steps.back().tick - 1) break; // (consecutive only)
        m.steps.push_back({ it->tick, it->in.held, it->in.pressed, it->fx, it->fy, it->fz });
    }
    NetWriter w; NetWrite(w, m);
    Net_Send(c.conn, w.bytes.data(), (uint32_t)w.bytes.size(), false);
}

bool NetClient_OtherState(const NetClient& c, uint32_t session, MsgPlayerState& out) {
    auto it = c.others.find(session);
    if (it == c.others.end() || it->second.empty()) return false;
    const auto& q = it->second;
    // kInterpDelay behind the newest state, in the player's own steps (its playback clock, NetClient_Update).
    auto v = c.viewTicks.find(session);
    const double t = v != c.viewTicks.end() ? v->second : (double)q.back().ackInput - kInterpDelay * kSimHz;
    if (t <= q.front().ackInput) { out = q.front(); return true; }
    for (size_t i = 1; i < q.size(); ++i) {
        if (q[i].ackInput < t) continue;
        const MsgPlayerState& a = q[i - 1];
        const MsgPlayerState& b = q[i];
        out = b; // flags, posture, the hands: the state being moved toward
        if (a.body != b.body) return true; // (a new body: no blending across frames)
        const double k = b.ackInput == a.ackInput ? 1.0 : (t - a.ackInput) / (double)(b.ackInput - a.ackInput);
        const float kf = (float)k;
        out.local = a.local + (b.local - a.local) * k;
        out.vx = a.vx + (b.vx - a.vx) * kf; out.vy = a.vy + (b.vy - a.vy) * kf; out.vz = a.vz + (b.vz - a.vz) * kf;
        float dyaw = b.yaw - a.yaw; // the short way round
        if (dyaw > PI) dyaw -= 2.0f * PI; else if (dyaw < -PI) dyaw += 2.0f * PI;
        out.yaw = a.yaw + dyaw * kf;
        out.pitch = a.pitch + (b.pitch - a.pitch) * kf;
        return true;
    }
    // Past the newest (the next state is late): carry on along its velocity for up to kExtrapolate, then hold. Holding at
    // once stopped the player and then jumped them ahead when the late state came in.
    out = q.back();
    const double ahead = std::min((t - q.back().ackInput) / kSimHz, kExtrapolate);
    if (ahead > 0.0) out.local = out.local + DVec3{ out.vx * ahead, out.vy * ahead, out.vz * ahead };
    return true;
}

bool NetClient_OtherView(const NetClient& c, uint32_t session, uint64_t& body, DVec3& local) {
    MsgPlayerState s;
    if (!NetClient_OtherState(c, session, s)) return false;
    body = s.body; local = s.local;
    return true;
}

void NetClient_FillAvatars(const NetClient& c, Scene& world, std::vector<AvatarState>& out) {
    if (c.state != NetClient::State::Joined) return;
    for (const auto& [session, q] : c.others) {
        if (session == c.welcome.session) continue;
        MsgPlayerState s;
        if (!NetClient_OtherState(c, session, s)) continue;
        Planetoid* b = Planetoid::findNet(s.body);
        if (!b) continue;
        AvatarState a;
        a.session = session;
        if (auto r = c.roster.find(session); r != c.roster.end()) a.name = r->second;
        a.position = b->toWorldD(s.local);
        const double len = s.local.length();
        a.up = len > 1e-6 ? b->dirToWorld(Vector3Normalize(s.local.toFloat())) : Vector3{ 0, 1, 0 };
        // The view from its body-local angles, then flattened onto the ground for the body; its tilt is the look pitch.
        const Vector3 viewLocal{ sinf(s.yaw) * cosf(s.pitch), sinf(s.pitch), cosf(s.yaw) * cosf(s.pitch) };
        const Vector3 view = b->dirToWorld(viewLocal);
        const Vector3 flat = Vector3Subtract(view, Vector3Scale(a.up, Vector3DotProduct(view, a.up)));
        if (Vector3Length(flat) > 1e-4f) a.facing = Vector3Normalize(flat);
        a.lookPitch = asinf(std::clamp(Vector3DotProduct(view, a.up), -1.0f, 1.0f));
        const Vector3 vl{ s.vx, s.vy, s.vz };
        a.velocity = Vector3RotateByQuaternion(vl, b->rotationQuat()); // (relative to the ground already)
        a.grounded = (s.flags & 1) != 0; a.flying = (s.flags & 2) != 0 || (s.flags & 16) != 0; a.swimming = (s.flags & 8) != 0;
        a.jetpack = (s.flags & 128) != 0;
        a.posture = (Player::Posture)std::min((s.flags >> 5) & 3, (int)Player::Posture::Prone);
        a.action = s.action; a.actionPhase = s.actionPhase;
        if (auto l = c.loadouts.find(session); l != c.loadouts.end()) {
            a.held = l->second.held; a.headlamp = l->second.headlamp; a.suitTier = l->second.suitTier;
            const uint32_t sc = l->second.suitColor;
            if (sc) a.suitColor = Color{ (unsigned char)(sc >> 16), (unsigned char)(sc >> 8), (unsigned char)sc, 255 };
        }
        out.push_back(a);
    }
}

void NetClient_SetLoadout(NetClient& c, const std::string& held, bool headlamp, int suitTier) {
    if (c.state != NetClient::State::Joined) return;
    MsgLoadout l; l.held = held; l.headlamp = headlamp; l.suitTier = (uint8_t)std::clamp(suitTier, 0, 255);
    l.suitColor = c.suitColor & 0xFFFFFFu;
    if (c.loadoutSent && l.held == c.loadout.held && l.headlamp == c.loadout.headlamp && l.suitTier == c.loadout.suitTier && l.suitColor == c.loadout.suitColor) return;
    c.loadout = l; c.loadoutSent = true;
    SendMsg(c, l);
}

NetClient* gNetClient = nullptr;

uint32_t NetClient_RequestEdit(NetClient& c, NetEditKind kind, const Planetoid& body, Vector3 hit, Vector3 centre, float radius, const std::string& item) {
    if (c.state != NetClient::State::Joined) return 0;
    MsgEditRequest q;
    q.seq = ++c.editSeq; q.kind = kind; q.body = body.netId;
    q.hit = DVec3(hit); q.center = DVec3(centre); q.radius = radius; q.item = item;
    SendMsg(c, q);
    return q.seq;
}

uint32_t NetClient_RequestPart(NetClient& c, const Planetoid& body, const Plot* plot, Vector3 anchor, const Part& part, bool remove, int partIndex) {
    if (c.state != NetClient::State::Joined) return 0;
    MsgPartRequest q;
    q.seq = ++c.partSeq; q.op = remove ? NetPartOp::Remove : NetPartOp::Place; q.body = body.netId;
    q.plot = plot ? plot->id : 0;
    q.anchor[0] = anchor.x; q.anchor[1] = anchor.y; q.anchor[2] = anchor.z;
    q.part = NetPart{ (uint8_t)part.type, part.origin.x, part.origin.y, part.origin.z, part.orient, part.material };
    q.partIndex = partIndex < 0 ? 0u : (uint32_t)partIndex;
    if (remove && !plot) return 0;
    SendMsg(c, q);
    return q.seq;
}

void NetClient_SendChat(NetClient& c, const std::string& text) {
    if (c.state != NetClient::State::Joined || text.empty()) return;
    SendMsg(c, MsgChat{ 0, text.substr(0, 256) });
}

// ---- Join sync (#776) helpers ----

// A body's plots as the save writes them (text; the same file mode both ways, so equal plots give equal bytes).
std::vector<uint8_t> Plots_Blob(const PlotSet& plots) {
    static const char* kTemp = "saves/.netsync_plots.tmp";
    if (!DirectoryExists("saves")) MakeDirectory("saves");
    std::vector<uint8_t> out;
    FILE* f = fopen(kTemp, "w");
    if (!f) return out;
    plots.save(f);
    fclose(f);
    if ((f = fopen(kTemp, "rb"))) {
        uint8_t buf[4096];
        size_t n;
        while ((n = fread(buf, 1, sizeof buf, f)) > 0) out.insert(out.end(), buf, buf + n);
        fclose(f);
    }
    std::remove(kTemp);
    return out;
}

bool Plots_Install(Planetoid& body, const std::vector<uint8_t>& blob) {
    static const char* kTemp = "saves/.netsync_plots_in.tmp";
    if (!DirectoryExists("saves")) MakeDirectory("saves");
    FILE* f = fopen(kTemp, "wb");
    if (!f) return false;
    fwrite(blob.data(), 1, blob.size(), f);
    fclose(f);
    auto set = std::make_unique<PlotSet>();
    bool ok = false;
    if ((f = fopen(kTemp, "r"))) { ok = set->load(f); fclose(f); }
    std::remove(kTemp);
    if (ok) body.plots = std::move(set);
    return ok;
}

uint64_t NetSync_TerrainFingerprint(Scene& world) {
    uint64_t h = 14695981039346656037ull;
    for (auto& o : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(o.get());
        if (!p || !p->edits) continue;
        for (TerrainDeltas::Key k : Regions_List(*p->edits)) {
            h = (h ^ p->netId) * 1099511628211ull;
            h = (h ^ k) * 1099511628211ull;
            h = (h ^ Region_Hash(Region_Encode(*p->edits, k))) * 1099511628211ull;
        }
    }
    return h;
}

// The client's cache (#776): its copy of the server's terrain edits, kept beside its client save (the regions as the
// save writes them), so a rejoin fetches only what changed. Keyed by the server's world (seed and bodies fingerprint):
// another world's cache is ignored. Bodies are numbered as the save numbers them.
static bool CachedBody(const Planetoid* p) { return p && p->ringRock < 0 && !p->vehicle; }
int NetClient_LoadCache(Scene& world, const std::string& save, uint32_t seed, uint64_t fingerprint) {
    FILE* f = fopen((save + ".synced").c_str(), "r");
    if (!f) return 0;
    unsigned s = 0; unsigned long long fp = 0;
    const bool same = fscanf(f, "%u %llu", &s, &fp) == 2 && s == seed && fp == (unsigned long long)fingerprint;
    fclose(f);
    if (!same) return 0;
    const std::string dir = save + ".regions";
    int index = 0, regions = 0;
    for (auto& o : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(o.get());
        if (!CachedBody(p)) continue;
        std::shared_ptr<const TerrainDeltas> d = LoadRegionsFor(index++, dir, 0);
        if (!d || d->empty()) continue;
        regions += (int)Regions_List(*d).size();
        p->edits = d;
        p->regionsSavedTo = dir; p->allRegionsDirty = false; p->dirtyRegions.clear();
        p->lodDirty = true;
    }
    return regions;
}
int NetClient_SaveCache(Scene& world, const std::string& save, uint32_t seed, uint64_t fingerprint) {
    if (!DirectoryExists("saves")) MakeDirectory("saves");
    const std::string dir = save + ".regions";
    if (!DirectoryExists(dir.c_str())) MakeDirectory(dir.c_str());
    int index = 0, written = 0;
    for (auto& o : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(o.get());
        if (!CachedBody(p) || !p->edits) { if (CachedBody(p)) ++index; continue; }
        written += SaveRegionsFor(*p->edits, p->dirtyRegions, p->allRegionsDirty, p->regionsSavedTo, index++, dir);
    }
    if (FILE* f = fopen((save + ".synced").c_str(), "w")) { fprintf(f, "%u %llu\n", seed, (unsigned long long)fingerprint); fclose(f); }
    return written;
}

int NetSync_ScatterStrokes(Scene& world, int strokes, uint32_t seed) {
    Planetoid* p = nullptr;
    for (auto& o : world.objects) if ((p = dynamic_cast<Planetoid*>(o.get()))) break;
    const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
    if (!pf || strokes <= 0) return 0;
    const Vector3 home = pf->findSpawnDir(p->seed);
    Vector3 e1 = Vector3CrossProduct({ 0, 1, 0 }, home);
    e1 = Vector3Normalize(Vector3Length(e1) < 1e-3f ? Vector3{ 1, 0, 0 } : e1);
    const Vector3 e2 = Vector3CrossProduct(home, e1);
    std::vector<TerrainEdit> edits;
    uint32_t s = seed * 2654435761u + 1u;
    auto next = [&]() { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return s; }; // xorshift: the same strokes everywhere
    for (int i = 0; i < strokes; ++i) {
        const float a = (float)(next() % 60001) / 100.0f - 300.0f, b = (float)(next() % 60001) / 100.0f - 300.0f;
        const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(e1, a / pf->radius), Vector3Scale(e2, b / pf->radius))));
        TerrainEdit e{}; e.center = Vector3Scale(d, pf->radius + pf->elevation(d)); e.radius = 1.5f + (float)(next() % 100) / 100.0f; e.carve = true;
        edits.push_back(e);
    }
    ApplyTerrainEdits(*p, edits.data(), edits.size());
    return (int)edits.size();
}

void NetClient_HitAnimal(NetClient& c, uint64_t netId) {
    if (c.state != NetClient::State::Joined || !netId) return;
    SendMsg(c, MsgAnimalHit{ netId });
}
void NetClient_UseAnimal(NetClient& c, uint64_t netId, int action) {
    if (c.state != NetClient::State::Joined || !netId || action < 1 || action > 2) return;
    SendMsg(c, MsgAnimalHit{ netId, (uint8_t)action });
}
