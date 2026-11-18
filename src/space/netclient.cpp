#include "netclient.h"
#include "hostserver.h"
#include "net.h"
#include "world.h"
#include "inventory.h"
#include "voxel.h"
#include "player.h"
#include "scene.h"
#include "simclock.h"
#include "worldtypes.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>

#include <raylib.h>
#include <cstdio>
#include <vector>

extern const char* kGameVersion; // "<branch> <commit>" (#501)

namespace {
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
        NetReader r(m.data.data() + (m.data.empty() ? 0 : 1), m.data.empty() ? 0 : m.data.size() - 1);
        switch (NetPeek(m.data)) {
            case NetMsg::Welcome: {
                MsgWelcome w;
                if (c.state != NetClient::State::AwaitingWelcome || !NetRead(r, w)) return;
                c.welcome = w; c.serverTick = w.serverTick; c.state = NetClient::State::Joined;
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
            case NetMsg::PlayerState: {
                MsgPlayerState s;
                if (c.state != NetClient::State::Joined || !NetRead(r, s)) return;
                if (s.session == c.welcome.session) {
                    if (!c.haveOwn || s.tick > c.own.tick) { c.own = s; c.haveOwn = true; } // (unreliable: may arrive out of order)
                    ++c.states;
                } else {
                    auto& q = c.others[s.session];
                    if (q.empty() || (s.tick > q.back().tick && s.ackInput > q.back().ackInput)) q.push_back(s); // (newer, and moved on: a state per step taken)
                    while (q.size() > 32) q.pop_front();
                }
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
    else { Net_Poll(); Net_Receive(in, 256); }
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
        if (auto l = c.loadouts.find(session); l != c.loadouts.end()) { a.held = l->second.held; a.headlamp = l->second.headlamp; a.suitTier = l->second.suitTier; }
        out.push_back(a);
    }
}

void NetClient_SetLoadout(NetClient& c, const std::string& held, bool headlamp, int suitTier) {
    if (c.state != NetClient::State::Joined) return;
    MsgLoadout l; l.held = held; l.headlamp = headlamp; l.suitTier = (uint8_t)std::clamp(suitTier, 0, 255);
    if (c.loadoutSent && l.held == c.loadout.held && l.headlamp == c.loadout.headlamp && l.suitTier == c.loadout.suitTier) return;
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

void NetClient_SendChat(NetClient& c, const std::string& text) {
    if (c.state != NetClient::State::Joined || text.empty()) return;
    SendMsg(c, MsgChat{ 0, text.substr(0, 256) });
}
