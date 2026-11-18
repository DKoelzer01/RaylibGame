#include "hostserver.h"
#include "netinterest.h"
#include "netguard.h"
#include "gputimer.h"
#include "events.h"
#include "joininfo.h"
#include "net.h"
#include "netproto.h"
#include "scene.h"
#include "simclock.h"
#include "survival.h"
#include "world.h"
#include "worldgenparams.h"
#include "player.h"
#include "netclient.h"
#include "terraintool.h"
#include "buildtool.h"
#include "plots.h"
#include "gamerng.h"
#include "tools.h"
#include "inventory.h"
#include <cmath>
#include <map>
#include "playerinput.h"
#include "worldtypes.h"
#include "voxel.h"
#include <raymath.h>
#include <deque>

#include <raylib.h>
#include <algorithm>
#include <functional>
#include <cstdio>
#include <set>
#include <string>
#include <unordered_map>

extern const char* kGameVersion; // "<branch> <commit>" (cmake/game_version.cmake, #501)

namespace {
    bool listening = false;
    int joins = 0;
    std::set<uint32_t> localClients;                   // tests: client connections opened in this process
    std::unordered_map<uint32_t, std::vector<NetMessage>> localInbox;
    double autosaveSeconds = 0.0, lastSave = -1.0;
    int autosaves = 0;
    std::string password;
    int maxPlayers = 8;

    // The connection lifecycle (#321). A peer is a connection; a session is a player, kept by name for the server's
    // run, so a reconnect under the same name gets its session (and later its PlayerSession) back.
    struct Peer { uint32_t session = 0; bool joined = false; double lastHeard = 0.0; NetGuard guard; }; // guard: #349
    std::unordered_map<uint32_t, Peer> peers;          // by connection
    struct Session {
        uint32_t id = 0; std::string name; uint32_t conn = 0; bool online = false;
        MsgLoadout loadout;          // what the player holds and wears (#750), kept across a reconnect
    };
    std::vector<Session> sessions;
    // #326 / #345: the dedicated server's hooks (main.cpp sets them in --server mode).
    std::string motd;
    std::function<bool(const std::string& name, const std::string& address, std::string& why)> joinCheck;
    std::function<void(uint32_t session, const std::string& name, const std::string& text)> chatHandler;
    std::function<void(const char* what, uint32_t session, const std::string& detail)> eventHook;
    uint32_t nextSession = 1;
    constexpr double kTimeoutSeconds = 15.0;           // nothing heard (pings come every second): gone

    // Players (#327): the server's own simulation of every session's player, from the inputs its client sends (one a sim
    // step; none waiting: the last keys stay held). Kept across a reconnect: the player is where they left.
    struct Remote {
        Player player;
        std::deque<MsgInput> inputs;
        uint32_t lastApplied = 0;   // the client tick of the last input stepped (the owner replays the ones after)
        uint32_t lastQueued = 0;    // the newest input tick received (redundant copies of older ones are skipped)
        uint64_t carried = 0;       // presses of inputs dropped from a full queue, for the next step
        double budget = 0.0;        // steps it may take (earned at real time, capped: a fast client can't run ahead)
        long long tick = 0;         // the server sim tick the budget has been earned to
        Vector3 forwardLocal{ 0, 0, 1 };
        uint8_t action = 0; float actionPhase = 0.0f; // the hands, from the newest input packet (relayed, #750)
        // Terrain edits (#329): the player's inventory as the server keeps it (yields go in, placing is paid from it; the
        // client mirrors it from InventoryDelta), and a stroke budget (a token per tool cooldown, a short burst allowed).
        Inventory inventory;
        double editTokens = 2.0, editTokensAt = -1.0;
        double partTokens = 4.0, partTokensAt = -1.0; // building (#330): placing and breaking, a few a second
    };
    std::unordered_map<uint32_t, Remote> remotes;      // by session
    long long lastSnapshot = 0;
    constexpr int kSnapshotEvery = 2;                  // sim steps (60 Hz) per snapshot: 30 Hz (ADR #309)
    constexpr size_t kMaxQueuedInputs = 120;           // 2 s of inputs waiting: more and the oldest go (keeping their presses)
    constexpr double kMaxBudget = 30.0;                // steps a player may catch up at once (0.5 s of hitch or burst)
    const int kMaxStepsPerFrame = getenv("RAYLIBGAME_REMOTE_STEPS_PER_FRAME") ? std::max(1, atoi(getenv("RAYLIBGAME_REMOTE_STEPS_PER_FRAME"))) : 2; // (#768) catch-up a frame

    Planetoid* StartBody(Scene& world) {
        for (auto& o : world.objects) if (auto* p = dynamic_cast<Planetoid*>(o.get())) return p;
        return nullptr;
    }
    void SpawnRemote(Scene& world, uint32_t session) {
        Planetoid* body = StartBody(world);
        const auto* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
        if (!body || !pf) return;
        // Beside the start, 3 m apart by session (never inside each other).
        const Vector3 start = pf->findSpawnDir(body->seed);
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, start));
        const Vector3 dir = Vector3Normalize(Vector3Add(start, Vector3Scale(east, 3.0f * (float)session / pf->radius)));
        Remote& r = remotes[session];
        r.player.spawnOn(world, *body, body->dirToWorld(dir));
        r.player.flyAllowed = false; // (#781) fly mode is noclip: never for a client's server-side player (admins later, #345)
        r.tick = SimTick();
    }
    // A player is stepped only by its own inputs, one step each, in order (as Source runs a player's usercmds): the server
    // walks exactly the path the client predicted, whatever the arrival times. No input waiting: the player holds still in
    // its body's frame (the ground or air anchor carries it with the planet) until the next arrives. A hitch or a burst is
    // caught up at once, but only as fast as real time allows (a budget earned per server tick, at most kMaxBudget), so
    // a client sending too many inputs can't run faster than anyone else. (Repeating the keys while an input was late,
    // then dropping one to catch up, put the server off the client's path after every long frame: metres, with jumps.)
    void StepRemotes(Scene& world) {
        for (auto& [session, r] : remotes) {
            r.budget = std::min(kMaxBudget, r.budget + (double)(SimTick() - r.tick));
            r.tick = SimTick();
            while (r.inputs.size() > kMaxQueuedInputs) { r.carried |= r.inputs.front().pressed; r.inputs.pop_front(); }
            // (#768) at most kMaxStepsPerFrame inputs a frame per player: after a slow frame the backlog is worked off over the
            // next frames, not all at once (all at once made the next frame slower still: the 31-40 ms catch-up bursts at 32 players).
            int stepsThisFrame = 0;
            while (r.budget >= 1.0 && !r.inputs.empty() && stepsThisFrame++ < kMaxStepsPerFrame) {
                r.budget -= 1.0;
                const MsgInput m = r.inputs.front();
                r.inputs.pop_front();
                PlayerInput in;
                in.enabled = true; in.held = m.held; in.pressed = m.pressed | r.carried;
                r.carried = 0;
                r.lastApplied = m.tick;
                r.forwardLocal = { m.fx, m.fy, m.fz };
                if (Planetoid* b = r.player.nearest ? r.player.nearest : StartBody(world)) r.player.forward = b->dirToWorld(r.forwardLocal);
                RenderOriginScope frame(world, r.player.position);
                r.player.update(world, in, SimDt());
            }
        }
    }
    MsgPlayerState StateOf(Scene& world, uint32_t session, const Remote& r) {
        MsgPlayerState s;
        Player::NetState ns;
        if (r.player.netState(world, ns)) s = NetClient_ToMsg(ns);
        else if (Planetoid* b = StartBody(world)) { s.body = b->netId; s.local = b->toLocalD(r.player.position); }
        s.session = session; s.tick = (uint32_t)SimTick(); s.ackInput = r.lastApplied;
        s.flags |= (r.player.flyMode ? 2 : 0) | (r.player.inWater ? 8 : 0) | (r.player.jetpackFiring ? 128 : 0);
        s.action = r.action; s.actionPhase = r.actionPhase;
        // The view, body-local, as spherical angles (yaw about the body's y axis from +z, pitch from its xz plane): the
        // others turn it back with the body (avatars, #328).
        s.yaw = atan2f(r.forwardLocal.x, r.forwardLocal.z);
        if (s.yaw < 0.0f) s.yaw += 2.0f * PI;
        s.pitch = asinf(std::clamp(r.forwardLocal.y, -1.0f, 1.0f));
        return s;
    }
    // #769 delta snapshots, per client: what each snapshot carried (to turn acks into bases), and per player the newest
    // state the client acknowledged (seq, state).
    constexpr uint32_t kDeltaRing = 64, kDeltaMaxAge = 48, kDeltaResync = 45; // snapshots (30 a second)
    struct SnapSent { uint32_t seq = 0; std::vector<PlayerStateQ> states; };
    struct DeltaState {
        uint32_t nextSeq = 1, lastAck = 0;
        SnapSent ring[kDeltaRing];
        std::unordered_map<uint32_t, std::pair<uint32_t, PlayerStateQ>> base;
    };
    std::unordered_map<uint32_t, DeltaState> deltas; // by connection

    void SendSnapshots(Scene& world) {
        if (SimTick() - lastSnapshot < kSnapshotEvery) return;
        lastSnapshot = SimTick();
        // Every online player's state, written once (#348: it was written again for each client); then each client gets
        // what interest management picks for it (#340: relevancy by distance and body, a byte budget; nothing starves).
        // RAYLIBGAME_NET_INTEREST=0: every state to every client, as before (A/B).
        static const bool interest = !(getenv("RAYLIBGAME_NET_INTEREST") && atoi(getenv("RAYLIBGAME_NET_INTEREST")) == 0);
        struct Out { uint32_t session; MsgPlayerState s; NetWriter w; };
        std::vector<Out> outs;
        for (const auto& [session, r] : remotes) {
            bool online = false;
            for (const Session& s : sessions) online = online || (s.id == session && s.online);
            if (!online) continue;
            Out o{ session, StateOf(world, session, r), {} };
            NetWrite(o.w, o.s);
            outs.push_back(std::move(o));
        }
        // #769: the chosen states go to each client in ONE snapshot message, each delta-coded against the newest state of that
        // player the client acknowledged (MsgInputs carries the acks) - else full. RAYLIBGAME_NET_DELTA=0: a message per state.
        static const bool delta = !(getenv("RAYLIBGAME_NET_DELTA") && atoi(getenv("RAYLIBGAME_NET_DELTA")) == 0);
        std::vector<PlayerStateQ> quant;
        if (delta) for (const Out& o : outs) quant.push_back(NetQuantize(o.s));
        std::vector<InterestItem> items;
        std::vector<uint32_t> chosen;
        for (const auto& [conn, p] : peers) {
            if (!p.joined) continue;
            const Out* self = nullptr;
            for (const Out& o : outs) if (o.session == p.session) self = &o;
            DeltaState& ds = deltas[conn];
            if (delta && ds.nextSeq - ds.lastAck > kDeltaResync) ds.base.clear(); // (no acks for a while: start over from full states)
            chosen.clear();
            if (!interest || !self) for (const Out& o : outs) chosen.push_back(o.session); // (no own state yet: everything)
            else {
                items.clear();
                for (const Out& o : outs) {
                    const uint32_t bytes = !delta ? (uint32_t)o.w.bytes.size() + 28 : ds.base.count(o.session) ? 12u : 34u; // (+ UDP/GNS overhead per message without deltas)
                    items.push_back({ o.session, o.s.body, o.s.local, bytes, o.session == p.session });
                }
                const Planetoid* vb = Planetoid::findNet(self->s.body);
                chosen = Interest_Select(conn, self->s.body, self->s.local, vb ? vb->scale : 1.0, items);
            }
            if (!delta) {
                for (uint32_t id : chosen)
                    for (const Out& o : outs) if (o.session == id) { Net_Send(conn, o.w.bytes.data(), (uint32_t)o.w.bytes.size(), false); break; } // unreliable: the next one replaces it
                continue;
            }
            MsgPlayerSnapshot snap; snap.seq = ds.nextSeq++; snap.tick = (uint32_t)SimTick();
            SnapSent& sent = ds.ring[snap.seq % kDeltaRing];
            sent.seq = snap.seq; sent.states.clear();
            for (uint32_t id : chosen)
                for (size_t k = 0; k < outs.size(); ++k) if (outs[k].session == id) {
                    MsgPlayerSnapshot::Entry e; e.s = quant[k];
                    auto b = ds.base.find(id);
                    if (b != ds.base.end() && snap.seq - b->second.first < kDeltaMaxAge) { e.baseSeq = b->second.first; e.base = b->second.second; }
                    snap.entries.push_back(e);
                    sent.states.push_back(quant[k]);
                    break;
                }
            NetWriter w; NetWrite(w, snap);
            Net_Send(conn, w.bytes.data(), (uint32_t)w.bytes.size(), false); // unreliable: acked back, and the next one replaces it
        }
    }
    // #769: what a client acknowledged moves its bases on (the newest state of each player it has).
    void OnSnapshotAck(uint32_t conn, uint32_t ackSnap, uint32_t ackMask) {
        auto it = deltas.find(conn);
        if (it == deltas.end() || ackSnap == 0) return;
        DeltaState& ds = it->second;
        ds.lastAck = std::max(ds.lastAck, ackSnap);
        for (int i = -1; i < 32; ++i) {
            if (i >= 0 && !((ackMask >> i) & 1u)) continue;
            const uint32_t seq = i < 0 ? ackSnap : ackSnap - 1 - (uint32_t)i;
            if (seq == 0 || seq > ackSnap) continue;
            const SnapSent& s = ds.ring[seq % kDeltaRing];
            if (s.seq != seq) continue;
            for (const PlayerStateQ& q : s.states) {
                auto& b = ds.base[q.session];
                if (seq > b.first) b = { seq, q };
            }
        }
    }

    void Send(uint32_t conn, const NetWriter& w) { Net_Send(conn, w.bytes.data(), (uint32_t)w.bytes.size(), true); }
    template <class M> void SendMsg(uint32_t conn, const M& m) { NetWriter w; NetWrite(w, m); Send(conn, w); }
    std::string BuildCommit() { const std::string v = kGameVersion; const size_t sp = v.rfind(' '); return sp == std::string::npos ? v : v.substr(sp + 1); }

    template <class M> void SendToJoined(const M& m) { NetWriter w; NetWrite(w, m); for (const auto& [c, p] : peers) if (p.joined) Send(c, w); }
    // The roster (#750): everyone online, to everyone joined, on each join and leave (reliable; the whole list each time:
    // a handful of names, and a client never has to merge).
    void SendRoster() {
        MsgRoster r;
        for (const Session& s : sessions) if (s.online) r.players.push_back({ s.id, s.name });
        SendToJoined(r);
    }

    void Drop(uint32_t conn, NetBye why) {
        auto it = peers.find(conn);
        if (it == peers.end()) return;
        for (Session& s : sessions) if (s.id == it->second.session && s.conn == conn) { s.online = false; s.conn = 0; }
        if (why != NetBye::Quit) SendMsg(conn, MsgBye{ why });
        printf("host: connection %u (session %u) %s\n", conn, it->second.session, NetByeText(why));
        if (eventHook && it->second.joined) { // (#345 audit)
            std::string name; for (const Session& s : sessions) if (s.id == it->second.session) name = s.name;
            eventHook("leave", it->second.session, name + " (" + NetByeText(why) + ")");
        }
        const bool wasJoined = it->second.joined;
        Interest_Forget(conn); // (#340) its accumulators
        deltas.erase(conn);    // (#769) its snapshot bases
        peers.erase(it);
        Net_Close(conn);
        if (wasJoined) SendRoster(); // the others see them go
    }
    void Refuse(uint32_t conn, NetReject why, const std::string& detail) {
        SendMsg(conn, MsgReject{ why, detail });
        printf("host: refused connection %u: %s%s%s\n", conn, NetRejectText(why), detail.empty() ? "" : ": ", detail.c_str());
        peers.erase(conn);
        Net_Close(conn); // (lingers: the reject goes out first)
    }

    void OnHello(Scene& world, uint32_t conn, const MsgHello& h) {
        if (h.protocol != kNetProtocol) { Refuse(conn, NetReject::Protocol, TextFormat("server %u, client %u", kNetProtocol, h.protocol)); return; }
        if (!h.build.empty() && h.build != BuildCommit()) { Refuse(conn, NetReject::Build, "server " + BuildCommit() + ", client " + h.build); return; }
        if (!password.empty() && h.password != password) { Refuse(conn, NetReject::Password, ""); return; }
        if (joinCheck) { // (#345) bans and the whitelist
            std::string why;
            if (!joinCheck(h.name, Net_RemoteAddress(conn), why)) {
                Refuse(conn, NetReject::Banned, why);
                if (eventHook) eventHook("refused", 0, h.name + " (" + why + ")");
                return;
            }
        }
        Session* s = nullptr;
        for (Session& x : sessions) if (x.name == h.name) s = &x;
        if (s && s->online) { Refuse(conn, NetReject::NameTaken, h.name); return; }
        const int online = (int)std::count_if(sessions.begin(), sessions.end(), [](const Session& x) { return x.online; });
        if (online >= maxPlayers) { Refuse(conn, NetReject::Full, TextFormat("%d players", maxPlayers)); return; }
        const bool reconnected = s != nullptr;
        if (!s) { sessions.push_back({ nextSession++, h.name, 0, false }); s = &sessions.back(); }
        s->conn = conn; s->online = true;
        Peer& p = peers[conn];
        p.session = s->id; p.joined = true; p.lastHeard = GetTime();
        const JoinInfo info = JoinInfo_FromWorld(world);
        MsgWelcome w;
        w.session = s->id; w.reconnected = reconnected; w.seed = info.seed; w.bodies = info.bodies; w.fingerprint = info.fingerprint;
        w.serverTick = (uint32_t)SimTick(); w.worldTime = WorldTime(); // #339: the client takes our clock
        w.worldgen = WorldGenParams_ToText(WorldGenParams_ForSave(WorldSavePath())); w.sandbox = gSandboxWorld; // #757: generate it as we did
        w.settings.difficulty = (uint8_t)std::clamp(worldDifficulty, 0, 2); w.settings.maxPlayers = (uint8_t)maxPlayers;
        w.motd = motd; // (#326)
        SendMsg(conn, w);
        const EventsShower sh = Events_GetShower();
        SendMsg(conn, MsgShower{ sh.start, sh.end, sh.id });
        if (!remotes.count(s->id)) SpawnRemote(world, s->id); // (a reconnect keeps its player where it was)
        else { // (#780) a reconnect keeps the player but its client counts ticks from 1 again: a fresh input stream
            Remote& rm = remotes[s->id];
            rm.inputs.clear(); rm.lastQueued = 0; rm.lastApplied = 0; rm.carried = 0; rm.budget = 0.0; rm.tick = SimTick();
        }
        // Who is here (everyone, the newcomer included), and what the others hold.
        SendRoster();
        for (const Session& x : sessions) if (x.online && x.id != s->id && x.loadout.session) SendMsg(conn, x.loadout);
        ++joins;
        printf("host: %s joined as session %u%s (connection %u, %d joins so far)\n", h.name.c_str(), s->id, reconnected ? " (reconnected)" : "", conn, joins);
        if (eventHook) eventHook("join", s->id, h.name + " from " + Net_RemoteAddress(conn)); // (#345 audit)
    }

    // ---- Terrain edits (#329) ----
    constexpr double kReachSlack = 1.5;                // m past the tool's reach (the player moved on while the request flew)
    constexpr double kEditBurst = 3.0;                 // strokes that may arrive together (a hitch, packet bunching)
    int editsApplied = 0;

    void SendInventoryDeltas(uint32_t conn, const Inventory& inv, const std::map<std::string, double>& before) {
        for (const auto& [item, was] : before) {
            const double d = inv.amount(item) - was;
            if (std::fabs(d) >= 0.005) SendMsg(conn, MsgInventoryDelta{ item, std::round(d * 100.0) / 100.0 });
        }
    }
    void RefuseEdit(uint32_t conn, uint32_t seq, NetEditRefusal why) {
        SendMsg(conn, MsgEditRefused{ seq, why });
        printf("host: edit %u refused: %s\n", seq, NetEditRefusalText(why));
    }
    // Validate a stroke as the tool would make it, from the server's own view of the player (where it stands, what it
    // holds), then apply it, pay or credit that player, and send it to everyone in this order.
    void OnEditRequest(Scene& world, uint32_t conn, uint32_t session, const MsgEditRequest& q) {
        Remote& rm = remotes[session];
        Planetoid* p = Planetoid::findNet(q.body);
        if (!p) { RefuseEdit(conn, q.seq, NetEditRefusal::Body); return; }
        const Vector3 hit = q.hit.toFloat(), centre = q.center.toFloat();
        // Reach: from the eye of the server's player (body-local), as far as the tool's raycast goes.
        Player::NetState ns;
        if (!rm.player.netState(world, ns) || ns.body != q.body) { RefuseEdit(conn, q.seq, NetEditRefusal::Reach); return; }
        const DVec3 up = ns.local * (1.0 / std::max(1e-9, ns.local.length()));
        const DVec3 eye = ns.local + up * (double)rm.player.eyeHeight;
        if ((q.hit - eye).length() * p->scale > TerrainToolSettings{}.reach + kReachSlack) { RefuseEdit(conn, q.seq, NetEditRefusal::Reach); return; }
        if (!p->field->editable(hit)) { RefuseEdit(conn, q.seq, NetEditRefusal::Mantle); return; }
        // The tool in the player's hands (its loadout; possession is #331's) and what it does to the aimed ground.
        std::string held;
        for (const Session& s : sessions) if (s.id == session) held = s.loadout.held;
        const ToolSpec tool = ToolOf(held);
        const Vector3 surfN = p->field->normal(hit, p->edits.get());
        const VoxelMaterial mat = p->materialAt(hit, surfN);
        const DigRule rule = DigRuleFor(tool, mat, TerrainTool_ItemAt(*p, hit, (int)mat));
        const float baseRadius = (rule.radius > 0.0f ? rule.radius : DigRuleFor(tool, VoxelMaterial::Regolith, nullptr).radius) / p->scale;
        const bool dig = q.kind == NetEditKind::Dig;
        if (dig && !rule.can) { RefuseEdit(conn, q.seq, NetEditRefusal::Tool); return; }
        // The stroke's shape: the tool's sphere (place: 1.3x, the soft blob), its centre where the tool puts it.
        const float maxRadius = baseRadius * (dig ? 1.0f : 1.3f) * 1.01f + 1e-4f;
        const float offset = Vector3Distance(hit, centre);
        if (q.radius > maxRadius || offset > q.radius * (dig ? 0.5f : 0.3f) + 0.05f) { RefuseEdit(conn, q.seq, NetEditRefusal::Shape); return; }
        // Rate: a token per cooldown (the tool's strokes a second, at least 2), up to a small burst.
        const double rate = std::max(rule.strokesPerSecond, 2.0f);
        if (rm.editTokensAt < 0.0) rm.editTokensAt = SimTime();
        rm.editTokens = std::min(kEditBurst, rm.editTokens + (SimTime() - rm.editTokensAt) * rate);
        rm.editTokensAt = SimTime();
        if (rm.editTokens < 0.999) { RefuseEdit(conn, q.seq, NetEditRefusal::Rate); return; }
        MsgTerrainEdit out;
        out.body = q.body; out.center = q.center; out.radius = q.radius; out.session = session; out.seq = q.seq;
        TerrainEdit e{};
        e.center = centre; e.radius = q.radius;
        std::map<std::string, double> before;
        if (dig) {
            const bool softOnly = tool.type != ToolType::Pickaxe;
            const DigYield got = TerrainTool_DigYield(*p, centre, q.radius, hit, softOnly, 1);
            for (const auto& [item, amt] : got.items) { before.emplace(item, rm.inventory.amount(item)); }
            if (got.gem) before.emplace(got.gem, rm.inventory.amount(got.gem));
            if (got.seed) before.emplace(got.seed, rm.inventory.amount(got.seed));
            for (const auto& [item, amt] : got.items) rm.inventory.addFractional(item, amt);
            if (got.gem) rm.inventory.add(got.gem, 1);
            if (got.seed) rm.inventory.add(got.seed, 1);
            e.carve = true; e.softOnly = softOnly;
            out.flags = 1 | (softOnly ? 8 : 0);
        } else {
            const int paint = TerrainTool_FillerMaterial(q.item);
            if (paint < 0) { RefuseEdit(conn, q.seq, NetEditRefusal::Shape); return; }
            before.emplace(q.item, rm.inventory.amount(q.item));
            if (!rm.inventory.removeAmount(q.item, 1.0)) { RefuseEdit(conn, q.seq, NetEditRefusal::Cost); return; }
            e.strength = 0.6f; e.paint = paint;
            out.strength = e.strength; out.paint = (int8_t)paint;
        }
        rm.editTokens -= 1.0;
        ApplyTerrainEdit(*p, e);
        SendToJoined(out); // everyone, the asker too, in the server's order
        SendInventoryDeltas(conn, rm.inventory, before);
        ++editsApplied;
    }

    // ---- Building (#330) ----
    constexpr double kPartRate = 4.0, kPartBurst = 4.0; // parts a second (a break takes a 0.45 s hold; placing is a click)
    int partsApplied = 0;

    void RefusePart(uint32_t conn, uint32_t seq, NetPartRefusal why) {
        SendMsg(conn, MsgPartRefused{ seq, why });
        printf("host: part request %u refused: %s\n", seq, NetPartRefusalText(why));
    }
    Part ToPart(const NetPart& n) { Part p{}; p.type = (PartType)n.type; p.origin = { n.x, n.y, n.z }; p.orient = n.orient; p.material = n.material; return p; }
    NetPart ToNet(const Part& p) { return NetPart{ (uint8_t)p.type, p.origin.x, p.origin.y, p.origin.z, p.orient, p.material }; }
    // Is a body-local point within the player's build reach (from the eye of the server's player)?
    bool InBuildReach(Scene& world, Remote& rm, uint64_t body, Vector3 local, double scale) {
        Player::NetState ns;
        if (!rm.player.netState(world, ns) || ns.body != body) return false;
        const DVec3 up = ns.local * (1.0 / std::max(1e-9, ns.local.length()));
        const DVec3 eye = ns.local + up * (double)rm.player.eyeHeight;
        return (DVec3(local) - eye).length() * scale <= BuildToolState{}.reach + kReachSlack + 1.0; // (+1 m: a part's far cells)
    }
    void OnPartRequest(Scene& world, uint32_t conn, uint32_t session, const MsgPartRequest& q) {
        Remote& rm = remotes[session];
        Planetoid* p = Planetoid::findNet(q.body);
        if (!p) { RefusePart(conn, q.seq, NetPartRefusal::Body); return; }
        if (rm.partTokensAt < 0.0) rm.partTokensAt = SimTime();
        rm.partTokens = std::min(kPartBurst, rm.partTokens + (SimTime() - rm.partTokensAt) * kPartRate);
        rm.partTokensAt = SimTime();
        if (rm.partTokens < 0.999) { RefusePart(conn, q.seq, NetPartRefusal::Rate); return; }
        std::map<std::string, double> before;
        if (q.op == NetPartOp::Remove) {
            Plot* plot = p->plots ? p->plots->byId(q.plot) : nullptr;
            if (!plot) { RefusePart(conn, q.seq, NetPartRefusal::NoPlot); return; }
            if (q.partIndex >= plot->parts.size() || plot->parts[q.partIndex].type == PartType::Count) { RefusePart(conn, q.seq, NetPartRefusal::NoPart); return; }
            const Part removed = plot->parts[q.partIndex];
            const Vector3 at = plot->toBody({ (removed.origin.x + 0.5f) * kCellSize, (removed.origin.y + 0.5f) * kCellSize, (removed.origin.z + 0.5f) * kCellSize });
            if (!InBuildReach(world, rm, q.body, at, p->scale)) { RefusePart(conn, q.seq, NetPartRefusal::Reach); return; }
            if (!plot->remove((int)q.partIndex)) { RefusePart(conn, q.seq, NetPartRefusal::NoPart); return; }
            // The refund (a meteorite breaks into ore, rolled here: the authority rolls, #316).
            std::vector<std::pair<std::string, int>> give;
            if (removed.type == PartType::Meteorite) {
                give.push_back({ "nickel_iron_ore", GameRng_Int(6, 10) });
                if (GameRng_Int(0, 99) < 35) give.push_back({ "platinum_ore", 2 });
                if (GameRng_Int(0, 99) < 12) give.push_back({ "iridium_ore", 1 });
            } else give.push_back({ PartItem(removed), PartCostCount(removed.type) });
            for (const auto& [item, n] : give) { before.emplace(item, rm.inventory.amount(item)); rm.inventory.add(item, n); }
            MsgPlotEdit out; out.session = session; out.seq = q.seq; out.op = NetPartOp::Remove; out.body = q.body; out.plot = q.plot; out.partIndex = q.partIndex;
            SendToJoined(out);
        } else {
            const Part part = ToPart(q.part);
            if (q.part.type >= (uint8_t)PartType::Count || part.type == PartType::Meteorite || (IsShape(part.type) && part.material >= BuildMaterialCount())) {
                RefusePart(conn, q.seq, NetPartRefusal::Kind); return;
            }
            Plot* plot = nullptr;
            Vector3 anchor{ q.anchor[0], q.anchor[1], q.anchor[2] };
            if (q.plot) {
                plot = p->plots ? p->plots->byId(q.plot) : nullptr;
                if (!plot) { RefusePart(conn, q.seq, NetPartRefusal::NoPlot); return; }
                if (!plot->canPlace(part)) { RefusePart(conn, q.seq, NetPartRefusal::Blocked); return; }
            }
            // Reach: the part's origin cell (in the plot, or from the new plot's anchor).
            const Vector3 cellCentre{ (part.origin.x + 0.5f) * kCellSize, (part.origin.y + 0.5f) * kCellSize, (part.origin.z + 0.5f) * kCellSize };
            const Vector3 at = plot ? plot->toBody(cellCentre) : anchor; // (a new plot: the part sits on its anchor)
            if (!InBuildReach(world, rm, q.body, at, p->scale)) { RefusePart(conn, q.seq, NetPartRefusal::Reach); return; }
            const char* item = PartItem(part);
            const int cost = PartCostCount(part.type);
            if (rm.inventory.count(item) < cost) { RefusePart(conn, q.seq, NetPartRefusal::Cost); return; }
            if (!plot) { // a new plot, made here first (the clients make the same one from the same id and anchor)
                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                plot = p->plots->plotFor(anchor, 0.0f, true);
                MsgPlotEdit made; made.session = session; made.seq = q.seq; made.op = NetPartOp::NewPlot; made.body = q.body; made.plot = plot->id;
                made.anchor[0] = anchor.x; made.anchor[1] = anchor.y; made.anchor[2] = anchor.z;
                SendToJoined(made);
            }
            before.emplace(item, rm.inventory.amount(item));
            rm.inventory.remove(item, cost);
            BuildTool_ApplyPart(*p, *plot, part);
            MsgPlotEdit out; out.session = session; out.seq = q.seq; out.op = NetPartOp::Place; out.body = q.body; out.plot = plot->id; out.part = q.part;
            SendToJoined(out);
        }
        rm.partTokens -= 1.0;
        SendInventoryDeltas(conn, rm.inventory, before);
        ++partsApplied;
    }

    // (#349) the guard's verdicts: a kick (with why, in the log and the audit hook), and violations found after reading.
    int guardDropped = 0, guardKicks = 0;
    void GuardKick(uint32_t conn, const std::string& why) {
        auto it = peers.find(conn);
        std::string name; if (it != peers.end()) for (const Session& s : sessions) if (s.id == it->second.session) name = s.name;
        printf("host: kicked connection %u (%s) for flooding: %s\n", conn, name.empty() ? "not joined" : name.c_str(), why.c_str());
        if (eventHook) eventHook("kicked", it != peers.end() ? it->second.session : 0, (name.empty() ? std::string("connection") : name) + " (flooding: " + why + ")");
        if (it != peers.end() && it->second.joined) SendMsg(conn, MsgChat{ 0, "You were kicked: too many invalid messages" });
        ++guardKicks;
        Drop(conn, NetBye::Kicked);
    }
    void GuardViolation(uint32_t conn, const char* why) {
        auto it = peers.find(conn);
        if (it == peers.end()) return;
        ++guardDropped;
        if (NetGuard_Violation(it->second.guard, GetTime()) == GuardVerdict::Kick) GuardKick(conn, why);
    }
    void OnMessage(Scene& world, const NetMessage& m) {
        if (!peers.count(m.conn) && !Net_IsConnected(m.conn)) return; // (#349: the rest of a batch from a connection just kicked)
        auto it = peers.find(m.conn);
        if (it == peers.end()) it = peers.emplace(m.conn, Peer{}).first; // first word on a new connection
        it->second.lastHeard = GetTime();
        { // (#349) size, kind, form and rate first: a bad message is dropped; a client that keeps sending them is kicked
            std::string why;
            const GuardVerdict v = NetGuard_Admit(it->second.guard, NetPeek(m.data), m.data.size(), NetWellFormed(m.data), GetTime(), why);
            if (v != GuardVerdict::Accept) {
                ++guardDropped;
                if (v == GuardVerdict::Kick) GuardKick(m.conn, why);
                return;
            }
        }
        NetReader r(m.data.data() + (m.data.empty() ? 0 : 1), m.data.empty() ? 0 : m.data.size() - 1);
        switch (NetPeek(m.data)) {
            case NetMsg::Hello: {
                MsgHello h;
                if (it->second.joined) return; // already in
                if (!NetRead(r, h)) { Refuse(m.conn, NetReject::Malformed, ""); return; }
                OnHello(world, m.conn, h);
                return;
            }
            case NetMsg::Ping: {
                MsgPing p;
                if (!it->second.joined || !NetRead(r, p)) return;
                SendMsg(m.conn, MsgPong{ p.seq, p.clientTime, (uint32_t)SimTick(), WorldTime() }); // #339: and follows it
                return;
            }
            case NetMsg::Input: {
                MsgInputs in;
                if (!it->second.joined || !NetRead(r, in)) return;
                OnSnapshotAck(m.conn, in.ackSnap, in.ackMask); // (#769)
                Remote& rm = remotes[it->second.session];
                if (in.steps[0].tick > rm.lastQueued) { rm.action = in.action; rm.actionPhase = in.actionPhase; } // (not an old packet)
                // (#349) only real actions, a unit view direction (else the last one): nothing a client sends turns into NaNs
                // or a key the game doesn't have; a jump of the tick counter far ahead (to freeze its own later inputs, or a
                // bug) is refused.
                constexpr uint64_t kActMask = ((int)Act::Count < 64 ? (1ull << (int)Act::Count) - 1 : ~0ull);
                if (in.steps[0].tick > rm.lastQueued + 3600) { GuardViolation(m.conn, "input tick jumped ahead"); return; }
                for (auto s = in.steps.rbegin(); s != in.steps.rend(); ++s) // oldest first; the ones already here are repeats
                    if (s->tick > rm.lastQueued) {
                        MsgInput q = *s;
                        q.held &= kActMask; q.pressed &= kActMask;
                        const float len = sqrtf(q.fx * q.fx + q.fy * q.fy + q.fz * q.fz);
                        if (!(len > 0.5f && len < 1.5f)) { q.fx = rm.forwardLocal.x; q.fy = rm.forwardLocal.y; q.fz = rm.forwardLocal.z; }
                        else { q.fx /= len; q.fy /= len; q.fz /= len; }
                        rm.inputs.push_back(q); rm.lastQueued = q.tick;
                    }
                return;
            }
            case NetMsg::Chat: { // (#326) a player's line: the server's handler (commands, logging), else to everyone
                MsgChat c;
                if (!it->second.joined || !NetRead(r, c) || c.text.empty()) return;
                c.from = it->second.session; // (a client speaks only for itself)
                std::string name; for (const Session& s : sessions) if (s.id == c.from) name = s.name;
                if (chatHandler) chatHandler(c.from, name, c.text);
                else SendToJoined(c);
                return;
            }
            case NetMsg::Loadout: {
                MsgLoadout l;
                if (!it->second.joined || !NetRead(r, l)) return;
                l.session = it->second.session; // (a client speaks only for itself)
                for (Session& s : sessions) if (s.id == l.session) s.loadout = l;
                SendToJoined(l);
                return;
            }
            case NetMsg::EditRequest: {
                MsgEditRequest q;
                if (!it->second.joined || !NetRead(r, q)) return;
                OnEditRequest(world, m.conn, it->second.session, q);
                return;
            }
            case NetMsg::PartRequest: {
                MsgPartRequest q;
                if (!it->second.joined || !NetRead(r, q)) return;
                OnPartRequest(world, m.conn, it->second.session, q);
                return;
            }
            case NetMsg::Bye: { MsgBye b; if (NetRead(r, b)) Drop(m.conn, NetBye::Quit); return; }
            default: return; // not joined, or not a message for the server yet: dropped
        }
    }
}

bool Host_Start(uint16_t port, std::string& status) {
    gNetRole = NetRole::Host; // the authority either way: the world runs and rolls here
    lastSave = GetTime();
    if (!Net_Available()) { status = "built without ENABLE_GNS: the world runs, nobody can connect"; return false; }
    std::string err;
    if (!Net_Init(err)) { status = "network init failed: " + err; return false; }
    if (!Net_Host(port)) { status = "could not listen on port " + std::to_string(port); return false; }
    listening = true;
    status = "listening on port " + std::to_string(port);
    return true;
}

void Host_SetAutosaveMinutes(double minutes) { autosaveSeconds = minutes > 0.0 ? minutes * 60.0 : 0.0; }
void Host_SetPassword(const std::string& p) { password = p; }
void Host_SetMaxPlayers(int n) { maxPlayers = std::clamp(n, 1, 64); }

void Host_Update(Scene& world) {
    // Autosave (user, 2026-11-15: every 30 minutes by default; --autosave-minutes N, 0 = off). Wall-clock time: a
    // server's saves protect against crashes, whatever the game clock does.
    if (autosaveSeconds > 0.0 && lastSave >= 0.0 && GetTime() - lastSave >= autosaveSeconds) {
        lastSave = GetTime();
        const bool ok = SaveWorld(world);
        if (ok) ++autosaves;
        printf("host: autosave %d to %s %s\n", autosaves, WorldSavePath(), ok ? "done" : "FAILED");
    }
    if (!listening) return;
    Net_Poll();
    std::vector<NetMessage> in;
    Net_Receive(in, 256);
    for (NetMessage& m : in) {
        if (localClients.count(m.conn)) { localInbox[m.conn].push_back(std::move(m)); continue; }
        OnMessage(world, m);
    }
    // Gone: closed under us (the transport noticed), or silent past the timeout.
    std::vector<uint32_t> gone, silent;
    for (const auto& [conn, p] : peers) {
        if (!Net_IsConnected(conn)) gone.push_back(conn);
        else if (GetTime() - p.lastHeard > kTimeoutSeconds) silent.push_back(conn);
    }
    for (uint32_t c : gone) Drop(c, NetBye::Quit);
    for (uint32_t c : silent) Drop(c, NetBye::TimedOut);
    GpuTimer_CpuMark("host: receive, handshakes"); // (#348)
    StepRemotes(world);
    GpuTimer_CpuMark("host: remote players");
    SendSnapshots(world);
    GpuTimer_CpuMark("host: snapshots");
}

bool Host_SaveOnShutdown(Scene& world) {
    const bool ok = SaveWorld(world);
    printf("host: saved on shutdown to %s %s\n", WorldSavePath(), ok ? "done" : "FAILED");
    return ok;
}

void Host_Stop() {
    if (listening) {
        // Everyone hears why (#321): the server closed (each close lingers until the bye is out).
        std::vector<uint32_t> conns;
        for (const auto& [conn, p] : peers) conns.push_back(conn);
        for (uint32_t c : conns) Drop(c, NetBye::ServerClosed);
        Net_Shutdown();
    }
    listening = false;
    localClients.clear();
    localInbox.clear();
    peers.clear();
    sessions.clear();
    nextSession = 1;
    remotes.clear();
    gNetRole = NetRole::Offline;
}

bool Host_Listening() { return listening; }
int Host_Joins() { return joins; }
int Host_Autosaves() { return autosaves; }
int Host_OnlinePlayers() { return (int)std::count_if(sessions.begin(), sessions.end(), [](const Session& s) { return s.online; }); }
void Host_SetLocalClient(uint32_t conn) { localClients.insert(conn); }
std::vector<NetMessage> Host_TakeLocalClientMessages(uint32_t conn) {
    std::vector<NetMessage> out;
    auto it = localInbox.find(conn);
    if (it != localInbox.end()) out.swap(it->second);
    return out;
}
bool Host_PlayerLocal(Scene& world, uint32_t session, uint64_t& body, DVec3& local) {
    auto it = remotes.find(session);
    if (it == remotes.end()) return false;
    const MsgPlayerState s = StateOf(world, session, it->second);
    body = s.body; local = s.local;
    return true;
}
int Host_EditsApplied() { return editsApplied; }
double Host_InventoryAmount(uint32_t session, const std::string& item) {
    auto it = remotes.find(session);
    return it == remotes.end() ? 0.0 : it->second.inventory.amount(item);
}
int Host_PartsApplied() { return partsApplied; }
bool Host_Give(uint32_t session, const std::string& item, int count) {
    auto it = remotes.find(session);
    if (it == remotes.end() || count <= 0 || !itemDB.find(item)) return false;
    std::map<std::string, double> before{ { item, it->second.inventory.amount(item) } };
    it->second.inventory.add(item, count);
    for (const Session& s : sessions) if (s.id == session && s.online) SendInventoryDeltas(s.conn, it->second.inventory, before);
    return true;
}

// ---- The dedicated server's console and admin (#326, #345) ----
std::vector<HostPlayer> Host_ListPlayers() {
    std::vector<HostPlayer> out;
    for (const Session& s : sessions) out.push_back({ s.id, s.name, s.online, s.online ? Net_RemoteAddress(s.conn) : std::string() });
    return out;
}
bool Host_Kick(const std::string& name, const std::string& reason) {
    for (const Session& s : sessions)
        if (s.online && s.conn && s.name == name) {
            if (!reason.empty()) SendMsg(s.conn, MsgChat{ 0, "You were kicked: " + reason });
            Drop(s.conn, NetBye::Kicked);
            return true;
        }
    return false;
}
void Host_Say(const std::string& text, uint32_t fromSession, uint32_t toSession) {
    const MsgChat c{ fromSession, text };
    if (toSession == 0) { SendToJoined(c); return; }
    for (const Session& s : sessions) if (s.id == toSession && s.online && s.conn) SendMsg(s.conn, c);
}
void Host_SetMotd(const std::string& text) { motd = text; }
void Host_SetJoinCheck(std::function<bool(const std::string&, const std::string&, std::string&)> check) { joinCheck = std::move(check); }
void Host_SetChatHandler(std::function<void(uint32_t, const std::string&, const std::string&)> handler) { chatHandler = std::move(handler); }
void Host_SetEventHook(std::function<void(const char*, uint32_t, const std::string&)> hook) { eventHook = std::move(hook); }
// Tests (#349): messages the guard dropped, and connections it kicked.
int Host_GuardDropped() { return guardDropped; }
int Host_GuardKicks() { return guardKicks; }
