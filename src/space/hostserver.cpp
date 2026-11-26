#include "hostserver.h"
#include "headless.h"
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
#include "machineact.h" // (#1684) the dungeon action codes
#include "plots.h"
#include "factory.h"
#include "gamerng.h"
#include "tools.h"
#include "animals.h"
#include "footsteps.h" // (#1224) players' step noise for creature hearing
#include "suit.h"
#include "combat.h"
#include "teamquests.h"
#include "navigation.h"
#include "recipebook.h"
#include "craftplan.h"
#include "skills.h"
#include "playerrecord.h"
#include "inventory.h"
#include "farming.h"
#include "claims.h"
#include "serveradmin.h"
#include "trees.h"
#include "structures.h" // (#1641 slice 2c) breakable props
#include <cmath>
#include <map>
#include <ctime>
#include <future>
#include <chrono>
#include <filesystem>
#include "playerinput.h"
#include "worldtypes.h"
#include "voxel.h"
#include <raymath.h>
#include <deque>

#include <raylib.h>
#include <algorithm>
#include <functional>
#include "gamesave.h"
#include "playerrecord.h"
#include <cstdio>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>

extern const char* kGameVersion; // "<branch> <commit>" (cmake/game_version.cmake, #501)
#include "utils/logger.h"
#include <cstdarg>

namespace {
    // (#1916) A connection's life (connect, refusal, join, join sync, kick, disconnect): on stdout as "host: ..." as before
    // (the console, Pterodactyl, tools/test_server.sh) and in log.txt under [Host], so a join can be read back afterwards.
#if defined(__GNUC__)
    __attribute__((format(gnu_printf, 1, 2)))
#endif
    void HostLog(const char* fmt, ...) {
        char b[1024];
        va_list a; va_start(a, fmt); vsnprintf(b, sizeof b, fmt, a); va_end(a);
        printf("host: %s\n", b);
        LOG_INFO("Host", "%s", b);
    }
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
    struct Peer { uint32_t session = 0; bool joined = false; double lastHeard = 0.0, firstHeard = 0.0; NetGuard guard; // guard, firstHeard: #349
                  bool inGame = false; }; // (#1033) its first input came: the client is past its world load
    std::unordered_map<uint32_t, Peer> peers;          // by connection
    struct Session {
        uint32_t id = 0; std::string name; uint32_t conn = 0; bool online = false;
        std::string uuid;            // (#346) the player's id: their record; a reconnect is by it, not the name
        MsgLoadout loadout;          // what the player holds and wears (#750), kept across a reconnect
    };
    std::vector<Session> sessions;
    // #326 / #345: the dedicated server's hooks (main.cpp sets them in --server mode).
    std::string motd;
    std::function<bool(const std::string& name, const std::string& address, const std::string& uuid, std::string& why)> joinCheck; // (#1868: and the id)
    std::function<void(uint32_t session, const std::string& name, const std::string& text)> chatHandler;
    std::function<void(const char* what, uint32_t session, const std::string& detail)> eventHook;
    uint32_t nextSession = 1;
    constexpr double kTimeoutSeconds = 15.0;           // nothing heard (pings come every second): gone
    // (#1033) A client that joined but hasn't sent an input yet is loading the world: one long frame there (the spawn load
    // took ~37 s on the runner) sends no pings. It gets longer before it counts as gone (a closed connection is still
    // noticed at once by the transport). 120 s: the user's choice (2026-11-20), room for slow machines.
    constexpr double kLoadingTimeoutSeconds = 120.0;

    // Players (#327): the server's own simulation of every session's player, from the inputs its client sends (one a sim
    // step; none waiting: the last keys stay held). Kept across a reconnect: the player is where they left.
    struct Remote {
        Player player;
        // (#844) The player's body-local state after its last step (or its spawn): a player is stepped only by its inputs,
        // and in a gap (a client loading the world, a stall) its world-frame velocity doesn't turn with the planet, so the
        // first step after it landed on ground that had moved on (a 36 m/s "fall" at --time-scale 60). Re-applied before a
        // step after a gap, so the player rides its body through it. restSpawns: a respawn or teleport makes it stale.
        Player::NetState rest{};
        bool restValid = false;
        int restSpawns = 0;
        long long restTick = 0;
        std::deque<MsgInput> inputs;
        uint32_t lastApplied = 0;   // the client tick of the last input stepped (the owner replays the ones after)
        uint32_t lastQueued = 0;    // the newest input tick received (redundant copies of older ones are skipped)
        uint64_t lastHeld = 0;      // (#334) the held actions of the last input stepped
        uint64_t carried = 0;       // presses of inputs dropped from a full queue, for the next step
        double budget = 0.0;        // steps it may take (earned at real time, capped: a fast client can't run ahead)
        long long tick = 0;         // the server sim tick the budget has been earned to
        Vector3 forwardLocal{ 0, 0, 1 };
        uint8_t action = 0; float actionPhase = 0.0f; // the hands, from the newest input packet (relayed, #750)
        // Terrain edits (#329): the player's inventory as the server keeps it (yields go in, placing is paid from it; the
        // client mirrors it from InventoryDelta), and a stroke budget (a token per tool cooldown, a short burst allowed).
        Inventory inventory;
        // #337: its survival, simulated here with each input stepped (oxygen, power, food, temperature, falls, hits, death and
        // respawn), and its suit's numbers (the starting suit until suits are replicated, #331).
        SurvivalState survival;
        SuitStats suitStats = Suit_StatsOf(SuitState{});
        bool deathHandled = false;   // (#804) this death's drops are done (a grave, losses, or none)
        bool inNight = false; int deathsAtDusk = 0; // (#338) a night survived (dark fell and day came back without a death)
        std::unordered_set<uint64_t> tilesAlone; // (#338 tests) the map tiles this player revealed itself, this session (logged on leaving)
        double editTokens = 2.0, editTokensAt = -1.0;
        double partTokens = 4.0, partTokensAt = -1.0; // building (#330): placing and breaking, a few a second
        PlayerRecordMeta record;    // (#346) their record's id, name, times, and the lines the server keeps as text
        bool creative = false;        // (#823) /gamemode creative (ops): building, digging and placing cost nothing; it may fly
        SkillState skills;           // (#920) XP and perks, counted here from what the server sees (strokes); saved in the record
        bool skillsDirty = true;     // (#920) to send (SkillState): at the join, after XP changes (at most every 2 s; a level-up at once)
        double skillsSentAt = -1e9;
        double joinedAt = 0.0;      // (play time, added at the next write)
        double lastFellAt = -1e9;   // (#335) its last plant felled (the fell rate)
        uint64_t farmHash = 0;      // (#335) the farm fields last sent to it
    };
    std::unordered_map<uint32_t, Remote> remotes;      // by session
    // (#1868) The item a player really holds: its loadout's, when the inventory the server keeps for it has one; else hands ("").
    // Read where it counts (digging, felling, hitting), so an item that left the inventory after the loadout came stops too.
    // Slots are matched by id (not Inventory::count: tools that carry item data are left out of counts). Creative players may
    // take any item from the Catalog, so their loadout stands as sent.
    bool InventoryHolds(const Inventory& inv, const std::string& id) {
        for (const InvSlot& s : inv.grid()) if (!s.empty() && s.id == id) return true;
        return false;
    }
    std::string HeldItem(uint32_t session) {
        std::string held;
        for (const Session& s : sessions) if (s.id == session) held = s.loadout.held;
        if (held.empty()) return held;
        const auto rit = remotes.find(session);
        if (rit == remotes.end()) return "";
        return rit->second.creative || InventoryHolds(rit->second.inventory, held) ? held : std::string();
    }
    int loadoutsRefused = 0; // (#1868) loadouts naming an item the player doesn't have (tests)
    void RemoteXp(uint32_t session, const char* event, double n = 1.0); // (#920) XP for a remote player's event (below)
    void TeamXp(uint32_t session, const char* event, double n);         // (#920) and teammates' share of it (below)
    // (#920) A remote player's perk multipliers by stat name (Skills_Stat on the skills the server keeps for them).
    std::function<float(const char*)> StatsOf(const Remote& rm) {
        const SkillState* s = &rm.skills;
        return [s](const char* name) { return (float)Skills_Stat(*s, name); };
    }
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
    constexpr long long kRestGapTicks = 6; // (#844) 0.1 s without a step counts as a gap (steps come 2 a frame at most)
    bool RestGapOff() { static const bool off = getenv("RAYLIBGAME_REST_GAP") && atoi(getenv("RAYLIBGAME_REST_GAP")) == 0; return off; } // A/B
    void StepRemotes(Scene& world) {
        for (auto& [session, r] : remotes) {
            // (#337) a dead player's client sends no inputs: its death clock (respawn after 3 s) runs here on the sim clock
            if (r.survival.dead) Survival_StepFor(world, SimDt(), PlayerInput::None(), r.survival, r.player, r.inventory, r.suitStats, false);
            r.budget = std::min(kMaxBudget, r.budget + (double)(SimTick() - r.tick));
            r.tick = SimTick();
            while (r.inputs.size() > kMaxQueuedInputs) { r.carried |= r.inputs.front().pressed; r.inputs.pop_front(); }
            // (#768) at most kMaxStepsPerFrame inputs a frame per player: after a slow frame the backlog is worked off over the
            // next frames, not all at once (all at once made the next frame slower still: the 31-40 ms catch-up bursts at 32 players).
            if (!r.survival.dead && (!r.restValid || r.restSpawns != r.player.spawns)) { // (#844) just spawned or moved: its state now
                RenderOriginScope frame(world, r.player.position);
                r.restValid = r.player.netState(world, r.rest);
                r.restSpawns = r.player.spawns; r.restTick = SimTick();
            }
            int stepsThisFrame = 0;
            while (r.budget >= 1.0 && !r.inputs.empty() && stepsThisFrame++ < kMaxStepsPerFrame) {
                r.budget -= 1.0;
                const MsgInput m = r.inputs.front();
                r.inputs.pop_front();
                PlayerInput in;
                in.enabled = true; in.held = m.held; in.pressed = m.pressed | r.carried;
                r.lastHeld = m.held; // (#334: a block held when an animal bites)
                r.carried = 0;
                r.lastApplied = m.tick;
                r.forwardLocal = { m.fx, m.fy, m.fz };
                if (Planetoid* b = r.player.nearest ? r.player.nearest : StartBody(world)) r.player.forward = b->dirToWorld(r.forwardLocal);
                RenderOriginScope frame(world, r.player.position);
                if (!r.survival.dead && r.restValid && r.restSpawns == r.player.spawns && SimTick() - r.restTick > kRestGapTicks && !RestGapOff())
                    r.player.applyNetState(world, r.rest); // (#844) back where its body has carried it, moving with it
                // (#1087) the player's own suit, as single player applies it every frame: the starting suit has no jetpack
                // (no thrust, an empty tank); without this a remote player kept the Player defaults, a full, working jetpack
                Suit_ApplyStatsTo(r.player, r.suitStats, (float)Skills_Stat(r.skills, "walk_speed"));
                r.survival.invulnerable = r.creative; // (#1561) creative: no damage, no death (this player's own mode)
                if (!r.survival.dead) { r.player.update(world, in, SimDt()); if (r.player.landingSpeed > 0.0f) Survival_OnLandingFor(r.survival, r.suitStats, r.player.landingSpeed, false); } // (dead: it lies still)
                if (!r.survival.dead) { r.restValid = r.player.netState(world, r.rest); r.restSpawns = r.player.spawns; r.restTick = SimTick(); }
                Survival_StepFor(world, SimDt(), in, r.survival, r.player, r.inventory, r.suitStats, false); // (#337)
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
    void SendInventoryDeltas(uint32_t conn, const Inventory& inv, const std::map<std::string, double>& before); // (below)
    template <class M> void SendMsg(uint32_t conn, const M& m); // (below)
    template <class M> void SendToJoined(const M& m); // (below)
    // ---- #334 animals ----
    // The animals near each client, 10 a second: on its player's body, within kAnimalRange, nearest first, at most
    // kAnimalsPerClient; and the ones it was sent last time that are gone (died, or out of its range).
    constexpr int kAnimalEvery = 6;          // sim steps (60 Hz): 10 a second
    constexpr float kAnimalRange = 200.0f;   // m (they draw to 160)
    constexpr size_t kAnimalsPerClient = 64;
    long long lastAnimals = 0;
    std::unordered_map<uint32_t, std::unordered_set<uint64_t>> animalsSent; // by connection
    std::unordered_map<uint32_t, double> lastAnimalHit;                     // by session (the cooldown)
    int animalHits = 0, animalHitsRefused = 0, animalHitsRewound = 0; // (rewound: #1108, allowed only by the rewind)
    bool RemoteLocal(const Remote& r, uint64_t& body, Vector3& local, double& scale) {
        Planetoid* b = r.player.nearest;
        if (!b) return false;
        body = b->netId; local = b->toLocalD(r.player.position).toFloat(); scale = b->scale;
        return true;
    }
    // (#1343, #1364) Creature shots launched and ended, and creature events, to each player near them (80 m), each tick
    // there are any: the clients fly the shots and play the creatures' sounds and captions.
    void SendCreatureFx(Scene&) {
        static MsgCreatureFx all;
        Animals_TakeFx(all);
        if (all.shots.empty() && all.ended.empty() && all.events.empty()) return;
        constexpr float kFxRange = 80.0f;
        for (const auto& [conn, p] : peers) {
            if (!p.joined) continue;
            auto rit = remotes.find(p.session);
            uint64_t body = 0; Vector3 at{}; double scale = 1.0;
            if (rit == remotes.end() || !RemoteLocal(rit->second, body, at, scale)) continue;
            MsgCreatureFx m;
            for (const MsgCreatureFx::Shot& s : all.shots)
                if (s.body == body && Vector3Distance(s.pos.toFloat(), at) * (float)scale <= kFxRange) m.shots.push_back(s);
            m.ended = all.ended; // (ids it never had are ignored there)
            for (const MsgCreatureFx::Event& e : all.events)
                if (e.body == body && Vector3Distance(e.pos.toFloat(), at) * (float)scale <= kFxRange) m.events.push_back(e);
            if (m.shots.empty() && m.ended.empty() && m.events.empty()) continue;
            NetWriter w; NetWrite(w, m);
            Net_Send(conn, w.bytes.data(), (uint32_t)w.bytes.size(), false);
        }
        all.shots.clear(); all.ended.clear(); all.events.clear();
    }
    void SendAnimals(Scene& world) {
        if (SimTick() - lastAnimals < kAnimalEvery) return;
        lastAnimals = SimTick();
        static std::vector<AnimalNet> all; Animals_CollectNet(world, all);
        for (const auto& [conn, p] : peers) {
            if (!p.joined) continue;
            auto rit = remotes.find(p.session);
            uint64_t body = 0; Vector3 at{}; double scale = 1.0;
            if (rit == remotes.end() || !RemoteLocal(rit->second, body, at, scale)) continue;
            std::vector<std::pair<float, size_t>> near;
            for (size_t i = 0; i < all.size(); ++i) {
                if (all[i].body != body) continue;
                const float d = Vector3Distance(all[i].pos, at) * (float)scale;
                if (d <= kAnimalRange) near.push_back({ d, i });
            }
            std::sort(near.begin(), near.end());
            if (near.size() > kAnimalsPerClient) near.resize(kAnimalsPerClient);
            MsgAnimals m; m.tick = (uint32_t)SimTick();
            std::unordered_set<uint64_t> now;
            for (const auto& [d, i] : near) {
                const AnimalNet& a = all[i];
                m.entries.push_back({ a.netId, a.species, a.body, DVec3(a.pos), a.heading, a.state, a.health, a.age, a.flags });
                now.insert(a.netId);
            }
            std::unordered_set<uint64_t>& was = animalsSent[conn];
            for (uint64_t id : was) if (!now.count(id)) m.removed.push_back(id);
            was.swap(now);
            if (m.entries.empty() && m.removed.empty()) continue;
            NetWriter w; NetWrite(w, m);
            Net_Send(conn, w.bytes.data(), (uint32_t)w.bytes.size(), false);
        }
    }
    // Animals live where players are (#318's interest centres: herds spawn and are simulated around each).
    void SetInterestCentres(Scene&) {
        gInterestCentres.clear();
        for (const Session& s : sessions) {
            if (!s.online) continue;
            auto rit = remotes.find(s.id);
            if (rit == remotes.end() || !rit->second.player.nearest) continue;
            const Planetoid* b = rit->second.player.nearest;
            gInterestCentres.push_back({ (uint64_t)b->id, b->toLocalD(rit->second.player.position).toFloat() });
        }
    }
    // A client's swing at an animal: the animal must be on its player's body and within reach of the server's player (+ the
    // distance it may have moved since the client saw it), and the player's swing cooldown spent. The damage is the held
    // tool's (as in single player); the drops go to that player's inventory.
    constexpr float kAnimalReach = 4.5f; // m: the 3 m aim, the animal's size, and its motion during the round trip
    constexpr double kHitRewind = 0.2;   // (#1108) s: how far back a hit may be checked (the combat design's cap, #1002)
    constexpr double kAnimalHitCooldown = 0.25; // s
    // (#337) each player's survival numbers to that player, 5 a second (its HUD; its own simulation fills in between)
    constexpr int kVitalsEvery = 12; // sim steps
    long long lastVitals = 0;
    int itemsUsed = 0;
    int gravesMade = 0; // (#804)
    int crafts = 0; // (#338) recipes crafted for players
    int dropsMade = 0, creativeGives = 0, creativeRefused = 0; // (#865, #866)
    // ---- Graves (#804; the user's rules on #337, 2026-11-18) ----
    // A remote player's death leaves what death_drops says (by difficulty, keep, or a grave owned by that player); walking up
    // to a grave takes what that player may (Combat_CanLoot); every client gets the list when it changes; the owners' teams
    // online keep their graves from decaying (14 days away and a grave goes, with what it holds).
    DeathDrops deathDrops = DeathDrops::Difficulty;
    uint32_t gravesSent = 0;
    long long lastGravePickup = 0;
    double lastGraveDecay = -1e9;
    std::map<std::string, double> AllAmounts(Inventory& inv) {
        std::map<std::string, double> m;
        for (const InvSlot& s : inv.grid()) if (!s.empty()) m[s.id] = inv.amount(s.id);
        return m;
    }
    void SendChangedItems(uint32_t conn, Inventory& inv, std::map<std::string, double> before) {
        for (const InvSlot& s : inv.grid()) if (!s.empty() && !before.count(s.id)) before[s.id] = 0.0; // new ones (picked up)
        SendInventoryDeltas(conn, inv, before);
    }
    void SendGraves(uint32_t conn) {
        MsgGraves m;
        for (const GraveInfo& g : Combat_GraveList()) { MsgGraves::Entry e{ g.id, g.body, g.local, g.owner, g.name, g.style, g.diedAt }; e.dropped = g.dropped; e.items = g.items; m.graves.push_back(std::move(e)); } // (#865: dropped piles too)
        SendMsg(conn, m);
    }
    void StepGraves(Scene& world) {
        const bool pickupTick = SimTick() - lastGravePickup >= 6; // ten times a second
        if (pickupTick) lastGravePickup = SimTick();
        std::vector<std::string> online;
        for (const auto& [conn, p] : peers) {
            if (!p.joined) continue;
            auto rit = remotes.find(p.session);
            if (rit == remotes.end()) continue;
            Remote& r = rit->second;
            std::string name;
            for (const Session& s : sessions) if (s.id == p.session) name = s.name;
            online.push_back(name);
            if (r.survival.dead && !r.deathHandled) {
                r.deathHandled = true;
                // (#844) what killed them, for the log (a death on a peaceful world had no cause to go on)
                printf("host: %s died (death %d): %s%s; food %.0f, oxygen %.0f, %s\n", name.c_str(), r.survival.deaths, r.survival.lastHarm[0] ? r.survival.lastHarm : "unknown",
                       !strcmp(r.survival.lastHarm, "fall") ? TextFormat(" at %.1f m/s", r.survival.lastFallSpeed) : "", r.survival.food, r.survival.oxygen,
                       r.player.flyMode ? "flying" : r.player.grounded ? "on the ground" : "in the air");
                const auto before = AllAmounts(r.inventory);
                if (Combat_OnDeathFor(r.inventory, nullptr, r.player, name, deathDrops, false)) ++gravesMade;
                SendChangedItems(conn, r.inventory, before);
            } else if (!r.survival.dead) {
                r.deathHandled = false;
                if (pickupTick) {
                    const auto before = AllAmounts(r.inventory);
                    Player::NetState feet; // (body-local feet as its last step left them: its world position goes stale as the planet turns, #856)
                    const Planetoid* body = r.player.netState(world, feet) ? Planetoid::findNet(feet.body) : nullptr;
                    if (body && Combat_PickUpGravesFor(body->name, feet.local, r.inventory, name, false, false)) SendChangedItems(conn, r.inventory, before);
                }
            }
        }
        if (GetTime() - lastGraveDecay > 10.0) { lastGraveDecay = GetTime(); Combat_UpdateGraveOwners(online, (int64_t)time(nullptr)); }
        if (Combat_GraveVersion() != gravesSent) {
            gravesSent = Combat_GraveVersion();
            for (const auto& [conn, p] : peers) if (p.joined) SendGraves(conn);
        }
    }

    // ---- Team quests (#338, #816; the user's rules on #303: quests and tiers per team, rewards to the completing player) ----
    // The server counts what it decides into the player's team's book (items obtained, machines placed, eating, animals
    // slain and fed), evaluates each team's book once a second while someone of it is online, pays a finished node's
    // reward to the member whose event came last (Host_Give: deltas), tells the team, and sends the book to its members
    // when it changes (and on join).
    std::map<int, uint32_t> questContributor;   // team -> the session whose event was counted last
    std::map<int, std::string> questBookSent;    // team -> the book text its members have
    long long lastQuestTick = 0;
    int questRewards = 0;
    int TeamOfSession(uint32_t session) {
        for (const Session& s : sessions) if (s.id == session) return Teams_Of(s.uuid);
        return 0;
    }
    void TeamEvent(uint32_t session, const char* kind, const std::string& what, double n = 1.0) {
        const int team = TeamOfSession(session);
        if (!team) return;
        if (!strcmp(kind, "count")) TeamQuests_Count(team, what, n);
        else if (!strcmp(kind, "placed")) TeamQuests_NotePlaced(team, what);
        else if (!strcmp(kind, "obtained")) TeamQuests_Obtained(team, what, n);
        questContributor[team] = session;
    }
    void SendQuestBook(uint32_t conn, int team) {
        MsgQuestBook m; m.team = (uint32_t)team; m.text = TeamQuests_Text(team);
        SendMsg(conn, m);
    }
    // ---- Machines (#332): the factory runs here only; players act on machines by request and see views ----
    // Each player's open machine (watch) gets its full view 5 times a second; the machines working on the plots within
    // kMachineRange of each player go out (as part lists) twice a second when they change. Results of a player's action
    // (extract, take a slot, leftovers) reach its inventory through gFactoryDeliver, then its client by deltas.
    constexpr float kMachineRange = 150.0f;      // m: plots whose working machines a player hears of
    constexpr float kMachineReach = 8.0f;        // m: how far a player reaches to act on a machine
    long long lastMachineViews = 0, lastMachineStates = 0;
    int machineActions = 0, machineActionsRefused = 0;
    struct MachineWatch { Factory::MachineKey key{ 0, 0, -1 }; bool on = false; };
    std::map<uint32_t, MachineWatch> machineWatch;                                    // session -> its open machine
    std::map<std::pair<uint32_t, std::pair<uint64_t, uint32_t>>, std::vector<int>> machineStatesSent; // (session, (body, plot)) -> parts
    // (#854) the chests near each player: their views go to it (what crafting may draw on), once a second when they change
    constexpr float kChestRange = 10.0f;
    long long lastChestViews = 0;
    std::map<std::pair<uint32_t, Factory::MachineKey>, uint64_t> chestViewSent; // (session, chest) -> hash of what it was sent
    bool ClaimBlocks(uint32_t session, const Planetoid* b, Vector3 local, bool count = true); // (#343, below)
    std::string ClaimRefusal(const Planetoid* b, Vector3 local);
    std::vector<std::pair<Planetoid*, std::pair<Plot*, int>>> ChestsNear(uint32_t session, const Remote& r) { // (#343: not another team's claimed chests)
        std::vector<std::pair<Planetoid*, std::pair<Plot*, int>>> out;
        Planetoid* b = r.player.nearest;
        if (!b || !b->plots) return out;
        const Vector3 me = b->toLocalD(r.player.position).toFloat();
        for (auto& pl : b->plots->plots)
            for (int k = 0; k < (int)pl->parts.size(); ++k) {
                if (pl->parts[k].type != PartType::Chest) continue;
                const Int3 c = pl->parts[k].origin;
                if (Vector3Distance(pl->toBody({ (c.x + 1.0f) * kCellSize, (c.y + 1.0f) * kCellSize, (c.z + 1.0f) * kCellSize }), me) <= kChestRange && !ClaimBlocks(session, b, pl->toBody({ (c.x + 1.0f) * kCellSize, (c.y + 1.0f) * kCellSize, (c.z + 1.0f) * kCellSize }), false)) out.push_back({ b, { pl.get(), k } });
            }
        return out;
    }
    uint32_t ConnOf(uint32_t session); // (below)
    // (#854) crafts waiting for what they pulled from chests to arrive (Factory::Take delivers next tick)
    struct PendingCraft { uint32_t session = 0; int recipe = 0, count = 0; double until = 0.0; };
    std::vector<PendingCraft> pendingCrafts;
    void StepMachines(Scene& world) {
        // The factory ticks the plots near every player at full rate and publishes what they see.
        std::vector<Vector3> near;
        for (auto& [sid, r] : remotes) near.push_back((r.player.position - world.renderOrigin).toFloat());
        Factory::SetNearPoints(near);
        Factory::SetCameraNear(!Headless()); // (#1302) a dedicated server's camera is nobody: only its players make plots near
        std::vector<Factory::MachineKey> keys;
        for (const auto& [sid, w] : machineWatch) if (w.on) keys.push_back(w.key);
        for (const auto& [sid, r] : remotes) for (const auto& [b, pc] : ChestsNear(sid, r)) keys.push_back({ b->id, pc.first->id, pc.second }); // (#854)
        Factory::WatchMany(keys);
        if (SimTick() - lastChestViews >= 60) { // (#854) 1 Hz
            lastChestViews = SimTick();
            for (const auto& [conn, p] : peers) {
                if (!p.joined) continue;
                auto rit = remotes.find(p.session);
                if (rit == remotes.end()) continue;
                for (const auto& [b, pc] : ChestsNear(p.session, rit->second)) {
                    const Factory::MachineKey key{ b->id, pc.first->id, pc.second };
                    MsgMachineView v; v.body = b->netId; v.plot = pc.first->id; v.part = pc.second; v.view = Factory::Watched(key);
                    if (!v.view.valid) continue;
                    uint64_t h = 1469598103934665603ull;
                    for (const auto& s : v.view.inSlots) { for (char ch : s.first) { h ^= (uint8_t)ch; h *= 1099511628211ull; } h ^= (uint64_t)s.second; h *= 1099511628211ull; }
                    uint64_t& sent = chestViewSent[{ p.session, key }];
                    if (sent == h) continue;
                    sent = h;
                    SendMsg(conn, v);
                }
            }
        }
        // (#854) crafts that pulled from chests: done once the items are here, refused after 3 s
        for (size_t i = 0; i < pendingCrafts.size();) {
            PendingCraft& pc = pendingCrafts[i];
            auto rit = remotes.find(pc.session);
            const uint32_t conn = ConnOf(pc.session);
            const Recipe& rc = itemDB.recipes()[pc.recipe];
            if (rit == remotes.end()) { pendingCrafts.erase(pendingCrafts.begin() + i); continue; }
            Inventory& inv = rit->second.inventory;
            if (inv.canCraft(rc)) {
                const auto before = AllAmounts(inv);
                int made = 0;
                while (made < pc.count && inv.canCraft(rc) && inv.craft(rc)) ++made;
                if (conn) SendChangedItems(conn, inv, before);
                crafts += made;
                if (made >= pc.count || GetTime() > pc.until) { pendingCrafts.erase(pendingCrafts.begin() + i); continue; }
                pc.count -= made;
            } else if (GetTime() > pc.until) {
                if (conn) SendMsg(conn, MsgChat{ 0, std::string("Missing items for ") + rc.output + " (the chests near you didn't have enough)" });
                pendingCrafts.erase(pendingCrafts.begin() + i);
                continue;
            }
            ++i;
        }
        if (SimTick() - lastMachineViews >= 12) { // 5 Hz
            lastMachineViews = SimTick();
            for (const auto& [conn, p] : peers) {
                if (!p.joined) continue;
                auto it = machineWatch.find(p.session);
                if (it == machineWatch.end() || !it->second.on) continue;
                const auto& [bid, plot, part] = it->second.key;
                Planetoid* b = Planetoid::find(bid);
                if (!b) continue;
                MsgMachineView v; v.body = b->netId; v.plot = plot; v.part = part; v.view = Factory::Watched(it->second.key);
                if (v.view.valid) SendMsg(conn, v);
            }
        }
        if (SimTick() - lastMachineStates >= 30) { // 2 Hz
            lastMachineStates = SimTick();
            for (const auto& [conn, p] : peers) {
                if (!p.joined) continue;
                auto rit = remotes.find(p.session);
                if (rit == remotes.end() || !rit->second.player.nearest || !rit->second.player.nearest->plots) continue;
                Planetoid* b = rit->second.player.nearest;
                const Vector3 me = b->toLocalD(rit->second.player.position).toFloat();
                for (const auto& pl : b->plots->plots) {
                    const Vector3 o = { (float)pl->origin.x, (float)pl->origin.y, (float)pl->origin.z };
                    if (Vector3Distance(o, me) > kMachineRange) continue;
                    std::vector<int> parts = Factory::WorkingParts(b->id, pl->id);
                    std::sort(parts.begin(), parts.end());
                    auto& sent = machineStatesSent[{ p.session, { b->id, pl->id } }];
                    if (parts == sent) continue;
                    sent = parts;
                    MsgMachineStates s; s.body = b->netId; s.plot = pl->id; s.working.assign(parts.begin(), parts.end());
                    SendMsg(conn, s);
                }
            }
        }
    }
    // (#1684) The part changes of a dungeon edit, as the clients apply them (in order: the removes and places keep every plot's indices equal).
    void SendPartEdits(uint64_t bodyNet, uint32_t plotId, const std::vector<PartEdit>& edits) {
        for (const PartEdit& e : edits) {
            MsgPlotEdit out; out.session = 0; out.seq = 0; out.body = bodyNet; out.plot = plotId;
            if (e.remove) { out.op = NetPartOp::Remove; out.partIndex = (uint32_t)e.index; }
            else { out.op = NetPartOp::Place; out.part = NetPart{ (uint8_t)e.part.type, e.part.origin.x, e.part.origin.y, e.part.origin.z, e.part.orient, e.part.material }; }
            SendToJoined(out);
        }
    }
    // A vault door (a dungeon_key from the requester's inventory) or a lever (#1684), asked for as a machine action.
    void OnDungeonAction(uint32_t conn, Remote& rm, Planetoid& b, Plot& plot, const MsgMachineAction& a) {
        auto refuse = [&](const char* why) { ++machineActionsRefused; SendMsg(conn, MsgChat{ 0, why }); };
        const Part& pt = plot.parts[a.part];
        const bool vault = a.action == kMachineActUnlockVault;
        if (vault ? !GateLocked(pt) : (pt.type != PartType::Lever || pt.material != 0)) { if (!vault) refuse("The lever is already down"); return; } // (a locked door that is open by now: nothing to say)
        const Int3 c = pt.origin;
        const Vector3 at = plot.toBody({ (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize });
        if (rm.player.nearest != &b || Vector3Distance(at, b.toLocalD(rm.player.position).toFloat()) > kMachineReach) { refuse("That is out of reach"); return; }
        std::vector<PartEdit> edits;
        const auto before = AllAmounts(rm.inventory);
        if (vault) {
            if (!rm.inventory.remove("dungeon_key", 1)) { refuse("Locked: it needs a dungeon key"); return; }
            Dungeon_UnlockVault(plot, a.part, &edits);
            SendMsg(conn, MsgChat{ 0, "The vault door unlocks" });
        } else {
            const int gates = Dungeon_PullLever(plot, a.part, &edits);
            SendMsg(conn, MsgChat{ 0, gates <= 0 ? "The lever clicks, and nothing happens" : "Somewhere, a gate grinds open" });
        }
        SendPartEdits(b.netId, plot.id, edits);
        SendChangedItems(conn, rm.inventory, before);
        ++machineActions;
    }
    void OnMachineAction(uint32_t conn, uint32_t session, const MsgMachineAction& a) {
        auto rit = remotes.find(session);
        Planetoid* b = Planetoid::findNet(a.body);
        Plot* plot = b && b->plots ? b->plots->byId(a.plot) : nullptr;
        auto refuse = [&](const char* why) { ++machineActionsRefused; SendMsg(conn, MsgChat{ 0, why }); };
        if (rit == remotes.end() || rit->second.survival.dead) return;
        if (a.action == kMachineActUnlockVault || a.action == kMachineActPullLever) { // (#1684) a dungeon's door or lever: not a machine
            if (!plot || a.part < 0 || a.part >= (int)plot->parts.size()) { refuse("That isn't there any more"); return; }
            const Vector3 at = plot->toBody({ (plot->parts[a.part].origin.x + 0.5f) * kCellSize, (plot->parts[a.part].origin.y + 0.5f) * kCellSize, (plot->parts[a.part].origin.z + 0.5f) * kCellSize });
            if (ClaimBlocks(session, b, at)) { const std::string why = ClaimRefusal(b, at); refuse(why.c_str()); return; } // (#343)
            OnDungeonAction(conn, rit->second, *b, *plot, a);
            return;
        }
        if (!plot || a.part < 0 || a.part >= (int)plot->parts.size() || !IsMachine(plot->parts[a.part].type)) { refuse("That machine isn't there any more"); return; }
        const Factory::MachineKey key{ b->id, plot->id, a.part };
        if (a.action == 7) { machineWatch[session].on = false; return; }
        // Reach: the machine's origin cell from the player (body-local)
        const Int3 c = plot->parts[a.part].origin;
        const Vector3 at = plot->toBody({ (c.x + 1.0f) * kCellSize, (c.y + 1.0f) * kCellSize, (c.z + 1.0f) * kCellSize });
        if (ClaimBlocks(session, b, at)) { const std::string why = ClaimRefusal(b, at); refuse(why.c_str()); return; } // (#343) no opening another team's chests and machines
        if (rit->second.player.nearest != b || Vector3Distance(at, b->toLocalD(rit->second.player.position).toFloat()) > kMachineReach) { if (a.action != 6) refuse("That machine is out of reach"); return; }
        if (a.action == 6) { machineWatch[session] = { key, true }; return; }
        Inventory& inv = rit->second.inventory;
        const PartType type = plot->parts[a.part].type;
        const auto before = AllAmounts(inv);
        ++machineActions;
        switch (a.action) {
            case 1: case 3: case 8: { // one item into the machine (or its fuel slot; 8: one storage slot of a chest, #862)
                const int n = std::min<int>(a.count, inv.count(a.item));
                const std::vector<std::string> fuels = a.action == 3 ? Factory::FuelsOf(type) : std::vector<std::string>{};
                const bool takes = a.action == 3 ? std::find(fuels.begin(), fuels.end(), a.item) != fuels.end() : Factory::Accepts(type, a.item);
                if (n <= 0 || !takes) { refuse("It doesn't take that"); return; }
                inv.remove(a.item, n);
                if (a.action == 3) Factory::InsertFuel(b->id, plot->id, a.part, a.item, n, session);
                else if (a.action == 8) Factory::InsertSlot(b->id, plot->id, a.part, a.item, n, a.slot, session);
                else Factory::Insert(b->id, plot->id, a.part, a.item, n, session);
                break;
            }
            case 2: { // everything it takes (I): not tools; a chest not what is on the hotbar (as the local I)
                const auto all = inv.all();
                for (const auto& [item, n] : all) {
                    if (n <= 0 || !Factory::Accepts(type, item) || KindOf(item) == ItemKind::Tool) continue;
                    if (inv.remove(item, n)) Factory::Insert(b->id, plot->id, a.part, item, n, session);
                }
                break;
            }
            case 4: Factory::Extract(b->id, plot->id, a.part, session); break;                       // O: the outputs (delivered next tick)
            case 9: // (#862) some of an item out of a chest into the inventory (the scroll wheel)
                if (type != PartType::Chest || a.count <= 0 || a.item.empty()) { refuse("Only a chest gives that way"); return; }
                Factory::Take(b->id, plot->id, a.part, a.item, a.count, session);
                break;
            case 5: // one slot of its screen; with a count (#854), that much of it (its client puts it on its cursor)
                if (a.count > 0) Factory::TakeSlotToCursor(b->id, plot->id, a.part, a.slot, a.which, a.count, session);
                else Factory::TakeSlot(b->id, plot->id, a.part, a.slot, a.which, session);
                break;
            case 10: // (#1103) light (count 1) or put out a campfire
                if (type != PartType::Campfire) { refuse("Only a campfire lights that way"); return; }
                Factory::SetFire(b->id, plot->id, a.part, a.count != 0);
                break;
            default: return;
        }
        SendChangedItems(conn, inv, before);
    }

    // ---- Plants and farming (#335) ----
    // Plants are generated from the seed (stable ids); the server's removals (id -> world time) are the world's: a fell is a
    // request checked here (standing in a loaded tile, within reach, not faster than the tool can fell it), the drops go to
    // the feller's inventory, and every removal goes to every client (MsgPlants), which regrows it on the shared clock.
    // Removals made any other way here (a host's own player) are found by diffing each body's copy-on-write removals.
    // Farm fields live here too: sowing and harvesting are requests; each client gets the fields near it when they change.
    constexpr float kPlantReach = 14.0f;  // m (the tool's 12 m reach, plus movement while the request travels)
    constexpr float kFarmRange = 300.0f;  // m: fields sent to a client around its player (crops draw to 120 m)
    std::unordered_map<uint64_t, std::shared_ptr<const PlantRemovals>> plantsSeen; // by body id: the removals last sent
    long long lastFarmTick = 0;
    int plantsFelled = 0, plantsRefused = 0, farmActions = 0, farmListsSent = 0;
    void SendPlants(uint32_t conn, uint64_t bodyNet, uint32_t session, const std::vector<std::pair<uint64_t, double>>& list) {
        for (size_t i = 0; i < list.size(); i += 4096) {
            MsgPlants m; m.body = bodyNet; m.session = session;
            m.removed.assign(list.begin() + i, list.begin() + std::min(list.size(), i + 4096));
            if (conn) SendMsg(conn, m); else SendToJoined(m);
        }
    }
    void SendAllPlants(Scene& world, uint32_t conn) { // at a join: every body's removals so far
        for (auto& o : world.objects)
            if (auto* p = dynamic_cast<Planetoid*>(o.get()); p && p->trees) {
                std::vector<std::pair<uint64_t, double>> list(p->trees->removals().begin(), p->trees->removals().end());
                if (!list.empty()) SendPlants(conn, p->netId, 0, list);
            }
    }
    void StepPlants(Scene& world) {
        for (auto& o : world.objects) {
            auto* p = dynamic_cast<Planetoid*>(o.get());
            if (!p || !p->trees) continue;
            const auto now = p->trees->removalsShared();
            auto& seen = plantsSeen[p->id];
            if (seen == now) continue;
            std::vector<std::pair<uint64_t, double>> fresh;
            for (const auto& [id, t] : *now) { if (!seen) { fresh.push_back({ id, t }); continue; } auto it = seen->find(id); if (it == seen->end() || it->second != t) fresh.push_back({ id, t }); }
            seen = now;
            if (!fresh.empty()) SendPlants(0, p->netId, 0, fresh);
        }
        // Farm fields near each player, 4 times a second when they changed (growth to 1/200: a crop of an hour every ~18 s)
        if (SimTick() - lastFarmTick < 15) return;
        lastFarmTick = SimTick();
        for (const auto& [conn, peer] : peers) {
            if (!peer.joined) continue;
            auto rit = remotes.find(peer.session);
            if (rit == remotes.end() || !rit->second.player.nearest) continue;
            Planetoid& b = *rit->second.player.nearest;
            MsgFarmFields m; m.body = b.netId;
            Farm_FieldsNear(b, b.toLocalD(rit->second.player.position).toFloat(), kFarmRange, m.fields);
            uint64_t h = 1469598103934665603ull ^ b.id;
            auto mix = [&h](uint64_t v) { h = (h ^ v) * 1099511628211ull; };
            for (const NetFarmField& f : m.fields) { mix((uint64_t)llroundf(f.x * 20.0f)); mix((uint64_t)llroundf(f.y * 20.0f)); mix((uint64_t)llroundf(f.z * 20.0f)); mix(std::hash<std::string>{}(f.crop)); mix((uint64_t)lroundf(f.growth * 200.0f)); }
            if (h == rit->second.farmHash) continue;
            rit->second.farmHash = h;
            SendMsg(conn, m);
            ++farmListsSent;
        }
    }
    void OnPlantAction(uint32_t conn, uint32_t session, const MsgPlantAction& a) {
        auto rit = remotes.find(session);
        Planetoid* b = Planetoid::findNet(a.body);
        if (rit == remotes.end() || rit->second.survival.dead || !b) return;
        Remote& rm = rit->second;
        auto refuse = [&](const char* why) { ++plantsRefused; if (why) SendMsg(conn, MsgChat{ 0, why }); };
        const Vector3 me = rm.player.nearest == b ? b->toLocalD(rm.player.position).toFloat() : Vector3{ 1e30f, 0, 0 };
        const auto before = AllAmounts(rm.inventory);
        if (a.action == 1) { // fell a plant
            Plant pl;
            // (#1107) every refusal is logged with its reason: plants near a fresh spawn couldn't be felled, silently
            auto refuseLog = [&](const char* reason, const char* tell) {
                printf("host: plant %016llx refused for session %u: %s\n", (unsigned long long)a.plant, session, reason);
                refuse(tell);
            };
            if (!b->trees) { refuseLog("the body has no trees here", nullptr); return; }
            if (!b->trees->findPlant(a.plant, pl)) {
                if (b->trees->removals().count(a.plant)) { refuseLog("already felled (its removal is on its way)", nullptr); return; } // (another player's)
                // Not in the server's tile: the tile isn't loaded here, or a clearing (a structure, a spring) took it on the
                // server only. Tell the client it's gone there, so the two worlds agree instead of the player swinging at nothing.
                refuseLog("not in the server's tile (not loaded, or cleared on the server only)", "That plant isn't there on the server");
                SendPlants(conn, b->netId, 0, { { a.plant, WorldTime() } });
                return;
            }
            if (const float d = Vector3Distance(pl.base, me); d > kPlantReach * b->scale) { refuseLog(TextFormat("out of reach (%.1f m)", d), "That's out of reach"); return; }
            if (ClaimBlocks(session, b, pl.base)) { const std::string why = ClaimRefusal(b, pl.base); refuseLog(why.c_str(), why.c_str()); return; } // (#343)
            const std::string held = HeldItem(session); // (#1868) only if the player has it
            const int hits = TerrainTool_HitsToFell(pl, ToolOf(held));
            const double minGap = 0.5 * 0.3 * (hits - 1); // half the swings it takes (latency, a held-on beat)
            if (SimTime() - rm.lastFellAt < minGap) { refuseLog(TextFormat("too soon (%.2f s of %.2f s, %d swings with '%s')", SimTime() - rm.lastFellAt, minGap, hits, held.c_str()), nullptr); return; }
            rm.lastFellAt = SimTime();
            PlantDrops drops = TerrainTool_PlantDrops(*b, pl);
            TerrainTool_ScalePlantDrops(drops, pl, StatsOf(rm)); // (#920) the player's Woodcutting and Foraging perks
            for (const auto& [item, n] : drops.items) rm.inventory.add(item, n);
            b->trees->remove(a.plant);
            plantsSeen[b->id] = b->trees->removalsShared(); // (sent here, with who felled it)
            SendPlants(0, b->netId, session, { { a.plant, b->trees->removals().at(a.plant) } });
            ++plantsFelled;
            if (pl.tree()) { RemoteXp(session, "chop", (double)hits); RemoteXp(session, "tree_felled"); } // (#920) Woodcutting; Foraging
            else RemoteXp(session, pl.stone() ? "loose_stone" : "bush");
        } else { // sow / harvest at a point of farmland ground (#1085)
            const Vector3 ground{ a.x, a.y, a.z };
            if (Vector3Distance(ground, me) > kMachineReach) { refuse("That farmland is out of reach"); return; }
            if (ClaimBlocks(session, b, ground)) { const std::string why = ClaimRefusal(b, ground); refuse(why.c_str()); return; } // (#343)
            std::string why;
            if (a.action == 2) {
                if (!Farm_SowFor(rm.inventory, *b, ground, a.seed, &why)) { refuse(why.empty() ? nullptr : why.c_str()); return; }
                TeamEvent(session, "count", "plant");
                RemoteXp(session, "plant"); // (#920) Farming
            } else {
                if (Farm_HarvestFor(rm.inventory, *b, ground, &why, StatsOf(rm)) <= 0) { refuse(why.empty() ? nullptr : why.c_str()); return; } // (#920) their Farming perks
                TeamEvent(session, "count", "harvest");
                RemoteXp(session, "harvest"); // (#920) Farming
            }
            for (auto& [s, r] : remotes) r.farmHash = 0; // (everyone near sees it at once)
            lastFarmTick = 0;
            ++farmActions;
        }
        SendChangedItems(conn, rm.inventory, before);
    }

    void StepTeamQuests(Scene& world) {
        if (SimTick() - lastQuestTick < 60) return; // once a second
        lastQuestTick = SimTick();
        // Nights survived, as single player counts them (survival.cpp): per player, for their team
        const Vector3 sunDir = Vector3Normalize(world.sunDirection);
        for (auto& [sid, r] : remotes) {
            if (!r.player.nearest) continue;
            bool dark = Vector3DotProduct(sunDir, r.player.up) < -0.1f;
            // tests: RAYLIBGAME_HOST_NIGHT="a,b": night for everyone from a to b seconds after the first join (a real night
            // here is half an hour of game time; a sped-up world makes streaming chase the spin)
            if (const char* nt = getenv("RAYLIBGAME_HOST_NIGHT")) {
                static double firstJoin = -1.0;
                if (firstJoin < 0.0) firstJoin = GetTime();
                double a = 0.0, b = 0.0;
                if (sscanf(nt, "%lf,%lf", &a, &b) == 2) { const double t = GetTime() - firstJoin; dark = t >= a && t < b; }
            }
            if (dark && !r.inNight) { r.inNight = true; r.deathsAtDusk = r.survival.deaths; }
            else if (!dark && r.inNight) { r.inNight = false; if (r.survival.deaths == r.deathsAtDusk) { TeamEvent(sid, "count", "nights"); RemoteXp(sid, "night"); } }
            // (#338 slice 3) exploring and visiting, from the server's own view of the player: the map tiles around it (the
            // local map's tiles and radius) for its team, and a moon or the ring within two radii
            if (const int team = TeamOfSession(sid)) {
                Planetoid* b = r.player.nearest;
                if (auto pf = std::dynamic_pointer_cast<const PlanetField>(b->field); pf && !b->flat) {
                    const DVec3 l = b->toLocalD(r.player.position);
                    if (l.length() > 1.0) {
                        TeamQuests_Explore(team, b->name, Vector3Normalize(l.toFloat()), pf->radius, r.record.uuid); // (#1886) the player's own tiles too
                        Nav_RevealTiles(Vector3Normalize(l.toFloat()), pf->radius, r.tilesAlone); // (tests: what it revealed alone)
                    }
                }
                if ((b->kind == BodyKind::Moon || b->kind == BodyKind::Asteroid) && b->toLocalD(r.player.position).length() < 2.0 * b->radius * b->scale)
                    TeamQuests_Visit(team, b->kind == BodyKind::Moon);
            }
        }
        std::map<int, std::vector<std::pair<uint32_t, uint32_t>>> online; // team -> (session, conn)
        for (const auto& [conn, p] : peers) if (p.joined) if (const int t = TeamOfSession(p.session)) online[t].push_back({ p.session, conn });
        for (auto& [team, members] : online) {
            const uint32_t payee = questContributor.count(team) ? questContributor[team] : members.front().first;
            TeamQuests_Update(team,
                [&](const std::string& item, int n) { if (Host_Give(payee, item, n)) ++questRewards; }, // (#816) deltas to that player
                [&](const std::string& text) { for (const auto& [s, c] : members) SendMsg(c, MsgChat{ 0, "Quest: " + text }); });
            const std::string text = TeamQuests_Text(team);
            if (text != questBookSent[team]) {
                questBookSent[team] = text;
                for (const auto& [s, c] : members) SendQuestBook(c, team);
            }
        }
    }

    void SendVitals() {
        if (SimTick() - lastVitals < kVitalsEvery) return;
        lastVitals = SimTick();
        for (const auto& [conn, p] : peers) {
            if (!p.joined) continue;
            auto rit = remotes.find(p.session);
            if (rit == remotes.end()) continue;
            const SurvivalState& sv = rit->second.survival;
            MsgVitals v;
            v.health = sv.health; v.food = sv.food; v.oxygen = sv.oxygen; v.power = sv.power; v.temperature = sv.temperature;
            v.flags = (uint8_t)((sv.dead ? 1 : 0) | (sv.breathable ? 2 : 0) | (sv.inVacuum ? 4 : 0) | (sv.underwater ? 8 : 0) | (sv.inLifeSupport ? 16 : 0) |
                                (SimTime() < sv.sickUntil ? 32 : 0) | (sv.headlamp ? 64 : 0) | (rit->second.creative ? 128 : 0)); // 128: creative (#823)
            v.deaths = (uint16_t)std::min(sv.deaths, 65535);
            NetWriter w; NetWrite(w, v);
            Net_Send(conn, w.bytes.data(), (uint32_t)w.bytes.size(), false);
        }
    }
    uint32_t ConnOf(uint32_t session); // (below)
    template <class M> void SendMsg(uint32_t conn, const M& m); // (below)
    int skillStatesSent = 0, skillSpends = 0, skillSpendsRefused = 0; // (#920 tests)
    // (#920) Each player's skills, to that player, when they changed: at most every 2 s (a level-up or the join at once).
    void SendSkills() {
        for (const auto& [conn, p] : peers) {
            if (!p.joined) continue;
            auto rit = remotes.find(p.session);
            if (rit == remotes.end() || !rit->second.skillsDirty || GetTime() - rit->second.skillsSentAt < 2.0) continue;
            SendMsg(conn, MsgSkillState{ Skills_SaveText(rit->second.skills) });
            rit->second.skillsDirty = false; rit->second.skillsSentAt = GetTime();
            ++skillStatesSent;
        }
    }
    // (#920) XP for a remote player's event (skills.txt's xp lines), at the world's rate; creative earns none. Sent with the
    // next SkillState (a level-up at once).
    void RemoteXp(uint32_t session, const char* event, double n) {
        auto it = remotes.find(session);
        if (it == remotes.end() || it->second.creative) return;
        const int ups = Skills_Event(it->second.skills, event, n, worldXpRate);
        it->second.skillsDirty = true;
        if (ups > 0) it->second.skillsSentAt = -1e9;
        TeamXp(session, event, n);
    }
    // (#920, the user's answer on #167) Teammates near the earner get a share of the XP (Skills_TeamShare: 10%, within
    // teamXpRange). n: the earner's event count after the anti-farming share.
    int teamXpGiven = 0; // tests: shares handed out
    Scene* hostWorld = nullptr; // (#920) the world Host_Update runs, for body-local positions where no world is passed in
    void TeamXp(uint32_t session, const char* event, double n) {
        const auto me = remotes.find(session);
        const int team = TeamOfSession(session);
        if (me == remotes.end() || team == 0 || !hostWorld) return;
        // Distances body-local, on the same body: a player not stepped for a while (no inputs) keeps a world position the
        // planet has long turned away from (#844), but its body-local anchor is right (Player::netState).
        Player::NetState mine;
        if (!me->second.player.netState(*hostWorld, mine)) return;
        const Planetoid* body = Planetoid::findNet(mine.body);
        if (!body) return;
        for (auto& [id, r] : remotes) {
            if (id == session || r.creative || TeamOfSession(id) != team) continue;
            Player::NetState theirs;
            if (!r.player.netState(*hostWorld, theirs) || theirs.body != mine.body) continue;
            const double share = Skills_TeamShare(n, (theirs.local - mine.local).length() * body->scale);
            if (share <= 0.0) continue;
            if (Skills_Event(r.skills, event, share, worldXpRate) > 0) r.skillsSentAt = -1e9;
            r.skillsDirty = true;
            ++teamXpGiven;
        }
    }
    // (#920) A stroke's or fill's cell for the anti-farming share: the body-local point on a 1 m grid (as skillplay's).
    uint64_t SkillCell(Vector3 local) {
        const auto q = [](float v) { return (uint64_t)(uint32_t)(int32_t)std::floor(v) & 0x1FFFFFull; };
        return q(local.x) | (q(local.y) << 21) | (q(local.z) << 42);
    }
    int animalBites = 0;
    // (#334 slice 2) the remote players for the animals to mind, each frame; a bite on one:
    // the server's copy takes it (Combat_DamageFor: health, the shove, a held block facing it cuts both - #337); its client gets the hit's feel.
    void SetAnimalTargets(Scene& world) {
        gAnimalTargets.clear();
        for (const Session& s : sessions) {
            if (!s.online) continue;
            auto rit = remotes.find(s.id);
            if (rit == remotes.end() || !rit->second.player.nearest) continue;
            const Planetoid* b = rit->second.player.nearest;
            gAnimalTargets.push_back({ b->netId, b->toLocalD(rit->second.player.position).toFloat(), s.id, !rit->second.player.flyMode && !rit->second.survival.dead, HeldItem(s.id) }); // (#337: the dead aren't prey; #1868 held: only if had)
        }
        // (#1223/#1224, Plan 12 testing) RAYLIBGAME_JOIN_CREATURE=1 (tools/test_join.sh): once two players have been playing for
        // 3 s, a thornback (a creature on the brain) by the first; it hunts and bites them; from 6 s on the server hits it for 10
        // every second until it falls. Both clients report the health they saw and the death (the replication of one creature).
        // (#1421) "Playing": each has walked 3 m from where it first stood - the test's players walk once their loading screens
        // are gone. Online was not enough (on the runner the join sync took 35-67 s), nor were inputs (they flow during the ~20 s
        // spawn-area loading screen after the sync): the creature lived and fell before the clients were in the world.
        static const bool joinCreature = getenv("RAYLIBGAME_JOIN_CREATURE") && getenv("RAYLIBGAME_JOIN_CREATURE")[0] == '1';
        if (joinCreature) {
            static double bothAt = -1.0, spawnedAt = -1.0, lastHit = 0.0; static uint64_t net = 0; static bool fell = false;
            static std::unordered_map<uint32_t, Vector3> firstAt; static std::unordered_set<uint32_t> walked;
            const double now = GetTime();
            for (const auto& kv : remotes) {
                uint64_t body = 0; Vector3 at{}; double scale = 1.0;
                if (kv.second.lastQueued == 0 || !RemoteLocal(kv.second, body, at, scale)) continue;
                auto it = firstAt.find(kv.first);
                if (it == firstAt.end()) firstAt[kv.first] = at;
                else if (Vector3Distance(it->second, at) * (float)scale > 3.0f) walked.insert(kv.first);
            }
            const int playing = (int)walked.size();
            if (playing >= 2 && gAnimalTargets.size() >= 2 && bothAt < 0.0) {
                bothAt = now; Animals_SetThreats(true); Animals_SetNaturalSpawns(false);
                printf("host: creature: both players in the world\n");
            }
            if (bothAt >= 0.0 && spawnedAt < 0.0 && now - bothAt > 3.0 && !gAnimalTargets.empty()) {
                const AnimalTarget& t = gAnimalTargets[0];
                Planetoid* b = Planetoid::findNet(t.body);
                const Vector3 up = Vector3Normalize(t.local);
                const Vector3 side = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const int before = Animals_Count();
                if (b) Animals_Spawn(world, *b, "thornback", Vector3Add(t.local, Vector3Scale(side, 5.0f)), 1, false, 0.01f);
                if (Animals_Count() > before) { net = Animals_NetId(Animals_Count() - 1); printf("host: creature: a thornback (%llu) by session %u\n", (unsigned long long)net, t.session); }
                spawnedAt = now; lastHit = now;
            }
            // (#1421) The test is the creature's replication, not the players' survival: while it lives, neither test player
            // goes below half health. On the runner a ram threw bob into the air and the hit killed him; he spent the
            // creature's whole life on the respawn screen and never saw it.
            if (net && !fell) for (auto& kv : remotes) kv.second.survival.health = std::max(kv.second.survival.health, 50.0f);
            if (net && !fell && spawnedAt >= 0.0 && now - spawnedAt > 6.0 && now - lastHit >= 1.0) {
                lastHit = now;
                const int i = Animals_FindNet(net);
                if (i < 0) { fell = true; printf("host: creature: fell (bites on players %d)\n", animalBites); }
                else { DamageEvent e; e.amount = 10.0f; Animals_HitEvent(i, e); if (Animals_FindNet(net) < 0) { fell = true; printf("host: creature: fell (bites on players %d)\n", animalBites); } else printf("host: creature: health %.0f\n", Animals_HealthOf(i)); }
            }
        }
        if (!gAnimalDamageRemote) gAnimalDamageRemote = [](uint32_t session, const DamageEvent& e) {
            auto rit = remotes.find(session);
            if (rit == remotes.end()) return;
            Player& pl = rit->second.player;
            rit->second.survival.blocking = (rit->second.lastHeld >> (int)Act::Secondary) & 1ull; // (its held input; its stamina is simulated here too)
            Combat_DamageFor(rit->second.survival, pl, e); // (#337: the server's copy takes the damage and the shove)
            const uint32_t conn = ConnOf(session);
            if (!conn || !pl.nearest) return;
            MsgDamage m; m.amount = e.amount; m.kind = (uint8_t)e.kind; m.hasSource = e.hasSource; m.knockback = e.knockback; m.sourceId = e.sourceId;
            m.body = pl.nearest->netId; m.local = pl.nearest->toLocalD(e.source);
            SendMsg(conn, m);
            ++animalBites;
        };
    }
    uint32_t ConnOf(uint32_t session) { for (const Session& s : sessions) if (s.id == session && s.online) return s.conn; return 0; }
    void OnAnimalHit(Scene& world, uint32_t session, uint64_t netId, int action = 0) {
        auto rit = remotes.find(session);
        uint64_t body = 0; Vector3 at{}; double scale = 1.0;
        const int i = Animals_FindNet(netId);
        if (rit == remotes.end() || i < 0 || !RemoteLocal(rit->second, body, at, scale)) { ++animalHitsRefused; return; }
        const Planetoid* b = rit->second.player.nearest;
        if (!b || Animals_Body(i) != b->name) { ++animalHitsRefused; return; }
        // (#1108) lag compensation: within reach of where the animal was in the last kHitRewind (the player swung at what
        // they saw: a round trip and their interpolation behind); a fleeing animal's hits were refused at the edge of reach
        const float nowD = Vector3Distance(Animals_Position(i), at) * (float)scale, rewD = Animals_HitDistance(i, at, kHitRewind) * (float)scale;
        if (rewD > kAnimalReach) { ++animalHitsRefused; return; }
        if (nowD > kAnimalReach) ++animalHitsRewound; // (allowed only by the rewind)
        if (ClaimBlocks(session, b, Animals_Position(i))) { ++animalHitsRefused; if (const uint32_t conn = ConnOf(session)) SendMsg(conn, MsgChat{ 0, ClaimRefusal(b, Animals_Position(i)) + ": its animals are theirs" }); return; } // (#343)
        double& last = lastAnimalHit[session];
        if (GetTime() - last < kAnimalHitCooldown) { ++animalHitsRefused; return; }
        last = GetTime();
        if (action != 0) { // (#334 slice 2) feed it, or take its product: from and into that player's inventory as the server keeps it
            std::string why; std::map<std::string, double> before;
            const bool done = Animals_UseNet(netId, action, rit->second.inventory, why, before);
            if (const uint32_t conn = ConnOf(session)) {
                if (done) SendInventoryDeltas(conn, rit->second.inventory, before);
                if (done && action == 1) TeamEvent(session, "count", "feed"); // (#338)
                else SendMsg(conn, MsgChat{ 0, (action == 1 ? "It won't eat: " : "Nothing to take: ") + why });
            }
            return;
        }
        const std::string held = HeldItem(session); // (#1868) only if the player has it
        const ToolSpec t = ToolOf(held);
        const float damage = t.type == ToolType::Hands ? 2.0f : 3.0f + 3.0f * t.tier;
        std::vector<std::pair<std::string, int>> drops;
        ++animalHits;
        if (Animals_HitNet(netId, damage, drops)) {
            TeamEvent(session, "count", "slain"); // (#338) a kill counts for the hunter's team
            // The drops into the hunter's inventory as the server keeps it, and the change to its client in one go (UI/UX's pattern).
            Remote& rm = rit->second;
            std::map<std::string, double> before;
            for (const auto& [item, n] : drops) before[item] = rm.inventory.amount(item);
            std::string got;
            for (const auto& [item, n] : drops) { rm.inventory.add(item, n); got += TextFormat(" +%d %s", n, item.c_str()); }
            for (const Session& s : sessions) if (s.id == session && s.online && s.conn) SendInventoryDeltas(s.conn, rm.inventory, before);
            printf("host: session %u felled animal %llx:%s\n", session, (unsigned long long)netId, got.c_str());
        }
        (void)world;
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

    void ForgetSync(uint32_t conn); // (below, with the join sync)
    bool SaveRecord(uint32_t session); // (#346, below)
    void StartSync(Scene& world, uint32_t conn, uint32_t session);
    void Drop(uint32_t conn, NetBye why) {
        auto it = peers.find(conn);
        if (it == peers.end()) return;
        for (Session& s : sessions) if (s.id == it->second.session && s.conn == conn) { s.online = false; s.conn = 0; }
        if (why != NetBye::Quit) SendMsg(conn, MsgBye{ why });
        HostLog("connection %u (session %u) %s", conn, it->second.session, NetByeText(why));
        if (eventHook && it->second.joined) { // (#345 audit)
            std::string name; for (const Session& s : sessions) if (s.id == it->second.session) name = s.name;
            eventHook("leave", it->second.session, name + " (" + NetByeText(why) + ")");
        }
        const bool wasJoined = it->second.joined;
        const uint32_t session = it->second.session;
        Interest_Forget(conn); // (#340) its accumulators
        ForgetSync(conn);      // (#776) its join sync, if it was still running
        deltas.erase(conn);    // (#769) its snapshot bases
        animalsSent.erase(conn); // (#334) the animals it was sent
        peers.erase(it);
        Net_Close(conn);
        // (#346) their record written, and the player gone from the world until they come back (unless the session is
        // already back on another connection)
        bool offline = true;
        for (const Session& s : sessions) if (s.id == session && s.online) offline = false;
        if (wasJoined && offline) {
            if (auto rit = remotes.find(session); rit != remotes.end() && !rit->second.tilesAlone.empty())
                printf("host: session %u revealed %zu map tiles itself; its team explored %.0f\n", session, rit->second.tilesAlone.size(), TeamQuests_Explored(TeamOfSession(session)));
            SaveRecord(session); remotes.erase(session);
        }
        if (wasJoined) SendRoster(); // the others see them go
    }
    void Refuse(uint32_t conn, NetReject why, const std::string& detail) {
        SendMsg(conn, MsgReject{ why, detail });
        HostLog("refused connection %u: %s%s%s", conn, NetRejectText(why), detail.empty() ? "" : ": ", detail.c_str());
        peers.erase(conn);
        Net_Close(conn); // (lingers: the reject goes out first)
    }

    // (#349) names and chat text from clients
    bool NameAcceptable(const std::string& n) {
        if (n.empty() || n.size() > 24 || n.front() == ' ' || n.back() == ' ') return false;
        for (unsigned char ch : n) if (ch < 0x20 || ch == 0x7F) return false; // (UTF-8 bytes >= 0x80 are fine: names in any script)
        return true;
    }
    bool SameName(const std::string& a, const std::string& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) if (tolower((unsigned char)a[i]) != tolower((unsigned char)b[i])) return false;
        return true;
    }
    // Player records (#346): a name is claimed by whoever has played under it here (their record, or a session this run).
    bool NameClaimedByOther(const std::string& name, const std::string& uuid) {
        for (const Session& x : sessions) if (x.uuid != uuid && SameName(x.name, name)) return true;
        for (const PlayerRecordMeta& m : PlayerRecord_List()) if (m.uuid != uuid && SameName(m.name, name)) return true;
        return false;
    }
    // A joining player's record, if they have one: their place, inventory and survival come back (the client is sent the
    // inventory); a first-timer starts at the world's start (SpawnRemote) with nothing.
    bool LoadRecord(Scene& world, uint32_t conn, const Session& s) { // true: their record was loaded
        Remote& rm = remotes[s.id];
        rm.record = PlayerRecordMeta{};
        rm.skills = SkillState{}; rm.skillsDirty = true; rm.skillsSentAt = -1e9; // (#920) (sent with the next update)
        const std::string path = PlayerRecord_Path(s.uuid);
        bool loaded = false;
        if (FileExists(path.c_str()) && PlayerRecord_Read(path, rm.player, rm.inventory, rm.survival, rm.record, &rm.skills)) {
            loaded = true;
            rm.player.flyAllowed = false; // (#781) never for a client's server-side player
            // (#856) anchored where the record put them: until their first step (their client still loading) the server reports
            // the anchor, and SpawnRemote's was the world's start
            if (Planetoid* b = rm.player.nearest ? rm.player.nearest : StartBody(world)) rm.player.anchorTo(*b);
            std::map<std::string, double> before;
            for (const InvSlot& c : rm.inventory.grid()) if (!c.empty()) before[c.id] = 0.0; // the client starts from nothing
            SendInventoryDeltas(conn, rm.inventory, before);
        }
        rm.record.uuid = s.uuid; rm.record.name = s.name;
        if (rm.record.firstSeen <= 0.0) rm.record.firstSeen = (double)time(nullptr);
        rm.joinedAt = GetTime();
        Teams_Join(s.uuid, s.name, 0); // a team of one to start with (#604 forms bigger ones)
        HostLog("%s (%s): %s", s.name.c_str(), s.uuid.c_str(), loaded ? "their record loaded" : "a new player");
        return loaded;
    }
    bool SaveRecord(uint32_t session) {
        auto it = remotes.find(session);
        if (it == remotes.end() || !PlayerId_Valid(it->second.record.uuid)) return false;
        Remote& rm = it->second;
        const double now = GetTime();
        const bool ok = PlayerRecord_Write(PlayerRecord_Path(rm.record.uuid), rm.player, rm.inventory, rm.survival, rm.record, now - rm.joinedAt, &rm.skills);
        if (ok) { rm.record.playTime += now - rm.joinedAt; rm.joinedAt = now; }
        else HostLog("could not write the record of %s", rm.record.name.c_str());
        return ok;
    }
    int SaveAllRecords() { int n = 0; for (const auto& [id, rm] : remotes) n += SaveRecord(id); return n; }

    // ---- Teams (#604): /team commands and team chat, from any player; the server decides, members hear TeamState ----
    // Defaults posted on #604: at most 8 a team; invites to players online, for 2 minutes; joining takes the team's quest
    // progress (a team of one left behind goes); leaving or being kicked starts a team of one with a copy of it; the owner
    // leaving passes ownership to the earliest member.
    struct TeamInvite { int team = 0; std::string to; double until = 0.0; };
    std::vector<TeamInvite> teamInvites;
    double nextTeamStates = 0.0;
    constexpr int kMaxTeam = 8;
    constexpr double kInviteSeconds = 120.0;
    const Session* SessionOfUuid(const std::string& uuid) { for (const Session& s : sessions) if (s.uuid == uuid && s.online) return &s; return nullptr; }
    const Session* OnlineNamed(const std::string& name) { for (const Session& s : sessions) if (s.online && SameName(s.name, name)) return &s; return nullptr; }
    std::string NameOfUuid(const std::string& uuid) {
        for (const Session& s : sessions) if (s.uuid == uuid) return s.name;
        for (const PlayerRecordMeta& m : PlayerRecord_List()) if (m.uuid == uuid) return m.name;
        return "a player";
    }
    void Tell(const Session& s, const std::string& text) { SendMsg(s.conn, MsgChat{ 0, text }); }
    void TellTeam(int team, const std::string& text) {
        if (const TeamRecord* t = Teams_Find(team)) for (const std::string& u : t->members) if (const Session* s = SessionOfUuid(u)) Tell(*s, text);
    }
    // (#828) An offline player's "last seen": seconds since their record file was last written (on leaving, and by autosaves
    // while they played); 0 when there is none.
    uint32_t SecondsSinceRecord(const std::string& uuid) {
        std::error_code ec;
        const auto when = std::filesystem::last_write_time(PlayerRecord_Path(uuid), ec);
        if (ec) return 0;
        const long long s = std::chrono::duration_cast<std::chrono::seconds>(std::filesystem::file_time_type::clock::now() - when).count();
        return s > 0 ? (uint32_t)std::min<long long>(s, 0xFFFFFFFFll) : 1;
    }
    void SendTeamState(const Session& s) {
        MsgTeamState m;
        m.team = Teams_Of(s.uuid);
        if (const TeamRecord* t = Teams_Find(m.team)) {
            m.name = t->name;
            const std::string owner = Teams_OwnerOf(m.team);
            for (const std::string& u : t->members) {
                MsgTeamState::Member x;
                x.name = NameOfUuid(u); x.owner = u == owner; x.you = u == s.uuid;
                if (const Session* o = SessionOfUuid(u)) { x.online = true; x.session = o->id; auto it = remotes.find(o->id); x.health = it != remotes.end() ? it->second.survival.health : 0.0f; }
                else x.lastSeen = SecondsSinceRecord(u); // (#828) their record is written when they leave
                m.members.push_back(x);
            }
            m.claimBudget = Claims_Budget(TeamQuests_BodiesReached(m.team)); // (#828)
            for (const int a : Claims_AlliesOf(m.team)) if (const TeamRecord* at = Teams_Find(a)) m.allies.push_back(at->name);
        }
        for (const TeamInvite& i : teamInvites) if (i.to == s.uuid && GetTime() < i.until) if (const TeamRecord* t = Teams_Find(i.team)) m.invite = t->name;
        SendMsg(s.conn, m);
    }
    void SendTeamStates() {
        if (GetTime() < nextTeamStates) return;
        nextTeamStates = GetTime() + 1.0; // (health changes: once a second)
        teamInvites.erase(std::remove_if(teamInvites.begin(), teamInvites.end(), [](const TeamInvite& i) { return GetTime() >= i.until || !Teams_Find(i.team); }), teamInvites.end());
        for (const Session& s : sessions) if (s.online) SendTeamState(s);
    }
    // Into another team. A joiner who was alone keeps their progress: their book is merged into the team's (#861, the user
    // on #604), then it goes with the emptied team. One who moves straight from a team that still has members brings nothing
    // across: that team keeps its book (leaving or a kick is SplitOff: a team of one with a copy, unchanged), except the
    // map they explored themselves, which comes with them (the user, on #1886).
    void JoinTeam(const Session& s, int team) {
        const int from = Teams_Of(s.uuid);
        if (from) TeamQuests_Book(from); // (loaded from its record now: Teams_Move removes an emptied team's record)
        const int emptied = Teams_Move(s.uuid, team);
        if (emptied) { if (emptied == from) TeamQuests_Merge(from, team); TeamQuests_Forget(emptied); Claims_MoveTeam(emptied, team); } // (#343) its claims go with them
        if (Teams_Of(s.uuid) == team) TeamQuests_CarryPlayer(s.uuid, team); // (#1886)
        SendQuestBook(s.conn, team);
    }
    // Out of their team into a team of one with a copy of its progress (leave, kick). Online or not.
    void SplitOff(const std::string& uuid, const std::string& name) {
        const int from = Teams_Of(uuid);
        const int solo = Teams_CreateSolo(uuid, name);
        if (from) TeamQuests_Copy(from, solo);
        if (const Session* s = SessionOfUuid(uuid)) SendQuestBook(s->conn, solo);
    }
    // A chat line that is a team command or team chat: handled (true), with a reply to the sender.
    bool TeamCommand(uint32_t session, const std::string& text) {
        const Session* me = nullptr;
        for (const Session& s : sessions) if (s.id == session && s.online) me = &s;
        if (!me) return false;
        const bool teamChat = text.rfind("/t ", 0) == 0;
        if (!teamChat && text != "/team" && text.rfind("/team ", 0) != 0) return false;
        const int mine = Teams_Of(me->uuid);
        TeamRecord* team = Teams_Find(mine);
        const bool owner = team && Teams_OwnerOf(mine) == me->uuid;
        if (teamChat) {
            const std::string said = text.substr(3);
            if (!said.empty()) TellTeam(mine, "[team] " + me->name + ": " + said);
            return true;
        }
        std::istringstream in(text.substr(5));
        std::string cmd, arg;
        in >> cmd; std::getline(in, arg);
        while (!arg.empty() && arg[0] == ' ') arg.erase(0, 1);
        nextTeamStates = 0.0; // (whatever changes, members hear it at once)
        if (cmd.empty()) {
            std::string list;
            if (team) for (const std::string& u : team->members) list += (list.empty() ? "" : ", ") + NameOfUuid(u) + (u == Teams_OwnerOf(mine) ? " (owner)" : "");
            Tell(*me, TextFormat("Your team: %s: %s. /team invite NAME, accept, decline, leave, kick NAME, name TEXT; /t TEXT says to your team", team ? team->name.c_str() : "none", list.c_str()));
        } else if (cmd == "invite") {
            const Session* to = OnlineNamed(arg);
            if (!to) Tell(*me, arg + " isn't online");
            else if (!team) Tell(*me, "You have no team");
            else if (Teams_Of(to->uuid) == mine) Tell(*me, to->name + " is in your team already");
            else if ((int)team->members.size() >= kMaxTeam) Tell(*me, TextFormat("A team holds %d at most", kMaxTeam));
            else {
                teamInvites.erase(std::remove_if(teamInvites.begin(), teamInvites.end(), [&](const TeamInvite& i) { return i.to == to->uuid && i.team == mine; }), teamInvites.end());
                teamInvites.push_back({ mine, to->uuid, GetTime() + kInviteSeconds });
                { const TeamRecord* theirs = Teams_Find(Teams_Of(to->uuid)); // (#861) alone: their progress joins the team's; in a team: that team keeps it
                  Tell(*to, me->name + " invites you to " + team->name + (theirs && theirs->members.size() > 1 ? ": /team accept (you share that team's quest progress; yours stays with your current team), or /team decline" : ": /team accept (your quest progress is added to the team's), or /team decline")); }
                Tell(*me, "Invited " + to->name);
            }
        } else if (cmd == "accept") {
            const TeamInvite* inv = nullptr;
            for (const TeamInvite& i : teamInvites) if (i.to == me->uuid && GetTime() < i.until && Teams_Find(i.team)) inv = &i;
            if (!inv) Tell(*me, "No invite waiting");
            else if ((int)Teams_Find(inv->team)->members.size() >= kMaxTeam) Tell(*me, "That team is full");
            else {
                const int to = inv->team;
                teamInvites.erase(std::remove_if(teamInvites.begin(), teamInvites.end(), [&](const TeamInvite& i) { return i.to == me->uuid; }), teamInvites.end());
                JoinTeam(*me, to);
                TellTeam(to, me->name + " joined the team");
            }
        } else if (cmd == "decline") {
            teamInvites.erase(std::remove_if(teamInvites.begin(), teamInvites.end(), [&](const TeamInvite& i) { return i.to == me->uuid; }), teamInvites.end());
            Tell(*me, "Invite declined");
        } else if (cmd == "leave") {
            if (!team || team->members.size() <= 1) Tell(*me, "You're a team of one already");
            else { const std::string name = team->name; SplitOff(me->uuid, me->name); TellTeam(mine, me->name + " left the team"); Tell(*me, "You left " + name + ": a team of one now, with its progress so far"); }
        } else if (cmd == "kick") {
            std::string target;
            if (team) for (const std::string& u : team->members) if (SameName(NameOfUuid(u), arg)) target = u;
            if (!owner) Tell(*me, "Only the team's owner can kick");
            else if (target.empty() || target == me->uuid) Tell(*me, arg + " isn't in your team");
            else {
                SplitOff(target, NameOfUuid(target));
                TellTeam(mine, NameOfUuid(target) + " was removed from the team");
                if (const Session* s = SessionOfUuid(target)) Tell(*s, "You were removed from the team: a team of one now, with its progress so far");
            }
        } else if (cmd == "name") {
            if (!owner) Tell(*me, "Only the team's owner can rename it");
            else if (arg.empty() || !NameAcceptable(arg)) Tell(*me, "/team name TEXT (1-24 characters)");
            else { Teams_SetName(mine, arg); TellTeam(mine, "The team is now " + arg); }
        } else Tell(*me, "Unknown: /team " + cmd + " (/team lists the commands)");
        return true;
    }

    // ---- Land claims (#343): a team's claims protect what's inside from other teams (ops may); /claim, /unclaim ----
    int claimsRefused = 0, claimsMade = 0;
    double nextClaimStep = 0.0;
    std::unordered_map<uint32_t, uint32_t> claimInside; // by session: the claim it stood in at the last step (enter/leave lines)
    int TeamOfName(const std::string& name) { // (#343) a player online or on record
        for (const Session& s : sessions) if (SameName(s.name, name)) return Teams_Of(s.uuid);
        for (const PlayerRecordMeta& m : PlayerRecord_List()) if (SameName(m.name, name)) return Teams_Of(m.uuid);
        return 0;
    }
    std::string TeamName(int team) { const TeamRecord* t = Teams_Find(team); return t ? t->name : "another team"; }
    bool IsOpSession(uint32_t session) { for (const Session& s : sessions) if (s.id == session) return Admin_RoleOf(s.name, s.uuid) != AdminRole::Player; return false; } // (#1086: by id too)
    bool ClaimBlocks(uint32_t session, const Planetoid* b, Vector3 local, bool count) {
        if (!b) return false;
        const Claim* c = Claims_At(b->name, local);
        if (!c || Claims_Allows(c->team, TeamOfSession(session)) || IsOpSession(session)) return false; // its team and allies may
        if (count) ++claimsRefused;
        return true;
    }
    std::string ClaimRefusal(const Planetoid* b, Vector3 local) {
        const Claim* c = b ? Claims_At(b->name, local) : nullptr;
        return c ? TeamName(c->team) + " has claimed this land" : "This land is claimed";
    }
    int TeamNamed(const std::string& name) { for (const TeamRecord& t : Teams_All()) if (SameName(t.name, name)) return t.id; return 0; }
    // Players claim with banners (#343, the user 2026-11-18); /claim and /unclaim stay for ops. Everyone: /claims (the team's
    // claims, budget and allies), /claim info, and the owner: /claim ally TEAM, /claim unally TEAM.
    bool ClaimCommand(uint32_t session, const std::string& text) {
        if (text != "/claim" && text != "/unclaim" && text != "/claims" && text.rfind("/claim ", 0) != 0) return false;
        const Session* me = nullptr;
        for (const Session& s : sessions) if (s.id == session && s.online) me = &s;
        auto rit = remotes.find(session);
        if (!me || rit == remotes.end()) return true;
        const int team = Teams_Of(me->uuid);
        const bool op = IsOpSession(session);
        Planetoid* b = rit->second.player.nearest;
        const Vector3 here = b ? b->toLocalD(rit->second.player.position).toFloat() : Vector3{};
        const int budget = Claims_Budget(TeamQuests_BodiesReached(team));
        if (text == "/claims" || text == "/claim list") {
            const std::vector<Claim> mine = Claims_OfTeam(team);
            std::string list, al;
            for (const Claim& c : mine) {
                list += TextFormat("%s#%u (%.0f m) on %s", list.empty() ? "" : ", ", c.id, c.radius * 2.0f, c.body.c_str());
                if (b && c.body == b->name) list += TextFormat(" %.0f m away", Vector3Distance(c.centre, here));
            }
            for (int a : Claims_AlliesOf(team)) al += (al.empty() ? "" : ", ") + TeamName(a);
            Tell(*me, TextFormat("Your team's claims: %d of %d%s%s; allies: %s. Place a claim banner to claim land; break it to unclaim", (int)mine.size(), budget,
                                 list.empty() ? "" : ": ", list.c_str(), al.empty() ? "none (/claim ally TEAM)" : al.c_str()));
            return true;
        }
        if (text.rfind("/claim ally ", 0) == 0 || text.rfind("/claim unally ", 0) == 0) {
            const bool on = text.rfind("/claim ally ", 0) == 0;
            const std::string name = text.substr(on ? 12 : 14);
            const int other = TeamNamed(name);
            if (Teams_OwnerOf(team) != me->uuid && !op) Tell(*me, "Only your team's owner can change its allies");
            else if (!other || other == team) Tell(*me, "No other team named " + name + " (/team lists yours; team names are what /team name set)");
            else {
                Claims_SetAlly(team, other, on);
                TellTeam(team, (on ? TeamName(other) + " may now build and open things on our claims" : TeamName(other) + " is no longer allowed on our claims"));
                TellTeam(other, (on ? TeamName(team) + " lets your team build and open things on their claims" : TeamName(team) + " no longer lets your team onto their claims"));
            }
            return true;
        }
        if (!b) { Tell(*me, "Stand on a body"); return true; }
        if (text == "/claim info") {
            const Claim* c = Claims_At(b->name, here);
            Tell(*me, c ? TextFormat("Claim #%u of %s (%.0f m across)%s", c->id, TeamName(c->team).c_str(), c->radius * 2.0f, Claims_Allows(c->team, team) ? ": you may build here" : "") : std::string("Unclaimed land"));
            return true;
        }
        if (!op) { Tell(*me, "Claim land with a claim banner (craft it, place it; break it to unclaim). /claims, /claim info, /claim ally TEAM, /claim unally TEAM"); return true; }
        if (text == "/unclaim") { // (ops) the claim they stand in, any team's
            const Claim* c = Claims_At(b->name, here);
            if (!c) Tell(*me, "No claim here");
            else { const int owner = c->team; const uint32_t id = c->id; Claims_Remove(id); TellTeam(owner, "An op removed claim #" + std::to_string(id)); Tell(*me, TextFormat("Removed claim #%u", id)); }
            return true;
        }
        if (text != "/claim") { Tell(*me, "/claim (ops: here, 100 m), /unclaim (ops), /claims, /claim info, /claim ally TEAM, /claim unally TEAM"); return true; }
        std::string why;
        const uint32_t id = Claims_Make(team, std::max(budget, 99), b->name, here, kClaimRadius, (int64_t)time(nullptr), why); // (an op's: over the budget)
        if (!id) { Tell(*me, why); return true; }
        ++claimsMade;
        TellTeam(team, TextFormat("%s claimed the land around them: claim #%u, %.0f m across", me->name.c_str(), id, kClaimRadius * 2.0f));
        claimInside.erase(session);
        return true;
    }
    // The claims list (NetMsg 41) to a client: every claim, with whether it is the client's team's and whether it may build
    // there; sent again when the claims change or the client's team does.
    std::unordered_map<uint32_t, std::pair<uint32_t, int>> claimsSent; // by connection: the version and team last sent
    int claimListsSent = 0;
    void SendClaims(const Session& s) {
        const int team = Teams_Of(s.uuid);
        auto& sent = claimsSent[s.conn];
        if (sent.first == Claims_Version() && sent.second == team) return;
        sent = { Claims_Version(), team };
        MsgClaims m;
        for (const Claim& c : Claims_All()) {
            Planetoid* b = nullptr;
            for (Planetoid* p : Planetoid::All()) if (p->name == c.body) { b = p; break; }
            if (!b) continue;
            NetClaim n; n.id = c.id; n.team = c.team; n.teamName = TeamName(c.team);
            const TeamRecord* t = Teams_Find(c.team);
            n.colour = Claims_TeamColour(c.team, t ? t->colour : 0);
            n.body = b->netId; n.centre[0] = c.centre.x; n.centre[1] = c.centre.y; n.centre[2] = c.centre.z; n.radius = c.radius;
            n.flags = (uint8_t)((c.team == team ? 1 : 0) | (Claims_Allows(c.team, team) ? 2 : 0));
            m.claims.push_back(std::move(n));
        }
        SendMsg(s.conn, m);
        ++claimListsSent;
    }
    // Once a second: decay (teams with a member online stay fresh), and a line to a player entering or leaving a claim.
    void StepClaims() {
        if (GetTime() < nextClaimStep) return;
        nextClaimStep = GetTime() + 1.0;
        std::vector<int> online;
        for (const Session& s : sessions) if (s.online) online.push_back(Teams_Of(s.uuid));
        if (const int gone = Claims_Update(online, (int64_t)time(nullptr))) printf("host: %d claim(s) decayed (their teams away %d days)\n", gone, (int)(kClaimDecaySeconds / 86400));
        for (const Session& s : sessions) {
            if (!s.online) continue;
            SendClaims(s); // (when they or its team changed)
            auto rit = remotes.find(s.id);
            Planetoid* b = rit != remotes.end() ? rit->second.player.nearest : nullptr;
            const Claim* c = b ? Claims_At(b->name, b->toLocalD(rit->second.player.position).toFloat()) : nullptr;
            const uint32_t now = c ? c->id : 0;
            auto& was = claimInside[s.id];
            if (now == was) continue;
            if (now) Tell(s, c->team == Teams_Of(s.uuid) ? std::string("Entering your team's claim") : Claims_Allows(c->team, Teams_Of(s.uuid)) ? "Entering " + TeamName(c->team) + "'s claim (allies: you may build here)" : "Entering " + TeamName(c->team) + "'s claim: you can't build, dig or open anything here");
            else if (was) Tell(s, "Leaving the claim");
            was = now;
        }
    }
    std::string CleanChat(std::string t) {
        for (char& ch : t) if ((unsigned char)ch < 0x20 || (unsigned char)ch == 0x7F) ch = ' ';
        return t;
    }
    // (#349) connections: at most kMaxPerAddress from one address at once (loopback exempt: the host's own machine and the
    // tests), and a connection must say Hello within kHelloSeconds.
    const int kMaxPerAddress = getenv("RAYLIBGAME_MAX_PER_ADDRESS") ? std::max(1, atoi(getenv("RAYLIBGAME_MAX_PER_ADDRESS"))) : 8;
    // (#1868) and at most kMaxPending connections that haven't joined yet (from any non-loopback address), counted from the
    // moment they connect (a silent one used to count only after its first message); the tests shorten the Hello wait
    const int kMaxPending = getenv("RAYLIBGAME_MAX_PENDING") ? std::max(1, atoi(getenv("RAYLIBGAME_MAX_PENDING"))) : 16;
    const double kHelloSeconds = getenv("RAYLIBGAME_HELLO_SECONDS") ? std::max(0.5, atof(getenv("RAYLIBGAME_HELLO_SECONDS"))) : 10.0;
    int pendingRefused = 0, silentDropped = 0; // (#1868) tests
    std::string AddressHost(const std::string& a) { const size_t c = a.rfind(':'); return c == std::string::npos ? a : a.substr(0, c); }
    // (loopback is counted too when RAYLIBGAME_MAX_PER_ADDRESS is set: the tests check the limit on one machine)
    bool Loopback(const std::string& host) { if (getenv("RAYLIBGAME_MAX_PER_ADDRESS")) return false; return host.empty() || host == "127.0.0.1" || host == "::1" || host == "[::1]" || host.rfind("[::ffff:127.", 0) == 0; }
    void OnHello(Scene& world, uint32_t conn, const MsgHello& h) {
        if (h.protocol != kNetProtocol) { Refuse(conn, NetReject::Protocol, TextFormat("server %u, client %u", kNetProtocol, h.protocol)); return; }
        if (!h.build.empty() && h.build != BuildCommit()) { Refuse(conn, NetReject::Build, "server " + BuildCommit() + ", client " + h.build); return; }
        if (!password.empty() && h.password != password) { Refuse(conn, NetReject::Password, ""); return; }
        if (joinCheck) { // (#345) bans and the whitelist
            std::string why;
            if (!joinCheck(h.name, Net_RemoteAddress(conn), h.uuid, why)) {
                Refuse(conn, NetReject::Banned, why);
                if (eventHook) eventHook("refused", 0, h.name + " (" + why + ")");
                return;
            }
        }
        Session* s = nullptr;
        // (#349) a name is 1-24 bytes, no control characters, no leading or trailing space (it goes into logs, the audit log,
        // nameplates and chat: a newline in it would forge lines). Names are one identity whatever their case (the admin lists
        // compare without case: "ALICE" got alice's op role in a session of its own).
        if (!NameAcceptable(h.name)) { Refuse(conn, NetReject::Malformed, "names are 1-24 characters, no control characters or edge spaces"); return; }
        // (#346) a player is their id; a name is theirs on this server once they've played under it: another id asking for
        // it is refused (not renamed), while the same id may come back under a new name.
        if (!PlayerId_Valid(h.uuid)) { Refuse(conn, NetReject::Malformed, "no player id"); return; }
        if (NameClaimedByOther(h.name, h.uuid)) { Refuse(conn, NetReject::NameTaken, h.name + " belongs to another player on this server"); return; }
        for (Session& x : sessions) if (x.uuid == h.uuid) s = &x;
        if (s && s->online) { Refuse(conn, NetReject::NameTaken, h.name + " is already playing"); return; }
        const int online = (int)std::count_if(sessions.begin(), sessions.end(), [](const Session& x) { return x.online; });
        if (online >= maxPlayers) { Refuse(conn, NetReject::Full, TextFormat("%d players", maxPlayers)); return; }
        const bool reconnected = s != nullptr;
        if (!s) { sessions.push_back({ nextSession++, h.name, 0, false }); s = &sessions.back(); s->uuid = h.uuid; }
        s->name = h.name; // (a new name of theirs)
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
        bool hadRecord = true; // (a reconnect within the run keeps its player: not a new one)
        if (!remotes.count(s->id)) { // (#346) at the world's start, or where their record left them, with what they carried
            SpawnRemote(world, s->id);
            hadRecord = LoadRecord(world, conn, *s);
        } else { // (#780) a reconnect keeps the player but its client counts ticks from 1 again: a fresh input stream
            Remote& rm = remotes[s->id];
            rm.inputs.clear(); rm.lastQueued = 0; rm.lastApplied = 0; rm.carried = 0; rm.budget = 0.0; rm.tick = SimTick();
            rm.farmHash = 0; // (#335) its fields again
        }
        // Who is here (everyone, the newcomer included), and what the others hold.
        SendRoster();
        for (const Session& x : sessions) if (x.online && x.id != s->id && x.loadout.session) { MsgLoadout l = x.loadout; l.held = HeldItem(x.id); SendMsg(conn, l); } // (#1868) as had
        StartSync(world, conn, s->id); // the world's edits so far (#776); later ones follow as live edits
        SendGraves(conn); // (#804) the graves so far
        SendAllPlants(world, conn); // (#335) the plants felled so far
        claimsSent.erase(conn); // (#343) the claims list goes with the next claims step
        ++joins;
        HostLog("%s joined as session %u%s (connection %u, %d joins so far)", h.name.c_str(), s->id, reconnected ? " (reconnected)" : "", conn, joins);
        if (eventHook) eventHook("join", s->id, h.name + " from " + Net_RemoteAddress(conn)); // (#345 audit)
        if (!hadRecord) { // (#804) a new player's starting kit, as single player's (the client's own inventory follows this one); (#346) never again for one with a record
            Host_Give(s->id, "telescope", 1);
            if (worldDifficulty < 2) Host_Give(s->id, "map", 1); // hard: find your own way
            if (worldDifficulty == 0) Host_Give(s->id, "berries", 6);
        }
        if (const char* d = getenv("RAYLIBGAME_HOST_DIFFICULTY")) worldDifficulty = std::clamp(atoi(d), 0, 2); // tests (#338): e.g. peaceful, no night threats
        if (const char* hunger = getenv("RAYLIBGAME_HOST_HUNGER")) remotes[s->id].survival.food = (float)atof(hunger); // tests (#337): a hungry joiner
        if (const char* give = getenv("RAYLIBGAME_HOST_GIVE")) { // tests (#337): "item count[,item count...]" to each player as they join
            for (const char* g = give; g && *g; ) {
                char item[64] = ""; int n = 0;
                if (sscanf(g, " %63[^ ,] %d", item, &n) == 2 && n > 0) printf("host: test give %s %d to session %u: %s\n", item, n, s->id, Host_Give(s->id, item, n) ? "given" : "unknown item");
                g = strchr(g, ','); if (g) ++g;
            }
        }
        // (#338) tests: RAYLIBGAME_HOST_TEAM="alice,bob": those players share one team (the first one's; uuids from the names)
        if (const char* team = getenv("RAYLIBGAME_HOST_TEAM")) {
            std::vector<std::string> names; std::string cur;
            for (const char* c = team; ; ++c) { if (*c == ',' || !*c) { if (!cur.empty()) names.push_back(cur); cur.clear(); if (!*c) break; } else cur += *c; }
            if (!names.empty() && std::find(names.begin(), names.end(), s->name) != names.end()) {
                const std::string lead = PlayerId_FromName(names[0]);
                const int t = Teams_Join(lead, names[0], 0);
                if (s->uuid != lead) Teams_Move(s->uuid, t);
                printf("host: test team: %s in team %d\n", s->name.c_str(), Teams_Of(s->uuid));
            }
        }
        // (#338) what this player obtains from now on counts for their team's book (the kit and a loaded record don't)
        {
            const uint32_t sid = s->id;
            remotes[sid].inventory.onAdded = [sid](const std::string& id, double n) { TeamEvent(sid, "obtained", id, n); };
        }
        if (const int t = Teams_Of(s->uuid)) SendQuestBook(conn, t); // (#338) their team's book
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
        if (ClaimBlocks(session, p, hit)) { RefuseEdit(conn, q.seq, NetEditRefusal::Claim); return; } // (#343)
        // The tool in the player's hands (its loadout, when its inventory here has one: #1868) and what it does to the ground.
        const std::string held = HeldItem(session);
        const ToolSpec tool = ToolOf(held);
        const Vector3 surfN = p->field->normal(hit, p->edits.get());
        const VoxelMaterial mat = p->materialAt(hit, surfN);
        const DigRule rule = DigRuleFor(tool, mat, TerrainTool_ItemAt(*p, hit, (int)mat));
        // (#1043) The player's own perks scale the stroke as their client does (Mining on rock and ore, Digging on soil).
        const auto stat = StatsOf(rm);
        const DigPerkScale perks = TerrainTool_DigPerkScale(tool.type == ToolType::Pickaxe, (int)mat, stat);
        const float baseRadius = (rule.radius > 0.0f ? rule.radius : DigRuleFor(tool, VoxelMaterial::Regolith, nullptr).radius) * perks.radius / p->scale;
        const bool dig = q.kind == NetEditKind::Dig, till = q.kind == NetEditKind::Till;
        if (dig && !rule.can) { RefuseEdit(conn, q.seq, NetEditRefusal::Tool); return; }
        // (#1085) a till: a shovel or hoe in hand, soil under the aim, the tool's brush at most, centred on the aim
        if (till && (TillRadius(tool) <= 0.0f || !IsTillable(mat))) { RefuseEdit(conn, q.seq, NetEditRefusal::Tool); return; }
        // The stroke's shape: the tool's sphere (place: 1.3x, the soft blob), its centre where the tool puts it.
        const float maxRadius = till ? TillRadius(tool) / p->scale * 1.01f + 1e-4f : baseRadius * (dig ? 1.0f : 1.3f) * 1.01f + 1e-4f;
        const float offset = Vector3Distance(hit, centre);
        if (q.radius > maxRadius || offset > (till ? 0.0f : q.radius * (dig ? 0.5f : 0.3f)) + 0.05f) { RefuseEdit(conn, q.seq, NetEditRefusal::Shape); return; }
        // Rate: a token per cooldown (the tool's strokes a second, at least 2, with the perks), up to a small burst.
        const double rate = std::max(rule.strokesPerSecond, 2.0f) * perks.rate;
        if (rm.editTokensAt < 0.0) rm.editTokensAt = SimTime();
        rm.editTokens = std::min(kEditBurst, rm.editTokens + (SimTime() - rm.editTokensAt) * rate);
        rm.editTokensAt = SimTime();
        if (rm.editTokens < 0.999) { // (#1059) the numbers behind a rate refusal, for the log
            printf("host: edit %u rate: %.2f tokens at sim %.2f s, wall %.2f s (%.2f strokes/s, tool '%s')\n", q.seq, rm.editTokens, SimTime(), GetTime(), rate, held.c_str());
            RefuseEdit(conn, q.seq, NetEditRefusal::Rate); return;
        }
        MsgTerrainEdit out;
        out.body = q.body; out.center = q.center; out.radius = q.radius; out.session = session; out.seq = q.seq;
        TerrainEdit e{};
        e.center = centre; e.radius = q.radius;
        std::map<std::string, double> before;
        if (dig) {
            const bool softOnly = tool.type != ToolType::Pickaxe;
            DigYield got = TerrainTool_DigYield(*p, centre, q.radius, hit, softOnly, stat("gem_chance")); // (#920) the player's perks
            TerrainTool_ScaleYield(got, stat("ore_yield"), stat("stone_yield"));
            for (const auto& [item, amt] : got.items) { before.emplace(item, rm.inventory.amount(item)); }
            if (got.gem) before.emplace(got.gem, rm.inventory.amount(got.gem));
            if (got.seed) before.emplace(got.seed, rm.inventory.amount(got.seed));
            for (const auto& [item, amt] : got.items) rm.inventory.addFractional(item, amt);
            if (got.gem) rm.inventory.add(got.gem, 1);
            if (got.seed) rm.inventory.add(got.seed, 1);
            e.carve = true; e.softOnly = softOnly;
            out.flags = 1 | (softOnly ? 8 : 0);
            if (!rm.creative) { // (#920) Mining or Digging XP for the stroke, as a game run here counts it (#883); ground the
                // player filled and dug out again soon pays little (the anti-farming share)
                bool rock = false;
                for (const auto& [item, amt] : got.items) rock = rock || ItemMatches("any_stone", item);
                const char* ev = got.ore ? "ore_stroke" : (!softOnly && rock) ? "rock_stroke" : "soil_stroke";
                const double share = Skills_DigShare(rm.skills, SkillCell(hit), SimTime());
                int ups = Skills_Event(rm.skills, ev, share, worldXpRate);
                if (got.gem) ups += Skills_Event(rm.skills, "gem", 1.0, worldXpRate);
                rm.skillsDirty = true;
                if (ups > 0) rm.skillsSentAt = -1e9; // (a level-up goes out at once)
                TeamXp(session, ev, share);
                if (got.gem) TeamXp(session, "gem", 1.0);
            }
        } else if (till) { // (#1085) the soil there becomes farmland: no cost, no yield
            e.till = true;
            out.flags = 16;
            TeamEvent(session, "count", "till"); RemoteXp(session, "till"); // (#335) the quest's count, Farming XP
            printf("host: session %u tilled (edit %u, radius %.2f m)\n", session, q.seq, q.radius); // (#1107) an accepted till, to trace it
        } else {
            const int paint = TerrainTool_FillerMaterial(q.item);
            if (paint < 0) { RefuseEdit(conn, q.seq, NetEditRefusal::Shape); return; }
            before.emplace(q.item, rm.inventory.amount(q.item));
            e.strength = 0.6f; e.paint = paint;
            const float cost = TerrainTool_PlaceCost(*p, e); // (#1097) by the ground it adds, as the client reckons it
            if (cost < 0.05f) { RefuseEdit(conn, q.seq, NetEditRefusal::Shape); return; }
            if (!rm.creative && !rm.inventory.removeAmount(q.item, cost)) { RefuseEdit(conn, q.seq, NetEditRefusal::Cost); return; } // (#823) free in creative
            for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) // (#920) our ground now
                Skills_NoteFill(rm.skills, SkillCell(Vector3Add(hit, Vector3{ (float)dx, (float)dy, (float)dz })), SimTime());
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
            if (PartUnbreakable(plot->parts[q.partIndex])) { RefusePart(conn, q.seq, NetPartRefusal::Blocked); return; } // (#1641) a vault door, a portcullis or a lever is opened or pulled, never broken
            const Part removed = plot->parts[q.partIndex];
            const Vector3 at = plot->toBody({ (removed.origin.x + 0.5f) * kCellSize, (removed.origin.y + 0.5f) * kCellSize, (removed.origin.z + 0.5f) * kCellSize });
            if (!InBuildReach(world, rm, q.body, at, p->scale)) { RefusePart(conn, q.seq, NetPartRefusal::Reach); return; }
            if (ClaimBlocks(session, p, at)) { RefusePart(conn, q.seq, NetPartRefusal::Claim); return; } // (#343)
            // (#343) a claim banner: only its team (or an op) may break it, and breaking it unclaims
            const uint32_t bannerClaim = ClaimBannerTier(removed.type) >= 0 ? Claims_OfBanner(p->name, plot->id, (int)q.partIndex) : 0;
            int bannerTeam = 0;
            for (const Claim& c : Claims_All()) if (c.id == bannerClaim) bannerTeam = c.team;
            if (bannerClaim && bannerTeam != TeamOfSession(session) && !IsOpSession(session)) { RefusePart(conn, q.seq, NetPartRefusal::Claim); return; }
            if (IsMachine(plot->parts[q.partIndex].type)) plot->refundTo = session; // (#854) its contents go to this player (the factory's next layout)
            if (!plot->remove((int)q.partIndex)) { RefusePart(conn, q.seq, NetPartRefusal::NoPart); return; }
            // The refund (a meteorite breaks into ore, rolled here: the authority rolls, #316).
            std::vector<std::pair<std::string, int>> give;
            if (removed.type == PartType::Meteorite) {
                give.push_back({ "nickel_iron_ore", GameRng_Int(6, 10) });
                if (GameRng_Int(0, 99) < 35) give.push_back({ "platinum_ore", 2 });
                if (GameRng_Int(0, 99) < 12) give.push_back({ "iridium_ore", 1 });
            } else {
                give.push_back({ PartItem(removed), PartCostCount(removed.type) });
                Structures_BreakLoot(*p, q.plot, (int)q.partIndex, removed, give); // (#1641 slice 2c) a dungeon urn or crate: its loot too
            }
            for (const auto& [item, n] : give) { before.emplace(item, rm.inventory.amount(item)); rm.inventory.add(item, n); }
            MsgPlotEdit out; out.session = session; out.seq = q.seq; out.op = NetPartOp::Remove; out.body = q.body; out.plot = q.plot; out.partIndex = q.partIndex;
            SendToJoined(out);
            if (bannerClaim) { Claims_Remove(bannerClaim); TellTeam(bannerTeam, TextFormat("Claim #%u is gone: its banner was broken", bannerClaim)); }
        } else {
            const Part part = ToPart(q.part);
            if (q.part.type >= (uint8_t)PartType::Count || part.type == PartType::Meteorite || part.type == PartType::PressurePlate || part.type == PartType::GasVent || part.type == PartType::Lever || (IsShape(part.type) && part.material >= BuildMaterialCount())) {
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
            if (ClaimBlocks(session, p, at)) { RefusePart(conn, q.seq, NetPartRefusal::Claim); return; } // (#343)
            const char* item = PartItem(part);
            if (part.type == PartType::Farmland) { RefusePart(conn, q.seq, NetPartRefusal::Blocked); return; } // (#1085) farmland is tilled ground now (a till edit), not a part
            const int cost = PartCostCount(part.type);
            if (!rm.creative && rm.inventory.count(item) < cost) { RefusePart(conn, q.seq, NetPartRefusal::Cost); return; } // (#823) free in creative
            uint32_t bannerClaim = 0; // (#343) a claim banner claims the land around it for the placer's team (or is refused: budget, overlap)
            if (const int tier = ClaimBannerTier(part.type); tier >= 0) {
                const int team = TeamOfSession(session);
                std::string why;
                bannerClaim = Claims_Make(team, Claims_Budget(TeamQuests_BodiesReached(team)), p->name, at, kClaimBannerAcross[tier] * 0.5f, (int64_t)time(nullptr), why);
                if (!bannerClaim) { if (const uint32_t c = ConnOf(session)) SendMsg(c, MsgChat{ 0, why }); RefusePart(conn, q.seq, NetPartRefusal::Claim); return; }
            }
            if (!plot) { // a new plot, made here first (the clients make the same one from the same id and anchor)
                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                plot = p->plots->plotFor(anchor, 0.0f, true);
                MsgPlotEdit made; made.session = session; made.seq = q.seq; made.op = NetPartOp::NewPlot; made.body = q.body; made.plot = plot->id;
                made.anchor[0] = anchor.x; made.anchor[1] = anchor.y; made.anchor[2] = anchor.z;
                SendToJoined(made);
            }
            before.emplace(item, rm.inventory.amount(item));
            if (!rm.creative) rm.inventory.remove(item, cost);
            BuildTool_ApplyPart(*p, *plot, part);
            if (bannerClaim) {
                Claims_SetBanner(bannerClaim, plot->id, plot->partAt(part.origin));
                const int team = TeamOfSession(session);
                for (const Session& s : sessions) if (s.id == session) TellTeam(team, TextFormat("%s raised a claim banner: claim #%u, %.0f m across", s.name.c_str(), bannerClaim, kClaimBannerAcross[ClaimBannerTier(part.type)]));
            }
            if (IsMachine(part.type) || IsRocketPart(part.type)) TeamEvent(session, "placed", PartCostItem(part.type)); // (#338) "place" tasks
            MsgPlotEdit out; out.session = session; out.seq = q.seq; out.op = NetPartOp::Place; out.body = q.body; out.plot = plot->id; out.part = q.part;
            SendToJoined(out);
        }
        rm.partTokens -= 1.0;
        SendInventoryDeltas(conn, rm.inventory, before);
        ++partsApplied;
    }

    // ---- Join-in-progress sync (#776) ----
    // At a join: every body's edited regions (encoded now from the immutable delta snapshot: later edits reach the
    // client as live messages after this) and its plots, hashed, in an index the client answers with what it lacks.
    // Those then stream nearest the player first, a byte budget a frame (RAYLIBGAME_SYNC_KBPS, default 1024 kB/s).
    struct SyncItem { NetSyncEntry e; std::vector<uint8_t> bytes; bool packed = false; double dist = 0.0; }; // bytes: plain, then DEFLATEd when its turn comes
    struct SyncJob {
        uint32_t snapshot = 0, session = 0;
        std::vector<SyncItem> items;
        bool wantIn = false;
        std::vector<uint32_t> order; size_t cur = 0; uint32_t offset = 0; // what to send, the one being sent, how far
        uint64_t sent = 0; double started = 0.0, budget = 0.0, lastAt = 0.0;
        int wanted = 0;
        double statusAt = 0.0; // (#1836) when its transport queue was last printed
        std::future<std::vector<SyncItem>> regions; // the worker encoding and hashing the regions
        bool indexed = false;
    };
    std::unordered_map<uint32_t, SyncJob> syncs;   // by connection
    constexpr uint32_t kSyncMaxPending = 1u << 20; // bytes in flight on a connection before the sync waits
    void ForgetSync(uint32_t conn) { syncs.erase(conn); }
    uint32_t nextSnapshot = 1;
    int syncsDone = 0;
    uint64_t syncBytesSent = 0;

    double SyncRate() { // bytes a second
        static double rate = [] { const char* s = getenv("RAYLIBGAME_SYNC_KBPS"); const double k = s ? atof(s) : 1024.0; return std::max(16.0, k) * 1024.0; }();
        return rate;
    }
    // The hash and size are of an entry as stored (a client compares its own). The bytes go DEFLATE-compressed (RGN5's
    // run-length chunks shrink ~4x more), each when its turn comes: compressing them all at the join held the server's
    // frame ~0.9 s (3000 strokes). Encoding and hashing the regions for the index runs on a worker (it held the frame
    // ~130 ms): the delta snapshots are immutable, so the thread reads them safely while play goes on.
    SyncItem MakeItem(NetSyncKind kind, uint64_t body, uint64_t key, std::vector<uint8_t>&& plain) {
        SyncItem it; it.e.kind = kind; it.e.body = body; it.e.key = key;
        it.e.hash = Region_Hash(plain); it.e.size = (uint32_t)plain.size();
        it.bytes = std::move(plain);
        return it;
    }
    void StartSync(Scene& world, uint32_t conn, uint32_t session) {
        SyncJob job;
        job.snapshot = nextSnapshot++; job.session = session; job.started = job.lastAt = GetTime();
        std::vector<std::pair<uint64_t, std::shared_ptr<const TerrainDeltas>>> snaps;
        for (auto& o : world.objects) {
            Planetoid* p = dynamic_cast<Planetoid*>(o.get());
            if (!p || !p->netId) continue;
            if (p->plots && p->plots->partCount() > 0) job.items.push_back(MakeItem(NetSyncKind::Plots, p->netId, 0, Plots_Blob(*p->plots))); // (small: here)
            if (p->edits && !p->edits->empty()) snaps.push_back({ p->netId, p->edits }); // immutable: later edits make a new one
        }
        job.regions = std::async(std::launch::async, [snaps]() {
            std::vector<SyncItem> out;
            for (const auto& [body, snap] : snaps)
                for (TerrainDeltas::Key k : Regions_List(*snap)) out.push_back(MakeItem(NetSyncKind::Region, body, k, Region_Encode(*snap, k)));
            return out;
        });
        syncs[conn] = std::move(job);
    }
    // The worker's regions are in: the index goes out, in pages.
    void SendIndex(uint32_t conn, SyncJob& job) {
        std::vector<SyncItem> regions = job.regions.get();
        for (SyncItem& r : regions) job.items.push_back(std::move(r));
        job.indexed = true;
        const size_t n = job.items.size();
        size_t raw = 0;
        for (const SyncItem& s : job.items) raw += s.bytes.size();
        const uint16_t pages = (uint16_t)std::max<size_t>(1, (n + kSyncIndexPage - 1) / kSyncIndexPage);
        for (uint16_t pg = 0; pg < pages; ++pg) {
            MsgSyncIndex x; x.snapshot = job.snapshot; x.page = pg; x.pages = pages; x.total = (uint32_t)n;
            for (size_t i = (size_t)pg * kSyncIndexPage; i < n && i < (size_t)(pg + 1) * kSyncIndexPage; ++i) x.entries.push_back(job.items[i].e);
            SendMsg(conn, x);
        }
        HostLog("join sync %u for session %u: %zu entries, %.1f MB as stored, indexed in %.0f ms (on a worker)", job.snapshot, job.session, n, raw / 1048576.0, (GetTime() - job.started) * 1000.0);
    }
    void OnSyncWant(Scene& world, uint32_t conn, const MsgSyncWant& w) {
        auto it = syncs.find(conn);
        if (it == syncs.end() || it->second.wantIn || w.snapshot != it->second.snapshot || w.total != it->second.items.size()) return;
        SyncJob& job = it->second;
        job.wantIn = true;
        // Nearest the player first (plots before regions: small, and what players care about most).
        Remote& rm = remotes[job.session];
        Player::NetState ns;
        const bool have = rm.player.netState(world, ns);
        const double side = kRegionChunks * kDeltaChunk * kVoxelSize;
        for (uint32_t i = 0; i < job.items.size(); ++i) {
            if (!(w.bits[i / 8] & (1u << (i % 8)))) continue;
            SyncItem& s = job.items[i];
            if (s.e.kind == NetSyncKind::Plots) s.dist = -1.0;
            else {
                int rx, ry, rz;
                TerrainDeltas::Unpack(s.e.key, rx, ry, rz);
                const DVec3 c{ (rx + 0.5) * side, (ry + 0.5) * side, (rz + 0.5) * side };
                s.dist = have && ns.body == s.e.body ? (c - ns.local).length() : 1e12;
            }
            job.order.push_back(i);
        }
        std::stable_sort(job.order.begin(), job.order.end(), [&](uint32_t a, uint32_t b) { return job.items[a].dist < job.items[b].dist; });
        job.wanted = (int)job.order.size();
        job.lastAt = GetTime();
        // Entries not wanted: free them now.
        std::vector<bool> keep(job.items.size(), false);
        for (uint32_t i : job.order) keep[i] = true;
        for (size_t i = 0; i < job.items.size(); ++i) if (!keep[i]) { job.items[i].bytes.clear(); job.items[i].bytes.shrink_to_fit(); }
    }
    void StepSyncs(Scene& world) {
        for (auto it = syncs.begin(); it != syncs.end();) {
            SyncJob& job = it->second;
            if (!job.indexed) { // the worker still encoding? (else: the index goes out now)
                if (job.regions.wait_for(std::chrono::seconds(0)) == std::future_status::ready) SendIndex(it->first, job);
                ++it; continue;
            }
            if (!job.wantIn) { ++it; continue; }
            if (GetTime() - job.statusAt >= 1.0) { // (#1836, #1421) the joining link's transport queue, once a second, to see a join stall
                job.statusAt = GetTime();
                NetQueueStatus q;
                if (Net_QueueStatus(it->first, q))
                    HostLog("join sync %u: %zu of %zu entries, %.2f MB sent; queue %d B reliable, %d B unreliable, %d B unacked, ping %d ms; %d refused sends",
                           job.snapshot, job.cur, job.order.size(), job.sent / 1048576.0, q.pendingReliable, q.pendingUnreliable, q.sentUnacked, q.pingMs, Net_SendFailures());
            }
            if (Net_PendingBytes(it->first) > kSyncMaxPending) { job.lastAt = GetTime(); ++it; continue; } // the link is behind: wait
            const double now = GetTime();
            job.budget = std::min(job.budget + (now - job.lastAt) * SyncRate(), SyncRate() * 0.25); // (a quarter second's worth at most)
            job.lastAt = now;
            while (job.cur < job.order.size() && job.budget > 0.0) {
                SyncItem& s = job.items[job.order[job.cur]];
                if (!s.packed) { // its turn: compress it now (a few ms for a busy region)
                    int z = 0;
                    unsigned char* c = s.bytes.empty() ? nullptr : CompressData(s.bytes.data(), (int)s.bytes.size(), &z);
                    if (c) { s.bytes.assign(c, c + z); MemFree(c); }
                    s.packed = true;
                }
                const uint32_t n = (uint32_t)std::min<size_t>(kSyncPiece, s.bytes.size() - job.offset);
                MsgSyncData d; d.snapshot = job.snapshot; d.entry = job.order[job.cur]; d.offset = job.offset; d.size = (uint32_t)s.bytes.size();
                d.bytes.assign(s.bytes.begin() + job.offset, s.bytes.begin() + job.offset + n);
                SendMsg(it->first, d);
                job.offset += n; job.sent += n; job.budget -= n; syncBytesSent += n;
                if (job.offset >= s.bytes.size()) { s.bytes.clear(); s.bytes.shrink_to_fit(); ++job.cur; job.offset = 0; }
            }
            if (job.cur >= job.order.size()) {
                SendMsg(it->first, MsgSyncDone{ job.snapshot, (uint32_t)job.wanted, job.sent });
                HostLog("join sync %u done: %d of %zu entries wanted, %.2f MB in %.2f s", job.snapshot, job.wanted, job.items.size(), job.sent / 1048576.0, GetTime() - job.started);
                if (getenv("RAYLIBGAME_SYNC_FINGERPRINT")) // (tests: a client's after its sync must equal this; costs an encode of every region)
                    HostLog("terrain fingerprint %016llx", (unsigned long long)NetSync_TerrainFingerprint(world));
                ++syncsDone;
                it = syncs.erase(it);
            } else ++it;
        }
    }

    // (#349) the guard's verdicts: a kick (with why, in the log and the audit hook), and violations found after reading.
    int guardDropped = 0, guardKicks = 0;
    void GuardKick(uint32_t conn, const std::string& why) {
        auto it = peers.find(conn);
        std::string name; if (it != peers.end()) for (const Session& s : sessions) if (s.id == it->second.session) name = s.name;
        HostLog("kicked connection %u (%s) for flooding: %s", conn, name.empty() ? "not joined" : name.c_str(), why.c_str());
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
    // (#349, #1868) A new connection's admission: 0 admit, 1 too many from its address, 2 too many not joined yet. Loopback
    // (the host's machine, the tests) is exempt.
    int AdmitVerdict(bool loopback, int sameAddress, int pending) {
        if (loopback) return 0;
        if (sameAddress >= kMaxPerAddress) return 1;
        return pending >= kMaxPending ? 2 : 0;
    }
    // Admit a new connection as a peer (its Hello clock starts now), or refuse it. Called when it connects (#1868) and, for a
    // connection the transport didn't report, at its first message.
    bool AdmitPeer(uint32_t conn) {
        const std::string host = AddressHost(Net_RemoteAddress(conn));
        const bool loop = Loopback(host);
        int same = 0, pending = 0;
        if (!loop)
            for (const auto& [c, p] : peers) {
                const std::string h = AddressHost(Net_RemoteAddress(c));
                if (h == host) ++same;
                if (!p.joined && !Loopback(h)) ++pending;
            }
        if (const int v = AdmitVerdict(loop, same, pending)) {
            const char* why = v == 1 ? "too many connections from your address" : "too many players connecting at once: try again shortly";
            HostLog("refused connection %u from %s: %s (%d from that address, %d not joined)", conn, host.c_str(), why, same, pending);
            if (eventHook) eventHook("refused", 0, "connection from " + host + (v == 1 ? " (too many from one address)" : " (too many connecting)"));
            Refuse(conn, NetReject::Full, why);
            ++pendingRefused;
            return false;
        }
        Peer& p = peers[conn];
        p.firstHeard = p.lastHeard = GetTime();
        HostLog("connection %u from %s", conn, host.empty() ? "?" : host.c_str());
        return true;
    }
    void OnMessage(Scene& world, const NetMessage& m) {
        if (!peers.count(m.conn) && !Net_IsConnected(m.conn)) return; // (#349: the rest of a batch from a connection just kicked)
        auto it = peers.find(m.conn);
        if (it == peers.end()) { // first word on a connection not admitted yet
            if (!AdmitPeer(m.conn)) return;
            it = peers.find(m.conn);
        }
        it->second.lastHeard = GetTime();
        if (NetPeek(m.data) == NetMsg::Input) it->second.inGame = true; // (#1033)
        { // (#349) size, kind, form and rate first: a bad message is dropped; a client that keeps sending them is kicked
            std::string why;
            const GuardVerdict v = NetGuard_Admit(it->second.guard, NetPeek(m.data), m.data.size(), NetWellFormed(m.data), GetTime(), why);
            if (v != GuardVerdict::Accept) {
                // (#1107) the first drops are logged: a player's request lost here was silent ("can't interact")
                if (guardDropped < 20) printf("host: dropped a %s from session %u (%zu bytes): %s\n", NetMsgName(NetPeek(m.data)), it->second.session, m.data.size(), why.c_str());
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
            case NetMsg::UseItem: { // (#337) eat or use an item: from the player's inventory as the server keeps it, into its vitals
                MsgUseItem u;
                if (!it->second.joined || !NetRead(r, u)) return;
                auto rit = remotes.find(it->second.session);
                if (rit == remotes.end() || rit->second.survival.dead) return;
                std::map<std::string, double> before{ { u.item, rit->second.inventory.amount(u.item) } };
                if (Survival_UseItemFor(u.item, rit->second.survival, rit->second.inventory, rit->second.suitStats, false)) {
                    SendInventoryDeltas(m.conn, rit->second.inventory, before);
                    if (Survival_FoodValue(u.item) > 0.0f) { TeamEvent(it->second.session, "count", "eat"); RemoteXp(it->second.session, "eat"); } // (#338) "do eat" tasks
                }
                else SendMsg(m.conn, MsgChat{ 0, "Can't use " + u.item + " now" });
                ++itemsUsed;
                return;
            }
            case NetMsg::MachineAction: { // (#332) insert, extract, take a slot, fuel; watch / unwatch its screen
                MsgMachineAction a;
                if (!it->second.joined || !NetRead(r, a)) return;
                OnMachineAction(m.conn, it->second.session, a);
                return;
            }
            case NetMsg::PlantAction: { // (#335) fell a plant, sow or harvest farmland
                MsgPlantAction a;
                if (!it->second.joined || !NetRead(r, a)) return;
                OnPlantAction(m.conn, it->second.session, a);
                return;
            }
            case NetMsg::DropItems: { // (#865) some of an item from the inventory this server keeps, as a pile at the player's feet
                MsgDropItems q;
                if (!it->second.joined || !NetRead(r, q)) return;
                auto rit = remotes.find(it->second.session);
                if (rit == remotes.end() || rit->second.survival.dead) return;
                Inventory& inv = rit->second.inventory;
                const int64_t have = (int64_t)llround(inv.amount(q.item) * 100.0);
                const int64_t centi = std::min(q.centi, have);
                if (centi <= 0) return; // (it has none: the client's view was stale; nothing to send)
                const auto before = AllAmounts(inv);
                if (!inv.removeAmount(q.item, centi / 100.0)) return;
                std::string name;
                for (const Session& s : sessions) if (s.id == it->second.session) name = s.name;
                Player::NetState feet; // (its body and body-local feet as the last step left them: right even before a first step, #856)
                const Planetoid* body = rit->second.player.netState(world, feet) ? Planetoid::findNet(feet.body) : nullptr;
                if (body && Combat_DropFor(body->name, feet.local, InvSlot{ q.item, centi }, name)) ++dropsMade;
                else inv.addAmount(q.item, centi / 100.0); // (not on a body: it stays)
                SendChangedItems(m.conn, inv, before);
                return;
            }
            case NetMsg::SkillSpend: { // (#920) a point on a perk: the server's copy of the skills decides (points, tier, ranks, needs)
                MsgSkillSpend q;
                if (!it->second.joined || !NetRead(r, q)) return;
                auto rit = remotes.find(it->second.session);
                if (rit == remotes.end()) return;
                if (!Skills_Buy(rit->second.skills, q.perk)) {
                    ++skillSpendsRefused;
                    const PerkDef* p = Skills_FindPerk(q.perk);
                    SendMsg(m.conn, MsgChat{ 0, std::string("Can't learn ") + (p ? p->name : q.perk) + " now" });
                    return;
                }
                ++skillSpends;
                rit->second.skillsDirty = true; rit->second.skillsSentAt = -1e9; // (the new state at once)
                return;
            }
            case NetMsg::CreativeGive: { // (#866) a creative player's Catalog: free items (or "" clears its inventory); refused otherwise
                MsgCreativeGive q;
                if (!it->second.joined || !NetRead(r, q)) return;
                auto rit = remotes.find(it->second.session);
                if (rit == remotes.end()) return;
                if (!rit->second.creative) { ++creativeRefused; SendMsg(m.conn, MsgChat{ 0, "Only creative players take from the Catalog" }); return; }
                Inventory& inv = rit->second.inventory;
                const auto before = AllAmounts(inv);
                if (q.item.empty()) inv.clear();
                else {
                    const ItemDef* d = itemDB.find(q.item);
                    if (!d || IsWildcardItem(q.item)) return;
                    inv.addAmount(q.item, (double)std::min<uint32_t>(q.count, (uint32_t)StackLimit(q.item)));
                }
                ++creativeGives;
                SendChangedItems(m.conn, inv, before);
                return;
            }
            case NetMsg::CraftRequest: { // (#338) crafting on the server: the recipe unlocked for the player's team, the inputs in its inventory
                MsgCraftRequest q;
                if (!it->second.joined || !NetRead(r, q)) return;
                auto rit = remotes.find(it->second.session);
                if (rit == remotes.end() || rit->second.survival.dead) return;
                const auto& recipes = itemDB.recipes();
                if (q.recipe >= recipes.size() || recipes[q.recipe].output != q.output) { SendMsg(m.conn, MsgChat{ 0, "Can't craft that here (unknown recipe)" }); return; }
                const Recipe& rc = recipes[q.recipe];
                bool unlocked = true;
                if (const int team = TeamOfSession(it->second.session)) {
                    QuestBookState& book = TeamQuests_Book(team);
                    QuestBook_Select(&book); unlocked = RecipeBook_Unlocked(rc); QuestBook_Select(nullptr);
                }
                if (rit->second.creative) unlocked = true; // (#823) creative: every recipe, free
                if (!unlocked) { SendMsg(m.conn, MsgChat{ 0, "Your team hasn't unlocked " + rc.output + " yet" }); return; }
                Inventory& inv = rit->second.inventory;
                const auto before = AllAmounts(inv);
                if (q.recursive && !rit->second.creative) {
                    // (#946) Recursive crafting (#928) planned here, whole: the player's known recipes as the server sees them
                    // (unlocked for its team, every ingredient held once by this inventory), its carried items; then the
                    // steps in order, intermediates first. A step that fails stops the rest (what was made stays).
                    const int team = TeamOfSession(it->second.session);
                    if (team) QuestBook_Select(&TeamQuests_Book(team));
                    std::vector<char> knownK(recipes.size(), 0);
                    for (size_t k = 0; k < recipes.size(); ++k) {
                        bool ok = !team || RecipeBook_Unlocked(recipes[k]);
                        for (const auto& [in, n] : recipes[k].inputs) {
                            bool held = inv.obtainedTotal(in) > 0.0 || inv.count(in) > 0;
                            if (!held) for (const ItemDef& d : itemDB.items()) if (ItemMatches(in, d.id) && (inv.obtainedTotal(d.id) > 0.0 || inv.count(d.id) > 0)) { held = true; break; }
                            ok = ok && held;
                        }
                        knownK[k] = ok;
                    }
                    if (team) QuestBook_Select(nullptr);
                    knownK[q.recipe] = 1; // (the target itself: unlocked, checked above)
                    std::set<std::string> machineMade;
                    for (const auto& ri : Factory::RecipeList()) for (const auto& o : ri.out) machineMade.insert(o.first);
                    std::map<std::string, int> carried;
                    for (const auto& kv : inv.all()) carried[kv.first] = kv.second;
                    const CraftPlan plan = CraftPlan_Make(recipes, [&](int k) { return k >= 0 && k < (int)knownK.size() && knownK[k]; }, carried, (int)q.recipe, q.count, machineMade);
                    if (!plan.ok) { SendMsg(m.conn, MsgChat{ 0, std::string("Missing items for ") + rc.output }); return; }
                    int made = 0;
                    bool stopped = false;
                    for (const CraftStep& s : plan.steps) {
                        for (int k = 0; k < s.times && !stopped; ++k) {
                            if (!inv.craft(recipes[s.recipe])) stopped = true;
                            else if (s.recipe == (int)q.recipe) ++made;
                        }
                        if (stopped) break;
                    }
                    SendChangedItems(m.conn, inv, before);
                    if (made < q.count) SendMsg(m.conn, MsgChat{ 0, TextFormat("Crafted %d of %d %s: a step ran short", made, (int)q.count, rc.output.c_str()) });
                    crafts += made;
                    return;
                }
                int made = 0;
                if (rit->second.creative) { for (; made < q.count; ++made) inv.addAmount(rc.output, rc.outputCount); } // (#823) creative: free
                else while (made < q.count && inv.canCraft(rc) && inv.craft(rc)) ++made; // the output counts for the team (inventory.onAdded)
                if (made > 0) SendChangedItems(m.conn, inv, before);
                if (made < q.count && q.useChests) { // (#854) what the inventory lacks, from the chests near the player; the craft waits for it
                    const int rest = q.count - made;
                    std::vector<std::pair<std::string, int>> missing;
                    for (const auto& [want, n] : rc.inputs) { const int lack = n * rest - inv.count(want); if (lack > 0) missing.push_back({ want, lack }); }
                    std::vector<std::tuple<Planetoid*, Plot*, int, std::string, int>> takes;
                    for (const auto& [b, pc] : ChestsNear(it->second.session, rit->second))
                        for (const auto& [id, have] : Factory::ChestContents(b->id, pc.first->id, pc.second))
                            for (auto& [want, lack] : missing)
                                if (lack > 0 && ItemMatches(want, id)) { const int k = std::min(lack, have); takes.push_back({ b, pc.first, pc.second, id, k }); lack -= k; break; }
                    bool enough = true;
                    for (const auto& [want, lack] : missing) enough = enough && lack <= 0;
                    if (enough) {
                        for (const auto& [b, pl, part, id, k] : takes) Factory::Take(b->id, pl->id, part, id, k, it->second.session); // into its inventory next tick
                        pendingCrafts.push_back({ it->second.session, (int)q.recipe, rest, GetTime() + 3.0 });
                        crafts += made;
                        return;
                    }
                }
                if (made < q.count) SendMsg(m.conn, MsgChat{ 0, made ? TextFormat("Crafted %d of %d: not enough for more", made, (int)q.count) : std::string("Missing items for ") + rc.output });
                crafts += made;
                return;
            }
            case NetMsg::BedUse: { // (#800) a right-click on a bed: within reach, a bed there: that player's spawn
                MsgBedUse b;
                if (!it->second.joined || !NetRead(r, b)) return;
                auto rit = remotes.find(it->second.session);
                if (rit == remotes.end() || rit->second.survival.dead) return;
                Planetoid* body = nullptr;
                for (Planetoid* p : Planetoid::All()) if (p->name == b.body) { body = p; break; }
                Plot* plot = body && body->plots ? body->plots->byId(b.plot) : nullptr;
                const int part = plot ? plot->partAt({ b.cell[0], b.cell[1], b.cell[2] }) : -1;
                bool ok = false;
                if (part >= 0) {
                    const Vector3 bl = plot->toBody({ (b.cell[0] + 1.0f) * kCellSize, b.cell[1] * kCellSize, (b.cell[2] + 1.0f) * kCellSize });
                    const DVec3 at = body->toWorldD(DVec3{ bl.x, bl.y, bl.z });
                    if (ClaimBlocks(it->second.session, body, bl)) { SendMsg(m.conn, MsgChat{ 0, ClaimRefusal(body, bl) + ": not your bed" }); return; } // (#343)
                    ok = (at - rit->second.player.position).length() <= 8.0
                      && Survival_SetBedFor(rit->second.survival, *body, *plot, part, false);
                }
                SendMsg(m.conn, MsgChat{ 0, ok ? "Spawn set: you will wake at this bed" : "That bed is out of reach" });
                return;
            }
            case NetMsg::AnimalHit: { // (#334) a swing at an animal: checked and applied here
                MsgAnimalHit h;
                if (!it->second.joined || !NetRead(r, h)) return;
                OnAnimalHit(world, it->second.session, h.netId, h.action);
                return;
            }
            case NetMsg::Chat: { // (#326) a player's line: the server's handler (commands, logging), else to everyone
                MsgChat c;
                if (!it->second.joined || !NetRead(r, c) || c.text.empty()) return;
                c.from = it->second.session; // (a client speaks only for itself)
                c.text = CleanChat(c.text); // (#349) no control characters into logs and chat
                std::string name; for (const Session& s : sessions) if (s.id == c.from) name = s.name;
                if (TeamCommand(c.from, c.text)) return; // (#604) /team and /t, from anyone
                if (ClaimCommand(c.from, c.text)) return; // (#343) /claim, /unclaim, /claims
                if (chatHandler) chatHandler(c.from, name, c.text);
                else SendToJoined(c);
                return;
            }
            case NetMsg::Loadout: {
                MsgLoadout l;
                if (!it->second.joined || !NetRead(r, l)) return;
                l.session = it->second.session; // (a client speaks only for itself)
                for (Session& s : sessions) if (s.id == l.session) s.loadout = l; // kept as sent: HeldItem checks it at each use (#1868),
                // so a loadout that came a moment before its item (a craft, a take) works once the item is there
                MsgLoadout shown = l;
                if (!shown.held.empty() && HeldItem(l.session).empty()) { shown.held.clear(); ++loadoutsRefused; } // others see bare hands
                SendToJoined(shown);
                return;
            }
            case NetMsg::EditRequest: {
                MsgEditRequest q;
                if (!it->second.joined) return;
                if (!NetRead(r, q)) { printf("host: an edit request from session %u didn't read (%zu bytes)\n", it->second.session, m.data.size()); return; } // (#1107)
                OnEditRequest(world, m.conn, it->second.session, q);
                return;
            }
            case NetMsg::SyncWant: {
                MsgSyncWant w;
                if (!it->second.joined || !NetRead(r, w)) return;
                OnSyncWant(world, m.conn, w);
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
    // (#332) what a player takes from a machine reaches its inventory as the server keeps it, then its client by deltas
    gFactoryDeliver = [](uint32_t to, const std::string& item, double count, float heat) { // (#1102: an amount)
        auto rit = remotes.find(to);
        if (rit == remotes.end() || count <= 0) return;
        const std::map<std::string, double> before{ { item, rit->second.inventory.amount(item) } };
        rit->second.inventory.addAmount(item, count, heat);
        if (const uint32_t conn = ConnOf(to)) SendInventoryDeltas(conn, rit->second.inventory, before);
    };
    // (#343) graves: the owner's team by name (#604's teams), and their team's claims keep a grave closed to others
    gSameTeam = [](const std::string& a, const std::string& b) { const int ta = TeamOfName(a); return ta != 0 && ta == TeamOfName(b); };
    gInOwnerTeamClaim = [](const std::string& owner, const std::string& body, const DVec3& local) {
        const int team = TeamOfName(owner);
        const Claim* c = team ? Claims_At(body, local.toFloat()) : nullptr;
        return c && c->team == team;
    };
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
void Host_SetDeathDrops(DeathDrops d) { deathDrops = d; } // (#804)
int Host_GravesMade() { return gravesMade; }
int Host_DropsMade() { return dropsMade; }
int Host_CreativeGives() { return creativeGives; }
int Host_CreativeRefused() { return creativeRefused; }
int Host_LoadoutsRefused() { return loadoutsRefused; }
std::string Host_HeldItem(uint32_t session) { return HeldItem(session); }
int Host_ConnectionsRefused() { return pendingRefused; }
int Host_SilentDropped() { return silentDropped; }

void Host_Update(Scene& world) {
    hostWorld = &world;
    // Autosave (user, 2026-11-15: every 30 minutes by default; --autosave-minutes N, 0 = off). Wall-clock time: a
    // server's saves protect against crashes, whatever the game clock does.
    if (autosaveSeconds > 0.0 && lastSave >= 0.0 && GetTime() - lastSave >= autosaveSeconds) {
        lastSave = GetTime();
        const int records = SaveAllRecords(); // (#346) everyone online, with the world
        TeamQuests_StoreAll(); // (#338) the team books into their teams' records first
        const bool ok = SaveWorld(world);
        if (ok) ++autosaves;
        printf("host: autosave %d to %s %s (%d player records)\n", autosaves, WorldSavePath(), ok ? "done" : "FAILED", records);
    }
    if (!listening) return;
    Net_Poll();
    { // (#1868) connections admitted (or refused) as they come up, before any message: a silent one counts and times out too
        static std::vector<uint32_t> fresh;
        fresh.clear();
        Net_TakeNewIncoming(fresh);
        for (uint32_t c : fresh) if (!peers.count(c) && Net_IsConnected(c)) AdmitPeer(c);
    }
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
        else if (GetTime() - p.lastHeard > Host_SilenceLimit(p.joined, p.inGame)) silent.push_back(conn);
        else if (!p.joined && GetTime() - p.firstHeard > kHelloSeconds) { silent.push_back(conn); ++silentDropped; } // (#349) never said Hello
    }
    for (uint32_t c : gone) Drop(c, NetBye::Quit);
    for (uint32_t c : silent) Drop(c, NetBye::TimedOut);
    GpuTimer_CpuMark("host: receive, handshakes"); // (#348)
    StepRemotes(world);
    StepSyncs(world);
    GpuTimer_CpuMark("host: remote players");
    SendSnapshots(world);
    SendAnimals(world);   // (#334) and the animals near each player
    SendCreatureFx(world); // (#1343, #1364) their shots and calls
    SendVitals();         // (#337) and each player's survival numbers, to that player
    SendSkills();         // (#920) and their skills, when they changed
    SendTeamStates(); // (#604) each player's team, once a second and when it changes
    StepGraves(world);         // (#804) deaths' drops, pickups, decay and the graves list to everyone
    StepTeamQuests(world); // (#338) team books: evaluate, pay, tell, send
    StepMachines(world);  // (#332) machines: near points, watched views, working states
    StepPlants(world);    // (#335) plant removals to everyone, farm fields near each player
    StepClaims();         // (#343) claim decay, enter/leave lines
    SetInterestCentres(world); // (#334) and animals live around every player (#318's centres)
    {   // (#1224) every connected player's footsteps (no sound on the server): what creatures hear (Footsteps_RecentNoises)
        static std::vector<std::pair<uint32_t, const Player*>> walkers;
        walkers.clear();
        for (const Session& s : sessions) {
            if (!s.online) continue;
            auto rit = remotes.find(s.id);
            if (rit != remotes.end() && rit->second.player.nearest) walkers.push_back({ s.id, &rit->second.player });
        }
        Footsteps_UpdateServer(world, walkers, GetFrameTime());
    }
    SetAnimalTargets(world);   // (#334 slice 2) and mind them: flee, follow, hunt, bite
    GpuTimer_CpuMark("host: snapshots");
}

bool Host_SaveOnShutdown(Scene& world) {
    const int records = SaveAllRecords(); // (#346) everyone online, with the world
    TeamQuests_StoreAll(); // (#338)
    const bool ok = SaveWorld(world);
    printf("host: saved on shutdown to %s %s (%d player records)\n", WorldSavePath(), ok ? "done" : "FAILED", records);
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
// (#1171) An op's /tp: the remote player to a body-local spot, at rest (hovering when flying). Their game follows: its next
// state from us is far from its prediction, so it goes back to ours (NetClient_Step's reconcile).
bool Host_TeleportPlayer(Scene& world, uint32_t session, uint64_t body, const DVec3& local, bool hover) {
    auto it = remotes.find(session);
    if (it == remotes.end()) return false;
    Player::NetState s;
    if (!it->second.player.netState(world, s)) return false;
    s.body = body; s.local = local; s.vel = {}; s.grounded = false; s.hovering = hover;
    if (hover && it->second.creative) it->second.player.flyMode = true;
    return it->second.player.applyNetState(world, s);
}
int Host_EditsApplied() { return editsApplied; }
double Host_InventoryAmount(uint32_t session, const std::string& item) {
    auto it = remotes.find(session);
    return it == remotes.end() ? 0.0 : it->second.inventory.amount(item);
}
int Host_PartsApplied() { return partsApplied; }
void Host_BroadcastPartEdits(const Planetoid& body, const Plot& plot, const std::vector<PartEdit>& edits) { SendPartEdits(body.netId, plot.id, edits); } // (#1684)
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
    for (const Session& s : sessions) out.push_back({ s.id, s.name, s.online, s.online ? Net_RemoteAddress(s.conn) : std::string(), s.uuid });
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
void Host_SetJoinCheck(std::function<bool(const std::string&, const std::string&, const std::string&, std::string&)> check) { joinCheck = std::move(check); }
void Host_SetChatHandler(std::function<void(uint32_t, const std::string&, const std::string&)> handler) { chatHandler = std::move(handler); }
void Host_SetEventHook(std::function<void(const char*, uint32_t, const std::string&)> hook) { eventHook = std::move(hook); }
int Host_SyncsDone() { return syncsDone; }
uint64_t Host_SyncBytesSent() { return syncBytesSent; }
// Tests (#349): messages the guard dropped, and connections it kicked.
int Host_GuardDropped() { return guardDropped; }
int Host_GuardKicks() { return guardKicks; }

// ---- Tests (#346): records across a restart ----
// (#865 tests) remote player `who` stands where `at` stands (its body-local feet and anchor)
bool Host_TestPutPlayerAt(Scene& world, uint32_t who, uint32_t at) {
    auto iw = remotes.find(who), ia = remotes.find(at);
    if (iw == remotes.end() || ia == remotes.end()) return false;
    Player::NetState s;
    return ia->second.player.netState(world, s) && iw->second.player.applyNetState(world, s);
}
// (#865 tests) how far apart two remote players stand, body-local (their world positions go stale between steps); -1: unknown
double Host_PlayersApart(Scene& world, uint32_t a, uint32_t b) {
    auto ia = remotes.find(a), ib = remotes.find(b);
    Player::NetState sa, sb;
    if (ia == remotes.end() || ib == remotes.end() || !ia->second.player.netState(world, sa) || !ib->second.player.netState(world, sb) || sa.body != sb.body) return -1.0;
    return (sa.local - sb.local).length();
}
std::string Host_SessionUuid(uint32_t session) { for (const Session& s : sessions) if (s.id == session) return s.uuid; return std::string(); }
int Host_ItemCount(uint32_t session, const std::string& item) {
    auto it = remotes.find(session);
    return it == remotes.end() ? -1 : it->second.inventory.count(item);
}
// (#920 tests) XP into a remote player's skill (sent with the next SkillState); their food (so eating is not refused as full);
// the perk spends the server allowed and refused.
double Host_SilenceLimit(bool joined, bool inGame) { return joined && !inGame ? kLoadingTimeoutSeconds : kTimeoutSeconds; }
bool Host_TestSkillXp(uint32_t session, const std::string& skill, double xp) {
    auto it = remotes.find(session);
    if (it == remotes.end()) return false;
    Skills_AddXp(it->second.skills, skill, xp, 1.0);
    it->second.skillsDirty = true; it->second.skillsSentAt = -1e9;
    return true;
}
bool Host_TestSetFood(uint32_t session, float food) { auto it = remotes.find(session); if (it == remotes.end()) return false; it->second.survival.food = food; return true; }
bool Host_TestPerk(uint32_t session, const std::string& perk, int rank) {
    auto it = remotes.find(session);
    if (it == remotes.end()) return false;
    it->second.skills.ranks[perk] = rank;
    it->second.skillsDirty = true; it->second.skillsSentAt = -1e9;
    return true;
}
bool Host_TestTeamUp(uint32_t a, uint32_t b) {
    const Session *sa = nullptr, *sb = nullptr;
    for (const Session& s : sessions) { if (s.id == a) sa = &s; if (s.id == b) sb = &s; }
    if (!sa || !sb) return false;
    const int team = Teams_Join(sa->uuid, sa->name, 0xFFFFFFFFu);
    Teams_Move(sb->uuid, team);
    return team != 0 && Teams_Of(sb->uuid) == team;
}
int Host_TeamXpGiven() { return teamXpGiven; }
int Host_SkillSpends() { return skillSpends; }
int Host_SkillSpendsRefused() { return skillSpendsRefused; }
bool Host_TestGive(uint32_t session, const std::string& item, int count) { // (tests) count < 0 takes
    auto it = remotes.find(session);
    if (it == remotes.end()) return false;
    return count >= 0 ? (it->second.inventory.add(item, count), true) : it->second.inventory.remove(item, -count);
}
bool Host_ForgetPlayers() {
    for (const Session& s : sessions) if (s.online) return false;
    sessions.clear(); remotes.clear();
    return true;
}

// ---- Game mode (#823): /gamemode by an op ----
bool Host_SetCreative(const std::string& name, bool on) {
    for (const Session& s : sessions) {
        if (!s.online || !SameName(s.name, name)) continue;
        auto it = remotes.find(s.id);
        if (it == remotes.end()) return false;
        it->second.creative = on;
        it->second.player.flyAllowed = on; // (#781: flying only for creative players, whom ops set)
        it->second.player.doubleTapFly = on; // (#1562) as the client predicts it
        if (!on) it->second.player.flyMode = false;
        lastVitals = 0; // (its client hears at once, with the next vitals)
        return true;
    }
    return false;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
TEST(HostTimeout, ALoadingClientGetsLongerThanOneInTheGame) { // (#1033)
    EXPECT_DOUBLE_EQ(Host_SilenceLimit(true, true), 15.0) << "in the game: pings every second, 15 s of silence is gone";
    EXPECT_GE(Host_SilenceLimit(true, false), 60.0) << "joined, still loading: one long frame (~37 s on the runner) must not drop it";
    EXPECT_DOUBLE_EQ(Host_SilenceLimit(false, false), 15.0) << "not joined yet: the usual limit (Hello has its own)";
}

// (#1868) The held item counts only when the server's copy of the player's inventory has it (creative players excepted), and
// stops counting once it leaves the inventory.
TEST(HostLoadout, TheHeldItemCountsOnlyWhenThePlayerHasIt) {
    const uint32_t id = 0xF1868u;
    Session s;
    s.id = id;
    s.loadout.session = id;
    s.loadout.held = "unit_pickaxe";
    sessions.push_back(s);
    Remote& r = remotes[id];
    EXPECT_EQ(HeldItem(id), "") << "claimed, not in the inventory: bare hands";
    r.inventory.grid()[3] = InvSlot{ "unit_pickaxe", 100, 0.0f };
    EXPECT_EQ(HeldItem(id), "unit_pickaxe") << "in the inventory";
    r.inventory.grid()[3] = InvSlot{};
    EXPECT_EQ(HeldItem(id), "") << "it left the inventory after the loadout came";
    r.creative = true;
    EXPECT_EQ(HeldItem(id), "unit_pickaxe") << "creative: any Catalog item";
    EXPECT_EQ(HeldItem(id + 1), "") << "no such session";
    remotes.erase(id);
    sessions.pop_back();
}

// (#1868) Admission: the per-address cap and the cap on connections that haven't joined yet; loopback exempt.
TEST(HostAdmission, PerAddressAndPendingCaps) {
    if (getenv("RAYLIBGAME_MAX_PER_ADDRESS") || getenv("RAYLIBGAME_MAX_PENDING")) GTEST_SKIP() << "limits overridden by the environment";
    EXPECT_EQ(AdmitVerdict(false, 0, 0), 0);
    EXPECT_EQ(AdmitVerdict(false, 7, 15), 0) << "under both caps";
    EXPECT_EQ(AdmitVerdict(false, 8, 0), 1) << "the 9th from one address";
    EXPECT_EQ(AdmitVerdict(false, 0, 16), 2) << "the 17th not joined yet, from anywhere";
    EXPECT_EQ(AdmitVerdict(true, 100, 100), 0) << "loopback: the host's own machine and the tests";
}
#endif
