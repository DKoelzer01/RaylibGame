// net/late-join-lag (#347; the runner's join flakes #1421, #797): two players join late, at once, on a bad link. The game
// hosts (--host; needs a GameNetworkingSockets build, else one SKIP check passes) a world with 3000 dig strokes over 600 m
// around the start; bob joins and walks, and keeps walking while dave and frank join together with nothing cached (each
// fetches every edit), and bob digs once mid-sync. Every send goes through net.cpp's simulated link: 80 ms each way,
// 20 ms jitter, 3 % loss. Checked: both late joiners are in the game and synced (every entry, hashes right) with bob's
// dig applied after their sync; bob's prediction held while they joined; and no send was refused. A refused send prints
// its message type and result (Net_SendFailureSummary, #1836): that line is what a failing run is read for.
#include "scenario.h"
#include "hostserver.h"
#include "netclient.h"
#include "net.h"
#include "player.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <memory>
#include <string>
#include <thread>

namespace {
// A frame at most every 1/60 s of wall time: scenarios step a fixed 1/60 s a frame as fast as the machine goes (thousands
// a second headless), and a client sending an input a frame that fast is a flood the server kicks. A game is real time.
void RealTime() {
    static double next = 0.0;
    const double now = GetTime();
    if (next > now) std::this_thread::sleep_for(std::chrono::duration<double>(next - now));
    next = std::max(now, next) + 1.0 / 60.0;
}
constexpr uint16_t kPort = 27954;
constexpr int kStrokes = 3000;
struct State {
    NetClient bob, dave, frank;
    std::unique_ptr<Player> bobSelf;
    long long stepped = 0;
    double t0 = 0.0, tDone = -1.0;
    int strokes = 0, correctionsAtJoin = 0;
    double worstAtJoin = 0.0;
    bool skip = false, dug = false;
} S;

// bob's steps up to now: walking forward (he turns a little now and then, so he stays on the edited ground).
void StepBob(ScenarioCtx& c) {
    if (!S.bobSelf || S.bob.state != NetClient::State::Joined) return;
    for (; S.stepped < SimTick(); ++S.stepped) {
        PlayerInput in; in.enabled = true; in.hold(Act::Forward);
        if (S.stepped % 240 < 30) in.hold(Act::Right);
        NetClient_Step(S.bob, c.world, *S.bobSelf, in);
    }
}
const char* StateName(const NetClient& n) { return NetClient_StateName(n.state); }
}

SCENARIO("net/late-join-lag") {
    s.tags = { "multiplayer", "net" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-late-join-lag.sav", "--port", std::to_string(kPort) };
    s.env = { { "RAYLIBGAME_NET_LAG_MS", "80" }, { "RAYLIBGAME_NET_JITTER_MS", "20" }, { "RAYLIBGAME_NET_LOSS", "3" }, { "RAYLIBGAME_NET_SEED", "347" } };
    s.timeout = 240.0;
    s.realTime = true; // (#1972) joins over real sockets
    s.step("edits, then bob joins and walks", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!Host_Listening()) { // (the harness has no SKIP verdict: one passing check, said plainly)
            c.check("SKIP: not listening (a build without GameNetworkingSockets)", true);
            S.skip = true; return true;
        }
        if (S.bob.state == NetClient::State::Idle) {
            S.strokes = NetSync_ScatterStrokes(c.world, kStrokes, 347);
            if (!c.check("the server's world has the edits", S.strokes > kStrokes / 2, TextFormat("%d strokes", S.strokes))) return true;
            NetClient_Connect(S.bob, TextFormat("127.0.0.1:%d", kPort), "bob");
            S.t0 = c.clock;
        }
        RealTime(); NetClient_Update(S.bob);
        if (!S.bobSelf && S.bob.state == NetClient::State::Joined && S.bob.haveOwn && !S.bob.syncing)
            if (Planetoid* b = Planetoid::findNet(S.bob.own.body)) {
                S.bobSelf = std::make_unique<Player>();
                S.bobSelf->spawnOn(c.world, *b, b->dirToWorld(Vector3Normalize(S.bob.own.local.toFloat())));
                S.bobSelf->applyNetState(c.world, NetClient_FromMsg(S.bob.own));
                S.bobSelf->flyMode = false;
                S.stepped = SimTick(); S.t0 = c.clock;
            }
        StepBob(c);
        c.waitingFor(TextFormat("bob in the game (%s, syncing %d)", StateName(S.bob), (int)S.bob.syncing));
        if (!S.bobSelf) {
            if (c.clock - S.t0 > 50.0) c.fail(TextFormat("bob not in the game after 50 s (%s, %d of %d entries)", StateName(S.bob), S.bob.syncEntriesDone, S.bob.syncEntriesWanted));
            return false;
        }
        if (c.clock - S.t0 < 2.0) return false; // (walking first)
        S.correctionsAtJoin = S.bob.corrections; S.worstAtJoin = S.bob.worstError; S.bob.worstError = 0.0;
        S.dave.syncAssumeEmpty = S.frank.syncAssumeEmpty = true; // nothing cached: everything comes over
        // His dig mid-sync must be one the server takes: a pickaxe stroke. The server takes a held tool only from the player's
        // inventory (#1868), or anything from a creative player: bob is made creative (his walk is the same).
        c.check("the server makes bob creative (/gamemode)", Host_SetCreative("bob", true));
        NetClient_SetLoadout(S.bob, "stone_pickaxe", false, 0);
        const std::string addr = TextFormat("127.0.0.1:%d", kPort);
        NetClient_Connect(S.dave, addr.c_str(), "dave");
        NetClient_Connect(S.frank, addr.c_str(), "frank");
        return true;
    });
    s.step("dave and frank join together, late", 120.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        RealTime(); for (NetClient* n : { &S.bob, &S.dave, &S.frank }) NetClient_Update(*n);
        StepBob(c);
        if (S.bob.state != NetClient::State::Joined) { // (#1972) why, and the host's side: the runner keeps no host log
            c.fail(TextFormat("bob lost the server (%s: %s); host: %d guard kicks, %d messages dropped by the guard",
                              StateName(S.bob), NetByeText(S.bob.bye), Host_GuardKicks(), Host_GuardDropped()));
            return true;
        }
        // bob digs at his feet once both syncs are under way: they must hold it until their sync is done.
        if (!S.dug && S.dave.syncing && S.frank.syncing && S.dave.syncEntriesDone > 0 && S.frank.syncEntriesDone > 0) {
            uint64_t bid = 0; DVec3 feet{};
            if (Host_PlayerLocal(c.world, S.bob.welcome.session, bid, feet))
                if (Planetoid* b = Planetoid::findNet(bid)) {
                    const Vector3 up = Vector3Normalize(feet.toFloat());
                    const Vector3 hit = Vector3Subtract(feet.toFloat(), Vector3Scale(up, 0.05f));
                    S.dug = NetClient_RequestEdit(S.bob, NetEditKind::Dig, *b, hit, Vector3Subtract(hit, Vector3Scale(up, 0.2f)), 0.4f, "") != 0;
                }
        }
        auto done = [](const NetClient& n) { return n.state == NetClient::State::Joined && !n.syncing && n.syncDoneAt >= 0.0; };
        c.waitingFor(TextFormat("dave %s %d/%d, frank %s %d/%d", StateName(S.dave), S.dave.syncEntriesDone, S.dave.syncEntriesWanted,
                                StateName(S.frank), S.frank.syncEntriesDone, S.frank.syncEntriesWanted));
        const bool closed = S.dave.state == NetClient::State::Closed || S.dave.state == NetClient::State::Rejected ||
                            S.frank.state == NetClient::State::Closed || S.frank.state == NetClient::State::Rejected;
        if (closed || c.clock > 110.0) {
            const std::string why = Net_SendFailureSummary();
            c.fail(TextFormat("late join failed: dave %s (%d/%d), frank %s (%d/%d) after %.0f s; refused sends: %s", StateName(S.dave),
                              S.dave.syncEntriesDone, S.dave.syncEntriesWanted, StateName(S.frank), S.frank.syncEntriesDone,
                              S.frank.syncEntriesWanted, c.clock, why.empty() ? "none" : why.c_str()));
            return true;
        }
        if (!(done(S.dave) && done(S.frank))) return false;
        if (S.tDone < 0.0) S.tDone = c.clock;
        return c.clock - S.tDone > 2.0; // (a moment after: bob's dig, held during the syncs, is applied)
    });
    s.step("checks", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        for (NetClient* n : { &S.bob, &S.dave, &S.frank }) NetClient_Update(*n);
        for (const NetClient* n : { &S.dave, &S.frank }) {
            const char* who = n == &S.dave ? "dave" : "frank";
            c.note(TextFormat("%s: %d of %zu entries, %.2f MB in %.1f s, %d bad hashes, %d live edits held during the sync, %d edits seen; roster %zu",
                              who, n->syncEntriesDone, n->syncIndex.size(), n->syncBytes / 1048576.0, n->syncDoneAt - n->syncStartedAt,
                              n->syncHashBad, n->syncQueued, n->editsSeen, n->roster.size()));
            c.check(TextFormat("%s got every entry", who), n->syncIndex.size() > 10 && n->syncEntriesWanted == (int)n->syncIndex.size() && n->syncEntriesDone == n->syncEntriesWanted,
                    TextFormat("%d of %d wanted, %zu in the index", n->syncEntriesDone, n->syncEntriesWanted, n->syncIndex.size()));
            c.eq(TextFormat("%s: entries with a wrong hash", who), n->syncHashBad, 0);
            c.check(TextFormat("%s sees all three players", who), n->roster.size() >= 3, TextFormat("roster %zu", n->roster.size()));
            if (S.dug) c.check(TextFormat("%s: bob's dig mid-sync held, then applied", who), n->syncQueued >= 1 && n->editsSeen >= 1,
                               TextFormat("held %d, seen %d", n->syncQueued, n->editsSeen));
        }
        c.check("bob dug mid-sync (the server took his stroke)", S.dug);
        const int corrections = S.bob.corrections - S.correctionsAtJoin;
        c.note(TextFormat("bob, while they joined: %d own states, %d corrections, worst prediction error %.3f m (%.3f m before)",
                          S.bob.states, corrections, S.bob.worstError, S.worstAtJoin));
        c.lt("bob's prediction held while they joined (m worst error at an acked input)", S.bob.worstError, 0.5);
        const std::string why = Net_SendFailureSummary();
        c.eq("sends the transport refused", Net_SendFailures(), 0);
        if (!why.empty()) c.note("refused sends: " + why);
        for (NetClient* n : { &S.bob, &S.dave, &S.frank }) NetClient_Disconnect(*n);
        return true;
    });
}
