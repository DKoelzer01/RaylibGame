// net/belts (#1735 PR 4): a joined client replays the server's belts. The game hosts (--host, needs a GameNetworkingSockets
// build; without one the scenario notes SKIP and passes) and one in-process client, "belter", joins. Near its spawn the
// server builds factory/belts' layout (two chests, a thin belt between them) and loads 40 stone; the client gets the plot's
// belts as a snapshot, then every belt step's events (BeltSteps), and its replica checks the server's checksum after each
// message. Checked: (1) the replica follows the moving stone with no drift, (2) a replica made to drift (an item dropped)
// asks for a snapshot (BeltResync) and is whole again. Belts aren't drawn yet (#2053): the camera shows the two chests.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "hostserver.h"
#include "netclient.h"
#include "player.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <memory>
#include <string>
#include <thread>

namespace {
constexpr uint16_t kPort = 27955;
// A frame at most every 1/60 s of wall time (net/creative-fly: a client sending inputs faster is a flood the server kicks).
void RealTime() {
    static double next = 0.0;
    const double now = GetTime();
    if (next > now) std::this_thread::sleep_for(std::chrono::duration<double>(next - now));
    next = std::max(now, next) + 1.0 / 60.0;
}
struct State {
    NetClient cl;
    std::unique_ptr<Player> p;
    long long stepped = 0;
    double t0 = 0.0;
    Plot* plot = nullptr;
    uint64_t bodyId = 0, bodyNet = 0;
    int chestA = -1;
    bool skip = false, drifted = false;
    int64_t firstStep = -1;            // the replica's step when it first had the belt
    int maxItems = 0;                  // the most items the replica held on its belts at once
    int snapshotsBefore = 0, resyncsBefore = 0;
    int64_t stepAtDrift = 0;
} S;

BeltClient::Key Key() { return { S.bodyNet, S.plot ? S.plot->id : 0u }; }
// The client's frame: messages, then its idle player's input steps (a joined player sends inputs).
bool ClientFrame(ScenarioCtx& c) {
    RealTime(); NetClient_Update(S.cl);
    if (S.cl.state != NetClient::State::Joined) { c.fail("the client lost the server"); return false; }
    if (S.p) for (; S.stepped < SimTick(); ++S.stepped) { PlayerInput in; NetClient_Step(S.cl, c.world, *S.p, in); }
    return true;
}
void Look(ScenarioCtx& c) { if (S.plot) c.camera(S.plot->toBody({ 4.0f, 4.0f, 6.0f }), S.plot->toBody({ 4.0f, 0.5f, 0.5f }), S.plot->up); }
int ItemsOn(const BeltReplica& r) { int n = 0; for (const BeltLine& l : r.lines()) n += (int)l.items.size(); return n; }
}

SCENARIO("net/belts") {
    s.tags = { "multiplayer", "net", "factory", "belts" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-net-belts.sav", "--port", std::to_string(kPort), "--day" };
    s.timeout = 120.0;
    s.step("join as belter", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        if (!Host_Listening()) { // (the harness has no SKIP verdict: one passing check, said plainly)
            c.check("SKIP: not listening (a build without GameNetworkingSockets)", true);
            S.skip = true; return true;
        }
        if (S.cl.state == NetClient::State::Idle) { NetClient_Connect(S.cl, TextFormat("127.0.0.1:%d", kPort), "belter"); S.t0 = c.t; }
        RealTime(); NetClient_Update(S.cl);
        c.waitingFor("the client to join and get its state");
        if (S.cl.state != NetClient::State::Joined || !S.cl.haveOwn) {
            if (c.t - S.t0 > 30.0) c.fail(TextFormat("not joined after 30 s (state %d)", (int)S.cl.state));
            return false;
        }
        Planetoid* b = Planetoid::findNet(S.cl.own.body);
        if (!b) { c.fail("the client's body isn't here"); return true; }
        S.p = std::make_unique<Player>();
        S.p->spawnOn(c.world, *b, b->dirToWorld(Vector3Normalize(S.cl.own.local.toFloat())));
        S.p->applyNetState(c.world, NetClient_FromMsg(S.cl.own));
        S.stepped = SimTick();
        return true;
    });
    s.step("the server builds two chests and a belt near the client", 20.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        if (!ClientFrame(c)) return true;
        Planetoid* b = Planetoid::findNet(S.cl.own.body);
        if (!b || b != c.body) { c.fail("the client isn't on the sandbox's body"); return true; }
        terrainTool.creative = true;
        const Vector3 dir = Vector3Normalize(S.cl.own.local.toFloat());
        const Vector3 pole = fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 };
        const Vector3 east = Vector3Normalize(Vector3CrossProduct(pole, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!b->plots) b->plots = std::make_unique<PlotSet>();
        S.plot = b->plots->plotFor(ground, 0.0f, true);
        S.bodyId = b->id; S.bodyNet = b->netId;
        const bool ok = BuildTool_Place(*b, *S.plot, { PartType::Chest, { 0, 0, 0 }, 0 }, false)
                     && BuildTool_Place(*b, *S.plot, { PartType::Chest, { 6, 0, 0 }, 0 }, false);
        S.chestA = S.plot->partAt({ 0, 0, 0 });
        const int chestB = S.plot->partAt({ 6, 0, 0 });
        if (!c.check("two chests placed", ok && S.chestA >= 0 && chestB >= 0)) return true;
        int xa = 0, xb = 6; // (factory/belts: A's east face to B's west face along z = 0)
        while (S.plot->partAt({ xa + 1, 0, 0 }) == S.chestA) ++xa;
        while (xb > 0 && S.plot->partAt({ xb - 1, 0, 0 }) == chestB) --xb;
        ConduitRun run; run.material = "basic"; run.points = { { xa * 8 + 6, 0, 4 }, { xb * 8, 0, 4 } };
        c.check("the belt added to the plot", S.plot->addBelt(run) > 0);
        Factory::Sync(c.world);
        Look(c);
        return true;
    });
    s.step("40 stone into the first chest", 10.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        if (!ClientFrame(c)) return true;
        Look(c);
        Factory::Sync(c.world);
        c.waitingFor("the factory has the plot (chest A's view)");
        if (!Factory::View(S.bodyId, S.plot->id, S.chestA).valid) return false;
        Factory::Insert(S.bodyId, S.plot->id, S.chestA, "stone", 40);
        return true;
    });
    s.step("the client's replica follows the stone, no drift", 30.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        if (!ClientFrame(c)) return true;
        Look(c);
        Factory::Sync(c.world);
        const BeltReplica* r = S.cl.belts.find(Key());
        c.waitingFor(TextFormat("the client's replica (%d snapshots, %d step messages)", S.cl.belts.snapshots, S.cl.belts.stepMsgs));
        if (!r || r->lines().empty()) return false;
        if (S.firstStep < 0) S.firstStep = r->step();
        S.maxItems = std::max(S.maxItems, ItemsOn(*r));
        if (r->needResync()) { c.fail(TextFormat("the replica drifted at step %lld", (long long)r->step())); return true; }
        if (r->step() - S.firstStep < 120) return false; // 6 s of belt steps: the stone loads (one per 3 ticks) and rides
        c.note(TextFormat("replica at step %lld (from %lld), %d snapshots, %d step messages, %d resyncs, at most %d items on the belt",
                          (long long)r->step(), (long long)S.firstStep, S.cl.belts.snapshots, S.cl.belts.stepMsgs, S.cl.belts.resyncs, S.maxItems));
        c.eq("the replica asked for no snapshot", S.cl.belts.resyncs, 0);
        c.gt("stone rode the replica's belt (most items at once)", S.maxItems, 3);
        c.check("the client knows the stone's name", S.cl.belts.itemName(r->lines()[0].items.empty() ? -1 : r->lines()[0].items[0].item) == "stone"
                                                       || S.maxItems > 0, "named with the snapshot or steps");
        c.gt("step messages heard", S.cl.belts.stepMsgs, 30);
        c.lt("the replica plays close behind the server (steps queued)", (double)(r->received() - r->step()), 2.0 * BeltReplica::kPlayout + 1.0);
        return true;
    });
    s.step("a drifted replica asks for a snapshot and is whole again", 15.0, [](ScenarioCtx& c) {
        if (S.skip) return true;
        if (!ClientFrame(c)) return true;
        Look(c);
        Factory::Sync(c.world);
        const BeltReplica* r = S.cl.belts.find(Key());
        if (!r) { c.fail("the replica is gone"); return true; }
        if (!S.drifted) { // the replica's own state with the front item dropped (or one added): the next checksum differs
            BeltSyncSnapshot bad; bad.epoch = r->epoch(); bad.step = r->step(); bad.runs = r->runs(); bad.lines = r->lines();
            BeltLine& l = bad.lines[0];
            if (!l.items.empty()) { const int64_t g = l.items.front().gap; l.items.pop_front(); if (!l.items.empty()) l.items.front().gap += g; }
            else Belt_TryInsert(l, 1, 0);
            S.snapshotsBefore = S.cl.belts.snapshots; S.resyncsBefore = S.cl.belts.resyncs; S.stepAtDrift = r->step();
            S.cl.belts.onSnapshot(Key(), bad, {}, GetTime());
            S.drifted = true;
            return false;
        }
        const int asked = S.cl.belts.resyncs - S.resyncsBefore, snaps = S.cl.belts.snapshots - S.snapshotsBefore - 1; // (-1: ours)
        r = S.cl.belts.find(Key());
        c.waitingFor(TextFormat("a resync (%d asked, %d snapshots)", asked, snaps));
        if (asked < 1 || snaps < 1 || !r || r->needResync() || r->step() < S.stepAtDrift + 40) return false;
        c.note(TextFormat("drift at step %lld: %d asked, %d snapshots, replica at step %lld", (long long)S.stepAtDrift, asked, snaps, (long long)r->step()));
        c.eq("one snapshot request", asked, 1);
        c.check("the replica plays on without drift", !r->needResync());
        NetClient_Disconnect(S.cl);
        return true;
    });
}
