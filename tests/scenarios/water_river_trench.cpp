// water/river-trench (#1504): what a dig into a river costs the fluid simulation. The river nearest the start (seed 42) gets
// (one at least 1.5 m deep: past its head, the same river whatever the valley widths) a dead-end trench, 15 m long and ~2.4 m wide, cut from its bank straight out, its floor 1.5 m under the river's level:
// normal play (a pond dug off a river). For 40 s: the fluid step time (mean and worst of 1 s windows), active and
// steady bricks, and whether the water sleeps. Run with RAYLIBGAME_RIVER_VALLEY=0 and =1 (#637's valley-filling widths)
// for the same river narrow and wide. Reported only (a measurement); tagged local (timing).
#include "scenario.h"
#include "fluid.h"
#include "rivers.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
struct TState {
    bool dug = false;
    double t0 = 0.0, lastT = 0.0, winStart = 0.0;
    float half = 0.0f;
    double sumStep = 0.0; int nStep = 0; double worstWin = 0.0, winSum = 0.0; int winN = 0;
    int maxActive = 0, maxSteady = 0, asleepAt = -1;
};
TState T;
}

SCENARIO("water/river-trench") {
    s.tags = { "water", "local" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;
    s.record = false;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 1.0 && WorldIdle(); });
    s.step("dig a trench off the nearest river", 5.0, [](ScenarioCtx& c) {
        T = TState{};
        const RiverNetwork* net = c.field->rivers();
        if (!net) { c.fail("no rivers"); return true; }
        // The nearest river point to the start: a spiral of 10 m rings out to 3 km.
        const Vector3 home = c.spawnDirLocal();
        const Vector3 ax = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), ay = Vector3CrossProduct(ax, home);
        const float R = c.field->radius;
        Vector3 hit{}; bool found = false;
        for (float r = 0.0f; r < 3000.0f && !found; r += 10.0f)
            for (int k = 0; k < 64 && !found; ++k) {
                const float a = 2.0f * PI * k / 64.0f;
                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(ax, r * cosf(a) / R), Vector3Scale(ay, r * sinf(a) / R))));
                if (c.field->waterSourceAt(d).kind == 1 && net->nearest(d, 60.0f).depth >= 1.5f) { hit = d; found = true; } // past its head: depth (from the discharge alone) picks the same river with valley widths off or on
            }
        if (!found) { c.fail("no river within 3 km"); return true; }
        const RiverNetwork::Hit rh = net->nearest(hit, 60.0f);
        const PlanetField::WaterSource ws = c.field->waterSourceAt(hit);
        T.half = rh.halfWidth;
        // Across the river, away from its centre: from the bank out 15 m.
        const Vector3 side = Vector3Normalize(Vector3CrossProduct(rh.down, hit));
        std::vector<TerrainEdit> cut;
        const float floorR = R + ws.level - 1.5f + 1.2f; // sphere centres 1.2 m over the floor
        for (float o = rh.halfWidth - 1.0f; o <= rh.halfWidth + 15.0f; o += 0.8f) {
            const Vector3 d = Vector3Normalize(Vector3Add(hit, Vector3Scale(side, o / R)));
            cut.push_back({ Vector3Scale(d, floorR), 1.2f, -1.0f, true });
            cut.push_back({ Vector3Scale(d, floorR + 1.5f), 1.2f, -1.0f, true });
        }
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.note(TextFormat("river at %.0f m: half width %.1f m, level %.2f m; trench from its bank 15 m out, floor 1.5 m under the level (%s)",
                          Vector3Distance(Vector3Scale(home, R), Vector3Scale(hit, R)), rh.halfWidth, ws.level,
                          getenv("RAYLIBGAME_RIVER_VALLEY") && getenv("RAYLIBGAME_RIVER_VALLEY")[0] == '0' ? "valley widths off" : "valley widths on"));
        const Vector3 eye = Vector3Add(Vector3Scale(hit, R + ws.level + 25.0f), Vector3Scale(side, 25.0f));
        c.camera(eye, Vector3Scale(hit, R + ws.level), hit);
        T.dug = true; T.t0 = T.winStart = c.T;
        return true;
    });
    s.step("40 s of water", 60.0, [](ScenarioCtx& c) {
        const Fluid::Stats fs = Fluid::GetStats();
        T.sumStep += fs.stepMs; ++T.nStep;
        T.winSum += fs.stepMs; ++T.winN;
        if (c.T - T.winStart >= 1.0) { T.worstWin = std::max(T.worstWin, T.winSum / std::max(1, T.winN)); T.winSum = 0.0; T.winN = 0; T.winStart = c.T; }
        T.maxActive = std::max(T.maxActive, fs.activeBricks);
        T.maxSteady = std::max(T.maxSteady, fs.steadyBricks);
        if (T.asleepAt < 0 && c.T - T.t0 > 2.0 && Fluid::Asleep()) T.asleepAt = (int)(c.T - T.t0);
        c.waitingFor("40 s of water after the dig");
        if (c.T - T.t0 < 40.0) return false;
        c.note(TextFormat("fluid step: mean %.2f ms, worst 1 s window %.2f ms; bricks: up to %d active, %d steady; asleep %s",
                          T.sumStep / std::max(1, T.nStep), T.worstWin, T.maxActive, T.maxSteady,
                          T.asleepAt >= 0 ? TextFormat("after %d s", T.asleepAt) : "never (40 s)"));
        c.metric("fluid step mean", T.sumStep / std::max(1, T.nStep), "ms");
        c.metric("fluid step worst 1 s", T.worstWin, "ms");
        c.metric("river half width", T.half, "m", false);
        c.check("dug the trench", T.dug);
        return true;
    });
}
