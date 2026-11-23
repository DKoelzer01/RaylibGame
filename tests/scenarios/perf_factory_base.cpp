// perf/factory-base (#1690): a budget scene. Sandbox; one plot of 300 working machines (75 groups as in perf/factory-scale: a
// generator on a cable to two macerators and a chest on their conveyor, all fed), rendered from a camera that circles the
// base 70 m up. 5 s of warm-up not measured, then 20 s timed: frame mean and p95, and the factory tick's mean and p99 (the
// scene's own cost). Limits live in tools/ci/perf_budgets.txt (runner only).
#include "perf_budget.h"
#include "factory.h"
#include "plots.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
constexpr int kGroups = 75, kPerRow = 10;
struct St { BudgetFrames frames; Vector3 dir{}, east{}, north{}; float ground = 0.0f; double t0 = 0.0; int placed = 0; } S;
void Orbit(ScenarioCtx& c, double t) {
    const float a = (float)(t * 0.15);
    const Vector3 centre = Vector3Add(Vector3Scale(S.dir, S.ground), Vector3Add(Vector3Scale(S.east, 50.0f), Vector3Scale(S.north, 15.0f)));
    const Vector3 off = Vector3Add(Vector3Scale(S.east, 60.0f * cosf(a)), Vector3Scale(S.north, 60.0f * sinf(a)));
    c.camera(Vector3Add(Vector3Add(centre, off), Vector3Scale(S.dir, 70.0f)), centre, S.dir);
}
}

SCENARIO("perf/factory-base") {
    s.tags = { "perf", "budget" }; // (not "factory": the runner runs that tag as part of its scenarios)
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    s.record = false;
    s.step("300 machines in working chains", 20.0, [](ScenarioCtx& c) {
        if (c.T < 1.0 || !WorldIdle()) return false;
        S.dir = c.spawnDirLocal();
        S.east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, S.dir));
        S.north = Vector3CrossProduct(S.dir, S.east);
        S.ground = c.field->radius + c.field->elevation(S.dir);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(Vector3Scale(S.dir, S.ground + 0.05f), 0.0f, true);
        for (int g = 0; g < kGroups; ++g) {
            const int bx = (g % kPerRow) * 12, bz = (g / kPerRow) * 5;
            const int gen = plot->place({ PartType::Generator, { bx, 0, bz }, 0 });
            const int m1 = plot->place({ PartType::Macerator, { bx + 3, 0, bz }, 0 });
            const int m2 = plot->place({ PartType::Macerator, { bx + 6, 0, bz }, 0 });
            const int ch = plot->place({ PartType::Chest, { bx + 9, 0, bz }, 0 });
            for (int x = 0; x < 8; ++x) plot->place({ PartType::Cable, { bx + x, 0, bz - 1 }, 0 });
            for (int x = 3; x < 11; ++x) plot->place({ PartType::Conveyor, { bx + x, 0, bz + 2 }, 0 });
            if (gen < 0 || m1 < 0 || m2 < 0 || ch < 0) continue;
            S.placed += 4;
            Factory::Sync(c.world); // (the layout first, so the seeds find the machines)
            Factory::Seed(c.body->id, plot->id, gen, "coal", 64);
            Factory::Seed(c.body->id, plot->id, m1, "stone", 64);
            Factory::Seed(c.body->id, plot->id, m2, "stone", 64);
        }
        c.metric("machines", S.placed, "", false);
        c.check("the machines placed", S.placed >= kGroups * 4 * 9 / 10, TextFormat("%d of %d", S.placed, kGroups * 4));
        Orbit(c, 0.0);
        S.t0 = c.T;
        return true;
    });
    s.step("warm-up", kBudgetWarmup + 5.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        Orbit(c, c.T - S.t0);
        return c.T - S.t0 > kBudgetWarmup;
    });
    s.step("measure", kBudgetMeasure + 10.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        Orbit(c, c.T - S.t0);
        S.frames.frame();
        c.waitingFor("the timed window");
        if (c.t <= kBudgetMeasure) return false;
        S.frames.report(c);
        const FactoryStats st = Factory::Stats();
        c.metric("factory tick mean", st.meanTickMs, "ms");
        c.metric("factory tick p99", st.p99TickMs, "ms");
        c.metric("machines working", st.activeMachines, "", false);
        c.gt("the fed chains work", st.activeMachines, st.machines * 0.4);
        return true;
    });
}
