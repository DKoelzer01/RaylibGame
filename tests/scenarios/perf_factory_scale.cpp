// perf/factory-scale (#1207, #997; the factory doc's section 4): N machines, all near (default 5000; RAYLIBGAME_FACTORY_N),
// in working chains on plots of ~400 machines: each group a generator on a cable to two macerators (stone -> sand) and a
// chest on their conveyor. After a 10 s warm-up the factory's own last-minute ring gives the tick's mean, p99 and worst
// and its phases (commands, plots, views). Reported only (perf policy); the doc's target, < 5 ms mean and p99 < 8 ms at
// 5000 all near, is for scale 2 (parallel plots). Quiet machine only.
#include "scenario.h"
#include "factory.h"
#include "plots.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

namespace {
constexpr int kMachinesPerPlot = 400;     // (guess, the doc's)
constexpr int kGroupMachines = 4;         // generator, 2 macerators, chest
constexpr float kPlotSpacing = 80.0f;     // m between plot origins (all within 1000 m of the camera: near)
constexpr double kWarmup = 10.0, kMeasure = 60.0; // s of factory time
int target = 5000;
double t0 = 0.0;
long long made0 = 0;
float idleShare = 0.0f; // RAYLIBGAME_FACTORY_IDLE: the share of groups starved (idle machines: the sleeping-machines case)
uint64_t body0 = 0; uint32_t plot0 = 0; int mac0 = -1, chest0 = -1; // the first group (the views checked at the end)
}

SCENARIO("perf/factory-scale") {
    s.tags = { "perf", "factory" };
    s.world = ScenarioWorld::Sandbox;
    s.record = false;
    s.timeout = 150.0;
    s.step("N machines in working chains on plots of 400", 20.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        if (const char* n = getenv("RAYLIBGAME_FACTORY_N")) target = std::max(kGroupMachines, atoi(n));
        if (const char* f = getenv("RAYLIBGAME_FACTORY_IDLE")) idleShare = std::clamp((float)atof(f), 0.0f, 1.0f); // (starved groups: coal, no stone)
        const int groups = target / kGroupMachines, perPlot = kMachinesPerPlot / kGroupMachines;
        const int plots = (groups + perPlot - 1) / perPlot, side = (int)ceilf(sqrtf((float)plots));
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir)), north = Vector3CrossProduct(dir, east);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        int placed = 0;
        for (int pi = 0; pi < plots; ++pi) {
            const float ex = ((pi % side) - (side - 1) * 0.5f) * kPlotSpacing, nz = ((pi / side) - (side - 1) * 0.5f) * kPlotSpacing;
            const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Add(Vector3Scale(east, ex), Vector3Scale(north, nz))));
            Plot* plot = c.body->plots->plotFor(Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f), 0.0f, true);
            const int here = std::min(perPlot, groups - pi * perPlot), perRow = 10;
            for (int g = 0; g < here; ++g) {
                const int bx = (g % perRow) * 12, bz = (g / perRow) * 5;
                const int gen = plot->place({ PartType::Generator, { bx, 0, bz }, 0 });
                const int m1 = plot->place({ PartType::Macerator, { bx + 3, 0, bz }, 0 });
                const int m2 = plot->place({ PartType::Macerator, { bx + 6, 0, bz }, 0 });
                const int ch = plot->place({ PartType::Chest, { bx + 9, 0, bz }, 0 });
                for (int x = 0; x < 8; ++x) plot->place({ PartType::Cable, { bx + x, 0, bz - 1 }, 0 });
                for (int x = 3; x < 11; ++x) plot->place({ PartType::Conveyor, { bx + x, 0, bz + 2 }, 0 });
                if (gen < 0 || m1 < 0 || m2 < 0 || ch < 0) continue;
                placed += 4;
                if (chest0 < 0) { body0 = c.body->id; plot0 = plot->id; mac0 = m1; chest0 = ch; }
                Factory::Sync(c.world); // (the layout first, so the seeds find the machines)
                Factory::Seed(c.body->id, plot->id, gen, "coal", 64);
                const int gi = pi * perPlot + g; // (#1207) a share of the groups starved: their macerators idle, their generator unloaded
                if (std::floor((gi + 1) * idleShare) > std::floor(gi * idleShare)) continue;
                Factory::Seed(c.body->id, plot->id, m1, "stone", 64);
                Factory::Seed(c.body->id, plot->id, m2, "stone", 64);
            }
        }
        c.metric("machines", placed, "");
        c.metric("plots", plots, "");
        c.check("the machines placed", placed >= target * 9 / 10, TextFormat("%d of %d", placed, target));
        c.camera(Vector3Scale(dir, c.field->radius + c.field->elevation(dir) + 120.0f), Vector3Scale(dir, c.field->radius), north); // the field from above
        t0 = c.T;
        return true;
    });
    s.step("warm up", kWarmup + 5.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        if (c.T - t0 < kWarmup) return false;
        made0 = Factory::Stats().itemsMade; t0 = c.T;
        return true;
    });
    s.step("a minute of factory ticks", kMeasure + 10.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        c.waitingFor(TextFormat("%.0f / %.0f s", c.T - t0, kMeasure));
        if (c.T - t0 < kMeasure + 1.0) return false;
        const FactoryStats st = Factory::Stats();
        c.metric("factory tick mean", st.meanTickMs, "ms");
        c.metric("factory tick p99", st.p99TickMs, "ms");
        c.metric("factory tick worst", st.worstRecentMs, "ms");
        c.metric("commands", st.cmdMs, "ms");
        c.metric("plots", st.tickOnlyMs, "ms");
        c.metric("views", st.viewMs, "ms");
        c.metric("machines working", st.activeMachines, "", false);
        c.metric("items made in the minute", (double)(st.itemsMade - made0), "", false);
        c.gt("ticks measured", st.recentTicks, 1000.0);
        c.gt("the fed chains work", st.activeMachines, st.machines * 0.4 * (1.0 - idleShare));
        c.metric("idle share (starved groups)", idleShare, "", false);
        // (#1207) the near-plot views (republished only on change) still show the chest filling and the macerator working
        int sand = 0;
        for (const auto& [item, n] : Factory::ChestContents(body0, plot0, chest0)) if (item == "sand") sand += n;
        c.gt("the first chest shows the sand its conveyor brought", sand, 0.0);
        c.check("the first macerator shows working", Factory::Working(body0, plot0, mac0));
        c.note(TextFormat("doc target (scale 2): mean < 5 ms, p99 < 8 ms at 5000 all near; now %.2f / %.2f ms", st.meanTickMs, st.p99TickMs));
        return true;
    });
}
