// visual/belts (#2053): a thin belt drawn with the items riding it. An S-shaped run (two corners, one each way) from chest
// A's east face to chest B's west face, laid in code (Plot::addBelt), and stone, iron ore and coal loaded into A: the
// factory publishes each item's position every tick (Factory::BeltItems), the renderer draws the slab, rails, scrolling
// chevrons, the items and their icons (Belts_LastStats). The camera looks down on the run from its open side.
#include "scenario.h"
#include "beltmodels.h"
#include "buildtool.h"
#include "factory.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int chestA = -1, chestB = -1;

void Look(ScenarioCtx& c) { c.camera(plot->toBody({ 3.4f, 1.1f, 2.4f }), plot->toBody({ 2.5f, 0.0f, 1.2f }), plot->up); }
int Items(Factory::BeltPlotView& v) { int n = 0; for (const auto& l : v.lines) n += (int)l.pos.size(); return n; }
int64_t Sum(Factory::BeltPlotView& v) { int64_t s = 0; for (const auto& l : v.lines) for (int64_t p : l.pos) s += p; return s; }
}

SCENARIO("visual/belts") {
    s.tags = { "factory", "belts", "visual" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("two chests and an S-shaped belt between them", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("2 s of play and an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        const bool ok = BuildTool_Place(*c.body, *plot, { PartType::Chest, { 0, 0, 0 }, 0 }, false)
                     && BuildTool_Place(*c.body, *plot, { PartType::Chest, { 8, 0, 2 }, 0 }, false);
        chestA = plot->partAt({ 0, 0, 0 }); chestB = plot->partAt({ 8, 0, 2 });
        if (!c.check("two chests placed", ok && chestA >= 0 && chestB >= 0)) return true;
        int xa = 0, xb = 8;
        while (plot->partAt({ xa + 1, 0, 0 }) == chestA) ++xa;
        while (xb > 0 && plot->partAt({ xb - 1, 0, 2 }) == chestB) --xb;
        // A's east face (cell row z 0: slot z 4) east to x 5 cells, a left turn north to row z 2 (slot z 20), then east into
        // B's west face: one corner turning each way.
        ConduitRun run; run.material = "basic"; run.points = { { xa * 8 + 6, 0, 4 }, { 40, 0, 4 }, { 40, 0, 20 }, { xb * 8, 0, 20 } };
        c.note(TextFormat("belt from x %d to %d", xa * 8 + 6, xb * 8));
        c.check("the belt added to the plot", plot->addBelt(run) > 0);
        Factory::Sync(c.world);
        Look(c);
        return true;
    });
    s.step("stone, iron ore and coal into the first chest", 10.0, [](ScenarioCtx& c) {
        Look(c);
        Factory::Sync(c.world);
        c.waitingFor("the factory has the plot (chest A's view)");
        if (!Factory::View(bodyId, plot->id, chestA).valid) return false;
        Factory::Insert(bodyId, plot->id, chestA, "stone", 24);
        Factory::Insert(bodyId, plot->id, chestA, "iron_ore", 24);
        Factory::Insert(bodyId, plot->id, chestA, "coal", 300);
        return true;
    });
    s.step("items ride the belt, drawn", 15.0, [](ScenarioCtx& c) {
        Look(c);
        Factory::Sync(c.world);
        Factory::BeltPlotView v;
        const bool got = Factory::BeltItems(bodyId, plot->id, v);
        c.waitingFor(TextFormat("%d items on the belt", got ? Items(v) : -1));
        if (c.t < 4.0 || !got || Items(v) < 10) return c.t > 12.0 ? (c.check("items on the belt", false), true) : false;
        const BeltDrawStats& st = Belts_LastStats();
        c.eq("one run drawn", st.runs, 1);
        c.eq("three segments", st.segments, 3);
        c.gt("chevrons", st.chevrons, 10);
        c.gt("items drawn", st.items, 9);
        c.gt("icons on them", st.icons, 9);
        c.eq("the factory names the items", (int)v.names.empty(), 0);
        const double ahead = Factory::TickClock() - (double)v.tick;
        c.check(TextFormat("the tick clock is 0 .. 2 ticks past the view (%.2f)", ahead), ahead >= 0.0 && ahead <= 2.0);
        return true;
    });
    // Long enough for a capture of the moving belt (record from ~7 s after idle).
    s.step("and keep moving", 20.0, [](ScenarioCtx& c) {
        static long long tick0 = -1;
        static int64_t sum0 = -1;
        Look(c);
        Factory::Sync(c.world);
        Factory::BeltPlotView v;
        if (!Factory::BeltItems(bodyId, plot->id, v) || v.lines.empty() || v.lines[0].pos.empty()) return c.t > 18.0 ? (c.check("items on the belt", false), true) : false;
        if (tick0 < 0) { tick0 = v.tick; sum0 = Sum(v); }
        if (c.t < 15.0) return false;
        c.note(TextFormat("ticks %lld -> %lld, %d items, positions summed %lld then %lld (1/64 cell)", tick0, v.tick, Items(v), (long long)sum0, (long long)Sum(v)));
        c.gt("the factory ticked on", (double)(v.tick - tick0), 100.0);
        c.check("the items moved", Sum(v) != sum0);
        return true;
    });
}
