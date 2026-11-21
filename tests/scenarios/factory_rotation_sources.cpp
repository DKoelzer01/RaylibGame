// factory/rotation-sources (#1270, #159): a water wheel and a windmill in the sandbox. A trench is dug along the plot's x with a
// 1 m^3/s spring at its head and a pit at its end; a water wheel (axis z, 1.5 m across) stands with its bottom in the trench,
// and once the stream runs it turns with the water (simulated water: Fluid::FlowAt), on its own shaft. A windmill on a mast
// of shafts turns in a set wind (RotationSources_SetWindOverride: the weather's wind changes) at full speed, at part speed in
// a light wind, and stands in a calm. The camera looks across the trench at both.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "fluid.h"
#include "plots.h"
#include "rotationplot.h"
#include "rotationsources.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int wheel = -1, wheelShaft = -1, mill = -1;
double lastNote = -10.0;
Vector3 At(float x, float y, float z = 0.0f) { return plot->toBody({ x, y, z }); }
bool Place(ScenarioCtx& c, PartType t, Int3 at, int orient = 0) { return BuildTool_Place(*c.body, *plot, RotationPlot_PartAt(t, at, orient), false); }
RotationView View(int part) { return Factory::RotationAt(bodyId, plot->id, part); }
}

SCENARIO("factory/rotation-sources") {
    s.tags = { "factory", "machines", "water" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 180.0;
    s.step("a trench with a spring, a water wheel in it, a windmill on a mast", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        terrainTool.creative = true;
        const Vector3 d = Vector3Normalize(c.spawnDirLocal());
        const Vector3 ground = Vector3Scale(d, c.field->radius + c.field->elevation(d) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        // The trench (as fluid/flow-meter digs it): ~2.8 m wide, ~1.6 m deep, along x from -10 to 10 m; a pit at 14 m.
        std::vector<TerrainEdit> cut;
        for (float x = -10.0f; x <= 10.0f; x += 0.8f) { cut.push_back({ At(x, -1.0f), 1.4f, -1.0f, true }); cut.push_back({ At(x, 0.2f), 1.4f, -1.0f, true }); }
        for (int k = 0; k < 4; ++k) cut.push_back({ At(14.0f, -3.0f - 2.5f * k), 4.0f, -1.0f, true });
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        Fluid::AddSource(*c.body, At(-9.5f, -2.2f), 0, 1.0f, 1.5f);
        // The wheel: cells x -1..1 (x -0.5..1 m), y -4..-2 (-2..-0.5 m: its bottom in the stream), z 0; its hub at (0, -3, 0).
        // A shaft on the hub out to z 1 (in the trench's side: parts stand in the ground).
        const bool ok = Place(c, PartType::WaterWheel, { -1, -4, 0 }) && Place(c, PartType::Shaft, { 0, -3, 1 })
            // The windmill: on a mast of 8 vertical shafts at x 8 m, z 6 m; its hub (axis z) on top.
            && Place(c, PartType::Shaft, { 16, 0, 12 }, RotationPlot_OrientFor({ 0, 1, 0 })) && Place(c, PartType::Windmill, { 16, 8, 12 });
        for (int y = 1; y < 8; ++y) Place(c, PartType::Shaft, { 16, y, 12 }, RotationPlot_OrientFor({ 0, 1, 0 }));
        if (!c.check("the wheel, its shaft, the mast and the windmill placed", ok)) return true;
        wheel = plot->partAt({ 0, -3, 0 }); wheelShaft = plot->partAt({ 0, -3, 1 }); mill = plot->partAt({ 16, 8, 12 });
        RotationSources_SetWindOverride(true, 0.0f); // calm until the wind step
        c.note("trench and pit dug, a 1 m^3/s spring; the wheel and the windmill placed");
        return true;
    });
    s.step("the stream reaches the wheel and turns it", 90.0, [](ScenarioCtx& c) {
        const Vector3 up = Vector3Normalize(At(0.0f, 0.0f));
        c.camera(Vector3Add(At(1.5f, 3.0f), Vector3Scale(plot->north, -7.0f)), At(2.0f, 0.5f, 1.5f), up);
        Factory::Sync(c.world);
        const RotationView w = View(wheel);
        if (c.t - lastNote >= 10.0) { lastNote = c.t; c.note(TextFormat("t=%d s: wheel %.0f RPM, flow %.2f m^3/s along x", (int)c.t, w.rpm, Vector3DotProduct(Fluid::FlowAt(*c.body, At(0.25f, -1.6f), 1.0f), plot->east))); }
        c.waitingFor(TextFormat("wheel %.0f RPM", w.rpm));
        if (!w.valid || fabsf(w.rpm) < 1.0f) return false;
        if (c.t < 40.0) return false; // let the stream settle
        c.gt("the wheel turns with the stream (RPM)", fabsf(w.rpm), 1.0);
        const RotationView sh = View(wheelShaft);
        c.near("its shaft turns with it", sh.rpm, w.rpm, 0.01);
        c.gt("it offers stress capacity (SU)", w.capacity, 0.0);
        return true;
    });
    s.step("a full wind turns the windmill at full speed", 10.0, [](ScenarioCtx& c) {
        RotationSources_SetWindOverride(true, 10.0f);
        Factory::Sync(c.world);
        const RotationView m = View(mill);
        if (c.t < 2.0 || !m.valid) return false;
        c.near("full wind: 8 RPM", fabsf(m.rpm), 8.0, 0.01);
        c.eq("four sails: 1024 SU", m.capacity, 1024);
        return true;
    });
    s.step("a light wind turns it slower", 10.0, [](ScenarioCtx& c) {
        RotationSources_SetWindOverride(true, 4.0f);
        Factory::Sync(c.world);
        const RotationView m = View(mill);
        if (c.t < 2.0 || !m.valid) return false;
        c.near("4 m/s: 4 RPM", fabsf(m.rpm), 4.0, 0.01);
        return true;
    });
    s.step("a calm stops it", 10.0, [](ScenarioCtx& c) {
        RotationSources_SetWindOverride(true, 1.0f);
        Factory::Sync(c.world);
        const RotationView m = View(mill);
        if (c.t < 2.0 || !m.valid) return false;
        c.eq("under the cut-in: it stands", m.rpm, 0.0);
        RotationSources_SetWindOverride(false);
        return true;
    });
}
