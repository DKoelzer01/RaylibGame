// determinism/state-hash (#1015): per-subsystem state hashes against goldens (tests/scenarios/golden_hashes.txt), in the
// sandbox. In one frame (so no frame timing gets in): a trench dug with a pit at its end (terrain edits), a hand crank,
// a gearbox, a shaft and a millstone with stone; then 200 factory ticks with the crank turned every 5 (WaitTicks), and the
// edits and the factory hashed. Then a spring at the trench's head runs for 15 s in fluid lockstep and the water is hashed.
// The same seed and steps must give the same hashes on every run and machine (the Linux runner and the Windows PC): a
// difference means the simulation doesn't repeat (thread timing, a float or libm result, undefined behaviour), unless the
// rules changed on purpose (then RAYLIBGAME_GOLDEN_UPDATE=1 and say why in the PR).
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "fluid.h"
#include "plots.h"
#include "rotationplot.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>
#include <vector>

namespace {
Plot* plot = nullptr;
Vector3 At(float x, float y, float z = 0.0f) { return plot->toBody({ x, y, z }); }
bool Place(ScenarioCtx& c, PartType t, Int3 at, int orient = 0) { return BuildTool_Place(*c.body, *plot, RotationPlot_PartAt(t, at, orient), false); }
}

SCENARIO("determinism/state-hash") {
    s.tags = { "factory", "water", "determinism" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.env = { { "RAYLIBGAME_FACTORY_LOCKSTEP", "1" }, { "RAYLIBGAME_FLUID_LOCKSTEP", "1" } };
    s.timeout = 120.0;
    s.step("dig, build, crank 200 factory ticks; hash the edits and the factory", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        terrainTool.creative = true;
        const Vector3 d = Vector3Normalize(c.spawnDirLocal());
        const Vector3 ground = Vector3Scale(d, c.field->radius + c.field->elevation(d) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        // The trench along x (-8..4 m, as factory/rotation-sources digs it) and a pit at 8 m; the machines stand beside it at z 4 m.
        std::vector<TerrainEdit> cut;
        for (float x = -8.0f; x <= 4.0f; x += 0.8f) { cut.push_back({ At(x, -1.0f), 1.4f, -1.0f, true }); cut.push_back({ At(x, 0.2f), 1.4f, -1.0f, true }); }
        for (int k = 0; k < 3; ++k) cut.push_back({ At(8.0f, -3.0f - 2.5f * k), 4.0f, -1.0f, true });
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        // Millstone (cells 0..1 high), a vertical shaft over it, a gearbox, and the crank on the gearbox's side (axis x).
        const int up = RotationPlot_OrientFor({ 0, 1, 0 }), alongX = RotationPlot_OrientFor({ 1, 0, 0 });
        const bool ok = Place(c, PartType::Millstone, { 0, 0, 8 }) && Place(c, PartType::Shaft, { 0, 2, 8 }, up)
            && Place(c, PartType::Gearbox, { 0, 3, 8 }) && Place(c, PartType::HandCrank, { -1, 3, 8 }, alongX);
        if (!c.check("the millstone, its shaft, the gearbox and the crank placed", ok)) return true;
        const int mill = plot->partAt({ 0, 0, 8 }), crank = plot->partAt({ -1, 3, 8 });
        Factory::Insert(c.body->id, plot->id, mill, "stone", 40);
        c.camera(At(-2.0f, 3.0f, -5.0f), At(1.0f, -0.5f, 2.0f), Vector3Normalize(At(0.0f, 0.0f))); // the trench and the machines
        Factory::Sync(c.world);
        for (int k = 0; k < 40; ++k) { Factory::Crank(c.body->id, plot->id, crank); Factory::WaitTicks(5); }
        const RotationView v = Factory::RotationAt(c.body->id, plot->id, crank);
        c.gt("the crank turned the line (RPM)", v.rpm, 0.0); // the hash covers a working line, not an idle one
        c.hashState("cranked", { "edits", "factory" });
        Fluid::AddSource(*c.body, At(-7.5f, -2.2f), 0, 1.0f, 1.5f); // 1 m^3/s at the trench's head
        return true;
    });
    s.step("15 s of a spring in the trench; hash the water", 40.0, [](ScenarioCtx& c) {
        c.waitingFor(TextFormat("t=%.1f s of 15", c.t));
        if (c.t < 15.0) return false;
        c.hashState("flowing", { "water" });
        return true;
    });
}
