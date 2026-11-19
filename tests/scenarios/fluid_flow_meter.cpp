// fluid/flow-meter (#870, moved from --test-flow-meter, #392): a spring of 1 m^3/s at the head of a 20 m trench that ends in a
// deep pit; once the stream runs, Fluid::FlowAt across the middle of the trench reads ~1 m^3/s along it. Steady flow (#674):
// on the GPU the trench's bricks turn steady; a pour at 45 s wakes them (fewer steady). PASS as the flag: 0.6-1.4 m^3/s
// along the trench at 60 s, and (GPU mode) steady bricks before the pour and fewer after.
#include "scenario.h"
#include "fluid.h"
#include "plots.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>

namespace {
struct State {
    Plot* plot = nullptr;
    int steadyBefore = 0, steadyAfter = 1 << 30;
    bool poured = false;
    double lastPrint = -10.0;
    Vector3 at(float x, float y) const { return plot->toBody({ x, y, 0.0f }); }
} S;
}

SCENARIO("fluid/flow-meter") {
    s.tags = { "fluid", "water", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 150.0;
    s.legacyFlag = "--seed 42 --test-flow-meter";
    // Lockstep fluids never mark bricks steady (fluid.cpp: a steady brick is read back one tick in several), and #674's check
    // is about steady bricks: this one runs the fluid as the game does.
    s.env = { { "RAYLIBGAME_FLUID_LOCKSTEP", "0" } };

    s.step("dig", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (!WorldIdle()) return false;
        const Vector3 d = c.spawnDirWorld(); // (as the flag: StartSpawnDir, the start's direction turned into the world frame)
        float gh = c.field->elevation(d) + 8.0f;
        while (gh > c.field->elevation(d) - 8.0f && c.body->densityAt(Vector3Scale(d, c.field->radius + gh)) < ISO_LEVEL) gh -= 0.1f;
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        S.plot = c.body->plots->plotFor(Vector3Scale(d, c.field->radius + gh), 0.0f, true);
        std::vector<TerrainEdit> cut;
        for (float x = -10.0f; x <= 10.0f; x += 0.8f) { cut.push_back({ S.at(x, -1.0f), 1.4f, -1.0f, true }); cut.push_back({ S.at(x, 0.2f), 1.4f, -1.0f, true }); }
        for (int k = 0; k < 4; ++k) cut.push_back({ S.at(14.0f, -3.0f - 2.5f * k), 4.0f, -1.0f, true }); // the pit: ~600 m^3
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        Fluid::AddSource(*c.body, S.at(-9.5f, -2.2f), 0, 1.0f, 1.5f);
        c.note("trench and pit dug, a 1 m^3/s spring at the trench's head");
        return true;
    });
    s.step("stream", 70.0, [](ScenarioCtx& c) {
        const Vector3 up = Vector3Normalize(S.at(0.0f, 0.0f));
        c.camera(Vector3Add(S.at(0.0f, 8.0f), Vector3Scale(S.plot->north, 8.0f)), S.at(0.0f, -1.5f), up);
        const Vector3 q = Fluid::FlowAt(*c.body, S.at(0.0f, -1.8f), 2.0f); // asked every frame: the readout smooths
        const Fluid::Stats fs = Fluid::GetStats();
        if (c.t - S.lastPrint >= 10.0) { S.lastPrint = c.t; c.note("t=" + std::to_string((int)c.t) + " s: " + std::to_string(Vector3DotProduct(q, S.plot->east)) + " m^3/s along the trench, " + std::to_string(fs.steadyBricks) + " steady bricks"); }
        if (c.t > 30.0 && c.t < 45.0) S.steadyBefore = std::max(S.steadyBefore, fs.steadyBricks);
        if (c.t > 45.0 && !S.poured) { S.poured = true; Fluid::Pour(*c.body, S.at(0.0f, -1.5f), 3.0f); }
        if (c.t > 45.5 && c.t < 48.0) S.steadyAfter = std::min(S.steadyAfter, fs.steadyBricks);
        c.waitingFor("60 s of stream");
        if (c.t <= 60.0) return false;
        const float along = Vector3DotProduct(q, S.plot->east);
        const bool gpu = Fluid::GpuMode();
        c.metric("flow along the trench", along, "m^3/s", false);
        c.metric("steady bricks before the pour", S.steadyBefore, "", false);
        c.metric("steady bricks after the pour", S.steadyAfter == (1 << 30) ? -1 : S.steadyAfter);
        c.between("flow along the trench (m^3/s) from a 1 m^3/s spring", along, 0.6, 1.4);
        if (gpu) c.check("steady bricks woken by a pour (#674)", S.steadyBefore > 0 && S.steadyAfter < S.steadyBefore,
                         std::to_string(S.steadyBefore) + " before, " + std::to_string(S.steadyAfter) + " after");
        else c.note("CPU mode: no steady bricks to check");
        return true;
    });
}
