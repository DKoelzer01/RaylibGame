// fluid/flow-meter (#870, moved from --test-flow-meter, #392): a spring of 1 m^3/s at the head of a 20 m trench that ends in a
// deep pit; once the stream runs, Fluid::FlowAt across the middle of the trench reads ~1 m^3/s along it. Steady flow (#674):
// on the GPU the trench's bricks turn steady; a pour into a steady one (found from 30 s on, #1131) wakes it. PASS: 0.6-1.4 m^3/s
// along the trench at 60 s, and (GPU mode) a steady brick found and woken by the pour.
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
    bool poured = false, woke = false;
    Vector3 pourAt{};
    double pourT = 0.0;
    double lastPrint = -10.0;
    // (#1407 diagnostics) the last window's counters: fluid ticks, GPU dispatches, flux reads, wall time, water in the pit
    long long ticks = -1, dispatches = 0, reads = 0; double wall = 0.0, pit = 0.0, simT = 0.0;
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
        // (#1407) Is the stream slow on the runner, or only the readout? Every 10 s: fluid ticks per sim and per wall second
        // (40 = the water at the game's clock), dispatches and flux reads per tick, and the water reaching the pit per sim and
        // per fluid second (a mass count FlowAt has no part in; the pour at 30-55 s adds 3 m^3 to it once).
        if (S.ticks < 0 || c.t - S.simT >= 10.0) {
            const double pit = Fluid::VolumeIn(*c.body, S.at(14.0f, -6.75f), S.plot->east, S.plot->north, S.plot->up, { 4.5f, 4.5f, 6.5f });
            if (S.ticks >= 0) {
                const double dt = c.t - S.simT, dw = GetTime() - S.wall, dk = (double)(fs.ticks - S.ticks);
                c.note(TextFormat("t=%d s: %.1f fluid ticks per sim s, %.1f per wall s; %.2f dispatches and %.2f flux reads a tick; the pit gained %.3f m^3 per sim s, %.3f per fluid s",
                                  (int)c.t, dk / dt, dk / std::max(1e-6, dw), (fs.dispatches - S.dispatches) / std::max(1.0, dk),
                                  (fs.flowReads - S.reads) / std::max(1.0, dk), (pit - S.pit) / dt, (pit - S.pit) / std::max(1e-6, dk * 0.025)));
                c.metric("fluid ticks per sim second", dk / dt, "", false);
            }
            S.ticks = fs.ticks; S.dispatches = fs.dispatches; S.reads = fs.flowReads; S.wall = GetTime(); S.pit = pit; S.simT = c.t;
        }
        if (c.t > 30.0 && !S.poured) S.steadyBefore = std::max(S.steadyBefore, fs.steadyBricks);
        // (#1131) Pour into a brick that is steady, found along the stream from 30 s on (by 55 s at the latest), and check
        // that this brick wakes: the global count before and after depended on which bricks happened to be steady (the
        // pour's fixed spot was not always one of them: "5 before, 5 after").
        if (c.t > 30.0 && c.t < 55.0 && !S.poured && Fluid::GpuMode())
            for (float x = -8.0f; x <= 10.0f && !S.poured; x += 1.0f) {
                const Vector3 p = S.at(x, -1.5f);
                if (Fluid::SteadyAt(*c.body, p)) { S.poured = true; S.pourAt = p; S.pourT = c.t; Fluid::Pour(*c.body, p, 3.0f); c.note("poured into a steady brick at x = " + std::to_string((int)x) + " m"); }
            }
        if (S.poured && c.t > S.pourT + 0.3 && c.t < S.pourT + 3.0) {
            S.steadyAfter = std::min(S.steadyAfter, fs.steadyBricks);
            if (!Fluid::SteadyAt(*c.body, S.pourAt)) S.woke = true;
        }
        c.waitingFor("60 s of stream");
        if (c.t <= 60.0 || (S.poured && c.t < S.pourT + 3.0)) return false;
        const float along = Vector3DotProduct(q, S.plot->east);
        const bool gpu = Fluid::GpuMode();
        c.metric("flow along the trench", along, "m^3/s", false);
        c.metric("steady bricks before the pour", S.steadyBefore, "", false);
        c.metric("steady bricks after the pour", S.steadyAfter == (1 << 30) ? -1 : S.steadyAfter);
        c.between("flow along the trench (m^3/s) from a 1 m^3/s spring", along, 0.6, 1.4);
        if (gpu) {
            c.check("a steady brick along the stream (#674)", S.poured, "none steady between 30 and 55 s");
            c.check("the steady brick a pour lands in wakes (#674)", S.woke, std::to_string(S.steadyBefore) + " steady at most before, " +
                    std::to_string(S.steadyAfter == (1 << 30) ? -1 : S.steadyAfter) + " fewest after");
        } else c.note("CPU mode: no steady bricks to check");
        return true;
    });
}
