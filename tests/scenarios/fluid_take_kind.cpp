// fluid/take-kind: Fluid::Take of a non-water fluid removes it, in GPU mode, from a thin moving film. Magic's alchemy
// fizzle (#1398) took essence with Take(body, outlet, amount, 5 m, kind) every 0.5 s and the published amount never fell
// (0.300 of 0.30 m^3); a bucket taking oil goes the same way. 0.3 m^3 of crude oil is poured on the sandbox's flat ground
// (it spreads into a film); from 3 s on, Take(0.05 m^3, 5 m, oil) every 0.5 s for 6 s, as alchemy does. Checked: the
// published amount (FluidAt, what the take itself asks first) and the simulation's own (Fluid::VolumeIn) both fall, by about
// what was taken. The default (GPU) path; -cpu runs the CPU kernel.
#include "scenario.h"
#include "fluid.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
struct State {
    Vector3 home{}, ax{}, ay{}, spot{};
    int kind = -1;
    double nextTake = 0.0;
    float taken = 0.0f, pub0 = 0.0f, sim0 = 0.0f;
    int takes = 0;
} S;

double SimOil(ScenarioCtx& c) { // the oil the simulation holds within 6 m of the spot (its own bricks)
    const Vector3 up = Vector3Normalize(S.spot);
    return Fluid::VolumeIn(*c.body, S.spot, S.ax, S.ay, up, { 6.0f, 6.0f, 2.0f }); // all kinds: the sandbox has no other fluid
}

void Build(ScenarioDef& s) {
    s.tags = { "fluid" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("pour 0.3 m^3 of oil on flat ground", 2.0, [](ScenarioCtx& c) {
        S = State{};
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.spot = Vector3Scale(S.home, c.field->radius + c.field->elevation(S.home) + 0.3f);
        S.kind = Fluid::KindByName("crude_oil");
        if (S.kind < 0) { c.fail("no crude_oil fluid"); return true; }
        Fluid::Pour(*c.body, S.spot, 0.3f, S.kind);
        c.camera(Vector3Add(S.spot, Vector3Add(Vector3Scale(S.home, 4.0f), Vector3Scale(S.ay, -6.0f))), S.spot, S.home);
        return true;
    });
    s.step("take it back as alchemy does", 20.0, [](ScenarioCtx& c) {
        c.camera(Vector3Add(S.spot, Vector3Add(Vector3Scale(S.home, 4.0f), Vector3Scale(S.ay, -6.0f))), S.spot, S.home);
        c.waitingFor("3 s spreading, then 6 s of takes");
        if (c.t < 3.0) return false;
        if (S.takes == 0 && S.nextTake == 0.0) { S.pub0 = Fluid::FluidAt(*c.body, S.spot, 5.0f, S.kind); S.sim0 = (float)SimOil(c); S.nextTake = c.t; }
        if (c.t >= S.nextTake && c.t < 9.0) {
            S.nextTake = c.t + 0.5;
            S.taken += Fluid::Take(*c.body, S.spot, 0.05f, 5.0f, S.kind);
            ++S.takes;
        }
        if (c.t < 10.0) return false; // a second for the last takes to reach the GPU and come back
        const float pub1 = Fluid::FluidAt(*c.body, S.spot, 5.0f, S.kind);
        const float sim1 = (float)SimOil(c);
        c.note(TextFormat("%d takes asked for %.3f m^3; published oil %.3f -> %.3f m^3, the simulation's %.3f -> %.3f m^3 (%s)",
                          S.takes, S.taken, S.pub0, pub1, S.sim0, sim1, Fluid::GpuMode() ? "GPU" : "CPU"));
        c.metric("published oil removed", S.pub0 - pub1, "m^3", false);
        c.metric("simulated oil removed", S.sim0 - sim1, "m^3", false);
        c.gt("the takes asked for some oil (m^3)", S.taken, 0.05);
        c.gt("the published amount fell by most of what was taken (m^3)", S.pub0 - pub1, 0.7f * S.taken);
        c.gt("the simulation's own amount fell by most of what was taken (m^3)", S.sim0 - sim1, 0.7f * S.taken);
        return true;
    });
}
}

SCENARIO("fluid/take-kind") { Build(s); }
SCENARIO("fluid/take-kind-cpu") { s.env = { { "FLUID_CPU", "1" } }; Build(s); }
