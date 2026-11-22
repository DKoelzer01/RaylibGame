// Plan 11 F5 (#1457) gates for the rocket in body frames. Both build and board --test-rocket's rocket on the planet, launch
// it, then place it coasting (Rocket_TestSetState: engine off, autopilot off):
//   rocket/orbit-closure: a circular orbit 45 km over Terra, 10 orbits at 100x warp; the band of its distance from Terra's
//     centre (max - min) must stay within 10 m (a float global velocity summed every 1/120 s drifted; the frame path holds
//     the state in doubles).
//   rocket/soi-trip: a flyby of Luna, from 1.3 SOI radii out, passing ~2 radii from its centre: it enters Luna's SOI
//     (Centred(Terra) -> Centred(Luna)) and leaves again. Each frame, the global path's second difference (a smooth path:
//     ~a dt^2) is measured; at every frame switch it must stay under 1 mm (a switch re-expresses the same state).
// The -pose variants run with RAYLIBGAME_FRAMES=1; without it the switch counts are 0 (the global path has no frames).
#include "scenario.h"
#include "buildtool.h"
#include "inventory.h"
#include "player.h"
#include "plots.h"
#include "rocket.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
struct RState {
    bool built = false, placed = false;
    double t0 = 0.0;
    Planetoid* terra = nullptr; Planetoid* luna = nullptr;
    double rMin = 1e30, rMax = 0.0, period = 0.0;
    DVec3 c1{}, c2{}; int have = 0;        // the last two global positions (second differences)
    FrameRef lastFrame{};
    int switches = 0, frames = 0;
    double worstSwitch = 0.0, worstSmooth = 0.0, closest = 1e30;
    bool inLuna = false, leftLuna = false;
};
RState S;

double Mu(const Planetoid& b) { const double R = b.gravityRadius > 0.0f ? b.gravityRadius : b.radius * b.scale; return b.surfaceGravity * R * R; }

// --test-rocket's rocket, built on a plot 14 m from the start and boarded (OnPad).
bool BuildAndBoard(ScenarioCtx& c) {
    Planetoid& p = *c.body;
    const Vector3 home = c.spawnDirLocal();
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 }));
    const Vector3 pad = Vector3Normalize(Vector3Add(home, Vector3Scale(t1, 14.0f / c.field->radius)));
    if (!p.plots) p.plots = std::make_unique<PlotSet>();
    Plot* plot = p.plots->plotFor(Vector3Scale(pad, c.field->radius + c.field->elevation(pad)), 0.0f, true);
    bool ok = BuildTool_Place(p, *plot, { PartType::RocketEngine, { 0, 0, 0 }, 0 }, false);
    const Int3 legs[4] = { { -1, 0, -1 }, { 2, 0, -1 }, { -1, 0, 2 }, { 2, 0, 2 } };
    for (const Int3& l : legs) ok &= BuildTool_Place(p, *plot, { PartType::LandingLeg, l, 0 }, false);
    for (int k = 0; k < 4; ++k) ok &= BuildTool_Place(p, *plot, { PartType::FuelTank, { 0, 2 + 2 * k, 0 }, 0 }, false);
    ok &= BuildTool_Place(p, *plot, { PartType::Guidance, { 0, 10, 0 }, 0 }, false);
    ok &= BuildTool_Place(p, *plot, { PartType::Cockpit, { 0, 11, 0 }, 0 }, false);
    inventory.add("rocket_fuel", 60);
    Rocket_LoadFuel(*plot, 60);
    inventory.add("oxygen_canister", 1);
    Rocket_LoadAir(*plot, 1);
    return ok && Rocket_Board(c.world, p, *plot);
}

void Common(ScenarioDef& s) {
    s.tags = { "frames", "local" }; // local: rocket flights (the regression runs the rocket tests on the PC)
    s.world = ScenarioWorld::Planet;
    s.timeout = 420.0;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 1.0 && WorldIdle(); });
    s.step("build, board and launch", 30.0, [](ScenarioCtx& c) {
        if (!S.built) {
            S = RState{};
            for (auto& obj : c.world.objects)
                if (auto* q = dynamic_cast<Planetoid*>(obj.get())) { if (!S.terra && q->field) S.terra = q; if (q->kind == BodyKind::Moon && q->orbitParent) S.luna = q; }
            S.terra = c.body;
            if (!BuildAndBoard(c)) { c.fail("the rocket wasn't built or boarded"); return true; }
            S.built = true;
            Rocket_SetAutopilot(RocketAutopilot::Ascent, 45000.0); // launches it (state Flying)
            return false;
        }
        DVec3 com; Vector3 vel;
        c.waitingFor("the rocket flying");
        return c.t > 1.0 && Rocket_TestState(com, vel);
    });
}

// Each frame: the second difference of the global path, and frame switches.
void Track(ScenarioCtx& c, const DVec3& com, const FrameRef& f) {
    if (S.have >= 2) {
        const double dd = (com - S.c1 * 2.0 + S.c2).length();
        if (f != S.lastFrame) { ++S.switches; S.worstSwitch = std::max(S.worstSwitch, dd); }
        else S.worstSmooth = std::max(S.worstSmooth, dd);
    }
    S.c2 = S.c1; S.c1 = com; S.have = std::min(S.have + 1, 2);
    S.lastFrame = f;
    ++S.frames;
    (void)c;
}
}

static void OrbitClosure(ScenarioDef& s) {
    Common(s);
    s.step("10 orbits at 45 km, 100x warp", 300.0, [](ScenarioCtx& c) {
        Planetoid& t = *S.terra;
        const double R = t.radius * t.scale, r = R + 45000.0, mu = Mu(t);
        if (!S.placed) {
            const DVec3 u(0.0, 1.0, 0.0), e(1.0, 0.0, 0.0);      // over the pole, moving along X: a polar orbit
            const DVec3 pos = t.worldPos + u * r;
            const DVec3 v = DVec3(t.linearVelocity) + e * sqrt(mu / r);
            if (!Rocket_TestSetState(pos, v.toFloat())) { c.fail("not flying"); return true; }
            Rocket_SetWarp(6);
            S.period = 2.0 * PI * sqrt(r * r * r / mu);
            S.placed = true; S.t0 = c.T;
            return false;
        }
        DVec3 com; Vector3 vel; FrameRef f;
        if (!Rocket_TestState(com, vel, &f)) { c.fail("the rocket stopped flying (crashed or landed)"); return true; }
        if (c.T - S.t0 > 1.0) { const double d = (com - t.worldPos).length(); S.rMin = std::min(S.rMin, d); S.rMax = std::max(S.rMax, d); }
        Track(c, com, f);
        const double orbits = (c.T - S.t0) * Rocket_Warp() / S.period;
        c.waitingFor(TextFormat("10 orbits (%.1f so far)", orbits));
        if (orbits < 10.0) return false;
        Rocket_SetWarp(0);
        c.note(TextFormat("%.1f orbits of %.0f s at %.0f km: distance from the centre %.1f-%.1f m (band %.2f m); %d frame switches (%s)",
                          orbits, S.period, (r - R) / 1000.0, S.rMin, S.rMax, S.rMax - S.rMin, S.switches, Frames_On() ? "frames on" : "frames off"));
        c.metric("radius band", S.rMax - S.rMin, "m");
        c.lt("the orbit's radius band after 10 orbits (m)", S.rMax - S.rMin, 10.0);
        return true;
    });
}

static void SoiTrip(ScenarioDef& s) {
    Common(s);
    s.step("a flyby of Luna", 300.0, [](ScenarioCtx& c) {
        if (!S.luna) { c.fail("no moon"); return true; }
        Planetoid& m = *S.luna;
        const double soi = Frames_Soi(m), Rm = m.radius * m.scale;
        if (!S.placed) {
            // From 1.3 SOI radii out, aimed to pass 2 Luna radii from its centre at 300 m/s relative.
            const DVec3 toTerra = (S.terra->worldPos - m.worldPos) * (1.0 / (S.terra->worldPos - m.worldPos).length());
            const DVec3 side(-toTerra.z, 0.0, toTerra.x);
            const DVec3 start = m.worldPos + toTerra * (1.3 * soi) + side * (2.0 * Rm);
            const DVec3 v = DVec3(m.linearVelocity) - toTerra * 300.0;
            if (!Rocket_TestSetState(start, v.toFloat())) { c.fail("not flying"); return true; }
            Rocket_SetWarp(2); // 5x
            S.placed = true; S.t0 = c.T;
            return false;
        }
        DVec3 com; Vector3 vel; FrameRef f;
        if (!Rocket_TestState(com, vel, &f)) { c.fail("the rocket stopped flying (crashed into Luna?)"); return true; }
        Track(c, com, f);
        const double d = (com - m.worldPos).length();
        S.closest = std::min(S.closest, d);
        if (d < soi) S.inLuna = true;
        if (S.inLuna && d > 1.3 * soi) S.leftLuna = true;
        c.waitingFor(TextFormat("through Luna's SOI and out (%.0f km from it)", d / 1000.0));
        if (!S.leftLuna && c.T - S.t0 < 280.0) return false;
        Rocket_SetWarp(0);
        c.note(TextFormat("closest %.1f km (Luna R %.1f km, SOI %.1f km); %d frame switches; the path's second difference: worst %.4f mm at a switch, %.4f mm elsewhere (%s)",
                          S.closest / 1000.0, Rm / 1000.0, soi / 1000.0, S.switches, S.worstSwitch * 1000.0, S.worstSmooth * 1000.0, Frames_On() ? "frames on" : "frames off"));
        c.metric("worst jump at a switch", S.worstSwitch * 1000.0, "mm");
        c.check("entered Luna's SOI and left it", S.inLuna && S.leftLuna);
        if (Frames_On()) {
            c.check("the frame switched entering and leaving", S.switches >= 2, TextFormat("%d switches", S.switches));
            c.lt("the global path's jump at a frame switch (m)", S.worstSwitch, 0.001 + S.worstSmooth);
        }
        return true;
    });
}

SCENARIO("rocket/orbit-closure") { OrbitClosure(s); }
SCENARIO("rocket/orbit-closure-pose") { s.env = { { "RAYLIBGAME_FRAMES", "1" } }; OrbitClosure(s); }
SCENARIO("rocket/soi-trip") { SoiTrip(s); }
SCENARIO("rocket/soi-trip-pose") { s.env = { { "RAYLIBGAME_FRAMES", "1" } }; SoiTrip(s); }
