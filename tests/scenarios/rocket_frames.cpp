// Plan 11 F5 (#1457) gates for the rocket in body frames. Both build and board --test-rocket's rocket on the planet, launch
// it, then place it coasting (Rocket_TestSetState: engine off, autopilot off):
//   rocket/orbit-closure: a near-circular orbit 45 km over Terra, 10 orbits at 100x warp; the drift of its lowest and
//     highest points from the first orbit to the 10th (Luna's pull makes it ~400 m; the gate only guards against a blow-up).
//   rocket/soi-trip: a flyby of Luna, from 1.3 SOI radii out, passing ~2 radii from its centre: it enters Luna's SOI
//     (Centred(Terra) -> Centred(Luna)) and leaves again. Each frame, the global path's second difference (a smooth path:
//     ~a dt^2) is measured; at every frame switch it must stay under 1 mm (a switch re-expresses the same state).
#include "scenario.h"
#include "buildtool.h"
#include "inventory.h"
#include "player.h"
#include "plots.h"
#include "rocket.h"
#include "station.h"
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
    DVec3 c1{}, c2{}, dd1{}; int have = 0; // the last two global positions and the last second difference
    FrameRef lastFrame{};
    int switches = 0, frames = 0; bool afterSwitch = false;
    double worstSwitch = 0.0, worstSmooth = 0.0, closest = 1e30;
    // Orbit drift: the radius range over the first orbit and over the last.
    double firstMin = 1e30, firstMax = 0.0, lastMin = 1e30, lastMax = 0.0, doneMin = 0.0, doneMax = 0.0; int lastOrbit = -1;
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

// Each frame, the jump in the global path: the change of its second difference from the frame before (the third
// difference). A smooth path's second difference is a dt^2 (2-4 mm a frame at 5x warp), but it changes only by its jerk
// dt^3; a one-frame jump j shows as j here. Measured at frame switches (gated) and elsewhere (reported).
void Track(ScenarioCtx& c, const DVec3& com, const FrameRef& f) {
    const bool switched = S.have >= 1 && f != S.lastFrame;
    if (switched) ++S.switches;
    if (S.have >= 3) {
        const DVec3 dd = com - S.c1 * 2.0 + S.c2;
        const double jump = (dd - S.dd1).length();
        if (switched || S.afterSwitch) S.worstSwitch = std::max(S.worstSwitch, jump); // a jump shows on two frames
        else S.worstSmooth = std::max(S.worstSmooth, jump);
        S.dd1 = dd;
    } else if (S.have == 2) S.dd1 = com - S.c1 * 2.0 + S.c2;
    S.afterSwitch = switched;
    S.c2 = S.c1; S.c1 = com; S.have = std::min(S.have + 1, 3);
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
        Track(c, com, f);
        const double orbits = (c.T - S.t0) * Rocket_Warp() / S.period;
        if (c.T - S.t0 > 1.0) {
            const double d = (com - t.worldPos).length();
            S.rMin = std::min(S.rMin, d); S.rMax = std::max(S.rMax, d);
            const int k = (int)orbits;
            if (k == 0) { S.firstMin = std::min(S.firstMin, d); S.firstMax = std::max(S.firstMax, d); }
            if (k != S.lastOrbit) { if (S.lastOrbit >= 0) { S.doneMin = S.lastMin; S.doneMax = S.lastMax; } S.lastOrbit = k; S.lastMin = 1e30; S.lastMax = 0.0; } // doneMin/Max: the last whole orbit
            S.lastMin = std::min(S.lastMin, d); S.lastMax = std::max(S.lastMax, d);
        }
        c.waitingFor(TextFormat("10 orbits (%.1f so far)", orbits));
        if (orbits < 10.0) return false;
        Rocket_SetWarp(0);
        c.note(TextFormat("%.1f orbits of %.0f s at %.0f km: distance from the centre %.1f-%.1f m (band %.2f m); %d frame switches",
                          orbits, S.period, (r - R) / 1000.0, S.rMin, S.rMax, S.rMax - S.rMin, S.switches));
        // The start isn't exactly circular (the other bodies pull, the ring), so the band is reported; the gate is the drift:
        // the 10th orbit's lowest and highest points against the first's.
        const double drift = std::max(fabs(S.doneMin - S.firstMin), fabs(S.doneMax - S.firstMax));
        c.note(TextFormat("first orbit %.1f-%.1f m, last %.1f-%.1f m from the centre: drift %.2f m", S.firstMin, S.firstMax, S.doneMin, S.doneMax, drift));
        c.metric("radius band", S.rMax - S.rMin, "m", true);
        c.metric("orbit drift over 10 orbits", drift, "m");
        // ~400 m on seed 42 (frames on 397.76; the global path, before F8b, 399.05): Luna's pull precesses the orbit over its
        // ~4 h, so the plan's 10 m bound (two bodies only) doesn't apply. This guards against a blow-up.
        c.lt("the orbit's drift over 10 orbits (m; physics, ~400 m with Luna)", drift, 600.0);
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
        c.note(TextFormat("closest %.1f km (Luna R %.1f km, SOI %.1f km); %d frame switches; the path's jump (third difference): worst %.4f mm at a switch, %.4f mm elsewhere",
                          S.closest / 1000.0, Rm / 1000.0, soi / 1000.0, S.switches, S.worstSwitch * 1000.0, S.worstSmooth * 1000.0));
        c.metric("worst jump at a switch", S.worstSwitch * 1000.0, "mm");
        c.check("entered Luna's SOI and left it", S.inLuna && S.leftLuna);
        c.check("the frame switched entering and leaving", S.switches >= 2, TextFormat("%d switches", S.switches));
        c.lt("the global path's jump at a frame switch (m)", S.worstSwitch, 0.001);
        return true;
    });
}

// rocket/station-hold (F5 stations): a station created 20 m beside the coasting rocket on its 45 km circular orbit, at
// the same velocity; for 30 s both fall freely, so their separation must hold (the tidal divergence over 20 m is far under
// a centimetre). Inside the station's small SOI a framed rocket flies in Fixed(station): if that frame's own acceleration
// is missing (stations fly Kepler elements, not the frames' orbit rule), the rocket feels the free fall twice and drifts
// metres away. Gate: the separation changes by under 1 m.
static void StationHold(ScenarioDef& s) {
    Common(s);
    s.step("coast 30 s beside a station", 60.0, [](ScenarioCtx& c) {
        Planetoid& t = *S.terra;
        const double R = t.radius * t.scale, r = R + 45000.0, mu = Mu(t);
        static Planetoid* st = nullptr;
        static double sep0 = 0.0;
        if (!S.placed) {
            const DVec3 u(0.0, 1.0, 0.0), e(1.0, 0.0, 0.0);
            const DVec3 pos = t.worldPos + u * r;
            const DVec3 v = DVec3(t.linearVelocity) + e * sqrt(mu / r);
            if (!Rocket_TestSetState(pos, v.toFloat())) { c.fail("not flying"); return true; }
            st = Station_Create(c.world, t, pos + DVec3(0.0, 0.0, 20.0), v.toFloat(), { 0, 0, 0, 1 }, 4.0f);
            if (!st) { c.fail("no station"); return true; }
            sep0 = 20.0;
            S.placed = true; S.t0 = c.T;
            return false;
        }
        DVec3 com; Vector3 vel; FrameRef f;
        if (!Rocket_TestState(com, vel, &f)) { c.fail("the rocket stopped flying"); return true; }
        Track(c, com, f);
        const double sep = (com - st->worldPos).length();
        S.rMax = std::max(S.rMax, fabs(sep - sep0));
        { // a jump or a drift? the separation at 1 s and 5 s, and the path's jumps (Track: at frame switches and elsewhere)
            static bool n1 = false, n5 = false;
            if (c.T - S.t0 < 0.5) n1 = n5 = false;
            if (!n1 && c.T - S.t0 >= 1.0) { n1 = true; c.note(TextFormat("at 1 s: %.3f m apart; %d frame switches so far", sep, S.switches)); }
            if (!n5 && c.T - S.t0 >= 5.0) { n5 = true; c.note(TextFormat("at 5 s: %.3f m apart", sep)); }
        }
        c.waitingFor(TextFormat("30 s beside the station (%.2f m apart)", sep));
        if (c.T - S.t0 < 30.0) return false;
        c.note(TextFormat("after 30 s: %.3f m apart (from %.1f), the separation changed by up to %.3f m; frame %s(%s)", sep, sep0, S.rMax,
                          f.kind == FrameKind::Fixed ? "Fixed" : "Centred", f.body == st->netId ? "the station" : f.body == t.netId ? "Terra" : "?"));
        { // Which frame the SOI rule gives here, and the station's own numbers (stations as Fixed(station) frames, F5)
            const std::vector<Planetoid*> all = Planetoid::All(); // (a copy: one vector for begin and end)
            const FrameRef pick = Frames_Choose(all, com, f);
            c.note(TextFormat("Frames_Choose here: %s of netId %llu (the station is netId %llu: SOI %.1f m, radius %.1f, kind %d, in All(): %s)",
                              pick.kind == FrameKind::Fixed ? "Fixed" : "Centred", (unsigned long long)pick.body, (unsigned long long)st->netId,
                              Frames_Soi(*st), st->radius * st->scale, (int)st->kind,
                              std::find(all.begin(), all.end(), st) != all.end() ? "yes" : "NO"));
        }
        c.note(TextFormat("the path's jump (third difference): worst %.4f mm at a frame switch (%d switches), %.4f mm elsewhere", S.worstSwitch * 1000.0, S.switches, S.worstSmooth * 1000.0));
        c.metric("separation change", S.rMax, "m");
        c.lt("the separation from a free-falling station changed by (m)", S.rMax, 1.0);
        return true;
    });
}

// rocket/frame-save (Plan 11 F7a, #2074), frames on: coasting in a circular orbit 45 km over Terra, the world is saved and
// loaded mid-flight. The rocket's frame state (its frame, double position and velocity) must come back to 1e-6 m and
// 1e-6 m/s (the save's rocket_frame line), and the orbit must go on as before: after 10 s more its radius is within 5 m
// of the orbit's.
static void FrameSave(ScenarioDef& s) {
    Common(s);
    s.step("coast, save and load mid-flight", 120.0, [](ScenarioCtx& c) {
        static uint64_t terraNet = 0;
        if (!S.placed) terraNet = S.terra->netId;
        Planetoid* tp = Planetoid::findNet(terraNet); // (a load makes the bodies again: S.terra dangles after it)
        if (!tp) { c.fail("no Terra after the load"); return true; }
        Planetoid& t = *tp;
        const double R = t.radius * t.scale, r = R + 45000.0, mu = Mu(t);
        static FramePos before{};
        static bool loaded = false;
        if (!S.placed) {
            const DVec3 u(0.0, 1.0, 0.0), e(1.0, 0.0, 0.0);
            const DVec3 pos = t.worldPos + u * r;
            const DVec3 v = DVec3(t.linearVelocity) + e * sqrt(mu / r);
            if (!Rocket_TestSetState(pos, v.toFloat())) { c.fail("not flying"); return true; }
            loaded = false;
            S.placed = true; S.t0 = c.T;
            return false;
        }
        DVec3 com; Vector3 vel; FramePos fp;
        if (!Rocket_TestState(com, vel, nullptr, &fp)) { c.fail("the rocket stopped flying"); return true; }
        if (!loaded) {
            c.waitingFor("3 s of coasting on the frame path");
            if (c.T - S.t0 < 3.0 || fp.frame.body == 0) return false;
            before = fp;
            if (!c.check("the save", SaveWorld(c.world))) return true;
            if (!c.check("the load", LoadWorldSave(c.world), LastLoadError())) return true;
            FramePos after;
            if (!c.check("flying after the load", Rocket_TestState(com, vel, nullptr, &after))) return true;
            c.check("the same frame after the load", after.frame == before.frame,
                    TextFormat("netId %llu kind %d, was %llu kind %d", (unsigned long long)after.frame.body, (int)after.frame.kind, (unsigned long long)before.frame.body, (int)before.frame.kind));
            c.lt("frame position moved by the load (m)", (after.pos - before.pos).length(), 1e-6);
            c.lt("frame velocity changed by the load (m/s)", (after.vel - before.vel).length(), 1e-6);
            loaded = true; S.t0 = c.T;
            return false;
        }
        c.waitingFor("10 s of coasting after the load");
        if (c.T - S.t0 < 10.0) return false;
        const double rNow = (com - t.worldPos).length();
        S.terra = tp;
        c.note(TextFormat("10 s after the load: %.3f m from Terra's centre (the orbit's radius %.3f)", rNow, r));
        c.lt("the orbit's radius after the load is off by (m)", fabs(rNow - r), 5.0);
        return true;
    });
}

// (F8b #2113: the -pose twins folded in; frames are the only path.)
SCENARIO("rocket/frame-save") { FrameSave(s); s.tags = { "frames", "saves", "local" }; }
SCENARIO("rocket/orbit-closure") { OrbitClosure(s); }
SCENARIO("rocket/soi-trip") { SoiTrip(s); }
SCENARIO("rocket/station-hold") { StationHold(s); }
