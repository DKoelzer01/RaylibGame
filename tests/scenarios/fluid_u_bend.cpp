// fluid/u-bend (#1168): communicating vessels. Two pits 4 m deep (4 x 2 m each) joined by a tunnel 1 m high and 2 m wide
// under 3 m of solid roof, 16 m long; 80 m^3 is poured into the near pit. Hydrostatic pressure must push the water through
// the tunnel and up the far pit until both stand at the same level: below the tunnel's top the pits and the tunnel hold
// 16 + 32 = 48 m^3, the other 32 m^3 stand 2 m deep over 16 m^2 of pits, so both levels would end 1 m below the ground.
// The carved space holds ~12 m^3 more than those boxes (measured, #1168: with the volume conserved, the levels settle at
// ~-1.75 m), so the gates are what holds whatever the carve's exact size: the two legs level, the tunnel full, no water lost
// or made; the level against the box arithmetic is reported. The playtest (2026-11-20, bfbec4c): "water cannot flow up a
// U-bend" - it stopped at the bottom or the near leg.
#include "scenario.h"
#include "fluid.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
constexpr float kPour = 80.0f;      // m^3
constexpr float kFloor = -4.0f;     // m under the ground: the pits' and the tunnel's floor
constexpr float kTunnelTop = -3.0f; // the tunnel's roof underside
constexpr float kExpected = kTunnelTop + (kPour - 48.0f) / 16.0f; // -1 m: the common level
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f, poured = 0.0f;
    double vol0 = 0.0; // the simulated volume once the pour is done (m^3)
    double nextPour = 0.0;
    int k = 0;
    Vector3 pt(float x, float y, float h) const {
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    TerrainEdit box(float x0, float x1, float y0, float y1, float h0b, float h1b) const {
        TerrainEdit e{};
        e.center = pt(0.5f * (x0 + x1), 0.5f * (y0 + y1), 0.5f * (h0b + h1b));
        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 0.5f * (y1 - y0), 0.5f * (h1b - h0b) };
        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        return e;
    }
    // The water level (m against the ground) in a pit at x: the surface above a point just over the floor.
    float level(const Planetoid& p, float x) const {
        const float s = Fluid::SurfaceAbove(p, pt(x, 0.0f, kFloor + 0.05f), 6.0f);
        return s > 0.0f ? kFloor + 0.05f + s : kFloor;
    }
    // Side on, a little above the ground: both pits and the line of the tunnel between them.
    void camera(ScenarioCtx& c) const { c.camera(pt(0.0f, -18.0f, 5.0f), pt(0.0f, 0.0f, -2.0f), home); }
} S;
}

static void Build(ScenarioDef& s) {
    s.tags = { "fluid", "water", "physics-reference" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 300.0;

    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("dig", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        std::vector<TerrainEdit> cut = {
            S.box(-12.0f, -8.0f, -1.0f, 1.0f, kFloor, 1.0f),       // the near pit (open to the sky)
            S.box(-8.0f, 8.0f, -1.0f, 1.0f, kFloor, kTunnelTop),   // the tunnel under 3 m of roof
            S.box(8.0f, 12.0f, -1.0f, 1.0f, kFloor, 1.0f),         // the far pit
        };
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.note("two pits joined under 3 m of roof (cells " + std::to_string(WaterField::CellSize()) + " m); expected common level " + std::to_string(kExpected) + " m");
        S.camera(c);
        return true;
    });
    s.step("pour into the near pit", 120.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.waitingFor("80 m^3 poured");
        if (c.t <= 1.0) return false;
        if (c.t >= S.nextPour) {
            S.nextPour = c.t + 0.1;
            const float x = -11.5f + 1.0f * (float)(S.k % 4), y = (S.k / 4) % 2 ? 0.5f : -0.5f;
            ++S.k;
            S.poured += Fluid::Pour(*c.body, S.pt(x, y, -0.5f), 0.25f);
        }
        if (S.poured < kPour) return false;
        S.vol0 = Fluid::GetStats().volume; // the sandbox has no other water: the simulated volume is what was poured
        return true;
    });
    s.step("level out", 150.0, [](ScenarioCtx& c) {
        S.camera(c);
        c.waitingFor("the water asleep (or 120 s)");
        if (!((c.t > 10.0 && Fluid::Asleep()) || c.t > 120.0)) return false;
        const float nearL = S.level(*c.body, -10.0f), farL = S.level(*c.body, 10.0f);
        const float tunnel = Fluid::WaterAt(*c.body, S.pt(0.0f, 0.0f, 0.5f * (kFloor + kTunnelTop)), 0.5f);
        const float lost = Fluid::PourDropped();
        const double vol1 = Fluid::GetStats().volume;
        c.note("simulated volume " + std::to_string(S.vol0) + " m^3 after the pour, " + std::to_string(vol1) + " m^3 at the end");
        c.metric("simulated volume after the pour", S.vol0, "m^3", false);
        c.metric("simulated volume at the end", vol1, "m^3", false);
        c.near("no water lost or made while it levels (m^3)", vol1, S.vol0, 0.5);
        c.note("near pit " + std::to_string(nearL) + " m, far pit " + std::to_string(farL) + " m (expected " + std::to_string(kExpected) +
               "), tunnel middle " + std::to_string(tunnel) + " m^3 in a 1 m box, poured " + std::to_string(S.poured) + " m^3 (dropped " + std::to_string(lost) + ")");
        c.metric("near pit level", nearL, "m", false);
        c.metric("far pit level", farL, "m", false);
        c.metric("level difference", std::fabs(nearL - farL), "m");
        c.gt("the tunnel is full of water (m^3 in its middle)", tunnel, 0.5);
        c.lt("the two pits stand level (m apart)", std::fabs(nearL - farL), 0.15);
        c.metric("far pit level against the box arithmetic's -1 m", farL - kExpected, "m", false);
        return true;
    });
}

SCENARIO("fluid/u-bend") { Build(s); }
// The CPU kernel without the pass: the baseline the pressure variant is compared with.
SCENARIO("fluid/u-bend-cpu") { s.env = { { "FLUID_CPU", "1" }, { "RAYLIBGAME_FLUID_PRESSURE", "0" } }; Build(s); }
// The pressure pass (#1168) on the CPU kernel, where it is built so far (fluid_step.comp follows): the same U-bend.
SCENARIO("fluid/u-bend-pressure") { s.env = { { "FLUID_CPU", "1" }, { "RAYLIBGAME_FLUID_PRESSURE", "1" } }; Build(s); }
