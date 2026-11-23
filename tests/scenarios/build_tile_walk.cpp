// build/tile-walk (#1732, the user: the player bounces off placed plot tiles instead of stepping onto them or sliding along):
// in the sandbox, a floor of one-cell stone-brick blocks (3 cells wide, 36 long), a 1-tile (0.5 m) step at its far end (a
// second layer, 2 cells deep, across the width) and a 2-cell wall along one long side. The player walks (W, scripted through
// the real input path) straight along the floor and up the step, then diagonally into the wall along it. Checks: no catch on
// the floor's seams, the step is climbed, no bounce (never moving backwards), and the wall is slid along.
// (The floor is longer than the 6 cells of the issue's sketch: 4 s of walking is ~16-20 m.)
#include "scenario.h"
#include "player.h"
#include "plots.h"
#include "playerinput.h"
#include "world.h"
#include "buildtool.h"
#include "terraintool.h"
#include <raymath.h>
#include <memory>
#include <algorithm>
#include <cmath>

namespace {
constexpr uint64_t kW = 1ull << (int)Act::Forward;
constexpr int kLen = 36, kWid = 3, kStepDeep = 2;   // floor cells along x and across z; the step's depth in cells
Plot* plot = nullptr;
Int3 base{};                                         // the floor's low corner (cell); y = the floor's row
Vector3 upL{}, xL{}, zL{}, diagL{};                  // body-local: up, the walk (+x), +z (toward the wall), the diagonal
double floorH = 0.0;                                 // the floor top's height (dot with up)
double steady = 0.0, flatMin = 1e9, allMin = 1e9, tReach = -1.0, stepRise = 0.0;
double wallMin = 1e9, wallAllMin = 1e9, wallX0 = 0.0, wallT0 = -1.0;
int phase = 0;

Vector3 Feet(ScenarioCtx& c) { return c.body->toLocalD(player.position).toFloat(); }
Vector3 AtM(float x, float y, float z) { return plot->toBody({ x, y, z }); }  // plot metres -> body-local
double Height(ScenarioCtx& c) { return Vector3DotProduct(Feet(c), upL); }
double Along(ScenarioCtx& c) { return plot->toPlot(Feet(c)).x; }              // metres along the floor (plot x)
// Speed relative to the ground along a body-local direction, from the feet's body-local motion since the last call (doubles:
// with the frames (Plan 11) player.velocity is already ground-relative near a body, and subtracting surfaceVelocity took the
// spin off twice). 0 on the first call or when no time passed.
DVec3 prevFeet{}; double prevT = -1.0, lastSpeed = 0.0;
double SpeedAlong(ScenarioCtx& c, Vector3 dirL) {
    const DVec3 now = c.body->toLocalD(player.position);
    if (prevT >= 0.0 && c.T > prevT) {
        const DVec3 d = now - prevFeet;
        lastSpeed = (d.x * dirL.x + d.y * dirL.y + d.z * dirL.z) / (c.T - prevT);
    }
    prevFeet = now; prevT = c.T;
    return lastSpeed;
}
void Hold(bool on) { if (on) gFakeHeld |= kW; else gFakeHeld &= ~kW; }
void Face(ScenarioCtx& c, Vector3 dirL) { player.forward = c.body->dirToWorld(dirL); }
void FollowCam(ScenarioCtx& c, Vector3 dirL) { // 4 m behind the player along the walk, 2 m up, looking ahead
    const Vector3 feet = Feet(c);
    const Vector3 eye = Vector3Add(Vector3Subtract(feet, Vector3Scale(dirL, 4.0f)), Vector3Scale(upL, 2.0f));
    c.world.camera.position = c.body->toWorld(eye);
    c.world.camera.target = c.body->toWorld(Vector3Add(Vector3Add(feet, Vector3Scale(dirL, 2.0f)), Vector3Scale(upL, 0.4f)));
    c.world.camera.up = c.body->dirToWorld(upL);
}
// Put the player on the floor at plot metres (x, z), standing 0.3 m above it (it drops).
void PlaceOnFloor(ScenarioCtx& c, float x, float z) {
    const float top = (base.y + 1) * kCellSize;
    player.flyMode = false;
    player.placeAt(c.world, *c.body, c.body->toWorld(AtM(x, top + 0.3f, z)));
    prevT = -1.0; lastSpeed = 0.0; // a teleport is not a speed
}
}

void BuildTileWalk(ScenarioDef& s, ScenarioWorld world);
SCENARIO("build/tile-walk") { BuildTileWalk(s, ScenarioWorld::Sandbox); }
// On the planet, 60 km from its centre, where a plot's positions pass through float transforms (#1732's playtest).
SCENARIO("build/tile-walk-planet") { BuildTileWalk(s, ScenarioWorld::Planet); }

void BuildTileWalk(ScenarioDef& s, ScenarioWorld world) {
    s.tags = { "build", "physics" };
    s.world = world;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("build the tile floor, the step and the wall", 20.0, [](ScenarioCtx& c) {
        player.flyMode = false; // (a scenario starts the player flying over the start)
        if (c.T < 1.0 || !player.grounded || !WorldIdle()) return false;
        phase = 0; flatMin = allMin = wallMin = wallAllMin = 1e9; steady = 0.0; tReach = -1.0; wallT0 = -1.0;
        const Vector3 feet = Feet(c);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.y += 2; // (clear of the ground)
        const int mat = std::max(0, BuildMaterialForItem("stone_brick"));
        int floor = 0, step = 0, wall = 0;
        for (int x = 0; x < kLen; ++x)
            for (int z = 0; z < kWid; ++z) {
                floor += plot->place(Part{ PartType::Cube, { base.x + x, base.y, base.z + z }, 0, (uint8_t)mat }) >= 0;
                if (x >= kLen - kStepDeep) step += plot->place(Part{ PartType::Cube, { base.x + x, base.y + 1, base.z + z }, 0, (uint8_t)mat }) >= 0;
            }
        for (int x = 0; x < kLen; ++x)
            for (int y = 0; y < 3; ++y) // (the floor's row and the two cells above it)
                wall += plot->place(Part{ PartType::Cube, { base.x + x, base.y + y, base.z + kWid }, 0, (uint8_t)mat }) >= 0;
        c.eq("floor cells placed", floor, kLen * kWid);
        c.eq("step cells placed", step, kStepDeep * kWid);
        c.eq("wall cells placed", wall, kLen * 3);
        upL = Vector3Normalize(Vector3Subtract(AtM(0, 1, 0), AtM(0, 0, 0)));
        xL = Vector3Normalize(Vector3Subtract(AtM(1, 0, 0), AtM(0, 0, 0)));
        zL = Vector3Normalize(Vector3Subtract(AtM(0, 0, 1), AtM(0, 0, 0)));
        diagL = Vector3Normalize(Vector3Add(xL, zL));
        floorH = Vector3DotProduct(AtM(0, (base.y + 1) * kCellSize, 0), upL);
        gScenarioInput = true;
        gScenarioCamera = true;
        PlaceOnFloor(c, (base.x + 2) * kCellSize, (base.z + kWid * 0.5f) * kCellSize);
        return true;
    });
    s.step("settle on the floor", 8.0, [](ScenarioCtx& c) {
        Face(c, xL); FollowCam(c, xL);
        return c.t > 1.0 && player.grounded;
    });
    s.step("walk the floor and up the step (W, ~4 s)", 20.0, [](ScenarioCtx& c) {
        Face(c, xL); FollowCam(c, xL);
        Hold(true);
        const double v = SpeedAlong(c, xL), along = Along(c);
        const double stepFront = (base.x + kLen - kStepDeep) * kCellSize;
        if (tReach < 0.0) { // on the flat floor
            if (c.t >= 1.0 && along < stepFront - 0.6) { steady = std::max(steady, v); flatMin = std::min(flatMin, v); }
            if (along >= stepFront - 0.3) tReach = c.t;
        }
        if (c.t > 0.2) allMin = std::min(allMin, v);
        if (tReach >= 0.0) {
            stepRise = std::max(stepRise, Height(c) - floorH);
            if (c.t - tReach < 2.0) return false;
        } else return false;
        Hold(false);
        c.gt("steady walking speed on the floor (m/s)", steady, 1.0);
        c.metric("walk: min speed over the seams (m/s)", flatMin, "m/s", false);
        c.metric("walk: steady speed (m/s)", steady, "m/s", false);
        c.check("no catch on seams: speed stays >= 70% of steady across the floor", flatMin >= 0.7 * steady,
                TextFormat("min %.2f m/s of steady %.2f m/s", flatMin, steady));
        c.metric("step: rise within 2 s of reaching it (m)", stepRise, "m", false);
        c.gt("the step is climbed: feet rise (m, want >= 0.4)", stepRise, 0.4);
        c.metric("walk: min speed along the walk (m/s)", allMin, "m/s", false);
        c.gt("no bounce: speed along the walk never negative (m/s)", allMin, -0.05);
        c.note(TextFormat("walk: steady %.2f m/s, min over seams %.2f, min overall %.2f, step rise %.2f m", steady, flatMin, allMin, stepRise));
        return true;
    });
    s.step("walk diagonally into the wall, along it (2 s)", 20.0, [](ScenarioCtx& c) {
        if (phase == 0) { // back on the floor, near the wall, early in it
            Hold(false);
            PlaceOnFloor(c, (base.x + 2) * kCellSize, (base.z + kWid) * kCellSize - 0.35f);
            phase = 1; wallT0 = -1.0;
            return false;
        }
        Face(c, diagL); FollowCam(c, diagL);
        if (phase == 1) { if (c.t > 1.0 && player.grounded) { phase = 2; wallT0 = c.t; wallX0 = Along(c); } return false; }
        Hold(true);
        const double v = SpeedAlong(c, xL), t = c.t - wallT0;
        if (t > 0.4) wallMin = std::min(wallMin, v);
        wallAllMin = std::min(wallAllMin, v);
        if (t < 2.0) return false;
        Hold(false);
        const double need = 0.5 * std::max(steady, 1.0);
        c.metric("wall: min speed along the wall (m/s)", wallMin, "m/s", false);
        c.metric("wall: distance slid along it (m)", Along(c) - wallX0, "m", false);
        c.check("slides along the wall: speed along it >= 50% of walking speed", wallMin >= need,
                TextFormat("min %.2f m/s, want >= %.2f", wallMin, need));
        c.gt("no reversal at the wall: speed along it never negative (m/s)", wallAllMin, -0.05);
        c.note(TextFormat("wall: min along %.2f m/s, min overall %.2f, slid %.2f m", wallMin, wallAllMin, Along(c) - wallX0));
        gFakeHeld = 0; gScenarioInput = false; gScenarioCamera = false;
        return true;
    });
}
