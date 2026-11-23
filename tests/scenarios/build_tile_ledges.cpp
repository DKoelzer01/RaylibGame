// build/tile-ledges (#1732, the user: "walking on placed tiles is fine until you walk into another tile, which you will bounce
// off of"): build/tile-walk climbs a one-cube (0.5 m) step and slides along a wall without a bounce, so this walks head-on into
// every other ledge a plot floor can have. Lanes side by side in the sandbox, each a floor of stone-brick cubes (4 cells wide,
// 16 long) with a ledge across its far end: a slab (0.25 m), a cube (0.5 m), a cube and a slab (0.75 m), 1 m blocks, and the
// blocks again jumped onto. The player walks (W through the real input path) into each. Checks per lane: never pushed back
// (speed along the walk >= -0.05 m/s), no pop (the feet never rise more than 0.15 m above where they end up), and the low
// ledges are climbed, the 1 m one stops the walk and the jump lands on it.
#include "scenario.h"
#include "player.h"
#include "plots.h"
#include "playerinput.h"
#include "world.h"
#include "buildtool.h"
#include <raymath.h>
#include <memory>
#include <algorithm>
#include <cmath>

namespace {
constexpr uint64_t kW = 1ull << (int)Act::Forward, kJump = 1ull << (int)Act::Jump;
constexpr int kLen = 16, kWid = 4, kLedgeAt = 10, kPitch = 6; // floor cells along x, across z; the ledge's first cell; lane spacing (z)
enum Ledge { Slab, Cube, CubeSlab, Block1m, Block1mJump, kLanes };
const char* kName[kLanes] = { "slab 0.25 m", "cube 0.5 m", "cube+slab 0.75 m", "block 1 m", "block 1 m, jumped" };
const float kRise[kLanes] = { 0.25f, 0.5f, 0.75f, 1.0f, 1.0f };   // the ledge's height over the floor (m)
const bool kClimb[kLanes] = { true, true, false, false, true };    // walked (or jumped) onto, or a wall that stops the walk
Plot* plot = nullptr;
Int3 base{};
Vector3 upL{}, xL{};
double floorH = 0.0;
int lane = 0, phase = 0;
double t0 = 0.0, vMin = 1e9, hMax = -1e9, vUpMax = -1e9, prevH = 0.0, prevHt = -1.0;
bool jumped = false;

Vector3 Feet(ScenarioCtx& c) { return c.body->toLocalD(player.position).toFloat(); }
Vector3 AtM(float x, float y, float z) { return plot->toBody({ x, y, z }); }
double Height(ScenarioCtx& c) { return Vector3DotProduct(Feet(c), upL) - floorH; } // above the floor's top (m)
double Along(ScenarioCtx& c) { return plot->toPlot(Feet(c)).x / kCellSize - base.x; } // cells along the lane
DVec3 prevFeet{}; double prevT = -1.0, lastSpeed = 0.0;
double SpeedAlong(ScenarioCtx& c) { // ground-relative, from the feet's body-local motion (as build/tile-walk)
    const DVec3 now = c.body->toLocalD(player.position);
    if (prevT >= 0.0 && c.T > prevT) { const DVec3 d = now - prevFeet; lastSpeed = (d.x * xL.x + d.y * xL.y + d.z * xL.z) / (c.T - prevT); }
    prevFeet = now; prevT = c.T;
    return lastSpeed;
}
void FollowCam(ScenarioCtx& c) { // beside and behind the player, so the ledge and the feet are both in view
    const Vector3 feet = Feet(c);
    const Vector3 side = Vector3CrossProduct(upL, xL);
    const Vector3 eye = Vector3Add(Vector3Add(Vector3Subtract(feet, Vector3Scale(xL, 3.0f)), Vector3Scale(upL, 1.6f)), Vector3Scale(side, 2.0f));
    c.world.camera.position = c.body->toWorld(eye);
    c.world.camera.target = c.body->toWorld(Vector3Add(Vector3Add(feet, Vector3Scale(xL, 1.5f)), Vector3Scale(upL, 0.6f)));
    c.world.camera.up = c.body->dirToWorld(upL);
}
int LaneZ(int l) { return base.z + l * kPitch; }
}

SCENARIO("build/tile-ledges") {
    s.tags = { "build", "physics" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 150.0;
    s.step("build the lanes", 20.0, [](ScenarioCtx& c) {
        player.flyMode = false;
        if (c.T < 1.0 || !player.grounded || !WorldIdle()) return false;
        const Vector3 feet = Feet(c);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.y += 2;
        base.z -= kPitch * kLanes / 2;
        const uint8_t mat = (uint8_t)std::max(0, BuildMaterialForItem("stone_brick"));
        int placed = 0, want = 0;
        auto put = [&](PartType t, int x, int y, int z) { ++want; placed += plot->place(Part{ t, { x, y, z }, 0, mat }) >= 0; };
        for (int l = 0; l < kLanes; ++l) {
            const int z0 = LaneZ(l);
            for (int x = 0; x < kLen; ++x)
                for (int z = 0; z < kWid; ++z) put(PartType::Cube, base.x + x, base.y, z0 + z);
            for (int x = kLedgeAt; x < kLen; ++x)
                for (int z = 0; z < kWid; ++z) {
                    if (l == Slab) put(PartType::Slab, base.x + x, base.y + 1, z0 + z);
                    if (l == Cube || l == CubeSlab) put(PartType::Cube, base.x + x, base.y + 1, z0 + z);
                    if (l == CubeSlab) put(PartType::Slab, base.x + x, base.y + 2, z0 + z);
                    if ((l == Block1m || l == Block1mJump) && (x - kLedgeAt) % 2 == 0 && z % 2 == 0) put(PartType::Block, base.x + x, base.y + 1, z0 + z);
                }
        }
        c.eq("parts placed", placed, want);
        upL = Vector3Normalize(Vector3Subtract(AtM(0, 1, 0), AtM(0, 0, 0)));
        xL = Vector3Normalize(Vector3Subtract(AtM(1, 0, 0), AtM(0, 0, 0)));
        floorH = Vector3DotProduct(AtM(0, (base.y + 1) * kCellSize, 0), upL);
        gScenarioInput = true;
        gScenarioCamera = true;
        lane = 0; phase = 0;
        return true;
    });
    s.step("walk into each ledge", 120.0, [](ScenarioCtx& c) {
        if (lane >= kLanes) { gFakeHeld = 0; gScenarioInput = false; gScenarioCamera = false; return true; }
        player.forward = c.body->dirToWorld(xL);
        if (phase == 0) { // on the lane's floor, 6 cells before the ledge, standing 0.3 m above it
            gFakeHeld = 0;
            player.flyMode = false;
            player.placeAt(c.world, *c.body, c.body->toWorld(AtM((base.x + kLedgeAt - 6) * kCellSize, (base.y + 1) * kCellSize + 0.3f, (LaneZ(lane) + kWid * 0.5f) * kCellSize)));
            prevT = -1.0; lastSpeed = 0.0; vMin = 1e9; hMax = -1e9; vUpMax = -1e9; prevHt = -1.0; jumped = false;
            phase = 1; t0 = c.t;
            return false;
        }
        FollowCam(c);
        if (phase == 1) { if (c.t - t0 > 1.0 && player.grounded) { phase = 2; t0 = c.t; } return false; }
        const double t = c.t - t0, v = SpeedAlong(c), h = Height(c), along = Along(c);
        gFakeHeld = kW;
        if (lane == Block1mJump && !jumped && along >= kLedgeAt - 1.6) { gFakeHeld |= kJump; jumped = true; } // ~0.8 m short of it
        if (t > 0.5) vMin = std::min(vMin, v);
        hMax = std::max(hMax, h);
        if (prevHt >= 0.0 && c.T > prevHt && lane != Block1mJump) vUpMax = std::max(vUpMax, (h - prevH) / (c.T - prevHt));
        prevH = h; prevHt = c.T;
        if (t < 3.5) return false;
        gFakeHeld = 0;
        const char* n = kName[lane];
        const double hEnd = h, rise = kRise[lane];
        c.metric(TextFormat("%s: min speed along the walk (m/s)", n), vMin, "m/s", false);
        c.metric(TextFormat("%s: highest feet over the floor (m)", n), hMax, "m", true);
        c.gt(TextFormat("%s: never pushed back (m/s)", n), vMin, -0.05);
        c.check(TextFormat("%s: no pop above where the feet end up", n), hMax <= hEnd + 0.15,
                TextFormat("highest %.2f m, end %.2f m", hMax, hEnd));
        if (kClimb[lane]) c.check(TextFormat("%s: climbed", n), hEnd >= rise - 0.08 && along > kLedgeAt + 0.5, TextFormat("feet %.2f m over the floor at cell %.1f", hEnd, along));
        else c.check(TextFormat("%s: stops the walk at its face", n), hEnd < 0.1 && along < kLedgeAt, TextFormat("feet %.2f m over the floor at cell %.1f", hEnd, along));
        c.note(TextFormat("%s: min along %.2f m/s, highest %.2f m, end %.2f m at cell %.1f, max rise rate %.2f m/s", n, vMin, hMax, hEnd, along, vUpMax));
        ++lane; phase = 0;
        return false;
    });
}
