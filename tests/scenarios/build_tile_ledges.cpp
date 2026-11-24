// build/tile-ledges (#1732, the user: on a flat pad of plot tiles, walking INTO something standing on it - another tile, a
// machine, a chest - bounces the player in an odd direction instead of just blocking): lanes side by side, each a floor of
// stone-brick cubes (4 cells wide, 16 long) with an obstruction across or on its far end - a slab (0.25 m), a cube (0.5 m), a
// cube and a slab (0.75 m), 1 m blocks, the blocks jumped onto, a chest, a furnace - and a wide lane whose 1 m wall is walked
// into at ~31 degrees. The player walks (W through the real input path) into each. The velocity is the feet's motion in the
// plot's axes (ground-relative). Checks: against a wall the velocity stays in its plane - never pushed back (>= -0.05 m/s),
// no sideways kick head-on (|side| <= 0.15 m/s), no upward kick (<= 0.3 m/s) - the angled walk keeps sliding along the wall;
// a slab is stepped onto without a launch (the feet never pass its top); the jump lands on the blocks. A 0.5 m cube is a wall.
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
constexpr int kLen = 16, kWid = 4, kLedgeAt = 10, kPitch = 6, kObliqueWid = 14; // cells: lane length, width, the ledge's first cell, lane spacing
enum Ledge { Slab, Cube, CubeSlab, Block1m, Block1mJump, Chest, Furnace, Oblique, kLanes };
const char* kName[kLanes] = { "slab 0.25 m", "cube 0.5 m", "cube+slab 0.75 m", "block 1 m", "block 1 m, jumped", "chest", "furnace", "1 m wall at 31 deg" };
Plot* plot = nullptr;
Int3 base{};
Vector3 upL{}, xL{}, zL{}, walkL{};
double floorH = 0.0;
int lane = 0, phase = 0;
double t0 = 0.0, tHit = -1.0, backMin = 1e9, sideMax = 0.0, upMax = -1e9, slideMin = 1e9, hMax = -1e9;
bool jumped = false;
DVec3 prevFeet{}; double prevT = -1.0; Vector3 vel{}; // the feet's motion in plot axes (x along the lane, y up, z across), m/s

Vector3 Feet(ScenarioCtx& c) { return c.body->toLocalD(player.position).toFloat(); }
Vector3 AtM(float x, float y, float z) { return plot->toBody({ x, y, z }); }
double Height(ScenarioCtx& c) { return Vector3DotProduct(Feet(c), upL) - floorH; } // above the floor's top (m)
double Along(ScenarioCtx& c) { return plot->toPlot(Feet(c)).x / kCellSize - base.x; } // cells along the lane
void Measure(ScenarioCtx& c) { // doubles: 60 km out a float position steps ~4 mm, a 0.25 m/s error at 60 Hz
    const DVec3 now = c.body->toLocalD(player.position);
    if (prevT >= 0.0 && c.T > prevT) {
        const DVec3 d = now - prevFeet; const double dt = c.T - prevT;
        auto dot = [&](Vector3 a) { return (float)((d.x * a.x + d.y * a.y + d.z * a.z) / dt); };
        vel = { dot(xL), dot(upL), dot(zL) };
    }
    prevFeet = now; prevT = c.T;
}
void FollowCam(ScenarioCtx& c) { // behind and beside the player, so the obstruction and the feet are both in view
    const Vector3 feet = Feet(c);
    const Vector3 side = Vector3CrossProduct(upL, xL);
    const Vector3 eye = Vector3Add(Vector3Add(Vector3Subtract(feet, Vector3Scale(xL, 3.0f)), Vector3Scale(upL, 1.6f)), Vector3Scale(side, 2.0f));
    c.world.camera.position = c.body->toWorld(eye);
    c.world.camera.target = c.body->toWorld(Vector3Add(Vector3Add(feet, Vector3Scale(xL, 1.5f)), Vector3Scale(upL, 0.6f)));
    c.world.camera.up = c.body->dirToWorld(upL);
}
int LaneZ(int l) { return base.z + l * kPitch; }
}

void BuildTileLedges(ScenarioDef& s, ScenarioWorld world);
SCENARIO("build/tile-ledges") { BuildTileLedges(s, ScenarioWorld::Sandbox); }
// On the planet, 60 km out on a spinning body (#1732's playtest): the frame carry and grounded flicker matter there.
SCENARIO("build/tile-ledges-planet") { BuildTileLedges(s, ScenarioWorld::Planet); }

void BuildTileLedges(ScenarioDef& s, ScenarioWorld world) {
    s.tags = { "world", "build", "physics" };
    s.world = world;
    s.args = { "--day" };
    s.timeout = 150.0;
    s.step("build the lanes", 20.0, [](ScenarioCtx& c) {
        player.flyMode = false;
        if (c.T < 1.0 || !player.grounded || !WorldIdle()) return false;
        const Vector3 feet = Feet(c);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.y += 3; // (clear of the ground)
        base.z -= kPitch * kLanes / 2;
        const uint8_t mat = (uint8_t)std::max(0, BuildMaterialForItem("stone_brick"));
        int placed = 0, want = 0;
        auto put = [&](PartType t, int x, int y, int z) { ++want; placed += plot->place(Part{ t, { x, y, z }, 0, mat }) >= 0; };
        for (int l = 0; l < kLanes; ++l) {
            const int z0 = LaneZ(l), wid = l == Oblique ? kObliqueWid : kWid;
            for (int x = 0; x < kLen; ++x)
                for (int z = 0; z < wid; ++z) put(PartType::Cube, base.x + x, base.y, z0 + z);
            if (l == Chest || l == Furnace) { put(l == Chest ? PartType::Chest : PartType::Furnace, base.x + kLedgeAt, base.y + 1, z0 + 1); continue; }
            for (int x = kLedgeAt; x < kLen; ++x)
                for (int z = 0; z < wid; ++z) {
                    if (l == Slab) put(PartType::Slab, base.x + x, base.y + 1, z0 + z);
                    if (l == Cube || l == CubeSlab || l == Oblique) put(PartType::Cube, base.x + x, base.y + 1, z0 + z);
                    if (l == CubeSlab) put(PartType::Slab, base.x + x, base.y + 2, z0 + z);
                    if (l == Oblique) put(PartType::Cube, base.x + x, base.y + 2, z0 + z);
                    if ((l == Block1m || l == Block1mJump) && (x - kLedgeAt) % 2 == 0 && z % 2 == 0) put(PartType::Block, base.x + x, base.y + 1, z0 + z);
                }
        }
        c.eq("parts placed", placed, want);
        upL = Vector3Normalize(Vector3Subtract(AtM(0, 1, 0), AtM(0, 0, 0)));
        xL = Vector3Normalize(Vector3Subtract(AtM(1, 0, 0), AtM(0, 0, 0)));
        zL = Vector3Normalize(Vector3Subtract(AtM(0, 0, 1), AtM(0, 0, 0)));
        floorH = Vector3DotProduct(AtM(0, (base.y + 1) * kCellSize, 0), upL);
        gScenarioInput = true;
        gScenarioCamera = true;
        lane = 0; phase = 0;
        return true;
    });
    s.step("walk into each obstruction", 140.0, [](ScenarioCtx& c) {
        if (lane >= kLanes) { gFakeHeld = 0; gScenarioInput = false; gScenarioCamera = false; return true; }
        walkL = lane == Oblique ? Vector3Normalize(Vector3Add(xL, Vector3Scale(zL, 0.6f))) : xL; // (atan 0.6 = 31 degrees)
        player.forward = c.body->dirToWorld(walkL);
        if (phase == 0) { // on the lane's floor, 6 cells before the obstruction, standing 0.3 m above it
            gFakeHeld = 0;
            player.flyMode = false;
            const float z = (LaneZ(lane) + (lane == Oblique ? 1.0f : kWid * 0.5f)) * kCellSize;
            player.placeAt(c.world, *c.body, c.body->toWorld(AtM((base.x + kLedgeAt - 6) * kCellSize, (base.y + 1) * kCellSize + 0.3f, z)));
            prevT = -1.0; vel = {}; tHit = -1.0; backMin = 1e9; sideMax = 0.0; upMax = -1e9; slideMin = 1e9; hMax = -1e9; jumped = false;
            phase = 1; t0 = c.t;
            return false;
        }
        FollowCam(c);
        if (phase == 1) { if (c.t - t0 > 1.0 && player.grounded) { phase = 2; t0 = c.t; } return false; }
        const double t = c.t - t0, h = Height(c), along = Along(c);
        Measure(c);
        gFakeHeld = kW;
        if (lane == Block1mJump && !jumped && along >= kLedgeAt - 1.6) { gFakeHeld |= kJump; jumped = true; } // ~0.8 m short of it
        hMax = std::max(hMax, h);
        if (lane == Cube && along > kLedgeAt - 2.0 && t < 2.0)
            TraceLog(LOG_INFO, "ledgetrace t %.3f along %.4f h %.4f vx %.2f vy %.2f vz %.2f grounded %d pv_up %.2f gn_up %.3f", t, along, h, vel.x, vel.y, vel.z,
                     (int)player.grounded, Vector3DotProduct(player.velocity, c.body->dirToWorld(upL)), Vector3DotProduct(player.groundNormal, c.body->dirToWorld(upL)));
        // At the obstruction: from 0.1 s after the feet come within a radius (0.45 m) of its face, while they stay below its top.
        if (tHit < 0.0 && along >= kLedgeAt - 1.0) tHit = t;
        if (tHit >= 0.0 && t > tHit + 0.1 && h < 0.2) {
            backMin = std::min(backMin, (double)vel.x);
            sideMax = std::max(sideMax, (double)fabsf(vel.z));
            upMax = std::max(upMax, (double)vel.y);
            slideMin = std::min(slideMin, (double)vel.z);
        }
        if (lane == Slab) upMax = std::max(upMax, (double)vel.y); // (stepped onto: the rise rate is the launch)
        const bool slidOff = lane == Oblique && tHit >= 0.0 && t > tHit + 1.2; // (~3.5 m along its 7 m wall)
        if (t < 3.0 && along < kLedgeAt + 3.0 && !slidOff) return false;
        gFakeHeld = 0;
        const char* n = kName[lane];
        c.metric(TextFormat("%s: min speed along the walk at it (m/s)", n), backMin, "m/s", false);
        c.metric(TextFormat("%s: max upward speed at it (m/s)", n), upMax, "m/s", true);
        if (lane == Slab) {
            c.check("slab: stepped onto", h > 0.2 && along > kLedgeAt + 1.0, TextFormat("feet %.2f m over the floor at cell %.1f", h, along));
            c.check("slab: no launch (the feet never pass its top by > 5 cm)", hMax <= 0.25 + 0.05, TextFormat("highest %.2f m, rising at up to %.2f m/s", hMax, upMax));
        } else if (lane == Block1mJump) {
            c.check("block 1 m, jumped: lands on it", h > 0.9 && along > kLedgeAt + 0.5, TextFormat("feet %.2f m over the floor at cell %.1f", h, along));
        } else {
            c.gt(TextFormat("%s: never pushed back (m/s)", n), backMin, -0.05);
            c.check(TextFormat("%s: no upward kick (<= 0.3 m/s)", n), upMax <= 0.3, TextFormat("%.2f m/s", upMax));
            if (lane == Oblique) {
                c.metric("1 m wall at 31 deg: min slide along it (m/s)", slideMin, "m/s", false);
                c.check("1 m wall at 31 deg: keeps sliding along it (>= 1 m/s)", slideMin >= 1.0, TextFormat("%.2f m/s", slideMin));
            } else c.check(TextFormat("%s: no sideways kick (<= 0.15 m/s)", n), sideMax <= 0.15, TextFormat("%.2f m/s", sideMax));
        }
        c.note(TextFormat("%s: at it: min along %.2f m/s, max |side| %.2f, min side %.2f, max up %.2f; highest %.2f m, end %.2f m at cell %.1f",
                          n, backMin, sideMax, slideMin, upMax, hMax, h, along));
        ++lane; phase = 0;
        return false;
    });
}
