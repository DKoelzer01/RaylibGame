// water/diagonal-corner (#1773, playtest: "Water will slightly visibly clip through diagonal corners, but wont flow through"):
// in the sandbox, a floating 5 x 5 m floor of stone-brick cubes (plot cells, 0.5 m) with a 2 x 2 m pen on it, walled 1 m high
// on all four sides, whose two inner walls meet only diagonally at one corner: the cells beside the corner are wall, the
// corner cell itself is floor outside the pen. Water poured into the pen (0.5 m^3/s for 4 s) fills it about 0.5 m deep.
// Checked while it settles, over the floor outside that corner (a 0.75 m square at the corner, on a 5 cm grid): the DRAWN
// surface (WaterField::ChunkSampler, what the water mesh is built from) shows no water under the pen's level, and the
// simulation put none there (so a drawn sliver is the mesher's, not a leak).
#include "scenario.h"
#include "fluid.h"
#include "waterfield.h"
#include "plots.h"
#include "player.h"
#include "buildtool.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <memory>

namespace {
constexpr int kFloor = 10; // floor cells across (5 m)
constexpr int kPen = 5;    // the pen's walls at cells 0 and kPen (inside: 1..kPen-1); the open corner cell is (kPen, kPen)
struct State {
    Plot* plot = nullptr;
    Int3 base{};
    double nextPour = 0.0;
    float poured = 0.0f;
    int checks = 0, samples = 0, drawnWet = 0, simWet = 0, penWet = 0;
    float worst = -1.0f; // the highest drawn water outside the corner, over the floor (m)
    Vector3 at(float x, float y, float z) const { return plot->toBody({ x, y, z }); } // plot metres -> body-local
    float floorTop() const { return (base.y + 1) * kCellSize; }
} S;
float CornerX() { return (S.base.x + kPen) * kCellSize; } // the open corner cell's near edges (plot metres)
float CornerZ() { return (S.base.z + kPen) * kCellSize; }
void CornerCamera(ScenarioCtx& c) { // outside the open corner, over it, looking at it and into the pen
    c.world.camera.position = c.body->toWorld(S.at(CornerX() + 2.0f, S.floorTop() + 2.2f, CornerZ() + 2.0f));
    c.world.camera.target = c.body->toWorld(S.at(CornerX() - 0.3f, S.floorTop() + 0.3f, CornerZ() - 0.3f));
    c.world.camera.up = c.body->dirToWorld(Vector3Normalize(Vector3Subtract(S.at(0, 1, 0), S.at(0, 0, 0))));
}
}

SCENARIO("water/diagonal-corner") {
    s.tags = { "fluid", "water", "build" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("build the floor and the pen", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        S.plot = c.body->plots->plotFor(feet, 0.0f, true);
        S.base = S.plot->cellAt(S.plot->toPlot(feet));
        S.base.x += 4; S.base.y += 2; S.base.z -= kFloor / 2; // ahead of the player, 1 m up (clear of the ground)
        const uint8_t mat = (uint8_t)std::max(0, BuildMaterialForItem("stone_brick"));
        int placed = 0, want = 0;
        auto put = [&](int x, int y, int z) { ++want; placed += S.plot->place(Part{ PartType::Cube, { S.base.x + x, S.base.y + y, S.base.z + z }, 0, mat }) >= 0; };
        for (int x = 0; x < kFloor; ++x)
            for (int z = 0; z < kFloor; ++z) put(x, 0, z);
        for (int y = 1; y <= 2; ++y)
            for (int i = 0; i <= kPen; ++i) {
                put(i, y, 0); if (i > 0) put(0, y, i);                  // the outer walls
                if (i > 0 && i < kPen) { put(kPen, y, i); put(i, y, kPen); } // the inner walls, not at (kPen, kPen)
            }
        c.eq("cubes placed", placed, want);
        gScenarioCamera = true;
        CornerCamera(c);
        return true;
    });
    s.step("fill the pen and look past the corner (12 s)", 40.0, [](ScenarioCtx& c) {
        CornerCamera(c);
        const float mid = (S.base.x + kPen * 0.5f + 0.5f) * kCellSize, midz = (S.base.z + kPen * 0.5f + 0.5f) * kCellSize;
        if (c.t >= S.nextPour && c.t < 4.0) { S.nextPour = c.t + 0.1; S.poured += Fluid::Pour(*c.body, S.at(mid, S.floorTop() + 0.3f, midz), 0.05f); }
        if (c.t < 3.0 || c.t > 12.0) return c.t > 12.0;
        if ((int)(c.t * 60.0) % 30 != 0) return false; // twice a second
        ++S.checks;
        const WaterField::SnapshotPtr snap = WaterField::Current();
        const Vector3 lo = S.at(CornerX() - 1.5f, S.floorTop() - 0.5f, CornerZ() - 1.5f);
        const Vector3 hi = S.at(CornerX() + 1.5f, S.floorTop() + 1.5f, CornerZ() + 1.5f);
        WaterField::ChunkSampler cs(snap.get(), c.body->id, Vector3Min(lo, hi), Vector3Max(lo, hi));
        auto drawnTop = [&](float x, float z) { // the highest wet point of the vertical from 0.9 m over the floor down to 2 cm over it
            for (float h = 0.9f; h >= 0.02f; h -= 0.01f) {
                float wd;
                if (cs.density(S.at(x, S.floorTop() + h, z), 0, wd) && wd >= ISO_LEVEL) return h;
            }
            return -1.0f;
        };
        // The pen holds water (its corner cell, inside).
        if (drawnTop(CornerX() - 0.25f, CornerZ() - 0.25f) > 0.0f) ++S.penWet;
        for (float dx = 0.02f; dx < 0.75f; dx += 0.05f)
            for (float dz = 0.02f; dz < 0.75f; dz += 0.05f) {
                ++S.samples;
                const float x = CornerX() + dx, z = CornerZ() + dz;
                const float top = drawnTop(x, z);
                if (top > 0.0f) { ++S.drawnWet; S.worst = std::max(S.worst, top); }
                for (float h = 0.05f; h < 0.9f; h += 0.1f)
                    if (WaterField::CellFill(*snap, c.body->id, S.at(x, S.floorTop() + h, z), 0) > 0) { ++S.simWet; break; }
            }
        return false;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        gScenarioCamera = false;
        c.note(TextFormat("%.1f m^3 poured; pen wet in %d of %d checks; outside the corner %d samples: %d drawn wet (highest %.3f m over the "
                          "floor), %d with simulated water",
                          S.poured, S.penWet, S.checks, S.samples, S.drawnWet, S.worst, S.simWet));
        c.metric("drawn water outside the corner", S.drawnWet, "samples");
        c.metric("simulated water outside the corner", S.simWet, "samples");
        c.gt("the pen held water", S.penWet, S.checks / 2);
        c.eq("the simulation keeps the water in the pen", S.simWet, 0);
        c.eq("no water drawn past the diagonal corner", S.drawnWet, 0);
        return true;
    });
}
