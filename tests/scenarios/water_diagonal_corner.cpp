// water/diagonal-corner (#1773, playtest: "Water will slightly visibly clip through diagonal corners, but wont flow through"):
// in the sandbox, a floating 5 x 5 m floor of stone-brick cubes (plot cells, 0.5 m) with a 2 x 2 m pen on it, walled 1 m high
// on all four sides, whose two inner walls meet only diagonally at one corner: the cells beside the corner are wall, the
// corner cell itself is floor outside the pen. Water poured into the pen (0.5 m^3/s for 4 s) fills it about 0.45 m deep.
// The fluid site's cells (0.25 m) are turned against the plot here (its frame comes from the first water near it), so pen cells
// reach past the corner: the drawn water poked out there as a wedge (capture 2026-11-24 diag-corner-1773-shift).
// Checked while it settles: the water MESH (WaterVertexProbe: every vertex the water builds put near the corner) has no vertex
// over the floor outside the corner, it has some in the pen, and the simulation put no water outside the corner either (so
// a drawn wedge is the mesher's, not a leak).
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
#include <vector>

namespace {
constexpr int kFloor = 10; // floor cells across (5 m)
constexpr int kPen = 5;    // the pen's walls at cells 0 and kPen (inside: 1..kPen-1); the open corner cell is (kPen, kPen)
struct State {
    Plot* plot = nullptr;
    Int3 base{};
    double nextPour = 0.0;
    float poured = 0.0f;
    int checks = 0, simWet = 0, penVerts = 0, pastVerts = 0, verts = 0;
    float worst = 0.0f; // the farthest a vertex lay past the corner (m, the nearer of its two distances from the walls' faces)
    Vector3 at(float x, float y, float z) const { return plot->toBody({ x, y, z }); } // plot metres -> body-local
    float floorTop() const { return (base.y + 1) * kCellSize; }
} S;
float CornerX() { return (S.base.x + kPen) * kCellSize; } // the open corner cell's near edges (plot metres)
float CornerZ() { return (S.base.z + kPen) * kCellSize; }
void CornerCamera(ScenarioCtx& c) { // close outside the open corner, looking into it (the wedge is a few cm)
    c.world.camera.position = c.body->toWorld(S.at(CornerX() + 0.8f, S.floorTop() + 0.7f, CornerZ() + 0.8f));
    c.world.camera.target = c.body->toWorld(S.at(CornerX(), S.floorTop() + 0.2f, CornerZ()));
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
        S.base.x += 4; S.base.y += 2; S.base.z -= kFloor / 2 + 3; // ahead of the player, 1 m up (clear of the ground)
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
        // The probe: a box (body-local) around the corner, 0.6 m each way.
        Vector3 lo{ 1e30f, 1e30f, 1e30f }, hi{ -1e30f, -1e30f, -1e30f };
        for (int k = 0; k < 8; ++k) {
            const Vector3 q = S.at(CornerX() + (k & 1 ? 0.6f : -0.6f), S.floorTop() + (k & 2 ? 1.0f : -0.2f), CornerZ() + (k & 4 ? 0.6f : -0.6f));
            lo = Vector3Min(lo, q); hi = Vector3Max(hi, q);
        }
        WaterVertexProbe_Set(c.body->id, lo, hi);
        gScenarioCamera = true;
        CornerCamera(c);
        return true;
    });
    s.step("fill the pen and look past the corner (12 s)", 40.0, [](ScenarioCtx& c) {
        CornerCamera(c);
        const float mid = (S.base.x + kPen * 0.5f + 0.5f) * kCellSize, midz = (S.base.z + kPen * 0.5f + 0.5f) * kCellSize;
        if (c.t >= S.nextPour && c.t < 4.0) { S.nextPour = c.t + 0.1; S.poured += Fluid::Pour(*c.body, S.at(mid, S.floorTop() + 0.3f, midz), 0.05f); }
        if (c.t > 12.0) return true;
        if (c.t < 3.0 || (int)(c.t * 60.0) % 30 != 0) return false; // twice a second
        ++S.checks;
        std::vector<Vector3> got;
        WaterVertexProbe_Take(got);
        for (const Vector3& v : got) {
            const Vector3 q = S.plot->toPlot(v);
            const float h = q.y - S.floorTop(), dx = q.x - CornerX(), dz = q.z - CornerZ();
            if (h < 0.01f || h > 1.0f) continue; // in the floor (hidden) or over the walls
            ++S.verts;
            if (dx < -0.05f && dz < -0.05f) ++S.penVerts;
            if (dx > 0.005f && dz > 0.005f) { ++S.pastVerts; S.worst = std::max(S.worst, std::min(dx, dz)); }
        }
        // The simulation's cells over the floor outside the corner (5 cm grid, 0.3 m square, up to 0.9 m over the floor).
        const WaterField::SnapshotPtr snap = WaterField::Current();
        for (float dx = 0.025f; dx < 0.3f; dx += 0.05f)
            for (float dz = 0.025f; dz < 0.3f; dz += 0.05f)
                for (float h = 0.05f; h < 0.9f; h += 0.1f)
                    if (WaterField::CellFill(*snap, c.body->id, S.at(CornerX() + dx, S.floorTop() + h, CornerZ() + dz), 0) > 0) { ++S.simWet; break; }
        return false;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        gScenarioCamera = false;
        WaterVertexProbe_Set(0, { 1, 1, 1 }, { 0, 0, 0 }); // (an empty box: records nothing)
        c.note(TextFormat("%.1f m^3 poured; %d checks; water mesh vertices by the corner %d: %d in the pen, %d past the corner (farthest "
                          "%.3f m); %d samples outside the corner with simulated water",
                          S.poured, S.checks, S.verts, S.penVerts, S.pastVerts, S.worst, S.simWet));
        c.metric("water vertices past the corner", S.pastVerts, "vertices");
        c.metric("simulated water outside the corner", S.simWet, "samples");
        c.gt("the pen's water was drawn by the corner", S.penVerts, 0);
        c.eq("the simulation keeps the water in the pen", S.simWet, 0);
        c.eq("no water drawn past the diagonal corner", S.pastVerts, 0);
        return true;
    });
}
