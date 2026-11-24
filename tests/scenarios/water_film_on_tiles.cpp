// water/film-on-tiles (#1772, playtest: "water pouring over the top of placed tiles doesn't show as a thin film on the tile, it
// is slightly clipped into the surface of the tile"): in the sandbox, a floating 4 x 4 m pad of stone-brick cubes (plot cells,
// 0.5 m) 1 m over the ground, and water poured onto its middle (0.5 m^3/s) runs across it and off its edges. Over the pad, on
// a 0.25 m grid, the DRAWN surface (WaterField::ChunkSampler, what the water mesh is built from) is found by scanning each
// point's vertical down from 0.3 m over the tiles' top. Checked while it pours: where the water stands over the pad, it is
// drawn above the tiles' top face (a film), not inside the tiles. The fluid cells (0.25 m) don't line up with the plot's top,
// so a cell straddling it stays open; its water stood at the cell's floor, under the top.
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
constexpr int kPad = 8;   // cells across (4 m)
struct State {
    Plot* plot = nullptr;
    Int3 base{};
    double nextPour = 0.0;
    float poured = 0.0f;
    int samples = 0, wet = 0, inside = 0;
    float worst = 1e9f; // the lowest drawn surface over the top among the wet points (m)
    Vector3 at(float x, float y, float z) const { return plot->toBody({ x, y, z }); } // plot metres -> body-local
    float topY() const { return (base.y + 1) * kCellSize; }
} S;
void Camera(ScenarioCtx& c) { // over a corner of the pad, looking at its middle
    const float mid = (S.base.x + kPad * 0.5f) * kCellSize, midz = (S.base.z + kPad * 0.5f) * kCellSize;
    c.world.camera.position = c.body->toWorld(S.at(mid - 3.5f, S.topY() + 2.5f, midz - 3.5f));
    c.world.camera.target = c.body->toWorld(S.at(mid, S.topY(), midz));
    c.world.camera.up = c.body->dirToWorld(Vector3Normalize(Vector3Subtract(S.at(0, 1, 0), S.at(0, 0, 0))));
}
}

SCENARIO("water/film-on-tiles") {
    s.tags = { "fluid", "water", "build" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("build the pad", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        S.plot = c.body->plots->plotFor(feet, 0.0f, true);
        S.base = S.plot->cellAt(S.plot->toPlot(feet));
        S.base.x += 4; S.base.y += 2; S.base.z -= kPad / 2; // ahead of the player, 1 m up (clear of the ground)
        const uint8_t mat = (uint8_t)std::max(0, BuildMaterialForItem("stone_brick"));
        int placed = 0;
        for (int x = 0; x < kPad; ++x)
            for (int z = 0; z < kPad; ++z) placed += S.plot->place(Part{ PartType::Cube, { S.base.x + x, S.base.y, S.base.z + z }, 0, mat }) >= 0;
        c.eq("pad cubes placed", placed, kPad * kPad);
        gScenarioCamera = true;
        Camera(c);
        return true;
    });
    s.step("pour onto the pad and measure the drawn film (8 s)", 40.0, [](ScenarioCtx& c) {
        Camera(c);
        const float mid = (S.base.x + kPad * 0.5f) * kCellSize, midz = (S.base.z + kPad * 0.5f) * kCellSize;
        if (c.t >= S.nextPour && c.t < 8.0) { S.nextPour = c.t + 0.1; S.poured += Fluid::Pour(*c.body, S.at(mid, S.topY() + 0.3f, midz), 0.05f); }
        if (c.t < 2.0 || c.t > 8.0) return c.t > 8.0;
        if ((int)(c.t * 60.0) % 30 != 0) return false; // twice a second
        const WaterField::SnapshotPtr snap = WaterField::Current();
        const Vector3 lo = S.at(S.base.x * kCellSize - 1.0f, S.topY() - 1.0f, S.base.z * kCellSize - 1.0f);
        const Vector3 hi = S.at((S.base.x + kPad) * kCellSize + 1.0f, S.topY() + 1.0f, (S.base.z + kPad) * kCellSize + 1.0f);
        WaterField::ChunkSampler cs(snap.get(), c.body->id, Vector3Min(lo, hi), Vector3Max(lo, hi));
        // Points over the pad's inner part (0.25 m from its edges: the edge itself is where water falls off).
        for (float x = S.base.x * kCellSize + 0.25f; x <= (S.base.x + kPad) * kCellSize - 0.25f; x += 0.25f)
            for (float z = S.base.z * kCellSize + 0.25f; z <= (S.base.z + kPad) * kCellSize - 0.25f; z += 0.25f) {
                ++S.samples;
                // The drawn surface: the highest point of the vertical (0.3 m over the top down to 0.25 m under it) where the
                // water's density is at or over ISO.
                float surf = -1e9f;
                for (float h = 0.3f; h >= -0.25f; h -= 0.005f) {
                    float wd;
                    if (cs.density(S.at(x, S.topY() + h, z), 0, wd) && wd >= ISO_LEVEL) { surf = h; break; }
                }
                if (surf < -1e8f) continue; // dry here
                ++S.wet;
                S.worst = std::min(S.worst, surf);
                if (surf < 0.005f) ++S.inside;
            }
        return false;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        gScenarioCamera = false;
        c.note(TextFormat("%.1f m^3 poured; %d samples over the pad, %d wet, %d drawn inside the tiles (< 5 mm over their top); lowest wet surface %.3f m",
                          S.poured, S.samples, S.wet, S.inside, S.worst));
        c.metric("wet samples drawn inside the tiles", S.inside, "samples");
        c.metric("lowest drawn surface over the tiles' top (m)", S.worst, "m", false);
        c.gt("water lay on the pad (wet samples)", S.wet, 20);
        c.eq("no wet sample drawn inside the tiles (surface >= 5 mm over their top)", S.inside, 0);
        return true;
    });
}
