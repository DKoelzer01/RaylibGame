// perf/blueprint-place (#1167, the user: "Game freezes for a moment when placing blueprints"): a 12-part and a 60-part
// blueprint of stone-brick cubes placed half sunk in the sandbox ground (their lower layer clears the terrain under it,
// as a foundation does), creative. Reported: the placement call's wall time and the worst frame over the next second
// (the terrain re-mesh lands there), against a baseline second at the same spot, the world idle (#1372). Metrics only (timing is the local machine's: compare A/B, RAYLIBGAME_PROFILE=1
// for the zones "blueprint place" and "blueprint foundations").
#include "scenario.h"
#include "blueprint.h"
#include "buildtool.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <chrono>

namespace {
Blueprint Slab(int w, int d, int h) { // w x d cubes, h layers
    Blueprint bp; bp.name = TextFormat("slab %dx%dx%d", w, d, h);
    const int mat = std::max(0, BuildMaterialForItem("stone_brick"));
    for (int y = 0; y < h; ++y) for (int x = 0; x < w; ++x) for (int z = 0; z < d; ++z)
        bp.parts.push_back({ PartType::Cube, { x, y, z }, 0, (uint8_t)mat });
    return bp;
}
// drain: (#1183) frames the queued foundations took; baseline: (#1372) the worst frame over the same second at the same spot,
// the world idle, nothing placed (the first version placed while the world still streamed in: its "worst frame" was mostly
// that streaming)
struct Run { int parts = 0; double placeMs = 0.0, worstFrameMs = 0.0, baselineMs = 0.0; int placed = 0, drainFrames = 0, pendingAfter = 0; };
Run runs[2];
int stage = 0; // per blueprint k: 4k aim (the camera, then the world idle), 4k+1 baseline, 4k+2 place, 4k+3 the second after
int idleFrames = 0;
double phaseUntil = 0.0;
std::chrono::steady_clock::time_point lastFrame;
Vector3 groundOf[2], atOf[2], eastOf[2];
}

SCENARIO("perf/blueprint-place") {
    s.tags = { "perf", "build" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("place a 12-part and a 60-part blueprint half in the ground", 60.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        const auto now = std::chrono::steady_clock::now();
        const double frameMs = std::chrono::duration<double, std::milli>(now - lastFrame).count();
        lastFrame = now;
        if (stage >= 8) {
            for (const Run& r : runs) {
                c.eq(TextFormat("%d-part blueprint: every part placed", r.parts), r.placed, r.parts);
                c.metric(TextFormat("%d-part blueprint: the placement call", r.parts), r.placeMs, "ms");
                c.metric(TextFormat("%d-part blueprint: the worst frame in the next second", r.parts), r.worstFrameMs, "ms");
                c.metric(TextFormat("%d-part blueprint: the worst frame in a second without it (baseline)", r.parts), r.baselineMs, "ms");
                c.metric(TextFormat("%d-part blueprint: frames until its foundations were all laid", r.parts), r.drainFrames, "frames");
                c.eq(TextFormat("%d-part blueprint: its foundations all laid within the second", r.parts), r.pendingAfter, 0);
            }
            return true;
        }
        const int k = stage / 4, phase = stage % 4;
        Run& r = runs[k];
        if (phase == 0) { // the camera on the spot; wait until the world has streamed what it shows
            terrainTool.creative = true; // (free: the cost isn't what is measured)
            const Vector3 dir = c.spawnDirLocal();
            eastOf[k] = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
            // A fresh plot 8 m (then 20 m) ahead, its floor 0.75 m under the ground: the lower layer is in the ground.
            atOf[k] = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(eastOf[k], k == 0 ? 8.0f : 20.0f)));
            groundOf[k] = Vector3Scale(atOf[k], c.field->radius + c.field->elevation(atOf[k]) - 0.75f);
            c.camera(Vector3Add(groundOf[k], Vector3Add(Vector3Scale(atOf[k], 6.0f), Vector3Scale(eastOf[k], -8.0f))), groundOf[k], atOf[k]); // the build in view
            idleFrames = WorldIdle() ? idleFrames + 1 : 0;
            if (idleFrames < 10) return false; // (the step's own timeout ends a world that never settles)
            idleFrames = 0; phaseUntil = c.T + 1.0; ++stage;
            return false;
        }
        if (phase == 1 || phase == 3) { // a second at the spot: its worst frame (3: after the placement)
            if (c.T > phaseUntil - 0.95) (phase == 1 ? r.baselineMs : r.worstFrameMs) = std::max(phase == 1 ? r.baselineMs : r.worstFrameMs, frameMs);
            if (phase == 3 && BuildTool_PendingFoundations() > 0) ++r.drainFrames; // (#1183) still being laid, a share a frame
            if (c.T < phaseUntil) return false;
            if (phase == 3) r.pendingAfter = BuildTool_PendingFoundations();
            ++stage;
            return false;
        }
        // phase 2: place it
        const Blueprint bp = k == 0 ? Slab(3, 2, 2) : Slab(6, 5, 2); // 12 and 60 cubes
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(groundOf[k], 0.0f, true);
        std::string missing;
        const auto t0 = std::chrono::steady_clock::now();
        r.placed = Blueprint_Place(*c.body, *plot, bp, { 0, 0, 0 }, true, missing);
        r.placeMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        r.parts = (int)bp.parts.size();
        phaseUntil = c.T + 1.0;
        lastFrame = std::chrono::steady_clock::now();
        ++stage;
        return false;
    });
}
