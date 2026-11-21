// perf/blueprint-place (#1167, the user: "Game freezes for a moment when placing blueprints"): a 12-part and a 60-part
// blueprint of stone-brick cubes placed half sunk in the sandbox ground (their lower layer clears the terrain under it,
// as a foundation does), creative. Reported: the placement call's wall time and the worst frame over the next second
// (the terrain re-mesh lands there). Metrics only (timing is the local machine's: compare A/B, RAYLIBGAME_PROFILE=1
// for the zones "blueprint place" and "blueprint foundations").
#include "scenario.h"
#include "blueprint.h"
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
struct Run { int parts = 0; double placeMs = 0.0, worstFrameMs = 0.0; int placed = 0; };
Run runs[2];
int stage = 0;
double settleUntil = 0.0;
std::chrono::steady_clock::time_point lastFrame;
}

SCENARIO("perf/blueprint-place") {
    s.tags = { "perf", "build" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("place a 12-part and a 60-part blueprint half in the ground", 40.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        const auto now = std::chrono::steady_clock::now();
        const double frameMs = std::chrono::duration<double, std::milli>(now - lastFrame).count();
        lastFrame = now;
        if (stage == 1 || stage == 3) { // the second after a placement: its worst frame
            Run& r = runs[stage / 2];
            if (c.T > settleUntil - 0.95) r.worstFrameMs = std::max(r.worstFrameMs, frameMs);
            if (c.T < settleUntil) return false;
            ++stage;
            return false;
        }
        if (stage >= 4) {
            for (const Run& r : runs) {
                c.eq(TextFormat("%d-part blueprint: every part placed", r.parts), r.placed, r.parts);
                c.metric(TextFormat("%d-part blueprint: the placement call", r.parts), r.placeMs, "ms");
                c.metric(TextFormat("%d-part blueprint: the worst frame in the next second", r.parts), r.worstFrameMs, "ms");
            }
            return true;
        }
        terrainTool.creative = true; // (free: the cost isn't what is measured)
        const int k = stage / 2;
        const Blueprint bp = k == 0 ? Slab(3, 2, 2) : Slab(6, 5, 2); // 12 and 60 cubes
        // A fresh plot 8 m (then 20 m) ahead, its floor 0.75 m under the ground: the lower layer is in the ground.
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, k == 0 ? 8.0f : 20.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) - 0.75f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(ground, 0.0f, true);
        c.camera(Vector3Add(ground, Vector3Add(Vector3Scale(at, 6.0f), Vector3Scale(east, -8.0f))), ground, at); // the build in view
        std::string missing;
        const auto t0 = std::chrono::steady_clock::now();
        runs[k].placed = Blueprint_Place(*c.body, *plot, bp, { 0, 0, 0 }, true, missing);
        runs[k].placeMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        runs[k].parts = (int)bp.parts.size();
        settleUntil = c.T + 1.0;
        lastFrame = std::chrono::steady_clock::now();
        ++stage;
        return false;
    });
}
