// visual/arena-compact (#1603): the terrain arena gives memory back by compacting (every chunk's vertices and indices copied
// packed into smaller buffers, its table entry's offsets rewritten). A wrong offset scrambles that chunk's triangles, so: the
// view settles, the screen is read, the arena is compacted now (TerrainArena_CompactForTest), and a few frames later the
// screen is read again: it must be the same picture (a small share of pixels may move with the wind, the water and clouds).
#include "scenario.h"
#include "terrainarena.h"
#include "world.h"
#include <cstdlib>
#include <string>

namespace {
struct State { Image before{}; long long capV = 0, capI = 0; int frame0 = 0; } S;
}

SCENARIO("visual/arena-compact") {
    s.tags = { "visual" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 120.0;
    s.step("settle", 90.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 3.0 && WorldIdle(); });
    s.step("read the screen, compact", 2.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        S.before = LoadImageFromScreen();
        const ArenaStats a = TerrainArena_Stats();
        S.capV = a.vertexCapacity; S.capI = a.indexCapacity;
        c.check("compacted", TerrainArena_CompactForTest());
        const ArenaStats b = TerrainArena_Stats();
        c.metric("arena vertex capacity before (M)", S.capV / 1e6, "M");
        c.metric("arena vertex capacity after (M)", b.vertexCapacity / 1e6, "M");
        c.metric("arena index capacity before (M)", S.capI / 1e6, "M");
        c.metric("arena index capacity after (M)", b.indexCapacity / 1e6, "M");
        c.check("nothing lost (vertices in use)", b.vertices == a.vertices, std::to_string(a.vertices) + " -> " + std::to_string(b.vertices));
        S.frame0 = (int)c.frame;
        return true;
    });
    s.step("the same picture", 3.0, [](ScenarioCtx& c) {
        if (c.frame < S.frame0 + 4) return false;
        Image after = LoadImageFromScreen();
        long long n = 0, diff = 0;
        if (after.width == S.before.width && after.height == S.before.height) {
            const Color* a = (const Color*)S.before.data; const Color* b = (const Color*)after.data;
            for (int i = 0; i < after.width * after.height; i += 3) {
                ++n;
                if (abs(a[i].r - b[i].r) + abs(a[i].g - b[i].g) + abs(a[i].b - b[i].b) > 24) ++diff;
            }
        }
        const double pct = n ? 100.0 * diff / n : 100.0;
        c.metric("pixels changed by the compaction (%)", pct, "%");
        c.check("the terrain draws the same after compacting (< 2% of pixels changed)", n > 0 && pct < 2.0, std::to_string(pct) + " %");
        UnloadImage(after); UnloadImage(S.before);
        return true;
    });
}
