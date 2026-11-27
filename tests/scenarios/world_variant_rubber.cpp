// world/variant-rubber (#2012): trees felled in the old-growth, mangrove and acacia variants give rubber a third of the time
// (biomes.txt `rubber 33`, the user's ruling 7 on #1964); a plain mixed forest gives none, the jungle always. On the planet at
// seed 42, for each variant: its nearest deep patch (the --bench-biome search), then 300 medium trees at points within 150 m
// that read as that variant, each with its own id, through TerrainTool_PlantDrops (what a chop, and the host, give).
#include "scenario.h"
#include "biomevariants.h"
#include "terraintool.h"
#include "trees.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

SCENARIO("world/variant-rubber") {
    s.tags = { "world", "foliage", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;
    s.record = false; // (one frame of counting: nothing to watch)
    s.step("rubber from trees in each variant", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const BiomeVariantTable& t = BiomeVariants();
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        struct Want { const char* key; double lo, hi; };
        const Want wants[] = { { "forest.old_growth", 0.25, 0.42 }, { "swamp.mangrove", 0.25, 0.42 }, { "savanna.acacia", 0.25, 0.42 },
                               { "forest.mixed", 0.0, 0.0 } };
        for (const Want& w : wants) {
            const int v = t.find(w.key);
            if (!c.check(std::string(w.key) + " is in biomes.txt", v >= 0)) continue;
            const Vector3 at = c.field->nearestVariantDir(home, v);
            if (c.field->variantAt(at, c.field->elevation(at)) != v) { c.note(std::string("no patch of ") + w.key + " found (skipped)"); continue; }
            const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(at, { 0, 1, 0 })), e2 = Vector3CrossProduct(e1, at);
            int trees = 0, rubber = 0;
            for (int k = 0; k < 4000 && trees < 300; ++k) {
                const float a = k * 2.39996f, r = 150.0f * sqrtf((k + 0.5f) / 4000.0f); // a sunflower spiral over 150 m
                const Vector3 d = Vector3Normalize(Vector3Add(at, Vector3Scale(Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a))), r / c.field->radius)));
                const float el = c.field->elevation(d);
                if (c.field->variantAt(d, el) != v) continue;
                Plant p{};
                p.id = 0x5EED0000ull + (uint64_t)k * 7919u; p.base = Vector3Scale(d, c.field->radius + el); p.up = d;
                p.height = 5.2f; p.radius = 0.3f; p.scale = 1.0f; p.variant = 0; // a medium oak
                ++trees;
                for (const auto& [item, n] : TerrainTool_PlantDrops(*c.body, p).items) if (item == "rubber") rubber += n;
            }
            const double share = trees ? (double)rubber / trees : 0.0;
            c.note(TextFormat("%s: %d trees, %d rubber (%.0f%%)", w.key, trees, rubber, share * 100.0));
            c.gt(std::string(w.key) + ": trees sampled", trees, 99.0);
            if (w.hi <= 0.0) c.eq(std::string(w.key) + ": no rubber", rubber, 0);
            else { c.gt(std::string(w.key) + ": rubber share", share, w.lo); c.lt(std::string(w.key) + ": rubber share", share, w.hi); }
        }
        // The jungle type: every tree, whatever its variant.
        const Vector3 j = c.field->nearestBiomeDir(home, Biome::Jungle);
        if (c.field->biomeAt(j, c.field->elevation(j)) == Biome::Jungle) {
            Plant p{};
            p.id = 12345; p.base = Vector3Scale(j, c.field->radius + c.field->elevation(j)); p.up = j;
            p.height = 5.2f; p.radius = 0.3f; p.scale = 1.0f; p.variant = 0;
            int n = 0;
            for (const auto& [item, k] : TerrainTool_PlantDrops(*c.body, p).items) if (item == "rubber") n += k;
            c.eq("a jungle tree gives one rubber", n, 1);
        } else c.note("no jungle patch found (skipped)");
        return true;
    });
}
