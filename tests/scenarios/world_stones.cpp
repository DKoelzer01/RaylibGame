// world/stones-per-biome (#893): loose stones (tree-system variant 9: the tier-0 quest's first 4 stone) lie in every biome
// but on the mountaintops. On the planet at seed 42, for each land biome: the nearest deep patch of it (the --bench-biome
// search, PlanetField::nearestBiomeDir), and the scatter as generated within 200 m of it (Trees_CountScatter: no edits or removals).
// Each biome must have at least 4 stones there (alpine is only reported: its high ground is bare by design). The trees and
// bushes are reported too. Before #893 the swamp (clay and mud underfoot) and the desert (sand) had none.
#include "scenario.h"
#include "trees.h"
#include "voxel.h"
#include "world.h"
#include <cstdio>

SCENARIO("world/stones-per-biome") {
    s.tags = { "world", "foliage", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;
    s.record = false; // (one frame of counting: nothing to watch)
    s.step("count the scatter in each biome", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const unsigned seed = c.body->trees ? c.field->foliageSeed : 0u;
        for (int b = (int)Biome::Jungle; b <= (int)Biome::Alpine; ++b) {
            const char* name = BiomeName((Biome)b);
            const Vector3 dir = c.field->nearestBiomeDir(c.field->findSpawnDir(c.body->seed), (Biome)b); // (not gSpawnBiome: it caches one per field)
            const float e = c.field->elevation(dir);
            if (c.field->biomeAt(dir, e) != (Biome)b) { c.note(std::string("no patch of ") + name + " found (skipped)"); continue; }
            int trees = 0, bushes = 0, stones = 0;
            Trees_CountScatter(*c.field, seed, dir, 200.0f, trees, bushes, stones);
            c.metric(std::string(name) + " trees within 200 m", trees, "", false);
            c.metric(std::string(name) + " bushes within 200 m", bushes, "", false);
            if ((Biome)b == Biome::Alpine) c.metric(std::string(name) + " stones within 200 m", stones, "", false);
            else c.gt(std::string(name) + " stones within 200 m", stones, 3.5);
        }
        return true;
    });
}
