// worldgen/metrics-1 (#1882): the terrain's numbers for one planet, the scenario world's own (seed 42's Terra): land %, biome
// shares, slope histogram, flat buildable share, bake hash, the spawn. Report-only; fails only if a land biome vanished.
// The PR-tier variant (#1871 runs it on terrain and erosion changes); worldgen/metrics-8 runs 8 planets in the full tier.
#include "worldgen_metrics.h"

SCENARIO("worldgen/metrics-1") {
    s.tags = { "worldgen" };
    worldgen_metrics::Define(s, 1);
}
