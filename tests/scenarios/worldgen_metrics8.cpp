// worldgen/metrics-8 (#1882): worldgen/metrics-1 over 8 planets (the scenario world's own, then terrain seeds 1-7): ~8
// erosion bakes, so it runs in the #1871 full tier (nightly, develop, manual staging runs). Report-only; fails only if a
// land biome appears on none of them. See worldgen_metrics.h.
#include "worldgen_metrics.h"

SCENARIO("worldgen/metrics-8") {
    s.tags = { "worldgen" };
    worldgen_metrics::Define(s, 8);
}
