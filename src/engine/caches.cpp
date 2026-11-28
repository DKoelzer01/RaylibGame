#include "caches.h"

// (#2072) A stub until Progression's caches.txt and roll land: one common parts roll, so the open flow can be built and
// tested with the real signature.
std::vector<CacheRoll> Caches_Roll(uint64_t, uint64_t, CacheKind, int, const CacheTeam&) {
    CacheRoll r;
    r.kind = CacheReward::Parts;
    r.rarity = CacheRarity::Common;
    r.items = { { "stone", 4 } };
    r.text = "4 stone";
    return { r };
}

void Caches_ApplyTeam(int, const std::vector<CacheRoll>&) {} // the stub rolls no team-side rewards
