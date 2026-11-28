#ifndef CACHES_H
#define CACHES_H

// Ruin cache rewards (#2072): what a sealed cache gives when it is opened. Pure and deterministic: the same seed, cache,
// kind, band and team progress always roll the same. Progression owns the pool (resources/data/caches.txt), alternate
// recipes and the research boost; UI/UX's open flow (ruincache.*) calls the roll on the host and applies parts and
// blueprints. Until caches.txt lands, the roll is a stub: one common parts roll.

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

enum class CacheKind : uint8_t { Small, Vault, DeepVault };
enum class CacheReward : uint8_t { Parts, Blueprint, Alt, Boost };
enum class CacheRarity : uint8_t { Common, Uncommon, Rare };

struct CacheTeam {                          // the opening team's progress, as the host holds it
    int tier = 0;                           // its highest open ladder tier (#1271)
    std::vector<std::string> researched;    // research node ids done (any order)
    std::vector<std::string> blueprints;    // team-library blueprint ids it has
    std::vector<std::string> alts;          // alternate recipe ids it has
};

struct CacheRoll {
    CacheReward kind = CacheReward::Parts;
    CacheRarity rarity = CacheRarity::Common;
    std::vector<std::pair<std::string, int>> items; // Parts: item, count
    std::string id;                         // Blueprint: preset id; Alt: alternate recipe id; Boost: research node id
    int units = 0;                          // Boost: units added (25% of the node, at least 5)
    std::string text;                       // the result screen's line, e.g. "Alternate recipe: Packed coke (coke oven)"
};

// band: 0 surface or shallow, 1 middle (250-500 m), 2 deep (500 m+), from the cache cell's depth.
std::vector<CacheRoll> Caches_Roll(uint64_t worldSeed, uint64_t cacheId, CacheKind kind, int band, const CacheTeam& team);
// Host, after the roll: the team-side rewards (alternate recipes into the teamstore, research units). The open flow applies
// parts (inventory, overflow dropped) and blueprints (the team library) itself.
void Caches_ApplyTeam(int teamId, const std::vector<CacheRoll>& rolls);

#endif
