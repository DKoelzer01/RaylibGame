#ifndef CACHES_H
#define CACHES_H

// Ruin cache rewards (#2072): what a sealed cache gives when it is opened. Pure and deterministic: the same seed, cache,
// kind, band and team progress always roll the same. Progression owns the pool (resources/data/caches.txt), alternate
// recipes and the research boost; UI/UX's open flow (ruincache.*) calls the roll on the host and applies parts and
// blueprints.

#include <cstdint>
#include <set>
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
    struct Node { std::string id, name; int tier = 0, units = 1; bool core = false; };
    std::vector<Node> frontier;             // research nodes it could start now (not done, every prerequisite done)
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
// Host, after the roll: the team-side rewards (alternate recipes into the team's list, research units). The open flow
// applies parts (inventory, overflow dropped) and blueprints (the team library) itself. (Defined in src/space: research.)
void Caches_ApplyTeam(int teamId, const std::vector<CacheRoll>& rolls);

// The pool (caches.txt), loaded on first use; a test loads its own text. Errors name the line ("" none).
bool Caches_LoadText(const std::string& text, std::string& error);
bool Caches_Load(const char* path = "resources/data/caches.txt");
struct CacheBlueprint { std::string id, name; int tier = 0; };
struct CacheAlt { std::string id, name; int tier = 0; };
const std::vector<CacheBlueprint>& Caches_Blueprints();
const std::vector<CacheAlt>& Caches_Alts();
const CacheAlt* Caches_FindAlt(const std::string& id);
// Every item a parts line names, for the data check against items.txt.
std::vector<std::string> Caches_PartItems();

// A team's alternate recipes (the team store's "alt_recipes" section: one id a line). The factory offers an alternate
// recipe to a machine only when its plot's team has it.
const std::set<std::string>& Caches_TeamAlts(int team);
bool Caches_TeamHasAlt(int team, const std::string& id);
void Caches_AddTeamAlt(int team, const std::string& id);
std::string Caches_TeamAltsText(int team);
void Caches_SetTeamAltsText(int team, const std::string& text);
void Caches_ResetTeams();

#endif
