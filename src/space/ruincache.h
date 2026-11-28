#pragma once

// Ruin caches (#2072): sealed containers that ruin sets place (a `cache <kind> x y z` line in resources/structures/*.st).
// They are never broken: aimed at, E held for kCacheHoldSeconds opens one if the opener has what its kind needs
// (CacheNeed). The authority (single player, a host's own player, or the host for a joined client as a machine action)
// rolls the rewards (Caches_Roll, caches.h) and applies them; the cache's part material then carries kCacheOpened, so the
// opened state saves with the plot and reaches every peer as a part edit: a cache opens once per world (user, 2026-11-27).

#include <cstdint>
#include <string>
#include <vector>
#include "caches.h"
#include "plots.h"

class Planetoid;
class Inventory;
struct PartEdit;

constexpr uint8_t kCacheOpened = 0x10;      // a RuinCache part's material: kind in the low bits, this bit once opened
constexpr float kCacheHoldSeconds = 3.0f;   // GUESS (the user confirms): how long E is held to open one
constexpr float kCacheReach = 4.0f;         // m: the host's reach check (as machine actions)

inline CacheKind CacheKindOf(const Part& p) { const int k = p.material & 3; return k > 2 ? CacheKind::Small : (CacheKind)k; }
inline bool CacheOpened(const Part& p) { return (p.material & kCacheOpened) != 0; }
const char* CacheKindName(CacheKind k);     // "sealed cache", "sealed vault", "deep vault"

// What a kind needs (GUESS, the user confirms): a small cache a power cell (spent), a vault a crowbar (kept), a deep vault
// both. The line the prompt shows when the opener lacks it ("" when they have it).
std::string CacheMissing(const Inventory& inv, CacheKind k);
std::string CacheNeedText(CacheKind k);     // "needs a power cell" ...

// The cache's id: the world's seed and where the cache is on its body (a ruin is built in the same place on every peer).
uint64_t CacheId(const Planetoid& body, const Plot& plot, const Part& p);
// The opening team's progress, as the roll takes it: its research, and its tier as the highest tier of its done nodes.
CacheTeam CacheTeamOf(int team);

// The authority opens it: checks it is sealed and the need, spends it, rolls, applies the parts to `inv` (what doesn't fit
// goes to `overflow`) and the team rewards, and marks the part opened (the edits for the clients). False + why: nothing done.
bool Cache_Open(Planetoid& body, Plot& plot, int part, Inventory& inv, int team, std::vector<CacheRoll>& rolls,
                std::vector<std::pair<std::string, int>>& overflow, std::vector<PartEdit>* edits, std::string* why);

// The opener's screen: the rolls' lines for a few seconds, and the hold ring while E is held (progress 0..1, <0 none).
void RuinCache_ShowResult(CacheKind kind, const std::vector<std::string>& lines);
void RuinCache_SetHold(float progress);
void RuinCache_DrawHUD();
const std::vector<std::string>& RuinCache_LastResult(); // tests: the lines last shown
