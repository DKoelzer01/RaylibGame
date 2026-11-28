#pragma once

// What opens a ruin cache (#2072): its opener, one per cache (user, 2026-11-27: "power, a specific part, or a tool; later
// also a key dropped by an enemy or hidden in the ruin"). The catalogue is data, resources/data/cache_openers.txt:
//   opener <n> <id> power <EU> <cap EU/t> <min V>   charged from a real grid or a battery bank (not the suit battery)
//   opener <n> <id> tool <item>                     the item is needed and kept (a crowbar)
//   opener <n> <id> part <item>                     the item is spent (a key)
//   default <small|vault|deep_vault> <id>           the opener of a cache whose part names none
// n (1-7) is what a RuinCache part's material keeps in bits 5-7 (0: its kind's default), so a saved or sent part keeps its
// opener whatever the file's order; n never changes once used. A structure set names one with `cache <kind> x y z opener <id>`.
//
// The power side is Progression's (the factory charges the cache's lock from a cable or a battery bank up to the cap, at the
// voltage or more): it reads eu, capEUt and minVolts from here and reports the charge through CacheOpeners_SetPower.

#include <cstdint>
#include <functional>
#include <string>
#include <vector>
#include "caches.h"
#include "plots.h"

class Planetoid;
class Plot;

enum class OpenerKind : uint8_t { Power, Tool, Part };

struct CacheOpener {
    int index = 0;                 // 1-7, the material bits
    std::string id;
    OpenerKind kind = OpenerKind::Power;
    long long eu = 0;              // Power: EU the lock takes to open
    int capEUt = 0;                // Power: the most it draws per factory tick
    int minVolts = 0;              // Power: a network below this can't charge it
    std::string item;              // Tool / Part: the item id
};

constexpr int kCacheOpenerShift = 5;      // a RuinCache part's material: bits 5-7 the opener index (0: the kind's default)
constexpr uint8_t kCacheOpenerMask = 0xE0;
inline uint8_t CacheOpenerBits(int index) { return (uint8_t)((index & 7) << kCacheOpenerShift); }
constexpr uint8_t kCacheOpened = 0x10;    // a RuinCache part's material: kind in the low bits, this bit once opened
inline bool CacheOpened(const Part& p) { return (p.material & kCacheOpened) != 0; }

// The catalogue (resources/data/cache_openers.txt, read on first use). Parse replaces it from text (unit tests); false and
// a reason on a line it can't read (the lines before it stay).
bool CacheOpeners_Parse(const std::string& text, std::string* error = nullptr);
void CacheOpeners_Reset(); // forget it: the file is read again on next use (after a test's own catalogue)
const CacheOpener* CacheOpener_Find(const std::string& id);
const CacheOpener* CacheOpener_ByIndex(int index);
// The opener of a cache part: its material's index, else its kind's default; null only with an empty catalogue.
const CacheOpener* CacheOpenerOf(const Part& p);

// The prompt's need ("needs 2,000 EU at 8 V or more (a battery bank or a cable)", "needs a Crowbar (kept)").
std::string CacheOpener_NeedText(const CacheOpener& o);

// The power side (Progression, the factory): charged = how full the cache's lock is (0..1); spend empties it when the cache
// opens. Until it is set, a power opener reads 0 and its cache stays sealed.
struct CachePowerHooks {
    std::function<float(const Planetoid&, const Plot&, int part)> charged;
    std::function<void(Planetoid&, Plot&, int part)> spend;
};
void CacheOpeners_SetPower(CachePowerHooks hooks);
float CacheOpener_Charge(const Planetoid& body, const Plot& plot, int part); // 0 without the hooks
void CacheOpener_Spend(Planetoid& body, Plot& plot, int part);
