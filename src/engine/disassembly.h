#pragma once

// Disassembly (#1427 ladder L7c): what an assembler set to Disassemble gives back for an item, from
// resources/data/disassembly.txt and the forward recipes (hand recipes and assembler recipes, by item name).
// - A made item: its first recipe's direct inputs, per craft of that recipe (outCount items). Inputs no forward recipe
//   makes are raw (amount x yield, fractions carried); the others are components (a chance per unit).
// - A composed tool or weapon (composed.h): the raw materials of its head, handle and binding (part lines).
// - Tiers: voltage tiers (0 LV, -1 one under ...), from the ladder tier (recipe_unlocks / research, via the hook) less
//   the file's lv_ladder_tier, else the ingredients' highest.
// Built once on the main thread before the factory runs; read-only afterwards (the sim thread reads it).
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

struct DisassemblyRecipe {
    std::string output;
    int outCount = 1;
    std::vector<std::pair<std::string, int>> in;
    int ticks = 0, eu = 0; // 0: a hand recipe (the file's hand line)
};
struct DisassemblyPart { std::string item; int count = 0; bool raw = true; };
struct DisassemblyEntry {
    int tier = 0;          // voltage tier
    int outCount = 1;      // items one craft of the forward recipe makes: each taken gives back parts / outCount
    int ticks = 0, eu = 0; // the forward recipe's (before overclocking)
    std::vector<DisassemblyPart> parts;
};

bool Disassembly_LoadText(const std::string& text);                       // disassembly.txt's lines (tests)
bool Disassembly_Load(const char* path = "resources/data/disassembly.txt");
// The index from the forward recipes (the first per output wins: hand recipes first, then the assembler's).
// ladderTier: an item's ladder tier from its unlock or research node, -1 if it has none.
void Disassembly_Build(const std::vector<DisassemblyRecipe>& recipes, const std::function<int(const std::string&)>& ladderTier);
const DisassemblyEntry* Disassembly_Find(const std::string& item); // nullptr: not taken apart (raw, excluded, unknown)
std::vector<std::string> Disassembly_Items();                       // every item the index takes apart
// The game's ladder tier of an item (recipe_unlocks "tier N", else its research node's tier; -1 none). Set before the
// factory loads its recipes (they build the index).
inline std::function<int(const std::string& item)> gDisassemblyLadderTier;
int Disassembly_ItemTier(const std::string& item);                  // voltage tier (raw items: their ladder tier, else the lowest)
// A composed item (class id and record) as an entry: its parts' raw materials at 100%, its tier the parts' highest, a
// hand recipe's time. False if it is not one or gives nothing.
bool Disassembly_Composed(const std::string& cls, const std::string& data, DisassemblyEntry& out);
bool Disassembly_Excluded(const std::string& item); // the file's excludes and the game's (foods: gDisassemblyExclude)
inline std::function<bool(const std::string& item)> gDisassemblyExclude;
// The parts of one item taken at `yieldPct`: raw parts as hundredths (count x yield / outCount, rounded down), components
// as rolls: `rolls` units each kept with chance yield / (100 x outCount). rng: a 32-bit state stepped per roll.
struct DisassemblyGive { std::string item; int centi = 0; };
void Disassembly_Give(const DisassemblyEntry& e, int yieldPct, uint32_t& rng, std::vector<DisassemblyGive>& out);
uint32_t Disassembly_Next(uint32_t& rng); // xorshift32 (a zero state is reseeded)
