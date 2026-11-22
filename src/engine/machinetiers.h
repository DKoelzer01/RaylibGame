#pragma once

// Tiered machines (#1427, ladder L7a): resources/data/machine_tiers.txt.
// - Voltage tiers (LV 32 V, MV 128 V, HV 512 V, ...), numbered 0, 1, 2, ... in the file's order.
// - Families: a machine (its machines.txt name) with the tier it starts at and the tiers it has variants for. Each variant is
//   its own item (the base machine's item at its first tier; "<base>_<tier>" or a name from the file above it) placing the
//   same part type, with the variant's step above the first tier in Part::material (machines use no build material).
// - A recipe has a tier (machine_recipes.txt's sixth column; empty: the machine's first tier). A machine runs the recipes
//   of its tier and below; one N tiers under it overclocks N steps (min(N, cap)): time / speed^N, EU a tick x eu^N.
// - The disassembly yields (L7c) by the item's tier against the assembler's.
#include <string>
#include <vector>

#include "plots.h"

struct MachineOverclock { int speed = 2, eu = 3, cap = 0; }; // cap 0: no cap

bool MachineTiers_Load(); // resources/data/machine_tiers.txt (once; true when read)
int MachineTiers_Count();                          // voltage tiers
const char* MachineTiers_Name(int tier);           // "LV", "MV", ... ("" out of range)
int MachineTiers_Volts(int tier);                  // 32, 128, ... (0 out of range)
int MachineTiers_ByName(const std::string& name);  // "lv" / "LV" -> 0; -1 if unknown
int MachineTiers_ForVolts(int volts);              // the tier whose voltage that is (the highest at or under it); -1 under LV
const MachineOverclock& MachineTiers_Overclock();

// Families: is this part type tiered, its first tier, its highest variant step (0: only its base).
bool MachineTiers_Tiered(PartType t);
int MachineTiers_First(PartType t);                // -1 if not tiered
int MachineTiers_MaxStep(PartType t);              // 0 if not tiered
// A part's voltage tier: a tiered machine's first tier + its step (Part::material); else its machines.txt voltage's tier.
int MachineTiers_PartTier(PartType t, int step, int baseVolts);
// The item of a variant (step 0: the base machine's own item); "" if there is none.
std::string MachineTiers_VariantItem(PartType t, int step);
// An item placing a machine variant: its part type and step; false for anything else (base items included: step 0).
bool MachineTiers_ItemPart(const std::string& item, PartType& type, int& step);
// The variant items in items.txt's format (ItemDB::loadExtraText): their names, colours and descriptions from the base's.
std::string MachineTiers_ItemsText();

// Overclock steps for a recipe of tier `recipeTier` in a machine of tier `machineTier` (0 at its tier; -1: too high to run).
int MachineTiers_Steps(int machineTier, int recipeTier);
int MachineTiers_Ticks(int ticks, int steps);      // max(1, ticks / speed^steps)
int MachineTiers_Eu(int eu, int steps);            // eu x eu-multiplier^steps

// (L7c) The disassembly yield in percent for an item `itemTier` taken apart by an assembler of `asmTier`.
int MachineTiers_YieldPercent(int asmTier, int itemTier);
