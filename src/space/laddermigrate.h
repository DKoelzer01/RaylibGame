#pragma once

// Save migration to the new tier ladder (#1272, ladder step L7; the factory and progression doc, section 7.7, decided
// 2026-11-20): nobody loses progress. A team whose research has no "ladder_migrated" mark is migrated once, when the world
// loads, and marked:
//  - its quest book's tier maps to the new ladder: 0, 1 the same; old 2 (Bronze and steam) becomes new 3 if steel was made,
//    else new 2 - with its milestone granted (new 3) when a boiler or forge hammer was placed, so a player mid-way through
//    the old tier 2 is never sent back to build a water wheel; 3-7 the same (tiers 8 and up are untouched);
//  - everything below the new tier is granted: its tasks, and from tier 3 every node (core and branch) and milestone line;
//  - inside the new tier: every node whose machine or recipe the team already placed or made, so nothing built is locked.
// The 8 V steam machines that are gone (forge hammer, alloy smelter, steam macerator) are dealt with by one swappable function
// (LadderMigrate_RemovedMachineHook): by default each becomes a chest in its place, holding the machine itself and what it held.

#include <functional>
#include <set>
#include <string>

struct QuestBookState;
struct ResearchState;
class Planetoid;
class Plot;

int LadderMigrate_MapTier(int oldTier, bool steelMade, bool steamPlaced);

struct LadderMigrateInput {
    int oldTier = 0;
    bool steelMade = false;                              // steel_ingot obtained
    std::set<std::string> placed;                        // machine items placed
    std::function<bool(const std::string&)> made;        // an item obtained (ever)
};
struct LadderMigratePlan {
    int newTier = 0;
    std::set<std::string> nodes;                         // research nodes granted
    std::set<std::pair<int, int>> milestoneLines;        // (tier, line) granted
};
LadderMigratePlan LadderMigrate_Plan(const LadderMigrateInput& in);
// Apply to a team's quest book and research: tier, tasks below it done (taskCount(tier) tasks each), nodes and milestone lines
// granted, the mark set. False (and nothing changed) when the research already carries the mark.
bool LadderMigrate_Apply(QuestBookState& book, ResearchState& research, const LadderMigratePlan& plan, const std::function<int(int)>& taskCount);

// The whole world after a load (main thread, from the save's load-done hook): every team's book and research, and every plot's
// removed machines (through the hook below). Returns how many teams were migrated.
int LadderMigrate_World();
bool LadderMigrate_IsRemovedMachine(int partType);       // a PartType the ladder removed (8 V)
// A removed machine on a plot, after a load: the default turns it into a chest (same footprint, same part index, so the saved
// contents load into it) and puts the machine's own item in it. Swap it for another rule (#1271 question 3).
inline std::function<void(Planetoid& body, Plot& plot, int part)> LadderMigrate_RemovedMachineHook;
void LadderMigrate_RefundIntoChest(Planetoid& body, Plot& plot, int part);
