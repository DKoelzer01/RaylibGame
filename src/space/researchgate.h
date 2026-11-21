#pragma once

// Research gating and tier milestones (#1271, ladder step L1; the factory and progression doc, section 7.6).
//
// INERT in this step: nothing is gated unless the gate is switched on (RAYLIBGAME_LADDER=1, or ResearchGate_SetOn in tests).
// Step L6 switches it on, in the same series as the tier rewrite and the save migration (L7), so no build gates recipes
// without migrating saves.
//
// Gating: a recipe or a machine that a research node unlocks (research.txt "unlocks recipe:<output> machine:<item>") is
// allowed once the team has that node done. Things no node names are not gated here (tiers 0-2 unlock with their tier's
// milestone: step L5). Hooks: RecipeBook_Unlocked (the crafting screens and the server's craft check) and BuildTool_Place.
//
// Tier milestones (tiers 3-7) come from tiers.txt, as indented lines under a tier:
//     milestone | deliver <item> <n>            | text     (consumed at a lab)
//     milestone | build lab_node                | text     (a boiler and steam engine feeding a lab that completes a node)
//     milestone | build throughput <item> <per minute> <minutes> [or throughput ...] | text   (sustained, from the statistics)
//     milestone | build rocket_on_pad           | text
//     milestone | find <find id>                | text     (coal_seam, crystal_sample, vent_sample, seafloor_survey,
//                                                           moon_survey, sounding_rocket; a ruin schematic replaces a find
//                                                           whose feature the world lacks: step L6)
// Every line is required; " or " separates alternatives within a line. A tier's milestone is met when all of its core nodes
// are done and every line is (the team's research state records the lines done, in any order).

#include "plots.h"
#include "research.h"
#include <string>
#include <vector>

bool ResearchGate_On();
void ResearchGate_SetOn(bool on);

// The node whose unlocks name it, or nullptr (not gated by research).
const ResearchNode* ResearchGate_NodeForRecipe(const std::string& output);
const ResearchNode* ResearchGate_NodeForPart(PartType t); // by its item (PartCostItem) or its machine name
// Allowed: the gate is off, no node names it, or its node is done.
bool ResearchGate_RecipeAllowed(const ResearchState& team, const std::string& output);
bool ResearchGate_PartAllowed(const ResearchState& team, PartType t);
// The research that gates now: a selected team's (a server checking a player's craft), else the local player's team.
void ResearchGate_SelectTeam(int team); // 0: the local player's again
const ResearchState& ResearchGate_Team();

// ---- Tier milestones ----
struct MilestoneOption {
    std::string kind;          // deliver, build, find
    std::string what;          // deliver: the item; build: lab_node | throughput | rocket_on_pad; find: the find id
    std::string item;          // build throughput: the item made
    int count = 0;             // deliver: how many; build throughput: per minute
    int minutes = 0;           // build throughput: the window
};
struct MilestoneLine { std::vector<MilestoneOption> options; std::string text; }; // any one option meets the line
struct TierMilestone { int tier = 0; std::vector<MilestoneLine> lines; };

bool ResearchGate_LoadMilestonesText(const std::string& text, std::string& error);
bool ResearchGate_LoadMilestones(const char* path = "resources/data/tiers.txt");
const TierMilestone* ResearchGate_Milestone(int tier); // nullptr: none (tiers 0-2 end on tasks)
// Checks: tiers 3-7 each have a milestone with at least one delivery, one build and one find; known build kinds and find
// ids; counts, rates and windows over zero; no milestone outside tiers 3-7. The problems, one a line ("" none).
std::string ResearchGate_ValidateMilestones();
bool ResearchGate_CoreDone(const ResearchState& s, int tier);     // every core node of the tier done
bool ResearchGate_MilestoneMet(const ResearchState& s, int tier); // core done and every milestone line done
