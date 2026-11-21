#pragma once

// Research (#1030; the factory and progression doc, section 7.6): from tier 3 a tree of nodes researched in labs with packs.
// This is the data model: the nodes (resources/data/research.txt), a team's research state (done, queue, progress, insight)
// kept as save lines in the team record, the join rule and what a player may see. Labs, packs, the Research tab, tier
// milestones and the recipe gating come in their own steps (doc section 8); nothing is gated by research yet.

#include <map>
#include <set>
#include <string>
#include <vector>

struct ResearchNode {
    std::string id, name;
    int tier = 3;                      // 3 Steam .. 7 Rocketry
    bool core = true;                  // required for its tier's milestone; else a branch on `line`
    std::string line;                  // a branch's line: combat, logistics, efficiency, magic, mining, exploration
    int units = 1;                     // research units; each takes one of every pack in `packs` and unitSeconds in one lab
    float unitSeconds = 10.0f;
    std::vector<std::string> packs;    // pack items consumed per unit
    std::vector<std::string> needs;    // prerequisite nodes
    std::vector<std::string> unlocks;  // "recipe:<output>", "machine:<part>"
    std::string tome;                  // a magic cross link: the tome page it needs ("" none)
};

// The node table: loaded once from resources/data/research.txt (or a test's text).
bool Research_LoadText(const std::string& text, std::string& error);
bool Research_Load(const char* path = "resources/data/research.txt");
const std::vector<ResearchNode>& Research_Nodes();
const ResearchNode* Research_Find(const std::string& id);
// Checks: unique ids, every prerequisite exists, no cycles, a node needs nothing of a later tier, its packs are its tier's
// kind and earlier ones, no core node needs a tome page. The problems, one a line ("" none).
std::string Research_Validate();
const char* Research_PackOfTier(int tier); // "steam_pack" (3) .. "rocket_pack" (7); "" outside

// A team's research.
struct ResearchState {
    std::set<std::string> done;
    std::vector<std::string> queue;          // in order; labs take the first they have the packs for
    std::map<std::string, int> progress;     // units finished of a node not yet done (kept when it leaves the queue)
    double insight = 0.0;                    // activity XP waiting to be drawn on (up to ~25% of a node, decided)
    std::set<std::pair<int, int>> milestoneDone; // (#1271) tier milestone lines met: (tier, line index in tiers.txt order)
    bool ladderMigrated = false;             // (#1272) an old save's tiers and research were mapped to the new ladder (once)
    std::string Text() const;                // as save lines
    bool FromText(const std::string& text);  // false: nothing parsed
    bool Done(const std::string& node) const { return done.count(node) > 0; }
};
bool Research_CanStart(const ResearchState& s, const ResearchNode& n); // not done, and every prerequisite done
// What the team may see of a node (decided): its current tier in full, the next tier's names only (no costs), later hidden.
// Done nodes are always shown.
enum class ResearchSight { Hidden, NameOnly, Full };
ResearchSight Research_Sight(const ResearchState& s, const ResearchNode& n, int teamTier);
// A player joining a team: the team gets what the joiner had (decided: the team's research, and leaving keeps yours):
// done is the union, progress the larger, insight the larger; the team's queue stays.
void Research_Merge(const ResearchState& from, ResearchState& into);
// (#1193) Labs: the first queued node a lab of `maxTier` can run (startable, its tier at most the lab's; nullptr none), and
// `units` research units of `tier` added to the first such queued node of that tier. A node with all its units is done
// and leaves the queue; the names of the nodes done are returned.
const ResearchNode* Research_HeadFor(const ResearchState& s, int maxTier);
std::vector<std::string> Research_AddUnits(ResearchState& s, int tier, int units);

// The team record's block (saved with the world as team_research <id> ... team_research_end).
ResearchState& Research_Team(int team);
void Research_Reset();
std::string Research_TeamText(int team);       // what the team record keeps
void Research_SetTeamText(int team, const std::string& text);
