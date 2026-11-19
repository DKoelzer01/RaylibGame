#pragma once

// Skills (#167, #883): use-based XP. Doing a thing earns its skill XP; each level gives a point to spend in that skill's
// perk tree (resources/data/skills.txt). Perks multiply stats (+10% a rank, capped by ranks), stacking with the suit's
// modules; a few are abilities. Skills never gate tech and have no cap. The user's answers on #167 (2026-11-18): respec with a
// craftable item, death costs some XP toward the next level (never a level), an XP-rate world setting, 10% team sharing.
//
// The model holds no UI and no game state beyond one player's SkillState: the game feeds events (Skills_Event) and reads
// multipliers (Skills_Stat); the Skills tab draws from it; saves and player records keep it as text lines.

#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

struct SkillDef { std::string id, name; std::map<std::string, double> xp; }; // xp: per event
struct PerkDef {
    std::string id, skill, name, text;
    int tier = 1, cost = 1, ranks = 1;
    std::vector<std::string> needs;                       // perks it needs (a rank of each)
    std::vector<std::pair<std::string, double>> effects;  // stat, per rank
    std::string ability;                                  // "" none
};

// One player's skills: XP per skill and perk ranks. Levels and points follow from them.
struct SkillState {
    std::map<std::string, double> xp;                     // by skill id
    std::map<std::string, int> ranks;                     // by perk id
    // Anti-farming (not saved): ground this player filled (a cell key) and when, so filling a hole and digging it out again
    // pays little (digging on into the same spot, a hole going down, pays in full).
    std::unordered_map<uint64_t, double> filledCells;
};

// The tables (resources/data/skills.txt), loaded once. Tests may load their own text.
bool Skills_Load(const char* path = "resources/data/skills.txt");
bool Skills_LoadText(const std::string& text);
const std::vector<SkillDef>& Skills_All();
const std::vector<PerkDef>& Skills_Perks();
const SkillDef* Skills_Find(const std::string& id);
const PerkDef* Skills_FindPerk(const std::string& id);

// The curve: going from level n-1 to n costs 100 x n^1.4 XP (level 10 ~2,500; 30 ~11,700); no cap.
double Skills_XpForLevel(int n);                          // the step to reach n (n >= 1)
int Skills_LevelOf(double xp);                            // the level a total XP gives
double Skills_Progress(double xp);                        // 0..1 toward the next level
constexpr int kSkillTierLevel[3] = { 5, 15, 30 };         // tiers open at these skill levels

// XP. rate: the world's XP-rate setting (0.5 / 1 / 2 / 5). Returns the levels gained (for a toast).
int Skills_AddXp(SkillState& s, const std::string& skill, double xp, double rate = 1.0);
// An event the game counted (rock_stroke, ore_stroke, gem, chop, ...): every skill that pays for it gains its XP x n.
// Returns the levels gained across skills; levelled (if given) gets each skill that levelled.
int Skills_Event(SkillState& s, const std::string& event, double n = 1.0, double rate = 1.0, std::vector<std::string>* levelled = nullptr);
// Cells are any stable key (a 1 m cell's coordinates packed); now in seconds. A fill notes its cell; a dig stroke pays
// 0.1 of its XP in a cell this player filled within the last 10 minutes, else all of it.
void Skills_NoteFill(SkillState& s, uint64_t cell, double now);
double Skills_DigShare(const SkillState& s, uint64_t cell, double now);
// Death (Q2): loses `share` of the progress toward the next level in every skill, never a level.
void Skills_OnDeath(SkillState& s, double share = 0.25);

int Skills_Level(const SkillState& s, const std::string& skill);
int Skills_PointsSpent(const SkillState& s, const std::string& skill);
int Skills_PointsFree(const SkillState& s, const std::string& skill); // its level, less what its tree has spent

// Perks. CanBuy explains why not ("" when it can).
std::string Skills_CanBuy(const SkillState& s, const std::string& perk);
bool Skills_Buy(SkillState& s, const std::string& perk);
void Skills_Respec(SkillState& s, const std::string& skill); // every point of that tree back (the cost is the caller's: an item)

// A stat's multiplier: 1 + the sum over perks of ranks x per rank (negative for wear: 1 - 0.1 x ranks, never below 0.1).
// Stats in slice 1: rock_speed, ore_speed, ore_yield, stone_yield, gem_chance, pick_wear, rock_radius. (#919) Slice 2a:
// soil_speed, soil_radius, shovel_wear (Digging); chop_speed, log_yield, axe_wear (Woodcutting); fibre_yield, stone_find
// (Foraging); crop_yield, seed_yield (Farming); break_speed (Building); walk_speed, fall_safe (Athletics); hunger, cold
// (Survival). All act for the local player only until the server reads them (#920).
double Skills_Stat(const SkillState& s, const std::string& stat);

// The world's XP rate (Q6: 0.5, 1, 2 or 5; saved with the world as "xp_rate", new worlds take the New World screen's).
inline float worldXpRate = 1.0f;

// Saves and player records: "skill_xp <id> <xp>" and "perk <id> <rank>" lines.
std::string Skills_SaveText(const SkillState& s);
bool Skills_LoadLine(SkillState& s, const char* line);     // true when it was one of ours
