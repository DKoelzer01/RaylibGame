#ifndef SURVIVAL_H
#define SURVIVAL_H

#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <map>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include <raylib.h>

class Scene;

// Suit survival loop: oxygen drains outside the ship's life-support field, power drains with
// recharging the jetpack and the headlamp, health drops when suffocating or from hard landings. Death
// respawns you at the ship.
struct SurvivalState {
    float oxygen = 100.0f, power = 100.0f, health = 100.0f;   // 0..100
    float fallDamageTaken = 0.0f;   // health lost to hard landings so far (not saved; tests tell it from other damage, #439)
    const char* lastHarm = "";      // (#844) what last took health: fall, hit, vacuum, suffocation, starvation, heat (not saved; logs)
    float lastFallSpeed = 0.0f;     // (#844) m/s into the ground at the last harmful landing (logs)
    float oxygenDrainPerSec = 0.35f;                          // ~5 minutes per tank
    float powerDrainPerSec = 0.05f;                           // life support baseline
    float jetpackPowerPerTank = 25.0f;                        // suit power to recharge an empty jetpack tank
    float lastJetpackFuel = 1.0f;                             // fuel at the end of the last update (recharge tracking)
    float suffocationDamagePerSec = 6.0f;
    float lifeSupportRefillPerSec = 12.0f;                    // inside / near the ship
    float headlampPowerPerSec = 0.15f;                        // ~11 min of lamp on a full suit
    bool headlamp = false;                                    // L toggles; goes out at 0 power
    // Hunger: the food bar empties in ~45 min of play (faster running or on the jetpack, and in the
    // cold); at 0 health falls slowly, above 80 it heals.
    float food = 100.0f;
    float foodDrainPerSec = 100.0f / (45.0f * 60.0f);
    float starvationDamagePerSec = 0.5f;
    float wellFedHealPerSec = 0.25f;
    // The air around you (deg C): latitude, height, day and night, shelter underground; NaN in
    // vacuum (the suit regulates). Below 5 C food drains faster; above 35 C slow damage.
    float temperature = 20.0f;
    float sickUntil = 0.0f;                                   // raw meat: food drains faster for a while
    bool inLifeSupport = false;
    float breathableAltitude = 15000.0f;                      // a planet's air is breathable below this (m above sea level)
    float airRefillPerSec = 4.0f;                             // suit oxygen top-up from breathable air
    bool breathable = false;
    bool inGas = false;                                       // in a dungeon's poison gas (#719): not breathable, oxygen drains x3
    bool underwater = false;                                  // the head is in water (oxygen drains faster)
    bool inVacuum = false;                                    // outside any breathable air (space, the moon, ring rocks)
    float vacuumTime = 0.0f;                                  // seconds in vacuum in a suit not rated for it (60 s, then it hurts)
    // Stamina (#165, combat.cpp): sprinting, swings and blocked hits spend it; it comes back after a short rest.
    float stamina = 100.0f;                                   // 0..100
    double staminaUsedAt = -1e9;                              // sim time of the last spend (regen waits a moment)
    bool blocking = false;                                    // holding a block this step
    float hitDamageTaken = 0.0f;                              // health lost to hits so far (not saved; tests)
    double lastHitAt = -1e9;                                  // sim time of the last hit taken (the camera shake, hit-stop)
    float lastHitStrength = 0.0f;                             // 0..1 (20 damage = 1)
    bool dead = false;
    double deathTime = 0.0;
    int deaths = 0;
    // The bed this player set as their spawn (#800): body name, plot id and the bed's first cell. Respawn goes there while a
    // bed still stands in that cell with room to stand on it; otherwise at the start (and the bed is forgotten).
    bool hasBed = false;
    std::string bedBody;
    uint32_t bedPlot = 0;
    int bedCell[3] = { 0, 0, 0 };
};
extern SurvivalState& survival; // the local PlayerSession's (#313)

// The world's difficulty (saved with it; new worlds take settings.difficulty): 0 peaceful (no night
// threats or hostile animals, hunger x0.6, extra food to start), 1 normal, 2 hard (no Map in the
// starting kit, hunger x1.25, threats bite harder and come more often).
extern int worldDifficulty;
float Difficulty_HungerScale();
float Difficulty_ThreatScale(); // 0 peaceful, 1 normal, 1.5 hard
const char* Difficulty_Name(int d);

struct PlayerInput;
void Survival_Update(Scene& world, float dt, const PlayerInput& in);
// Consumables (the hotbar uses them): oxygen canister, power cell (charges the ship while piloting).
bool Survival_IsUsable(const std::string& id);
bool Survival_UseItem(const std::string& id);
// #337: a client process's item use goes to the server (set by the join); a server applies one for a remote player.
inline std::function<void(const std::string& id)> gSurvivalUseRequest;
struct SuitStats; class Inventory;
bool Survival_UseItemFor(const std::string& id, SurvivalState& survival, Inventory& inventory, const SuitStats& stats, bool local);
bool Survival_IsFood(const std::string& id);  // resources/data/foods.txt
float Survival_FoodValue(const std::string& id);
void Survival_OnLanding(float impactSpeed);   // fall damage
void Survival_Respawn(Scene& world);
// #337: the same for a given player (a server's remote players: local false - no UI, sounds, quests, rocket or ship; the weather
// and temperature where that player is). Survival_Update / _Respawn / _OnLanding are these for the local player.
struct Player; class Inventory; struct SuitStats;
void Survival_StepFor(Scene& world, float dt, const PlayerInput& in, SurvivalState& survival, Player& player, Inventory& inventory, const SuitStats& stats, bool local);
void Survival_RespawnFor(Scene& world, SurvivalState& survival, Player& player, bool local);
// Beds (#800). SetBedFor: the part at (body, plot, part) is a bed: make it this player's spawn (false: not a bed). BedSpawnFor: where
// the player would wake on their bed (render-space feet on top of it), false when it is gone or something stands on it.
class Planetoid; class Plot;
bool Survival_SetBedFor(SurvivalState& survival, const Planetoid& body, const Plot& plot, int part, bool local);
bool Survival_BedSpawnFor(const SurvivalState& survival, Planetoid*& body, Vector3& renderFeet);
// A joined game (#800): right-clicking a bed asks the server (main.cpp sets it); unset, the bed is set locally.
inline std::function<void(const Planetoid& body, const Plot& plot, int part)> gBedUseRequest;
void Survival_OnLandingFor(SurvivalState& survival, const SuitStats& stats, float impactSpeed, bool local);
void Survival_DrawHUD();
const char* Survival_BiomeName();   // the biome where the player stands ("" off-planet)

// Lightweight objective tracker for the vertical slice.
void Objectives_Update(Scene& world);
void Objectives_DrawHUD();
// Quest book (tiers.txt): current tier, advanced when its gate item is made or placed.
// One quest book's progress (#338): single player has one (the local game's); a server keeps one per team (quests and tiers
// belong to teams, user 2026-11-16 on #303). Every QuestBook_* / QuestTree_* call works on the selected book.
struct QuestBookState {
    int tier = 0;
    std::set<std::string> placed;                    // items placed ("place" tasks, tier gates)
    std::set<int> visitedBodies;
    bool reachedOrbit = false, visitedMoon = false, visitedRing = false;
    std::map<std::string, double> counts;            // events for "do" tasks: eat, feed, penned, till, plant, harvest, ...
    std::set<std::pair<int, int>> tasksDone;         // main line: (tier, task)
    std::set<std::string> treeDone;                  // side lines: node ids (questtree.cpp)
    std::map<int, double> tierSeconds;               // (#1199) the tier clock: seconds of play in each tier (sim time: not loading or paused)
    // "get" tasks count what was obtained. The local book reads the local inventory's totals (saved with them, as before);
    // a team book keeps its own, the sum of what its members obtained on the server.
    bool ownObtained = false;
    std::unordered_map<std::string, double> obtained;
};
QuestBookState& QuestBook_LocalState();
void QuestBook_Select(QuestBookState* book);         // nullptr: the local one again
QuestBookState& QuestBook_Selected();
// Rewards and announcements (#816): unset, a finished side-line node pays the local inventory and toasts; a server sets
// these while it evaluates a team's book (the completing player is paid through Host_Give; the team hears the news).
inline std::function<void(const std::string& item, int count)> gQuestReward;
inline std::function<void(const std::string& text)> gQuestNews;
void QuestBook_Update();
void QuestBook_NotePlaced(const std::string& item);
// Events counted for "do <event> <n>" tasks (eat, feed, penned, till, plant, harvest, waypoint, modules, explored).
void QuestBook_Count(const std::string& event, double n = 1.0);
void QuestBook_AtLeast(const std::string& event, double value); // a running maximum (penned animals, modules fitted)
double QuestBook_CountOf(const std::string& event);
// Whether a task spec ("get iron_ingot 3", "place kiln", "do eat 1", ...) is met now (the quest tree).
bool QuestBook_Met(const std::string& kind, const std::string& target, double n);
int QuestBook_Tier();
void QuestBook_SetTier(int tier);
// (#1199) The tier clock: call once a frame with the sim time (negative while loading); the time since the last call (none after a load or a jump
// back) goes to the current tier. QuestBook_TierSeconds: that tier's seconds so far.
void QuestBook_TickClock(double simTime);
double QuestBook_TierSeconds(int tier);
int QuestBook_TierReport(const char* savePath); // (#1199) --tier-report: prints each quest book's hours per tier
void QuestBook_Reset();                         // a new game: tier 0, nothing done
void QuestBook_Save(FILE* f);                   // tier, done tasks, placed machines, obtained totals
bool QuestBook_LoadLine(const char* tok, FILE* f);
// The tiers and their tasks (HUD, the Quests tab).
int QuestBook_TierCount();
const char* QuestBook_TierName(int tier);
const char* QuestBook_TierText(int tier);
int QuestBook_TaskCount(int tier);
const char* QuestBook_TaskText(int tier, int task);
bool QuestBook_TaskDone(int tier, int task);
// Pinned quests (#825): up to kMaxPins tasks shown on the HUD; the player's own (their record), not the team's.
constexpr int kMaxPins = 3;
bool QuestBook_Pin(int tier, int task, bool on);   // false: already kMaxPins pinned (or no such task)
bool QuestBook_Pinned(int tier, int task);
const std::vector<std::pair<int, int>>& QuestBook_Pins();
void QuestBook_ClearPins();
void QuestBook_SavePins(FILE* f);
bool QuestBook_LoadPinLine(const char* tok, FILE* f);
void Objectives_MarkVisitedBody(int bodyIndex);

#endif
