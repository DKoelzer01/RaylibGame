#ifndef SURVIVAL_H
#define SURVIVAL_H

#include <cstdio>
#include <string>

class Scene;

// Suit survival loop: oxygen drains outside the ship's life-support field, power drains with
// recharging the jetpack and the headlamp, health drops when suffocating or from hard landings. Death
// respawns you at the ship.
struct SurvivalState {
    float oxygen = 100.0f, power = 100.0f, health = 100.0f;   // 0..100
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
    bool underwater = false;                                  // the head is in water (oxygen drains faster)
    bool inVacuum = false;                                    // outside any breathable air (space, the moon, ring rocks)
    float vacuumTime = 0.0f;                                  // seconds in vacuum in a suit not rated for it (60 s, then it hurts)
    bool dead = false;
    double deathTime = 0.0;
    int deaths = 0;
};
extern SurvivalState survival;

// The world's difficulty (saved with it; new worlds take settings.difficulty): 0 peaceful (no night
// threats or hostile animals, hunger x0.6, extra food to start), 1 normal, 2 hard (no Map in the
// starting kit, hunger x1.25, threats bite harder and come more often).
extern int worldDifficulty;
float Difficulty_HungerScale();
float Difficulty_ThreatScale(); // 0 peaceful, 1 normal, 1.5 hard
const char* Difficulty_Name(int d);

void Survival_Update(Scene& world, float dt, bool inputEnabled);
// Consumables (the hotbar uses them): oxygen canister, power cell (charges the ship while piloting).
bool Survival_IsUsable(const std::string& id);
bool Survival_UseItem(const std::string& id);
bool Survival_IsFood(const std::string& id);  // resources/data/foods.txt
float Survival_FoodValue(const std::string& id);
void Survival_OnLanding(float impactSpeed);   // fall damage
void Survival_Respawn(Scene& world);
void Survival_DrawHUD();
const char* Survival_BiomeName();   // the biome where the player stands ("" off-planet)

// Lightweight objective tracker for the vertical slice.
void Objectives_Update(Scene& world);
void Objectives_DrawHUD();
// Quest book (tiers.txt): current tier, advanced when its gate item is made or placed.
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
void Objectives_MarkVisitedBody(int bodyIndex);

#endif
