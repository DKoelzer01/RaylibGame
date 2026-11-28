#ifndef SUIT_H
#define SUIT_H

// Environmental suit (plan 2, milestone 13). The player wears one suit (resources/data/suits.txt:
// basic, field, explorer) with 2-6 module slots; modules (resources/data/suit_modules.txt) each
// answer a hazard or a chore. Suits and modules are tool items (never stack); the Suit tab of the
// inventory fits them. Suit_Apply() turns the worn setup into the player's and survival's numbers
// every frame.

#include <raylib.h>
#include <cstdio>
#include <string>
#include <vector>

class Scene;

struct SuitStats {
    int slots = 2;
    float oxygen = 1.0f, power = 1.0f;   // capacity multipliers (the bars drain that much slower)
    float jetFuel = 1.0f, jetThrust = 1.0f, walk = 1.0f, lamp = 1.0f, food = 1.0f;
    float underwater = 2.0f;             // oxygen drain multiplier underwater (rebreather: 0.5)
    float cold = 1.0f;                   // share of the cold's extra hunger felt (insulation 0.4, regulator 0)
    bool heatProof = false;              // thermal regulator: no heat damage
    // (#1472, #1465) the heat rating: the depth bands' heat (0 above 250 m, 1 to 500 m, 2 below) this suit stands. The best
    // module wins, no stacking: a cooling vest 1, the thermal regulator 2 with suit power (1 without).
    int heat = 0;
    bool grip = false, scanner = false, servo = false;
    bool gasSensor = false;              // (#1472) cave gas warnings from 15 m, with a direction
    bool hover = false;                  // flight pack (#163): hover mode (H)
    bool powered = false;                // a battery pack fitted: suit power at all (#217: the basic suit has none)
    bool jetpack = false;                // a thruster pack fitted: a jetpack at all (#217: flight comes at MV)
    bool vacuum = false;                 // rated for vacuum (the explorer suit)
};

// One player's suit (#313, in their PlayerSession): what is worn and fitted, its numbers, and the ore scanner.
struct SuitState {
    std::string worn = "basic_suit";
    std::vector<std::string> fitted = std::vector<std::string>(2);
    SuitStats stats;
    // The ore scanner: the nearest deposit within 64 m, looked for once a second.
    std::string scanKind;
    float scanDist = -1.0f;
    Vector3 scanLocal = { 0, 0, 0 };
    double lastScan = -10.0;
};

const SuitStats& Suit_Stats();
SuitStats Suit_StatsOf(const SuitState& suit); // (#337) any player's (a server's remote players)
int Suit_HeatRating(const SuitState& suit);    // (#1472) SuitStats::heat of any player's suit: the depth heat it stands
// (#1840) Packs fitted modules forward into `slots` slots; returns the ones that don't fit (each once), for the inventory.
std::vector<std::string> Suit_PackSlots(std::vector<std::string>& fitted, int slots);
const std::string& Suit_Worn();                       // the worn suit's item id
const std::vector<std::string>& Suit_Modules();       // one per slot, "" empty
bool Suit_IsSuit(const std::string& id);
bool Suit_IsModule(const std::string& id);
const char* Suit_ModuleText(const std::string& id);   // its effect, for the UI

bool Suit_Wear(const std::string& suitId);            // from the inventory; the old suit goes back, modules that don't fit too
bool Suit_Fit(const std::string& moduleId);           // into the first free slot, from the inventory
bool Suit_Unfit(int slot);                            // back to the inventory

void Suit_Apply(Scene& world, float dt);              // stats into player/survival; the ore scanner
// (#1087) A suit's stats into any player's movement: walk speed (x walkMul: Athletics), the jetpack (none without a
// thruster pack and power, #217: no thrust, an empty tank), hover, grip. Suit_Apply for the local player; a server for
// each remote player before its step (it had left them the Player defaults: a full, working jetpack).
struct Player;
void Suit_ApplyStatsTo(Player& p, const SuitStats& st, float walkMul = 1.0f);
void Suit_DrawHUD();                                  // the scanner's ping
bool Suit_FitInto(int slot, const std::string& moduleId); // (#822) into that slot (a fitted one comes out first)
void Suit_DrawCarried();                              // (#822) carried suits and modules, Wear / Fit (ImGui)
void Suit_Reset();                                    // a new game: the basic suit, empty slots
void Suit_Save(FILE* f);
bool Suit_LoadLine(const char* tok, FILE* f);

#endif
