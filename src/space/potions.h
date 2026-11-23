#pragma once

// Potions (#1254, alchemy v1; the magic doc sections 4 and 6), behind RAYLIBGAME_MAGIC. A cauldron brews them (alchemy.h);
// they are drunk from the hotbar (its key, like the oxygen canister). What each does comes from alchemy.txt: heal (health at once),
// night vision (a light at the eye), water breathing (no oxygen lost under water), cold resistance (the cold costs no extra
// food), mining haste (faster strokes), and the cleansing draught (burden down by its amount, once per in-game day). A timed
// effect drunk again starts its time over. The local player's for now; their times and the draught's day are saved with
// the magic lines.

#include <map>
#include <string>

struct MagicPlayerState;
struct SurvivalState;
struct SuitStats;
struct PotionDef;

struct PotionState {
    std::map<int, float> left;   // PotionEffect -> seconds left
    long long cleansedDay = -1;  // the in-game day the cleansing draught last worked
};

// The rules (pure). Drinking: false (nothing used) when it can't work now (the draught again the same day, healing at full
// health); `why` says why.
bool Potion_Drink(const PotionDef& p, PotionState& ps, MagicPlayerState& magic, SurvivalState& survival, long long day, std::string* why = nullptr);
void Potion_Tick(PotionState& ps, float dt);
float Potion_Left(const PotionState& ps, int effect);
// A suit's figures under the local player's timed effects: water breathing (no extra drain under water: 0), cold (none).
SuitStats Potion_Suit(const PotionState& ps, const SuitStats& suit);
float Potion_StrokeScale(const PotionState& ps); // mining haste's factor (1 without it)
std::string Potion_SaveText(const PotionState& ps);
bool Potion_LoadLine(PotionState& ps, const char* line); // "magic_potion ..." lines

// The game (the local player).
PotionState& Potion_Local();
bool Potion_IsPotion(const std::string& item);
bool Potion_UseLocal(const std::string& item);    // drink one from the inventory (true: drunk)
void Potion_Update(float dt);                      // the timers; mining haste's hook
bool Potion_NightVision();                         // a light at the eye this frame
long long Potion_Day();                            // the in-game day now (the nearest body's day length)
