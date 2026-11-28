#pragma once

// Strain and burden (#1255, magic step 8; the magic doc section 6). Strain is short-term: big spells, transmutation by hand
// and rituals add it, it drains over minutes (faster asleep), and at the top comes a backlash (damage that never kills,
// an emptied pool). Burden lasts (Thaumcraft's warp; tome flavour: Onus): a quarter of the strain gained while over half
// strain, and the strongest acts directly. It barely drains by itself (1 a day without casting), dying never clears it
// (decided), and it is cured at a shrine, by the cleansing draught or a purification ritual. Symptoms rise in bands.
// Never unfair: symptoms never touch builds, machines, items or other players, never kill, and pause while loading,
// paused, in a menu, asleep, in a fight or inside a shrine's ward. Behind RAYLIBGAME_MAGIC.

#include <string>

struct MagicPlayerState;

enum class BurdenBand { Clear, Unsettled, Marked, Frayed, Overburdened };

BurdenBand Burden_Band(float burden);              // 0-19 Clear .. 80-100 Overburdened
const char* Burden_BandName(BurdenBand b);
const char* Burden_BandSymptoms(BurdenBand b);     // what it does (the HUD tooltip)
float Burden_RegenFactor(float burden);            // pool regeneration (Unsettled and up: -10%)
float Burden_CostFactor(float burden);             // mana costs (Marked and up: +10%)
float Burden_MisfireChance(float burden, int spellParts); // big spells (5+ parts) at Frayed and up
float Burden_BacklashAt(float burden);             // the strain a backlash comes at (100; Overburdened: 85)

// Strain. Returns true if it set off a backlash (the caller applies the damage: Burden_BacklashHealth).
bool Strain_Add(MagicPlayerState& s, float amount);
void Strain_Drain(MagicPlayerState& s, float dt, bool asleep);     // over minutes; asleep 3x; Frayed and up 25% slower
float Burden_BacklashHealth(float health, float maxHealth);         // health after a backlash: never below a quarter

// Burden.
void Burden_Add(MagicPlayerState& s, float amount);                // a strong act directly (1-5 each); the charm cuts it
void Burden_DayPassed(MagicPlayerState& s, bool castToday);        // 1 a day without casting
void Burden_Cure(MagicPlayerState& s, float amount);               // a night at a shrine (10), the draught (15), the ritual (50)
bool Burden_LatePagesShown(const MagicPlayerState& s);             // the optional Void and Rot mastery pages (burden 20+)
