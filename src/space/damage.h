#ifndef DAMAGE_H
#define DAMAGE_H

// Combat step 1 (#1215; combat doc section 14): damage types, resistances and the five status effects, one path for players,
// creatures, turrets and traps. Pure functions over small state structs (no globals of a player), so a server resolves any
// target the same way. The numbers are named constants, starting values to tune in play.
//
// Damage taken = amount x (1 - resistance[type]) - flat armour, never below 10% of the amount; poison's ticks ignore armour.
// Status effects (decided 2026-11-21):
//   burn   (heat)                     3/s for 4 s; water puts it out
//   chill  (cold)                     -30% move and attack speed for 3 s; three stacks freeze for 1 s
//   bleed  (slash heavies, back hits) 2/s for 6 s per stack, stacks to 3
//   poison (bites, ooze)              1.5/s for 10 s, ignores armour; a herb cures it
//   shock  (energy)                   interrupts a wind-up; x2 damage in water
// RAYLIBGAME_COMBAT1=0: off (types, resistances and effects ignored: today's damage).

#include <cstdint>

enum class DamageType : uint8_t { Blunt, Slash, Pierce, Heat, Cold, Energy, Count };
enum class StatusKind : uint8_t { None, Burn, Chill, Bleed, Poison, Shock };

const char* DamageType_Name(DamageType t);
const char* StatusKind_Name(StatusKind s);

namespace StatusTune {
    constexpr float kMinDamageShare = 0.10f;   // a hit never does less than this share of its amount
    constexpr float kBurnDps = 3.0f, kBurnTime = 4.0f;
    constexpr float kChillSlow = 0.30f, kChillTime = 3.0f, kFreezeTime = 1.0f;
    constexpr int kChillFreezeStacks = 3;
    constexpr float kBleedDps = 2.0f, kBleedTime = 6.0f;
    constexpr int kBleedMaxStacks = 3;
    constexpr float kPoisonDps = 1.5f, kPoisonTime = 10.0f;
    constexpr float kShockWaterFactor = 2.0f;
}

bool Combat1_On(); // RAYLIBGAME_COMBAT1=0: off (A/B)

struct Resistances {
    float r[(int)DamageType::Count] = {}; // 0..1 per type (1: immune)
    float armour = 0.0f;                  // flat, subtracted after the resistance
};

// The damage a hit of `amount` of `type` does through `res` (the 10% floor; ignoresArmour: poison).
float Damage_Apply(float amount, DamageType type, const Resistances& res, bool ignoresArmour = false);

struct StatusState {
    float burn = 0.0f;                 // seconds left
    float chill = 0.0f; int chillStacks = 0; float frozen = 0.0f;
    float bleed = 0.0f; int bleedStacks = 0;
    float poison = 0.0f;
    bool interrupted = false;          // a shock this step (Status_Step clears it; a wind-up checks it, #1216)
    bool any() const { return burn > 0.0f || chill > 0.0f || frozen > 0.0f || bleed > 0.0f || poison > 0.0f; }
};

// The effect a hit brings (explicit, else from its type: heat burn, cold chill, energy shock, a slash heavy or a back hit
// bleed; a bite poisons).
StatusKind Status_FromHit(DamageType type, StatusKind explicitStatus, bool heavy, bool fromBehind, bool bite);
// Applies an effect: refreshes its time, adds a stack where it stacks (three chills freeze).
void Status_Apply(StatusState& s, StatusKind k);
// One step: the damage the effects deal this step. Each rate goes through the damage formula per second (burn as heat,
// bleed as slash, both with armour; poison ignores armour), then x dt, so flat armour never wipes out a tick. Water puts
// burn out. Clears `interrupted` (a shock's wind-up interrupt lasts one step).
float Status_Step(StatusState& s, float dt, bool inWater, const Resistances& res);
// Move and attack speed multipliers: chilled 0.7, frozen 0.
float Status_SpeedFactor(const StatusState& s);
void Status_CurePoison(StatusState& s);

#endif
