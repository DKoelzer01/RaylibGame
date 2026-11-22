#ifndef MELEE_H
#define MELEE_H

// Combat step 2 (#1216): melee swings in phases on the 60 Hz sim clock (combat doc 2.2-2.3).
// Every attack has a wind-up (cancellable into a block or dodge in its first half), a hit window (the weapon's arc is
// swept each step, one hit per target per swing) and a recovery (a dodge may cancel its last third). Tap LMB for a light
// attack; hold it through the light wind-up for a heavy (wind-up x1.6, damage x1.8, poise x2, stamina x1.7). Two lights
// chain (the second starts a third of the way into the first's recovery); no longer combos. Heavy and two-handed
// wind-ups have hyper-armour (a shock doesn't interrupt them). RAYLIBGAME_COMBAT2=0: the old one-click hit.
//
// This file is the logic only (no world, no input devices): animals.cpp drives it with this frame's aim.

#include "damage.h"

#include <raylib.h>
#include <cstdint>
#include <string>
#include <vector>

bool Combat2_On();

enum class WeaponClass : uint8_t { Tools, Knife, Sword, Axe, Mace, Spear, TwoHanded, Count };

struct WeaponSpec {
    const char* name;
    float windUp, hit, recovery; // s, light (the doc's table)
    float reach;                 // m
    DamageType type;
    float stamina;               // a light swing's cost
    float poise;                 // poise damage of a light hit (#1217 uses it)
    float arcDeg;                // the swept arc (a spear thrusts: narrow)
    bool twoHanded;              // hyper-armour on every wind-up
};

namespace MeleeTune {
constexpr float kHeavyWindUp = 1.6f, kHeavyDamage = 1.8f, kHeavyPoise = 2.0f, kHeavyStamina = 1.7f;
constexpr float kChainAt = 1.0f / 3.0f;     // the second light starts this far into the first's recovery
constexpr float kDodgeCancelFrom = 2.0f / 3.0f; // a dodge may cancel the recovery from here
constexpr float kBladeStart = 0.3f;         // m in front of the eye: the blade's near end
constexpr float kSweepStepDeg = 4.0f;       // blade samples at most this far apart across a step's sweep
constexpr float kHitStop = 0.07f;           // s, the view's pause on contact (reduce motion: none)
constexpr float kLockRange = 12.0f, kLockKeep = 16.0f, kLockConeDeg = 35.0f, kLockTurnTime = 0.12f;
}

const WeaponSpec& Weapon_Spec(WeaponClass c);
// The class an item fights as: knives, swords, axes (the axe tool too), maces/hammers/clubs, spears, two-handers by
// their id; anything else (hands, pickaxe, shovel, hoe) is Tools.
WeaponClass Weapon_ClassOf(const std::string& itemId);
DamageType Weapon_TypeOf(const std::string& itemId); // the class's type; a pickaxe pierces
// (#1219) A light hit's base damage with this item: a composed copy's from its parts (data: its record); otherwise hands 2,
// a tool 3 + 3 x tier (as before).
float Weapon_BaseDamage(const std::string& itemId, const std::string& data);

enum class MeleePhase : uint8_t { Idle, WindUp, Hit, Recovery };

struct MeleeState {
    MeleePhase phase = MeleePhase::Idle;
    float t = 0.0f;              // s into the phase
    bool heavy = false;
    int chain = 0;               // 0 the first light, 1 the chained second
    bool queued = false;         // a light pressed during the hit window or recovery: chains
    float sweep = 0.0f;          // 0..1 across the arc in the hit window
    std::vector<uint64_t> struck; // targets hit by this swing (one hit each)
    int swings = 0, heavies = 0; // counters (tests)
};

struct MeleeInput {
    bool pressed = false;        // LMB went down
    bool held = false;           // LMB is down
    bool block = false;          // a block started (cancels the first half of a wind-up)
    bool dodge = false;          // a dodge started (#1218: the first half of a wind-up, the last third of a recovery)
    bool interrupted = false;    // shocked this step
    float stamina = 1e9f;        // what's left
};

struct MeleeStep {
    bool started = false;        // a swing began (its stamina is in cost)
    bool becameHeavy = false;
    bool cancelled = false;      // by a block, dodge or interrupt
    bool tired = false;          // pressed without the stamina for it
    bool sweeping = false;       // the arc from sweepFrom to sweepTo (0..1) is live this step
    float sweepFrom = 0.0f, sweepTo = 0.0f;
    bool finished = false;       // back to idle after a recovery
    float cost = 0.0f;           // stamina spent this step
};

MeleeStep Melee_Step(MeleeState& s, const WeaponSpec& w, const MeleeInput& in, float dt);
bool Melee_HyperArmour(const MeleeState& s, const WeaponSpec& w);
float Melee_WindUpTime(const MeleeState& s, const WeaponSpec& w);
// The swing's angle (radians, + = to the right of the view) at sweep fraction f: right to left, the chained light back.
float Melee_SweepAngle(const MeleeState& s, const WeaponSpec& w, float f);
float Melee_Damage(float base, const MeleeState& s);
float Melee_Poise(const MeleeState& s, const WeaponSpec& w);
bool Melee_Struck(const MeleeState& s, uint64_t id);

// Geometry. The closest distance between segments p0-p1 and q0-q1 (and the point on the first).
float Melee_SegmentDistance(Vector3 p0, Vector3 p1, Vector3 q0, Vector3 q1, Vector3* onP = nullptr);
// Does the blade (eye + dir x [kBladeStart, reach], dir = forward turned about up) cross the capsule q0-q1 radius r anywhere
// between the angles a0 and a1? The contact point (on the blade) and its distance from the eye come back.
bool Melee_SweepCapsule(Vector3 eye, Vector3 forward, Vector3 up, float reach, float a0, float a1,
                        Vector3 q0, Vector3 q1, float r, Vector3* contact = nullptr, float* dist = nullptr);
// Soft lock-on: forward turned toward `toTarget` by the share of the angle a frame of dt closes (horizontal and pitch).
Vector3 Melee_SoftTurn(Vector3 forward, Vector3 up, Vector3 toTarget, float dt);

#endif
