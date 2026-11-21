#ifndef DEFENSE_H
#define DEFENSE_H

// Combat step 4 (#1218, combat doc section 3): stamina use and the winded state, shields, the parry and the dodge.
// - Stamina (one pool of 100): attacks, blocks, dodges (20) and sprinting (8/s) spend it. Run dry, the player is winded for
//   1.5 s: no heavies, no dodge, walk x0.7.
// - Shields: a buckler (stops 55% within +-45 deg, parry window 220 ms), the round wooden shield (75%, +-60, 150 ms; today's)
//   and a tower shield (95%, +-75, 100 ms, walk x0.7 while raised).
// - Parry: a block raised within the shield's window before a parryable melee hit or bite stops it completely, costs no
//   stamina and staggers the attacker (2 s, #1217).
// - Dodge: double-tap a direction (on by default, Settings), or the dodge key (unbound by default: Alt is free look). A 0.45 s
//   step of 3 m with 0.2 s of invulnerability for 20 stamina; on the ground only, not while winded. It cancels the first
//   half of a melee wind-up and the last third of a recovery (#1216).
// RAYLIBGAME_COMBAT4=0: off (today: no parry, dodge, sprint cost or winded state; the wooden shield as before).
// This file is the logic only; combat.cpp and player.cpp use it.

#include <string>

bool Combat4_On();

namespace DefenseTune {
constexpr float kDodgeTime = 0.45f, kDodgeDistance = 3.0f, kDodgeInvuln = 0.2f, kDodgeStamina = 20.0f;
constexpr double kDoubleTap = 0.25;          // s between two presses of a direction
constexpr float kWindedTime = 1.5f, kWindedWalk = 0.7f;
constexpr float kSprintStamina = 8.0f;       // per second sprinting on foot (guess)
}

struct ShieldSpec {
    const char* name;
    float stops;        // share of a blocked hit's damage
    float arcDeg;       // half-angle of the guarded arc around the level view
    float parryWindow;  // s before a hit in which raising the shield parries it
    float raisedWalk;   // walk speed factor while it is raised
};
// The shield an item is (nullptr: not a shield): a buckler, a tower shield, else a shield is the round one.
const ShieldSpec* Shield_Of(const std::string& itemId);

// Was the block raised within the shield's parry window before `now`?
bool Parry_InWindow(double blockStartedAt, double now, const ShieldSpec& s);

// Double-tap detection: a press of direction dir (0 forward, 1 back, 2 left, 3 right) at `now`. Returns dir when it is the
// second press of that direction within kDoubleTap (the pair is used up), else -1.
struct DoubleTap {
    double last[4] = { -1e9, -1e9, -1e9, -1e9 };
    int press(int dir, double now);
};

#endif
