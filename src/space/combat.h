#ifndef COMBAT_H
#define COMBAT_H

// Combat (#165): hits on players resolved in one place, stamina, blocking and knockback.
// Combat_Damage is the only way a hit lands on a player. It takes everything it needs from the event and the
// target's session (no input, no camera, no globals of the local player), so a server can call it for any
// connected player (rolls are server-side, #309).

#include "dvec3.h"
#include "netid.h"

#include <cstdio>

struct PlayerSession;
struct PlayerInput;
class Scene;

enum class DamageKind { Melee, Bite, Projectile, Other };

struct DamageEvent {
    float amount = 0.0f;          // health points before any block
    DamageKind kind = DamageKind::Melee;
    bool hasSource = false;       // a source position: a block must face it, knockback pushes away from it
    DVec3 source{};               // the attacker, world space (double: planet scale)
    float knockback = 0.0f;       // m/s away from the source along the ground, before a block
    NetId sourceId = 0;           // who hit (0: the world)
};

struct DamageResult {
    float dealt = 0.0f;           // health taken
    bool blocked = false;         // a block took part of it
    bool guardBroken = false;     // the block ran out of stamina
    bool killed = false;          // this hit took the last health
};

DamageResult Combat_Damage(PlayerSession& target, const DamageEvent& e);

// One sim step of a player's combat state: the block (RMB with a shield in the off hand), stamina back at rest.
void Combat_Step(PlayerSession& s, float dt, const PlayerInput& in);

// Hit feel, the local view only (never the simulation): a short camera shake after applyToCamera, and whether the held
// tool's animation pauses (hit-stop). Reduce motion turns both off.
void Combat_ShakeCamera(Scene& world, const PlayerSession& s, bool reduceMotion);
bool Combat_HitStop(const PlayerSession& s, bool reduceMotion);

// Death losses by the world's difficulty (#165): peaceful, tools wear a tenth; normal, also a tenth of each resource
// stack; hard, everything goes into a grave at the death spot (a waypoint; walk up to it to take it all back).
void Combat_OnDeath(PlayerSession& s);
int Combat_GraveCount();
void Combat_DrawGraves(Scene& world);
void Combat_Reset();                          // a new world: no graves
void Combat_Save(FILE* f);
bool Combat_LoadLine(const char* tok, FILE* f);

// --test-combat: called once a frame once the world is idle; true when finished (it prints PASS or FAIL).
bool Combat_TestFrame(Scene& world);

#endif
