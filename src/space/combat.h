#ifndef COMBAT_H
#define COMBAT_H

// Combat (#165): hits on players resolved in one place, stamina, blocking and knockback.
// Combat_Damage is the only way a hit lands on a player. It takes everything it needs from the event and the
// target's session (no input, no camera, no globals of the local player), so a server can call it for any
// connected player (rolls are server-side, #309).

#include "dvec3.h"
#include "netid.h"
#include "damage.h" // (#1215) damage types and status effects

#include <cstdint>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

struct PlayerSession;
struct PlayerInput;
struct Player;
class Scene;

enum class DamageKind { Melee, Bite, Projectile, Other };

struct DamageEvent {
    float amount = 0.0f;          // health points before any block
    DamageKind kind = DamageKind::Melee;
    bool hasSource = false;       // a source position: a block must face it, knockback pushes away from it
    DVec3 source{};               // the attacker, world space (double: planet scale)
    float knockback = 0.0f;       // m/s away from the source along the ground, before a block
    NetId sourceId = 0;           // who hit (0: the world)
    // (#1215, combat step 1) What the hit is made of. The defaults (blunt, no effect) leave today's damage unchanged.
    DamageType type = DamageType::Blunt;
    StatusKind status = StatusKind::None; // an explicit effect (None: from the type, the heavy / back flags and a bite)
    float poise = 0.0f;           // poise damage (#1217)
    uint8_t faction = 0;          // the attacker's faction (friendly fire, #1223)
    bool parryable = true;        // a melee hit a parry can turn (#1218)
    bool heavy = false;           // a heavy attack (a slash heavy bleeds)
    bool fromBehind = false;      // landed on the target's back (a slash bleeds)
    bool fromFront = false;       // (#1227) came at the target's face (a front shield can take it)
};

struct DamageResult {
    float dealt = 0.0f;           // health taken
    bool blocked = false;         // a block took part of it
    bool guardBroken = false;     // the block ran out of stamina
    bool killed = false;          // this hit took the last health
    bool staggered = false;       // (#1217) it broke the player's poise
    bool parried = false;         // (#1218) a parry stopped it whole (the attacker staggers)
    bool dodged = false;          // (#1218) a dodge's invulnerability let it pass
};

DamageResult Combat_Damage(PlayerSession& target, const DamageEvent& e);
struct SurvivalState;
DamageResult Combat_DamageFor(SurvivalState& survival, Player& player, const DamageEvent& e); // (#337) any player's state (a server's remote players)

// One sim step of a player's combat state: the block (RMB with a shield in the off hand), stamina back at rest.
void Combat_Step(PlayerSession& s, float dt, const PlayerInput& in);
void Combat_CurePoison(PlayerSession& s); // (#1215) what bitterroot does (#1374: eating it cures through Survival_UseItemFor)
inline bool gDoubleTapDodge = true; // (#1218) double-tapping a direction dodges (Settings, on by default)

// Hit feel, the local view only (never the simulation): a short camera shake after applyToCamera, and whether the held
// tool's animation pauses (hit-stop). Reduce motion turns both off.
void Combat_ShakeCamera(Scene& world, const PlayerSession& s, bool reduceMotion);
bool Combat_HitStop(const PlayerSession& s, bool reduceMotion);

// Death losses by the world's difficulty (#165): peaceful, tools wear a tenth; normal, also a tenth of each resource
// stack; hard, everything goes into a grave at the death spot (a waypoint; walk up to it to take it all back).
void Combat_OnDeath(PlayerSession& s);
int Combat_GraveCount();
// (#828) Whether the waypoint `name` on `body` is a grave; `label` then reads as the map shows it ("bob's grave 2, open to
// all in 12 min" in multiplayer, its name alone in single player, where anyone may take it).
bool Combat_GraveMark(const std::string& body, const std::string& name, std::string& label);
// Drops a stack at the player's feet as a pile (#564), picked up by walking back to it after leaving; false: nowhere.
bool Combat_Drop(PlayerSession& s, const struct InvSlot& stack);
// (#865) A pile at `at` (body-local) on that body: a server dropping a remote player's stack (who: its name, the pile's owner;
// only their walk arms it). A joined client's
// Combat_Drop sends gDropRequest instead (the item, hundredths of a unit): the server makes the pile and the deltas.
bool Combat_DropFor(const std::string& body, const DVec3& at, const struct InvSlot& stack, const std::string& who); // at: body-local feet
inline std::function<void(const std::string& item, int64_t centi)> gDropRequest;
int Combat_DroppedPiles();
// (#1225) Cave scavengers steal dropped piles: the nearest within radius (body-local), and taking it whole.
bool Combat_PileNear(const std::string& body, Vector3 at, float radius, Vector3& pileAt);
bool Combat_TakePileNear(const std::string& body, Vector3 at, float radius, std::vector<struct InvSlot>& out);
void Combat_DrawGraves(Scene& world);
void Combat_Reset();                          // a new world: no graves
void Combat_Save(FILE* f);
bool Combat_LoadLine(const char* tok, FILE* f);

// ---- Graves in multiplayer (#804; the user's rules, 2026-11-18 on #337) ----
// death_drops (a server setting, --death-drops): by the world's difficulty (single player's rules, #165), keep everything,
// or always a grave.
enum class DeathDrops : uint8_t { Difficulty, Keep, Grave };
class Inventory; struct ToolWearState;
// One player's losses at death (Combat_OnDeath is this for the local player). owner: the player's name ("" in single player:
// anyone may loot). wear: their tool wear (nullptr: none kept, a server's remote player). True when a grave was made.
bool Combat_OnDeathFor(Inventory& inv, ToolWearState* wear, const Player& p, const std::string& owner, DeathDrops mode, bool local);
// One player walking up to graves takes what they may (Combat_Step does this for the local player). True when anything moved.
// body, me: the player's body and body-local feet (a server: from its netState, right between its steps, #856).
bool Combat_PickUpGravesFor(const std::string& body, const DVec3& me, Inventory& inv, const std::string& who, bool dead, bool local);
// Who may loot a grave: the owner and their team at once; anyone else 15 minutes after the death, never inside the owner's
// team claim. Times are wall-clock seconds (graves outlive server restarts). An owner "" (single player): anyone.
bool Combat_CanLoot(const std::string& owner, int64_t diedAt, const std::string& who, int64_t now, bool inOwnerClaim);
constexpr int64_t kGraveOpenToAllSeconds = 15 * 60;
constexpr int64_t kGraveDecaySeconds = 14 * 24 * 3600; // the claims clock (#343): 14 days without the owner's team online
// Teams and claims arrive with #343; until then the same name is the only team, and no grave lies in a claim.
inline std::function<bool(const std::string& a, const std::string& b)> gSameTeam;
inline std::function<bool(const std::string& owner, const std::string& body, const DVec3& local)> gInOwnerTeamClaim;
bool Combat_SameTeam(const std::string& a, const std::string& b);
// Decay (server): the names online now keep their team's graves fresh; a grave whose owner's team has been away 14 days is
// removed with what it holds. Returns the graves removed.
int Combat_UpdateGraveOwners(const std::vector<std::string>& online, int64_t now);
// Replication (#804, NetMsg 31): the graves (not dropped piles) as the server lists them, and a client's mirror of them
// (drawn; its own and its team's get a compass waypoint and a toast when new). Version: bumps on every change.
struct GraveInfo { uint32_t id = 0; std::string body; DVec3 local{}; std::string owner, name; uint8_t style = 0; int64_t diedAt = 0;
                   bool dropped = false; std::vector<std::string> items; }; // (#865) a dropped pile and up to 3 of its item ids
std::vector<GraveInfo> Combat_GraveList();
int Combat_TestGraveOffHand(size_t index); // (#1856 tests) the off hand's item index in the index-th grave or pile (-1 none)
void Combat_SetGraveMirror(const std::vector<GraveInfo>& list, const std::string& me);
uint32_t Combat_GraveVersion();
uint8_t Combat_GraveStyle(const std::string& body, const DVec3& local, uint32_t id); // 0 cross, 1 headstone, 2 monolith
void Combat_TestAddGrave(const std::string& body, const DVec3& local, uint8_t style, const std::string& owner); // --shot-graves

// --test-combat: called once a frame once the world is idle; true when finished (it prints PASS or FAIL).
bool Combat_TestFrame(Scene& world);

#endif
