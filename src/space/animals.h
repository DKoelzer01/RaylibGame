#ifndef ANIMALS_H
#define ANIMALS_H

// Animals (plan 2, milestone 12): passive placeholder animals - coloured boxes with a head - from
// resources/data/animals.txt. Herds are placed from the world seed in 256 m tiles near the player
// (once seen, an animal is saved); only animals within 300 m are simulated.
//
// Behaviour: wander and graze, rest at night, flee when hit (wild ones also when the player comes
// close). They walk on the terrain surface, avoid water and steep ground, and treat build parts as
// walls, so a ring of blocks or beams is a pen. Look at one: I feeds it its food (it then follows
// you for 30 s - lead it into a pen; two fed adults near each other breed a young one that grows up
// over a day), O collects its product (wool with shears, milk with a bucket, eggs), LMB attacks
// (meat and hide when it falls).

#include <cstdint>
#include <functional>
#include <map>
#include <cstdio>
#include <string>
#include <utility>
#include <vector>
#include <raylib.h>

class Scene;
class Planetoid;

void Animals_Update(Scene& world, float dt, float worldDt, bool inputEnabled);
void Animals_Draw3D(Scene& world);
void Animals_DrawHUD();
struct PlayerInput;
void Animals_Aim(Scene& world, const PlayerInput& in); // this frame's aim and actions (after the camera is placed, before the terrain tool)
bool Animals_Aiming();                    // the crosshair is on an animal (the terrain tool holds off)
// (#1216) Melee in phases (melee.h): a swing under way (the terrain tool holds off too), its pose for the held tool (phase:
// 1 wind-up, 2 hit window, 3 recovery; f its progress; charge 0..1 the heavy's extra pull; back: the chained light), the soft
// lock-on (the view eases toward it; 0 none), hits landed (tests).
bool Animals_MeleeBusy();
bool Animals_MeleePose(int& phase, float& f, float& charge, bool& back);
void Animals_LockOnView(Scene& world, float dt);
uint64_t Animals_LockOn();
int Animals_MeleeHits();
const struct MeleeState& Animals_Melee();
void Animals_TestPin(int index, Vector3 bodyLocal, Vector3 faceLocal); // (#1216 tests) a dummy held in place
void Animals_Reset();
void Animals_Save(FILE* f);
bool Animals_LoadLine(const char* tok, FILE* f);

// Tests.
int Animals_Count();
void Animals_PathStats(int& plans, int& found, double& avgMs, double& worstMs); // route planning (#253): counts and cost
uint64_t Animals_NetId(int index);                              // its stable id (netid.h, #317)
int Animals_SimulatedCount();                                   // simulated in the last update: near an interest centre (#318)
int Animals_SimulatedNear(const Planetoid& body, Vector3 local, float range); // of those, within range of a body-local point
int Animals_Spawn(Scene& world, Planetoid& body, const std::string& species, Vector3 bodyLocal, int count, bool young = false, float spread = 18.0f);
// Every species as (id, display name), in animals.txt order: the Creative tab's spawner (#529).
std::vector<std::pair<std::string, std::string>> Animals_SpeciesList();
double Animals_BiteWindUp();            // s from a hostile's warning to its strike (#165)
int Animals_Nearest(Vector3 bodyLocal, float range);          // index or -1
Vector3 Animals_Position(int index);                           // body-local
bool Animals_Hit(int index, float damage);                     // true if it fell
// (#1215) A typed hit (damage.h): through the species' resistances, with its status effect; true if it fell. dealtOut: the damage.
bool Animals_HitEvent(int index, const struct DamageEvent& e, float* dealtOut = nullptr);
const struct StatusState* Animals_Status(int index);         // (#1215 tests) its status effects (nullptr: no such animal)
bool Animals_Feed(int index);                                  // with the player's food for it
int Animals_YoungCount();
double Animals_LastUpdateMs();
const std::string& Animals_Body(int index);                   // the body's name
bool Animals_Wild(int index);
bool Animals_Hostile(int index);
bool Animals_IsFeed(const std::string& item);                   // some animal eats it
void Animals_SetThreats(bool on);
void Animals_SetNaturalSpawns(bool on); // herds placed around the player as it moves (tests turn it off to measure their own animals)
void Animals_RefreshNow();                                     // tests: rescan campfires, allow a night spawn, recheck pens
int Animals_NightCount();                              // hostile behaviour and night spawns (bench runs: off)
int Animals_PennedCount();                                   // near the player, last checked (every 2 s)

// #334 multiplayer. The server runs the animals (around every player: interest centres, #318) and sends each client the
// ones near it; a client (gAnimalsReplicated) decides nothing - no AI, no spawns - and eases each copy toward the state it
// last heard, dropping any not heard of for 3 s. A client's swing at an animal is a request (gAnimalHitRequest); the server
// checks it, applies the damage, and credits the drops to that player.
inline bool gAnimalsReplicated = false;
inline std::function<void(uint64_t netId)> gAnimalHitRequest;
// A client's feed (action 1) and take-the-product (2: shear, milk, eggs) go to the server too (slice 2).
inline std::function<void(uint64_t netId, int action)> gAnimalUseRequest;
// The server's other players, for the animals to mind (slice 2): flee, follow their food, hunt, bite. The host lists them each
// frame; a bite on one goes to gAnimalDamageRemote (its client applies it: shield and health are its own until #337).
struct AnimalTarget { uint64_t body = 0; Vector3 local{}; uint32_t session = 0; bool vulnerable = true; std::string held; };
inline std::vector<AnimalTarget> gAnimalTargets;
struct DamageEvent;
class Inventory;
inline std::function<void(uint32_t session, const DamageEvent& e)> gAnimalDamageRemote;
// Feed (its food from `inv`) or take its product (into `inv`; the tool must be in `inv`) for a remote player; false with why.
bool Animals_UseNet(uint64_t netId, int action, Inventory& inv, std::string& why, std::map<std::string, double>& before); // before: the items it may change, as they were
struct AnimalNet {
    uint64_t netId = 0; uint8_t species = 0; uint64_t body = 0; // the body's NetId
    Vector3 pos{}; float heading = 0.0f;                         // body-local, radians
    uint8_t state = 0; float health = 0.0f, age = 1.0f;
    uint8_t flags = 0;                                           // 1 hostile, 2 angry, 4 winding up a bite
};
void Animals_CollectNet(Scene& world, std::vector<AnimalNet>& out);                       // the server: every animal now
void Animals_ApplyNet(const std::vector<AnimalNet>& in, const std::vector<uint64_t>& removed); // a client: states heard, and removals
int Animals_FindNet(uint64_t netId);                                                      // index, -1 none
bool Animals_HitNet(uint64_t netId, float damage, std::vector<std::pair<std::string, int>>& drops); // the server; true if it fell (drops filled)
float Animals_HealthOf(int index);

#endif
