#ifndef CREATURES_H
#define CREATURES_H

// Creature framework (#1223, Plan 12 step 1): what makes an animal a creature. Species keep their body, model, health
// and drops in animals.txt; resources/data/creatures.txt adds, per species, a faction, senses, behaviour, an attack,
// resistances, spawn rules and extra loot, and resources/data/factions.txt says which factions hunt, fear or ignore
// which. A species with a creatures.txt row is driven by the brain below (idle, wander, alert, chase, attack, retreat,
// flee, sleep); the rest stay plain animals.
//
// This file is the logic only (no world): data parsing, the faction table, senses (sight scaled by light, hearing from
// noise events) and the brain, a pure step function animals.cpp feeds with what the creature perceives.

#include "damage.h"

#include <raylib.h>
#include <cstdint>
#include <functional>
#include <string>
#include <utility>
#include <vector>

// ---- Factions ----
enum class Stance : uint8_t { Neutral, Hostile, Fears, Friendly };

struct FactionTable {
    std::vector<std::string> names;
    std::vector<Stance> m;                 // names.size()^2, row = the one deciding, column = the other
    int find(const std::string& name) const;
    int add(const std::string& name);      // the index (existing or new)
    Stance stance(int a, int b) const;     // a toward b (Neutral when unknown; a faction is Friendly to itself)
    void set(int a, int b, Stance s);
};
// Lines "faction | hostile: a, b | fears: c | friendly: d" ('#' comments). Factions named only as targets exist too.
bool Factions_Parse(const std::string& text, FactionTable& out, std::string* error = nullptr);
const FactionTable& Factions();            // resources/data/factions.txt (loaded once)
int Faction_Player();                      // the "player" faction's index

// ---- Species data ----
enum class Aggro : uint8_t { Passive, Defensive, Territorial, Predator };
enum class LightPref : uint8_t { None, Avoids, Drawn };
enum class Habitat : uint8_t { Surface, Night, Cave, Mouth, Dungeon };

struct CreatureSpec {
    std::string id, faction;
    int factionIdx = -1;
    // Senses: sight reach in full light (m), times darkSight..1 with the light at what it looks at; hearing reach (m) for a
    // noise of loudness 1.
    float sight = 25.0f, darkSight = 0.35f, hearing = 15.0f;
    LightPref light = LightPref::None;
    // Behaviour.
    Aggro aggro = Aggro::Passive;
    float aggroRange = 12.0f;  // m: a hostile target this near (and seen) is chased
    float chaseRange = 30.0f;  // m: further than this it gives up
    float angryTime = 20.0f;   // s a hit keeps a defensive or territorial creature after its attacker
    float retreatAt = 0.0f;    // health share under which it retreats (0: never)
    float threat = 1.0f;       // how dangerous it looks to others (fear and target choice)
    bool sleepsByDay = false;  // a night creature: rests in daylight
    bool fireScares = false;   // keeps 14 m from fires; a player by a fire within 10 m is safe from it
    // Its attack (the combat steps' timing: telegraphed wind-up, the strike, a cooldown).
    float damage = 8.0f, reach = 1.0f, windUp = 0.5f, cooldown = 1.2f, knockback = 6.0f, poise = 0.0f;
    DamageType type = DamageType::Blunt;
    bool bite = false;         // a bite (poisons through Status_FromHit)
    Resistances resist;
    // Spawning.
    Habitat habitat = Habitat::Surface;
    float depthMin = 0.0f, depthMax = 1e9f; // m under the surface (caves)
    float darkMax = 1.0f;      // spawns only where the light is at most this (0..1)
    int capPerChunk = 2;       // alive at once per 64 m chunk
    float respawn = 600.0f;    // s for a chunk to refill one, after it was cleared
    bool nests = false;        // larger caves hold a nest: destroyed, the chunk stops respawning
    std::vector<std::pair<std::string, int>> loot; // extra drops on top of animals.txt's
    // Sound names per event (sound_idle= ... sound_death=; empty: none). The audio side (#1246) maps names to sounds; this
    // framework only names them and raises Creature_Event.
    std::string sounds[5];
};

// "id | faction | key=value ..." (see creatures.txt). False with a reason on an unknown key or a bad value.
bool Creatures_ParseLine(const std::string& line, CreatureSpec& out, std::string* error = nullptr);
const CreatureSpec* Creature_Spec(const std::string& speciesId); // resources/data/creatures.txt; nullptr: a plain animal
const std::vector<CreatureSpec>& Creature_All();
bool Creatures_On(); // RAYLIBGAME_CREATURES=0: the old hostile/night animal code (A/B)

// Events for sound and effects (#1246 plugs in as data): idle now and then, alert (a noise or a sighting), the attack's
// telegraph (its wind-up), hurt, death. The hook gets the species, the event, where, and the row's sound name for it.
enum class CreatureEvent : uint8_t { Idle, Alert, Telegraph, Hurt, Death };
const char* CreatureEvent_Name(CreatureEvent e);
inline std::function<void(const std::string& species, CreatureEvent e, uint64_t body, Vector3 local, const std::string& sound)> gCreatureEvent;
void Creature_Event(const CreatureSpec& c, CreatureEvent e, uint64_t body, Vector3 local); // calls the hook (counts for tests)
int Creature_EventCount(CreatureEvent e);

// ---- Senses ----
// Light (0..1) at what is looked at: the sky's light there, or a carried light. Seen within sight x (darkSight..1).
float Sense_SightRange(const CreatureSpec& c, float lightAtTarget);
bool Sense_Sees(const CreatureSpec& c, float dist, float lightAtTarget, bool lineOfSight);
bool Sense_Hears(const CreatureSpec& c, float dist, float loudness);

// Noise events (body-local, kept 1.5 s): footsteps and sprinting, digging, combat, falls. Loudness 1 ~ a sprint; 2 a hit.
struct NoiseEvent { uint64_t body = 0; Vector3 pos{}; float loudness = 1.0f; double at = 0.0; uint64_t source = 0; };
namespace NoiseTune { constexpr double kKeep = 1.5; constexpr float kWalk = 0.4f, kSprint = 1.0f, kDig = 1.2f, kHit = 2.0f; }
void Noise_Emit(uint64_t body, Vector3 local, float loudness, double now, uint64_t source = 0);
const std::vector<NoiseEvent>& Noise_Recent(double now); // drops the old ones
void Noise_Clear();

// ---- The brain ----
enum class Mood : uint8_t { Idle, Wander, Alert, Chase, Attack, Retreat, Flee, Sleep };
const char* Mood_Name(Mood m);

// One thing it perceives (another creature or a player), as animals.cpp found it.
struct Percept {
    uint64_t id = 0;           // a creature's net id, or a player's (any stable id)
    Vector3 pos{};             // body-local
    float dist = 1e9f;
    int faction = -1;
    float threat = 1.0f;
    bool seen = false;         // within sight (light, line of sight)
    bool safe = false;         // a player by a fire it fears: untouchable
    bool vulnerable = true;    // a player inside a cockpit is not
};

struct BrainInput {
    double now = 0.0;
    float health01 = 1.0f;
    bool daylight = true;
    bool peaceful = false;     // no hostility at all (decided: peaceful has no hostile creatures)
    bool nearFire = false;     // a fire within 14 m (fireScares: it backs off)
    Vector3 fire{};
    Vector3 self{};
    std::vector<Percept> others;
    bool hurt = false; uint64_t hurtBy = 0; Vector3 hurtFrom{}; // hit this step, and by whom
    bool heardNoise = false; Vector3 noiseAt{};
};

struct BrainState {
    Mood mood = Mood::Idle;
    float timer = 0.0f;        // s left in idle / wander / alert / retreat / flee
    uint64_t target = 0;       // the one chased or attacked
    double angryUntil = -1.0;  // defensive / territorial: after a hit
    uint64_t angryAt = 0;
    double windUpAt = -1.0;    // an attack winding up since (-1 none)
    double lastStrike = -1e9;
    Vector3 goal{};            // alert: where the noise was; flee / retreat: what it runs from
};

struct BrainOut {
    Mood mood = Mood::Idle;
    bool move = false;         // toward goal (chase, alert) or away from it (flee, retreat)
    bool away = false;
    Vector3 goal{};
    float speed = 0.0f;        // 0..1 of its run (walk ~0.3)
    bool windUpStarted = false;// the telegraph (a sound, a toast)
    bool strike = false;       // the attack lands now (if the target is still in reach + 0.4 m)
    uint64_t target = 0;
};

BrainOut Brain_Step(BrainState& s, const CreatureSpec& c, const BrainInput& in, float dt);

#endif
