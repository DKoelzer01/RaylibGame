#include "fluid.h"
#include "factory.h"
#include "animals.h"
#include "gamerng.h"
#include "netid.h"
#include "playerinput.h"
#include "instancebatch.h"
#include "audio.h"
#include "combat.h"
#include "melee.h"
#include "creatures.h" // (#1223)
#include "subtitles.h" // (#1635) the telegraph captions
#include "profiler.h" // (#1649) zones in the update
#include "jobsystem.h" // (#1650 slice 3) the probes on the pool
#include "footsteps.h" // (#1240) step noise for creature hearing
#include "avatar.h" // (#1225) humanoid creatures
#include "locallights.h"
#include "toolwear.h" // (#1226) an ooze corrodes tools
#include "claims.h"
#include "headless.h"
#include "net.h"
#include "netproto.h" // (#1343, #1364) creature shots and events to clients
#include "playersession.h"
#include "buildtool.h"
#include "gameui.h"
#include "hotbar.h"
#include "rocket.h"
#include "weather.h"
#include "inventory.h"
#include "player.h"
#include "plots.h"
#include "survival.h"
#include "terraintool.h"
#include "tools.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <atomic>
#include <functional>
#include <memory>
#include <mutex>
#include <thread>
#include <chrono>
#include <sstream>
#include <cmath>
#include <cstring>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <vector>

namespace {

// ---------------------------------------------------------------------------------------------
// Species (resources/data/animals.txt)
// ---------------------------------------------------------------------------------------------

struct Species {
    std::string id, name;
    Vector3 size{ 1, 1, 1 };                         // length, height, width (m)
    Color body = WHITE, head = GRAY;
    float walk = 1.0f, run = 4.0f, health = 10.0f;
    Resistances resist;                              // (#1215) per damage type and flat armour (species data later, #1223)
    std::vector<std::string> diet;
    std::vector<std::pair<std::string, int>> drops;
    std::string product, productTool;                // "-" none / no tool needed
    float regrow = 600.0f;                           // world seconds until the product is ready again
    bool wild = false;                               // flees the player at 25 m
    bool hostile = false;                            // charges the player within 12 m (by day too)
    bool night = false;                              // a night hunter: spawns near the player after dark, leaves at dawn
    int herdMin = 2, herdMax = 5;
    float latMin = 0.0f, latMax = 0.8f;
    std::vector<Biome> biomes;         // Plan 5 B: where herds live (empty: by latitude)
    bool biped = false, wool = false, horns = false, spikes = false, beak = false, tail = false, ears = false, longBody = false; // model (Plan 5 A)
    bool nest = false;                                // (#1223) a creature nest: a still mound, no legs or head
    bool hexa = false;                                // (#1224) six legs (the tunnel lurker)
    bool humanoid = false;                            // (#1225) drawn with the avatar rig (cave scavengers)
    std::string held;                                 // (#1227) a humanoid's item in hand (model held:<item>; default a torch)
    bool plankShield = false;                         // (#1227) a humanoid with a plank shield on its arm (model plankshield)
    bool blob = false;                                // (#1226) a mineral ooze: a wobbling heap, the colour of its ore
    bool turret = false;                              // (#1228) a ruin sentinel: a pedestal, a turning head and a barrel; its eye glows as it charges
};
std::vector<Species> species;
bool speciesLoaded = false;

std::vector<std::string> Split(const std::string& s, char sep) {
    std::vector<std::string> out;
    size_t start = 0, at;
    while ((at = s.find(sep, start)) != std::string::npos) { out.push_back(s.substr(start, at - start)); start = at + 1; }
    out.push_back(s.substr(start));
    for (auto& c : out) { const size_t a = c.find_first_not_of(" \t"), b = c.find_last_not_of(" \t\r\n"); c = a == std::string::npos ? "" : c.substr(a, b - a + 1); }
    return out;
}

void LoadSpecies() {
    if (speciesLoaded) return;
    speciesLoaded = true;
    FILE* f = fopen("resources/data/animals.txt", "r");
    if (!f) return;
    char buf[1024];
    while (fgets(buf, sizeof buf, f)) {
        if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
        const auto c = Split(buf, '|');
        if (c.size() < 13) continue;
        Species s;
        s.id = c[0]; s.name = c[1];
        sscanf(c[2].c_str(), "%f %f %f", &s.size.x, &s.size.y, &s.size.z);
        int r, g, b;
        if (sscanf(c[3].c_str(), "%d,%d,%d", &r, &g, &b) == 3) s.body = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 };
        if (sscanf(c[4].c_str(), "%d,%d,%d", &r, &g, &b) == 3) s.head = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 };
        sscanf(c[5].c_str(), "%f %f", &s.walk, &s.run);
        s.health = (float)atof(c[6].c_str());
        for (const auto& d : Split(c[7], ',')) if (!d.empty()) s.diet.push_back(d);
        for (const auto& d : Split(c[8], ',')) { char id[64] = {}; int n = 1; if (sscanf(d.c_str(), "%63s %d", id, &n) >= 1) s.drops.push_back({ id, n }); }
        char prod[64] = {}, tool[64] = {};
        if (sscanf(c[9].c_str(), "%63s %63s %f", prod, tool, &s.regrow) >= 1) { s.product = prod; s.productTool = tool; }
        s.wild = c[10] == "wild";
        s.hostile = c[10] == "hostile" || c[10] == "night";
        s.night = c[10] == "night";
        sscanf(c[11].c_str(), "%d %d", &s.herdMin, &s.herdMax);
        sscanf(c[12].c_str(), "%f %f", &s.latMin, &s.latMax);
        if (c.size() > 13) for (const auto& bn : Split(c[13], ',')) for (int b = 0; b < (int)Biome::Count; ++b) if (bn == BiomeName((Biome)b)) s.biomes.push_back((Biome)b);
        if (c.size() > 14) { // the model (Plan 5 A): quad or biped, then features
            std::istringstream ms(c[14]); std::string w;
            while (ms >> w) { if (w == "biped") s.biped = true; else if (w == "wool") s.wool = true; else if (w == "horns") s.horns = true; else if (w == "spikes") s.spikes = true; else if (w == "beak") s.beak = true; else if (w == "tail") s.tail = true; else if (w == "ears") s.ears = true; else if (w == "long") s.longBody = true; else if (w == "nest") s.nest = true; else if (w == "hexa") s.hexa = true; else if (w == "humanoid") s.humanoid = true; else if (w.rfind("held:", 0) == 0) s.held = w.substr(5); else if (w == "plankshield") s.plankShield = true; else if (w == "blob") s.blob = true; else if (w == "turret") s.turret = true; }
        }
        species.push_back(s);
    }
    fclose(f);
    // (#1635) A creature's telegraph sound is captioned with its own wind-up ("Lost miner raises its pick").
    for (const Species& s : species)
        if (const CreatureSpec* cs = Creature_Spec(s.id); cs && !cs->sayWindup.empty())
            Subtitles_SetCaption(cs->sounds[(int)CreatureEvent::Telegraph], s.name + " " + cs->sayWindup, true);
}
int SpeciesIndex(const std::string& id) { LoadSpecies(); for (size_t i = 0; i < species.size(); ++i) if (species[i].id == id) return (int)i; return -1; }

// ---------------------------------------------------------------------------------------------
// Animals
// ---------------------------------------------------------------------------------------------

enum class State : uint8_t { Graze, Walk, Flee, Follow, Rest, Hunt };

struct Animal {
    int species = 0;
    uint64_t netId = 0;          // stable, saved (netid.h, #317)
    DVec3 prevPos;               // pos and heading before the last simulation step: drawing blends them by SimAlpha (#315)
    float prevHeading = 0.0f;
    bool prevValid = false;
    int simFrame = -1;           // the update (frameNo) that last simulated it (#318 tests)
    std::string body;            // the body's name (saves); its pointer is looked up by id each frame
    DVec3 pos;                   // body-local metres, on the surface; doubles: a float has ~4 mm steps 60 km out, so it moved in visible jumps (#1698)
    float heading = 0.0f;        // radians from east toward north, in the local tangent frame
    float health = 10.0f;
    float age = 1.0f;            // 0 newborn .. 1 adult
    bool fed = false;
    double followUntil = -1.0;   // world time
    double lastBred = -1e9;
    double productReady = 0.0;   // world time the product (wool, milk, eggs) is ready
    State state = State::Graze;
    float timer = 0.0f, fleeTimer = 0.0f;
    Vector3 fleeFrom{};
    uint32_t rng = 1;
    double angryUntil = -1.0, lastBite = -1.0; // hostile: charging until (world time); the last bite
    double windUpAt = -1.0;      // an attack winding up since (simClock; #165), -1 none
    float windUpHeading = 0.0f;  // (#1227) its heading when the wind-up began: held until the strike
    Vector3 home{};              // (#1227) a cave creature's spawn floor (body-local): where it goes back to
    bool hasHome = false;
    StatusState status;          // (#1215) burn, chill, bleed, poison
    PoiseState poise;            // (#1217) worn down by hits; broken, it staggers (set up on its first hit)
    BrainState brain;            // (#1223) a creature's mind (creatures.txt rows only)
    bool hurtPending = false;    // (#1223) hit since its last update: by hurtBy (a player: kPlayerPercept | session), from hurtFrom
    uint64_t hurtBy = 0;
    Vector3 hurtFrom{};
    double heardPlayerAt = -1e9; // (#1224) sim time it last heard a player (a blind hunter tracks by this)
    bool clinging = false;       // (#1224) on a cave ceiling (pos), waiting over clingFloor
    bool dropOnSteps = false;    // (#1224) this drop was set off by prey walking under it (not hunting)
    float dropIn = -1.0f, clingWait = 0.0f; // the drop's telegraph left (s, -1 none); s it has heard prey not under it
    Vector3 clingFloor{};
    uint64_t group = 0;          // (#1225) its group (the first member's net id; 0: alone)
    double lastDecideAt = -1.0;  // (#1650) simClock of its last brain decision (AI tiers); -1: none yet
    double roofCheckAt = 0.0;    // (#1650) a ceiling clinger's next look for a ceiling (simClock)
    BrainOut lastBo;             // (#1650) that decision, reused until the next
    bool thrower = false;        // (#1225) the group's thrower: keeps its distance and throws
    double throwWindAt = -1.0, lastThrow = -1e9; // (#1225) a throw winding up since (simClock), the last throw
    std::vector<InvSlot> carried; // (#1225) what a thief took (dropped when it falls)
    std::string ore;             // (#1226) an ooze: the ore it ate (its colour, its drop); set from where it lives
    int generation = 0;          // (#1226) an ooze: 0 whole, 1-2 split off
    float daylight = 0.0f;                     // night hunters: seconds of daylight seen (they leave)
    // Plan 5 A: needs (0..1) that drive the choice between grazing, wandering and resting; the herd; a
    // chase target (predators); the gait for the walk animation.
    float hunger = 0.3f, fatigue = 0.2f;
    int chase = -1;                             // predators: the animal index being chased (this update)
    Vector3 herdCentre{}; int herdSize = 0;
    float gait = 0.0f, moveSpeed = 0.0f;        // leg cycle (radians) and the last speed (m/s)
    float stuckTime = 0.0f;                     // seconds without progress while trying to move
    float settleIn = 0.0f;                      // seconds to the next ground check (#223)
    float detourSide = 0.0f;                    // +1 / -1: the side it goes round an obstacle on (0: none)
    float detourHeading = 0.0f, detourTime = 0.0f;   // the detour being followed, and for how much longer (s)
    std::vector<Vector3> path;                  // a planned route (#253): body-local ground points, next first
    Vector3 pathGoal{};                         // where the route was planned to (replanned when the goal moves on)
    float pathAge = 0.0f, planWait = 0.0f;      // seconds on this route; no new plan before this runs out (s)
    // #334 a client's copy (replicated mode): the server's latest state of it, eased toward; when it was last heard of.
    Vector3 netPos{}; float netHeading = 0.0f; double netSeen = -1.0;
    // (#1108) where it was lately (sim time, body-local; the last kTrailSeconds): a server checks a joined player's hit against
    // where the animal was when that player saw it, not where it is by the time the hit arrives
    std::vector<std::pair<double, Vector3>> trail;
};
constexpr double kTrailSeconds = 0.3;
std::vector<Animal> animals;

// (#1650 slice 2) AI tiers: how often a creature's brain decides, by its distance to the nearest player (GUESS ranges and
// rates, from the #1650 design: engaged every update, near 5 Hz, far 1 Hz). RAYLIBGAME_AI_TIERS=0: every update (A/B).
namespace AiTier {
    constexpr float kEngagedRange = 15.0f, kNearRange = 60.0f; // m
    constexpr double kNearPeriod = 0.2, kFarPeriod = 1.0;      // s between decisions
    inline bool On() { static const bool on = !(getenv("RAYLIBGAME_AI_TIERS") && getenv("RAYLIBGAME_AI_TIERS")[0] == '0'); return on; }
}
long long gAiReused = 0; // decisions reused (tests)
// (#1650) Who is near whom: the creatures in 40 m cells (body-local), so perceiving others looks at the 27 cells around a
// creature instead of every creature (that scan was O(n^2) a frame). Rebuilt once an update, and again after a removal
// shifted the indices (gAiGridDirty).
bool gAiGridDirty = true;
namespace AiGrid {
    constexpr float kCell = 40.0f;
    std::unordered_map<uint64_t, std::vector<int>> cells;
    uint64_t Key(const std::string& body, int x, int y, int z) {
        return (((uint64_t)(uint32_t)(x & 0x1FFFFF)) | ((uint64_t)(uint32_t)(y & 0x1FFFFF) << 21) | ((uint64_t)(uint32_t)(z & 0x1FFFFF) << 42)) ^ (std::hash<std::string>()(body) * 0x9E3779B97F4A7C15ull);
    }
    void Build() {
        for (auto& kv : cells) kv.second.clear();
        for (size_t i = 0; i < animals.size(); ++i) {
            const Vector3 p = animals[i].pos.toFloat();
            cells[Key(animals[i].body, (int)floorf(p.x / kCell), (int)floorf(p.y / kCell), (int)floorf(p.z / kCell))].push_back((int)i);
        }
        if (cells.size() > 4096) { std::unordered_map<uint64_t, std::vector<int>> keep; for (auto& kv : cells) if (!kv.second.empty()) keep.emplace(kv.first, std::move(kv.second)); cells.swap(keep); }
        gAiGridDirty = false;
    }
    template <class F> void Near(const Animal& a, F f) {
        if (gAiGridDirty) Build();
        const int cx = (int)floorf(a.pos.x / kCell), cy = (int)floorf(a.pos.y / kCell), cz = (int)floorf(a.pos.z / kCell);
        for (int dz = -1; dz <= 1; ++dz) for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            auto it = cells.find(Key(a.body, cx + dx, cy + dy, cz + dz));
            if (it == cells.end()) continue;
            for (int j : it->second) if (j < (int)animals.size()) f(j);
        }
    }
}
constexpr uint64_t kPlayerPercept = kPlayerSource; // (#1223) a player's id among a creature's percepts (| its session)
int IndexOfNet(uint64_t id) { for (size_t i = 0; i < animals.size(); ++i) if (animals[i].netId == id) return (int)i; return -1; }
// (#1223) Its faction: its creatures.txt row's; a plain animal is livestock (calm) or wild.
int AnimalFaction(const Animal& a) {
    const Species& s = species[a.species];
    if (const CreatureSpec* c = Creature_Spec(s.id)) return c->factionIdx;
    return Factions().find(s.wild ? "wild" : "livestock");
}
struct CreatureHit { uint64_t target, attacker; Vector3 from; DamageEvent e; }; // (#1223) one creature's strike on another
std::vector<CreatureHit> creatureHits;
int shieldBlocks = 0; // (#1227) hits a front shield took (tests)
// (#1343, #1364) A server's creature shots and events since the host last took them (Animals_TakeFx), for its clients.
std::vector<MsgCreatureFx::Shot> fxShots;
std::vector<uint32_t> fxEnded;
std::vector<MsgCreatureFx::Event> fxEvents;
void FxEvent(const Animal& a, const CreatureSpec& c, CreatureEvent e, uint64_t bodyId, Vector3 pos) {
    Creature_Event(c, e, bodyId, pos); // here (sounds through gCreatureEvent)
    if (gNetRole != NetRole::Host || fxEvents.size() >= 256) return;
    if (Planetoid* b = Planetoid::find(bodyId)) fxEvents.push_back({ a.netId, (uint8_t)e, (uint8_t)a.species, b->netId, DVec3(pos) });
}
std::unordered_map<std::string, std::unordered_set<uint64_t>> spawnedTiles; // body -> herd tiles decided
int aimed = -1;
static bool hitLatch = false; // LMB went to an animal: the terrain tool waits for its release (#222)
double lastSpawnCheck = -1.0, lastHit = -1.0, lastPenCheck = -10.0, lastFireScan = -10.0, lastNightSpawn = -10.0;
double simClock = 0.0; // seconds of animal steps (the sum of Animals_Update's dt): the bite wind-up and cooldown (#165), so tests that step animals many times in a frame still see time pass
std::vector<Vector3> fires; // light sources near the players (body-local): placed lights, and (#560) players carrying one
// (#560, the user: dusk stalkers fear any light source, not only campfires) Placed: campfires, torches, lamp blocks, and fire
// machines while they burn. Carried: a torch, fibre torch or lantern in hand (or the local player's off hand).
// RAYLIBGAME_STALKER_ANY_LIGHT=0: campfires only, as before (A/B).
bool AnyLightScares() { static const bool on = !(getenv("RAYLIBGAME_STALKER_ANY_LIGHT") && atoi(getenv("RAYLIBGAME_STALKER_ANY_LIGHT")) == 0); return on; }
bool IsCarriedLight(const std::string& id) { return id == "torch" || id == "fibre_torch" || id == "lantern"; }
bool threatsOn = true;     // off in bench runs except --test-threats
bool naturalSpawns = true;  // herds placed around the player (SpawnHerds); tests that measure their own animals turn it off
bool NearFire(Vector3 p, float r, Vector3& at) { for (const Vector3& f : fires) if (Vector3Distance(f, p) < r) { at = f; return true; } return false; }
int pennedNow = 0;
double lastUpdateMs = 0.0;
int frameNo = 0;
int simulatedNow = 0; // animals within reach of an interest centre in the last update (#318)

float Rand01(uint32_t& s) { s = s * 1664525u + 1013904223u; return (s >> 8) * (1.0f / 16777216.0f); }

bool TangentFrame(Vector3 up, Vector3& east, Vector3& north) {
    const Vector3 c = Vector3CrossProduct({ 0, 1, 0 }, up);
    if (Vector3Length(c) < 1e-3f) { east = { 1, 0, 0 }; north = Vector3Normalize(Vector3CrossProduct(up, east)); return false; }
    east = Vector3Normalize(c);
    north = Vector3CrossProduct(up, east);
    return true;
}

Planetoid* BodyNamed(Scene& world, const std::string& name) {
    for (auto& o : world.objects) if (auto* p = dynamic_cast<Planetoid*>(o.get()); p && p->name == name) return p;
    return nullptr;
}
uint64_t BodyIdOf(const std::string& name) { for (Planetoid* p : Planetoid::All()) if (p->name == name) return p->id; return 0; } // (#1223) creature events
// Its facing (body-local, level): heading from east toward north.
Vector3 AnimalForward(const Animal& a) {
    Vector3 east, north;
    TangentFrame(Vector3Normalize(a.pos.toFloat()), east, north);
    return Vector3Add(Vector3Scale(east, cosf(a.heading)), Vector3Scale(north, sinf(a.heading)));
}

// (#1226) An ooze's ore: the deposit where it lives (iron ore where there is none), set the first time it is asked.
const std::string& OozeOre(Animal& a) {
    if (a.ore.empty()) {
        const char* k = nullptr;
        for (Planetoid* p : Planetoid::All()) if (p->name == a.body) { k = OreKindOf(*p, a.pos.toFloat()); break; }
        a.ore = k && *k ? k : "iron_ore";
    }
    return a.ore;
}

// Solid enough to block an animal (a wall, an overhang): ~0.4 m into rock. Thin near-iso films hover just over the
// ground in places (density within ~0.06 of iso); they must not stop a walk.
constexpr float kRockMargin = 0.1f;
constexpr double kBiteWindUp = 0.5; // s from the warning to the strike (#165): time to step back or raise a shield

// The surface point under a body-local direction, and whether an animal may stand there. The ground is the
// density field (#223: the height map ignored the 3D band, crags, tors, digs and builds, so animals walked
// through them): searched down from `above` metres over radius `ref` to `below` under it, in 0.25 m steps
// with the crossing interpolated. ref < 0 (spawning): around the height map, the whole 3D band. Fails where
// the start is inside rock (a wall or overhang) or nothing solid is in reach (a cliff or a hole).
bool Ground(const Planetoid& body, const PlanetField& pf, Vector3 dir, Vector3& at, float& elevation,
            float ref = -1.0f, float above = 0.0f, float below = 0.0f) {
    if (ref < 0.0f) { ref = pf.radius + pf.elevation(dir); above = 30.0f; below = 16.0f; }
    auto dens = [&](float r) { return body.densityAt(Vector3Scale(dir, r)); };
    float prevR = ref + above, prevD = dens(prevR);
    if (prevD >= ISO_LEVEL + kRockMargin) return false;
    bool found = false;
    for (float r = prevR - 0.25f; r >= ref - below - 1e-4f; r -= 0.25f) {
        const float d = dens(r);
        // Air above, solid here (a start inside a thin film looks on down). The air must be real: 0.3 m over the crossing
        // still open, not a hairline seam in rock where the density dips a hair under the iso level (#373: animals
        // were set into such seams, 'inside the ground' for seconds, where the 3D terrain's density hovers near it).
        if (d >= ISO_LEVEL && prevD < ISO_LEVEL && dens(r + (prevR - r) * 0.5f + 0.3f) < ISO_LEVEL) {
            const float t = std::clamp((ISO_LEVEL - prevD) / std::max(d - prevD, 1e-6f), 0.0f, 1.0f);
            at = Vector3Scale(dir, prevR + (r - prevR) * t);
            found = true;
            break;
        }
        prevR = r; prevD = d;
    }
    if (!found) return false;
    elevation = Vector3Length(at) - pf.radius;
    // (#1223) Under the surface (a cave floor, 3 m or more down): no sea or river there, whatever the height.
    if (pf.elevation(dir) - elevation > 3.0f) return true;
    if (elevation < 0.3f) return false; // the sea
    if (pf.inlandWaterNear(dir, 0.0f)) {
        const float w = pf.waterElevation(dir);
        if (w != PlanetField::kNoWater && w > elevation + 0.3f && Fluid::GeneratedWaterValid(body.id, Vector3Scale(dir, pf.radius + w - 0.2f)))
            return false; // a river or pond deeper than 0.3 m (not one the simulation drained, #360)
    }
    return true;
}

bool BlockedByParts(const Planetoid& body, Vector3 at, Vector3 up, float height) {
    if (!body.plots) return false;
    for (const auto& pl : body.plots->plots)
        for (float h = 0.3f; h < height + 0.3f; h += 0.5f)
            if (pl->partAt(pl->cellAt(pl->toPlot(Vector3Add(at, Vector3Scale(up, h))))) >= 0) return true;
    return false;
}

// Every 0.2 m along a step (a long step, far away or at a low frame rate, must not hop a wall).
bool BlockedAlong(const Planetoid& body, Vector3 from, Vector3 to, float height) {
    if (!body.plots || body.plots->plots.empty()) return false;
    const int n = std::max(1, (int)ceilf(Vector3Distance(from, to) / 0.2f));
    for (int k = 1; k <= n; ++k) {
        const Vector3 at = Vector3Lerp(from, to, (float)k / n);
        if (BlockedByParts(body, at, Vector3Normalize(at), height)) return true;
    }
    return false;
}

// Pathfinding (#253): A* on a 1 m grid in the tangent plane around the animal (+-kPathCells), each cell's ground
// found with the same probe and limits as a step (ledge up/down from the cell it is reached from, headroom, build
// parts). Cells are probed lazily. One search runs at a time and is resumed across updates: at most kPathFrameMs a
// frame and kPathTotalMs in all (on eroded terrain a probe costs more, and getting out of a pocket first leads away
// from the goal, so a single 3 ms slice found nothing). Without a full route it returns the part of the way that gets
// closest, if that is real progress.
constexpr int kPathCells = 40;          // grid half-width (m): 81 x 81 cells
constexpr int kPathMaxExpand = 3000;
constexpr double kPathFrameMs = 1.0, kPathTotalMs = 60.0; // (#1649: 3 ms slices ran to 4-5 ms and were the creature step's p99)
int gPathPlans = 0, gPathFound = 0;     // tests: plans made and routes found
double gPathMs = 0.0, gPathWorstMs = 0.0; // total ms over all searches; the worst single frame's slice

struct PathSearch {
    static constexpr int W = 2 * kPathCells + 1;
    bool active = false;
    int owner = -1;
    long long lastFrame = 0;
    const Planetoid* body = nullptr;
    const PlanetField* pf = nullptr;
    Vector3 from{}, east{}, north{};
    float height = 1.0f, reachUp = 0.35f, reachDown = 0.35f;
    int gx = 0, gy = 0, start = 0, goalI = 0, best = 0, expanded = 0;
    float bestH = 0.0f;
    double spentMs = 0.0;
    std::vector<float> cost, failR;
    std::vector<int> parent;
    std::vector<Vector3> ground;
    std::vector<uint8_t> state; // 0 unseen, 1 ground known, 2 probe failed (from failR), 3 closed
    struct Node { float f; int i; bool operator<(const Node& o) const { return f > o.f; } };
    std::priority_queue<Node> open;
    int idx(int x, int y) const { return (y + kPathCells) * W + (x + kPathCells); }
    float heur(int x, int y) const { return sqrtf((float)((x - gx) * (x - gx) + (y - gy) * (y - gy))); }
};
PathSearch gSearch;

void PathBegin(const Planetoid& body, const PlanetField& pf, Vector3 from, Vector3 goal, float height, float reachUp, float reachDown, int owner, long long frame) {
    PathSearch& S = gSearch;
    ++gPathPlans;
    S.active = true; S.owner = owner; S.lastFrame = frame;
    S.body = &body; S.pf = &pf; S.from = from; S.height = height; S.reachUp = reachUp; S.reachDown = reachDown;
    TangentFrame(Vector3Normalize(from), S.east, S.north);
    const Vector3 g = Vector3Subtract(goal, from);
    S.gx = std::clamp((int)lroundf(Vector3DotProduct(g, S.east)), -kPathCells, kPathCells);
    S.gy = std::clamp((int)lroundf(Vector3DotProduct(g, S.north)), -kPathCells, kPathCells);
    constexpr int N = PathSearch::W * PathSearch::W;
    S.cost.assign(N, 1e30f); S.parent.assign(N, -1); S.ground.assign(N, Vector3{}); S.state.assign(N, 0); S.failR.assign(N, 0.0f);
    S.open = {};
    S.start = S.idx(0, 0); S.goalI = S.idx(S.gx, S.gy);
    S.cost[S.start] = 0.0f; S.ground[S.start] = from; S.state[S.start] = 1;
    S.open.push({ S.heur(0, 0), S.start });
    S.best = S.start; S.bestH = S.heur(0, 0); S.expanded = 0; S.spentMs = 0.0;
}

// One slice of the search: false while it goes on; true when finished (`out`: the route, empty if none worth taking).
bool PathStep(std::vector<Vector3>& out, long long frame) {
    PathSearch& S = gSearch;
    const auto t0 = std::chrono::steady_clock::now();
    S.lastFrame = frame;
    const Planetoid& body = *S.body;
    const PlanetField& pf = *S.pf;
    bool done = false;
    while (true) {
        if (S.open.empty() || S.expanded >= kPathMaxExpand) { done = true; break; }
        if ((S.expanded & 3) == 3 && std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() > kPathFrameMs) break; // (#1649) checked every 4: an expansion probes 8 cells
        const PathSearch::Node n = S.open.top(); S.open.pop();
        if (S.state[n.i] == 3) continue;
        S.state[n.i] = 3;
        ++S.expanded;
        const int x = n.i % PathSearch::W - kPathCells, y = n.i / PathSearch::W - kPathCells;
        const float h = S.heur(x, y);
        if (h < S.bestH) { S.bestH = h; S.best = n.i; }
        if (n.i == S.goalI) { done = true; break; }
        const Vector3 at = S.ground[n.i];
        const float r = Vector3Length(at);
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
                if (!dx && !dy) continue;
                const int nx = x + dx, ny = y + dy;
                if (nx < -kPathCells || nx > kPathCells || ny < -kPathCells || ny > kPathCells) continue;
                const int j = S.idx(nx, ny);
                if (S.state[j] == 3) continue;
                Vector3 gj; float e;
                // A cell's ground is probed once and reused from other neighbours within their ledge reach; a failed
                // probe is retried only from a height over 0.8 m different (a cell under a cliff may be reached from below).
                if (S.state[j] == 1) {
                    gj = S.ground[j];
                    const float dz = Vector3Length(gj) - r;
                    if (dz > S.reachUp || dz < -S.reachDown || BlockedAlong(body, at, gj, S.height)) continue;
                } else {
                    if (S.state[j] == 2 && fabsf(S.failR[j] - r) < 0.8f) continue;
                    const Vector3 dir = Vector3Normalize(Vector3Add(S.from, Vector3Add(Vector3Scale(S.east, (float)nx), Vector3Scale(S.north, (float)ny))));
                    if (!Ground(body, pf, dir, gj, e, r, S.reachUp, S.reachDown) || body.densityAt(Vector3Add(gj, Vector3Scale(dir, S.height * 0.8f))) >= ISO_LEVEL + kRockMargin ||
                        BlockedAlong(body, at, gj, S.height)) { S.state[j] = 2; S.failR[j] = r; continue; }
                    S.ground[j] = gj; S.state[j] = 1;
                }
                const float c = S.cost[n.i] + Vector3Distance(at, gj) + 0.5f * fabsf(Vector3Length(gj) - r);
                if (c < S.cost[j]) { S.cost[j] = c; S.parent[j] = n.i; S.open.push({ c + S.heur(nx, ny), j }); }
            }
    }
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    S.spentMs += ms; gPathMs += ms; gPathWorstMs = std::max(gPathWorstMs, ms);
    if (!done && S.spentMs < kPathTotalMs) return false;
    out.clear();
    const bool reached = S.best == S.goalI || S.bestH < 1.5f;
    // Partial routes only when they get at least 4 m closer: otherwise the local detour does as well.
    if (S.best != S.start && (reached || S.heur(0, 0) - S.bestH > 4.0f)) {
        for (int i = S.best; i != S.start && i >= 0; i = S.parent[i]) out.push_back(S.ground[i]);
        std::reverse(out.begin(), out.end());
    }
    if (reached) ++gPathFound;
    S.active = false;
    return true;
}

// An animal at `pos` is penned when build parts stand within 30 m in all eight directions, from 1.7 m under its feet
// to 2.3 m over them: on sloping ground a fence 8 m away stands a metre or more above or below the animal (#477: the
// old 0.3-1.3 m band at its own height missed fences on eroded terrain).
bool IsPenned(const Planetoid& body, Vector3 pos) {
    if (!body.plots || body.plots->plots.empty()) return false;
    const Vector3 up = Vector3Normalize(pos);
    Vector3 east, north;
    TangentFrame(up, east, north);
    const float r = Vector3Length(pos);
    for (int k = 0; k < 8; ++k) {
        const float ang = k * PI / 4.0f;
        const Vector3 dirT = Vector3Add(Vector3Scale(east, cosf(ang)), Vector3Scale(north, sinf(ang)));
        bool hit = false;
        for (float s = 0.25f; s <= 30.0f && !hit; s += 0.25f) {
            const Vector3 at = Vector3Scale(Vector3Normalize(Vector3Add(pos, Vector3Scale(dirT, s))), r - 2.0f);
            hit = BlockedByParts(body, at, Vector3Normalize(at), 4.0f);
        }
        if (!hit) return false;
    }
    return true;
}

constexpr float kTileSize = 256.0f;
uint64_t TileKey(Vector3 dir, float radius, Vector3& centreDir) {
    static const Vector3 kN[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
    static const Vector3 kU[6] = { { 0, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 0 } };
    static const Vector3 kV[6] = { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 } };
    int face = 0; float best = -2.0f;
    for (int f = 0; f < 6; ++f) { const float d = Vector3DotProduct(dir, kN[f]); if (d > best) { best = d; face = f; } }
    const int n = std::max(1, (int)ceilf(0.5f * PI * radius / kTileSize));
    const float da = 0.5f * PI / n;
    const float a = atanf(Vector3DotProduct(dir, kU[face]) / best), b = atanf(Vector3DotProduct(dir, kV[face]) / best);
    const int i = std::clamp((int)floorf((a + 0.25f * PI) / da), 0, n - 1), j = std::clamp((int)floorf((b + 0.25f * PI) / da), 0, n - 1);
    const float ca = -0.25f * PI + (i + 0.5f) * da, cb = -0.25f * PI + (j + 0.5f) * da;
    centreDir = Vector3Normalize(Vector3Add(kN[face], Vector3Add(Vector3Scale(kU[face], tanf(ca)), Vector3Scale(kV[face], tanf(cb)))));
    return ((uint64_t)face << 48) | ((uint64_t)i << 24) | (uint64_t)j;
}
uint32_t Hash(uint64_t k, uint32_t seed) {
    uint32_t h = (uint32_t)(k ^ (k >> 32)) ^ seed * 2654435761u;
    h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
    return h;
}

// New herds in tiles near the player the first time they come within 300 m (planets with life).
Vector3 StartDir(const Planetoid& body, const PlanetField& pf) {
    static std::unordered_map<std::string, Vector3> cache;
    auto it = cache.find(body.name);
    if (it != cache.end()) return it->second;
    return cache[body.name] = pf.findSpawnDir(body.seed);
}

void SpawnHerds(Scene& world, Planetoid& body, const PlanetField& pf, Vector3 up) {
    if (!pf.foliage || species.empty()) return;
    auto& done = spawnedTiles[body.name];
    Vector3 east, north;
    TangentFrame(up, east, north);
    for (float a = -300.0f; a <= 300.0f; a += 128.0f)
        for (float b = -300.0f; b <= 300.0f; b += 128.0f) {
            const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(up, pf.radius), Vector3Add(Vector3Scale(east, a), Vector3Scale(north, b))));
            Vector3 centre;
            const uint64_t key = TileKey(d, pf.radius, centre);
            if (!done.insert(key).second) continue;
            uint32_t h = Hash(key, body.seed ^ 0xa11a7u);
            if (Rand01(h) > 0.35f) continue; // about a third of the tiles hold a herd
            const float lat = fabsf(centre.y);
            const Biome tileBiome = pf.biomeAt(centre, pf.elevation(centre, 16.0f)); // Plan 5 B: herds by biome
            auto livesHere = [&](const Species& sp) {
                if (sp.biomes.empty()) return lat >= sp.latMin && lat <= sp.latMax;
                return std::find(sp.biomes.begin(), sp.biomes.end(), tileBiome) != sp.biomes.end();
            };
            std::vector<int> fit;
            for (size_t i = 0; i < species.size(); ++i)
                if (livesHere(species[i]) && species[i].herdMax > 0 && !species[i].night && !(species[i].hostile && worldDifficulty == 0) &&
                    (!body.profile || body.profile->allowsCreature(species[i].id))) // the world file's creatures (#1274)
                    fit.push_back((int)i);
            if (fit.empty()) continue;
            const Species& s = species[fit[(size_t)(Rand01(h) * fit.size()) % fit.size()]];
            if (s.hostile && Vector3Distance(Vector3Scale(centre, pf.radius), Vector3Scale(StartDir(body, pf), pf.radius)) < 500.0f) continue; // a safe start
            Vector3 at; float e;
            if (!Ground(body, pf, centre, at, e) || e > 900.0f) continue; // no herds at sea or above the tree line
            const int n = s.herdMin + (int)(Rand01(h) * (s.herdMax - s.herdMin + 1));
            Animals_Spawn(world, body, s.id, at, n);
        }
}

// toGround (#1223): killed by another creature, its drops lie where it fell (a pile anyone picks up), no quest counts.
void Kill(int i, std::vector<std::pair<std::string, int>>* dropsTo = nullptr, bool toGround = false) { // dropsTo (#334): a server's, for the hitting player
    const Animal& a = animals[i];
    const Species& s = species[a.species];
    const CreatureSpec* cs = Creature_Spec(s.id);
    if (cs) FxEvent(a, *cs, CreatureEvent::Death, BodyIdOf(a.body), a.pos.toFloat());
    if (!dropsTo && !toGround) QuestBook_Count("slain"); // (a server's kill isn't its own player's quest)
    if (!dropsTo && !toGround && s.night) QuestBook_Count("stalker");
    std::vector<std::pair<std::string, int>> all = s.drops;
    if (cs) all.insert(all.end(), cs->loot.begin(), cs->loot.end()); // (#1223) the creature's extra loot
    if (cs) { uint32_t r = a.rng ^ 0x9e3779b9u; for (const CreatureSpec::Rare& rr : cs->rare) if (Rand01(r) < rr.chance) all.push_back({ rr.item, rr.n }); } // (#1227) sometimes
    if (cs && cs->split) all.push_back({ OozeOre(animals[i]), a.generation == 0 ? 3 : 1 }); // (#1226) the ore an ooze ate
    std::string got;
    for (const auto& d : all) {
        const int n = std::max(1, (int)lroundf(d.second * (0.4f + 0.6f * a.age)));
        if (!itemDB.find(d.first)) continue;
        if (toGround) { InvSlot st; st.id = d.first; st.centi = (int64_t)n * 100; Combat_DropFor(a.body, a.pos, st, ""); continue; }
        if (dropsTo) { dropsTo->push_back({ d.first, n }); continue; }
        inventory.add(d.first, n); got += TextFormat("  +%d %s", n, d.first.c_str());
    }
    for (const InvSlot& st : a.carried) { // (#1225) what a thief took comes back
        const int n = std::max(1, (int)((st.centi + 50) / 100));
        if (toGround) { Combat_DropFor(a.body, a.pos, st, ""); continue; }
        if (dropsTo) { dropsTo->push_back({ st.id, n }); continue; }
        inventory.add(st.id, n); got += TextFormat("  +%d %s", n, st.id.c_str());
    }
    if (!dropsTo && !toGround) GameUI_Toast(TextFormat("The %s fell.%s", s.name.c_str(), got.c_str()));
    animals.erase(animals.begin() + i); gAiGridDirty = true; // (#1650) indices shifted
    if (aimed == i) aimed = -1;
}

// (#1223) The resistances a hit meets: its creatures.txt row's, else the species'.
const Resistances& ResistOf(const Animal& a) {
    if (const CreatureSpec* cs = Creature_Spec(species[a.species].id)) return cs->resist;
    return species[a.species].resist;
}

// A typed hit's damage on an animal (#1215 types and status, #1217 poise): what to take off its health.
float TypedAmount(Animal& a, const DamageEvent& e) {
    float amount = e.amount;
    // (#1227) A front shield (the lost miner's plank shield) takes most of a non-blunt hit that comes at its face.
    if (const CreatureSpec* sh = Creature_Spec(species[a.species].id); sh && sh->frontShield > 0.0f && e.fromFront && e.type != DamageType::Blunt) {
        amount *= 1.0f - sh->frontShield;
        ++shieldBlocks;
    }
    // (#1217) Poise: the hit wears it (its poise damage; none given: its amount); broken, the creature staggers (its attack
    // cancelled), and the first hit while it staggers does x1.5. Small creatures stagger on any heavy hit.
    if (Combat3_On()) {
        const Species& sp = species[a.species];
        if (a.poise.max <= 0.0f) Poise_Init(a.poise, Poise_CreatureMax(sp.health));
        const bool small = std::max(sp.size.x, std::max(sp.size.y, sp.size.z)) * (0.6f + 0.4f * a.age) <= PoiseTune::kSmallSize;
        const PoiseHit ph = Poise_Hit(a.poise, e.poise > 0.0f ? e.poise : e.amount, e.heavy, small);
        amount *= ph.damageScale;
        if (ph.staggered) { a.windUpAt = -1.0; a.brain.windUpAt = -1.0; }
    }
    if (Combat1_On()) {
        amount = Damage_Apply(amount, e.type, ResistOf(a));
        Status_Apply(a.status, Status_FromHit(e.type, e.status, e.heavy, e.fromBehind, e.kind == DamageKind::Bite));
        if (const CreatureSpec* im = Creature_Spec(species[a.species].id)) { // (#1228) a construct neither bleeds nor is poisoned
            if (im->immuneBleed) { a.status.bleed = 0.0f; a.status.bleedStacks = 0; }
            if (im->immunePoison) a.status.poison = 0.0f;
        }
    }
    return amount;
}


// (#1226) An ooze cut (a slashing or piercing hit) splits in two: the half it loses becomes a second ooze beside it, of
// the same ore (two generations at most). Blunt hits don't split it, nor does anything while it burns.
int oozeSplits = 0;
void MaybeSplit(int i, DamageType type) {
    Animal& a = animals[i];
    const CreatureSpec* cs = Creature_Spec(species[a.species].id);
    if (!cs || !cs->split || (type != DamageType::Slash && type != DamageType::Pierce) || a.status.burn > 0.0f || a.generation >= 2 || a.health <= 6.0f) return;
    OozeOre(a);
    Animal b = a;
    b.health = a.health * 0.5f;
    a.health -= b.health;
    b.generation = ++a.generation;
    const Vector3 up = Vector3Normalize(a.pos.toFloat());
    const Vector3 side = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    b.pos = b.prevPos = a.pos + DVec3(side) * 0.9;
    b.netId = NetId_Allocate(NetKind::Animal);
    b.brain = BrainState{}; b.poise = PoiseState{}; b.hurtPending = false; b.carried.clear();
    animals.push_back(b);
    ++oozeSplits;
}

// (#1223) One creature's strike on another: typed damage; the struck one learns who hit it (a creature turns on it, a plain
// animal runs); a kill leaves its loot on the ground.
void HitByCreature(int i, const DamageEvent& e, uint64_t attacker, Vector3 from) {
    Animal& a = animals[i];
    DamageEvent ef = e; // (#1227) at its face or not, from where the attacker stood
    { const CreatureSpec* fs = Creature_Spec(species[a.species].id);
      ef.fromFront = Vector3DotProduct(AnimalForward(a), Vector3Normalize(Vector3Subtract(from, a.pos.toFloat()))) > cosf((fs ? fs->shieldArcDeg : 60.0f) * DEG2RAD); }
    a.health -= TypedAmount(a, ef);
    if (const CreatureSpec* cs = Creature_Spec(species[a.species].id)) {
        a.hurtPending = true; a.hurtBy = attacker; a.hurtFrom = from;
        FxEvent(a, *cs, CreatureEvent::Hurt, BodyIdOf(a.body), a.pos.toFloat());
    } else { a.fleeTimer = 5.0f; a.fleeFrom = from; a.state = State::Flee; }
    if (a.health <= 0.0f) Kill(i, nullptr, true);
    else MaybeSplit(i, e.type); // (#1226)
}

} // namespace

// ---------------------------------------------------------------------------------------------

void Animals_Reset() { Animals_ResetCaveSpawns(); animals.clear(); spawnedTiles.clear(); aimed = -1; gSearch.active = false; gSearch.body = nullptr; } // (#902: the search held the old world's body)
int Animals_Count() { return (int)animals.size(); }
Vector3 Animals_Position(int i) { return i >= 0 && i < (int)animals.size() ? animals[i].pos.toFloat() : Vector3{ 0, 0, 0 }; }
double Animals_LastUpdateMs() { return lastUpdateMs; }
int Animals_PennedCount() { return pennedNow; }
void Animals_SetThreats(bool on) { threatsOn = on; }
void Animals_SetNaturalSpawns(bool on) { naturalSpawns = on; }
void Animals_RefreshNow() { lastFireScan = -100.0; lastNightSpawn = -100.0; lastPenCheck = -100.0; }
int Animals_NightCount() { int n = 0; for (const Animal& a : animals) n += species[a.species].night ? 1 : 0; return n; }
bool Animals_IsFeed(const std::string& id) { LoadSpecies(); for (const Species& s : species) for (const std::string& d : s.diet) if (d == id) return true; return false; }
const std::string& Animals_Body(int i) { static const std::string none; return i >= 0 && i < (int)animals.size() ? animals[i].body : none; }
bool Animals_Wild(int i) { return i >= 0 && i < (int)animals.size() && species[animals[i].species].wild; }
bool Animals_Hostile(int i) { return i >= 0 && i < (int)animals.size() && species[animals[i].species].hostile; }
int Animals_AimedIndex() { return aimed >= 0 && aimed < (int)animals.size() ? aimed : -1; }
std::pair<std::string, std::string> Animals_SpeciesOf(int i) {
    LoadSpecies();
    if (i < 0 || i >= (int)animals.size() || animals[i].species < 0 || animals[i].species >= (int)species.size()) return {};
    const Species& s = species[animals[i].species];
    return { s.id, s.name };
}
bool Animals_Aiming() { return aimed >= 0 || hitLatch || Animals_MeleeBusy(); }
int Animals_YoungCount() { int n = 0; for (const Animal& a : animals) n += a.age < 1.0f ? 1 : 0; return n; }

double Animals_BiteWindUp() { return kBiteWindUp; }

std::vector<std::pair<std::string, std::string>> Animals_SpeciesList() {
    LoadSpecies();
    std::vector<std::pair<std::string, std::string>> out;
    for (const Species& s : species) out.push_back({ s.id, s.name });
    return out;
}

int Animals_Spawn(Scene& world, Planetoid& body, const std::string& id, Vector3 local, int count, bool young, float spread) {
    (void)world;
    const int si = SpeciesIndex(id);
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (si < 0 || !pf) return 0;
    const Vector3 up = Vector3Normalize(local);
    Vector3 east, north;
    TangentFrame(up, east, north);
    uint32_t h = Hash((uint64_t)(local.x * 7.0f) ^ ((uint64_t)(local.z * 13.0f) << 20), (uint32_t)animals.size() + 17u);
    int made = 0;
    for (int k = 0; k < count; ++k) {
        const float r = spread >= 18.0f ? 4.0f + 14.0f * Rand01(h) : spread * Rand01(h), t = 2.0f * PI * Rand01(h);
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(up, pf->radius), Vector3Add(Vector3Scale(east, r * cosf(t)), Vector3Scale(north, r * sinf(t)))));
        Vector3 at; float e;
        if (!Ground(body, *pf, d, at, e)) continue;
        Animal a;
        a.species = si;
        a.body = body.name;
        a.pos = DVec3(at);
        a.heading = 2.0f * PI * Rand01(h);
        a.health = species[si].health;
        a.age = young ? 0.0f : 1.0f;
        a.rng = h ^ (uint32_t)(animals.size() * 7919u);
        a.timer = 1.0f + 4.0f * Rand01(h);
        a.settleIn = Rand01(h); // first ground check spread over a second: a herd spawned together doesn't all probe in one update (#476)
        a.productReady = WorldTime();
        a.netId = NetId_Allocate(NetKind::Animal);
        animals.push_back(a);
        ++made;
    }
    return made;
}

int Animals_Nearest(Vector3 local, float range) {
    int best = -1;
    float bestD = range;
    for (size_t i = 0; i < animals.size(); ++i) {
        const float d = Vector3Distance(animals[i].pos.toFloat(), local);
        if (d < bestD) { bestD = d; best = (int)i; }
    }
    return best;
}

bool Animals_Hit(int i, float damage) {
    if (i < 0 || i >= (int)animals.size()) return false;
    Animal& a = animals[i];
    a.health -= damage;
    if (const CreatureSpec* cs = Creature_Spec(species[a.species].id)) {
        // (#1223) A creature: its brain hears who hit it (the local player) and decides (turn on them, retreat, run).
        a.hurtPending = true; a.hurtBy = kPlayerPercept; a.hurtFrom = a.pos.toFloat();
        FxEvent(a, *cs, CreatureEvent::Hurt, BodyIdOf(a.body), a.pos.toFloat());
    } else {
    // Hostile ones fight back until badly hurt; the rest run.
    a.fleeTimer = !species[a.species].hostile || a.health < 0.35f * species[a.species].health ? 5.0f : 0.0f;
    if (species[a.species].hostile && a.fleeTimer <= 0.0f) a.angryUntil = WorldTime() + 20.0;
    a.fleeFrom = a.pos.toFloat(); // refined to the player's position in the update
    a.state = State::Flee;
    }
    Audio_Play(Sfx::Dig, 0.5f, 1.6f);
    if (a.health <= 0.0f) { Kill(i); return true; }
    return false;
}

// (#1215) A typed hit: through the species' resistances and armour (the 10% floor), its status effect, then the hit itself.
// RAYLIBGAME_COMBAT1=0: the amount as it is (today).
bool Animals_HitEvent(int i, const DamageEvent& e, float* dealtOut) {
    if (i < 0 || i >= (int)animals.size()) return false;
    const float amount = TypedAmount(animals[i], e); // (#1215, #1217; #1223: a creature's own resistances)
    if (dealtOut) *dealtOut = amount;
    const bool fell = Animals_Hit(i, amount);
    if (!fell) MaybeSplit(i, e.type); // (#1226) an ooze cut in two
    return fell;
}

const StatusState* Animals_Status(int i) { return i >= 0 && i < (int)animals.size() ? &animals[i].status : nullptr; } // (#1215 tests)
const PoiseState* Animals_Poise(int i) { return i >= 0 && i < (int)animals.size() ? &animals[i].poise : nullptr; } // (#1217 tests)

bool Animals_Feed(int i) {
    if (i < 0 || i >= (int)animals.size()) return false;
    Animal& a = animals[i];
    const Species& s = species[a.species];
    for (const std::string& food : s.diet)
        if (inventory.count(food) > 0 && inventory.remove(food, 1)) {
            a.fed = true;
            a.followUntil = WorldTime() + 30.0;
            QuestBook_Count("feed");
            a.fleeTimer = 0.0f;
            Audio_Play(Sfx::Pickup, 0.4f, 1.3f);
            return true;
        }
    std::string want;
    for (const std::string& f : s.diet) want += (want.empty() ? "" : " or ") + f;
    GameUI_Toast(TextFormat("The %s eats %s", s.name.c_str(), want.c_str()));
    return false;
}

// ---- (#1225) Thrown projectiles and humanoid creatures ----
// A creature's throw (a rock or spear) flies on a ballistic arc under the planet's pull (~9.8 m/s^2 along the local down),
// hits the first player whose body it passes within 0.45 m of (through Combat_Damage: a raised shield facing it blocks; a
// dodge passes through; never parried), and stops in rock. Gone after 6 s.
namespace {
struct Projectile { std::string body; Vector3 pos, vel; float damage; uint64_t owner; double born; bool spear; bool dart = false; uint32_t id = 0; }; // (#1343) id: a server's shots carry it to clients // (#1228) dart: a sentinel's bolt (pierce, fast, flat)
std::vector<Projectile> projectiles;
int projectileHits = 0, projectilesThrown = 0;

void ThrowAt(const std::string& body, Vector3 from, Vector3 target, float damage, uint64_t owner, bool spear, float speed = 16.0f, bool dart = false) {
    Vector3 vel;
    Throw_Aim(from, target, Vector3Normalize(from), speed, 9.81f, vel);
    static uint32_t nextShot = 1;
    const uint32_t id = nextShot++;
    projectiles.push_back({ body, from, vel, damage, owner, simClock, spear, dart, id });
    if (gNetRole == NetRole::Host && fxShots.size() < 128) // (#1343) its clients fly it too
        for (Planetoid* b : Planetoid::All()) if (b->name == body) { fxShots.push_back({ id, b->netId, DVec3(from), vel, (uint8_t)(dart ? 2 : spear ? 1 : 0), owner }); break; }
    ++projectilesThrown;
}

void StepProjectiles(Scene& world, float dt) {
    if (projectiles.empty() || dt <= 0.0f) return;
    Planetoid* here = player.nearest;
    const Vector3 feet = here ? here->toLocalOf(player.position, world.renderOrigin) : Vector3{};
    for (size_t k = 0; k < projectiles.size();) {
        Projectile& p = projectiles[k];
        Planetoid* body = BodyNamed(world, p.body);
        bool gone = !body || simClock - p.born > 6.0;
        if (!gone) {
            const Vector3 up = Vector3Normalize(p.pos);
            const Vector3 before = p.pos;
            p.vel = Vector3Add(p.vel, Vector3Scale(up, -9.81f * dt));
            p.pos = Vector3Add(p.pos, Vector3Scale(p.vel, dt));
            if (body->densityAt(p.pos) >= ISO_LEVEL) gone = true; // into rock or the ground
            // The local player: their body from the feet up 1.7 m.
            if (!gone && here == body && !survival.dead && !player.flyMode && !gAnimalsReplicated) { // (a client only flies them: hits are the server's)
                const Vector3 pu = Vector3Normalize(feet);
                if (Melee_SegmentDistance(before, p.pos, feet, Vector3Add(feet, Vector3Scale(pu, 1.7f))) < 0.45f) {
                    DamageEvent e;
                    e.amount = p.damage * Difficulty_ThreatScale();
                    e.kind = DamageKind::Projectile; e.type = (p.spear || p.dart) ? DamageType::Pierce : DamageType::Blunt;
                    e.parryable = false; e.hasSource = true;
                    e.source = body->toWorldD(DVec3(Vector3Subtract(p.pos, Vector3Scale(Vector3Normalize(p.vel), 2.0f)))); // where it came from
                    e.knockback = 2.0f; e.sourceId = p.owner;
                    const DamageResult r = Combat_Damage(LocalSession(), e);
                    GameUI_Toast(p.dart ? (r.dodged ? "You dodge a dart" : r.blocked ? "Your shield takes a dart" : "A dart hits you!") : r.dodged ? "You dodge a thrown stone" : r.blocked ? "Your shield takes a thrown stone" : p.spear ? "A spear hits you!" : "A thrown stone hits you!");
                    ++projectileHits;
                    gone = true;
                }
            }
            // Remote players (a server): their feet as the animals see them.
            if (!gone)
                for (const AnimalTarget& t : gAnimalTargets) {
                    if (Planetoid::findNet(t.body) != body || !t.vulnerable) continue;
                    const Vector3 tu = Vector3Normalize(t.local);
                    if (Melee_SegmentDistance(before, p.pos, t.local, Vector3Add(t.local, Vector3Scale(tu, 1.7f))) >= 0.45f) continue;
                    DamageEvent e; e.amount = p.damage * Difficulty_ThreatScale(); e.kind = DamageKind::Projectile; e.parryable = false;
                    e.hasSource = true; e.source = body->toWorldD(DVec3(before)); e.knockback = 2.0f; e.sourceId = p.owner;
                    if (gAnimalDamageRemote) gAnimalDamageRemote(t.session, e);
                    ++projectileHits;
                    gone = true;
                    break;
                }
        }
        if (gone) { if (gNetRole == NetRole::Host && fxEnded.size() < 128) fxEnded.push_back(p.id); projectiles[k] = projectiles.back(); projectiles.pop_back(); }
        else ++k;
    }
}

void DrawProjectiles(Scene& world) {
    for (const Projectile& p : projectiles) {
        Planetoid* body = BodyNamed(world, p.body);
        if (!body) continue;
        const Vector3 w = body->toWorld(p.pos);
        if (Vector3Distance(w, world.camera.position) > 120.0f) continue;
        if (p.dart) { // (#1228) a glowing bolt along its flight
            const Vector3 d = Vector3Normalize(body->dirToWorld(p.vel));
            for (int s = 0; s < 3; ++s) IB_Cube(Vector3Add(w, Vector3Scale(d, -0.12f * s)), { 0.05f, 0.05f, 0.05f }, s == 0 ? Color{ 255, 120, 80, 255 } : Color{ 160, 160, 170, 255 });
            continue;
        }
        if (!p.spear) { IB_Cube(w, { 0.16f, 0.14f, 0.16f }, Color{ 120, 115, 105, 255 }); continue; }
        const Vector3 d = Vector3Normalize(body->dirToWorld(p.vel));
        for (int s = -2; s <= 2; ++s) IB_Cube(Vector3Add(w, Vector3Scale(d, 0.3f * s)), { 0.05f, 0.05f, 0.05f }, s == 2 ? Color{ 90, 90, 95, 255 } : Color{ 120, 90, 60, 255 });
    }
}
} // namespace


// (#1343, #1364) The host: what happened since the last call (shots launched, shots ended, creature events), moved out.
void Animals_TakeFx(MsgCreatureFx& out) {
    out.shots.swap(fxShots); fxShots.clear();
    out.ended.swap(fxEnded); fxEnded.clear();
    out.events.swap(fxEvents); fxEvents.clear();
}

// (#1343, #1364) A client: the server's shots fly here too (a flight only: hits are the server's), ended ones go, and creature
// events reach this machine's hook (the sounds and captions, #1363).
void Animals_ApplyFx(const MsgCreatureFx& m) {
    LoadSpecies();
    for (const MsgCreatureFx::Shot& s : m.shots) {
        const Planetoid* b = Planetoid::findNet(s.body);
        if (!b) continue;
        projectiles.push_back({ b->name, s.pos.toFloat(), s.vel, 0.0f, s.owner, simClock, s.kind == 1, s.kind == 2, s.id });
    }
    for (uint32_t id : m.ended)
        for (size_t k = 0; k < projectiles.size(); ++k)
            if (projectiles[k].id == id) { projectiles[k] = projectiles.back(); projectiles.pop_back(); break; }
    for (const MsgCreatureFx::Event& e : m.events) {
        if (e.species >= species.size()) continue;
        const CreatureSpec* cs = Creature_Spec(species[e.species].id);
        const Planetoid* b = Planetoid::findNet(e.body);
        if (cs && b) Creature_Event(*cs, (CreatureEvent)e.kind, b->id, e.pos.toFloat());
    }
}

int Animals_ProjectilesThrown() { return projectilesThrown; }
int Animals_ProjectileHits() { return projectileHits; }

namespace {
// (#1698) A body-local double point in render space: through the body's double path, not its float toWorld (at planet scale a
// float local has ~4 mm steps, and the animal stepped with it).
Vector3 RenderPos(const Scene& world, const Planetoid& body, const DVec3& local) { return (body.toWorldD(local) - world.renderOrigin).toFloat(); }
}

// (#1225) Torches carried by creatures (torch=1): a warm light each, for the frame's local lights (render space).
void Animals_AddLights(Scene& world) {
    for (const Animal& a : animals) {
        const CreatureSpec* cs = Creature_Spec(species[a.species].id);
        if (cs && cs->ranged && cs->dart && a.throwWindAt >= 0.0) { // (#1228) a sentinel charging: its eye lights the room red
            if (Planetoid* sb = BodyNamed(world, a.body)) {
                const float ch = std::clamp((float)((simClock - a.throwWindAt) / std::max(0.1f, cs->throwWindUp)), 0.0f, 1.0f);
                LocalLight L; L.position = RenderPos(world, *sb, a.pos + DVec3(Vector3Normalize(a.pos.toFloat())) * 1.0); L.range = 4.0f + 4.0f * ch;
                L.color = { 1.0f, 0.3f + 0.5f * ch, 0.2f + 0.5f * ch }; L.intensity = 0.6f + 1.2f * ch;
                if (Vector3Distance(L.position, world.camera.position) < 70.0f) world.localLights.push_back(L);
            }
            continue;
        }
        if (!cs || !cs->torch || a.thrower) continue; // (the thrower's hands hold stones)
        Planetoid* body = BodyNamed(world, a.body);
        if (!body) continue;
        const Vector3 w = RenderPos(world, *body, a.pos + DVec3(Vector3Normalize(a.pos.toFloat())) * 1.5);
        if (Vector3Distance(w, world.camera.position) > 70.0f) continue;
        LocalLight L;
        L.position = w; L.range = 9.0f; L.color = { 1.0f, 0.62f, 0.32f }; L.intensity = 1.3f;
        world.localLights.push_back(L);
    }
}

// ---- (#1223) Cave creatures: spawning ----
// Every 2 s around each interest centre (players) on a planet with graph caves: the cave segments within 64 m give spawn
// sites (their floors, found down from the axis), each with its depth, light and whether an entrance is near. Each cave
// species (creatures.txt habitat=cave or mouth) fills each 64 m chunk where it fits up to its cap on the first visit, then
// adds one every respawn seconds while below it; a species with nests gets a nest in a chunk's large chamber, and once that
// nest is destroyed the chunk stays clear (decided). None on peaceful, none inside claims, none within 24 m of the player.
namespace {
struct ChunkId { int x, y, z, species; bool operator<(const ChunkId& o) const { return std::tie(x, y, z, species) < std::tie(o.x, o.y, o.z, o.species); } };
std::map<ChunkId, SpawnChunk> spawnChunks;
double lastCaveSpawn = -1e9;
int caveSpawned = 0, caveSites = 0;
int lurkerDrops = 0; // (#1224) drops that landed on a player (tests)
int thefts = 0;      // (#1225) piles taken by thieves (tests)

ChunkId ChunkOf(Vector3 p, int si) {
    return { (int)floorf(p.x / SpawnTune::kChunk), (int)floorf(p.y / SpawnTune::kChunk), (int)floorf(p.z / SpawnTune::kChunk), si };
}

// (#1633) Populations: a cave creature (cave or mouth habitat, and the nests) belongs to the chunk of its home, wherever it
// has wandered (one that chased the player out was no longer counted, so its chunk spawned another every respawn: a save
// gathered 294 around the spawn). Around each player a species also has an area budget, and creatures far from every
// player are let go (their chunk fills again on the next visit).
bool CavePopulation(int si) {
    if (si < 0 || si >= (int)species.size()) return false;
    if (species[si].id == "creature_nest") return true;
    const CreatureSpec* c = Creature_Spec(species[si].id);
    return c && (c->habitat == Habitat::Cave || c->habitat == Habitat::Mouth);
}
ChunkId HomeChunk(const Animal& a) { return ChunkOf(a.hasHome ? a.home : a.pos.toFloat(), a.species); }
// (#1717) RAYLIBGAME_CAVE_SETTLE=0: the old settle (a cave creature off its feet could be put on the surface; A/B)
bool CaveSettleFix() { static const bool on = !(getenv("RAYLIBGAME_CAVE_SETTLE") && getenv("RAYLIBGAME_CAVE_SETTLE")[0] == '0'); return on; }
namespace PopTune {
    constexpr float kAreaRange = 128.0f;  // the area budget's radius around a player (GUESS)
    constexpr int kAreaPerSpecies = 8;    // at most this many of one cave species within it (GUESS: four chunks' worth)
    constexpr float kReleaseRange = 320.0f; // farther than this from every player (past the 300 m simulation): let go
}
bool gPopulationCleanup = false; // set by a load: trim the loaded cave creatures once (Animals_Update)

// (#1227) No rock within 40 m over a body-local point (2 m steps): under the open sky.
// (#1227) The creature facing rules (wind-up holds facing, a turn rate, strikes only in front); RAYLIBGAME_CREATURE_TURN=0: off (A/B).
static int CreatureTurnMask() { static const int m = getenv("RAYLIBGAME_CREATURE_TURN") ? atoi(getenv("RAYLIBGAME_CREATURE_TURN")) : 7; return m; } // bits: 1 the wind-up hold, 2 the turn rate, 4 strikes in front only
static bool CreatureTurnRules() { return CreatureTurnMask() != 0; }
// (#1649) Column probes (open sky over a point, no roof within reach) cached per 1 m cell for a few seconds: a crowd asked
// them every step (20-79 density samples each; a lurker idling in the open searched for a ceiling every step).
namespace ColumnCache {
    constexpr double kKeep = 3.0; // s (sim time): a dug or built column is seen within this
    std::unordered_map<uint64_t, std::pair<bool, double>> cells;
    uint64_t Key(const Planetoid& body, Vector3 p, uint64_t kind) {
        const uint64_t x = (uint64_t)(int64_t)floorf(p.x) & 0xFFFFF, y = (uint64_t)(int64_t)floorf(p.y) & 0xFFFFF, z = (uint64_t)(int64_t)floorf(p.z) & 0xFFFFF;
        return (x | (y << 20) | (z << 40) | (kind << 60)) ^ ((uint64_t)(uintptr_t)&body * 0x9E3779B97F4A7C15ull);
    }
    std::mutex m; // (#1650 slice 3) the step's probes run on worker threads too: the map under a lock, the probe outside it
    template <class F> bool Get(const Planetoid& body, Vector3 p, uint64_t kind, F probe) {
        const double now = SimTime();
        const uint64_t key = Key(body, p, kind);
        {
            std::lock_guard<std::mutex> lk(m);
            if (cells.size() > 200000) cells.clear();
            auto it = cells.find(key);
            if (it != cells.end() && now - it->second.second <= kKeep) return it->second.first;
        }
        const bool v = probe();
        std::lock_guard<std::mutex> lk(m);
        cells[key] = { v, now };
        return v;
    }
}
static bool OpenSkyOver(const Planetoid& body, Vector3 p) {
    return ColumnCache::Get(body, p, 0, [&] {
        const Vector3 up = Vector3Normalize(p);
        for (float h = 2.0f; h <= 40.0f; h += 2.0f) if (body.densityAt(Vector3Add(p, Vector3Scale(up, h))) >= ISO_LEVEL) return false;
        return true;
    });
}

// (#1650 slice 3) The step in three passes: the update loop decides and records each mover's step (a MoveReq); the
// terrain probes for all of them run in parallel (each request reads only the terrain, the parts and its own record);
// then the results are applied in the loop's order on the main thread (path plans, positions, the stuck rule). Same
// results whatever the threads: a request never reads another's. RAYLIBGAME_AI_THREADS=0: the probes on one thread.
struct MoveReq {
    uint64_t net = 0;
    const Planetoid* body = nullptr; const PlanetField* field = nullptr;
    Vector3 east{}, north{}, curPos{}; DVec3 pos{};
    float stepLen = 0.0f, height = 0.0f, reachUp = 0.0f, reachDown = 0.0f, speed = 0.0f, step = 0.0f, heading = 0.0f, detourSide = 0.0f, sizeX = 0.0f;
    bool leashed = false, outside = false, steer = false, closeHunt = false; // closeHunt (#1717): sidesteps only, no flip
    Vector3 goalPt{};
    uint32_t rng = 1;
};
struct MoveRes {
    bool ok0 = false; Vector3 at0{}; DVec3 next0{};          // the step straight on
    bool detourTried = false, detourOk = false;              // the detour (steering and blocked), as the serial code ran it
    float detourSideInit = 0.0f, detourTurn = 0.0f, detourSide = 0.0f; uint32_t rngAfter = 1; Vector3 atD{}; DVec3 nextD{};
};
// One step's probe along `heading` (the old tryStep): ground within reach, room for the body, no part in the way, a leashed
// cave creature not stepping out under the open sky.
bool ProbeStep(const MoveReq& r, float heading, Vector3& at, DVec3& nextPos) {
    const Vector3 fwd = Vector3Add(Vector3Scale(r.east, cosf(heading)), Vector3Scale(r.north, sinf(heading)));
    const DVec3 target = r.pos + DVec3(fwd) * (double)r.stepLen; // the stride is added in doubles (#1698)
    const Vector3 nd = Vector3Normalize(target.toFloat());
    float e;
    if (!Ground(*r.body, *r.field, nd, at, e, Vector3Length(r.curPos), r.reachUp, r.reachDown)) return false;
    nextPos = target * (Vector3Length(at) / target.length()); // the probe gives the height; the sideways position stays double
    return r.body->densityAt(Vector3Add(at, Vector3Scale(nd, r.height * 0.8f))) < ISO_LEVEL + kRockMargin &&
           !BlockedAlong(*r.body, r.curPos, at, r.height) &&
           !(r.leashed && !r.outside && OpenSkyOver(*r.body, at)); // (#1227) a cave creature doesn't step out under the sky
}
// The straight step, and the detour a blocked steering creature would take (its side from a copy of its random state: the
// real one moves only if the detour is used, pass 3).
void ProbeMove(const MoveReq& r, MoveRes& o) {
    o.ok0 = ProbeStep(r, r.heading, o.at0, o.next0);
    if (o.ok0 || !r.steer) return;
    o.detourTried = true;
    uint32_t rng = r.rng;
    float side0 = r.detourSide;
    if (side0 == 0.0f) side0 = Rand01(rng) < 0.5f ? 1.0f : -1.0f;
    o.detourSideInit = side0; o.rngAfter = rng;
    for (float turn : { 0.5f, 1.0f, 1.6f, 2.2f, 3.1f }) {
        if (r.closeHunt && turn > 1.05f) break; // (#1717) close in on its prey: sidesteps only
        for (float side : { side0, -side0 })
            if (!o.detourOk && ProbeStep(r, r.heading + side * turn, o.atD, o.nextD)) { o.detourOk = true; o.detourTurn = turn; o.detourSide = side; }
    }
}
// Work items 0..n-1 on the main thread and up to 6 pool workers: the main thread takes items too, so it never waits on a
// worker still busy with a chunk build, only on items already started. Late helpers find nothing left and return.
template <class F> void AiParallelFor(int n, F f) {
    static const bool on = !(getenv("RAYLIBGAME_AI_THREADS") && getenv("RAYLIBGAME_AI_THREADS")[0] == '0');
    if (!on || n < 24) { for (int k = 0; k < n; ++k) f(k); return; }
    struct Shared { std::atomic<int> next{ 0 }, done{ 0 }; int n = 0; std::function<void(int)> fn; };
    auto sh = std::make_shared<Shared>();
    sh->n = n; sh->fn = f;
    constexpr int kChunk = 8;
    auto work = [sh] {
        for (int k; (k = sh->next.fetch_add(kChunk)) < sh->n;) {
            const int e = std::min(k + kChunk, sh->n);
            for (int j = k; j < e; ++j) sh->fn(j);
            sh->done.fetch_add(e - k);
        }
    };
    const int helpers = std::min(6, n / kChunk - 1);
    for (int h = 0; h < helpers; ++h) WorldJobs().submit(-1e30f, work);
    work();
    while (sh->done.load() < n) std::this_thread::yield();
}

// An animal placed at a body-local point as it is (a cave floor), not dropped from above as Animals_Spawn does.
uint64_t SpawnAt(Planetoid& body, int si, Vector3 at, uint32_t seed) {
    Animal a;
    a.species = si;
    a.body = body.name;
    a.pos = a.prevPos = DVec3(at);
    a.home = at; a.hasHome = true; // (#1227) the leash's anchor
    a.heading = 2.0f * PI * Rand01(seed);
    a.health = species[si].health;
    a.age = 1.0f;
    a.rng = seed ^ (uint32_t)(animals.size() * 7919u);
    a.timer = 1.0f + 4.0f * Rand01(seed);
    a.settleIn = 1.0f;
    a.productReady = WorldTime();
    a.netId = NetId_Allocate(NetKind::Animal);
    animals.push_back(a);
    return a.netId;
}

void SpawnCaveCreatures(Planetoid& body, const PlanetField& f, Vector3 centre, Vector3 player) {
    std::vector<std::pair<int, const CreatureSpec*>> kinds;
    for (const CreatureSpec& c : Creature_All())
        if ((c.habitat == Habitat::Cave || c.habitat == Habitat::Mouth) && (!body.profile || body.profile->allowsCreature(c.id))) { // (#1274) the world file's creatures
            const int si = SpeciesIndex(c.id); if (si >= 0) kinds.push_back({ si, &c });
        }
    if (kinds.empty()) return;
    const float R = SpawnTune::kChunk;
    // Within 64 m sideways of the player, down to 200 m under them (caves run deep under a surface player) and 30 m over.
    const float kDown = 200.0f, kUp = 30.0f;
    const Vector3 cu = Vector3Normalize(centre);
    const Vector3 lo = Vector3SubtractValue(centre, kDown), hi = Vector3AddValue(centre, kDown);
    std::vector<PlanetField::CaveSpan> spans;
    std::vector<PlanetField::CaveEntrance> ents;
    f.caveSpansNear(lo, hi, spans);
    f.caveEntrancesNear(Vector3SubtractValue(lo, 20.0f), Vector3AddValue(hi, 20.0f), ents);
    std::vector<SpawnSite> sites;
    int nSolid = 0, nNoFloor = 0, nNear = 0, nClaim = 0, nFar = 0, nSky = 0;
    for (const PlanetField::CaveSpan& sp : spans) {
        if (sp.kind == 2) continue; // rubble is solid
        const float len = Vector3Distance(sp.a, sp.b);
        const int n = std::max(1, (int)(len / 12.0f));
        for (int k = 0; k < n; ++k) {
            const float t = (k + 0.5f) / n;
            const Vector3 axis = Vector3Lerp(sp.a, sp.b, t);
            { const Vector3 rel = Vector3Subtract(axis, centre); const float v = Vector3DotProduct(rel, cu);
              if (Vector3Length(Vector3Subtract(rel, Vector3Scale(cu, v))) > R || v < -kDown || v > kUp) { ++nFar; continue; } }
            const float r = sp.ra + (sp.rb - sp.ra) * t;
            const Vector3 up = Vector3Normalize(sp.up);
            if (body.densityAt(axis) >= ISO_LEVEL) { ++nSolid; continue; }          // the axis itself in rock (a squeezed span): no room
            Vector3 floor{}; bool found = false;
            for (float d = 0.25f; d <= r + 2.0f; d += 0.25f) {
                const Vector3 q = Vector3Subtract(axis, Vector3Scale(up, d));
                if (body.densityAt(q) >= ISO_LEVEL) { floor = Vector3Add(q, Vector3Scale(up, 0.25f)); found = true; break; }
            }
            if (!found || body.densityAt(Vector3Add(floor, Vector3Scale(up, 1.0f))) >= ISO_LEVEL) { ++nNoFloor; continue; } // no headroom
            if (Vector3Distance(floor, player) < SpawnTune::kMinFromPlayer) { ++nNear; continue; }
            if (Claims_At(body.name, floor)) { ++nClaim; continue; } // never inside claimed land
            SpawnSite s;
            s.floor = floor; s.up = up; s.kind = sp.kind; s.radius = r;
            const Vector3 dir = Vector3Normalize(floor);
            s.depth = f.radius + f.elevation(dir) - Vector3Length(floor);
            for (const PlanetField::CaveEntrance& e : ents) if (Vector3Distance(e.at, floor) < 20.0f) { s.nearMouth = true; break; }
            // (#1592) Under rock? The depth is the smooth field surface over the floor, so a span through an outcrop, an eroded
            // valley or an opened roof put "deep" sites on open ground: lurkers on an Alpine outcrop, slimes on the grass,
            // the Lost Miner walking the surface by day. Up from the floor in 2 m steps to 80 m: no rock, open sky.
            s.openSky = true;
            for (float h = 2.0f; h <= 80.0f; h += 2.0f)
                if (body.densityAt(Vector3Add(floor, Vector3Scale(up, h))) >= ISO_LEVEL) { s.openSky = false; break; }
            if (s.openSky) ++nSky;
            s.light = s.nearMouth ? 0.5f : s.depth < 4.0f || s.openSky ? 0.6f : 0.0f;
            s.wet = f.caveWetColumn(dir, f.elevation(dir)) > 0.3f; // (#1226) wet karst: by the cave water, the sumps
            sites.push_back(s);
        }
    }
    caveSites = (int)sites.size();
    if (getenv("RAYLIBGAME_CAVE_SPAWN_LOG")) printf("cave spawn: spans %d, far %d, solid %d, no floor %d, near %d, claimed %d, sites %d (%d under open sky)\n", (int)spans.size(), nFar, nSolid, nNoFloor, nNear, nClaim, caveSites, nSky);
    if (sites.empty()) return;
    const int nestSi = SpeciesIndex("creature_nest");
    const double now = SimTime();
    for (const auto& kv : kinds) {
        const int si = kv.first; const CreatureSpec& c = *kv.second;
        std::map<ChunkId, std::vector<const SpawnSite*>> byChunk;
        for (const SpawnSite& s : sites) if (Spawn_Fits(c, s)) byChunk[ChunkOf(s.floor, si)].push_back(&s);
        for (auto& ck : byChunk) {
            SpawnChunk& ch = spawnChunks[ck.first];
            int alive = 0, inArea = 0;
            for (const Animal& a : animals) {
                if (a.species != si || a.body != body.name) continue;
                const ChunkId id = HomeChunk(a); // (#1633) by its home, wherever it wandered
                alive += !(id < ck.first) && !(ck.first < id);
                inArea += Vector3Distance(a.pos.toFloat(), centre) < PopTune::kAreaRange ? 1 : 0;
            }
            if (c.nests && nestSi >= 0) {
                if (!ch.hasNest) // the nest in the chunk's large chamber (none without one)
                    for (const SpawnSite* s : ck.second)
                        if (s->kind == 3 && s->radius >= SpawnTune::kNestChamber) { ch.nest = SpawnAt(body, nestSi, s->floor, Hash((uint64_t)ck.first.x * 73856093u ^ (uint64_t)ck.first.z, 3u)); ch.hasNest = true; break; }
                if (ch.hasNest && IndexOfNet(ch.nest) < 0) ch.nestDestroyed = true;
            }
            int due = Spawn_Due(c, ch, alive, now);
            due = std::min(due, std::max(0, PopTune::kAreaPerSpecies - inArea)); // (#1633) the area budget around the player
            inArea += due;
            uint32_t h = Hash(((uint64_t)(uint32_t)ck.first.x << 32) ^ (uint32_t)ck.first.y ^ ((uint64_t)ck.first.z << 16), (uint32_t)(now * 10.0) + (uint32_t)si);
            for (int k = 0; k < due;) {
                const SpawnSite* s = ck.second[std::min((int)ck.second.size() - 1, (int)(Rand01(h) * ck.second.size()))];
                // (#1225) A group species comes as a group at one site: up to `group` together (within what is due), the first
                // its thrower when the species has one.
                const int n = std::max(1, std::min(c.group, due - k));
                uint64_t first = 0;
                for (int m = 0; m < n; ++m) {
                    const uint64_t id = SpawnAt(body, si, s->floor, h);
                    if (const int j = IndexOfNet(id); j >= 0 && c.group > 1) {
                        if (!first) first = id;
                        animals[j].group = first;
                        animals[j].thrower = c.thrower && m == 0;
                    }
                    ++caveSpawned;
                }
                k += n;
            }
        }
    }
}
} // namespace

int Animals_CaveSpawned() { return caveSpawned; }
bool Animals_Clinging(int i) { return i >= 0 && i < (int)animals.size() && animals[i].clinging; }
Vector3 Animals_ClingFloor(int i) { return i >= 0 && i < (int)animals.size() ? animals[i].clingFloor : Vector3{ 0, 0, 0 }; }
int Animals_LurkerDrops() { return lurkerDrops; }
int Animals_Thefts() { return thefts; }
int Animals_ShieldBlocks() { return shieldBlocks; }
const char* Animals_MoodName(int i) { return i >= 0 && i < (int)animals.size() ? Mood_Name(animals[i].brain.mood) : "-"; }
int Animals_OozeSplits() { return oozeSplits; }
int Animals_Generation(int i) { return i >= 0 && i < (int)animals.size() ? animals[i].generation : -1; }
const std::string& Animals_Ore(int i) { static const std::string none; return i >= 0 && i < (int)animals.size() ? animals[i].ore : none; }
bool Animals_Thrower(int i) { return i >= 0 && i < (int)animals.size() && animals[i].thrower; }
void Animals_TestSetGroup(int i, uint64_t group, bool thrower) { if (i >= 0 && i < (int)animals.size()) { animals[i].group = group; animals[i].thrower = thrower; } }
void Animals_TestProvokeAll(int part, int parts) { for (size_t i = 0; i < animals.size(); ++i) if (parts <= 1 || (int)(i % (size_t)parts) == part) { Animal& a = animals[i]; a.hurtPending = true; a.hurtBy = kPlayerPercept; a.angryUntil = simClock + 30.0; } } // (#1649; #1650: one part of them a call)
float Animals_Heading(int i) { return i >= 0 && i < (int)animals.size() ? animals[i].heading : 0.0f; } // (#1717)
void Animals_TestSetHeading(int i, float h) { if (i >= 0 && i < (int)animals.size()) { animals[i].heading = animals[i].prevHeading = h; } } // (#1717)
double Animals_Clock() { return simClock; }
int Animals_CaveSites() { return caveSites; }
void Animals_CaveSpawnNow() { lastCaveSpawn = -1e9; }
void Animals_ResetCaveSpawns() { spawnChunks.clear(); caveSpawned = 0; caveSites = 0; lastCaveSpawn = -1e9; }

void Animals_Update(Scene& world, float dt, float worldDt, bool inputEnabled) {
    const auto t0 = std::chrono::steady_clock::now();
    simClock += dt; // the attack wind-up and the bite cooldown run on it (#165): tests step animals many times a frame
    LoadSpecies();
    ++frameNo;
    simulatedNow = 0; gAiGridDirty = true; // (#1650) positions moved: the grid is rebuilt on its first use this update
    if (gPopulationCleanup && !gAnimalsReplicated) { // (#1633) once after a load: cave creatures from before homes were saved
        gPopulationCleanup = false;                  // (wanderers with no way back) go, and each chunk keeps its cap
        std::map<ChunkId, int> kept;
        int dropped = 0;
        for (size_t i = 0; i < animals.size();) {
            const Animal& a = animals[i];
            bool drop = false;
            if (CavePopulation(a.species)) {
                const CreatureSpec* c = Creature_Spec(species[a.species].id);
                drop = !a.hasHome || ++kept[HomeChunk(a)] > (c ? c->capPerChunk : 1);
            }
            if (drop) { animals.erase(animals.begin() + (long)i); ++dropped; } else ++i;
        }
        if (dropped) { aimed = -1; printf("animals: let %d cave creatures go after the load (#1633: no home saved, or over their chunk's cap)\n", dropped); }
    }
    for (Animal& a : animals) { a.prevPos = a.pos; a.prevHeading = a.heading; a.prevValid = true; } // the pose drawing blends from (#315)
    if (gAnimalsReplicated) { // #334: a client - the server runs the animals; each eases toward its latest state, nothing decides here
        const double now = GetTime();
        const float k = 1.0f - expf(-dt / 0.08f); // ~80 ms to close a gap (states come 10 a second)
        for (size_t i = 0; i < animals.size();) {
            Animal& a = animals[i];
            if (now - a.netSeen > 3.0) { animals.erase(animals.begin() + i); if (aimed == (int)i) aimed = -1; continue; } // out of view (or the removal was lost)
            const DVec3 before = a.pos;
            a.pos += (DVec3(a.netPos) - a.pos) * (double)k;
            float dh = a.netHeading - a.heading; while (dh > PI) dh -= 2 * PI; while (dh < -PI) dh += 2 * PI;
            a.heading += dh * k;
            const float step = (float)(a.pos - before).length();
            a.moveSpeed = dt > 0.0f ? step / dt : 0.0f;
            if (step > 1e-4f) a.gait += step / std::max(0.25f, species[a.species].size.x * 0.45f) * PI; // the walk cycle, as the AI drives it
            ++i;
        }
        StepProjectiles(world, dt); // (#1343) the server's shots, flown here (no hits: those are the server's)
        lastUpdateMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return;
    }
    Planetoid* here = player.nearest;
    const PlanetField* pf = here ? dynamic_cast<const PlanetField*>(here->field.get()) : nullptr;
    Vector3 playerLocal{};
    if (here) playerLocal = here->toLocalOf(player.position, world.renderOrigin);
    // (#1223) Players' noise, for creatures' hearing. Footsteps and landings come from the footsteps module (#1240: each step's
    // loudness by surface, pace and posture - crouched and prone are quiet); 1 there (walking on stone) is 0.4 here (a walk is
    // heard within 0.4 x a creature's hearing). Each dig or swing stroke of the local player is loud (melee hits: MeleeFrame).
    {
        static double lastStepTime = -1e9;
        static int lastSwings = 0;
        for (const StepNoise& n : Footsteps_RecentNoises()) {
            if (n.time <= lastStepTime) continue;
            lastStepTime = n.time;
            Planetoid* b = nullptr;
            if (n.who == 0) b = here;
            else for (const AnimalTarget& t : gAnimalTargets) if (t.session == n.who) { b = Planetoid::findNet(t.body); break; }
            if (!b) continue;
            Noise_Emit(b->id, b->toLocalD(n.pos).toFloat(), n.loudness * 0.4f, simClock, kPlayerPercept | n.who);
        }
        if (here && !player.flyMode && !survival.dead && terrainTool.swings != lastSwings) {
            lastSwings = terrainTool.swings;
            Noise_Emit(here->id, playerLocal, NoiseTune::kDig, simClock, kPlayerPercept);
        }
    }
    // Interest centres (#318): the local player and every other centre (other players on a server). Herds spawn and animals
    // are simulated around all of them; player-specific behaviour (fleeing the player, following food in hand, night
    // hunters, campfires) still looks at the local player only.
    static std::vector<InterestCentre> centres; // (reused: no allocation per update)
    centres.clear();
    if (here) centres.push_back({ here->id, playerLocal });
    for (const InterestCentre& c : gInterestCentres) centres.push_back(c);
    auto nearestCentre = [&](const Planetoid* body, Vector3 at) {
        float best = 1e9f;
        for (const InterestCentre& c : centres) if (body && c.body == body->id) best = std::min(best, Vector3Distance(at, c.local));
        return best;
    };
    if (naturalSpawns && ((here && pf && !player.flyMode) || !gInterestCentres.empty()) && SimTime() - lastSpawnCheck > 1.0) { // naturalSpawns: m56 (#373) // (single-player: as before)
        lastSpawnCheck = SimTime();
        for (size_t k = 0; k < centres.size(); ++k) {
            const InterestCentre& c = centres[k];
            Planetoid* b = Planetoid::find(c.body);
            const PlanetField* cf = b ? dynamic_cast<const PlanetField*>(b->field.get()) : nullptr;
            if (!cf || (here && k == 0 && player.flyMode)) continue; // (the local player flying: no herds under the camera)
            const float d = Vector3Length(c.local);
            if (d > 1.0f && d - cf->radius < 400.0f) SpawnHerds(world, *b, *cf, Vector3Scale(c.local, 1.0f / d));
        }
    }
    const double now = WorldTime();
    // #334 the players animals react to: the local one (not a dedicated server's idle stand-in) and every remote player the
    // host lists (gAnimalTargets). Each animal minds the nearest on its body: flees it, follows its food, hunts it, bites it.
    struct Tgt { Planetoid* body; Vector3 local; uint32_t session; bool vulnerable; std::string held; bool lit; };
    static std::vector<Tgt> targets; targets.clear();
    if (here && !(gHeadless && gNetRole == NetRole::Host)) {
        // (#560) a light in hand, or in the off hand (as the game lights it: main.cpp's held torch)
        const std::string inHand = Hotbar_Selected() >= 0 ? Hotbar_Slot(Hotbar_Selected()) : std::string();
        const std::string& off = inventory.offHand().id;
        const bool lit = (IsCarriedLight(inHand) && inventory.count(inHand) > 0) || (IsCarriedLight(off) && inventory.count(off) > 0);
        targets.push_back({ here, playerLocal, 0, !survival.dead && !player.flyMode && !Rocket_Boarded(), Hotbar_HeldFood(), lit });
    }
    for (const AnimalTarget& t : gAnimalTargets) if (Planetoid* b = Planetoid::findNet(t.body)) targets.push_back({ b, t.local, t.session, t.vulnerable, t.held, IsCarriedLight(t.held) }); // (remote: the main hand only)
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    // Light sources near the players (placed ones refreshed every second; carried ones every step): night hunters keep away.
    static std::vector<Vector3> placedLights;
    if (!targets.empty() && SimTime() - lastFireScan > 1.0) { // (#334: near every player, on each one's body)
        lastFireScan = SimTime();
        placedLights.clear();
        static int lampMat = -2;
        if (lampMat == -2) lampMat = BuildMaterialForItem("lamp");
        static std::vector<std::pair<Int3, CellShape>> fp;
        std::vector<Planetoid*> scanned;
        for (const Tgt& t : targets) {
            if (!t.body->plots || std::find(scanned.begin(), scanned.end(), t.body) != scanned.end()) continue;
            scanned.push_back(t.body);
            for (const auto& pl : t.body->plots->plots)
                for (size_t i = 0; i < pl->parts.size(); ++i) {
                    const Part& p = pl->parts[i];
                    bool light = p.type == PartType::Campfire;
                    if (!light && AnyLightScares()) {
                        const bool fireMachine = p.type == PartType::Kiln || p.type == PartType::CokeOven || p.type == PartType::BlastFurnace ||
                                                 p.type == PartType::Furnace || p.type == PartType::AlloySmelter || p.type == PartType::ArcFurnace;
                        light = p.type == PartType::Torch || (lampMat >= 0 && IsShape(p.type) && p.material == lampMat) ||
                                (fireMachine && Factory::Working(t.body->id, pl->id, (int)i));
                    }
                    if (!light) continue;
                    Plot::Footprint(p, fp);
                    if (fp.empty()) continue;
                    Vector3 sum{};
                    for (auto& cell : fp) sum = Vector3Add(sum, { (cell.first.x + 0.5f) * kCellSize, (cell.first.y + 0.5f) * kCellSize, (cell.first.z + 0.5f) * kCellSize });
                    const Vector3 c = pl->toBody(Vector3Scale(sum, 1.0f / fp.size()));
                    bool nearOne = false;
                    for (const Tgt& u : targets) nearOne = nearOne || (u.body == t.body && Vector3Distance(c, u.local) < 200.0f);
                    if (nearOne) placedLights.push_back(c);
                }
        }
    }
    fires = placedLights;
    if (AnyLightScares()) for (const Tgt& t : targets) if (t.lit) fires.push_back(t.local); // (#560) a carried light moves with its player
    // Night threats (not on peaceful): every 12 s after dark, a dusk stalker comes out of the dark 45-60 m
    // away, up to 2 near the player (3 on hard) - not underground, flying, or by a campfire.
    if (!targets.empty() && threatsOn && Difficulty_ThreatScale() > 0.0f && SimTime() - lastNightSpawn > 12.0) {
        lastNightSpawn = SimTime();
        int hunter = -1;
        for (size_t k = 0; k < species.size(); ++k) if (species[k].night) { hunter = (int)k; break; }
        for (const Tgt& t : targets) { // (#334: around every player out in the dark, as around the local one)
        const PlanetField* tf = dynamic_cast<const PlanetField*>(t.body->field.get());
        if (!tf || !t.vulnerable || hunter < 0) continue;
        if (t.body->profile && !t.body->profile->allowsCreature(species[hunter].id)) continue; // the world file's creatures (#1274)
        const Vector3 up = Vector3Normalize(t.local);
        const bool dark = Vector3DotProduct(Vector3Normalize(world.sunDirection), t.body->dirToWorld(up)) < -0.1f;
        const bool surface = Vector3Length(t.local) > tf->radius + tf->elevation(up) - 6.0f;
        Vector3 fire;
        if (dark && surface && !NearFire(t.local, 10.0f, fire)) {
            int near = 0;
            for (const Animal& a : animals) near += species[a.species].night && a.body == t.body->name && Vector3Distance(a.pos.toFloat(), t.local) < 100.0f ? 1 : 0;
            if (near < (worldDifficulty == 2 ? 3 : 2)) {
                Vector3 east, north;
                TangentFrame(up, east, north);
                for (int attempt = 0; attempt < 8; ++attempt) { // a dry spot (not water, not by a fire)
                const float ang = 2.0f * PI * GameRng_Unit(), r = 45.0f + 15.0f * GameRng_Unit();
                const Vector3 d = Vector3Normalize(Vector3Add(t.local, Vector3Add(Vector3Scale(east, r * cosf(ang)), Vector3Scale(north, r * sinf(ang)))));
                const Vector3 at = Vector3Scale(d, tf->radius + tf->elevation(d));
                if (!NearFire(at, 14.0f, fire) && Animals_Spawn(world, *t.body, species[hunter].id, at, 1, false, 0.5f) > 0) break;
                }
            }
        }
        }
    }
    // (#1223) Cave creatures around each interest centre, every 2 s (none on peaceful, or with threats off: benches).
    if (Creatures_On() && threatsOn && Difficulty_ThreatScale() > 0.0f && SimTime() - lastCaveSpawn > SpawnTune::kPass) {
        lastCaveSpawn = SimTime();
        for (const InterestCentre& c : centres) {
            Planetoid* b = Planetoid::find(c.body);
            const PlanetField* cf = b ? dynamic_cast<const PlanetField*>(b->field.get()) : nullptr;
            if (cf) SpawnCaveCreatures(*b, *cf, c.local, c.local);
        }
    }
    // Ecology (Plan 5 A), twice a second over the animals near the player: predators (the hostile kinds)
    // pick the nearest prey within 30 m to chase (chase only: no kills yet, per the user); prey that spots a
    // predator within 18 m bolts; each herd's centre (same kind within 40 m) for cohesion.
    static double ecoAt = 0.0;
    if (here && SimTime() > ecoAt) {
        ecoAt = SimTime() + 0.5; PROFILE_SCOPE("creatures: ecology"); // (#1649)
        std::vector<int> near;
        for (size_t i = 0; i < animals.size(); ++i)
            if (animals[i].body == here->name && nearestCentre(here, animals[i].pos.toFloat()) < 300.0f) near.push_back((int)i);
        for (int i : near) {
            Animal& a = animals[i];
            const Species& s = species[a.species];
            a.chase = -1;
            Vector3 centre{}; int herd = 0;
            float bestPrey = 30.0f, bestThreat = 18.0f;
            for (int j : near) {
                if (j == i) continue;
                const Animal& b = animals[j];
                const Vector3 bp = b.pos.toFloat();
                const float d = Vector3Distance(a.pos.toFloat(), bp);
                const Species& t = species[b.species];
                if (b.species == a.species && d < 40.0f) { centre = Vector3Add(centre, bp); ++herd; }
                if (s.hostile && !t.hostile && d < bestPrey && threatsOn) { bestPrey = d; a.chase = j; }
                if (!s.hostile && t.hostile && d < bestThreat && threatsOn && !Creature_Spec(s.id)) /* (#1225: a creature's fears are its faction's) */ { bestThreat = d; a.fleeTimer = std::max(a.fleeTimer, 3.0f); a.fleeFrom = bp; }
            }
            a.herdCentre = herd ? Vector3Scale(centre, 1.0f / herd) : a.pos.toFloat();
            a.herdSize = herd;
        }
    }
    static const bool pathsOn = !(getenv("RAYLIBGAME_ANIMAL_PATHS") && getenv("RAYLIBGAME_ANIMAL_PATHS")[0] == '0'); // A/B
    int plansLeft = pathsOn ? 1 : 0; // route plans this update (#253)
    std::vector<MoveReq> moveReqs; // (#1650 slice 3) the steps to probe (pass 2) and apply (pass 3) after the loop
    moveReqs.reserve(animals.size());
    std::vector<int> gone;
    std::unordered_map<std::string, Planetoid*> bodies;
    const auto tLoop0 = std::chrono::steady_clock::now(); // (#1650)
    double accSettle = 0.0, accBrain = 0.0; // (#1650) the loop's settle searches and brain time (the slow-update report)
    double worstBrain = 0.0; std::string worstBrainWho; // (#1650) the slowest single brain block this update, and whose
    struct AccTimer { double& acc; std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now(); ~AccTimer() { acc += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(); } };
    for (size_t i = 0; i < animals.size(); ++i) {
        Animal& a = animals[i];
        const Species& s = species[a.species];
        auto bit = bodies.find(a.body);
        Planetoid* body = bit != bodies.end() ? bit->second : (bodies[a.body] = BodyNamed(world, a.body));
        if (!body) continue;
        const auto* field = dynamic_cast<const PlanetField*>(body->field.get());
        if (!field) continue;
        // Growing up and products run on world time wherever the animal is.
        if (a.age < 1.0f) a.age = std::min(1.0f, a.age + worldDt / 3600.0f);
        // Only near an interest centre (the player, other players: #318) is it simulated; past 80 m a quarter as often.
        const float dist = nearestCentre(body, a.pos.toFloat());
        // (#334) the nearest player on its body: what it flees, follows, hunts (none: it keeps to itself)
        const Tgt* tg = nullptr; float pdist = 1e9f;
        for (const Tgt& t : targets) if (t.body == body) { const float d = Vector3Distance(a.pos.toFloat(), t.local); if (d < pdist) { pdist = d; tg = &t; } }
        const bool sameBody = tg != nullptr;
        const Vector3 tgLocal = tg ? tg->local : Vector3{};
        if (dist > PopTune::kReleaseRange && CavePopulation(a.species) && naturalSpawns) { // (#1633) far from every player: let go
            auto ch = spawnChunks.find(HomeChunk(a));
            if (ch != spawnChunks.end() && !ch->second.nestDestroyed) ch->second.filled = false; // the next visit fills it again
            gone.push_back((int)i);
            continue;
        }
        if (dist > 300.0f) continue;
        ++simulatedNow; a.simFrame = frameNo;
        float step = dt;
        if (dist > 80.0f) { if ((frameNo + (int)i) % 4) continue; step = dt * 4.0f; }
        else if (dist > 40.0f) { if ((frameNo + (int)i) % 2) continue; step = dt * 2.0f; } // (#1649) 40-80 m: every other update
        if (Combat1_On() && a.status.any()) { // (#1215) burn, bleed and poison over time (creatures on land: no water to douse)
            a.health -= Status_Step(a.status, step, false, ResistOf(a)); // (#1223: a creature's own resistances)
            if (a.health <= 0.0f) { Kill((int)i); --i; continue; }
        }
        if ((a.settleIn -= step) <= 0.0f && !a.clinging) { AccTimer accSettleT{ accSettle }; // (#1224: not while on a ceiling)
            a.settleIn = 1.0f;
            // Settle (#223): ground dug out from under it drops it onto what is left (up to 3 m); inside rock
            // (an older save, ground built over it) or over a deeper hole, it goes to the surface found from above.
            const Vector3 p0 = a.pos.toFloat();
            const Vector3 d = Vector3Normalize(p0);
            // Still standing on the ground (solid just under its feet, open 0.3 m over them): nothing to do. Two density
            // samples instead of a probe of 15 (and its fallback of ~180 on rough ground): the eroded field costs more a
            // sample, and every animal checked once a second (#476).
            const bool standing = body->densityAt(Vector3Subtract(p0, Vector3Scale(d, 0.15f))) >= ISO_LEVEL &&
                                  body->densityAt(Vector3Add(p0, Vector3Scale(d, 0.3f))) < ISO_LEVEL;
            Vector3 g; float ge;
            // (#1717 / #1227) A cave creature never takes the surface fallback (the field's whole band finds the ground at the TOP): a
            // miner on a sloped cave floor failed the standing test, was put on the surface over its cave and stood there by day,
            // its home straight below (nothing to walk to). No ground near it: back to its home.
            if (!standing) {
                if (Ground(*body, *field, d, g, ge, Vector3Length(p0), 0.6f, 3.0f)) a.pos = DVec3(g);
                else if (CavePopulation(a.species) && CaveSettleFix()) { if (a.hasHome) a.pos = DVec3(a.home); }
                else if (Ground(*body, *field, d, g, ge)) a.pos = DVec3(g);
            }
        }
        const Vector3 up = Vector3Normalize(a.pos.toFloat());
        Vector3 east, north;
        TangentFrame(up, east, north);
        // (#1217) Staggered: its attack is off and it stands where it is until the stagger ends; poise comes back after a rest.
        if (Combat3_On() && a.poise.max > 0.0f) {
            Poise_Step(a.poise, step);
            if (Poise_Staggered(a.poise)) { a.windUpAt = -1.0; a.moveSpeed = 0.0f; continue; } // (#1218) its state kept: it was set to Flee here and read "fleeing"
        }
        // Decide.
        float speed = 0.0f;
        Vector3 toward{};
        bool steer = false;
        const bool dark = Vector3DotProduct(sun, body->dirToWorld(up)) < -0.1f;
        // Night hunters: in daylight they run off, and are gone once out of sight (or after 40 s).
        if (s.night && !dark) {
            a.daylight += step;
            if (a.daylight > 40.0f || (pdist > 70.0f && a.daylight > 5.0f)) { gone.push_back((int)i); continue; }
            a.fleeTimer = std::max(a.fleeTimer, 1.0f);
            a.fleeFrom = sameBody ? tgLocal : a.pos.toFloat();
        }
        // Hostile ones (not on peaceful): thornbacks charge within 12 m (or for 20 s after being hit),
        // dusk stalkers hunt within 70 m after dark. A campfire within 10 m keeps the player safe, and
        // stalkers will not come within 14 m of one.
        bool hunting = false;
        bool goingTo = false, ranged = false; // (#1225) going for a pile (a thief); throwing from a distance (a thrower)
        // (#1223) A creature (a creatures.txt row) is driven by the brain (creatures.h): what it perceives (the nearest player,
        // creatures of factions it hunts or fears), whether it was hurt, noise it heard, fire, daylight. A plain animal, or
        // RAYLIBGAME_CREATURES=0: the code below it (the brain's thornback and stalker rows reproduce that code).
        const CreatureSpec* cs = Creature_Spec(s.id);
        if (cs && cs->split) OozeOre(a); // (#1226) an ooze takes the colour of the ore where it lives
        BrainOut bo;
        int preyIdx = -1;                      // another creature it goes for (infighting)
        Vector3 huntAt = tgLocal; float huntDist = pdist;
        if (cs) {
            PROFILE_SCOPE("creatures: brain"); AccTimer accBrainT{ accBrain }; // (#1649)
            struct WorstT { double& w; std::string& who; const Animal& a; const char* sp; std::chrono::steady_clock::time_point t0 = std::chrono::steady_clock::now();
                ~WorstT() { const double d = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(); if (d > w) { w = d; who = std::string(sp) + " mood " + Mood_Name(a.brain.mood) + (a.clinging ? " clinging" : ""); } } } worstT{ worstBrain, worstBrainWho, a, s.id.c_str() };
            const Mood before = a.brain.mood;
            const Vector3 selfPos = a.pos.toFloat(); // (#1698) where it stands, in floats for the senses
            // (#1650 slice 2) AI tiers: engaged (near a player, winding up, attacking, just hit) decide every update as before;
            // near ones every kNearPeriod, far ones every kFarPeriod, reusing the last decision between (its one-shot events
            // cleared). The motor below still runs at the creature's step rate.
            const bool engaged = a.lastDecideAt < 0.0 || pdist < AiTier::kEngagedRange || a.brain.windUpAt >= 0.0 || a.hurtPending ||
                                 a.brain.mood == Mood::Attack || !AiTier::On();
            const double period = engaged ? 0.0 : pdist < AiTier::kNearRange ? AiTier::kNearPeriod : AiTier::kFarPeriod;
            // (#1650) Each creature on its own phase (0.75-1.25 of the period, from its id): a crowd made together decided
            // together, all in one update (4-7 ms spikes) and nothing between.
            const double phase = 0.75 + 0.5 * (double)((a.netId * 0x9E3779B97F4A7C15ull) >> 40) / (double)(1ull << 24);
            const bool decide = simClock - a.lastDecideAt >= period * phase;
            BrainInput bi; // (empty between decisions: a hit makes it engaged, so bi.hurt is always fresh)
            if (decide) {
            const float brainDt = engaged ? step : (float)std::min(simClock - a.lastDecideAt, 2.0);
            bi.now = simClock;
            bi.health01 = s.health > 0.0f ? a.health / s.health : 1.0f;
            bi.daylight = !dark;
            bi.peaceful = !threatsOn || Difficulty_ThreatScale() <= 0.0f;
            bi.self = selfPos;
            Vector3 fire{};
            bi.nearFire = cs->fireScares && NearFire(selfPos, 14.0f, fire);
            bi.fire = fire;
            if (cs->wakeByLight) { Vector3 lf; bi.lightNear = NearFire(selfPos, cs->wakeRange, lf) || (sameBody && tg->lit && pdist < cs->wakeRange); } // (#1227)
            if (sameBody) {
                Percept p;
                p.id = kPlayerPercept | tg->session; p.pos = tgLocal; p.dist = pdist; p.faction = Faction_Player();
                p.seen = Sense_Sees(*cs, pdist, dark && !tg->lit ? 0.0f : 1.0f, true);
                Vector3 f2; p.safe = cs->fireSafe && NearFire(tgLocal, 10.0f, f2); // (#1225: not for torch-bearers)
                p.vulnerable = tg->vulnerable;
                bi.others.push_back(p);
            }
            AiGrid::Near(a, [&](int j) { // other creatures it has a stance toward, within its sight (#1650: the 40 m cells around it)
                if (j == (int)i || animals[j].body != a.body) return;
                const Vector3 jp = animals[j].pos.toFloat();
                const float d = Vector3Distance(jp, selfPos);
                if (d > std::min(40.0f, cs->sight)) return;
                const int fj = AnimalFaction(animals[j]);
                if (Factions().stance(cs->factionIdx, fj) == Stance::Neutral) return;
                Percept p; p.id = animals[j].netId; p.pos = jp; p.dist = d; p.faction = fj; p.seen = Sense_Sees(*cs, d, dark ? 0.0f : 1.0f, true);
                bi.others.push_back(p);
            });
            if (a.hurtPending) { bi.hurt = true; bi.hurtBy = a.hurtBy; bi.hurtFrom = (a.hurtBy & kPlayerPercept) ? tgLocal : a.hurtFrom; a.hurtPending = false; }
            for (const NoiseEvent& n : Noise_Recent(simClock))
                if (n.body == body->id && n.source != a.netId && Sense_Hears(*cs, Vector3Distance(n.pos, selfPos), n.loudness)) {
                    bi.heardNoise = true; bi.noiseAt = n.pos;
                    if (n.source & kPlayerPercept) a.heardPlayerAt = simClock; // (#1224) the player made it
                }
            for (Percept& p : bi.others) if (p.id & kPlayerPercept) p.heard = simClock - a.heardPlayerAt < 4.0;
            bo = Brain_Step(a.brain, *cs, bi, brainDt);
            a.lastBo = bo; a.lastDecideAt = simClock;
            } else { bo = a.lastBo; bo.windUpStarted = bo.strike = false; ++gAiReused; } // (#1650) between decisions
            if (bo.mood == Mood::Idle && Rand01(a.rng) < step / 12.0f) FxEvent(a, *cs, CreatureEvent::Idle, body->id, selfPos); // (#1363) an idle call now and then (~12 s)
            if (getenv("RAYLIBGAME_CREATURE_LOG") && (frameNo % 60) == 0) { const Percept* pp = nullptr; for (const Percept& q : bi.others) if (q.id & kPlayerPercept) pp = &q;
                printf("creature %s %llu: mood %s, player %.1f m seen %d heard %d safe %d, flee %.1f, health %.0f\n", s.id.c_str(), (unsigned long long)a.netId, Mood_Name(bo.mood), pp ? pp->dist : -1.0f, pp ? (int)pp->seen : -1, pp ? (int)pp->heard : -1, pp ? (int)pp->safe : -1, a.fleeTimer, a.health); }
            if ((bo.mood == Mood::Alert || bo.mood == Mood::Chase) && (before == Mood::Idle || before == Mood::Sleep || before == Mood::Wander))
                FxEvent(a, *cs, CreatureEvent::Alert, body->id, selfPos);
            if (bo.mood == Mood::Flee || bo.mood == Mood::Retreat) { a.fleeTimer = std::max(a.fleeTimer, a.brain.timer); a.fleeFrom = bo.goal; }
            if ((bo.mood == Mood::Chase || bo.mood == Mood::Attack) && a.fleeTimer <= 0.0f) {
                if (bo.target & kPlayerPercept) hunting = sameBody;
                else if ((preyIdx = IndexOfNet(bo.target)) >= 0) { hunting = true; huntAt = animals[preyIdx].pos.toFloat(); huntDist = Vector3Distance(huntAt, selfPos); }
            }
            if (bo.mood == Mood::Alert && a.fleeTimer <= 0.0f) { // goes to look where the noise was
                a.state = State::Walk;
                toward = Vector3Subtract(bo.goal, selfPos); steer = true; speed = s.walk;
            }
            // (#1224) Ceiling clingers (the tunnel lurker): idle in a cave with a ceiling 2.2-20 m over it (chambers run tall), it climbs up and waits
            // there. Prey heard or seen right under it (within dropRange sideways): a rattle (the telegraph), and 0.35 s later it
            // drops onto them (damage x dropScale, a heavy hit). Prey it hunts but that stays clear for 3 s: it drops anyway
            // and gives chase. Hit or staggered: it falls.
            if (cs->ceiling) {
                const Vector3 aup = Vector3Normalize(selfPos);
                if (getenv("RAYLIBGAME_CREATURE_LOG") && (frameNo % 60) == 0) { float r = -1.0f; for (float h = 0.5f; h <= 20.0f; h += 0.25f) if (body->densityAt(Vector3Add(selfPos, Vector3Scale(aup, h))) >= ISO_LEVEL) { r = h; break; }
                    printf("creature %llu: mood %s, hunting %d, flee %.1f, clinging %d, since bite %.1f, roof %.2f m\n", (unsigned long long)a.netId, Mood_Name(bo.mood), (int)hunting, a.fleeTimer, (int)a.clinging, simClock - a.lastBite, r); }
                if (a.clinging) {
                    // (#1224 playtest: "walking under/near doesn't trigger them to drop, no audio cue") It never heard the
                    // player: it hears from the ceiling (up to 20 m over the floor) and a quiet walk carries a few metres, so
                    // the drop (and its rattle) never started. Prey right under it (within its drop range sideways, 3 m of the
                    // floor) now sets it off as well: it feels the steps.
                    bool preyUnder = false;
                    if (sameBody && preyIdx < 0) {
                        const Vector3 d = Vector3Subtract(tgLocal, a.clingFloor);
                        const float vert = Vector3DotProduct(d, aup);
                        preyUnder = Vector3Length(Vector3Subtract(d, Vector3Scale(aup, vert))) < cs->dropRange && fabsf(vert) < 3.0f;
                    }
                    const bool knocked = bi.hurt || (Combat3_On() && Poise_Staggered(a.poise));
                    auto fall = [&]() { a.pos = a.prevPos = DVec3(a.clingFloor); a.clinging = false; a.dropIn = -1.0f; a.clingWait = 0.0f; };
                    if (knocked) fall();
                    else if (a.dropIn >= 0.0f) {
                        a.dropIn -= step;
                        if (a.dropIn <= 0.0f) {
                            const Vector3 below = a.clingFloor;
                            const bool onSteps = a.dropOnSteps; a.dropOnSteps = false;
                            fall();
                            if (sameBody && preyIdx < 0 && (hunting || onSteps)) { // (#1224) or set off by steps under it
                                // It drops onto its prey: within its drop range sideways of where it hung (and not far above or
                                // below that floor), it lands beside them and the drop hits.
                                const Vector3 d = Vector3Subtract(tgLocal, below);
                                const float vert = Vector3DotProduct(d, aup);
                                const Vector3 flat = Vector3Subtract(d, Vector3Scale(aup, vert));
                                const float side = Vector3Length(flat);
                                if (getenv("RAYLIBGAME_CREATURE_LOG")) printf("creature %llu: drop, prey %.2f m sideways, %.2f m up from the floor below\n", (unsigned long long)a.netId, side, vert);
                                if (side < cs->dropRange + 0.5f && fabsf(vert) < 3.0f) {
                                    if (side > 0.9f) a.pos = a.prevPos = DVec3(below) + DVec3(flat) * (double)((side - 0.9f) / side); // beside them
                                    DamageEvent e;
                                    e.amount = cs->damage * cs->dropScale * Difficulty_ThreatScale();
                                    e.poise = e.amount * PoiseTune::kAttackPoisePerDamage;
                                    e.heavy = true; e.kind = DamageKind::Melee; e.type = cs->type;
                                    e.hasSource = true; e.source = body->toWorldD(DVec3(below)); e.knockback = cs->knockback; e.sourceId = a.netId;
                                    e.faction = (uint8_t)std::max(0, cs->factionIdx);
                                    if (tg->session != 0) { if (gAnimalDamageRemote) gAnimalDamageRemote(tg->session, e); }
                                    else {
                                        const DamageResult r = Combat_Damage(LocalSession(), e);
                                        GameUI_Toast(r.dodged ? TextFormat("The %s drops - you roll clear", s.name.c_str()) : TextFormat("The %s drops on you!", s.name.c_str()));
                                    }
                                    ++lurkerDrops;
                                    a.lastBite = simClock;
                                    a.brain.lastStrike = simClock; // its cooldown starts with the drop
                                }
                            }
                        }
                    } else if ((hunting && preyIdx < 0) || preyUnder) {
                        const Vector3 d = Vector3Subtract(tgLocal, a.clingFloor);
                        const float side = Vector3Length(Vector3Subtract(d, Vector3Scale(aup, Vector3DotProduct(d, aup))));
                        a.clingWait += step;
                        if (side < cs->dropRange || (hunting && a.clingWait > 3.0f)) {
                            a.dropIn = 0.35f;
                            a.dropOnSteps = preyUnder && !hunting;
                            FxEvent(a, *cs, CreatureEvent::Telegraph, body->id, selfPos);
                            if (tg->session == 0 && side < cs->dropRange) { Audio_Play(Sfx::Error, 0.4f, 1.8f); GameUI_Toast("Something rattles overhead!"); }
                        }
                    }
                    if (a.clinging) { a.moveSpeed = 0.0f; a.state = State::Rest; continue; } // hanging on: no walking
                } else if (!hunting && a.fleeTimer <= 0.0f && (bo.mood == Mood::Idle || bo.mood == Mood::Sleep) && simClock - a.lastBite > 5.0 && simClock >= a.roofCheckAt) {
                    // (#1650) at most about once a second (its own phase): a dozen idle lurkers searching 79 samples each in one
                    // update was a 3-5 ms spike in a crowd
                    a.roofCheckAt = simClock + 0.75 + 0.5 * (double)((a.netId * 0x9E3779B97F4A7C15ull) >> 40) / (double)(1ull << 24);
                    float up = 0.0f; bool roof = false;
                    auto findRoof = [&] { for (float h = 0.5f; h <= 20.0f; h += 0.25f) if (body->densityAt(Vector3Add(selfPos, Vector3Scale(aup, h))) >= ISO_LEVEL) { up = h; roof = true; break; } return roof; };
                    // (#1649) "no roof here" is cached per cell (in the open it searched 79 samples every step); a roof is found exactly
                    if (ColumnCache::Get(*body, selfPos, 1, findRoof) && !roof) findRoof();
                    if (roof && up >= 2.2f) {
                        a.clingFloor = selfPos;
                        a.pos = a.prevPos = a.pos + DVec3(aup) * (double)(up - 0.3f);
                        a.clinging = true; a.clingWait = 0.0f; a.dropIn = -1.0f;
                        a.moveSpeed = 0.0f; a.state = State::Rest;
                        continue;
                    }
                }
            }
            // (#1225) Thieves: a dropped pile within 12 m (and nothing carried yet): it goes for it, takes it whole and runs
            // (12 s away from the player). It carries it until it falls (then it drops with its loot).
            if (cs->thief && a.carried.empty() && a.fleeTimer <= 0.0f) {
                Vector3 pile{};
                const Vector3 thiefAt = a.pos.toFloat();
                if (Combat_PileNear(a.body, thiefAt, 12.0f, pile)) {
                    if (Vector3Distance(pile, thiefAt) < 1.4f) {
                        if (Combat_TakePileNear(a.body, thiefAt, 1.6f, a.carried) && !a.carried.empty()) {
                            a.fleeTimer = 12.0f; a.fleeFrom = sameBody ? tgLocal : thiefAt; hunting = false;
                            ++thefts;
                            if (sameBody && tg->session == 0 && pdist < 40.0f)
                                GameUI_Toast(TextFormat("A %s grabs your %s and runs!", s.name.c_str(), a.carried[0].id.c_str()));
                        }
                    } else { hunting = true; goingTo = true; preyIdx = -1; huntAt = pile; huntDist = Vector3Distance(pile, thiefAt); }
                }
            }
            // (#1225) The group's thrower keeps its distance from the player it hunts and throws instead of closing in.
            ranged = hunting && !goingTo && preyIdx < 0 && ((a.thrower && cs->thrower) || cs->ranged); // (#1228: every sentinel shoots)
        } else if (s.hostile && threatsOn && Difficulty_ThreatScale() > 0.0f && sameBody && a.fleeTimer <= 0.0f && tg->vulnerable) { // not inside a cockpit
            Vector3 fire;
            const bool playerSafe = NearFire(tgLocal, 10.0f, fire);
            if (s.night && NearFire(a.pos.toFloat(), 14.0f, fire)) { a.fleeTimer = 3.0f; a.fleeFrom = fire; }
            else if (s.night) hunting = dark && pdist < 70.0f && !playerSafe;
            else hunting = (pdist < 12.0f || now < a.angryUntil) && pdist < 30.0f && !playerSafe;
        }
        if (!hunting) a.windUpAt = -1.0; // a wind-up broken off (the prey gone or safe) never lands later
        const Vector3 curPos = a.pos.toFloat(); // where it stands for the rest of the step (a.pos changes only when it moves, at the end)
        if (hunting) {
            a.state = State::Hunt;
            toward = Vector3Subtract(huntAt, a.pos.toFloat());
            steer = true;
            speed = s.run * (s.night ? 0.8f : 1.0f);
            if (cs && cs->burst) speed *= fmodf((float)simClock + (a.rng % 97) * 0.01f, cs->burstOn + cs->burstOff) < cs->burstOn ? 1.0f : 0.33f; // (#1224) bursts
            const float reach = cs ? cs->reach : 0.5f * s.size.x + 1.0f;
            const bool atPlayer = preyIdx < 0;
            // A telegraphed attack (#165): in reach it winds up (a growl, a warning) and strikes kBiteWindUp later; the
            // hit goes through Combat_Damage (a shield block, knockback away from it). Out of reach by then: it missed.
            bool windUp = false, strike = false;
            if (cs) { // (#1223) the brain's timing (the same numbers for the surface hunters)
                if (bo.mood == Mood::Attack) speed = 0.0f;
                windUp = bo.windUpStarted; strike = bo.strike;
            } else {
                if (pdist < reach) { speed = 0.0f; windUp = a.windUpAt < 0.0 && simClock - a.lastBite > 1.2; } // (on simClock: tests step animals many times in a frame)
                strike = a.windUpAt >= 0.0 && simClock - a.windUpAt >= kBiteWindUp;
            }
            // (#1225) A thief going for a pile walks to it (no attack). A thrower: too close (under throwMin), it backs off; in
            // range, it stands and throws (a telegraphed wind-up, then a stone on an arc at the player's chest); further, it
            // closes in. Neither melees.
            if (goingTo) { a.state = State::Walk; windUp = strike = false; }
            else if (ranged) {
                windUp = strike = false;
                if (pdist < cs->throwMin) { toward = Vector3Subtract(a.pos.toFloat(), tgLocal); speed = s.walk * 1.5f; a.throwWindAt = -1.0; }
                else if (pdist <= cs->throwRange) {
                    speed = 0.0f;
                    if (a.throwWindAt < 0.0 && simClock - a.lastThrow > cs->throwCooldown) {
                        a.throwWindAt = simClock;
                        FxEvent(a, *cs, CreatureEvent::Telegraph, body->id, a.pos.toFloat());
                        if (tg->session == 0) GameUI_Toast(TextFormat(cs->dart ? "The %s charges up!" : "The %s winds up a throw!", s.name.c_str())); // (#1228: a sentinel's charge)
                    }
                    if (a.throwWindAt >= 0.0 && simClock - a.throwWindAt >= cs->throwWindUp) {
                        const Vector3 throwAt = a.pos.toFloat();
                        const Vector3 aup = Vector3Normalize(throwAt), tup = Vector3Normalize(tgLocal);
                        ThrowAt(a.body, Vector3Add(throwAt, Vector3Scale(aup, cs->dart ? 1.1f : 1.6f)), Vector3Add(tgLocal, Vector3Scale(tup, 1.1f)), cs->throwDamage, a.netId, false, cs->dart ? cs->dartSpeed : 16.0f, cs->dart);
                        a.throwWindAt = -1.0; a.lastThrow = simClock;
                    }
                } else a.throwWindAt = -1.0;
            }
            if (windUp) {
                a.windUpAt = simClock; a.windUpHeading = a.heading; // (#1227) it faces this way until the strike
                if (cs) FxEvent(a, *cs, CreatureEvent::Telegraph, body->id, curPos);
                if (atPlayer && tg->session == 0) { Audio_Play(Sfx::Error, 0.5f, 0.35f);
                const DVec3 at = body->toWorldD(DVec3(a.pos)); // (#1635) a caption, not a toast: "<name> raises its pick", with its direction
                Subtitles_Event(s.name + " " + (cs && !cs->sayWindup.empty() ? cs->sayWindup : std::string(s.night ? "lunges" : "lowers its head")), true, &at, Audio_GetListener()); } // (a remote player sees the wind-up in the animal's flags)
            }
            if (strike) {
                a.windUpAt = -1.0;
                a.lastBite = simClock;
                // (#1227) A creature's strike lands only in front of it (within 75 degrees of where it faces, GUESS): with its
                // facing held through the wind-up, stepping round it is the counter.
                bool inFront = true;
                if (cs && (CreatureTurnMask() & 4)) {
                    const Vector3 upA = Vector3Normalize(curPos);
                    Vector3 eA, nA; TangentFrame(upA, eA, nA);
                    const Vector3 fwdA = Vector3Add(Vector3Scale(eA, cosf(a.heading)), Vector3Scale(nA, sinf(a.heading)));
                    Vector3 toT = Vector3Subtract(huntAt, curPos);
                    toT = Vector3Subtract(toT, Vector3Scale(upA, Vector3DotProduct(toT, upA)));
                    if (Vector3Length(toT) > 0.3f) inFront = Vector3DotProduct(Vector3Normalize(toT), fwdA) >= cosf(75.0f * DEG2RAD);
                }
                if (huntDist < reach + 0.4f && inFront) {
                    DamageEvent e;
                    e.amount = (cs ? cs->damage : 8.0f) * Difficulty_ThreatScale();
                    e.poise = cs && cs->poise > 0.0f ? cs->poise : e.amount * PoiseTune::kAttackPoisePerDamage; // (#1217) a strong bite or ram: 20
                    e.kind = (cs ? cs->bite : s.night) ? DamageKind::Bite : DamageKind::Melee;
                    if (cs) e.type = cs->type;
                    e.hasSource = true;
                    e.source = body->toWorldD(a.pos);
                    e.knockback = cs ? cs->knockback : s.night ? 3.0f : 6.0f; // a ram shoves harder than a bite
                    e.sourceId = a.netId;
                    if (cs) e.faction = (uint8_t)std::max(0, cs->factionIdx);
                    if (!atPlayer) creatureHits.push_back({ animals[preyIdx].netId, a.netId, curPos, e }); // (#1223) applied after the loop
                    else if (tg->session != 0) { // (#334) a remote player: its client applies the hit (shield, health), the server the shove
                        if (gAnimalDamageRemote) gAnimalDamageRemote(tg->session, e);
                    } else {
                        const DamageResult r = Combat_Damage(LocalSession(), e);
                        if (cs && cs->corrode > 0 && !r.dodged && !r.parried && ToolOf(terrainTool.heldTool).type != ToolType::Hands) { // (#1226) its touch eats the tool
                            ToolWear_OnStroke(terrainTool.heldTool, cs->corrode);
                            GameUI_Toast(TextFormat("The %s corrodes your %s", s.name.c_str(), terrainTool.heldTool.c_str()));
                        }
                        if (r.parried) { // (#1218) a parry staggers the attacker (2 s, #1217)
                            if (a.poise.max <= 0.0f) Poise_Init(a.poise, Poise_CreatureMax(s.health));
                            Poise_Stagger(a.poise, PoiseTune::kParryStaggerTime);
                        }
                        Audio_Play(Sfx::Error, 0.7f, r.blocked ? 0.9f : 0.5f);
                        if (r.parried || r.dodged || r.blocked) // what you did about it stays a toast
                            GameUI_Toast(r.parried ? TextFormat("Parried! The %s staggers", s.name.c_str()) : r.dodged ? TextFormat("You dodge the %s", s.name.c_str()) : TextFormat("You block the %s", s.name.c_str()));
                        else { // (#1635) the hit: a caption with its direction, not a toast
                            const DVec3 at = body->toWorldD(DVec3(a.pos));
                            Subtitles_Event(s.name + " " + (cs && !cs->sayHit.empty() ? cs->sayHit : std::string(s.night ? "bites you" : "hits you")), true, &at, Audio_GetListener());
                        }
                    }
                } else if (atPlayer && tg->session == 0) { const DVec3 at = body->toWorldD(DVec3(a.pos)); Subtitles_Event(s.name + " misses", false, &at, Audio_GetListener()); } // (#1635)
            }
        } else if (a.fleeTimer > 0.0f) {
            a.fleeTimer -= step;
            a.state = State::Flee;
            if (sameBody) a.fleeFrom = tgLocal;
            toward = Vector3Subtract(curPos, a.fleeFrom);
            steer = true;
            speed = s.run;
        } else if (sameBody && pdist > 2.5f && (now < a.followUntil || (pdist < 16.0f && !tg->held.empty() && std::find(s.diet.begin(), s.diet.end(), tg->held) != s.diet.end()))) {
            // Fed (30 s), or the player holds its food: it follows (wild ones too: the food wins).
            a.state = State::Follow;
            toward = Vector3Subtract(tgLocal, curPos);
            steer = true;
            speed = std::min(s.run, s.walk * 2.0f);
        } else if (s.wild && sameBody && pdist < 25.0f) {
            a.fleeTimer = 4.0f;
            a.fleeFrom = tgLocal;
            continue;
        } else if (s.hostile && a.chase >= 0 && a.chase < (int)animals.size() && animals[a.chase].body == a.body &&
                   Vector3Distance(animals[a.chase].pos.toFloat(), curPos) < 35.0f && !(s.night && !dark)) {
            // A predator after prey (Plan 5 A; chase only): it runs at it; close in, the prey bolts.
            a.state = State::Hunt;
            const Vector3 preyPos = animals[a.chase].pos.toFloat();
            toward = Vector3Subtract(preyPos, curPos);
            steer = true;
            speed = s.run * 0.75f;
            Animal& prey = animals[a.chase];
            if (Vector3Distance(preyPos, curPos) < 6.0f) { prey.fleeTimer = std::max(prey.fleeTimer, 3.0f); prey.fleeFrom = curPos; }
        } else if (dark || (sameBody && Weather_StormHere())) { // night, or a storm: they shelter and rest
            a.state = State::Rest; // night
            a.fatigue = std::max(0.0f, a.fatigue - step / 90.0f);
        } else {
            // Needs (Plan 5 A): hunger grows (grazing eats it down), so does tiredness (resting sleeps it off);
            // when the current activity runs out, the strongest pull wins - graze, wander (more so when the
            // herd has drifted off: it heads back to it), or lie down for a rest.
            a.hunger = std::min(1.0f, a.hunger + step / 240.0f);
            a.fatigue = std::min(1.0f, a.fatigue + step / (a.state == State::Walk ? 300.0f : 600.0f));
            if (a.state == State::Graze) a.hunger = std::max(0.0f, a.hunger - step / 25.0f);
            if (a.state == State::Rest) a.fatigue = std::max(0.0f, a.fatigue - step / 40.0f);
            const float herdAway = a.herdSize > 0 ? Vector3Distance(a.herdCentre, curPos) : 0.0f;
            a.timer -= step;
            if (a.timer <= 0.0f) {
                const float graze = a.hunger * 1.2f + 0.1f * Rand01(a.rng);
                const float walk = 0.35f + (herdAway > 12.0f ? 0.6f : 0.0f) + 0.25f * Rand01(a.rng);
                const float rest = a.fatigue * 1.1f + 0.05f * Rand01(a.rng);
                a.state = graze >= walk && graze >= rest ? State::Graze : walk >= rest ? State::Walk : State::Rest;
                if (a.state == State::Walk) a.heading += (Rand01(a.rng) - 0.5f) * 2.2f;
                a.timer = a.state == State::Walk ? 2.0f + 5.0f * Rand01(a.rng) : a.state == State::Rest ? 8.0f + 10.0f * Rand01(a.rng) : 2.0f + 6.0f * Rand01(a.rng);
            }
            if (a.state == State::Walk) {
                speed = s.walk;
                if (herdAway > 12.0f) { toward = Vector3Subtract(a.herdCentre, curPos); steer = true; } // back to the herd
            }
        }
        if (steer) {
            const float x = Vector3DotProduct(toward, east), y = Vector3DotProduct(toward, north);
            if (x * x + y * y > 1e-6f) a.heading = atan2f(y, x);
        }
        // A planned route (#253): head for its next point while the goal stays put (within 6 m) and the route is fresh.
        // Fleeing, the goal is a point 25 m on away from the threat.
        a.planWait = std::max(0.0f, a.planWait - step);
        const Vector3 goalPt = !steer ? curPos : a.state == State::Flee && Vector3Length(toward) > 1e-3f
            ? Vector3Add(curPos, Vector3Scale(Vector3Normalize(toward), 25.0f)) : Vector3Add(curPos, toward);
        if (!a.path.empty()) {
            a.pathAge += step;
            if (!steer || speed <= 0.0f || a.pathAge > 10.0f || Vector3Distance(a.pathGoal, goalPt) > 6.0f) a.path.clear();
            while (!a.path.empty() && Vector3Distance(a.path.front(), curPos) < 0.7f) a.path.erase(a.path.begin());
            if (!a.path.empty()) {
                const Vector3 w = Vector3Subtract(a.path.front(), curPos);
                a.heading = atan2f(Vector3DotProduct(w, north), Vector3DotProduct(w, east));
                a.detourTime = 0.0f;
            }
        }
        // (#1227 playtest: "perfectly tracks player, so getting a hit from behind is very hard ... when it starts winding up an
        // attack, it should lock its looking angle") Creatures: a wound-up attack holds the facing it began with (step aside or
        // behind it then), and hunting they turn at most kCreatureTurnRate (GUESS), not straight at the player every step.
        // (#1227, user 2026-11-22) Cave species stay underground: they give up a chase at the cave mouth (a step out under the
        // open sky fails, tryStep below). Only the Lost Miner may come out near a mouth at night; by day, or a cave creature
        // caught outside, it walks back to where it spawned.
        PROFILE_SCOPE("creatures: move"); // (#1649) the leash, the steps, the path plan (to the end of this animal)
        bool leashed = false, outside = false;
        if (cs && cs->habitat == Habitat::Cave && a.hasHome) {
            leashed = true;
            if (s.id == "lost_miner" && dark && field) {
                std::vector<PlanetField::CaveEntrance> ents;
                field->caveEntrancesNear(Vector3SubtractValue(curPos, 30.0f), Vector3AddValue(curPos, 30.0f), ents);
                if (!ents.empty()) leashed = false; // out at night near a mouth: allowed
            }
            outside = leashed && OpenSkyOver(*body, curPos);
            if (outside) { // back under rock: toward its home, at a walk, the chase dropped
                Vector3 h = Vector3Subtract(a.home, curPos);
                h = Vector3Subtract(h, Vector3Scale(up, Vector3DotProduct(h, up)));
                if (Vector3Length(h) > 0.5f) {
                    a.heading = atan2f(Vector3DotProduct(h, north), Vector3DotProduct(h, east));
                    a.windUpAt = -1.0; a.state = State::Walk; steer = true; speed = std::max(speed, s.walk);
                } else if (Vector3Distance(a.home, curPos) > 3.0f) { // (#1717) its home straight below (stranded over its cave by
                    a.pos = a.prevPos = DVec3(a.home);                  // the old surface fallback): back down to it
                    a.windUpAt = -1.0; a.moveSpeed = 0.0f;
                }
            }
        }
        // (#1717) Close in on its prey (the turn-rate rule's case): the detours and the stuck rule turned it up to ~180 degrees at
        // once, past the rule; it then turned back at 180 deg/s, met the same obstacle and flipped again (the user: "turning
        // slowly one way, flipping 180, turning back"). Close in it only sidesteps (up to ~57 degrees) or stands facing its prey.
        static const bool closeHuntOn = !(getenv("RAYLIBGAME_CLOSE_HUNT") && getenv("RAYLIBGAME_CLOSE_HUNT")[0] == '0'); // A/B: 0 the old turns
        const bool closeHunt = closeHuntOn && cs && CreatureTurnRules() && (CreatureTurnMask() & 2) && (a.state == State::Hunt || bo.mood == Mood::Alert) && !goingTo && huntDist < 6.0f; // (#1717: Alert too: it snapped to the noise, then turned back)
        if (cs && CreatureTurnRules()) {
            if ((CreatureTurnMask() & 1) && a.windUpAt >= 0.0 && a.state == State::Hunt) a.heading = a.windUpHeading; // (an attack under way: not a thief turned to a pile)
            else if ((CreatureTurnMask() & 2) && (a.state == State::Hunt || bo.mood == Mood::Alert) && !goingTo && huntDist < 6.0f) { // (#1717: Alert too) // (close in, where getting behind it matters; beyond, it paths and turns freely: a thief at a distance missed its pile)
                constexpr float kCreatureTurnRate = 3.14159265f; // rad/s: 180 degrees a second
                float d = a.heading - a.prevHeading;
                d = atan2f(sinf(d), cosf(d));
                const float maxTurn = kCreatureTurnRate * step;
                if (fabsf(d) > maxTurn) a.heading = a.prevHeading + (d > 0.0f ? maxTurn : -maxTurn);
            }
        }
        if (Combat1_On()) speed *= Status_SpeedFactor(a.status); // (#1215) chilled 0.7, frozen 0
        if (speed <= 0.0f) { a.moveSpeed = 0.0f; continue; }
        // Move along the surface; turn away from water, cliffs and build parts.
        const float stepLen = speed * step * (0.6f + 0.4f * a.age);
        const float height = s.size.y * (0.6f + 0.4f * a.age);
        // The ground ahead within a 40 degree climb or drop of where it stands (a steeper step is a wall or a
        // cliff), with room for its body above (no walking under overhangs into rock).
        // A small ledge (up to 0.35 m, less for small kinds) can be stepped up or hopped down whatever the stride:
        // the 0.25 m terrain lattice and digs make steps the old height map never had.
        const float reach = std::max(stepLen * 0.84f + 0.05f, std::min(0.35f, height * 0.45f));
        // Going somewhere (chasing, fleeing, following, back to the herd) it leaps ledges up to 1.3x its height (1.2 m at
        // most) and hops down 1.5 m; grazing and wandering it keeps to small steps. Taller cliffs still stop it.
        const float reachUp = steer ? std::max(reach, std::min(1.2f, height * 1.3f)) : reach;
        const float reachDown = steer ? std::max(reach, 1.5f) : reach;
        // Heading somewhere (the player, the herd, away from a threat) and blocked: go round (#223: real terrain has
        // ledges, pits and basins the height map hid). The detour is held for a while - following the obstacle on one
        // side, turning as far as backing off - before aiming at the target again, so a pocket or a long ledge
        // doesn't trap it pushing straight into the same edge.
        if (closeHunt) a.detourTime = 0.0f; // (#1717) no held detour close in: it turned it past the turn-rate rule
        else if (steer && a.detourTime > 0.0f) { a.detourTime -= step; a.heading = a.detourHeading; }
        else if (a.detourTime > 0.0f) a.detourTime = 0.0f;
        // (#1650 slice 3) The step's probes run after this loop for every mover at once (pass 2), and what follows them in
        // the loop's order (pass 3, below the loop).
        MoveReq mr;
        mr.net = a.netId; mr.body = body; mr.field = field; mr.east = east; mr.north = north; mr.curPos = curPos; mr.pos = a.pos;
        mr.stepLen = stepLen; mr.height = height; mr.reachUp = reachUp; mr.reachDown = reachDown; mr.speed = speed; mr.step = step;
        mr.heading = a.heading; mr.detourSide = a.detourSide; mr.sizeX = s.size.x; mr.leashed = leashed; mr.outside = outside;
        mr.steer = steer; mr.goalPt = goalPt; mr.rng = a.rng; mr.closeHunt = closeHunt;
        moveReqs.push_back(mr);
    }
    const auto tProbe0 = std::chrono::steady_clock::now(); // (#1650)
    // (#1650 slice 3) Pass 2: the probes, in parallel.
    std::vector<MoveRes> moveRes(moveReqs.size());
    { PROFILE_SCOPE("creatures: probes"); AiParallelFor((int)moveReqs.size(), [&](int k) { ProbeMove(moveReqs[k], moveRes[k]); }); }
    // Pass 3: applied in the loop's order on this thread (the path search is one shared search).
    std::unordered_map<uint64_t, int> indexOfNet;
    if (!moveReqs.empty()) { indexOfNet.reserve(animals.size()); for (size_t k = 0; k < animals.size(); ++k) indexOfNet[animals[k].netId] = (int)k; }
    const auto tApply0 = std::chrono::steady_clock::now(); // (#1650)
    PROFILE_SCOPE("creatures: apply");
    for (size_t k = 0; k < moveReqs.size(); ++k) {
        const MoveReq& mr = moveReqs[k];
        const MoveRes& res = moveRes[k];
        auto itI = indexOfNet.find(mr.net);
        if (itI == indexOfNet.end()) continue;
        const int i = itI->second;
        Animal& a = animals[i];
        const PlanetField* field = mr.field;
        const Planetoid* body = mr.body;
        const Vector3 curPos = mr.curPos, east = mr.east, north = mr.north, goalPt = mr.goalPt;
        const bool steer = mr.steer;
        const float step = mr.step, speed = mr.speed, stepLen = mr.stepLen, height = mr.height, reachUp = mr.reachUp, reachDown = mr.reachDown;
        Vector3 at = res.at0;
        DVec3 nextPos = res.next0;
        bool headingMoved = false; // a plan turned it: the detour from pass 2 was for the old heading
        bool ok = res.ok0;
        if (ok && a.detourTime <= 0.0f) a.detourSide = 0.0f;
        if (!ok && steer && !a.path.empty()) { a.path.clear(); a.planWait = 1.0f; } // the route's step failed (a thin wall between cells)
        // Blocked on the way somewhere: plan a route round (one plan a frame for all animals, then a 3 s wait after a
        // failed one); the local detour below covers the time between.
        // One search at a time, resumed across updates: its owner continues it (standing while it thinks, a few frames);
        // a search nobody has touched for 30 updates (its animal got going, or was removed) is dropped.
        if (gSearch.active && frameNo - gSearch.lastFrame > 30) gSearch.active = false;
        const bool mySearch = gSearch.active && gSearch.owner == (int)i && Vector3Distance(gSearch.from, curPos) < 3.0f;
        if (!ok && steer && a.path.empty() && (mySearch || (a.planWait <= 0.0f && !gSearch.active)) && plansLeft > 0) {
            --plansLeft; PROFILE_SCOPE("creatures: path"); // (#1649)
            if (!mySearch) PathBegin(*body, *field, curPos, goalPt, height, std::max(reachUp, 0.35f), std::max(reachDown, 0.35f), (int)i, frameNo);
            std::vector<Vector3> route;
            if (!PathStep(route, frameNo)) { a.moveSpeed = 0.0f; continue; } // still searching: wait for it
            if (!route.empty()) {
                a.path = std::move(route);
                a.pathGoal = goalPt; a.pathAge = 0.0f; a.planWait = 1.0f;
                const Vector3 w = Vector3Subtract(a.path.front(), curPos);
                a.heading = atan2f(Vector3DotProduct(w, north), Vector3DotProduct(w, east));
                a.detourTime = 0.0f;
                MoveReq turned = mr; turned.heading = a.heading;
                ok = ProbeStep(turned, a.heading, at, nextPos); // (one creature a frame at most: on this thread)
                headingMoved = true;
            } else a.planWait = 3.0f;
        }
        if (!ok && steer) {
            if (!headingMoved && res.detourTried) { // the detour pass 2 probed, for this heading and side
                a.detourSide = res.detourSideInit; a.rng = res.rngAfter;
                if (res.detourOk) {
                    ok = true; at = res.atD; nextPos = res.nextD;
                    a.heading += res.detourSide * res.detourTurn;
                    a.detourSide = res.detourSide;
                    a.detourHeading = a.heading;
                    a.detourTime = 1.5f;
                }
            } else { // a plan turned it this frame: the detour from the new heading, here
                if (a.detourSide == 0.0f) a.detourSide = Rand01(a.rng) < 0.5f ? 1.0f : -1.0f;
                MoveReq turned = mr; turned.heading = a.heading;
                for (float turn : { 0.5f, 1.0f, 1.6f, 2.2f, 3.1f }) if (!mr.closeHunt || turn < 1.05f) // (#1717) sidesteps only close in
                    for (float side : { a.detourSide, -a.detourSide })
                        if (!ok && ProbeStep(turned, a.heading + side * turn, at, nextPos)) {
                            ok = true;
                            a.heading += side * turn;
                            a.detourSide = side;
                            a.detourHeading = a.heading;
                            a.detourTime = 1.5f;
                        }
            }
        }
        if (ok) {
            a.pos = nextPos;
            a.moveSpeed = speed;
            a.gait += stepLen / std::max(0.25f, mr.sizeX * 0.45f) * PI; // one stride per ~0.9 body lengths
            a.stuckTime = 0.0f;
        } else {
            // Blocked (water, a cliff, a wall): turn away; stuck a while, turn harder and pick a new activity.
            a.moveSpeed = 0.0f;
            a.stuckTime += step;
            if (mr.closeHunt && a.stuckTime <= 4.0f) continue; // (#1717) close in on its prey: it stands facing it, no flip
            a.heading += PI * (0.6f + 0.6f * Rand01(a.rng)) * (a.stuckTime > 2.0f ? 1.5f : 1.0f);
            a.timer = 1.0f + Rand01(a.rng);
            if (a.stuckTime > 4.0f) { a.state = State::Graze; a.stuckTime = 0.0f; }
        }
        // (#1717) Close in on its prey, its facing keeps the turn rate through this pass too: a sidestep or a planned route turned
        // it up to ~75 degrees in one update. It may step sideways, but it turns at 180 deg/s.
        if (mr.closeHunt) {
            const float d = atan2f(sinf(a.heading - mr.heading), cosf(a.heading - mr.heading)), m = PI * step;
            if (fabsf(d) > m) a.heading = mr.heading + (d > 0.0f ? m : -m);
        }
    }
    const auto tAfter0 = std::chrono::steady_clock::now(); // (#1650)
    // Penned: build parts within 30 m in all eight directions around an animal near the player
    // (checked every 2 s; the quest's "pen two animals").
    // (#1650) A pass over the animals starts every 2 s and is spread over the next updates, 8 animals each: at once it cost
    // 6-7 ms in a crowd (eight 30 m walls probed for each of 200 animals): the creature step's p99.
    {
        static size_t penCursor = 0; static int penCount = 0; static bool penPass = false;
        if (here && !penPass && SimTime() - lastPenCheck > 2.0) { lastPenCheck = SimTime(); penPass = true; penCursor = 0; penCount = 0; }
        if (penPass) {
            if (!here) penPass = false;
            for (int k = 0; penPass && k < 8; ++k, ++penCursor) {
                if (penCursor >= animals.size()) { penPass = false; pennedNow = penCount; QuestBook_AtLeast("penned", penCount); break; }
                const Animal& a = animals[penCursor];
                const Vector3 ap = a.pos.toFloat();
                if (a.body != here->name || Vector3Distance(ap, playerLocal) > 120.0f) continue;
                if (IsPenned(*here, ap)) ++penCount;
            }
        }
    }
    const auto tGone0 = std::chrono::steady_clock::now(); // (#1650)
    // Night hunters gone at daylight.
    for (int k = (int)gone.size() - 1; k >= 0; --k) animals.erase(animals.begin() + gone[k]);
    if (!gone.empty()) aimed = -1;
    // (#1223) Creatures' strikes on each other, gathered in the loop (a kill shifts the list), found again by id.
    for (const CreatureHit& h : creatureHits) { const int t = IndexOfNet(h.target); if (t >= 0) HitByCreature(t, h.e, h.attacker, h.from); }
    creatureHits.clear();
    const auto tProj0 = std::chrono::steady_clock::now(); // (#1650)
    StepProjectiles(world, dt); // (#1225) thrown stones and spears in flight
    const auto tBreed0 = std::chrono::steady_clock::now(); // (#1650)
    // Breeding: two fed adults of a kind within 6 m make a young one (once every 5 min of world time).
    for (size_t i = 0; i < animals.size(); ++i) {
        Animal& a = animals[i];
        if (!a.fed || a.age < 1.0f || now - a.lastBred < 300.0) continue;
        for (size_t j = i + 1; j < animals.size(); ++j) {
            Animal& b = animals[j];
            if (b.species != a.species || b.body != a.body || !b.fed || b.age < 1.0f || now - b.lastBred < 300.0) continue;
            if ((a.pos - b.pos).length() > 6.0) continue;
            a.fed = b.fed = false;
            a.lastBred = b.lastBred = now;
            Animal young = a;
            young.pos = young.prevPos = (a.pos + b.pos) * 0.5;
            young.age = 0.0f;
            young.health = species[a.species].health * 0.5f;
            young.followUntil = -1.0;
            young.rng = a.rng ^ 0x9e37u;
            young.netId = NetId_Allocate(NetKind::Animal);
            const bool seen = here && here->name == a.body; // before the push: it may move the vector (a, b)
            const std::string kind = species[a.species].name;
            animals.push_back(young);
            if (seen) GameUI_Toast(TextFormat("A young %s was born", kind.c_str()));
            QuestBook_Count("bred");
            break;
        }
    }
    const auto tTrail0 = std::chrono::steady_clock::now(); // (#1650)
    { // (#1108) the trail: where each animal was, for the server's rewound hit checks
        const double now = SimTime();
        for (Animal& a : animals) {
            if (a.trail.empty() || a.trail.back().first < now) a.trail.push_back({ now, a.pos.toFloat() });
            size_t drop = 0;
            while (drop < a.trail.size() && a.trail[drop].first < now - kTrailSeconds) ++drop;
            if (drop) a.trail.erase(a.trail.begin(), a.trail.begin() + (long)drop);
        }
    }
    lastUpdateMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    { // (#1650) RAYLIBGAME_CREATURE_SLOW=<ms>: an update slower than that prints where its time went
        static const double slowMs = getenv("RAYLIBGAME_CREATURE_SLOW") ? atof(getenv("RAYLIBGAME_CREATURE_SLOW")) : 0.0;
        if (slowMs > 0.0 && lastUpdateMs > slowMs) { auto ms = [](std::chrono::steady_clock::time_point a, std::chrono::steady_clock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
            printf("creatures: slow update %.2f ms: before the loop %.2f, loop %.2f (settle %.2f, brain %.2f; %zu movers), probes %.2f, apply %.2f, pens %.2f, removals+hits %.2f, projectiles %.2f, breeding %.2f, trails %.2f\n",
                   lastUpdateMs, ms(t0, tLoop0), ms(tLoop0, tProbe0), accSettle, accBrain, moveReqs.size(), ms(tProbe0, tApply0), ms(tApply0, tAfter0),
                   ms(tAfter0, tGone0), ms(tGone0, tProj0), ms(tProj0, tBreed0), ms(tBreed0, tTrail0), ms(tTrail0, std::chrono::steady_clock::now()));
            printf("creatures:   slowest brain %.2f ms: %s\n", worstBrain, worstBrainWho.c_str()); }
    }
}

float Animals_RewindDistance(const std::vector<std::pair<double, Vector3>>& trail, Vector3 current, Vector3 from, double now, double window) {
    float best = Vector3Distance(current, from);
    for (const auto& [t, p] : trail) if (t >= now - window) best = std::min(best, Vector3Distance(p, from));
    return best;
}
float Animals_HitDistance(int i, Vector3 from, double window) {
    if (i < 0 || i >= (int)animals.size()) return 1e30f;
    return Animals_RewindDistance(animals[i].trail, animals[i].pos.toFloat(), from, SimTime(), window);
}

// ---- (#1216) Melee in phases (combat step 2, melee.h) ----
// LMB at a creature (aimed at, locked on, or inside the light swing's arc and reach) starts a swing; the rest of LMB stays
// the terrain tool's. The swing runs on sim time with this frame's view: each step of its hit window sweeps the blade
// across the arc against every creature's body capsule, once per creature per swing; terrain and parts in front of a
// creature shield it. Middle mouse toggles a soft lock-on (the view eases toward the creature). RAYLIBGAME_COMBAT2=0: the
// old one-click hit.
namespace {
MeleeState gMelee;
double gMeleeClock = -1.0;
int gMeleeShocks = 0, gMeleeStaggers = 0, gMeleeDodges = 0;
uint64_t gLockNet = 0;              // the soft lock-on's creature (0: none)
Vector2 gLockScreen{ -1.0f, -1.0f }; // where it is on screen (the HUD marker)
WeaponClass gMeleeClass = WeaponClass::Tools;
WeaponSpec gMeleeSpec = Weapon_Spec(WeaponClass::Tools); // (#1219) this swing's spec (the held copy's timing)
int gMeleeHits = 0;                 // hits landed (tests)



// Its body as a capsule along its length (body-local), its radius half its height or width.
void AnimalCapsule(const Animal& a, Vector3& q0, Vector3& q1, float& r) {
    const Species& s = species[a.species];
    const float scale = 0.6f + 0.4f * a.age;
    const Vector3 p = a.pos.toFloat();
    const Vector3 up = Vector3Normalize(p), fwd = AnimalForward(a);
    r = 0.5f * std::max(s.size.x, s.size.y) * scale;
    const Vector3 c = Vector3Add(p, Vector3Scale(up, 0.5f * s.size.y * scale));
    const float hl = std::max(0.0f, 0.5f * s.size.z * scale - r);
    q0 = Vector3Subtract(c, Vector3Scale(fwd, hl)); q1 = Vector3Add(c, Vector3Scale(fwd, hl));
}

Vector3 AnimalCentre(const Animal& a) {
    const Species& s = species[a.species];
    const Vector3 p = a.pos.toFloat();
    return Vector3Add(p, Vector3Scale(Vector3Normalize(p), 0.5f * s.size.y * (0.6f + 0.4f * a.age)));
}

// Terrain or a part between the eye and a point (body-local) on `here`.
bool Shielded(Scene& world, Planetoid& here, Vector3 eye, Vector3 at) {
    const Vector3 d = Vector3Subtract(at, eye);
    const float len = Vector3Length(d);
    if (len < 0.2f) return false;
    const Vector3 dir = Vector3Scale(d, 1.0f / len);
    TerrainHit th;
    if (RaycastTerrain(world, Ray{ here.toWorld(eye), here.dirToWorld(dir) }, len, th) && th.distance < len - 0.15f) return true;
    if (here.plots)
        for (auto& pl : here.plots->plots) { Int3 c, f; float dd; if (pl->raycast(eye, dir, len - 0.15f, c, f, dd)) return true; }
    return false;
}

void MeleeFrame(Scene& world, const PlayerInput& in, Planetoid* here, Vector3 dirW) {
    const double now = SimTime();
    const float dt = gMeleeClock < 0.0 ? 0.0f : (float)std::clamp(now - gMeleeClock, 0.0, 0.1);
    gMeleeClock = now;
    if (gMelee.phase == MeleePhase::Idle) { // (#1219) the held copy's own timing: its weight against its class's iron one
        gMeleeClass = Weapon_ClassOf(terrainTool.heldTool);
        gMeleeSpec = Weapon_SpecOf(terrainTool.heldTool, gToolDataOf ? gToolDataOf(terrainTool.heldTool) : std::string());
    }
    const WeaponSpec& w = gMeleeSpec;
    const bool live = here && !player.flyMode && in.enabled && !buildTool.active && !survival.dead;
    if (!here || survival.dead || player.flyMode) { gMelee = MeleeState{}; gLockNet = 0; gLockScreen = { -1.0f, -1.0f }; return; }
    const Vector3 eye = here->toLocal(world.camera.position);
    const Vector3 fwd = here->dirToLocal(dirW);
    const Vector3 up = Vector3Length(eye) > 1.0f ? Vector3Normalize(eye) : Vector3{ 0, 1, 0 };

    // The soft lock-on: middle mouse takes the creature nearest the view's centre (within 35 degrees and 12 m), again lets go.
    if (live && in.hit(Act::LockOn)) {
        if (gLockNet) gLockNet = 0;
        else {
            float best = cosf(MeleeTune::kLockConeDeg * DEG2RAD);
            for (const Animal& a : animals) {
                if (a.body != here->name) continue;
                const Vector3 to = Vector3Subtract(AnimalCentre(a), eye);
                const float d = Vector3Length(to);
                if (d > MeleeTune::kLockRange || d < 0.1f) continue;
                const float c = Vector3DotProduct(Vector3Scale(to, 1.0f / d), fwd);
                if (c > best && !Shielded(world, *here, eye, AnimalCentre(a))) { best = c; gLockNet = a.netId; }
            }
            if (!gLockNet) GameUI_Toast("Nothing to lock on to");
        }
    }
    int lockIdx = gLockNet ? IndexOfNet(gLockNet) : -1;
    if (gLockNet && (lockIdx < 0 || animals[lockIdx].body != here->name ||
                     Vector3Distance(AnimalCentre(animals[lockIdx]), eye) > MeleeTune::kLockKeep)) { gLockNet = 0; lockIdx = -1; }
    gLockScreen = lockIdx >= 0 ? GetWorldToScreen(here->toWorld(AnimalCentre(animals[lockIdx])), world.camera) : Vector2{ -1.0f, -1.0f };

    // A creature the light swing would reach (from the view as it is now).
    auto inArc = [&]() {
        const float half = 0.5f * w.arcDeg * DEG2RAD;
        for (const Animal& a : animals) {
            if (a.body != here->name || Vector3Distance(a.pos.toFloat(), eye) > w.reach + 4.0f) continue;
            Vector3 q0, q1, hit; float r;
            AnimalCapsule(a, q0, q1, r);
            if (Melee_SweepCapsule(eye, fwd, up, w.reach, half, -half, q0, q1, r, &hit) && !Shielded(world, *here, eye, hit)) return true;
        }
        return false;
    };
    // (#1639) A weapon (knife, sword, mace, spear, two-handed) always swings: what it hits is the sweep's reach against the
    // hitboxes below, not the crosshair. Tools and axes still swing only at a creature (else LMB mines, digs or fells).
    const bool weaponInHand = terrainTool.heldTool.rfind("weapon_", 0) == 0 && gMeleeClass != WeaponClass::Axe; // (composed weapons; not parts or machines named like them)
    MeleeInput mi;
    mi.pressed = live && in.hit(Act::Primary) && !gPrimaryLatch.latched &&
                 (weaponInHand || gMelee.phase != MeleePhase::Idle || aimed >= 0 || lockIdx >= 0 || inArc());
    mi.held = live && in.down(Act::Primary) && survival.windedLeft <= 0.0f; // (#1218) winded: no heavies
    mi.block = survival.blocking;
    mi.dodge = survival.dodges != gMeleeDodges; // (#1218) a dodge started
    gMeleeDodges = survival.dodges;
    mi.interrupted = survival.shocks != gMeleeShocks || survival.staggers != gMeleeStaggers; // a shock or (#1217) a stagger
    gMeleeShocks = survival.shocks; gMeleeStaggers = survival.staggers;
    mi.stamina = survival.stamina;
    const MeleeStep o = Melee_Step(gMelee, w, mi, dt);
    survival.hyperArmour = Combat3_On() && Melee_HyperArmour(gMelee, w); // (#1217) hits don't stagger this wind-up
    if (o.cost > 0.0f) { survival.stamina = std::max(0.0f, survival.stamina - o.cost); survival.staminaUsedAt = now; }
    if (o.tired) GameUI_Toast("Too tired to swing");
    if (o.started) { hitLatch = true; ++terrainTool.swings; Audio_Play(Sfx::Dig, 0.25f, 1.9f); }
    if (!o.sweeping) return;

    // The hit window: the blade from this step's start angle to its end, against every creature on this body.
    const float a0 = Melee_SweepAngle(gMelee, w, o.sweepFrom), a1 = Melee_SweepAngle(gMelee, w, o.sweepTo);
    struct Hit { int i; Vector3 at; };
    std::vector<Hit> hits;
    for (size_t i = 0; i < animals.size(); ++i) {
        const Animal& a = animals[i];
        if (a.body != here->name || Vector3Distance(a.pos.toFloat(), eye) > w.reach + 4.0f || Melee_Struck(gMelee, a.netId)) continue;
        Vector3 q0, q1, at; float r;
        AnimalCapsule(a, q0, q1, r);
        if (!Melee_SweepCapsule(eye, fwd, up, w.reach, a0, a1, q0, q1, r, &at)) continue;
        if (Shielded(world, *here, eye, at)) continue; // terrain or a part takes the blow
        hits.push_back({ (int)i, at });
    }
    const float base = Weapon_BaseDamage(terrainTool.heldTool, gToolDataOf ? gToolDataOf(terrainTool.heldTool) : std::string()); // (#1219) a composed copy: from its head
    std::sort(hits.begin(), hits.end(), [](const Hit& x, const Hit& y) { return x.i > y.i; }); // a kill removes it: high first
    for (const Hit& h : hits) {
        const Animal& a = animals[h.i];
        gMelee.struck.push_back(a.netId);
        DamageEvent e;
        e.amount = Melee_Damage(base, gMelee);
        e.kind = DamageKind::Melee;
        e.type = Weapon_TypeOf(terrainTool.heldTool);
        e.heavy = gMelee.heavy;
        e.poise = Melee_Poise(gMelee, w);
        { const GemEffect ge = Weapon_GemEffect(terrainTool.heldTool, gToolDataOf ? gToolDataOf(terrainTool.heldTool) : std::string()); // (#1219) set gems
          if (ge.status != StatusKind::None) e.status = ge.status;
          e.poise *= ge.poiseScale; }
        e.fromBehind = Vector3DotProduct(AnimalForward(a), Vector3Normalize(Vector3Subtract(h.at, eye))) > 0.5f; // its back to the blade
        { const CreatureSpec* fs = Creature_Spec(species[a.species].id); // (#1227) at its face: within its shield's arc of its facing
          e.fromFront = Vector3DotProduct(AnimalForward(a), Vector3Normalize(Vector3Subtract(eye, a.pos.toFloat()))) > cosf((fs ? fs->shieldArcDeg : 60.0f) * DEG2RAD); }
        survival.lastDealtAt = now;
        ++gMeleeHits;
        Noise_Emit(here->id, eye, NoiseTune::kHit, simClock, kPlayerPercept); // (#1224) combat is loud
        if (gAnimalsReplicated) { if (gAnimalHitRequest) gAnimalHitRequest(a.netId, (e.heavy ? 1 : 0) | (e.fromBehind ? 2 : 0)); } // (#334) a client: the server decides
        else {
            const CreatureSpec* hs = Creature_Spec(species[animals[h.i].species].id);
            if (hs && hs->corrode > 0 && ToolOf(terrainTool.heldTool).type != ToolType::Hands) ToolWear_OnStroke(terrainTool.heldTool, std::max(1, hs->corrode / 3)); // (#1226) a tool that strikes an ooze corrodes
            Animals_HitEvent(h.i, e);
        }
    }
}
} // namespace

bool Animals_MeleePose(int& phase, float& f, float& charge, bool& back) {
    if (!Combat2_On() || gMelee.phase == MeleePhase::Idle) return false;
    const WeaponSpec& w = gMeleeSpec; // (#1219) the pose follows the swing's own timing
    phase = (int)gMelee.phase;
    back = gMelee.chain == 1;
    charge = gMelee.heavy ? 1.0f : 0.0f;
    switch (gMelee.phase) {
        case MeleePhase::WindUp:
            f = std::min(1.0f, gMelee.t / w.windUp);
            if (gMelee.heavy) charge = std::clamp((gMelee.t - w.windUp) / (w.windUp * (MeleeTune::kHeavyWindUp - 1.0f)), 0.0f, 1.0f);
            break;
        case MeleePhase::Hit: f = gMelee.sweep; break;
        default: f = std::min(1.0f, gMelee.t / w.recovery); break;
    }
    return true;
}

static Vector3 gLockPrevTo{}; static uint64_t gLockPrevNet = 0; // (#1216) the bearing last frame, and to which lock
void Animals_LockOnView(Scene& world, float dt) {
    const int i = gLockNet ? IndexOfNet(gLockNet) : -1;
    Planetoid* here = player.nearest;
    if (i < 0 || !here || player.flyMode) { gLockPrevNet = 0; return; }
    const Vector3 to = Vector3Subtract(here->toWorld(AnimalCentre(animals[i])), world.camera.position);
    // (#1216 playtest: "does not lock to the enemy hard enough when strafing sideways") The ease alone trails a target whose
    // bearing keeps turning (strafing at 5 m/s, 3 m away: ~1.7 rad/s x 0.12 s = ~11 degrees behind). First the view turns
    // with the bearing's own change since last frame, then the ease takes out what is left of the offset.
    Vector3& prevTo = gLockPrevTo; uint64_t& prevLock = gLockPrevNet;
    if (prevLock == gLockNet && Vector3Length(prevTo) > 1e-3f && Vector3Length(to) > 1e-3f) {
        const Vector3 turned = Vector3RotateByQuaternion(player.forward, QuaternionFromVector3ToVector3(Vector3Normalize(prevTo), Vector3Normalize(to)));
        if (fabsf(Vector3DotProduct(Vector3Normalize(turned), Vector3Normalize(player.up))) < 0.985f) player.forward = Vector3Normalize(turned); // the look's pitch clamp
    }
    prevTo = to; prevLock = gLockNet;
    player.forward = Melee_SoftTurn(player.forward, player.up, to, dt);
    // (#1218) A sideways dodge just started (more across the view than along it): it circles the locked creature, the way
    // the dodge was going (player.cpp follows the circle each step).
    if (player.dodgeLeft > 0.0f && player.dodgeOrbitSign == 0.0f && !player.dodgeOrbitDone) {
        player.dodgeOrbitDone = true;
        const Vector3 n = Vector3Normalize(player.up);
        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(n, Vector3DotProduct(player.forward, n))));
        if (fabsf(Vector3DotProduct(Vector3Normalize(player.dodgeDir), flat)) < 0.5f) {
            const DVec3 c = here->toWorldD(DVec3(AnimalCentre(animals[i])));
            const Vector3 rel = (player.position - c).toFloat();
            const Vector3 tangent = Vector3CrossProduct(n, Vector3Subtract(rel, Vector3Scale(n, Vector3DotProduct(rel, n))));
            if (Vector3Length(tangent) > 1e-3f) { player.dodgeOrbit = c; player.dodgeOrbitSign = Vector3DotProduct(tangent, player.dodgeDir) >= 0.0f ? 1.0f : -1.0f; }
        }
    }
    if (player.dodgeLeft <= 0.0f) player.dodgeOrbitDone = false;
}

uint64_t Animals_LockOn() { return gLockNet; }
// (#1218) The lock-on creature's centre in the world (doubles), for a dodge that circles it; false: no lock (or it's gone).
bool Animals_LockOnWorld(Scene& world, DVec3& out) {
    if (!gLockNet) return false;
    const int i = Animals_FindNet(gLockNet);
    if (i < 0) return false;
    Planetoid* b = BodyNamed(world, animals[i].body);
    if (!b) return false;
    out = b->toWorldD(DVec3(AnimalCentre(animals[i])));
    return true;
}
int Animals_MeleeHits() { return gMeleeHits; }
const MeleeState& Animals_Melee() { return gMelee; }
bool Animals_MeleeBusy() { return Combat2_On() && gMelee.phase != MeleePhase::Idle; }

// Aim along the crosshair and the player's actions on the aimed animal (#222). Called after the camera is placed
// for this frame (Animals_Update ran on last frame's), right before the terrain tool, which holds off while
// Animals_Aiming(). The nearer of an animal and the ground wins: an animal behind a slope or a rock is not aimed.
// Reach matches the tool's. After a hit LMB stays with the animals until released, so a held button doesn't dig
// the ground behind one that bolted.
void Animals_Aim(Scene& world, const PlayerInput& in) {
    Planetoid* here = player.nearest;
    const double now = WorldTime();
    if (!in.down(Act::Primary)) hitLatch = false;
    aimed = -1;
    if (here && !player.flyMode) {
        const Vector3 dirW = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        const Vector3 o = here->toLocal(world.camera.position);
        const Vector3 d = here->dirToLocal(dirW);
        float best = std::max(5.0f, terrainTool.reach);
        for (size_t i = 0; i < animals.size(); ++i) {
            const Animal& a = animals[i];
            if (a.body != here->name) continue;
            const Species& s = species[a.species];
            const float scale = 0.6f + 0.4f * a.age;
            const Vector3 ap = a.pos.toFloat();
            const Vector3 c = Vector3Add(ap, Vector3Scale(Vector3Normalize(ap), s.size.y * 0.5f * scale));
            const Vector3 oc = Vector3Subtract(c, o);
            const float t = Vector3DotProduct(oc, d);
            // The body's sphere (half its length or height), a little wider with distance: easier to land a hit.
            const float r = 0.5f * std::max(s.size.z, s.size.y) * scale + 0.1f + 0.03f * std::max(t, 0.0f);
            if (t < 0.0f || t > best) continue;
            if (Vector3Length(Vector3Subtract(oc, Vector3Scale(d, t))) > r) continue;
            best = t; aimed = (int)i;
        }
        // Ground in front of it: the terrain is what's aimed at.
        TerrainHit th;
        if (aimed >= 0 && RaycastTerrain(world, Ray{ world.camera.position, dirW }, best, th) && th.distance < best - 0.3f) aimed = -1;
    }
    // (#1216) Melee in phases: LMB at a creature goes to the swing (the old one-click hit below with RAYLIBGAME_COMBAT2=0).
    if (Combat2_On()) MeleeFrame(world, in, here, Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position)));
    if (aimed >= 0 && in.enabled && !buildTool.active) {
        Animal& a = animals[aimed];
        const Species& s = species[a.species];
        if (!Combat2_On() && in.hit(Act::Primary) && !gPrimaryLatch.latched && SimTime() - lastHit > 0.35) { // (#896: not while eating; #895: not on a used-up hold)
            lastHit = SimTime();
            hitLatch = true;
            const ToolSpec t = ToolOf(terrainTool.heldTool);
            // A swing at a creature costs stamina (#165: 4 with bare hands, 8 with a tool); with too little left, no swing.
            const float cost = t.type == ToolType::Hands ? 4.0f : 8.0f;
            if (survival.stamina < cost) GameUI_Toast("Too tired to swing");
            else {
                survival.stamina -= cost;
                survival.staminaUsedAt = SimTime();
                ++terrainTool.swings;
                if (gAnimalsReplicated) { if (gAnimalHitRequest) gAnimalHitRequest(animals[aimed].netId, 0); } // (#334) a client: the server decides
                else Animals_Hit(aimed, t.type == ToolType::Hands ? 2.0f : 3.0f + 3.0f * t.tier);
            }
        } else if (in.hit(Act::Insert) && gAnimalsReplicated) {
            if (gAnimalUseRequest) gAnimalUseRequest(a.netId, 1); // (#334) a client: the server feeds it from our inventory
        } else if (in.hit(Act::Extract) && gAnimalsReplicated) {
            if (gAnimalUseRequest) gAnimalUseRequest(a.netId, 2);
        } else if (in.hit(Act::Insert)) {
            Animals_Feed(aimed);
        } else if (in.hit(Act::Extract)) {
            if (s.product.empty() || s.product == "-") GameUI_Toast(TextFormat("A %s gives nothing but meat", s.name.c_str()));
            else if (a.age < 1.0f) GameUI_Toast("Too young");
            else if (now < a.productReady) GameUI_Toast(TextFormat("Not yet (%.0f s)", a.productReady - now));
            else if (s.productTool != "-" && inventory.count(s.productTool) <= 0) GameUI_Toast(TextFormat("You need %s", s.productTool.c_str()));
            else {
                inventory.add(s.product, s.product == "wool" ? 2 : 1);
                a.productReady = now + s.regrow;
                Audio_Play(Sfx::Pickup, 0.5f);
            }
        }
    }
}

// The animal model (Plan 5 A): built from boxes in the animal's frame (x right, y up, z forward), proportioned
// from animals.txt (length, height, width) and its model features. Legs swing on the gait (trot: diagonal
// pairs; bipeds alternate) with knees bending on the back swing; the head dips to graze; resting folds the
// legs and lowers the body; the tail sways. Box() draws one part rotated about its own top (a joint).
static inline void DrawCubeB(Vector3 p, float w, float h, float l, Color c) { IB_Cube(p, { w, h, l }, c); } // batched, lit (instancebatch.h)
template <typename ShadeFn>
static void DrawAnimalModel(const Species& s, const Animal& a, ShadeFn shade) {
    const Vector3 sz = s.size;                       // length (z), height (y), width (x)
    if (s.blob) { // (#1226) a mineral ooze: a wobbling heap, half its own green and half the colour of the ore it ate, with lumps of it
        const ItemDef* od = a.ore.empty() ? nullptr : itemDB.find(a.ore);
        const Color oc = od ? od->color : s.head;
        const Color mix = { (unsigned char)((s.body.r + oc.r) / 2), (unsigned char)((s.body.g + oc.g) / 2), (unsigned char)((s.body.b + oc.b) / 2), 255 };
        const float t = (float)GameClock(), w = 1.0f + 0.07f * sinf(t * 2.3f + (float)(a.rng % 17));
        const float g = a.generation == 0 ? 1.0f : a.generation == 1 ? 0.8f : 0.62f; // the split-off ones are smaller
        DrawCubeB({ 0, sz.y * 0.28f * w * g, 0 }, sz.z * g / w, sz.y * 0.56f * w * g, sz.x * g / w, shade(mix, 1.0f));
        DrawCubeB({ 0, sz.y * 0.6f * w * g, 0 }, sz.z * 0.68f * g, sz.y * 0.4f * g, sz.x * 0.68f * g, shade(mix, 1.12f));
        for (int k = 0; k < 3; ++k) DrawCubeB({ (k - 1) * sz.z * 0.28f * g, sz.y * (0.55f + 0.1f * (k & 1)) * g, (k == 1 ? 0.25f : -0.15f) * sz.x * g }, sz.z * 0.16f * g, sz.y * 0.16f * g, sz.z * 0.16f * g, shade(oc, 0.95f));
        return;
    }
    if (s.turret) { // (#1228) a ruin sentinel: stone pedestal, a bronze head on it facing its heading, a barrel, and an eye that
        // glows from dull red to white-hot as a shot charges (the telegraph you can see)
        const CreatureSpec* tc = Creature_Spec(s.id);
        const float charge = a.throwWindAt >= 0.0 && tc ? std::clamp((float)((simClock - a.throwWindAt) / std::max(0.1f, tc->throwWindUp)), 0.0f, 1.0f) : 0.0f;
        DrawCubeB({ 0, sz.y * 0.2f, 0 }, sz.z * 0.9f, sz.y * 0.4f, sz.x * 0.9f, shade(s.body, 0.9f));            // pedestal
        DrawCubeB({ 0, sz.y * 0.42f, 0 }, sz.z * 0.6f, sz.y * 0.08f, sz.x * 0.6f, shade(s.body, 1.1f));          // collar
        DrawCubeB({ 0, sz.y * 0.62f, 0 }, sz.z * 0.7f, sz.y * 0.32f, sz.x * 0.7f, shade(s.head, 1.0f));          // head
        DrawCubeB({ 0, sz.y * 0.64f, sz.x * 0.55f }, sz.z * 0.14f, sz.y * 0.12f, sz.x * 0.5f, shade(s.head, 0.7f)); // barrel
        const Color eye = { 255, (unsigned char)(60 + 195 * charge), (unsigned char)(40 + 200 * charge * charge), 255 };
        DrawCubeB({ 0, sz.y * 0.72f, sz.x * 0.36f }, sz.z * 0.22f, sz.y * 0.1f, sz.x * 0.04f, eye);              // the eye
        return;
    }
    if (s.nest) { // (#1223) a nest: a low mound of lumps, the egg-like pods on top (no legs, no head, never moves)
        DrawCubeB({ 0, sz.y * 0.3f, 0 }, sz.z, sz.y * 0.6f, sz.x, shade(s.body, 1.0f));
        DrawCubeB({ sz.z * 0.2f, sz.y * 0.55f, -sz.x * 0.15f }, sz.z * 0.55f, sz.y * 0.5f, sz.x * 0.5f, shade(s.body, 0.9f));
        for (int k = 0; k < 3; ++k) DrawCubeB({ (k - 1) * sz.z * 0.25f, sz.y * 0.8f, (k % 2 ? 0.15f : -0.1f) * sz.x }, sz.z * 0.18f, sz.y * 0.3f, sz.z * 0.18f, shade(s.head, 1.0f));
        return;
    }
    const float t = (float)GameClock();
    const bool resting = a.state == State::Rest, grazing = a.state == State::Graze && a.moveSpeed <= 0.0f;
    const float legH = sz.y * (s.biped ? 0.45f : 0.42f);
    const float lower = resting ? legH * 0.75f : 0.0f;              // folded legs: the body sinks
    const float swing = std::min(0.75f, a.moveSpeed * 0.22f);       // stride angle grows with speed
    auto box = [&](Vector3 at, Vector3 size, float pitchDeg, float yawDeg, Color c) {
        rlPushMatrix();
        rlTranslatef(at.x, at.y, at.z);
        if (yawDeg != 0.0f) rlRotatef(yawDeg, 0, 1, 0);
        if (pitchDeg != 0.0f) rlRotatef(pitchDeg, 1, 0, 0);
        DrawCubeB({ 0, -size.y * 0.5f, 0 }, size.x, size.y, size.z, c);  // hangs from its joint
        rlPopMatrix();
    };
    // Body: a torso with a deeper chest and a slightly narrower rump.
    const float bodyY = legH - lower, bodyH = sz.y - legH, len = sz.x * (s.longBody ? 1.1f : 1.0f);
    const float breathe = 1.0f + 0.015f * sinf(t * 1.7f + a.rng % 7);
    const Color bodyC = shade(s.body, 1.0f);
    DrawCubeB({ 0, bodyY + bodyH * 0.5f, len * 0.08f }, sz.z * breathe, bodyH * breathe, len * 0.62f, bodyC);
    DrawCubeB({ 0, bodyY + bodyH * 0.45f, -len * 0.28f }, sz.z * 0.9f, bodyH * 0.88f, len * 0.4f, shade(s.body, 0.93f));
    if (s.wool)   // fleece: puffy lumps over the back
        for (int k = 0; k < 3; ++k) DrawCubeB({ (k - 1) * sz.z * 0.18f, bodyY + bodyH * 0.92f, len * (0.2f - 0.22f * k) }, sz.z * 0.6f, bodyH * 0.35f, len * 0.3f, shade(s.body, 1.08f));
    if (s.spikes) // a ridge of spines down the back
        for (int k = 0; k < 5; ++k) box({ 0, bodyY + bodyH + 0.12f * sz.y, len * (0.3f - 0.15f * k) }, { sz.z * 0.08f, sz.y * 0.22f, len * 0.06f }, 180.0f - 20.0f, 0.0f, shade(s.head, 0.8f));
    // Neck and head: forward and up; down to the grass when grazing; a little bob while walking.
    // The neck hangs along -y from its base: a pitch of -125 deg points it forward and up, -55 forward and down.
    const float neckPitch = grazing ? -55.0f : resting ? -100.0f : -125.0f + 6.0f * sinf(a.gait * 2.0f) * (swing > 0.0f ? 1.0f : 0.0f);
    const Vector3 neckBase = { 0, bodyY + bodyH * 0.75f, len * 0.36f };
    const float neckLen = sz.y * (s.biped ? 0.35f : 0.5f);
    rlPushMatrix();
    rlTranslatef(neckBase.x, neckBase.y, neckBase.z);
    rlRotatef(neckPitch, 1, 0, 0);                 // the neck hangs "down" from its base, rotated up and forward
    DrawCubeB({ 0, -neckLen * 0.5f, 0 }, sz.z * 0.38f, neckLen, sz.z * 0.38f, shade(s.body, 0.97f));
    rlTranslatef(0, -neckLen, 0);
    rlRotatef(-neckPitch + (grazing ? 40.0f : 0.0f), 1, 0, 0); // the head held level (nose down when grazing)
    const float hs = std::min(sz.y, sz.z) * (s.biped ? 0.5f : 0.55f);
    const Color headC = shade(s.head, 1.0f);
    DrawCubeB({ 0, 0, hs * 0.25f }, hs * 0.8f, hs * 0.85f, hs * 1.05f, headC);
    DrawCubeB({ 0, -hs * 0.15f, hs * 0.85f }, hs * 0.55f, hs * 0.5f, hs * 0.5f, shade(s.head, 0.9f));            // snout
    if (s.beak) DrawCubeB({ 0, -hs * 0.05f, hs * 0.95f }, hs * 0.3f, hs * 0.22f, hs * 0.45f, shade({ 230, 170, 50, 255 }, 1.0f));
    DrawCubeB({ hs * 0.25f, hs * 0.15f, hs * 0.75f }, hs * 0.1f, hs * 0.1f, hs * 0.05f, shade({ 20, 20, 20, 255 }, 1.0f)); // eyes
    DrawCubeB({ -hs * 0.25f, hs * 0.15f, hs * 0.75f }, hs * 0.1f, hs * 0.1f, hs * 0.05f, shade({ 20, 20, 20, 255 }, 1.0f));
    if (s.ears) {
        const float flick = 8.0f * sinf(t * 2.3f + a.rng % 11);
        box({ hs * 0.32f, hs * 0.45f, 0 }, { hs * 0.18f, hs * 0.4f, hs * 0.08f }, 180.0f + 20.0f + flick, -20.0f, headC);
        box({ -hs * 0.32f, hs * 0.45f, 0 }, { hs * 0.18f, hs * 0.4f, hs * 0.08f }, 180.0f + 20.0f - flick, 20.0f, headC);
    }
    if (s.horns) {
        const Color hc = shade({ 225, 215, 190, 255 }, 1.0f);
        box({ hs * 0.3f, hs * 0.4f, hs * 0.1f }, { hs * 0.12f, hs * 0.55f, hs * 0.12f }, 180.0f + 35.0f, -35.0f, hc);
        box({ -hs * 0.3f, hs * 0.4f, hs * 0.1f }, { hs * 0.12f, hs * 0.55f, hs * 0.12f }, 180.0f + 35.0f, 35.0f, hc);
    }
    rlPopMatrix();
    // Tail: from the rump, swaying.
    if (s.tail) {
        const float sway = 18.0f * sinf(t * 2.5f + a.rng % 13) + (a.moveSpeed > 0.0f ? 10.0f * sinf(a.gait) : 0.0f);
        box({ 0, bodyY + bodyH * 0.8f, -len * 0.48f }, { sz.z * 0.12f, sz.y * (s.biped ? 0.3f : 0.45f), sz.z * 0.12f }, s.biped ? 140.0f : 35.0f, sway, shade(s.body, 0.85f));
    }
    // Legs: two segments, swinging from the hips; diagonal pairs in step (bipeds: alternating).
    // (#1224) Six legs (hexa): three pairs splayed out to the sides, an alternating tripod gait, the body slung low.
    const int legs = s.biped ? 2 : s.hexa ? 6 : 4;
    const Color legC = shade(s.body, 0.8f);
    for (int l = 0; l < legs; ++l) {
        const bool right = l & 1, front = s.biped ? true : s.hexa ? (l >> 1) == 0 : (l & 2) != 0;
        const int row = l >> 1; // hexa: 0 front, 1 middle, 2 back
        const float phase = s.biped ? (right ? PI : 0.0f) : s.hexa ? ((right != (row == 1)) ? PI : 0.0f) : ((right != front) ? PI : 0.0f);
        const float ang = resting ? 80.0f : swing * sinf(a.gait + phase) * RAD2DEG;
        const float knee = resting ? -150.0f : std::max(0.0f, -sinf(a.gait + phase)) * swing * 60.0f;
        const float hipZ = s.biped ? 0.0f : s.hexa ? (0.34f - 0.34f * row) * len : (front ? 0.3f : -0.33f) * len;
        const Vector3 hip = { (right ? 0.32f : -0.32f) * sz.z * (s.hexa ? 1.3f : 1.0f), bodyY + bodyH * 0.15f, hipZ };
        const float upper = legH * 0.55f, lowerLen = legH * 0.5f, w = sz.z * (s.biped ? 0.1f : s.hexa ? 0.1f : 0.16f);
        rlPushMatrix();
        rlTranslatef(hip.x, hip.y, hip.z);
        if (s.hexa) rlRotatef(right ? 40.0f : -40.0f, 0, 0, 1); // splayed out sideways
        rlRotatef(ang, 1, 0, 0);
        DrawCubeB({ 0, -upper * 0.5f, 0 }, w * 1.15f, upper, w * 1.15f, legC);
        rlTranslatef(0, -upper, 0);
        rlRotatef(front ? -knee : knee, 1, 0, 0);
        DrawCubeB({ 0, -lowerLen * 0.5f, 0 }, w, lowerLen, w, legC);
        DrawCubeB({ 0, -lowerLen, w * 0.15f }, w * 1.1f, w * 0.5f, w * 1.4f, shade(s.head, 0.6f));   // hoof / foot
        rlPopMatrix();
    }
}

void Animals_Draw3D(Scene& world) {
    DrawProjectiles(world); // (#1225) stones and spears in flight
    if (animals.empty()) return;
    static std::vector<AvatarState> humanoids; humanoids.clear(); // (#1225) drawn after the loop, through the avatar rig
    static std::vector<Matrix> shieldMats; shieldMats.clear(); // (#1227) plank shields, drawn after the avatars
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    std::unordered_map<std::string, Planetoid*> bodies;
    for (size_t i = 0; i < animals.size(); ++i) {
        const Animal& a = animals[i];
        auto bit = bodies.find(a.body);
        Planetoid* body = bit != bodies.end() ? bit->second : (bodies[a.body] = BodyNamed(world, a.body));
        if (!body) continue;
        // Between the last two simulation steps (#315): smooth at any frame rate (a jump over 2 m is a teleport).
        DVec3 at = a.pos;
        float heading = a.heading;
        if (a.prevValid && (a.prevPos - a.pos).length() < 2.0) {
            const double al = SimAlpha();
            at = a.prevPos + (a.pos - a.prevPos) * al;
            heading = a.prevHeading + remainderf(a.heading - a.prevHeading, 2.0f * PI) * (float)al;
        }
        const Vector3 base = RenderPos(world, *body, at); // (#1698)
        if (Vector3Distance(base, world.camera.position) > 160.0f) continue;
        const Species& s = species[a.species];
        const Vector3 lup = Vector3Normalize(at.toFloat());
        Vector3 east, north;
        TangentFrame(lup, east, north);
        const Vector3 f = body->dirToWorld(Vector3Add(Vector3Scale(east, cosf(heading)), Vector3Scale(north, sinf(heading))));
        const Vector3 u = body->dirToWorld(a.clinging ? Vector3Negate(lup) : lup); // (#1224) upside down on a ceiling
        const Vector3 r = Vector3CrossProduct(u, f); // right-handed (r, u, f)
        const float k = 0.6f + 0.4f * a.age;
        const float light = 1.0f; // lit by the instance batch shader (sun, shadows, headlamp)
        auto shade = [&](Color c, float m) { return Color{ (unsigned char)(c.r * light * m), (unsigned char)(c.g * light * m), (unsigned char)(c.b * light * m), 255 }; };
        if (s.humanoid) { // (#1225) a cave scavenger: the avatar rig in ragged colours, a torch in hand (the thrower: a stone)
            AvatarState av;
            av.session = 0x40000000u | (uint32_t)(a.netId & 0x3FFFFFFFu);
            av.position = body->toWorldD(at); // (absolute world, as player.position)
            av.up = u; av.facing = f;
            av.velocity = Vector3Scale(f, a.moveSpeed);
            av.held = a.thrower ? "stone" : !s.held.empty() ? s.held : "torch"; // (#1227: a miner's rusted pick)
            const bool winding = a.windUpAt >= 0.0 || a.throwWindAt >= 0.0;
            av.action = winding ? 1 : 0;
            av.actionPhase = winding ? (float)std::min(1.0, (simClock - std::max(a.windUpAt, a.throwWindAt)) / 0.6) : 0.0f;
            av.suitColor = s.body;
            humanoids.push_back(av);
            if (s.plankShield) shieldMats.push_back(Matrix{ r.x, u.x, f.x, base.x, r.y, u.y, f.y, base.y, r.z, u.z, f.z, base.z, 0, 0, 0, 1 }); // (#1227)
            continue;
        }
        Matrix m = { r.x, u.x, f.x, base.x, r.y, u.y, f.y, base.y, r.z, u.z, f.z, base.z, 0, 0, 0, 1 };
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        rlScalef(k, k, k);
        DrawAnimalModel(s, a, shade);
        if ((int)i == aimed) DrawCubeWires({ 0, s.size.y * 0.5f, 0 }, s.size.z * 1.15f, s.size.y * 1.1f, s.size.x * 1.15f, YELLOW);
        rlPopMatrix();
    }
    if (!humanoids.empty()) Avatar_Draw3D(world, humanoids); // (#1225)
    for (const Matrix& sm : shieldMats) { // (#1227) the lost miner's plank shield on its arm: three boards and a rusted band, facing out
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(sm));
        for (int k = 0; k < 3; ++k) DrawCubeB({ -0.3f + (k - 1) * 0.15f, 1.05f, 0.34f }, 0.14f, 0.7f, 0.05f, Color{ (unsigned char)(120 + 10 * k), 88, 56, 255 });
        DrawCubeB({ -0.3f, 1.2f, 0.37f }, 0.48f, 0.07f, 0.03f, Color{ 105, 70, 48, 255 });
        rlPopMatrix();
    }
}

void Animals_DrawHUD() {
    if (gLockNet && gLockScreen.x >= 0.0f) { // (#1216) the soft lock-on's marker: a diamond around the creature's centre
        const Vector2 p = gLockScreen; const float r = 14.0f; const Color c = Color{ 255, 200, 90, 230 };
        DrawLineEx({ p.x, p.y - r }, { p.x + r, p.y }, 2.0f, c); DrawLineEx({ p.x + r, p.y }, { p.x, p.y + r }, 2.0f, c);
        DrawLineEx({ p.x, p.y + r }, { p.x - r, p.y }, 2.0f, c); DrawLineEx({ p.x - r, p.y }, { p.x, p.y - r }, 2.0f, c);
    }
    if (aimed < 0 || aimed >= (int)animals.size()) return;
    const Animal& a = animals[aimed];
    const Species& s = species[a.species];
    const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
    DrawRectangle(x - 8, y - 6, 380, 74, Fade(BLACK, 0.55f));
    const char* state = Combat3_On() && a.poise.max > 0.0f && Poise_Staggered(a.poise) ? "staggered" : a.state == State::Hunt ? "hunting you" : a.state == State::Flee ? "fleeing" : a.state == State::Follow ? "follows you" : a.state == State::Rest ? "resting" : a.state == State::Walk ? "wandering" : "grazing";
    DrawText(TextFormat("%s%s  (%s%s)", a.age < 1.0f ? "young " : "", s.name.c_str(), state, a.fed ? ", fed" : ""), x, y, 20, Color{ 240, 220, 160, 255 });
    DrawRectangle(x, y + 26, 200, 6, Fade(RAYWHITE, 0.2f));
    DrawRectangle(x, y + 26, (int)(200 * std::max(0.0f, a.health) / s.health), 6, RED);
    std::string diet;
    for (const auto& d : s.diet) diet += (diet.empty() ? "" : "/") + d;
    DrawText(TextFormat("I feed (%s)   O %s   LMB attack%s", diet.c_str(), s.product.empty() || s.product == "-" ? "-" : s.product.c_str(), Combat2_On() ? " (hold: heavy)   MMB lock on" : ""), x, y + 40, 15, Fade(RAYWHITE, 0.75f));
}

void Animals_Save(FILE* f) {
    for (const Animal& a : animals) {
        fprintf(f, "animal %s %s %f %f %f %f %f %f %d %.2f %.2f\n", species[a.species].id.c_str(), a.body.c_str(), a.pos.x, a.pos.y, a.pos.z,
                a.heading, a.health, a.age, a.fed ? 1 : 0, a.lastBred, a.productReady);
        fprintf(f, "animal_id %llu\n", (unsigned long long)a.netId); // its own line: older builds skip it (#317)
        if (a.hasHome) fprintf(f, "animal_home %f %f %f\n", a.home.x, a.home.y, a.home.z); // (#1633) a cave creature's home (older builds skip it)
    }
    for (const auto& kv : spawnedTiles) {
        fprintf(f, "animal_tiles %s %zu", kv.first.c_str(), kv.second.size());
        for (uint64_t k : kv.second) fprintf(f, " %llu", (unsigned long long)k);
        fprintf(f, "\n");
    }
}

static int lastLoaded = -1; // the animal the next "animal_id" line belongs to
int Animals_SimulatedCount() { return simulatedNow; }
int Animals_SimulatedNear(const Planetoid& body, Vector3 local, float range) {
    int n = 0;
    for (const Animal& a : animals) n += a.body == body.name && a.simFrame == frameNo && Vector3Distance(a.pos.toFloat(), local) < range ? 1 : 0;
    return n;
}
uint64_t Animals_NetId(int index) { return index >= 0 && index < (int)animals.size() ? animals[index].netId : 0; }

bool Animals_LoadLine(const char* tok, FILE* f) {
    LoadSpecies();
    char id[64] = {}, body[128] = {};
    if (strcmp(tok, "animal") == 0) {
        Animal a;
        int fed = 0;
        if (fscanf(f, " %63s %127s %lf %lf %lf %f %f %f %d %lf %lf", id, body, &a.pos.x, &a.pos.y, &a.pos.z, &a.heading, &a.health, &a.age, &fed, &a.lastBred, &a.productReady) == 11) {
            a.species = SpeciesIndex(id);
            if (a.species >= 0) {
                a.body = body; a.fed = fed != 0; a.rng = (uint32_t)animals.size() * 7919u + 1u;
                a.netId = NetId_Allocate(NetKind::Animal); // a save from before #317; its "animal_id" line replaces it
                animals.push_back(a);
                lastLoaded = (int)animals.size() - 1;
                gPopulationCleanup = true; // (#1633)
            } else lastLoaded = -1;
        }
        return true;
    }
    if (strcmp(tok, "animal_id") == 0) {
        unsigned long long id = 0;
        if (fscanf(f, " %llu", &id) == 1 && id != 0 && lastLoaded >= 0 && lastLoaded < (int)animals.size()) {
            animals[lastLoaded].netId = (uint64_t)id;
            NetId_Note((uint64_t)id);
        }
        return true;
    }
    if (strcmp(tok, "animal_home") == 0) { // (#1633)
        Vector3 h{};
        if (fscanf(f, " %f %f %f", &h.x, &h.y, &h.z) == 3 && lastLoaded >= 0 && lastLoaded < (int)animals.size()) { animals[lastLoaded].home = h; animals[lastLoaded].hasHome = true; }
        return true;
    }
    if (strcmp(tok, "animal_tiles") == 0) {
        size_t n = 0;
        if (fscanf(f, " %127s %zu", body, &n) != 2) return true;
        auto& set = spawnedTiles[body];
        for (size_t i = 0; i < n; ++i) { unsigned long long k = 0; if (fscanf(f, " %llu", &k) != 1) break; set.insert((uint64_t)k); }
        return true;
    }
    return false;
}

void Animals_PathStats(int& plans, int& found, double& avgMs, double& worstMs) {
    plans = gPathPlans; found = gPathFound; avgMs = gPathPlans ? gPathMs / gPathPlans : 0.0; worstMs = gPathWorstMs;
}

// ---- #334 multiplayer: the server's animals, sent to clients; a client's copies ----
void Animals_CollectNet(Scene& world, std::vector<AnimalNet>& out) {
    LoadSpecies();
    out.clear();
    std::unordered_map<std::string, uint64_t> bodyIds; // name -> NetId, looked up once a call
    for (const Animal& a : animals) {
        auto it = bodyIds.find(a.body);
        if (it == bodyIds.end()) {
            uint64_t id = 0;
            if (Planetoid* p = BodyNamed(world, a.body)) id = p->netId;
            it = bodyIds.emplace(a.body, id).first;
        }
        if (!it->second) continue;
        const Species& s = species[a.species];
        AnimalNet n;
        n.netId = a.netId; n.species = (uint8_t)a.species; n.body = it->second; n.pos = a.pos.toFloat(); n.heading = a.heading;
        n.state = (uint8_t)a.state; n.health = a.health; n.age = a.age;
        n.flags = (uint8_t)((s.hostile ? 1 : 0) | (a.angryUntil > WorldTime() ? 2 : 0) | (a.windUpAt >= 0.0 ? 4 : 0));
        out.push_back(n);
    }
}

int Animals_FindNet(uint64_t netId) {
    for (size_t i = 0; i < animals.size(); ++i) if (animals[i].netId == netId) return (int)i;
    return -1;
}

void Animals_ApplyNet(const std::vector<AnimalNet>& in, const std::vector<uint64_t>& removed) {
    LoadSpecies();
    const double now = GetTime();
    for (uint64_t id : removed) { const int i = Animals_FindNet(id); if (i >= 0) { animals.erase(animals.begin() + i); if (aimed == i) aimed = -1; else if (aimed > i) --aimed; } }
    for (const AnimalNet& n : in) {
        if (n.species >= species.size()) continue;
        int i = Animals_FindNet(n.netId);
        if (i < 0) { // first sight: appear where it is
            Planetoid* b = Planetoid::findNet(n.body);
            if (!b) continue;
            Animal a;
            a.netId = n.netId; a.species = n.species; a.body = b->name; a.netPos = n.pos; a.pos = a.prevPos = DVec3(n.pos); a.heading = a.netHeading = a.prevHeading = n.heading;
            a.prevValid = true;
            animals.push_back(a);
            i = (int)animals.size() - 1;
        }
        Animal& a = animals[i];
        a.netPos = n.pos; a.netHeading = n.heading; a.netSeen = now;
        a.state = (State)std::min<uint8_t>(n.state, (uint8_t)State::Hunt); a.health = n.health; a.age = n.age;
        a.angryUntil = (n.flags & 2) ? WorldTime() + 1.0 : -1.0;
        a.windUpAt = (n.flags & 4) ? simClock : -1.0;
    }
}

bool Animals_HitNet(uint64_t netId, float damage, std::vector<std::pair<std::string, int>>& drops) {
    DamageEvent e; e.amount = damage; e.kind = DamageKind::Melee;
    return Animals_HitNetEvent(netId, e, drops, true);
}
bool Animals_HitNetEvent(uint64_t netId, const DamageEvent& e, std::vector<std::pair<std::string, int>>& drops, bool raw, float* dealtOut) {
    const int i = Animals_FindNet(netId);
    if (i < 0) return false;
    Animal& a = animals[i];
    const float damage = raw ? e.amount : TypedAmount(a, e); // (#1509) a client's hit as a local one: resistances, poise, its effect
    if (dealtOut) *dealtOut = damage;
    a.health -= damage;
    a.fleeTimer = !species[a.species].hostile || a.health < 0.35f * species[a.species].health ? 5.0f : 0.0f;
    if (species[a.species].hostile && a.fleeTimer <= 0.0f) a.angryUntil = WorldTime() + 20.0;
    a.fleeFrom = a.pos.toFloat();
    a.state = State::Flee;
    if (a.health <= 0.0f) { Kill(i, &drops); return true; }
    return false;
}
float Animals_HealthOf(int index) { return index >= 0 && index < (int)animals.size() ? animals[index].health : 0.0f; }

bool Animals_UseNet(uint64_t netId, int action, Inventory& inv, std::string& why, std::map<std::string, double>& before) {
    const int i = Animals_FindNet(netId);
    if (i < 0) { why = "it isn't there"; return false; }
    Animal& a = animals[i];
    const Species& s = species[a.species];
    const double now = WorldTime();
    for (const std::string& food : s.diet) before[food] = inv.amount(food);
    if (!s.product.empty() && s.product != "-") before[s.product] = inv.amount(s.product);
    if (action == 1) { // feed: one of its foods from the player's inventory
        for (const std::string& food : s.diet)
            if (inv.count(food) > 0 && inv.remove(food, 1)) { a.fed = true; a.followUntil = now + 30.0; a.fleeTimer = 0.0f; return true; }
        why = "no food it eats"; return false;
    }
    if (action == 2) { // its product (wool, milk, eggs)
        if (s.product.empty() || s.product == "-") { why = "nothing but meat"; return false; }
        if (a.age < 1.0f) { why = "too young"; return false; }
        if (now < a.productReady) { why = "not yet"; return false; }
        if (s.productTool != "-" && inv.count(s.productTool) <= 0) { why = "needs " + s.productTool; return false; }
        inv.add(s.product, s.product == "wool" ? 2 : 1);
        a.productReady = now + s.regrow;
        return true;
    }
    why = "unknown action"; return false;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
TEST(AnimalLagCompensation, AHitCountsWhereTheAnimalWasWithinTheRewind) { // (#1108)
    // A deer fleeing at 6 m/s along x, its trail a sample every sim step (60 Hz) for 0.3 s; the hunter stands at x = -4.4,
    // 4.4 m behind where the deer was 0.15 s ago, 5.3 m from where it is now (out of the 4.5 m reach).
    std::vector<std::pair<double, Vector3>> trail;
    const double now = 10.0;
    for (int k = 18; k >= 1; --k) { const double t = now - k / 60.0; trail.push_back({ t, Vector3{ (float)(6.0 * (t - now)), 0.0f, 0.0f } }); }
    const Vector3 current{ 0.0f, 0.0f, 0.0f }, hunter{ -5.3f, 0.0f, 0.0f };
    EXPECT_NEAR(Animals_RewindDistance(trail, current, hunter, now, 0.0), 5.3f, 1e-4f) << "no rewind: where it is now";
    EXPECT_LT(Animals_RewindDistance(trail, current, hunter, now, 0.2), 4.5f) << "within 0.2 s it was in reach";
    EXPECT_GT(Animals_RewindDistance(trail, current, Vector3{ -7.0f, 0.0f, 0.0f }, now, 0.2), 4.5f) << "the rewind is capped: further back doesn't count";
}
#endif
// (#1216 tests) A training dummy: the animal held at `pos` (body-local, on the ground) facing `face` (body-local), at rest.
void Animals_TestPin(int i, Vector3 pos, Vector3 face) {
    if (i < 0 || i >= (int)animals.size()) return;
    Animal& a = animals[i];
    Vector3 east, north;
    TangentFrame(Vector3Normalize(pos), east, north);
    a.pos = a.prevPos = DVec3(pos);
    a.heading = a.prevHeading = atan2f(Vector3DotProduct(face, north), Vector3DotProduct(face, east));
    a.state = State::Rest; a.fleeTimer = 0.0f; a.path.clear();
}
