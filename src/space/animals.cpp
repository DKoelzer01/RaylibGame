#include "animals.h"
#include "instancebatch.h"
#include "audio.h"
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
#include <chrono>
#include <sstream>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
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
            while (ms >> w) { if (w == "biped") s.biped = true; else if (w == "wool") s.wool = true; else if (w == "horns") s.horns = true; else if (w == "spikes") s.spikes = true; else if (w == "beak") s.beak = true; else if (w == "tail") s.tail = true; else if (w == "ears") s.ears = true; else if (w == "long") s.longBody = true; }
        }
        species.push_back(s);
    }
    fclose(f);
}
int SpeciesIndex(const std::string& id) { LoadSpecies(); for (size_t i = 0; i < species.size(); ++i) if (species[i].id == id) return (int)i; return -1; }

// ---------------------------------------------------------------------------------------------
// Animals
// ---------------------------------------------------------------------------------------------

enum class State : uint8_t { Graze, Walk, Flee, Follow, Rest, Hunt };

struct Animal {
    int species = 0;
    std::string body;            // the body's name (saves); its pointer is looked up by id each frame
    Vector3 pos{};               // body-local, on the surface
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
};
std::vector<Animal> animals;
std::unordered_map<std::string, std::unordered_set<uint64_t>> spawnedTiles; // body -> herd tiles decided
int aimed = -1;
double lastSpawnCheck = -1.0, lastHit = -1.0, lastPenCheck = -10.0, lastFireScan = -10.0, lastNightSpawn = -10.0;
std::vector<Vector3> fires; // campfires near the player (body-local)
bool threatsOn = true;     // off in bench runs except --test-threats
bool NearFire(Vector3 p, float r, Vector3& at) { for (const Vector3& f : fires) if (Vector3Distance(f, p) < r) { at = f; return true; } return false; }
int pennedNow = 0;
double lastUpdateMs = 0.0;
int frameNo = 0;

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

// Solid enough to block an animal (a wall, an overhang): ~0.4 m into rock. Thin near-iso films hover just over the
// ground in places (density within ~0.06 of iso); they must not stop a walk.
constexpr float kRockMargin = 0.1f;

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
        if (d >= ISO_LEVEL && prevD < ISO_LEVEL) { // air above, solid here (a start inside a thin film looks on down)
            const float t = std::clamp((ISO_LEVEL - prevD) / std::max(d - prevD, 1e-6f), 0.0f, 1.0f);
            at = Vector3Scale(dir, prevR + (r - prevR) * t);
            found = true;
            break;
        }
        prevR = r; prevD = d;
    }
    if (!found) return false;
    elevation = Vector3Length(at) - pf.radius;
    if (elevation < 0.3f) return false; // the sea
    if (pf.inlandWaterNear(dir, 0.0f)) {
        const float w = pf.waterElevation(dir);
        if (w != PlanetField::kNoWater && w > elevation + 0.3f) return false; // a river or pond deeper than 0.3 m
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

// An animal at `pos` is penned when build parts stand within 30 m in all eight directions.
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
            const Vector3 at = Vector3Scale(Vector3Normalize(Vector3Add(pos, Vector3Scale(dirT, s))), r);
            hit = BlockedByParts(body, at, Vector3Normalize(at), 1.0f);
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
                if (livesHere(species[i]) && species[i].herdMax > 0 && !species[i].night && !(species[i].hostile && worldDifficulty == 0))
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

void Kill(int i) {
    const Animal& a = animals[i];
    const Species& s = species[a.species];
    QuestBook_Count("slain");
    if (s.night) QuestBook_Count("stalker");
    std::string got;
    for (const auto& d : s.drops) {
        const int n = std::max(1, (int)lroundf(d.second * (0.4f + 0.6f * a.age)));
        if (itemDB.find(d.first)) { inventory.add(d.first, n); got += TextFormat("  +%d %s", n, d.first.c_str()); }
    }
    GameUI_Toast(TextFormat("The %s fell.%s", s.name.c_str(), got.c_str()));
    animals.erase(animals.begin() + i);
    if (aimed == i) aimed = -1;
}

} // namespace

// ---------------------------------------------------------------------------------------------

void Animals_Reset() { animals.clear(); spawnedTiles.clear(); aimed = -1; }
int Animals_Count() { return (int)animals.size(); }
Vector3 Animals_Position(int i) { return i >= 0 && i < (int)animals.size() ? animals[i].pos : Vector3{ 0, 0, 0 }; }
double Animals_LastUpdateMs() { return lastUpdateMs; }
int Animals_PennedCount() { return pennedNow; }
void Animals_SetThreats(bool on) { threatsOn = on; }
void Animals_RefreshNow() { lastFireScan = -100.0; lastNightSpawn = -100.0; lastPenCheck = -100.0; }
int Animals_NightCount() { int n = 0; for (const Animal& a : animals) n += species[a.species].night ? 1 : 0; return n; }
bool Animals_IsFeed(const std::string& id) { LoadSpecies(); for (const Species& s : species) for (const std::string& d : s.diet) if (d == id) return true; return false; }
const std::string& Animals_Body(int i) { static const std::string none; return i >= 0 && i < (int)animals.size() ? animals[i].body : none; }
bool Animals_Wild(int i) { return i >= 0 && i < (int)animals.size() && species[animals[i].species].wild; }
bool Animals_Hostile(int i) { return i >= 0 && i < (int)animals.size() && species[animals[i].species].hostile; }
bool Animals_Aiming() { return aimed >= 0; }
int Animals_YoungCount() { int n = 0; for (const Animal& a : animals) n += a.age < 1.0f ? 1 : 0; return n; }

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
        a.pos = at;
        a.heading = 2.0f * PI * Rand01(h);
        a.health = species[si].health;
        a.age = young ? 0.0f : 1.0f;
        a.rng = h ^ (uint32_t)(animals.size() * 7919u);
        a.timer = 1.0f + 4.0f * Rand01(h);
        a.productReady = WorldTime();
        animals.push_back(a);
        ++made;
    }
    return made;
}

int Animals_Nearest(Vector3 local, float range) {
    int best = -1;
    float bestD = range;
    for (size_t i = 0; i < animals.size(); ++i) {
        const float d = Vector3Distance(animals[i].pos, local);
        if (d < bestD) { bestD = d; best = (int)i; }
    }
    return best;
}

bool Animals_Hit(int i, float damage) {
    if (i < 0 || i >= (int)animals.size()) return false;
    Animal& a = animals[i];
    a.health -= damage;
    // Hostile ones fight back until badly hurt; the rest run.
    a.fleeTimer = !species[a.species].hostile || a.health < 0.35f * species[a.species].health ? 5.0f : 0.0f;
    if (species[a.species].hostile && a.fleeTimer <= 0.0f) a.angryUntil = WorldTime() + 20.0;
    a.fleeFrom = a.pos; // refined to the player's position in the update
    a.state = State::Flee;
    Audio_Play(Sfx::Dig, 0.5f, 1.6f);
    if (a.health <= 0.0f) { Kill(i); return true; }
    return false;
}

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

void Animals_Update(Scene& world, float dt, float worldDt, bool inputEnabled) {
    const auto t0 = std::chrono::steady_clock::now();
    LoadSpecies();
    ++frameNo;
    Planetoid* here = player.nearest;
    const PlanetField* pf = here ? dynamic_cast<const PlanetField*>(here->field.get()) : nullptr;
    Vector3 playerLocal{};
    if (here) playerLocal = here->toLocal((player.position - world.renderOrigin).toFloat());
    if (here && pf && !player.flyMode && GetTime() - lastSpawnCheck > 1.0) {
        lastSpawnCheck = GetTime();
        const float d = Vector3Length(playerLocal);
        if (d > 1.0f && d - pf->radius < 400.0f) SpawnHerds(world, *here, *pf, Vector3Scale(playerLocal, 1.0f / d));
    }
    const double now = WorldTime();
    const std::string held = Hotbar_HeldFood(); // its food in the player's hand: animals follow
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    // Campfires near the player (refreshed every second): night hunters keep away from them.
    if (here && GetTime() - lastFireScan > 1.0) {
        lastFireScan = GetTime();
        fires.clear();
        if (here->plots)
            for (const auto& pl : here->plots->plots)
                for (const Part& p : pl->parts)
                    if (p.type == PartType::Campfire) {
                        const Vector3 c = pl->toBody({ (p.origin.x + 1.0f) * kCellSize, (p.origin.y + 0.5f) * kCellSize, (p.origin.z + 1.0f) * kCellSize });
                        if (Vector3Distance(c, playerLocal) < 200.0f) fires.push_back(c);
                    }
    }
    // Night threats (not on peaceful): every 12 s after dark, a dusk stalker comes out of the dark 45-60 m
    // away, up to 2 near the player (3 on hard) - not underground, flying, or by a campfire.
    if (here && pf && threatsOn && !player.flyMode && !survival.dead && Difficulty_ThreatScale() > 0.0f && GetTime() - lastNightSpawn > 12.0) {
        lastNightSpawn = GetTime();
        const Vector3 up = Vector3Normalize(playerLocal);
        const bool dark = Vector3DotProduct(Vector3Normalize(world.sunDirection), here->dirToWorld(up)) < -0.1f;
        const bool surface = Vector3Length(playerLocal) > pf->radius + pf->elevation(up) - 6.0f;
        int hunter = -1;
        for (size_t k = 0; k < species.size(); ++k) if (species[k].night) { hunter = (int)k; break; }
        Vector3 fire;
        if (dark && surface && hunter >= 0 && !NearFire(playerLocal, 10.0f, fire)) {
            int near = 0;
            for (const Animal& a : animals) near += species[a.species].night && a.body == here->name && Vector3Distance(a.pos, playerLocal) < 100.0f ? 1 : 0;
            if (near < (worldDifficulty == 2 ? 3 : 2)) {
                Vector3 east, north;
                TangentFrame(up, east, north);
                for (int attempt = 0; attempt < 8; ++attempt) { // a dry spot (not water, not by a fire)
                const float ang = 2.0f * PI * (float)GetRandomValue(0, 1000) / 1000.0f, r = 45.0f + 15.0f * (float)GetRandomValue(0, 1000) / 1000.0f;
                const Vector3 d = Vector3Normalize(Vector3Add(playerLocal, Vector3Add(Vector3Scale(east, r * cosf(ang)), Vector3Scale(north, r * sinf(ang)))));
                const Vector3 at = Vector3Scale(d, pf->radius + pf->elevation(d));
                if (!NearFire(at, 14.0f, fire) && Animals_Spawn(world, *here, species[hunter].id, at, 1, false, 0.5f) > 0) break;
                }
            }
        }
    }
    // Ecology (Plan 5 A), twice a second over the animals near the player: predators (the hostile kinds)
    // pick the nearest prey within 30 m to chase (chase only: no kills yet, per the user); prey that spots a
    // predator within 18 m bolts; each herd's centre (same kind within 40 m) for cohesion.
    static double ecoAt = 0.0;
    if (here && GetTime() > ecoAt) {
        ecoAt = GetTime() + 0.5;
        std::vector<int> near;
        for (size_t i = 0; i < animals.size(); ++i)
            if (animals[i].body == here->name && Vector3Distance(animals[i].pos, playerLocal) < 300.0f) near.push_back((int)i);
        for (int i : near) {
            Animal& a = animals[i];
            const Species& s = species[a.species];
            a.chase = -1;
            Vector3 centre{}; int herd = 0;
            float bestPrey = 30.0f, bestThreat = 18.0f;
            for (int j : near) {
                if (j == i) continue;
                const Animal& b = animals[j];
                const float d = Vector3Distance(a.pos, b.pos);
                const Species& t = species[b.species];
                if (b.species == a.species && d < 40.0f) { centre = Vector3Add(centre, b.pos); ++herd; }
                if (s.hostile && !t.hostile && d < bestPrey && threatsOn) { bestPrey = d; a.chase = j; }
                if (!s.hostile && t.hostile && d < bestThreat && threatsOn) { bestThreat = d; a.fleeTimer = std::max(a.fleeTimer, 3.0f); a.fleeFrom = b.pos; }
            }
            a.herdCentre = herd ? Vector3Scale(centre, 1.0f / herd) : a.pos;
            a.herdSize = herd;
        }
    }
    std::vector<int> gone;
    std::unordered_map<std::string, Planetoid*> bodies;
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
        // Only near the player is it simulated; past 80 m a quarter as often.
        const bool sameBody = body == here;
        const float dist = sameBody ? Vector3Distance(a.pos, playerLocal) : 1e9f;
        if (dist > 300.0f) continue;
        float step = dt;
        if (dist > 80.0f) { if ((frameNo + (int)i) % 4) continue; step = dt * 4.0f; }
        if ((a.settleIn -= step) <= 0.0f) {
            a.settleIn = 1.0f;
            // Settle (#223): ground dug out from under it drops it onto what is left (up to 3 m); inside rock
            // (an older save, ground built over it) or over a deeper hole, it goes to the surface found from above.
            const Vector3 d = Vector3Normalize(a.pos);
            Vector3 g; float ge;
            if (Ground(*body, *field, d, g, ge, Vector3Length(a.pos), 0.6f, 3.0f) || Ground(*body, *field, d, g, ge)) a.pos = g;
        }
        const Vector3 up = Vector3Normalize(a.pos);
        Vector3 east, north;
        TangentFrame(up, east, north);
        // Decide.
        float speed = 0.0f;
        Vector3 toward{};
        bool steer = false;
        const bool dark = Vector3DotProduct(sun, body->dirToWorld(up)) < -0.1f;
        // Night hunters: in daylight they run off, and are gone once out of sight (or after 40 s).
        if (s.night && !dark) {
            a.daylight += step;
            if (a.daylight > 40.0f || (dist > 70.0f && a.daylight > 5.0f)) { gone.push_back((int)i); continue; }
            a.fleeTimer = std::max(a.fleeTimer, 1.0f);
            a.fleeFrom = sameBody ? playerLocal : a.pos;
        }
        // Hostile ones (not on peaceful): thornbacks charge within 12 m (or for 20 s after being hit),
        // dusk stalkers hunt within 70 m after dark. A campfire within 10 m keeps the player safe, and
        // stalkers will not come within 14 m of one.
        bool hunting = false;
        if (s.hostile && threatsOn && Difficulty_ThreatScale() > 0.0f && sameBody && a.fleeTimer <= 0.0f && !survival.dead && !player.flyMode && !Rocket_Boarded()) { // not inside a cockpit
            Vector3 fire;
            const bool playerSafe = NearFire(playerLocal, 10.0f, fire);
            if (s.night && NearFire(a.pos, 14.0f, fire)) { a.fleeTimer = 3.0f; a.fleeFrom = fire; }
            else if (s.night) hunting = dark && dist < 70.0f && !playerSafe;
            else hunting = (dist < 12.0f || now < a.angryUntil) && dist < 30.0f && !playerSafe;
        }
        if (hunting) {
            a.state = State::Hunt;
            toward = Vector3Subtract(playerLocal, a.pos);
            steer = true;
            speed = s.run * (s.night ? 0.8f : 1.0f);
            const float reach = 0.5f * s.size.x + 1.0f;
            if (dist < reach) {
                speed = 0.0f;
                if (GetTime() - a.lastBite > 1.2) {
                    a.lastBite = GetTime();
                    survival.health -= 8.0f * Difficulty_ThreatScale();
                    Audio_Play(Sfx::Error, 0.7f, 0.5f);
                    GameUI_Toast(TextFormat("The %s %s you!", s.name.c_str(), s.night ? "bites" : "rams"));
                }
            }
        } else if (a.fleeTimer > 0.0f) {
            a.fleeTimer -= step;
            a.state = State::Flee;
            if (sameBody) a.fleeFrom = playerLocal;
            toward = Vector3Subtract(a.pos, a.fleeFrom);
            steer = true;
            speed = s.run;
        } else if (sameBody && dist > 2.5f && (now < a.followUntil || (dist < 16.0f && !held.empty() && std::find(s.diet.begin(), s.diet.end(), held) != s.diet.end()))) {
            // Fed (30 s), or the player holds its food: it follows (wild ones too: the food wins).
            a.state = State::Follow;
            toward = Vector3Subtract(playerLocal, a.pos);
            steer = true;
            speed = std::min(s.run, s.walk * 2.0f);
        } else if (s.wild && sameBody && dist < 25.0f) {
            a.fleeTimer = 4.0f;
            a.fleeFrom = playerLocal;
            continue;
        } else if (s.hostile && a.chase >= 0 && a.chase < (int)animals.size() && animals[a.chase].body == a.body &&
                   Vector3Distance(animals[a.chase].pos, a.pos) < 35.0f && !(s.night && !dark)) {
            // A predator after prey (Plan 5 A; chase only): it runs at it; close in, the prey bolts.
            a.state = State::Hunt;
            toward = Vector3Subtract(animals[a.chase].pos, a.pos);
            steer = true;
            speed = s.run * 0.75f;
            Animal& prey = animals[a.chase];
            if (Vector3Distance(prey.pos, a.pos) < 6.0f) { prey.fleeTimer = std::max(prey.fleeTimer, 3.0f); prey.fleeFrom = a.pos; }
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
            const float herdAway = a.herdSize > 0 ? Vector3Distance(a.herdCentre, a.pos) : 0.0f;
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
                if (herdAway > 12.0f) { toward = Vector3Subtract(a.herdCentre, a.pos); steer = true; } // back to the herd
            }
        }
        if (steer) {
            const float x = Vector3DotProduct(toward, east), y = Vector3DotProduct(toward, north);
            if (x * x + y * y > 1e-6f) a.heading = atan2f(y, x);
        }
        if (speed <= 0.0f) { a.moveSpeed = 0.0f; continue; }
        // Move along the surface; turn away from water, cliffs and build parts.
        const float stepLen = speed * step * (0.6f + 0.4f * a.age);
        Vector3 at; float e;
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
        auto tryStep = [&](float heading) {
            const Vector3 fwd = Vector3Add(Vector3Scale(east, cosf(heading)), Vector3Scale(north, sinf(heading)));
            const Vector3 nd = Vector3Normalize(Vector3Add(a.pos, Vector3Scale(fwd, stepLen)));
            return Ground(*body, *field, nd, at, e, Vector3Length(a.pos), reachUp, reachDown) &&
                   body->densityAt(Vector3Add(at, Vector3Scale(nd, height * 0.8f))) < ISO_LEVEL + kRockMargin &&
                   !BlockedAlong(*body, a.pos, at, height);
        };
        // Heading somewhere (the player, the herd, away from a threat) and blocked: go round (#223: real terrain has
        // ledges, pits and basins the height map hid). The detour is held for a while - following the obstacle on one
        // side, turning as far as backing off - before aiming at the target again, so a pocket or a long ledge
        // doesn't trap it pushing straight into the same edge.
        if (steer && a.detourTime > 0.0f) { a.detourTime -= step; a.heading = a.detourHeading; }
        else if (a.detourTime > 0.0f) a.detourTime = 0.0f;
        bool ok = tryStep(a.heading);
        if (ok && a.detourTime <= 0.0f) a.detourSide = 0.0f;
        if (!ok && steer) {
            if (a.detourSide == 0.0f) a.detourSide = Rand01(a.rng) < 0.5f ? 1.0f : -1.0f;
            for (float turn : { 0.5f, 1.0f, 1.6f, 2.2f, 3.1f })
                for (float side : { a.detourSide, -a.detourSide })
                    if (!ok && tryStep(a.heading + side * turn)) {
                        ok = true;
                        a.heading += side * turn;
                        a.detourSide = side;
                        a.detourHeading = a.heading;
                        a.detourTime = 1.5f;
                    }
        }
        if (ok) {
            a.pos = at;
            a.moveSpeed = speed;
            a.gait += stepLen / std::max(0.25f, s.size.x * 0.45f) * PI; // one stride per ~0.9 body lengths
            a.stuckTime = 0.0f;
        } else {
            // Blocked (water, a cliff, a wall): turn away; stuck a while, turn harder and pick a new activity.
            a.moveSpeed = 0.0f;
            a.stuckTime += step;
            a.heading += PI * (0.6f + 0.6f * Rand01(a.rng)) * (a.stuckTime > 2.0f ? 1.5f : 1.0f);
            a.timer = 1.0f + Rand01(a.rng);
            if (a.stuckTime > 4.0f) { a.state = State::Graze; a.stuckTime = 0.0f; }
        }
    }
    // Penned: build parts within 30 m in all eight directions around an animal near the player
    // (checked every 2 s; the quest's "pen two animals").
    if (here && GetTime() - lastPenCheck > 2.0) {
        lastPenCheck = GetTime();
        int n = 0;
        for (const Animal& a : animals) {
            if (a.body != here->name || Vector3Distance(a.pos, playerLocal) > 120.0f) continue;
            if (IsPenned(*here, a.pos)) ++n;
        }
        pennedNow = n;
        QuestBook_AtLeast("penned", n);
    }
    // Night hunters gone at daylight.
    for (int k = (int)gone.size() - 1; k >= 0; --k) animals.erase(animals.begin() + gone[k]);
    if (!gone.empty()) aimed = -1;
    // Breeding: two fed adults of a kind within 6 m make a young one (once every 5 min of world time).
    for (size_t i = 0; i < animals.size(); ++i) {
        Animal& a = animals[i];
        if (!a.fed || a.age < 1.0f || now - a.lastBred < 300.0) continue;
        for (size_t j = i + 1; j < animals.size(); ++j) {
            Animal& b = animals[j];
            if (b.species != a.species || b.body != a.body || !b.fed || b.age < 1.0f || now - b.lastBred < 300.0) continue;
            if (Vector3Distance(a.pos, b.pos) > 6.0f) continue;
            a.fed = b.fed = false;
            a.lastBred = b.lastBred = now;
            Animal young = a;
            young.pos = Vector3Scale(Vector3Add(a.pos, b.pos), 0.5f);
            young.age = 0.0f;
            young.health = species[a.species].health * 0.5f;
            young.followUntil = -1.0;
            young.rng = a.rng ^ 0x9e37u;
            const bool seen = here && here->name == a.body; // before the push: it may move the vector (a, b)
            const std::string kind = species[a.species].name;
            animals.push_back(young);
            if (seen) GameUI_Toast(TextFormat("A young %s was born", kind.c_str()));
            QuestBook_Count("bred");
            break;
        }
    }
    // Aim (within 5 m along the crosshair) and the player's actions.
    aimed = -1;
    if (here && !player.flyMode) {
        const Vector3 o = here->toLocal(world.camera.position);
        const Vector3 d = here->dirToLocal(Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position)));
        float best = 5.0f;
        for (size_t i = 0; i < animals.size(); ++i) {
            const Animal& a = animals[i];
            if (a.body != here->name) continue;
            const Species& s = species[a.species];
            const float scale = 0.6f + 0.4f * a.age;
            const Vector3 c = Vector3Add(a.pos, Vector3Scale(Vector3Normalize(a.pos), s.size.y * 0.5f * scale));
            const float r = 0.5f * std::max(s.size.x, s.size.y) * scale + 0.1f;
            const Vector3 oc = Vector3Subtract(c, o);
            const float t = Vector3DotProduct(oc, d);
            if (t < 0.0f || t > best) continue;
            if (Vector3Length(Vector3Subtract(oc, Vector3Scale(d, t))) > r) continue;
            best = t; aimed = (int)i;
        }
    }
    if (aimed >= 0 && inputEnabled && !buildTool.active) {
        Animal& a = animals[aimed];
        const Species& s = species[a.species];
        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && GetTime() - lastHit > 0.35) {
            lastHit = GetTime();
            const ToolSpec t = ToolOf(terrainTool.heldTool);
            ++terrainTool.swings;
            Animals_Hit(aimed, t.type == ToolType::Hands ? 2.0f : 3.0f + 3.0f * t.tier);
        } else if (IsKeyPressed(KEY_I)) {
            Animals_Feed(aimed);
        } else if (IsKeyPressed(KEY_O)) {
            if (s.product.empty() || s.product == "-") GameUI_Toast(TextFormat("A %s gives nothing but meat", s.name.c_str()));
            else if (a.age < 1.0f) GameUI_Toast("Too young");
            else if (now < a.productReady) GameUI_Toast(TextFormat("Not yet (%.0f s)", a.productReady - now));
            else if (s.productTool != "-" && inventory.count(s.productTool) <= 0) GameUI_Toast(TextFormat("You need %s", s.productTool.c_str()));
            else {
                inventory.add(s.product, s.product == "wool" ? 2 : 1);
                a.productReady = now + s.regrow;
                GameUI_Toast(TextFormat("+%s", s.product.c_str()));
                Audio_Play(Sfx::Pickup, 0.5f);
            }
        }
    }
    lastUpdateMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// The animal model (Plan 5 A): built from boxes in the animal's frame (x right, y up, z forward), proportioned
// from animals.txt (length, height, width) and its model features. Legs swing on the gait (trot: diagonal
// pairs; bipeds alternate) with knees bending on the back swing; the head dips to graze; resting folds the
// legs and lowers the body; the tail sways. Box() draws one part rotated about its own top (a joint).
static inline void DrawCubeB(Vector3 p, float w, float h, float l, Color c) { IB_Cube(p, { w, h, l }, c); } // batched, lit (instancebatch.h)
template <typename ShadeFn>
static void DrawAnimalModel(const Species& s, const Animal& a, ShadeFn shade) {
    const Vector3 sz = s.size;                       // length (z), height (y), width (x)
    const float t = (float)GetTime();
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
    const int legs = s.biped ? 2 : 4;
    const Color legC = shade(s.body, 0.8f);
    for (int l = 0; l < legs; ++l) {
        const bool right = l & 1, front = s.biped ? true : (l & 2) != 0;
        const float phase = s.biped ? (right ? PI : 0.0f) : ((right != front) ? PI : 0.0f);
        const float ang = resting ? 80.0f : swing * sinf(a.gait + phase) * RAD2DEG;
        const float knee = resting ? -150.0f : std::max(0.0f, -sinf(a.gait + phase)) * swing * 60.0f;
        const Vector3 hip = { (right ? 0.32f : -0.32f) * sz.z, bodyY + bodyH * 0.15f, s.biped ? 0.0f : (front ? 0.3f : -0.33f) * len };
        const float upper = legH * 0.55f, lowerLen = legH * 0.5f, w = sz.z * (s.biped ? 0.1f : 0.16f);
        rlPushMatrix();
        rlTranslatef(hip.x, hip.y, hip.z);
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
    if (animals.empty()) return;
    const Vector3 sun = Vector3Normalize(world.sunDirection);
    std::unordered_map<std::string, Planetoid*> bodies;
    for (size_t i = 0; i < animals.size(); ++i) {
        const Animal& a = animals[i];
        auto bit = bodies.find(a.body);
        Planetoid* body = bit != bodies.end() ? bit->second : (bodies[a.body] = BodyNamed(world, a.body));
        if (!body) continue;
        const Vector3 base = body->toWorld(a.pos);
        if (Vector3Distance(base, world.camera.position) > 160.0f) continue;
        const Species& s = species[a.species];
        const Vector3 lup = Vector3Normalize(a.pos);
        Vector3 east, north;
        TangentFrame(lup, east, north);
        const Vector3 f = body->dirToWorld(Vector3Add(Vector3Scale(east, cosf(a.heading)), Vector3Scale(north, sinf(a.heading))));
        const Vector3 u = body->dirToWorld(lup);
        const Vector3 r = Vector3CrossProduct(u, f); // right-handed (r, u, f)
        const float k = 0.6f + 0.4f * a.age;
        const float light = 1.0f; // lit by the instance batch shader (sun, shadows, headlamp)
        auto shade = [&](Color c, float m) { return Color{ (unsigned char)(c.r * light * m), (unsigned char)(c.g * light * m), (unsigned char)(c.b * light * m), 255 }; };
        Matrix m = { r.x, u.x, f.x, base.x, r.y, u.y, f.y, base.y, r.z, u.z, f.z, base.z, 0, 0, 0, 1 };
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        rlScalef(k, k, k);
        DrawAnimalModel(s, a, shade);
        if ((int)i == aimed) DrawCubeWires({ 0, s.size.y * 0.5f, 0 }, s.size.z * 1.15f, s.size.y * 1.1f, s.size.x * 1.15f, YELLOW);
        rlPopMatrix();
    }
}

void Animals_DrawHUD() {
    if (aimed < 0 || aimed >= (int)animals.size()) return;
    const Animal& a = animals[aimed];
    const Species& s = species[a.species];
    const int x = GetScreenWidth() / 2 + 20, y = GetScreenHeight() / 2 + 40;
    DrawRectangle(x - 8, y - 6, 380, 74, Fade(BLACK, 0.55f));
    const char* state = a.state == State::Hunt ? "hunting you" : a.state == State::Flee ? "fleeing" : a.state == State::Follow ? "follows you" : a.state == State::Rest ? "resting" : a.state == State::Walk ? "wandering" : "grazing";
    DrawText(TextFormat("%s%s  (%s%s)", a.age < 1.0f ? "young " : "", s.name.c_str(), state, a.fed ? ", fed" : ""), x, y, 20, Color{ 240, 220, 160, 255 });
    DrawRectangle(x, y + 26, 200, 6, Fade(RAYWHITE, 0.2f));
    DrawRectangle(x, y + 26, (int)(200 * std::max(0.0f, a.health) / s.health), 6, RED);
    std::string diet;
    for (const auto& d : s.diet) diet += (diet.empty() ? "" : "/") + d;
    DrawText(TextFormat("I feed (%s)   O %s   LMB attack", diet.c_str(), s.product.empty() || s.product == "-" ? "-" : s.product.c_str()), x, y + 40, 15, Fade(RAYWHITE, 0.75f));
}

void Animals_Save(FILE* f) {
    for (const Animal& a : animals)
        fprintf(f, "animal %s %s %f %f %f %f %f %f %d %.2f %.2f\n", species[a.species].id.c_str(), a.body.c_str(), a.pos.x, a.pos.y, a.pos.z,
                a.heading, a.health, a.age, a.fed ? 1 : 0, a.lastBred, a.productReady);
    for (const auto& kv : spawnedTiles) {
        fprintf(f, "animal_tiles %s %zu", kv.first.c_str(), kv.second.size());
        for (uint64_t k : kv.second) fprintf(f, " %llu", (unsigned long long)k);
        fprintf(f, "\n");
    }
}

bool Animals_LoadLine(const char* tok, FILE* f) {
    LoadSpecies();
    char id[64] = {}, body[128] = {};
    if (strcmp(tok, "animal") == 0) {
        Animal a;
        int fed = 0;
        if (fscanf(f, " %63s %127s %f %f %f %f %f %f %d %lf %lf", id, body, &a.pos.x, &a.pos.y, &a.pos.z, &a.heading, &a.health, &a.age, &fed, &a.lastBred, &a.productReady) == 11) {
            a.species = SpeciesIndex(id);
            if (a.species >= 0) { a.body = body; a.fed = fed != 0; a.rng = (uint32_t)animals.size() * 7919u + 1u; animals.push_back(a); }
        }
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
