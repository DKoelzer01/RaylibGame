#include "animals.h"
#include "audio.h"
#include "buildtool.h"
#include "gameui.h"
#include "hotbar.h"
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

// The surface point under a body-local direction, and whether an animal may stand there.
bool Ground(const PlanetField& pf, Vector3 dir, Vector3& at, float& elevation) {
    elevation = pf.elevation(dir);
    at = Vector3Scale(dir, pf.radius + elevation);
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
            std::vector<int> fit;
            for (size_t i = 0; i < species.size(); ++i)
                if (lat >= species[i].latMin && lat <= species[i].latMax && species[i].herdMax > 0 && !species[i].night && !(species[i].hostile && worldDifficulty == 0))
                    fit.push_back((int)i);
            if (fit.empty()) continue;
            const Species& s = species[fit[(size_t)(Rand01(h) * fit.size()) % fit.size()]];
            if (s.hostile && Vector3Distance(Vector3Scale(centre, pf.radius), Vector3Scale(StartDir(body, pf), pf.radius)) < 500.0f) continue; // a safe start
            Vector3 at; float e;
            if (!Ground(pf, centre, at, e) || e > 900.0f) continue; // no herds at sea or above the tree line
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
        if (!Ground(*pf, d, at, e)) continue;
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
        if (s.hostile && threatsOn && Difficulty_ThreatScale() > 0.0f && sameBody && a.fleeTimer <= 0.0f && !survival.dead && !player.flyMode) {
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
        } else if (Vector3DotProduct(sun, body->dirToWorld(up)) < -0.1f) {
            a.state = State::Rest; // night
        } else {
            a.timer -= step;
            if (a.timer <= 0.0f) {
                const bool walk = Rand01(a.rng) < 0.55f;
                a.state = walk ? State::Walk : State::Graze;
                if (walk) a.heading += (Rand01(a.rng) - 0.5f) * 2.2f;
                a.timer = walk ? 2.0f + 5.0f * Rand01(a.rng) : 2.0f + 6.0f * Rand01(a.rng);
            }
            if (a.state == State::Walk) speed = s.walk;
        }
        if (steer) {
            const float x = Vector3DotProduct(toward, east), y = Vector3DotProduct(toward, north);
            if (x * x + y * y > 1e-6f) a.heading = atan2f(y, x);
        }
        if (speed <= 0.0f) continue;
        // Move along the surface; turn away from water, cliffs and build parts.
        const Vector3 fwd = Vector3Add(Vector3Scale(east, cosf(a.heading)), Vector3Scale(north, sinf(a.heading)));
        const float stepLen = speed * step * (0.6f + 0.4f * a.age);
        const Vector3 nd = Vector3Normalize(Vector3Add(a.pos, Vector3Scale(fwd, stepLen)));
        Vector3 at; float e;
        const float height = s.size.y * (0.6f + 0.4f * a.age);
        const bool ok = Ground(*field, nd, at, e) &&
                        fabsf(Vector3Length(at) - Vector3Length(a.pos)) <= stepLen * 0.84f + 0.05f && // < 40 degrees
                        !BlockedAlong(*body, a.pos, at, height);
        if (ok) a.pos = at;
        else { a.heading += PI * (0.6f + 0.6f * Rand01(a.rng)); a.timer = 1.0f + Rand01(a.rng); }
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
            animals.push_back(young);
            if (here && here->name == a.body) GameUI_Toast(TextFormat("A young %s was born", species[a.species].name.c_str()));
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
        const float light = 0.45f + 0.55f * std::clamp(Vector3DotProduct(sun, u) * 2.0f + 0.5f, 0.0f, 1.0f);
        auto shade = [&](Color c, float m) { return Color{ (unsigned char)(c.r * light * m), (unsigned char)(c.g * light * m), (unsigned char)(c.b * light * m), 255 }; };
        Matrix m = { r.x, u.x, f.x, base.x, r.y, u.y, f.y, base.y, r.z, u.z, f.z, base.z, 0, 0, 0, 1 };
        rlPushMatrix();
        rlMultMatrixf(MatrixToFloat(m));
        rlScalef(k, k, k);
        const Vector3 sz = s.size; // length (forward), height, width
        const float legH = sz.y * 0.35f;
        DrawCube({ 0, legH + (sz.y - legH) * 0.5f, 0 }, sz.z, sz.y - legH, sz.x, shade(s.body, 1.0f));                     // body
        DrawCubeWires({ 0, legH + (sz.y - legH) * 0.5f, 0 }, sz.z, sz.y - legH, sz.x, shade(s.body, 0.6f));
        const float hs = std::min(sz.y, sz.z) * 0.55f;
        DrawCube({ 0, sz.y * 0.8f, sz.x * 0.5f + hs * 0.4f }, hs * 0.8f, hs, hs, shade(s.head, 1.0f));                // head
        for (int l = 0; l < 4; ++l)                                                                                       // legs
            DrawCube({ (l & 1 ? 0.3f : -0.3f) * sz.z, legH * 0.5f, (l & 2 ? 0.35f : -0.35f) * sz.x }, sz.z * 0.15f, legH, sz.z * 0.15f, shade(s.body, 0.75f));
        if ((int)i == aimed) DrawCubeWires({ 0, sz.y * 0.5f, 0 }, sz.z * 1.15f, sz.y * 1.1f, sz.x * 1.15f, YELLOW);
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
