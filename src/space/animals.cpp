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
#include "headless.h"
#include "net.h"
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
#include <chrono>
#include <sstream>
#include <cmath>
#include <cstring>
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
    uint64_t netId = 0;          // stable, saved (netid.h, #317)
    Vector3 prevPos{};           // pos and heading before the last simulation step: drawing blends them by SimAlpha (#315)
    float prevHeading = 0.0f;
    bool prevValid = false;
    int simFrame = -1;           // the update (frameNo) that last simulated it (#318 tests)
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
    double windUpAt = -1.0;      // an attack winding up since (simClock; #165), -1 none
    StatusState status;          // (#1215) burn, chill, bleed, poison
    PoiseState poise;            // (#1217) worn down by hits; broken, it staggers (set up on its first hit)
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
};
std::vector<Animal> animals;
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
constexpr double kPathFrameMs = 3.0, kPathTotalMs = 60.0;
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
        if ((S.expanded & 15) == 15 && std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() > kPathFrameMs) break;
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

void Kill(int i, std::vector<std::pair<std::string, int>>* dropsTo = nullptr) { // dropsTo (#334): a server's, for the hitting player
    const Animal& a = animals[i];
    const Species& s = species[a.species];
    if (!dropsTo) QuestBook_Count("slain"); // (a server's kill isn't its own player's quest)
    if (!dropsTo && s.night) QuestBook_Count("stalker");
    std::string got;
    for (const auto& d : s.drops) {
        const int n = std::max(1, (int)lroundf(d.second * (0.4f + 0.6f * a.age)));
        if (!itemDB.find(d.first)) continue;
        if (dropsTo) { dropsTo->push_back({ d.first, n }); continue; }
        inventory.add(d.first, n); got += TextFormat("  +%d %s", n, d.first.c_str());
    }
    if (!dropsTo) GameUI_Toast(TextFormat("The %s fell.%s", s.name.c_str(), got.c_str()));
    animals.erase(animals.begin() + i);
    if (aimed == i) aimed = -1;
}

} // namespace

// ---------------------------------------------------------------------------------------------

void Animals_Reset() { animals.clear(); spawnedTiles.clear(); aimed = -1; gSearch.active = false; gSearch.body = nullptr; } // (#902: the search held the old world's body)
int Animals_Count() { return (int)animals.size(); }
Vector3 Animals_Position(int i) { return i >= 0 && i < (int)animals.size() ? animals[i].pos : Vector3{ 0, 0, 0 }; }
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
bool Animals_Aiming() { return aimed >= 0 || hitLatch || Animals_MeleeBusy(); }
int Animals_AimedIndex() { return aimed >= 0 && aimed < (int)animals.size() ? aimed : -1; }
std::pair<std::string, std::string> Animals_SpeciesOf(int i) {
    LoadSpecies();
    if (i < 0 || i >= (int)animals.size() || animals[i].species < 0 || animals[i].species >= (int)species.size()) return {};
    const Species& s = species[animals[i].species];
    return { s.id, s.name };
}
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
        a.pos = at;
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

// (#1215) A typed hit: through the species' resistances and armour (the 10% floor), its status effect, then the hit itself.
// RAYLIBGAME_COMBAT1=0: the amount as it is (today).
bool Animals_HitEvent(int i, const DamageEvent& e, float* dealtOut) {
    if (i < 0 || i >= (int)animals.size()) return false;
    Animal& a = animals[i];
    float amount = e.amount;
    // (#1217) Poise: the hit wears it (its poise damage; none given: its amount); broken, the creature staggers (its attack
    // cancelled), and the first hit while it staggers does x1.5. Small creatures stagger on any heavy hit.
    if (Combat3_On()) {
        const Species& sp = species[a.species];
        if (a.poise.max <= 0.0f) Poise_Init(a.poise, Poise_CreatureMax(sp.health));
        const bool small = std::max(sp.size.x, std::max(sp.size.y, sp.size.z)) * (0.6f + 0.4f * a.age) <= PoiseTune::kSmallSize;
        const PoiseHit ph = Poise_Hit(a.poise, e.poise > 0.0f ? e.poise : e.amount, e.heavy, small);
        amount *= ph.damageScale;
        if (ph.staggered) a.windUpAt = -1.0;
    }
    if (Combat1_On()) {
        amount = Damage_Apply(amount, e.type, species[a.species].resist);
        Status_Apply(a.status, Status_FromHit(e.type, e.status, e.heavy, e.fromBehind, e.kind == DamageKind::Bite));
    }
    if (dealtOut) *dealtOut = amount;
    return Animals_Hit(i, amount);
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

void Animals_Update(Scene& world, float dt, float worldDt, bool inputEnabled) {
    const auto t0 = std::chrono::steady_clock::now();
    simClock += dt; // the attack wind-up and the bite cooldown run on it (#165): tests step animals many times a frame
    LoadSpecies();
    ++frameNo;
    simulatedNow = 0;
    for (Animal& a : animals) { a.prevPos = a.pos; a.prevHeading = a.heading; a.prevValid = true; } // the pose drawing blends from (#315)
    if (gAnimalsReplicated) { // #334: a client - the server runs the animals; each eases toward its latest state, nothing decides here
        const double now = GetTime();
        const float k = 1.0f - expf(-dt / 0.08f); // ~80 ms to close a gap (states come 10 a second)
        for (size_t i = 0; i < animals.size();) {
            Animal& a = animals[i];
            if (now - a.netSeen > 3.0) { animals.erase(animals.begin() + i); if (aimed == (int)i) aimed = -1; continue; } // out of view (or the removal was lost)
            const Vector3 before = a.pos;
            a.pos = Vector3Lerp(a.pos, a.netPos, k);
            float dh = a.netHeading - a.heading; while (dh > PI) dh -= 2 * PI; while (dh < -PI) dh += 2 * PI;
            a.heading += dh * k;
            const float step = Vector3Distance(before, a.pos);
            a.moveSpeed = dt > 0.0f ? step / dt : 0.0f;
            if (step > 1e-4f) a.gait += step / std::max(0.25f, species[a.species].size.x * 0.45f) * PI; // the walk cycle, as the AI drives it
            ++i;
        }
        lastUpdateMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
        return;
    }
    Planetoid* here = player.nearest;
    const PlanetField* pf = here ? dynamic_cast<const PlanetField*>(here->field.get()) : nullptr;
    Vector3 playerLocal{};
    if (here) playerLocal = here->toLocal((player.position - world.renderOrigin).toFloat());
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
            for (const Animal& a : animals) near += species[a.species].night && a.body == t.body->name && Vector3Distance(a.pos, t.local) < 100.0f ? 1 : 0;
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
    // Ecology (Plan 5 A), twice a second over the animals near the player: predators (the hostile kinds)
    // pick the nearest prey within 30 m to chase (chase only: no kills yet, per the user); prey that spots a
    // predator within 18 m bolts; each herd's centre (same kind within 40 m) for cohesion.
    static double ecoAt = 0.0;
    if (here && SimTime() > ecoAt) {
        ecoAt = SimTime() + 0.5;
        std::vector<int> near;
        for (size_t i = 0; i < animals.size(); ++i)
            if (animals[i].body == here->name && nearestCentre(here, animals[i].pos) < 300.0f) near.push_back((int)i);
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
    static const bool pathsOn = !(getenv("RAYLIBGAME_ANIMAL_PATHS") && getenv("RAYLIBGAME_ANIMAL_PATHS")[0] == '0'); // A/B
    int plansLeft = pathsOn ? 1 : 0; // route plans this update (#253)
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
        // Only near an interest centre (the player, other players: #318) is it simulated; past 80 m a quarter as often.
        const float dist = nearestCentre(body, a.pos);
        // (#334) the nearest player on its body: what it flees, follows, hunts (none: it keeps to itself)
        const Tgt* tg = nullptr; float pdist = 1e9f;
        for (const Tgt& t : targets) if (t.body == body) { const float d = Vector3Distance(a.pos, t.local); if (d < pdist) { pdist = d; tg = &t; } }
        const bool sameBody = tg != nullptr;
        const Vector3 tgLocal = tg ? tg->local : Vector3{};
        if (dist > 300.0f) continue;
        ++simulatedNow; a.simFrame = frameNo;
        float step = dt;
        if (dist > 80.0f) { if ((frameNo + (int)i) % 4) continue; step = dt * 4.0f; }
        if (Combat1_On() && a.status.any()) { // (#1215) burn, bleed and poison over time (creatures on land: no water to douse)
            a.health -= Status_Step(a.status, step, false, s.resist);
            if (a.health <= 0.0f) { Kill((int)i); --i; continue; }
        }
        if ((a.settleIn -= step) <= 0.0f) {
            a.settleIn = 1.0f;
            // Settle (#223): ground dug out from under it drops it onto what is left (up to 3 m); inside rock
            // (an older save, ground built over it) or over a deeper hole, it goes to the surface found from above.
            const Vector3 d = Vector3Normalize(a.pos);
            // Still standing on the ground (solid just under its feet, open 0.3 m over them): nothing to do. Two density
            // samples instead of a probe of 15 (and its fallback of ~180 on rough ground): the eroded field costs more a
            // sample, and every animal checked once a second (#476).
            const bool standing = body->densityAt(Vector3Subtract(a.pos, Vector3Scale(d, 0.15f))) >= ISO_LEVEL &&
                                  body->densityAt(Vector3Add(a.pos, Vector3Scale(d, 0.3f))) < ISO_LEVEL;
            Vector3 g; float ge;
            if (!standing && (Ground(*body, *field, d, g, ge, Vector3Length(a.pos), 0.6f, 3.0f) || Ground(*body, *field, d, g, ge))) a.pos = g;
        }
        const Vector3 up = Vector3Normalize(a.pos);
        Vector3 east, north;
        TangentFrame(up, east, north);
        // (#1217) Staggered: its attack is off and it stands where it is until the stagger ends; poise comes back after a rest.
        if (Combat3_On() && a.poise.max > 0.0f) {
            Poise_Step(a.poise, step);
            if (Poise_Staggered(a.poise)) { a.windUpAt = -1.0; a.moveSpeed = 0.0f; a.state = State::Flee; continue; }
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
            a.fleeFrom = sameBody ? tgLocal : a.pos;
        }
        // Hostile ones (not on peaceful): thornbacks charge within 12 m (or for 20 s after being hit),
        // dusk stalkers hunt within 70 m after dark. A campfire within 10 m keeps the player safe, and
        // stalkers will not come within 14 m of one.
        bool hunting = false;
        if (s.hostile && threatsOn && Difficulty_ThreatScale() > 0.0f && sameBody && a.fleeTimer <= 0.0f && tg->vulnerable) { // not inside a cockpit
            Vector3 fire;
            const bool playerSafe = NearFire(tgLocal, 10.0f, fire);
            if (s.night && NearFire(a.pos, 14.0f, fire)) { a.fleeTimer = 3.0f; a.fleeFrom = fire; }
            else if (s.night) hunting = dark && pdist < 70.0f && !playerSafe;
            else hunting = (pdist < 12.0f || now < a.angryUntil) && pdist < 30.0f && !playerSafe;
        }
        if (!hunting) a.windUpAt = -1.0; // a wind-up broken off (the prey gone or safe) never lands later
        if (hunting) {
            a.state = State::Hunt;
            toward = Vector3Subtract(tgLocal, a.pos);
            steer = true;
            speed = s.run * (s.night ? 0.8f : 1.0f);
            const float reach = 0.5f * s.size.x + 1.0f;
            // A telegraphed attack (#165): in reach it winds up (a growl, a warning) and strikes kBiteWindUp later; the
            // hit goes through Combat_Damage (a shield block, knockback away from it). Out of reach by then: it missed.
            if (pdist < reach) {
                speed = 0.0f;
                if (a.windUpAt < 0.0 && simClock - a.lastBite > 1.2) { // (on simClock: tests step animals many times in a frame)
                    a.windUpAt = simClock;
                    if (tg->session == 0) { Audio_Play(Sfx::Error, 0.5f, 0.35f);
                    GameUI_Toast(TextFormat("The %s %s!", s.name.c_str(), s.night ? "lunges" : "lowers its head")); } // (a remote player sees the wind-up in the animal's flags)
                }
            }
            if (a.windUpAt >= 0.0 && simClock - a.windUpAt >= kBiteWindUp) {
                a.windUpAt = -1.0;
                a.lastBite = simClock;
                if (pdist < reach + 0.4f) {
                    DamageEvent e;
                    e.amount = 8.0f * Difficulty_ThreatScale();
                    e.poise = e.amount * PoiseTune::kAttackPoisePerDamage; // (#1217) a strong bite or ram: 20
                    e.kind = s.night ? DamageKind::Bite : DamageKind::Melee;
                    e.hasSource = true;
                    e.source = body->toWorldD(DVec3(a.pos));
                    e.knockback = s.night ? 3.0f : 6.0f; // a ram shoves harder than a bite
                    e.sourceId = a.netId;
                    if (tg->session != 0) { // (#334) a remote player: its client applies the hit (shield, health), the server the shove
                        if (gAnimalDamageRemote) gAnimalDamageRemote(tg->session, e);
                    } else {
                        const DamageResult r = Combat_Damage(LocalSession(), e);
                        Audio_Play(Sfx::Error, 0.7f, r.blocked ? 0.9f : 0.5f);
                        GameUI_Toast(r.blocked ? TextFormat("You block the %s", s.name.c_str())
                                               : TextFormat("The %s %s you!", s.name.c_str(), s.night ? "bites" : "rams"));
                    }
                } else if (tg->session == 0) GameUI_Toast(TextFormat("The %s misses", s.name.c_str()));
            }
        } else if (a.fleeTimer > 0.0f) {
            a.fleeTimer -= step;
            a.state = State::Flee;
            if (sameBody) a.fleeFrom = tgLocal;
            toward = Vector3Subtract(a.pos, a.fleeFrom);
            steer = true;
            speed = s.run;
        } else if (sameBody && pdist > 2.5f && (now < a.followUntil || (pdist < 16.0f && !tg->held.empty() && std::find(s.diet.begin(), s.diet.end(), tg->held) != s.diet.end()))) {
            // Fed (30 s), or the player holds its food: it follows (wild ones too: the food wins).
            a.state = State::Follow;
            toward = Vector3Subtract(tgLocal, a.pos);
            steer = true;
            speed = std::min(s.run, s.walk * 2.0f);
        } else if (s.wild && sameBody && pdist < 25.0f) {
            a.fleeTimer = 4.0f;
            a.fleeFrom = tgLocal;
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
        // A planned route (#253): head for its next point while the goal stays put (within 6 m) and the route is fresh.
        // Fleeing, the goal is a point 25 m on away from the threat.
        a.planWait = std::max(0.0f, a.planWait - step);
        const Vector3 goalPt = !steer ? a.pos : a.state == State::Flee && Vector3Length(toward) > 1e-3f
            ? Vector3Add(a.pos, Vector3Scale(Vector3Normalize(toward), 25.0f)) : Vector3Add(a.pos, toward);
        if (!a.path.empty()) {
            a.pathAge += step;
            if (!steer || speed <= 0.0f || a.pathAge > 10.0f || Vector3Distance(a.pathGoal, goalPt) > 6.0f) a.path.clear();
            while (!a.path.empty() && Vector3Distance(a.path.front(), a.pos) < 0.7f) a.path.erase(a.path.begin());
            if (!a.path.empty()) {
                const Vector3 w = Vector3Subtract(a.path.front(), a.pos);
                a.heading = atan2f(Vector3DotProduct(w, north), Vector3DotProduct(w, east));
                a.detourTime = 0.0f;
            }
        }
        if (Combat1_On()) speed *= Status_SpeedFactor(a.status); // (#1215) chilled 0.7, frozen 0
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
        if (!ok && steer && !a.path.empty()) { a.path.clear(); a.planWait = 1.0f; } // the route's step failed (a thin wall between cells)
        // Blocked on the way somewhere: plan a route round (one plan a frame for all animals, then a 3 s wait after a
        // failed one); the local detour below covers the time between.
        // One search at a time, resumed across updates: its owner continues it (standing while it thinks, a few frames);
        // a search nobody has touched for 30 updates (its animal got going, or was removed) is dropped.
        if (gSearch.active && frameNo - gSearch.lastFrame > 30) gSearch.active = false;
        const bool mySearch = gSearch.active && gSearch.owner == (int)i && Vector3Distance(gSearch.from, a.pos) < 3.0f;
        if (!ok && steer && a.path.empty() && (mySearch || (a.planWait <= 0.0f && !gSearch.active)) && plansLeft > 0) {
            --plansLeft;
            if (!mySearch) PathBegin(*body, *field, a.pos, goalPt, height, std::max(reachUp, 0.35f), std::max(reachDown, 0.35f), (int)i, frameNo);
            std::vector<Vector3> route;
            if (!PathStep(route, frameNo)) { a.moveSpeed = 0.0f; continue; } // still searching: wait for it
            if (!route.empty()) {
                a.path = std::move(route);
                a.pathGoal = goalPt; a.pathAge = 0.0f; a.planWait = 1.0f;
                const Vector3 w = Vector3Subtract(a.path.front(), a.pos);
                a.heading = atan2f(Vector3DotProduct(w, north), Vector3DotProduct(w, east));
                a.detourTime = 0.0f;
                ok = tryStep(a.heading);
            } else a.planWait = 3.0f;
        }
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
    if (here && SimTime() - lastPenCheck > 2.0) {
        lastPenCheck = SimTime();
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
            young.netId = NetId_Allocate(NetKind::Animal);
            const bool seen = here && here->name == a.body; // before the push: it may move the vector (a, b)
            const std::string kind = species[a.species].name;
            animals.push_back(young);
            if (seen) GameUI_Toast(TextFormat("A young %s was born", kind.c_str()));
            QuestBook_Count("bred");
            break;
        }
    }
    lastUpdateMs = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
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
int gMeleeShocks = 0, gMeleeStaggers = 0;
uint64_t gLockNet = 0;              // the soft lock-on's creature (0: none)
Vector2 gLockScreen{ -1.0f, -1.0f }; // where it is on screen (the HUD marker)
WeaponClass gMeleeClass = WeaponClass::Tools;
int gMeleeHits = 0;                 // hits landed (tests)

int IndexOfNet(uint64_t id) { for (size_t i = 0; i < animals.size(); ++i) if (animals[i].netId == id) return (int)i; return -1; }

Vector3 AnimalForward(const Animal& a) {
    Vector3 east, north;
    TangentFrame(Vector3Normalize(a.pos), east, north);
    return Vector3Add(Vector3Scale(east, cosf(a.heading)), Vector3Scale(north, sinf(a.heading)));
}

// Its body as a capsule along its length (body-local), its radius half its height or width.
void AnimalCapsule(const Animal& a, Vector3& q0, Vector3& q1, float& r) {
    const Species& s = species[a.species];
    const float scale = 0.6f + 0.4f * a.age;
    const Vector3 up = Vector3Normalize(a.pos), fwd = AnimalForward(a);
    r = 0.5f * std::max(s.size.x, s.size.y) * scale;
    const Vector3 c = Vector3Add(a.pos, Vector3Scale(up, 0.5f * s.size.y * scale));
    const float hl = std::max(0.0f, 0.5f * s.size.z * scale - r);
    q0 = Vector3Subtract(c, Vector3Scale(fwd, hl)); q1 = Vector3Add(c, Vector3Scale(fwd, hl));
}

Vector3 AnimalCentre(const Animal& a) {
    const Species& s = species[a.species];
    return Vector3Add(a.pos, Vector3Scale(Vector3Normalize(a.pos), 0.5f * s.size.y * (0.6f + 0.4f * a.age)));
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
    if (gMelee.phase == MeleePhase::Idle) gMeleeClass = Weapon_ClassOf(terrainTool.heldTool);
    const WeaponSpec& w = Weapon_Spec(gMeleeClass);
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
            if (a.body != here->name || Vector3Distance(a.pos, eye) > w.reach + 4.0f) continue;
            Vector3 q0, q1, hit; float r;
            AnimalCapsule(a, q0, q1, r);
            if (Melee_SweepCapsule(eye, fwd, up, w.reach, half, -half, q0, q1, r, &hit) && !Shielded(world, *here, eye, hit)) return true;
        }
        return false;
    };
    MeleeInput mi;
    mi.pressed = live && in.hit(Act::Primary) && !gPrimaryLatch.latched &&
                 (gMelee.phase != MeleePhase::Idle || aimed >= 0 || lockIdx >= 0 || inArc());
    mi.held = live && in.down(Act::Primary);
    mi.block = survival.blocking;
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
        if (a.body != here->name || Vector3Distance(a.pos, eye) > w.reach + 4.0f || Melee_Struck(gMelee, a.netId)) continue;
        Vector3 q0, q1, at; float r;
        AnimalCapsule(a, q0, q1, r);
        if (!Melee_SweepCapsule(eye, fwd, up, w.reach, a0, a1, q0, q1, r, &at)) continue;
        if (Shielded(world, *here, eye, at)) continue; // terrain or a part takes the blow
        hits.push_back({ (int)i, at });
    }
    const ToolSpec t = ToolOf(terrainTool.heldTool);
    const float base = t.type == ToolType::Hands ? 2.0f : 3.0f + 3.0f * t.tier; // today's damage by tool (weapons: #1219)
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
        e.fromBehind = Vector3DotProduct(AnimalForward(a), Vector3Normalize(Vector3Subtract(h.at, eye))) > 0.5f; // its back to the blade
        survival.lastDealtAt = now;
        ++gMeleeHits;
        if (gAnimalsReplicated) { if (gAnimalHitRequest) gAnimalHitRequest(a.netId); } // (#334) a client: the server decides
        else Animals_HitEvent(h.i, e);
    }
}
} // namespace

bool Animals_MeleePose(int& phase, float& f, float& charge, bool& back) {
    if (!Combat2_On() || gMelee.phase == MeleePhase::Idle) return false;
    const WeaponSpec& w = Weapon_Spec(gMeleeClass);
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

void Animals_LockOnView(Scene& world, float dt) {
    const int i = gLockNet ? IndexOfNet(gLockNet) : -1;
    Planetoid* here = player.nearest;
    if (i < 0 || !here || player.flyMode) return;
    const Vector3 to = Vector3Subtract(here->toWorld(AnimalCentre(animals[i])), world.camera.position);
    player.forward = Melee_SoftTurn(player.forward, player.up, to, dt);
}

uint64_t Animals_LockOn() { return gLockNet; }
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
            const Vector3 c = Vector3Add(a.pos, Vector3Scale(Vector3Normalize(a.pos), s.size.y * 0.5f * scale));
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
                if (gAnimalsReplicated) { if (gAnimalHitRequest) gAnimalHitRequest(animals[aimed].netId); } // (#334) a client: the server decides
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
        // Between the last two simulation steps (#315): smooth at any frame rate (a jump over 2 m is a teleport).
        Vector3 at = a.pos;
        float heading = a.heading;
        if (a.prevValid && Vector3Distance(a.prevPos, a.pos) < 2.0f) {
            const float al = SimAlpha();
            at = Vector3Lerp(a.prevPos, a.pos, al);
            heading = a.prevHeading + remainderf(a.heading - a.prevHeading, 2.0f * PI) * al;
        }
        const Vector3 base = body->toWorld(at);
        if (Vector3Distance(base, world.camera.position) > 160.0f) continue;
        const Species& s = species[a.species];
        const Vector3 lup = Vector3Normalize(at);
        Vector3 east, north;
        TangentFrame(lup, east, north);
        const Vector3 f = body->dirToWorld(Vector3Add(Vector3Scale(east, cosf(heading)), Vector3Scale(north, sinf(heading))));
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
    const char* state = a.state == State::Hunt ? "hunting you" : a.state == State::Flee ? "fleeing" : a.state == State::Follow ? "follows you" : a.state == State::Rest ? "resting" : a.state == State::Walk ? "wandering" : "grazing";
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
    for (const Animal& a : animals) n += a.body == body.name && a.simFrame == frameNo && Vector3Distance(a.pos, local) < range ? 1 : 0;
    return n;
}
uint64_t Animals_NetId(int index) { return index >= 0 && index < (int)animals.size() ? animals[index].netId : 0; }

bool Animals_LoadLine(const char* tok, FILE* f) {
    LoadSpecies();
    char id[64] = {}, body[128] = {};
    if (strcmp(tok, "animal") == 0) {
        Animal a;
        int fed = 0;
        if (fscanf(f, " %63s %127s %f %f %f %f %f %f %d %lf %lf", id, body, &a.pos.x, &a.pos.y, &a.pos.z, &a.heading, &a.health, &a.age, &fed, &a.lastBred, &a.productReady) == 11) {
            a.species = SpeciesIndex(id);
            if (a.species >= 0) {
                a.body = body; a.fed = fed != 0; a.rng = (uint32_t)animals.size() * 7919u + 1u;
                a.netId = NetId_Allocate(NetKind::Animal); // a save from before #317; its "animal_id" line replaces it
                animals.push_back(a);
                lastLoaded = (int)animals.size() - 1;
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
        n.netId = a.netId; n.species = (uint8_t)a.species; n.body = it->second; n.pos = a.pos; n.heading = a.heading;
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
            a.netId = n.netId; a.species = n.species; a.body = b->name; a.pos = a.netPos = a.prevPos = n.pos; a.heading = a.netHeading = a.prevHeading = n.heading;
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
    const int i = Animals_FindNet(netId);
    if (i < 0) return false;
    Animal& a = animals[i];
    a.health -= damage;
    a.fleeTimer = !species[a.species].hostile || a.health < 0.35f * species[a.species].health ? 5.0f : 0.0f;
    if (species[a.species].hostile && a.fleeTimer <= 0.0f) a.angryUntil = WorldTime() + 20.0;
    a.fleeFrom = a.pos;
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

// (#1216 tests) A training dummy: the animal held at `pos` (body-local, on the ground) facing `face` (body-local), at rest.
void Animals_TestPin(int i, Vector3 pos, Vector3 face) {
    if (i < 0 || i >= (int)animals.size()) return;
    Animal& a = animals[i];
    Vector3 east, north;
    TangentFrame(Vector3Normalize(pos), east, north);
    a.pos = a.prevPos = pos;
    a.heading = a.prevHeading = atan2f(Vector3DotProduct(face, north), Vector3DotProduct(face, east));
    a.state = State::Rest; a.fleeTimer = 0.0f; a.path.clear();
}
