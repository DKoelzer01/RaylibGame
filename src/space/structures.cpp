#include "structures.h"
#include "magicdiscovery.h" // (#1250) First Sight pages only in chests near a ley line
#include "buildtool.h"
#include "factory.h"
#include "inventory.h"
#include "player.h"
#include "plots.h"
#include "scene.h"
#include "trees.h"
#include "animals.h"
#include "navigation.h"
#include "net.h"
#include "gameui.h"
#include "telescope.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

namespace fs = std::filesystem;

namespace {

struct LootEntry { std::string item; int lo = 1, hi = 1; float chance = 1.0f; };
std::map<std::string, std::vector<LootEntry>> lootTables;

struct Loot { Int3 cell; std::string table; };
struct Set {
    std::string name, kind;
    std::vector<int> biomes;       // empty: any
    float weight = 1.0f, bury = 0.0f, decay = 0.0f;
    std::vector<Part> parts;       // structure cells (x east, y up, z north before turning), lowest corner 0
    std::vector<Loot> loots;       // chests (also in parts)
    Int3 size{ 1, 1, 1 };          // cells
    std::vector<std::pair<Vector3, float>> carves; // wonders: spheres cut from the ground (metres from the site centre at ground level; radius)
    int yOffset = 0;               // the lowest part cell before normalising (parts may sit underground: wonders)
    struct Guard { std::string species; int count = 1; float radius = 6.0f; };
    std::vector<Guard> guards;     // animals set down round the site when it is built ("guard <species> <n> [radius m]")
    float dungeon = 0.0f;          // the odds of a dungeon under it (#719: "dungeon <chance>"; RAYLIBGAME_DUNGEONS=0 off)
};
std::vector<Set> sets;
bool loaded = false;

// Sites: one per 512 m tile at most (#164). Near the start (design 2026-11-17): none within kSpawnKeepOut, and two
// guaranteed sites - a camp and a ruin - between kRingNear and kRingFar (the "ring" sites, tiles kRingTile + 0/1).
constexpr float kTile = 512.0f, kBuildRange = 250.0f, kLootRange = 60.0f, kSiteOdds = 0.35f;
constexpr float kSpawnKeepOut = 150.0f, kRingNear = 300.0f, kRingFar = 800.0f, kFoundRange = 60.0f;
constexpr float kWaterMargin = 8.0f, kMaxRise = 2.5f; // m: clear of water this far round the footprint; rise over it
constexpr uint64_t kRingTile = 7ull << 48;             // past the six cube faces
constexpr int kSpotTries = 4;                          // spots tried in a tile before it is left empty
std::unordered_map<std::string, std::unordered_set<uint64_t>> builtTiles;   // body name -> tiles decided (built or skipped)
// Built sites (#164 discovery): where they are and whether the player has found them (a waypoint then).
struct BuiltSite { std::string body; uint64_t tile = 0; Vector3 dir{}; std::string name, kind; bool found = false; };
std::vector<BuiltSite> builtSites;
std::unordered_map<std::string, std::unordered_set<uint64_t>> foundTiles; // body name -> tiles whose site was found (saved)
std::set<std::string> foundKinds; // kinds of site found (saved): the Help tab lists only these with "Discovered only"
struct PendingLoot { std::string body; uint32_t plot; int part; std::string table; uint32_t seed; };
std::vector<PendingLoot> pending;
std::vector<std::pair<uint32_t, int>> testFilled; // tests: the chests Structures_Test built (plot, part)
double lastCheck = -10.0;
const char* lastSkip = "";  // why the last site was not built (tests)
int sitesBuilt = 0;

uint32_t Hash(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33;
    return (uint32_t)x;
}
float Unit(uint32_t& s) { s = Hash(s + 0x9e3779b97f4a7c15ull); return (s >> 8) * (1.0f / 16777216.0f); }

PartType TypeNamed(const std::string& n) {
    static const std::pair<const char*, PartType> kNames[] = {
        { "cube", PartType::Cube }, { "block", PartType::Block }, { "slope", PartType::Slope }, { "corner", PartType::Corner },
        { "slab", PartType::Slab }, { "panel", PartType::Panel }, { "pane", PartType::Pane }, { "post", PartType::Post },
        { "chest", PartType::Chest }, { "fence", PartType::Fence }, { "gate", PartType::Gate }, { "campfire", PartType::Campfire },
        { "fuel_tank", PartType::FuelTank }, { "rocket_engine", PartType::RocketEngine }, { "cabin", PartType::Cabin },
        { "landing_leg", PartType::LandingLeg }, { "torch", PartType::Torch },
        { "farmland", PartType::Farmland }, { "generator", PartType::Generator }, { "battery", PartType::Battery }, // (#720: salvage)
    };
    for (auto& k : kNames) if (n == k.first) return k.second;
    return PartType::Count;
}
int CellsOf(PartType t) { return t == PartType::Block ? 2 : 1; } // the step for fill/walls

void Footprint(const Part& p, std::vector<std::pair<Int3, CellShape>>& fp) { Plot::Footprint(p, fp); }

void LoadLoot() {
    FILE* f = fopen("resources/data/loot.txt", "r");
    if (!f) return;
    char line[1024];
    while (fgets(line, sizeof line, f)) {
        line[strcspn(line, "\r\n")] = 0; // CRLF checkouts: a trailing CR would spoil the last word
        if (line[0] == '#' || line[0] == 0) continue;
        std::stringstream ss(line);
        std::string cell, table;
        bool first = true;
        while (std::getline(ss, cell, '|')) {
            std::stringstream cs(cell);
            if (first) { cs >> table; first = false; continue; }
            LootEntry e;
            std::string range;
            if (!(cs >> e.item >> range)) continue;
            if (sscanf(range.c_str(), "%d-%d", &e.lo, &e.hi) < 1) continue;
            e.hi = std::max(e.hi, e.lo);
            cs >> e.chance;
            if (itemDB.find(e.item)) lootTables[table].push_back(e); // unknown items are left out
        }
    }
    fclose(f);
}

bool LoadSet(const std::string& path, Set& s) {
    FILE* f = fopen(path.c_str(), "r");
    if (!f) return false;
    s.name = fs::path(path).stem().string();
    char line[512];
    auto add = [&](PartType t, Int3 o, int orient, int mat) { s.parts.push_back({ t, o, (uint8_t)orient, (uint8_t)std::max(0, mat) }); };
    auto matOf = [](const std::string& m) { return IsShape(PartType::Cube) ? BuildMaterialForItem(m.c_str()) : 0; };
    while (fgets(line, sizeof line, f)) {
        line[strcspn(line, "\r\n")] = 0; // CRLF checkouts: a trailing CR would spoil the last word ("stone\r" is no material)
        if (line[0] == '#') continue;
        std::stringstream ss(line);
        std::string cmd;
        if (!(ss >> cmd)) continue;
        if (cmd == "name") { std::getline(ss >> std::ws, s.name); while (!s.name.empty() && (s.name.back() == '\r' || s.name.back() == '\n')) s.name.pop_back(); }
        else if (cmd == "kind") ss >> s.kind;
        else if (cmd == "weight") ss >> s.weight;
        else if (cmd == "bury") ss >> s.bury;
        else if (cmd == "decay") ss >> s.decay;
        else if (cmd == "biomes") {
            std::string b;
            static const char* kBiomes[] = { "Ocean", "Jungle", "Savanna", "Desert", "Badlands", "Forest", "Meadow", "Swamp", "Taiga", "Tundra", "Alpine", "Polar" };
            while (ss >> b) { if (b == "any") { s.biomes.clear(); break; } for (int i = 0; i < 12; ++i) if (b == kBiomes[i]) s.biomes.push_back(i); }
        } else if (cmd == "part") {
            std::string t, m; int x, y, z, o;
            if (ss >> t >> x >> y >> z >> o >> m && TypeNamed(t) != PartType::Count)
                add(TypeNamed(t), { x, y, z }, o, IsShape(TypeNamed(t)) ? matOf(m) : 0);
        } else if (cmd == "fill" || cmd == "walls") {
            std::string t, m; int x0, y0, z0, x1, y1, z1;
            if (!(ss >> t >> x0 >> y0 >> z0 >> x1 >> y1 >> z1 >> m) || TypeNamed(t) == PartType::Count) continue;
            const PartType pt = TypeNamed(t);
            const int st = CellsOf(pt);
            for (int y = y0; y + st - 1 <= y1; y += st)
                for (int z = z0; z + st - 1 <= z1; z += st)
                    for (int x = x0; x + st - 1 <= x1; x += st) {
                        if (cmd == "walls" && !(x == x0 || x + st - 1 == x1 || z == z0 || z + st - 1 == z1)) continue;
                        add(pt, { x, y, z }, 0, IsShape(pt) ? matOf(m) : 0);
                    }
        } else if (cmd == "clear") {
            int x0, y0, z0, x1, y1, z1;
            if (!(ss >> x0 >> y0 >> z0 >> x1 >> y1 >> z1)) continue;
            std::vector<std::pair<Int3, CellShape>> fp;
            s.parts.erase(std::remove_if(s.parts.begin(), s.parts.end(), [&](const Part& p) {
                Footprint(p, fp);
                for (auto& c : fp) if (c.first.x >= x0 && c.first.x <= x1 && c.first.y >= y0 && c.first.y <= y1 && c.first.z >= z0 && c.first.z <= z1) return true;
                return false;
            }), s.parts.end());
        } else if (cmd == "dungeon") {
            ss >> s.dungeon;
        } else if (cmd == "guard") {
            Set::Guard g; if (ss >> g.species >> g.count) { float r; if (ss >> r) g.radius = r; s.guards.push_back(g); }
        } else if (cmd == "carve") {
            Vector3 c; float r;
            if (ss >> c.x >> c.y >> c.z >> r) s.carves.push_back({ c, r });
        } else if (cmd == "loot") {
            Loot L; if (!(ss >> L.cell.x >> L.cell.y >> L.cell.z >> L.table)) continue;
            s.loots.push_back(L);
            add(PartType::Chest, L.cell, 0, 0);
        }
    }
    fclose(f);
    if (s.parts.empty()) return false;
    // Size and lowest corner 0 (footprints, so turned parts count).
    std::vector<std::pair<Int3, CellShape>> fp;
    Int3 lo{ 1 << 20, 1 << 20, 1 << 20 }, hi{ -(1 << 20), -(1 << 20), -(1 << 20) };
    for (const Part& p : s.parts) { Footprint(p, fp); for (auto& c : fp) { lo = { std::min(lo.x, c.first.x), std::min(lo.y, c.first.y), std::min(lo.z, c.first.z) }; hi = { std::max(hi.x, c.first.x), std::max(hi.y, c.first.y), std::max(hi.z, c.first.z) }; } }
    for (Part& p : s.parts) p.origin = { p.origin.x - lo.x, p.origin.y - lo.y, p.origin.z - lo.z };
    for (Loot& L : s.loots) L.cell = { L.cell.x - lo.x, L.cell.y - lo.y, L.cell.z - lo.z };
    s.size = { hi.x - lo.x + 1, hi.y - lo.y + 1, hi.z - lo.z + 1 };
    s.yOffset = lo.y;
    return true;
}

void Ensure() {
    if (loaded) return;
    loaded = true;
    LoadLoot();
    std::error_code ec;
    std::vector<std::string> files;
    if (fs::exists("resources/structures", ec))
        for (const auto& e : fs::directory_iterator("resources/structures", ec))
            if (e.is_regular_file() && e.path().extension() == ".st") files.push_back(e.path().string());
    std::sort(files.begin(), files.end());
    for (const std::string& f : files) { Set s; if (LoadSet(f, s)) sets.push_back(std::move(s)); }
}

// Turning a part a quarter turn about up (plot y): the origin cell turns as a vector (x, y, z) -> (z, y, -x) and the
// orientation picks up the same turn (RotateOrientation pre-multiplies), so every footprint cell turns rigidly.
Part Turned(Part p, int quarters) {
    for (int k = 0; k < (quarters & 3); ++k) {
        p.origin = { p.origin.z, p.origin.y, -p.origin.x };
        p.orient = (uint8_t)RotateOrientation(p.orient, 1);
    }
    return p;
}
Int3 TurnedCell(Int3 c, int quarters) { for (int k = 0; k < (quarters & 3); ++k) c = { c.z, c.y, -c.x }; return c; }

// The biome's stone (#164: biome materials): sandstone in the dry south, granite in the cold, limestone in the wet.
int BiomeStone(int mat, Biome b) {
    static const int stone = BuildMaterialForItem("stone_brick");
    if (mat != stone) return mat;
    const char* m = nullptr;
    switch (b) {
        case Biome::Desert: case Biome::Badlands: case Biome::Savanna: m = "sandstone_brick"; break;
        case Biome::Tundra: case Biome::Alpine: case Biome::Taiga: case Biome::Polar: m = "granite_brick"; break;
        case Biome::Jungle: case Biome::Swamp: m = "limestone_brick"; break;
        default: break;
    }
    const int r = m ? BuildMaterialForItem(m) : -1;
    return r >= 0 ? r : mat;
}

// Tiles: the cube-sphere face of the direction, then a grid of ~kTile metres on it.
uint64_t TileOf(Vector3 dir, float radius, Vector3* centre) {
    const Vector3 a = { fabsf(dir.x), fabsf(dir.y), fabsf(dir.z) };
    int face; float u, v, m;
    if (a.x >= a.y && a.x >= a.z) { face = dir.x > 0 ? 0 : 1; m = a.x; u = dir.y / m; v = dir.z / m; }
    else if (a.y >= a.z) { face = dir.y > 0 ? 2 : 3; m = a.y; u = dir.x / m; v = dir.z / m; }
    else { face = dir.z > 0 ? 4 : 5; m = a.z; u = dir.x / m; v = dir.y / m; }
    const int n = std::max(1, (int)ceilf(0.5f * PI * radius / kTile));
    const int i = std::clamp((int)floorf((u + 1.0f) * 0.5f * n), 0, n - 1), j = std::clamp((int)floorf((v + 1.0f) * 0.5f * n), 0, n - 1);
    if (centre) {
        const float cu = (i + 0.5f) / n * 2.0f - 1.0f, cv = (j + 0.5f) / n * 2.0f - 1.0f;
        Vector3 c;
        switch (face) {
            case 0: c = { 1, cu, cv }; break; case 1: c = { -1, cu, cv }; break;
            case 2: c = { cu, 1, cv }; break; case 3: c = { cu, -1, cv }; break;
            case 4: c = { cu, cv, 1 }; break; default: c = { cu, cv, -1 }; break;
        }
        *centre = Vector3Normalize(c);
    }
    return ((uint64_t)face << 48) | ((uint64_t)i << 24) | (uint64_t)j;
}

struct Site { bool any = false; int set = -1; Vector3 dir{}; int quarters = 0; uint32_t seed = 0; };

// Metres along the ground between two directions on a body of this radius.
float GroundDistance(Vector3 a, Vector3 b, float radius) { return acosf(std::clamp(Vector3DotProduct(a, b), -1.0f, 1.0f)) * radius; }
// Tangent axes at a direction (the same everywhere this file needs them).
void TangentAxes(Vector3 dir, Vector3& t1, Vector3& t2) {
    t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    t2 = Vector3CrossProduct(t1, dir);
}
bool Wet(const PlanetField& pf, Vector3 d) { // the sea, or a river or lake at the surface there
    const float e = pf.elevation(d);
    if (e < 1.0f) return true;
    if (!pf.inlandWaterNear(d, 0.0f)) return false;
    const float w = pf.waterElevation(d);
    return w != PlanetField::kNoWater && w > e - 0.5f;
}
// Whether a set's footprint suits a spot (#164 design): dry at the centre, the corners and the edge midpoints of the
// footprint grown by kWaterMargin, and at most kMaxRise between the lowest and highest of the footprint's own points.
bool SpotOk(const PlanetField& pf, Vector3 dir, const Set& S, const char** why) {
    Vector3 t1, t2; TangentAxes(dir, t1, t2);
    const float half = 0.5f * kCellSize * std::max(S.size.x, S.size.z);
    float lo = 1e9f, hi = -1e9f;
    for (int k = 0; k < 9; ++k) {
        const float u = (float)(k % 3 - 1), v = (float)(k / 3 - 1);
        auto at = [&](float r) { return Vector3Normalize(Vector3Add(dir, Vector3Add(Vector3Scale(t1, u * r / pf.radius), Vector3Scale(t2, v * r / pf.radius)))); };
        if (Wet(pf, at(half + kWaterMargin))) { if (why) *why = "water nearby"; return false; }
        const float e = pf.elevation(at(half));
        lo = std::min(lo, e); hi = std::max(hi, e);
    }
    if (hi - lo > kMaxRise) { if (why) *why = "too steep"; return false; }
    return true;
}
// The body's start (where a new game spawns): the keep-out and the ring sites are measured from it.
Vector3 StartOf(const PlanetField& pf, unsigned int bodySeed) {
    static std::unordered_map<uint64_t, Vector3> cache;
    const uint64_t key = (uint64_t)bodySeed * 0x9e3779b97f4a7c15ull ^ (uint64_t)pf.foliageSeed;
    auto it = cache.find(key);
    if (it != cache.end()) return it->second;
    return cache[key] = pf.findSpawnDir(bodySeed);
}
bool Allows(const Set& st, int biome, const char* kind, bool unguarded = false) {
    return (!kind || st.kind == kind) && (!unguarded || st.guards.empty()) && (st.biomes.empty() || std::find(st.biomes.begin(), st.biomes.end(), biome) != st.biomes.end());
}
// A set the biome allows (of `kind` if given), by weight; -1 none.
int PickSet(int biome, const char* kind, uint32_t& h, bool unguarded = false) {
    float total = 0.0f;
    for (const Set& st : sets) if (Allows(st, biome, kind, unguarded)) total += st.weight;
    if (total <= 0.0f) return -1;
    float pick = Unit(h) * total;
    for (int i = 0; i < (int)sets.size(); ++i)
        if (Allows(sets[i], biome, kind, unguarded) && (pick -= sets[i].weight) <= 0.0f) return i;
    return -1;
}
// The two ring sites near the start (#164): a camp (k 0) and a ruin (k 1) somewhere between kRingNear and kRingFar on
// suitable ground, the bearing and distance drawn from the seed (tries until one suits; none: no ring site).
Site RingSite(const PlanetField& pf, Vector3 start, int k) {
    // Cached per world and start: every tile's check asks for both ring sites.
    static std::unordered_map<uint64_t, Site> cache;
    const uint64_t key = ((uint64_t)pf.foliageSeed * 0x9e3779b97f4a7c15ull) ^ ((uint64_t)(uint32_t)(start.x * 1e6f) * 31 + (uint64_t)(uint32_t)(start.z * 1e6f)) ^ (uint64_t)k;
    if (auto it = cache.find(key); it != cache.end()) return it->second;
    Site s;
    static const char* kKinds[2] = { "camp", "ruin" };
    uint32_t h = Hash(((uint64_t)pf.foliageSeed << 8) ^ 0x51e5u ^ (uint64_t)(k + 1) * 0x9e3779b97f4a7c15ull);
    Vector3 t1, t2; TangentAxes(start, t1, t2);
    for (int tries = 0; tries < 24; ++tries) {
        const float ang = Unit(h) * 2.0f * PI, dist = kRingNear + Unit(h) * (kRingFar - kRingNear);
        const Vector3 d = Vector3Normalize(Vector3Add(start, Vector3Add(Vector3Scale(t1, cosf(ang) * dist / pf.radius), Vector3Scale(t2, sinf(ang) * dist / pf.radius))));
        if (Wet(pf, d)) continue;
        const int set = PickSet((int)pf.biomeAt(d, pf.elevation(d)), kKinds[k], h, true); // no guards this near the start
        if (set < 0 || !SpotOk(pf, d, sets[set], nullptr)) continue;
        if (k == 1) { const Site camp = RingSite(pf, start, 0); if (camp.any && GroundDistance(camp.dir, d, pf.radius) < kSpawnKeepOut) continue; }
        s.any = true; s.set = set; s.dir = d; s.quarters = (int)(Unit(h) * 4.0f) & 3; s.seed = h;
        return cache[key] = s;
    }
    return cache[key] = s;
}
// The tile's site (if it has one): a set the biome there allows, chosen by weight, on a spot in the tile's middle 60%
// that suits its footprint (up to kSpotTries spots). Nothing within kSpawnKeepOut of the start or of a ring site.
// Ring tiles (kRingTile + k) are the ring sites.
Site SiteOf(const PlanetField& pf, uint64_t tile, Vector3 centreDir, float radius, Vector3 start) {
    if (tile >= kRingTile) return RingSite(pf, start, (int)(tile - kRingTile));
    Site s;
    uint32_t h = Hash(tile * 0x9e3779b97f4a7c15ull ^ (uint64_t)pf.foliageSeed * 0x632be59bd9b4e019ull ^ 0x5717c7u);
    if (Unit(h) > kSiteOdds || sets.empty()) return s;
    Vector3 t1, t2; TangentAxes(centreDir, t1, t2);
    for (int tries = 0; tries < kSpotTries; ++tries) {
        const float a = (Unit(h) - 0.5f) * 0.6f * kTile, b = (Unit(h) - 0.5f) * 0.6f * kTile;
        const Vector3 d = Vector3Normalize(Vector3Add(centreDir, Vector3Add(Vector3Scale(t1, a / radius), Vector3Scale(t2, b / radius))));
        if (Wet(pf, d) || GroundDistance(d, start, radius) < kSpawnKeepOut) continue;
        bool nearRing = false;
        for (int k = 0; k < 2 && !nearRing; ++k) { const Site r = RingSite(pf, start, k); nearRing = r.any && GroundDistance(r.dir, d, radius) < kSpawnKeepOut; }
        if (nearRing) continue;
        const int set = PickSet((int)pf.biomeAt(d, pf.elevation(d)), nullptr, h);
        if (set < 0 || !SpotOk(pf, d, sets[set], nullptr)) continue;
        s.any = true; s.set = set; s.dir = d; s.quarters = (int)(Unit(h) * 4.0f) & 3; s.seed = h;
        return s;
    }
    return s;
}

// A site's parts still to place (#164: a few a frame, so building a site is no hitch): the plot is made and the ground
// shaped when the site is decided; RunJob places parts and their foundations until done or out of time, and queues
// each chest's loot as it goes down.
struct BuildJob {
    uint64_t body = 0; uint32_t plot = 0; uint32_t seed = 0;
    std::vector<std::pair<Part, std::string>> parts; // the part, and its loot table for a chest ("" none)
    size_t next = 0;
    std::vector<TerrainEdit> edits;                  // carved first (a dungeon's rooms and corridors, #719)
    size_t nextEdit = 0;
};
std::vector<BuildJob> jobs;
int lastPlaced = 0; // parts the last RunJob placed (all of them, for a forced build)
bool RunJob(Planetoid& body, BuildJob& job, double budgetMs) { // true when the job is done
    Plot* plot = body.plots ? body.plots->byId(job.plot) : nullptr;
    if (!plot) return true; // the plot is gone (a reload, a world change): nothing to finish
    const auto t0 = std::chrono::steady_clock::now();
    auto spent = [&]() { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count(); };
    lastPlaced = 0;
    while (job.nextEdit < job.edits.size()) {
        ApplyTerrainEdit(body, job.edits[job.nextEdit++]);
        if (spent() > budgetMs) return false;
    }
    while (job.next < job.parts.size()) {
        const auto& [p, table] = job.parts[job.next++];
        const int idx = plot->place(p);
        if (idx >= 0) {
            BuildTool_LayFoundation(body, *plot, p);
            ++lastPlaced;
            if (!table.empty()) pending.push_back({ body.name, plot->id, idx, table, Hash(job.seed + idx) });
        }
        if (spent() > budgetMs) break;
    }
    return job.next >= job.parts.size();
}
// Each frame: the oldest job gets a little time (kJobBudgetMs).
constexpr double kJobBudgetMs = 1.5;
void PlaceJobs() {
    while (!jobs.empty()) {
        Planetoid* b = Planetoid::find(jobs.front().body);
        if (b && !RunJob(*b, jobs.front(), kJobBudgetMs)) return; // out of time: more next frame
        jobs.erase(jobs.begin());
        return;
    }
}

// Dungeons (#719, per the #164 design): under some ruins, a ramp beside the ruin goes down 8 m to the first of 2-3 floors
// 10 m apart (the deepest ~29 m down), joined by stairs. Each floor is 2-4 rooms joined in a tree by corridors that leave
// the rooms square through doorways with a door (two stacked gates). Rooms come in four shapes: plain, a pillared hall, a
// vaulted room (a high middle) and a crypt (alcoves in the walls). Each room has a torch; the last room of a floor a chest
// (the deepest a dungeon_rare one, with the dungeon-only relic now and then), the others now and then a common chest. Traps
// (slice 2): pressure plates in some corridors (a falling block or darts) and gas vents in some rooms (dungeontraps.cpp).
// Everything is terrain carves and plot parts on the site's plot (saved as edits and parts, nothing new).
// RAYLIBGAME_DUNGEONS=0 turns them off.
bool DungeonsOn() { static const bool on = !(getenv("RAYLIBGAME_DUNGEONS") && atoi(getenv("RAYLIBGAME_DUNGEONS")) == 0); return on; }
struct DungeonPlan {
    int rooms = 0, floors = 0, doors = 0, plates = 0, darts = 0, vents = 0, shapes = 0; // shapes: a bit per room shape used
    Vector3 entrance{}; float floorY = 0.0f, deepest = 0.0f;
    std::vector<Vector3> roomCentres;               // 1 m over each room's floor
    std::vector<std::pair<int, int>> links;
    std::vector<Vector3> corridorMids, stairPoints, rampPoints; // 1 m over the floor: must be open
};
DungeonPlan lastDungeon; // tests: the last dungeon planned (plot metres)
bool forceDungeon = false; // tests: every ruin with a dungeon chance gets one
bool noDungeons = false;   // tests: none (--test-structures checks the surface sets and their chests)
struct DungeonArgs { float floorY = 0.0f, half = 0.0f; int stone = 0; };
DungeonArgs lastDungeonArgs; // tests: what the last site would plan its dungeon with (built with noDungeons, then planned by hand)
void PlanDungeon(BuildJob& job, Planetoid& body, const PlanetField& pf, const Plot& plot, float floorY, float ruinHalfX, uint32_t seed, int stoneMat) {
    // Rock over every room and corridor (at least kCover): where the ground falls away from the ruin a room came up
    // under a thin crust and trees showed through its ceiling. A floor that fits fewer than 2 rooms ends the dungeon there,
    // and a dungeon whose first floor fits none is dropped.
    constexpr float kCover = 3.0f, kDepth = 8.0f, kFloorGap = 10.0f, kRoomH = 3.0f, kCorrHalfW = 0.75f, kCorrH = 2.4f;
    constexpr float kStepDown = 0.35f, kStepRun = 0.5f, kStub = 1.2f; // stairs at 35 deg; a doorway's passage out of the wall
    auto cover = [&](Vector3 plotTop) { const Vector3 t = plot.toBody(plotTop); return pf.radius + pf.elevation(Vector3Normalize(t)) - Vector3Length(t); };
    auto snap = [](float y) { return floorf(y / kCellSize) * kCellSize; };
    auto cellOf = [](float x, float y, float z) { return Int3{ (int)floorf(x / kCellSize), (int)floorf(y / kCellSize), (int)floorf(z / kCellSize) }; };
    std::vector<TerrainEdit> cuts;
    std::vector<std::pair<Part, std::string>> parts;
    struct Room { Vector3 c; float hx, hz; int shape; int floor; uint8_t walls = 0; }; // walls: a bit per wall with a door or stairs (+x, -x, +z, -z)
    static const Vector2 kTemplates[] = { { 3.0f, 3.0f }, { 5.0f, 2.5f }, { 2.5f, 2.5f }, { 4.0f, 4.0f }, { 2.5f, 4.5f }, { 4.5f, 3.5f } }; // half sizes (m)
    constexpr int kTemplateCount = 6;
    uint32_t h = Hash((uint64_t)seed * 0x9e3779b97f4a7c15ull ^ 0xd06e0ull);
    auto box = [&](Vector3 c, Vector3 half) {
        TerrainEdit e{};
        e.carve = true; e.box = true;
        e.center = plot.toBody(c);
        e.half = half; e.ax = plot.east; e.ay = plot.up; e.az = plot.north;
        e.radius = Vector3Length(half);
        cuts.push_back(e);
    };
    auto covered = [&](const Room& r, float flr) { // the ceiling's middle and corners under kCover of rock
        const float top = flr + kRoomH + 1.5f;     // (a vaulted room's middle is higher)
        for (int k = 0; k < 5; ++k) {
            const float u = k == 0 ? 0.0f : (k & 1 ? 1.0f : -1.0f), v = k == 0 ? 0.0f : (k & 2 ? 1.0f : -1.0f);
            if (cover({ r.c.x + u * (r.hx + 1.0f), top, r.c.z + v * (r.hz + 1.0f) }) < kCover) return false; // +1 m: alcoves
        }
        return true;
    };
    // A part whose footprint's lowest cell lands on `lo` (gates turn about their origin cell, so place by the footprint).
    auto placeAt = [&](PartType t, uint8_t orient, Int3 lo, uint8_t material, const std::string& table) {
        Part p{ t, { 0, 0, 0 }, orient, material };
        std::vector<std::pair<Int3, CellShape>> fp;
        Plot::Footprint(p, fp);
        Int3 mn = fp.empty() ? Int3{ 0, 0, 0 } : fp[0].first;
        for (const auto& [c, s] : fp) { mn.x = std::min(mn.x, c.x); mn.y = std::min(mn.y, c.y); mn.z = std::min(mn.z, c.z); }
        p.origin = { lo.x - mn.x, lo.y - mn.y, lo.z - mn.z };
        parts.push_back({ p, table });
    };
    // A doorway in room r's wall `wall` (0 +x, 1 -x, 2 +z, 3 -z) at lateral `at`: a 1 x 2 m passage kStub out of the wall,
    // a door in it. Returns the point kStub out, where the corridor starts.
    auto doorway = [&](Room& r, int wall, float at, float flr) {
        r.walls |= (uint8_t)(1 << wall);
        const bool xWall = wall < 2;
        const float sign = (wall & 1) ? -1.0f : 1.0f;
        const float wallPos = xWall ? r.c.x + sign * r.hx : r.c.z + sign * r.hz;
        const float lim = (xWall ? r.hz : r.hx) - 0.8f;
        const float lat0 = floorf((std::clamp(at, (xWall ? r.c.z : r.c.x) - lim, (xWall ? r.c.z : r.c.x) + lim) - 0.5f) / kCellSize + 0.5f) * kCellSize; // the door's low edge, on a cell
        const float mid = wallPos + sign * (kStub - 0.3f) * 0.5f, depth = (kStub + 0.3f) * 0.5f;
        if (xWall) box({ mid, flr + 1.0f, lat0 + 0.5f }, { depth, 1.0f, 0.5f });
        else box({ lat0 + 0.5f, flr + 1.0f, mid }, { 0.5f, 1.0f, depth });
        // The door: two gates, one over the other, in the wall's plane.
        const uint8_t orient = xWall ? 0 : (uint8_t)RotateOrientation(0, 1);
        const Int3 lo = xWall ? cellOf(wallPos, flr + 0.1f, lat0 + 0.01f) : cellOf(lat0 + 0.01f, flr + 0.1f, wallPos);
        for (int k = 0; k < 2; ++k) placeAt(PartType::Gate, orient, { lo.x, lo.y + 2 * k, lo.z }, 0, "");
        ++lastDungeon.doors;
        const float out = wallPos + sign * kStub;
        return xWall ? Vector3{ out, flr, lat0 + 0.5f } : Vector3{ lat0 + 0.5f, flr, out };
    };
    // The ramp: from just outside the ruin (+x), 0.5 m forward and 0.35 m down a step (35 deg), 2.4 m of headroom.
    lastDungeon = DungeonPlan{};
    const float x0 = ruinHalfX + 2.0f;
    const int steps = (int)ceilf(kDepth / kStepDown);
    for (int i = 0; i <= steps; ++i) {
        const float y = floorY - kStepDown * i;
        box({ x0 + kStepRun * i, y + 1.2f, 0.0f }, { 0.4f, 1.2f, kCorrHalfW });
        if (i > 1 && i < steps - 1 && i % 4 == 0) lastDungeon.rampPoints.push_back({ x0 + kStepRun * i, y + 1.0f, 0.0f });
    }
    float flr = snap(floorY - kStepDown * steps);
    Vector3 foot = { x0 + kStepRun * steps, flr, 0.0f }; // where this floor's first room starts
    Vector3 footDir = { 1, 0, 0 };
    std::vector<Room> rooms;
    std::vector<std::pair<int, int>> links;
    const int floors = 2 + (Unit(h) < 0.4f ? 1 : 0);
    int floorsBuilt = 0;
    for (int fl = 0; fl < floors; ++fl) {
        const size_t first = rooms.size();
        // The first room at the foot of the ramp or stairs, along their direction.
        {
            const int ti = (int)(Unit(h) * kTemplateCount) % kTemplateCount;
            const Vector2 t = kTemplates[ti];
            const float reach = fabsf(footDir.x) > 0.5f ? t.x : t.y;
            Room r{ Vector3Add(foot, Vector3Scale(footDir, reach)), t.x, t.y, (int)(Unit(h) * 4.0f) & 3, fl };
            r.c.y = flr;
            if (!covered(r, flr)) break; // the ground falls away here: no (deeper) floor
            rooms.push_back(r);
        }
        const int n = 2 + (int)(Unit(h) * 3.0f) % 3;
        for (int tries = 0; (int)(rooms.size() - first) < n && tries < 60; ++tries) {
            const size_t fi = first + (size_t)(Unit(h) * (rooms.size() - first)) % (rooms.size() - first);
            const Vector2 t = kTemplates[(int)(Unit(h) * kTemplateCount) % kTemplateCount];
            const int wall = (int)(Unit(h) * 4.0f) & 3;
            if (rooms[fi].walls & (1 << wall)) continue;
            const Room& from = rooms[fi];
            // Off one of its walls: square out, shifted sideways a little (the corridor turns once).
            const float gap = 4.0f + Unit(h) * 5.0f, side = (Unit(h) - 0.5f) * 6.0f, sign = (wall & 1) ? -1.0f : 1.0f;
            Vector3 c = from.c;
            if (wall < 2) { c.x += sign * (from.hx + gap + t.x); c.z += side; } else { c.z += sign * (from.hz + gap + t.y); c.x += side; }
            Room r{ c, t.x, t.y, (int)(Unit(h) * 4.0f) & 3, fl };
            bool clear = fl > 0 || fabsf(r.c.x) > ruinHalfX + 4.0f + r.hx || fabsf(r.c.z) > ruinHalfX + 4.0f + r.hz; // not under the ruin (or its ramp's top)
            for (size_t k = first; k < rooms.size(); ++k) { const Room& o = rooms[k]; clear = clear && (fabsf(o.c.x - r.c.x) > o.hx + r.hx + 3.0f || fabsf(o.c.z - r.c.z) > o.hz + r.hz + 3.0f); }
            if (fl == 0) clear = clear && (fabsf(r.c.z) > r.hz + 2.0f || r.c.x - r.hx > foot.x + 1.0f || r.c.x + r.hx < x0 - 1.0f); // off the ramp
            if (!clear || !covered(r, flr) || cover({ 0.5f * (from.c.x + r.c.x), flr + kCorrH, 0.5f * (from.c.z + r.c.z) }) < kCover) continue;
            // The corridor: out of a doorway in `from`, into one in `r` (the wall facing back), boxes every 0.6 m between.
            const int back = wall ^ 1;
            const size_t ri = rooms.size();
            rooms.push_back(r);
            const Vector3 a = doorway(rooms[fi], wall, wall < 2 ? r.c.z : r.c.x, flr);
            const Vector3 b = doorway(rooms[ri], back, wall < 2 ? rooms[fi].c.z : rooms[fi].c.x, flr);
            const float len = Vector2Length({ b.x - a.x, b.z - a.z });
            for (float s = 0.0f; s <= len + 0.01f; s += 0.6f) { const float k = len > 0.0f ? s / len : 0.0f; box({ a.x + (b.x - a.x) * k, flr + kCorrH * 0.5f, a.z + (b.z - a.z) * k }, { kCorrHalfW, kCorrH * 0.5f, kCorrHalfW }); }
            const Vector3 m = { 0.5f * (a.x + b.x), flr + 1.0f, 0.5f * (a.z + b.z) };
            lastDungeon.corridorMids.push_back(m);
            links.push_back({ (int)fi, (int)ri });
            // A pressure plate in the middle of some corridors: a falling block or darts.
            if (len > 2.5f && Unit(h) < 0.45f) {
                const uint8_t kind = Unit(h) < 0.55f ? 0 : 1;
                placeAt(PartType::PressurePlate, 0, cellOf(m.x, flr + 0.1f, m.z), kind, "");
                ++lastDungeon.plates; lastDungeon.darts += kind;
            }
        }
        if (rooms.size() - first < 2 && fl == 0) { lastDungeon = DungeonPlan{}; return; } // too little room under good cover
        ++floorsBuilt;
        lastDungeon.deepest = flr;
        if (fl + 1 == floors || rooms.size() - first < 2) break; // (a lone room deeper down: the stairs end in it)
        // Stairs down from a wall of this floor's last room without a door, 10 m to the next floor.
        bool stairs = false;
        for (int k = 0; k < 4 && !stairs; ++k) {
            Room& r = rooms.back();
            if (r.walls & (1 << k)) continue;
            const float sign = (k & 1) ? -1.0f : 1.0f;
            const Vector3 dir = k < 2 ? Vector3{ sign, 0, 0 } : Vector3{ 0, 0, sign };
            const Vector3 start = { r.c.x + dir.x * r.hx, flr, r.c.z + dir.z * r.hz };
            const int n2 = (int)ceilf(kFloorGap / kStepDown);
            // The stairwell's first 4 m (above the room floors under it) must miss this floor's rooms.
            bool clear = true;
            for (size_t q = first; q + 1 < rooms.size(); ++q) {
                const Room& o = rooms[q];
                const Vector3 e = Vector3Add(start, Vector3Scale(dir, 4.0f));
                const float lox = std::min(start.x, e.x) - kCorrHalfW, hix = std::max(start.x, e.x) + kCorrHalfW, loz = std::min(start.z, e.z) - kCorrHalfW, hiz = std::max(start.z, e.z) + kCorrHalfW;
                clear = clear && (hix < o.c.x - o.hx - 1.0f || lox > o.c.x + o.hx + 1.0f || hiz < o.c.z - o.hz - 1.0f || loz > o.c.z + o.hz + 1.0f);
            }
            const Vector3 end = Vector3Add(start, Vector3Scale(dir, kStepRun * n2));
            if (!clear || cover({ end.x, flr - kFloorGap + kCorrH, end.z }) < kCover) continue;
            r.walls |= (uint8_t)(1 << k);
            for (int i = 0; i <= n2; ++i) {
                const Vector3 p = Vector3Add(start, Vector3Scale(dir, kStepRun * i + 0.25f));
                const float y = flr - kStepDown * i;
                box({ p.x, y + 1.2f, p.z }, k < 2 ? Vector3{ 0.4f, 1.2f, kCorrHalfW } : Vector3{ kCorrHalfW, 1.2f, 0.4f });
                if (i > 1 && i < n2 - 1 && i % 5 == 0) lastDungeon.stairPoints.push_back({ p.x, y + 1.0f, p.z });
            }
            { // a torch on the 10th step (3.5 m down: a whole number of cells), at the stairwell's side
                const Vector3 p = Vector3Add(start, Vector3Scale(dir, kStepRun * 10 + 0.25f));
                parts.push_back({ Part{ PartType::Torch, cellOf(p.x - dir.z * 0.5f, flr - 10 * kStepDown + 0.1f, p.z + dir.x * 0.5f), 0, 0 }, "" });
            }
            foot = end; footDir = dir;
            flr = snap(flr - kStepDown * n2);
            stairs = true;
        }
        if (!stairs) break;
    }
    if (floorsBuilt == 0) { lastDungeon = DungeonPlan{}; return; }
    lastDungeon.rooms = (int)rooms.size(); lastDungeon.floors = floorsBuilt;
    lastDungeon.entrance = { x0, floorY, 0.0f }; lastDungeon.floorY = rooms[0].c.y; lastDungeon.links = links;
    // The rooms: carved by shape, a torch each, chests and gas.
    for (size_t i = 0; i < rooms.size(); ++i) {
        const Room& r = rooms[i];
        const float rf = r.c.y;
        box({ r.c.x, rf + kRoomH * 0.5f, r.c.z }, { r.hx, kRoomH * 0.5f, r.hz });
        lastDungeon.shapes |= 1 << r.shape;
        if (r.shape == 1 && r.hx >= 2.5f && r.hz >= 2.5f) { // a pillared hall: four stone pillars, floor to ceiling
            for (int k = 0; k < 4; ++k) {
                const float px = r.c.x + (k & 1 ? 0.5f : -0.5f) * r.hx, pz = r.c.z + (k & 2 ? 0.5f : -0.5f) * r.hz;
                const Int3 c = cellOf(px, rf + 0.1f, pz);
                for (int y = 0; y < (int)(kRoomH / kCellSize); ++y) parts.push_back({ Part{ PartType::Cube, { c.x, c.y + y, c.z }, 0, (uint8_t)stoneMat }, "" });
            }
        } else if (r.shape == 2) { // vaulted: the middle 1.5 m higher
            box({ r.c.x, rf + kRoomH + 0.75f, r.c.z }, { r.hx * 0.6f, 0.75f, r.hz * 0.6f });
        } else if (r.shape == 3) { // a crypt: an alcove 1 m deep in each wall, off its middle (doors take the middle)
            for (int k = 0; k < 4; ++k) {
                const float s = (k & 1) ? -1.0f : 1.0f;
                if (k < 2) box({ r.c.x + s * (r.hx + 0.5f), rf + 1.0f, r.c.z - 0.5f * r.hz }, { 0.5f, 1.0f, 0.6f });
                else box({ r.c.x + 0.5f * r.hx, rf + 1.0f, r.c.z + s * (r.hz + 0.5f) }, { 0.6f, 1.0f, 0.5f });
            }
        }
        lastDungeon.roomCentres.push_back({ r.c.x, rf + 1.0f, r.c.z });
        if (body.trees) body.trees->addClearing(plot.toBody({ r.c.x, rf + kRoomH, r.c.z }), std::max(r.hx, r.hz) + 2.0f); // (roots)
        // A torch in a corner (light: rooms 8 m down are black), and a chest against the far wall now and then: always in
        // a floor's last room, of rare loot in the deepest.
        parts.push_back({ Part{ PartType::Torch, cellOf(r.c.x + r.hx - 0.4f, rf + 0.1f, r.c.z + r.hz - 0.4f), 0, 0 }, "" });
        const bool deepest = i + 1 == rooms.size(), lastOfFloor = deepest || rooms[i + 1].floor != r.floor;
        if (lastOfFloor || Unit(h) < 0.4f)
            parts.push_back({ Part{ PartType::Chest, cellOf(r.c.x - r.hx + 0.6f, rf + 0.1f, r.c.z - 0.5f), 0, 0 }, deepest ? "dungeon_rare" : lastOfFloor ? "ruin_rare" : "ruin_common" });
        // Gas in some rooms past the first (more of them deeper); never the deepest, its chest is reward enough.
        if (i > 0 && !deepest && Unit(h) < (r.floor == 0 ? 0.2f : 0.35f)) {
            parts.push_back({ Part{ PartType::GasVent, cellOf(r.c.x + 0.25f, rf + 0.1f, r.c.z + 0.25f), 0, 0 }, "" });
            ++lastDungeon.vents;
        }
    }
    job.edits.insert(job.edits.end(), cuts.begin(), cuts.end());
    job.parts.insert(job.parts.end(), parts.begin(), parts.end());
}

// Build one site (its tile id for the discovery record; world for its guards, null: none): -1 if the spot doesn't
// suit (water, too steep, someone's plot nearby), else the parts placed.
int BuildSite(Scene* world, Planetoid& body, const PlanetField& pf, uint64_t tile, const Site& site, bool force) {
    const Set& S = sets[site.set];
    const float e = pf.elevation(site.dir);
    if (e < 1.0f) { lastSkip = "the sea"; return -1; }
    if (pf.inlandWaterNear(site.dir, 0.0f)) {
        const float w = pf.waterElevation(site.dir);
        if (w != PlanetField::kNoWater && w > e - 0.5f) { lastSkip = "a river or lake"; return -1; }
    }
    const Vector3 ground = Vector3Scale(site.dir, pf.radius + e);
    // Steepness over the footprint's corners; the structure stands on the lowest of them.
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(site.dir, fabsf(site.dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(t1, site.dir);
    const float half = 0.5f * kCellSize * std::max(S.size.x, S.size.z);
    float lo = e, hi = e;
    for (int k = 0; k < 4; ++k) {
        const Vector3 d = Vector3Normalize(Vector3Add(site.dir, Vector3Add(Vector3Scale(t1, (k & 1 ? half : -half) / pf.radius), Vector3Scale(t2, (k & 2 ? half : -half) / pf.radius))));
        const float ek = pf.elevation(d);
        lo = std::min(lo, ek); hi = std::max(hi, ek);
    }
    (void)hi; // (the footprint's rise is SpotOk's)
    if (const char* why = ""; !force && !SpotOk(pf, site.dir, S, &why)) { lastSkip = why; return -1; }
    if (!body.plots) body.plots = std::make_unique<PlotSet>();
    if (!force && body.plots->plotFor(ground, 40.0f, false)) { lastSkip = "a build nearby"; return -1; }
    Plot* plot = body.plots->plotFor(ground, 0.0f, true);
    if (!plot) return -1;
    const Biome biome = pf.biomeAt(site.dir, e);
    // Turned, centred on the site, its floor on the lowest corner's ground, then sunk `bury`.
    std::vector<Part> parts;
    for (const Part& p : S.parts) parts.push_back(Turned(p, site.quarters));
    std::vector<std::pair<Int3, CellShape>> fp;
    Int3 mn{ 1 << 20, 0, 1 << 20 }, mx{ -(1 << 20), 0, -(1 << 20) };
    for (const Part& p : parts) { Footprint(p, fp); for (auto& c : fp) { mn.x = std::min(mn.x, c.first.x); mn.z = std::min(mn.z, c.first.z); mx.x = std::max(mx.x, c.first.x); mx.z = std::max(mx.z, c.first.z); } }
    const int dyGround = (int)floorf((lo - e) / kCellSize) - (int)lroundf(S.bury / kCellSize);
    const int dy = dyGround + S.yOffset; // parts below the ground (wonders) keep their depth
    // Wonders shape the ground first (#164): each carve turned with the site, from the site's centre at ground level.
    for (const auto& c : S.carves) {
        Vector3 o = c.first;
        for (int k = 0; k < (site.quarters & 3); ++k) o = { o.z, o.y, -o.x };
        TerrainEdit ce{};
        ce.center = plot->toBody({ o.x, o.y + dyGround * kCellSize, o.z });
        ce.radius = c.second;
        ce.carve = true;
        ApplyTerrainEdit(body, ce);
    }
    const Int3 shift = { -(mn.x + mx.x + 1) / 2, dy, -(mn.z + mx.z + 1) / 2 };
    // The ground pad (#164 design): over the footprint and a metre round it, ground more than 0.5 m above the floor is
    // cut away (up to 2 m: a structure on a slope stands in a level cut, not half inside the hill). Flat ground is left
    // alone (cutting down to the floor sank flat sites into a pit); the parts clear their own cells, and the ground under
    // them comes from their foundations. Wonders keep their shaped ground, buried sets the ground over their sunk base.
    if (S.carves.empty() && S.yOffset == 0) {
        const float floorY = (float)(dyGround + (int)lroundf(S.bury / kCellSize)) * kCellSize; // the ground's floor, before burying
        const float x0 = (mn.x + shift.x) * kCellSize - 1.0f, x1 = (mx.x + 1 + shift.x) * kCellSize + 1.0f;
        const float z0 = (mn.z + shift.z) * kCellSize - 1.0f, z1 = (mx.z + 1 + shift.z) * kCellSize + 1.0f;
        TerrainEdit pad{};
        pad.carve = true; pad.box = true;
        pad.center = plot->toBody({ 0.5f * (x0 + x1), floorY + 1.25f, 0.5f * (z0 + z1) }); // floor + 0.5 .. floor + 2 m
        pad.half = { 0.5f * (x1 - x0), 0.75f, 0.5f * (z1 - z0) };
        pad.ax = plot->east; pad.ay = plot->up; pad.az = plot->north;
        pad.radius = Vector3Length(pad.half);
        ApplyTerrainEdit(body, pad);
    }
    // The parts to place (decayed, moved, in the biome's stone), each chest with its loot table. They go down in a build
    // job, a few a frame (PlaceJobs): placing a site's 60-180 parts with their foundations at once was a 0.4 s frame.
    BuildJob job;
    job.body = body.id; job.plot = plot->id; job.seed = site.seed;
    uint32_t h = site.seed ^ 0xdecaf;
    const float topY = (float)std::max(1, S.size.y);
    for (Part p : parts) {
        const bool chest = p.type == PartType::Chest;
        const float rise = (p.origin.y + 0.5f) / topY;
        if (!chest && Unit(h) < S.decay * (0.35f + 1.3f * rise)) continue;   // decay: more near the top
        p.origin = { p.origin.x + shift.x, p.origin.y + shift.y, p.origin.z + shift.z };
        if (IsShape(p.type)) p.material = (uint8_t)BiomeStone(p.material, biome);
        std::string table;
        if (chest)
            for (const Loot& L : S.loots) {
                const Int3 c = TurnedCell(L.cell, site.quarters);
                const Int3 at = { c.x + shift.x, c.y + shift.y, c.z + shift.z };
                // The chest's origin is its cell after turning (Turned moved the origin the same way).
                if (at.x == p.origin.x && at.y == p.origin.y && at.z == p.origin.z) { table = L.table; break; }
            }
        job.parts.push_back({ p, table });
    }
    lastDungeonArgs = { (float)(dyGround + (int)lroundf(S.bury / kCellSize)) * kCellSize, half, BiomeStone(BuildMaterialForItem("stone_brick"), biome) };
    uint32_t dh = Hash((uint64_t)site.seed ^ 0xd0d0ull);
    if (S.dungeon > 0.0f && DungeonsOn() && !noDungeons && (forceDungeon || Unit(dh) < S.dungeon)) // a dungeon under it (#719)
        PlanDungeon(job, body, pf, *plot, lastDungeonArgs.floorY, lastDungeonArgs.half, site.seed, lastDungeonArgs.stone);
    const int planned = (int)job.parts.size();
    jobs.push_back(std::move(job));
    if (force) { RunJob(body, jobs.back(), 1e9); jobs.pop_back(); } // tests (forced sites): all at once
    if (body.trees) body.trees->addClearing(ground, half * 1.42f + 2.0f);
    // Guards (#164: "guard <species> <n> [radius]"): ordinary animals set down round the site, saved like any other.
    if (world) for (const Set::Guard& g : S.guards) Animals_Spawn(*world, body, g.species, ground, g.count, false, g.radius);
    builtSites.push_back({ body.name, tile, site.dir, S.name, S.kind, false });
    return force ? lastPlaced : planned;
}

void FillLoot(Planetoid& body, const PendingLoot& L) {
    auto it = lootTables.find(L.table);
    if (it == lootTables.end()) return;
    uint32_t h = L.seed;
    for (const LootEntry& e : it->second) {
        if (Unit(h) > e.chance) continue;
        const int n = e.lo + (int)(Unit(h) * (e.hi - e.lo + 1));
        if (n > 0 && MagicDiscovery_IsPage(e.item)) { // (#1250) a First Sight page only lies near a ley line (magic items come last: the rolls before stay)
            Vector3 at{};
            bool placed = false;
            if (body.plots) for (auto& p : body.plots->plots)
                if (p->id == L.plot && L.part >= 0 && L.part < (int)p->parts.size()) {
                    const Int3 o = p->parts[L.part].origin;
                    at = p->toBody({ (o.x + 1) * kCellSize, (o.y + 1) * kCellSize, (o.z + 1) * kCellSize });
                    placed = true;
                }
            if (!placed || !MagicDiscovery_LootAllowed(body, at, e.item)) continue;
        }
        if (n > 0) Factory::Insert(body.id, L.plot, L.part, e.item, std::min(n, e.hi));
    }
}

} // namespace

int Structures_SetCount() { Ensure(); return (int)sets.size(); }

void Structures_Reset() { builtTiles.clear(); pending.clear(); builtSites.clear(); foundTiles.clear(); foundKinds.clear(); jobs.clear(); lastCheck = -10.0; }

namespace {
// Sites round one interest centre (a body-local point near the ground): every tile within build range, and the two
// ring sites near the start, built once each.
void BuildAround(Scene& world, Planetoid& body, Vector3 me) {
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) return;
    const float r = Vector3Length(me);
    if (r < 1.0f || r - pf->radius > 400.0f) return;
    const Vector3 up = Vector3Scale(me, 1.0f / r);
    const Vector3 start = StartOf(*pf, body.seed);
    Vector3 t1, t2; TangentAxes(up, t1, t2);
    auto& done = builtTiles[body.name];
    auto consider = [&](uint64_t tile, Vector3 centre) {
        if (done.count(tile)) return;
        const Site site = SiteOf(*pf, tile, centre, pf->radius, start);
        if (!site.any) { done.insert(tile); return; }
        if (Vector3Distance(Vector3Scale(site.dir, pf->radius + pf->elevation(site.dir)), me) > kBuildRange) return; // not yet
        done.insert(tile);
        if (BuildSite(&world, body, *pf, tile, site, false) >= 0) ++sitesBuilt;
    };
    for (int k = 0; k < 25; ++k) { // a 5 x 5 grid 125 m apart: every 512 m tile within build range gets a sample
        const float a = (k % 5 - 2) * 0.5f * kBuildRange / pf->radius, b = (k / 5 - 2) * 0.5f * kBuildRange / pf->radius;
        Vector3 centre;
        const uint64_t tile = TileOf(Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, a), Vector3Scale(t2, b)))), pf->radius, &centre);
        consider(tile, centre);
    }
    if (GroundDistance(up, start, pf->radius) < kRingFar + kBuildRange)
        for (int k = 0; k < 2; ++k) consider(kRingTile + (uint64_t)k, start);
}

// Discovery (#164): a site is found by coming within kFoundRange of it, or by centring it in the telescope (the ground
// the view rests on within 60 m of a site, up to 3 km out). Found sites get a waypoint; undiscovered ones stay hidden.
void Discover(Planetoid& body, uint64_t tile, Vector3 dir, const std::string& name, const std::string& kind) {
    if (!foundTiles[body.name].insert(tile).second) return;
    if (!kind.empty()) foundKinds.insert(kind); // the Help tab's Places (#164)
    Nav_AddWaypoint(body.name, dir, name, false);
    GameUI_Toast(TextFormat("Found: %s (marked on the map and compass)", name.c_str()));
}
void DiscoverNear(Scene& world, Planetoid& body, Vector3 me) {
    for (const BuiltSite& b : builtSites)
        if (b.body == body.name && !foundTiles[body.name].count(b.tile)) {
            const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
            const Vector3 at = Vector3Scale(b.dir, (pf ? pf->radius + pf->elevation(b.dir) : Vector3Length(me)));
            if (Vector3Distance(at, me) < kFoundRange) Discover(body, b.tile, b.dir, b.name, b.kind);
        }
    if (!telescope.active) return;
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) return;
    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
    TerrainHit th;
    if (!RaycastTerrain(world, Ray{ world.camera.position, dir }, 3000.0f, th) || th.planetoid != &body) return;
    const Vector3 hd = Vector3Normalize(th.local);
    Vector3 centre;
    const uint64_t tile = TileOf(hd, pf->radius, &centre);
    const Vector3 start = StartOf(*pf, body.seed);
    for (uint64_t t : { tile, kRingTile, kRingTile + 1 }) {
        const Site s = SiteOf(*pf, t, centre, pf->radius, start);
        if (s.any && GroundDistance(s.dir, hd, pf->radius) < kFoundRange) Discover(body, t, s.dir, sets[s.set].name, sets[s.set].kind);
    }
}
}

bool structuresOn = true; // tests that count every part on the planet turn world structures off (Structures_SetEnabled)
void Structures_SetEnabled(bool on) { structuresOn = on; }

void Structures_Update(Scene& world) {
    if (!structuresOn) return;
    PlaceJobs(); // every frame: a site's parts go down a few at a time
    if (GetTime() - lastCheck < 1.0) return;
    lastCheck = GetTime();
    Ensure();
    if (sets.empty()) return;
    Planetoid* body = player.nearest;
    const bool onFoot = body && !player.flyMode;
    const Vector3 me = body ? body->toLocalOf(player.position, world.renderOrigin) : Vector3{};
    // Only the authority generates (#164 design, MP): a client gets the structures as plot and terrain edits. It builds
    // round the local player and every interest centre (#318: a server's players).
    if (Net_IsAuthority()) {
        if (onFoot) BuildAround(world, *body, me);
        for (const InterestCentre& ic : gInterestCentres)
            if (Planetoid* b = Planetoid::find(ic.body)) BuildAround(world, *b, ic.local);
    }
    if (!onFoot) return;
    DiscoverNear(world, *body, me);
    // Chests near the player get their loot (the factory knows the chest by now: it was built at least a second ago).
    for (size_t i = 0; i < pending.size();) {
        const PendingLoot& L = pending[i];
        Plot* pl = nullptr;
        if (L.body == body->name && body->plots) for (auto& p : body->plots->plots) if (p->id == L.plot) { pl = p.get(); break; }
        if (!pl || L.part < 0 || L.part >= (int)pl->parts.size() || pl->parts[L.part].type != PartType::Chest) {
            if (L.body == body->name) { pending.erase(pending.begin() + i); continue; } // broken before it was filled
            ++i; continue;
        }
        const Int3 o = pl->parts[L.part].origin;
        if (Vector3Distance(pl->toBody({ (o.x + 1) * kCellSize, (o.y + 1) * kCellSize, (o.z + 1) * kCellSize }), me) > kLootRange) { ++i; continue; }
        FillLoot(*body, L);
        pending.erase(pending.begin() + i);
    }
}

void Structures_Save(FILE* f) {
    // Sites still going down are finished first: a half-built site would never be resumed after a load.
    for (BuildJob& j : jobs) if (Planetoid* b = Planetoid::find(j.body)) RunJob(*b, j, 1e9);
    jobs.clear();
    for (const auto& kv : builtTiles) {
        fprintf(f, "structure_tiles %s %zu", kv.first.c_str(), kv.second.size());
        for (uint64_t k : kv.second) fprintf(f, " %llu", (unsigned long long)k);
        fprintf(f, "\n");
    }
    for (const PendingLoot& L : pending) fprintf(f, "structure_loot %s %u %d %s %u\n", L.body.c_str(), L.plot, L.part, L.table.c_str(), L.seed);
    if (!foundKinds.empty()) { fprintf(f, "structure_kinds %zu", foundKinds.size()); for (const auto& k : foundKinds) fprintf(f, " %s", k.c_str()); fprintf(f, "\n"); }
    // Built sites (where, which set) and the ones found (#164 discovery); the name is last (it has spaces).
    for (const BuiltSite& b : builtSites)
        fprintf(f, "structure_site %s %llu %.7f %.7f %.7f %s %s\n", b.body.c_str(), (unsigned long long)b.tile, b.dir.x, b.dir.y, b.dir.z, b.kind.empty() ? "-" : b.kind.c_str(), b.name.c_str());
    for (const auto& kv : foundTiles) {
        fprintf(f, "structure_found %s %zu", kv.first.c_str(), kv.second.size());
        for (uint64_t k : kv.second) fprintf(f, " %llu", (unsigned long long)k);
        fprintf(f, "\n");
    }
}

bool Structures_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "structure_tiles") == 0) {
        char name[128]; size_t n = 0;
        if (fscanf(f, " %127s %zu", name, &n) != 2) return true;
        auto& s = builtTiles[name];
        for (size_t i = 0; i < n; ++i) { unsigned long long k; if (fscanf(f, " %llu", &k) == 1) s.insert(k); }
        return true;
    }
    if (strcmp(tok, "structure_site") == 0) {
        char name[128], kind[32], title[128] = {}; unsigned long long t = 0; Vector3 d{};
        if (fscanf(f, " %127s %llu %f %f %f %31s %127[^\n]", name, &t, &d.x, &d.y, &d.z, kind, title) == 7) {
            std::string tt = title; while (!tt.empty() && (tt.back() == '\r' || tt.back() == ' ')) tt.pop_back();
            builtSites.push_back({ name, (uint64_t)t, d, tt, strcmp(kind, "-") ? kind : "", false });
        }
        return true;
    }
    if (strcmp(tok, "structure_kinds") == 0) {
        size_t n = 0; char k[32];
        if (fscanf(f, " %zu", &n) == 1) for (size_t i = 0; i < n; ++i) if (fscanf(f, " %31s", k) == 1) foundKinds.insert(k);
        return true;
    }
    if (strcmp(tok, "structure_found") == 0) {
        char name[128]; size_t n = 0;
        if (fscanf(f, " %127s %zu", name, &n) != 2) return true;
        auto& s = foundTiles[name];
        for (size_t i = 0; i < n; ++i) { unsigned long long k; if (fscanf(f, " %llu", &k) == 1) s.insert(k); }
        return true;
    }
    if (strcmp(tok, "structure_loot") == 0) {
        char name[128], table[64]; unsigned plot = 0, seed = 0; int part = 0;
        if (fscanf(f, " %127s %u %d %63s %u", name, &plot, &part, table, &seed) == 5) pending.push_back({ name, plot, part, table, seed });
        return true;
    }
    return false;
}

bool Structures_NearestSite(const Planetoid& body, Vector3 local, float range, Vector3& at, std::string& set) {
    Ensure();
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) return false;
    const Vector3 up = Vector3Normalize(local);
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(t1, up);
    float best = range; bool found = false;
    std::unordered_set<uint64_t> seen;
    const int R = (int)ceilf(range / kTile) + 1;
    for (int a = -R; a <= R; ++a)
        for (int b = -R; b <= R; ++b) {
            Vector3 centre;
            const uint64_t tile = TileOf(Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, a * kTile * 0.5f / pf->radius), Vector3Scale(t2, b * kTile * 0.5f / pf->radius)))), pf->radius, &centre);
            if (!seen.insert(tile).second) continue;
            const Site s = SiteOf(*pf, tile, centre, pf->radius, StartOf(*pf, body.seed));
            if (!s.any) continue;
            const Vector3 p = Vector3Scale(s.dir, pf->radius + pf->elevation(s.dir));
            const float d = Vector3Distance(p, local);
            if (d < best) { best = d; at = p; set = sets[s.set].name; found = true; }
        }
    for (int k = 0; k < 2; ++k) { // and the ring sites near the start
        const Site s = RingSite(*pf, StartOf(*pf, body.seed), k);
        if (!s.any) continue;
        const Vector3 p = Vector3Scale(s.dir, pf->radius + pf->elevation(s.dir));
        const float d = Vector3Distance(p, local);
        if (d < best) { best = d; at = p; set = sets[s.set].name; found = true; }
    }
    return found;
}

int Structures_SitesNear(const Planetoid& body, Vector3 local, float range, std::vector<SiteInfo>& out) {
    out.clear();
    Ensure();
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf) return 0;
    const Vector3 up = Vector3Normalize(local);
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(t1, up);
    std::unordered_set<uint64_t> seen;
    auto add = [&](const Site& s) {
        if (!s.any) return;
        const Vector3 p = Vector3Scale(s.dir, pf->radius + pf->elevation(s.dir));
        const float d = Vector3Distance(p, local);
        if (d > range) return;
        const Set& S = sets[s.set];
        uint32_t dh = Hash((uint64_t)s.seed ^ 0xd0d0ull);
        const bool dungeon = S.dungeon > 0.0f && DungeonsOn() && Unit(dh) < S.dungeon; // (as Build's roll)
        out.push_back({ p, S.name, S.kind, dungeon, d });
    };
    const int R = (int)ceilf(range / kTile) + 1;
    for (int a = -R; a <= R; ++a)
        for (int b = -R; b <= R; ++b) {
            Vector3 centre;
            const uint64_t tile = TileOf(Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, a * kTile * 0.5f / pf->radius), Vector3Scale(t2, b * kTile * 0.5f / pf->radius)))), pf->radius, &centre);
            if (!seen.insert(tile).second) continue;
            add(SiteOf(*pf, tile, centre, pf->radius, StartOf(*pf, body.seed)));
        }
    for (int k = 0; k < 2; ++k) add(RingSite(*pf, StartOf(*pf, body.seed), k)); // and the ring sites near the start
    std::sort(out.begin(), out.end(), [](const SiteInfo& x, const SiteInfo& y) { return x.distance < y.distance; });
    out.erase(std::unique(out.begin(), out.end(), [](const SiteInfo& x, const SiteInfo& y) { return Vector3Distance(x.at, y.at) < 1.0f; }), out.end());
    return (int)out.size();
}

bool Structures_Test(Scene& world, Planetoid& body, Vector3 nearLocal) {
    // Every set built 30 m from the player on each side (forced: slope and nearby plots ignored), each turned a
    // different way; then the loot filled (the chests are within reach). Checks: parts placed (allowing for decay),
    // one chest each with pending loot, turned footprints the same shape.
    Ensure();
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf || sets.empty()) { printf("bench: structures: no sets or no planet -> FAIL\n"); return false; }
    const Vector3 up = Vector3Normalize(nearLocal);
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    noDungeons = true; // the surface sets only (--test-dungeon builds one)
    bool ok = true;
    const size_t pending0 = pending.size();
    for (int i = 0; i < (int)sets.size(); ++i) {
        Site s; s.any = true; s.set = i; s.quarters = i & 3; s.seed = 1234u + i;
        const Vector3 t2 = Vector3CrossProduct(t1, up); // four sides of the player, 30 m out (chests report only near the player)
        const Vector3 side = i % 4 == 0 ? t1 : i % 4 == 1 ? Vector3Negate(t1) : i % 4 == 2 ? t2 : Vector3Negate(t2); // (was i == 0.. 3: sets 5 on all went to one side, on top of each other)
        s.dir = Vector3Normalize(Vector3Add(up, Vector3Scale(side, (30.0f + 22.0f * (i / 4)) / pf->radius))); // 22 m apart: the camps and the farmstead are ~12 m across
        const int placed = BuildSite(&world, body, *pf, kRingTile + 100 + (uint64_t)i, s, true); // (a test tile id: never a real one)
        // Turning keeps the shape: the footprint cells of the turned set are the turned cells of the set.
        std::vector<std::pair<Int3, CellShape>> fp;
        std::unordered_set<Int3> a, b;
        // Cells turn about the corner, (x, y, z) -> (z, y, -x-1); origins turn as vectors (Turned).
        for (const Part& p : sets[i].parts) {
            Footprint(p, fp);
            for (auto& c : fp) { Int3 q = c.first; for (int k = 0; k < s.quarters; ++k) q = { q.z, q.y, -q.x - 1 }; a.insert(q); }
            Footprint(Turned(p, s.quarters), fp);
            for (auto& c : fp) b.insert(c.first);
        }
        const bool rigid = a == b;
        const int expect = (int)sets[i].parts.size();
        const bool good = placed >= expect * (1.0f - sets[i].decay * 1.7f) - 1 && rigid;
        printf("bench: structures: %s (%s): %d of %d parts placed (decay %.2f), turned %d quarters, rigid %d\n", sets[i].name.c_str(), sets[i].kind.c_str(), placed, expect, sets[i].decay, s.quarters, (int)rigid);
        ok = ok && good;
    }
    noDungeons = false;
    for (size_t k = pending0; k < pending.size(); ++k) testFilled.push_back({ pending[k].plot, pending[k].part });
    const int chests = (int)(pending.size() - pending0);
    int lootCount = 0;
    for (size_t k = pending0; k < pending.size(); ++k) lootCount += (int)lootTables[pending[k].table].size();
    printf("bench: structures: %d sets, %d chests with loot waiting (%d table entries), %zu loot tables\n", (int)sets.size(), chests, lootCount, lootTables.size());
    (void)world;
    return ok && chests == (int)sets.size() && lootCount > 0;
}


void Structures_TestFill(Planetoid& body) { // fill every pending chest on this body now (the test's approach)
    Ensure();
    for (size_t i = 0; i < pending.size();) {
        if (pending[i].body != body.name) { ++i; continue; }
        FillLoot(body, pending[i]);

        pending.erase(pending.begin() + i);
    }
}
int Structures_TestCheck(Planetoid& body, int& items) { // chests filled by the test that hold something now
    int n = 0; items = 0;
    for (auto& c : testFilled) {
        const auto contents = Factory::ChestContents(body.id, c.first, c.second);
        int k = 0; for (auto& s : contents) k += s.second;
        if (k > 0) ++n;
        printf("bench: structures: chest plot %u part %d holds %d (%zu kinds)\n", c.first, c.second, k, contents.size());
        items += k;
    }
    return n;
}
std::string Structures_PlacesText(bool discoveredOnly) {
    // The Help tab's Places (#164): one paragraph per kind of site, only the kinds found when "Discovered only" is on.
    static const std::pair<const char*, const char*> kKinds[] = {
        { "camp", "Camps: lean-tos round a cold fire, left by whoever passed before. A chest of supplies; at a hunter's camp, hides and meat, and the thornbacks that drove the hunters off." },
        { "ruin", "Ruins: what is left of an earlier people - watchtowers, house walls, farmsteads, desert shrines. Built in the stone of their biome, half fallen, with a chest inside. Their stone and planks can be broken for materials." },
        { "wreck", "Crash sites: wreckage of a lander, scattered over a scorched patch. Tech parts in its chest: plates, wire, circuits, sometimes a power cell." },
        { "wonder", "Wonders: stone circles, crystal geodes, sinkholes. Odd places where the ground holds something - gems, quartz, a long fall." },
    };
    Ensure();
    std::string out;
    int shown = 0;
    for (const auto& k : kKinds) {
        if (discoveredOnly && !foundKinds.count(k.first)) continue;
        out += k.second; out += "\n\n"; ++shown;
    }
    if (discoveredOnly && shown < (int)(sizeof kKinds / sizeof kKinds[0])) out += TextFormat("(%d kinds of place not found yet.)", (int)(sizeof kKinds / sizeof kKinds[0]) - shown);
    return out;
}

const char* Structures_LastSkip() { return lastSkip; }

bool Structures_TestDungeon(Scene& world, Planetoid& body, Vector3 nearLocal, Vector3& eye, Vector3& target) {
    // A ruined watchtower with its dungeon forced, 40 m from the player (#719), built at once. Checks: 2-3 floors down to
    // ~30 m, 4-12 rooms; the air (density under ISO) at every room's middle, every corridor's middle, along the ramp and the
    // stairs; a door (two gates) in every doorway, all placed; the traps planned all placed; one dungeon_rare chest waiting.
    // Tries a few seeds for one with both kinds of trap. eye/target (body-local): a view across the first room by its torch.
    Ensure();
    if (!DungeonsOn()) { printf("bench: dungeon: RAYLIBGAME_DUNGEONS=0: off -> SKIP\n"); return true; }
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    int set = -1;
    for (int i = 0; i < (int)sets.size(); ++i) if (sets[i].name == "Ruined watchtower") set = i;
    if (!pf || set < 0) { printf("bench: dungeon: no watchtower set -> FAIL\n"); return false; }
    const Vector3 up = Vector3Normalize(nearLocal);
    Vector3 t1, t2; TangentAxes(up, t1, t2);
    Site s; s.any = true; s.set = set; s.quarters = 0;
    s.dir = Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 40.0f / pf->radius)));
    // The ruin first (no dungeon), then the first dungeon seed (of a few, planned dry: cheap) with a falling block, darts
    // and gas, built under it.
    s.seed = 4242u;
    noDungeons = true;
    const size_t pending0 = pending.size();
    const int placed = BuildSite(&world, body, *pf, kRingTile + 200, s, true);
    noDungeons = false;
    Plot* plot = body.plots ? body.plots->plotFor(Vector3Scale(s.dir, pf->radius + pf->elevation(s.dir)), 15.0f, false) : nullptr;
    if (placed < 0 || !plot) { printf("bench: dungeon: the ruin not built (%s) -> FAIL\n", lastSkip); return false; }
    const DungeonArgs A = lastDungeonArgs;
    uint32_t pick = 4242u;
    for (uint32_t k = 0; k < 40; ++k) {
        const uint32_t sd = 4242u + 7919u * k;
        BuildJob dry;
        PlanDungeon(dry, body, *pf, *plot, A.floorY, A.half, sd, A.stone);
        if (lastDungeon.rooms > 0 && lastDungeon.vents > 0 && lastDungeon.darts > 0 && lastDungeon.darts < lastDungeon.plates) { pick = sd; break; }
    }
    BuildJob job;
    job.body = body.id; job.plot = plot->id; job.seed = pick;
    PlanDungeon(job, body, *pf, *plot, A.floorY, A.half, pick, A.stone);
    if (lastDungeon.rooms == 0) { printf("bench: dungeon: none fits under the ruin -> FAIL\n"); return false; }
    RunJob(body, job, 1e9);
    auto air = [&](Vector3 plotPos) { return body.densityAt(plot->toBody(plotPos)) < ISO_LEVEL; };
    const DungeonPlan& D = lastDungeon;
    int roomsOpen = 0, linksOpen = 0, rampOpen = 0, stairsOpen = 0;
    for (const Vector3& c : D.roomCentres) roomsOpen += air(c);
    for (const Vector3& m : D.corridorMids) linksOpen += air(m);
    for (const Vector3& p : D.rampPoints) rampOpen += air(p);
    for (const Vector3& p : D.stairPoints) stairsOpen += air(p);
    int gates = 0, plates = 0, vents = 0;
    for (const Part& p : plot->parts) { gates += p.type == PartType::Gate; plates += p.type == PartType::PressurePlate; vents += p.type == PartType::GasVent; }
    int rare = 0;
    for (size_t k = pending0; k < pending.size(); ++k) rare += pending[k].table == "dungeon_rare";
    const float depth = D.entrance.y - D.deepest;
    const bool ok = D.floors >= 2 && D.floors <= 3 && depth > 18.0f && depth < 31.0f && D.rooms >= 4 && D.rooms <= 12 && roomsOpen == D.rooms &&
                    linksOpen == (int)D.corridorMids.size() && rampOpen == (int)D.rampPoints.size() && !D.stairPoints.empty() && stairsOpen == (int)D.stairPoints.size() &&
                    D.doors == 2 * (int)D.links.size() && gates == 2 * D.doors && plates == D.plates && vents == D.vents && rare == 1;
    printf("bench: dungeon: seed %u: %d floors, %.1f m down at the deepest; %d rooms (%d open, shapes 0x%x), %d corridors (%d open at the middle), ramp open %d/%d, stairs open %d/%d; "
           "%d doors (%d gates placed); %d plates (%d darts; %d placed), %d gas vents (%d placed); %d dungeon_rare chest -> %s\n",
           pick, D.floors, depth, D.rooms, roomsOpen, D.shapes, (int)D.corridorMids.size(), linksOpen, rampOpen, (int)D.rampPoints.size(), stairsOpen, (int)D.stairPoints.size(),
           D.doors, gates, D.plates, D.darts, plates, D.vents, vents, rare, ok ? "PASS" : "FAIL");
    const Vector3 c0 = D.roomCentres[0];
    eye = plot->toBody({ c0.x - 1.5f, c0.y + 0.6f, c0.z - 1.5f });
    target = plot->toBody({ c0.x + 3.0f, c0.y - 0.6f, c0.z + 3.0f });
    return ok;
}

bool Structures_TestPlacement(const Planetoid& body, float range, uint64_t& hash) {
    const auto clk0 = std::chrono::steady_clock::now();
    // The world's own placement round the start (#164 design), checked independently of SpotOk: every site within
    // `range` on dry ground (25 samples over its footprint grown by the water margin) with at most 3 m of rise over 25
    // samples of the footprint, none within the keep-out of the start, and the ring's camp and ruin between kRingNear
    // and kRingFar. `hash` folds every site's tile, set and position (run twice: placement is deterministic).
    Ensure();
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    if (!pf || sets.empty()) return false;
    const Vector3 start = StartOf(*pf, body.seed);
    Vector3 t1, t2; TangentAxes(start, t1, t2);
    std::unordered_set<uint64_t> seen;
    std::vector<std::pair<uint64_t, Site>> found;
    const int R = (int)ceilf(range / (kTile * 0.5f)) + 1;
    for (int a = -R; a <= R; ++a)
        for (int b = -R; b <= R; ++b) {
            Vector3 centre;
            const uint64_t tile = TileOf(Vector3Normalize(Vector3Add(start, Vector3Add(Vector3Scale(t1, a * kTile * 0.5f / pf->radius), Vector3Scale(t2, b * kTile * 0.5f / pf->radius)))), pf->radius, &centre);
            if (!seen.insert(tile).second) continue;
            const Site s = SiteOf(*pf, tile, centre, pf->radius, start);
            if (s.any && GroundDistance(s.dir, start, pf->radius) <= range) found.push_back({ tile, s });
        }
    for (int k = 0; k < 2; ++k) { const Site s = RingSite(*pf, start, k); if (s.any) found.push_back({ kRingTile + (uint64_t)k, s }); }
    bool ok = true;
    hash = 1469598103934665603ull;
    int wet = 0, steep = 0, close = 0;
    for (const auto& [tile, s] : found) {
        const Set& S = sets[s.set];
        Vector3 u1, u2; TangentAxes(s.dir, u1, u2);
        const float half = 0.5f * kCellSize * std::max(S.size.x, S.size.z);
        float lo = 1e9f, hi = -1e9f; bool dry = true;
        for (int i = 0; i < 25; ++i) {
            const float u = (float)(i % 5 - 2) * 0.5f, v = (float)(i / 5 - 2) * 0.5f;
            auto at = [&](float r) { return Vector3Normalize(Vector3Add(s.dir, Vector3Add(Vector3Scale(u1, u * r / pf->radius), Vector3Scale(u2, v * r / pf->radius)))); };
            dry = dry && !Wet(*pf, at(half + kWaterMargin));
            const float e = pf->elevation(at(half));
            lo = std::min(lo, e); hi = std::max(hi, e);
        }
        const float d = GroundDistance(s.dir, start, pf->radius);
        if (!dry) ++wet;
        if (hi - lo > 3.0f) ++steep;
        if (d < kSpawnKeepOut) ++close;
        if (!dry || hi - lo > 3.0f || d < kSpawnKeepOut)
            printf("bench: structures: placement: %s %.0f m from the start: %s%s%s\n", S.name.c_str(), d, dry ? "" : "wet ", hi - lo > 3.0f ? TextFormat("rise %.1f m ", hi - lo) : "", d < kSpawnKeepOut ? "too close" : "");
        for (uint64_t x : { tile, (uint64_t)s.set, (uint64_t)(uint32_t)(s.dir.x * 1e6f), (uint64_t)(uint32_t)(s.dir.z * 1e6f) }) hash = (hash ^ x) * 1099511628211ull;
    }
    ok = wet == 0 && steep == 0 && close == 0;
    const Site camp = RingSite(*pf, start, 0), ruin = RingSite(*pf, start, 1);
    const float dc = camp.any ? GroundDistance(camp.dir, start, pf->radius) : -1.0f, dr = ruin.any ? GroundDistance(ruin.dir, start, pf->radius) : -1.0f;
    const bool ring = camp.any && ruin.any && sets[camp.set].kind == "camp" && sets[ruin.set].kind == "ruin" &&
                      dc >= kRingNear - 1.0f && dc <= kRingFar + 1.0f && dr >= kRingNear - 1.0f && dr <= kRingFar + 1.0f;
    printf("bench: structures: placement: scanned in %.0f ms\n", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - clk0).count());
    printf("bench: structures: placement: %d sites within %.0f m of the start (%d wet, %d steep, %d inside %.0f m); near the start a %s at %.0f m and a %s at %.0f m; layout hash %016llx -> %s\n",
           (int)found.size(), range, wet, steep, close, kSpawnKeepOut, camp.any ? sets[camp.set].name.c_str() : "(no camp)", dc,
           ruin.any ? sets[ruin.set].name.c_str() : "(no ruin)", dr, (unsigned long long)hash, ok && ring ? "ok" : "WRONG");
    return ok && ring;
}
int Structures_SitesBuilt() { return sitesBuilt; }
