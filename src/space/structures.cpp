#include "structures.h"
#include "buildtool.h"
#include "factory.h"
#include "inventory.h"
#include "player.h"
#include "plots.h"
#include "scene.h"
#include "trees.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <map>
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
};
std::vector<Set> sets;
bool loaded = false;

// Sites: one per 512 m tile at most (#164).
constexpr float kTile = 512.0f, kBuildRange = 250.0f, kLootRange = 60.0f, kSiteOdds = 0.35f;
std::unordered_map<std::string, std::unordered_set<uint64_t>> builtTiles;   // body name -> tiles decided (built or skipped)
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
    static const int stone = BuildMaterialForItem("stone");
    if (mat != stone) return mat;
    const char* m = nullptr;
    switch (b) {
        case Biome::Desert: case Biome::Badlands: case Biome::Savanna: m = "sandstone"; break;
        case Biome::Tundra: case Biome::Alpine: case Biome::Taiga: case Biome::Polar: m = "granite"; break;
        case Biome::Jungle: case Biome::Swamp: m = "limestone"; break;
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
// The tile's site (if it has one): a set the biome there allows, chosen by weight; a spot in the tile's middle 60%.
Site SiteOf(const PlanetField& pf, uint64_t tile, Vector3 centreDir, float radius) {
    Site s;
    uint32_t h = Hash(tile * 0x9e3779b97f4a7c15ull ^ (uint64_t)pf.foliageSeed * 0x632be59bd9b4e019ull ^ 0x5717c7u);
    if (Unit(h) > kSiteOdds || sets.empty()) return s;
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(centreDir, fabsf(centreDir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(t1, centreDir);
    const float a = (Unit(h) - 0.5f) * 0.6f * kTile, b = (Unit(h) - 0.5f) * 0.6f * kTile;
    s.dir = Vector3Normalize(Vector3Add(centreDir, Vector3Add(Vector3Scale(t1, a / radius), Vector3Scale(t2, b / radius))));
    const float e = pf.elevation(s.dir);
    if (e < 1.0f) return s;                                                 // in the sea: no site
    if (pf.inlandWaterNear(s.dir, 0.0f)) {                                  // or in a river or lake
        const float w = pf.waterElevation(s.dir);
        if (w != PlanetField::kNoWater && w > e - 0.5f) return s;
    }
    const int biome = (int)pf.biomeAt(s.dir, e);
    float total = 0.0f;
    for (const Set& st : sets) if (st.biomes.empty() || std::find(st.biomes.begin(), st.biomes.end(), biome) != st.biomes.end()) total += st.weight;
    if (total <= 0.0f) return s;
    float pick = Unit(h) * total;
    for (int i = 0; i < (int)sets.size(); ++i) {
        const Set& st = sets[i];
        if (!(st.biomes.empty() || std::find(st.biomes.begin(), st.biomes.end(), biome) != st.biomes.end())) continue;
        if ((pick -= st.weight) <= 0.0f) { s.set = i; break; }
    }
    if (s.set < 0) return s;
    s.quarters = (int)(Unit(h) * 4.0f) & 3;
    s.seed = h;
    s.any = true;
    return s;
}

// Build one site: false if the spot doesn't suit (water, too steep, someone's plot nearby).
int BuildSite(Planetoid& body, const PlanetField& pf, const Site& site, bool force) {
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
    if (!force && hi - lo > 2.5f) { lastSkip = "too steep"; return -1; }
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
    uint32_t h = site.seed ^ 0xdecaf;
    int placed = 0;
    const float topY = (float)std::max(1, S.size.y);
    for (Part p : parts) {
        const bool chest = p.type == PartType::Chest;
        const float rise = (p.origin.y + 0.5f) / topY;
        if (!chest && Unit(h) < S.decay * (0.35f + 1.3f * rise)) continue;   // decay: more near the top
        p.origin = { p.origin.x + shift.x, p.origin.y + shift.y, p.origin.z + shift.z };
        if (IsShape(p.type)) p.material = (uint8_t)BiomeStone(p.material, biome);
        const int idx = plot->place(p);
        if (idx < 0) continue;
        BuildTool_LayFoundation(body, *plot, p);
        ++placed;
        if (chest) {
            for (const Loot& L : S.loots) {
                const Int3 c = TurnedCell(L.cell, site.quarters);
                const Int3 at = { c.x + shift.x, c.y + shift.y, c.z + shift.z };
                // The chest's origin is its cell after turning (Turned moved the origin the same way).
                if (at.x == p.origin.x && at.y == p.origin.y && at.z == p.origin.z) { pending.push_back({ body.name, plot->id, idx, L.table, Hash(site.seed + idx) }); break; }
            }
        }
    }
    if (body.trees) body.trees->addClearing(ground, half * 1.42f + 2.0f);
    return placed;
}

void FillLoot(Planetoid& body, const PendingLoot& L) {
    auto it = lootTables.find(L.table);
    if (it == lootTables.end()) return;
    uint32_t h = L.seed;
    for (const LootEntry& e : it->second) {
        if (Unit(h) > e.chance) continue;
        const int n = e.lo + (int)(Unit(h) * (e.hi - e.lo + 1));
        if (n > 0) Factory::Insert(body.id, L.plot, L.part, e.item, std::min(n, e.hi));
    }
}

} // namespace

int Structures_SetCount() { Ensure(); return (int)sets.size(); }

void Structures_Reset() { builtTiles.clear(); pending.clear(); lastCheck = -10.0; }

void Structures_Update(Scene& world) {
    if (GetTime() - lastCheck < 1.0) return;
    lastCheck = GetTime();
    Ensure();
    Planetoid* body = player.nearest;
    if (!body || player.flyMode) return;
    const auto* pf = dynamic_cast<const PlanetField*>(body->field.get());
    if (!pf || sets.empty()) return;
    const Vector3 me = body->toLocal((player.position - world.renderOrigin).toFloat());
    const float r = Vector3Length(me);
    if (r < 1.0f || r - pf->radius > 400.0f) return;
    const Vector3 up = Vector3Scale(me, 1.0f / r);
    // The tiles round the player, sampled on a grid out to the build range.
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(t1, up);
    auto& done = builtTiles[body->name];
    for (int k = 0; k < 25; ++k) { // a 5 x 5 grid 125 m apart: every 512 m tile within build range gets a sample
        const float a = (k % 5 - 2) * 0.5f * kBuildRange / pf->radius, b = (k / 5 - 2) * 0.5f * kBuildRange / pf->radius;
        Vector3 centre;
        const uint64_t tile = TileOf(Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, a), Vector3Scale(t2, b)))), pf->radius, &centre);
        if (done.count(tile)) continue;
        const Site site = SiteOf(*pf, tile, centre, pf->radius);
        if (!site.any) { done.insert(tile); continue; }
        if (Vector3Distance(Vector3Scale(site.dir, pf->radius + pf->elevation(site.dir)), me) > kBuildRange) continue; // not yet
        done.insert(tile);
        if (BuildSite(*body, *pf, site, false) >= 0) ++sitesBuilt;
    }
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
    for (const auto& kv : builtTiles) {
        fprintf(f, "structure_tiles %s %zu", kv.first.c_str(), kv.second.size());
        for (uint64_t k : kv.second) fprintf(f, " %llu", (unsigned long long)k);
        fprintf(f, "\n");
    }
    for (const PendingLoot& L : pending) fprintf(f, "structure_loot %s %u %d %s %u\n", L.body.c_str(), L.plot, L.part, L.table.c_str(), L.seed);
}

bool Structures_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "structure_tiles") == 0) {
        char name[128]; size_t n = 0;
        if (fscanf(f, " %127s %zu", name, &n) != 2) return true;
        auto& s = builtTiles[name];
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
            const Site s = SiteOf(*pf, tile, centre, pf->radius);
            if (!s.any) continue;
            const Vector3 p = Vector3Scale(s.dir, pf->radius + pf->elevation(s.dir));
            const float d = Vector3Distance(p, local);
            if (d < best) { best = d; at = p; set = sets[s.set].name; found = true; }
        }
    return found;
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
    bool ok = true;
    const size_t pending0 = pending.size();
    for (int i = 0; i < (int)sets.size(); ++i) {
        Site s; s.any = true; s.set = i; s.quarters = i & 3; s.seed = 1234u + i;
        const Vector3 t2 = Vector3CrossProduct(t1, up); // four sides of the player, 30 m out (chests report only near the player)
        const Vector3 side = i == 0 ? t1 : i == 1 ? Vector3Negate(t1) : i == 2 ? t2 : Vector3Negate(t2);
        s.dir = Vector3Normalize(Vector3Add(up, Vector3Scale(side, (30.0f + 12.0f * (i / 4)) / pf->radius)));
        const int placed = BuildSite(body, *pf, s, true);
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
const char* Structures_LastSkip() { return lastSkip; }
int Structures_SitesBuilt() { return sitesBuilt; }
