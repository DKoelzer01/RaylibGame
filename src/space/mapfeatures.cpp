#include "mapfeatures.h"
#include "jobsystem.h"
#include "plots.h"
#include "profiler.h"
#include "waterfield.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>

MapKind MapKindOf(PartType t) {
    switch (t) {
        case PartType::Conveyor: case PartType::Cable: case PartType::Pipe: case PartType::Outlet:
        case PartType::Shaft: case PartType::Cogwheel: case PartType::LargeCogwheel: case PartType::Gearbox:
        case PartType::HandCrank: case PartType::StressGauge:
            return MapKind::Transport;
        case PartType::Chest: case PartType::EssenceTank: case PartType::FuelTank:
            return MapKind::Storage;
        case PartType::Generator: case PartType::MvGenerator: case PartType::HvGenerator: case PartType::Battery:
        case PartType::Boiler: case PartType::WaterWheel: case PartType::Windmill: case PartType::FusionReactor:
        case PartType::ManaDynamo: case PartType::SteamEngine: case PartType::Dynamo: case PartType::Motor:
            return MapKind::Power;
        case PartType::Farmland: case PartType::Snare:
            return MapKind::Farm;
        case PartType::Fence: case PartType::Gate: case PartType::Ladder:
            return MapKind::Structure;
        case PartType::MbController: case PartType::MbCasing: case PartType::MbHatch:
            return MapKind::Machine; // a multiblock reads as one machine
        default: break;
    }
    if (IsMachine(t)) return MapKind::Machine;
    if (IsShape(t)) return MapKind::Structure;
    return t == PartType::Count ? MapKind::None : MapKind::Other; // rocket parts, torches, beds, banners...
}

const char* MapKindName(MapKind k) {
    static const char* names[(int)MapKind::Count] = { "none", "structure", "machine", "transport", "storage", "power", "farm", "other" };
    return names[(int)k];
}

namespace {

// Which part a column shows: machines over transport over the rest; the higher part on a tie.
int Rank(MapKind k) {
    switch (k) {
        case MapKind::Machine: case MapKind::Power: case MapKind::Storage: return 3;
        case MapKind::Transport: return 2;
        case MapKind::Farm: case MapKind::Other: return 1;
        default: return 0;
    }
}

struct Snapshot {
    uint32_t id = 0, version = 0;
    DVec3 origin;
    Vector3 east{}, up{}, north{};
    std::vector<Part> parts;
};

std::shared_ptr<const PlotMap> Build(const Snapshot& s) {
    auto m = std::make_shared<PlotMap>();
    m->id = s.id;
    m->version = s.version;
    m->origin = s.origin;
    m->east = s.east;
    m->north = s.north;
    struct Cell { int x, y, z; int part; MapKind kind; };
    std::vector<Cell> cells;
    std::vector<std::pair<Int3, CellShape>> fp;
    int x0 = INT32_MAX, z0 = INT32_MAX, x1 = INT32_MIN, z1 = INT32_MIN;
    for (size_t i = 0; i < s.parts.size(); ++i) {
        const Part& p = s.parts[i];
        const MapKind k = MapKindOf(p.type);
        if (k == MapKind::None) continue;
        fp.clear();
        Plot::Footprint(p, fp);
        if (fp.empty()) continue;
        float cx = 0.0f, cy = 0.0f, cz = 0.0f;
        int ex0 = INT32_MAX, ex1 = INT32_MIN, ez0 = INT32_MAX, ez1 = INT32_MIN;
        for (const auto& c : fp) {
            cells.push_back({ c.first.x, c.first.y, c.first.z, (int)i, k });
            x0 = std::min(x0, c.first.x); x1 = std::max(x1, c.first.x);
            z0 = std::min(z0, c.first.z); z1 = std::max(z1, c.first.z);
            ex0 = std::min(ex0, c.first.x); ex1 = std::max(ex1, c.first.x);
            ez0 = std::min(ez0, c.first.z); ez1 = std::max(ez1, c.first.z);
            cx += c.first.x + 0.5f; cy += c.first.y + 0.5f; cz += c.first.z + 0.5f;
        }
        if (IsMachine(p.type) || k == MapKind::Power || k == MapKind::Storage) {
            const float n = (float)fp.size();
            const DVec3 at = s.origin + Vector3Scale(s.east, cx / n * kCellSize) + Vector3Scale(s.up, cy / n * kCellSize) + Vector3Scale(s.north, cz / n * kCellSize);
            const char* item = PartItem(p);
            if (item && *item) m->icons.push_back({ item, Vector3Normalize(at.toFloat()), std::max(ex1 - ex0 + 1, ez1 - ez0 + 1) * kCellSize });
        }
    }
    if (cells.empty()) return m;
    const int w = x1 - x0 + 1, h = z1 - z0 + 1;
    m->x0 = x0;
    m->z0 = z0;
    PlotMap::Level l0{ w, h, std::vector<uint8_t>((size_t)w * h, 0) };
    m->part.assign((size_t)w * h, -1);
    std::vector<int> top((size_t)w * h, INT32_MIN);
    for (const Cell& c : cells) {
        const size_t i = (size_t)(c.z - z0) * w + (c.x - x0);
        const MapKind was = (MapKind)l0.kind[i];
        const bool better = was == MapKind::None || Rank(c.kind) > Rank(was) || (Rank(c.kind) == Rank(was) && c.y > top[i]);
        if (!better) continue;
        l0.kind[i] = (uint8_t)c.kind;
        m->part[i] = c.part;
        top[i] = c.y;
    }
    // Coarser levels: per block the highest-ranked kind among its columns, so a machine on a floor or a conveyor row still
    // shows zoomed out (any part at all keeps the block filled).
    using Counts = std::array<uint32_t, (int)MapKind::Count>;
    std::vector<Counts> cnt((size_t)w * h, Counts{});
    for (size_t i = 0; i < cnt.size(); ++i) cnt[i][l0.kind[i]] = 1;
    m->levels.push_back(std::move(l0));
    int lw = w, lh = h;
    while (lw > 1 || lh > 1) {
        const int nw = (lw + 1) / 2, nh = (lh + 1) / 2;
        std::vector<Counts> next((size_t)nw * nh, Counts{});
        PlotMap::Level l{ nw, nh, std::vector<uint8_t>((size_t)nw * nh, 0) };
        for (int z = 0; z < lh; ++z)
            for (int x = 0; x < lw; ++x) {
                Counts& d = next[(size_t)(z / 2) * nw + x / 2];
                const Counts& c = cnt[(size_t)z * lw + x];
                for (int k = 0; k < (int)MapKind::Count; ++k) d[k] += c[k];
            }
        for (size_t i = 0; i < next.size(); ++i) { // the highest-ranked kind present, the most common one on a tie
            int best = 0;
            for (int k = 1; k < (int)MapKind::Count; ++k) {
                if (!next[i][k]) continue;
                if (!best || Rank((MapKind)k) > Rank((MapKind)best) || (Rank((MapKind)k) == Rank((MapKind)best) && next[i][k] > next[i][best])) best = k;
            }
            l.kind[i] = (uint8_t)best;
        }
        cnt.swap(next);
        m->levels.push_back(std::move(l));
        lw = nw;
        lh = nh;
    }
    const float hx = 0.5f * w * kCellSize, hz = 0.5f * h * kCellSize;
    m->centre = s.origin + Vector3Scale(s.east, (x0 * kCellSize) + hx) + Vector3Scale(s.north, (z0 * kCellSize) + hz);
    m->radius = sqrtf(hx * hx + hz * hz);
    return m;
}

Snapshot Take(const Plot& p) {
    Snapshot s;
    s.id = p.id;
    s.version = p.version;
    s.origin = p.origin;
    s.east = p.east;
    s.up = p.up;
    s.north = p.north;
    s.parts = p.parts;
    return s;
}

struct Slot {
    std::atomic<bool> done{ false };
    std::shared_ptr<const PlotMap> map;
};
struct Entry {
    uint32_t builtVersion = UINT32_MAX;          // the version `map` or the build in flight is of
    std::shared_ptr<const PlotMap> map;
    std::shared_ptr<Slot> building;
};
std::unordered_map<uint64_t, std::unordered_map<uint32_t, Entry>> cache; // body id -> plot id -> footprint (main thread)

// resources/data/map_colours.txt: "kind | r,g,b". The defaults are the Okabe-Ito set, told apart in every colour vision.
std::array<Color, (int)MapKind::Count> colours;
bool coloursLoaded = false;
void LoadColours() {
    coloursLoaded = true;
    colours = { Color{ 0, 0, 0, 0 }, { 176, 176, 170, 255 }, { 230, 159, 0, 255 }, { 240, 228, 66, 255 },
                { 204, 121, 167, 255 }, { 86, 180, 233, 255 }, { 0, 158, 115, 255 }, { 213, 94, 0, 255 } };
    std::ifstream f("resources/data/map_colours.txt");
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        char name[32] = {};
        int r = 0, g = 0, b = 0;
        if (sscanf(line.c_str(), " %31[a-z_] | %d , %d , %d", name, &r, &g, &b) != 4) continue;
        for (int k = 1; k < (int)MapKind::Count; ++k)
            if (strcmp(name, MapKindName((MapKind)k)) == 0) colours[k] = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 };
    }
}

// Plot cells (continuous) of a body-local point.
inline void ToCells(const PlotMap& m, const DVec3& p, double& u, double& v) {
    const DVec3 d = p - m.origin;
    u = (d.x * m.east.x + d.y * m.east.y + d.z * m.east.z) / kCellSize;
    v = (d.x * m.north.x + d.y * m.north.y + d.z * m.north.z) / kCellSize;
}

inline int PartAt(const PlotMap& m, double u, double v) {
    const int x = (int)floor(u) - m.x0, z = (int)floor(v) - m.z0;
    const PlotMap::Level& l = m.levels[0];
    if (x < 0 || z < 0 || x >= l.w || z >= l.h) return -1;
    return m.part[(size_t)z * l.w + x];
}

} // namespace

bool MapFeatures_Update(const Planetoid& body) {
    PROFILE_SCOPE("nav: plots");
    if (!coloursLoaded) LoadColours();
    bool landed = false;
    auto& plots = cache[body.id];
    if (!body.plots) {
        landed = !plots.empty();
        plots.clear();
        return landed;
    }
    // Plots that went away.
    for (auto it = plots.begin(); it != plots.end();) {
        bool alive = false;
        for (const auto& p : body.plots->plots) if (p && p->id == it->first) { alive = true; break; }
        if (alive) { ++it; continue; }
        landed |= it->second.map != nullptr;
        it = plots.erase(it);
    }
    for (const auto& p : body.plots->plots) {
        if (!p) continue;
        Entry& e = plots[p->id];
        if (e.building && e.building->done.load()) {
            e.map = e.building->map;
            e.building.reset();
            landed = true;
        }
        if (e.building || e.builtVersion == p->version) continue;
        // One build in flight per plot: a plot changing every frame (the player building) re-builds once the last one ends.
        auto slot = std::make_shared<Slot>();
        e.building = slot;
        e.builtVersion = p->version;
        auto snap = std::make_shared<Snapshot>(Take(*p));
        JobSystem::tagNext(JobTag::NavMap);
        WorldJobs().submit(250.0f, [slot, snap] { slot->map = Build(*snap); slot->done.store(true); });
    }
    return landed;
}

std::shared_ptr<const MapPlots> MapFeatures_PlotsNear(const Planetoid& body, const DVec3& centre, float radius) {
    auto it = cache.find(body.id);
    if (it == cache.end()) return nullptr;
    std::shared_ptr<MapPlots> s;
    for (const auto& kv : it->second) {
        const auto& m = kv.second.map;
        if (!m || m->levels.empty()) continue;
        if ((m->centre - centre).length() > radius + m->radius) continue;
        if (!s) {
            s = std::make_shared<MapPlots>();
            for (int k = 0; k < (int)MapKind::Count; ++k) s->colour[k] = colours[k];
        }
        s->plots.push_back(m);
    }
    return s;
}

bool MapFeatures_Sample(const MapPlots& s, const DVec3& p, Vector3 east, Vector3 north, float mpp, Color& out) {
    for (const auto& mp : s.plots) {
        const PlotMap& m = *mp;
        if ((p - m.centre).length() > m.radius + mpp) continue;
        double u, v;
        ToCells(m, p, u, v);
        // The finest level whose blocks are at least a pixel across, so a row of single cells never falls between pixels.
        int L = 0;
        while (L + 1 < (int)m.levels.size() && (float)(1 << L) * kCellSize < mpp) ++L;
        const PlotMap::Level& l = m.levels[L];
        const int x = (int)floor((u - m.x0) / (1 << L)), z = (int)floor((v - m.z0) / (1 << L));
        if (x < 0 || z < 0 || x >= l.w || z >= l.h) continue;
        const MapKind k = (MapKind)l.kind[(size_t)z * l.w + x];
        if (k == MapKind::None) continue;
        Color c = s.colour[(int)k];
        // Zoomed in (a cell 4+ pixels across; at 2, nearly every pixel touches a cell edge): a dark edge where the pixel's
        // corners land on another part.
        if (L == 0 && mpp <= 0.25f * kCellSize) {
            const int self = PartAt(m, u, v);
            const float h = 0.5f * mpp;
            for (int i = 0; i < 4; ++i) {
                const float a = (i & 1) ? h : -h, b = (i & 2) ? h : -h;
                double cu, cv;
                ToCells(m, p + Vector3Add(Vector3Scale(east, a), Vector3Scale(north, b)), cu, cv);
                if (PartAt(m, cu, cv) != self) { c = { (unsigned char)(c.r * 2 / 5), (unsigned char)(c.g * 2 / 5), (unsigned char)(c.b * 2 / 5), 255 }; break; }
            }
        }
        out = c;
        return true;
    }
    return false;
}

std::shared_ptr<const PlotMap> MapFeatures_BuildNow(const Plot& plot) {
    if (!coloursLoaded) LoadColours();
    return Build(Take(plot));
}

Color MapFeatures_Colour(MapKind k) {
    if (!coloursLoaded) LoadColours();
    return colours[(int)k];
}

// ---------------------------------------------------------------------------------------------
// Terrain edits (#1795 part 2)
// ---------------------------------------------------------------------------------------------

namespace {

constexpr float kEditTile = 32.0f; // m: an edit tile's side on the cube-sphere grid
constexpr int kEditCols = 64;      // columns along a tile's side (0.5 m)
const Vector3 kFaceN[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
const Vector3 kFaceU[6] = { { 0, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 0 } };
const Vector3 kFaceV[6] = { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 } };

int TilesPerFace(float R) { return std::max(1, (int)ceilf(0.5f * PI * R / kEditTile)); }

// A direction's tile (face, i, j) and its place in the tile (0..1 along each side), as navigation's TileKey splits faces.
uint64_t EditKey(Vector3 dir, float R, float* fu = nullptr, float* fv = nullptr) {
    int face = 0;
    float best = -2.0f;
    for (int f = 0; f < 6; ++f) { const float d = Vector3DotProduct(dir, kFaceN[f]); if (d > best) { best = d; face = f; } }
    const int n = TilesPerFace(R);
    const float da = 0.5f * PI / n;
    const float a = (atanf(Vector3DotProduct(dir, kFaceU[face]) / best) + 0.25f * PI) / da, b = (atanf(Vector3DotProduct(dir, kFaceV[face]) / best) + 0.25f * PI) / da;
    const int i = std::clamp((int)floorf(a), 0, n - 1), j = std::clamp((int)floorf(b), 0, n - 1);
    if (fu) *fu = std::clamp(a - i, 0.0f, 0.99999f);
    if (fv) *fv = std::clamp(b - j, 0.0f, 0.99999f);
    return ((uint64_t)face << 48) | ((uint64_t)i << 24) | (uint64_t)j;
}

Vector3 EditColumnDir(uint64_t key, float R, int cu, int cv) {
    const int face = (int)(key >> 48), i = (int)((key >> 24) & 0xFFFFFF), j = (int)(key & 0xFFFFFF);
    const float da = 0.5f * PI / TilesPerFace(R);
    const float a = -0.25f * PI + (i + (cu + 0.5f) / kEditCols) * da, b = -0.25f * PI + (j + (cv + 0.5f) / kEditCols) * da;
    return Vector3Normalize(Vector3Add(kFaceN[face], Vector3Add(Vector3Scale(kFaceU[face], tanf(a)), Vector3Scale(kFaceV[face], tanf(b)))));
}

// Per tile, an order-free hash of the delta chunks reaching it (their keys and the snapshot's chunk pointers: an edit makes
// a new chunk, untouched chunks are shared between snapshots). A chunk marks the tiles under its 8 corners.
std::unordered_map<uint64_t, uint64_t> ScanTiles(const TerrainDeltas& d, float R) {
    std::unordered_map<uint64_t, uint64_t> h;
    const float side = kDeltaChunk * kVoxelSize;
    for (const auto& kv : d.chunks) {
        int cx, cy, cz;
        TerrainDeltas::Unpack(kv.first, cx, cy, cz);
        uint64_t v = kv.first * 0x9E3779B97F4A7C15ull ^ (uint64_t)(uintptr_t)kv.second.get();
        v ^= v >> 29; v *= 0xBF58476D1CE4E5B9ull; v ^= v >> 32;
        uint64_t seen[8];
        int n = 0;
        for (int c = 0; c < 8; ++c) {
            const Vector3 p = { (cx + (c & 1)) * side, (cy + ((c >> 1) & 1)) * side, (cz + (c >> 2)) * side };
            const float len = Vector3Length(p);
            if (len < 1.0f) continue;
            const uint64_t k = EditKey(Vector3Scale(p, 1.0f / len), R);
            bool dup = false;
            for (int q = 0; q < n; ++q) dup |= seen[q] == k;
            if (dup) continue;
            seen[n++] = k;
            h[k] ^= v;
        }
    }
    return h;
}

// The edited surface down each column: from just over the generated ground, down through air (a pit) or up through solid
// (a fill) in 1 m steps to the first crossing, then 6 halvings (~1.5 cm).
std::shared_ptr<const EditTile> SampleTile(uint64_t key, const DensityField& field, const PlanetField& pf, const TerrainDeltas* d) {
    auto t = std::make_shared<EditTile>();
    t->key = key;
    const float R = pf.radius;
    t->elev.resize((size_t)kEditCols * kEditCols);
    t->change.resize(t->elev.size());
    for (int cv = 0; cv < kEditCols; ++cv)
        for (int cu = 0; cu < kEditCols; ++cu) {
            const Vector3 dir = EditColumnDir(key, R, cu, cv);
            const float g = pf.elevation(dir, 0.5f);
            auto surface = [&](const TerrainDeltas* deltas) {
                auto solid = [&](float e) { return field.density(Vector3Scale(dir, R + e), deltas) >= ISO_LEVEL; };
                float lo, hi;
                float e = g + 1.0f;
                if (solid(e)) { // filled over the ground: climb to the air
                    int k = 0;
                    while (k < 24 && solid(e + 1.0f)) { e += 1.0f; ++k; }
                    lo = e; hi = e + 1.0f;
                } else {        // the ground or a pit: descend to the solid
                    int k = 0;
                    while (k < 40 && !solid(e - 1.0f)) { e -= 1.0f; ++k; }
                    hi = e; lo = e - 1.0f;
                }
                for (int s = 0; s < 6; ++s) { const float m = 0.5f * (lo + hi); if (solid(m)) lo = m; else hi = m; }
                return 0.5f * (lo + hi);
            };
            // The change against the same search without the edits (the generator's elevation() and the density surface
            // differ by up to ~0.3 m where the 3D band shapes the ground).
            const size_t i = (size_t)cv * kEditCols + cu;
            t->elev[i] = surface(d);
            t->change[i] = t->elev[i] - surface(nullptr);
        }
    return t;
}

struct ScanSlot { std::atomic<bool> done{ false }; std::unordered_map<uint64_t, uint64_t> hashes; };
struct TileSlot { std::atomic<bool> done{ false }; std::shared_ptr<const EditTile> tile; };
struct EditTileState {
    uint64_t hash = 0, sampledHash = 0;
    double sampledAt = -1e18;
    std::shared_ptr<TileSlot> building;
    std::shared_ptr<const EditTile> tile;
};
struct EditState {
    const TerrainDeltas* scanned = nullptr; // the snapshot the last scan read
    double scanAt = -1e18;
    std::shared_ptr<ScanSlot> scanning;
    std::unordered_map<uint64_t, EditTileState> tiles;
    std::shared_ptr<const MapEdits> set;
};
std::unordered_map<uint64_t, EditState> editStates; // body id -> its edit tiles (main thread)

} // namespace

bool MapFeatures_UpdateEdits(const Planetoid& body, double now) {
    PROFILE_SCOPE("nav: edits");
    auto pf = std::dynamic_pointer_cast<const PlanetField>(body.field);
    if (!pf) return false;
    EditState& st = editStates[body.id];
    const std::shared_ptr<const TerrainDeltas> snap = body.edits;
    // A scan of a newer snapshot, at most every 2 s.
    if (st.scanning && st.scanning->done.load()) {
        for (auto& kv : st.scanning->hashes) st.tiles[kv.first].hash = kv.second;
        for (auto& kv : st.tiles) if (!st.scanning->hashes.count(kv.first)) kv.second.hash = 0; // its edits went (a reload)
        st.scanning.reset();
    }
    if (!st.scanning && snap && snap.get() != st.scanned && now - st.scanAt >= 2.0) {
        auto slot = std::make_shared<ScanSlot>();
        st.scanning = slot;
        st.scanned = snap.get();
        st.scanAt = now;
        const float R = pf->radius;
        JobSystem::tagNext(JobTag::NavMap);
        WorldJobs().submit(150.0f, [slot, snap, R] { slot->hashes = ScanTiles(*snap, R); slot->done.store(true); });
    }
    // Re-sample the dirty tiles (a changed hash), each at most every kEditResampleS seconds.
    bool landed = false;
    for (auto& kv : st.tiles) {
        EditTileState& t = kv.second;
        if (t.building && t.building->done.load()) {
            t.tile = t.building->tile;
            t.building.reset();
            landed = true;
        }
        if (t.building || t.hash == t.sampledHash || now - t.sampledAt < kEditResampleS || !snap) continue;
        auto slot = std::make_shared<TileSlot>();
        t.building = slot;
        t.sampledHash = t.hash;
        t.sampledAt = now;
        const uint64_t key = kv.first;
        auto field = body.field;
        JobSystem::tagNext(JobTag::NavMap);
        WorldJobs().submit(150.0f, [slot, key, field, pf, snap] { slot->tile = SampleTile(key, *field, *pf, snap.get()); slot->done.store(true); });
    }
    if (landed || (!st.set && !st.tiles.empty())) {
        auto set = std::make_shared<MapEdits>();
        set->radius = pf->radius;
        for (const auto& kv : st.tiles) if (kv.second.tile) set->tiles[kv.first] = kv.second.tile;
        st.set = set;
    }
    return landed;
}

std::shared_ptr<const MapEdits> MapFeatures_Edits(const Planetoid& body) {
    auto it = editStates.find(body.id);
    if (it == editStates.end() || !it->second.set || it->second.set->tiles.empty()) return nullptr;
    return it->second.set;
}

bool MapFeatures_EditedAt(const MapEdits& e, Vector3 dir, float& elev, float& change) {
    float fu, fv;
    const uint64_t key = EditKey(dir, e.radius, &fu, &fv);
    const auto it = e.tiles.find(key);
    if (it == e.tiles.end()) return false;
    // Bilinear between column centres (a nearest column drew 0.5 m steps, and their hillshade a grid, when zoomed in);
    // clamped at the tile's edge.
    const float x = std::clamp(fu * kEditCols - 0.5f, 0.0f, kEditCols - 1.0f), y = std::clamp(fv * kEditCols - 0.5f, 0.0f, kEditCols - 1.0f);
    const int x0 = std::min((int)x, kEditCols - 2), y0 = std::min((int)y, kEditCols - 2);
    const float tx = x - x0, ty = y - y0;
    auto lerp2 = [&](const std::vector<float>& a) {
        const float* r0 = &a[(size_t)y0 * kEditCols + x0];
        const float* r1 = r0 + kEditCols;
        return (r0[0] * (1 - tx) + r0[1] * tx) * (1 - ty) + (r1[0] * (1 - tx) + r1[1] * tx) * ty;
    };
    elev = lerp2(it->second->elev);
    change = lerp2(it->second->change);
    return true;
}

// ---------------------------------------------------------------------------------------------
// Simulated water (#1795 part 3)
// ---------------------------------------------------------------------------------------------

namespace {

inline int64_t ColKey(int i, int j) { return ((int64_t)(i + (1 << 20)) << 21) | (int64_t)(j + (1 << 20)); }

std::shared_ptr<const MapWater> SampleWater(const WaterField::Snapshot& snap, uint64_t body) {
    auto w = std::make_shared<MapWater>();
    const float cs = WaterField::CellSize();
    for (const WaterField::Site& s : snap.sites) {
        if (s.body != body || s.bricks.empty()) continue;
        MapWaterSite ms;
        ms.origin = s.origin; ms.e = s.e; ms.n = s.n; ms.u = s.u; ms.lo = s.lo; ms.hi = s.hi;
        for (const auto& kv : s.bricks) {
            const int bi = (int)(kv.first >> 22) - 1024, bj = (int)((kv.first >> 11) & 2047) - 1024, bk = (int)(kv.first & 2047) - 1024;
            const WaterField::Brick& b = *kv.second;
            for (int y = 0; y < WaterField::kBrick; ++y)
                for (int x = 0; x < WaterField::kBrick; ++x) {
                    float top = -1e9f, depth = 0.0f;
                    uint8_t kind = 0;
                    for (int z = WaterField::kBrick - 1; z >= 0; --z) {
                        const int i = WaterField::CellIdx(x, y, z);
                        if (!b.fill[i] || (b.flags[i] & WaterField::kSolid)) continue;
                        const float t = bk * WaterField::kBrick + z + (b.shut[i] + b.fill[i]) / 255.0f; // water stands on the cell's ground share (#454)
                        if (t > top) { top = t; kind = (uint8_t)WaterField::CellKind(b.flags[i]); }
                        depth += b.fill[i] / 255.0f * cs;
                    }
                    if (depth <= 0.0f) continue;
                    MapWaterSite::Col& c = ms.cols[ColKey(bi * WaterField::kBrick + x, bj * WaterField::kBrick + y)];
                    if (top > c.top) { c.top = top; c.kind = kind; }
                    c.depth += depth; // (bricks stacked on one column add up)
                }
        }
        for (auto it = ms.cols.begin(); it != ms.cols.end();) it = it->second.depth < 0.1f ? ms.cols.erase(it) : std::next(it); // films
        if (!ms.cols.empty()) w->sites.push_back(std::move(ms));
    }
    return w;
}

struct WaterSlot { std::atomic<bool> done{ false }; std::shared_ptr<const MapWater> water; };
struct WaterState {
    uint64_t seq = UINT64_MAX; // the snapshot sampled last
    double at = -1e18;
    std::shared_ptr<WaterSlot> building;
    std::shared_ptr<const MapWater> water;
};
std::unordered_map<uint64_t, WaterState> waterStates; // body id -> its sampled water (main thread)

} // namespace

bool MapFeatures_UpdateWater(const Planetoid& body, double now) {
    PROFILE_SCOPE("nav: water");
    WaterState& st = waterStates[body.id];
    bool landed = false;
    if (st.building && st.building->done.load()) {
        landed = (st.water && !st.water->sites.empty()) || !st.building->water->sites.empty();
        st.water = st.building->water;
        st.building.reset();
    }
    WaterField::SnapshotPtr snap = WaterField::Current();
    if (st.building || !snap || snap->seq == st.seq) return landed;
    // The first water at once, then at most every kWaterResampleS seconds.
    const bool first = !st.water || st.water->sites.empty();
    if (!first && now - st.at < kWaterResampleS) return landed;
    if (first) {
        bool any = false;
        for (const auto& s : snap->sites) any |= s.body == body.id && !s.bricks.empty();
        if (!any) { st.seq = snap->seq; return landed; }
    }
    auto slot = std::make_shared<WaterSlot>();
    st.building = slot;
    st.seq = snap->seq;
    st.at = now;
    const uint64_t id = body.id;
    JobSystem::tagNext(JobTag::NavMap);
    WorldJobs().submit(150.0f, [slot, snap, id] { slot->water = SampleWater(*snap, id); slot->done.store(true); });
    return landed;
}

std::shared_ptr<const MapWater> MapFeatures_Water(const Planetoid& body) {
    auto it = waterStates.find(body.id);
    if (it == waterStates.end() || !it->second.water || it->second.water->sites.empty()) return nullptr;
    return it->second.water;
}

bool MapFeatures_WaterAt(const MapWater& w, Vector3 p, float& depth, int& kind) {
    const float cs = WaterField::CellSize();
    for (const MapWaterSite& s : w.sites) {
        if (p.x < s.lo.x - 2.0f || p.y < s.lo.y - 2.0f || p.z < s.lo.z - 2.0f || p.x > s.hi.x + 2.0f || p.y > s.hi.y + 2.0f || p.z > s.hi.z + 2.0f) continue;
        const Vector3 d = Vector3Subtract(p, s.origin);
        const float fx = Vector3DotProduct(d, s.e) / cs, fy = Vector3DotProduct(d, s.n) / cs, fz = Vector3DotProduct(d, s.u) / cs;
        const auto it = s.cols.find(ColKey((int)floorf(fx), (int)floorf(fy)));
        if (it == s.cols.end()) continue;
        const float over = (it->second.top - fz) * cs;
        if (over < -0.3f) continue; // water below this ground (another level of a cave, or ground that rose)
        depth = std::max(over, 0.0f);
        kind = it->second.kind;
        return true;
    }
    return false;
}

void MapFeatures_Reset() { cache.clear(); editStates.clear(); waterStates.clear(); }

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
// (#1795) A plot at a 30 degree heading on a 60 km planet: a floor of cubes, a furnace on it and a row of conveyors. Each
// shows in its kind's colour at its true place and heading, thin rows stay visible zoomed out, and part edges are dark.
TEST(MapFeatures, PlotFootprintAtItsHeadingInKindColours) {
    const double R = 60000.0;
    const Vector3 up = Vector3Normalize({ 0.3f, 0.8f, 0.5f });
    const Vector3 e0 = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, up)), n0 = Vector3CrossProduct(up, e0);
    const float a = 30.0f * DEG2RAD;
    Plot plot;
    plot.up = up;
    plot.east = Vector3Add(Vector3Scale(e0, cosf(a)), Vector3Scale(n0, sinf(a)));
    plot.north = Vector3CrossProduct(up, plot.east);
    plot.origin = DVec3(up) * R;
    for (int x = 0; x < 10; ++x)
        for (int z = -1; z < 3; ++z) ASSERT_GE(plot.place({ PartType::Cube, { x, 0, z }, 0, 0 }), 0);
    const int furnace = plot.place({ PartType::Furnace, { 0, 1, 0 }, 0, 0 });
    ASSERT_GE(furnace, 0);
    for (int x = 4; x < 10; ++x) ASSERT_GE(plot.place({ PartType::Conveyor, { x, 1, 0 }, 0, 0 }), 0);

    auto map = MapFeatures_BuildNow(plot);
    ASSERT_FALSE(map->levels.empty());
    MapPlots s;
    s.plots.push_back(map);
    for (int k = 0; k < (int)MapKind::Count; ++k) s.colour[k] = colours[k];
    // A point at plot cell (u, v) (continuous), seen in an image whose axes are the unturned east and north.
    auto at = [&](float u, float v) { return plot.origin + Vector3Scale(plot.east, u * kCellSize) + Vector3Scale(plot.north, v * kCellSize) + Vector3Scale(up, 0.6f); };
    auto kindAt = [&](float u, float v, float mpp) {
        Color c;
        if (!MapFeatures_Sample(s, at(u, v), e0, n0, mpp, c)) return MapKind::None;
        for (int k = 1; k < (int)MapKind::Count; ++k) if (c.r == s.colour[k].r && c.g == s.colour[k].g && c.b == s.colour[k].b) return (MapKind)k;
        return MapKind::Count; // an outline (darkened)
    };
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(plot.parts[furnace], fp);
    const float fu = fp[0].first.x + 0.5f, fv = fp[0].first.z + 0.5f;
    EXPECT_EQ(kindAt(fu, fv, 0.1f), MapKind::Machine) << "the furnace shows over the floor";
    EXPECT_EQ(kindAt(6.5f, 0.5f, 0.1f), MapKind::Transport) << "a conveyor";
    EXPECT_EQ(kindAt(6.5f, 2.5f, 0.1f), MapKind::Structure) << "the bare floor";
    EXPECT_EQ(kindAt(12.5f, 0.5f, 0.1f), MapKind::None) << "past the plot's east end";
    EXPECT_EQ(kindAt(6.5f, -3.0f, 0.1f), MapKind::None) << "south of it";
    // The heading: 4.5 m along the plot's own east is on the plot, 4.5 m along the image's east falls off its south side.
    {
        Color c;
        EXPECT_TRUE(MapFeatures_Sample(s, plot.origin + Vector3Scale(plot.east, 4.5f) + Vector3Scale(plot.north, 0.25f), e0, n0, 0.1f, c));
        EXPECT_FALSE(MapFeatures_Sample(s, plot.origin + Vector3Scale(e0, 4.5f) + Vector3Scale(n0, 0.25f), e0, n0, 0.1f, c));
    }
    // Zoomed in, the edge between two parts is dark; the middle of a part is not.
    EXPECT_EQ(kindAt(6.02f, 0.5f, 0.1f), MapKind::Count) << "an outline at the edge of two conveyors";
    EXPECT_EQ(kindAt(6.5f, 0.5f, 0.1f), MapKind::Transport);
    EXPECT_EQ(kindAt(6.02f, 0.5f, 0.25f), MapKind::Transport) << "no outlines at 2 pixels a cell (every pixel would touch an edge)";
    // Zoomed out, a block shows its highest-ranked part: the furnace and the conveyors over the floor they stand on.
    EXPECT_EQ(kindAt(fu, fv, 2.0f), MapKind::Machine) << "the furnace's 2 m block";
    EXPECT_EQ(kindAt(5.5f, 2.5f, 2.0f), MapKind::Transport) << "a block of floor with a conveyor row through it";
    // Zoomed out (4 m a pixel, 8 cells a block): the plot is still drawn wherever it has parts.
    EXPECT_NE(kindAt(6.5f, 0.5f, 4.0f), MapKind::None);
    EXPECT_EQ(kindAt(30.0f, 0.5f, 4.0f), MapKind::None);
    // A machine gets an icon at its centre.
    ASSERT_EQ(map->icons.size(), 1u);
    EXPECT_EQ(map->icons[0].item, PartItem(plot.parts[furnace]));
    // Kinds.
    EXPECT_EQ(MapKindOf(PartType::Chest), MapKind::Storage);
    EXPECT_EQ(MapKindOf(PartType::Generator), MapKind::Power);
    EXPECT_EQ(MapKindOf(PartType::Shaft), MapKind::Transport);
    EXPECT_EQ(MapKindOf(PartType::Farmland), MapKind::Farm);
    EXPECT_EQ(MapKindOf(PartType::Slope), MapKind::Structure);
    EXPECT_EQ(MapKindOf(PartType::MbCasing), MapKind::Machine);
    EXPECT_EQ(MapKindOf(PartType::Count), MapKind::None);
}

// (#1795 part 2) An edit tile's columns: each column's direction falls back in its own tile and column (random tiles, a
// cube edge and corner), and neighbouring columns are about 0.5 m apart on a 60 km planet.
TEST(MapFeatures, EditTileColumnsMapBackToThemselves) {
    const float R = 60000.0f;
    std::vector<Vector3> dirs = { Vector3Normalize({ 1.0f, 1.0f, 1.0f }), Vector3Normalize({ 1.0f, 1.0f, 0.2f }), { 0, 1, 0 } };
    for (int k = 0; k < 6; ++k) dirs.push_back(Vector3Normalize({ sinf(k * 1.7f), cosf(k * 2.3f), sinf(k * 0.9f + 1.0f) }));
    for (const Vector3 d : dirs) {
        const uint64_t key = EditKey(d, R);
        int bad = 0;
        float worstGap = 0.0f;
        for (int cv = 0; cv < kEditCols; cv += 7)
            for (int cu = 0; cu < kEditCols; cu += 7) {
                const Vector3 c = EditColumnDir(key, R, cu, cv);
                float fu, fv;
                if (EditKey(c, R, &fu, &fv) != key || (int)(fu * kEditCols) != cu || (int)(fv * kEditCols) != cv) ++bad;
                if (cu + 1 < kEditCols) worstGap = std::max(worstGap, Vector3Distance(c, EditColumnDir(key, R, cu + 1, cv)) * R);
            }
        EXPECT_EQ(bad, 0) << "columns that map elsewhere, near " << d.x << " " << d.y << " " << d.z;
        EXPECT_LT(worstGap, 0.6f);
        EXPECT_GT(worstGap, 0.3f);
    }
}
#endif
