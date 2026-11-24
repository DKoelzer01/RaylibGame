#include "mapfeatures.h"
#include "jobsystem.h"
#include "plots.h"
#include "profiler.h"
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
    // Coarser levels: per block the most common kind of its columns (any part at all keeps the block filled).
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
        for (size_t i = 0; i < next.size(); ++i) {
            int best = 0;
            for (int k = 1; k < (int)MapKind::Count; ++k) if (next[i][k] > (best ? next[i][best] : 0u)) best = k;
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
        // Zoomed in (a cell 2+ pixels across): a dark edge where the pixel's corners land on another part.
        if (L == 0 && mpp <= 0.5f * kCellSize) {
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

void MapFeatures_Reset() { cache.clear(); }

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
#endif
