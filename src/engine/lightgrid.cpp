#include "lightgrid.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <tuple>

static_assert(LightGrid::kBrick == 16, "BrickOf uses >> 4 and & 15");

namespace {
struct Move { int dx, dy, dz, cost; };
// The 26 neighbours with chamfer costs (8 per metre: face 8, edge 8*sqrt2 ~ 11, corner 8*sqrt3 ~ 14).
struct Moves {
    Move m[26];
    Moves() {
        int n = 0;
        for (int dz = -1; dz <= 1; ++dz)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dx = -1; dx <= 1; ++dx) {
                    const int axes = (dx != 0) + (dy != 0) + (dz != 0);
                    if (axes == 0) continue;
                    m[n++] = { dx, dy, dz, axes == 1 ? 8 : axes == 2 ? 11 : 14 };
                }
    }
};
const Moves kMoves;

uint64_t CellKey(int x, int y, int z) {
    const uint64_t m = (1u << 21) - 1;
    return ((uint64_t)(x & m)) | ((uint64_t)(y & m) << 21) | ((uint64_t)(z & m) << 42);
}
// The brighter of two arrivals: by level, then by colour, so every order of work ends on the same grid.
bool Brighter(const LightCell& a, const LightCell& b) {
    return std::tie(a.level, a.r, a.g, a.b) > std::tie(b.level, b.r, b.g, b.b);
}
}

LightGrid::LightGrid(OpacityFn opacity) : opacity_(std::move(opacity)) {}

LightCell LightGrid::Emission(float radiusM, float r, float g, float b) {
    LightCell c;
    const float top = std::max(r, std::max(g, b));
    if (top <= 0.0f || radiusM <= 0.0f) return c;
    auto hue = [&](float v) { return (uint8_t)std::lround(std::min(std::max(v / top, 0.0f), 1.0f) * 255.0f); };
    c.r = hue(r); c.g = hue(g); c.b = hue(b);
    c.level = (uint8_t)std::lround(std::min(radiusM, 31.0f) * kStep); // 31 m * 8 = 248 <= 255
    return c;
}

uint64_t LightGrid::Key(int bx, int by, int bz) { return CellKey(bx, by, bz); }

LightGrid::Brick* LightGrid::Find(int bx, int by, int bz) const {
    if (last_ && last_->bx == bx && last_->by == by && last_->bz == bz) return const_cast<Brick*>(last_);
    auto it = bricks_.find(Key(bx, by, bz));
    if (it == bricks_.end()) return nullptr;
    last_ = &it->second;
    return const_cast<Brick*>(last_);
}

LightGrid::Brick& LightGrid::Get(int bx, int by, int bz) {
    if (Brick* b = Find(bx, by, bz)) return *b;
    Brick& b = bricks_[Key(bx, by, bz)];
    b.bx = bx; b.by = by; b.bz = bz;
    std::memset(b.opaque, 0, sizeof(b.opaque)); // b.light starts unlit (LightCell's initializers)
    if (opacity_) opacity_(bx, by, bz, b.opaque);
    last_ = &b;
    return b;
}

LightGrid::Brick& LightGrid::BrickOf(int x, int y, int z, int& index) {
    index = (x & 15) + (y & 15) * kBrick + (z & 15) * kBrick * kBrick;
    return Get(x >> 4, y >> 4, z >> 4);
}

LightCell LightGrid::At(int x, int y, int z) const {
    const Brick* b = Find(x >> 4, y >> 4, z >> 4);
    return b ? b->light[(x & 15) + (y & 15) * kBrick + (z & 15) * kBrick * kBrick] : LightCell{};
}

bool LightGrid::OpaqueAt(int x, int y, int z) const {
    const Brick* b = Find(x >> 4, y >> 4, z >> 4);
    return b && b->opaque[(x & 15) + (y & 15) * kBrick + (z & 15) * kBrick * kBrick];
}

bool LightGrid::Opaque(int x, int y, int z) {
    int i;
    return BrickOf(x, y, z, i).opaque[i] != 0;
}

// A diagonal step needs every cell of the box it spans open, so light does not slip between two solid cells that only
// touch at an edge or a corner.
bool LightGrid::CanStep(int x, int y, int z, int dx, int dy, int dz) {
    for (int oz = 0; oz <= (dz != 0); ++oz)
        for (int oy = 0; oy <= (dy != 0); ++oy)
            for (int ox = 0; ox <= (dx != 0); ++ox) {
                if (!ox && !oy && !oz) continue;
                if (Opaque(x + ox * dx, y + oy * dy, z + oz * dz)) return false;
            }
    return true;
}

LightCell LightGrid::EmissionAt(int x, int y, int z) const {
    auto it = sources_.find(CellKey(x, y, z));
    LightCell e;
    if (it == sources_.end()) return e;
    for (const LightCell& s : it->second) if (Brighter(s, e)) e = s;
    return e;
}

void LightGrid::Darken(int x, int y, int z) {
    int i;
    Brick& b = BrickOf(x, y, z, i);
    const uint8_t old = b.light[i].level;
    if (old == 0) return;
    b.light[i] = {};
    MarkDirty(b);
    removeQ_.push_back({ x, y, z, old });
}

void LightGrid::Relight(int x, int y, int z) {
    const LightCell e = EmissionAt(x, y, z);
    int i;
    Brick& b = BrickOf(x, y, z, i);
    if (Brighter(e, b.light[i])) { b.light[i] = e; MarkDirty(b); }
    addQ_.push_back({ x, y, z });
}

void LightGrid::AddSource(int x, int y, int z, LightCell light) {
    if (light.level == 0) return;
    sources_[CellKey(x, y, z)].push_back(light);
    Relight(x, y, z);
}

void LightGrid::RemoveSource(int x, int y, int z, LightCell light) {
    auto it = sources_.find(CellKey(x, y, z));
    if (it == sources_.end()) return;
    auto& v = it->second;
    auto f = std::find(v.begin(), v.end(), light);
    if (f == v.end()) return;
    v.erase(f);
    if (v.empty()) sources_.erase(it);
    Darken(x, y, z);
    Relight(x, y, z); // what is left there (a second source in the cell); the removal pass refills from its edge too
}

void LightGrid::ClearSources() {
    std::vector<uint64_t> keys;
    for (auto& kv : sources_) keys.push_back(kv.first);
    auto coord = [](uint64_t k, int shift) { return (int)((int64_t)(k >> shift << 43) >> 43); }; // sign-extend 21 bits
    for (uint64_t k : keys) {
        sources_.erase(k);
        Darken(coord(k, 0), coord(k, 21), coord(k, 42));
    }
}

void LightGrid::SetOpaque(int x, int y, int z, bool opaque) {
    int i;
    Brick& b = BrickOf(x, y, z, i);
    if ((b.opaque[i] != 0) == opaque) return;
    b.opaque[i] = opaque ? 1 : 0;
    if (opaque) {
        // Light through this cell (and diagonal steps past it) goes; the removal pass reaches everything lit that way
        // (a cell lit past it is in its 26-neighbourhood and dimmer than it was).
        Darken(x, y, z);
        if (EmissionAt(x, y, z).level) Relight(x, y, z); // a source keeps shining from inside its block
        return;
    }
    // Opened: every lit neighbour spreads again (that also covers diagonal steps this cell used to block: both of
    // their ends are among its neighbours).
    for (const Move& m : kMoves.m)
        if (At(x + m.dx, y + m.dy, z + m.dz).level) addQ_.push_back({ x + m.dx, y + m.dy, z + m.dz });
    if (EmissionAt(x, y, z).level) Relight(x, y, z);
}

void LightGrid::RefreshOpacity(int minX, int minY, int minZ, int maxX, int maxY, int maxZ) {
    std::vector<uint8_t> fresh(kBrickCells);
    for (int bz = minZ >> 4; bz <= maxZ >> 4; ++bz)
        for (int by = minY >> 4; by <= maxY >> 4; ++by)
            for (int bx = minX >> 4; bx <= maxX >> 4; ++bx) {
                if (!Find(bx, by, bz) || !opacity_) continue;
                opacity_(bx, by, bz, fresh.data());
                for (int z = std::max(minZ, bz * kBrick); z <= std::min(maxZ, bz * kBrick + kBrick - 1); ++z)
                    for (int y = std::max(minY, by * kBrick); y <= std::min(maxY, by * kBrick + kBrick - 1); ++y)
                        for (int x = std::max(minX, bx * kBrick); x <= std::min(maxX, bx * kBrick + kBrick - 1); ++x)
                            SetOpaque(x, y, z, fresh[(x & 15) + (y & 15) * kBrick + (z & 15) * kBrick * kBrick] != 0);
            }
}

void LightGrid::TakeDirty(std::vector<int>& out) {
    out.clear();
    for (Brick* b : dirty_) { b->dirty = false; out.push_back(b->bx); out.push_back(b->by); out.push_back(b->bz); }
    dirty_.clear();
}

bool LightGrid::Update(int budget) {
    int done = 0;
    // Removal first (Minecraft's two-queue relight): a neighbour dimmer than the light that left could have been lit by
    // it, so it goes dark and passes the removal on; a neighbour at least as bright is lit from elsewhere and refills
    // the dark region in the add pass. Source cells that go dark get their own emission back at once.
    while (removeHead_ < removeQ_.size()) {
        if (budget > 0 && done >= budget) return false;
        const RemoveItem it = removeQ_[removeHead_++];
        ++done;
        for (const Move& m : kMoves.m) {
            const int nx = it.x + m.dx, ny = it.y + m.dy, nz = it.z + m.dz;
            Brick* b = Find(nx >> 4, ny >> 4, nz >> 4);
            if (!b) continue;
            LightCell& l = b->light[(nx & 15) + (ny & 15) * kBrick + (nz & 15) * kBrick * kBrick];
            const uint8_t ln = l.level;
            if (ln == 0) continue;
            if (ln < it.old) {
                l = {};
                MarkDirty(*b);
                removeQ_.push_back({ nx, ny, nz, ln });
                const LightCell e = EmissionAt(nx, ny, nz);
                if (e.level) { l = e; addQ_.push_back({ nx, ny, nz }); }
            } else {
                addQ_.push_back({ nx, ny, nz });
            }
        }
    }
    removeQ_.clear();
    removeHead_ = 0;
    while (addHead_ < addQ_.size()) {
        if (budget > 0 && done >= budget) return false;
        const Cell c = addQ_[addHead_++];
        ++done;
        const LightCell l = At(c.x, c.y, c.z);
        for (const Move& m : kMoves.m) {
            if (l.level <= m.cost) continue;
            const int nx = c.x + m.dx, ny = c.y + m.dy, nz = c.z + m.dz;
            if (!CanStep(c.x, c.y, c.z, m.dx, m.dy, m.dz)) continue;
            LightCell arrive = l;
            arrive.level = (uint8_t)(l.level - m.cost);
            int i;
            Brick& tb = BrickOf(nx, ny, nz, i);
            if (!Brighter(arrive, tb.light[i])) continue;
            tb.light[i] = arrive;
            MarkDirty(tb);
            addQ_.push_back({ nx, ny, nz });
        }
    }
    addQ_.clear();
    addHead_ = 0;
    return true;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
#include <random>
#include <unordered_set>

namespace {
// A test world: opaque cells in a set, read by the grid's OpacityFn.
struct SolidSet {
    std::unordered_set<uint64_t> cells;
    bool Has(int x, int y, int z) const { return cells.count(CellKey(x, y, z)) != 0; }
    void Set(int x, int y, int z, bool on) { if (on) cells.insert(CellKey(x, y, z)); else cells.erase(CellKey(x, y, z)); }
    LightGrid::OpacityFn Fn() const {
        return [this](int bx, int by, int bz, uint8_t* out) {
            for (int z = 0; z < 16; ++z)
                for (int y = 0; y < 16; ++y)
                    for (int x = 0; x < 16; ++x) out[x + y * 16 + z * 256] = Has(bx * 16 + x, by * 16 + y, bz * 16 + z) ? 1 : 0;
        };
    }
};
struct Src { int x, y, z; LightCell l; };

// The grid built from scratch for this world and these sources must match the incrementally updated one everywhere.
void ExpectMatchesRebuild(const LightGrid& g, const SolidSet& world, const std::vector<Src>& sources, int r) {
    LightGrid fresh(world.Fn());
    for (const Src& s : sources) fresh.AddSource(s.x, s.y, s.z, s.l);
    fresh.Update();
    int bad = 0;
    for (int z = -r; z <= r; ++z)
        for (int y = -r; y <= r; ++y)
            for (int x = -r; x <= r; ++x) {
                const LightCell a = g.At(x, y, z), b = fresh.At(x, y, z);
                if (a != b && bad++ < 5)
                    ADD_FAILURE() << "cell " << x << "," << y << "," << z << ": incremental " << (int)a.level << " (" << (int)a.r << " "
                                  << (int)a.g << " " << (int)a.b << ") rebuild " << (int)b.level << " (" << (int)b.r << " " << (int)b.g << " " << (int)b.b << ")";
            }
    EXPECT_EQ(bad, 0);
}
}

TEST(LightGrid, OpenSpaceFalloffIsChamferedKeepsItsColourAndEndsAtTheReach) {
    SolidSet world;
    LightGrid g(world.Fn());
    const LightCell torch = LightGrid::Emission(12.0f, 1.0f, 0.5f, 0.0f);
    EXPECT_EQ(torch.level, 96); EXPECT_EQ(torch.r, 255); EXPECT_EQ(torch.g, 128); EXPECT_EQ(torch.b, 0);
    g.AddSource(0, 0, 0, torch);
    EXPECT_TRUE(g.Busy());
    EXPECT_TRUE(g.Update());
    EXPECT_FALSE(g.Busy());
    EXPECT_EQ(g.At(0, 0, 0).level, 96);
    EXPECT_EQ(g.At(1, 0, 0).level, 88);   // face 8
    EXPECT_EQ(g.At(1, 1, 0).level, 85);   // edge 11
    EXPECT_EQ(g.At(1, 1, 1).level, 82);   // corner 14
    EXPECT_EQ(g.At(-3, 0, 0).level, 72);  // the same in every direction, negative cells included
    EXPECT_EQ(g.At(0, 0, -11).level, 8);
    EXPECT_EQ(g.At(0, 0, 12).level, 0);   // the reach: 12 m
    EXPECT_EQ(g.At(0, 0, 12).r, 0);       // unlit cells carry no colour
    EXPECT_EQ(g.At(0, 10, 0).g, 128);     // the colour is the same far out
    // close to a sphere: 8 m along an axis and along the (1,1,0) diagonal light about the same
    EXPECT_NEAR((int)g.At(8, 0, 0).level, (int)g.At(6, 6, 0).level, 8);
    // the lit bricks are reported once for upload: 12 m round (0,0,0) spans bricks -1..0 on each axis
    std::vector<int> dirty;
    g.TakeDirty(dirty);
    EXPECT_EQ(dirty.size(), 8u * 3u);
    g.TakeDirty(dirty);
    EXPECT_TRUE(dirty.empty());
    EXPECT_NE(g.BrickLight(-1, -1, -1), nullptr);
    EXPECT_EQ(g.BrickLight(5, 5, 5), nullptr);
}

TEST(LightGrid, TheBrighterLightWinsAndTiesDoNotDependOnOrder) {
    SolidSet world;
    const LightCell red = LightGrid::Emission(10.0f, 1, 0, 0), blue = LightGrid::Emission(10.0f, 0, 0, 1);
    LightGrid a(world.Fn()), b(world.Fn());
    a.AddSource(-5, 0, 0, red); a.AddSource(5, 0, 0, blue); a.Update();
    b.AddSource(5, 0, 0, blue); b.Update(); b.AddSource(-5, 0, 0, red); b.Update();
    EXPECT_EQ(a.At(-2, 0, 0).r, 255);     // nearer the red light: red
    EXPECT_EQ(a.At(2, 0, 0).b, 255);      // nearer the blue: blue
    for (int x = -10; x <= 10; ++x) EXPECT_EQ(a.At(x, 3, 0), b.At(x, 3, 0)) << x; // the midline tie: the same either way
    EXPECT_EQ(a.At(0, 0, 0).level, 40);   // overlapping lights do not add
}

TEST(LightGrid, AWallStopsLightAndItGoesAroundTheEdge) {
    SolidSet world;
    for (int z = -20; z <= 20; ++z)
        for (int y = -20; y <= 20; ++y) world.Set(2, y, z, true);  // a wall at x = 2...
    for (int z = -20; z <= 20; ++z) world.Set(2, 4, z, false);     // ...with a slot at y = 4
    LightGrid g(world.Fn());
    g.AddSource(0, 0, 0, LightGrid::Emission(20.0f, 1, 1, 1));
    g.Update();
    EXPECT_EQ(g.At(2, 0, 0).level, 0);    // inside the wall
    EXPECT_GT(g.At(3, 4, 0).level, 0);    // through the slot
    // behind the wall, straight across: only light that came round through the slot, so dimmer than in front
    EXPECT_LT(g.At(3, 0, 0).level, g.At(1, 0, 0).level - 4 * LightGrid::kStep);
    // no light leaks through a diagonal wall: cells on x + y = 1 touch only at their edges, and still split (0,0) from (1,1)
    SolidSet edge;
    for (int z = -20; z <= 20; ++z)
        for (int i = -20; i <= 20; ++i) edge.Set(i, 1 - i, z, true);
    LightGrid e(edge.Fn());
    e.AddSource(0, 0, 0, LightGrid::Emission(10.0f, 1, 1, 1));
    e.Update();
    EXPECT_EQ(e.At(1, 1, 0).level, 0);
}

TEST(LightGrid, RemovingASourceDarksItsLightAndKeepsTheOther) {
    SolidSet world;
    LightGrid g(world.Fn());
    const LightCell a = LightGrid::Emission(12.0f, 1, 0.6f, 0.3f), b = LightGrid::Emission(20.0f, 0.4f, 0.8f, 1);
    g.AddSource(0, 0, 0, a);
    g.AddSource(9, 2, 0, b);
    g.Update();
    g.RemoveSource(0, 0, 0, a);
    g.Update();
    ExpectMatchesRebuild(g, world, { { 9, 2, 0, b } }, 32);
    g.RemoveSource(9, 2, 0, b);
    g.Update();
    ExpectMatchesRebuild(g, world, {}, 32);
    // two sources in one cell: removing one leaves the other
    g.AddSource(1, 1, 1, a);
    g.AddSource(1, 1, 1, b);
    g.Update();
    g.RemoveSource(1, 1, 1, b);
    g.Update();
    ExpectMatchesRebuild(g, world, { { 1, 1, 1, a } }, 32);
}

TEST(LightGrid, DigAndFillMatchARebuild) {
    // A random cave-ish world, random lights, then random digs and fills: after each batch the incremental grid equals
    // one built from scratch (the fixed point is unique, so any correct relight lands on it).
    std::mt19937 rng(1739);
    SolidSet world;
    std::uniform_int_distribution<int> c(-14, 14);
    for (int i = 0; i < 6000; ++i) world.Set(c(rng), c(rng), c(rng), true);
    std::vector<Src> sources;
    LightGrid g(world.Fn());
    const float radii[] = { 6, 12, 20, 30 };
    for (int i = 0; i < 6; ++i) {
        Src s{ c(rng), c(rng), c(rng), LightGrid::Emission(radii[i % 4], (i % 3 == 0) ? 1.0f : 0.5f, 0.8f, (i % 2) ? 1.0f : 0.2f) };
        sources.push_back(s);
        g.AddSource(s.x, s.y, s.z, s.l);
    }
    g.Update();
    ExpectMatchesRebuild(g, world, sources, 48);
    for (int round = 0; round < 6; ++round) {
        for (int k = 0; k < 40; ++k) {
            const int x = c(rng), y = c(rng), z = c(rng);
            const bool fill = (rng() & 1) != 0;
            world.Set(x, y, z, fill);
            g.RefreshOpacity(x, y, z, x, y, z);
        }
        if (round == 3) { g.RemoveSource(sources[1].x, sources[1].y, sources[1].z, sources[1].l); sources.erase(sources.begin() + 1); }
        // spread over many small updates, as the game does per frame
        int frames = 0;
        while (!g.Update(500)) ++frames;
        EXPECT_LT(frames, 2000);
        ExpectMatchesRebuild(g, world, sources, 48);
    }
}

TEST(LightGrid, ASourceInsideASolidCellStillShinesOut) {
    SolidSet world;
    world.Set(0, 0, 0, true); // a lamp block, or a torch the opacity test calls solid
    LightGrid g(world.Fn());
    g.AddSource(0, 0, 0, LightGrid::Emission(6.0f, 1, 1, 1));
    g.Update();
    EXPECT_EQ(g.At(0, 0, 0).level, 48);
    EXPECT_EQ(g.At(0, 2, 0).level, 32);
    g.ClearSources();
    g.Update();
    ExpectMatchesRebuild(g, world, {}, 16);
}
#endif
