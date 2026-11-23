#include "lightgrid.h"

#include <algorithm>
#include <cmath>
#include <cstring>

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
uint8_t Channel(const LightRGB& c, int ch) { return ch == 0 ? c.r : ch == 1 ? c.g : c.b; }
uint8_t& Channel(LightRGB& c, int ch) { return ch == 0 ? c.r : ch == 1 ? c.g : c.b; }
}

LightGrid::LightGrid(OpacityFn opacity) : opacity_(std::move(opacity)) {}

LightRGB LightGrid::Emission(float radiusM, float r, float g, float b) {
    const float reach = std::min(std::max(radiusM, 0.0f), 31.0f) * kStep; // 31 m * 8 = 248 <= 255
    auto lvl = [&](float c) { return (uint8_t)std::lround(std::min(std::max(c, 0.0f), 1.0f) * reach); };
    return { lvl(r), lvl(g), lvl(b) };
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
    std::memset(b.opaque, 0, sizeof(b.opaque)); // b.light starts at 0 (LightRGB's initializers)
    if (opacity_) opacity_(bx, by, bz, b.opaque);
    last_ = &b;
    return b;
}

LightGrid::Brick& LightGrid::BrickOf(int x, int y, int z, int& index) {
    index = (x & 15) + (y & 15) * kBrick + (z & 15) * kBrick * kBrick;
    return Get(x >> 4, y >> 4, z >> 4);
}

LightRGB LightGrid::At(int x, int y, int z) const {
    const Brick* b = Find(x >> 4, y >> 4, z >> 4);
    return b ? b->light[(x & 15) + (y & 15) * kBrick + (z & 15) * kBrick * kBrick] : LightRGB{};
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

LightRGB LightGrid::EmissionAt(int x, int y, int z) const {
    auto it = sources_.find(CellKey(x, y, z));
    LightRGB e;
    if (it == sources_.end()) return e;
    for (const LightRGB& s : it->second) { e.r = std::max(e.r, s.r); e.g = std::max(e.g, s.g); e.b = std::max(e.b, s.b); }
    return e;
}

void LightGrid::Darken(int x, int y, int z) {
    int i;
    Brick& b = BrickOf(x, y, z, i);
    for (int ch = 0; ch < 3; ++ch) {
        const uint8_t old = Channel(b.light[i], ch);
        if (old == 0) continue;
        Channel(b.light[i], ch) = 0;
        removeQ_.push_back({ x, y, z, (uint8_t)ch, old });
    }
}

void LightGrid::Relight(int x, int y, int z) {
    const LightRGB e = EmissionAt(x, y, z);
    int i;
    Brick& b = BrickOf(x, y, z, i);
    LightRGB& l = b.light[i];
    l.r = std::max(l.r, e.r); l.g = std::max(l.g, e.g); l.b = std::max(l.b, e.b);
    addQ_.push_back({ x, y, z });
}

void LightGrid::AddSource(int x, int y, int z, LightRGB level) {
    sources_[CellKey(x, y, z)].push_back(level);
    Relight(x, y, z);
}

void LightGrid::RemoveSource(int x, int y, int z, LightRGB level) {
    auto it = sources_.find(CellKey(x, y, z));
    if (it == sources_.end()) return;
    auto& v = it->second;
    auto f = std::find(v.begin(), v.end(), level);
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
        const LightRGB e = EmissionAt(x, y, z);
        if (e.r || e.g || e.b) Relight(x, y, z); // a source keeps shining from inside its block
        return;
    }
    // Opened: every lit neighbour spreads again (that also covers diagonal steps this cell used to block: both of
    // their ends are among its neighbours).
    for (const Move& m : kMoves.m) {
        const LightRGB l = At(x + m.dx, y + m.dy, z + m.dz);
        if (l.r || l.g || l.b) addQ_.push_back({ x + m.dx, y + m.dy, z + m.dz });
    }
    const LightRGB e = EmissionAt(x, y, z);
    if (e.r || e.g || e.b) Relight(x, y, z);
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
            LightRGB& l = b->light[(nx & 15) + (ny & 15) * kBrick + (nz & 15) * kBrick * kBrick];
            const uint8_t ln = Channel(l, it.channel);
            if (ln == 0) continue;
            if (ln < it.old) {
                Channel(l, it.channel) = 0;
                removeQ_.push_back({ nx, ny, nz, it.channel, ln });
                const uint8_t e = Channel(EmissionAt(nx, ny, nz), it.channel);
                if (e) { Channel(l, it.channel) = e; addQ_.push_back({ nx, ny, nz }); }
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
        const LightRGB l = At(c.x, c.y, c.z);
        const int top = std::max(l.r, std::max(l.g, l.b));
        for (const Move& m : kMoves.m) {
            if (top <= m.cost) continue;
            const int nx = c.x + m.dx, ny = c.y + m.dy, nz = c.z + m.dz;
            if (!CanStep(c.x, c.y, c.z, m.dx, m.dy, m.dz)) continue;
            int i;
            LightRGB& t = BrickOf(nx, ny, nz, i).light[i];
            bool raised = false;
            for (int ch = 0; ch < 3; ++ch) {
                const int nv = (int)Channel(l, ch) - m.cost;
                if (nv > (int)Channel(t, ch)) { Channel(t, ch) = (uint8_t)nv; raised = true; }
            }
            if (raised) addQ_.push_back({ nx, ny, nz });
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
struct Src { int x, y, z; LightRGB l; };

// The grid built from scratch for this world and these sources must match the incrementally updated one everywhere.
void ExpectMatchesRebuild(const LightGrid& g, const SolidSet& world, const std::vector<Src>& sources, int r) {
    LightGrid fresh(world.Fn());
    for (const Src& s : sources) fresh.AddSource(s.x, s.y, s.z, s.l);
    fresh.Update();
    int bad = 0;
    for (int z = -r; z <= r; ++z)
        for (int y = -r; y <= r; ++y)
            for (int x = -r; x <= r; ++x)
                if (g.At(x, y, z) != fresh.At(x, y, z) && bad++ < 5)
                    ADD_FAILURE() << "cell " << x << "," << y << "," << z << ": incremental " << (int)g.At(x, y, z).r << "/"
                                  << (int)g.At(x, y, z).g << "/" << (int)g.At(x, y, z).b << " rebuild " << (int)fresh.At(x, y, z).r
                                  << "/" << (int)fresh.At(x, y, z).g << "/" << (int)fresh.At(x, y, z).b;
    EXPECT_EQ(bad, 0);
}
}

TEST(LightGrid, OpenSpaceFalloffIsChamferedAndEndsAtTheReach) {
    SolidSet world;
    LightGrid g(world.Fn());
    const LightRGB torch = LightGrid::Emission(12.0f, 1.0f, 0.5f, 0.0f);
    EXPECT_EQ(torch.r, 96); EXPECT_EQ(torch.g, 48); EXPECT_EQ(torch.b, 0);
    g.AddSource(0, 0, 0, torch);
    EXPECT_TRUE(g.Busy());
    EXPECT_TRUE(g.Update());
    EXPECT_FALSE(g.Busy());
    EXPECT_EQ(g.At(0, 0, 0).r, 96);
    EXPECT_EQ(g.At(1, 0, 0).r, 88);   // face 8
    EXPECT_EQ(g.At(1, 1, 0).r, 85);   // edge 11
    EXPECT_EQ(g.At(1, 1, 1).r, 82);   // corner 14
    EXPECT_EQ(g.At(-3, 0, 0).r, 72);  // the same in every direction, negative cells included
    EXPECT_EQ(g.At(0, 0, -11).r, 8);
    EXPECT_EQ(g.At(0, 0, 12).r, 0);   // the reach: 12 m
    EXPECT_EQ(g.At(6, 0, 0).g, 0);    // green's reach is half
    EXPECT_EQ(g.At(5, 0, 0).g, 8);
    EXPECT_EQ(g.At(2, 0, 0).b, 0);
    // close to a sphere: 8 m along an axis and along the (1,1,0) diagonal light about the same
    EXPECT_NEAR((int)g.At(8, 0, 0).r, (int)g.At(6, 6, 0).r, 8);
}

TEST(LightGrid, AWallStopsLightAndItGoesAroundTheEdge) {
    SolidSet world;
    for (int z = -20; z <= 20; ++z)
        for (int y = -20; y <= 20; ++y) world.Set(2, y, z, true);  // a wall at x = 2...
    for (int z = -20; z <= 20; ++z) world.Set(2, 4, z, false);     // ...with a slot at y = 4
    LightGrid g(world.Fn());
    g.AddSource(0, 0, 0, LightGrid::Emission(20.0f, 1, 1, 1));
    g.Update();
    EXPECT_EQ(g.At(2, 0, 0).r, 0);    // inside the wall
    EXPECT_GT(g.At(3, 4, 0).r, 0);    // through the slot
    // behind the wall, straight across: only light that came round through the slot, so dimmer than in front
    EXPECT_LT(g.At(3, 0, 0).r, g.At(1, 0, 0).r - 4 * LightGrid::kStep);
    // no light leaks through a diagonal wall: cells on x + y = 1 touch only at their edges, and still split (0,0) from (1,1)
    SolidSet edge;
    for (int z = -20; z <= 20; ++z)
        for (int i = -20; i <= 20; ++i) edge.Set(i, 1 - i, z, true);
    LightGrid e(edge.Fn());
    e.AddSource(0, 0, 0, LightGrid::Emission(10.0f, 1, 1, 1));
    e.Update();
    EXPECT_EQ(e.At(1, 1, 0).r, 0);
}

TEST(LightGrid, RemovingASourceDarksItsLightAndKeepsTheOther) {
    SolidSet world;
    LightGrid g(world.Fn());
    const LightRGB a = LightGrid::Emission(12.0f, 1, 0.6f, 0.3f), b = LightGrid::Emission(20.0f, 0.4f, 0.8f, 1);
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
    EXPECT_EQ(g.At(0, 0, 0).r, 48);
    EXPECT_EQ(g.At(0, 2, 0).r, 32);
    g.ClearSources();
    g.Update();
    ExpectMatchesRebuild(g, world, {}, 16);
}
#endif
