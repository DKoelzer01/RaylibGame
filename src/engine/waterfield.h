#ifndef WATERFIELD_H
#define WATERFIELD_H

// One voxel water field (#360, Plan 8 item 3).
//
// Water anywhere = generated water (the sea, rivers, lakes and ponds: open cells under their basin's level,
// PlanetField::waterSourceAt, computed on demand and never stored) overridden by the simulation's regions
// (fluid.h): a published brick of 16^3 one-metre cells is the truth for its cells - water, no water (a drained
// pond) or solid - and everything outside the bricks is generated.
//
// The simulation publishes an immutable Snapshot every tick it changes something; chunk builds (worker
// threads), queries and the water-only chunk rebuilds read it without locks. The chunk mesher draws every
// fluid surface from it at every level of detail (MeshChunkWater, voxel.cpp): there are no per-brick meshes.

#include <raylib.h>
#include <array>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

namespace WaterField {

constexpr int kBrick = 16, kBrick3 = kBrick * kBrick * kBrick;
static_assert(kBrick == 16, "ChunkSampler::brickAt / cellAt index with >> 4 and & 15");
// Cell flags (fluid.cpp): solid; a source (infinite generated water, held at its level); inside generated
// infinite water; open only because of an edit; a sink (a source that also takes water in and drops what is
// over its level: the downstream end of a released river reach, #361); the fluid kind in bits 5-7.
constexpr uint8_t kSolid = 1, kSource = 2, kStatic = 4, kDug = 8, kSink = 16;
constexpr int kKindShift = 5;
constexpr uint8_t kKindMask = 0xE0;
inline int CellKind(uint8_t fl) { return fl >> kKindShift; }
inline int CellIdx(int x, int y, int z) { return x + y * kBrick + z * kBrick * kBrick; }
inline int64_t BrickKey(int i, int j, int k) { return ((int64_t)(i + 1024) << 22) | ((int64_t)(j + 1024) << 11) | (int64_t)(k + 1024); }
inline int FloorDiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }

struct Brick {
    std::array<uint8_t, kBrick3> fill{}, flags{};
    std::array<uint8_t, kBrick3> shut{}; // the share of an open cell that is ground, at its bottom (#454): holds 255 - shut, water stands at shut + fill
    uint8_t kinds = 0; // bit k: some cell holds fluid of kind k
};

// A simulated region: a tangent frame (east, north, up) of 1 m cells anchored at `origin`.
struct Site {
    uint64_t body = 0;
    Vector3 origin{}, e{}, n{}, u{};
    Vector3 lo{ 1e30f, 1e30f, 1e30f }, hi{ -1e30f, -1e30f, -1e30f }; // body-local bounds of its bricks
    std::unordered_map<int64_t, std::shared_ptr<const Brick>> bricks;

    // Site coordinates of a body-local point (and its height in cells, continuous).
    void frame(Vector3 local, float& fx, float& fy, float& fz) const {
        const Vector3 d = { local.x - origin.x, local.y - origin.y, local.z - origin.z };
        fx = d.x * e.x + d.y * e.y + d.z * e.z;
        fy = d.x * n.x + d.y * n.y + d.z * n.z;
        fz = d.x * u.x + d.y * u.y + d.z * u.z;
    }
    const Brick* brick(int bi, int bj, int bk) const {
        auto it = bricks.find(BrickKey(bi, bj, bk));
        return it == bricks.end() ? nullptr : it->second.get();
    }
    // A cell by site coordinates: false if no brick holds it.
    bool cell(int ci, int cj, int ck, uint8_t& fill, uint8_t& flags) const {
        const int bi = FloorDiv(ci, kBrick), bj = FloorDiv(cj, kBrick), bk = FloorDiv(ck, kBrick);
        const Brick* b = brick(bi, bj, bk);
        if (!b) return false;
        const int i = CellIdx(ci - bi * kBrick, cj - bj * kBrick, ck - bk * kBrick);
        fill = b->fill[i]; flags = b->flags[i];
        return true;
    }
    // The same with the cell's ground share (#454).
    bool cell(int ci, int cj, int ck, uint8_t& fill, uint8_t& flags, uint8_t& shut) const {
        const int bi = FloorDiv(ci, kBrick), bj = FloorDiv(cj, kBrick), bk = FloorDiv(ck, kBrick);
        const Brick* b = brick(bi, bj, bk);
        if (!b) return false;
        const int i = CellIdx(ci - bi * kBrick, cj - bj * kBrick, ck - bk * kBrick);
        fill = b->fill[i]; flags = b->flags[i]; shut = b->shut[i];
        return true;
    }
    // Body-local corners of a brick (its bounds), for dirty regions and the site's bounds.
    void brickBounds(int bi, int bj, int bk, Vector3& blo, Vector3& bhi) const;
    void cellBounds(int i0, int j0, int k0, int i1, int j1, int k1, Vector3& blo, Vector3& bhi) const; // cells [i0, i1) x ...
};

struct Snapshot {
    uint64_t seq = 0;
    std::vector<Site> sites;
};
using SnapshotPtr = std::shared_ptr<const Snapshot>;

// The latest published snapshot (never null; empty before the first publish).
SnapshotPtr Current();
// The simulation publishes a new snapshot with the body-local boxes whose water changed (the chunks there redraw).
struct Dirty { uint64_t body; Vector3 lo, hi; int mag = 255; }; // mag: the largest fill change in it (255: new, solidity, kind)
void Publish(SnapshotPtr snap, std::vector<Dirty>&& dirty);
void Clear();                         // a new world / a load
void TakeDirty(std::vector<Dirty>& out); // main thread: boxes since the last call

// A point's state in the published field: -1 not in a simulated cell (ask the generated water), else the
// fill (0-255) of fluid kind `kind` (other kinds and solid count as 0). `frac`: the point's height in its cell.
int CellFill(const Snapshot& s, uint64_t body, Vector3 local, int kind, float* frac = nullptr, uint8_t* flagsOut = nullptr);

// Chunk meshing: the fluid surface of the published cells around one chunk as a density (ISO-based, the same
// gradient as the generated water so both meet as one surface).
class ChunkSampler {
public:
    // probeOnly: just whether any() (no dense caches, stops at the first site; for the main thread).
    ChunkSampler(const Snapshot* snap, uint64_t body, Vector3 lo, Vector3 hi, bool probeOnly = false);
    bool any() const { return !sites.empty(); }  // some published brick reaches the box
    uint8_t kinds() const { return kindMask; }    // fluid kinds present (bit k)
    // Water density (ISO + gradient * (surface - height)) of kind `kind` at a body-local point; false where no
    // published cell covers it (the generated water decides there).
    bool density(Vector3 p, int kind, float& wd);
    bool covered(Vector3 p) const;
    // Metres of fluid in the column under a covered point (up to 8), for the shader's depth tint.
    float depthBelow(Vector3 p, int kind);
private:
    std::vector<const Site*> sites;
    uint8_t kindMask = 0;
    std::unordered_map<int64_t, float> memo; // (site, cell, kind) -> surface height
    // Per site, over the chunk's box in site cells (+2 margin): water surfaces (NaN: not yet asked) and which bricks are
    // published (#423: the blended density asks four columns a sample; hashing every one doubled a water rebuild).
    struct Dense { int c0[3] = {}, n[3] = {}, b0[3] = {}, bn[3] = {}; std::vector<float> surf; std::vector<const Brick*> brick; };
    std::vector<Dense> dense;
    bool blend = true; // blend columns (#423): only in fine chunks (samples finer than a cell); coarse ones gain nothing (#381)
    const Brick* brickPtr(int si, int bi, int bj, int bk) const; // by brick coordinates, through the dense cache
    bool brickAt(int si, int ci, int cj, int ck) const { return brickPtr(si, ci >> 4, cj >> 4, ck >> 4) != nullptr; }
    bool cellAt(int si, int ci, int cj, int ck, uint8_t& fill, uint8_t& flags) const; // Site::cell through the cache
    bool cellAt(int si, int ci, int cj, int ck, uint8_t& fill, uint8_t& flags, uint8_t& shut) const; // with the ground share (#454)
    float surface(int si, int ci, int cj, int ck, int kind);
    float openSurface(int si, int ci, int cj, int ck, int kind, int depth);
};

} // namespace WaterField

#endif
