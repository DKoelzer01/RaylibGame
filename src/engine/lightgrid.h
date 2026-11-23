#ifndef LIGHTGRID_H
#define LIGHTGRID_H

// Static light-level grid (#1739): torches, lamps and other placed lights flood a 1 m grid through open cells and stop at
// solid ones, like Minecraft's block light, so any number of them can light a cave (the point-light path keeps only the
// nearest kMaxLocalLights). Held and moving lights stay point lights.
//
// Cells are 1 m cubes in the body-local frame: cell (x, y, z) covers [x, x+1) etc. Each holds a level (0-255) and the
// colour of the light that brought it (r, g, b scaled so the largest is 255). A source emits `Emission(radius, colour)`:
// level = 8 per metre of reach. Light loses kStep per metre of travel (chamfer steps: face 8, edge 11, corner 14, so the lit
// shape is close to a sphere, not Minecraft's diamond) and keeps its colour, so a torch is the same hue near and far. A
// cell keeps the brightest arrival (ties broken by colour, so the result does not depend on the order of the work);
// overlapping lights do not add, as in Minecraft. 0-255 instead of the 0-15 written on the issue: 1 m cells and a 30 m
// reach need finer steps.
//
// Storage is sparse: 16^3 bricks created when light reaches them; each brick also caches which of its cells are opaque
// (filled once from the caller's OpacityFn, e.g. terrain density at the cell centres plus solid parts). Changes are
// queued and applied by Update(budget) a bounded number of cells at a time (Minecraft-style removal then refill), so
// digging through a lit cave spreads its relight over frames. Nothing is saved: the game rebuilds it from the terrain
// and the parts, and in multiplayer each machine builds the same grid from the same data.
//
// Not thread-safe: one owner (the main thread) calls everything.

#include <cstdint>
#include <cstddef>
#include <functional>
#include <unordered_map>
#include <vector>

struct LightCell {
    uint8_t r = 0, g = 0, b = 0;    // colour, the largest channel 255 (0 0 0 where unlit)
    uint8_t level = 0;              // 8 per metre of reach left
    bool operator==(const LightCell& o) const { return r == o.r && g == o.g && b == o.b && level == o.level; }
    bool operator!=(const LightCell& o) const { return !(*this == o); }
};

class LightGrid {
public:
    static constexpr int kBrick = 16;                   // cells per brick edge
    static constexpr int kBrickCells = kBrick * kBrick * kBrick;
    static constexpr int kStep = 8;                     // level lost per metre (face step); edge 11, corner 14

    // Fills out[kBrickCells] (x fastest, then y, then z) for the brick whose first cell is (bx, by, bz) * kBrick:
    // 1 = opaque (light stops), 0 = open.
    using OpacityFn = std::function<void(int bx, int by, int bz, uint8_t* out)>;

    explicit LightGrid(OpacityFn opacity);

    // A source of this reach (m, capped at 31) and linear colour (0..1 per channel; only the ratios count).
    static LightCell Emission(float radiusM, float r, float g, float b);

    // Sources: several may share a cell (the brightest wins). Remove takes the same value that was added.
    void AddSource(int x, int y, int z, LightCell light);
    void RemoveSource(int x, int y, int z, LightCell light);
    void ClearSources();                                // removes every source (queued like any removal)

    // The caller changed what is solid inside [min, max] (inclusive cells): re-read the opacity of the bricks it
    // touches and relight the cells that changed. Bricks not created yet are read when light first reaches them.
    void RefreshOpacity(int minX, int minY, int minZ, int maxX, int maxY, int maxZ);
    // One cell's opacity set directly (tests, and parts that know their own cells).
    void SetOpaque(int x, int y, int z, bool opaque);

    // Applies queued work, at most `budget` cells (<= 0: everything). Returns true when nothing is left queued.
    bool Update(int budget = 0);
    bool Busy() const { return removeHead_ < removeQ_.size() || addHead_ < addQ_.size(); }

    LightCell At(int x, int y, int z) const;            // unlit where no brick exists
    bool OpaqueAt(int x, int y, int z) const;           // false where no brick exists (not read yet)
    size_t BrickCount() const { return bricks_.size(); }
    size_t MemoryBytes() const { return bricks_.size() * sizeof(Brick); }

    // One brick's cells (kBrickCells, x fastest), or null where none exists.
    const LightCell* BrickLight(int bx, int by, int bz) const { const Brick* b = Find(bx, by, bz); return b ? b->light : nullptr; }
    // Bricks whose light changed since the last call (brick coordinates, 3 ints each); clears the list.
    void TakeDirty(std::vector<int>& out);

private:
    struct Brick {
        int bx, by, bz;
        bool dirty = false;                              // light changed since TakeDirty
        LightCell light[kBrickCells];
        uint8_t opaque[kBrickCells];
    };
    struct RemoveItem { int x, y, z; uint8_t old; };
    struct Cell { int x, y, z; };

    static uint64_t Key(int bx, int by, int bz);
    Brick* Find(int bx, int by, int bz) const;
    Brick& Get(int bx, int by, int bz);                  // creates (and reads its opacity) when missing
    Brick& BrickOf(int x, int y, int z, int& index);
    bool Opaque(int x, int y, int z);                    // creates the brick to know
    bool CanStep(int x, int y, int z, int dx, int dy, int dz); // the target and the cells the step passes are open
    LightCell EmissionAt(int x, int y, int z) const;
    void Darken(int x, int y, int z);                    // queue removal of a cell's light
    void Relight(int x, int y, int z);                   // the cell's emission back, and queue it to spread
    void MarkDirty(Brick& b) { if (!b.dirty) { b.dirty = true; dirty_.push_back(&b); } }

    OpacityFn opacity_;
    std::unordered_map<uint64_t, Brick> bricks_;
    mutable const Brick* last_ = nullptr;                // lookup cache (pointers into the node-based map are stable)
    std::unordered_map<uint64_t, std::vector<LightCell>> sources_; // by cell key
    std::vector<RemoveItem> removeQ_;
    size_t removeHead_ = 0;
    std::vector<Cell> addQ_;
    size_t addHead_ = 0;
    std::vector<Brick*> dirty_;
};

#endif
