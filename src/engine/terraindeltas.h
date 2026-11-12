#ifndef TERRAINDELTAS_H
#define TERRAINDELTAS_H

// Player terrain edits as density deltas (milestone 2, design doc section 3).
//
// Every brush stroke is folded straight into a lattice of density changes at 1-voxel spacing,
// stored in 32^3 chunks of 16-bit values (1/64 density steps). Density = base field + the
// trilinear delta, so sampling costs the same whether the world has ten strokes or a million (the
// old brush list was scanned on every sample). A TerrainDeltas value is an immutable snapshot
// shared with the worker threads; an edit builds a new snapshot that shares every untouched chunk.

#include <raylib.h>
#include <cstdint>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class DensityField;
struct TerrainEdit;

constexpr int kDeltaChunk = 32;                 // lattice samples per chunk edge
constexpr int kRegionChunks = 16;               // chunks per region file edge
// Density per stored unit: +-512 in int16, about 2 km of rock (0.25 a metre) - enough to dig from a
// mountain top to the mantle (1/256 saturated at ~510 m: an invisible floor). ~6 cm of surface precision.
constexpr float kDeltaStep = 1.0f / 64.0f;

struct DeltaChunk {
    static constexpr int kN = kDeltaChunk * kDeltaChunk * kDeltaChunk;
    int16_t d[kN] = {};
    // Placed ground's material per lattice point: 0 = the field's own, else VoxelMaterial + 1
    // (gravel placed stays gravel instead of taking on the dirt or grass the ground has there).
    // Allocated on first use: most chunks are only dug (the 0.25 m lattice made chunks 64x more).
    std::unique_ptr<uint8_t[]> mats;
    DeltaChunk() = default;
    DeltaChunk(const DeltaChunk& o) {
        memcpy(d, o.d, sizeof d);
        if (o.mats) { mats.reset(new uint8_t[kN]); memcpy(mats.get(), o.mats.get(), kN); }
    }
    DeltaChunk& operator=(const DeltaChunk&) = delete;
    uint8_t mat(int i) const { return mats ? mats[i] : 0; }
    void setMat(int i, uint8_t v) { if (!mats) { if (!v) return; mats.reset(new uint8_t[kN]()); } mats[i] = v; }
};

class TerrainDeltas {
public:
    using Key = uint64_t;                       // packed chunk coordinates
    std::unordered_map<Key, std::shared_ptr<const DeltaChunk>> chunks;
    size_t strokes = 0;                         // strokes folded in (stats, saves)

    static Key KeyOf(int cx, int cy, int cz);
    static void Unpack(Key k, int& cx, int& cy, int& cz);
    static Key RegionOf(Key chunk);             // packed region coordinates of a chunk

    bool empty() const { return chunks.empty(); }
    float at(int x, int y, int z) const;        // lattice delta (0 where never edited)
    float sample(Vector3 p) const;              // trilinear between lattice points
    // The placed material at a point: the nearest lattice point's, else any of the cell's corners
    // (a surface point lies between a placed solid point and air); -1 where nothing was placed.
    int materialAt(Vector3 p) const;
    // Chunks overlapping a body-local box (a build's quick test; empty -> no edits there).
    void chunksIn(Vector3 lo, Vector3 hi, std::vector<std::pair<Key, const DeltaChunk*>>& out) const;
};

// Builds the next snapshot: copies the chunk map once, and each chunk only the first time a stroke
// touches it. `touched` collects the chunk keys changed (for region saves).
class TerrainDeltasBuilder {
public:
    explicit TerrainDeltasBuilder(const TerrainDeltas& from) : next(std::make_shared<TerrainDeltas>(from)) {}
    // Fold a brush stroke in (smoothstep falloff, as the brush always had); samples the field
    // refuses (a planet's mantle) are skipped.
    void apply(const TerrainEdit& e, const DensityField& field, std::unordered_set<TerrainDeltas::Key>* touched = nullptr);
    void setChunk(TerrainDeltas::Key k, std::shared_ptr<const DeltaChunk> c) { next->chunks[k] = std::move(c); owned.erase(k); }
    std::shared_ptr<const TerrainDeltas> finish() { return next; }
    TerrainDeltas& state() { return *next; }

private:
    std::shared_ptr<TerrainDeltas> next;
    std::unordered_map<TerrainDeltas::Key, DeltaChunk*> owned; // chunks already copied in this build
};

#endif
