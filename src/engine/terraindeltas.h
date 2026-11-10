#ifndef TERRAINDELTAS_H
#define TERRAINDELTAS_H

// Player terrain edits as density deltas (milestone 2, design doc section 3).
//
// Every brush stroke is folded straight into a lattice of density changes at 1-voxel spacing,
// stored in 32^3 chunks of 16-bit values (1/256 density steps). Density = base field + the
// trilinear delta, so sampling costs the same whether the world has ten strokes or a million (the
// old brush list was scanned on every sample). A TerrainDeltas value is an immutable snapshot
// shared with the worker threads; an edit builds a new snapshot that shares every untouched chunk.

#include <raylib.h>
#include <cstdint>
#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class DensityField;
struct TerrainEdit;

constexpr int kDeltaChunk = 32;                 // lattice samples per chunk edge
constexpr int kRegionChunks = 16;               // chunks per region file edge
constexpr float kDeltaStep = 1.0f / 256.0f;     // density per stored unit

struct DeltaChunk {
    int16_t d[kDeltaChunk * kDeltaChunk * kDeltaChunk] = {};
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
