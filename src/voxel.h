#ifndef VOXEL_H
#define VOXEL_H

#include <raylib.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <vector>

#include "SimplexNoise.h"

extern const int CHUNK_SIZE;          // cells per chunk edge (samples per edge = CHUNK_SIZE + 1)
constexpr float ISO_LEVEL = 0.5f;     // density >= ISO_LEVEL is solid

struct Int3 {
    int x, y, z;
    bool operator==(const Int3& o) const { return x == o.x && y == o.y && z == o.z; }
};

namespace std {
    template<> struct hash<Int3> {
        size_t operator()(const Int3& k) const {
            // Large primes; much better spread than xor-shifting std::hash<int> (identity on MSVC/GCC).
            return (size_t)((int64_t)k.x * 73856093LL ^ (int64_t)k.y * 19349663LL ^ (int64_t)k.z * 83492791LL);
        }
    };
}

// Voxel materials. Meshes carry them as per-vertex blend weights in the vertex colour
// (R = rock, G = regolith, B = ice, A = ore), so they interpolate smoothly across triangles.
enum class VoxelMaterial : uint8_t { Rock = 0, Regolith = 1, Ice = 2, Ore = 3, Count };
const char* VoxelMaterialName(VoxelMaterial m);

// A spherical density brush in planetoid-local space. strength < 0 digs, > 0 adds.
struct TerrainEdit {
    Vector3 center;
    float radius;
    float strength;
};
using EditList = std::vector<TerrainEdit>;

// Immutable description of one planetoid's terrain. Shared read-only with worker threads.
struct PlanetoidField {
    SimplexNoise shape;    // overall form
    SimplexNoise detail;   // ore veins / ice patches
    float radius = 50.0f;  // density fades to 0 at this distance from the centre
    float iceAmount = 0.0f;// 0..1: how much of the body is ice

    PlanetoidField(unsigned int seed, float radius);

    // Base density (no edits) at planetoid-local position p.
    float baseDensity(Vector3 p) const;
    // Density including edits (the list should already be filtered to ones near p).
    float density(Vector3 p, const TerrainEdit* edits, size_t editCount) const;
    // Outward surface normal: negative density gradient, by central differences.
    Vector3 normal(Vector3 p, const TerrainEdit* edits, size_t editCount) const;
    VoxelMaterial material(Vector3 p, float densityAtP) const;
    // Densities above this are impossible without edits; used to skip far chunks.
    float maxSurfaceRadius() const { return radius; }

private:
    PlanetoidField(std::mt19937 rng, float radius);
};

// CPU-side result of meshing one chunk (produced on a worker thread).
struct ChunkMeshData {
    uint64_t planetoidId = 0;
    Int3 coord{};
    uint32_t version = 0;
    std::vector<float> positions;          // xyz, chunk-local
    std::vector<float> normals;            // xyz
    std::vector<unsigned char> colors;     // rgba material weights
    std::vector<unsigned short> indices;
    bool overflow = false;                 // > 65535 vertices: mesh dropped
    float buildMs = 0.0f;
};

// Samples the field for chunk `coord` (offset coord * CHUNK_SIZE) and runs marching cubes.
// Thread-safe: touches only its arguments.
void BuildChunkMesh(const PlanetoidField& field, const EditList& allEdits, ChunkMeshData& out);

// Edits whose influence (radius + margin) overlaps the axis-aligned box [lo, hi].
void FilterEdits(const EditList& all, Vector3 lo, Vector3 hi, float margin, EditList& out);

#endif
