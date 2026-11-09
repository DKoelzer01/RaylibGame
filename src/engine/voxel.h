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
constexpr float DENSITY_GRADIENT = 0.25f; // density change per voxel across the surface

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


// Octree LOD node. A node at `level` covers CHUNK_SIZE * 2^level voxels per axis, meshed with
// CHUNK_SIZE cells (stride 2^level), so every node costs about the same to build and draw.
// (x, y, z) are in units of the node size at that level.
struct NodeKey {
    int level, x, y, z;
    bool operator==(const NodeKey& o) const { return level == o.level && x == o.x && y == o.y && z == o.z; }
    int stride() const { return 1 << level; }
    float size() const { return (float)(CHUNK_SIZE << level); }
    Vector3 origin() const { const float s = size(); return { x * s, y * s, z * s }; }
    // True if this node's box contains `o`'s box (o is the same node or a descendant).
    bool contains(const NodeKey& o) const {
        if (o.level > level) return false;
        const int shift = level - o.level;
        return (o.x >> shift) == x && (o.y >> shift) == y && (o.z >> shift) == z;
    }
    bool overlaps(const NodeKey& o) const { return contains(o) || o.contains(*this); }
};

namespace std {
    template<> struct hash<NodeKey> {
        size_t operator()(const NodeKey& k) const {
            return (size_t)((int64_t)k.x * 73856093LL ^ (int64_t)k.y * 19349663LL ^ (int64_t)k.z * 83492791LL ^ (int64_t)k.level * 2654435761LL);
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

// A terrain density field: density >= ISO_LEVEL is solid. Implementations must be immutable and
// thread-safe (they are shared read-only with the worker threads that build chunk meshes).
class DensityField {
public:
    virtual ~DensityField() = default;
    // Density without player edits, in body-local space.
    virtual float baseDensity(Vector3 p) const = 0;
    // `base` is the unedited density at p (lets materials know "original surface" vs "dug out").
    virtual VoxelMaterial material(Vector3 p, float base) const = 0;
    // Used when the gradient vanishes (e.g. deep inside a flat region).
    virtual Vector3 fallbackUp(Vector3 p) const = 0;

    // Density including edits (the list should already be filtered to ones near p).
    float density(Vector3 p, const TerrainEdit* edits, size_t editCount) const;
    // Outward surface normal: negative density gradient, by central differences.
    Vector3 normal(Vector3 p, const TerrainEdit* edits, size_t editCount, float step = 0.5f) const;
};

// Asteroid / small body: a displaced sphere centred on the local origin.
class PlanetoidField : public DensityField {
public:
    SimplexNoise shape;    // overall form
    SimplexNoise detail;   // ore veins / ice patches
    float radius = 50.0f;  // the surface never reaches past this distance from the centre
    static constexpr float MEAN_SURFACE = 0.72f; // mean surface radius as a fraction of `radius`
    float iceAmount = 0.0f;// 0..1: how much of the body is ice

    PlanetoidField(unsigned int seed, float radius);
    float baseDensity(Vector3 p) const override;
    VoxelMaterial material(Vector3 p, float base) const override;
    Vector3 fallbackUp(Vector3 p) const override;

private:
    PlanetoidField(std::mt19937 rng, float radius);
};

// Open-world terrain: +Y up, rolling hills and mountains from a heightmap, carved by 3D cave
// tunnels, with grass/dirt near the surface, snow on high ground and crystal veins underground.
// (Materials reuse the VoxelMaterial slots: Regolith = grass/soil, Ice = snow, Ore = crystal.)
class TerrainField : public DensityField {
public:
    SimplexNoise hills, detail, caves, veins;
    float minHeight = -48.0f, maxHeight = 140.0f; // the surface stays inside this band
    float snowLine = 70.0f;

    explicit TerrainField(unsigned int seed);
    float surfaceHeight(float x, float z) const;
    float baseDensity(Vector3 p) const override;
    VoxelMaterial material(Vector3 p, float base) const override;
    Vector3 fallbackUp(Vector3) const override { return { 0, 1, 0 }; }
};

// CPU-side result of meshing one chunk (produced on a worker thread).
struct ChunkMeshData {
    uint64_t planetoidId = 0;
    NodeKey key{};
    uint32_t version = 0;
    std::vector<float> positions;          // xyz, relative to key.origin(), in voxels
    std::vector<float> normals;            // xyz
    std::vector<unsigned char> colors;     // rgba material weights
    std::vector<unsigned short> indices;
    bool overflow = false;                 // > 65535 vertices: mesh dropped
    float buildMs = 0.0f;
};

// Samples the field over node `key` (CHUNK_SIZE+1 cells at the node stride, i.e. one cell of
// overlap with the +x/+y/+z neighbours to hide cracks between LOD levels) and runs marching cubes.
// Thread-safe: touches only its arguments.
void BuildChunkMesh(const DensityField& field, const EditList& allEdits, ChunkMeshData& out);

// Edits whose influence (radius + margin) overlaps the axis-aligned box [lo, hi].
void FilterEdits(const EditList& all, Vector3 lo, Vector3 hi, float margin, EditList& out);

#endif
