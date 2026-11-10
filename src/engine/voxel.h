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
// Eight materials. Meshes carry them as per-vertex blend weights: 0-3 in the vertex colour (RGBA),
// 4-5 in texcoord (xy), 6-7 in texcoord2 (xy), so they interpolate smoothly across triangles.
// Asteroids use 0-3; planets all eight (Regolith is rendered as dirt there).
enum class VoxelMaterial : uint8_t { Rock = 0, Regolith = 1, Ice = 2, Ore = 3, Sand = 4, Gravel = 5, Grass = 6, Snow = 7, Count };
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
    // False where player edits have no effect (a planet's mantle).
    virtual bool editable(Vector3) const { return true; }
    // Upper bound on the surface's distance from the local origin (LOD root size, spawn raycasts).
    virtual float maxRadius() const = 0;
    // Lower bound: no surface lies entirely inside this radius (a planet's mantle). 0 = unknown.
    virtual float minRadius() const { return 0.0f; }
    // Radius of a sphere that is always opaque (nothing can be dug below it): terrain hidden behind
    // it from the camera is culled. 0 = no occluder.
    virtual float occluderRadius() const { return 0.0f; }
    // Sea level: a water surface at this radius wherever it meets open air (0 = no sea).
    virtual float seaLevelRadius() const { return 0.0f; }

    // Density for samples `footprint` apart (chunk meshing at a LOD stride): fields may leave out
    // detail finer than that. Defaults to baseDensity; footprint 0 = full detail.
    virtual float baseDensityLod(Vector3 p, float /*footprint*/) const { return baseDensity(p); }
    // A whole S^3 block of baseDensityLod samples at origin + (x, y, z) * stride, x fastest (chunk
    // meshing). Fields can share work between samples (a planet evaluates its height map once per
    // column); the default just loops.
    virtual void baseDensityBlock(Vector3 origin, float stride, int S, float footprint, float* out) const;
    // Material of a mesh vertex sampled at that footprint (so coarse vertices, which sit on the
    // simplified surface, are judged against it). Defaults to material().
    // `normal` = the vertex's outward surface normal (slope-dependent materials).
    virtual VoxelMaterial materialLod(Vector3 p, float base, float /*footprint*/, Vector3 /*normal*/) const { return material(p, base); }
    // materialLod plus how exposed the vertex is to the sky (1 open, 0 buried: caves, deep pits), which
    // scales its sunlight. Fields that know their surface override it (computed together, once).
    virtual VoxelMaterial materialSky(Vector3 p, float base, float footprint, Vector3 normal, float& sky) const {
        sky = 1.0f;
        return materialLod(p, base, footprint, normal);
    }

    // Density including edits (the list should already be filtered to ones near p).
    float density(Vector3 p, const TerrainEdit* edits, size_t editCount, float footprint = 0.0f) const;
    // Adds the edits' contributions to a base density at p (nothing where !editable(p)).
    float applyEdits(Vector3 p, float base, const TerrainEdit* edits, size_t editCount) const;
    // Outward surface normal: negative density gradient, by central differences.
    Vector3 normal(Vector3 p, const TerrainEdit* edits, size_t editCount, float step = 0.5f, float footprint = 0.0f) const;
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
    float maxRadius() const override { return radius; }

private:
    PlanetoidField(std::mt19937 rng, float radius);
};

// Planet or moon: a sphere of `radius` (sea level) with a height map over it (continents, ocean
// floor, ridged mountains, hills, river valleys), 3D caves, ore veins that get richer with depth,
// and an unmineable mantle 1.5% of the radius below sea level. Same uniform density gradient as
// asteroids. Moons: no oceans or snow, craters instead of ridges.
class PlanetField : public DensityField {
public:
    float radius;                 // sea level (m)
    float mantleRadius;           // edits stop below this
    bool moon;
    float minElevation, maxElevation; // the surface stays inside [radius + min, radius + max]
    SimplexNoise continents, ridges, hills, valleys, detail, cavesA, cavesB, veins, climate, caverns, ravines, ravineMask;
    static constexpr float kTunnelLevel = 0.096f; // worm tunnel where both cave fields are within this of zero
    static constexpr float kCavernLevel = 0.68f;  // cavern where the cavern field exceeds this (rare)

    PlanetField(unsigned int seed, float radius, bool moon);
    // Surface height above sea level along a unit direction; footprint > 0 leaves out detail finer
    // than samples that far apart (see SimplexNoise::fractalLimited).
    float elevation(Vector3 dir, float footprint = 0.0f) const;
    // A body-local direction to a good start spot: temperate latitude, flat grassland a little above
    // the sea, solid ground underfoot (no cave mouth). Deterministic for a seed.
    Vector3 findSpawnDir(unsigned int seed) const;
    float baseDensity(Vector3 p) const override { return baseDensityLod(p, 0.0f); }
    float baseDensityLod(Vector3 p, float footprint) const override;
    void baseDensityBlock(Vector3 origin, float stride, int S, float footprint, float* out) const override;
    VoxelMaterial material(Vector3 p, float base) const override { return materialLod(p, base, 0.0f, fallbackUp(p)); }
    VoxelMaterial materialLod(Vector3 p, float base, float footprint, Vector3 normal) const override {
        float sky;
        return materialSky(p, base, footprint, normal, sky);
    }
    VoxelMaterial materialSky(Vector3 p, float base, float footprint, Vector3 normal, float& sky) const override;
    Vector3 fallbackUp(Vector3 p) const override;
    bool editable(Vector3 p) const override;
    float maxRadius() const override { return radius + maxElevation; }
    float minRadius() const override { return mantleRadius - 1.0f; } // digging stops at the mantle
    float occluderRadius() const override { return mantleRadius; }
    float seaLevelRadius() const override { return moon ? 0.0f : radius; }

private:
    float densityNear(Vector3 p, float d, float surface, float footprint) const; // after the height map
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
    float maxRadius() const override { return 1e9f; } // unbounded horizontally
};

// CPU-side result of meshing one chunk (produced on a worker thread).
struct ChunkMeshData {
    uint64_t planetoidId = 0;
    NodeKey key{};
    uint32_t version = 0;
    std::vector<float> positions;          // xyz, relative to key.origin(), in voxels
    std::vector<float> normals;            // xyz
    std::vector<unsigned char> colors;     // material weights 0-3 (rgba)
    std::vector<float> texcoords;          // material weights 4-5
    std::vector<float> texcoords2;         // material weights 6-7
    std::vector<float> tangents;           // (0, 0, 0, sky exposure) per vertex
    std::vector<float> waterPositions;     // sea surface where it meets open air (same frame as positions)
    std::vector<float> waterNormals;
    std::vector<unsigned short> waterIndices;
    std::vector<unsigned short> indices;
    bool overflow = false;                 // > 65535 vertices: mesh dropped
    float buildMs = 0.0f;
    float sampleMs = 0.0f;                 // part of buildMs spent sampling the density field
};

// Samples the field over node `key` (CHUNK_SIZE+1 cells at the node stride, i.e. one cell of
// overlap with the +x/+y/+z neighbours to hide cracks between LOD levels) and runs marching cubes.
// Thread-safe: touches only its arguments.
void BuildChunkMesh(const DensityField& field, const EditList& allEdits, ChunkMeshData& out);

// Edits whose influence (radius + margin) overlaps the axis-aligned box [lo, hi].
void FilterEdits(const EditList& all, Vector3 lo, Vector3 hi, float margin, EditList& out);

#endif
