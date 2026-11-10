#ifndef TREES_H
#define TREES_H

// Trees and bushes (design doc section 11, foliage step 2).
//
// Shapes: a handful of variants built procedurally once from connected rectangular prisms (a
// trunk of stacked tapering segments, branches, boxy leaf blocks; conifers as stacked leaf slabs;
// bushes as clusters of leaf blocks), drawn instanced with a per-instance scale, turn and up.
//
// Placement: independent of the terrain mesh (trees must show far beyond the finest chunks), on a
// cube-sphere grid of 64 m tiles. A tile is built on a worker from the seed alone: 4 m cells, a
// low-frequency forest field decides trees and bushes, and a candidate survives on grass or dirt,
// above the beach and below the tree line, on gentle slopes, with solid ground underfoot and open
// air above (terrain edits included). Same seed, same forests; nothing is stored.
//
// Drawing: per body, one instance buffer per variant over all loaded tiles, rebuilt (a few times a
// second at most) when tiles arrive or leave, positioned relative to an anchor near the camera.
// Trees draw to ~1 km, bushes to ~300 m, both with a dithered fade.
//
// Collision: each tree trunk is a vertical capsule for the player and the ship (bushes are soft).
//
// Interaction (step 3): every plant has a stable id (its tile and cell). Chopping or harvesting
// records the id with the world time; it vanishes at once, its tile rebuilds without it, and after
// kRegrowDelay it comes back as a sapling that reaches full size kRegrowTime later. Removals are
// saved with the world.

#include <raylib.h>
#include <atomic>
#include <cstdint>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <vector>
#include "dvec3.h"
#include "voxel.h"

class PlanetField;
class Planetoid;
struct RenderContext;

constexpr int kTreeVariants = 9;            // 0-3 broadleaf, 4-5 conifer, 6-8 bush
constexpr int kFirstBush = 6;
constexpr float kTreeTileSize = 64.0f;      // m
constexpr float kTreeDrawDistance = 1000.0f;
constexpr float kBushDrawDistance = 300.0f;
constexpr float kTreeNearDistance = 200.0f; // full meshes inside this, simple ones beyond

struct TreeTileKey {
    int face, i, j;
    bool operator==(const TreeTileKey& o) const { return face == o.face && i == o.i && j == o.j; }
};
struct TreeTileKeyHash {
    size_t operator()(const TreeTileKey& k) const { return ((size_t)k.face * 73856093u) ^ ((size_t)k.i * 19349663u) ^ ((size_t)k.j * 83492791u); }
};

constexpr double kRegrowDelay = 600.0;      // s of world time before a sapling appears
constexpr double kRegrowTime = 1200.0;      // s from sapling to full size

// A plant for picking (chop / harvest), in body-local space.
struct Plant {
    uint64_t id;
    Vector3 base, up;
    float height, radius;  // trees: the trunk; bushes: the whole bush
    float scale;
    int variant;
    bool tree() const { return variant < kFirstBush; }
};
struct PlantHit {
    Plant plant;
    float distance;        // along the ray
};
using PlantRemovals = std::unordered_map<uint64_t, double>; // id -> world time it was removed

// A trunk for collision, in body-local space (metres).
struct Trunk {
    Vector3 base;
    Vector3 up;
    float height, radius;
};

struct TreeTile {
    bool ready = false;
    uint32_t version = 0;
    Vector3 centre{};                          // body-local, on the surface
    std::vector<float> inst[kTreeVariants];    // 8 floats each: body-local position, scale, up, hash
    std::vector<Trunk> trunks;
    std::vector<Plant> plants;                 // every tree and bush (picking)
    std::vector<uint64_t> ids[kTreeVariants];  // per instance, parallel to inst
    bool regrowing = false;                    // has saplings: rebuild now and then as they grow
    double builtAt = 0.0;                      // world time the build used
    bool pending = false;                      // a rebuild is queued
};

class TreeSystem {
public:
    TreeSystem(std::shared_ptr<const PlanetField> field, unsigned int seed);
    ~TreeSystem();

    // Main thread, once a frame: request tiles near the camera (body-local), drop far ones, take
    // finished tiles, rebuild the draw batches when needed.
    void update(Planetoid& body, Vector3 camLocal);
    // Instanced draw (tree shader bound by the scene; its view-projection set).
    void draw(const Planetoid& body, const RenderContext& ctx);
    // A terrain edit: rebuild tiles it reaches (trees lose their ground, new ground may grow some).
    void onEdit(Vector3 localCentre, float radius);
    // Trunks within `range` of a body-local point.
    void trunksNear(Vector3 local, float range, std::vector<Trunk>& out) const;
    // The nearest plant a body-local ray hits within maxDistance.
    bool raycast(Vector3 localOrigin, Vector3 localDir, float maxDistance, PlantHit& hit) const;
    // Chop / harvest: gone now, regrows later.
    void remove(uint64_t id);
    const PlantRemovals& removals() const { return *removed; }
    void setRemovals(const PlantRemovals& r);  // from a save: rebuilds every tile
    // Keep an area free of plants (the start: where the player stands and the ship parks).
    void addClearing(Vector3 localCentre, float radius);

    int loadedTiles() const { return (int)tiles.size(); }
    int instanceCount() const { return drawnInstances; }

private:
    struct Shared;                             // results queue shared with worker jobs
    std::shared_ptr<Shared> shared;
    std::shared_ptr<const PlanetField> field;
    unsigned int seed;
    float radius;
    std::unordered_map<TreeTileKey, TreeTile, TreeTileKeyHash> tiles;
    std::shared_ptr<const EditList> edits;     // the body's edits as last seen (tiles built with them)
    std::shared_ptr<const PlantRemovals> removed; // copy-on-write, shared with tile builds
    std::shared_ptr<const std::vector<Vector4>> clearings; // xyz centre, w radius (body-local)

    // Draw batches.
    DVec3 anchor{};                            // body-local origin of the batch positions
    bool batchDirty = true;
    double lastBatch = -1.0;
    // Two detail levels: [0..kTreeVariants) full meshes near the camera, then simple far meshes.
    unsigned int vao[2 * kTreeVariants] = {}, vbo[2 * kTreeVariants] = {};
    int count[2 * kTreeVariants] = {};
    Vector3 lastWantCam = { 1e30f, 0, 0 };    // tile set recomputed after the camera moves 16 m
    int drawnInstances = 0;

    void request(const TreeTileKey& key, TreeTile& tile, float priority);
    void rebuildBatches(Vector3 camLocal);
};

// Tile helpers (exposed for tests).
TreeTileKey TreeTileOf(Vector3 dir, float radius);
Vector3 TreeTileCentreDir(const TreeTileKey& key, float radius);
// Fraction of the ground under tree crowns around a body-local direction (the same forest field
// and odds as the placement), for painting distant forests into the terrain.
float TreeCanopy(Vector3 dir, float radius, float elevation, unsigned int seed);

#endif
