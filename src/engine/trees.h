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
// Trees draw to kTreeDrawDistance (2 km), bushes to ~300 m, both with a dithered fade; past the trees
// the terrain shader paints their canopy in (lighting.fs, treeDistance).
//
// Collision: each tree trunk is a vertical capsule for the player and the ship (bushes are soft).
//
// Interaction (step 3): every plant has a stable id (its tile and cell). Chopping or harvesting
// records the id with the world time; it vanishes at once, its tile rebuilds without it, and after
// kRegrowDelay it comes back as a sapling that reaches full size kRegrowTime later. Removals are
// saved with the world.

#include <raylib.h>
#include <atomic>
#include <cmath>
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

// Tree variants are species x 4 size classes (#274: small, medium, large, extra large, each built by its own rules):
// 0-3 oak, 4-7 birch, 8-11 maple (broadleaf), 12-15 pine, 16-19 spruce
// (conifers); then 20-22 bushes and 23 loose stones.
constexpr int kTreeVariants = 24;
constexpr int kFirstBush = 20;
constexpr int kTreeClasses = 4;              // small, medium, large, extra large (old-growth giants)
constexpr int kStoneVariant = 23;           // loose stones on the ground: picked up for stone (the first stone tools)
// Wood species: what a tree drops when felled (its log) and how it is drawn (bark and leaf tints
// over the shader's latitude colours). An area has a dominant species (60%) and two others (20% each).
enum class WoodSpecies : uint8_t { Oak, Birch, Maple, Pine, Spruce, Count };
struct WoodSpeciesDef { const char* log; const char* name; bool conifer; Vector3 bark, leaf; };
const WoodSpeciesDef& WoodSpeciesInfo(WoodSpecies s);
inline WoodSpecies SpeciesOfVariant(int v) { return (WoodSpecies)(v / kTreeClasses); }  // trees only (v < kFirstBush)
inline int SizeClassOfVariant(int v) { return v % kTreeClasses; }                          // 0 small .. 3 extra large

constexpr float kTreeTileSize = 64.0f;      // m
constexpr float kTreeInterestRange = 192.0f; // m: tiles kept around other interest centres (other players, #318)
inline float kTreeDrawDistance = 2000.0f;    // trees (simple meshes past kTreeNearDistance); --tree-dist
constexpr float kBushDrawDistance = 300.0f;
constexpr float kTreeNearDistance = 200.0f; // full meshes inside this, simple ones beyond
// Shadow casters (#204 perf): trees within these of the camera go into the far / near cascade passes (boxes +-200 m and
// +-32 m, plus long low-sun shadows and the 40 m the camera may move between batch rebuilds).
constexpr float kTreeShadowFarReach = 355.0f, kTreeShadowNearReach = 115.0f; // from the sun axis through the camera: box half-diagonal + 40 m drift + 30 m tree
constexpr int kTreeSets = 3; // draw set 0: all; 1: far-cascade casters; 2: near-cascade casters

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
    bool stone() const { return variant == kStoneVariant; }
    // A tree trunk's volume relative to a medium broadleaf at scale 1 (0.3 m radius, ~5.2 m): chop hits and logs (#274).
    float trunkVolume() const { return tree() ? radius * radius * height / (0.09f * 5.2f) : 0.0f; }
    // Chop work and yield grow a little slower than the volume (volume^0.7): a medium oak ~14 hits / 5 logs, a giant ~55 / 26.
    float chopSize() const { return std::pow(trunkVolume(), 0.7f); }
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

// (#893 tests) the scatter around dir as generated (no edits or removals), counted within range m: trees, bushes, loose stones.
class PlanetField;
void Trees_CountScatter(const PlanetField& f, uint32_t seed, Vector3 dir, float range, int& trees, int& bushes, int& stones);
class TreeSystem {
public:
    TreeSystem(std::shared_ptr<const PlanetField> field, unsigned int seed);
    ~TreeSystem();

    // Main thread, once a frame: request tiles near the camera (body-local), drop far ones, take
    // finished tiles, rebuild the draw batches when needed.
    void update(Planetoid& body, Vector3 camLocal);
    // Other interest centres on this body (body-local; #318: other players on a server): their nearby tiles stay loaded.
    void setInterest(std::vector<Vector3> centresLocal) { interest = std::move(centresLocal); }
    int readyTilesNear(Vector3 local, float range) const; // built tiles whose centre is within range (arc, m)
    // Instanced draw (tree shader bound by the scene; its view-projection set).
    void draw(const Planetoid& body, const RenderContext& ctx, int onlyLod = -1, int set = 0); // onlyLod 1: simple far meshes only (the mid pass); set: 0 all, 1/2 the far/near shadow casters
    // The shadow cascades' axis (body-local direction toward the sun): the caster sets are trees near it. Rebuilds the
    // sets when it turns more than ~1 degree.
    void setShadowSun(Vector3 dirLocal) { if (dirLocal.x * shadowSunLocal.x + dirLocal.y * shadowSunLocal.y + dirLocal.z * shadowSunLocal.z < 0.99985f) { shadowSunLocal = dirLocal; batchDirty = true; } }
    // A terrain edit: rebuild tiles it reaches (trees lose their ground, new ground may grow some).
    void onEdit(Vector3 localCentre, float radius, std::shared_ptr<const TerrainDeltas> edits);
    // Trunks within `range` of a body-local point.
    void trunksNear(Vector3 local, float range, std::vector<Trunk>& out) const;
    void plantsNear(Vector3 local, float range, std::vector<Plant>& out) const; // trees, bushes and stones (tests)
    // The nearest loaded plant (tree or bush) to a body-local point.
    bool nearestPlant(Vector3 local, float range, bool tree, Plant& out, bool stone = false, int sizeClass = -1) const; // stone: loose stones (else bushes when !tree); sizeClass: trees of that class only (#274)
    // The nearest plant a body-local ray hits within maxDistance.
    bool raycast(Vector3 localOrigin, Vector3 localDir, float maxDistance, PlantHit& hit) const;
    // Chop / harvest: gone now, regrows later.
    void remove(uint64_t id, double when = -1.0); // when: the world time it was removed (< 0: now; a joined client: the server's, #335)
    bool findPlant(uint64_t id, Plant& out) const; // a plant standing in a loaded tile (#335: the server checks a chop)
    const PlantRemovals& removals() const { return *removed; }
    std::shared_ptr<const PlantRemovals> removalsShared() const { return removed; } // (#335: the server diffs it; copy-on-write)
    void mergeRemovals(const std::vector<std::pair<uint64_t, double>>& list); // (#335) a joined client: the server's, many at once
    void setRemovals(const PlantRemovals& r);  // from a save: rebuilds every tile
    // Keep an area free of plants (the start: where the player stands and the ship parks).
    void addClearing(Vector3 localCentre, float radius);

    int loadedTiles() const { return (int)tiles.size(); }
    int instanceCount() const { return drawnInstances; }
    // Tests: trees per species in loaded tiles within `range` of a body-local point.
    void speciesCounts(Vector3 local, float range, int out[(int)WoodSpecies::Count], int* classes = nullptr) const; // classes: kTreeClasses counts (#274)

private:
    struct Shared;                             // results queue shared with worker jobs
    std::shared_ptr<Shared> shared;
    std::shared_ptr<const PlanetField> field;
    unsigned int seed;
    float radius;
    std::unordered_map<TreeTileKey, TreeTile, TreeTileKeyHash> tiles;
    std::shared_ptr<const TerrainDeltas> edits; // the body's edits as last seen (tiles built with them)
    std::shared_ptr<const PlantRemovals> removed; // copy-on-write, shared with tile builds
    std::shared_ptr<const std::vector<Vector4>> clearings; // xyz centre, w radius (body-local)

    // Draw batches.
    DVec3 anchor{};                            // body-local origin of the batch positions
    bool batchDirty = true;
    double lastBatch = -1.0;
    // Two detail levels: [0..kTreeVariants) full meshes near the camera, then simple far meshes.
    unsigned int vao[kTreeSets][2 * kTreeVariants] = {}, vbo[kTreeSets][2 * kTreeVariants] = {};
    int count[kTreeSets][2 * kTreeVariants] = {};
    Vector3 shadowSunLocal = { 0, 0, 0 };     // (0: unknown, the caster sets use distance from the camera)
    // Main-view GPU culling (#292): per slot of set 0, the visible instances (compacted) and its VAO; one command buffer.
    unsigned int culledVao[2 * kTreeVariants] = {}, culledVbo[2 * kTreeVariants] = {}, cmdBuf = 0;
    Vector3 lastWantCam = { 1e30f, 0, 0 };    // tile set recomputed after the camera moves 16 m
    std::vector<Vector3> interest, lastInterest; // other interest centres (#318) and where they were at the last recompute
    int drawnInstances = 0;

    void request(const TreeTileKey& key, TreeTile& tile, float priority);
    void rebuildBatches(Vector3 camLocal);
};

// Tile helpers (exposed for tests).
TreeTileKey TreeTileOf(Vector3 dir, float radius);
Vector3 TreeTileCentreDir(const TreeTileKey& key, float radius);
// Fraction of the ground under tree crowns around a body-local direction (the same forest field
// and odds as the placement), for painting distant forests into the terrain.
class PlanetField;
float TreeCanopy(Vector3 dir, float radius, float elevation, unsigned int seed, const PlanetField* field = nullptr); // field: by its biomes (Plan 5 B)

#endif
