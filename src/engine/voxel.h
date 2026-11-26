#ifndef VOXEL_H
#define VOXEL_H

#include <raylib.h>
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <vector>
#include <mutex>
#include <atomic>
#include <unordered_set>

#include "SimplexNoise.h"
#include "strata.h"
#include "worldgenparams.h"

extern const int CHUNK_SIZE;          // cells per chunk edge (samples per edge = CHUNK_SIZE + 1)
// The finest terrain cell (body-local units, metres on planets): the edit lattice spacing and the
// octree's level-0 stride. 0.25 m makes the dig spheres the tools show (0.35-1.2 m) exactly what is
// carved; it was 1 m (a 0.4 m sphere moved the nearest lattice point, so holes depended on the grid).
#ifdef RAYLIBGAME_VOXEL_1M   // a comparison build at the old 1 m lattice (tests only)
constexpr float kVoxelSize = 1.0f;
constexpr int kVoxelShift = 0;
#else
constexpr float kVoxelSize = 0.25f;
constexpr int kVoxelShift = 2;                  // log2(1 / kVoxelSize): levels finer than the old 1 m cell
#endif
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


// Octree LOD node. A node at `level` covers CHUNK_SIZE * 2^level voxels (kVoxelSize each) per axis, meshed with
// CHUNK_SIZE cells (stride 2^level), so every node costs about the same to build and draw.
// (x, y, z) are in units of the node size at that level.
struct NodeKey {
    int level, x, y, z;
    bool operator==(const NodeKey& o) const { return level == o.level && x == o.x && y == o.y && z == o.z; }
    float stride() const { return kVoxelSize * (float)(1 << level); }   // metres between samples
    int latticeStride() const { return 1 << level; }                      // edit-lattice points between samples
    float size() const { return kVoxelSize * (float)(CHUNK_SIZE << level); }
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
// Ids 0-7 are the legacy weights (per-chunk draws); later ones ride the arena's per-chunk palette (#480: up to 256).
enum class VoxelMaterial : uint8_t { Rock = 0, Regolith = 1, Ice = 2, Ore = 3, Sand = 4, Gravel = 5, Grass = 6, Snow = 7, Clay = 8,
                                    Granite = 9, Basalt = 10, Sandstone = 11, Limestone = 12, Marble = 13,
                                    Pebbles = 14, Silt = 15, Mud = 16, Farmland = 17, Count }; // 9-13: rock strata (#587); 14-16: river beds, floodplains, wet shores (#595); 17: tilled dirt (#1085)
enum class Biome : uint8_t { Ocean = 0, Jungle, Savanna, Desert, Badlands, Forest, Meadow, Swamp, Taiga, Tundra, Alpine, Polar, Count };
// Bench override (--bench-biome): findSpawnDir returns the nearest deep patch of this biome (-1: off).
inline int gSpawnBiome = -1;
// Bench override (--bench-pond): findSpawnDir returns the rim of the nearest pond (#251 crack checks).
inline bool gSpawnPond = false;
// Erosion bake (Plan 7, #355): RAYLIBGAME_EROSION=1 bakes the whole planet at world creation (erosion.h); off until the plan's gate.
inline bool gErosionPrototype = true; // erosion-baked terrain by default (user 2026-11-15, #362); RAYLIBGAME_EROSION=0: the old noise terrain
// The 3D surface band (cliffs, overhangs: terrain3D). Off while the erosion bake is being tested (user, 2026-11-13: keep
// it out of the erosion comparisons); RAYLIBGAME_3D=1 turns it back on.
inline bool gTerrain3D = true;
inline bool gGullies = true; // erosion-filter gullies on the baked height (#357); RAYLIBGAME_GULLIES=0 for A/B
inline bool gCliffFaces = true; // one face per steep resistant layer (#588, Plan 9 step 5); RAYLIBGAME_CLIFF_FACES=0 for A/B
inline bool gTerrainAO = true; // per-vertex ambient occlusion in chunk meshing (#132); RAYLIBGAME_AO=0 for A/B
inline bool gSkyProbes = true; // sky light under overhangs from short rays (#369); RAYLIBGAME_SKY_PROBES=0 for A/B
class PlanetErosion;
const char* BiomeName(Biome b);
const char* VoxelMaterialName(VoxelMaterial m);

// A spherical density brush in planetoid-local space. strength < 0 digs, > 0 adds.
struct TerrainEdit {
    Vector3 center;
    float radius;
    float strength;
    // Carve: subtract exactly this sphere (CSG, like the caves): density becomes at most
    // ISO + gradient * (distance - radius), a hard edge where the surface meets it. Strength unused.
    bool carve = false;
    bool softOnly = false; // carve: only soft ground (dirt, sand, gravel, snow); rock stays (no pickaxe)
    // Fill (#228): add exactly this sphere (CSG union, sharp-edged): density becomes at least ISO - gradient * (distance -
    // radius); `paint` colours what it adds. Strength unused.
    bool fill = false;
    int paint = -1;        // adding: the VoxelMaterial the added ground is made of (-1: the field's own)
    // Till (#1085): no change to the ground's shape; soil (grass, dirt, clay, silt, mud) within `radius` of the centre and
    // within half a metre under the surface (and the lattice points just over it) becomes Farmland. A free brush.
    bool till = false;
    // Carve an oriented box instead of the sphere (#210: a part's cells, exactly): half extents along the unit axes
    // ax/ay/az (body-local), centred on `center`. Set `radius` to its bounding reach (the half extents' length): the
    // dirty-chunk, tree and fluid bounds use it.
    bool box = false;
    Vector3 half{ 0, 0, 0 }, ax{ 1, 0, 0 }, ay{ 0, 1, 0 }, az{ 0, 0, 1 };
};
#include "terraindeltas.h"

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

    // Density including the player's edits (base + trilinear delta; deltas may be null).
    float density(Vector3 p, const TerrainDeltas* deltas, float footprint = 0.0f) const;
    // Outward surface normal: negative density gradient, by central differences.
    Vector3 normal(Vector3 p, const TerrainDeltas* deltas, float step = 0.5f, float footprint = 0.0f) const;
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
class RiverNetwork;
class PlanetField : public DensityField {
public:
    float radius;                 // sea level (m)
    float mantleRadius;           // edits stop below this
    bool moon;
    float minElevation, maxElevation; // the surface stays inside [radius + min, radius + max]
    WorldGenParams gen;           // the world's generation parameters (#504): defaults = the terrain as before
    uint32_t caveEpoch = 0;       // #1055: unique per field and land (renewed by setGenParams / setSandbox): part of the cave cells' key
    bool caveBaked() const;       // #1055: whether the erosion bake is visible here (part of the cave cells' key)
    // Before any sampling or the erosion bake: the parameters, and the elevation bounds they need (taller ranges).
    void setGenParams(const WorldGenParams& p);
    bool foliage = false;         // grows trees (planets with life); the terrain paints distant canopy
    // A sandbox (#662, --sandbox): a flat plane kSandboxHeight over the radius - no height map, rivers, ponds, caves, rocks
    // or sea - for sim benches that shape their own ground with edits. Before any sampling.
    static constexpr float kSandboxHeight = 50.0f;
    bool sandbox = false;
    void setSandbox();
    unsigned int foliageSeed = 0; // the tree system's seed (TreeCanopy matches its forests)
    unsigned int oreSeed = 0;     // ore deposits (oreKind)
    SimplexNoise continents, ridges, hills, valleys, detail, cavesA, cavesB, veins, climate, caverns, ravines, ravineMask, streams, rocks, rockMask, rugged, cliffs, warp3, tempNoise, moistNoise, biomeJitter, caveMouths;
    static constexpr float kTunnelLevel = 0.096f; // worm tunnel where both cave fields are within this of zero
    static constexpr float kCavernLevel = 0.68f;  // cavern where the cavern field exceeds this (rare)

    PlanetField(unsigned int seed, float radius, bool moon);
    // Surface height above sea level along a unit direction; footprint > 0 leaves out detail finer
    // than samples that far apart (see SimplexNoise::fractalLimited).
    float elevation(Vector3 dir, float footprint = 0.0f) const;
    std::string debugColumn(Vector3 dir) const; // tests: the height terms at dir (land, base, rivers, lakes, ponds)
    float elevationBase(Vector3 dir, float footprint = 0.0f) const; // land + detail + biome landforms: no rivers, lakes, ponds (Plan 6: GPU port check)
    // The terrain noise fields in a fixed order for the GPU port (terrain.glsl): frequency + seed per field.
    void gpuNoiseParams(std::vector<float>& out) const;
    // Inland water (planets): rivers of the drainage network (rivers.h) in channels they carve, lakes
    // in filled depressions, and
    // ponds (12-38 m, flat surface, a low embankment where the ground falls away) in some
    // cells of a ~380 m grid on gentle inland ground. The surface's height above sea level here, or
    // kNoWater. The sea is separate (sea level).
    static constexpr float kNoWater = -1e9f;
    float waterElevation(Vector3 dir, float footprint = 0.0f) const;
    // The ground beside a river before the river shapes it (#640): elevation() without the channel and with the gullies
    // and cliff faces faded out, as elevation() does near a river. The bake-native network holds its levels under it.
    float riverBankHeight(Vector3 dir) const;
    // The generated water's level for the chunk water lattice (#609): waterElevation, plus a narrow stream's level out to
    // about half a lattice step beside it where the ground stands over the water (the terrain hides it there), so a
    // lattice of nodes 4 samples apart does not miss a stream narrower than its spacing.
    float waterElevationMesh(Vector3 dir, float footprint) const;
    // Whether inland water may lie within `extent` metres of a direction (cheap; meshing skips the
    // per-sample work elsewhere).
    bool inlandWaterNear(Vector3 dir, float extent, float footprint = 0.0f) const;
    // Ore at a point `depth` m under the original surface: an index into kOreKinds, or -1.
    //   Deposits: a 96 m grid; about a third of its cells hold one lumpy blob 7-18 m across of a single
    //   kind (coal, iron, copper, tin; bauxite below 20 m, ilmenite below 60 m) - a find worth mining.
    //   Scattered: small pockets of the common ores (iron, copper, tin) and coal everywhere, most near the surface.
    static constexpr int kOreKindCount = 6;
    static const char* OreKindName(int kind);   // "iron_ore", "copper_ore", "tin_ore", "coal", "bauxite_ore", "ilmenite_ore"
    int oreKind(Vector3 p, float depth) const;
    int depositNear(Vector3 p, float reach) const; // the kind of a deposit within `reach` m, or -1 (miners)
    // The ore deposit whose cell holds p (#207: a vein mined out): its centre, radius (the blob reaches ~1.35x it with its
    // lumpiness) and kind; false if the cell has none.
    bool depositAt(Vector3 p, Vector3& centre, float& radius, int& kind) const;
    // The kind of stone that rock is here (M21, building variety): an item id - marble in rare deep
    // pockets, basalt deep down and under the sea, granite in the mountains, sandstone on low coastal
    // land, limestone elsewhere (in broad regions).
    const char* stoneKind(Vector3 p) const;
    // Rock strata (#587): the column here (from region-scale fields, as the erosion bake sees it) and the rock at a
    // body-local point (its layer at the point's height, in the bake's simulation heights).
    StrataColumn strataColumn(Vector3 dir) const;
    RockType rockTypeAt(Vector3 p) const;
    // Tests (#561): true on this thread = the bake is ignored (erosion() answers nullptr), so land() / elevationBase()
    // give the noise terrain the GPU point kernel ports (it has no bake).
    static inline thread_local bool tIgnoreBake = false;
    // True while land() answers for the erosion bake's input (#544: canyons go in before the bake).
    static inline thread_local bool tBakeInput = false;
    float canyonFault(float sx, float sy, float sz) const; // #544: 0..1 along the canyon lines
    float canyonLine(float sx, float sy, float sz) const;  // #544: |the ravine field| (0 on a canyon line; the notch is 0.02 wide)
    float faultUplift(Vector3 s) const { return Strata_FaultUplift(s, fieldSeed); } // #590: metres (s = surface point)
    float foldWeight(Vector3 s) const { return Strata_FoldWeight(s, fieldSeed); }   // #591: 0..1
    float volcano(Vector3 s) const { return gen.volcanoes ? Strata_Volcano(s, radius, fieldSeed) : 0.0f; } // #593: metres
    float karstWeight(Vector3 s) const { return gen.karst ? Strata_KarstWeight(s, fieldSeed) : 0.0f; } // #592: 0..1
    unsigned int terrainSeed() const { return fieldSeed; } // #1248: the mana field is seeded from it
    float volcanoWeight(Vector3 s) const { float w = 0.0f; if (gen.volcanoes) Strata_Volcano(s, radius, fieldSeed, &w); return w; } // #1063: 0..1 toward a vent (0 off volcanoes)
    // #588: the cliff face of a steep resistant layer at surface point s (m), height h (final m), the bake's gradient.
    float cliffFace(Vector3 s, float h, Vector3 grad, float fp) const;
    // On a river or pond bed (or its wet shore): digging there also turns up clay.
    bool wetBed(Vector3 local) const;
    // The static water at a direction and its kind (flowing water, fluid.h): 0 none, 1 river, 2 lake,
    // 3 pond; level above sea level. Ponds carry an id and their centre and radius.
    struct WaterSource { int kind = 0; float level = kNoWater; uint64_t pond = 0; Vector3 pondDir{}; float pondR = 0.0f, pondDepth = 4.0f; };
    // Small lakes of the network (#360): finite like ponds. Their id (bit 63 set: never a pond cell id),
    // or 0 when the lake is large (infinite) or there is no bake.
    uint64_t drainableLake(int lakeId) const;
    WaterSource waterSourceAt(Vector3 dir) const;
    // Tests: the bake's ponds (#487: 309 ponds on seed 42, none near the start).
    std::vector<Vector3> bakePondsByDistance(Vector3 dir) const; // their centres, nearest first (empty without a bake)    // A body-local direction to a good start spot: temperate latitude, flat grassland a little above
    // the sea, solid ground underfoot (no cave mouth). Deterministic for a seed.
    Vector3 findSpawnDir(unsigned int seed) const;
    Vector3 nearestBiomeDir(Vector3 home, Biome b) const; // nearest deep patch of b (home if none)
    // Biomes (Plan 5 B): climate from latitude, height and regional noise (2-5 km patches); the biome
    // from it (Whittaker-style: temperature x moisture). biomeAt jitters the climate a little at a fine
    // scale so borders interleave instead of running as hard lines; tints use the smooth climate.
    struct Climate { float temp = 15.0f, moist = 0.5f; }; // temp in C (annual mean at the surface), moist 0..1
    Climate climateAt(Vector3 dir, float elev) const;
    Biome biomeAt(Vector3 dir, float elev) const;
    Biome biomeOf(const Climate& c, float elev) const;
    // Biome landforms (Plan 5 B) added to the height before rivers carve it: desert dunes, badlands
    // mesas (flat tops on steep risers), swamps flattened toward the water table.
    float biomeShape(Vector3 dir, float h, float fp) const;
    float ruggedness(Vector3 dir, float elev) const;
    float floatingIsland(Vector3 p, float d) const;  // Plan 5 C: metres inside a floating island (< 0 outside)
    int floatingIslandsNear(Vector3 dir, float range, std::vector<Vector3>& out) const;
    bool islandTopAt(Vector3 dir, float& topRadius) const; // a floating island over this column: its top (radius from the centre) // Plan 5 C: 0 flat (meadows, swamps) .. 1 dramatic (badlands, mountains)
    float baseDensity(Vector3 p) const override { return baseDensityLod(p, 0.0f); }
    float baseDensityLod(Vector3 p, float footprint) const override;
    void baseDensityBlock(Vector3 origin, float stride, int S, float footprint, float* out) const override;
    VoxelMaterial material(Vector3 p, float base) const override { return materialLod(p, base, 0.0f, fallbackUp(p)); }
    VoxelMaterial materialLod(Vector3 p, float base, float footprint, Vector3 normal) const override {
        float sky;
        return materialSky(p, base, footprint, normal, sky);
    }
    VoxelMaterial materialSky(Vector3 p, float base, float footprint, Vector3 normal, float& sky) const override { return materialSkyHinted(p, base, footprint, normal, sky, nullptr, nullptr); }
    Vector3 fallbackUp(Vector3 p) const override;
    bool editable(Vector3 p) const override;
    // The river drainage network (rivers.h), built on first use; nullptr on moons.
    const RiverNetwork* rivers() const;
    float maxRadius() const override { return radius + maxElevation; }
    float minRadius() const override { return mantleRadius - 1.0f; } // digging stops at the mantle
    float occluderRadius() const override { return mantleRadius; }
    float seaLevelRadius() const override { return moon || sandbox ? 0.0f : radius; }

    // Rock outcrops and boulders near the surface: metres inside the nearest one (< 0 outside).
    float rockFeature(Vector3 p, float d, float surface, float footprint) const;
    // 3D terrain: metres the surface is raised (+) or carved (-) here by the 3D noise band (0 outside it).
    struct Column3D { float mask = 0.0f, f = 0.0f, wet = 0.0f, river = 0.0f; }; // river (#1438): riverChannel, 1 in a river channel and its bank, fading to 0; // the 3D band's per-column strength (Plan 5 C). (The cliff lip over a column went with the terraces, #543.) // wet (#971): caveWetColumn, the cave roofs' keeper
    // With the height-map elevation and 3D-band column at p already known (gridSurface, from the chunk's block grid):
    // chunk meshing (#148) skips re-evaluating the height map at every vertex (~80% of the per-vertex cost).
    VoxelMaterial materialSkyHinted(Vector3 p, float base, float footprint, Vector3 normal, float& sky, const float* elevHint, const Column3D* colHint) const;
    Column3D terrain3DColumn(Vector3 dir, float elev) const;
    // (#971) 0..1: a wet lowland (swamp, or a river's floodplain), where cave roofs never collapse into trenches (#549's mouths
    // stay closed). Once per height-map grid node (cached with the column) and for single samples in a mouth region only.
    float caveWetColumn(Vector3 dir, float elev) const;
    // (#1438) 1 in a river's channel and 3 m beyond its edge, fading to 0 by 15 m: the density's detail noise and rock features
    // fade out there, so they don't refill the channel over its water (rivers seen as dry sand).
    float riverChannel(Vector3 dir) const;
    // (#615 step 1) All the caves in one place: metres from p to the nearest cave wall (< 0 inside a cave; 1e9 where no cave
    // reaches). densityNear carves where it is under its meshing reach; terrain.glsl's caveDistance is the same, bit-exact.
    // hit (tools: --cave-stats): what made the distance.
    struct CaveHit { int kind = 0; float mouth = 0.0f; bool collapsed = false; }; // kind: 0 none, 1 tunnel, 2 cavern
    float caveDistance(Vector3 p, float d, float depth, float surface, float fp, const Column3D* col, CaveHit* hit = nullptr, float* solid = nullptr) const; // solid (#615 step 3): metres inside rubble aprons (densityNear adds them)
    // (#615 step 2) The cave graph: passages generated per 256 m cell (sea-level projection) from the seed, cached. A block
    // gets the segments whose capsules can reach it, 5 vec4 each: (a, ra), (b, rb), (up, flat), box lo, box hi. False
    // when more than maxSegs (the GPU's cap: the block takes the CPU path). RAYLIBGAME_CAVES=1: graph tunnels instead of
    // the noise tunnels (and their #549 mouths); 0 (the default until natural mouths, step 3) keeps today's caves.
    bool blockCaveSegs(Vector3 origin, float stride, int S, float footprint, std::vector<float>& out, int maxSegs) const;
    static bool CavesGraph();
    // (#615 step 3) The nearest hillside cave mouth to a direction within range (m): its axis point at the slope face and the
    // outward (downhill) direction, body-local. Tests and captures.
    bool nearestCaveMouth(Vector3 dir, float range, Vector3& at, Vector3& outward, int kind = -1) const; // kind: -1 any entrance, 0 hillside, 1 karst sink, 2 resurgence (#1052), 3 a chamber centre (#1058)
    void caveMouthsList(Vector3 dir, float range, int kind, std::vector<Vector3>& out) const; // #1116: the entries of one kind within range (tests)
    struct CaveShaft { Vector3 top, up; float radius, depth; };          // #1067: a shaft open to the sky (karst sink, tube skylight)
    void shaftsNear(Vector3 lo, Vector3 hi, std::vector<CaveShaft>& out) const; // those whose shaft can reach the box [lo, hi]
    int caveMouthsNear(Vector3 dir, float range) const; // (#615 step 3) how many mouths within range (m) of a direction (stats)
    // (#1223) For creature spawning: the cave segments (axis a-b, radii, the up of their floor, kind: 0 passage, 1 mouth, 2 rubble
    // (solid), 3 chamber or lava tube, 4 slot) and the entrances (kind as nearestCaveMouth's) whose axis passes through the box.
    struct CaveSpan { Vector3 a, b, up; float ra, rb; int kind; };
    struct CaveEntrance { Vector3 at, out; int kind; };
    void caveSpansNear(Vector3 lo, Vector3 hi, std::vector<CaveSpan>& out) const;
    void caveEntrancesNear(Vector3 lo, Vector3 hi, std::vector<CaveEntrance>& out) const;
    Vector3 uphillAt(Vector3 dir) const; // the eroded height's gradient (body frame, m per m; zero without the bake): benches face cliffs
    float terrain3D(Vector3 p, float d, float surface, float footprint, const Column3D* col = nullptr) const;

    // baseDensityBlock in two halves (Plan 6 step 4, #147): the height map on a gnomonic grid tangent at the
    // block centre (rivers, lakes, ponds included: CPU), then every sample from it (densityNear: on the CPU
    // here, or on the GPU, terrain_density.comp, with identical results). prepareBlock returns false for
    // blocks that take the per-sample path instead (wholly outside the surface band, or too large).
    struct BlockGrid {
        Vector3 n{}, t1{}, t2{};
        float umin = 0, vmin = 0, su = 1, sv = 1;
        int nu = 0, nv = 0;
        std::vector<float> h;            // elevation per node
        std::vector<Column3D> cols;      // 3D band strength per node
    };
    bool prepareBlock(Vector3 origin, float stride, int S, float footprint, BlockGrid& g) const;
    void sampleBlock(Vector3 origin, float stride, int S, float footprint, const BlockGrid& g, float* out) const;
    // The elevation and 3D-band column at p as sampleBlock interpolates them from the grid (false: p outside the surface band).
    bool gridSurface(const BlockGrid& g, Vector3 p, float& elev, Column3D& col) const;
    // The floating islands (Plan 5 C) whose 4 km cells a block's samples may fall in, for the GPU: per island
    // (cell x, y, z, ok), (dir xyz, waist), (r, depth, 0, 0). False when too many cells (> maxCells).
    bool blockIslands(Vector3 origin, float stride, int S, std::vector<float>& out, int maxCells) const;
    float mantleRadiusValue() const { return mantleRadius; }
    float minElevationValue() const { return minElevation; }
    float maxElevationValue() const { return maxElevation; }
    bool isMoon() const { return moon; }
private:
    struct RiverState { std::once_flag once; std::unique_ptr<RiverNetwork> net; ~RiverState(); };
    // Height-grid nodes (prepareBlock) on the cube faces' global lattices, shared by the blocks that use them (#251):
    // stacked blocks of a column read the same nodes. Sharded; a shard is dropped when full.
    struct GridNodeCache;
    std::shared_ptr<GridNodeCache> gridCache;
    // Erosion bake (#352 prototype): built on first use, around the start (fieldSeed's findSpawnDir).
    struct ErosionState;
    std::shared_ptr<ErosionState> erosionState;
    unsigned int fieldSeed = 0;
    const PlanetErosion* erosion() const;
public:
    // Bake erosion around the start (#352 prototype, gErosionPrototype): once, on the main thread, before the field is
    // shared with the workers.
    void bakeErosion(const std::string& cachePath); // loads the bake from cachePath when it matches, else bakes and saves it
    float profileLandHeight(Vector3 dir, float footprint) const; // land(dir, fp).h, for --profile-density
    float testBroadLand(Vector3 dir, float footprint) const { return broadLand(land(dir, footprint), dir, footprint); } // (#1856 M11 tests) rivers' smooth land
    // Base density (no edits) at p, |p| = d, in a column whose surface elevation and Column3D are known (a block grid's):
    // skips the height map. The sky probe's up ray uses it (#387).
    float densityInColumn(Vector3 p, float d, float elev, float footprint, const Column3D* col) const { return densityNear(p, d, radius + elev, footprint, col); }
private:
    std::shared_ptr<RiverState> riverState;
    unsigned int riverSeed = 0;
    float densityNear(Vector3 p, float d, float surface, float footprint, const Column3D* col = nullptr) const; // after the height map
    struct Land { float h, c, mountain, valley, valleyWeight, detail, stream, ridge; float relief = 0.0f; }; // relief: what the gullies and cliff faces added to the bake (#547: rivers keep to the bake's own floor)
    Land land(Vector3 dir, float footprint) const;     // the height map before rivers and ponds
    struct Pond { bool exists; Vector3 dir; float r, level, depth; };
    Pond pondOfCell(int cx, int cy, int cz) const;       // deterministic per cell (cached per thread)
    float riverWeight(const Land& l) const;
    // A river here: its weight (0: none), water level (smooth: from the land without its fine hills
    // and detail, so it does not follow every bump) and width factor (wider toward the sea).
    struct River { float w, level, wf; };
    River river(const Land& l, Vector3 dir, float footprint) const;
    River stream(const Land& l, Vector3 dir, float footprint) const; // the same for streams (wf unused)
    float broadLand(const Land& l, Vector3 dir, float smoothFootprint) const; // the land without fine hills and detail
};

// Gems (M22): the rare gem a stone kind holds (granite amethyst, marble ruby, basalt sapphire,
// limestone emerald), or nullptr; any rock may also give quartz.
const char* GemForStone(const char* stone);

// No terrain at all: a body that only carries a build plot (a rocket in flight, rocket.h).
class EmptyField : public DensityField {
public:
    explicit EmptyField(float radius) : r(radius) {}
    float baseDensity(Vector3) const override { return ISO_LEVEL - 100.0f; }
    VoxelMaterial material(Vector3, float) const override { return VoxelMaterial::Rock; }
    Vector3 fallbackUp(Vector3) const override { return { 0, 1, 0 }; }
    bool editable(Vector3) const override { return false; }
    float maxRadius() const override { return r; }
private:
    float r;
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

// Water (#360): the generated water a chunk may hold (GeneratedWaterOf, voxel.cpp): whether the node spans the sea's
// radius, and whether inland water may lie in it (its levels: the water surface map, watersurface.h, #639).
struct ChunkWaterGen { bool seaSpan = false, inland = false; };
// What a water-only rebuild of a chunk reuses (kept only for chunks next to simulated water): the terrain density with
// edits (S^3), the generated levels, and the block's height grid.
struct ChunkWaterCache {
    std::vector<float> dens;
    ChunkWaterGen gen;
    std::shared_ptr<const PlanetField::BlockGrid> grid; // shared with the build that made it (not a copy)
    long long counted = 0; // bytes it adds to ChunkWaterCacheBytes
    ChunkWaterCache();
    ~ChunkWaterCache();
    ChunkWaterCache(const ChunkWaterCache&) = delete;
    ChunkWaterCache& operator=(const ChunkWaterCache&) = delete;
};
// Bytes held by all chunk water caches (benches; past kChunkWaterCacheCap no new ones are kept: those chunks redraw
// their water with a full build instead).
long long ChunkWaterCacheBytes();
// Water mesh phase times (ns, summed over workers): sampling, marching cubes, shore clip, vertex output (#381),
void ChunkWaterPhases(long long out[6]); // + exact density samples, all samples
// (#1773) Test probe: water builds on `body` record their triangles (3 body-local corners each) with a corner inside [lo, hi];
// Take hands over (and clears) what they recorded since. Scenarios check the drawn water mesh with it.
void WaterTriangleProbe_Set(uint64_t body, Vector3 lo, Vector3 hi);
void WaterTriangleProbe_Take(std::vector<Vector3>& out);
void ChunkWaterFlowStats(long long out[3]); // moving simulated water vertices: all, steep, steep drawn running down their surface (#677)
constexpr long long kChunkWaterCacheCap = 128ll << 20;

// CPU-side result of meshing one chunk (produced on a worker thread).
// The terrain arena's packed vertices (#302, #480), filled on the mesh worker (TerrainArena_Pack): the attributes, the
// chunk's material palette, and, when triangles on material borders were split, its own positions and indices.
struct ArenaPacked {
    std::vector<unsigned char> attr;          // terrainarena.cpp Attr per vertex
    std::vector<float> positions;             // empty: the mesh's own (no triangle was split)
    std::vector<unsigned short> indices;
    unsigned char palette[16] = {};           // palette slot -> VoxelMaterial id
    bool valid = false;
};

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
    std::vector<float> occlusion;          // per vertex: 0 open .. 1 enclosed, short-range ambient occlusion (#132; near LODs)
    std::vector<float> waterPositions;     // sea surface where it meets open air (same frame as positions)
    std::vector<float> waterNormals;
    std::vector<float> waterDepth;         // per water vertex: (depth to the floor in m, fluid kind: fluids.txt)
    std::vector<unsigned char> waterFlow;  // per water vertex (#361): the current, body-local (rgb = dir * 0.5 + 0.5, a = speed / 4 m/s)
    std::vector<unsigned short> waterIndices;
    std::vector<float> grass;              // grass tufts (foliage.h: kGrassFloats each), finest level only
    std::vector<unsigned short> indices;
    std::vector<unsigned char> matIds;     // #480: per vertex, its VoxelMaterial id (0-255; the arena packs a per-chunk palette of up to 16)
    ArenaPacked arena;                     // #302, #480: the terrain arena's packed vertices, filled on the worker
    bool overflow = false;                 // > 65535 vertices: mesh dropped
    float buildMs = 0.0f;
    float sampleMs = 0.0f;                 // part of buildMs spent sampling the density field
    float priority = 0.0f;                 // as requested: < 0 edits (urgent), else distance / (1 + level) (#258 upload order)
    // Water (#360): the snapshot the water was meshed from (newer water wins), a water-only rebuild, and the cache
    // a water-only rebuild needs (set when the chunk touches simulated water).
    uint64_t waterSeq = 0;
    bool waterOnly = false;
    std::shared_ptr<const ChunkWaterCache> waterCache;
};

// Samples the field over node `key` (CHUNK_SIZE+1 cells at the node stride, i.e. one cell of
// overlap with the +x/+y/+z neighbours to hide cracks between LOD levels) and runs marching cubes.
// Thread-safe: touches only its arguments.
void BuildChunkMesh(const DensityField& field, const TerrainDeltas& deltas, ChunkMeshData& out);
// The density block a chunk build samples (E^3 at origin + i * stride, the given footprint), and the build's
// second half from that block (edits, marching cubes, materials): Plan 6 step 4 samples the block on the GPU.
struct ChunkBlock { Vector3 origin{}; float stride = 1.0f, footprint = 1.0f; int E = 0; };
ChunkBlock ChunkBlockOf(const NodeKey& key);
// grid: the block's height grid when it was sampled from one (PlanetField::prepareBlock): materials read the surface from it.
// gridShared: the same grid when the caller holds it in a shared_ptr (the water cache keeps it instead of a copy).
void BuildChunkMeshFromBase(const DensityField& field, const TerrainDeltas& deltas, ChunkMeshData& out, const float* base, float sampledMs, const PlanetField::BlockGrid* grid = nullptr,
                            std::shared_ptr<const PlanetField::BlockGrid> gridShared = nullptr);
// The water surface alone from a chunk's water cache and the latest published water (WaterField::Current): flowing
// water redraws without rebuilding the terrain (#360). out.key / planetoidId set by the caller.
void BuildChunkWater(const DensityField& field, const ChunkWaterCache& cache, ChunkMeshData& out);


#endif
