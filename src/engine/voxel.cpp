#include <functional>
#include <unordered_map>
#include <filesystem>
#include <cstring>
#include "voxel.h"
#include "envflag.h"
#include "utils/logger.h"
#include "erosion.h"
#include "canyons.h"
#include <thread>
#include "rivers.h"
#include "foliage.h"
#include "terrainarena.h"
#include "trees.h"
#include "cubemarch.h"
#include "waterfield.h"
#include "watersurface.h"

#include <cstring>
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>
#include <cstdlib>

const char* VoxelMaterialName(VoxelMaterial m) {
    switch (m) {
        case VoxelMaterial::Rock: return "rock";
        case VoxelMaterial::Regolith: return "regolith";
        case VoxelMaterial::Ice: return "ice";
        case VoxelMaterial::Ore: return "ore";
        case VoxelMaterial::Sand: return "sand";
        case VoxelMaterial::Gravel: return "gravel";
        case VoxelMaterial::Grass: return "grass";
        case VoxelMaterial::Snow: return "snow";
        case VoxelMaterial::Clay: return "clay";
        case VoxelMaterial::Granite: return "granite";
        case VoxelMaterial::Basalt: return "basalt";
        case VoxelMaterial::Sandstone: return "sandstone";
        case VoxelMaterial::Limestone: return "limestone";
        case VoxelMaterial::Marble: return "marble";
        case VoxelMaterial::Pebbles: return "pebbles";
        case VoxelMaterial::Silt: return "silt";
        case VoxelMaterial::Mud: return "mud";
        case VoxelMaterial::Farmland: return "farmland";
        default: return "?";
    }
}

// ---------------------------------------------------------------------------------------------
// PlanetoidField
// ---------------------------------------------------------------------------------------------

static SimplexNoise MakeShapeNoise(std::mt19937& rng, float radius) {
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    const float frequencyJitter = unit(rng) - 0.5f; // -0.5..0.5, same range as the original generator
    const float sx = unit(rng) * 10000.0f, sy = unit(rng) * 10000.0f, sz = unit(rng) * 10000.0f;
    return SimplexNoise((1.5f + frequencyJitter) / radius, 1.0f, 2.0f, 0.5f, sx, sy, sz);
}

static SimplexNoise MakeDetailNoise(std::mt19937& rng) {
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    return SimplexNoise(0.045f, 1.0f, 2.0f, 0.5f, unit(rng) * 10000.0f, unit(rng) * 10000.0f, unit(rng) * 10000.0f);
}

PlanetoidField::PlanetoidField(unsigned int seed, float radius_)
    : PlanetoidField(std::mt19937(seed), radius_) {}

PlanetoidField::PlanetoidField(std::mt19937 rng, float radius_)
    : shape(MakeShapeNoise(rng, radius_)), detail(MakeDetailNoise(rng)), radius(radius_) {
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    iceAmount = (unit(rng) < 0.4f) ? 0.15f + 0.45f * unit(rng) : 0.0f; // ~40% of bodies have ice
}

float PlanetoidField::baseDensity(Vector3 p) const {
    // Displaced sphere: density rises DENSITY_GRADIENT per unit toward the centre, and 3D noise
    // pushes the surface in/out by up to ~25% of the radius (3D, so overhangs are possible).
    // A uniform gradient keeps surfaces crisp and makes dig strokes carve predictable depths;
    // the earlier (noise+1)*(1-d/R) field ramped ~0.02/unit and left mushy near-iso regions.
    const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    const float displacement = radius * 0.25f * shape.fractal(4, p.x, p.y, p.z); // sampled in local space
    return ISO_LEVEL + DENSITY_GRADIENT * (radius * MEAN_SURFACE + displacement - d);
}

void DensityField::baseDensityBlock(Vector3 o, float stride, int S, float fp, float* out) const {
    size_t i = 0;
    for (int z = 0; z < S; ++z)
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x, ++i) out[i] = baseDensityLod({ o.x + x * stride, o.y + y * stride, o.z + z * stride }, fp);
}

float DensityField::density(Vector3 p, const TerrainDeltas* deltas, float footprint) const {
    const float base = footprint > 0.0f ? baseDensityLod(p, footprint) : baseDensity(p);
    return deltas ? base + deltas->sample(p) : base; // edits never reach a planet's mantle (folding skips it)
}

Vector3 DensityField::normal(Vector3 p, const TerrainDeltas* deltas, float h, float fp) const {
    const Vector3 g = {
        density({ p.x + h, p.y, p.z }, deltas, fp) - density({ p.x - h, p.y, p.z }, deltas, fp),
        density({ p.x, p.y + h, p.z }, deltas, fp) - density({ p.x, p.y - h, p.z }, deltas, fp),
        density({ p.x, p.y, p.z + h }, deltas, fp) - density({ p.x, p.y, p.z - h }, deltas, fp),
    };
    const float len = sqrtf(g.x * g.x + g.y * g.y + g.z * g.z);
    if (len < 1e-6f) {
        return fallbackUp(p);
    }
    return { -g.x / len, -g.y / len, -g.z / len };
}

// `base` is the unedited density at p: near ISO_LEVEL means "original surface", higher means
// deeper inside (exposed by digging), lower means material that was added by the player.
VoxelMaterial PlanetoidField::material(Vector3 p, float base) const {
    // Thin ore sheets: where the detail noise crosses zero.
    if (fabsf(detail.fractal(2, p.x, p.y, p.z)) < 0.035f) return VoxelMaterial::Ore;
    if (iceAmount > 0.0f && detail.fractal(2, p.x + 1000.0f, p.y, p.z) > 0.55f - iceAmount) return VoxelMaterial::Ice;
    if (base < ISO_LEVEL + 1.0f * DENSITY_GRADIENT && detail.fractal(1, p.x + 2000.0f, p.y, p.z) > -0.25f) return VoxelMaterial::Regolith;
    return VoxelMaterial::Rock;
}

Vector3 PlanetoidField::fallbackUp(Vector3 p) const {
    const float pl = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    return pl > 1e-6f ? Vector3{ p.x / pl, p.y / pl, p.z / pl } : Vector3{ 0, 1, 0 };
}

// ---------------------------------------------------------------------------------------------
// TerrainField (open world, +Y up)
// ---------------------------------------------------------------------------------------------

// #532: the seed offsets were drawn in [0, 10000) METRES, so a field with a long wavelength moved by only offset x
// frequency noise units between seeds (continents on a 60 km planet: at most 0.42 units, while the planet spans 5):
// every seed sampled nearly the same patch of the noise, and worlds looked alike. Now the span is max(10000 m, U /
// frequency): at least U = 64 noise units for every field. RAYLIBGAME_SEED_SPAN=U overrides it (0: the old 10 km).
float gSeedSpanUnits = -1.0f; // < 0: read the environment once
static float SeedSpan(float frequency) {
    if (gSeedSpanUnits < 0.0f) { const char* e = getenv("RAYLIBGAME_SEED_SPAN"); gSeedSpanUnits = e ? (float)atof(e) : 64.0f; }
    return gSeedSpanUnits > 0.0f ? std::max(10000.0f, gSeedSpanUnits / frequency) : 10000.0f;
}

float NoiseSeedSpanUnits() { SeedSpan(1.0f); return gSeedSpanUnits; }

// #1055: a process-wide counter for PlanetField::caveEpoch (never 0, so a default key never matches a field's).
static uint32_t NextCaveEpoch() { static std::atomic<uint32_t> n{ 0 }; return ++n; }

static SimplexNoise SeededNoise(std::mt19937& rng, float frequency) {
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    const float span = SeedSpan(frequency);
    return SimplexNoise(frequency, 1.0f, 2.0f, 0.5f, unit(rng) * span, unit(rng) * span, unit(rng) * span);
}

TerrainField::TerrainField(unsigned int seed) {
    std::mt19937 rng(seed);
    hills = SeededNoise(rng, 0.0025f);  // ~400 m features: valleys and mountains
    detail = SeededNoise(rng, 0.02f);   // ~50 m lumps
    caves = SeededNoise(rng, 0.012f);   // tunnel network
    veins = SeededNoise(rng, 0.05f);    // crystal veins
}

float TerrainField::surfaceHeight(float x, float z) const {
    // (2D fractal in this noise library ignores the seed, so sample the 3D one on a plane.)
    const float h = hills.fractal(5, x, 0.0f, z);         // -1..1
    const float mountains = std::max(0.0f, h) * std::max(0.0f, h) * 140.0f; // sharpen peaks
    return 6.0f + h * 22.0f + mountains + detail.fractal(3, x, 100.0f, z) * 4.0f;
}

float TerrainField::baseDensity(Vector3 p) const {
    const float h = surfaceHeight(p.x, p.z);
    float d = ISO_LEVEL + DENSITY_GRADIENT * (h - p.y);
    // A little 3D variation for overhangs and crags.
    d += DENSITY_GRADIENT * 2.5f * detail.fractal(2, p.x, p.y * 1.5f, p.z);
    // Caves: tunnels where the cave noise crosses zero, kept a few metres below the surface.
    const float depth = h - p.y;
    if (depth > 4.0f) {
        const float c = fabsf(caves.fractal(3, p.x, p.y * 1.6f, p.z));
        const float tunnel = 1.0f - c / 0.07f;
        if (tunnel > 0.0f) d -= tunnel * std::min(1.0f, (depth - 4.0f) / 6.0f) * 3.0f;
    }
    return d;
}

VoxelMaterial TerrainField::material(Vector3 p, float base) const {
    const float depth = (base - ISO_LEVEL) / DENSITY_GRADIENT; // ~metres below the original surface
    if (depth > 6.0f && fabsf(veins.fractal(2, p.x, p.y, p.z)) < 0.03f) return VoxelMaterial::Ore; // crystal
    if (depth < 2.5f) {
        if (p.y > snowLine + detail.fractal(1, p.x, 300.0f, p.z) * 12.0f) return VoxelMaterial::Ice;  // snow
        return VoxelMaterial::Regolith;                                                               // grass / soil
    }
    return VoxelMaterial::Rock;
}

// ---------------------------------------------------------------------------------------------
// PlanetField (planets and moons)
// ---------------------------------------------------------------------------------------------

constexpr float kBand3DBelow = 16.0f, kBand3DAbove = 50.0f; // the 3D terrain band around the height map (Plan 5 C: up to 50 m up)
static float Smoothstep(float a, float b, float x) {
    const float t = std::clamp((x - a) / (b - a), 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}
// A sine-like wave from arithmetic only (two parabolas per period, within ~5% of sin): the density field avoids
// sinf / powf so its GPU port (resources/compute/terrain.glsl) computes exactly the same.
static float WaveSin(float x) {
    const float u = x * 0.15915494f, f = u - floorf(u);
    return f < 0.5f ? 16.0f * f * (0.5f - f) : -16.0f * (f - 0.5f) * (1.0f - f);
}

struct PlanetField::GridNodeCache {
    struct Key {
        long long i, j; uint32_t fp; int face;
        bool operator==(const Key& o) const { return i == o.i && j == o.j && fp == o.fp && face == o.face; }
    };
    struct Hash {
        size_t operator()(const Key& k) const {
            uint64_t h = (uint64_t)k.i * 0x9E3779B97F4A7C15ull ^ (uint64_t)k.j * 0xC2B2AE3D27D4EB4Full ^ ((uint64_t)k.fp << 3 | (uint64_t)k.face) * 0x165667B19E3779F9ull;
            h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ull; h ^= h >> 29;
            return (size_t)h;
        }
    };
    struct Node { float h; Column3D col; };
    static constexpr int kShards = 64;
    static constexpr size_t kShardCap = 1u << 14; // ~1M nodes in all (~60 MB)
    struct Shard { std::mutex m; std::unordered_map<Key, Node, Hash> map; };
    Shard shards[kShards];
    Shard& shardOf(const Key& k) { return shards[Hash()(k) >> 58]; }
};

PlanetField::PlanetField(unsigned int seed, float radius_, bool moon_)
    : radius(radius_), mantleRadius(radius_ * (1.0f - 0.015f)), moon(moon_) {
    oreSeed = seed ^ 0x0e5e0e5u;
    riverSeed = seed;
    riverState = std::make_shared<RiverState>();
    gridCache = std::make_shared<GridNodeCache>();
    erosionState = std::make_shared<ErosionState>();
    fieldSeed = seed;
    caveEpoch = NextCaveEpoch();      // #1055
    std::mt19937 rng(seed);
    // Large features scale with the body (a moon gets proportionally sized basins and ridges);
    // hills, detail, caves and veins are absolute sizes, so they feel the same everywhere.
    continents = SeededNoise(rng, 2.5f / radius);  // ~24 km continents on a 60 km planet
    ridges = SeededNoise(rng, 8.5f / radius);      // mountain ranges
    hills = SeededNoise(rng, 1.0f / 900.0f);
    valleys = SeededNoise(rng, 6.5f / radius);     // river valleys: where this crosses zero
    detail = SeededNoise(rng, 1.0f / 70.0f);
    cavesA = SeededNoise(rng, 1.0f / 150.0f);      // tunnels run where both cave fields cross zero
    cavesB = SeededNoise(rng, 1.0f / 150.0f);
    veins = SeededNoise(rng, 1.0f / 30.0f);
    climate = SeededNoise(rng, 1.0f / 15000.0f);
    // Added later: drawn after the others so existing seeds keep their continents and mountains.
    caverns = SeededNoise(rng, 1.0f / 220.0f);       // rare underground chambers
    ravines = SeededNoise(rng, 5.0f / radius);       // narrow canyons where this crosses zero
    ravineMask = SeededNoise(rng, 3.0f / radius);    // ...only in some regions
    streams = SeededNoise(rng, 1.0f / 900.0f);        // small streams where this crosses zero (added last)
    rocks = SeededNoise(rng, 1.0f / 22.0f);           // outcrop fins and crags (playtest 4)
    rockMask = SeededNoise(rng, 1.0f / 340.0f);       // ...only in rocky patches
    rugged = SeededNoise(rng, 1.0f / 1600.0f);        // 3D terrain (plan item 6): where hills turn to cliffs
    cliffs = SeededNoise(rng, 1.0f / 45.0f);          // ...their shape (ridged, stretched upward)
    warp3 = SeededNoise(rng, 1.0f / 70.0f);           // ...and its domain warp (leaning crests, undercuts)
    tempNoise = SeededNoise(rng, 1.0f / 4200.0f);     // biomes (Plan 5 B): regional warm and cool patches
    moistNoise = SeededNoise(rng, 1.0f / 3200.0f);    // ...wet and dry ones
    biomeJitter = SeededNoise(rng, 1.0f / 90.0f);     // ...and a fine wobble that interleaves their borders
    caveMouths = SeededNoise(rng, 1.0f / 700.0f);     // Plan 5 C: regions where cave tunnels break the surface
    // Bounds of elevation() below (each noise term is in [-1, 1]).
    if (moon) { minElevation = -90.0f; maxElevation = 340.0f; }
    else { minElevation = -740.0f; maxElevation = 1690.0f + 50.0f + 15.0f; } // + the 3D terrain band + biome landforms (dunes, mesas)
}

void PlanetField::setGenParams(const WorldGenParams& p) {
    caveEpoch = NextCaveEpoch();      // #1055: the land changes: cave cells built before must not be reused
    gen = p;
    if (moon) return;
    // land(): 200 m of gentle land + the ranges (1400 m x mountains) + hills (85 m x roughness) + detail; the default
    // parameters give the bounds above exactly.
    if (p.isDefault()) { maxElevation = 1690.0f + 50.0f + 15.0f; return; }
    maxElevation = 200.0f + 1400.0f * std::max(p.mountains, 0.0f) + 90.0f * std::max(p.roughness, 1.0f) + 50.0f + 15.0f;
}

void PlanetField::setSandbox() {
    caveEpoch = NextCaveEpoch();      // #1055
    if (moon) return;
    sandbox = true;
    gen.erosion = false;               // no bake
    minElevation = kSandboxHeight - 2.0f; maxElevation = kSandboxHeight + 2.0f;
}

// River valleys: 70 m deep, with a flat floodplain floor (valley-noise units < kValleyFloor, room for
// the widest river) and sides rising over the next 0.04 units. (A V-shaped floor rose ~27 m across a
// river's width and tilted its water.)
static constexpr float kValleyFloor = 0.035f, kValleyEdge = 0.075f;
static float ValleyCarve(float v) {
    if (v >= kValleyEdge) return 0.0f;
    const float t = 1.0f - std::max(v - kValleyFloor, 0.0f) / (kValleyEdge - kValleyFloor);
    return t * t * 70.0f;
}


struct PlanetField::ErosionState {
    std::atomic<bool> active{ false }; PlanetErosion bake; std::once_flag pondOnce; std::unordered_map<uint64_t, int> pondByCell; // pondByCell: 260 m pond cell -> the bake's pond there (#360)
    CanyonGrade grade; bool gradeOn = false; // #544: the canyons' graded floors, while the bake reads its input (set before its threads start)
};

bool PlanetField::caveBaked() const { return erosion() != nullptr; }
const PlanetErosion* PlanetField::erosion() const {
    if (tIgnoreBake) return nullptr;
    return erosionState && erosionState->active.load(std::memory_order_acquire) ? &erosionState->bake : nullptr;
}

void PlanetField::bakeErosion(const std::string& cachePath) {
    if (moon || !erosionState) return;
    // The whole planet (Plan 7 step 1, #355), on the raw land; then switched on. Loaded from cachePath when a bake of
    // this seed and radius is there (a world load reads it instead of baking again). Any river network built before
    // came from the raw land: start a fresh one, built on first use from the eroded land.
    PlanetErosion& b = erosionState->bake;
    ErosionParams p;
    auto envf = [](const char* n, float v) { const char* e = getenv(n); return e ? (float)atof(e) : v; }; // tuning
    p.k = envf("RAYLIBGAME_EROSION_K", p.k); p.steps = (int)envf("RAYLIBGAME_EROSION_STEPS", (float)p.steps);
    p.cell = envf("RAYLIBGAME_EROSION_CELL", p.cell); p.uplift = envf("RAYLIBGAME_EROSION_UPLIFT", p.uplift);
    p.start = envf("RAYLIBGAME_EROSION_START", p.start);
    p.glacial = envf("RAYLIBGAME_GLACIAL", p.glacial); p.glacialLine = envf("RAYLIBGAME_GLACIAL_LINE", p.glacialLine);
    p.wind = envf("RAYLIBGAME_WIND", p.wind); // #370: aeolian erosion (0: off, A/B)
    p.refine = (int)envf("RAYLIBGAME_EROSION_REFINE", (float)p.refine); p.refineSteps = (int)envf("RAYLIBGAME_EROSION_REFINE_STEPS", (float)p.refineSteps);
    if (!gen.isDefault()) { p.k *= gen.erosionStrength; p.uplift *= gen.mountains; p.landKey = gen.key(); } // #504: the world's parameters (defaults: the key and bake as before)
    // The cache must match the seed, radius, parameters and bake code version (#355 review): tuning runs and code
    // changes bake again instead of loading a stale bake.
    const auto tl0 = std::chrono::steady_clock::now();
    const bool loaded = !cachePath.empty() && b.load(cachePath, radius, fieldSeed, ErosionParamsKey(p));
    const float loadMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - tl0).count();
    if (!loaded) {
        int threads = std::max(1u, std::thread::hardware_concurrency());
        if (const char* e = getenv("RAYLIBGAME_BAKE_THREADS")) threads = std::max(1, atoi(e)); // #605: the determinism check bakes with other counts
        // Canyon floors graded to drain end to end (#544, #547) before the bake reads its input; RAYLIBGAME_CANYON_GRADE=0: off (A/B).
        static const bool grade = !(getenv("RAYLIBGAME_CANYON_GRADE") && getenv("RAYLIBGAME_CANYON_GRADE")[0] == '0');
        if (grade) {
            CanyonGrade& g = erosionState->grade;
            g.build(radius, 100.0f, [this](Vector3 d) { return canyonFault(d.x * radius, d.y * radius, d.z * radius); },
                    [this](Vector3 d) { tBakeInput = true; const float h = land(d, 50.0f).h; tBakeInput = false; return h; }, threads);
            erosionState->gradeOn = true;
            LOG_INFO("Erosion", "canyons graded (#544): %d canyons (%d drain both ways from a divide), %.0f km, %d corridor cells of 100 m; floors cut %.1f m on average (up to %.0f m), filled %.1f m (up to %.0f m); %.2f s\n",
                     g.canyons, g.divided, g.lengthKm, g.corridorCells, g.meanCut, g.maxCut, g.meanFill, g.maxFill, g.seconds);
        }
        b.build([this](Vector3 d, float fp) { tBakeInput = true; const float h = land(d, fp).h; tBakeInput = false; return h; }, radius, fieldSeed, p, threads, [this](Vector3 d) { return strataColumn(d); }); // + the rock strata (#587)
        erosionState->gradeOn = false;
        erosionState->grade = CanyonGrade(); // only the bake's input needs it
        const auto ts0 = std::chrono::steady_clock::now();
        if (!cachePath.empty()) {
            // The save's folder may not exist yet (a first bench run): the cache silently went unwritten and every run baked again.
            std::error_code ec;
            const std::filesystem::path dirPath = std::filesystem::path(cachePath).parent_path();
            if (!dirPath.empty()) std::filesystem::create_directories(dirPath, ec);
            if (!b.save(cachePath)) LOG_WARN("Erosion", "could not write the bake cache %s\n", cachePath.c_str());
        }
        LOG_INFO("Erosion", "cache saved in %.1f s\n", std::chrono::duration<float>(std::chrono::steady_clock::now() - ts0).count());
        LOG_INFO("Erosion", "baked the planet, %d cells at %.0f m + %d at %.0f m, in %.1f s (coarse %.1f s, faces %.1f s): cut up to %.0f m, raised up to %.0f m, hash %016llx\n",
               b.cells(), p.cell, b.fineCells(), b.fineCell(), b.buildSeconds, b.coarseSeconds, b.refineSeconds, b.maxCut, b.maxRaise, (unsigned long long)b.hash());
        LOG_INFO("Erosion", "bake profile (#605): %s\n", b.profile.c_str());
    } else {
        LOG_INFO("Erosion", "loaded %s (%d cells, hash %016llx) in %.0f ms\n", cachePath.c_str(), b.cells(), (unsigned long long)b.hash(), loadMs);
    }
    LOG_INFO("Erosion", "%zu ponds from the bake (#360)\n", b.ponds().size());
    erosionState->active.store(true, std::memory_order_release);
    riverState = std::make_shared<RiverState>();
}
float PlanetField::profileLandHeight(Vector3 dir, float fp) const { return land(dir, fp).h; }
Vector3 PlanetField::uphillAt(Vector3 dir) const {
    const PlanetErosion* eroded = erosion();
    return eroded ? eroded->slopeAt(dir, 1.0f).grad : Vector3{ 0, 0, 0 };
}
PlanetField::Land PlanetField::land(Vector3 dir, float fp) const {
    const float sx = dir.x * radius, sy = dir.y * radius, sz = dir.z * radius; // the surface point (m)
    Land l{ 0, 0, 0, 1, 0, 0, 1, 0 };
    if (moon) {
        // Grey highlands and lowland plains, rough ridged uplands, no oceans.
        const float base = continents.fractalLimited(4, fp, sx, sy, sz) * 60.0f;
        const float r = 1.0f - fabsf(ridges.fractalLimited(4, fp, sx, sy, sz));
        l.h = base + r * r * 250.0f + hills.fractalLimited(3, fp, sx, sy, sz) * 25.0f + detail.fractalLimited(2, fp, sx, sy, sz) * 2.5f;
        return l;
    }
    const float c = continents.fractalLimited(4, fp, sx, sy, sz) - 0.15f * gen.ocean; // < 0: ocean, > 0: land (#504: ocean shifts it)
    const float base = c < 0.0f ? c * 650.0f : c * 200.0f;        // ocean floor down to -650 m, gentle land
    const float landMask = Smoothstep(-0.02f, 0.06f, c);
    const float mountainMask = Smoothstep(0.08f, 0.35f, c);      // ranges only well inland
    // With the erosion bake (Plan 7, #356) the baked height replaces the noise relief (base, ranges, hills), and gullies
    // steered by its slope add the 5-100 m detail (#357). l.ridge then comes from the height (the rivers' and ponds'
    // mountain mask, until terrain-first water replaces them).
    const PlanetErosion* eroded = erosion();
    float h;
    if (eroded) {
        const PlanetErosion::Slope s = eroded->slopeAt(dir, fp);
        h = s.h;
        if (gGullies && fp < 450.0f) h += (l.relief = ErosionGullies({ sx, sy, sz }, dir, s.grad, fp, fieldSeed, Smoothstep(2.0f, 15.0f, s.h)));
        { const float cf = cliffFace({ sx, sy, sz }, h, s.grad, fp); h += cf; l.relief += cf; } // one face per steep resistant layer (#588)
        { // inside the elevation bound (less its 65 m allowance): the bake and the faces can lift formations (#593)
            const float top = maxElevation - 65.0f, knee = top - 200.0f;
            if (h > knee) h = knee + 200.0f * (1.0f - expf(-(h - knee) / 200.0f));
        }
        l.ridge = std::max(0.0f, h - 150.0f);
    } else {
        h = base;
        if (mountainMask > 0.0f) {
            const float r = 1.0f - fabsf(ridges.fractalLimited(5, fp, sx, sy, sz)); // ridged: sharp crests
            l.ridge = r * r * r * 1400.0f * mountainMask * gen.mountains; // the height the ranges add here
            h += l.ridge;
        }
        h += hills.fractalLimited(3, fp, sx, sy, sz) * (25.0f + 60.0f * mountainMask) * gen.roughness;
    }
    // (River valleys are carved by the drainage network's rivers in elevation(); the old valley-noise
    // valleys and streams are gone, so l.valley and l.stream stay at 1: none.)
    const float valleyWeight = landMask * std::clamp(2.0f - fp / 64.0f, 0.0f, 1.0f);
    l.valley = 1.0f;
    // Ravines: rare narrow canyons (~60 m wide, up to ~90 m deep, steep sqrt-profile walls) in some
    // land regions, faded out where samples are too far apart to resolve them.
    // With the bake (#544, Plan 9) they are canyons in its input instead: a ~200 m notch along the same lines (wide enough
    // for the 50 m grid) on soft rock (strataColumn), which the rivers find and cut, so they come out eroded and wet.
    // The baked terrain gets no post-bake cut.
    if (tBakeInput) {
        const float notch = canyonFault(sx, sy, sz);
        if (notch > 0.0f) h -= notch * 60.0f * landMask;
        // Fault scarps (#590): a few long faults lift one side 50-300 m (scaled by "Cliffs and canyons"); the bake wears
        // the step into an escarpment with spurs and gullies (hard rock on the raised side: Strata_Column).
        h += Strata_FaultUplift({ sx, sy, sz }, fieldSeed) * landMask * std::min(gen.cliffs, 1.5f);
        // Cuestas and hogbacks (#591): in fold regions the resistant layers stand up to 80 m proud where they crop out
        // (the differential erosion the bake's short run can't do alone; it rounds and gullies them).
        // Karst (#592): limestone towers and dolines in karst regions (New World "Karst"); the bake rounds them and the
        // cliff profile (#588) gives the hard limestone its faces.
        // Volcanoes (#593): cones 1-3 km high with craters or calderas, on land or rising from the sea as islands; their
        // basalt flows over soft ash (strataColumn) erode into radial gullies and lava-capped mesas.
        if (gen.volcanoes) h += Strata_Volcano({ sx, sy, sz }, radius, fieldSeed);
        if (gen.karst && h > 0.0f)
            if (const float kw = Strata_KarstWeight({ sx, sy, sz }, fieldSeed); kw > 0.0f)
                h += kw * landMask * Strata_KarstRelief({ sx, sy, sz }, fieldSeed);
        if (const float fw = Strata_FoldWeight({ sx, sy, sz }, fieldSeed); fw > 0.0f && h > 0.0f)
            h += 80.0f * fw * landMask * std::min(gen.cliffs, 1.5f) * Strata_OutcropBump(strataColumn(dir), h, fieldSeed);
        // Graded canyon floors (#544, #547; canyons.h): along each canyon line the floor falls toward its outlet, cutting
        // through the ridges it crosses and filling its dips, so the bake's rivers run along it end to end. The walls
        // widen with the depth of the cut (soft rock holds ~35 deg), and past the corridors the plain notch stays.
        if (erosionState && erosionState->gradeOn && landMask > 0.0f) {
            float cov = 0.0f;
            const float floorH = erosionState->grade.floorAt(dir, cov);
            if (cov > 0.0f) {
                const float hl = h + notch * 60.0f * landMask; // the land without the notch
                const float D = std::max(0.0f, hl - floorH);
                const float r0 = 0.02f * std::max(1.0f, D / 75.0f);
                const float t = std::clamp(1.0f - canyonLine(sx, sy, sz) / r0, 0.0f, 1.0f);
                const float hg = hl - D * t * t * (3.0f - 2.0f * t);
                h += (hg - h) * cov * landMask;
            }
        }
        // Keep the formations inside the field's elevation bound (maxElevation less its 65 m for detail and landforms):
        // a volcano or a fault on high ground compresses smoothly into the top 200 m instead of poking above the band
        // the chunk selection trusts (#593: seed 42 reached 1759 m against 1755).
        const float top = maxElevation - 65.0f, knee = top - 200.0f;
        if (h > knee) h = knee + 200.0f * (1.0f - expf(-(h - knee) / 200.0f));
    }
    const float ravineWeight = eroded || tBakeInput ? 0.0f : landMask * std::clamp(2.0f - fp / 24.0f, 0.0f, 1.0f);
    if (ravineWeight > 0.0f) {
        const float region = Smoothstep(0.25f, 0.4f, ravineMask.fractal(2, sx, sy, sz));
        if (region > 0.0f) {
            const float r = fabsf(ravines.fractal(2, sx, sy, sz));
            if (r < 0.006f) h -= sqrtf(1.0f - r / 0.006f) * 90.0f * region * ravineWeight;
        }
    }
    l.detail = detail.fractalLimited(2, fp, sx, sy, sz) * 2.5f * gen.roughness;
    // Polar ice shelves: near the poles the sea is frozen, so low ground is raised to an ice sheet a
    // couple of metres above sea level (the material rules make it ice).
    const float iceCap = Smoothstep(0.86f, 0.92f, fabsf(dir.y) + 0.03f * climate.fractalLimited(2, fp, sx, sy, sz));
    float hd = h + l.detail;
    if (iceCap > 0.0f && hd < 2.0f) { hd += (2.0f - hd) * iceCap; h = hd - l.detail; }
    l.h = h; // without the fine detail (elevation() adds it)
    l.c = c;
    l.mountain = mountainMask;
    l.valleyWeight = valleyWeight * (1.0f - iceCap);
    return l;
}

namespace {
    constexpr float kRiverHalf = 0.011f;   // channel half-width in valley-noise units (~15 m)
    constexpr float kRiverZone = 0.022f;   // where the banks are smoothed (detail damped)
    constexpr float kRiverDepth = 3.2f;    // channel depth at its centre (m)
    constexpr float kRiverDrop = 0.9f;     // water surface below the banks (m)
    constexpr float kPondCell = 260.0f;    // m
    constexpr float kStreamHalf = 0.029f;  // stream half-width in stream-noise units (~4 m)
    constexpr float kStreamZone = 0.06f;
    constexpr float kStreamDepth = 1.6f;
    constexpr float kStreamDrop = 0.5f;
    uint32_t PondHash(int x, int y, int z, uint32_t seed) {
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
        h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
        return h;
    }
    float PondUnit(uint32_t& h) { h = h * 1664525u + 1013904223u; return (h >> 8) * (1.0f / 16777216.0f); }
    uint64_t PondId(int x, int y, int z) { return ((uint64_t)(x & 0x1FFFFF) << 42) | ((uint64_t)(y & 0x1FFFFF) << 21) | (uint64_t)(z & 0x1FFFFF); }
    constexpr float kPondReach = 70.0f;    // a pond's effect (embankment) ends within this of its centre
    // The pond cells whose pond could reach within `reach` of surface point s (usually just one).
    template <class F> void ForPondCellsNear(Vector3 s, float reach, F&& f) {
        const int x0 = (int)floorf((s.x - reach) / kPondCell), x1 = (int)floorf((s.x + reach) / kPondCell);
        const int y0 = (int)floorf((s.y - reach) / kPondCell), y1 = (int)floorf((s.y + reach) / kPondCell);
        const int z0 = (int)floorf((s.z - reach) / kPondCell), z1 = (int)floorf((s.z + reach) / kPondCell);
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) f(x, y, z);
    }
}

float PlanetField::riverWeight(const Land& l) const {
    // Rivers on land valleys, not in the mountains, above the sea (they meet it at the coast). With the bake the
    // network follows its drainage, which starts high in the ranges (#359): no mountain or altitude cut.
    if (erosion()) return l.valleyWeight * Smoothstep(1.0f, 6.0f, l.h);
    return l.valleyWeight * (1.0f - Smoothstep(80.0f, 250.0f, l.ridge)) * Smoothstep(1.0f, 6.0f, l.h) * (1.0f - Smoothstep(700.0f, 900.0f, l.h));
}

float PlanetField::broadLand(const Land& l, Vector3 dir, float smoothFp) const {
    if (const PlanetErosion* e = erosion()) return e->heightAt(dir, smoothFp);
    const float sx = dir.x * radius, sy = dir.y * radius, sz = dir.z * radius;
    const float base = l.c < 0.0f ? l.c * 650.0f : l.c * 200.0f;
    return base + l.ridge + hills.fractalLimited(3, smoothFp, sx, sy, sz) * (25.0f + 60.0f * l.mountain) - ValleyCarve(l.valley) * l.valleyWeight;
}

PlanetField::River PlanetField::stream(const Land& l, Vector3 dir, float fp) const {
    River r{ riverWeight(l) * std::clamp(2.0f - fp / 12.0f, 0.0f, 1.0f), kNoWater, 1.0f };
    if (r.w <= 0.0f || l.stream >= kStreamZone) { r.w = 0.0f; return r; }
    // Smoothed like rivers, less so (streams follow the land more closely): without the hills'
    // finest octave; faded out where it would need a levee or a deep cut.
    r.level = std::max(broadLand(l, dir, std::max(fp, 32.0f)) - kStreamDrop, 0.3f);
    const float m = r.level - l.h;
    r.w *= (1.0f - Smoothstep(2.5f, 4.0f, m)) * (1.0f - Smoothstep(5.0f, 8.0f, -m));
    return r;
}

PlanetField::River PlanetField::river(const Land& l, Vector3 dir, float fp) const {
    River r{ riverWeight(l), kNoWater, 1.0f };
    if (r.w <= 0.0f) return r;
    // Wider toward the sea: lowland rivers have gathered their tributaries.
    r.wf = 1.0f + 1.5f * (1.0f - Smoothstep(5.0f, 250.0f, l.h));
    if (l.valley >= kRiverZone * r.wf) { r.w = 0.0f; return r; }
    // The water level from the land with only its broad hills (no finer octaves, no detail) and the
    // valley: smooth along the river instead of following every bump of the valley floor.
    r.level = std::max(broadLand(l, dir, std::max(fp, 96.0f)) - kRiverDrop, 0.3f);
    // Where that would stand well above the valley floor (a tall levee) or deep below it (a gorge),
    // the river fades out rather than build walls.
    const float m = r.level - l.h;
    r.w *= (1.0f - Smoothstep(4.0f, 8.0f, m)) * (1.0f - Smoothstep(14.0f, 22.0f, -m));
    return r;
}

PlanetField::Pond PlanetField::pondOfCell(int cx, int cy, int cz) const {
    struct Entry { int x, y, z; const PlanetField* f; Pond p; };
    thread_local Entry cache[512];
    const uint32_t slot = PondHash(cx, cy, cz, 77u) & 511u;
    Entry& e = cache[slot];
    if (e.f == this && e.x == cx && e.y == cy && e.z == cz) return e.p;
    Pond p{ false, { 0, 0, 0 }, 0, 0, 0 };
    uint32_t h = PondHash(cx, cy, cz, (uint32_t)(radius * 7.0f) ^ 0x9e37u);
    if (const PlanetErosion* ero = moon ? nullptr : erosion()) {
        // Ponds from the bake (#360): the small closed depressions its 50 m refinement left, by the 260 m cell their
        // centre lies in (the largest when two share one). The level is held 0.3 m under the ground sampled around the
        // rim (the runtime ground has detail and gullies the bake has not); the bowl only cuts (elevation()).
        std::call_once(erosionState->pondOnce, [this, ero] {
            auto& idx = erosionState->pondByCell;
            const auto& ponds = ero->ponds();
            for (int i = 0; i < (int)ponds.size(); ++i) {
                const Vector3 s = Vector3Scale(ponds[i].dir, radius);
                const uint64_t key = PondId((int)floorf(s.x / kPondCell), (int)floorf(s.y / kPondCell), (int)floorf(s.z / kPondCell));
                auto it = idx.find(key);
                if (it == idx.end() || ponds[it->second].r < ponds[i].r) idx[key] = i;
            }
        });
        const auto it = erosionState->pondByCell.find(PondId(cx, cy, cz));
        if (it != erosionState->pondByCell.end()) {
            const PlanetErosion::Pond& bp = ero->ponds()[it->second];
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(bp.dir, fabsf(bp.dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const Vector3 t2 = Vector3CrossProduct(t1, bp.dir);
            float rim = 1e30f;
            for (int k = 0; k < 8; ++k) {
                const float a = k * (PI / 4.0f), s = bp.r * 1.15f / radius;
                rim = std::min(rim, elevationBase(Vector3Normalize(Vector3Add(bp.dir, Vector3Add(Vector3Scale(t1, cosf(a) * s), Vector3Scale(t2, sinf(a) * s)))), 8.0f));
            }
            const float level = std::min(bp.level, rim - 0.3f);
            const float bottom = elevationBase(bp.dir, 8.0f);
            // The bowl cuts down to the level where the middle stands higher; a middle far below it is a ravine or a
            // gully the bake has not (one pond's middle lay 129 m under its level): no pond there.
            // Ponds stand apart from the network's water: none in a lake's region or by a river (one 16 m from a river
            // inside a lake's region never settled once drained).
            bool apart = true;
            if (const RiverNetwork* net = rivers()) {
                float lv, lo;
                const RiverNetwork::Hit rh = net->nearest(bp.dir, bp.r + 40.0f);
                apart = !net->lakeAt(bp.dir, lv, lo) && !(rh.found && rh.dist < bp.r + rh.halfWidth + 20.0f);
            }
            if (apart && level > 1.0f && level - bottom > -1.0f && level - bottom < 6.0f) {
                p.exists = true; p.dir = bp.dir; p.r = bp.r; p.level = level;
                p.depth = std::clamp(bp.depth, 1.2f, 4.0f);
            }
        }
        e = { cx, cy, cz, this, p };
        return p;
    }
    if (!moon && PondUnit(h) < 0.4f) {
        const Vector3 c = { (cx + 0.2f + 0.6f * PondUnit(h)) * kPondCell, (cy + 0.2f + 0.6f * PondUnit(h)) * kPondCell, (cz + 0.2f + 0.6f * PondUnit(h)) * kPondCell };
        const float len = Vector3Length(c);
        if (len > 1.0f) {
            p.dir = Vector3Scale(c, 1.0f / len);
            p.r = 12.0f + 26.0f * PondUnit(h);
            p.depth = 1.5f + 2.5f * PondUnit(h);
            // The cell's point projects to the sphere: only the cell the projected point lies in
            // owns the pond (each sphere point has one owner cell, so ponds are not duplicated).
            const Vector3 s = Vector3Scale(p.dir, radius);
            const bool owned = (int)floorf(s.x / kPondCell) == cx && (int)floorf(s.y / kPondCell) == cy && (int)floorf(s.z / kPondCell) == cz;
            if (owned) {
                const Land l = land(p.dir, 32.0f);
                if (l.c > 0.06f && l.ridge < 60.0f && l.h > 6.0f && l.h < 700.0f && l.valley > 0.06f && fabsf(p.dir.y) < 0.8f) {
                    // Gentle ground only: the height a radius away on two sides.
                    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(p.dir, fabsf(p.dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                    const Vector3 t2 = Vector3CrossProduct(t1, p.dir);
                    const float a = p.r / radius;
                    const float h1 = land(Vector3Normalize(Vector3Add(p.dir, Vector3Scale(t1, a))), 32.0f).h;
                    const float h2 = land(Vector3Normalize(Vector3Add(p.dir, Vector3Scale(t2, a))), 32.0f).h;
                    if (fabsf(h1 - l.h) < 0.22f * p.r && fabsf(h2 - l.h) < 0.22f * p.r) {
                        p.exists = true;
                        p.level = l.h - 0.3f;
                    }
                }
            }
        }
    }
    e = { cx, cy, cz, this, p };
    return p;
}

float PlanetField::riverBankHeight(Vector3 dir) const {
    if (sandbox) return kSandboxHeight;
    constexpr float fp = 2.0f; // the banks at the detail a player sees
    const Land l = land(dir, fp);
    if (moon) return l.h;
    const float h = l.h - l.relief + l.detail;
    return h + biomeShape(dir, h, fp);
}

float PlanetField::elevationBase(Vector3 dir, float fp) const {
    if (sandbox) return kSandboxHeight;
    const Land l = land(dir, fp);
    if (moon) return l.h;
    const float h = l.h + l.detail;
    return h + biomeShape(dir, h, fp);
}

void PlanetField::gpuNoiseParams(std::vector<float>& out) const {
    // Order = the N_* indices in resources/compute/terrain.glsl.
    const SimplexNoise* fields[] = { &continents, &ridges, &hills, &valleys, &detail, &cavesA, &cavesB, &veins, &climate,
                                     &caverns, &ravines, &ravineMask, &streams, &rocks, &rockMask, &rugged, &cliffs, &warp3,
                                     &tempNoise, &moistNoise, &biomeJitter, &caveMouths };
    out.clear();
    for (const SimplexNoise* f : fields) { float p[4]; f->gpuParams(p); out.insert(out.end(), p, p + 4); }
}

float PlanetField::elevation(Vector3 dir, float fp) const {
    if (sandbox) return kSandboxHeight; // #662: a flat plane
    const Land l = land(dir, fp);
    if (moon) return l.h;
    float h = l.h + l.detail;
    h += biomeShape(dir, h, fp);
    // Rivers of the drainage network: a valley (a floodplain 2.5 channel half-widths across at the
    // water, then sides rising 0.3 m per m) cut into higher ground, low banks where the ground dips
    // below the water, and the channel bed. The valley stays resolved further out than the channel.
    // With the bake, water only cuts (#359): no lake bank, river bank or raised bed (they made pedestals and levees).
    const bool carveOnly = erosion() != nullptr;
    const RiverNetwork* netAll = rivers();
    if (netAll && netAll->bakeNative()) {
        // Bake-native rivers (#547, #608): the bake already wore the valleys and channels; only the channel itself is
        // shaped under the water (the lowest of every piece's bed and bank near, faded out where samples are too far
        // apart to resolve it). No valley carve, no banks, nothing raised.
        // The gullies and cliff faces are detail drawn from the bake's slope, not its drainage: by a river they fade out,
        // so the valley floor it runs on is the bake's own (they stood the water as walls over their troughs).
        if (l.relief != 0.0f) {
            const RiverNetwork::Hit rh = netAll->nearest(dir, 160.0f);
            if (rh.found) h -= l.relief * (1.0f - Smoothstep(rh.halfWidth + 20.0f, rh.halfWidth + 140.0f, rh.dist));
        }
        const float cw = std::clamp(2.0f - fp / 16.0f, 0.0f, 1.0f) * Smoothstep(-2.0f, 0.5f, l.h);
        if (cw > 0.0f) { const float c = netAll->carve(dir, h); h += (c - h) * cw; }
    }
    else if (const RiverNetwork* net = netAll) {
        // Lakes (rivers.h): the level is held under the land around the basin, so the shore is the
        // terrain's own. Between the samples that checked it the ground may still dip lower; a low
        // bank in the outer part of the water region (at most 4 m of fill, 0.3 m above the water)
        // keeps the water in there, and fades out beyond the region.
        float lv = 0.0f, lout = 0.0f;
        const bool nearLake = net->lakeAt(dir, lv, lout);
        const bool inLake = nearLake && lout < RiverNetwork::kLakeMargin;
        if (!carveOnly && nearLake && lout > 0.0f && h < lv + 0.3f) {
            const float M = RiverNetwork::kLakeMargin;
            const float k = Smoothstep(0.3f * M, 0.8f * M, lout) * (1.0f - Smoothstep(1.2f * M, 2.0f * M, lout));
            h += std::min(lv + 0.3f - h, 4.0f) * k;
        }
        const RiverNetwork::Hit rh = net->nearest(dir, net->reach());
        // In a lake, only a river arriving at the lake's level shapes the ground (its valley and bed
        // lead into the lake; no banks, which stood up out of the water as sand ridges); the outflow
        // and any other river leave the lake alone.
        if (rh.found && (!inLake || fabsf(rh.level - lv) < 1.0f)) {
            const float d = rh.dist, half = rh.halfWidth;
            // With the bake the valleys are the terrain's own: only the channel and its near banks are cut (#359).
            const float flat = 2.5f * half, reach = std::min(flat + (carveOnly ? 60.0f : 200.0f), net->reach());
            const float allowed = rh.level + 0.5f + std::max(0.0f, d - flat) * 0.3f;
            const float vw = (1.0f - Smoothstep(reach * 0.7f, reach, d)) * std::clamp(2.0f - fp / 512.0f, 0.0f, 1.0f);
            if (h > allowed) h += (allowed - h) * vw;
            const float cw = std::clamp(2.0f - fp / (4.0f * half), 0.0f, 1.0f) * Smoothstep(-2.0f, 0.5f, l.h); // not out at sea
            const float bank = rh.level + 0.4f;
            if (!carveOnly && !inLake && h < bank && d < half + 12.0f) h += (bank - h) * (1.0f - Smoothstep(half, half + 12.0f, d)) * cw;
            if (d < half) {
                const float t = 1.0f - (d / half) * (d / half);
                const float bed = rh.level + 0.3f - (rh.depth + 0.3f) * t;
                if (!carveOnly || bed < h) h += (bed - h) * cw;
            }
        }
    }
    // Ponds: a bowl to the pond's level inside its radius, an embankment just outside it where the
    // ground falls away. Faded out where samples are too far apart to resolve them.
    const float pw = std::clamp(2.0f - fp / 48.0f, 0.0f, 1.0f);
    if (pw > 0.0f) {
        ForPondCellsNear(Vector3Scale(dir, radius), kPondReach, [&](int x, int y, int z) {
            const Pond p = pondOfCell(x, y, z);
            if (!p.exists) return;
            const float dist = Vector3Distance(dir, p.dir) * radius;
            if (dist >= 1.8f * p.r) return;
            float target;
            if (dist < p.r) {
                const float q = dist / p.r;
                target = p.level + 0.4f - (p.depth + 0.4f) * (1.0f - q * q);
            } else {
                target = std::max(h, p.level + 0.4f);                               // an embankment 0.4 m above the water...
                target = target + (h - target) * Smoothstep(p.r, 1.8f * p.r, dist); // ...blending out to the ground
            }
            if (!carveOnly || target < h) h += (target - h) * pw; // with the bake: the bowl only cuts, no embankment (#360)
        });
    }
    return h;
}

float PlanetField::biomeShape(Vector3 dir, float h, float fp) const {
    if (moon || h < 2.0f) return 0.0f;
    const Climate c = climateAt(dir, h);
    auto smooth = [](float a, float b, float x) { const float t = std::clamp((x - a) / (b - a), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); };
    const float hot = smooth(16.0f, 20.0f, c.temp);
    const float desert = hot * smooth(0.27f, 0.2f, c.moist);
    const float badlands = hot * smooth(0.2f, 0.27f, c.moist) * (1.0f - smooth(0.36f, 0.42f, c.moist));
    const float swamp = smooth(0.7f, 0.82f, c.moist) * smooth(9.0f, 15.0f, c.temp) * (1.0f - smooth(60.0f, 110.0f, h));
    if (desert + badlands + swamp < 1e-3f) return 0.0f;
    const float x = dir.x * radius, y = dir.y * radius, z = dir.z * radius;
    float dh = 0.0f;
    if (desert > 0.0f && fp < 48.0f) {
        // Dunes: ~55 m crests across the prevailing wind, wandering, their height swelling and fading
        // (up to ~9 m); gone where samples are too far apart to show them. Fade near the coast.
        const float warp = 2.2f * biomeJitter.fractal(2, x * 0.25f, y * 0.25f, z * 0.25f);
        const float u = (x * 0.8f + z * 0.6f) / 55.0f + warp;
        // Crest profile: a smoothed triangle wave, squared (close to (0.5 + 0.5 sin)^2.2, which it replaced:
        // only floor and arithmetic, so the GPU port (terrain.glsl) computes exactly the same).
        const float fu = u - floorf(u), tri = 1.0f - fabsf(2.0f * fu - 1.0f), sm = tri * tri * (3.0f - 2.0f * tri);
        const float crest = sm * sm;
        const float swell = 0.55f + 0.45f * moistNoise.fractal(1, x * 6.0f, y * 6.0f, z * 6.0f);
        dh += desert * crest * 9.0f * swell * smooth(4.0f, 20.0f, h) * std::clamp(2.0f - fp / 24.0f, 0.0f, 1.0f);
    }
    // (The badlands' 14 m mesa steps went with the cliff terraces, #543: mesas come from caprock strata in the bake, #587.)
    if (swamp > 0.0f) {
        // Swamps: pulled down toward a low, gently rolling floor 2-5 m above the sea.
        const float floorH = 3.5f + 1.5f * biomeJitter.fractal(1, x, y, z);
        if (h > floorH) dh += (floorH - h) * 0.85f * swamp;
    }
    return dh;
}

float PlanetField::waterElevation(Vector3 dir, float fp) const {
    if (moon || sandbox) return kNoWater;
    float w = kNoWater;
    if (const RiverNetwork* net = rivers()) {
        // A lake's water is its level over its whole region (the terrain decides the shore); rivers
        // in it take its level. Lakes are big: kept out to coarse levels of detail.
        float lake;
        int lid = -1; float lo = 0.0f;
        if (fp < 1024.0f && net->lakeAt(dir, lake, lo, lid) && lo < RiverNetwork::kLakeMargin) w = lake;
        else if (fp < 128.0f) {
            const RiverNetwork::Hit rh = net->nearest(dir, 60.0f);
            if (rh.found && rh.dist < rh.halfWidth * 1.05f) w = rh.level;
        }
    }
    if (fp >= 96.0f) return w; // ponds are not resolved (their bowls fade out by here)
    ForPondCellsNear(Vector3Scale(dir, radius), kPondReach, [&](int x, int y, int z) {
        const Pond p = pondOfCell(x, y, z);
        if (p.exists && Vector3Distance(dir, p.dir) * radius < p.r * 1.05f) w = std::max(w, p.level);
    });
    return w;
}

float PlanetField::waterElevationMesh(Vector3 dir, float fp) const {
    const float w = waterElevation(dir, fp);
    if (w != kNoWater || fp < 0.2f || fp > 8.0f) return w;
    const RiverNetwork* net = rivers();
    if (!net || !net->bakeNative()) return w;
    // The lattice puts a node every 4 samples: a stream narrower than that falls between them (the water at a river's head
    // was never drawn, #609). Out to half a step from the centreline, a sample whose ground stands over the water takes it.
    const float ext = 2.0f * fp + 0.5f;
    const RiverNetwork::Hit rh = net->nearest(dir, ext + 1.0f);
    if (!rh.found || rh.dist >= ext || rh.dist < rh.halfWidth * 1.05f) return w;
    float lv, lo;
    if (net->lakeAt(dir, lv, lo) &&lo < 2.0f * RiverNetwork::kLakeMargin) return w; // lakes keep their own shore
    return elevation(dir, fp) > rh.level - 0.1f ? rh.level : w;
}

bool PlanetField::inlandWaterNear(Vector3 dir, float extent, float fp) const {
    if (moon || sandbox) return false;
    const float sx = dir.x * radius, sy = dir.y * radius, sz = dir.z * radius;
    // Each kind only where samples resolve it (as in elevation() and waterElevation()).
    if (const RiverNetwork* net = rivers()) {
        if (fp < 1024.0f && net->lakeNear(dir, extent)) return true;
        if (fp < 256.0f && net->riverNear(dir, extent)) return true;
    }
    if (fp >= 96.0f) return false;
    bool near = false;
    ForPondCellsNear({ sx, sy, sz }, kPondReach + extent, [&](int x, int y, int z) {
        if (near) return;
        const Pond p = pondOfCell(x, y, z);
        if (p.exists && Vector3Distance(dir, p.dir) * radius < p.r * 1.05f + extent) near = true;
    });
    return near;
}

PlanetField::WaterSource PlanetField::waterSourceAt(Vector3 dir) const {
    WaterSource s;
    if (moon || sandbox) return s;
    if (const RiverNetwork* net = rivers()) {
        float lake;
        int lid = -1; float lo = 0.0f;
        if (net->lakeAt(dir, lake, lo, lid) && lo < RiverNetwork::kLakeMargin) {
            if (const uint64_t cid = drainableLake(lid)) {
                // A small lake (#360): finite water (the simulation drains it like a pond once a region of it is simulated).
                {
                    const RiverNetwork::LakeInfo& li = net->lakeInfo(lid);
                    s.kind = 3; s.level = lake; s.pond = cid; s.pondDir = li.dir; s.pondR = li.r; s.pondDepth = li.depth + 1.0f;
                }
            } else { s.kind = 2; s.level = lake; }
        }
        else {
            const RiverNetwork::Hit rh = net->nearest(dir, 60.0f);
            if (rh.found && rh.dist < rh.halfWidth * 1.05f) { s.kind = 1; s.level = rh.level; }
        }
    }
    ForPondCellsNear(Vector3Scale(dir, radius), kPondReach, [&](int x, int y, int z) {
        const Pond p = pondOfCell(x, y, z);
        if (!p.exists || p.level <= s.level || Vector3Distance(dir, p.dir) * radius >= p.r * 1.05f) return;
        const uint64_t id = PondId(x, y, z);
        s.kind = 3; s.level = p.level; s.pond = id; s.pondDir = p.dir; s.pondR = p.r;
    });
    return s;
}

uint64_t PlanetField::drainableLake(int lakeId) const {
    // Small lakes are finite (user, 2026-11-13: large lakes cannot be drained): at most 2 of the bake's 200 m cells
    // (~8 ha). Only with the bake: there the random ponds are gone and these replace them.
    constexpr int kDrainableCells = 0; // off: a 2-cell lake held ~486,000 m^3 and did not settle (#360); ponds come from the bake instead
    const RiverNetwork* net = erosion() ? rivers() : nullptr;
    if (!net || lakeId < 0 || net->lakeInfo(lakeId).cells > kDrainableCells) return 0;
    return (1ull << 63) | (uint64_t)lakeId;
}


PlanetField::RiverState::~RiverState() = default;

const RiverNetwork* PlanetField::rivers() const {
    static const bool noRivers = EnvFlag("RAYLIBGAME_NO_RIVERS"); // profiling (#365): no river or lake overlay
    if (moon || sandbox || !riverState || noRivers) return nullptr;
    std::call_once(riverState->once, [this] {
        RiverNetwork::Params rp;
        if (erosion()) {
            // Terrain-first rivers (Plan 7 step 5, #359): on the bake's own 200 m grid, from ~2 km^2 of drainage (upland
            // streams), barely meandering (the baked valleys already wind), each level held under the full ground
            // along its curve (rivers only cut now).
            rp.cell = 200.0f; rp.riverArea = 2.0e6f; rp.meander = 0.05f; rp.maxReach = 180.0f; rp.lakeDepth = 2.5f;
            rp.fineHeight = [this](Vector3 d) { return elevationBase(d, 8.0f); };
        }
        const auto t0 = std::chrono::steady_clock::now();
        const PlanetErosion* ero = erosion();
        if (ero && ero->hasFlow()) {
            // Bake-native rivers (#547): the bake's own final drainage on its 50 m grid; streams from 1 km^2.
            rp.riverArea = 1.0e6f; rp.lakeDepth = 1.5f;
            // Runoff (#637): moister climates and a busier water cycle put more water into every km^2 drained (wider, deeper rivers).
            rp.runoff = [this](Vector3 d) { return (0.5f + climateAt(d, 100.0f).moist) * gen.waterCycle; };
            if (!EnvFlag("RAYLIBGAME_RIVER_BANKS_OLD")) rp.bankHeight = [this](Vector3 d) { return riverBankHeight(d); }; // #640 (A/B: =1 off)
            // #903: the biomes' regional reshaping (swamps pulled ~24 m down) moves the water with the ground; at a coarse
            // footprint only the regional part counts (dunes are bumps the channel cuts through). A/B: =1 the old levels.
            if (!EnvFlag("RAYLIBGAME_RIVER_SHIFT_OLD")) rp.regionalShift = [this](Vector3 d, float h) { return biomeShape(d, h, 512.0f); };
            riverState->net = std::make_unique<RiverNetwork>(*ero, radius, rp);
        } else
        // With the bake, the network sees the baked surface itself: land() adds the gullies (+-75 m at this footprint),
        // whose troughs read as closed basins: phantom lakes 130 m under the ground and canyons cut to them (#360).
        riverState->net = std::make_unique<RiverNetwork>(radius, riverSeed, [this, ero](Vector3 d) {
            return ero ? ero->heightAt(d, 256.0f) : land(d, 256.0f).h; // (ice shelves: the network's own rule)
        }, rp);
        LOG_INFO("Rivers", "network built in %.2f s (%d cells of %.0f m)\n", std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count(),
               riverState->net->stats().cells, riverState->net->stats().cellSize);
    });
    return riverState->net.get();
}

bool PlanetField::wetBed(Vector3 local) const {
    const float d = Vector3Length(local);
    if (moon || d < 1.0f) return false;
    const Vector3 dir = Vector3Scale(local, 1.0f / d);
    const float w = waterElevation(dir);
    if (w == kNoWater) return false;
    const float e = elevation(dir);
    return e < w + 0.6f && d > radius + e - 4.0f;
}

float PlanetField::baseDensityLod(Vector3 p, float fp) const {
    const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    if (sandbox) return ISO_LEVEL + DENSITY_GRADIENT * (radius + kSandboxHeight - d); // #662: a flat plane, nothing in it
    // Far above or far below the surface band only the sign matters (no surface to mesh), so skip
    // the height map: most samples of a chunk that does not touch the surface take these paths.
    const float outer = radius + maxElevation, inner = radius + minElevation;
    if (d > outer + 4.0f) return ISO_LEVEL + DENSITY_GRADIENT * (outer - d);
    const float surface = (d < inner - 4.0f || d < 1e-3f) ? inner : radius + elevation({ p.x / d, p.y / d, p.z / d }, fp);
    return densityNear(p, d, surface, fp);
}

namespace {
    constexpr float kBoulderCell = 24.0f;  // m: at most one boulder per cell
    constexpr float kBoulderMax = 5.0f;    // m: largest boulder radius
    struct Boulder { bool exists; Vector3 dir; float r, lift; };
}

// Rock outcrops and boulders (playtest 4: "all terrain is round"). Both are unions with the ground,
// so they mine and collide like any rock, and both fade out where samples are too coarse to show them.
// - Outcrops: tors in rocky patches (a 340 m mask, more on high ground): columns with vertical cliff
//   sides, 3-12 m tall, their tops flaring into a lip (small overhangs), strata on the walls.
// - Boulders: at most one per 24 m cell (owned by the cell its sea-level point lies in, like ponds),
//   1.2-5 m across, sunk 10-40% into the ground, lumpy; none under water or on beaches.
// Returns metres inside the nearest feature (< 0 outside it).
float PlanetField::rockFeature(Vector3 p, float d, float surface, float fp) const {
    if (moon || d < 1e-3f) return -1e9f;
    const float w = std::clamp(2.0f - fp / 4.0f, 0.0f, 1.0f);
    const float h = d - surface;                      // metres above the height-map ground
    if (w <= 0.0f || h > 14.0f || h < -3.0f) return -1e9f; // tops reach 12 m: a lower cap cut crowns flat with a density jump (#242)
    const float elev = surface - radius;
    if (elev < 4.0f) return -1e9f;                    // sea, beaches, shallows
    float best = -1e9f;
    // Outcrops.
    const float mask = Smoothstep(0.5f, 0.75f, rockMask.fractal(2, p.x, p.y, p.z) + 0.2f * Smoothstep(150.0f, 600.0f, elev));
    const Vector3 dir = { p.x / d, p.y / d, p.z / d };
    const Vector3 s = Vector3Scale(dir, radius);
    if (mask > 0.0f) {
        // Tors: the footprint is noise at the sea-level point, the same all the way up a column, so
        // the sides are vertical cliffs; each tor's top height varies, its top flares a little
        // (an overhang lip) and rough strata break the lower walls.
        const float n = rocks.fractal(2, s.x, s.y, s.z);
        const float top = (3.0f + 9.0f * mask) * (0.6f + 0.4f * rocks.noise(s.x * 0.01f, s.y * 0.01f, s.z * 0.01f));
        const float f = std::clamp(h / std::max(top, 0.5f), 0.0f, 1.0f);
        const float thr = 0.3f - 0.06f * f * f + 0.04f * WaveSin(h * 2.1f + n * 9.0f) * (1.0f - Smoothstep(0.35f, 0.7f, f));  // flared top; strata on the walls only (on the crown they made ridges)
        const float side = (n - thr) * 14.0f;                   // ~metres in from the wall
        best = std::max(best, std::min(side, top - h) * w);
    }
    // Boulders.
    const float reach = kBoulderMax * 1.3f;
    const int x0 = (int)floorf((s.x - reach) / kBoulderCell), x1 = (int)floorf((s.x + reach) / kBoulderCell);
    const int y0 = (int)floorf((s.y - reach) / kBoulderCell), y1 = (int)floorf((s.y + reach) / kBoulderCell);
    const int z0 = (int)floorf((s.z - reach) / kBoulderCell), z1 = (int)floorf((s.z + reach) / kBoulderCell);
    for (int cz = z0; cz <= z1; ++cz)
        for (int cy = y0; cy <= y1; ++cy)
            for (int cx = x0; cx <= x1; ++cx) {
                struct Entry { int x, y, z; const PlanetField* f; Boulder b; };
                thread_local Entry cache[1024];
                Entry& e = cache[PondHash(cx, cy, cz, 91u) & 1023u];
                if (!(e.f == this && e.x == cx && e.y == cy && e.z == cz)) {
                    Boulder b{ false, { 0, 0, 0 }, 0, 0 };
                    uint32_t hh = PondHash(cx, cy, cz, (uint32_t)(radius * 3.0f) ^ 0xb01du);
                    if (PondUnit(hh) < 0.16f) {
                        const Vector3 c = { (cx + PondUnit(hh)) * kBoulderCell, (cy + PondUnit(hh)) * kBoulderCell, (cz + PondUnit(hh)) * kBoulderCell };
                        const float len = Vector3Length(c);
                        if (len > 1.0f) {
                            b.dir = Vector3Scale(c, 1.0f / len);
                            const Vector3 q = Vector3Scale(b.dir, radius);
                            if ((int)floorf(q.x / kBoulderCell) == cx && (int)floorf(q.y / kBoulderCell) == cy && (int)floorf(q.z / kBoulderCell) == cz) {
                                const float u = PondUnit(hh);
                                b.r = 1.2f + (kBoulderMax - 1.2f) * u * u;  // mostly small
                                b.lift = b.r * (0.1f + 0.3f * PondUnit(hh));
                                b.exists = true;
                            }
                        }
                    }
                    e = { cx, cy, cz, this, b };
                }
                const Boulder& b = e.b;
                if (!b.exists) continue;
                const float horiz = Vector3Distance(dir, b.dir) * radius;
                if (horiz > b.r * 1.3f) continue;
                const float vert = (h - b.lift) / 0.75f;      // a little squat
                const float lump = 0.25f * b.r * rocks.fractal(2, p.x * 1.7f, p.y * 1.7f, p.z * 1.7f);
                best = std::max(best, (b.r + lump - sqrtf(horiz * horiz + vert * vert)) * w);
            }
    return best;
}

// 3D terrain (plan item 6): a band of full 3D noise from 16 m under the height map to 26 m over it.
// Where a regional mask (~1.6 km patches, more of it the higher the land) makes the land rugged,
// smooth noise stretched upward (features ~2.5x taller than wide: cliff faces), edged with its ridged
// form, and domain-warped
// (crests lean: overhangs, undercuts) raises ridges and carves alcoves, up to 12 m either way.
// The beaches and the low land by the sea stay smooth. Coarse levels of detail keep its low octaves,
// so far silhouettes match near ones. Returns metres added to the surface here (+ solid, - carved).
// How dramatic the 3D terrain is here (Plan 5 C), 0..1, from a cheap one-octave climate (the 3D term runs per
// density sample) blended smoothly so biome borders don't step: badlands and the high mountains most, taiga and
// jungle (karst) a fair amount, forests some, savanna and deserts a little, meadows and swamps not at all (the
// flat building land).
float PlanetField::ruggedness(Vector3 dir, float elev) const {
    const float x = dir.x * radius, y = dir.y * radius, z = dir.z * radius;
    const float t = 30.0f - 44.0f * (fabsf(dir.y) * sqrtf(fabsf(dir.y))) - 6.5f * std::max(elev, 0.0f) / 1000.0f + 7.0f * tempNoise.fractal(1, x, y, z);
    const float m = std::clamp(0.5f + 0.85f * moistNoise.fractal(1, x, y, z) + 0.08f * (1.0f - Smoothstep(10.0f, 150.0f, elev)), 0.0f, 1.0f);
    const float hot = Smoothstep(16.0f, 20.0f, t), cold = 1.0f - Smoothstep(5.0f, 9.0f, t), temperate = (1.0f - hot) * (1.0f - cold);
    float f = 0.35f;  // forests and other temperate land: some
    f += (0.15f - f) * hot * Smoothstep(0.26f, 0.2f, m);                                        // desert
    f += (1.0f - f) * hot * Smoothstep(0.18f, 0.25f, m) * (1.0f - Smoothstep(0.38f, 0.45f, m)); // badlands
    f += (0.25f - f) * hot * Smoothstep(0.38f, 0.45f, m) * (1.0f - Smoothstep(0.55f, 0.62f, m)); // savanna
    f += (0.55f - f) * hot * Smoothstep(0.55f, 0.62f, m);                                        // jungle
    f += (0.0f - f) * temperate * Smoothstep(0.37f, 0.29f, m);                                   // meadow
    f += (0.6f - f) * cold * Smoothstep(0.25f, 0.35f, m);                                        // taiga
    f += (0.0f - f) * Smoothstep(0.68f, 0.78f, m) * (1.0f - Smoothstep(60.0f, 110.0f, elev));    // swamp
    f += (0.9f - f) * Smoothstep(600.0f, 850.0f, elev);                                          // alpine
    return std::clamp(f, 0.0f, 1.0f);
}

// The column part of the 3D band (Plan 5 C): where and how strongly, from the surface point of a column -
// regional rugged patches, wider and stronger in dramatic biomes, none in meadows and swamps. Depends only on
// the column, so blocks compute it once per height-grid node (baseDensityBlock) instead of per sample.
PlanetField::Column3D PlanetField::terrain3DColumn(Vector3 dir, float elev) const {
    Column3D c;
    if (moon || elev < 6.0f || !gTerrain3D) return c;
    const Vector3 s = Vector3Scale(dir, radius + elev);
    const float rn = rugged.fractal(2, s.x, s.y, s.z);
    if (rn < -0.25f) return c;
    const float f = ruggedness(dir, elev);
    c.mask = Smoothstep(-0.25f * f, 0.4f - 0.15f * f, rn) * Smoothstep(6.0f, 60.0f, elev) * Smoothstep(0.0f, 0.2f, f);
    c.f = f;
    return c;
}

// (#971, user 2026-11-19) Wet lowlands keep their cave roofs: in swamps (the climate's swamp rule, smoothed: moist > ~0.74
// under ~90 m) and on rivers' floodplains (within ~120 m of a real river's edge, under ~8 m over its level) a tunnel just
// under the ground stays closed instead of collapsing into a trench (#549). RAYLIBGAME_WET_CAVE_ROOF=0: off.
float PlanetField::caveWetColumn(Vector3 dir, float elev) const {
    static const bool on = !(getenv("RAYLIBGAME_WET_CAVE_ROOF") && getenv("RAYLIBGAME_WET_CAVE_ROOF")[0] == '0');
    if (!on || moon || elev <= 0.0f) return 0.0f;
    const Climate cl = climateAt(dir, elev);
    const float swamp = Smoothstep(0.70f, 0.78f, cl.moist) * (1.0f - Smoothstep(75.0f, 95.0f, elev));
    float flood = 0.0f;
    if (const RiverNetwork* net = rivers()) {
        const RiverNetwork::Hit rh = net->nearest(dir, 160.0f);
        if (rh.found) {
            const float edge = std::max(0.0f, rh.dist - rh.halfWidth), over = elev - rh.level;
            flood = (1.0f - Smoothstep(60.0f, 120.0f, edge)) * (1.0f - Smoothstep(4.0f, 8.0f, over));
        }
    }
    return std::max(swamp, flood);
}

float PlanetField::terrain3D(Vector3 p, float d, float surface, float fp, const Column3D* col) const {
    if (moon || d < 1e-3f) return 0.0f;
    const float h = d - surface;
    if (h < -kBand3DBelow || h > kBand3DAbove) return 0.0f;
    const Vector3 dir = { p.x / d, p.y / d, p.z / d };
    const Column3D c = col ? *col : terrain3DColumn(dir, surface - radius);
    if (c.mask <= 0.0f) return 0.0f;
    const float amp = 16.0f * (0.5f + 1.5f * c.f), top = 20.0f + 26.0f * c.f; // up to ~32 m of rise, shelves up to ~46 m
    Vector3 q = Vector3Subtract(p, Vector3Scale(dir, h * 0.6f)); // vertical stretch
    const float w = 9.0f * warp3.fractalLimited(2, fp, q.x, q.y, q.z);
    q = { q.x + w, q.y - 0.7f * w, q.z + 0.4f * w };
    const float n = cliffs.fractalLimited(3, fp, q.x, q.y, q.z);
    // Smooth warped noise gives the body (overhangs as thick as they are wide); a little of its ridged
    // form sharpens the rock edges. (All-ridged made paper-thin shelves high up.)
    const float shape = n + 0.3f * (1.0f - 2.0f * fabsf(n)) + 0.15f;
    const float fade = (1.0f - Smoothstep(3.0f, top, h)) * Smoothstep(-16.0f, -8.0f, h); // shelves up the band, not floating sheets
    // Raised more than carved: hills and cliffs rise, alcoves are their undercuts (deep carving read as
    // pits in the plains), so dips stop at 3 m.
    return std::max(amp * c.mask * shape * fade, -3.0f);
}

// Floating islands (Plan 5 C, user: big enough for a small/medium factory): one in ~8% of the 4 km cells on
// the sea-level sphere, kept inside its cell. 100-250 m across with a lobed outline; a grassy top with hills
// (~50 m features, up to ~14 m high) terraced into small cliffs in places; a stone underside 40-80 m deep,
// lumpy and tapering; hung 30-60 m above the ground (land below 700 m) or the sea. Trees grow on the tops
// (trees.cpp asks islandTopAt). The ground height under each island is cached per thread.
namespace {
    constexpr float kIslandCell = 4000.0f;
    struct IslandInfo { Vector3 dir; float waist, r, depth; bool ok; };
}
static bool IslandInCell(const PlanetField& pf, int cx, int cy, int cz, IslandInfo& out) {
    uint32_t h = (uint32_t)cx * 73856093u ^ (uint32_t)cy * 19349663u ^ (uint32_t)cz * 83492791u ^ (pf.oreSeed + 0x15a7du) * 2654435761u;
    h ^= h >> 15; h *= 2246822519u; h ^= h >> 13;
    if ((h & 1023) > 82) return false;
    // A direct-mapped cache of plain data: nothing to destroy when a worker thread ends (#485: a thread_local
    // unordered_map here crashed in its destructor at thread exit under MinGW).
    struct Slot { int64_t key; IslandInfo info; };
    constexpr int kSlots = 4096;
    thread_local Slot cache[kSlots];
    thread_local const PlanetField* cacheOwner = nullptr;
    if (cacheOwner != &pf) { for (Slot& s : cache) s.key = -1; cacheOwner = &pf; } // one field's cells (#180: the seed was XORed into the cell bits)
    const int64_t key = ((int64_t)(cx + 32768) << 32) | ((int64_t)(uint16_t)(cy + 32768) << 16) | (int64_t)(uint16_t)(cz + 32768);
    Slot& slot = cache[(uint32_t)((uint64_t)key * 0x9e3779b97f4a7c15ull >> 52) & (kSlots - 1)];
    if (slot.key != key) {
        IslandInfo is{};
        auto u = [&](int shift) { return ((h >> shift) & 255) / 255.0f; };
        const Vector3 c = { (cx + 0.15f + 0.7f * u(10)) * kIslandCell, (cy + 0.15f + 0.7f * u(18)) * kIslandCell, (cz + 0.15f + 0.7f * u(24)) * kIslandCell };
        is.dir = Vector3Normalize(c);
        const Vector3 own = Vector3Scale(is.dir, pf.radius);
        is.ok = (int)floorf(own.x / kIslandCell) == cx && (int)floorf(own.y / kIslandCell) == cy && (int)floorf(own.z / kIslandCell) == cz;
        const float e = is.ok ? pf.elevation(is.dir, 64.0f) : 0.0f;
        is.ok = is.ok && e < 700.0f && fabsf(is.dir.y) < 0.85f;
        is.r = 50.0f + 75.0f * u(5);                          // 100-250 m across
        is.depth = 40.0f + 40.0f * u(2);                       // stone underside 40-80 m deep
        is.waist = pf.radius + std::max(e, 0.0f) + is.depth + 30.0f + 30.0f * u(14); // its lowest point 30-60 m up
        slot.key = key; slot.info = is;
    }
    out = slot.info;
    return out.ok;
}

// The top surface over an island point (metres above the waist), q = 0 centre .. 1 rim, p for the noise.
static float IslandTop(const PlanetField& pf, Vector3 p, float q) {
    const float hills = pf.detail.fractal(2, p.x * 1.4f, p.y * 1.4f, p.z * 1.4f);         // ~50 m hills
    float t = 4.0f + 10.0f * std::max(hills, -0.2f) * (1.0f - q * q);                     // fall away toward the rim
    const float step = 5.0f, terr = t / step, fr = terr - floorf(terr);
    const float terraced = (floorf(terr) + Smoothstep(0.65f, 0.9f, fr)) * step;          // small cliffs
    const float cliffMask = Smoothstep(0.0f, 0.3f, pf.rugged.fractal(1, p.x * 8.0f, p.y * 8.0f, p.z * 8.0f));
    t += (terraced - t) * cliffMask;
    return t + 2.0f * (1.0f - q);                                                         // a dome overall
}

float PlanetField::floatingIsland(Vector3 p, float d) const {
    const Vector3 s = Vector3Scale(p, radius / d);
    IslandInfo is;
    if (!IslandInCell(*this, (int)floorf(s.x / kIslandCell), (int)floorf(s.y / kIslandCell), (int)floorf(s.z / kIslandCell), is)) return -1e9f;
    const float along = Vector3DotProduct(p, is.dir) - is.waist;           // up from the island's waist
    if (along > 30.0f || along < -is.depth - 6.0f) return -1e9f;
    const float lateral = Vector3Length(Vector3Subtract(p, Vector3Scale(is.dir, Vector3DotProduct(p, is.dir))));
    if (lateral > is.r * 1.5f) return -1e9f;
    // Irregular: a lobed outline (~100 m and ~30 m lobes), a lumpy underside of uneven depth.
    const float wobble = 0.35f * detail.fractal(2, p.x * 0.8f, p.y * 0.8f, p.z * 0.8f) + 0.1f * detail.fractal(2, p.x * 2.5f, p.y * 2.5f, p.z * 2.5f);
    const float lump = detail.fractal(2, p.x * 2.0f + 31.0f, p.y * 2.0f, p.z * 2.0f);
    const float rr = is.r * (1.0f + wobble), q = std::clamp(lateral / rr, 0.0f, 1.0f);
    const float top = IslandTop(*this, p, q);
    const float bottom = -is.depth * (0.75f + 0.25f * lump) * ((1.0f - q) * sqrtf(sqrtf(1.0f - q))) - 2.0f;
    return std::min({ top - along, along - bottom, rr - lateral });
}

bool PlanetField::islandTopAt(Vector3 dir, float& topRadius) const {
    if (moon) return false;
    const Vector3 s = Vector3Scale(dir, radius);
    IslandInfo is;
    if (!IslandInCell(*this, (int)floorf(s.x / kIslandCell), (int)floorf(s.y / kIslandCell), (int)floorf(s.z / kIslandCell), is)) return false;
    const Vector3 p = Vector3Scale(dir, is.waist);
    const float lateral = Vector3Length(Vector3Subtract(p, Vector3Scale(is.dir, Vector3DotProduct(p, is.dir))));
    const float wobble = 0.35f * detail.fractal(2, p.x * 0.8f, p.y * 0.8f, p.z * 0.8f) + 0.1f * detail.fractal(2, p.x * 2.5f, p.y * 2.5f, p.z * 2.5f);
    const float q = lateral / (is.r * (1.0f + wobble));
    if (q > 0.92f) return false;
    topRadius = is.waist + IslandTop(*this, p, q);
    return true;
}

// Floating islands within range of a direction (tests, maps): their centres (body-local, at the waist).
int PlanetField::floatingIslandsNear(Vector3 dir, float range, std::vector<Vector3>& out) const {
    out.clear();
    if (moon) return 0;
    const Vector3 s = Vector3Scale(Vector3Normalize(dir), radius);
    const int R = (int)ceilf(range / kIslandCell) + 1;
    const int cx0 = (int)floorf(s.x / kIslandCell), cy0 = (int)floorf(s.y / kIslandCell), cz0 = (int)floorf(s.z / kIslandCell);
    for (int dz = -R; dz <= R; ++dz)
        for (int dy = -R; dy <= R; ++dy)
            for (int dx = -R; dx <= R; ++dx) {
                IslandInfo is;
                if (!IslandInCell(*this, cx0 + dx, cy0 + dy, cz0 + dz, is)) continue;
                if (Vector3Distance(Vector3Scale(is.dir, radius), s) <= range) out.push_back(Vector3Scale(is.dir, is.waist));
            }
    return (int)out.size();
}

// ---- (#615 step 2) The cave graph ----
// Passages are capsule chains generated per 256 m cell of the sea-level projection (as floating islands, in 3D cells of
// the projected point), from the field's seed only: the GPU never repeats the random choices, it gets the segments.
bool PlanetField::CavesGraph() {
    static const bool on = !(getenv("RAYLIBGAME_CAVES") && getenv("RAYLIBGAME_CAVES")[0] == '0'); // on by default since step 3 (#1022)
    return on;
}
namespace {
constexpr float kCaveCell = 256.0f;
constexpr int kCaveCellReach = 3;      // a passage and its branches run at most ~700 m: cells within 3 of a point can reach it
struct CaveSeg { Vector3 a, b, up; float ra, rb, flat; int kind = 0; }; // kind: 0 passage, 1 mouth (rougher, rubble floor), 2 rubble apron (solid)
struct CaveMouth { Vector3 at, out; int kind = 0; float radius = 0.0f, depth = 0.0f; }; // radius, depth: shafts (kind 1, #1067) // (#615 step 3) an entrance: its axis point at the face, outward; kind 0 hillside, 1 karst sink (#1052: at the doline floor, out = up), 2 resurgence; 3 a chamber's centre (#1058: not an entrance)
struct CaveCellData { std::vector<CaveSeg> segs; std::vector<CaveMouth> mouths; };
// (#1055) A cell depends on more than (seed, cell): on the field's land (its gen params: caveEpoch, renewed by setGenParams /
// setSandbox and unique per field) and on whether the erosion bake was visible when it was built. Both are in the key
// (equality only: the hash and so the cell's random numbers stay (seed, cell)), so a read before the bake is attached, or
// an earlier world with the same seed in this process, can no longer leave cells that later reads reuse.
struct CaveCellKey { unsigned seed; int x, y, z; uint32_t epoch = 0; bool baked = false;
    bool operator==(const CaveCellKey& o) const { return seed == o.seed && x == o.x && y == o.y && z == o.z && epoch == o.epoch && baked == o.baked; } };
struct CaveCellHash { size_t operator()(const CaveCellKey& k) const { return ((size_t)k.seed * 2654435761u) ^ ((size_t)k.x * 73856093u) ^ ((size_t)k.y * 19349663u) ^ ((size_t)k.z * 83492791u); } };
struct CaveCache {
    struct Shard { std::mutex m; std::unordered_map<CaveCellKey, std::shared_ptr<const CaveCellData>, CaveCellHash> map; };
    Shard shards[32];
};
CaveCache& Caves() { static CaveCache c; return c; }
// The cave field's own per-point list (single samples: collision, tests); blocks set theirs (sampleBlock).
thread_local const std::vector<float>* tlCaveSegs = nullptr;
struct CaveSegScope {
    const std::vector<float>* prev;
    explicit CaveSegScope(const std::vector<float>* s) : prev(tlCaveSegs) { tlCaveSegs = s; }
    ~CaveSegScope() { tlCaveSegs = prev; }
};
} // namespace

// One cell's passages: background tunnels (step 2), 0-2 per land cell, each a gently turning random walk of 6-16 segments
// 15-30 m long, 30-150 m deep and never shallower than 25 m + its radius (so they never reach the surface), the radius
// wandering between 2 and 8 m, flattened along the local vertical (height = flat x width).
// (#1058) RAYLIBGAME_CAVE_CHAMBERS=0: no chambers (cost A/B; a skipped chamber draws no numbers, so the passages walked after
// it in that cell come out differently: compare volumes and costs, not hashes).
static bool CaveChambersOn() { static const bool on = !(getenv("RAYLIBGAME_CAVE_CHAMBERS") && getenv("RAYLIBGAME_CAVE_CHAMBERS")[0] == '0'); return on; }
static std::shared_ptr<const CaveCellData> CaveCell(const PlanetField& f, unsigned seed, int cx, int cy, int cz) {
    const CaveCellKey key{ seed, cx, cy, cz, f.caveEpoch, f.caveBaked() };
    CaveCache::Shard& sh = Caves().shards[CaveCellHash{}(key) % 32];
    {
        std::lock_guard<std::mutex> lock(sh.m);
        if (auto it = sh.map.find(key); it != sh.map.end()) return it->second;
    }
    auto data = std::make_shared<CaveCellData>();
    std::vector<CaveSeg>* segs = &data->segs;
    std::mt19937 rng((unsigned)(CaveCellHash{}(key) ^ 0x9E3779B9u));
    std::uniform_real_distribution<float> U(0.0f, 1.0f);
    const float R = f.radius;
    // A passage: a gently turning random walk of segments 15-30 m long, kept 25 m + its radius deep (never at the surface),
    // the radius wandering in [rmin, 8] m. Main passages branch (25% of their vertices) into shorter, narrower side passages,
    // so the tunnels form networks (#615, #999: spaghetti and noodle caves). Branches at 35% of main vertices.
    auto minDepth = [](float rad) { return 25.0f + 1.3f * rad; };
    // (#1058, Plan 10 step 5) A chamber at a deep passage vertex c (depth m under the ground; 2% of deep background vertices,
    // 10% of karst ones): an ellipsoid flattened
    // along the local vertical (half-width 15-45 m, mostly the smaller, height x0.3-0.45: wide and low, its roof domed),
    // a level floor 0.4 r below its axis (kind 3), a breakdown floor of 3-8 sunk boulders (solid, r 2-6 m) and 0-2
    // columns floor to roof (solid, r 1.5-4 m). Skipped where its dome would come within 10 m of the surface.
    auto chamber = [&](Vector3 c, Vector3 head, float depth, std::mt19937& g) {
        std::uniform_real_distribution<float> Ug(0.0f, 1.0f);
        const float u0 = Ug(g), u1 = Ug(g), u2 = Ug(g), u3 = Ug(g);
        const float r = 15.0f + 30.0f * u0 * u0, flat = 0.3f + 0.15f * u1;
        if (depth - r * flat * 2.3f < 10.0f) return;
        const Vector3 up = Vector3Normalize(c);
        const Vector3 hd = Vector3Normalize(Vector3Subtract(head, Vector3Scale(up, Vector3DotProduct(head, up))));
        const Vector3 sd = Vector3CrossProduct(up, hd);
        auto add = [&](Vector3 a, Vector3 b, float ra, float rb, float fl, int kind) {
            CaveSeg sg; sg.a = a; sg.b = b; sg.ra = ra; sg.rb = rb; sg.flat = fl; sg.up = up; sg.kind = kind; segs->push_back(sg);
        };
        add(Vector3Subtract(c, Vector3Scale(hd, 0.2f * r)), Vector3Add(c, Vector3Scale(hd, 0.2f * r)), r, r, flat, 3);
        data->mouths.push_back({ c, up, 3 });                            // (for tests: --test-caves RAYLIBGAME_CAVES_CHAMBER=1)
        const float floorH = -0.4f * r * flat;
        auto onFloor = [&](float dist, float ang, float h) {
            return Vector3Add(c, Vector3Add(Vector3Add(Vector3Scale(hd, dist * cosf(ang)), Vector3Scale(sd, dist * sinf(ang))), Vector3Scale(up, h)));
        };
        const int boulders = 3 + (int)(6.0f * u2);
        for (int i = 0; i < boulders; ++i) {
            const float ang = 2.0f * PI * Ug(g), dist = 0.75f * r * sqrtf(Ug(g)), br = 2.0f + 4.0f * Ug(g);
            const Vector3 bc = onFloor(dist, ang, floorH - 0.35f * br);
            add(Vector3Subtract(bc, Vector3Scale(sd, 0.3f)), Vector3Add(bc, Vector3Scale(sd, 0.3f)), br, br, 0.6f, 2);
        }
        const int columns = (int)(3.0f * u3);
        for (int i = 0; i < columns; ++i) {
            const float ang = 2.0f * PI * Ug(g), dist = 0.6f * r * sqrtf(Ug(g)), cr = 1.5f + 2.5f * Ug(g);
            add(onFloor(dist, ang, floorH - 1.0f), onFloor(dist, ang, 1.4f * r * flat), cr, cr, 1.0f, 2);
        }
    };
    std::function<void(Vector3, Vector3, float, float, float, int, int)> walk = [&](Vector3 pos, Vector3 head, float depth, float r, float flat, int n, int level) {
        Vector3 dir = Vector3Normalize(pos);
        const float rmin = level == 0 ? 2.0f : 1.2f, rmax = level == 0 ? 8.0f : 4.0f;
        for (int k = 0; k < n; ++k) {
            const float L = 15.0f + 15.0f * U(rng);
            const float turn = (U(rng) - 0.5f) * 1.2f;                      // up to ~35 degrees a segment
            const Vector3 side = Vector3CrossProduct(dir, head);
            head = Vector3Normalize(Vector3Add(Vector3Scale(head, cosf(turn)), Vector3Scale(side, sinf(turn))));
            const float r2 = std::clamp(r + (U(rng) - 0.5f) * 3.0f, rmin, rmax);
            const float depth2 = std::clamp(depth + (U(rng) - 0.5f) * 0.4f * L, minDepth(r2), 150.0f);
            const Vector3 dir2 = Vector3Normalize(Vector3Add(Vector3Scale(dir, R), Vector3Scale(head, L)));
            const float e2 = f.elevation(dir2);
            if (e2 < 0.0f) break;                                          // not under the sea
            const Vector3 pos2 = Vector3Scale(dir2, R + e2 - depth2);
            // The middle must stay deep too (the ground between may dip): else the passage ends here.
            const Vector3 mid = Vector3Scale(Vector3Add(pos, pos2), 0.5f);
            const float ml = Vector3Length(mid);
            if (R + f.elevation(Vector3Scale(mid, 1.0f / ml)) - ml < minDepth(std::max(r, r2))) break;
            CaveSeg sg;
            sg.a = pos; sg.b = pos2; sg.ra = r; sg.rb = r2; sg.flat = flat;
            sg.up = Vector3Normalize(Vector3Add(pos, pos2));
            segs->push_back(sg);
            pos = pos2; dir = dir2; r = r2; depth = depth2;
            if (level == 0 && depth2 >= 40.0f && U(rng) < 0.02f && CaveChambersOn()) chamber(pos2, head, depth2, rng); // (#1058) a chamber at a deep vertex
            head = Vector3Normalize(Vector3Subtract(head, Vector3Scale(dir, Vector3DotProduct(head, dir)))); // back into the new tangent plane
            if (level == 0 && k + 1 < n && U(rng) < 0.35f) {               // a side passage
                const float bt = (U(rng) < 0.5f ? -1.0f : 1.0f) * (1.0f + 0.7f * U(rng)); // 57-97 degrees off
                const Vector3 s2 = Vector3CrossProduct(dir, head);
                const Vector3 bh = Vector3Normalize(Vector3Add(Vector3Scale(head, cosf(bt)), Vector3Scale(s2, sinf(bt))));
                walk(pos, bh, depth, std::min(r, 1.5f + 2.5f * U(rng)), std::clamp(flat + (U(rng) - 0.5f) * 0.2f, 0.45f, 0.95f), 3 + (int)(U(rng) * 5.0f), 1);
            }
        }
    };
    // (#615) Abundance (the New World "Caves" slider, step 11; RAYLIBGAME_CAVE_ABUNDANCE for now): 0 none, 1 as tuned, 2 ~3x.
    // The first starts draw the same numbers at any abundance, so more caves add passages and keep the others.
    static const float abundance = getenv("RAYLIBGAME_CAVE_ABUNDANCE") ? std::clamp((float)atof(getenv("RAYLIBGAME_CAVE_ABUNDANCE")), 0.0f, 2.0f) : 1.0f;
    const int starts = (int)lroundf(4.0f * (abundance <= 1.0f ? abundance : 1.0f + 2.0f * (abundance - 1.0f)));
    for (int tries = 0; tries < starts; ++tries) {
        const Vector3 s = { (cx + U(rng)) * kCaveCell, (cy + U(rng)) * kCaveCell, (cz + U(rng)) * kCaveCell };
        const float sl = Vector3Length(s);
        if (sl < 1.0f) continue;
        const Vector3 dir = Vector3Scale(s, 1.0f / sl);
        const Vector3 onR = Vector3Scale(dir, R); // the start must project into this cell
        if ((int)floorf(onR.x / kCaveCell) != cx || (int)floorf(onR.y / kCaveCell) != cy || (int)floorf(onR.z / kCaveCell) != cz) continue;
        if (f.elevation(dir) < 5.0f) continue;     // land only
        if (U(rng) < 0.1f) continue;               // not every start
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const float ang = U(rng) * 2.0f * PI;
        const Vector3 head = Vector3Normalize(Vector3Add(Vector3Scale(t1, cosf(ang)), Vector3Scale(Vector3CrossProduct(dir, t1), sinf(ang))));
        const float r = 2.0f + 6.0f * U(rng);
        const float depth = std::max(25.0f + 95.0f * U(rng), minDepth(r));
        walk(Vector3Scale(dir, R + f.elevation(dir) - depth), head, depth, r, 0.5f + 0.4f * U(rng), 6 + (int)(U(rng) * 11.0f), 0);
    }
    // (#615 step 3) Natural mouths: at most one per cell, on a slope of ~20 degrees or more. The entrance runs level into the
    // hill from the slope face (its axis set so the flat floor meets the ground there, its opening ~1.6x the passage, wider
    // than tall), then descends ~1 m every 8 m with the hill above as its roof (3 m+ thick, or there is no mouth here),
    // and once 25 m+ down it goes on as an ordinary branching passage. Caves reach the surface only here (#999: entrances,
    // not accidents). Its own random numbers: the same mouths at any abundance, more of them with more caves.
    {
        std::mt19937 rm((unsigned)(CaveCellHash{}(key) ^ 0x5BD1E995u));
        std::uniform_real_distribution<float> UM(0.0f, 1.0f);
        const float odds = std::min(0.95f, 0.18f * (abundance <= 1.0f ? abundance : 1.0f + 2.0f * (abundance - 1.0f)));
        const bool cellHasMouth = UM(rm) < odds;                           // one roll a cell; the tries only look for a face
        for (int tries = 0; cellHasMouth && tries < 10 && data->mouths.empty(); ++tries) {
            const Vector3 s = { (cx + UM(rm)) * kCaveCell, (cy + UM(rm)) * kCaveCell, (cz + UM(rm)) * kCaveCell };
            const float coin = UM(rm), rr = UM(rm), rf = UM(rm), rd = UM(rm), rn = UM(rm);
            const float sl = Vector3Length(s);
            (void)coin;
            if (sl < 1.0f) continue;
            Vector3 dir = Vector3Scale(s, 1.0f / sl);
            const Vector3 onR = Vector3Scale(dir, R);
            if ((int)floorf(onR.x / kCaveCell) != cx || (int)floorf(onR.y / kCaveCell) != cy || (int)floorf(onR.z / kCaveCell) != cz) continue;
            const float e0 = f.elevation(dir);
            if (e0 < 2.0f) continue;
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const Vector3 t2 = Vector3CrossProduct(dir, t1);
            auto at = [&](Vector3 d0, Vector3 off) { return Vector3Normalize(Vector3Add(Vector3Scale(d0, R), off)); };
            const float g1 = (f.elevation(at(dir, Vector3Scale(t1, 4.0f))) - f.elevation(at(dir, Vector3Scale(t1, -4.0f)))) / 8.0f;
            const float g2 = (f.elevation(at(dir, Vector3Scale(t2, 4.0f))) - f.elevation(at(dir, Vector3Scale(t2, -4.0f)))) / 8.0f;
            const float slope = sqrtf(g1 * g1 + g2 * g2);
            if (slope < 0.36f || slope > 3.0f) continue;                  // ~20-72 degrees
            Vector3 uphill = Vector3Normalize(Vector3Add(Vector3Scale(t1, g1), Vector3Scale(t2, g2)));
            // A real face, not a dimple in a field: the ground rises 5 m+ within 10 m behind the mouth (2 m+ within 5 m),
            // and the opening is sized to fit under it (its roof 1.6 x rMouth x flat over the ground at the face, plus a
            // 1.5 m lip), so the passage enters a hillside and is seen from the front with rock overhead.
            const float rise5 = f.elevation(at(dir, Vector3Scale(uphill, 5.0f))) - e0;
            const float rise10 = f.elevation(at(dir, Vector3Scale(uphill, 10.0f))) - e0;
            const float fall4 = e0 - f.elevation(at(dir, Vector3Scale(uphill, -4.0f)));
            if (rise10 < 5.0f || rise5 < 2.0f || fall4 < 0.5f) continue;
            const float flat = 0.55f + 0.15f * rf;
            const float r = std::min(2.5f + 2.5f * rr, (rise10 - 1.5f) / (1.6f * 1.6f * flat));
            if (r < 2.0f) continue;                                         // the face is too low for a passable mouth
            const float rMouth = 1.6f * r;
            { // (#1022) no tors or boulders (rock over the height map) in or before the opening: they hid it from the front
                const Vector3 side = Vector3Normalize(Vector3CrossProduct(dir, uphill));
                bool blocked = false;
                for (int sx = -1; sx <= 1 && !blocked; ++sx)
                    for (float out = -1.0f; out <= 8.0f && !blocked; out += 3.0f) {
                        const Vector3 dq = at(dir, Vector3Add(Vector3Scale(side, 0.6f * rMouth * sx), Vector3Scale(uphill, -out)));
                        const float sq = R + f.elevation(dq);
                        for (float hq = 0.5f; hq <= 2.5f && !blocked; hq += 1.0f) {
                            const Vector3 pq = Vector3Scale(dq, sq + hq);
                            blocked = f.rockFeature(pq, sq + hq, sq, 0.0f) > -0.5f;
                        }
                    }
                if (blocked) continue;
            }
            const float axisH = e0 + 0.6f * rMouth * flat;                 // the flat floor at the ground at the face
            const float drop = 0.8f + 0.8f * rd;                           // m per 8 m segment
            std::vector<CaveSeg> ent;
            int kindNext = 1;                                              // the mouth's first two segments are rougher, with a rubble floor
            auto seg = [&](Vector3 a, Vector3 b, float ra, float rb) {
                CaveSeg sg; sg.a = a; sg.b = b; sg.ra = ra; sg.rb = rb; sg.flat = flat; sg.up = Vector3Normalize(Vector3Add(a, b)); sg.kind = kindNext; ent.push_back(sg);
            };
            const Vector3 outside = Vector3Scale(at(dir, Vector3Scale(uphill, -rMouth)), R + axisH);
            Vector3 pos = Vector3Scale(dir, R + axisH);
            seg(outside, pos, rMouth, rMouth);
            CaveSeg lip; // (#1022) the overhanging lip: a thin slab of rock following the hillside down to the top of the
            {            // opening and 2 m out past the face (its rounded end included), ragged-edged (solid, kind 2): its inner end
                         // lies 1.2 m under the slope 6 m+ in, its outer end at the roof, so it reads as the hill hanging over
                         // The hill meets roof height some metres in from the face (the face slopes): the lip hangs 2 m out
                         // from there, never in the open air in front of the slope.
                const float roofH = axisH + 0.8f * rMouth * flat;
                float dRoof = 0.0f;
                while (dRoof < 15.0f && f.elevation(at(dir, Vector3Scale(uphill, dRoof))) < roofH + 0.5f) dRoof += 0.5f;
                const float cap = 0.8f * rMouth;
                const float bOff = dRoof - 2.0f + cap;                       // the outer end, so the cap reaches 2 m out
                const float aOff = bOff + std::max(2.0f, cap);
                const Vector3 aDir = at(dir, Vector3Scale(uphill, aOff));
                lip.a = Vector3Scale(aDir, R + f.elevation(aDir) - 1.2f);
                lip.b = Vector3Scale(at(dir, Vector3Scale(uphill, bOff)), R + roofH);
                lip.ra = lip.rb = cap; lip.flat = 0.15f; lip.up = Vector3Normalize(Vector3Add(lip.a, lip.b)); lip.kind = 2;
            }
            float h = axisH, depth = 0.0f;
            bool ok = true;
            int k = 0;
            for (; k < 14; ++k) {
                const Vector3 dir2 = at(dir, Vector3Scale(uphill, 8.0f));
                h -= drop;
                const float e2 = f.elevation(dir2);
                depth = e2 - h;
                const float rb = k == 0 ? r + 0.4f * (rMouth - r) : r;
                const float roof = depth - rb * flat * 2.3f;                // rock over the passage's highest point (swell, bulges, roof pockets)
                if (roof < (k == 0 ? 1.0f : 3.0f)) { ok = false; break; }   // the hill is too thin over it: no mouth here
                const Vector3 pos2 = Vector3Scale(dir2, R + h);
                kindNext = k == 0 ? 1 : 0;
                seg(pos, pos2, k == 0 ? rMouth : (k == 1 ? r + 0.4f * (rMouth - r) : r), rb);
                pos = pos2; dir = dir2;
                uphill = Vector3Normalize(Vector3Subtract(uphill, Vector3Scale(dir, Vector3DotProduct(uphill, dir))));
                if (depth >= minDepth(r)) break;
            }
            if (!ok || k < 1) continue;
            segs->insert(segs->end(), ent.begin(), ent.end());
            { // a rubble apron: a low, lumpy mound of fallen blocks on the slope just outside and below the opening
                const Vector3 od = Vector3Normalize(ent.front().a);
                const Vector3 apDir = at(od, Vector3Scale(Vector3Normalize(Vector3Subtract(ent.front().a, ent.front().b)), 0.2f * rMouth)); // (#1022: low talus, walked over)
                const Vector3 c = Vector3Scale(apDir, R + f.elevation(apDir) - 0.26f * rMouth);
                const Vector3 across = Vector3Normalize(Vector3CrossProduct(apDir, Vector3Subtract(ent.front().a, ent.front().b))); // along the slope's contour
                CaveSeg ap; ap.a = Vector3Subtract(c, Vector3Scale(across, 0.5f * rMouth));
                ap.b = Vector3Add(c, Vector3Scale(across, 0.5f * rMouth));
                ap.ra = ap.rb = 0.65f * rMouth; ap.flat = 0.3f; ap.up = apDir; ap.kind = 2;
                segs->push_back(ap);
            }
            segs->push_back(lip);
            data->mouths.push_back({ Vector3Scale(dir, 0.0f), Vector3Negate(uphill) });
            data->mouths.back().at = ent.front().b;                        // the axis point at the slope face
            data->mouths.back().out = Vector3Normalize(Vector3Subtract(ent.front().a, ent.front().b));
            if (depth >= minDepth(r)) walk(pos, uphill, depth, r, flat, 4 + (int)(rn * 8.0f), 0); // on into the network
        }
    }
    // (#1052, Plan 10 step 4) Karst river caves: in karst country (#592, New World "Karst"; weight > 0.3) a doline whose
    // centre falls in this cell may be a sink: a round shaft from its floor down 10-30 m with a rubble cone at its foot, and
    // from there a river passage following the drainage (the 240 m-scale downhill, not the bowl's own), its axis falling
    // 1-5% of the way: keyhole sections first (a tube with a narrow slot above it, made of two stacked round tubes), then
    // round tubes; level side passages at 40% of its vertices make mazes. Where the roof thins out after 100 m+ it has met
    // a valley side: it opens there (a resurgence). At most ~640 m, inside the cells' reach (kCaveCellReach). Each doline
    // draws its own numbers (seeded by its id): the same sinks at any abundance, more of them with more.
    {
        const Vector3 clo = { cx * kCaveCell, cy * kCaveCell, cz * kCaveCell };
        std::vector<KarstDoline> dolines;
        Strata_KarstDolines(clo, Vector3AddValue(clo, kCaveCell), R, seed, dolines);
        const float sinkOdds = std::min(0.95f, 0.6f * (abundance <= 1.0f ? abundance : 1.0f + 2.0f * (abundance - 1.0f)));
        auto at = [&](Vector3 d0, Vector3 off) { return Vector3Normalize(Vector3Add(Vector3Scale(d0, R), off)); };
        auto frame = [](Vector3 d0, Vector3& a1, Vector3& a2) {
            a1 = Vector3Normalize(Vector3CrossProduct(d0, fabsf(d0.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            a2 = Vector3CrossProduct(d0, a1);
        };
        auto downhill = [&](Vector3 d0, Vector3 fallback) {               // the drainage direction at the 240 m scale
            Vector3 a1, a2; frame(d0, a1, a2);
            const float g1 = f.elevation(at(d0, Vector3Scale(a1, 120.0f))) - f.elevation(at(d0, Vector3Scale(a1, -120.0f)));
            const float g2 = f.elevation(at(d0, Vector3Scale(a2, 120.0f))) - f.elevation(at(d0, Vector3Scale(a2, -120.0f)));
            const Vector3 g = Vector3Add(Vector3Scale(a1, -g1), Vector3Scale(a2, -g2));
            const float gl = Vector3Length(g);
            return gl > 1.0f ? Vector3Scale(g, 1.0f / gl) : fallback;
        };
        auto push = [&](Vector3 a, Vector3 b, float ra, float rb, float flat, int kind) {
            CaveSeg sg; sg.a = a; sg.b = b; sg.ra = ra; sg.rb = rb; sg.flat = flat; sg.up = Vector3Normalize(Vector3Add(a, b)); sg.kind = kind;
            segs->push_back(sg);
        };
        for (const KarstDoline& dl : dolines) {
            if (f.karstWeight(Vector3Scale(dl.dir, R)) < 0.3f) continue;
            const float e0 = f.elevation(dl.dir);
            if (e0 < 3.0f) continue;                                         // a flooded hollow
            std::mt19937 rk(dl.id ^ 0x6B4A57u);
            std::uniform_real_distribution<float> UK(0.0f, 1.0f);
            float rnd[8];
            for (float& v : rnd) v = UK(rk);
            if (rnd[0] >= sinkOdds) continue;
            const Vector3 up = dl.dir;
            Vector3 a1, a2; frame(up, a1, a2);
            // The shaft: round, its walls vertical, from 3 m over the bowl's floor down to the passage; the cone on its floor.
            const float sr = 3.0f + 3.0f * rnd[1], sdepth = 10.0f + 20.0f * rnd[2];
            const Vector3 bot = Vector3Scale(up, R + e0 - sdepth);
            push(Vector3Scale(up, R + e0 + 3.0f), bot, sr, 0.85f * sr, 1.0f, 0);
            const Vector3 coneC = Vector3Scale(up, R + e0 - sdepth - 0.45f * sr);
            push(Vector3Add(coneC, Vector3Scale(a1, 0.3f)), Vector3Subtract(coneC, Vector3Scale(a1, 0.3f)), 0.7f * sr, 0.7f * sr, 0.35f, 2);
            data->mouths.push_back({ Vector3Scale(up, R + e0), up, 1, sr, sdepth });
            // The river passage.
            const float slope = 0.01f + 0.04f * rnd[3];
            const int n = 20 + (int)(12.0f * rnd[4]);
            const float keyLen = (0.25f + 0.35f * rnd[5]) * n * 20.0f;
            const float ang = rnd[6] * 2.0f * PI;
            const Vector3 head0 = downhill(up, Vector3Add(Vector3Scale(a1, cosf(ang)), Vector3Scale(a2, sinf(ang))));
            std::function<void(Vector3, Vector3, Vector3, float, float, float, int, int, float)> kwalk =
                [&](Vector3 pos, Vector3 dir, Vector3 head, float h, float r, float fall, int steps, int level, float keyRun) {
                float run = 0.0f;
                for (int k = 0; k < steps; ++k) {
                    const float L = 15.0f + 10.0f * UK(rk);
                    const float turn = (UK(rk) - 0.5f) * 0.7f;
                    head = Vector3Normalize(Vector3Add(Vector3Scale(head, cosf(turn)), Vector3Scale(Vector3CrossProduct(dir, head), sinf(turn))));
                    if (level == 0) head = Vector3Normalize(Vector3Add(Vector3Scale(head, 0.7f), Vector3Scale(downhill(dir, head), 0.3f)));
                    const float r2 = std::clamp(r + (UK(rk) - 0.5f) * 1.5f, level ? 1.2f : 2.0f, level ? 3.0f : 5.0f);
                    const float h2 = h - fall * L;
                    const Vector3 dir2 = at(dir, Vector3Scale(head, L));
                    const float e2 = f.elevation(dir2);
                    if (e2 < 0.0f) break;                                    // not under the sea
                    const bool key = run < keyRun;
                    const float topR = key ? 2.25f * r2 : 0.95f * r2;        // the section's top over the axis
                    const Vector3 pos2 = Vector3Scale(dir2, R + h2);
                    if (e2 - h2 - 1.5f * topR < 3.0f) {                      // the roof is thinning out (pockets x1.5)
                        if (level == 0 && run > 100.0f) {                    // a valley side: a resurgence
                            const Vector3 fwd = Vector3Normalize(Vector3Subtract(pos2, pos));
                            push(pos, Vector3Add(pos2, Vector3Scale(fwd, 1.6f * r2)), r, 1.6f * r2, 0.8f, 1);
                            data->mouths.push_back({ pos2, fwd, 2 });
                        }
                        break;
                    }
                    push(pos, pos2, r, r2, key ? 0.95f : 0.9f, 0);
                    if (key) {                                               // the keyhole's slot: two stacked round tubes
                        const Vector3 u1 = Vector3Normalize(pos), u2 = Vector3Normalize(pos2);
                        push(Vector3Add(pos, Vector3Scale(u1, 1.0f * r)), Vector3Add(pos2, Vector3Scale(u2, 1.0f * r2)), 0.45f * r, 0.45f * r2, 1.0f, 0);
                        push(Vector3Add(pos, Vector3Scale(u1, 1.8f * r)), Vector3Add(pos2, Vector3Scale(u2, 1.8f * r2)), 0.45f * r, 0.45f * r2, 1.0f, 0);
                    }
                    pos = pos2; dir = dir2; r = r2; h = h2; run += L;
                    if (level == 0 && e2 - h2 >= 40.0f && UK(rk) < 0.10f && CaveChambersOn()) chamber(pos2, head, e2 - h2, rk); // (#1058) karst: more chambers
                    head = Vector3Normalize(Vector3Subtract(head, Vector3Scale(dir, Vector3DotProduct(head, dir))));
                    if (level == 0 && k + 1 < steps && UK(rk) < 0.4f) {     // a maze branch, level and narrower
                        const float bt = (UK(rk) < 0.5f ? -1.0f : 1.0f) * (1.0f + 0.6f * UK(rk));
                        const Vector3 bh = Vector3Normalize(Vector3Add(Vector3Scale(head, cosf(bt)), Vector3Scale(Vector3CrossProduct(dir, head), sinf(bt))));
                        kwalk(pos, dir, bh, h, std::min(r, 1.5f + 1.5f * UK(rk)), 0.0f, 3 + (int)(4.0f * UK(rk)), 1, 0.0f);
                    }
                }
            };
            kwalk(bot, up, head0, e0 - sdepth, 2.0f + 3.0f * rnd[7], slope, n, 0, keyLen);
        }
    }
    // (#1063, Plan 10 step 6) Lava tubes: a cell on a volcano's flank (weight 0.3-0.85: below the crater, above the toe;
    // New World "Volcanoes") starts one at 1% x abundance (its own numbers): a sinuous tube walking down the 120 m-scale
    // slope 5-25 m under the ground, radius 2.5-6.5 m and height x0.6-0.8 (5-13 m wide, 3-10 m high), with the level floor
    // and arched roof of kind 3, thin ledges along one wall (half its segments, solid), skylights (round shafts to the
    // surface at 12% of its segments under 15 m; entrances of kind 1), and an exit where its roof thins out (kind 2).
    // At most ~640 m (the cells' reach).
    if (f.gen.volcanoes) {
        std::mt19937 rt((unsigned)(CaveCellHash{}(key) ^ 0x1A7A70BEu));
        std::uniform_real_distribution<float> UT(0.0f, 1.0f);
        const float tubeOdds = std::min(0.95f, 0.01f * (abundance <= 1.0f ? abundance : 1.0f + 2.0f * (abundance - 1.0f)));
        const float coin = UT(rt), sx = UT(rt), sy = UT(rt), sz = UT(rt);
        const Vector3 s = { (cx + sx) * kCaveCell, (cy + sy) * kCaveCell, (cz + sz) * kCaveCell };
        const float sl = Vector3Length(s);
        float vw = 0.0f;
        if (coin < tubeOdds && sl > 1.0f) Strata_Volcano(Vector3Scale(s, R / sl), R, seed, &vw);
        Vector3 dir = sl > 1.0f ? Vector3Scale(s, 1.0f / sl) : Vector3{ 0, 1, 0 };
        const Vector3 onR = Vector3Scale(dir, R);
        const bool inCell = (int)floorf(onR.x / kCaveCell) == cx && (int)floorf(onR.y / kCaveCell) == cy && (int)floorf(onR.z / kCaveCell) == cz;
        if (coin < tubeOdds && inCell && vw > 0.3f && vw < 0.85f && f.elevation(dir) > 5.0f) {
            auto at = [&](Vector3 d0, Vector3 off) { return Vector3Normalize(Vector3Add(Vector3Scale(d0, R), off)); };
            auto push = [&](Vector3 a, Vector3 b, float ra, float rb, float flat, int kind) {
                CaveSeg sg; sg.a = a; sg.b = b; sg.ra = ra; sg.rb = rb; sg.flat = flat; sg.up = Vector3Normalize(Vector3Add(a, b)); sg.kind = kind;
                segs->push_back(sg);
            };
            auto downhill = [&](Vector3 d0, Vector3 fallback) {
                const Vector3 a1 = Vector3Normalize(Vector3CrossProduct(d0, fabsf(d0.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const Vector3 a2 = Vector3CrossProduct(d0, a1);
                const float g1 = f.elevation(at(d0, Vector3Scale(a1, 60.0f))) - f.elevation(at(d0, Vector3Scale(a1, -60.0f)));
                const float g2 = f.elevation(at(d0, Vector3Scale(a2, 60.0f))) - f.elevation(at(d0, Vector3Scale(a2, -60.0f)));
                const Vector3 g = Vector3Add(Vector3Scale(a1, -g1), Vector3Scale(a2, -g2));
                const float gl = Vector3Length(g);
                return gl > 1.0f ? Vector3Scale(g, 1.0f / gl) : fallback;
            };
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const float ang = 2.0f * PI * UT(rt);
            Vector3 head = downhill(dir, Vector3Normalize(Vector3Add(Vector3Scale(t1, cosf(ang)), Vector3Scale(Vector3CrossProduct(dir, t1), sinf(ang)))));
            float r = 2.5f + 4.0f * UT(rt), depth = 6.0f + 10.0f * UT(rt);
            const float flat = 0.6f + 0.2f * UT(rt);
            const int n = 20 + (int)(12.0f * UT(rt));
            Vector3 pos = Vector3Scale(dir, R + f.elevation(dir) - depth);
            for (int k = 0; k < n; ++k) {
                const float L = 15.0f + 10.0f * UT(rt);
                const float turn = (UT(rt) - 0.5f) * 0.5f;                     // sinuous
                head = Vector3Normalize(Vector3Add(Vector3Scale(head, cosf(turn)), Vector3Scale(Vector3CrossProduct(dir, head), sinf(turn))));
                head = Vector3Normalize(Vector3Add(Vector3Scale(head, 0.6f), Vector3Scale(downhill(dir, head), 0.4f)));
                const float r2 = std::clamp(r + (UT(rt) - 0.5f) * 1.5f, 2.5f, 6.5f);
                const float d2 = std::clamp(depth + (UT(rt) - 0.5f) * 3.0f, 5.0f, 25.0f);
                const Vector3 dir2 = at(dir, Vector3Scale(head, L));
                const float e2 = f.elevation(dir2);
                if (e2 < 0.0f) break;
                const Vector3 pos2 = Vector3Scale(dir2, R + e2 - d2);
                const float ledgeU = UT(rt), skyU = UT(rt), sideU = UT(rt);
                // The roof must hold over the arch (pockets x1.5): on the lower flank where it thins out, the tube opens.
                if (d2 - 1.5f * r2 * flat < 3.0f) {
                    if (k > 3) {
                        const Vector3 fwd = Vector3Normalize(Vector3Subtract(pos2, pos));
                        push(pos, Vector3Add(pos2, Vector3Scale(fwd, 1.6f * r2)), r, 1.6f * r2, 0.8f, 1);
                        data->mouths.push_back({ pos2, fwd, 2 });
                    }
                    break;
                }
                push(pos, pos2, r, r2, flat, 3);
                if (ledgeU < 0.5f) {                                           // a thin shelf along one wall, mid-height
                    const Vector3 up = Vector3Normalize(Vector3Add(pos, pos2));
                    Vector3 side = Vector3CrossProduct(up, Vector3Subtract(pos2, pos));
                    side = Vector3Scale(Vector3Normalize(side), sideU < 0.5f ? -1.0f : 1.0f);
                    push(Vector3Add(pos, Vector3Scale(side, 0.85f * r)), Vector3Add(pos2, Vector3Scale(side, 0.85f * r2)), 0.22f * r, 0.22f * r2, 0.35f, 2);
                }
                if (skyU < 0.12f && d2 < 15.0f && k > 0) {                    // a skylight: a round shaft to the surface
                    const float sr = std::max(1.5f, 0.35f * r2);
                    push(Vector3Scale(dir2, R + e2 + 3.0f), pos2, sr, sr, 1.0f, 0);
                    data->mouths.push_back({ Vector3Scale(dir2, R + e2), dir2, 1, sr, d2 });
                }
                pos = pos2; dir = dir2; r = r2; depth = d2;
                head = Vector3Normalize(Vector3Subtract(head, Vector3Scale(dir, Vector3DotProduct(head, dir))));
            }
        }
    }
    // (#1181) Sinkholes: hillside mouths need a real face, so on flat or gently rolling ground (a swamp or meadow start)
    // there were none for kilometres. A collapse pit: a round shaft 2-4 m across at a gentle spot of the cell's ground,
    // flaring into a funnel at the top, 8-16 m down to a rubble cone, and a passage from its foot wandering level-ish
    // under the ground (it ends where its roof would thin under 4 m). ~1 in 16 surface cells at 1x (~1 per km^2), its own
    // random numbers (the same sinks at any abundance, more with more); none in karst (the dolines are its sinks), the sea
    // or water. Mouth kind 1 (a shaft), so #1067's light pool falls under it. Built from the same segment kinds as the
    // others: caveDistance and terrain.glsl are unchanged. RAYLIBGAME_CAVE_SINKHOLES=0: off.
    static const bool sinkholesOn = !(getenv("RAYLIBGAME_CAVE_SINKHOLES") && getenv("RAYLIBGAME_CAVE_SINKHOLES")[0] == '0');
    if (sinkholesOn) {
        std::mt19937 rs((unsigned)(CaveCellHash{}(key) ^ 0x51A4C0DEu));
        std::uniform_real_distribution<float> US(0.0f, 1.0f);
        const float sinkOdds = std::min(0.95f, 0.065f * (abundance <= 1.0f ? abundance : 1.0f + 2.0f * (abundance - 1.0f)));
        const float roll = US(rs);
        const Vector3 cc = { (cx + 0.5f) * kCaveCell, (cy + 0.5f) * kCaveCell, (cz + 0.5f) * kCaveCell };
        const float cl = Vector3Length(cc);
        if (roll < sinkOdds && cl > 1.0f) {
            const Vector3 d0 = Vector3Scale(cc, 1.0f / cl);
            const Vector3 a1 = Vector3Normalize(Vector3CrossProduct(d0, fabsf(d0.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const Vector3 a2 = Vector3CrossProduct(d0, a1);
            auto at = [&](Vector3 base, Vector3 off) { return Vector3Normalize(Vector3Add(Vector3Scale(base, R), off)); };
            auto push = [&](Vector3 a, Vector3 b, float ra, float rb, float flat, int kind) {
                CaveSeg sg; sg.a = a; sg.b = b; sg.ra = ra; sg.rb = rb; sg.flat = flat; sg.up = Vector3Normalize(Vector3Add(a, b)); sg.kind = kind;
                segs->push_back(sg);
            };
            for (int tries = 0; tries < 8; ++tries) {
                const float u = (US(rs) - 0.5f) * 1.4f * kCaveCell, v = (US(rs) - 0.5f) * 1.4f * kCaveCell;
                const Vector3 dir = at(d0, Vector3Add(Vector3Scale(a1, u), Vector3Scale(a2, v)));
                const Vector3 onR = Vector3Scale(dir, R);
                if ((int)floorf(onR.x / kCaveCell) != cx || (int)floorf(onR.y / kCaveCell) != cy || (int)floorf(onR.z / kCaveCell) != cz) continue;
                const float e0 = f.elevation(dir);
                if (e0 < 3.0f) continue;                                               // the sea, a shore, a flooded hollow
                if (f.karstWeight(onR) >= 0.3f) continue;                              // karst has its dolines
                if (f.inlandWaterNear(dir, 20.0f)) continue;                           // not into a river or a lake
                const float g1 = (f.elevation(at(dir, Vector3Scale(a1, 6.0f))) - f.elevation(at(dir, Vector3Scale(a1, -6.0f)))) / 12.0f;
                const float g2 = (f.elevation(at(dir, Vector3Scale(a2, 6.0f))) - f.elevation(at(dir, Vector3Scale(a2, -6.0f)))) / 12.0f;
                if (g1 * g1 + g2 * g2 > 0.18f * 0.18f) continue;                       // gentle ground (~10 degrees): hills have mouths
                const float sr = 1.0f + 1.0f * US(rs), sdepth = 8.0f + 8.0f * US(rs);
                const Vector3 top = Vector3Scale(dir, R + e0 + 2.0f), rim = Vector3Scale(dir, R + e0 - 1.5f);
                const Vector3 bot = Vector3Scale(dir, R + e0 - sdepth);
                push(top, rim, 2.2f * sr, 1.15f * sr, 1.0f, 0);                       // the funnel: a collapsed lip
                push(rim, bot, 1.15f * sr, sr, 1.0f, 0);                               // the shaft
                const Vector3 coneC = Vector3Scale(dir, R + e0 - sdepth - 0.45f * sr);
                push(Vector3Add(coneC, Vector3Scale(a1, 0.3f)), Vector3Subtract(coneC, Vector3Scale(a1, 0.3f)), 0.8f * sr, 0.8f * sr, 0.35f, 2);
                data->mouths.push_back({ Vector3Scale(dir, R + e0), dir, 1, sr, sdepth });
                // The passage from the foot: 6-12 segments, turning gently, falling ~1 m in 15, radius 1.8-3.2 m.
                const float ang = US(rs) * 2.0f * PI;
                Vector3 head = Vector3Normalize(Vector3Add(Vector3Scale(a1, cosf(ang)), Vector3Scale(a2, sinf(ang))));
                Vector3 pos = Vector3Scale(dir, R + e0 - sdepth + 0.6f * sr), pdir = dir;
                float r = 1.6f + 0.6f * US(rs), h = e0 - sdepth + 0.6f * sr;
                const int n = 6 + (int)(7.0f * US(rs));
                for (int k = 0; k < n; ++k) {
                    const float L = 12.0f + 8.0f * US(rs);
                    const float turn = (US(rs) - 0.5f) * 0.8f;
                    head = Vector3Normalize(Vector3Add(Vector3Scale(head, cosf(turn)), Vector3Scale(Vector3CrossProduct(pdir, head), sinf(turn))));
                    const float r2 = std::clamp(r + (US(rs) - 0.5f) * 1.2f, 1.8f, 3.2f);
                    const float h2 = h - L / 15.0f;
                    const Vector3 dir2 = at(pdir, Vector3Scale(head, L));
                    const float e2 = f.elevation(dir2);
                    if (e2 - h2 - 1.5f * r2 < 4.0f) break;                             // the roof would thin out: it ends here
                    const Vector3 pos2 = Vector3Scale(dir2, R + h2);
                    push(pos, pos2, r, r2, 0.9f, 0);
                    pos = pos2; pdir = dir2; r = r2; h = h2;
                    head = Vector3Normalize(Vector3Subtract(head, Vector3Scale(pdir, Vector3DotProduct(head, pdir))));
                }
                break;
            }
        }
    }
    std::lock_guard<std::mutex> lock(sh.m);
    if (sh.map.size() > 4096) sh.map.clear();
    return sh.map.emplace(key, data).first->second;
}

// Segments of the cells within reach of a box (body-local), whose capsule (radius x 1.3 for the wall roughness) comes
// within `reach` of it. (seed: fieldSeed; the cells by the box corners' sea-level projections, +-kCaveCellReach.)
static void GatherCaveSegs(const PlanetField& f, unsigned seed, Vector3 lo, Vector3 hi, float reach, std::vector<float>& out) {
    out.clear();
    int clo[3] = { INT32_MAX, INT32_MAX, INT32_MAX }, chi[3] = { INT32_MIN, INT32_MIN, INT32_MIN };
    for (int k = 0; k < 8; ++k) {
        const Vector3 q = { (k & 1) ? hi.x : lo.x, (k & 2) ? hi.y : lo.y, (k & 4) ? hi.z : lo.z };
        const float d = Vector3Length(q);
        if (d < 1.0f) return;
        const Vector3 s = Vector3Scale(q, f.radius / d);
        const int c[3] = { (int)floorf(s.x / kCaveCell), (int)floorf(s.y / kCaveCell), (int)floorf(s.z / kCaveCell) };
        for (int a = 0; a < 3; ++a) { clo[a] = std::min(clo[a], c[a] - kCaveCellReach); chi[a] = std::max(chi[a], c[a] + kCaveCellReach); }
    }
    for (int z = clo[2]; z <= chi[2]; ++z)
        for (int y = clo[1]; y <= chi[1]; ++y)
            for (int x = clo[0]; x <= chi[0]; ++x) {
                const auto cell = CaveCell(f, seed, x, y, z);
                for (const CaveSeg& s : cell->segs) {
                    // The box beyond which the shape can't reach: the skeleton at most 1.3 x 1.4 r (swell, the large-scale term), its
                    // centre wandering 0.4 of that sideways (2.55 r in all; the warp moves only the noise; solids 1.3 r); plus the reach.
                    const float rmax = std::max(s.ra, s.rb);
                    const float rr = 3.0f * rmax + 2.0f + reach + std::max(0.0f, 0.8f * rmax - 6.4f); // (#1058: big chambers keep the single-sample bound)
                    const Vector3 slo = Vector3SubtractValue(Vector3Min(s.a, s.b), rr), shi = Vector3AddValue(Vector3Max(s.a, s.b), rr);
                    if (shi.x < lo.x || slo.x > hi.x || shi.y < lo.y || slo.y > hi.y || shi.z < lo.z || slo.z > hi.z) continue;
                    Vector3 side = Vector3CrossProduct(s.up, Vector3Subtract(s.b, s.a));     // across the passage, level
                    const float sl = Vector3Length(side);
                    side = sl > 1e-4f ? Vector3Scale(side, 1.0f / sl) : Vector3{ 0, 0, 0 };
                    const float e[24] = { s.a.x, s.a.y, s.a.z, s.ra, s.b.x, s.b.y, s.b.z, s.rb, s.up.x, s.up.y, s.up.z, s.flat,
                                          slo.x, slo.y, slo.z, (float)s.kind, shi.x, shi.y, shi.z, 0.0f, side.x, side.y, side.z, 0.0f };
                    out.insert(out.end(), e, e + 24);
                }
            }
}

bool PlanetField::nearestCaveMouth(Vector3 dir, float range, Vector3& at, Vector3& outward, int kind) const {
    if (!CavesGraph() || moon) return false;
    const Vector3 s = Vector3Scale(Vector3Normalize(dir), radius);
    const int R = (int)ceilf(range / kCaveCell) + 1;
    const int c0[3] = { (int)floorf(s.x / kCaveCell), (int)floorf(s.y / kCaveCell), (int)floorf(s.z / kCaveCell) };
    float best = range;
    bool found = false;
    for (int z = -R; z <= R; ++z)
        for (int y = -R; y <= R; ++y)
            for (int x = -R; x <= R; ++x) {
                const Vector3 cc = { (c0[0] + x + 0.5f) * kCaveCell, (c0[1] + y + 0.5f) * kCaveCell, (c0[2] + z + 0.5f) * kCaveCell };
                if (fabsf(Vector3Length(cc) - radius) > 2.0f * kCaveCell) continue; // only cells on the surface shell
                const auto cell = CaveCell(*this, fieldSeed, c0[0] + x, c0[1] + y, c0[2] + z);
                for (const CaveMouth& m : cell->mouths) {
                    if (kind >= 0 ? m.kind != kind : m.kind >= 3) continue; // "any" = entrances (not chamber centres, #1058)
                    const float dd = Vector3Distance(Vector3Scale(Vector3Normalize(m.at), radius), s);
                    if (dd < best) { best = dd; at = m.at; outward = m.out; found = true; }
                }
            }
    return found;
}

void PlanetField::shaftsNear(Vector3 lo, Vector3 hi, std::vector<CaveShaft>& out) const {
    // (#1067) The shaft entrances (kind 1) whose shaft (top to its depth, plus the light pool under its foot) can reach
    // the box: from the cells around the box's corners projected onto the sea-level sphere.
    out.clear();
    if (!CavesGraph() || moon) return;
    int clo[3] = { INT32_MAX, INT32_MAX, INT32_MAX }, chi[3] = { INT32_MIN, INT32_MIN, INT32_MIN };
    for (int k = 0; k < 8; ++k) {
        const Vector3 q = { (k & 1) ? hi.x : lo.x, (k & 2) ? hi.y : lo.y, (k & 4) ? hi.z : lo.z };
        const float d = Vector3Length(q);
        if (d < 1.0f) return;
        const Vector3 s = Vector3Scale(q, radius / d);
        const int c[3] = { (int)floorf(s.x / kCaveCell), (int)floorf(s.y / kCaveCell), (int)floorf(s.z / kCaveCell) };
        for (int a = 0; a < 3; ++a) { clo[a] = std::min(clo[a], c[a] - 1); chi[a] = std::max(chi[a], c[a] + 1); }
    }
    for (int z = clo[2]; z <= chi[2]; ++z)
        for (int y = clo[1]; y <= chi[1]; ++y)
            for (int x = clo[0]; x <= chi[0]; ++x)
                for (const CaveMouth& m : CaveCell(*this, fieldSeed, x, y, z)->mouths) {
                    if (m.kind != 1 || m.radius <= 0.0f) continue;
                    const float reach = 3.0f * m.radius + 12.0f;               // the pool's reach around the shaft
                    const Vector3 foot = Vector3Subtract(m.at, Vector3Scale(m.out, m.depth + 12.0f));
                    const Vector3 blo = Vector3SubtractValue(Vector3Min(m.at, foot), reach), bhi = Vector3AddValue(Vector3Max(m.at, foot), reach);
                    if (bhi.x < lo.x || blo.x > hi.x || bhi.y < lo.y || blo.y > hi.y || bhi.z < lo.z || blo.z > hi.z) continue;
                    out.push_back({ m.at, m.out, m.radius, m.depth });
                }
}

int PlanetField::caveMouthsNear(Vector3 dir, float range) const {
    if (!CavesGraph() || moon) return 0;
    const Vector3 s = Vector3Scale(Vector3Normalize(dir), radius);
    const int R = (int)ceilf(range / kCaveCell) + 1;
    const int c0[3] = { (int)floorf(s.x / kCaveCell), (int)floorf(s.y / kCaveCell), (int)floorf(s.z / kCaveCell) };
    int n = 0;
    for (int z = -R; z <= R; ++z)
        for (int y = -R; y <= R; ++y)
            for (int x = -R; x <= R; ++x) {
                const Vector3 cc = { (c0[0] + x + 0.5f) * kCaveCell, (c0[1] + y + 0.5f) * kCaveCell, (c0[2] + z + 0.5f) * kCaveCell };
                if (fabsf(Vector3Length(cc) - radius) > 2.0f * kCaveCell) continue;
                for (const CaveMouth& m : CaveCell(*this, fieldSeed, c0[0] + x, c0[1] + y, c0[2] + z)->mouths)
                    n += m.kind < 3 && Vector3Distance(Vector3Scale(Vector3Normalize(m.at), radius), s) <= range; // entrances only (#1058)
            }
    return n;
}

bool PlanetField::blockCaveSegs(Vector3 o, float stride, int S, float fp, std::vector<float>& out, int maxSegs) const {
    out.clear();
    if (!CavesGraph() || moon) return true;
    if (2.0f - fp / 8.0f <= 0.0f) return true;            // caveWeight 0 at this footprint: no caves to carve
    const float ext = stride * (float)(S - 1);
    const float reach = 2.0f * std::max(fp, 1.0f) + 2.0f; // densityNear's meshing reach
    GatherCaveSegs(*this, fieldSeed, o, Vector3AddValue(o, ext), reach, out);
    return (int)(out.size() / 24) <= maxSegs;
}

// The graph tunnels' distance at q (metres to the nearest wall), and `solid` (metres inside a rubble apron, < 0 outside).
// Organic shapes (#615 step 3, #1022; user 2026-11-19: "define the actual path and diameter, then apply heavy 3D noise
// inwards"): the skeleton is the graph segment's bore (radius x a +-30% swell and a 10 m outline term, a gently flattened
// ellipse whose centre wanders sideways, a floor 20-70% down that tilts across the passage, roof pockets). On it, heavy
// domain-warped 3D noise (2 octaves from ~4 m plus a gentle ~2 m one) only ever grows rock inward from the wall: none
// inside half the bore, so a clear core always runs through, and up to ~0.6 of the bore near the wall, so walls, ledges and
// overhangs come out of the noise. Mouths: a larger outline term, taller on one side, a level floor with low rubble. Solids
// (kind 2): rubble aprons and the lip over a mouth. A segment whose box doesn't hold q is skipped and the noises are looked
// up only once one does. The same operations in the same order as terrain.glsl's caveGraphDistance (bit-exact).
struct CaveNoises { const SimplexNoise &a, &b, &det, &cav, &mou; };
// RAYLIBGAME_CAVE_WALLS=0 (#1022 A/B): the cave walls as first merged (v22); default 1, the larger outcrops. The GPU gets it
// per block (c.w).
static const int gCaveWalls = getenv("RAYLIBGAME_CAVE_WALLS") && getenv("RAYLIBGAME_CAVE_WALLS")[0] == '0' ? 0 : 1;
static float CaveGraphDistance(const CaveNoises& N, Vector3 q, const float* segs, int n, float& solid, float shrink = 0.0f) {
    solid = -1e9f;
    bool any = false;
    float swell = 0.0f, lev = 0.0f, tilt = 0.0f, pocket = 0.0f, rub = 0.0f, shift = 0.0f, big = 0.0f;
    float m = 1e9f;
    for (int i = 0; i < n; ++i) {
        const float* s = segs + 24 * i;
        if (q.x < s[12] + shrink || q.y < s[13] + shrink || q.z < s[14] + shrink || q.x > s[16] - shrink || q.y > s[17] - shrink || q.z > s[18] - shrink) continue;
        if (!any) {
            any = true;
            // The noises carry their own base frequencies (a, b 1/150 m, det 1/70, cav 1/220, mou 1/700), so each factor
            // below is (the wanted features per metre) x (1 / that frequency): 0.09/m on a 1/150 noise is 13.5, and so on.
            // (#1022: before this, every factor was 70-700x too low, and all these terms were constant over a passage.)
            const float wx = q.x + 3.0f * N.a.fractal(1, q.x * 13.5f + 17.0f, q.y * 13.5f, q.z * 13.5f);   // warp: 3 m at ~11 m
            const float wy = q.y + 3.0f * N.a.fractal(1, q.x * 13.5f, q.y * 13.5f + 29.0f, q.z * 13.5f);
            const float wz = q.z + 3.0f * N.a.fractal(1, q.x * 13.5f, q.y * 13.5f, q.z * 13.5f + 41.0f);
            swell = 1.0f + 0.3f * N.b.fractal(1, wx * 15.0f, wy * 15.0f, wz * 15.0f);                     // the bore, +-30% over ~10 m
            lev = 0.45f + 0.25f * N.cav.fractal(1, wx * 13.2f, wy * 13.2f, wz * 13.2f);    // floor level, ~17 m
            tilt = 0.6f * N.mou.fractal(1, wx * 35.0f + 7.0f, wy * 35.0f, wz * 35.0f);       // floor tilt, ~20 m
            pocket = std::max(0.0f, N.det.fractal(1, wx * 14.0f + 53.0f, wy * 14.0f, wz * 14.0f)); // roof pockets, ~5 m
            rub = N.det.fractal(1, wx * 35.0f + 71.0f, wy * 35.0f, wz * 35.0f);             // rubble, ~2 m
            shift = 0.4f * N.b.fractal(1, wx * 6.0f + 13.0f, wy * 6.0f, wz * 6.0f);         // the bore wanders to one side, ~25 m
            big = N.cav.fractal(1, wx * 22.0f + 91.0f, wy * 22.0f, wz * 22.0f);             // (#1022) the large-scale outline, ~10 m
        }
        const float abx = s[4] - s[0], aby = s[5] - s[1], abz = s[6] - s[2];
        const float apx = q.x - s[0], apy = q.y - s[1], apz = q.z - s[2];
        const float num = apx * abx + apy * aby + apz * abz;
        const float den = abx * abx + aby * aby + abz * abz;
        const float t = std::clamp(num / den, 0.0f, 1.0f);
        const float hx = apx - abx * t, hy = apy - aby * t, hz = apz - abz * t;
        const float hv = hx * s[8] + hy * s[9] + hz * s[10];
        const float hhx = hx - s[8] * hv, hhy = hy - s[9] * hv, hhz = hz - s[10] * hv;
        const float hs = hhx * s[20] + hhy * s[21] + hhz * s[22];
        const float r = s[3] + (s[7] - s[3]) * t;
        const float kind = s[15];
        const bool mouth = kind > 0.5f && kind < 1.5f, hall = kind > 2.5f; // kinds: 0 passage, 1 mouth, 2 solid, 3 chamber (#1058)
        if (kind > 1.5f && kind < 2.5f) { // solid inside: a rubble apron, the lip over a mouth, a chamber's boulder or column
            const float va = hv / s[11];                                  // squashed: aprons 0.35, the lip over a mouth 0.2
            const float la2 = hhx * hhx + hhy * hhy + hhz * hhz + va * va;
            const float la = sqrtf(la2);
            const float rr = r * (1.0f + 0.3f * rub);
            solid = std::max(solid, rr - la);
            continue;
        }
        const float bk = mouth ? 0.4f : 0.2f;                       // the large-scale term: mouths' noise-cut outline
        const float rough = 1.0f + bk * big;
        const float rs = swell * rough;
        const float rw = r * std::max(0.7f, rs);                           // the skeleton: the clear bore (never pinched shut)
        const float pk = hv > 0.0f ? pocket : 0.0f;
        const float vfac = s[11] * (1.0f + 0.5f * pk);                    // roof pockets
        const float hp = hs - shift * rw;
        const float hvT = mouth ? hv - 0.5f * hp : hv;              // mouths: taller on one side
        const float v2 = hvT / vfac;
        const float l2 = std::max(0.0f, hhx * hhx + hhy * hhy + hhz * hhz - hs * hs + hp * hp) + v2 * v2;
        const float len = sqrtf(l2);
        const float rubble = mouth ? 0.35f * std::max(0.0f, rub) : 0.0f; // mouths: rubble mounds on the floor (<= ~0.3 m)
        const float levK = mouth ? 0.45f : lev;                      // mouths: a predictable, level floor (no tilt), so
        const float tiltK = mouth || hall ? 0.0f : tilt;                     //   the opening meets the ground and stays open
        const float flo = hall ? 0.4f * r : levK * rw;                     // chambers: a level floor 0.4 r down, so boulders sit on it
        const float floorCut = (0.0f - hv) - flo * s[11] + 0.5f * tiltK * hp + rubble;
        // (#1022, user) Heavy, jagged 3D noise grows the wall inward: none inside half the bore (a clear core; the floor
        // stays the floor), up to ~0.6 of it near the wall. Voids behind the faces come from the field falling again
        // further out along a ray from the axis, i.e. from the noise's radial change: the noise is looked up at a point
        // squeezed radially toward the skeleton's wall (radial distance x 0.1 around it), so it changes at full rate
        // across the wall (jagged faces) but 10x slower outward (solid rock behind them, no shells or floating pieces).
        // Only in the shell where it can matter (inside half the bore it is 0; 40 m out the distance is never used).
        float intrude = 0.0f;
        if (len > 0.5f * rw && len < rw + 40.0f) {
            const float hl2 = hx * hx + hy * hy + hz * hz;
            const float hl = sqrtf(hl2);
            if (hl > 1e-3f) {
                const float sc = 0.1f + 0.9f * rw / hl;
                const float px = s[0] + abx * t + hx * sc, py = s[1] + aby * t + hy * sc, pz = s[2] + abz * t + hz * sc;
                // (#1022, user 2026-11-20) Large irregular outcrops: rock forms ~9 m (2 octaves) and creases ~3.2 m (|n|: hard
                // edges), weights 2.9 / 0.8. RAYLIBGAME_CAVE_WALLS=0: the walls as first merged (v22: ~3 m forms in 3 octaves,
                // ~1 m creases in 2, 2.6 / 1.4).
                const bool v22 = gCaveWalls == 0;
                const int o3 = v22 ? 3 : 2, o1 = v22 ? 2 : 1;
                const float f3 = v22 ? 24.5f : 7.8f, f1 = v22 ? 135.0f : 47.0f;
                const float g3 = v22 ? 2.6f : 2.9f, g1 = v22 ? 1.4f : 0.8f;
                const float r3 = N.det.fractal(o3, px * f3 + 101.0f, py * f3, pz * f3);
                const float r1 = N.b.fractal(o1, px * f1 + 7.0f, py * f1, pz * f1);
                const float wn = std::max(0.0f, 0.35f + g3 * r3 + g1 * (0.3f - fabsf(r1)));
                const float fa = (len - 0.5f * rw) / (0.4f * rw);
                const float fc = std::clamp(fa, 0.0f, 1.0f);
                const float fade = fc * fc * (3.0f - 2.0f * fc);
                const float ra = hall ? std::min(rw, 8.0f) : rw;            // chambers: outcrops at most passage-sized (else tens of metres)
                intrude = (mouth ? 0.65f : 0.55f) * ra * fade * wn;
            }
        }
        m = std::min(m, std::max(len - rw + intrude, floorCut));
    }
    return m;
}

// (#615 step 1) The caves: worm tunnels (with #549's mouth regions and their collapsing thin roofs, kept closed in wet
// lowlands by #971) and caverns. The same operations in the same order as terrain.glsl's caveDistance (bit-exact).
float PlanetField::caveDistance(Vector3 p, float d, float depth, float surface, float fp, const Column3D* col, CaveHit* hit, float* solidOut) const {
    if (solidOut) *solidOut = -1e9f;
    // Caves, from a few metres down to 20 m above the mantle (which stays intact). Each cave shape is
    // its own distance field in metres, combined with min(): rock density grows with depth, so
    // subtracting a fixed amount would only carve shallow caves.
    // - Worm tunnels (~14 m radius) where both cave fields cross zero: at 1/150 m the fields change
    //   ~0.01 per metre, so (t - kTunnelLevel) * 100 ~ metres from the tunnel wall.
    // - Caverns: rare chambers tens of metres across where a low-frequency field peaks, 30 m+ down.
    // Coarse LOD samples can't resolve them, so both fade out there.
    const float caveWeight = std::clamp(2.0f - fp / 8.0f, 0.0f, 1.0f);
    if (CavesGraph()) { // (#615 step 2) graph tunnels (never at the surface: no mouths yet) and the caverns
        float m = 1e9f;
        if (caveWeight > 0.0f && depth > -12.0f && d > mantleRadius + 20.0f) { // (#615 step 3: up to the surface at mouths)
            if (depth > -12.0f) {
                const float* segData = tlCaveSegs ? tlCaveSegs->data() : nullptr;
                int n = tlCaveSegs ? (int)(tlCaveSegs->size() / 24) : 0;
                if (!tlCaveSegs) { // a single sample: this 256 m cell's list, kept per thread (reach 36 m: any footprint that carves),
                    // binned (#1037) into 8^3 bins of 32 m, each a contiguous copy of the segments whose box overlaps it, so a
                    // sample tests only its bin's few instead of the whole cell's. The same result: CaveGraphDistance uses only
                    // segments whose box holds q, and every one of those overlaps q's bin (min / max don't depend on order).
                    constexpr int kBins = 8;
                    constexpr float kBin = kCaveCell / kBins;
                    struct CaveBins { std::vector<float> own, binned; std::vector<int> binStart; }; // (#1133, #485: MinGW frees thread_local storage before running its destructor at thread exit: never destroyed)
                    thread_local CaveBins* bins = nullptr; if (!bins) bins = new CaveBins(); auto& own = bins->own; auto& binned = bins->binned; auto& binStart = bins->binStart; // binStart: kBins^3 + 1 offsets (in segments) into binned
                    thread_local int ok = 0, kx = 0, ky = 0, kz = 0; thread_local unsigned kseed = 0, kepoch = 0; thread_local bool kbaked = false; thread_local const void* kf = nullptr;
                    const int cx = (int)floorf(p.x / kCaveCell), cy = (int)floorf(p.y / kCaveCell), cz = (int)floorf(p.z / kCaveCell);
                    const Vector3 lo = { cx * kCaveCell, cy * kCaveCell, cz * kCaveCell };
                    if (!ok || cx != kx || cy != ky || cz != kz || kseed != fieldSeed || kf != this || kepoch != caveEpoch || kbaked != (erosion() != nullptr)) { // (#1055: + epoch, bake)
                        GatherCaveSegs(*this, fieldSeed, lo, Vector3AddValue(lo, kCaveCell), 36.0f, own);
                        const int ns = (int)(own.size() / 24);
                        auto binRange = [&](const float* s, int a, int& b0, int& b1) {
                            const float l = a == 0 ? lo.x : a == 1 ? lo.y : lo.z;
                            b0 = std::clamp((int)floorf((s[12 + a] - l) / kBin), 0, kBins - 1);
                            b1 = std::clamp((int)floorf((s[16 + a] - l) / kBin), 0, kBins - 1);
                        };
                        std::vector<int> count(kBins * kBins * kBins + 1, 0);
                        for (int pass = 0; pass < 2; ++pass) {                  // count, then fill
                            if (pass == 1) {
                                binStart.assign(kBins * kBins * kBins + 1, 0);
                                for (int b = 0; b < kBins * kBins * kBins; ++b) binStart[b + 1] = binStart[b] + count[b];
                                binned.assign((size_t)binStart.back() * 24, 0.0f);
                                std::fill(count.begin(), count.end(), 0);
                            }
                            for (int i = 0; i < ns; ++i) {
                                const float* s = own.data() + 24 * i;
                                int x0, x1, y0, y1, z0, z1; binRange(s, 0, x0, x1); binRange(s, 1, y0, y1); binRange(s, 2, z0, z1);
                                for (int z = z0; z <= z1; ++z) for (int y = y0; y <= y1; ++y) for (int x = x0; x <= x1; ++x) {
                                    const int b = (z * kBins + y) * kBins + x;
                                    if (pass == 1) std::copy(s, s + 24, binned.begin() + (size_t)(binStart[b] + count[b]) * 24);
                                    ++count[b];
                                }
                            }
                        }
                        ok = 1; kx = cx; ky = cy; kz = cz; kseed = fieldSeed; kf = this; kepoch = caveEpoch; kbaked = erosion() != nullptr;
                    }
                    const int bx = std::clamp((int)floorf((p.x - lo.x) / kBin), 0, kBins - 1);
                    const int by = std::clamp((int)floorf((p.y - lo.y) / kBin), 0, kBins - 1);
                    const int bz = std::clamp((int)floorf((p.z - lo.z) / kBin), 0, kBins - 1);
                    const int b = (bz * kBins + by) * kBins + bx;
                    segData = binned.data() + (size_t)binStart[b] * 24;
                    n = binStart[b + 1] - binStart[b];
                }
                if (n > 0) {
                    float sol = -1e9f;
                    // (#1037) A single sample's list holds boxes grown by 36 m (any footprint that carves); a segment matters to
                    // this sample's density only within the meshing reach (2 x fp + 2, at least 4 m: the solids' -4 cut). Out of a
                    // segment's shape box by d, its distance is >= 0.59 d - 2.57 and its solid < -4 (bores flat <= 1, so at most
                    // x1.5 tall with roof pockets; radius <= 8 m; mouths' lopsided height x0.89; #1052 corrected the #1037 bound,
                    // which assumed x1.05), so testing the box shrunk to 1.7 x that reach + 4.5 m changes no density but skips
                    // most noise lookups in rock near caves. Blocks: 0.
                    const float keep = std::max(2.0f * std::max(fp, 1.0f) + 2.0f, 4.0f) * 1.7f + 4.5f;
                    const float shrink = tlCaveSegs ? 0.0f : std::max(0.0f, 36.0f - keep);
                    m = CaveGraphDistance(CaveNoises{ cavesA, cavesB, detail, caverns, caveMouths }, p, segData, n, sol, shrink) + (1.0f - caveWeight) * 20.0f;
                    if (solidOut && sol > -1e8f) *solidOut = sol - (1.0f - caveWeight) * 20.0f;
                    if (hit && m < 1e8f) hit->kind = 1;
                }
            }
            // (#1058) No noise caverns with graph caves: chambers in the graph replace them (RAYLIBGAME_CAVES=0 keeps them).
        }
        return m;
    }
    // Cave mouths (Plan 5 C): in some regions (~15%) tunnels may reach up through the surface, so caves are found
    // by walking; elsewhere they stay 5 m+ down. Only looked up for shallow samples.
    float caveTop = 5.0f;
    float mouth = 0.0f; // #549: in mouth regions tunnels open at full width where they meet the surface (no grille)
    if (depth < 14.0f && depth > -14.0f && caveWeight > 0.0f && !moon) {
        mouth = Smoothstep(0.3f, 0.45f, caveMouths.fractal(1, p.x, p.y, p.z));
        // (#971) wet lowlands keep their roofs: no mouth there (the GPU reads the same interpolated column, terrain.glsl)
        if (mouth > 0.0f) mouth = mouth * (1.0f - (col ? col->wet : caveWetColumn(Vector3Scale(p, 1.0f / d), surface - radius)));
        caveTop = 5.0f - 7.5f * mouth;
    }
    if (hit) hit->mouth = mouth;
    // #549 (playtest 8): a tunnel running just under the ground left a roof thin enough for the surface noise to punch
    // rows of small holes through it (the grille), and tor slabs bridged the mouth. In mouth regions a roof thinner than
    // kRoof collapses: every sample over the tunnel whose roof (its depth + its distance to the tunnel) is that thin is
    // carved. The roof depth is about the same all down a column, so whole columns open, with clean edges; up to the
    // tops of tors (14 m) so slabs over a collapsed roof go too.
    constexpr float kRoof = 4.0f;
    const float gate = mouth > 0.0f ? -14.0f : caveTop;
    float m = 1e9f;                                            // metres to the nearest cave wall (< 0 inside)
    if (caveWeight > 0.0f && depth > gate && d > mantleRadius + 20.0f) {
        const float near = 2.0f * std::max(fp, 1.0f) + 2.0f; // farther than this from a cave = plain rock for meshing
        const float reach = mouth > 0.0f ? std::max(near, kRoof + 1.0f - depth) : near; // a thin roof's samples too
        // A tunnel needs both fields near zero: if the first alone rules one out, skip the second.
        const float a = fabsf(cavesA.fractal(2, p.x, p.y, p.z));
        if ((a - kTunnelLevel) * 100.0f <= reach) {
            float top = std::clamp((depth - caveTop) / 8.0f, 0.0f, 1.0f);
            top = top + (1.0f - top) * mouth;
            const float fade = top * std::min(1.0f, (d - mantleRadius - 20.0f) / 10.0f);
            const float t = std::max(a, fabsf(cavesB.fractal(2, p.x, p.y, p.z)));
            m = (t - kTunnelLevel) * 100.0f + (1.0f - fade * caveWeight) * 20.0f;
            if (hit) hit->kind = 1;
            if (mouth > 0.0f && m > 0.0f) {
                const float roof = depth + m;
                const float c = mouth * std::clamp((kRoof + 1.0f - roof) * 0.5f, 0.0f, 1.0f);
                m = m + (-2.0f - m) * c;
                if (hit && c > 0.0f) hit->collapsed = true;
            }
        }
        if (depth > 30.0f && d > mantleRadius + 40.0f) {
            const float fade = std::min(1.0f, (depth - 30.0f) / 20.0f) * std::min(1.0f, (d - mantleRadius - 40.0f) / 20.0f);
            const float c = caverns.fractal(2, p.x, p.y * 1.4f, p.z); // a little flattened: wider than tall
            const float mc = (kCavernLevel - c) * 120.0f + (1.0f - fade * caveWeight) * 30.0f;
            if (hit && mc < m) hit->kind = 2;
            m = std::min(m, mc);
        }
    }
    return m;
}

// Everything after the height map, shared by single samples and blocks so both agree exactly.
float PlanetField::densityNear(Vector3 p, float d, float surface, float fp, const Column3D* col) const {
    float v = ISO_LEVEL + DENSITY_GRADIENT * (surface - d);
    const float depth = surface - d;
    if (depth < kBand3DBelow && depth > -kBand3DAbove) v += DENSITY_GRADIENT * terrain3D(p, d, surface, fp, col); // hills, cliffs, overhangs
    if (fabsf(depth) < 6.0f) v += DENSITY_GRADIENT * 1.5f * detail.fractalLimited(2, fp, p.x, p.y * 1.3f, p.z); // crags, overhangs
    if (depth < 3.0f && depth > -14.0f) {
        const float m = rockFeature(p, d, surface, fp);
        if (m > -4.0f) v = std::max(v, ISO_LEVEL + DENSITY_GRADIENT * m);
    }
    // Caves (#615 step 1: one shared caveDistance, bit-exact with terrain.glsl's): carved where under the meshing reach.
    {
        float solid = -1e9f;
        const float m = caveDistance(p, d, depth, surface, fp, col, nullptr, &solid);
        const float near = 2.0f * std::max(fp, 1.0f) + 2.0f; // farther than this from a cave = plain rock for meshing
        if (m < near) v = std::min(v, ISO_LEVEL + DENSITY_GRADIENT * m);
        if (solid > -4.0f) v = std::max(v, ISO_LEVEL + DENSITY_GRADIENT * solid); // (#615 step 3) rubble aprons, after the carve
    }
    // Floating islands (Plan 5 C, rare): only for samples well above the ground.
    if (depth < -15.0f && !moon) { const float m = floatingIsland(p, d); if (m > -4.0f) v = std::max(v, ISO_LEVEL + DENSITY_GRADIENT * m); }
    return v;
}

void PlanetField::baseDensityBlock(Vector3 o, float stride, int S, float fp, float* out) const {
    BlockGrid g;
    if (!prepareBlock(o, stride, S, fp, g)) { DensityField::baseDensityBlock(o, stride, S, fp, out); return; }
    sampleBlock(o, stride, S, fp, g, out);
}

bool PlanetField::blockIslands(Vector3 o, float stride, int S, std::vector<float>& out, int maxCells) const {
    out.clear();
    if (moon) return true;
    // The cells of the samples' sea-level points: bounded by the cells of the block corners' projections
    // (a block is small next to a 4 km cell, so its projection is close to the corners' box).
    const float ext = stride * (float)(S - 1);
    int lo[3] = { INT32_MAX, INT32_MAX, INT32_MAX }, hi[3] = { INT32_MIN, INT32_MIN, INT32_MIN };
    for (int k = 0; k < 8; ++k) {
        const Vector3 q = { o.x + ((k & 1) ? ext : 0.0f), o.y + ((k & 2) ? ext : 0.0f), o.z + ((k & 4) ? ext : 0.0f) };
        const float d = Vector3Length(q);
        if (d < 1.0f) return false;
        const Vector3 s = Vector3Scale(q, radius / d);
        const int c[3] = { (int)floorf(s.x / kIslandCell), (int)floorf(s.y / kIslandCell), (int)floorf(s.z / kIslandCell) };
        for (int a = 0; a < 3; ++a) { lo[a] = std::min(lo[a], c[a] - 1); hi[a] = std::max(hi[a], c[a] + 1); }
    }
    if ((hi[0] - lo[0] + 1) * (hi[1] - lo[1] + 1) * (hi[2] - lo[2] + 1) > maxCells) return false;
    for (int z = lo[2]; z <= hi[2]; ++z)
        for (int y = lo[1]; y <= hi[1]; ++y)
            for (int x = lo[0]; x <= hi[0]; ++x) {
                IslandInfo is;
                if (!IslandInCell(*this, x, y, z, is)) continue;
                const float e[12] = { (float)x, (float)y, (float)z, 1.0f, is.dir.x, is.dir.y, is.dir.z, is.waist, is.r, is.depth, 0.0f, 0.0f };
                out.insert(out.end(), e, e + 12);
            }
    return true;
}

bool PlanetField::prepareBlock(Vector3 o, float stride, int S, float fp, BlockGrid& g) const {
    if (sandbox) return false; // #662: per sample (the plane is cheap; the GPU path is off too)
    const float ext = stride * (float)(S - 1);
    const Vector3 c = { o.x + 0.5f * ext, o.y + 0.5f * ext, o.z + 0.5f * ext };
    const float dc = sqrtf(c.x * c.x + c.y * c.y + c.z * c.z), halfDiag = 0.8661f * ext;
    const float outer = radius + maxElevation, inner = radius + minElevation;
    // Blocks wholly outside the surface band take the cheap per-sample paths anyway; blocks too large
    // for a tangent-plane grid (near the centre) fall back to the per-sample loop too.
    if (dc - halfDiag > outer + 4.0f || dc + halfDiag < inner - 4.0f || halfDiag > 0.5f * dc) return false;
    // Height map on a gnomonic grid tangent at the block centre's direction: every sample along the
    // same line of sight from the planet centre shares one height, so a 35^3 block needs ~35^2
    // height evaluations instead of 35^3. Bilinear interpolation at one stride spacing is accurate
    // to centimetres (the finest height detail kept at this footprint is >= 2 strides long).
    // The plane is the cube face of the block's direction and the nodes sit on that face's global lattice (#251):
    // neighbouring blocks of a level share their nodes, so the interpolated surface is one function across them.
    // A grid tangent at each block's own centre interpolated a shared border from different nodes, and where the
    // height map kinks (pond bowls and banks, river beds) the borders disagreed: cracks along chunk seams. Only
    // the 12 cube edges still join two lattices.
    const int ax = fabsf(c.x) >= fabsf(c.y) && fabsf(c.x) >= fabsf(c.z) ? 0 : fabsf(c.y) >= fabsf(c.z) ? 1 : 2;
    const float sgn = (ax == 0 ? c.x : ax == 1 ? c.y : c.z) >= 0.0f ? 1.0f : -1.0f;
    Vector3 n{}, t1{}, t2{};
    if (ax == 0)      { n = { sgn, 0, 0 }; t1 = { 0, 1, 0 }; t2 = { 0, 0, 1 }; }
    else if (ax == 1) { n = { 0, sgn, 0 }; t1 = { 0, 0, 1 }; t2 = { 1, 0, 0 }; }
    else              { n = { 0, 0, sgn }; t1 = { 1, 0, 0 }; t2 = { 0, 1, 0 }; }
    float umin = 1e30f, umax = -1e30f, vmin = 1e30f, vmax = -1e30f;
    for (int k = 0; k < 8; ++k) {
        const Vector3 q = { o.x + ((k & 1) ? ext : 0.0f), o.y + ((k & 2) ? ext : 0.0f), o.z + ((k & 4) ? ext : 0.0f) };
        const float s = Vector3DotProduct(q, n);
        if (s <= 1e-3f) return false; // reaches past the face's hemisphere (only huge blocks): the per-sample path
        const float u = Vector3DotProduct(q, t1) / s, v = Vector3DotProduct(q, t2) / s;
        umin = std::min(umin, u); umax = std::max(umax, u); vmin = std::min(vmin, v); vmax = std::max(vmax, v);
    }
    // Spaced at the footprint (>= 1 m) at the face centre, not the stride: nodes finer than 1 m (for edits)
    // interpolate a height map that has no finer detail anyway. Gnomonic spacing tightens toward the face edges.
    const double du = (double)std::max(stride, fp) / radius;
    const long long i0 = (long long)floor(umin / du) - 1, j0 = (long long)floor(vmin / du) - 1;
    const int nu = (int)std::min<long long>(256, (long long)ceil(umax / du) + 2 - i0), nv = (int)std::min<long long>(256, (long long)ceil(vmax / du) + 2 - j0);
    umin = (float)(i0 * du); vmin = (float)(j0 * du);
    const float su = (float)du, sv = (float)du;
    std::vector<float>& h = g.h;
    std::vector<Column3D>& cols = g.cols;
    h.assign((size_t)nu * nv, 0.0f);
    cols.assign((size_t)nu * nv, Column3D{});
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) {
            // Shared with the other blocks on this lattice (stacked blocks of a column read the same nodes).
            GridNodeCache::Key key{ i0 + i, j0 + j, 0u, ax * 2 + (sgn > 0.0f ? 1 : 0) };
            memcpy(&key.fp, &fp, sizeof(float));
            GridNodeCache::Shard& sh = gridCache->shardOf(key);
            {
                std::lock_guard<std::mutex> lock(sh.m);
                if (auto it = sh.map.find(key); it != sh.map.end()) { h[(size_t)j * nu + i] = it->second.h; cols[(size_t)j * nu + i] = it->second.col; continue; }
            }
            const float u = (float)((i0 + i) * du), v = (float)((j0 + j) * du); // from the global index: identical in every block
            const Vector3 gdir = Vector3Normalize(Vector3Add(n, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v))));
            const float e = elevation(gdir, fp);
            Column3D col = terrain3DColumn(gdir, e); // the 3D band's column strength, once per node
            col.wet = caveWetColumn(gdir, e); // (#971)
            h[(size_t)j * nu + i] = e;
            cols[(size_t)j * nu + i] = col;
            std::lock_guard<std::mutex> lock(sh.m);
            if (sh.map.size() >= GridNodeCache::kShardCap) sh.map.clear();
            sh.map.emplace(key, GridNodeCache::Node{ e, col });
        }
    g.n = n; g.t1 = t1; g.t2 = t2;
    g.umin = umin; g.vmin = vmin; g.su = su; g.sv = sv; g.nu = nu; g.nv = nv;
    return true;
}


bool PlanetField::gridSurface(const BlockGrid& g, Vector3 p, float& elev, Column3D& col) const {
    const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    if (g.nu < 2 || g.nv < 2 || d < radius + minElevation - 4.0f || d > radius + maxElevation + 4.0f) return false;
    const float s = Vector3DotProduct(p, g.n);
    if (s <= 0.0f) return false;
    const float fu = std::clamp((Vector3DotProduct(p, g.t1) / s - g.umin) / g.su, 0.0f, (float)(g.nu - 1) - 1e-3f);
    const float fv = std::clamp((Vector3DotProduct(p, g.t2) / s - g.vmin) / g.sv, 0.0f, (float)(g.nv - 1) - 1e-3f);
    const int iu = (int)fu, iv = (int)fv;
    const float au = fu - iu, av = fv - iv;
    const float* r0 = &g.h[(size_t)iv * g.nu + iu];
    const float* r1 = r0 + g.nu;
    elev = (r0[0] * (1 - au) + r0[1] * au) * (1 - av) + (r1[0] * (1 - au) + r1[1] * au) * av;
    const Column3D* c0 = &g.cols[(size_t)iv * g.nu + iu];
    const Column3D* c1 = c0 + g.nu;
    col.mask = (c0[0].mask * (1 - au) + c0[1].mask * au) * (1 - av) + (c1[0].mask * (1 - au) + c1[1].mask * au) * av;
    col.wet = (c0[0].wet * (1 - au) + c0[1].wet * au) * (1 - av) + (c1[0].wet * (1 - au) + c1[1].wet * au) * av; // (#971) as the GPU
    col.f = (c0[0].f * (1 - au) + c0[1].f * au) * (1 - av) + (c1[0].f * (1 - au) + c1[1].f * au) * av;
    return true;
}

void PlanetField::sampleBlock(Vector3 o, float stride, int S, float fp, const BlockGrid& g, float* out) const {
    std::vector<float> caveSegs; // (#615 step 2) the block's cave segments: the same list the GPU gets (blockCaveSegs)
    if (CavesGraph()) blockCaveSegs(o, stride, S, fp, caveSegs, INT32_MAX);
    CaveSegScope caveScope(CavesGraph() ? &caveSegs : nullptr);
    const float outer = radius + maxElevation, inner = radius + minElevation;
    const Vector3 n = g.n, t1 = g.t1, t2 = g.t2;
    const float umin = g.umin, vmin = g.vmin, su = g.su, sv = g.sv;
    const int nu = g.nu, nv = g.nv;
    const std::vector<float>& h = g.h;
    const std::vector<Column3D>& cols = g.cols;
    size_t idx = 0;
    for (int z = 0; z < S; ++z)
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x, ++idx) {
                const Vector3 p = { o.x + x * stride, o.y + y * stride, o.z + z * stride };
                const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
                if (d > outer + 4.0f) { out[idx] = ISO_LEVEL + DENSITY_GRADIENT * (outer - d); continue; }
                float surface = inner;
                Column3D col;
                if (d >= inner - 4.0f) {
                    const float s = Vector3DotProduct(p, n);
                    const float fu = std::clamp((Vector3DotProduct(p, t1) / s - umin) / su, 0.0f, (float)(nu - 1) - 1e-3f);
                    const float fv = std::clamp((Vector3DotProduct(p, t2) / s - vmin) / sv, 0.0f, (float)(nv - 1) - 1e-3f);
                    const int iu = (int)fu, iv = (int)fv;
                    const float au = fu - iu, av = fv - iv;
                    const float* r0 = &h[(size_t)iv * nu + iu];
                    const float* r1 = r0 + nu;
                    surface = radius + (r0[0] * (1 - au) + r0[1] * au) * (1 - av) + (r1[0] * (1 - au) + r1[1] * au) * av;
                    const Column3D* c0 = &cols[(size_t)iv * nu + iu];
                    const Column3D* c1 = c0 + nu;
                    col.mask = (c0[0].mask * (1 - au) + c0[1].mask * au) * (1 - av) + (c1[0].mask * (1 - au) + c1[1].mask * au) * av;
                    col.wet = (c0[0].wet * (1 - au) + c0[1].wet * au) * (1 - av) + (c1[0].wet * (1 - au) + c1[1].wet * au) * av; // (#971) as the GPU
                    col.f = (c0[0].f * (1 - au) + c0[1].f * au) * (1 - av) + (c1[0].f * (1 - au) + c1[1].f * au) * av;
                }
                out[idx] = densityNear(p, d, surface, fp, &col);
            }
}

// ---------------------------------------------------------------------------------------------
// Biomes (Plan 5 B)
// ---------------------------------------------------------------------------------------------

const char* BiomeName(Biome b) {
    static const char* names[(int)Biome::Count] = { "Ocean", "Jungle", "Savanna", "Desert", "Badlands", "Forest", "Meadow", "Swamp", "Taiga", "Tundra", "Alpine", "Polar" };
    return (int)b < (int)Biome::Count ? names[(int)b] : "?";
}

PlanetField::Climate PlanetField::climateAt(Vector3 dir, float elev) const {
    Climate c;
    if (moon) { c.temp = -20.0f; c.moist = 0.0f; return c; }
    const float lat = fabsf(dir.y);
    const float bs = 1.0f / std::max(gen.biomeSize, 0.1f); // #504: larger biomes = the climate noise stretched
    const float x = dir.x * radius * bs, y = dir.y * radius * bs, z = dir.z * radius * bs;
    // Warm at the equator, cold at the poles; 6.5 C colder per km of height; regional patches +-7 C.
    c.temp = 30.0f - 44.0f * (lat * sqrtf(lat)) - 6.5f * std::max(elev, 0.0f) / 1000.0f + 7.0f * tempNoise.fractal(2, x, y, z);
    // Moisture: regional wet and dry patches, wetter in the lowlands near the sea.
    c.moist = 0.5f + 0.85f * moistNoise.fractal(3, x, y, z) + 0.08f * (1.0f - Smoothstep(10.0f, 150.0f, elev)); // noise sits mostly in +-0.4
    c.moist += 0.15f * (gen.waterCycle - 1.0f); // #637: a busier water cycle, wetter climates (biomes follow)
    c.moist = std::clamp(c.moist, 0.0f, 1.0f);
    return c;
}

Biome PlanetField::biomeOf(const Climate& c, float elev) const {
    if (moon) return Biome::Polar;
    if (elev < 0.0f) return Biome::Ocean;
    if (c.temp < -12.0f) return Biome::Polar;
    if (elev > 850.0f) return Biome::Alpine;
    if (c.temp < -2.0f) return Biome::Tundra;
    if (c.temp < 7.0f) return c.moist > 0.3f ? Biome::Taiga : Biome::Tundra;
    if (c.temp < 19.0f) {
        if (c.moist < 0.33f) return Biome::Meadow;
        if (c.moist > 0.74f && elev < 90.0f) return Biome::Swamp;
        return Biome::Forest;
    }
    if (c.moist < 0.24f) return Biome::Desert;
    if (c.moist < 0.37f) return Biome::Badlands;
    if (c.moist < 0.58f) return Biome::Savanna;
    if (c.moist > 0.8f && elev < 70.0f) return Biome::Swamp;
    return Biome::Jungle;
}

Biome PlanetField::biomeAt(Vector3 dir, float elev) const {
    Climate c = climateAt(dir, elev);
    // A fine wobble (~90 m) on the climate: neighbouring biomes interleave over ~100-150 m at a border.
    const float j = biomeJitter.fractal(2, dir.x * radius, dir.y * radius, dir.z * radius);
    c.temp += 1.2f * j;
    c.moist += 0.035f * j;
    return biomeOf(c, elev);
}

Vector3 PlanetField::findSpawnDir(unsigned int seed) const {
    std::mt19937 rng(seed ^ 0x5eedu);
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    Vector3 best = { 0, 1, 0 };
    float bestScore = -1e30f;
    for (int attempt = 0; attempt < 4000; ++attempt) {
        // Temperate band: 25-45 degrees north or south.
        const float lat = (25.0f + 20.0f * unit(rng)) * DEG2RAD * (unit(rng) < 0.5f ? 1.0f : -1.0f);
        const float lon = unit(rng) * 2.0f * PI;
        const Vector3 dir = { cosf(lat) * cosf(lon), sinf(lat), cosf(lat) * sinf(lon) };
        const float e = elevation(dir);
        if (moon || e < 8.0f || e > 400.0f) continue;
        // Slope from the height a few metres away in two tangent directions.
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, { 0, 1, 0 }));
        const Vector3 t2 = Vector3CrossProduct(t1, dir);
        const float step = 4.0f, ang = step / (radius + e);
        const float e1 = elevation(Vector3Normalize(Vector3Add(dir, Vector3Scale(t1, ang))));
        const float e2 = elevation(Vector3Normalize(Vector3Add(dir, Vector3Scale(t2, ang))));
        const Vector3 normal = Vector3Normalize(Vector3Subtract(dir, Vector3Add(Vector3Scale(t1, (e1 - e) / step), Vector3Scale(t2, (e2 - e) / step))));
        const float slope = 1.0f - Vector3DotProduct(normal, dir);
        if (slope > 0.05f) continue;
        const Vector3 surface = Vector3Scale(dir, radius + e - 0.3f);
        if (materialLod(surface, 1.0f, 0.0f, normal) != VoxelMaterial::Grass) continue;
        // Solid underfoot (caves and ravines carve the density, not the height map).
        bool solid = true;
        for (float depth : { 1.0f, 3.0f, 6.0f, 12.0f })
            if (baseDensity(Vector3Scale(dir, radius + e - depth)) < ISO_LEVEL) { solid = false; break; }
        if (!solid) continue;
        // Prefer flat, low ground: good enough is good enough.
        // Start in a temperate forest or meadow (Plan 5 B); other biomes only when none is found.
        const Biome b = biomeAt(dir, e);
        const bool temperate = b == Biome::Forest || b == Biome::Meadow;
        const float score = -slope * 100.0f - e * 0.01f - (temperate ? 0.0f : 1000.0f);
        if (score > bestScore) { bestScore = score; best = dir; }
        if (slope < 0.01f && temperate) break;
    }
    if (gSpawnPond && !moon) { // bench override: the rim of the nearest pond within 6 km (cached per field)
        static std::mutex m; static std::unordered_map<const PlanetField*, Vector3> cache;
        std::lock_guard<std::mutex> lock(m);
        auto it = cache.find(this);
        if (it == cache.end()) {
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(best, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, best);
            float bestDist = 1e30f; Vector3 rim = best;
            for (int pass = 0; pass < 2 && bestDist > 1e29f; ++pass) // 6 km at 24 m; then 20 km at 60 m (#360)
                for (float e = pass ? -20000.0f : -6000.0f; e <= (pass ? 20000.0f : 6000.0f); e += pass ? 60.0f : 24.0f)
                for (float nn = pass ? -20000.0f : -6000.0f; nn <= (pass ? 20000.0f : 6000.0f); nn += pass ? 60.0f : 24.0f) {
                    const float dist = sqrtf(e * e + nn * nn);
                    if (dist >= bestDist) continue;
                    const Vector3 d = Vector3Normalize(Vector3Add(best, Vector3Add(Vector3Scale(t1, e / radius), Vector3Scale(t2, nn / radius))));
                    const WaterSource ws = waterSourceAt(d);
                    if (ws.kind != 3) continue;
                    bestDist = dist;
                    const Vector3 out = Vector3Normalize(Vector3CrossProduct(ws.pondDir, { 0, 1, 0 }));
                    rim = getenv("RAYLIBGAME_POND_CENTRE") ? ws.pondDir : Vector3Normalize(Vector3Add(ws.pondDir, Vector3Scale(out, ws.pondR / radius))); // (centre: #367 checks)
                }
            it = cache.emplace(this, rim).first;
        }
        return it->second;
    }
    if (gSpawnBiome >= 0 && gSpawnBiome < (int)Biome::Count) { // bench override: start in that biome (cached per field)
        static std::mutex m; static std::unordered_map<const PlanetField*, Vector3> cache;
        std::lock_guard<std::mutex> lock(m);
        auto it = cache.find(this);
        if (it == cache.end()) it = cache.emplace(this, nearestBiomeDir(best, (Biome)gSpawnBiome)).first;
        return it->second;
    }
    return best;
}

// The nearest deep patch of biome b (the same biome 200 m out on four sides) to home, from 60k random directions;
// home if none. Used by --bench-biome (gSpawnBiome).
Vector3 PlanetField::nearestBiomeDir(Vector3 home, Biome b) const {
    std::mt19937 rng(11);
    std::uniform_real_distribution<float> u(-1.0f, 1.0f);
    float best = -2.0f;
    Vector3 out = home;
    for (int i = 0; i < 60000; ++i) {
        Vector3 d = { u(rng), u(rng), u(rng) };
        if (Vector3LengthSqr(d) > 1.0f || Vector3LengthSqr(d) < 1e-4f) continue;
        d = Vector3Normalize(d);
        const float close = Vector3DotProduct(d, home);
        if (close <= best) continue;
        const float e = elevation(d, 32.0f);
        if (e < 6.0f || biomeAt(d, e) != b) continue;
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(d, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, d);
        bool deep = true;
        for (int s = 0; s < 4 && deep; ++s) {
            const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(s < 2 ? t1 : t2, (s & 1 ? 200.0f : -200.0f) / radius)));
            deep = biomeAt(q, elevation(q, 32.0f)) == b;
        }
        if (deep) { best = close; out = d; }
    }
    return out;
}

// The surface material of a stratum's rock (#587; soft rock is plain rock).
static VoxelMaterial RockMaterialOf(RockType r) {
    switch (r) {
        case RockType::Granite: return VoxelMaterial::Granite;
        case RockType::Basalt: return VoxelMaterial::Basalt;
        case RockType::Sandstone: return VoxelMaterial::Sandstone;
        case RockType::Limestone: return VoxelMaterial::Limestone;
        case RockType::Marble: return VoxelMaterial::Marble;
        default: return VoxelMaterial::Rock;
    }
}

VoxelMaterial PlanetField::materialSkyHinted(Vector3 p, float base, float fp, Vector3 normal, float& sky, const float* elevHint, const Column3D* colHint) const {
    sky = 1.0f;
    const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    if (d < mantleRadius + 0.5f || d < 1e-3f) { sky = 0.0f; return VoxelMaterial::Rock; } // mantle: always deep underground
    // Depth below the original surface, from the height map (the density is no measure of depth
    // here: cave tunnels take the minimum with it, so cave walls stay rock rather than topsoil).
    // `base` (the unedited density) tells undug ground (near the iso level) from dug-out ground.
    const Vector3 up = { p.x / d, p.y / d, p.z / d };
    const float surfaceElevation = elevHint ? *elevHint : elevation(up, fp);
    // Against the surface this LOD meshes, moved by the 3D terrain band (a carved alcove 10 m under the
    // height map is open air, not buried rock: it keeps its sky light and surface materials).
    const float depth = radius + surfaceElevation - d + terrain3D(p, d, radius + surfaceElevation, fp, colHint);
    // A coarse LOD mesh strays from the height map it approximates by up to a couple of samples, so
    // its "depth" must allow for that, or whole patches of surface read as buried (black, no sun) and
    // take underground materials (seen from orbit as black triangles on peaks and rings along the grid).
    // On steep ground a vertex a sample off sideways stands up to (footprint x gradient) under the height map: the slack
    // grows as the surface tilts (#541: steep rock walls at mid range read as half-buried - dark blotches, no sky).
    const float upDot = std::max(0.25f, Vector3DotProduct(normal, up));
    const float slack = fp > 1.0f ? 2.0f * fp / upDot : 0.0f;
    sky = std::clamp(1.0f - (depth - 2.0f - slack) / 8.0f, 0.0f, 1.0f); // buried more than ~10 m: no sky (caves, deep pits)
    if (moon) return depth < 2.5f + slack ? VoxelMaterial::Regolith : VoxelMaterial::Rock; // grey dust over rock (asteroid palette)
    // Bare rock is its stratum's (#587): granite, basalt, sandstone, limestone, marble, or plain (soft) rock.
    auto rockHere = [&] { return RockMaterialOf(rockTypeAt(p)); };
    if (depth < -15.0f) return Vector3DotProduct(normal, up) > 0.35f ? VoxelMaterial::Grass : rockHere(); // floating islands: turf on top, rock beneath
    // Outcrops and boulders are bare rock (snow still caps them up high, below).
    if (depth < 2.5f + slack && surfaceElevation < 1050.0f && rockFeature(p, d, radius + surfaceElevation, fp) > -0.35f) return rockHere();
    if (depth < 2.5f + slack) {
        // Surface biome from height above sea level, latitude, slope and a climate wobble. The top
        // ~0.6 m is the surface material; below it a soil layer (dirt under grass, sand under sand,
        // gravel under scree), then rock.
        const float wobble = climate.fractal(2, p.x, p.y, p.z);
        const float lat = fabsf(up.y);                                     // 0 at the equator, 1 at the poles
        const float slope = 1.0f - std::clamp(Vector3DotProduct(normal, up), -1.0f, 1.0f); // 0 flat .. 1 vertical
        const float e = surfaceElevation;
        const bool skin = depth < 0.6f + slack;
        // Snow (#548, user 2026-11-16): today's altitude, varied along the slopes - shaded (poleward-facing) slopes hold it
        // ~120 m lower, sunny ones lose it ~120 m higher, steep faces shed it (up to +300 m), and a fine noise makes the
        // edge ragged. Ice only in the coldest hollows: high, flattish and shaded (cirque floors), plus the polar caps.
        Vector3 pole = { -up.x * up.y, 1.0f - up.y * up.y, -up.z * up.y };   // the tangent toward the nearer pole
        const float pl = Vector3Length(pole);
        pole = pl > 1e-4f ? Vector3Scale(pole, (up.y >= 0.0f ? 1.0f : -1.0f) / pl) : Vector3{ 0, 0, 0 };
        const float shade = std::clamp(Vector3DotProduct(normal, pole) * 3.0f, -1.0f, 1.0f); // +1 facing the pole
        const float ragged = 45.0f * climate.fractal(2, p.x * 25.0f, p.y * 25.0f, p.z * 25.0f);
        const float snowLine = 1050.0f + 150.0f * wobble - 120.0f * shade + std::min(300.0f, 900.0f * slope) + ragged;
        if (lat > 0.93f + 0.03f * wobble) return VoxelMaterial::Ice;                                         // polar cap
        if (e > snowLine + 250.0f && slope < 0.06f && shade > -0.2f) return VoxelMaterial::Ice;              // cold hollows
        if (lat > 0.82f + 0.08f * wobble + 0.02f * ragged / 45.0f || e > snowLine) return skin ? VoxelMaterial::Snow : rockHere(); // rock under the snow (ice only in the hollows)
        // Inland water (#609, playtest 8: "rivers should replace existing terrain materials in the shape that already exists"):
        // a material swap on the ground as it lies, by how wet it is - its height over the water, its distance from the
        // water's edge (scaled to the river: a trickle wets almost nothing) and how fast the reach runs. No band of a fixed
        // width: every edge is a threshold on a smooth wetness dithered by a fine (~1.5 m) and a coarse (~15 m) noise, so
        // the per-vertex materials blend into ragged, mottled edges instead of a line. Beds by flow speed (#595): river
        // stones in fast reaches, sand with gravel bars at a moderate pace, silt and sand in slow ones, clay under the skin
        // (#206: digging a river bed finds clay); lakes and ponds keep clay with gravel bars in the deeper middle.
        constexpr float kWetReach = 60.0f; // m: floodplains and wet shores at most this far out
        // Under sea level only a bed with inland water over it (#940: a river reaching the coast through a swamp, its bed
        // a metre under the sea's level, read as beach sand - no clay to dig); the sea floor by a river mouth stays the sea's.
        auto inlandBedBelowSea = [&] { if (e <= -12.0f || !inlandWaterNear(up, 0.0f)) return false; const float w0 = waterElevation(up, fp); return w0 != kNoWater && w0 > 0.0f && e < w0; };
        if ((e > 0.0f || inlandBedBelowSea()) && inlandWaterNear(up, kWetReach)) {
            const RiverNetwork::NaturalScope natural; // a river dried by a dam keeps its bed and banks (#656)
            const RiverNetwork* net = rivers();
            const RiverNetwork::Hit rh = net && fp < 128.0f ? net->nearest(up, kWetReach + 50.0f) : RiverNetwork::Hit{};
            const float w = waterElevation(up, fp);
            const float fine = rocks.fractal(2, p.x * 14.0f + 97.0f, p.y * 14.0f, p.z * 14.0f);   // ~1.5 m mottling
            const float patch = rocks.fractal(2, p.x * 1.5f + 311.0f, p.y * 1.5f, p.z * 1.5f);    // ~15 m patches
            const float mid = rocks.fractal(2, p.x * 3.7f + 53.0f, p.y * 3.7f, p.z * 3.7f);       // ~6 m
            // The edges' threshold noise: mostly the ~6 m one, so an edge is a smooth wandering curve at the vertices (the
            // 1.5 m one alone flipped single vertices: stair steps); the shader's height blend frays it at the pixel scale.
            const float dither = 0.25f * fine + 0.5f * mid + 0.25f * patch;
            const float flat = Vector3DotProduct(normal, up);
            // A bed's skin is thicker than the turf's: a coarse mesh across a rapid's step strays up to ~1 m under the height map, and
            // its lower vertices read as buried clay (round terracotta patches in the bed). Digging a bed still finds clay.
            const bool bedSkin = depth < 1.2f + slack;
            // Speed class of the reach: 0 fast (stones), 1 moderate (sand, gravel), 2 slow (silt, mud).
            const int speed = !rh.found ? 2 : rh.slope > 0.015f ? 0 : rh.slope > 0.003f ? 1 : 2;
            // The ground right at the waterline: cobbles by fast water, sand at a moderate pace, mud with sand by slow water.
            auto waterline = [&]() {
                if (!bedSkin) return speed == 0 ? VoxelMaterial::Gravel : VoxelMaterial::Clay;
                if (speed == 0) return fine > 0.4f ? VoxelMaterial::Gravel : VoxelMaterial::Pebbles;
                if (speed == 1) return fine > 0.3f ? VoxelMaterial::Gravel : VoxelMaterial::Sand;
                return patch + 0.5f * fine > 0.15f ? VoxelMaterial::Sand : VoxelMaterial::Mud;
            };
            // The river's size: 0 for a trickle (half width 0.8 m at a head), 1 from ~6 m.
            const float size = rh.found ? std::clamp((rh.halfWidth - 0.8f) / 5.2f, 0.0f, 1.0f) : 1.0f;
            // How much beach the bank has here (0..1): patches along the river (grass meets the water in some, a strand in
            // others), more by slow water, less by fast; on gentle ground the same height band spreads wider by itself.
            const float beach = std::clamp(0.45f + 1.3f * patch + 0.5f * mid + (speed == 2 ? 0.25f : speed == 0 ? -0.25f : 0.0f), 0.0f, 1.0f);
            bool keepCover = false; // a shallow margin that keeps the ground's own cover (turf to the water's edge)
            // In a lake's water a river's channel is the lake's floor (#634: river stones and silt showed on lake floors where a
            // stream enters or leaves; its pieces cross the lake at the lake's level).
            float lakeLv = 0.0f, lakeOut = 0.0f;
            const bool inLake = net && fp < 1024.0f && net->lakeAt(up, lakeLv, lakeOut) && lakeOut < RiverNetwork::kLakeMargin;
            if (w != kNoWater && e < w) { // under the water (ground just over it is the bank's, below: its beach varies)
                const float deep = w - e; // the water's depth over this ground
                const bool river = !inLake && rh.found && rh.dist < rh.halfWidth * 1.05f && fabsf(rh.level - w) < 1.5f; // (a pond in a river stands a little over it)
                // The shallow margin (up to ~0.45 m deep, by the beach and ragged) is the waterline's ground - by a small
                // stream mostly the meadow's own turf, by a river where it has no beach - and the bed beyond it by the speed.
                if (deep < 0.06f + 0.4f * beach * (0.6f + dither)) {
                    keepCover = river && skin && (dither + 0.5f > 0.3f + 1.4f * size || beach < 0.2f + 0.3f * dither);
                    if (!keepCover) return waterline();
                } else {
                    if (!bedSkin) return VoxelMaterial::Clay;
                    if (river && speed == 0) return fine > 0.45f ? VoxelMaterial::Gravel : VoxelMaterial::Pebbles;
                    if (river && speed == 1) return deep > 0.6f && patch + 0.4f * fine > 0.2f ? VoxelMaterial::Gravel : VoxelMaterial::Sand;
                    if (river) return patch + 0.4f * fine > 0.3f ? VoxelMaterial::Sand : VoxelMaterial::Silt;
                    return deep > 1.4f && patch > 0.25f ? VoxelMaterial::Gravel : VoxelMaterial::Clay;
                }
            }
            if (!keepCover && !inLake && rh.found && rh.dist > rh.halfWidth * 0.9f) {
                // Beside a river.
                const float edge = std::max(0.0f, rh.dist - rh.halfWidth); // m out from the water's edge
                const float over = e - rh.level;                            // m over the water
                // The waterline: a thin wet margin, 0.3 m by a trickle to ~2.5 m by a river (none to ~3.5 m by the beach), only
                // near the water's own height.
                const float wl = (0.3f + 2.2f * size) * (0.1f + 1.3f * beach);
                const float wet = (1.0f - Smoothstep(0.0f, wl, edge)) * (1.0f - Smoothstep(0.1f, (0.5f + 0.4f * size) * (0.3f + 1.2f * beach), over));
                if (over > -0.5f && wet > 0.5f + 0.7f * dither) return waterline();
                // Floodplain silt: only by real rivers, on flat ground barely over the water, in ragged patches thinning
                // out with the distance and the height; none by streams.
                const float reach = size * size * (rh.halfWidth + 15.0f);
                const float high = 0.3f + 1.2f * size;
                const float plain = (1.0f - Smoothstep(0.3f * reach, reach, edge)) * (1.0f - Smoothstep(0.4f * high, high, over));
                if (skin && flat > 0.97f && reach > 1.0f && over > -0.3f && plain > 0.55f + 0.5f * patch + 0.4f * fine) return VoxelMaterial::Silt;
            }
            // Lake and pond shores: wet mud and sand by the height over the level, ragged.
            if (skin && flat > 0.9f) {
                float lv, lo;
                if (net && fp < 1024.0f && net->lakeAt(up, lv, lo) && lo < 2.0f * RiverNetwork::kLakeMargin && e > lv - 0.5f) {
                    const float wet = 1.0f - Smoothstep(0.1f, 0.9f, e - lv);
                    if (wet > 0.5f + 0.7f * dither) return patch + 0.5f * fine > 0.2f ? VoxelMaterial::Sand : VoxelMaterial::Mud;
                }
                float pondWet = 0.0f;
                if (fp < 96.0f)
                    ForPondCellsNear(Vector3Scale(up, radius), kPondReach, [&](int x, int y, int z) {
                        const Pond pd = pondOfCell(x, y, z);
                        if (!pd.exists) return;
                        const float dist = Vector3Distance(up, pd.dir) * radius;
                        if (dist > pd.r * 0.9f && dist < pd.r + 6.0f && e > pd.level - 0.4f)
                            pondWet = std::max(pondWet, 1.0f - Smoothstep(0.1f, 0.8f, e - pd.level));
                    });
                if (pondWet > 0.5f + 0.7f * dither) return VoxelMaterial::Mud;
            }
        }
        if (slope > 0.45f) return rockHere();                                                               // cliffs: bare stone
        if (e < -80.0f) return VoxelMaterial::Gravel;                                                        // deep sea floor
        // Beaches, shallows: a ragged edge (#609: the smooth climate wobble alone drew the sand's edge as a line).
        if (e < 7.5f + 4.0f * wobble && e < 6.0f + 4.0f * wobble + 1.5f * rocks.fractal(2, p.x * 6.0f + 53.0f, p.y * 6.0f, p.z * 6.0f)) return VoxelMaterial::Sand;
        // Biomes (Plan 5 B). The terrain shader tints each by the climate (red sand in the badlands,
        // dry grass on the savanna, dark mud in swamps).
        const Biome biome = biomeAt(up, e);
        if (biome == Biome::Desert) return VoxelMaterial::Sand;                                              // dunes all the way down
        if (biome == Biome::Badlands) return slope > 0.18f ? rockHere() : skin ? VoxelMaterial::Sand : rockHere(); // red sand over banded rock
        if (slope > 0.3f || e > 800.0f + 100.0f * wobble) return VoxelMaterial::Gravel;                     // scree, alpine
        if (biome == Biome::Tundra) return skin ? VoxelMaterial::Gravel : VoxelMaterial::Regolith;          // tundra: stony ground over dirt
        if (biome == Biome::Swamp) { // mud flats between the reeds: clay (#206), and clay under the swamp's turf
            if (skin && rocks.fractal(1, p.x * 0.12f, p.y * 0.12f, p.z * 0.12f) > 0.15f + 0.2f * rocks.fractal(2, p.x * 14.0f + 97.0f, p.y * 14.0f, p.z * 14.0f)) return VoxelMaterial::Clay; // ragged edges (#609)
            if (!skin) return VoxelMaterial::Clay;
        }
        // Grassland. Grass is a state of the generated surface (#1085): the band within ~1 m of the density surface, and the
        // height map's 0.6 m skin, read as turf (the mesh strays from the height map, which used to bare patches of dirt); an
        // edit that exposes or makes ground writes bare dirt over it in the material channel (TerrainDeltas::apply), so a dig
        // never shows turf. (The skin stays: the spawn search, trees and other height-map samplers read it - without it seed
        // 42's start moved to frozen tundra.)
        const bool natural = base < ISO_LEVEL + DENSITY_GRADIENT * (1.0f + slack);
        if (natural && slope < 0.3f) return VoxelMaterial::Grass;
        if (skin) return VoxelMaterial::Grass;
        // Rare small clay pockets in the dirt (#206: a dry start can still find clay by digging), a metre or two across.
        if (rocks.fractal(2, p.x * 0.35f + 731.0f, p.y * 0.35f, p.z * 0.35f) > 0.42f) return VoxelMaterial::Clay;
        return VoxelMaterial::Regolith;                                                                     // grassland over dirt
    }
    // Ore: deposits of one kind, and scattered pockets of the common ores (oreKind).
    if (depth > 3.0f && oreKind(p, depth) >= 0) return VoxelMaterial::Ore;
    return rockHere();
}

Vector3 PlanetField::fallbackUp(Vector3 p) const {
    const float pl = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    return pl > 1e-6f ? Vector3{ p.x / pl, p.y / pl, p.z / pl } : Vector3{ 0, 1, 0 };
}

bool PlanetField::editable(Vector3 p) const {
    return p.x * p.x + p.y * p.y + p.z * p.z > mantleRadius * mantleRadius;
}

// ---------------------------------------------------------------------------------------------
// Chunk meshing (worker thread)
// ---------------------------------------------------------------------------------------------

ChunkBlock ChunkBlockOf(const NodeKey& key) {
    // Sample one extra layer on every side (E = S + 2 per axis): normals are central differences of
    // this grid, so they need no extra field evaluations, and neighbouring chunks see the same
    // samples on their shared face (consistent lighting across seams).
    ChunkBlock b;
    b.stride = (float)key.stride();
    b.E = CHUNK_SIZE + 4;
    b.origin = Vector3AddValue(key.origin(), -b.stride);
    // Footprint at least 1 m: the generated terrain keeps no detail finer than it had on the 1 m lattice
    // (finer nodes exist for edits: the dig spheres), and sub-metre octaves would only cost.
    b.footprint = std::max(b.stride, 1.0f);
    return b;
}

void BuildChunkMesh(const DensityField& field, const TerrainDeltas& deltas, ChunkMeshData& out) {
    const auto t0 = std::chrono::steady_clock::now();
    const ChunkBlock b = ChunkBlockOf(out.key);
    std::vector<float> base((size_t)b.E * b.E * b.E);
    // Planets: sample through the block's height grid (as baseDensityBlock does) and keep it for the materials (#148).
    PlanetField::BlockGrid grid;
    const auto* pf = dynamic_cast<const PlanetField*>(&field);
    const bool gridded = pf && pf->prepareBlock(b.origin, b.stride, b.E, b.footprint, grid);
    if (gridded) pf->sampleBlock(b.origin, b.stride, b.E, b.footprint, grid, base.data());
    else field.baseDensityBlock(b.origin, b.stride, b.E, b.footprint, base.data());
    BuildChunkMeshFromBase(field, deltas, out, base.data(), std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count(), gridded ? &grid : nullptr);
}

// Sky light under overhangs (#369, Plan 8 item 2): the depth under the height map darkens caves and pits, but ground
// under rock that stands over it (a cliff lip, a 3D-band shelf, a tor, a dug-out hollow) read as open sky. Five short
// rays (up, and four tilted 45 degrees) through the block's own samples (edits included): the open share sets the light,
// the up ray counting double, so a shelter's mouth is dim and a few metres in is black. The up ray continues past the
// block through the full field to 12 m (so block borders don't show).
// #541: the rays start half a sample off the surface along the vertex normal and read the block trilinearly. They used to
// start on the vertex and take the nearest sample: on a steep wall a point a metre above a vertex rounded to the solid
// sample behind it about half the time once samples were 1-2 m apart, so mid-range cliffs read as roofed (dark blotches),
// up to the stride where the probe stops (a ring).
static float OverheadOpenness(const DensityField& field, const TerrainDeltas& deltas, const float* ext, int E, Vector3 eo, float stride,
                              Vector3 p, Vector3 n, const PlanetField::Column3D* col, const float* elev, bool verifyOutside = false) {
    const float d = Vector3Length(p);
    if (d < 1e-3f) return 1.0f;
    const Vector3 up = Vector3Scale(p, 1.0f / d);
    p = Vector3Add(p, Vector3Scale(n, 0.5f * stride)); // clear of the vertex's own surface cell
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(up, t1);
    const Vector3 dirs[5] = { up,
        Vector3Normalize(Vector3Add(up, t1)), Vector3Normalize(Vector3Subtract(up, t1)),
        Vector3Normalize(Vector3Add(up, t2)), Vector3Normalize(Vector3Subtract(up, t2)) };
    // Past the block only where something can stand over the ground there: the 3D band (no column info:
    // assume so). On plain height-map ground nothing natural is above a surface vertex beyond its own block.
    const bool mayRoof = !col || col->mask > 0.0f;
    auto blockSolid = [&](Vector3 q, bool& inside) {
        const float fx = (q.x - eo.x) / stride, fy = (q.y - eo.y) / stride, fz = (q.z - eo.z) / stride;
        inside = fx >= 0.0f && fy >= 0.0f && fz >= 0.0f && fx <= E - 1.0f && fy <= E - 1.0f && fz <= E - 1.0f;
        if (!inside) return false;
        const int x0 = std::min((int)fx, E - 2), y0 = std::min((int)fy, E - 2), z0 = std::min((int)fz, E - 2);
        const float tx = fx - x0, ty = fy - y0, tz = fz - z0;
        auto A = [&](int x, int y, int z) { return ext[((size_t)z * E + y) * E + x]; };
        const float c00 = A(x0, y0, z0) + (A(x0 + 1, y0, z0) - A(x0, y0, z0)) * tx;
        const float c10 = A(x0, y0 + 1, z0) + (A(x0 + 1, y0 + 1, z0) - A(x0, y0 + 1, z0)) * tx;
        const float c01 = A(x0, y0, z0 + 1) + (A(x0 + 1, y0, z0 + 1) - A(x0, y0, z0 + 1)) * tx;
        const float c11 = A(x0, y0 + 1, z0 + 1) + (A(x0 + 1, y0 + 1, z0 + 1) - A(x0, y0 + 1, z0 + 1)) * tx;
        const float c0 = c00 + (c10 - c00) * ty, c1 = c01 + (c11 - c01) * ty;
        return c0 + (c1 - c0) * tz >= ISO_LEVEL;
    };
    float open = 0.0f;
    for (int r = 0; r < 5; ++r) {
        bool blocked = false, inside = true;
        // Eight distances, closer together near the vertex (a low roof matters most), from clear of its own surface.
        static const float kDist[8] = { 1.0f, 1.5f, 2.0f, 3.0f, 4.0f, 6.0f, 8.0f, 12.0f };
        int k = 0;
        for (; k < 8 && !blocked; ++k) {
            if (blockSolid(Vector3Add(p, Vector3Scale(dirs[r], kDist[k])), inside)) blocked = true;
            if (!inside) break;
        }
        if (!blocked && !inside && r == 0 && (mayRoof || verifyOutside)) { // the up ray beyond the block: the rest from the field (#898: also for a buried vertex)
            // Straight up stays in this vertex's column, so with the column known (the block's grid) the surface is too:
            // densityNear skips the height map, which made full density() calls here the bulk of a build in the 3D band
            // (no erosion: "Preparing spawn" 8 s -> 40 s, #387).
            const auto* pf = (col && elev) ? dynamic_cast<const PlanetField*>(&field) : nullptr;
            const float fp = std::max(stride, 1.0f);
            for (; k < 8 && !blocked; k += 2) {
                const Vector3 q = Vector3Add(p, Vector3Scale(up, kDist[k]));
                const float v = pf ? pf->densityInColumn(q, d + kDist[k], *elev, fp, col) + deltas.sample(q)
                                   : field.density(q, &deltas, fp);
                blocked = v >= ISO_LEVEL;
            }
        }
        if (!blocked && !inside && verifyOutside && r != 0) {
            // (#898) a buried vertex: a ray that left the block is open only if the ground beyond really is (else cave walls
            // at a block's edge would read as lit), from the field itself.
            const float fp = std::max(stride, 1.0f);
            for (; k < 8 && !blocked; k += 2) blocked = field.density(Vector3Add(p, Vector3Scale(dirs[r], kDist[k])), &deltas, fp) >= ISO_LEVEL; // (every other: the cost)
        }
        if (!blocked) open += r == 0 ? 3.0f : 1.0f; // (#899) the sky overhead weighs most
    }
    // (#899) Open ground and ground by one wall stay 1; under an overhang (only the up ray blocked) 0.71, in an alcove 0.54. With
    // the old weights (up 2, x1.5) an overhang open at the sides read as full sky.
    return std::clamp(open / 7.0f * 1.25f, 0.0f, 1.0f);
}

// ---------------------------------------------------------------------------------------------
// Water (#360): one surface for the generated water and the simulation's published cells
// ---------------------------------------------------------------------------------------------

// The generated water a chunk may hold (never stored, from the field): whether it spans the sea's radius, and whether
// inland water (rivers, lakes, ponds) may lie in it. Its levels come from the water surface map (watersurface.h, #639).
static ChunkWaterGen GeneratedWaterOf(const DensityField& field, Vector3 origin, int cells, float stride) {
    ChunkWaterGen g;
    const float sea = field.seaLevelRadius();
    if (sea <= 0.0f) return g;
    const Vector3 lo = origin, hi = Vector3AddValue(origin, cells * stride);
    const Vector3 nearest = { std::clamp(0.0f, lo.x, hi.x), std::clamp(0.0f, lo.y, hi.y), std::clamp(0.0f, lo.z, hi.z) };
    const float dmin = Vector3Length(nearest);
    const float dmax = Vector3Length({ std::max(fabsf(lo.x), fabsf(hi.x)), std::max(fabsf(lo.y), fabsf(hi.y)), std::max(fabsf(lo.z), fabsf(hi.z)) });
    const PlanetField* pf = dynamic_cast<const PlanetField*>(&field);
    const Vector3 mid = Vector3Scale(Vector3Add(lo, hi), 0.5f);
    const float extent = Vector3Distance(lo, hi) * 0.5f;
    g.seaSpan = dmin < sea && dmax > sea;
    g.inland = pf && !pf->moon && dmax > sea && dmin < sea + pf->maxElevation && Vector3Length(mid) > 1.0f &&
               pf->inlandWaterNear(Vector3Normalize(mid), extent, stride);
    return g;
}

// Every fluid surface of a chunk from the one water field: per kind, marching cubes on a density that is the
// published simulation cells' surface where a brick covers the sample (WaterField::ChunkSampler) and the generated
// water's level elsewhere - the same plane where nothing moved, so the two meet without a seam. Triangles inside
// the terrain (the chunk's own density) are dropped, so the water meets the shore at the chunk's resolution; the
// generated water also stays out of caves under it (only over the original bed), the simulated water does not
// (cave water, dug channels). The vertex's second texcoord carries the fluid kind (water.fs colours it).
static constexpr float kSkin = 0.15f; // m: water shallower than this over the ground is not drawn (MeshChunkWater)
static constexpr float kSkinMoving = 0.02f; // m: the same for moving water (#649: films are lifted ~5 cm, waterfield.cpp)
static std::atomic<long long> gWaterPhaseNs[4];
static std::atomic<long long> gWaterSamples[2]; // exact water density samples, all samples (#381)
void ChunkWaterPhases(long long out[6]) { for (int i = 0; i < 4; ++i) out[i] = gWaterPhaseNs[i].load(); out[4] = gWaterSamples[0].load(); out[5] = gWaterSamples[1].load(); }
// Moving simulated water's vertices (#677, --test-flow-shade): all, steep ones (surface 40+ degrees off level), steep ones
// drawn with a current of 1 m/s or more along their surface.
static std::atomic<long long> gSimFlowVerts[3];
// RAYLIBGAME_FLOW_SURFACE=0: the current stays level and simulated water keeps its gradient normal (before #677/#696), for A/B.
static const bool kFlowSurface = [] { const char* e = getenv("RAYLIBGAME_FLOW_SURFACE"); return !(e && e[0] == '0'); }();
void ChunkWaterFlowStats(long long out[3]) { for (int i = 0; i < 3; ++i) out[i] = gSimFlowVerts[i].load(); }
// RAYLIBGAME_WATER_COARSE=0: sample every water density point (no coarse-first skipping), for A/B.
static const bool kWaterCoarse = [] { const char* e = getenv("RAYLIBGAME_WATER_COARSE"); return !(e && e[0] == '0'); }();
static long long NowNs() { return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static void MeshChunkWater(const DensityField& field, ChunkMeshData& out, const float* dens, const PlanetField::BlockGrid* grid,
                           const ChunkWaterGen& gen, WaterField::ChunkSampler& sim) {
    const int cells = CHUNK_SIZE + 1, S = cells + 1;
    const float stride = (float)out.key.stride();
    const Vector3 origin = out.key.origin();
    const float sea = field.seaLevelRadius();
    const bool genAny = sea > 0.0f && (gen.seaSpan || gen.inland);
    const uint8_t kinds = sim.any() ? sim.kinds() : 0;
    if (!genAny && !kinds) return;
    const PlanetField* pf = dynamic_cast<const PlanetField*>(&field);
    const RiverNetwork* rivers = pf ? pf->rivers() : nullptr; // the current on river water (#361)
    // Inland levels from the water surface map at this chunk's sample spacing (#639): one shared 2D map, so a stream
    // narrower than the old per-chunk lattice (every 4 samples) is still drawn and shores don't grow teeth between nodes.
    std::unique_ptr<WaterSurface::Reader> surf;
    if (gen.inland && pf) surf = std::make_unique<WaterSurface::Reader>(*pf, stride);
    auto levelAt = [&](int x, int y, int z) {
        if (!surf) return sea;
        return sea + surf->drawLevel(Vector3Normalize({ origin.x + x * stride, origin.y + y * stride, origin.z + z * stride }));
    };
    // Trilinear over an S^3 grid in sample coordinates.
    auto trilinear = [&](const float* f, Vector3 g) {
        const int x0 = std::clamp((int)floorf(g.x), 0, S - 2), y0 = std::clamp((int)floorf(g.y), 0, S - 2), z0 = std::clamp((int)floorf(g.z), 0, S - 2);
        const float fx = std::clamp(g.x - x0, 0.0f, 1.0f), fy = std::clamp(g.y - y0, 0.0f, 1.0f), fz = std::clamp(g.z - z0, 0.0f, 1.0f);
        auto d = [&](int x, int y, int z) { return f[((size_t)z * S + y) * S + x]; };
        const float c00 = d(x0, y0, z0) * (1 - fx) + d(x0 + 1, y0, z0) * fx, c10 = d(x0, y0 + 1, z0) * (1 - fx) + d(x0 + 1, y0 + 1, z0) * fx;
        const float c01 = d(x0, y0, z0 + 1) * (1 - fx) + d(x0 + 1, y0, z0 + 1) * fx, c11 = d(x0, y0 + 1, z0 + 1) * (1 - fx) + d(x0 + 1, y0 + 1, z0 + 1) * fx;
        return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz;
    };
    // Height-map elevation under a body point: from the block's grid when there is one (#148), as the terrain.
    auto groundAt = [&](Vector3 point, Vector3 dir) {
        float e; PlanetField::Column3D col;
        if (grid && pf && pf->gridSurface(*grid, point, e, col)) return e;
        return pf->elevation(dir, stride);
    };
    std::vector<float> wd((size_t)S * S * S);
    std::vector<Vector3> wv;
    std::vector<int> wi;
    std::vector<int> cache((size_t)S * S * S * 3);
    // Coarse-first sampling (#381), level-0 chunks only: the coarse lattice lines and, per axis, the blocks a sample lies in
    // (one between two coarse lines, two on a shared line), the corner lines of the first and its fraction between them
    // (0 or 1 on a shared line, where the blocks agree).
    const bool coarse = kWaterCoarse && S >= 4 && out.key.level == 0;
    std::vector<int> cc, A0, A1, Lo, Hi, work;
    std::vector<float> T;
    std::vector<char> mixed, done; // per block: holds a surface (1) or was taken in (2); per sample: exact
    int NB = 0;
    if (coarse) {
        for (int v = 0; v < S - 1; v += 2) cc.push_back(v);
        cc.push_back(S - 1);
        NB = (int)cc.size() - 1;
        auto isC = [&](int v) { return (v & 1) == 0 || v == S - 1; };
        A0.resize(S); A1.resize(S); Lo.resize(S); Hi.resize(S); T.resize(S);
        for (int v = 0; v < S; ++v) {
            if (!isC(v)) A0[v] = A1[v] = (v - 1) / 2;
            else if (v == S - 1) A0[v] = A1[v] = NB - 1;
            else { A0[v] = std::max(0, v / 2 - 1); A1[v] = std::min(NB - 1, v / 2); }
            Lo[v] = cc[A0[v]]; Hi[v] = cc[A0[v] + 1];
            T[v] = (float)(v - Lo[v]) / (float)(Hi[v] - Lo[v]);
        }
        mixed.resize((size_t)NB * NB * NB);
        done.resize((size_t)S * S * S);
    }
    for (int kind = 0; kind < 8; ++kind) {
        const bool simKind = (kinds >> kind) & 1;
        const bool genKind = kind == 0 && genAny;
        if (!simKind && !genKind) continue;
        size_t k = 0;
        bool anyIn = false, anyOut = false;
        long long tPhase = NowNs();
        auto phase = [&](int i) { const long long t = NowNs(); gWaterPhaseNs[i] += t - tPhase; tPhase = t; };
        auto exact = [&](int x, int y, int z) {
            const Vector3 p = { origin.x + x * stride, origin.y + y * stride, origin.z + z * stride };
            float w;
            if (!(simKind && sim.density(p, kind, w)))
                w = genKind ? ISO_LEVEL + DENSITY_GRADIENT * (levelAt(x, y, z) - sqrtf(p.x * p.x + p.y * p.y + p.z * p.z)) : ISO_LEVEL - DENSITY_GRADIENT;
            return w;
        };
        long long exactCount = 0;
        if (!coarse) {
            for (int z = 0; z < S; ++z)
                for (int y = 0; y < S; ++y)
                    for (int x = 0; x < S; ++x, ++k) wd[k] = exact(x, y, z);
            exactCount = (long long)S * S * S;
        } else {
            // Coarse first (#381): the samples on every other lattice line (and the last), then the 2-sample blocks between
            // them. A block whose 8 corners are all on one side of the surface, a quarter metre or more from it, holds no
            // surface; its other samples are interpolated, not asked. Blocks with a corner near or across the surface are
            // sampled in full, with a one-sample border (the normals read half a sample around their vertices), and a block
            // one of those exact samples shows on the other side of its corners (a feature thinner than a block) is taken in
            // as well, until none is left: every cube that makes triangles has exact corners. Fine chunks only: a cell is
            // 2 blocks across there, so a simulated column (1 m) can't hide between coarse lines (#381 review).
            auto at = [&](int x, int y, int z) -> float& { return wd[((size_t)z * S + y) * S + x]; };
            std::fill(done.begin(), done.end(), (char)0);
            auto ex = [&](int x, int y, int z) {
                const size_t i = ((size_t)z * S + y) * S + x;
                if (!done[i]) { wd[i] = exact(x, y, z); done[i] = 1; ++exactCount; }
                return wd[i];
            };
            for (int z : cc) for (int y : cc) for (int x : cc) ex(x, y, z);
            constexpr float kMargin = 0.25f * DENSITY_GRADIENT; // a quarter metre of the fluid's (or the sea's) surface
            std::fill(mixed.begin(), mixed.end(), (char)0);
            for (int bz = 0; bz < NB; ++bz)
                for (int by = 0; by < NB; ++by)
                    for (int bx = 0; bx < NB; ++bx) {
                        bool in = false, outS = false, nearS = false;
                        for (int c = 0; c < 8; ++c) {
                            const float w = at(cc[bx + (c & 1)], cc[by + ((c >> 1) & 1)], cc[bz + (c >> 2)]);
                            if (w >= ISO_LEVEL) in = true; else outS = true;
                            if (fabsf(w - ISO_LEVEL) < kMargin) nearS = true;
                        }
                        if ((in && outS) || nearS) { mixed[((size_t)bz * NB + by) * NB + bx] = 1; work.push_back(((bz * NB) + by) * NB + bx); }
                    }
            while (!work.empty()) {
                const int b = work.back();
                work.pop_back();
                const int bx = b % NB, by = (b / NB) % NB, bz = b / (NB * NB);
                for (int z = std::max(0, cc[bz] - 1); z <= std::min(S - 1, cc[bz + 1] + 1); ++z)
                    for (int y = std::max(0, cc[by] - 1); y <= std::min(S - 1, cc[by + 1] + 1); ++y)
                        for (int x = std::max(0, cc[bx] - 1); x <= std::min(S - 1, cc[bx + 1] + 1); ++x) {
                            if (done[((size_t)z * S + y) * S + x]) continue;
                            const bool wet = ex(x, y, z) >= ISO_LEVEL;
                            for (int qz = A0[z]; qz <= A1[z]; ++qz)
                                for (int qy = A0[y]; qy <= A1[y]; ++qy)
                                    for (int qx = A0[x]; qx <= A1[x]; ++qx) {
                                        char& m = mixed[((size_t)qz * NB + qy) * NB + qx];
                                        if (!m && wet != (at(cc[qx], cc[qy], cc[qz]) >= ISO_LEVEL)) { m = 2; work.push_back((qz * NB + qy) * NB + qx); }
                                    }
                        }
            }
            for (int z = 0; z < S; ++z)
                for (int y = 0; y < S; ++y)
                    for (int x = 0; x < S; ++x) {
                        if (done[((size_t)z * S + y) * S + x]) continue;
                        // Trilinear in the first block holding it (its corners all on one side: so is the result).
                        const float tx = T[x], ty = T[y], tz = T[z];
                        const int xa = Lo[x], xb = Hi[x], ya = Lo[y], yb = Hi[y], za = Lo[z], zb = Hi[z];
                        const float c00 = at(xa, ya, za) * (1 - tx) + at(xb, ya, za) * tx, c10 = at(xa, yb, za) * (1 - tx) + at(xb, yb, za) * tx;
                        const float c01 = at(xa, ya, zb) * (1 - tx) + at(xb, ya, zb) * tx, c11 = at(xa, yb, zb) * (1 - tx) + at(xb, yb, zb) * tx;
                        at(x, y, z) = (c00 * (1 - ty) + c10 * ty) * (1 - tz) + (c01 * (1 - ty) + c11 * ty) * tz;
                    }
        }
        gWaterSamples[0] += exactCount; gWaterSamples[1] += (long long)S * S * S;
        for (size_t i = 0; i < (size_t)S * S * S; ++i) { if (wd[i] >= ISO_LEVEL) anyIn = true; else anyOut = true; }
        phase(0);
        if (!anyIn || !anyOut) continue;
        const size_t v0 = out.waterPositions.size() / 3, i0 = out.waterIndices.size(); // this kind's start
        wv.clear(); wi.clear();
        std::fill(cache.begin(), cache.end(), -1);
        if (coarse) {
            // Coarse sampling: only the blocks that hold a surface can make triangles (a skipped block's cubes all lie on one
            // side), so march just their cubes. The same triangles as the full sweep, in another order.
            for (int bz = 0; bz < NB; ++bz)
                for (int by = 0; by < NB; ++by)
                    for (int bx = 0; bx < NB; ++bx) {
                        if (!mixed[((size_t)bz * NB + by) * NB + bx]) continue;
                        for (int z = cc[bz]; z < cc[bz + 1]; ++z)
                            for (int y = cc[by]; y < cc[by + 1]; ++y)
                                for (int x = cc[bx]; x < cc[bx + 1]; ++x) marchCube(x, y, z, wd, S, wv, wi, cache, ISO_LEVEL);
                    }
        } else {
            for (int z = 0; z < cells; ++z)
                for (int y = 0; y < cells; ++y)
                    for (int x = 0; x < cells; ++x)
                        marchCube(x, y, z, wd, S, wv, wi, cache, ISO_LEVEL);
        }
        phase(1);
        // The shore (#423): a triangle crossing the ground keeps only its part out of the ground, cut where its edges
        // cross the terrain's iso level (linear in the terrain density), so the waterline is a continuous curve instead
        // of the edge of whichever whole triangles stayed out (the cull below dropped a triangle by its centre).
        {
            std::vector<int> clipped;
            clipped.reserve(wi.size());
            auto cut = [&](int a, int b, float da, float db) {
                const float t = da / (da - db); // da, db: terrain density - ISO, of opposite signs
                wv.push_back(Vector3Lerp(wv[a], wv[b], t));
                return (int)wv.size() - 1;
            };
            // The terrain density at each marched vertex, once (a vertex is shared by ~6 triangles; #381).
            const size_t nv = wv.size();
            std::vector<float> vd(nv);
            for (size_t i = 0; i < nv; ++i) vd[i] = trilinear(dens, wv[i]) - ISO_LEVEL;
            for (size_t t = 0; t + 2 < wi.size(); t += 3) {
                int v[3] = { wi[t], wi[t + 1], wi[t + 2] };
                float d[3];
                int in = 0;
                for (int e = 0; e < 3; ++e) { d[e] = vd[v[e]]; in += d[e] >= 0.0f ? 1 : 0; }
                if (in == 0) { clipped.insert(clipped.end(), { v[0], v[1], v[2] }); continue; }
                if (in == 3) continue;
                // Rotate (keeping the winding) so that: one in the ground -> it is first; two -> the one out is last.
                for (int r = 0; r < 3; ++r) {
                    const bool ok = in == 1 ? d[0] >= 0.0f : d[2] < 0.0f;
                    if (ok) break;
                    std::swap(v[0], v[1]); std::swap(v[1], v[2]); std::swap(d[0], d[1]); std::swap(d[1], d[2]);
                }
                if (in == 1) { // X in, Y Z out: the quad P Y Z Q
                    const int p = cut(v[0], v[1], d[0], d[1]), q = cut(v[0], v[2], d[0], d[2]);
                    clipped.insert(clipped.end(), { p, v[1], v[2], p, v[2], q });
                } else {       // X Y in, Z out: the triangle P Z Q
                    const int p = cut(v[1], v[2], d[1], d[2]), q = cut(v[2], v[0], d[2], d[0]);
                    clipped.insert(clipped.end(), { p, v[2], q });
                }
            }
            wi.swap(clipped);
        }
        phase(2);
        std::vector<int> remap(wv.size(), -1);
        for (size_t t = 0; t + 2 < wi.size(); t += 3) {
            const Vector3 c = Vector3Scale(Vector3Add(Vector3Add(wv[wi[t]], wv[wi[t + 1]]), wv[wi[t + 2]]), 1.0f / 3.0f);
            if (trilinear(dens, c) >= ISO_LEVEL) continue; // under land
            const Vector3 wc = Vector3Add(origin, Vector3Scale(c, stride));
            const bool covered = simKind && sim.covered(wc);
            if (covered) {
                // Simulated water as a skin over the ground (under ~15 cm deep): not drawn. Its plane at the ground's own
                // height broke up over every ripple of the bed into a grid of spikes (the films a drained pond leaves).
                // Generated shores keep their shallow band (#382 review). Moving water keeps a thinner skin (#649): a film running
                // over the ground is drawn (it was hidden, so the water seemed to skip to where it pooled).
                const float skin = sim.moving(wc) ? kSkinMoving : kSkin;
                const Vector3 down = Vector3Scale(Vector3Normalize(wc), -skin / stride);
                if (trilinear(dens, Vector3Add(c, down)) >= ISO_LEVEL) continue;
            } else {
                if (!genKind) continue;
                if (pf) {
                    // Underground: where the original ground stands over the water (a cave or tunnel under a river, or
                    // below sea level) the generated water surface is not open water. Rivers used to line the cave walls.
                    const float wl = Vector3Length(wc);
                    if (wl > 1.0f && groundAt(wc, Vector3Scale(wc, 1.0f / wl)) > wl - sea + (stride <= 2.0f ? 0.15f : 1.0f)) continue; // near: only over the original bed
                    // Steep: generated water lies (nearly) flat - a river's rapids slope at most ~0.4. A steep surface is the
                    // level lattice interpolated across a waterfall's drop: blobs in the rock behind the fall, and from far
                    // LODs (their lattice is coarse) a giant slanted sheet (#613). The fall's own sheet draws the water there.
                    const Vector3 tn = Vector3CrossProduct(Vector3Subtract(wv[wi[t + 1]], wv[wi[t]]), Vector3Subtract(wv[wi[t + 2]], wv[wi[t]]));
                    const float tl = Vector3Length(tn);
                    if (wl > 1.0f && tl > 1e-9f && fabsf(Vector3DotProduct(tn, wc)) < 0.7f * tl * wl) continue;
                }
            }
            for (int e = 0; e < 3; ++e) {
                int& r = remap[wi[t + e]];
                if (r < 0) {
                    r = (int)(out.waterPositions.size() / 3);
                    const Vector3 g = wv[wi[t + e]];
                    const Vector3 lp = Vector3Scale(g, stride);
                    const Vector3 bp = Vector3Add(origin, lp);
                    const Vector3 radial = Vector3Normalize(bp);
                    Vector3 n = radial; // the sea is a sphere: radial normals
                    Vector3 simFlow{}; float simSpeed = 0.0f; // simulated water's current (#677)
                    float depth = 50.0f;
                    if (simKind && sim.covered(bp)) {
                        // Moving fluid: the surface's own normal (falls, steps), the fluid under it for the tint.
                        const float h = 0.5f;
                        // Central differences kept inside the sample grid, divided by the step actually taken: at a chunk's
                        // lower faces the lookup clamped to the edge and halved that axis's difference, which tilted the
                        // normal of a surface not aligned with the axes - a shading line along chunk borders.
                        const float lo = 0.0f, hiG = (float)(S - 1);
                        auto d = [&](float c) { return std::pair<float, float>{ std::max(lo, c - h), std::min(hiG, c + h) }; };
                        const auto [x0, x1] = d(g.x); const auto [y0, y1] = d(g.y); const auto [z0, z1] = d(g.z);
                        const Vector3 grad = { (trilinear(wd.data(), { x1, g.y, g.z }) - trilinear(wd.data(), { x0, g.y, g.z })) / std::max(1e-3f, x1 - x0),
                                               (trilinear(wd.data(), { g.x, y1, g.z }) - trilinear(wd.data(), { g.x, y0, g.z })) / std::max(1e-3f, y1 - y0),
                                               (trilinear(wd.data(), { g.x, g.y, z1 }) - trilinear(wd.data(), { g.x, g.y, z0 })) / std::max(1e-3f, z1 - z0) };
                        const float gl = Vector3Length(grad);
                        if (gl > 1e-6f) n = Vector3Scale(grad, -1.0f / gl);
                        depth = std::max(0.3f, sim.depthBelow(bp, kind));
                        // Its current (#677): down the simulated surface's slope at Manning's speed, v = depth^(2/3) sqrt(slope) /
                        // n (n 0.03, at most 4 m/s): a channel runs, a level pool stands still - awake, steady or asleep alike.
                        Vector3 down; float grade;
                        Vector3 solverDir; float solverSpeed;
                        if (sim.simFlow(bp, solverDir, solverSpeed)) { // the solver's own flow, published per column (#695)
                            simFlow = Vector3Subtract(solverDir, Vector3Scale(radial, Vector3DotProduct(solverDir, radial)));
                            const float fl = Vector3Length(simFlow);
                            simFlow = fl > 1e-4f ? Vector3Scale(simFlow, 1.0f / fl) : Vector3{ 0, 0, 0 };
                            simSpeed = solverSpeed;
                        } else if (sim.slope(bp, kind, down, grade) && grade > 1e-4f) {
                            simFlow = Vector3Subtract(down, Vector3Scale(radial, Vector3DotProduct(down, radial)));
                            const float fl = Vector3Length(simFlow);
                            if (fl > 1e-4f) {
                                simFlow = Vector3Scale(simFlow, 1.0f / fl);
                                simSpeed = std::min(4.0f, powf(std::max(0.05f, sim.depthBelow(bp, kind)), 2.0f / 3.0f) * sqrtf(grade) / 0.03f);
                            }
                        }
                        // Along the water's own surface (#677): the solver's flow is level (east, north), so on a step, a chute or
                        // a falling sheet it pointed into the surface and water.fs had nothing left to draw. Moving water on a
                        // slope runs down it: tilt the current into the surface plane, towards its downhill as it steepens, and a
                        // sheet sliding down a step or a fall is fast (3 m/s at vertical) even where its column's flow is small.
                        const float upN = Vector3DotProduct(n, radial); // 1: level, 0: a vertical sheet, < 0: an underside
                        if (kFlowSurface && upN < 0.97f && upN > -0.3f && sim.moving(bp)) {
                            const Vector3 along = Vector3Subtract(simFlow, Vector3Scale(n, Vector3DotProduct(simFlow, n)));
                            Vector3 downhill = Vector3Subtract(Vector3Scale(n, upN), radial);
                            const float dl = Vector3Length(downhill);
                            if (dl > 1e-4f) {
                                downhill = Vector3Scale(downhill, 1.0f / dl);
                                const float steep = std::clamp((0.97f - upN) / 0.6f, 0.0f, 1.0f);
                                const Vector3 f = Vector3Add(Vector3Scale(along, 1.0f - steep), Vector3Scale(downhill, steep));
                                const float fl = Vector3Length(f);
                                if (fl > 1e-4f) { simFlow = Vector3Scale(f, 1.0f / fl); simSpeed = std::min(4.0f, std::max(simSpeed, 3.0f * steep)); }
                            }
                        }
                        if (sim.moving(bp)) {
                            ++gSimFlowVerts[0];
                            if (upN < 0.77f && upN > -0.3f) {
                                ++gSimFlowVerts[1];
                                const float alongN = fabsf(Vector3DotProduct(simFlow, n)); // drawn along the surface, not into it
                                if (simSpeed >= 1.0f && alongN < 0.2f) ++gSimFlowVerts[2];
                            }
                        }
                        // Shade near-level water with the radial normal, as the generated water is (#696): the fill grid's
                        // gradient tilts from cell to cell over a resting or gently running surface (its levels step), and
                        // at grazing angles that tilt cost most of the Fresnel sky reflection - disturbed water turned clear
                        // and dull along a hard line where it met the generated water. Steps and falls keep their own normal.
                        if (kFlowSurface && upN > 0.5f) {
                            const float t = std::clamp((upN - 0.5f) / 0.4f, 0.0f, 1.0f), w = t * t * (3.0f - 2.0f * t);
                            n = Vector3Normalize(Vector3Lerp(n, radial, w));
                        }
                        // Over simulated sea or lake cells: no shallower than the generated water around them (a lighter
                        // square at the edge of a site otherwise: the column count stops at 8 m, #382 review).
                        if (genKind && pf) depth = std::max(depth, Vector3Length(bp) - sea - groundAt(bp, radial));
                    } else if (pf) {
                        // Depth to the sea floor from the height map (the shader grades clarity and colour by it).
                        const float surfaceE = Vector3Length(bp) - sea; // the water surface's height above sea level
                        depth = std::max(0.0f, surfaceE - groundAt(bp, radial));
                    }
                    out.waterPositions.insert(out.waterPositions.end(), { lp.x, lp.y, lp.z });
                    out.waterNormals.insert(out.waterNormals.end(), { n.x, n.y, n.z });
                    out.waterDepth.insert(out.waterDepth.end(), { depth, (float)kind });
                    // The current (#361): generated river water flows downstream, ~1 m/s mid-channel, slower to the banks.
                    // Lakes and the sea stand still here; simulated water runs down its surface (#677, above).
                    Vector3 flow = simFlow; float speed = simSpeed;
                    // (#1121, playtest: settled rivers lost their moving-water look) Simulated water whose solver current has
                    // died down (its bricks settled or asleep: < 5 cm/s) is still a river where a river runs: it takes the river's
                    // current as the generated water does. Ponds and lakes have no river there and stay still.
                    // RAYLIBGAME_SETTLED_RIVER_FLOW=0: as before.
                    static const bool settledRiver = !(getenv("RAYLIBGAME_SETTLED_RIVER_FLOW") && getenv("RAYLIBGAME_SETTLED_RIVER_FLOW")[0] == '0');
                    const bool simHere = simKind && sim.covered(bp);
                    if (kind == 0 && gen.inland && rivers && (!simHere || (settledRiver && simSpeed < 0.05f))) {
                        const RiverNetwork::Hit hr = rivers->nearest(radial, 60.0f);
                        if (hr.found && hr.dist < hr.halfWidth + 3.0f) {
                            flow = Vector3Subtract(hr.down, Vector3Scale(radial, Vector3DotProduct(hr.down, radial))); // along the surface
                            const float fl = Vector3Length(flow);
                            if (fl > 1e-4f) { flow = Vector3Scale(flow, 1.0f / fl); speed = 1.0f * sqrtf(std::max(0.0f, 1.0f - hr.dist / (hr.halfWidth + 3.0f))); }
                        }
                    }
                    auto enc = [](float c) { return (unsigned char)std::clamp((int)lroundf((c * 0.5f + 0.5f) * 255.0f), 0, 255); };
                    out.waterFlow.insert(out.waterFlow.end(), { enc(flow.x), enc(flow.y), enc(flow.z), (unsigned char)std::clamp((int)lroundf(speed / 4.0f * 255.0f), 0, 255) });
                }
                out.waterIndices.push_back((unsigned short)r);
            }
        }
        phase(3);
        if (out.waterPositions.size() / 3 > 65535) { out.waterPositions.resize(v0 * 3); out.waterNormals.resize(v0 * 3); out.waterDepth.resize(v0 * 2); out.waterFlow.resize(v0 * 4); out.waterIndices.resize(i0); }
    }
    // (16-bit indices: a kind that would pass 65535 vertices is dropped alone, never the water before it.)
}

void BuildChunkWater(const DensityField& field, const ChunkWaterCache& cache, ChunkMeshData& out) {
    const auto t0 = std::chrono::steady_clock::now();
    const WaterField::SnapshotPtr snap = WaterField::Current();
    const Vector3 origin = out.key.origin();
    WaterField::ChunkSampler sim(snap.get(), out.planetoidId, origin, Vector3AddValue(origin, (CHUNK_SIZE + 1) * (float)out.key.stride()));
    out.waterSeq = snap->seq;
    out.waterOnly = true;
    MeshChunkWater(field, out, cache.dens.data(), cache.grid.get(), cache.gen, sim);
    out.buildMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

static std::atomic<long long> gWaterCacheBytes{ 0 };
ChunkWaterCache::ChunkWaterCache() = default;
ChunkWaterCache::~ChunkWaterCache() { gWaterCacheBytes -= counted; }
long long ChunkWaterCacheBytes() { return gWaterCacheBytes.load(); }

// #679 A/B: RAYLIBGAME_DELTA_POINT=1 = coarse nodes point-sample the edit lattice (the old way).
static const bool gDeltaPointSample = [] { const char* e = getenv("RAYLIBGAME_DELTA_POINT"); return e && e[0] == '1'; }();

// Player edits on a coarse node (#679). A carve writes its exact shape only within 1 m of it (the lattice there holds
// the carved density); farther out the field is the unedited ground again, a step of up to the dig's depth in density.
// A coarse sample a stride under a dug floor read that step, and a sample over it the ground the dig took away, so
// marching cubes put the floor anywhere between them: a flat dug floor came out lumpy by up to half a stride (a
// honeycomb under shallow water, the floor plane cutting the sample lattice at a slant) and a dig shallower than a
// stride could vanish between two samples.
// Along a vertical column the unedited density is linear (gradient times the height under the height map), so the
// edit at a point q above a sample p says what the field at p would be if q's edit reached down to it:
// base(p) + delta(q) (the floor-plane value). Each coarse sample takes, within two strides along its column (the
// lattice axis nearest the vertical): the most negative such dig delta - from above only while p is solid (a dig
// below solid ground must not open it), from either side once p is air (it only gets more air) - and the most
// positive fill delta, mirrored. Off-vertical columns pay the base's drift (1 - cos) per step. A sample keeps its own
// delta when nothing in reach goes further. The result is exact for flat floors and walls a stride from a sample
// (the samples bracketing them read the floor plane), and unedited ground is never touched.
static void AddDeltasColumnFiltered(const DensityField& field, const std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>>& near,
                                    float* ext, int E, int st, Vector3 eo) {
    std::unordered_map<TerrainDeltas::Key, const DeltaChunk*> chunks;
    chunks.reserve(near.size() * 2);
    for (const auto& c : near) chunks.emplace(c.first, c.second);
    const float ext3 = (E - 1) * st * kVoxelSize;
    Vector3 up = dynamic_cast<const TerrainField*>(&field) ? Vector3{ 0, 1, 0 }
                                                           : Vector3Normalize(Vector3AddValue(eo, 0.5f * ext3));
    const float ua[3] = { up.x, up.y, up.z };
    const int a = fabsf(ua[0]) >= fabsf(ua[1]) && fabsf(ua[0]) >= fabsf(ua[2]) ? 0 : fabsf(ua[1]) >= fabsf(ua[2]) ? 1 : 2;
    const int b = (a + 1) % 3, c = (a + 2) % 3;
    const int upSign = ua[a] >= 0.0f ? 1 : -1;
    const float penU = DENSITY_GRADIENT * kVoxelSize * (1.0f - fabsf(ua[a])) / kDeltaStep; // stored units per lattice step
    const int o[3] = { (int)lroundf(eo.x / kVoxelSize), (int)lroundf(eo.y / kVoxelSize), (int)lroundf(eo.z / kVoxelSize) };
    const int W = 2 * st;
    auto fdiv = [](int v) { return v >= 0 ? v / kDeltaChunk : -((-v + kDeltaChunk - 1) / kDeltaChunk); };
    struct Entry { int pos; int v; };
    std::vector<Entry> col;
    const int A0 = o[a] - W, A1 = o[a] + (E - 1) * st + W;
    const int strides[3] = { 1, E, E * E }; // ext index = (z * E + y) * E + x
    for (int ic = 0; ic < E; ++ic)
        for (int ib = 0; ib < E; ++ib) {
            int L[3];
            L[b] = o[b] + ib * st; L[c] = o[c] + ic * st;
            const int cb = fdiv(L[b]), cc = fdiv(L[c]);
            const int lb = L[b] - cb * kDeltaChunk, lc = L[c] - cc * kDeltaChunk;
            col.clear();
            for (int pos = A0; pos <= A1;) {
                const int ca = fdiv(pos);
                int C[3]; C[a] = ca; C[b] = cb; C[c] = cc;
                auto it = chunks.find(TerrainDeltas::KeyOf(C[0], C[1], C[2]));
                const int chunkEnd = (ca + 1) * kDeltaChunk; // first lattice index past this chunk along the axis
                if (it == chunks.end() || it->second->empty()) { pos = chunkEnd; continue; }
                const DeltaChunk* dc = it->second;
                for (; pos <= A1 && pos < chunkEnd;) {
                    int l[3]; l[a] = pos - ca * kDeltaChunk; l[b] = lb; l[c] = lc;
                    const int i = (l[2] * kDeltaChunk + l[1]) * kDeltaChunk + l[0];
                    if (!dc->slot[DeltaChunk::BrickOf(i)]) { pos += DeltaChunk::kBrick - (l[a] & (DeltaChunk::kBrick - 1)); continue; }
                    const int v = dc->get(i);
                    if (v) col.push_back({ pos, v });
                    ++pos;
                }
            }
            if (col.empty()) continue;
            size_t first = 0;
            for (int k = 0; k < E; ++k) {
                const int pk = o[a] + k * st;
                while (first < col.size() && col[first].pos < pk - W) ++first;
                int v0 = 0;
                float dneg = 1e30f, dpos = -1e30f;
                int idx3[3]; idx3[a] = k; idx3[b] = ib; idx3[c] = ic;
                float& e = ext[(size_t)idx3[0] * strides[0] + (size_t)idx3[1] * strides[1] + (size_t)idx3[2] * strides[2]];
                for (size_t q = first; q < col.size() && col[q].pos <= pk + W; ++q) if (col[q].pos == pk) { v0 = col[q].v; break; }
                const bool solid = e + v0 * kDeltaStep >= ISO_LEVEL;
                for (size_t q = first; q < col.size() && col[q].pos <= pk + W; ++q) {
                    const int off = (col[q].pos - pk) * upSign, dist = off >= 0 ? off : -off;
                    if (col[q].v < 0) { if (!solid || off >= 0) dneg = std::min(dneg, col[q].v + penU * dist); }
                    else if (solid || off <= 0) dpos = std::max(dpos, col[q].v - penU * dist);
                }
                float d = (float)v0;
                if (dneg < d) d = dneg;
                else if (dpos > d) d = dpos;
                e += d * kDeltaStep;
            }
        }
}

void BuildChunkMeshFromBase(const DensityField& field, const TerrainDeltas& deltas, ChunkMeshData& out, const float* base, float sampledMs, const PlanetField::BlockGrid* grid,
                            std::shared_ptr<const PlanetField::BlockGrid> gridShared) {
    const auto t0 = std::chrono::steady_clock::now() - std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<float, std::milli>(sampledMs));
    const int cells = CHUNK_SIZE + 1;        // +1 cell overlap with the positive neighbours (LOD crack cover)
    const int S = cells + 1;                 // samples per axis
    const float stride = (float)out.key.stride();
    const Vector3 origin = out.key.origin();
    // `base` (E = S + 2 per axis, ChunkBlockOf) keeps the unedited values for materials ("original surface"
    // vs "dug out").
    const int E = S + 2;
    const Vector3 eo = Vector3AddValue(origin, -stride);
    std::vector<float> ext(base, base + (size_t)E * E * E);
    // (#908) The plots near this chunk with building blocks in them: the sky probe below lets their roofs shade the ground.
    std::shared_ptr<const GrassBlockers> roofSnap;
    std::vector<const GrassBlockPlot*> roofPlots;
    static const bool roofSkyOn = !(getenv("RAYLIBGAME_ROOF_SKY") && getenv("RAYLIBGAME_ROOF_SKY")[0] == '0');
    if (roofSkyOn && stride <= 2.0f && (roofSnap = Grass_BlockersOf(&field))) {
        const Vector3 cc = Vector3AddValue(origin, 0.5f * out.key.size());
        const float reach = 0.8660254f * out.key.size() + 12.0f;
        for (const GrassBlockPlot& b : *roofSnap) {
            if (b.solid.empty()) continue;
            if ((DVec3{ cc.x, cc.y, cc.z } - b.origin).length() <= reach + b.radius) roofPlots.push_back(&b);
        }
    }
    // (#1067) The cave shafts open to the sky that reach this chunk (karst sinks, tube skylights): their light pool below.
    static const bool shaftSkyOn = !(getenv("RAYLIBGAME_SHAFT_SKY") && getenv("RAYLIBGAME_SHAFT_SKY")[0] == '0');
    std::vector<PlanetField::CaveShaft> shafts;
    if (shaftSkyOn && stride <= 1.0f)
        if (const auto* pfs = dynamic_cast<const PlanetField*>(&field)) pfs->shaftsNear(origin, Vector3AddValue(origin, out.key.size()), shafts);
    // Player edits: the delta lattice at the samples (lattice points at every stride); coarse nodes filter it along
    // their columns (AddDeltasColumnFiltered, #679), so big excavations show from afar with flat floors. Cost depends
    // on the edited chunks here, never on how many strokes the world has had.
    std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>> near;
    deltas.chunksIn(eo, Vector3AddValue(eo, (E - 1) * stride), near);
    const int stEdit = out.key.latticeStride();
    if (!near.empty() && stEdit > 1 && !gDeltaPointSample) {
        // Coarse nodes (#679): the column filter below, not a point sample. Needs the deltas a column filter reaches.
        near.clear();
        const float reach = 2.0f * stEdit * kVoxelSize;
        deltas.chunksIn(Vector3AddValue(eo, -reach), Vector3AddValue(eo, (E - 1) * stride + reach), near);
        AddDeltasColumnFiltered(field, near, ext.data(), E, stEdit, eo);
    } else if (!near.empty()) {
        // Lattice indices (kVoxelSize apart): a node of level L samples every 2^L-th lattice point.
        const int st = stEdit, ox = (int)lroundf(eo.x / kVoxelSize), oy = (int)lroundf(eo.y / kVoxelSize), oz = (int)lroundf(eo.z / kVoxelSize);
        auto chunkAt = [&](int cx, int cy, int cz) -> const DeltaChunk* {
            const TerrainDeltas::Key key = TerrainDeltas::KeyOf(cx, cy, cz);
            for (const auto& c : near) if (c.first == key) return c.second;
            return nullptr;
        };
        auto fdiv = [](int a) { return a >= 0 ? a / kDeltaChunk : -((-a + kDeltaChunk - 1) / kDeltaChunk); };
        size_t k = 0;
        for (int z = 0; z < E; ++z)
            for (int y = 0; y < E; ++y) {
                const int wy = oy + y * st, wz = oz + z * st, cy = fdiv(wy), cz = fdiv(wz);
                int lastCx = INT32_MIN;
                const DeltaChunk* c = nullptr;
                for (int x = 0; x < E; ++x, ++k) {
                    const int wx = ox + x * st, cx = fdiv(wx);
                    if (cx != lastCx) { lastCx = cx; c = chunkAt(cx, cy, cz); }
                    if (!c) continue;
                    const int lx = wx - cx * kDeltaChunk, ly = wy - cy * kDeltaChunk, lz = wz - cz * kDeltaChunk;
                    ext[k] += c->get((lz * kDeltaChunk + ly) * kDeltaChunk + lx) * kDeltaStep;
                }
            }
    }
    auto at = [&](const float* g, int x, int y, int z) { return g[((size_t)z * E + y) * E + x]; }; // extended coords

    std::vector<float> dens((size_t)S * S * S);
    bool anySolid = false, anyEmpty = false;
    size_t i = 0;
    for (int z = 0; z < S; ++z)
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x, ++i) {
                const float d = at(ext.data(), x + 1, y + 1, z + 1);
                dens[i] = d;
                if (d >= ISO_LEVEL) anySolid = true; else anyEmpty = true;
            }

    out.sampleMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
    // Density gradient at an inner sample (inner coords), by central differences on the extended grid.
    auto gradient = [&](int x, int y, int z) {
        return Vector3{ at(ext.data(), x + 2, y + 1, z + 1) - at(ext.data(), x, y + 1, z + 1),
                        at(ext.data(), x + 1, y + 2, z + 1) - at(ext.data(), x + 1, y, z + 1),
                        at(ext.data(), x + 1, y + 1, z + 2) - at(ext.data(), x + 1, y + 1, z) };
    };
    if (anySolid && anyEmpty) {
        std::vector<Vector3> verts;
        std::vector<int> idx;
        std::vector<int> edgeCache((size_t)S * S * S * 3, -1);
        verts.reserve(4096);
        idx.reserve(12288);
        for (int z = 0; z < cells; ++z)
            for (int y = 0; y < cells; ++y)
                for (int x = 0; x < cells; ++x)
                    marchCube(x, y, z, dens, S, verts, idx, edgeCache, ISO_LEVEL);

        if (verts.size() > 65535) {
            out.overflow = true; // raylib meshes use 16-bit indices
        } else if (!verts.empty()) {
            out.positions.resize(verts.size() * 3);
            out.normals.resize(verts.size() * 3);
            out.colors.assign(verts.size() * 4, 0);
            out.texcoords.assign(verts.size() * 2, 0.0f);
            out.texcoords2.assign(verts.size() * 2, 0.0f);
            out.tangents.assign(verts.size() * 4, 0.0f);
            out.occlusion.assign(verts.size(), 0.0f);
            out.matIds.assign(verts.size(), 0);
            const PlanetField* planetLife = dynamic_cast<const PlanetField*>(&field);
            if (planetLife && !planetLife->foliage) planetLife = nullptr;
            const PlanetField* gridField = grid ? dynamic_cast<const PlanetField*>(&field) : nullptr;
            for (size_t v = 0; v < verts.size(); ++v) {
                const Vector3 lp = Vector3Scale(verts[v], stride);     // node-relative, in voxels
                const Vector3 p = Vector3Add(origin, lp);
                // A marching-cubes vertex lies on a grid edge: interpolate the gradient (and the unedited
                // density) between the edge's two samples, at the vertex's position along it.
                const Vector3 g = verts[v];
                int ax = 0;
                float t = g.x - floorf(g.x);
                if (g.y - floorf(g.y) > t) { ax = 1; t = g.y - floorf(g.y); }
                if (g.z - floorf(g.z) > t) { ax = 2; t = g.z - floorf(g.z); }
                const int x0 = std::min((int)floorf(g.x), S - 1), y0 = std::min((int)floorf(g.y), S - 1), z0 = std::min((int)floorf(g.z), S - 1);
                const int x1 = std::min(x0 + (ax == 0 ? 1 : 0), S - 1), y1 = std::min(y0 + (ax == 1 ? 1 : 0), S - 1), z1 = std::min(z0 + (ax == 2 ? 1 : 0), S - 1);
                const Vector3 grad = Vector3Lerp(gradient(x0, y0, z0), gradient(x1, y1, z1), t);
                const float gl = Vector3Length(grad);
                const Vector3 n = gl > 1e-6f ? Vector3Scale(grad, -1.0f / gl) : field.fallbackUp(p);
                const float b = at(base, x0 + 1, y0 + 1, z0 + 1) * (1.0f - t) + at(base, x1 + 1, y1 + 1, z1 + 1) * t;
                float sky = 1.0f;
                // The surface under the vertex from the block's height grid (the one its density was sampled from), not
                // a fresh height-map evaluation per vertex (#148).
                float gElev = 0.0f;
                PlanetField::Column3D gCol;
                const bool hinted = gridField && gridField->gridSurface(*grid, p, gElev, gCol);
                // Dug ground (#679): the materials allow a coarse mesh to stray metres from the height map (slack ~ 2 x
                // the footprint), so a floor dug 3 m down still read as undug turf on coarse nodes: grass under the water
                // of a dug basin from afar. Where the edge's samples were dug by more than half a metre, judge the vertex
                // at a 1 m footprint (no slack): the edit, not the LOD, put the surface there.
                const bool dugHere = !near.empty() && !gDeltaPointSample && std::min(at(ext.data(), x0 + 1, y0 + 1, z0 + 1) - at(base, x0 + 1, y0 + 1, z0 + 1),
                                                              at(ext.data(), x1 + 1, y1 + 1, z1 + 1) - at(base, x1 + 1, y1 + 1, z1 + 1)) < -0.5f * DENSITY_GRADIENT;
                const float matFp = dugHere ? std::min(stride, 1.0f) : stride;
                VoxelMaterial m = hinted ? gridField->materialSkyHinted(p, b, matFp, n, sky, &gElev, &gCol) : field.materialSky(p, b, matFp, n, sky);
                if (!near.empty()) { const int placed = deltas.materialAt(p); if (placed >= 0) m = (VoxelMaterial)placed; } // placed ground keeps its material
                float roofOpen = 1.0f; // (#908) the built roofs' share of the sky here (ambient only: folded into the occlusion below)
                if (gSkyProbes && stride <= 2.0f) {
                    const PlanetField::Column3D* pc = hinted ? &gCol : nullptr; const float* pe = hinted ? &gElev : nullptr;
                    // (#898) Under the surface (the depth term < 1) and less than 40 m down: light comes in where the sky can be
                    // seen (a cave mouth, a shaft, a pit): at least 0.6 x (1 - depth / 40 m) of the open share (0.42 at 12 m: walls stay
                    // under the open ground and fall off deeper in), the rays checked past the block. The
                    // depth term alone cut to black ~10 m down, under an opening or not; a closed tunnel stays black.
                    const float depthHere = hinted ? gridField->radius + gElev - Vector3Length(p) : 1e9f;
                    static const bool mouthOn = !(getenv("RAYLIBGAME_CAVE_MOUTH_SKY") && getenv("RAYLIBGAME_CAVE_MOUTH_SKY")[0] == '0'); // A/B (#898): 0 = the depth term alone
                    if (mouthOn && sky < 0.999f && depthHere < 40.0f && stride <= 1.0f) sky = std::max(sky, 0.6f * (1.0f - depthHere / 40.0f)) * OverheadOpenness(field, deltas, ext.data(), E, eo, stride, p, n, pc, pe, true);
                    else if (sky > 0.02f) sky *= OverheadOpenness(field, deltas, ext.data(), E, eo, stride, p, n, pc, pe);
                    // (#1067, user on Plan 10: "a vertical pit or shaft open to the sky must be lit properly from outside all the way
                    // down") The probes above see only 12 m and a shaft's walls block their tilted rays, so a shaft read as buried.
                    // 1. A vertex (not a roof facing down) whose column is open all the way to the surface sees the sky straight up:
                    //    at least 0.85, at any depth (the column walked in 2 m steps; a closed cave stops at its first roof sample).
                    // 2. Under a shaft's foot (karst sink, tube skylight): a pool of up to 0.5, fading from the shaft's radius to 2.5
                    //    radii around its axis and from 4 to 12 m below its foot. RAYLIBGAME_SHAFT_SKY=0: off.
                    if (shaftSkyOn && stride <= 1.0f && hinted && depthHere > 2.0f && depthHere < 150.0f && sky < 0.85f) {
                        const Vector3 upv = Vector3Normalize(p);
                        if (Vector3DotProduct(n, upv) > -0.3f) {
                            const Vector3 p0 = Vector3Add(p, Vector3Scale(n, 0.5f * stride));
                            const float fpc = std::max(stride, 1.0f);
                            bool open = true;
                            for (float t = 1.0f; t < depthHere + 3.0f && open; t += 2.0f)
                                open = field.density(Vector3Add(p0, Vector3Scale(upv, t)), &deltas, fpc) < ISO_LEVEL;
                            if (open) sky = std::max(sky, 0.85f);
                        }
                        for (const PlanetField::CaveShaft& s : shafts) {
                            const Vector3 v = Vector3Subtract(p, s.top);
                            const float along = -Vector3DotProduct(v, s.up);      // metres below the shaft's top
                            if (along < 0.5f * s.depth || along > s.depth + 12.0f) continue;
                            const float horiz = Vector3Length(Vector3Add(v, Vector3Scale(s.up, along)));
                            const float a = std::clamp((horiz - s.radius) / (1.5f * s.radius), 0.0f, 1.0f);
                            const float b = std::clamp((along - s.depth - 4.0f) / 8.0f, 0.0f, 1.0f);
                            sky = std::max(sky, 0.5f * (1.0f - a * a * (3.0f - 2.0f * a)) * (1.0f - b * b * (3.0f - 2.0f * b)));
                        }
                    }
                    // (#908) Built roofs: the same rays (up, and four at 45 degrees; the up one weighs 3) through the solid cells of the
                    // plots near this chunk (the #830 snapshot): the ground under a shelter, a covered yard or a house without a floor
                    // no longer takes the full sky. RAYLIBGAME_ROOF_SKY=0: off.
                    if (sky > 0.02f && !roofPlots.empty()) {
                        constexpr float kCellSize = 0.5f; // plots.h (not included here)
                        float partsOpen = 1.0f;
                        for (const GrassBlockPlot* b : roofPlots) {
                            const DVec3 d = DVec3{ p.x, p.y, p.z } - b->origin;
                            const Vector3 l = { (float)d.x, (float)d.y, (float)d.z };
                            const Vector3 c = { Vector3DotProduct(l, b->east) / kCellSize, Vector3DotProduct(l, b->up) / kCellSize, Vector3DotProduct(l, b->north) / kCellSize };
                            if (Vector3Length(c) * kCellSize > b->radius + 12.0f) continue;
                            static const Vector3 dirs[5] = { { 0, 1, 0 }, { 0.7071f, 0.7071f, 0 }, { -0.7071f, 0.7071f, 0 }, { 0, 0.7071f, 0.7071f }, { 0, 0.7071f, -0.7071f } };
                            float open = 0.0f;
                            for (int r = 0; r < 5; ++r) {
                                bool blocked = false;
                                for (float s = 0.6f; s <= 24.0f && !blocked; s += 0.5f) { // cells: 12 m
                                    const Vector3 q = Vector3Add(c, Vector3Scale(dirs[r], s));
                                    blocked = b->solid.count(Int3{ (int)floorf(q.x), (int)floorf(q.y), (int)floorf(q.z) }) > 0;
                                }
                                if (!blocked) open += r == 0 ? 3.0f : 1.0f;
                            }
                            partsOpen = std::min(partsOpen, std::clamp(open / 7.0f * 1.25f, 0.0f, 1.0f));
                        }
                        roofOpen = partsOpen; // not the sky term: that scales the sun too, and a roof's shadow comes from the shadow maps (parts cast)
                    }
                }
                out.tangents[v * 4 + 3] = sky;
                // Ambient occlusion (#132): how much of the hemisphere over the vertex is solid within ~1.5 m, from the
                // block's own samples (edits included): 5 directions (the normal, and four 45 degrees off it) at 0.6 and
                // 1.4 m. Dug holes, crevices and the foot of walls darken; open ground stays 0. Near LODs only (the
                // samples are too coarse beyond 2 m spacing).
                if (gTerrainAO && stride <= 2.0f) {
                    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(n, fabsf(n.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                    const Vector3 t2 = Vector3CrossProduct(n, t1);
                    const Vector3 dirs[5] = { n, Vector3Normalize(Vector3Add(n, t1)), Vector3Normalize(Vector3Subtract(n, t1)),
                                              Vector3Normalize(Vector3Add(n, t2)), Vector3Normalize(Vector3Subtract(n, t2)) };
                    auto solidAt = [&](Vector3 e) { // trilinear on the extended block (coords in samples)
                        e.x = std::clamp(e.x, 0.0f, E - 1.001f); e.y = std::clamp(e.y, 0.0f, E - 1.001f); e.z = std::clamp(e.z, 0.0f, E - 1.001f);
                        const int x0 = (int)e.x, y0 = (int)e.y, z0 = (int)e.z;
                        const float fx = e.x - x0, fy = e.y - y0, fz = e.z - z0;
                        auto A = [&](int x, int y, int z) { return at(ext.data(), x, y, z); };
                        const float c00 = A(x0, y0, z0) * (1 - fx) + A(x0 + 1, y0, z0) * fx, c10 = A(x0, y0 + 1, z0) * (1 - fx) + A(x0 + 1, y0 + 1, z0) * fx;
                        const float c01 = A(x0, y0, z0 + 1) * (1 - fx) + A(x0 + 1, y0, z0 + 1) * fx, c11 = A(x0, y0 + 1, z0 + 1) * (1 - fx) + A(x0 + 1, y0 + 1, z0 + 1) * fx;
                        return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz >= ISO_LEVEL;
                    };
                    const Vector3 e0 = Vector3AddValue(g, 1.0f); // the vertex in extended-block samples
                    float occ = 0.0f;
                    for (int r = 0; r < 5; ++r)
                        for (float m : { 0.6f, 1.4f })
                            if (solidAt(Vector3Add(e0, Vector3Scale(dirs[r], m / stride)))) occ += (r == 0 ? 1.5f : 1.0f) * (m < 1.0f ? 1.3f : 0.7f);
                    out.occlusion[v] = std::clamp(occ / 7.0f, 0.0f, 1.0f);
                }
                if (roofOpen < 1.0f) out.occlusion[v] = std::max(out.occlusion[v], 1.0f - roofOpen); // (#908) also with RAYLIBGAME_AO=0
                if (m == VoxelMaterial::Ore) {
                    // Which ore (the shader colours each kind): 1 + kind, over 8.
                    if (const auto* pfo = dynamic_cast<const PlanetField*>(&field)) {
                        const float dl = Vector3Length(p);
                        const int k = dl > 1.0f ? pfo->oreKind(p, pfo->radius + (hinted ? gElev : pfo->elevation(Vector3Scale(p, 1.0f / dl), (float)stride)) - dl) : -1;
                        out.tangents[v * 4 + 2] = (k + 1) / 8.0f;
                    }
                }
                if (planetLife) { // under overhangs too: there was no climate where the sky was hidden, so no tint (#254)
                    // The climate (Plan 5 B) for the terrain shader's biome tints, packed in y: 2 + the
                    // temperature step (0-62 over -30..40 C) + moisture in the fraction (0 = no climate). The shader reads it
                    // flat (per triangle): interpolated, the packing swept through whole steps and striped the tint (#254).
                    const float d = Vector3Length(p);
                    const Vector3 up = Vector3Scale(p, 1.0f / d);
                    const PlanetField::Climate cl = planetLife->climateAt(up, d - planetLife->radius);
                    out.tangents[v * 4 + 1] = 2.0f + floorf(std::clamp((cl.temp + 30.0f) / 70.0f, 0.0f, 1.0f) * 62.0f) + std::min(cl.moist, 0.999f); // from 2 (#254)
                    // Distant forests: canopy cover (x), painted in where the real trees fade out.
                    if (sky > 0.9f && (m == VoxelMaterial::Grass || m == VoxelMaterial::Regolith)) {
                        const float steep = 1.0f - Vector3DotProduct(n, up);
                        out.tangents[v * 4 + 0] = TreeCanopy(up, planetLife->radius, d - planetLife->radius, planetLife->foliageSeed, planetLife) * std::clamp(1.0f - (steep - 0.1f) / 0.15f, 0.0f, 1.0f);
                    }
                }
                out.positions[v * 3 + 0] = lp.x;
                out.positions[v * 3 + 1] = lp.y;
                out.positions[v * 3 + 2] = lp.z;
                out.normals[v * 3 + 0] = n.x;
                out.normals[v * 3 + 1] = n.y;
                out.normals[v * 3 + 2] = n.z;
                out.matIds[v] = (unsigned char)m; // #480: the arena's per-chunk palette
                // The legacy one-hot weight (per-chunk draws, grass scattering): materials past the 8 channels show as
                // their nearest (clay as dirt).
                const int mi = m == VoxelMaterial::Pebbles ? (int)VoxelMaterial::Gravel                                              // river stones as gravel (#595)
                             : (m == VoxelMaterial::Clay || m == VoxelMaterial::Silt || m == VoxelMaterial::Mud || m == VoxelMaterial::Farmland) ? (int)VoxelMaterial::Regolith
                             : (m >= VoxelMaterial::Granite && m <= VoxelMaterial::Marble) ? (int)VoxelMaterial::Rock : (int)m; // rock strata (#587): as rock
                if (mi < 4) out.colors[v * 4 + mi] = 255;
                else if (mi < 6) out.texcoords[v * 2 + (mi - 4)] = 1.0f;
                else out.texcoords2[v * 2 + (mi - 6)] = 1.0f;
            }
            out.indices.assign(idx.begin(), idx.end());
        }
    }

    // Water (#360): the generated water and the simulation's published cells, one surface. Runs even when the node
    // holds no terrain surface (open ocean). Chunks next to simulated water keep what a water-only rebuild needs.
    {
        const WaterField::SnapshotPtr snap = WaterField::Current();
        WaterField::ChunkSampler sim(snap.get(), out.planetoidId, origin, Vector3AddValue(origin, cells * stride));
        out.waterSeq = snap->seq;
        ChunkWaterGen gen = GeneratedWaterOf(field, origin, cells, stride);
        MeshChunkWater(field, out, dens.data(), grid, gen, sim);
        if (sim.any() && gWaterCacheBytes.load() < kChunkWaterCacheCap) {
            auto wc = std::make_shared<ChunkWaterCache>();
            wc->dens = dens;
            wc->gen = std::move(gen);
            const bool copied = !gridShared && grid;
            // The grid the build already shares (GPU path), else a copy (the CPU path's grid lives on its stack).
            if (gridShared) wc->grid = std::move(gridShared);
            else if (grid) wc->grid = std::make_shared<const PlanetField::BlockGrid>(*grid);
            wc->counted = (long long)wc->dens.size() * (long long)sizeof(float) + (copied ? (long long)sizeof(PlanetField::BlockGrid) + (long long)grid->h.size() * 4 + (long long)grid->cols.size() * (long long)sizeof(PlanetField::Column3D) : 0);
            gWaterCacheBytes += wc->counted;
            out.waterCache = std::move(wc);
        }
    }
    ScatterGrass(field, out); // foliage: grass tufts on the finest level (in the build time)
    if (TerrainArena_PackOnWorkers() && !out.positions.empty()) TerrainArena_Pack(out, out.arena); // #302: off the main thread
    out.buildMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
}

// ---------------------------------------------------------------------------------------------
// Ore: deposits and scattered pockets
// ---------------------------------------------------------------------------------------------

namespace {
    constexpr float kDepositCell = 96.0f;
    uint32_t OreHash(int x, int y, int z, uint32_t seed) {
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
        h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
        return h;
    }
    float OreUnit(uint32_t& h) { h = h * 1664525u + 1013904223u; return (h >> 8) * (1.0f / 16777216.0f); }
    struct Deposit { bool exists; Vector3 centre; float radius; int kind; };
    Deposit DepositOf(int cx, int cy, int cz, uint32_t seed) {
        uint32_t h = OreHash(cx, cy, cz, seed);
        Deposit d{ false, {}, 0, -1 };
        if (OreUnit(h) > 0.35f) return d;
        d.exists = true;
        d.centre = { (cx + 0.25f + 0.5f * OreUnit(h)) * kDepositCell, (cy + 0.25f + 0.5f * OreUnit(h)) * kDepositCell, (cz + 0.25f + 0.5f * OreUnit(h)) * kDepositCell };
        d.radius = 3.5f + 5.5f * OreUnit(h); // 7-18 m across
        const float r = OreUnit(h);          // coal 22%, iron 22%, copper 20%, tin 12%, bauxite 14%, ilmenite 10%
        d.kind = r < 0.22f ? 3 : r < 0.44f ? 0 : r < 0.64f ? 1 : r < 0.76f ? 2 : r < 0.90f ? 4 : 5;
        return d;
    }
    float MinDepth(int kind) { return kind == 4 ? 20.0f : kind == 5 ? 60.0f : 4.0f; }
}

const char* PlanetField::OreKindName(int kind) {
    static const char* names[kOreKindCount] = { "iron_ore", "copper_ore", "tin_ore", "coal", "bauxite_ore", "ilmenite_ore" };
    return kind >= 0 && kind < kOreKindCount ? names[kind] : nullptr;
}

int PlanetField::oreKind(Vector3 p, float depth) const {
    if (moon || depth < 3.0f) return -1;
    // A deposit (its blob stays inside its own cell: centre at least 24 m from the edges, radius <= 9 m
    // plus the lumpiness).
    const int cx = (int)floorf(p.x / kDepositCell), cy = (int)floorf(p.y / kDepositCell), cz = (int)floorf(p.z / kDepositCell);
    const Deposit d = DepositOf(cx, cy, cz, oreSeed);
    if (d.exists && depth >= MinDepth(d.kind)) {
        const float k = Vector3Distance(p, d.centre) / d.radius + 0.35f * veins.fractal(2, p.x * 3.0f, p.y * 3.0f, p.z * 3.0f);
        if (k < 1.0f) return d.kind;
    }
    // Scattered pockets of the common ores and coal, more of them near the surface.
    const float pocket = veins.fractal(1, p.x * 4.0f, p.y * 4.0f, p.z * 4.0f);
    const float threshold = 0.72f + 0.12f * std::clamp((depth - 20.0f) / 200.0f, 0.0f, 1.0f) - 0.08f * (gen.resources - 1.0f); // #504: richer, more pockets
    if (pocket > threshold) {
        uint32_t h = OreHash((int)floorf(p.x / 8.0f), (int)floorf(p.y / 8.0f), (int)floorf(p.z / 8.0f), oreSeed ^ 0x5eedu);
        const float r = OreUnit(h);
        return r < 0.35f ? 0 : r < 0.60f ? 1 : r < 0.75f ? 2 : 3; // iron, copper, tin, and coal (the seams near the surface)
    }
    return -1;
}

bool PlanetField::depositAt(Vector3 p, Vector3& centre, float& radius, int& kind) const {
    if (moon) return false;
    const Deposit d = DepositOf((int)floorf(p.x / kDepositCell), (int)floorf(p.y / kDepositCell), (int)floorf(p.z / kDepositCell), oreSeed);
    if (!d.exists) return false;
    centre = d.centre; radius = d.radius; kind = d.kind;
    return true;
}

int PlanetField::depositNear(Vector3 p, float reach) const {
    const int x0 = (int)floorf((p.x - reach) / kDepositCell), x1 = (int)floorf((p.x + reach) / kDepositCell);
    const int y0 = (int)floorf((p.y - reach) / kDepositCell), y1 = (int)floorf((p.y + reach) / kDepositCell);
    const int z0 = (int)floorf((p.z - reach) / kDepositCell), z1 = (int)floorf((p.z + reach) / kDepositCell);
    int best = -1;
    float bestD = 1e30f;
    for (int z = z0; z <= z1; ++z)
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) {
                const Deposit d = DepositOf(x, y, z, oreSeed);
                if (!d.exists) continue;
                const float dist = Vector3Distance(p, d.centre) - d.radius;
                if (dist < reach && dist < bestD) { bestD = dist; best = d.kind; }
            }
    return best;
}

StrataColumn PlanetField::strataColumn(Vector3 dir) const {
    // The continent value at a fixed coarse footprint: the bake's cells and the mesher's vertices see the same province.
    const float sx = dir.x * radius, sy = dir.y * radius, sz = dir.z * radius;
    const float c = continents.fractalLimited(4, 2000.0f, sx, sy, sz) - 0.15f * gen.ocean;
    StrataColumn col = Strata_Column({ sx, sy, sz }, c, fieldSeed);
    col.hardShare = (uint8_t)std::clamp(col.hardShare * gen.cliffs, 0.0f, 255.0f); // New World "Cliffs and canyons" (#589)
    if (gen.karst && col.province == 0 && Strata_KarstWeight({ sx, sy, sz }, fieldSeed) > 0.3f) col.province = 3; // karst limestone (#592)
    if (canyonFault(sx, sy, sz) > 0.3f) { col.province = 0; col.hardShare = 0; } // canyon lines are shattered, soft rock (#544)
    if (gen.volcanoes) { float vw = 0.0f; if (Strata_Volcano({ sx, sy, sz }, radius, fieldSeed, &vw) > 0.0f) col.province = 2; } // lava over ash (#593)
    return col;
}

float PlanetField::canyonLine(float sx, float sy, float sz) const { return fabsf(ravines.fractal(2, sx, sy, sz)); }

float PlanetField::canyonFault(float sx, float sy, float sz) const {
    // The ravine lines (zero crossings of `ravines` in the `ravineMask` regions) as a 0..1 notch ~200 m wide (#544).
    const float region = Smoothstep(0.25f, 0.4f, ravineMask.fractal(2, sx, sy, sz));
    if (region <= 0.0f) return 0.0f;
    const float r = canyonLine(sx, sy, sz);
    if (r >= 0.02f) return 0.0f;
    const float t = 1.0f - r / 0.02f;
    return t * t * (3.0f - 2.0f * t) * region;
}

// The strata column at a body-local point: it varies over kilometres, so it is cached per thread on a 16 m grid (the
// mesher asks for every rock vertex, land() for every fine sample on steep rock).
static const StrataColumn& CachedColumn(const PlanetField& pf, Vector3 p, float d) {
    struct Entry { const PlanetField* f = nullptr; int32_t x = 0, y = 0, z = 0; StrataColumn c; };
    thread_local Entry cache[1024];
    const int32_t qx = (int32_t)floorf(p.x / 16.0f), qy = (int32_t)floorf(p.y / 16.0f), qz = (int32_t)floorf(p.z / 16.0f);
    Entry& e = cache[((uint32_t)qx * 73856093u ^ (uint32_t)qy * 19349663u ^ (uint32_t)qz * 83492791u) & 1023u];
    if (e.f != &pf || e.x != qx || e.y != qy || e.z != qz) { e.f = &pf; e.x = qx; e.y = qy; e.z = qz; e.c = pf.strataColumn(Vector3Scale(p, 1.0f / d)); }
    return e.c;
}

RockType PlanetField::rockTypeAt(Vector3 p) const {
    if (moon) return RockType::Basalt; // the maria: dark basalt plains
    const float d = Vector3Length(p);
    if (d < 1.0f) return RockType::Soft;
    const PlanetErosion* eroded = erosion();
    const float scale = eroded && eroded->heightScale > 0.0f ? eroded->heightScale : 1.0f;
    return Strata_RockAt(CachedColumn(*this, p, d), (d - radius) / scale, fieldSeed);
}

float PlanetField::cliffFace(Vector3 s, float h, Vector3 grad, float fp) const {
    // Plan 9 step 5 (#588): the bake's 50 m grid spreads a resistant layer's outcrop over a cell or two, never sheer.
    // Where the bake left the ground steep and the layer here is resistant, its height range is redistributed into one
    // face in the layer's middle half, a ledge above it and a talus apron below; the column keeps its mean height and
    // the soft layers between stay as the bake made them. One face per resistant layer (they are irregular: 40-150 m,
    // 35-60% of the layers, warped and tilted), not the regular treads of the old terraces (#543). Granite cores keep
    // the bake's shape. Fine footprints only (faded out by 32 m).
    if (!gCliffFaces || fp >= 32.0f || h < 2.0f) return 0.0f;
    const float steep = Smoothstep(0.6f, 1.2f, Vector3Length(grad));
    if (steep <= 0.0f) return 0.0f;
    const float d = Vector3Length(s);
    if (d < 1.0f) return 0.0f;
    const StrataColumn& c = CachedColumn(*this, s, d);
    if (c.province == 1) return 0.0f;
    const PlanetErosion* eroded = erosion();
    const float scale = eroded && eroded->heightScale > 0.0f ? eroded->heightScale : 1.0f;
    const float hs = h / scale;
    if (Strata_RockAt(c, hs, fieldSeed) == RockType::Soft) return 0.0f;
    const float thick = (c.province == 2 ? 0.6f : 1.0f) * c.thick;
    const float t = (hs + c.offset) / thick, f = t - floorf(t);
    const float face = Smoothstep(0.25f, 0.75f, f);
    return (face - f) * thick * scale * steep * std::clamp(2.0f - fp / 16.0f, 0.0f, 1.0f);
}

// The rock item at a point: its stratum's (#587: one data source with the bake and the surface material).
const char* PlanetField::stoneKind(Vector3 p) const { return Strata_Item(rockTypeAt(p)); }

const char* GemForStone(const char* stone) {
    if (!stone) return nullptr;
    if (strcmp(stone, "granite") == 0) return "amethyst";
    if (strcmp(stone, "marble") == 0) return "ruby";
    if (strcmp(stone, "basalt") == 0) return "sapphire";
    if (strcmp(stone, "limestone") == 0) return "emerald";
    return nullptr;
}


std::string PlanetField::debugColumn(Vector3 dir) const {
    char buf[512];
    const Land l = land(dir, 0.0f);
    int n = snprintf(buf, sizeof(buf), "land %.1f (fp 256: %.1f) base %.1f final %.1f", l.h, land(dir, 256.0f).h, elevationBase(dir, 0.0f), elevation(dir, 0.0f));
    if (const RiverNetwork* net = rivers()) {
        const RiverNetwork::Hit rh = net->nearest(dir, net->reach());
        if (rh.found) n += snprintf(buf + n, sizeof(buf) - n, "; river %.0f m away, level %.1f, half %.1f", rh.dist, rh.level, rh.halfWidth);
        float lv, lo; int id;
        if (net->lakeAt(dir, lv, lo, id)) n += snprintf(buf + n, sizeof(buf) - n, "; lake %d level %.1f outside %.0f", id, lv, lo);
    }
    return buf;
}

std::vector<Vector3> PlanetField::bakePondsByDistance(Vector3 dir) const {
    std::vector<Vector3> out;
    const PlanetErosion* ero = moon ? nullptr : erosion();
    if (!ero) return out;
    for (const PlanetErosion::Pond& p : ero->ponds()) out.push_back(p.dir);
    std::sort(out.begin(), out.end(), [&](const Vector3& a, const Vector3& b) { return Vector3DotProduct(dir, a) > Vector3DotProduct(dir, b); });
    return out;
}
