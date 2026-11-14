#include <unordered_map>
#include <cstring>
#include "voxel.h"
#include "erosion.h"
#include <thread>
#include "rivers.h"
#include "foliage.h"
#include "trees.h"
#include "cubemarch.h"

#include <cstring>
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>

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

static SimplexNoise SeededNoise(std::mt19937& rng, float frequency) {
    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    return SimplexNoise(frequency, 1.0f, 2.0f, 0.5f, unit(rng) * 10000.0f, unit(rng) * 10000.0f, unit(rng) * 10000.0f);
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

// River valleys: 70 m deep, with a flat floodplain floor (valley-noise units < kValleyFloor, room for
// the widest river) and sides rising over the next 0.04 units. (A V-shaped floor rose ~27 m across a
// river's width and tilted its water.)
static constexpr float kValleyFloor = 0.035f, kValleyEdge = 0.075f;
static float ValleyCarve(float v) {
    if (v >= kValleyEdge) return 0.0f;
    const float t = 1.0f - std::max(v - kValleyFloor, 0.0f) / (kValleyEdge - kValleyFloor);
    return t * t * 70.0f;
}


struct PlanetField::ErosionState { std::atomic<bool> active{ false }; PlanetErosion bake; };

const PlanetErosion* PlanetField::erosion() const {
    return erosionState && erosionState->active.load(std::memory_order_acquire) ? &erosionState->bake : nullptr;
}

void PlanetField::bakeErosion(const std::string& cachePath) {
    if (moon || !erosionState) return;
    // The whole planet (Plan 7 step 1, #355), on the raw land; then switched on. Loaded from cachePath when a bake of
    // this seed and radius is there (a world load reads it instead of baking again). Any river network built before
    // came from the raw land: start a fresh one, built on first use from the eroded land.
    PlanetErosion& b = erosionState->bake;
    const bool loaded = !cachePath.empty() && b.load(cachePath, radius, fieldSeed);
    if (!loaded) {
        ErosionParams p;
        auto envf = [](const char* n, float v) { const char* e = getenv(n); return e ? (float)atof(e) : v; }; // tuning
        p.k = envf("RAYLIBGAME_EROSION_K", p.k); p.steps = (int)envf("RAYLIBGAME_EROSION_STEPS", (float)p.steps);
        p.cell = envf("RAYLIBGAME_EROSION_CELL", p.cell); p.uplift = envf("RAYLIBGAME_EROSION_UPLIFT", p.uplift);
        p.start = envf("RAYLIBGAME_EROSION_START", p.start);
        p.refine = (int)envf("RAYLIBGAME_EROSION_REFINE", (float)p.refine); p.refineSteps = (int)envf("RAYLIBGAME_EROSION_REFINE_STEPS", (float)p.refineSteps);
        const int threads = std::max(1u, std::thread::hardware_concurrency());
        b.build([this](Vector3 d, float fp) { return land(d, fp).h; }, radius, fieldSeed, p, threads);
        if (!cachePath.empty()) b.save(cachePath);
        printf("erosion: baked the planet, %d cells at %.0f m + %d at %.0f m, in %.1f s (coarse %.1f s, faces %.1f s): cut up to %.0f m, raised up to %.0f m, hash %016llx\n",
               b.cells(), p.cell, b.fineCells(), b.fineCell(), b.buildSeconds, b.coarseSeconds, b.refineSeconds, b.maxCut, b.maxRaise, (unsigned long long)b.hash());
    } else {
        printf("erosion: loaded %s (%d cells, hash %016llx)\n", cachePath.c_str(), b.cells(), (unsigned long long)b.hash());
    }
    erosionState->active.store(true, std::memory_order_release);
    riverState = std::make_shared<RiverState>();
}
float PlanetField::profileLandHeight(Vector3 dir, float fp) const { return land(dir, fp).h; }
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
    const float c = continents.fractalLimited(4, fp, sx, sy, sz); // < 0: ocean, > 0: land
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
        if (gGullies && fp < 450.0f) h += ErosionGullies({ sx, sy, sz }, dir, s.grad, fp, fieldSeed, Smoothstep(2.0f, 15.0f, s.h));
        l.ridge = std::max(0.0f, h - 150.0f);
    } else {
        h = base;
        if (mountainMask > 0.0f) {
            const float r = 1.0f - fabsf(ridges.fractalLimited(5, fp, sx, sy, sz)); // ridged: sharp crests
            l.ridge = r * r * r * 1400.0f * mountainMask; // the height the ranges add here
            h += l.ridge;
        }
        h += hills.fractalLimited(3, fp, sx, sy, sz) * (25.0f + 60.0f * mountainMask);
    }
    // (River valleys are carved by the drainage network's rivers in elevation(); the old valley-noise
    // valleys and streams are gone, so l.valley and l.stream stay at 1: none.)
    const float valleyWeight = landMask * std::clamp(2.0f - fp / 64.0f, 0.0f, 1.0f);
    l.valley = 1.0f;
    // Ravines: rare narrow canyons (~60 m wide, up to ~90 m deep, steep sqrt-profile walls) in some
    // land regions, faded out where samples are too far apart to resolve them.
    const float ravineWeight = landMask * std::clamp(2.0f - fp / 24.0f, 0.0f, 1.0f);
    if (ravineWeight > 0.0f) {
        const float region = Smoothstep(0.25f, 0.4f, ravineMask.fractal(2, sx, sy, sz));
        if (region > 0.0f) {
            const float r = fabsf(ravines.fractal(2, sx, sy, sz));
            if (r < 0.006f) h -= sqrtf(1.0f - r / 0.006f) * 90.0f * region * ravineWeight;
        }
    }
    l.detail = detail.fractalLimited(2, fp, sx, sy, sz) * 2.5f;
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
    // Rivers on land valleys, not in the mountains, above the sea (they meet it at the coast).
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

float PlanetField::elevationBase(Vector3 dir, float fp) const {
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
    const Land l = land(dir, fp);
    if (moon) return l.h;
    float h = l.h + l.detail;
    h += biomeShape(dir, h, fp);
    // Rivers of the drainage network: a valley (a floodplain 2.5 channel half-widths across at the
    // water, then sides rising 0.3 m per m) cut into higher ground, low banks where the ground dips
    // below the water, and the channel bed. The valley stays resolved further out than the channel.
    if (const RiverNetwork* net = rivers()) {
        // Lakes (rivers.h): the level is held under the land around the basin, so the shore is the
        // terrain's own. Between the samples that checked it the ground may still dip lower; a low
        // bank in the outer part of the water region (at most 4 m of fill, 0.3 m above the water)
        // keeps the water in there, and fades out beyond the region.
        float lv = 0.0f, lout = 0.0f;
        const bool nearLake = net->lakeAt(dir, lv, lout);
        const bool inLake = nearLake && lout < RiverNetwork::kLakeMargin;
        if (nearLake && lout > 0.0f && h < lv + 0.3f) {
            const float M = RiverNetwork::kLakeMargin;
            const float k = Smoothstep(0.3f * M, 0.8f * M, lout) * (1.0f - Smoothstep(1.2f * M, 2.0f * M, lout));
            h += std::min(lv + 0.3f - h, 4.0f) * k;
        }
        const RiverNetwork::Hit rh = net->nearest(dir, RiverNetwork::kMaxReach);
        // In a lake, only a river arriving at the lake's level shapes the ground (its valley and bed
        // lead into the lake; no banks, which stood up out of the water as sand ridges); the outflow
        // and any other river leave the lake alone.
        if (rh.found && (!inLake || fabsf(rh.level - lv) < 1.0f)) {
            const float d = rh.dist, half = rh.halfWidth;
            const float flat = 2.5f * half, reach = std::min(flat + 200.0f, RiverNetwork::kMaxReach);
            const float allowed = rh.level + 0.5f + std::max(0.0f, d - flat) * 0.3f;
            const float vw = (1.0f - Smoothstep(reach * 0.7f, reach, d)) * std::clamp(2.0f - fp / 512.0f, 0.0f, 1.0f);
            if (h > allowed) h += (allowed - h) * vw;
            const float cw = std::clamp(2.0f - fp / (4.0f * half), 0.0f, 1.0f) * Smoothstep(-2.0f, 0.5f, l.h); // not out at sea
            const float bank = rh.level + 0.4f;
            if (!inLake && h < bank && d < half + 12.0f) h += (bank - h) * (1.0f - Smoothstep(half, half + 12.0f, d)) * cw;
            if (d < half) {
                const float t = 1.0f - (d / half) * (d / half);
                const float bed = rh.level + 0.3f - (rh.depth + 0.3f) * t;
                h += (bed - h) * cw;
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
            h += (target - h) * pw;
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
    if (badlands > 0.0f) {
        // Mesas: the height in 14 m steps, each step a flat top that rises sharply at its edge.
        const float step = 14.0f, t = h / step, f = t - floorf(t);
        const float terraced = (floorf(t) + smooth(0.72f, 0.92f, f)) * step;
        dh += (terraced - h) * badlands;
    }
    if (swamp > 0.0f) {
        // Swamps: pulled down toward a low, gently rolling floor 2-5 m above the sea.
        const float floorH = 3.5f + 1.5f * biomeJitter.fractal(1, x, y, z);
        if (h > floorH) dh += (floorH - h) * 0.85f * swamp;
    }
    return dh;
}

float PlanetField::waterElevation(Vector3 dir, float fp) const {
    if (moon) return kNoWater;
    float w = kNoWater;
    if (const RiverNetwork* net = rivers()) {
        // A lake's water is its level over its whole region (the terrain decides the shore); rivers
        // in it take its level. Lakes are big: kept out to coarse levels of detail.
        float lake;
        if (fp < 1024.0f && net->lakeAt(dir, lake)) w = lake;
        else if (fp < 128.0f) {
            const RiverNetwork::Hit rh = net->nearest(dir, 60.0f);
            if (rh.found && rh.dist < rh.halfWidth * 1.05f) w = rh.level;
        }
    }
    if (fp >= 96.0f) return w; // ponds are not resolved (their bowls fade out by here)
    ForPondCellsNear(Vector3Scale(dir, radius), kPondReach, [&](int x, int y, int z) {
        const Pond p = pondOfCell(x, y, z);
        if (p.exists && !pondClaimed(PondId(x, y, z)) && Vector3Distance(dir, p.dir) * radius < p.r * 1.05f) w = std::max(w, p.level);
    });
    return w;
}

bool PlanetField::inlandWaterNear(Vector3 dir, float extent, float fp) const {
    if (moon) return false;
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
        if (p.exists && !pondClaimed(PondId(x, y, z)) && Vector3Distance(dir, p.dir) * radius < p.r * 1.05f + extent) near = true;
    });
    return near;
}

PlanetField::WaterSource PlanetField::waterSourceAt(Vector3 dir) const {
    WaterSource s;
    if (moon) return s;
    if (const RiverNetwork* net = rivers()) {
        float lake;
        if (net->lakeAt(dir, lake)) { s.kind = 2; s.level = lake; }
        else {
            const RiverNetwork::Hit rh = net->nearest(dir, 60.0f);
            if (rh.found && rh.dist < rh.halfWidth * 1.05f) { s.kind = 1; s.level = rh.level; }
        }
    }
    ForPondCellsNear(Vector3Scale(dir, radius), kPondReach, [&](int x, int y, int z) {
        const Pond p = pondOfCell(x, y, z);
        if (!p.exists || p.level <= s.level || Vector3Distance(dir, p.dir) * radius >= p.r * 1.05f) return;
        const uint64_t id = PondId(x, y, z);
        if (pondClaimed(id)) return;
        s.kind = 3; s.level = p.level; s.pond = id; s.pondDir = p.dir; s.pondR = p.r;
    });
    return s;
}

void PlanetField::claimPond(uint64_t id) const {
    std::lock_guard<std::mutex> lock(pondState->m);
    if (pondState->claimed.insert(id).second) pondState->n++;
}
bool PlanetField::pondClaimed(uint64_t id) const {
    if (pondState->n.load(std::memory_order_relaxed) == 0) return false;
    std::lock_guard<std::mutex> lock(pondState->m);
    return pondState->claimed.count(id) > 0;
}
std::vector<uint64_t> PlanetField::claimedPonds() const {
    std::lock_guard<std::mutex> lock(pondState->m);
    return std::vector<uint64_t>(pondState->claimed.begin(), pondState->claimed.end());
}
void PlanetField::clearClaimedPonds() const {
    std::lock_guard<std::mutex> lock(pondState->m);
    pondState->claimed.clear();
    pondState->n = 0;
}

PlanetField::RiverState::~RiverState() = default;

const RiverNetwork* PlanetField::rivers() const {
    if (moon || !riverState) return nullptr;
    std::call_once(riverState->once, [this] {
        riverState->net = std::make_unique<RiverNetwork>(radius, riverSeed, [this](Vector3 d) { return land(d, 256.0f).h; });
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
    thread_local std::unordered_map<int64_t, IslandInfo> cache;
    const int64_t key = ((int64_t)(cx + 32768) << 32) ^ ((int64_t)(cy + 32768) << 16) ^ (int64_t)(cz + 32768) ^ ((int64_t)pf.oreSeed << 48);
    auto it = cache.find(key);
    if (it == cache.end()) {
        if (cache.size() > 4096) cache.clear();
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
        it = cache.emplace(key, is).first;
    }
    out = it->second;
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
    // Caves, from a few metres down to 20 m above the mantle (which stays intact). Each cave shape is
    // its own distance field in metres, combined with min(): rock density grows with depth, so
    // subtracting a fixed amount would only carve shallow caves.
    // - Worm tunnels (~14 m radius) where both cave fields cross zero: at 1/150 m the fields change
    //   ~0.01 per metre, so (t - kTunnelLevel) * 100 ~ metres from the tunnel wall.
    // - Caverns: rare chambers tens of metres across where a low-frequency field peaks, 30 m+ down.
    // Coarse LOD samples can't resolve them, so both fade out there.
    const float caveWeight = std::clamp(2.0f - fp / 8.0f, 0.0f, 1.0f);
    // Cave mouths (Plan 5 C): in some regions (~15%) tunnels may reach up through the surface, so caves are found
    // by walking; elsewhere they stay 5 m+ down. Only looked up for shallow samples.
    float caveTop = 5.0f;
    if (depth < 14.0f && depth > -3.0f && caveWeight > 0.0f && !moon) caveTop = 5.0f - 7.5f * Smoothstep(0.3f, 0.45f, caveMouths.fractal(1, p.x, p.y, p.z));
    if (caveWeight > 0.0f && depth > caveTop && d > mantleRadius + 20.0f) {
        const float near = 2.0f * std::max(fp, 1.0f) + 2.0f; // farther than this from a cave = plain rock for meshing
        float m = 1e9f;                                        // metres to the nearest cave wall (< 0 inside)
        // A tunnel needs both fields near zero: if the first alone rules one out, skip the second.
        const float a = fabsf(cavesA.fractal(2, p.x, p.y, p.z));
        if ((a - kTunnelLevel) * 100.0f <= near) {
            const float fade = std::min(1.0f, (depth - caveTop) / 8.0f) * std::min(1.0f, (d - mantleRadius - 20.0f) / 10.0f);
            const float t = std::max(a, fabsf(cavesB.fractal(2, p.x, p.y, p.z)));
            m = (t - kTunnelLevel) * 100.0f + (1.0f - fade * caveWeight) * 20.0f;
        }
        if (depth > 30.0f && d > mantleRadius + 40.0f) {
            const float fade = std::min(1.0f, (depth - 30.0f) / 20.0f) * std::min(1.0f, (d - mantleRadius - 40.0f) / 20.0f);
            const float c = caverns.fractal(2, p.x, p.y * 1.4f, p.z); // a little flattened: wider than tall
            m = std::min(m, (kCavernLevel - c) * 120.0f + (1.0f - fade * caveWeight) * 30.0f);
        }
        if (m < near) v = std::min(v, ISO_LEVEL + DENSITY_GRADIENT * m);
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
            const Column3D col = terrain3DColumn(gdir, e); // the 3D band's column strength, once per node
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
    col.f = (c0[0].f * (1 - au) + c0[1].f * au) * (1 - av) + (c1[0].f * (1 - au) + c1[1].f * au) * av;
    return true;
}

void PlanetField::sampleBlock(Vector3 o, float stride, int S, float fp, const BlockGrid& g, float* out) const {
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
    const float x = dir.x * radius, y = dir.y * radius, z = dir.z * radius;
    // Warm at the equator, cold at the poles; 6.5 C colder per km of height; regional patches +-7 C.
    c.temp = 30.0f - 44.0f * (lat * sqrtf(lat)) - 6.5f * std::max(elev, 0.0f) / 1000.0f + 7.0f * tempNoise.fractal(2, x, y, z);
    // Moisture: regional wet and dry patches, wetter in the lowlands near the sea.
    c.moist = 0.5f + 0.85f * moistNoise.fractal(3, x, y, z) + 0.08f * (1.0f - Smoothstep(10.0f, 150.0f, elev)); // noise sits mostly in +-0.4
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
            for (float e = -6000.0f; e <= 6000.0f; e += 24.0f)
                for (float nn = -6000.0f; nn <= 6000.0f; nn += 24.0f) {
                    const float dist = sqrtf(e * e + nn * nn);
                    if (dist >= bestDist) continue;
                    const Vector3 d = Vector3Normalize(Vector3Add(best, Vector3Add(Vector3Scale(t1, e / radius), Vector3Scale(t2, nn / radius))));
                    const WaterSource ws = waterSourceAt(d);
                    if (ws.kind != 3) continue;
                    bestDist = dist;
                    const Vector3 out = Vector3Normalize(Vector3CrossProduct(ws.pondDir, { 0, 1, 0 }));
                    rim = Vector3Normalize(Vector3Add(ws.pondDir, Vector3Scale(out, ws.pondR / radius)));
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
    const float slack = fp > 1.0f ? 2.0f * fp : 0.0f;
    sky = std::clamp(1.0f - (depth - 2.0f - slack) / 8.0f, 0.0f, 1.0f); // buried more than ~10 m: no sky (caves, deep pits)
    if (moon) return depth < 2.5f + slack ? VoxelMaterial::Regolith : VoxelMaterial::Rock; // grey dust over rock (asteroid palette)
    if (depth < -15.0f) return Vector3DotProduct(normal, up) > 0.35f ? VoxelMaterial::Grass : VoxelMaterial::Rock; // floating islands: turf on top, rock beneath
    // Outcrops and boulders are bare rock (snow still caps them up high, below).
    if (depth < 2.5f + slack && surfaceElevation < 1050.0f && rockFeature(p, d, radius + surfaceElevation, fp) > -0.35f) return VoxelMaterial::Rock;
    if (depth < 2.5f + slack) {
        // Surface biome from height above sea level, latitude, slope and a climate wobble. The top
        // ~0.6 m is the surface material; below it a soil layer (dirt under grass, sand under sand,
        // gravel under scree), then rock.
        const float wobble = climate.fractal(2, p.x, p.y, p.z);
        const float lat = fabsf(up.y);                                     // 0 at the equator, 1 at the poles
        const float slope = 1.0f - std::clamp(Vector3DotProduct(normal, up), -1.0f, 1.0f); // 0 flat .. 1 vertical
        const float e = surfaceElevation;
        const bool skin = depth < 0.6f + slack;
        if (lat > 0.93f + 0.03f * wobble || e > 1350.0f + 100.0f * wobble) return VoxelMaterial::Ice;          // polar cap, glaciers
        if (lat > 0.82f + 0.08f * wobble || e > 1050.0f + 150.0f * wobble) return skin ? VoxelMaterial::Snow : VoxelMaterial::Ice;
        // Inland water: sand on the beds and wet shores, gravel in the deeper middle.
        if (e > 0.0f && inlandWaterNear(up, 0.0f)) {
            const float w = waterElevation(up, fp);
            if (w != kNoWater && e < w + 0.6f) return w - e > 1.4f ? VoxelMaterial::Gravel : VoxelMaterial::Sand;
        }
        if (slope > 0.45f) return VoxelMaterial::Rock;                                                       // cliffs: bare stone
        if (e < -80.0f) return VoxelMaterial::Gravel;                                                        // deep sea floor
        if (e < 6.0f + 4.0f * wobble) return VoxelMaterial::Sand;                                            // beaches, shallows
        // Biomes (Plan 5 B). The terrain shader tints each by the climate (red sand in the badlands,
        // dry grass on the savanna, dark mud in swamps).
        const Biome biome = biomeAt(up, e);
        if (biome == Biome::Desert) return VoxelMaterial::Sand;                                              // dunes all the way down
        if (biome == Biome::Badlands) return slope > 0.18f ? VoxelMaterial::Rock : skin ? VoxelMaterial::Sand : VoxelMaterial::Rock; // red sand over banded rock
        if (slope > 0.3f || e > 800.0f + 100.0f * wobble) return VoxelMaterial::Gravel;                     // scree, alpine
        if (biome == Biome::Tundra) return skin ? VoxelMaterial::Gravel : VoxelMaterial::Regolith;          // tundra: stony ground over dirt
        if (biome == Biome::Swamp && skin && rocks.fractal(1, p.x * 0.12f, p.y * 0.12f, p.z * 0.12f) > 0.15f) return VoxelMaterial::Regolith; // mud flats between the reeds
        // Grassland. Undug ground facing up is turf all the way down to the soil layer (the mesh strays
        // from the height map, which used to bare patches of dirt); dug ground is dirt.
        const bool natural = base < ISO_LEVEL + DENSITY_GRADIENT * (1.0f + slack);
        if (natural && slope < 0.3f) return VoxelMaterial::Grass;
        return skin ? VoxelMaterial::Grass : VoxelMaterial::Regolith;                                        // grassland over dirt
    }
    // Ore: deposits of one kind, and scattered pockets of the common ores (oreKind).
    if (depth > 3.0f && oreKind(p, depth) >= 0) return VoxelMaterial::Ore;
    return VoxelMaterial::Rock;
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

void BuildChunkMeshFromBase(const DensityField& field, const TerrainDeltas& deltas, ChunkMeshData& out, const float* base, float sampledMs, const PlanetField::BlockGrid* grid) {
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
    // Player edits: the delta lattice at the samples (lattice points at every stride; coarse nodes
    // point-sample it, so big excavations show from afar). Cost depends on the edited chunks here,
    // never on how many strokes the world has had.
    std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>> near;
    deltas.chunksIn(eo, Vector3AddValue(eo, (E - 1) * stride), near);
    if (!near.empty()) {
        // Lattice indices (kVoxelSize apart): a node of level L samples every 2^L-th lattice point.
        const int st = out.key.latticeStride(), ox = (int)lroundf(eo.x / kVoxelSize), oy = (int)lroundf(eo.y / kVoxelSize), oz = (int)lroundf(eo.z / kVoxelSize);
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
                VoxelMaterial m = hinted ? gridField->materialSkyHinted(p, b, stride, n, sky, &gElev, &gCol) : field.materialSky(p, b, stride, n, sky);
                if (!near.empty()) { const int placed = deltas.materialAt(p); if (placed >= 0) m = (VoxelMaterial)placed; } // placed ground keeps its material
                out.tangents[v * 4 + 3] = sky;
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
                const int mi = (int)m; // one-hot weight in the material's channel (see VoxelMaterial)
                if (mi < 4) out.colors[v * 4 + mi] = 255;
                else if (mi < 6) out.texcoords[v * 2 + (mi - 4)] = 1.0f;
                else out.texcoords2[v * 2 + (mi - 6)] = 1.0f;
            }
            out.indices.assign(idx.begin(), idx.end());
        }
    }

    // Sea surface: marching cubes on the sea-level sphere over the same samples, keeping only the
    // triangles that are in open air (their centre is not inside rock), so it shows over oceans and
    // lakes but not through land. Runs even when the node holds no terrain surface (open ocean).
    const float sea = field.seaLevelRadius();
    if (sea > 0.0f) {
        const Vector3 lo = origin, hi = Vector3AddValue(origin, cells * stride);
        const Vector3 nearest = { std::clamp(0.0f, lo.x, hi.x), std::clamp(0.0f, lo.y, hi.y), std::clamp(0.0f, lo.z, hi.z) };
        const float dmin = Vector3Length(nearest);
        const float dmax = Vector3Length({ std::max(fabsf(lo.x), fabsf(hi.x)), std::max(fabsf(lo.y), fabsf(hi.y)), std::max(fabsf(lo.z), fabsf(hi.z)) });
        // Inland water (rivers, ponds) above the sea: the surface follows the local water level. It
        // is sampled on a lattice every 4 samples and interpolated (it changes slowly), and only
        // for nodes near some.
        const PlanetField* pf = dynamic_cast<const PlanetField*>(&field);
        const Vector3 mid = Vector3Scale(Vector3Add(lo, hi), 0.5f);
        const float extent = Vector3Distance(lo, hi) * 0.5f;
        const bool inland = pf && !pf->moon && dmax > sea && dmin < sea + pf->maxElevation && Vector3Length(mid) > 1.0f &&
                            pf->inlandWaterNear(Vector3Normalize(mid), extent, (float)stride);
        if ((dmin < sea && dmax > sea) || inland) {
            std::vector<float> wd((size_t)S * S * S);
            std::vector<float> level; // water surface radius on the lattice
            const int L = 4, LN = (S - 1 + L - 1) / L + 1;
            auto lat = [&](int i) { return std::min(i * L, S - 1); };
            if (inland) {
                level.resize((size_t)LN * LN * LN);
                // Raw levels on the lattice and one node beyond it on every side (the neighbouring
                // chunks' nodes: lattices of chunks at the same level line up), so a dry node on the
                // chunk's edge sees wet ones next door and both chunks' surfaces meet at the seam.
                const int LE = LN + 2;
                auto coord = [&](int i) { return i < 0 ? -L : i >= LN ? lat(LN - 1) + L : lat(i); };
                std::vector<float> raw((size_t)LE * LE * LE);
                auto R = [&](int x, int y, int z) -> float& { return raw[((size_t)(z + 1) * LE + (y + 1)) * LE + (x + 1)]; };
                auto P = [&](int x, int y, int z) { return Vector3{ origin.x + coord(x) * stride, origin.y + coord(y) * stride, origin.z + coord(z) * stride }; };
                for (int z = -1; z <= LN; ++z)
                    for (int y = -1; y <= LN; ++y)
                        for (int x = -1; x <= LN; ++x)
                            R(x, y, z) = pf->waterElevation(Vector3Normalize(P(x, y, z)), (float)stride);
                // A dry point next to wet ones takes the lower of their level and just under its own
                // ground, not sea level: between them the interpolated surface then stays level out to
                // where the shore rises, or runs just under low ground, instead of sloping down across
                // it (sheets of water at chunk corners). Dry points with no wet neighbour: sea level.
                for (int z = 0; z < LN; ++z)
                    for (int y = 0; y < LN; ++y)
                        for (int x = 0; x < LN; ++x) {
                            float w = R(x, y, z);
                            if (w == PlanetField::kNoWater) {
                                float wet = PlanetField::kNoWater;
                                for (int c = -1; c <= 1; ++c)
                                    for (int b = -1; b <= 1; ++b)
                                        for (int a = -1; a <= 1; ++a) wet = std::max(wet, R(x + a, y + b, z + c));
                                w = wet == PlanetField::kNoWater ? 0.0f : std::min(wet, pf->elevation(Vector3Normalize(P(x, y, z)), (float)stride) - 2.0f);
                            }
                            level[((size_t)z * LN + y) * LN + x] = sea + std::max(0.0f, w);
                        }
            }
            auto levelAt = [&](int x, int y, int z) {
                if (!inland) return sea;
                const int ix = std::min(x / L, LN - 2), iy = std::min(y / L, LN - 2), iz = std::min(z / L, LN - 2);
                const float fx = (float)(x - lat(ix)) / std::max(1, lat(ix + 1) - lat(ix)), fy = (float)(y - lat(iy)) / std::max(1, lat(iy + 1) - lat(iy)), fz = (float)(z - lat(iz)) / std::max(1, lat(iz + 1) - lat(iz));
                auto g = [&](int a, int b, int c) { return level[((size_t)(iz + c) * LN + (iy + b)) * LN + (ix + a)]; };
                const float c00 = g(0, 0, 0) * (1 - fx) + g(1, 0, 0) * fx, c10 = g(0, 1, 0) * (1 - fx) + g(1, 1, 0) * fx;
                const float c01 = g(0, 0, 1) * (1 - fx) + g(1, 0, 1) * fx, c11 = g(0, 1, 1) * (1 - fx) + g(1, 1, 1) * fx;
                return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz;
            };
            size_t k = 0;
            for (int z = 0; z < S; ++z)
                for (int y = 0; y < S; ++y)
                    for (int x = 0; x < S; ++x, ++k) {
                        const Vector3 p = { origin.x + x * stride, origin.y + y * stride, origin.z + z * stride };
                        wd[k] = ISO_LEVEL + DENSITY_GRADIENT * (levelAt(x, y, z) - sqrtf(p.x * p.x + p.y * p.y + p.z * p.z));
                    }
            std::vector<Vector3> wv;
            std::vector<int> wi;
            std::vector<int> cache((size_t)S * S * S * 3, -1);
            for (int z = 0; z < cells; ++z)
                for (int y = 0; y < cells; ++y)
                    for (int x = 0; x < cells; ++x)
                        marchCube(x, y, z, wd, S, wv, wi, cache, ISO_LEVEL);
            // Terrain density at a point in sample coordinates (trilinear over the terrain grid).
            auto terrainAt = [&](Vector3 g) {
                const int x0 = std::clamp((int)floorf(g.x), 0, S - 2), y0 = std::clamp((int)floorf(g.y), 0, S - 2), z0 = std::clamp((int)floorf(g.z), 0, S - 2);
                const float fx = g.x - x0, fy = g.y - y0, fz = g.z - z0;
                auto d = [&](int x, int y, int z) { return dens[((size_t)z * S + y) * S + x]; };
                const float c00 = d(x0, y0, z0) * (1 - fx) + d(x0 + 1, y0, z0) * fx, c10 = d(x0, y0 + 1, z0) * (1 - fx) + d(x0 + 1, y0 + 1, z0) * fx;
                const float c01 = d(x0, y0, z0 + 1) * (1 - fx) + d(x0 + 1, y0, z0 + 1) * fx, c11 = d(x0, y0 + 1, z0 + 1) * (1 - fx) + d(x0 + 1, y0 + 1, z0 + 1) * fx;
                return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz;
            };
            const PlanetField* planetField = dynamic_cast<const PlanetField*>(&field);
            // Height-map elevation under a body point: from the block's grid when there is one (#148), as the terrain.
            auto groundAt = [&](Vector3 point, Vector3 dir) {
                float e; PlanetField::Column3D col;
                if (grid && planetField && planetField->gridSurface(*grid, point, e, col)) return e;
                return planetField->elevation(dir, (float)stride);
            };
            std::vector<int> remap(wv.size(), -1);
            for (size_t t = 0; t + 2 < wi.size(); t += 3) {
                const Vector3 c = Vector3Scale(Vector3Add(Vector3Add(wv[wi[t]], wv[wi[t + 1]]), wv[wi[t + 2]]), 1.0f / 3.0f);
                if (terrainAt(c) >= ISO_LEVEL) continue; // under land
                if (planetField) {
                    // Underground: where the original ground stands over the water (a cave or tunnel under a
                    // river, or below sea level) the water surface is not open water. Rivers used to line the
                    // cave walls.
                    const Vector3 wc = Vector3Add(origin, Vector3Scale(c, stride));
                    const float wl = Vector3Length(wc);
                    if (wl > 1.0f && groundAt(wc, Vector3Scale(wc, 1.0f / wl)) > wl - sea + (stride <= 2.0f ? 0.15f : 1.0f)) continue; // near: only over the original bed (dug ground beside it is the simulation's)
                }
                for (int e = 0; e < 3; ++e) {
                    int& r = remap[wi[t + e]];
                    if (r < 0) {
                        r = (int)(out.waterPositions.size() / 3);
                        const Vector3 lp = Vector3Scale(wv[wi[t + e]], stride);
                        const Vector3 n = Vector3Normalize(Vector3Add(origin, lp)); // the sea is a sphere: radial normals
                        out.waterPositions.insert(out.waterPositions.end(), { lp.x, lp.y, lp.z });
                        out.waterNormals.insert(out.waterNormals.end(), { n.x, n.y, n.z });
                        // Depth to the sea floor from the height map (the shader grades clarity and colour by it).
                        const float surfaceE = Vector3Length(Vector3Add(origin, lp)) - sea; // the water surface's height above sea level
                        const float depth = planetField ? std::max(0.0f, surfaceE - groundAt(Vector3Add(origin, lp), n)) : 50.0f;
                        out.waterDepth.insert(out.waterDepth.end(), { depth, 0.0f });
                    }
                    out.waterIndices.push_back((unsigned short)r);
                }
            }
            if (out.waterPositions.size() / 3 > 65535) { out.waterPositions.clear(); out.waterNormals.clear(); out.waterDepth.clear(); out.waterIndices.clear(); }
        }
    }
    ScatterGrass(field, out); // foliage: grass tufts on the finest level (in the build time)
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
    const float threshold = 0.72f + 0.12f * std::clamp((depth - 20.0f) / 200.0f, 0.0f, 1.0f);
    if (pocket > threshold) {
        uint32_t h = OreHash((int)floorf(p.x / 8.0f), (int)floorf(p.y / 8.0f), (int)floorf(p.z / 8.0f), oreSeed ^ 0x5eedu);
        const float r = OreUnit(h);
        return r < 0.35f ? 0 : r < 0.60f ? 1 : r < 0.75f ? 2 : 3; // iron, copper, tin, and coal (the seams near the surface)
    }
    return -1;
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

const char* PlanetField::stoneKind(Vector3 p) const {
    if (moon) return "basalt"; // the maria: dark basalt plains
    const float d = Vector3Length(p);
    if (d < 1.0f) return "stone";
    const Vector3 dir = Vector3Scale(p, 1.0f / d);
    const float surface = radius + elevation(dir, 8.0f);
    const float depth = surface - d;
    // Marble: small pockets below 20 m (the vein noise at a coarse scale).
    if (depth > 20.0f && fabsf(veins.fractal(1, p.x * 0.2f, p.y * 0.2f, p.z * 0.2f)) < 0.05f) return "marble";
    if (depth > 60.0f) return "basalt";
    const Land l = land(dir, 32.0f);
    if (l.h < 0.0f) return "basalt";            // the sea floor
    if (l.mountain > 0.4f) return "granite";    // the ranges
    if (l.h < 25.0f) return "sandstone";        // low coastal land
    // Broad regions of granite and limestone elsewhere.
    return climate.fractal(1, p.x * 0.3f, p.y * 0.3f, p.z * 0.3f) > 0.15f ? "granite" : "limestone";
}

const char* GemForStone(const char* stone) {
    if (!stone) return nullptr;
    if (strcmp(stone, "granite") == 0) return "amethyst";
    if (strcmp(stone, "marble") == 0) return "ruby";
    if (strcmp(stone, "basalt") == 0) return "sapphire";
    if (strcmp(stone, "limestone") == 0) return "emerald";
    return nullptr;
}
