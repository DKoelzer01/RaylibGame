#include "voxel.h"
#include "rivers.h"
#include "foliage.h"
#include "trees.h"
#include "cubemarch.h"

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

static float Smoothstep(float a, float b, float x) {
    const float t = std::clamp((x - a) / (b - a), 0.0f, 1.0f);
    return t * t * (3.0f - 2.0f * t);
}

PlanetField::PlanetField(unsigned int seed, float radius_, bool moon_)
    : radius(radius_), mantleRadius(radius_ * (1.0f - 0.015f)), moon(moon_) {
    oreSeed = seed ^ 0x0e5e0e5u;
    riverSeed = seed;
    riverState = std::make_shared<RiverState>();
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
    // Bounds of elevation() below (each noise term is in [-1, 1]).
    if (moon) { minElevation = -90.0f; maxElevation = 340.0f; }
    else { minElevation = -740.0f; maxElevation = 1690.0f; }
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
    float h = base;
    if (mountainMask > 0.0f) {
        const float r = 1.0f - fabsf(ridges.fractalLimited(5, fp, sx, sy, sz)); // ridged: sharp crests
        l.ridge = r * r * r * 1400.0f * mountainMask; // the height the ranges add here
        h += l.ridge;
    }
    h += hills.fractalLimited(3, fp, sx, sy, sz) * (25.0f + 60.0f * mountainMask);
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

float PlanetField::elevation(Vector3 dir, float fp) const {
    const Land l = land(dir, fp);
    if (moon) return l.h;
    float h = l.h + l.detail;
    // Rivers of the drainage network: a valley (a floodplain 2.5 channel half-widths across at the
    // water, then sides rising 0.3 m per m) cut into higher ground, low banks where the ground dips
    // below the water, and the channel bed. The valley stays resolved further out than the channel.
    if (const RiverNetwork* net = rivers()) {
        // Next to a lake, the ground is banked up 0.3 m above its surface, fading out over 140 m, so
        // the water stays in (the outflow river below cuts through it).
        float lv, ld;
        if (net->lakeRim(dir, lv, ld) && h < lv + 0.3f) h += (lv + 0.3f - h) * (1.0f - Smoothstep(0.0f, 140.0f, ld));
        const RiverNetwork::Hit rh = net->nearest(dir, RiverNetwork::kMaxReach);
        if (rh.found) {
            const float d = rh.dist, half = rh.halfWidth;
            const float flat = 2.5f * half, reach = std::min(flat + 200.0f, RiverNetwork::kMaxReach);
            const float allowed = rh.level + 0.5f + std::max(0.0f, d - flat) * 0.3f;
            const float vw = (1.0f - Smoothstep(reach * 0.7f, reach, d)) * std::clamp(2.0f - fp / 512.0f, 0.0f, 1.0f);
            if (h > allowed) h += (allowed - h) * vw;
            const float cw = std::clamp(2.0f - fp / (4.0f * half), 0.0f, 1.0f) * Smoothstep(-2.0f, 0.5f, l.h); // not out at sea
            const float bank = rh.level + 0.4f;
            if (h < bank && d < half + 12.0f) h += (bank - h) * (1.0f - Smoothstep(half, half + 12.0f, d)) * cw;
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

float PlanetField::waterElevation(Vector3 dir, float fp) const {
    if (moon) return kNoWater;
    float w = kNoWater;
    if (const RiverNetwork* net = rivers()) {
        if (fp < 128.0f) {
            const RiverNetwork::Hit rh = net->nearest(dir, 60.0f);
            if (rh.found && rh.dist < rh.halfWidth * 1.05f) w = rh.level;
        }
        float lake;
        if (fp < 256.0f && net->lakeAt(dir, lake)) w = std::max(w, lake);
    }
    if (fp >= 96.0f) return w; // ponds are not resolved (their bowls fade out by here)
    ForPondCellsNear(Vector3Scale(dir, radius), kPondReach, [&](int x, int y, int z) {
        const Pond p = pondOfCell(x, y, z);
        if (p.exists && Vector3Distance(dir, p.dir) * radius < p.r * 1.05f) w = std::max(w, p.level);
    });
    return w;
}

bool PlanetField::inlandWaterNear(Vector3 dir, float extent, float fp) const {
    if (moon) return false;
    const float sx = dir.x * radius, sy = dir.y * radius, sz = dir.z * radius;
    // Each kind only where samples resolve it (as in elevation() and waterElevation()).
    if (fp < 256.0f) if (const RiverNetwork* net = rivers()) if (net->waterNear(dir, extent)) return true;
    if (fp >= 96.0f) return false;
    bool near = false;
    ForPondCellsNear({ sx, sy, sz }, kPondReach + extent, [&](int x, int y, int z) {
        if (near) return;
        const Pond p = pondOfCell(x, y, z);
        if (p.exists && Vector3Distance(dir, p.dir) * radius < p.r * 1.05f + extent) near = true;
    });
    return near;
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

// Everything after the height map, shared by single samples and blocks so both agree exactly.
float PlanetField::densityNear(Vector3 p, float d, float surface, float fp) const {
    float v = ISO_LEVEL + DENSITY_GRADIENT * (surface - d);
    const float depth = surface - d;
    if (fabsf(depth) < 6.0f) v += DENSITY_GRADIENT * 1.5f * detail.fractalLimited(2, fp, p.x, p.y * 1.3f, p.z); // crags, overhangs
    // Caves, from a few metres down to 20 m above the mantle (which stays intact). Each cave shape is
    // its own distance field in metres, combined with min(): rock density grows with depth, so
    // subtracting a fixed amount would only carve shallow caves.
    // - Worm tunnels (~14 m radius) where both cave fields cross zero: at 1/150 m the fields change
    //   ~0.01 per metre, so (t - kTunnelLevel) * 100 ~ metres from the tunnel wall.
    // - Caverns: rare chambers tens of metres across where a low-frequency field peaks, 30 m+ down.
    // Coarse LOD samples can't resolve them, so both fade out there.
    const float caveWeight = std::clamp(2.0f - fp / 8.0f, 0.0f, 1.0f);
    if (caveWeight > 0.0f && depth > 5.0f && d > mantleRadius + 20.0f) {
        const float near = 2.0f * std::max(fp, 1.0f) + 2.0f; // farther than this from a cave = plain rock for meshing
        float m = 1e9f;                                        // metres to the nearest cave wall (< 0 inside)
        // A tunnel needs both fields near zero: if the first alone rules one out, skip the second.
        const float a = fabsf(cavesA.fractal(2, p.x, p.y, p.z));
        if ((a - kTunnelLevel) * 100.0f <= near) {
            const float fade = std::min(1.0f, (depth - 5.0f) / 8.0f) * std::min(1.0f, (d - mantleRadius - 20.0f) / 10.0f);
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
    return v;
}

void PlanetField::baseDensityBlock(Vector3 o, float stride, int S, float fp, float* out) const {
    const float ext = stride * (float)(S - 1);
    const Vector3 c = { o.x + 0.5f * ext, o.y + 0.5f * ext, o.z + 0.5f * ext };
    const float dc = sqrtf(c.x * c.x + c.y * c.y + c.z * c.z), halfDiag = 0.8661f * ext;
    const float outer = radius + maxElevation, inner = radius + minElevation;
    // Blocks wholly outside the surface band take the cheap per-sample paths anyway; blocks too large
    // for a tangent-plane grid (near the centre) fall back to the per-sample loop too.
    if (dc - halfDiag > outer + 4.0f || dc + halfDiag < inner - 4.0f || halfDiag > 0.5f * dc) {
        DensityField::baseDensityBlock(o, stride, S, fp, out);
        return;
    }
    // Height map on a gnomonic grid tangent at the block centre's direction: every sample along the
    // same line of sight from the planet centre shares one height, so a 35^3 block needs ~35^2
    // height evaluations instead of 35^3. Bilinear interpolation at one stride spacing is accurate
    // to centimetres (the finest height detail kept at this footprint is >= 2 strides long).
    const Vector3 n = { c.x / dc, c.y / dc, c.z / dc };
    const Vector3 t1 = Vector3Normalize(fabsf(n.y) < 0.9f ? Vector3CrossProduct(n, { 0, 1, 0 }) : Vector3CrossProduct(n, { 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(n, t1);
    float umin = 1e30f, umax = -1e30f, vmin = 1e30f, vmax = -1e30f;
    for (int k = 0; k < 8; ++k) {
        const Vector3 q = { o.x + ((k & 1) ? ext : 0.0f), o.y + ((k & 2) ? ext : 0.0f), o.z + ((k & 4) ? ext : 0.0f) };
        const float s = Vector3DotProduct(q, n);
        const float u = Vector3DotProduct(q, t1) / s, v = Vector3DotProduct(q, t2) / s;
        umin = std::min(umin, u); umax = std::max(umax, u); vmin = std::min(vmin, v); vmax = std::max(vmax, v);
    }
    const float du = stride / radius;
    const int nu = std::min(256, (int)ceilf((umax - umin) / du) + 2), nv = std::min(256, (int)ceilf((vmax - vmin) / du) + 2);
    const float su = (umax - umin) / (float)(nu - 1), sv = (vmax - vmin) / (float)(nv - 1);
    std::vector<float> h((size_t)nu * nv);
    for (int j = 0; j < nv; ++j)
        for (int i = 0; i < nu; ++i) {
            const float u = umin + i * su, v = vmin + j * sv;
            h[(size_t)j * nu + i] = elevation(Vector3Normalize(Vector3Add(n, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v)))), fp);
        }
    size_t idx = 0;
    for (int z = 0; z < S; ++z)
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x, ++idx) {
                const Vector3 p = { o.x + x * stride, o.y + y * stride, o.z + z * stride };
                const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
                if (d > outer + 4.0f) { out[idx] = ISO_LEVEL + DENSITY_GRADIENT * (outer - d); continue; }
                float surface = inner;
                if (d >= inner - 4.0f) {
                    const float s = Vector3DotProduct(p, n);
                    const float fu = std::clamp((Vector3DotProduct(p, t1) / s - umin) / su, 0.0f, (float)(nu - 1) - 1e-3f);
                    const float fv = std::clamp((Vector3DotProduct(p, t2) / s - vmin) / sv, 0.0f, (float)(nv - 1) - 1e-3f);
                    const int iu = (int)fu, iv = (int)fv;
                    const float au = fu - iu, av = fv - iv;
                    const float* r0 = &h[(size_t)iv * nu + iu];
                    const float* r1 = r0 + nu;
                    surface = radius + (r0[0] * (1 - au) + r0[1] * au) * (1 - av) + (r1[0] * (1 - au) + r1[1] * au) * av;
                }
                out[idx] = densityNear(p, d, surface, fp);
            }
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
        const float score = -slope * 100.0f - e * 0.01f;
        if (score > bestScore) { bestScore = score; best = dir; }
        if (slope < 0.01f) break;
    }
    return best;
}

VoxelMaterial PlanetField::materialSky(Vector3 p, float base, float fp, Vector3 normal, float& sky) const {
    sky = 1.0f;
    const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    if (d < mantleRadius + 0.5f || d < 1e-3f) { sky = 0.0f; return VoxelMaterial::Rock; } // mantle: always deep underground
    // Depth below the original surface, from the height map (the density is no measure of depth
    // here: cave tunnels take the minimum with it, so cave walls stay rock rather than topsoil).
    // `base` (the unedited density) tells undug ground (near the iso level) from dug-out ground.
    const Vector3 up = { p.x / d, p.y / d, p.z / d };
    const float surfaceElevation = elevation(up, fp);
    const float depth = radius + surfaceElevation - d; // against the surface this LOD meshes
    // A coarse LOD mesh strays from the height map it approximates by up to a couple of samples, so
    // its "depth" must allow for that, or whole patches of surface read as buried (black, no sun) and
    // take underground materials (seen from orbit as black triangles on peaks and rings along the grid).
    const float slack = fp > 1.0f ? 2.0f * fp : 0.0f;
    sky = std::clamp(1.0f - (depth - 2.0f - slack) / 8.0f, 0.0f, 1.0f); // buried more than ~10 m: no sky (caves, deep pits)
    if (moon) return depth < 2.5f + slack ? VoxelMaterial::Regolith : VoxelMaterial::Rock; // grey dust over rock (asteroid palette)
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
        if (slope > 0.3f || e > 800.0f + 100.0f * wobble) return VoxelMaterial::Gravel;                     // scree, alpine
        if (lat > 0.62f + 0.06f * wobble) return skin ? VoxelMaterial::Gravel : VoxelMaterial::Regolith;    // tundra: stony ground over dirt
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

void BuildChunkMesh(const DensityField& field, const TerrainDeltas& deltas, ChunkMeshData& out) {
    const auto t0 = std::chrono::steady_clock::now();
    const int cells = CHUNK_SIZE + 1;        // +1 cell overlap with the positive neighbours (LOD crack cover)
    const int S = cells + 1;                 // samples per axis
    const float stride = (float)out.key.stride();
    const Vector3 origin = out.key.origin();

    // Sample one extra layer on every side (E = S + 2 per axis): normals are central differences of
    // this grid, so they need no extra field evaluations, and neighbouring chunks see the same
    // samples on their shared face (consistent lighting across seams). `base` keeps the unedited
    // values for materials ("original surface" vs "dug out").
    const int E = S + 2;
    std::vector<float> base((size_t)E * E * E);
    const Vector3 eo = Vector3AddValue(origin, -stride);
    field.baseDensityBlock(eo, stride, E, stride, base.data());
    std::vector<float> ext = base;
    // Player edits: the delta lattice at the samples (lattice points at every stride; coarse nodes
    // point-sample it, so big excavations show from afar). Cost depends on the edited chunks here,
    // never on how many strokes the world has had.
    std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>> near;
    deltas.chunksIn(eo, Vector3AddValue(eo, (E - 1) * stride), near);
    if (!near.empty()) {
        const int st = (int)stride, ox = (int)lroundf(eo.x), oy = (int)lroundf(eo.y), oz = (int)lroundf(eo.z);
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
                    ext[k] += c->d[(lz * kDeltaChunk + ly) * kDeltaChunk + lx] * kDeltaStep;
                }
            }
    }
    auto at = [&](const std::vector<float>& g, int x, int y, int z) { return g[((size_t)z * E + y) * E + x]; }; // extended coords

    std::vector<float> dens((size_t)S * S * S);
    bool anySolid = false, anyEmpty = false;
    size_t i = 0;
    for (int z = 0; z < S; ++z)
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x, ++i) {
                const float d = at(ext, x + 1, y + 1, z + 1);
                dens[i] = d;
                if (d >= ISO_LEVEL) anySolid = true; else anyEmpty = true;
            }

    out.sampleMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
    // Density gradient at an inner sample (inner coords), by central differences on the extended grid.
    auto gradient = [&](int x, int y, int z) {
        return Vector3{ at(ext, x + 2, y + 1, z + 1) - at(ext, x, y + 1, z + 1),
                        at(ext, x + 1, y + 2, z + 1) - at(ext, x + 1, y, z + 1),
                        at(ext, x + 1, y + 1, z + 2) - at(ext, x + 1, y + 1, z) };
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
                const VoxelMaterial m = field.materialSky(p, b, stride, n, sky);
                out.tangents[v * 4 + 3] = sky;
                if (m == VoxelMaterial::Ore) {
                    // Which ore (the shader colours each kind): 1 + kind, over 8.
                    if (const auto* pfo = dynamic_cast<const PlanetField*>(&field)) {
                        const float dl = Vector3Length(p);
                        const int k = dl > 1.0f ? pfo->oreKind(p, pfo->radius + pfo->elevation(Vector3Scale(p, 1.0f / dl), (float)stride) - dl) : -1;
                        out.tangents[v * 4 + 2] = (k + 1) / 8.0f;
                    }
                }
                if (planetLife && sky > 0.9f && (m == VoxelMaterial::Grass || m == VoxelMaterial::Regolith)) {
                    // Distant forests: canopy cover (x) and latitude (y) for the terrain shader, which
                    // paints them in where the real trees fade out.
                    const float d = Vector3Length(p);
                    const Vector3 up = Vector3Scale(p, 1.0f / d);
                    const float steep = 1.0f - Vector3DotProduct(n, up);
                    out.tangents[v * 4 + 0] = TreeCanopy(up, planetLife->radius, d - planetLife->radius, planetLife->foliageSeed) * std::clamp(1.0f - (steep - 0.1f) / 0.15f, 0.0f, 1.0f);
                    out.tangents[v * 4 + 1] = fabsf(up.y);
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
                for (int z = 0; z < LN; ++z)
                    for (int y = 0; y < LN; ++y)
                        for (int x = 0; x < LN; ++x) {
                            const Vector3 p = { origin.x + lat(x) * stride, origin.y + lat(y) * stride, origin.z + lat(z) * stride };
                            const float w = pf->waterElevation(Vector3Normalize(p), (float)stride);
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
                    if (wl > 1.0f && planetField->elevation(Vector3Scale(wc, 1.0f / wl), (float)stride) > wl - sea + 1.0f) continue;
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
                        const float depth = planetField ? std::max(0.0f, surfaceE - planetField->elevation(n, (float)stride)) : 50.0f;
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
