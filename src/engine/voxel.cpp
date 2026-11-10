#include "voxel.h"
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
    // Bounds of elevation() below (each noise term is in [-1, 1]).
    if (moon) { minElevation = -90.0f; maxElevation = 340.0f; }
    else { minElevation = -740.0f; maxElevation = 1690.0f; }
}

float PlanetField::elevation(Vector3 dir, float fp) const {
    const float sx = dir.x * radius, sy = dir.y * radius, sz = dir.z * radius; // the surface point (m)
    if (moon) {
        // Grey highlands and lowland plains, rough ridged uplands, no oceans.
        const float base = continents.fractalLimited(4, fp, sx, sy, sz) * 60.0f;
        const float r = 1.0f - fabsf(ridges.fractalLimited(4, fp, sx, sy, sz));
        return base + r * r * 250.0f + hills.fractalLimited(3, fp, sx, sy, sz) * 25.0f + detail.fractalLimited(2, fp, sx, sy, sz) * 2.5f;
    }
    const float c = continents.fractalLimited(4, fp, sx, sy, sz); // < 0: ocean, > 0: land
    const float base = c < 0.0f ? c * 650.0f : c * 200.0f;        // ocean floor down to -650 m, gentle land
    const float land = Smoothstep(-0.02f, 0.06f, c);
    const float mountainMask = Smoothstep(0.08f, 0.35f, c);      // ranges only well inland
    float h = base;
    if (mountainMask > 0.0f) {
        const float r = 1.0f - fabsf(ridges.fractalLimited(5, fp, sx, sy, sz)); // ridged: sharp crests
        h += r * r * r * 1400.0f * mountainMask;
    }
    h += hills.fractalLimited(3, fp, sx, sy, sz) * (25.0f + 60.0f * mountainMask);
    // River valleys (~150 m wide) fade out once samples are too far apart to resolve them.
    const float valleyWeight = land * std::clamp(2.0f - fp / 64.0f, 0.0f, 1.0f);
    if (valleyWeight > 0.0f) {
        const float v = fabsf(valleys.fractal(3, sx, sy, sz));
        if (v < 0.05f) { const float t = 1.0f - v / 0.05f; h -= t * t * 70.0f * valleyWeight; }
    }
    // Ravines: rare narrow canyons (~60 m wide, up to ~90 m deep, steep sqrt-profile walls) in some
    // land regions, faded out where samples are too far apart to resolve them.
    const float ravineWeight = land * std::clamp(2.0f - fp / 24.0f, 0.0f, 1.0f);
    if (ravineWeight > 0.0f) {
        const float region = Smoothstep(0.25f, 0.4f, ravineMask.fractal(2, sx, sy, sz));
        if (region > 0.0f) {
            const float r = fabsf(ravines.fractal(2, sx, sy, sz));
            if (r < 0.006f) h -= sqrtf(1.0f - r / 0.006f) * 90.0f * region * ravineWeight;
        }
    }
    h += detail.fractalLimited(2, fp, sx, sy, sz) * 2.5f;
    // Polar ice shelves: near the poles the sea is frozen, so low ground is raised to an ice sheet a
    // couple of metres above sea level (the material rules make it ice).
    const float iceCap = Smoothstep(0.86f, 0.92f, fabsf(dir.y) + 0.03f * climate.fractalLimited(2, fp, sx, sy, sz));
    if (iceCap > 0.0f && h < 2.0f) h += (2.0f - h) * iceCap;
    return h;
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
    (void)base;
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
        if (slope > 0.45f) return VoxelMaterial::Rock;                                                       // cliffs: bare stone
        if (e < -80.0f) return VoxelMaterial::Gravel;                                                        // deep sea floor
        if (e < 6.0f + 4.0f * wobble) return VoxelMaterial::Sand;                                            // beaches, shallows
        if (slope > 0.3f || e > 800.0f + 100.0f * wobble) return VoxelMaterial::Gravel;                     // scree, alpine
        if (lat > 0.62f + 0.06f * wobble) return VoxelMaterial::Regolith;                                    // tundra: bare dirt
        return skin ? VoxelMaterial::Grass : VoxelMaterial::Regolith;                                        // grassland over dirt
    }
    // Ore veins (thin sheets where the vein field crosses zero), richer with depth.
    if (depth > 8.0f) {
        const float thickness = depth < 150.0f ? 0.025f : depth < 500.0f ? 0.04f : 0.055f;
        if (fabsf(veins.fractal(2, p.x, p.y, p.z)) < thickness) return VoxelMaterial::Ore;
    }
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
        if (dmin < sea && dmax > sea) {
            std::vector<float> wd((size_t)S * S * S);
            size_t k = 0;
            for (int z = 0; z < S; ++z)
                for (int y = 0; y < S; ++y)
                    for (int x = 0; x < S; ++x, ++k) {
                        const Vector3 p = { origin.x + x * stride, origin.y + y * stride, origin.z + z * stride };
                        wd[k] = ISO_LEVEL + DENSITY_GRADIENT * (sea - sqrtf(p.x * p.x + p.y * p.y + p.z * p.z));
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
                for (int e = 0; e < 3; ++e) {
                    int& r = remap[wi[t + e]];
                    if (r < 0) {
                        r = (int)(out.waterPositions.size() / 3);
                        const Vector3 lp = Vector3Scale(wv[wi[t + e]], stride);
                        const Vector3 n = Vector3Normalize(Vector3Add(origin, lp)); // the sea is a sphere: radial normals
                        out.waterPositions.insert(out.waterPositions.end(), { lp.x, lp.y, lp.z });
                        out.waterNormals.insert(out.waterNormals.end(), { n.x, n.y, n.z });
                        // Depth to the sea floor from the height map (the shader grades clarity and colour by it).
                        const float depth = planetField ? std::max(0.0f, -planetField->elevation(n, (float)stride)) : 50.0f;
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
