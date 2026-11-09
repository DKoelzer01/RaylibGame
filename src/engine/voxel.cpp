#include "voxel.h"
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

float DensityField::density(Vector3 p, const TerrainEdit* edits, size_t editCount, float footprint) const {
    float v = footprint > 0.0f ? baseDensityLod(p, footprint) : baseDensity(p);
    if (editCount == 0 || !editable(p)) return v; // e.g. a planet's mantle ignores edits
    for (size_t i = 0; i < editCount; ++i) {
        const TerrainEdit& e = edits[i];
        const float dx = p.x - e.center.x, dy = p.y - e.center.y, dz = p.z - e.center.z;
        const float d2 = dx * dx + dy * dy + dz * dz;
        if (d2 >= e.radius * e.radius) continue;
        const float t = 1.0f - sqrtf(d2) / e.radius;
        v += e.strength * t * t * (3.0f - 2.0f * t); // smoothstep falloff: soft-edged craters/mounds
    }
    // Clamp so repeated digging/adding stays reversible in a few clicks.
    return v; // no clamp: the uniform gradient keeps values meaningful at any depth
}

Vector3 DensityField::normal(Vector3 p, const TerrainEdit* edits, size_t n, float h, float fp) const {

    const Vector3 g = {
        density({ p.x + h, p.y, p.z }, edits, n, fp) - density({ p.x - h, p.y, p.z }, edits, n, fp),
        density({ p.x, p.y + h, p.z }, edits, n, fp) - density({ p.x, p.y - h, p.z }, edits, n, fp),
        density({ p.x, p.y, p.z + h }, edits, n, fp) - density({ p.x, p.y, p.z - h }, edits, n, fp),
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
    cavesA = SeededNoise(rng, 1.0f / 60.0f);       // tunnels run where both cave fields cross zero
    cavesB = SeededNoise(rng, 1.0f / 60.0f);
    veins = SeededNoise(rng, 1.0f / 30.0f);
    climate = SeededNoise(rng, 1.0f / 15000.0f);
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
    return h + detail.fractalLimited(2, fp, sx, sy, sz) * 2.5f;
}

float PlanetField::baseDensityLod(Vector3 p, float fp) const {
    const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    // Far above or far below the surface band only the sign matters (no surface to mesh), so skip
    // the height map: most samples of a chunk that does not touch the surface take these paths.
    const float outer = radius + maxElevation, inner = radius + minElevation;
    if (d > outer + 4.0f) return ISO_LEVEL + DENSITY_GRADIENT * (outer - d);
    const float surface = (d < inner - 4.0f || d < 1e-3f) ? inner : radius + elevation({ p.x / d, p.y / d, p.z / d }, fp);
    float v = ISO_LEVEL + DENSITY_GRADIENT * (surface - d);
    const float depth = surface - d;
    if (fabsf(depth) < 6.0f) v += DENSITY_GRADIENT * 1.5f * detail.fractalLimited(2, fp, p.x, p.y * 1.3f, p.z); // crags, overhangs
    // Caves: worm tunnels where both cave fields cross zero, from a few metres down to 20 m above
    // the mantle (which stays intact). The tunnel is its own distance field (the noise changes
    // ~0.025 per metre at this frequency, so t * 40 ~ metres from the tunnel's core), combined with
    // min(): rock density grows with depth, so subtracting a fixed amount would only carve shallow caves.
    // Tunnels a few metres wide can't be resolved by coarse LOD samples, so they fade out there.
    const float caveWeight = std::clamp(2.0f - fp / 4.0f, 0.0f, 1.0f);
    if (caveWeight > 0.0f && depth > 5.0f && d > mantleRadius + 20.0f) {
        const float t = std::max(fabsf(cavesA.fractal(2, p.x, p.y, p.z)), fabsf(cavesB.fractal(2, p.x, p.y, p.z)));
        const float fade = std::min(1.0f, (depth - 5.0f) / 8.0f) * std::min(1.0f, (d - mantleRadius - 20.0f) / 10.0f);
        const float tunnel = ISO_LEVEL + DENSITY_GRADIENT * ((t - 0.06f) * 40.0f + (1.0f - fade * caveWeight) * 20.0f);
        v = std::min(v, tunnel);
    }
    return v;
}

VoxelMaterial PlanetField::materialLod(Vector3 p, float base, float fp) const {
    const float d = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
    if (d < mantleRadius + 0.5f || d < 1e-3f) return VoxelMaterial::Rock;   // mantle
    // Depth below the original surface, from the height map (the density is no measure of depth
    // here: cave tunnels take the minimum with it, so cave walls stay rock rather than topsoil).
    (void)base;
    const float depth = radius + elevation({ p.x / d, p.y / d, p.z / d }, fp) - d; // against the surface this LOD meshes
    if (depth < 2.5f) {
        if (!moon && d > 1e-3f) {
            // Snow at high latitudes and on high peaks, wobbled by a climate field.
            const float wobble = climate.fractal(2, p.x, p.y, p.z);
            if (fabsf(p.y) / d > 0.82f + 0.08f * wobble || d - radius > 1050.0f + 150.0f * wobble) return VoxelMaterial::Ice;
        }
        return VoxelMaterial::Regolith;
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

void FilterEdits(const EditList& all, Vector3 lo, Vector3 hi, float margin, EditList& out) {
    out.clear();
    for (const TerrainEdit& e : all) {
        const float cx = std::max(lo.x, std::min(e.center.x, hi.x));
        const float cy = std::max(lo.y, std::min(e.center.y, hi.y));
        const float cz = std::max(lo.z, std::min(e.center.z, hi.z));
        const float dx = e.center.x - cx, dy = e.center.y - cy, dz = e.center.z - cz;
        const float r = e.radius + margin;
        if (dx * dx + dy * dy + dz * dz <= r * r) out.push_back(e);
    }
}

void BuildChunkMesh(const DensityField& field, const EditList& allEdits, ChunkMeshData& out) {
    const auto t0 = std::chrono::steady_clock::now();
    const int cells = CHUNK_SIZE + 1;        // +1 cell overlap with the positive neighbours (LOD crack cover)
    const int S = cells + 1;                 // samples per axis
    const float stride = (float)out.key.stride();
    const Vector3 origin = out.key.origin();

    // Margin = one stride: normals sample density half a stride outside the node.
    EditList edits;
    FilterEdits(allEdits, origin, Vector3AddValue(origin, cells * stride), stride, edits);

    std::vector<float> dens((size_t)S * S * S);
    bool anySolid = false, anyEmpty = false;
    size_t i = 0;
    for (int z = 0; z < S; ++z)
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x, ++i) {
                const float d = field.density({ origin.x + x * stride, origin.y + y * stride, origin.z + z * stride }, edits.data(), edits.size(), stride);
                dens[i] = d;
                if (d >= ISO_LEVEL) anySolid = true; else anyEmpty = true;
            }

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
            for (size_t v = 0; v < verts.size(); ++v) {
                const Vector3 lp = Vector3Scale(verts[v], stride);     // node-relative, in voxels
                const Vector3 p = Vector3Add(origin, lp);
                // Coarser levels use a wider gradient step: smoother, cheaper-looking normals far away.
                const Vector3 n = field.normal(p, edits.data(), edits.size(), 0.5f * stride, stride);
                const VoxelMaterial m = field.materialLod(p, field.baseDensity(p), stride);
                out.positions[v * 3 + 0] = lp.x;
                out.positions[v * 3 + 1] = lp.y;
                out.positions[v * 3 + 2] = lp.z;
                out.normals[v * 3 + 0] = n.x;
                out.normals[v * 3 + 1] = n.y;
                out.normals[v * 3 + 2] = n.z;
                out.colors[v * 4 + (int)m] = 255;
            }
            out.indices.assign(idx.begin(), idx.end());
        }
    }
    out.buildMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
}
