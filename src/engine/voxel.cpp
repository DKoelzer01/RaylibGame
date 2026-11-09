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
    return ISO_LEVEL + DENSITY_GRADIENT * (radius * 0.72f + displacement - d);
}

float PlanetoidField::density(Vector3 p, const TerrainEdit* edits, size_t editCount) const {
    float v = baseDensity(p);
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

Vector3 PlanetoidField::normal(Vector3 p, const TerrainEdit* edits, size_t n, float h) const {

    const Vector3 g = {
        density({ p.x + h, p.y, p.z }, edits, n) - density({ p.x - h, p.y, p.z }, edits, n),
        density({ p.x, p.y + h, p.z }, edits, n) - density({ p.x, p.y - h, p.z }, edits, n),
        density({ p.x, p.y, p.z + h }, edits, n) - density({ p.x, p.y, p.z - h }, edits, n),
    };
    const float len = sqrtf(g.x * g.x + g.y * g.y + g.z * g.z);
    if (len < 1e-6f) {
        const float pl = sqrtf(p.x * p.x + p.y * p.y + p.z * p.z);
        return pl > 1e-6f ? Vector3{ p.x / pl, p.y / pl, p.z / pl } : Vector3{ 0, 1, 0 };
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

void BuildChunkMesh(const PlanetoidField& field, const EditList& allEdits, ChunkMeshData& out) {
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
                const float d = field.density({ origin.x + x * stride, origin.y + y * stride, origin.z + z * stride }, edits.data(), edits.size());
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
                const Vector3 n = field.normal(p, edits.data(), edits.size(), 0.5f * stride);
                const VoxelMaterial m = field.material(p, field.baseDensity(p));
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
