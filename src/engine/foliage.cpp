#include "foliage.h"
#include "voxel.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdint>

bool foliageEnabled = true;

namespace {
    // Small, fast integer hash (lowbias32) -> [0, 1).
    uint32_t Hash(uint32_t x) {
        x ^= x >> 16; x *= 0x7feb352dU;
        x ^= x >> 15; x *= 0x846ca68bU;
        x ^= x >> 16;
        return x;
    }
    float Unit(uint32_t& state) { state = Hash(state + 0x9e3779b9U); return (state >> 8) * (1.0f / 16777216.0f); }

    Vector3 Vert(const std::vector<float>& v, int i) { return { v[3 * i], v[3 * i + 1], v[3 * i + 2] }; }

    // The shared tuft meshes, one after another in one buffer: y runs 0..1 from root to tip (the
    // shader scales it and bends the tips with the wind).
    //   0: a sparse tuft, five blades fanned around the up axis;
    //   1: a clump, 14 blades from roots spread over ~12 cm;
    //   2: a dense clump, 24 blades over ~20 cm, the fullest.
    unsigned int tuftVbo = 0;
    int variantFirst[kGrassVariants] = {}, variantCount[kGrassVariants] = {};
    void AddBlade(std::vector<float>& v, float rx, float rz, float dirA, float lean, float height, float w) {
        const float c = cosf(dirA), s = sinf(dirA);
        // Root left, root right (across the lean direction), tip leaning out along it.
        const float pts[3][3] = { { rx - s * w, 0.0f, rz + c * w }, { rx + s * w, 0.0f, rz - c * w }, { rx + c * lean, height, rz + s * lean } };
        for (auto& p : pts) v.insert(v.end(), { p[0], p[1], p[2] });
    }
    void EnsureTuftMesh() {
        if (tuftVbo) return;
        std::vector<float> v;
        // Variant 0: the original tuft.
        variantFirst[0] = 0;
        for (int b = 0; b < 5; ++b) {
            const float a = b * 2.0f * PI / 5 + 0.4f * b;
            AddBlade(v, 0.0f, 0.0f, a, 0.18f + 0.07f * (b % 3), 1.0f - 0.1f * (b % 2), 0.035f);
        }
        // Clumps: blades from roots scattered in a disc, leaning outward (a little random), the
        // middle ones tallest. A fixed seed: every clump of a variant is the same mesh, turned and
        // scaled per instance.
        uint32_t state = 12345u;
        const int blades[3] = { 5, 14, 24 };
        const float spread[3] = { 0.0f, 0.12f, 0.2f };
        for (int var = 1; var < kGrassVariants; ++var) {
            variantFirst[var] = (int)(v.size() / 3);
            for (int b = 0; b < blades[var]; ++b) {
                const float ra = Unit(state) * 2.0f * PI, rr = spread[var] * sqrtf(Unit(state));
                const float rx = cosf(ra) * rr, rz = sinf(ra) * rr;
                const float out = rr / std::max(spread[var], 1e-3f);           // 0 centre .. 1 rim
                const float dirA = ra + (Unit(state) - 0.5f) * 1.2f;
                const float lean = 0.08f + 0.25f * out + 0.1f * Unit(state);
                const float height = (0.65f + 0.35f * Unit(state)) * (1.0f - 0.25f * out);
                AddBlade(v, rx, rz, dirA, lean, height, 0.028f + 0.012f * Unit(state));
            }
        }
        for (int var = 0; var < kGrassVariants; ++var)
            variantCount[var] = (var + 1 < kGrassVariants ? variantFirst[var + 1] : (int)(v.size() / 3)) - variantFirst[var];
        tuftVbo = rlLoadVertexBuffer(v.data(), (int)(v.size() * sizeof(float)), false);
    }
}

void ScatterGrass(const DensityField& field, ChunkMeshData& out) {
    // Nodes at 1 m resolution or finer (the old level 0 and the finer ones around the player).
    if (out.key.stride() > 1.0f || out.indices.empty() || out.texcoords2.empty() || out.tangents.empty()) return;
    const float sea = field.seaLevelRadius();
    if (sea <= 0.0f) return; // no seas, no grass (moons, asteroids)
    const Vector3 origin = out.key.origin();
    const float edge = out.key.size(); // keep the overlap cell to the neighbour (it draws its own)
    const uint32_t nodeSeed = Hash((uint32_t)out.key.x * 73856093U ^ (uint32_t)out.key.y * 19349663U ^ (uint32_t)out.key.z * 83492791U);
    const size_t tris = out.indices.size() / 3;
    std::vector<float> byVariant[kGrassVariants];
    for (size_t t = 0; t < tris; ++t) {
        const int i0 = out.indices[3 * t], i1 = out.indices[3 * t + 1], i2 = out.indices[3 * t + 2];
        // Grass is material 6 (texcoords2.x); sky exposure is tangent.w.
        const float g0 = out.texcoords2[2 * i0], g1 = out.texcoords2[2 * i1], g2 = out.texcoords2[2 * i2];
        if (g0 + g1 + g2 < 1.2f) continue;
        if (out.tangents[4 * i0 + 3] + out.tangents[4 * i1 + 3] + out.tangents[4 * i2 + 3] < 2.9f) continue;
        const Vector3 a = Vert(out.positions, i0), b = Vert(out.positions, i1), c = Vert(out.positions, i2);
        const Vector3 cr = Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a));
        const float area = 0.5f * Vector3Length(cr);
        if (area < 1e-4f) continue;
        uint32_t state = Hash(nodeSeed ^ (uint32_t)t * 2654435761U);
        // Expected count = area * density; the fractional part is a coin flip.
        const float want = area * kGrassPerSquareMetre;
        int n = (int)want;
        if (Unit(state) < want - n) ++n;
        for (int k = 0; k < n; ++k) {
            float u = Unit(state), v = Unit(state);
            if (u + v > 1.0f) { u = 1.0f - u; v = 1.0f - v; }
            const float w = 1.0f - u - v;
            const Vector3 p = Vector3Add(Vector3Add(Vector3Scale(a, w), Vector3Scale(b, u)), Vector3Scale(c, v));
            if (p.x >= edge || p.y >= edge || p.z >= edge) continue;
            if (w * g0 + u * g1 + v * g2 < 0.5f) continue;
            const Vector3 local = Vector3Add(origin, p);
            const float r = Vector3Length(local);
            if (r < sea + 0.3f) continue;
            const Vector3 up = Vector3Scale(local, 1.0f / r);
            const Vector3 n0 = Vert(out.normals, i0), n1 = Vert(out.normals, i1), n2 = Vert(out.normals, i2);
            const Vector3 nrm = Vector3Normalize(Vector3Add(Vector3Add(Vector3Scale(n0, w), Vector3Scale(n1, u)), Vector3Scale(n2, v)));
            if (Vector3DotProduct(nrm, up) < 0.82f) continue; // steeper than ~35 degrees: bare
            const float h = Unit(state);
            // Variant: a quarter sparse tufts, the rest clumps (taller), so open ground is covered.
            const float pick = Unit(state), su = Unit(state);
            const int var = pick < 0.25f ? 0 : pick < 0.7f ? 1 : 2;
            const float scale = var == 0 ? 0.22f + 0.35f * su * su   // mostly short, a few tall
                              : var == 1 ? 0.35f + 0.3f * su
                                         : 0.42f + 0.35f * su;
            // The variant rides in the scale (scale + 10 * variant); grouped by variant for one draw each.
            byVariant[var].insert(byVariant[var].end(), { p.x, p.y, p.z, scale + 10.0f * var, up.x, up.y, up.z, h });
        }
    }
    for (auto& v : byVariant) out.grass.insert(out.grass.end(), v.begin(), v.end());
}

void Foliage_Upload(FoliageGpu& g, const std::vector<float>& inst) {
    Foliage_Release(g);
    if (inst.empty()) return;
    EnsureTuftMesh();
    // Instances arrive grouped by variant (ScatterGrass); count each group from the scale's tens.
    const int n = (int)(inst.size() / kGrassFloats);
    for (int i = 0; i < n; ++i) ++g.count[std::clamp((int)(inst[(size_t)i * kGrassFloats + 3] / 10.0f), 0, kGrassVariants - 1)];
    g.vbo = rlLoadVertexBuffer(inst.data(), (int)(inst.size() * sizeof(float)), false);
    const int stride = kGrassFloats * (int)sizeof(float);
    int first = 0;
    for (int v = 0; v < kGrassVariants; ++v) {
        if (!g.count[v]) continue;
        // One VAO per variant: the shared mesh, and the instance buffer from this variant's first instance.
        g.vao[v] = rlLoadVertexArray();
        rlEnableVertexArray(g.vao[v]);
        rlEnableVertexBuffer(tuftVbo);
        rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0);
        rlEnableVertexAttribute(0);
        rlEnableVertexBuffer(g.vbo);
        rlSetVertexAttribute(10, 4, RL_FLOAT, false, stride, first * stride);                       // position + scale
        rlEnableVertexAttribute(10);
        rlSetVertexAttributeDivisor(10, 1);
        rlSetVertexAttribute(11, 4, RL_FLOAT, false, stride, first * stride + 4 * (int)sizeof(float)); // up + hash
        rlEnableVertexAttribute(11);
        rlSetVertexAttributeDivisor(11, 1);
        rlDisableVertexArray();
        first += g.count[v];
    }
    rlDisableVertexBuffer();
}

void Foliage_Release(FoliageGpu& g) {
    for (unsigned int& vao : g.vao) if (vao) rlUnloadVertexArray(vao);
    if (g.vbo) rlUnloadVertexBuffer(g.vbo);
    g = FoliageGpu{};
}

void Foliage_Draw(const FoliageGpu& g, Shader shader, int modelLoc, const Matrix& model) {
    if (!g.vbo) return;
    rlSetUniformMatrix(modelLoc, model);
    (void)shader;
    for (int v = 0; v < kGrassVariants; ++v) {
        if (!g.vao[v] || !g.count[v]) continue;
        rlEnableVertexArray(g.vao[v]);
        rlDrawVertexArrayInstanced(variantFirst[v], variantCount[v], g.count[v]);
    }
    rlDisableVertexArray();
}
