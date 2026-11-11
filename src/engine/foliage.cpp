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

    // The shared tuft: five tapered blades fanned around the up axis, leaning out. y runs 0..1 from
    // root to tip (the shader scales it and bends the tips with the wind).
    unsigned int tuftVbo = 0;
    int tuftVertexCount = 0;
    void EnsureTuftMesh() {
        if (tuftVbo) return;
        std::vector<float> v;
        const int blades = 5;
        for (int b = 0; b < blades; ++b) {
            const float a = b * 2.0f * PI / blades + 0.4f * b;
            const float c = cosf(a), s = sinf(a);
            const float lean = 0.18f + 0.07f * (b % 3);
            const float w = 0.035f;
            // Root left, root right, tip (leaning outward along the blade's direction).
            const float pts[3][3] = { { -s * w, 0.0f, c * w }, { s * w, 0.0f, -c * w }, { c * lean, 1.0f - 0.1f * (b % 2), s * lean } };
            for (auto& p : pts) v.insert(v.end(), { p[0], p[1], p[2] });
        }
        tuftVertexCount = (int)(v.size() / 3);
        tuftVbo = rlLoadVertexBuffer(v.data(), (int)(v.size() * sizeof(float)), false);
    }
}

void ScatterGrass(const DensityField& field, ChunkMeshData& out) {
    if (out.key.level != 0 || out.indices.empty() || out.texcoords2.empty() || out.tangents.empty()) return;
    const float sea = field.seaLevelRadius();
    if (sea <= 0.0f) return; // no seas, no grass (moons, asteroids)
    const Vector3 origin = out.key.origin();
    const float edge = (float)CHUNK_SIZE; // keep the overlap cell to the neighbour (it draws its own)
    const uint32_t nodeSeed = Hash((uint32_t)out.key.x * 73856093U ^ (uint32_t)out.key.y * 19349663U ^ (uint32_t)out.key.z * 83492791U);
    const size_t tris = out.indices.size() / 3;
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
            const float scale = 0.18f + 0.4f * Unit(state) * Unit(state); // mostly short, a few tall
            out.grass.insert(out.grass.end(), { p.x, p.y, p.z, scale, up.x, up.y, up.z, h });
        }
    }
}

void Foliage_Upload(FoliageGpu& g, const std::vector<float>& inst) {
    Foliage_Release(g);
    if (inst.empty()) return;
    EnsureTuftMesh();
    g.count = (int)(inst.size() / kGrassFloats);
    g.vao = rlLoadVertexArray();
    rlEnableVertexArray(g.vao);
    rlEnableVertexBuffer(tuftVbo);
    rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0);
    rlEnableVertexAttribute(0);
    g.vbo = rlLoadVertexBuffer(inst.data(), (int)(inst.size() * sizeof(float)), false); // binds it
    const int stride = kGrassFloats * (int)sizeof(float);
    rlSetVertexAttribute(10, 4, RL_FLOAT, false, stride, 0);                 // position + scale
    rlEnableVertexAttribute(10);
    rlSetVertexAttributeDivisor(10, 1);
    rlSetVertexAttribute(11, 4, RL_FLOAT, false, stride, 4 * sizeof(float)); // up + hash
    rlEnableVertexAttribute(11);
    rlSetVertexAttributeDivisor(11, 1);
    rlDisableVertexArray();
    rlDisableVertexBuffer();
}

void Foliage_Release(FoliageGpu& g) {
    if (g.vao) rlUnloadVertexArray(g.vao);
    if (g.vbo) rlUnloadVertexBuffer(g.vbo);
    g = FoliageGpu{};
}

void Foliage_Draw(const FoliageGpu& g, Shader shader, int modelLoc, const Matrix& model) {
    if (!g.vao || !g.count) return;
    rlSetUniformMatrix(modelLoc, model);
    (void)shader;
    rlEnableVertexArray(g.vao);
    rlDrawVertexArrayInstanced(0, tuftVertexCount, g.count);
    rlDisableVertexArray();
}
