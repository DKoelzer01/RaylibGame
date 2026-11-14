#include "trees.h"
#include "jobsystem.h"
#include "worldtypes.h"
#include "world.h"

#include <raymath.h>
#include <rlgl.h>
#include "gpucompute.h"
#include "external/glad.h"
#include <algorithm>
#include <cmath>
#include <mutex>
#include <random>

// Bark and leaf tints multiply the shader's bark brown and latitude leaf colour.
const WoodSpeciesDef& WoodSpeciesInfo(WoodSpecies s) {
    static const WoodSpeciesDef kDefs[(int)WoodSpecies::Count] = {
        { "oak_log",    "oak",    false, { 1.0f, 1.0f, 1.0f },  { 1.0f, 1.0f, 1.0f } },
        { "birch_log",  "birch",  false, { 3.7f, 5.5f, 8.6f },  { 1.45f, 1.35f, 0.8f } },  // white bark, light yellow-green
        { "maple_log",  "maple",  false, { 1.2f, 0.85f, 0.8f }, { 4.4f, 0.55f, 0.6f } },   // red-orange crowns
        { "pine_log",   "pine",   true,  { 1.15f, 1.0f, 0.9f }, { 1.0f, 1.0f, 1.0f } },
        { "spruce_log", "spruce", true,  { 0.75f, 0.7f, 0.7f }, { 0.7f, 0.85f, 1.7f } },   // dark blue-green
    };
    return kDefs[std::clamp((int)s, 0, (int)WoodSpecies::Count - 1)];
}

namespace {

// ---------------------------------------------------------------------------------------------
// Variant geometry: connected rectangular prisms
// ---------------------------------------------------------------------------------------------

struct VariantMesh {
    std::vector<float> pos, nrm;
    std::vector<unsigned char> col;   // r: 255 leaf / 0 bark, g: shade variation
    float trunkRadius = 0.0f, trunkHeight = 0.0f;
    unsigned int vboPos = 0, vboNrm = 0, vboCol = 0;
    int vertexCount = 0;
    float boundRadius = 0.0f;         // max vertex distance from the base (m at scale 1): the GPU view cull (#292)
};
VariantMesh gVariants[2][kTreeVariants]; // [0] full, [1] far (a few prisms)
unsigned char gStoneFlag = 0; // set while building the loose-stone mesh: its vertices draw grey (colour b)
std::once_flag gVariantsOnce;

void Emit(VariantMesh& m, Vector3 p, Vector3 n, bool leaf, unsigned char shade) {
    m.pos.insert(m.pos.end(), { p.x, p.y, p.z });
    m.nrm.insert(m.nrm.end(), { n.x, n.y, n.z });
    m.col.insert(m.col.end(), { (unsigned char)(leaf ? 255 : 0), shade, gStoneFlag, 255 });
}

void Quad(VariantMesh& m, Vector3 a, Vector3 b, Vector3 c, Vector3 d, Vector3 centre, bool leaf, unsigned char shade) {
    Vector3 n = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(d, a)));
    const Vector3 mid = Vector3Scale(Vector3Add(Vector3Add(a, b), Vector3Add(c, d)), 0.25f);
    if (Vector3DotProduct(n, Vector3Subtract(mid, centre)) < 0.0f) { n = Vector3Negate(n); std::swap(b, d); } // outward, CCW
    Emit(m, a, n, leaf, shade); Emit(m, b, n, leaf, shade); Emit(m, c, n, leaf, shade);
    Emit(m, a, n, leaf, shade); Emit(m, c, n, leaf, shade); Emit(m, d, n, leaf, shade);
}

// A prism from a to b with a square cross-section w0 wide at a, w1 at b, turned by `twist`.
void Prism(VariantMesh& m, Vector3 a, Vector3 b, float w0, float w1, float twist, bool leaf, unsigned char shade) {
    const Vector3 d = Vector3Normalize(Vector3Subtract(b, a));
    const Vector3 ref = fabsf(d.y) < 0.95f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 };
    Vector3 s = Vector3Normalize(Vector3CrossProduct(d, ref));
    Vector3 t = Vector3CrossProduct(d, s);
    const float c = cosf(twist), sn = sinf(twist);
    const Vector3 s2 = Vector3Add(Vector3Scale(s, c), Vector3Scale(t, sn));
    const Vector3 t2 = Vector3Subtract(Vector3Scale(t, c), Vector3Scale(s, sn));
    s = s2; t = t2;
    const float sx[4] = { 1, -1, -1, 1 }, ty[4] = { 1, 1, -1, -1 };
    Vector3 A[4], B[4];
    for (int k = 0; k < 4; ++k) {
        A[k] = Vector3Add(a, Vector3Add(Vector3Scale(s, sx[k] * w0 * 0.5f), Vector3Scale(t, ty[k] * w0 * 0.5f)));
        B[k] = Vector3Add(b, Vector3Add(Vector3Scale(s, sx[k] * w1 * 0.5f), Vector3Scale(t, ty[k] * w1 * 0.5f)));
    }
    const Vector3 centre = Vector3Scale(Vector3Add(a, b), 0.5f);
    for (int k = 0; k < 4; ++k) Quad(m, A[k], A[(k + 1) % 4], B[(k + 1) % 4], B[k], centre, leaf, shade);
    Quad(m, B[0], B[1], B[2], B[3], centre, leaf, shade); // top
    Quad(m, A[0], A[1], A[2], A[3], centre, leaf, shade); // bottom
}

float U(std::mt19937& r) { return std::uniform_real_distribution<float>(0.0f, 1.0f)(r); }
unsigned char Shade(std::mt19937& r) { return (unsigned char)(90 + (int)(U(r) * 140.0f)); }

// Shape knobs per species: height, crown size, trunk thickness, branch count.
struct BroadleafShape { float h = 1.0f, crown = 1.0f, trunk = 1.0f; int branches = 0; };
struct ConiferShape { float h = 1.0f, width = 1.0f; int layers = 0; };

void Broadleaf(VariantMesh& m, unsigned int seed, BroadleafShape sh = {}) {
    std::mt19937 r(seed);
    const float H = (4.2f + 1.6f * U(r)) * sh.h;
    const Vector3 lean = { (U(r) - 0.5f) * 0.6f, 0.0f, (U(r) - 0.5f) * 0.6f };
    auto trunkAt = [&](float y) { const float f = y / H; return Vector3{ lean.x * f * f, y, lean.z * f * f }; };
    // Trunk: three stacked tapering segments, slightly bent.
    const float ys[4] = { -0.4f, H * 0.35f, H * 0.7f, H };
    const float ws[4] = { 0.62f * sh.trunk, 0.45f * sh.trunk, 0.34f * sh.trunk, 0.24f * sh.trunk };
    for (int s = 0; s < 3; ++s) Prism(m, trunkAt(ys[s]), trunkAt(ys[s + 1]), ws[s], ws[s + 1], 0.2f * s, false, Shade(r));
    // Branches up and out, each ending in a leaf block.
    const int nb = 3 + (int)(U(r) * 3.0f) + sh.branches;
    const float yaw0 = U(r) * 2.0f * PI;
    for (int k = 0; k < nb; ++k) {
        const float y = H * (0.5f + 0.4f * U(r));
        const float yaw = yaw0 + k * 2.0f * PI / nb + (U(r) - 0.5f) * 0.6f;
        const float el = (35.0f + 25.0f * U(r)) * DEG2RAD;
        const float len = 1.5f + 1.3f * U(r);
        const Vector3 s = trunkAt(y);
        const Vector3 dir = { cosf(el) * cosf(yaw), sinf(el), cosf(el) * sinf(yaw) };
        const Vector3 e = Vector3Add(s, Vector3Scale(dir, len));
        Prism(m, s, e, 0.2f, 0.12f, U(r), false, Shade(r));
        const float c = (1.7f + 1.1f * U(r)) * sh.crown;
        Prism(m, Vector3Subtract(e, { 0, c * 0.4f, 0 }), Vector3Add(e, { 0, c * 0.6f, 0 }), c, c * 0.85f, U(r) * 1.5f, true, Shade(r));
    }
    // Crown.
    const float c = (2.5f + 0.9f * U(r)) * sh.crown;
    const Vector3 top = trunkAt(H);
    Prism(m, Vector3Subtract(top, { 0, 0.9f, 0 }), Vector3Add(top, { 0, c - 0.6f, 0 }), c, c * 0.8f, U(r) * 1.5f, true, Shade(r));
    m.trunkRadius = 0.3f * sh.trunk;
    m.trunkHeight = H;
}

void Conifer(VariantMesh& m, unsigned int seed, ConiferShape sh = {}) {
    std::mt19937 r(seed);
    const float H = (8.5f + 2.0f * U(r)) * sh.h;
    Prism(m, { 0, -0.4f, 0 }, { 0, H, 0 }, 0.5f, 0.12f, 0.0f, false, Shade(r));
    // Stacked leaf slabs, wide at the bottom, turned alternately.
    const int n = 6 + (int)(U(r) * 2.0f) + sh.layers;
    for (int k = 0; k < n; ++k) {
        const float f = (float)k / (n - 1);
        const float y = 1.6f + f * (H - 2.2f);
        const float w = (3.8f * (1.0f - f) + 0.8f) * sh.width;
        Prism(m, { 0, y, 0 }, { 0, y + 1.5f, 0 }, w, w * 0.5f, k * 0.785f + U(r) * 0.3f, true, Shade(r));
    }
    Prism(m, { 0, H - 0.2f, 0 }, { 0, H + 0.9f, 0 }, 0.7f, 0.15f, 0.4f, true, Shade(r));
    m.trunkRadius = 0.25f;
    m.trunkHeight = H;
}

void Bush(VariantMesh& m, unsigned int seed) {
    std::mt19937 r(seed);
    const int n = 2 + (int)(U(r) * 3.0f);
    for (int k = 0; k < n; ++k) {
        const float c = 0.7f + 0.6f * U(r);
        const Vector3 o = { (U(r) - 0.5f) * 0.9f, 0.15f + 0.4f * U(r), (U(r) - 0.5f) * 0.9f };
        Prism(m, Vector3Subtract(o, { 0, c * 0.45f, 0 }), Vector3Add(o, { 0, c * 0.55f, 0 }), c, c * 0.8f, U(r) * 1.5f, true, Shade(r));
    }
    m.trunkRadius = 0.0f;
    m.trunkHeight = 0.0f;
}

// Loose stones: two to four low, tilted lumps a hand or two across, lying on the ground.
void Stones(VariantMesh& m, unsigned int seed) {
    std::mt19937 r(seed);
    gStoneFlag = 255;
    const int n = 2 + (int)(U(r) * 3.0f);
    for (int k = 0; k < n; ++k) {
        const float c = 0.22f + 0.28f * U(r);
        const Vector3 o = { (U(r) - 0.5f) * 0.8f, 0.0f, (U(r) - 0.5f) * 0.8f };
        const Vector3 tilt = { (U(r) - 0.5f) * 0.3f, 0.0f, (U(r) - 0.5f) * 0.3f };
        Prism(m, Vector3Subtract(o, { 0, c * 0.25f, 0 }), Vector3Add(Vector3Add(o, tilt), { 0, c * 0.45f, 0 }), c, c * 0.7f, U(r) * 1.5f, false, Shade(r));
    }
    gStoneFlag = 0;
    m.trunkRadius = 0.0f;
    m.trunkHeight = 0.0f;
}

// Far meshes: the same silhouette from a few prisms (trunk + crown; trunk + two cones; one block).
void FarFrom(const VariantMesh& full, int v, VariantMesh& m) {
    float lo[3] = { 1e9f, 1e9f, 1e9f }, hi[3] = { -1e9f, -1e9f, -1e9f };
    for (size_t k = 0; k + 2 < full.pos.size(); k += 3)
        if (full.col[(k / 3) * 4] > 128) // leaf vertices: the crown's extent
            for (int a = 0; a < 3; ++a) { lo[a] = std::min(lo[a], full.pos[k + a]); hi[a] = std::max(hi[a], full.pos[k + a]); }
    const float w = std::max(hi[0] - lo[0], hi[2] - lo[2]) * 0.8f;
    const Vector3 c = { (lo[0] + hi[0]) * 0.5f, 0.0f, (lo[2] + hi[2]) * 0.5f };
    const unsigned char shade = 160;
    if (full.trunkHeight > 0.0f) Prism(m, { 0, -0.4f, 0 }, { c.x, lo[1] + 0.3f, c.z }, full.trunkRadius * 2.0f, full.trunkRadius * 1.4f, 0.0f, false, shade);
    if (v < kFirstBush && WoodSpeciesInfo(SpeciesOfVariant(v)).conifer) {
        const float mid = lo[1] + (hi[1] - lo[1]) * 0.55f;
        Prism(m, { c.x, lo[1], c.z }, { c.x, mid, c.z }, w, w * 0.6f, 0.785f, true, shade);
        Prism(m, { c.x, mid, c.z }, { c.x, hi[1], c.z }, w * 0.6f, 0.15f, 0.0f, true, shade);
    } else {
        Prism(m, { c.x, lo[1], c.z }, { c.x, hi[1], c.z }, w, w * 0.8f, 0.4f, true, shade);
    }
    m.trunkRadius = full.trunkRadius;
    m.trunkHeight = full.trunkHeight;
}

void BuildVariants() {
    for (int v = 0; v < kTreeVariants; ++v) {
        VariantMesh& m = gVariants[0][v];
        if (v < kFirstBush) {
            // Two shapes per species (oak broad and heavy, birch tall and slender, maple round and full;
            // pine tall and sparse, spruce dense and broad).
            static const BroadleafShape kBroad[3] = { { 1.0f, 1.25f, 1.3f, 2 }, { 1.35f, 0.65f, 0.7f, 0 }, { 1.05f, 1.1f, 1.0f, 1 } };
            static const ConiferShape kCone[2] = { { 1.2f, 0.75f, -1 }, { 0.9f, 1.1f, 2 } };
            const int sp = (int)SpeciesOfVariant(v);
            if (WoodSpeciesInfo((WoodSpecies)sp).conifer) Conifer(m, 2000u + v * 7919u, kCone[sp - (int)WoodSpecies::Pine]);
            else Broadleaf(m, 1000u + v * 7919u, kBroad[sp]);
        }
        else if (v == kStoneVariant) Stones(m, 4000u + v * 7919u);
        else Bush(m, 3000u + v * 7919u);
        m.vertexCount = (int)(m.pos.size() / 3);
        if (v == kStoneVariant) gVariants[1][v] = m; // already a few prisms (and no crown for FarFrom to fit)
        else FarFrom(m, v, gVariants[1][v]);
        gVariants[1][v].vertexCount = (int)(gVariants[1][v].pos.size() / 3);
    }
}
const VariantMesh& Variant(int v, int lod = 0) { std::call_once(gVariantsOnce, BuildVariants); return gVariants[lod][v]; }

void EnsureVariantGpu() {
    for (int lod = 0; lod < 2; ++lod)
        for (int v = 0; v < kTreeVariants; ++v) {
            VariantMesh& m = const_cast<VariantMesh&>(Variant(v, lod));
            if (m.vboPos) continue;
            m.vboPos = rlLoadVertexBuffer(m.pos.data(), (int)(m.pos.size() * sizeof(float)), false);
            m.vboNrm = rlLoadVertexBuffer(m.nrm.data(), (int)(m.nrm.size() * sizeof(float)), false);
            m.vboCol = rlLoadVertexBuffer(m.col.data(), (int)m.col.size(), false);
            for (size_t k = 0; k + 2 < m.pos.size(); k += 3) m.boundRadius = std::max(m.boundRadius, sqrtf(m.pos[k] * m.pos[k] + m.pos[k + 1] * m.pos[k + 1] + m.pos[k + 2] * m.pos[k + 2]));
        }
}

// ---------------------------------------------------------------------------------------------
// Cube-sphere tiles
// ---------------------------------------------------------------------------------------------

const Vector3 kN[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
const Vector3 kU[6] = { { 0, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 0 } };
const Vector3 kV[6] = { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 } };
constexpr int kCells = 16; // 4 m candidate cells per 64 m tile

int TilesPerFace(float radius) { return std::max(1, (int)ceilf(0.5f * PI * radius / kTreeTileSize)); }

Vector3 FaceDir(int face, float alpha, float beta) {
    return Vector3Normalize(Vector3Add(kN[face], Vector3Add(Vector3Scale(kU[face], tanf(alpha)), Vector3Scale(kV[face], tanf(beta)))));
}

uint32_t Hash(uint32_t x) {
    x ^= x >> 16; x *= 0x7feb352dU;
    x ^= x >> 15; x *= 0x846ca68bU;
    x ^= x >> 16;
    return x;
}
float Unit(uint32_t& s) { s = Hash(s + 0x9e3779b9U); return (s >> 8) * (1.0f / 16777216.0f); }

// Smooth value noise on a 3D lattice (forest patches), two octaves, in [0, 1].
float Lattice(int x, int y, int z, uint32_t seed) {
    return (Hash((uint32_t)x * 73856093U ^ (uint32_t)y * 19349663U ^ (uint32_t)z * 83492791U ^ seed) >> 8) * (1.0f / 16777216.0f);
}
float ValueNoise(Vector3 p, uint32_t seed) {
    const int x0 = (int)floorf(p.x), y0 = (int)floorf(p.y), z0 = (int)floorf(p.z);
    float fx = p.x - x0, fy = p.y - y0, fz = p.z - z0;
    fx = fx * fx * (3 - 2 * fx); fy = fy * fy * (3 - 2 * fy); fz = fz * fz * (3 - 2 * fz);
    auto L = [&](int dx, int dy, int dz) { return Lattice(x0 + dx, y0 + dy, z0 + dz, seed); };
    const float x00 = L(0, 0, 0) + (L(1, 0, 0) - L(0, 0, 0)) * fx, x10 = L(0, 1, 0) + (L(1, 1, 0) - L(0, 1, 0)) * fx;
    const float x01 = L(0, 0, 1) + (L(1, 0, 1) - L(0, 0, 1)) * fx, x11 = L(0, 1, 1) + (L(1, 1, 1) - L(0, 1, 1)) * fx;
    const float y0v = x00 + (x10 - x00) * fy, y1v = x01 + (x11 - x01) * fy;
    return y0v + (y1v - y0v) * fz;
}
float Forest(Vector3 local, uint32_t seed) {
    return 0.65f * ValueNoise(Vector3Scale(local, 1.0f / 700.0f), seed) + 0.35f * ValueNoise(Vector3Scale(local, 1.0f / 230.0f), seed ^ 0x51ed27u);
}
float Smooth(float a, float b, float x) { const float t = std::clamp((x - a) / (b - a), 0.0f, 1.0f); return t * t * (3 - 2 * t); }

// What grows in each biome (Plan 5 B): tree, bush and loose-stone density (multipliers), conifer or
// broadleaf woods, and tree size. Jungles and taiga are dense, savanna and meadows open with lone
// trees, deserts and badlands nearly bare (scrub and stones), tundra a few stunted conifers.
struct BiomeFlora { float trees, bushes, stones; bool cold; float size; float cover = 0.0f; }; // cover: least tree odds anywhere in it (closed canopy)
static BiomeFlora FloraOf(Biome b) {
    switch (b) {
        case Biome::Jungle:   return { 1.5f, 1.4f, 0.5f, false, 1.35f, 0.34f };
        case Biome::Savanna:  return { 0.08f, 0.45f, 1.0f, false, 1.15f };
        case Biome::Desert:   return { 0.0f, 0.06f, 2.0f, false, 1.0f };
        case Biome::Badlands: return { 0.02f, 0.12f, 2.5f, false, 0.8f };
        case Biome::Forest:   return { 1.0f, 1.0f, 1.0f, false, 1.0f };
        case Biome::Meadow:   return { 0.07f, 1.3f, 0.8f, false, 1.1f };
        case Biome::Swamp:    return { 0.45f, 1.5f, 0.3f, false, 0.9f };
        case Biome::Taiga:    return { 1.35f, 0.6f, 1.0f, true, 1.05f, 0.26f };
        case Biome::Tundra:   return { 0.04f, 0.35f, 1.6f, true, 0.6f };
        case Biome::Alpine:   return { 0.25f, 0.5f, 1.5f, true, 0.8f };
        default:              return { 0.0f, 0.0f, 0.0f, true, 1.0f };   // ocean, polar
    }
}
constexpr float kMaxTreeFactor = 1.5f, kMaxBushFactor = 1.5f, kMaxStoneFactor = 2.5f;

// Chance of a tree in a 4 m cell from the forest field, scaled by the biome's tree factor.
float TreeChanceFromForest(float forest, float biomeTrees) {
    return (0.55f * Smooth(0.48f, 0.72f, forest) + 0.015f) * biomeTrees;
}

uint64_t PlantId(const TreeTileKey& k, int cell) {
    return ((uint64_t)k.face << 56) | ((uint64_t)(k.i & 0xFFFFFF) << 32) | ((uint64_t)(k.j & 0xFFFFFF) << 8) | (uint64_t)(cell & 0xFF);
}
TreeTileKey TileOfPlant(uint64_t id) { return { (int)(id >> 56), (int)((id >> 32) & 0xFFFFFF), (int)((id >> 8) & 0xFFFFFF) }; }

void BuildTile(const PlanetField& f, const TerrainDeltas& edits, const PlantRemovals& removed, const std::vector<Vector4>& clearings, double now, uint32_t seed, const TreeTileKey& key, TreeTile& out) {
    const float R = f.radius;
    const int N = TilesPerFace(R);
    const float da = 0.5f * PI / N, dc = da / kCells;
    const float a0 = -0.25f * PI + key.i * da, b0 = -0.25f * PI + key.j * da;
    {
        const Vector3 cd = FaceDir(key.face, a0 + 0.5f * da, b0 + 0.5f * da);
        out.centre = Vector3Scale(cd, R + std::max(0.0f, f.elevation(cd, 16.0f)));
    }
    auto solid = [&](Vector3 p) { return f.density(p, &edits) >= ISO_LEVEL; };
    for (int ci = 0; ci < kCells; ++ci)
        for (int cj = 0; cj < kCells; ++cj) {
            uint32_t s = Hash(seed ^ Hash((uint32_t)key.face * 2654435761U ^ (uint32_t)(key.i * kCells + ci) * 40503U ^ (uint32_t)(key.j * kCells + cj) * 9973U));
            const float jx = Unit(s), jy = Unit(s), roll = Unit(s), h = Unit(s), sz = Unit(s);
            const Vector3 dir = FaceDir(key.face, a0 + (ci + jx) * dc, b0 + (cj + jy) * dc);
            const float forest = Forest(Vector3Scale(dir, R), seed);
            // First a cheap cut at the densest any biome allows, then the biome's own odds (it needs the
            // height, which is the expensive part).
            const float pTreeMax = std::max(TreeChanceFromForest(forest, kMaxTreeFactor), 0.34f);
            const float pBushBase = 0.05f + 0.3f * Smooth(0.35f, 0.6f, forest);
            const float pStoneBase = 0.025f;                              // loose stones: the first stone tools' stone
            if (roll >= pTreeMax + pBushBase * kMaxBushFactor + pStoneBase * kMaxStoneFactor) continue;
            const float lat = fabsf(dir.y);
            if (lat > 0.86f) continue;                                  // polar: nothing grows
            // Floating islands (Plan 5 C) have woods on their tops: the column's ground is the island's top.
            float islandTop = 0.0f;
            const bool onIsland = f.islandTopAt(dir, islandTop);
            const float e = onIsland ? islandTop - R : f.elevation(dir);
            const Biome biome = onIsland ? Biome::Forest : f.biomeAt(dir, e);
            const BiomeFlora flora = FloraOf(biome);
            const float pTree = std::max(TreeChanceFromForest(forest, flora.trees), flora.cover), pBush = pBushBase * flora.bushes, pStone = pStoneBase * flora.stones;
            const bool tree = roll < pTree;
            const bool stone = !tree && roll >= pTree + pBush && roll < pTree + pBush + pStone;
            if (!tree && !stone && roll >= pTree + pBush) continue;
            if (!onIsland && (e < 4.0f || e > (tree ? 900.0f : stone ? 1300.0f : 1100.0f))) continue;   // beach / sea, above the tree line
            Vector3 normal = dir;
            Vector3 ground = Vector3Scale(dir, R + e);
            if (onIsland) {
                // The island's real top (hills, terraces): the first solid going down; no trees on its cliffs.
                bool found = false;
                for (float hh = 6.0f; hh > -10.0f && !found; hh -= 0.25f)
                    if (solid(Vector3Scale(dir, R + e + hh))) { ground = Vector3Scale(dir, R + e + hh + 0.1f); found = true; }
                if (!found) continue;
            } else {
                // Slope from the height a few metres away.
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const Vector3 t2 = Vector3CrossProduct(t1, dir);
                const float step = 3.0f, ang = step / (R + e);
                const float e1 = f.elevation(Vector3Normalize(Vector3Add(dir, Vector3Scale(t1, ang))));
                const float e2 = f.elevation(Vector3Normalize(Vector3Add(dir, Vector3Scale(t2, ang))));
                normal = Vector3Normalize(Vector3Subtract(dir, Vector3Add(Vector3Scale(t1, (e1 - e) / step), Vector3Scale(t2, (e2 - e) / step))));
                const float slope = 1.0f - Vector3DotProduct(normal, dir);
                if (slope > (tree ? 0.12f : stone ? 0.3f : 0.2f)) continue;
                // The drawn surface is the density field (3D band, crags, tors), up to metres off the height map
                // (#34: bushes and stones floated or sank). Walk down from 3 m above it to the first solid sample,
                // at the footprint the finest chunks mesh with, and interpolate the crossing.
                auto dens = [&](float hh) { return f.density(Vector3Scale(dir, R + e + hh), &edits, 1.0f); };
                float prevH = 3.0f, prevD = dens(prevH);
                if (prevD >= ISO_LEVEL) continue;                    // buried under an overhang or a tor
                bool found = false;
                for (float hh = 2.75f; hh >= -4.0f; hh -= 0.25f) {
                    const float d = dens(hh);
                    if (d >= ISO_LEVEL) {
                        const float t = (ISO_LEVEL - prevD) / std::max(d - prevD, 1e-6f); // air above (prevD < ISO), solid here
                        ground = Vector3Scale(dir, R + e + prevH + (hh - prevH) * t);
                        found = true;
                        break;
                    }
                    prevH = hh; prevD = d;
                }
                if (!found) continue;
            }
            const VoxelMaterial m = f.materialLod(Vector3Subtract(ground, Vector3Scale(dir, 0.3f)), 1.0f, 0.0f, normal);
            if (m != VoxelMaterial::Grass && m != VoxelMaterial::Regolith && !(stone && m == VoxelMaterial::Gravel)) continue;
            // Solid ground underfoot and open air above (caves, ravines, edits).
            if (!solid(Vector3Subtract(ground, Vector3Scale(dir, 1.0f))) || !solid(Vector3Subtract(ground, Vector3Scale(dir, 4.0f)))) continue;
            if (solid(Vector3Add(ground, Vector3Scale(dir, 1.5f))) || (tree && solid(Vector3Add(ground, Vector3Scale(dir, 6.0f))))) continue;
            int v;
            float scale;
            if (tree) {
                // Broadleaf woods (oak, birch, maple) in the temperate band; conifer woods (pine, spruce,
                // birch) from ~40 degrees and up high, mixing at the edge. Each area (a ~400 m noise
                // patch) has a dominant species: 60% of its trees, the other two 20% each.
                const bool cold = flora.cold || e > 550.0f + 300.0f * (h - 0.5f); // conifers in cold biomes and up high (Plan 5 B)
                static const WoodSpecies kWarm[3] = { WoodSpecies::Oak, WoodSpecies::Maple, WoodSpecies::Birch };
                static const WoodSpecies kCold[3] = { WoodSpecies::Pine, WoodSpecies::Spruce, WoodSpecies::Birch };
                const float area = ValueNoise(Vector3Scale(Vector3Scale(dir, R), 1.0f / 400.0f), seed ^ 0x7a3c19u);
                const int dom = std::min(2, (int)(Smooth(0.3f, 0.7f, area) * 3.0f));
                const float pick = Unit(s);
                const int k = pick < 0.6f ? dom : pick < 0.8f ? (dom + 1) % 3 : (dom + 2) % 3;
                const WoodSpecies sp = (cold ? kCold : kWarm)[k];
                v = (int)sp * 2 + (int)(h * 2.0f) % 2;
                scale = (0.6f + 0.7f * sz + 0.35f * forest) * flora.size;
                if (sp == WoodSpecies::Birch) scale *= 0.85f; // slender
            } else if (stone) {
                v = kStoneVariant;
                scale = 0.8f + 0.6f * sz;
            } else {
                v = kFirstBush + (int)(h * 3.0f) % 3;
                scale = 0.7f + 0.7f * sz;
            }
            bool cleared = false;
            for (const Vector4& c : clearings) if (Vector3Distance(ground, { c.x, c.y, c.z }) < c.w) { cleared = true; break; }
            if (cleared) continue;
            const uint64_t id = PlantId(key, ci * kCells + cj);
            out.builtAt = now;
            if (auto it = removed.find(id); it != removed.end()) {
                // Chopped / harvested: gone for kRegrowDelay, then a sapling growing to full size.
                const double age = now - it->second - kRegrowDelay;
                if (age < 0.0) { out.regrowing = true; continue; }
                const float growth = (float)std::clamp(age / kRegrowTime, 0.15, 1.0);
                if (growth < 1.0f) out.regrowing = true;
                scale *= growth;
            }
            const Vector3 p = Vector3Subtract(ground, Vector3Scale(dir, 0.15f)); // sink the base a little
            // The climate (Plan 5 B) in the up vector's length, packed as for the terrain (trees.vs decodes it).
            const PlanetField::Climate cl = f.climateAt(dir, e);
            const float code = 2.0f + floorf(std::clamp((cl.temp + 30.0f) / 70.0f, 0.0f, 1.0f) * 62.0f) + std::min(cl.moist, 0.999f); // from 2: below 1.5 means no climate (#254)
            out.inst[v].insert(out.inst[v].end(), { p.x, p.y, p.z, scale, dir.x * code, dir.y * code, dir.z * code, h });
            out.ids[v].push_back(id);
            const VariantMesh& vm = Variant(v);
            if (vm.trunkRadius > 0.0f) {
                out.trunks.push_back({ p, dir, vm.trunkHeight * scale, vm.trunkRadius * scale });
                out.plants.push_back({ id, p, dir, vm.trunkHeight * scale, vm.trunkRadius * scale, scale, v });
            } else {
                out.plants.push_back({ id, p, dir, (stone ? 0.5f : 1.2f) * scale, 0.7f * scale, scale, v });
            }
        }
}

} // namespace

TreeTileKey TreeTileOf(Vector3 dir, float radius) {
    int face = 0;
    float best = -1.0f;
    for (int f = 0; f < 6; ++f) { const float d = Vector3DotProduct(dir, kN[f]); if (d > best) { best = d; face = f; } }
    const int N = TilesPerFace(radius);
    const float da = 0.5f * PI / N;
    const float a = atanf(Vector3DotProduct(dir, kU[face]) / best), b = atanf(Vector3DotProduct(dir, kV[face]) / best);
    return { face, std::clamp((int)floorf((a + 0.25f * PI) / da), 0, N - 1), std::clamp((int)floorf((b + 0.25f * PI) / da), 0, N - 1) };
}

Vector3 TreeTileCentreDir(const TreeTileKey& k, float radius) {
    const float da = 0.5f * PI / TilesPerFace(radius);
    return FaceDir(k.face, -0.25f * PI + (k.i + 0.5f) * da, -0.25f * PI + (k.j + 0.5f) * da);
}

// ---------------------------------------------------------------------------------------------
// TreeSystem
// ---------------------------------------------------------------------------------------------

struct TreeSystem::Shared {
    std::mutex mutex;
    std::vector<std::pair<TreeTileKey, TreeTile>> results;
    std::atomic<bool> alive{ true };
};

TreeSystem::TreeSystem(std::shared_ptr<const PlanetField> f, unsigned int s)
    : shared(std::make_shared<Shared>()), field(std::move(f)), seed(s), radius(field->radius), edits(std::make_shared<TerrainDeltas>()),
      removed(std::make_shared<PlantRemovals>()), clearings(std::make_shared<std::vector<Vector4>>()) {}

TreeSystem::~TreeSystem() {
    shared->alive = false;
    for (int s = 0; s < kTreeSets; ++s)
        for (int slot = 0; slot < 2 * kTreeVariants; ++slot) {
            if (vao[s][slot]) rlUnloadVertexArray(vao[s][slot]);
            if (vbo[s][slot]) rlUnloadVertexBuffer(vbo[s][slot]);
        }
    for (int slot = 0; slot < 2 * kTreeVariants; ++slot) {
        if (culledVao[slot]) rlUnloadVertexArray(culledVao[slot]);
        if (culledVbo[slot]) rlUnloadVertexBuffer(culledVbo[slot]);
    }
    if (cmdBuf) glDeleteBuffers(1, &cmdBuf);
}

void TreeSystem::request(const TreeTileKey& key, TreeTile& tile, float priority) {
    const uint32_t version = ++tile.version;
    tile.pending = true;
    auto sh = shared;
    auto fld = field;
    auto ed = edits;
    auto rm = removed;
    auto cl = clearings;
    const double now = WorldTime();
    const uint32_t sd = seed;
    WorldJobs().submit(priority, [sh, fld, ed, rm, cl, now, sd, key, version] {
        if (!sh->alive) return;
        TreeTile t;
        t.version = version;
        t.ready = true;
        BuildTile(*fld, *ed, *rm, *cl, now, sd, key, t);
        std::lock_guard<std::mutex> lock(sh->mutex);
        sh->results.emplace_back(key, std::move(t));
    });
}

void TreeSystem::onEdit(Vector3 c, float r, std::shared_ptr<const TerrainDeltas> snapshot) {
    edits = std::move(snapshot);
    for (auto& kv : tiles)
        if (Vector3Distance(kv.second.centre.x != 0.0f || kv.second.centre.y != 0.0f || kv.second.centre.z != 0.0f
                                ? kv.second.centre : Vector3Scale(TreeTileCentreDir(kv.first, radius), radius), c) < r + kTreeTileSize)
            request(kv.first, kv.second, 0.0f); // the old trees stay until the rebuild arrives
}

void TreeSystem::update(Planetoid& body, Vector3 camLocal) {
    if (body.edits != edits) {
        // Edits replaced without onEdit (a load): rebuild every tile against them.
        edits = body.edits;
        for (auto& kv : tiles) request(kv.first, kv.second, 0.0f);
    }
    {
        std::lock_guard<std::mutex> lock(shared->mutex);
        for (auto& r : shared->results) {
            auto it = tiles.find(r.first);
            if (it == tiles.end() || it->second.version != r.second.version) continue;
            it->second = std::move(r.second);
            batchDirty = true;
        }
        shared->results.clear();
    }
    // Saplings grow: rebuild their tiles every 20 s of world time.
    for (auto& kv : tiles)
        if (kv.second.ready && kv.second.regrowing && !kv.second.pending && WorldTime() - kv.second.builtAt > 20.0) request(kv.first, kv.second, 500.0f);
    // Tiles near the camera (arc distance on the sphere), nothing from high above.
    const float camR = Vector3Length(camLocal);
    const Vector3 camDir = Vector3Scale(camLocal, 1.0f / std::max(camR, 1.0f));
    const float altitude = camR - radius;
    const float range = kTreeDrawDistance + kTreeTileSize;
    std::vector<TreeTileKey> wanted;
    const bool recompute = Vector3Distance(camLocal, lastWantCam) > 16.0f;
    if (recompute) lastWantCam = camLocal;
    if (recompute && altitude < kTreeDrawDistance + 2000.0f) {
        const int N = TilesPerFace(radius);
        const float da = 0.5f * PI / N;
        const int K = (int)ceilf(range / kTreeTileSize) + 1;
        for (int f = 0; f < 6; ++f) {
            const float dn = Vector3DotProduct(camDir, kN[f]);
            if (dn < 0.2f) continue;
            const float a = std::clamp(atanf(Vector3DotProduct(camDir, kU[f]) / dn), -0.25f * PI, 0.25f * PI);
            const float b = std::clamp(atanf(Vector3DotProduct(camDir, kV[f]) / dn), -0.25f * PI, 0.25f * PI);
            const int ic = (int)floorf((a + 0.25f * PI) / da), jc = (int)floorf((b + 0.25f * PI) / da);
            for (int i = std::max(0, ic - K); i <= std::min(N - 1, ic + K); ++i)
                for (int j = std::max(0, jc - K); j <= std::min(N - 1, jc + K); ++j) {
                    const TreeTileKey key{ f, i, j };
                    const float arc = acosf(std::clamp(Vector3DotProduct(TreeTileCentreDir(key, radius), camDir), -1.0f, 1.0f)) * radius;
                    if (arc > range) continue;
                    wanted.push_back(key);
                    auto it = tiles.find(key);
                    if (it == tiles.end()) {
                        TreeTile& t = tiles[key];
                        request(key, t, arc);
                    }
                }
        }
    }
    // Drop tiles well outside the range (hysteresis so the edge doesn't churn).
    for (auto it = tiles.begin(); recompute && it != tiles.end();) {
        const float arc = acosf(std::clamp(Vector3DotProduct(TreeTileCentreDir(it->first, radius), camDir), -1.0f, 1.0f)) * radius;
        if (arc > range * 1.25f || altitude > kTreeDrawDistance + 2500.0f) { it = tiles.erase(it); batchDirty = true; }
        else ++it;
    }
    const double now = GetTime();
    const bool moved = Vector3Distance(camLocal, Vector3{ (float)anchor.x, (float)anchor.y, (float)anchor.z }) > 40.0f; // keep the near/far split current
    if ((batchDirty || moved) && now - lastBatch > 0.2) rebuildBatches(camLocal);
}

void TreeSystem::rebuildBatches(Vector3 camLocal) {
    lastBatch = GetTime();
    batchDirty = false;
    anchor = DVec3(camLocal);
    EnsureVariantGpu();
    drawnInstances = 0;
    // Set 0: everything drawn (to the draw distance). Sets 1 and 2: the shadow casters for the far (+-200 m) and
    // near (+-32 m) cascades, only the trees within reach of their boxes (with room for long low-sun shadows and
    // the 40 m the camera may move before the next rebuild): the shadow passes drew the whole 2 km forest twice,
    // ~1.3 ms a frame flying, for boxes a few hundred metres across.
    const float setReach[kTreeSets] = { 1e30f, kTreeShadowFarReach, kTreeShadowNearReach };
    std::vector<float> buf[kTreeSets][2];
    for (int v = 0; v < kTreeVariants; ++v) {
        const float reach = v < kFirstBush ? kTreeDrawDistance : kBushDrawDistance;
        for (auto& s : buf) { s[0].clear(); s[1].clear(); }
        for (auto& kv : tiles) {
            const TreeTile& t = kv.second;
            if (!t.ready || t.inst[v].empty() || Vector3Distance(t.centre, camLocal) > reach + kTreeTileSize) continue;
            for (size_t k = 0; k + 7 < t.inst[v].size(); k += 8) {
                const float x = (float)(t.inst[v][k] - anchor.x), y = (float)(t.inst[v][k + 1] - anchor.y), z = (float)(t.inst[v][k + 2] - anchor.z);
                const float d2 = x * x + y * y + z * z;
                if (d2 > reach * reach) continue;
                const int lod = d2 < kTreeNearDistance * kTreeNearDistance ? 0 : 1;
                const float inst[8] = { x, y, z, t.inst[v][k + 3], t.inst[v][k + 4], t.inst[v][k + 5], t.inst[v][k + 6], t.inst[v][k + 7] };
                // Casters: near the sun axis through the camera (the cascades' boxes are centred there and run along it);
                // without a sun direction yet, by distance from the camera.
                const float along = x * shadowSunLocal.x + y * shadowSunLocal.y + z * shadowSunLocal.z;
                const bool sunKnown = shadowSunLocal.x != 0.0f || shadowSunLocal.y != 0.0f || shadowSunLocal.z != 0.0f;
                const float perp2 = sunKnown ? std::max(0.0f, d2 - along * along) : d2;
                for (int s = 0; s < kTreeSets; ++s)
                    if (s == 0 || perp2 <= setReach[s] * setReach[s]) buf[s][lod].insert(buf[s][lod].end(), inst, inst + 8);
            }
        }
        for (int s = 0; s < kTreeSets; ++s)
            for (int lod = 0; lod < 2; ++lod) {
                const int slot = v + lod * kTreeVariants;
                unsigned int& va = vao[s][slot];
                unsigned int& vb = vbo[s][slot];
                if (va) { rlUnloadVertexArray(va); va = 0; }
                if (vb) { rlUnloadVertexBuffer(vb); vb = 0; }
                count[s][slot] = (int)(buf[s][lod].size() / 8);
                if (!count[s][slot]) continue;
                if (s == 0) drawnInstances += count[s][slot];
                const VariantMesh& m = Variant(v, lod);
                va = rlLoadVertexArray();
                rlEnableVertexArray(va);
                rlEnableVertexBuffer(m.vboPos); rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(0);
                rlEnableVertexBuffer(m.vboNrm); rlSetVertexAttribute(1, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(1);
                rlEnableVertexBuffer(m.vboCol); rlSetVertexAttribute(2, 4, RL_UNSIGNED_BYTE, true, 0, 0); rlEnableVertexAttribute(2);
                vb = rlLoadVertexBuffer(buf[s][lod].data(), (int)(buf[s][lod].size() * sizeof(float)), false);
                const int stride = 8 * (int)sizeof(float);
                rlSetVertexAttribute(10, 4, RL_FLOAT, false, stride, 0); rlEnableVertexAttribute(10); rlSetVertexAttributeDivisor(10, 1);
                rlSetVertexAttribute(11, 4, RL_FLOAT, false, stride, 4 * sizeof(float)); rlEnableVertexAttribute(11); rlSetVertexAttributeDivisor(11, 1);
                rlDisableVertexArray();
                rlDisableVertexBuffer();
                if (s == 0) {
                    // #292: room for the visible instances of this slot (the view cull writes them), and a VAO on it.
                    if (culledVao[slot]) { rlUnloadVertexArray(culledVao[slot]); culledVao[slot] = 0; }
                    if (culledVbo[slot]) { rlUnloadVertexBuffer(culledVbo[slot]); culledVbo[slot] = 0; }
                    culledVao[slot] = rlLoadVertexArray();
                    rlEnableVertexArray(culledVao[slot]);
                    rlEnableVertexBuffer(m.vboPos); rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(0);
                    rlEnableVertexBuffer(m.vboNrm); rlSetVertexAttribute(1, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(1);
                    rlEnableVertexBuffer(m.vboCol); rlSetVertexAttribute(2, 4, RL_UNSIGNED_BYTE, true, 0, 0); rlEnableVertexAttribute(2);
                    culledVbo[slot] = rlLoadVertexBuffer(nullptr, (int)(buf[s][lod].size() * sizeof(float)), true);
                    rlSetVertexAttribute(10, 4, RL_FLOAT, false, stride, 0); rlEnableVertexAttribute(10); rlSetVertexAttributeDivisor(10, 1);
                    rlSetVertexAttribute(11, 4, RL_FLOAT, false, stride, 4 * sizeof(float)); rlEnableVertexAttribute(11); rlSetVertexAttributeDivisor(11, 1);
                    rlDisableVertexArray();
                    rlDisableVertexBuffer();
                }
            }
    }
}

// #292: the main view's tree culling program (tree_cull.comp), loaded once; off with RAYLIBGAME_TREE_CULL=0.
namespace {
    ComputeProgram gTreeCull;
    int gTreeCullState = 0; // 0 untried, 1 ok, -1 unavailable
    GLint gTcN = -1, gTcCmd = -1, gTcModel = -1, gTcRadius = -1, gTcPlanes = -1;
    bool TreeCullReady() {
        if (gTreeCullState == 0) {
            const char* e = getenv("RAYLIBGAME_TREE_CULL");
            gTreeCullState = (e && e[0] == '0') || !Compute_Available() || !gTreeCull.load("resources/compute/tree_cull.comp") ? -1 : 1;
            if (gTreeCullState == 1) {
                const GLuint p = gTreeCull.id;
                gTcN = glGetUniformLocation(p, "n"); gTcCmd = glGetUniformLocation(p, "cmdIndex"); gTcModel = glGetUniformLocation(p, "model");
                gTcRadius = glGetUniformLocation(p, "radius"); gTcPlanes = glGetUniformLocation(p, "planes");
            }
        }
        return gTreeCullState == 1;
    }
}

void TreeSystem::draw(const Planetoid& body, const RenderContext& ctx, int onlyLod, int set) {
    set = std::clamp(set, 0, kTreeSets - 1);
    if (!drawnInstances || !ctx.treeShader) return;
    const Matrix model = body.localModel(anchor, ctx.renderOrigin);
    rlSetUniformMatrix(ctx.treeModelLoc, model);
    const Vector3 c = (body.worldPos - ctx.renderOrigin).toFloat();
    const float occluder[4] = { c.x, c.y, c.z, std::max(0.0f, body.radius - 60.0f) * body.scale };
    rlSetUniform(ctx.treeOccluderLoc, occluder, RL_SHADER_UNIFORM_VEC4, 1);
    // #292: the main view culls on the GPU first (all slots), then draws each slot's visible instances indirectly.
    const bool gpuCull = set == 0 && ctx.frustum && TreeCullReady();
    if (gpuCull) {
        constexpr int kCmds = 2 * kTreeVariants;
        GLuint cmds[kCmds][4];
        for (int slot = 0; slot < kCmds; ++slot) {
            const int v = slot % kTreeVariants, lod = slot / kTreeVariants;
            cmds[slot][0] = (GLuint)Variant(v, lod).vertexCount; cmds[slot][1] = 0; cmds[slot][2] = 0; cmds[slot][3] = 0;
        }
        if (!cmdBuf) { glGenBuffers(1, &cmdBuf); glBindBuffer(GL_DRAW_INDIRECT_BUFFER, cmdBuf); glBufferData(GL_DRAW_INDIRECT_BUFFER, sizeof(cmds), nullptr, GL_DYNAMIC_DRAW); }
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, cmdBuf);
        glBufferSubData(GL_DRAW_INDIRECT_BUFFER, 0, sizeof(cmds), cmds);
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
        glUseProgram(gTreeCull.id);
        glUniformMatrix4fv(gTcModel, 1, GL_FALSE, MatrixToFloat(model));
        glUniform4fv(gTcPlanes, 6, (const float*)ctx.frustum->planes);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, cmdBuf);
        for (int slot = 0; slot < kCmds; ++slot) {
            const int v = slot % kTreeVariants, lod = slot / kTreeVariants;
            if (!count[0][slot] || !culledVbo[slot] || (onlyLod >= 0 && lod != onlyLod)) continue;
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, vbo[0][slot]);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, culledVbo[slot]);
            glUniform1ui(gTcN, (GLuint)count[0][slot]);
            glUniform1ui(gTcCmd, (GLuint)slot);
            glUniform1f(gTcRadius, Variant(v, lod).boundRadius);
            glDispatchCompute((GLuint)((count[0][slot] + 63) / 64), 1, 1);
        }
        glMemoryBarrier(GL_COMMAND_BARRIER_BIT | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT);
        glUseProgram(ctx.treeShader->id);
    }
    for (int slot = 0; slot < 2 * kTreeVariants; ++slot) {
        if (!count[set][slot]) continue;
        const int v = slot % kTreeVariants, lod = slot / kTreeVariants;
        if (onlyLod >= 0 && lod != onlyLod) continue;
        const float fade = v < kFirstBush ? kTreeDrawDistance : kBushDrawDistance;
        rlSetUniform(ctx.treeFadeLoc, &fade, RL_SHADER_UNIFORM_FLOAT, 1);
        if (ctx.treeBarkLoc >= 0) {
            const Vector3 one = { 1.0f, 1.0f, 1.0f };
            const Vector3 bark = v < kFirstBush ? WoodSpeciesInfo(SpeciesOfVariant(v)).bark : one;
            const Vector3 leaf = v < kFirstBush ? WoodSpeciesInfo(SpeciesOfVariant(v)).leaf : one;
            rlSetUniform(ctx.treeBarkLoc, &bark, RL_SHADER_UNIFORM_VEC3, 1);
            rlSetUniform(ctx.treeLeafLoc, &leaf, RL_SHADER_UNIFORM_VEC3, 1);
        }
        if (gpuCull && culledVao[slot]) {
            glBindVertexArray(culledVao[slot]);
            glBindBuffer(GL_DRAW_INDIRECT_BUFFER, cmdBuf);
            glDrawArraysIndirect(GL_TRIANGLES, (const void*)(size_t)(slot * 4 * sizeof(GLuint)));
            glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
            glBindVertexArray(0);
        } else {
            rlEnableVertexArray(vao[set][slot]);
            rlDrawVertexArrayInstanced(0, Variant(v, lod).vertexCount, count[set][slot]);
            rlDisableVertexArray();
        }
        if (ctx.drawCalls) ++*ctx.drawCalls;
    }
}

// The tile under a body-local direction and its neighbours, across cube-face seams too: sampled a
// little under a tile away in eight directions (the same face's i+-1, j+-1 missed the trees just
// over a seam, so you walked through them and could not chop them).
static void TilesAround(Vector3 dir, float radius, TreeTileKey out[9], int& n) {
    n = 0;
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(t1, dir);
    const float a = 0.9f * kTreeTileSize / radius;
    for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj) {
            const TreeTileKey k = TreeTileOf(Vector3Normalize(Vector3Add(dir, Vector3Add(Vector3Scale(t1, di * a), Vector3Scale(t2, dj * a)))), radius);
            bool dup = false;
            for (int m = 0; m < n; ++m) dup |= out[m] == k;
            if (!dup) out[n++] = k;
        }
}

void TreeSystem::trunksNear(Vector3 local, float range, std::vector<Trunk>& out) const {
    const float r = Vector3Length(local);
    if (r < 1.0f) return;
    TreeTileKey keys[9];
    int n = 0;
    TilesAround(Vector3Scale(local, 1.0f / r), radius, keys, n);
    for (int m = 0; m < n; ++m) {
        auto it = tiles.find(keys[m]);
        if (it == tiles.end() || !it->second.ready) continue;
        for (const Trunk& t : it->second.trunks)
            if (Vector3Distance(t.base, local) < range + t.height) out.push_back(t);
    }
}

void TreeSystem::remove(uint64_t id) {
    auto next = std::make_shared<PlantRemovals>(*removed);
    (*next)[id] = WorldTime();
    removed = next;
    // Gone at once from the loaded tile (the next rebuild agrees), and its tile starts watching
    // for regrowth.
    auto it = tiles.find(TileOfPlant(id));
    if (it == tiles.end()) return;
    TreeTile& t = it->second;
    for (int v = 0; v < kTreeVariants; ++v)
        for (size_t k = 0; k < t.ids[v].size(); ++k)
            if (t.ids[v][k] == id) {
                t.ids[v].erase(t.ids[v].begin() + k);
                t.inst[v].erase(t.inst[v].begin() + k * 8, t.inst[v].begin() + k * 8 + 8);
                break;
            }
    t.plants.erase(std::remove_if(t.plants.begin(), t.plants.end(), [id](const Plant& p) { return p.id == id; }), t.plants.end());
    t.trunks.clear();
    for (const Plant& p : t.plants) if (p.tree()) t.trunks.push_back({ p.base, p.up, p.height, p.radius });
    t.regrowing = true;
    batchDirty = true;
    // A rebuild already queued was built with the old removals: it would bring the plant back (and,
    // not regrowing, keep it). Re-request so that result is dropped (stale version) for a fresh one.
    if (t.pending) request(it->first, t, 0.0f);
}

void TreeSystem::setRemovals(const PlantRemovals& r) {
    removed = std::make_shared<PlantRemovals>(r);
    for (auto& kv : tiles) request(kv.first, kv.second, 0.0f);
}

bool TreeSystem::raycast(Vector3 o, Vector3 d, float maxDistance, PlantHit& best) const {
    const float r = Vector3Length(o);
    if (r < 1.0f) return false;
    TreeTileKey keys[9];
    int nk = 0;
    TilesAround(Vector3Scale(o, 1.0f / r), radius, keys, nk);
    bool found = false;
    best.distance = maxDistance;
    for (int m = 0; m < nk; ++m) {
            auto it = tiles.find(keys[m]);
            if (it == tiles.end() || !it->second.ready) continue;
            for (const Plant& p : it->second.plants) {
                // Closest approach between the ray and the plant's axis segment (trunk / bush core).
                const float pick = p.tree() ? std::max(p.radius, 0.35f) : p.radius;
                const Vector3 a = p.tree() ? p.base : Vector3Add(p.base, Vector3Scale(p.up, 0.2f * p.scale));
                const float len = p.tree() ? p.height * 0.8f : 0.6f * p.scale;
                const Vector3 w = Vector3Subtract(o, a);
                const float b = Vector3DotProduct(d, p.up), dd = Vector3DotProduct(d, w), e = Vector3DotProduct(p.up, w);
                const float den = 1.0f - b * b;
                float t = den > 1e-6f ? (b * e - dd) / den : 0.0f;          // along the ray
                t = std::clamp(t, 0.0f, maxDistance);
                const float s = std::clamp(e + b * t, 0.0f, len);           // along the axis
                const Vector3 pr = Vector3Add(o, Vector3Scale(d, t)), pa = Vector3Add(a, Vector3Scale(p.up, s));
                if (Vector3Distance(pr, pa) > pick) continue;
                if (t < best.distance) { best.distance = t; best.plant = p; found = true; }
            }
        }
    return found;
}

float TreeCanopy(Vector3 dir, float radius, float elevation, unsigned int seed, const PlanetField* field) {
    const float lat = fabsf(dir.y);
    if (lat > 0.86f || elevation < 4.0f || elevation > 900.0f) return 0.0f;
    const float edge = std::min(Smooth(4.0f, 12.0f, elevation), 1.0f - Smooth(800.0f, 900.0f, elevation));
    // The biome's tree factor (Plan 5 B); without a field, the old latitude rule.
    const BiomeFlora fl = field ? FloraOf(field->biomeAt(dir, elevation)) : BiomeFlora{};
    if (field && fl.cover > 0.0f) return std::max(Smooth(0.04f, 0.3f, fl.cover), Smooth(0.04f, 0.3f, TreeChanceFromForest(Forest(Vector3Scale(dir, radius), seed), fl.trees))) * edge;
    const float trees = field ? fl.trees : (lat < 0.6f ? 1.3f - 0.3f * Smooth(0.2f, 0.4f, lat) : 1.0f - 0.55f * (lat - 0.6f) / 0.26f);
    return Smooth(0.04f, 0.3f, TreeChanceFromForest(Forest(Vector3Scale(dir, radius), seed), trees)) * edge; // forest presence
}

void TreeSystem::addClearing(Vector3 c, float r) {
    auto next = std::make_shared<std::vector<Vector4>>(*clearings);
    next->push_back({ c.x, c.y, c.z, r });
    clearings = next;
    for (auto& kv : tiles)
        if (kv.second.ready && Vector3Distance(kv.second.centre, c) < r + kTreeTileSize) request(kv.first, kv.second, 0.0f);
}

bool TreeSystem::nearestPlant(Vector3 local, float range, bool tree, Plant& out, bool stone) const {
    float best = range;
    bool found = false;
    for (const auto& kv : tiles) {
        if (!kv.second.ready || Vector3Distance(kv.second.centre, local) > range + kTreeTileSize) continue;
        for (const Plant& p : kv.second.plants) {
            if (p.tree() != tree || p.stone() != stone) continue;
            const float d = Vector3Distance(p.base, local);
            if (d < best) { best = d; out = p; found = true; }
        }
    }
    return found;
}

void TreeSystem::speciesCounts(Vector3 local, float range, int out[(int)WoodSpecies::Count]) const {
    for (int s = 0; s < (int)WoodSpecies::Count; ++s) out[s] = 0;
    for (const auto& kv : tiles) {
        const TreeTile& t = kv.second;
        if (!t.ready) continue;
        for (int v = 0; v < kFirstBush; ++v)
            for (size_t k = 0; k + 7 < t.inst[v].size(); k += 8)
                if (Vector3Distance({ t.inst[v][k], t.inst[v][k + 1], t.inst[v][k + 2] }, local) < range) ++out[(int)SpeciesOfVariant(v)];
    }
}
