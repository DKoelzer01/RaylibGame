#include "trees.h"
#include "jobsystem.h"
#include "worldtypes.h"
#include "world.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <mutex>
#include <random>

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
};
VariantMesh gVariants[2][kTreeVariants]; // [0] full, [1] far (a few prisms)
std::once_flag gVariantsOnce;

void Emit(VariantMesh& m, Vector3 p, Vector3 n, bool leaf, unsigned char shade) {
    m.pos.insert(m.pos.end(), { p.x, p.y, p.z });
    m.nrm.insert(m.nrm.end(), { n.x, n.y, n.z });
    m.col.insert(m.col.end(), { (unsigned char)(leaf ? 255 : 0), shade, 0, 255 });
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

void Broadleaf(VariantMesh& m, unsigned int seed) {
    std::mt19937 r(seed);
    const float H = 4.2f + 1.6f * U(r);
    const Vector3 lean = { (U(r) - 0.5f) * 0.6f, 0.0f, (U(r) - 0.5f) * 0.6f };
    auto trunkAt = [&](float y) { const float f = y / H; return Vector3{ lean.x * f * f, y, lean.z * f * f }; };
    // Trunk: three stacked tapering segments, slightly bent.
    const float ys[4] = { -0.4f, H * 0.35f, H * 0.7f, H };
    const float ws[4] = { 0.62f, 0.45f, 0.34f, 0.24f };
    for (int k = 0; k < 3; ++k) Prism(m, trunkAt(ys[k]), trunkAt(ys[k + 1]), ws[k], ws[k + 1], 0.2f * k, false, Shade(r));
    // Branches up and out, each ending in a leaf block.
    const int nb = 3 + (int)(U(r) * 3.0f);
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
        const float c = 1.7f + 1.1f * U(r);
        Prism(m, Vector3Subtract(e, { 0, c * 0.4f, 0 }), Vector3Add(e, { 0, c * 0.6f, 0 }), c, c * 0.85f, U(r) * 1.5f, true, Shade(r));
    }
    // Crown.
    const float c = 2.5f + 0.9f * U(r);
    const Vector3 top = trunkAt(H);
    Prism(m, Vector3Subtract(top, { 0, 0.9f, 0 }), Vector3Add(top, { 0, c - 0.6f, 0 }), c, c * 0.8f, U(r) * 1.5f, true, Shade(r));
    m.trunkRadius = 0.3f;
    m.trunkHeight = H;
}

void Conifer(VariantMesh& m, unsigned int seed) {
    std::mt19937 r(seed);
    const float H = 8.5f + 2.0f * U(r);
    Prism(m, { 0, -0.4f, 0 }, { 0, H, 0 }, 0.5f, 0.12f, 0.0f, false, Shade(r));
    // Stacked leaf slabs, wide at the bottom, turned alternately.
    const int n = 6 + (int)(U(r) * 2.0f);
    for (int k = 0; k < n; ++k) {
        const float f = (float)k / (n - 1);
        const float y = 1.6f + f * (H - 2.2f);
        const float w = 3.8f * (1.0f - f) + 0.8f;
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
    if (v >= 4 && v < kFirstBush) {
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
        if (v < 4) Broadleaf(m, 1000u + v * 7919u);
        else if (v < kFirstBush) Conifer(m, 2000u + v * 7919u);
        else Bush(m, 3000u + v * 7919u);
        m.vertexCount = (int)(m.pos.size() / 3);
        FarFrom(m, v, gVariants[1][v]);
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

void BuildTile(const PlanetField& f, const EditList& edits, uint32_t seed, const TreeTileKey& key, TreeTile& out) {
    const float R = f.radius;
    const int N = TilesPerFace(R);
    const float da = 0.5f * PI / N, dc = da / kCells;
    const float a0 = -0.25f * PI + key.i * da, b0 = -0.25f * PI + key.j * da;
    {
        const Vector3 cd = FaceDir(key.face, a0 + 0.5f * da, b0 + 0.5f * da);
        out.centre = Vector3Scale(cd, R + std::max(0.0f, f.elevation(cd, 16.0f)));
    }
    auto solid = [&](Vector3 p) { return f.density(p, edits.data(), edits.size()) >= ISO_LEVEL; };
    for (int ci = 0; ci < kCells; ++ci)
        for (int cj = 0; cj < kCells; ++cj) {
            uint32_t s = Hash(seed ^ Hash((uint32_t)key.face * 2654435761U ^ (uint32_t)(key.i * kCells + ci) * 40503U ^ (uint32_t)(key.j * kCells + cj) * 9973U));
            const float jx = Unit(s), jy = Unit(s), roll = Unit(s), h = Unit(s), sz = Unit(s);
            const Vector3 dir = FaceDir(key.face, a0 + (ci + jx) * dc, b0 + (cj + jy) * dc);
            const float forest = Forest(Vector3Scale(dir, R), seed);
            const float pTree = 0.55f * Smooth(0.48f, 0.72f, forest) + 0.015f;
            const float pBush = 0.05f + 0.3f * Smooth(0.35f, 0.6f, forest);
            const bool tree = roll < pTree;
            if (!tree && roll >= pTree + pBush) continue;
            const float lat = fabsf(dir.y);
            if (lat > 0.86f) continue;                                  // polar: nothing grows
            const float e = f.elevation(dir);
            if (e < 4.0f || e > (tree ? 900.0f : 1100.0f)) continue;   // beach / sea, above the tree line
            // Slope from the height a few metres away.
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const Vector3 t2 = Vector3CrossProduct(t1, dir);
            const float step = 3.0f, ang = step / (R + e);
            const float e1 = f.elevation(Vector3Normalize(Vector3Add(dir, Vector3Scale(t1, ang))));
            const float e2 = f.elevation(Vector3Normalize(Vector3Add(dir, Vector3Scale(t2, ang))));
            const Vector3 normal = Vector3Normalize(Vector3Subtract(dir, Vector3Add(Vector3Scale(t1, (e1 - e) / step), Vector3Scale(t2, (e2 - e) / step))));
            const float slope = 1.0f - Vector3DotProduct(normal, dir);
            if (slope > (tree ? 0.12f : 0.2f)) continue;
            const Vector3 ground = Vector3Scale(dir, R + e);
            const VoxelMaterial m = f.materialLod(Vector3Subtract(ground, Vector3Scale(dir, 0.3f)), 1.0f, 0.0f, normal);
            if (m != VoxelMaterial::Grass && m != VoxelMaterial::Regolith) continue;
            // Solid ground underfoot and open air above (caves, ravines, edits).
            if (!solid(Vector3Subtract(ground, Vector3Scale(dir, 1.0f))) || !solid(Vector3Subtract(ground, Vector3Scale(dir, 4.0f)))) continue;
            if (solid(Vector3Add(ground, Vector3Scale(dir, 1.5f))) || (tree && solid(Vector3Add(ground, Vector3Scale(dir, 6.0f))))) continue;
            int v;
            float scale;
            if (tree) {
                // Broadleaf in the temperate band, mixing into conifers from ~40 degrees and up high.
                const bool conifer = lat > 0.66f + 0.12f * (h - 0.5f) || e > 550.0f + 300.0f * (h - 0.5f);
                v = conifer ? 4 + (int)(h * 2.0f) % 2 : (int)(h * 4.0f) % 4;
                scale = 0.6f + 0.7f * sz + 0.35f * forest;
            } else {
                v = kFirstBush + (int)(h * 3.0f) % 3;
                scale = 0.7f + 0.7f * sz;
            }
            const Vector3 p = Vector3Subtract(ground, Vector3Scale(dir, 0.15f)); // sink the base a little
            out.inst[v].insert(out.inst[v].end(), { p.x, p.y, p.z, scale, dir.x, dir.y, dir.z, h });
            const VariantMesh& vm = Variant(v);
            if (vm.trunkRadius > 0.0f) out.trunks.push_back({ p, dir, vm.trunkHeight * scale, vm.trunkRadius * scale });
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
    : shared(std::make_shared<Shared>()), field(std::move(f)), seed(s), radius(field->radius), edits(std::make_shared<EditList>()) {}

TreeSystem::~TreeSystem() {
    shared->alive = false;
    for (int v = 0; v < kTreeVariants; ++v) {
        if (vao[v]) rlUnloadVertexArray(vao[v]);
        if (vao[v + kTreeVariants]) rlUnloadVertexArray(vao[v + kTreeVariants]);
        if (vbo[v + kTreeVariants]) rlUnloadVertexBuffer(vbo[v + kTreeVariants]);
        if (vbo[v]) rlUnloadVertexBuffer(vbo[v]);
    }
}

void TreeSystem::request(const TreeTileKey& key, TreeTile& tile, float priority) {
    const uint32_t version = ++tile.version;
    auto sh = shared;
    auto fld = field;
    auto ed = edits;
    const uint32_t sd = seed;
    WorldJobs().submit(priority, [sh, fld, ed, sd, key, version] {
        if (!sh->alive) return;
        TreeTile t;
        t.version = version;
        t.ready = true;
        BuildTile(*fld, *ed, sd, key, t);
        std::lock_guard<std::mutex> lock(sh->mutex);
        sh->results.emplace_back(key, std::move(t));
    });
}

void TreeSystem::onEdit(Vector3 c, float r) {
    for (auto& kv : tiles)
        if (Vector3Distance(kv.second.centre.x != 0.0f || kv.second.centre.y != 0.0f || kv.second.centre.z != 0.0f
                                ? kv.second.centre : Vector3Scale(TreeTileCentreDir(kv.first, radius), radius), c) < r + kTreeTileSize)
            request(kv.first, kv.second, 0.0f); // the old trees stay until the rebuild arrives
}

void TreeSystem::update(Planetoid& body, Vector3 camLocal) {
    if (body.edits != edits) {
        // Edits changed since last frame: rebuild the tiles each new edit reaches.
        auto old = edits;
        edits = body.edits;
        for (size_t k = old->size(); k < edits->size(); ++k) onEdit((*edits)[k].center, (*edits)[k].radius);
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
    std::vector<float> buf[2];
    for (int v = 0; v < kTreeVariants; ++v) {
        const float reach = v < kFirstBush ? kTreeDrawDistance : kBushDrawDistance;
        buf[0].clear(); buf[1].clear();
        for (auto& kv : tiles) {
            const TreeTile& t = kv.second;
            if (!t.ready || t.inst[v].empty() || Vector3Distance(t.centre, camLocal) > reach + kTreeTileSize) continue;
            for (size_t k = 0; k + 7 < t.inst[v].size(); k += 8) {
                const float x = (float)(t.inst[v][k] - anchor.x), y = (float)(t.inst[v][k + 1] - anchor.y), z = (float)(t.inst[v][k + 2] - anchor.z);
                const float d2 = x * x + y * y + z * z;
                if (d2 > reach * reach) continue;
                std::vector<float>& b = buf[d2 < kTreeNearDistance * kTreeNearDistance ? 0 : 1];
                b.insert(b.end(), { x, y, z, t.inst[v][k + 3], t.inst[v][k + 4], t.inst[v][k + 5], t.inst[v][k + 6], t.inst[v][k + 7] });
            }
        }
        for (int lod = 0; lod < 2; ++lod) {
            const int slot = v + lod * kTreeVariants;
            if (vao[slot]) { rlUnloadVertexArray(vao[slot]); vao[slot] = 0; }
            if (vbo[slot]) { rlUnloadVertexBuffer(vbo[slot]); vbo[slot] = 0; }
            count[slot] = (int)(buf[lod].size() / 8);
            if (!count[slot]) continue;
            drawnInstances += count[slot];
            const VariantMesh& m = Variant(v, lod);
            vao[slot] = rlLoadVertexArray();
            rlEnableVertexArray(vao[slot]);
            rlEnableVertexBuffer(m.vboPos); rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(0);
            rlEnableVertexBuffer(m.vboNrm); rlSetVertexAttribute(1, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(1);
            rlEnableVertexBuffer(m.vboCol); rlSetVertexAttribute(2, 4, RL_UNSIGNED_BYTE, true, 0, 0); rlEnableVertexAttribute(2);
            vbo[slot] = rlLoadVertexBuffer(buf[lod].data(), (int)(buf[lod].size() * sizeof(float)), false);
            const int stride = 8 * (int)sizeof(float);
            rlSetVertexAttribute(10, 4, RL_FLOAT, false, stride, 0); rlEnableVertexAttribute(10); rlSetVertexAttributeDivisor(10, 1);
            rlSetVertexAttribute(11, 4, RL_FLOAT, false, stride, 4 * sizeof(float)); rlEnableVertexAttribute(11); rlSetVertexAttributeDivisor(11, 1);
            rlDisableVertexArray();
            rlDisableVertexBuffer();
        }
    }
}

void TreeSystem::draw(const Planetoid& body, const RenderContext& ctx) {
    if (!drawnInstances || !ctx.treeShader) return;
    const Matrix model = body.localModel(anchor, ctx.renderOrigin);
    rlSetUniformMatrix(ctx.treeModelLoc, model);
    const Vector3 c = (body.worldPos - ctx.renderOrigin).toFloat();
    const float occluder[4] = { c.x, c.y, c.z, std::max(0.0f, body.radius - 60.0f) * body.scale };
    rlSetUniform(ctx.treeOccluderLoc, occluder, RL_SHADER_UNIFORM_VEC4, 1);
    for (int slot = 0; slot < 2 * kTreeVariants; ++slot) {
        if (!count[slot]) continue;
        const int v = slot % kTreeVariants, lod = slot / kTreeVariants;
        const float fade = v < kFirstBush ? kTreeDrawDistance : kBushDrawDistance;
        rlSetUniform(ctx.treeFadeLoc, &fade, RL_SHADER_UNIFORM_FLOAT, 1);
        rlEnableVertexArray(vao[slot]);
        rlDrawVertexArrayInstanced(0, Variant(v, lod).vertexCount, count[slot]);
        rlDisableVertexArray();
        if (ctx.drawCalls) ++*ctx.drawCalls;
    }
}

void TreeSystem::trunksNear(Vector3 local, float range, std::vector<Trunk>& out) const {
    const float r = Vector3Length(local);
    if (r < 1.0f) return;
    const TreeTileKey k = TreeTileOf(Vector3Scale(local, 1.0f / r), radius);
    for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj) {
            auto it = tiles.find({ k.face, k.i + di, k.j + dj });
            if (it == tiles.end() || !it->second.ready) continue;
            for (const Trunk& t : it->second.trunks)
                if (Vector3Distance(t.base, local) < range + t.height) out.push_back(t);
        }
}
