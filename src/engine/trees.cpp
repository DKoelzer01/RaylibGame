#include "trees.h"
#include "biomevariants.h"
#include "envflag.h"
#include "gpumem.h"
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
VariantMesh gVariants[2][kTreeMeshes]; // [0] full, [1] far (a few prisms)
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

// (#933 layer 2) Ambient occlusion baked into the vertex colour's alpha (255 = open; trees.fs scales the ambient term
// only). Leaves darken toward the crown's middle (the inner canopy sees little sky) and on their undersides; bark toward
// its foot (the ground's contact) and where it stands under the crown; loose stones at the ground. Analytic, from the
// variant's own geometry: no rays. RAYLIBGAME_FOLIAGE_AO=0: off (all 255).
void BakeFoliageAo(VariantMesh& m, bool stone) {
    static const bool on = !(getenv("RAYLIBGAME_FOLIAGE_AO") && getenv("RAYLIBGAME_FOLIAGE_AO")[0] == '0');
    const size_t nv = m.pos.size() / 3;
    if (!on || nv == 0) return;
    // The crown: the leaf vertices' middle and their reach from it (horizontal and vertical apart: crowns are wide).
    Vector3 c{}; int nl = 0; float top = 0.0f;
    for (size_t i = 0; i < nv; ++i) {
        top = std::max(top, m.pos[3 * i + 1]);
        if (m.col[4 * i] > 127) { c = Vector3Add(c, { m.pos[3 * i], m.pos[3 * i + 1], m.pos[3 * i + 2] }); ++nl; }
    }
    float rh = 0.0f, rv = 0.0f, lowLeaf = 1e9f;
    if (nl > 0) {
        c = Vector3Scale(c, 1.0f / nl);
        for (size_t i = 0; i < nv; ++i) if (m.col[4 * i] > 127) {
            const float dx = m.pos[3 * i] - c.x, dz = m.pos[3 * i + 2] - c.z;
            rh = std::max(rh, sqrtf(dx * dx + dz * dz)); rv = std::max(rv, fabsf(m.pos[3 * i + 1] - c.y));
            lowLeaf = std::min(lowLeaf, m.pos[3 * i + 1]);
        }
    }
    auto ss = [](float a, float b, float x) { const float t = std::clamp((x - a) / (b - a), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); };
    for (size_t i = 0; i < nv; ++i) {
        const Vector3 p = { m.pos[3 * i], m.pos[3 * i + 1], m.pos[3 * i + 2] };
        const float ny = m.nrm[3 * i + 1];
        float ao = 1.0f;
        if (stone) {
            ao = 0.6f + 0.4f * ss(0.0f, 0.25f, p.y);                                   // the stone's foot against the ground
        } else if (m.col[4 * i] > 127 && nl > 0) {
            const float dx = (p.x - c.x) / std::max(rh, 0.1f), dz = (p.z - c.z) / std::max(rh, 0.1f), dy = (p.y - c.y) / std::max(rv, 0.1f);
            const float d = sqrtf(dx * dx + dy * dy + dz * dz);                          // 0 in the crown's middle, ~1 at its edge
            ao = 0.5f + 0.5f * ss(0.15f, 0.95f, d);
            // A crown's lower half sees less sky than its top (the canopy above it): blocky crowns show only their outer
            // faces (all near the edge, d ~ 1), so this is the part of the AO that reads.
            const float h = std::clamp((p.y - (c.y - rv)) / std::max(2.0f * rv, 0.1f), 0.0f, 1.0f); // 0 bottom .. 1 top
            ao = std::min(ao, 0.6f + 0.4f * ss(0.05f, 0.85f, h));
            if (ny < -0.5f) ao *= 0.78f;                                                 // leaf undersides face the dark ground
        } else {
            ao = 0.55f + 0.45f * ss(-0.2f, 0.9f, p.y);                                   // the trunk's foot
            if (nl > 0 && p.y > lowLeaf - 0.5f) {                                        // bark up in the crown
                const float dx = p.x - c.x, dz = p.z - c.z;
                if (dx * dx + dz * dz < rh * rh) ao *= 0.72f;
            }
        }
        static const float strength = getenv("RAYLIBGAME_FOLIAGE_AO_STRENGTH") ? (float)atof(getenv("RAYLIBGAME_FOLIAGE_AO_STRENGTH")) : 2.0f; // (#933) 2 = the default (user, 2026-11-19); 1 = the first tuning
        ao = std::clamp(1.0f - strength * (1.0f - ao), 0.0f, 1.0f);
        m.col[4 * i + 3] = (unsigned char)std::clamp((int)lroundf(ao * 255.0f), 0, 255);
    }
}

float U(std::mt19937& r) { return std::uniform_real_distribution<float>(0.0f, 1.0f)(r); }
unsigned char Shade(std::mt19937& r) { return (unsigned char)(90 + (int)(U(r) * 140.0f)); }

// Shape knobs per species and size class (#274): height, crown size, trunk thickness, branch count and reach; old
// giants add buttress roots and a second, lower ring of branches.
struct BroadleafShape { float h = 1.0f, crown = 1.0f, trunk = 1.0f; int branches = 0; float reach = 1.0f; bool giant = false; };
struct ConiferShape { float h = 1.0f, width = 1.0f; int layers = 0; float trunk = 1.0f; bool giant = false; float k = 1.0f; }; // k: the class scale (slabs, tip)

// Buttress roots: flared fins from the trunk's foot out along the ground (old-growth giants).
void Buttress(VariantMesh& m, std::mt19937& r, float trunkW, float h) {
    const int n = 4 + (int)(U(r) * 2.0f);
    const float yaw0 = U(r) * 2.0f * PI;
    for (int k = 0; k < n; ++k) {
        const float yaw = yaw0 + k * 2.0f * PI / n + (U(r) - 0.5f) * 0.5f;
        const float out = trunkW * (1.3f + 0.6f * U(r));
        const Vector3 foot = { cosf(yaw) * out, -0.3f, sinf(yaw) * out };
        Prism(m, foot, { cosf(yaw) * trunkW * 0.2f, h * (0.6f + 0.3f * U(r)), sinf(yaw) * trunkW * 0.2f }, trunkW * 0.45f, trunkW * 0.2f, yaw, false, Shade(r));
    }
}

void Broadleaf(VariantMesh& m, unsigned int seed, BroadleafShape sh = {}) {
    std::mt19937 r(seed);
    const float H = (4.2f + 1.6f * U(r)) * sh.h;
    const Vector3 lean = { (U(r) - 0.5f) * 0.6f, 0.0f, (U(r) - 0.5f) * 0.6f };
    auto trunkAt = [&](float y) { const float f = y / H; return Vector3{ lean.x * f * f, y, lean.z * f * f }; };
    // Trunk: three stacked tapering segments, slightly bent.
    const float ys[4] = { -0.4f, H * 0.35f, H * 0.7f, H };
    const float ws[4] = { 0.62f * sh.trunk, 0.45f * sh.trunk, 0.34f * sh.trunk, 0.24f * sh.trunk };
    for (int s = 0; s < 3; ++s) Prism(m, trunkAt(ys[s]), trunkAt(ys[s + 1]), ws[s], ws[s + 1], 0.2f * s, false, Shade(r));
    if (sh.giant) Buttress(m, r, ws[0], H * 0.12f);
    // Branches up and out, each ending in a leaf block; giants add a lower ring, spreading wider.
    for (int tier = 0; tier < (sh.giant ? 2 : 1); ++tier) {
        const int nb = std::max(2, 3 + (int)(U(r) * 3.0f) + sh.branches - tier * 2);
        const float yaw0 = U(r) * 2.0f * PI;
        for (int k = 0; k < nb; ++k) {
            const float y = tier == 0 ? H * (0.5f + 0.4f * U(r)) : H * (0.3f + 0.15f * U(r));
            const float yaw = yaw0 + k * 2.0f * PI / nb + (U(r) - 0.5f) * 0.6f;
            const float el = (tier == 0 ? 35.0f + 25.0f * U(r) : 15.0f + 15.0f * U(r)) * DEG2RAD;
            const float len = (1.5f + 1.3f * U(r)) * sh.reach * (tier ? 1.2f : 1.0f);
            const Vector3 s = trunkAt(y);
            const Vector3 dir = { cosf(el) * cosf(yaw), sinf(el), cosf(el) * sinf(yaw) };
            const Vector3 e = Vector3Add(s, Vector3Scale(dir, len));
            Prism(m, s, e, 0.2f * std::sqrt(sh.trunk), 0.12f * std::sqrt(sh.trunk), U(r), false, Shade(r));
            const float c = (1.7f + 1.1f * U(r)) * sh.crown * (tier ? 0.8f : 1.0f);
            Prism(m, Vector3Subtract(e, { 0, c * 0.4f, 0 }), Vector3Add(e, { 0, c * 0.6f, 0 }), c, c * 0.85f, U(r) * 1.5f, true, Shade(r));
        }
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
    Prism(m, { 0, -0.4f, 0 }, { 0, H, 0 }, 0.5f * sh.trunk, 0.12f * sh.trunk, 0.0f, false, Shade(r));
    if (sh.giant) Buttress(m, r, 0.5f * sh.trunk, H * 0.08f);
    // Stacked leaf slabs, wide at the bottom, turned alternately; old giants keep a bare lower trunk.
    const int n = std::max(3, 6 + (int)(U(r) * 2.0f) + sh.layers);
    const float y0 = sh.giant ? H * 0.28f : std::min(1.6f, H * 0.2f);
    const float slab = 1.5f * sh.k;
    for (int k = 0; k < n; ++k) {
        const float f = (float)k / (n - 1);
        const float y = y0 + f * (H - y0 - 0.6f);
        const float w = (3.8f * (1.0f - f) + 0.8f) * sh.width;
        Prism(m, { 0, y, 0 }, { 0, y + slab, 0 }, w, w * 0.5f, k * 0.785f + U(r) * 0.3f, true, Shade(r));
    }
    Prism(m, { 0, H - 0.2f, 0 }, { 0, H + 0.9f * sh.k, 0 }, 0.7f * sh.k, 0.15f, 0.4f, true, Shade(r));
    m.trunkRadius = 0.25f * sh.trunk;
    m.trunkHeight = H;
}

} // namespace
bool TreeFormsOn() { static const bool on = !(getenv("RAYLIBGAME_TREE_FORMS") && !EnvFlag("RAYLIBGAME_TREE_FORMS")); return on; }
namespace {
// ---- #1995: shapes by species and by age, three seeded forms each (RAYLIBGAME_TREE_FORMS=0: the shapes above) ----
// Every species has its own structure (oak spreading and heavy, birch a slender column, maple round and full; pine tall
// with a bare trunk and a few flat whorls, spruce a dense narrow cone) and each size class its own form, not only its
// size: 0 sapling (a stem and a tuft), 1 young (few short branches, a narrow crown), 2 mature (the full form, forked
// branches), 3 old (thick, buttressed, wide and broken crowns with dead limbs). Numbers are GUESSes for playtests.
struct BroadSpecies { float h, crown, trunk, reach, elevLo, elevHi, yLo, yHi; int branches; float blockW, blockH; };
const BroadSpecies kBroadV2[3] = {
    { 1.00f, 1.25f, 1.30f, 1.30f, 15.0f, 38.0f, 0.42f, 0.78f, 2, 1.30f, 0.65f }, // oak: low, wide, flat-topped blocks
    { 1.40f, 0.62f, 0.65f, 0.55f, 55.0f, 75.0f, 0.35f, 0.95f, 2, 0.75f, 1.15f }, // birch: steep short twigs all up the stem
    { 1.05f, 1.05f, 1.00f, 0.90f, 32.0f, 55.0f, 0.55f, 0.85f, 1, 1.00f, 0.95f }, // maple: a round cluster
};
void BroadleafV2(VariantMesh& m, unsigned int seed, int sp, int cls) {
    std::mt19937 r(seed);
    const BroadSpecies& p = kBroadV2[sp];
    static const float kH[4] = { 0.4f, 0.8f, 1.15f, 1.6f }, kT[4] = { 0.35f, 0.75f, 1.25f, 2.1f }, kC[4] = { 0.5f, 0.8f, 1.1f, 1.4f };
    const float H = (4.2f + 1.6f * U(r)) * p.h * kH[cls];
    const float trunk = p.trunk * kT[cls];
    const Vector3 lean = { (U(r) - 0.5f) * 0.7f, 0.0f, (U(r) - 0.5f) * 0.7f };
    auto trunkAt = [&](float y) { const float f = y / H; return Vector3{ lean.x * f * f, y, lean.z * f * f }; };
    auto leafBlock = [&](Vector3 at, float c) {
        const float w = c * p.blockW, h = c * p.blockH;
        Prism(m, Vector3Subtract(at, { 0, h * 0.4f, 0 }), Vector3Add(at, { 0, h * 0.6f, 0 }), w, w * 0.85f, U(r) * 1.5f, true, Shade(r));
    };
    m.trunkRadius = 0.3f * trunk;
    m.trunkHeight = H;
    if (cls == 0) { // sapling: one thin stem and a tuft (birch: two)
        Prism(m, { 0, -0.4f, 0 }, trunkAt(H), 0.3f * trunk + 0.06f, 0.06f, 0.0f, false, Shade(r));
        const float c = (1.6f + 0.5f * U(r)) * p.crown * kC[0];
        if (sp == 1) { leafBlock(trunkAt(H * 0.62f), c * 0.8f); leafBlock(trunkAt(H * 0.95f), c * 0.7f); }
        else leafBlock(trunkAt(H * 0.92f), c);
        return;
    }
    const float ys[4] = { -0.4f, H * 0.35f, H * 0.7f, H };
    const float ws[4] = { 0.62f * trunk, 0.45f * trunk, 0.34f * trunk, 0.24f * trunk };
    for (int s = 0; s < 3; ++s) Prism(m, trunkAt(ys[s]), trunkAt(ys[s + 1]), ws[s], ws[s + 1], 0.2f * s, false, Shade(r));
    if (cls == 3) Buttress(m, r, ws[0], H * 0.12f);
    const int nb = std::max(2, 1 + p.branches + cls + (int)(U(r) * 2.0f));
    const float yaw0 = U(r) * 2.0f * PI;
    int dead = cls == 3 ? 1 + (int)(U(r) * 2.0f) : 0; // old trees carry a dead limb or two
    for (int k = 0; k < nb; ++k) {
        const float y = H * (p.yLo + (p.yHi - p.yLo) * (nb > 1 ? (float)k / (nb - 1) : 0.5f) * (0.85f + 0.15f * U(r)));
        const float yaw = yaw0 + k * 2.39996f + (U(r) - 0.5f) * 0.5f; // golden-angle spiral: no two branches line up
        const float el = (p.elevLo + (p.elevHi - p.elevLo) * U(r)) * DEG2RAD;
        const float len = (1.3f + 1.0f * U(r)) * p.reach * (0.6f + 0.25f * cls);
        const Vector3 s = trunkAt(y);
        const Vector3 dir = { cosf(el) * cosf(yaw), sinf(el), cosf(el) * sinf(yaw) };
        const Vector3 e = Vector3Add(s, Vector3Scale(dir, len));
        const float bw = 0.2f * std::sqrt(trunk);
        if (dead > 0 && k < nb - 1 && U(r) < 0.6f) { // bare, thinner, drooping at the tip
            --dead;
            const Vector3 tip = Vector3Add(e, { dir.x * len * 0.3f, -0.4f * len * 0.3f, dir.z * len * 0.3f });
            Prism(m, s, e, bw * 0.8f, bw * 0.5f, U(r), false, Shade(r));
            Prism(m, e, tip, bw * 0.5f, bw * 0.3f, U(r), false, Shade(r));
            continue;
        }
        Prism(m, s, e, bw, bw * 0.6f, U(r), false, Shade(r));
        const float c = (1.4f + 0.9f * U(r)) * p.crown * kC[cls];
        leafBlock(e, c);
        if (cls >= 2 && U(r) < 0.65f) { // a fork: a second, smaller block off the branch's end (two-level structure)
            const float yaw2 = yaw + (U(r) < 0.5f ? -0.9f : 0.9f);
            const Vector3 d2 = { cosf(el * 0.7f) * cosf(yaw2), sinf(el * 0.7f), cosf(el * 0.7f) * sinf(yaw2) };
            const Vector3 e2 = Vector3Add(Vector3Lerp(s, e, 0.6f), Vector3Scale(d2, len * 0.6f));
            Prism(m, Vector3Lerp(s, e, 0.6f), e2, bw * 0.6f, bw * 0.4f, U(r), false, Shade(r));
            leafBlock(e2, c * 0.7f);
        }
    }
    // The crown at the top: oak broad and flat, birch narrow, maple a cluster of three; old crowns are sometimes broken.
    if (cls == 3 && U(r) < 0.35f) return;
    const Vector3 top = trunkAt(H);
    const float c = (2.3f + 0.8f * U(r)) * p.crown * kC[cls];
    if (sp == 2) {
        for (int k = 0; k < 3; ++k) {
            const float a = yaw0 + k * 2.0944f;
            leafBlock(Vector3Add(top, { cosf(a) * c * 0.35f, c * 0.1f * (k - 1), sinf(a) * c * 0.35f }), c * 0.75f);
        }
    } else {
        leafBlock(Vector3Add(top, { 0, c * 0.2f, 0 }), c);
    }
}

struct ConeSpecies { float h, width, trunk, bare; int layers; float droop, offset; };
const ConeSpecies kConeV2[2] = {
    { 1.25f, 0.85f, 0.95f, 0.45f, 4, 0.0f, 0.35f }, // pine: a bare trunk, a few flat, uneven whorls high up
    { 0.95f, 0.95f, 1.00f, 0.08f, 8, 0.3f, 0.06f }, // spruce: dense, narrow, down to the ground
};
void ConiferV2(VariantMesh& m, unsigned int seed, int sp, int cls) {
    std::mt19937 r(seed);
    const ConeSpecies& p = kConeV2[sp];
    static const float kH[4] = { 0.35f, 0.75f, 1.3f, 1.8f }, kW[4] = { 0.45f, 0.8f, 1.15f, 1.35f }, kT[4] = { 0.45f, 0.8f, 1.4f, 2.1f }, kK[4] = { 0.55f, 0.85f, 1.15f, 1.4f };
    static const int kL[4] = { -1, 0, 1, 2 };
    const float H = (8.5f + 2.0f * U(r)) * p.h * kH[cls];
    const float tw = 0.5f * p.trunk * kT[cls];
    Prism(m, { 0, -0.4f, 0 }, { 0, H, 0 }, tw, std::max(0.06f, 0.12f * kT[cls]), 0.0f, false, Shade(r));
    if (cls == 3) Buttress(m, r, tw, H * 0.08f);
    m.trunkRadius = 0.5f * tw;
    m.trunkHeight = H;
    const int n = std::max(3, p.layers + kL[cls] + (int)(U(r) * 2.0f));
    const float bare = cls == 0 ? 0.15f : p.bare + (cls == 3 ? 0.12f : 0.0f);
    const float y0 = H * bare;
    const float slab = 1.5f * kK[cls] * (sp == 0 ? 0.7f : 1.0f);
    int gaps = cls == 3 && sp == 1 ? 1 : 0;
    for (int k = 0; k < n; ++k) {
        const float f = (float)k / (n - 1);
        const float y = y0 + f * (H - y0 - 0.6f);
        float w = (3.8f * (1.0f - f) + 0.8f) * p.width * kW[cls];
        if (sp == 0) w = (2.2f + 1.6f * U(r)) * kW[cls] * (1.0f - 0.5f * f); // pine whorls: about even, irregular
        if (gaps > 0 && k > 0 && k < n - 2 && U(r) < 0.3f) { --gaps; continue; } // a gap where a whorl broke off
        // Off-centre whorls (pine): each on a short branch, so the crown is ragged rather than a stack.
        const float a = U(r) * 2.0f * PI, off = p.offset * w * U(r);
        const Vector3 c = { cosf(a) * off, y, sinf(a) * off };
        if (off > 0.3f) Prism(m, { 0, y - 0.2f, 0 }, c, tw * 0.35f, tw * 0.25f, a, false, Shade(r));
        Prism(m, Vector3Subtract(c, { 0, slab * p.droop, 0 }), Vector3Add(c, { 0, slab, 0 }), w, w * 0.5f, k * 0.785f + U(r) * 0.3f, true, Shade(r));
    }
    if (!(cls == 3 && sp == 0)) // old pines end flat; the rest in a spire
        Prism(m, { 0, H - 0.2f, 0 }, { 0, H + 0.9f * kK[cls], 0 }, 0.7f * kK[cls], 0.15f, 0.4f, true, Shade(r));
}

// Bushes (#1995): three kinds - a round shrub, a low spreading mat and an upright clump.
void BushV2(VariantMesh& m, unsigned int seed, int kind) {
    std::mt19937 r(seed);
    const int n = (kind == 1 ? 3 : 2) + (int)(U(r) * 3.0f);
    for (int k = 0; k < n; ++k) {
        float c = 0.7f + 0.6f * U(r), h = c;
        Vector3 o = { (U(r) - 0.5f) * 0.9f, 0.15f + 0.4f * U(r), (U(r) - 0.5f) * 0.9f };
        if (kind == 1) { c *= 1.25f; h = c * 0.45f; o.x *= 1.6f; o.z *= 1.6f; o.y = 0.1f + 0.15f * U(r); }
        if (kind == 2) { c *= 0.7f; h = c * 1.6f; o.x *= 0.5f; o.z *= 0.5f; o.y = 0.2f + 0.5f * k; }
        Prism(m, Vector3Subtract(o, { 0, h * 0.45f, 0 }), Vector3Add(o, { 0, h * 0.55f, 0 }), c, c * 0.8f, U(r) * 1.5f, true, Shade(r));
    }
    if (kind == 2) Prism(m, { 0, -0.1f, 0 }, { 0, 0.4f, 0 }, 0.12f, 0.08f, 0.0f, false, Shade(r)); // a woody stem under the clump
    m.trunkRadius = 0.0f;
    m.trunkHeight = 0.0f;
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
    const bool forms = TreeFormsOn();
    for (int mesh = 0; mesh < kTreeMeshes; ++mesh) {
        const int v = VariantOfTreeMesh(mesh), form = mesh % kTreeForms;
        VariantMesh& m = gVariants[0][mesh];
        if (forms && v != kStoneVariant) {
            // #1995: form 0 keeps the variant's old seed; the others are their own trees.
            const unsigned int fs = (unsigned)form * 104729u;
            if (v < kFirstBush) {
                const int sp = (int)SpeciesOfVariant(v), cls = SizeClassOfVariant(v);
                if (WoodSpeciesInfo((WoodSpecies)sp).conifer) ConiferV2(m, 2000u + v * 7919u + fs, sp - (int)WoodSpecies::Pine, cls);
                else BroadleafV2(m, 1000u + v * 7919u + fs, sp, cls);
            } else {
                BushV2(m, 3000u + v * 7919u + fs, v - kFirstBush);
            }
        } else if (form > 0) {
            m = gVariants[0][TreeMeshOf(v, 0)]; // without forms (and for stones, already varied by their scatter) form 0 stands in
            m.vboPos = m.vboNrm = m.vboCol = 0;
        } else if (v < kFirstBush) {
            // The species' shape (oak broad and heavy, birch tall and slender, maple round and full; pine tall and
            // sparse, spruce dense and broad), then the size class's rules (#274): small = young and thin with few
            // branches; medium = the species as it was; large = taller, thicker, more and longer branches; extra
            // large = old-growth giants with buttress roots, a wide lower ring of branches (broadleaf) or a bare
            // lower trunk (conifers).
            static const BroadleafShape kBroad[3] = { { 1.0f, 1.25f, 1.3f, 2 }, { 1.35f, 0.65f, 0.7f, 0 }, { 1.05f, 1.1f, 1.0f, 1 } };
            static const ConiferShape kCone[2] = { { 1.2f, 0.75f, -1 }, { 0.9f, 1.1f, 2 } };
            const int sp = (int)SpeciesOfVariant(v), cls = SizeClassOfVariant(v);
            if (WoodSpeciesInfo((WoodSpecies)sp).conifer) {
                ConiferShape c = kCone[sp - (int)WoodSpecies::Pine];
                static const float kH[4] = { 0.5f, 1.0f, 1.35f, 1.8f }, kW[4] = { 0.55f, 1.0f, 1.2f, 1.35f }, kT[4] = { 0.6f, 1.0f, 1.4f, 2.1f }, kK[4] = { 0.6f, 1.0f, 1.2f, 1.4f };
                static const int kL[4] = { -3, 0, 2, 5 };
                c.h *= kH[cls]; c.width *= kW[cls]; c.trunk = kT[cls]; c.layers += kL[cls]; c.k = kK[cls]; c.giant = cls == 3;
                Conifer(m, 2000u + v * 7919u, c);
            } else {
                BroadleafShape b = kBroad[sp];
                static const float kH[4] = { 0.55f, 1.0f, 1.3f, 1.75f }, kC[4] = { 0.55f, 1.0f, 1.25f, 1.5f }, kT[4] = { 0.5f, 1.0f, 1.45f, 2.2f }, kR[4] = { 0.6f, 1.0f, 1.25f, 1.6f };
                static const int kB[4] = { -2, 0, 2, 5 };
                b.h *= kH[cls]; b.crown *= kC[cls]; b.trunk *= kT[cls]; b.reach = kR[cls]; b.branches += kB[cls]; b.giant = cls == 3;
                Broadleaf(m, 1000u + v * 7919u, b);
            }
        }
        else if (v == kStoneVariant) Stones(m, 4000u + v * 7919u);
        else Bush(m, 3000u + v * 7919u);
        m.vertexCount = (int)(m.pos.size() / 3);
        if (v == kStoneVariant) gVariants[1][mesh] = m; // already a few prisms (and no crown for FarFrom to fit)
        else FarFrom(m, v, gVariants[1][mesh]);
        gVariants[1][mesh].vertexCount = (int)(gVariants[1][mesh].pos.size() / 3);
        BakeFoliageAo(gVariants[0][mesh], v == kStoneVariant);
        BakeFoliageAo(gVariants[1][mesh], v == kStoneVariant);
    }
}
// By mesh id (TreeMeshOf, #1995).
const VariantMesh& Variant(int mesh, int lod = 0) { std::call_once(gVariantsOnce, BuildVariants); return gVariants[lod][mesh]; }
// Each variant's extent at scale 1 (#260 picking): height range along up and horizontal reach from the axis.
struct VariantBounds { float minY = 0.0f, maxY = 1.0f, reach = 0.5f; };
static const VariantBounds& BoundsOf(int v) {
    static VariantBounds b[kTreeVariants];
    static std::once_flag once;
    std::call_once(once, [] {
        for (int k = 0; k < kTreeVariants; ++k) { // over all its forms (#1995)
            VariantBounds vb{ 1e9f, -1e9f, 0.0f };
            for (int form = 0; form < kTreeForms; ++form) {
                const VariantMesh& m = Variant(TreeMeshOf(k, form), 0);
                for (size_t i = 0; i + 2 < m.pos.size(); i += 3) {
                    vb.minY = std::min(vb.minY, m.pos[i + 1]); vb.maxY = std::max(vb.maxY, m.pos[i + 1]);
                    vb.reach = std::max(vb.reach, sqrtf(m.pos[i] * m.pos[i] + m.pos[i + 2] * m.pos[i + 2]));
                }
            }
            if (vb.minY <= vb.maxY) b[k] = vb;
        }
    });
    return b[std::clamp(v, 0, kTreeVariants - 1)];
}

void EnsureVariantGpu() {
    for (int lod = 0; lod < 2; ++lod)
        for (int mesh = 0; mesh < kTreeMeshes; ++mesh) {
            VariantMesh& m = const_cast<VariantMesh&>(Variant(mesh, lod));
            if (m.vboPos) continue;
            GpuMemTag gm265_("trees"); // VRAM accounting (#302)
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
        case Biome::Jungle:   return { 1.5f, 1.4f, 0.8f, false, 1.35f, 0.34f }; // (#893: stones 0.5 -> 0.8)
        case Biome::Savanna:  return { 0.08f, 0.45f, 1.0f, false, 1.15f };
        case Biome::Desert:   return { 0.0f, 0.06f, 2.0f, false, 1.0f };
        case Biome::Badlands: return { 0.02f, 0.12f, 2.5f, false, 0.8f };
        case Biome::Forest:   return { 1.0f, 1.0f, 1.0f, false, 1.0f };
        case Biome::Meadow:   return { 0.07f, 1.3f, 0.8f, false, 1.1f };
        case Biome::Swamp:    return { 0.45f, 1.5f, 0.8f, false, 0.9f };          // (#893: stones 0.3 -> 0.8)
        case Biome::Taiga:    return { 1.35f, 0.6f, 1.0f, true, 1.05f, 0.26f };
        case Biome::Tundra:   return { 0.04f, 0.35f, 1.6f, true, 0.6f };
        case Biome::Alpine:   return { 0.25f, 0.5f, 1.5f, true, 0.8f };
        default:              return { 0.0f, 0.0f, 0.0f, true, 1.0f };   // ocean, polar
    }
}
// A biome variant's flora over its type's (#1978, biomes.txt): density and size multipliers, its own canopy floor.
static_assert(BiomeVariantDef::kSpecies == (int)WoodSpecies::Count, "biomes.txt species names follow WoodSpecies");
static BiomeFlora FloraOf(Biome b, const BiomeVariantDef* v) {
    BiomeFlora f = FloraOf(b);
    if (!v) return f;
    f.trees *= v->trees; f.bushes *= v->bushes; f.stones *= v->stones; f.size *= v->size;
    if (v->cover >= 0.0f) f.cover = v->cover;
    return f;
}
// The densest any variant allows: BuildTile's cheap first cut, so it never rejects a plant a variant would grow.
struct FloraMax { float trees = 0.0f, bushes = 0.0f, stones = 0.0f, cover = 0.0f; };
static const FloraMax& FloraMaxima() {
    static const FloraMax m = [] {
        FloraMax r;
        const BiomeVariantTable& t = BiomeVariants();
        for (const BiomeVariantDef& d : t.defs) {
            const BiomeFlora f = FloraOf(d.type, &d);
            r.trees = std::max(r.trees, f.trees); r.bushes = std::max(r.bushes, f.bushes);
            r.stones = std::max(r.stones, f.stones); r.cover = std::max(r.cover, f.cover);
        }
        const BiomeFlora island = FloraOf(Biome::Forest); // floating islands grow the plain forest
        r.trees = std::max(r.trees, island.trees); r.bushes = std::max(r.bushes, island.bushes); r.stones = std::max(r.stones, island.stones);
        return r;
    }();
    return m;
}

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
            const FloraMax& fmax = FloraMaxima();
            const float pTreeMax = std::max(TreeChanceFromForest(forest, fmax.trees), fmax.cover);
            const float pBushBase = 0.05f + 0.3f * Smooth(0.35f, 0.6f, forest);
            const float pStoneBase = 0.025f;                              // loose stones: the first stone tools' stone
            if (roll >= pTreeMax + pBushBase * fmax.bushes + pStoneBase * fmax.stones) continue;
            const float lat = fabsf(dir.y);
            if (lat > 0.86f) continue;                                  // polar: nothing grows
            // Floating islands (Plan 5 C) have woods on their tops: the column's ground is the island's top.
            float islandTop = 0.0f;
            const bool onIsland = f.islandTopAt(dir, islandTop);
            const float e = onIsland ? islandTop - R : f.elevation(dir);
            // The variant (#1978), read ~100 m off at random per plant so neighbouring variants' flora cross-fade over
            // ~200 m instead of meeting on a line (its own hash: the cell's draws below stay as they were).
            const BiomeVariantDef* vdef = nullptr;
            Biome biome = Biome::Forest;
            if (!onIsland) {
                uint32_t sv = Hash(s ^ 0x19780bu);
                const float ox = Unit(sv) * 2.0f - 1.0f, oy = Unit(sv) * 2.0f - 1.0f, k = 100.0f / R;
                const Vector3 u1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const Vector3 u2 = Vector3CrossProduct(u1, dir);
                const Vector3 uDir = Vector3Normalize(Vector3Add(dir, Vector3Add(Vector3Scale(u1, ox * k), Vector3Scale(u2, oy * k))));
                const int vid = f.variantAt(dir, e, uDir);
                const BiomeVariantTable& vt = BiomeVariants();
                if (vid >= 0 && vid < (int)vt.defs.size()) { vdef = &vt.defs[vid]; biome = vdef->type; }
                else biome = f.biomeAt(dir, e);
            }
            const BiomeFlora flora = FloraOf(biome, vdef);
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
            // Never under inland water (#189): a pond the fluid simulation took over, a lake or river bed - asked of the
            // generated water itself, not through the bed's material (the strata, #587, put trees in a taken-over pond).
            if (f.inlandWaterNear(dir, 0.0f)) {
                const float w = f.waterElevation(dir, 0.0f);
                if (w != PlanetField::kNoWater && Vector3Length(ground) - R < w + 0.2f) continue;
            }
            const VoxelMaterial m = f.materialLod(Vector3Subtract(ground, Vector3Scale(dir, 0.3f)), 1.0f, 0.0f, normal);
            // Trees and bushes root in turf or soil; loose stones lie on any natural ground but snow and ice (#893: the swamp's clay
            // and mud and the desert's sand had none, though their biomes call for them; the mountaintops stay bare by height).
            if (stone ? (m == VoxelMaterial::Snow || m == VoxelMaterial::Ice || m == VoxelMaterial::Ore) : (m != VoxelMaterial::Grass && m != VoxelMaterial::Regolith)) continue;
            // Solid ground underfoot and open air above (caves, ravines, edits).
            if (!solid(Vector3Subtract(ground, Vector3Scale(dir, 1.0f))) || !solid(Vector3Subtract(ground, Vector3Scale(dir, 4.0f)))) continue;
            if (solid(Vector3Add(ground, Vector3Scale(dir, 1.5f))) || (tree && solid(Vector3Add(ground, Vector3Scale(dir, 6.0f))))) continue;
            int v;
            float scale;
            if (tree) {
                // Broadleaf woods (oak, birch, maple) in the temperate band; conifer woods (pine, spruce,
                // birch) from ~40 degrees and up high, mixing at the edge. Each area (a ~400 m noise
                // patch) has a dominant species: 60% of its trees, the other two 20% each.
                const bool high = e > 550.0f + 300.0f * (h - 0.5f);
                const bool cold = flora.cold || high; // conifers in cold biomes and up high (Plan 5 B)
                static const WoodSpecies kWarm[3] = { WoodSpecies::Oak, WoodSpecies::Maple, WoodSpecies::Birch };
                static const WoodSpecies kCold[3] = { WoodSpecies::Pine, WoodSpecies::Spruce, WoodSpecies::Birch };
                const float area = ValueNoise(Vector3Scale(Vector3Scale(dir, R), 1.0f / 400.0f), seed ^ 0x7a3c19u);
                const float pick = Unit(s);
                WoodSpecies sp;
                if (vdef && vdef->hasSpecies() && (flora.cold || !high)) {
                    // The variant's own mix (#1978): its weights, with the area's dominant species (one of the listed,
                    // by the same ~400 m patch) counted double, so stands still form inside the mix.
                    int listed[BiomeVariantDef::kSpecies], n = 0;
                    for (int i = 0; i < BiomeVariantDef::kSpecies; ++i) if (vdef->species[i] > 0.0f) listed[n++] = i;
                    const int dom = listed[std::min(n - 1, (int)(Smooth(0.3f, 0.7f, area) * n))];
                    float w[BiomeVariantDef::kSpecies], total = 0.0f;
                    for (int i = 0; i < BiomeVariantDef::kSpecies; ++i) total += (w[i] = vdef->species[i] * (i == dom && n > 1 ? 2.0f : 1.0f));
                    float at = pick * total;
                    int k = dom;
                    for (int i = 0; i < BiomeVariantDef::kSpecies; ++i) if (w[i] > 0.0f) { k = i; if ((at -= w[i]) < 0.0f) break; }
                    sp = (WoodSpecies)k;
                } else {
                    const int dom = std::min(2, (int)(Smooth(0.3f, 0.7f, area) * 3.0f));
                    const int k = pick < 0.6f ? dom : pick < 0.8f ? (dom + 1) % 3 : (dom + 2) % 3;
                    sp = (cold ? kCold : kWarm)[k];
                }
                // Size class (#274): from the stand's age (a ~250 m noise patch) and the forest's density. Young
                // stands and open ground grow small trees; old dense forest holds the large ones and, rarely, giants.
                const float age = ValueNoise(Vector3Scale(Vector3Scale(dir, R), 1.0f / 250.0f), seed ^ 0x51ed27u);
                const float old = Smooth(0.45f, 0.85f, age) * Smooth(0.35f, 0.8f, forest);
                const float pc = Unit(s);
                const float pXL = 0.12f * old, pL = 0.08f + 0.3f * old, pS = 0.12f + 0.3f * (1.0f - Smooth(0.2f, 0.7f, age));
                const int cls = pc < pXL ? 3 : pc < pXL + pL ? 2 : pc < pXL + pL + pS ? 0 : 1;
                v = (int)sp * kTreeClasses + cls;
                scale = (0.8f + 0.35f * sz + 0.15f * forest) * flora.size; // the class carries the size; a little spread within it
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
            // #1995: one of the variant's forms, by its own hash (the cell's draws stay as they were); stones have one.
            const int form = TreeFormsOn() && v != kStoneVariant ? (int)(Hash(s ^ 0x1995a5u) % (uint32_t)kTreeForms) : 0;
            const int mesh = TreeMeshOf(v, form);
            out.inst[mesh].insert(out.inst[mesh].end(), { p.x, p.y, p.z, scale, dir.x * code, dir.y * code, dir.z * code, h });
            out.ids[mesh].push_back(id);
            const VariantMesh& vm = Variant(mesh);
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
        for (int slot = 0; slot < 2 * kTreeMeshes; ++slot) {
            if (vao[s][slot]) rlUnloadVertexArray(vao[s][slot]);
            if (vbo[s][slot]) rlUnloadVertexBuffer(vbo[s][slot]);
        }
    for (int slot = 0; slot < 2 * kTreeMeshes; ++slot) {
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
    JobSystem::tagNext(JobTag::Trees);
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
    // Other interest centres (#318: other players on a server): tiles within kTreeInterestRange of them stay loaded too, for
    // trunk collision and plant picking. Recomputed when the camera or any centre moves 16 m, or the set changes.
    bool interestMoved = interest.size() != lastInterest.size();
    for (size_t k = 0; !interestMoved && k < interest.size(); ++k) interestMoved = Vector3Distance(interest[k], lastInterest[k]) > 16.0f;
    const bool recompute = Vector3Distance(camLocal, lastWantCam) > 16.0f || interestMoved;
    if (recompute) { lastWantCam = camLocal; lastInterest = interest; }
    auto want = [&](Vector3 dir, float within) {
        const int N = TilesPerFace(radius);
        const float da = 0.5f * PI / N;
        const int K = (int)ceilf(within / kTreeTileSize) + 1;
        for (int f = 0; f < 6; ++f) {
            const float dn = Vector3DotProduct(dir, kN[f]);
            if (dn < 0.2f) continue;
            const float a = std::clamp(atanf(Vector3DotProduct(dir, kU[f]) / dn), -0.25f * PI, 0.25f * PI);
            const float b = std::clamp(atanf(Vector3DotProduct(dir, kV[f]) / dn), -0.25f * PI, 0.25f * PI);
            const int ic = (int)floorf((a + 0.25f * PI) / da), jc = (int)floorf((b + 0.25f * PI) / da);
            for (int i = std::max(0, ic - K); i <= std::min(N - 1, ic + K); ++i)
                for (int j = std::max(0, jc - K); j <= std::min(N - 1, jc + K); ++j) {
                    const TreeTileKey key{ f, i, j };
                    const float arc = acosf(std::clamp(Vector3DotProduct(TreeTileCentreDir(key, radius), dir), -1.0f, 1.0f)) * radius;
                    if (arc > within) continue;
                    auto it = tiles.find(key);
                    if (it == tiles.end()) {
                        TreeTile& t = tiles[key];
                        request(key, t, arc);
                    }
                }
        }
    };
    if (recompute && altitude < kTreeDrawDistance + 2000.0f) want(camDir, range);
    if (recompute) for (const Vector3& c : interest) want(Vector3Normalize(c), kTreeInterestRange);
    // Drop tiles well outside every range (hysteresis so the edge doesn't churn).
    for (auto it = tiles.begin(); recompute && it != tiles.end();) {
        const Vector3 tileDir = TreeTileCentreDir(it->first, radius);
        const float arc = acosf(std::clamp(Vector3DotProduct(tileDir, camDir), -1.0f, 1.0f)) * radius;
        bool keep = !(arc > range * 1.25f || altitude > kTreeDrawDistance + 2500.0f);
        for (const Vector3& c : interest)
            keep = keep || acosf(std::clamp(Vector3DotProduct(tileDir, Vector3Normalize(c)), -1.0f, 1.0f)) * radius < kTreeInterestRange * 1.5f;
        if (!keep) { it = tiles.erase(it); batchDirty = true; }
        else ++it;
    }
    const double now = GetTime();
    const bool moved = Vector3Distance(camLocal, Vector3{ (float)anchor.x, (float)anchor.y, (float)anchor.z }) > 40.0f; // keep the near/far split current
    if ((batchDirty || moved) && now - lastBatch > 0.2) rebuildBatches(camLocal);
}

int TreeSystem::readyTilesNear(Vector3 local, float range) const {
    const Vector3 dir = Vector3Normalize(local);
    int n = 0;
    for (const auto& kv : tiles)
        if (kv.second.ready && acosf(std::clamp(Vector3DotProduct(TreeTileCentreDir(kv.first, radius), dir), -1.0f, 1.0f)) * radius < range) ++n;
    return n;
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
    for (int mesh = 0; mesh < kTreeMeshes; ++mesh) {
        const int v = VariantOfTreeMesh(mesh);
        const float reach = v < kFirstBush ? kTreeDrawDistance : kBushDrawDistance;
        for (auto& s : buf) { s[0].clear(); s[1].clear(); }
        for (auto& kv : tiles) {
            const TreeTile& t = kv.second;
            if (!t.ready || t.inst[mesh].empty() || Vector3Distance(t.centre, camLocal) > reach + kTreeTileSize) continue;
            for (size_t k = 0; k + 7 < t.inst[mesh].size(); k += 8) {
                const float x = (float)(t.inst[mesh][k] - anchor.x), y = (float)(t.inst[mesh][k + 1] - anchor.y), z = (float)(t.inst[mesh][k + 2] - anchor.z);
                const float d2 = x * x + y * y + z * z;
                if (d2 > reach * reach) continue;
                const int lod = d2 < kTreeNearDistance * kTreeNearDistance ? 0 : 1;
                const float inst[8] = { x, y, z, t.inst[mesh][k + 3], t.inst[mesh][k + 4], t.inst[mesh][k + 5], t.inst[mesh][k + 6], t.inst[mesh][k + 7] };
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
                const int slot = mesh + lod * kTreeMeshes;
                unsigned int& va = vao[s][slot];
                unsigned int& vb = vbo[s][slot];
                if (va) { rlUnloadVertexArray(va); va = 0; }
                if (vb) { rlUnloadVertexBuffer(vb); vb = 0; }
                count[s][slot] = (int)(buf[s][lod].size() / 8);
                if (!count[s][slot]) continue;
                if (s == 0) drawnInstances += count[s][slot];
                const VariantMesh& m = Variant(mesh, lod);
                va = rlLoadVertexArray();
                rlEnableVertexArray(va);
                rlEnableVertexBuffer(m.vboPos); rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(0);
                rlEnableVertexBuffer(m.vboNrm); rlSetVertexAttribute(1, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(1);
                rlEnableVertexBuffer(m.vboCol); rlSetVertexAttribute(2, 4, RL_UNSIGNED_BYTE, true, 0, 0); rlEnableVertexAttribute(2);
                GpuMemTag gm667_("trees"); // VRAM accounting (#302)
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
                    GpuMemTag gm682_("trees"); // VRAM accounting (#302)
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
        constexpr int kCmds = 2 * kTreeMeshes;
        GLuint cmds[kCmds][4];
        for (int slot = 0; slot < kCmds; ++slot) {
            const int mesh = slot % kTreeMeshes, lod = slot / kTreeMeshes;
            cmds[slot][0] = (GLuint)Variant(mesh, lod).vertexCount; cmds[slot][1] = 0; cmds[slot][2] = 0; cmds[slot][3] = 0;
        }
        GpuMemTag gm728_("trees"); // VRAM accounting (#302)
        if (!cmdBuf) { glGenBuffers(1, &cmdBuf); glBindBuffer(GL_DRAW_INDIRECT_BUFFER, cmdBuf); glBufferData(GL_DRAW_INDIRECT_BUFFER, sizeof(cmds), nullptr, GL_DYNAMIC_DRAW); }
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, cmdBuf);
        glBufferSubData(GL_DRAW_INDIRECT_BUFFER, 0, sizeof(cmds), cmds);
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
        glUseProgram(gTreeCull.id);
        glUniformMatrix4fv(gTcModel, 1, GL_FALSE, MatrixToFloat(model));
        glUniform4fv(gTcPlanes, 6, (const float*)ctx.frustum->planes);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, cmdBuf);
        for (int slot = 0; slot < kCmds; ++slot) {
            const int mesh = slot % kTreeMeshes, lod = slot / kTreeMeshes;
            if (!count[0][slot] || !culledVbo[slot] || (onlyLod >= 0 && lod != onlyLod)) continue;
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, vbo[0][slot]);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, culledVbo[slot]);
            glUniform1ui(gTcN, (GLuint)count[0][slot]);
            glUniform1ui(gTcCmd, (GLuint)slot);
            // #745: the far-forest boost grows far trees up to 1.7x in trees.vs; pad the bound to match (RAYLIBGAME_FAR_FOREST=0: off).
            static const float farPad = [] { const char* e = getenv("RAYLIBGAME_FAR_FOREST"); return e && e[0] == '0' ? 1.0f : 1.7f; }(); // on by default (user, 2026-11-19); up to 1.7x
            constexpr float kVariationPad = 1.25f; // #1979: trees.vs stretches a tree up to 1.2x and leans it 4 degrees
            glUniform1f(gTcRadius, Variant(mesh, lod).boundRadius * farPad * kVariationPad);
            glDispatchCompute((GLuint)((count[0][slot] + 63) / 64), 1, 1);
        }
        glMemoryBarrier(GL_COMMAND_BARRIER_BIT | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT);
        glUseProgram(ctx.treeShader->id);
    }
    for (int slot = 0; slot < 2 * kTreeMeshes; ++slot) {
        if (!count[set][slot]) continue;
        const int mesh = slot % kTreeMeshes, lod = slot / kTreeMeshes, v = VariantOfTreeMesh(mesh);
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
        {   // #1995: which texture (trees.fs foliageKind); locations per program (the shader hot-reloads)
            static unsigned int prog = 0;
            static GLint kindLoc = -1, offLoc = -1;
            if (prog != ctx.treeShader->id) {
                prog = ctx.treeShader->id;
                kindLoc = glGetUniformLocation(prog, "foliageKind");
                offLoc = glGetUniformLocation(prog, "foliageTexOff");
            }
            static const float texOff = getenv("RAYLIBGAME_FOLIAGE_TEX") && !EnvFlag("RAYLIBGAME_FOLIAGE_TEX") ? 1.0f : 0.0f;
            const float kind = v < kFirstBush ? (float)(int)SpeciesOfVariant(v) : v == kStoneVariant ? 6.0f : 5.0f;
            if (kindLoc >= 0) glUniform1f(kindLoc, kind);
            if (offLoc >= 0) glUniform1f(offLoc, texOff);
        }
        if (gpuCull && culledVao[slot]) {
            glBindVertexArray(culledVao[slot]);
            glBindBuffer(GL_DRAW_INDIRECT_BUFFER, cmdBuf);
            glDrawArraysIndirect(GL_TRIANGLES, (const void*)(size_t)(slot * 4 * sizeof(GLuint)));
            glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
            glBindVertexArray(0);
        } else {
            rlEnableVertexArray(vao[set][slot]);
            rlDrawVertexArrayInstanced(0, Variant(mesh, lod).vertexCount, count[set][slot]);
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

// (#893 tests) The scatter around `dir` as generated (no edits, removals or clearings), counted within `range` m of the
// ground point there: trees, bushes and loose stones. Builds the tiles itself (the tree system need not have them).
void Trees_CountScatter(const PlanetField& f, uint32_t seed, Vector3 dir, float range, int& trees, int& bushes, int& stones) {
    trees = bushes = stones = 0;
    const Vector3 d = Vector3Normalize(dir);
    const Vector3 at = Vector3Scale(d, f.radius + f.elevation(d));
    TreeTileKey keys[9]; int n = 0;
    TilesAround(d, f.radius, keys, n);
    const TerrainDeltas none;
    const PlantRemovals noRemovals;
    const std::vector<Vector4> noClearings;
    for (int k = 0; k < n; ++k) {
        TreeTile t;
        BuildTile(f, none, noRemovals, noClearings, 0.0, seed, keys[k], t);
        for (const Plant& p : t.plants) {
            if (Vector3Distance(p.base, at) > range) continue;
            if (p.stone()) ++stones; else if (p.tree()) ++trees; else ++bushes;
        }
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

void TreeSystem::plantsNear(Vector3 local, float range, std::vector<Plant>& out) const {
    const float r = Vector3Length(local);
    if (r < 1.0f) return;
    TreeTileKey keys[9];
    int n = 0;
    TilesAround(Vector3Scale(local, 1.0f / r), radius, keys, n);
    for (int m = 0; m < n; ++m) {
        auto it = tiles.find(keys[m]);
        if (it == tiles.end() || !it->second.ready) continue;
        for (const Plant& p : it->second.plants)
            if (Vector3Distance(p.base, local) < range) out.push_back(p);
    }
}

bool TreeSystem::findPlant(uint64_t id, Plant& out) const {
    auto it = tiles.find(TileOfPlant(id));
    if (it == tiles.end() || !it->second.ready) return false;
    for (const Plant& p : it->second.plants) if (p.id == id) { out = p; return true; }
    return false;
}

void TreeSystem::remove(uint64_t id, double when) {
    auto next = std::make_shared<PlantRemovals>(*removed);
    (*next)[id] = when >= 0.0 ? when : WorldTime();
    removed = next;
    // Gone at once from the loaded tile (the next rebuild agrees), and its tile starts watching
    // for regrowth.
    auto it = tiles.find(TileOfPlant(id));
    if (it == tiles.end()) return;
    TreeTile& t = it->second;
    for (int v = 0; v < kTreeMeshes; ++v) // by mesh (#1995)
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

void TreeSystem::mergeRemovals(const std::vector<std::pair<uint64_t, double>>& list) {
    if (list.size() <= 8) { for (const auto& [id, t] : list) remove(id, t); return; }
    PlantRemovals next(*removed);
    for (const auto& [id, t] : list) next[id] = t;
    setRemovals(next);
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
                if (!p.tree()) {
                    // Bushes and stones (#260): a vertical cylinder fitted to the variant's mesh (its height range and
                    // horizontal reach, at the plant's scale). The fixed capsule picked stones in the air above them and
                    // missed bush edges and tops.
                    const VariantBounds& vb = BoundsOf(p.variant);
                    const float R = vb.reach * p.scale * 0.85f, y0 = vb.minY * p.scale, y1 = vb.maxY * p.scale;
                    const Vector3 w = Vector3Subtract(o, p.base);
                    const float oy = Vector3DotProduct(w, p.up), dy = Vector3DotProduct(d, p.up);
                    const Vector3 oh = Vector3Subtract(w, Vector3Scale(p.up, oy)), dh = Vector3Subtract(d, Vector3Scale(p.up, dy));
                    float tIn = 0.0f, tOut = maxDistance;
                    const float a = Vector3DotProduct(dh, dh), b = Vector3DotProduct(oh, dh), c = Vector3DotProduct(oh, oh) - R * R;
                    if (a < 1e-8f) { if (c > 0.0f) continue; }                 // straight up/down: inside the circle or not
                    else {
                        const float disc = b * b - a * c;
                        if (disc < 0.0f) continue;
                        const float sq = sqrtf(disc);
                        tIn = std::max(tIn, (-b - sq) / a); tOut = std::min(tOut, (-b + sq) / a);
                    }
                    if (fabsf(dy) < 1e-8f) { if (oy < y0 || oy > y1) continue; }
                    else {
                        float ta = (y0 - oy) / dy, tb = (y1 - oy) / dy;
                        if (ta > tb) std::swap(ta, tb);
                        tIn = std::max(tIn, ta); tOut = std::min(tOut, tb);
                    }
                    if (tIn > tOut) continue;
                    if (tIn < best.distance) { best.distance = tIn; best.plant = p; found = true; }
                    continue;
                }
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
    BiomeFlora fl{};
    if (field) { // the variant's flora (#1978), as BuildTile (without its border dither)
        const int vid = field->variantAt(dir, elevation);
        const BiomeVariantTable& vt = BiomeVariants();
        fl = vid >= 0 && vid < (int)vt.defs.size() ? FloraOf(vt.defs[vid].type, &vt.defs[vid]) : FloraOf(field->biomeAt(dir, elevation));
    }
    if (field && fl.cover > 0.0f) return std::max(Smooth(0.04f, 0.3f, fl.cover), Smooth(0.04f, 0.3f, TreeChanceFromForest(Forest(Vector3Scale(dir, radius), seed), fl.trees))) * edge;
    const float trees = field ? fl.trees : (lat < 0.6f ? 1.3f - 0.3f * Smooth(0.2f, 0.4f, lat) : 1.0f - 0.55f * (lat - 0.6f) / 0.26f);
    return Smooth(0.04f, 0.3f, TreeChanceFromForest(Forest(Vector3Scale(dir, radius), seed), trees)) * edge; // forest presence
}

// (#1856 M20) Whether a clearing reaches a tile: built or not - one still building has no centre yet (its key's is used),
// and its build captured the clearings as they were (it must be asked again, or its trees grow inside the new one).
bool TreeClearingReaches(const TreeTileKey& key, const TreeTile& tile, float planetRadius, Vector3 c, float r) {
    const bool hasCentre = tile.centre.x != 0.0f || tile.centre.y != 0.0f || tile.centre.z != 0.0f;
    const Vector3 centre = hasCentre ? tile.centre : Vector3Scale(TreeTileCentreDir(key, planetRadius), planetRadius);
    return Vector3Distance(centre, c) < r + kTreeTileSize;
}

void TreeSystem::addClearing(Vector3 c, float r) {
    auto next = std::make_shared<std::vector<Vector4>>(*clearings);
    next->push_back({ c.x, c.y, c.z, r });
    clearings = next;
    for (auto& kv : tiles)
        if (TreeClearingReaches(kv.first, kv.second, radius, c, r)) request(kv.first, kv.second, 0.0f);
}

bool TreeSystem::nearestPlant(Vector3 local, float range, bool tree, Plant& out, bool stone, int sizeClass) const {
    float best = range;
    bool found = false;
    for (const auto& kv : tiles) {
        if (!kv.second.ready || Vector3Distance(kv.second.centre, local) > range + kTreeTileSize) continue;
        for (const Plant& p : kv.second.plants) {
            if (p.tree() != tree || p.stone() != stone || (sizeClass >= 0 && (!p.tree() || SizeClassOfVariant(p.variant) != sizeClass))) continue;
            const float d = Vector3Distance(p.base, local);
            if (d < best) { best = d; out = p; found = true; }
        }
    }
    return found;
}

void TreeSystem::speciesCounts(Vector3 local, float range, int out[(int)WoodSpecies::Count], int* classes) const {
    for (int s = 0; s < (int)WoodSpecies::Count; ++s) out[s] = 0;
    if (classes) for (int c = 0; c < kTreeClasses; ++c) classes[c] = 0;
    for (const auto& kv : tiles) {
        const TreeTile& t = kv.second;
        if (!t.ready) continue;
        for (int mesh = 0; mesh < kFirstBush * kTreeForms; ++mesh) { const int v = VariantOfTreeMesh(mesh);
            for (size_t k = 0; k + 7 < t.inst[mesh].size(); k += 8)
                if (Vector3Distance({ t.inst[mesh][k], t.inst[mesh][k + 1], t.inst[mesh][k + 2] }, local) < range) { ++out[(int)SpeciesOfVariant(v)]; if (classes) ++classes[SizeClassOfVariant(v)]; } }
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
// (#1856 M20) A clearing reaches a tile still building (no centre yet: its key's is used), not only built ones.
TEST(Trees, ClearingReachesATileStillBuilding) {
    const float R = 60000.0f;
    const TreeTileKey key{ 0, 5, 5 };
    const Vector3 at = Vector3Scale(TreeTileCentreDir(key, R), R);
    TreeTile building; // ready false, centre zero: its first build is pending
    EXPECT_TRUE(TreeClearingReaches(key, building, R, Vector3Add(at, { 10.0f, 0.0f, 0.0f }), 20.0f));
    EXPECT_FALSE(TreeClearingReaches(key, building, R, Vector3Scale(at, -1.0f), 20.0f)); // the far side of the planet
    TreeTile built; built.ready = true; built.centre = at;
    EXPECT_TRUE(TreeClearingReaches(key, built, R, at, 5.0f));
}

// (#1995) Every mesh is built; a variant's forms are different trees; each species grows with its size class (sapling <
// young < mature < old, mean over the forms); trees have a trunk to collide with and chop, bushes none.
TEST(Trees, FormsAndSizeClasses) {
    if (!TreeFormsOn()) GTEST_SKIP() << "RAYLIBGAME_TREE_FORMS=0";
    for (int mesh = 0; mesh < kTreeMeshes; ++mesh) {
        EXPECT_GT(Variant(mesh, 0).vertexCount, 0) << mesh;
        EXPECT_GT(Variant(mesh, 1).vertexCount, 0) << mesh;
    }
    for (int v = 0; v < kTreeVariants; ++v) {
        if (v == kStoneVariant) continue;
        const VariantMesh &a = Variant(TreeMeshOf(v, 0)), &b = Variant(TreeMeshOf(v, 1)), &c = Variant(TreeMeshOf(v, 2));
        EXPECT_TRUE(a.pos != b.pos && b.pos != c.pos && a.pos != c.pos) << "variant " << v << ": forms are clones";
        if (v < kFirstBush) { EXPECT_GT(a.trunkRadius, 0.0f) << v; EXPECT_GT(a.trunkHeight, 0.5f) << v; }
        else EXPECT_EQ(a.trunkRadius, 0.0f) << v;
    }
    for (int sp = 0; sp < (int)WoodSpecies::Count; ++sp) {
        float prev = 0.0f;
        for (int cls = 0; cls < kTreeClasses; ++cls) {
            float h = 0.0f;
            for (int form = 0; form < kTreeForms; ++form) h += Variant(TreeMeshOf(sp * kTreeClasses + cls, form)).trunkHeight / kTreeForms;
            EXPECT_GT(h, prev) << WoodSpeciesInfo((WoodSpecies)sp).name << " class " << cls;
            prev = h;
        }
    }
}
#endif
