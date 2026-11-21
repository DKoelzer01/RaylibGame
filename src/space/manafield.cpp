#include "manafield.h"

#include "SimplexNoise.h"
#include "strata.h"
#include "voxel.h" // Biome

#include <raymath.h>
#include <algorithm>
#include <cmath>

// All numbers here are first guesses (the magic doc marks them so), to tune with the overlay in game.
namespace {
    constexpr float kLineScale[kManaBases] = { 3500.0f, 3000.0f, 4000.0f, 5000.0f }; // m: a base's contour noise (lines a few km apart)
    constexpr float kRegionScale = 12000.0f; // m: the kinds' regions
    constexpr float kFloor = 0.2f;           // every kind's strength anywhere: faint lines everywhere, nodes need more

    uint32_t Hash(uint32_t a, uint32_t b, uint32_t c, uint32_t seed) {
        uint32_t h = a * 73856093u ^ b * 19349663u ^ c * 83492791u ^ seed * 2654435761u;
        h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
        return h;
    }
    float Offset(unsigned int seed, uint32_t salt) { return (float)(Hash(salt, 7u, 11u, seed) & 0xFFFFu) * 1.5f; } // 0..98 km
    float Smooth(float e0, float e1, float x) { const float t = std::clamp((x - e0) / (e1 - e0), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); }

    // A tangent basis at a unit direction.
    void Tangents(Vector3 d, Vector3& e, Vector3& n) {
        e = Vector3CrossProduct({ 0.0f, 1.0f, 0.0f }, d);
        if (Vector3Length(e) < 1e-4f) e = Vector3CrossProduct({ 1.0f, 0.0f, 0.0f }, d);
        e = Vector3Normalize(e);
        n = Vector3CrossProduct(d, e);
    }
    Vector3 Step(Vector3 d, Vector3 e, Vector3 n, float x, float y, float R) { // a direction x m east and y m north of d
        return Vector3Normalize(Vector3Add(d, Vector3Add(Vector3Scale(e, x / R), Vector3Scale(n, y / R))));
    }

    // The cube-sphere tiles (equal-angle): face 0-2 = +X,+Y,+Z, 3-5 = -X,-Y,-Z; a, b = angles in [-pi/4, pi/4].
    constexpr float kQuarter = 0.78539816f;
    Vector3 FaceDir(int face, float a, float b) {
        const int m = face % 3;
        float c[3];
        c[m] = face < 3 ? 1.0f : -1.0f;
        c[(m + 1) % 3] = tanf(a);
        c[(m + 2) % 3] = tanf(b);
        return Vector3Normalize({ c[0], c[1], c[2] });
    }
    void DirFace(Vector3 d, int& face, float& a, float& b) {
        const float c[3] = { d.x, d.y, d.z };
        int m = 0;
        if (fabsf(c[1]) > fabsf(c[m])) m = 1;
        if (fabsf(c[2]) > fabsf(c[m])) m = 2;
        face = m + (c[m] < 0.0f ? 3 : 0);
        const float inv = 1.0f / std::max(fabsf(c[m]), 1e-6f);
        a = atanf(c[(m + 1) % 3] * inv);
        b = atanf(c[(m + 2) % 3] * inv);
    }
}

struct ManaField::Noise {
    SimplexNoise line[kManaBases], region[kManaBases];
};

ManaField::ManaField(const ManaFieldParams& params) : p(params), noise(std::make_shared<Noise>()) {
    for (int k = 0; k < kManaBases; ++k) {
        noise->line[k] = SimplexNoise(1.0f / kLineScale[k], 1.0f, 2.0f, 0.45f, Offset(p.seed, 100u + k), Offset(p.seed, 200u + k), Offset(p.seed, 300u + k));
        noise->region[k] = SimplexNoise(1.0f / kRegionScale, 1.0f, 2.0f, 0.5f, Offset(p.seed, 400u + k), Offset(p.seed, 500u + k), Offset(p.seed, 600u + k));
        families.push_back({ k, k });
    }
    families.push_back({ kManaEarth, -1 }); // Earth also runs along the faults (#590)
    if (!Mana_Loaded()) Mana_Load();
    static const char* kBaseIds[kManaBases] = { "earth", "water", "fire", "air" };
    for (int k = 0; k < kManaBases; ++k) kBase[k] = Mana_Find(kBaseIds[k]);
    for (int a = 0; a < kManaBases; ++a)
        for (int b = 0; b < kManaBases; ++b) kindOfPair[a][b] = a == b ? kBase[a] : Mana_Combine(kBase[a], kBase[b]);
    kAether = Mana_Aether();
    kFrost = Mana_Find("frost"); kStorm = Mana_Find("storm"); kShard = Mana_Find("shard"); kRot = Mana_Find("rot");
}

float ManaField::strength(int base, Vector3 dir) const {
    if (base < 0 || base >= kManaBases) return 0.0f;
    const Vector3 s = Vector3Scale(dir, p.radius);
    float w = Smooth(0.45f, 0.75f, 0.5f + 0.5f * noise->region[base].fractal(2, s.x, s.y, s.z));
    switch (base) {
        case kManaEarth: { float region = 0.0f; Strata_FaultValue(s, p.seed, &region); w += 0.8f * region; break; }
        case kManaWater: if (p.karst) w += 0.8f * Strata_KarstWeight(s, p.seed); break;
        case kManaFire: if (p.volcanoes) { float vw = 0.0f; Strata_Volcano(s, p.radius, p.seed, &vw); w += 0.9f * vw; } break;
        default: break;
    }
    return std::clamp(kFloor + (1.0f - kFloor) * w, 0.0f, 1.0f);
}

float ManaField::familyValue(int fam, Vector3 s) const {
    const Family& f = families[fam];
    if (f.noise >= 0) return noise->line[f.noise].fractal(3, s.x, s.y, s.z);
    float region = 0.0f;
    const float v = Strata_FaultValue(s, p.seed, &region);
    return region >= 0.3f ? v : NAN; // outside the fault regions there is no fault line
}

float ManaField::familyDistance(int fam, Vector3 dir) const {
    constexpr float eps = 8.0f;
    Vector3 e, n; Tangents(dir, e, n);
    const float g0 = familyValue(fam, Vector3Scale(dir, p.radius));
    const float gx = familyValue(fam, Vector3Scale(Step(dir, e, n, eps, 0.0f, p.radius), p.radius));
    const float gy = familyValue(fam, Vector3Scale(Step(dir, e, n, 0.0f, eps, p.radius), p.radius));
    if (std::isnan(g0) || std::isnan(gx) || std::isnan(gy)) return 1e9f;
    const float grad = sqrtf((gx - g0) * (gx - g0) + (gy - g0) * (gy - g0)) / eps;
    return grad > 1e-9f ? fabsf(g0) / grad : 1e9f;
}

float ManaField::lineDistance(int base, Vector3 dir) const {
    float d = 1e9f;
    for (int f = 0; f < (int)families.size(); ++f) if (families[f].base == base) d = std::min(d, familyDistance(f, dir));
    return d;
}

float ManaField::lineFlow(int base, Vector3 dir) const {
    const float d = lineDistance(base, dir);
    if (d > 20.0f * kLineWidth) return 0.0f;
    const float q = d / kLineWidth, h = d / (5.0f * kLineWidth);
    return strength(base, dir) * (expf(-q * q) + 0.25f * expf(-h * h)) / 1.25f; // a core and a faint halo
}

int ManaField::tilesPerFace() const { return std::max(1, (int)ceilf(2.0f * kQuarter * p.radius / kTileSize)); }

const std::vector<ManaNode>& ManaField::tile(int face, int ti, int tj) const {
    const uint64_t key = ((uint64_t)face << 40) | ((uint64_t)(uint32_t)ti << 20) | (uint64_t)(uint32_t)tj;
    {
        std::lock_guard<std::mutex> lock(mu);
        auto it = tiles.find(key);
        if (it != tiles.end()) return it->second;
    }
    // A grid of ~100 m squares over the tile; a crossing of two families of different kinds is a node, owned by the
    // square its root lies in (so neighbouring squares and tiles never both keep it).
    const int n = tilesPerFace();
    const float span = 2.0f * kQuarter / n;
    const float a0 = -kQuarter + ti * span, b0 = -kQuarter + tj * span;
    const int M = std::max(4, (int)ceilf(kTileSize / 100.0f));
    const float cell = span / M;
    const int F = (int)families.size();
    std::vector<float> g((size_t)(M + 1) * (M + 1) * F);
    auto at = [&](int i, int j, int f) -> float& { return g[((size_t)j * (M + 1) + i) * F + f]; };
    for (int j = 0; j <= M; ++j)
        for (int i = 0; i <= M; ++i) {
            const Vector3 s = Vector3Scale(FaceDir(face, a0 + i * cell, b0 + j * cell), p.radius);
            for (int f = 0; f < F; ++f) at(i, j, f) = familyValue(f, s);
        }
    auto changes = [&](int i, int j, int f) {
        const float v[4] = { at(i, j, f), at(i + 1, j, f), at(i, j + 1, f), at(i + 1, j + 1, f) };
        bool neg = false, pos = false;
        for (float x : v) { if (std::isnan(x)) return false; neg |= x < 0.0f; pos |= x >= 0.0f; }
        return neg && pos;
    };
    std::vector<ManaNode> out;
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < M; ++i)
            for (int f1 = 0; f1 < F; ++f1) {
                if (!changes(i, j, f1)) continue;
                for (int f2 = f1 + 1; f2 < F; ++f2) {
                    if (families[f1].base == families[f2].base || !changes(i, j, f2)) continue;
                    // Newton on (g1, g2) = 0 in the face angles from the square's centre.
                    float a = a0 + (i + 0.5f) * cell, b = b0 + (j + 0.5f) * cell;
                    const float h = cell * 0.05f;
                    bool ok = true;
                    for (int it = 0; it < 6 && ok; ++it) {
                        auto val = [&](float x, float y, int f) { return familyValue(f, Vector3Scale(FaceDir(face, x, y), p.radius)); };
                        const float u1 = val(a, b, f1), u2 = val(a, b, f2);
                        const float j11 = (val(a + h, b, f1) - u1) / h, j12 = (val(a, b + h, f1) - u1) / h;
                        const float j21 = (val(a + h, b, f2) - u2) / h, j22 = (val(a, b + h, f2) - u2) / h;
                        const float det = j11 * j22 - j12 * j21;
                        if (!(fabsf(det) > 1e-12f)) { ok = false; break; } // (also catches NaN: off a fault region)
                        a -= (u1 * j22 - u2 * j12) / det;
                        b -= (u2 * j11 - u1 * j21) / det;
                    }
                    if (!ok || a < a0 + i * cell || a >= a0 + (i + 1) * cell || b < b0 + j * cell || b >= b0 + (j + 1) * cell) continue;
                    ManaNode node;
                    node.dir = FaceDir(face, a, b);
                    node.a = std::min(families[f1].base, families[f2].base);
                    node.b = std::max(families[f1].base, families[f2].base);
                    node.strength = std::min(strength(node.a, node.dir), strength(node.b, node.dir));
                    if (node.strength < kNodeMinStrength) continue;
                    node.kind = kindOfPair[node.a][node.b];
                    if (node.b == kManaFire || node.a == kManaFire) // Fire meeting frozen ground: Fire + Frost cancel into Aether
                        if (p.meanTempC && p.meanTempC(node.dir) < -5.0f) node.kind = kAether;
                    node.id = Hash((uint32_t)(face * 4096 + ti), (uint32_t)tj, (uint32_t)((j * M + i) * 64 + f1 * 8 + f2), p.seed);
                    out.push_back(node);
                }
            }
    std::lock_guard<std::mutex> lock(mu);
    return tiles.emplace(key, std::move(out)).first->second; // (a race builds the same tile twice: the first stays)
}

void ManaField::nodesNear(Vector3 dir, float range, std::vector<ManaNode>& out) const {
    dir = Vector3Normalize(dir);
    const int n = tilesPerFace();
    const float span = 2.0f * kQuarter / n;
    // The tiles under the centre and rings of points around it, closer together than a tile.
    std::vector<uint64_t> seen;
    Vector3 e, nn; Tangents(dir, e, nn);
    auto visit = [&](Vector3 d) {
        int face; float a, b; DirFace(d, face, a, b);
        const int ti = std::clamp((int)floorf((a + kQuarter) / span), 0, n - 1), tj = std::clamp((int)floorf((b + kQuarter) / span), 0, n - 1);
        const uint64_t key = ((uint64_t)face << 40) | ((uint64_t)ti << 20) | (uint64_t)tj;
        if (std::find(seen.begin(), seen.end(), key) != seen.end()) return;
        seen.push_back(key);
        const float cosRange = cosf(std::min(range / p.radius, 3.14159f));
        for (const ManaNode& node : tile(face, ti, tj)) if (Vector3DotProduct(node.dir, dir) >= cosRange) out.push_back(node);
    };
    visit(dir);
    const float ringStep = 0.4f * kTileSize;
    const int rings = std::max(1, (int)ceilf(range / ringStep));
    for (int ring = 1; ring <= rings; ++ring) {
        const float r = range * ring / rings;
        const int pts = std::max(16, (int)ceilf(6.2831853f * r / ringStep));
        for (int k = 0; k < pts; ++k) {
            const float t = k * 6.2831853f / pts;
            const float th = std::min(r / p.radius, 3.14159f); // along a great circle
            visit(Vector3Add(Vector3Scale(dir, cosf(th)), Vector3Scale(Vector3Add(Vector3Scale(e, cosf(t)), Vector3Scale(nn, sinf(t))), sinf(th))));
        }
    }
}

ManaSample ManaField::density(Vector3 dir, const ManaEnv& env) const {
    dir = Vector3Normalize(dir);
    ManaSample out;
    auto add = [&](int kind, float v) { if (kind >= 0 && kind < kManaMaxKinds) out.v[kind] += v; };
    float flow[kManaBases];
    for (int k = 0; k < kManaBases; ++k) flow[k] = lineFlow(k, dir);
    // What the place adds to the bases (doc 3: deep stone, open water, height and wind).
    flow[kManaEarth] *= 1.0f + 0.5f * std::clamp(env.depth / 50.0f, 0.0f, 1.0f);
    flow[kManaWater] *= env.water ? 1.3f : 1.0f;
    flow[kManaAir] *= std::clamp(0.4f + env.elevation / 1500.0f, 0.4f, 1.6f) * (1.0f + 0.04f * std::max(0.0f, env.wind - 4.0f))
                    * (1.0f + env.storm) * (env.depth > 5.0f ? 0.2f : 1.0f);
    for (int k = 0; k < kManaBases; ++k) add(kBase[k], flow[k]);
    // Nodes: a glow of the node's kind, and half as much of its two lines' kinds.
    std::vector<ManaNode> nodes;
    nodesNear(dir, 3.0f * kNodeReach, nodes);
    float nodeGlow = 0.0f;
    for (const ManaNode& node : nodes) {
        const float d = acosf(std::clamp(Vector3DotProduct(dir, node.dir), -1.0f, 1.0f)) * p.radius;
        const float w = 1.5f * node.strength * expf(-(d / kNodeReach) * (d / kNodeReach));
        if (node.kind >= 0 && Mana_Kind(node.kind).slice) add(node.kind, w);
        add(kBase[node.a], 0.5f * w);
        add(kBase[node.b], 0.5f * w);
        nodeGlow += w;
    }
    // The slice's other kinds (doc 3.1 "found at").
    const Biome biome = (Biome)env.biome;
    const bool coldBiome = biome == Biome::Tundra || biome == Biome::Taiga || biome == Biome::Polar;
    add(kFrost, std::clamp((2.0f - env.tempC) / 25.0f, 0.0f, 1.0f) * (1.0f + 0.5f * env.night) + (coldBiome ? 0.2f : 0.0f));
    add(kStorm, env.storm * (0.3f + 0.6f * (flow[kManaAir] + flow[kManaFire])) + (biome == Biome::Alpine ? 0.25f : 0.0f));
    add(kShard, biome == Biome::Desert || biome == Biome::Badlands ? 0.35f : 0.0f);
    add(kRot, biome == Biome::Swamp ? 0.5f : biome == Biome::Jungle ? 0.1f : 0.0f);
    add(kAether, 0.03f + 0.3f * nodeGlow); // a trace everywhere, and wherever lines cross
    for (float& v : out.v) v = std::clamp(v, 0.0f, 2.0f);
    return out;
}

void ManaField::lineSegments(Vector3 dir, float range, float step, float minStrength, std::vector<ManaSegment>& out) const {
    dir = Vector3Normalize(dir);
    Vector3 e, n; Tangents(dir, e, n);
    const int H = std::max(1, (int)ceilf(range / step)), W = 2 * H + 1;
    const int F = (int)families.size();
    std::vector<Vector3> dirs((size_t)W * W);
    std::vector<float> g((size_t)W * W * F);
    for (int j = 0; j < W; ++j)
        for (int i = 0; i < W; ++i) {
            const Vector3 d = Step(dir, e, n, (i - H) * step, (j - H) * step, p.radius);
            dirs[(size_t)j * W + i] = d;
            for (int f = 0; f < F; ++f) g[((size_t)j * W + i) * F + f] = familyValue(f, Vector3Scale(d, p.radius));
        }
    auto G = [&](int i, int j, int f) { return g[((size_t)j * W + i) * F + f]; };
    auto D = [&](int i, int j) { return dirs[(size_t)j * W + i]; };
    auto cross = [&](int i0, int j0, int i1, int j1, int f) { // where the contour crosses the edge
        const float v0 = G(i0, j0, f), v1 = G(i1, j1, f);
        const float t = v0 / (v0 - v1);
        return Vector3Normalize(Vector3Lerp(D(i0, j0), D(i1, j1), t));
    };
    for (int f = 0; f < F; ++f) {
        const int base = families[f].base;
        for (int j = 0; j + 1 < W; ++j)
            for (int i = 0; i + 1 < W; ++i) {
                const float v[4] = { G(i, j, f), G(i + 1, j, f), G(i + 1, j + 1, f), G(i, j + 1, f) };
                if (std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]) || std::isnan(v[3])) continue;
                const int mask = (v[0] < 0) | (v[1] < 0) << 1 | (v[2] < 0) << 2 | (v[3] < 0) << 3;
                if (mask == 0 || mask == 15) continue;
                // The square's edges 0 (bottom) .. 3 (left) and where the contour crosses each.
                Vector3 pt[4]; bool has[4] = {};
                const int ei[4][4] = { { i, j, i + 1, j }, { i + 1, j, i + 1, j + 1 }, { i + 1, j + 1, i, j + 1 }, { i, j + 1, i, j } };
                for (int k = 0; k < 4; ++k) if ((v[k] < 0) != (v[(k + 1) % 4] < 0)) { has[k] = true; pt[k] = cross(ei[k][0], ei[k][1], ei[k][2], ei[k][3], f); }
                std::vector<int> hit;
                for (int k = 0; k < 4; ++k) if (has[k]) hit.push_back(k);
                const float s = strength(base, D(i, j));
                if (s < minStrength) continue;
                for (size_t k = 0; k + 1 < hit.size(); k += 2) out.push_back({ pt[hit[k]], pt[hit[k + 1]], s, base });
            }
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1248): the field is seeded, its lines and nodes sit where the contours say, its density is bounded ----
#include <gtest/gtest.h>
#include <cstdio>

namespace {
    ManaFieldParams TestParams(unsigned int seed) { ManaFieldParams p; p.seed = seed; p.radius = 47000.0f; return p; }
    Vector3 Dir(float x, float y, float z) { return Vector3Normalize({ x, y, z }); }
}

TEST(ManaField, TheSameSeedGivesTheSameFieldAndAnotherSeedAnother) {
    ASSERT_TRUE(Mana_Load());
    ManaField a(TestParams(42)), b(TestParams(42)), c(TestParams(7));
    const Vector3 here = Dir(0.3f, 0.8f, 0.5f);
    std::vector<ManaNode> na, nb, nc;
    a.nodesNear(here, 3000.0f, na);
    b.nodesNear(here, 3000.0f, nb);
    c.nodesNear(here, 3000.0f, nc);
    ASSERT_EQ(na.size(), nb.size());
    for (size_t i = 0; i < na.size(); ++i) {
        EXPECT_EQ(na[i].id, nb[i].id);
        EXPECT_EQ(na[i].kind, nb[i].kind);
        EXPECT_TRUE(na[i].dir.x == nb[i].dir.x && na[i].dir.y == nb[i].dir.y && na[i].dir.z == nb[i].dir.z) << "bit-identical directions";
    }
    int differ = 0;
    for (int k = 0; k < 64; ++k) {
        const Vector3 d = Dir(sinf(k * 1.3f), cosf(k * 0.7f), sinf(k * 2.1f + 1.0f));
        for (int base = 0; base < kManaBases; ++base) {
            EXPECT_EQ(a.lineDistance(base, d), b.lineDistance(base, d));
            EXPECT_EQ(a.strength(base, d), b.strength(base, d));
            differ += fabsf(a.lineDistance(base, d) - c.lineDistance(base, d)) > 1.0f;
        }
    }
    EXPECT_GT(differ, 128) << "another seed moves the lines";
    bool sameNodes = na.size() == nc.size();
    for (size_t i = 0; sameNodes && i < na.size(); ++i) sameNodes = na[i].id == nc[i].id;
    EXPECT_FALSE(sameNodes && !na.empty());
}

TEST(ManaField, LinesAreWhereTheirContoursAre) {
    ASSERT_TRUE(Mana_Load());
    ManaField f(TestParams(42));
    std::vector<ManaSegment> segs;
    const Vector3 here = Dir(-0.4f, 0.5f, 0.77f);
    f.lineSegments(here, 6000.0f, 50.0f, 0.0f, segs);
    ASSERT_FALSE(segs.empty()) << "lines a few km apart: some within 6 km";
    int checked = 0, onLine = 0;
    for (size_t i = 0; i < segs.size(); i += std::max<size_t>(1, segs.size() / 200)) {
        const Vector3 mid = Vector3Normalize(Vector3Add(segs[i].a, segs[i].b));
        ++checked;
        onLine += f.lineDistance(segs[i].base, mid) < 10.0f;
        EXPECT_GT(f.lineFlow(segs[i].base, mid), 0.0f);
    }
    EXPECT_GT(onLine, checked * 9 / 10) << "the overlay's pieces lie on the lines the density sees";
    // Away from every line the flow fades: a point 20+ line widths from a line has none of it.
    int far = 0;
    for (int k = 0; k < 400 && !far; ++k) {
        const Vector3 d = Dir(sinf(k * 0.37f), cosf(k * 0.53f), sinf(k * 0.11f + 2.0f));
        if (f.lineDistance(kManaFire, d) > 25.0f * ManaField::kLineWidth) { EXPECT_EQ(f.lineFlow(kManaFire, d), 0.0f); far = 1; }
    }
    EXPECT_TRUE(far);
}

TEST(ManaField, NodesSitOnCrossingsAndTakeTheirKindFromTheAlgebra) {
    ASSERT_TRUE(Mana_Load());
    ManaFieldParams p = TestParams(42);
    p.meanTempC = [](Vector3 d) { return 30.0f - 60.0f * fabsf(d.y); }; // frozen poles
    ManaField f(p);
    // A 40 x 40 km patch at mid latitude, and one near a pole.
    int count = 0, aetherNodes = 0;
    std::vector<uint32_t> ids;
    for (const Vector3 centre : { Dir(1.0f, 0.2f, 0.3f), Dir(0.1f, 1.0f, 0.15f) }) {
        std::vector<ManaNode> nodes;
        f.nodesNear(centre, 20000.0f, nodes);
        for (const ManaNode& n : nodes) {
            ids.push_back(n.id);
            ++count;
            EXPECT_LT(f.lineDistance(n.a, n.dir), 5.0f) << "on a line of its first kind";
            EXPECT_LT(f.lineDistance(n.b, n.dir), 5.0f) << "and of its second";
            EXPECT_NE(n.a, n.b);
            EXPECT_GE(n.strength, ManaField::kNodeMinStrength);
            const bool fire = n.a == kManaFire || n.b == kManaFire;
            if (n.kind == Mana_Aether()) { ++aetherNodes; EXPECT_TRUE(fire); EXPECT_LT(p.meanTempC(n.dir), -5.0f); }
            else {
                EXPECT_EQ(n.kind, Mana_Combine(Mana_Find(n.a == 0 ? "earth" : n.a == 1 ? "water" : n.a == 2 ? "fire" : "air"),
                                               Mana_Find(n.b == 0 ? "earth" : n.b == 1 ? "water" : n.b == 2 ? "fire" : "air")));
            }
        }
    }
    std::sort(ids.begin(), ids.end());
    EXPECT_EQ(std::unique(ids.begin(), ids.end()), ids.end()) << "each crossing is kept once, by the tile that owns it";
    const float perKm2 = count / (2.0f * 3.14159f * 20.0f * 20.0f);
    printf("mana field: %d nodes in two 20 km discs: %.2f per km^2 (%d Aether)\n", count, perKm2, aetherNodes);
    EXPECT_GT(perKm2, 0.02f) << "the doc: a few per 10 km^2";
    EXPECT_LT(perKm2, 1.5f);
}

TEST(ManaField, DensityIsBoundedAndFollowsThePlace) {
    ASSERT_TRUE(Mana_Load());
    ManaField f(TestParams(42));
    ManaEnv calm; calm.tempC = 20.0f;
    ManaEnv high = calm; high.elevation = 2000.0f; high.wind = 12.0f; high.storm = 1.0f;
    ManaEnv cold = calm; cold.tempC = -20.0f; cold.night = 1.0f; cold.biome = (int)Biome::Tundra;
    ManaEnv swamp = calm; swamp.biome = (int)Biome::Swamp;
    const int air = Mana_Find("air"), frost = Mana_Find("frost"), storm = Mana_Find("storm"), rot = Mana_Find("rot"), aether = Mana_Find("aether");
    float airCalm = 0.0f, airHigh = 0.0f;
    for (int k = 0; k < 200; ++k) {
        const Vector3 d = Dir(sinf(k * 0.91f), cosf(k * 0.29f), sinf(k * 0.47f + 0.5f));
        for (const ManaEnv* env : { &calm, &high, &cold, &swamp }) {
            const ManaSample s = f.density(d, *env);
            for (int i = 0; i < kManaMaxKinds; ++i) { EXPECT_GE(s.v[i], 0.0f); EXPECT_LE(s.v[i], 2.0f); }
            for (int i = 0; i < Mana_Count(); ++i) if (!Mana_Kind(i).slice) { EXPECT_EQ(s.v[i], 0.0f) << Mana_Kind(i).id << ": only the slice"; }
            EXPECT_GT(s.of(aether), 0.0f) << "a trace of Aether everywhere";
        }
        airCalm += f.density(d, calm).of(air);
        airHigh += f.density(d, high).of(air);
        EXPECT_EQ(f.density(d, calm).of(frost), 0.0f) << "no frost in warm air";
        EXPECT_GT(f.density(d, cold).of(frost), 0.5f);
        EXPECT_GT(f.density(d, high).of(storm), f.density(d, calm).of(storm));
        EXPECT_GT(f.density(d, swamp).of(rot), 0.4f);
    }
    EXPECT_GT(airHigh, 2.0f * airCalm) << "Air rises with height, wind and storms";
}
#endif
