#include "ring.h"
#include "bodies.h"
#include "player.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"
#include "SimplexNoise.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <random>

float RingRockGravity(float radius) { return std::clamp(0.8f + 0.006f * radius, 0.8f, 4.0f); }

namespace {
    constexpr int kVariants = 12;
    constexpr int kLevels = 2;            // 0: 320 triangles (near), 1: 80 (far)
    constexpr int kBins = 2048;
    constexpr double kRebinSeconds = 5.0;
    constexpr double kSmallReach = 40000.0;  // small rocks are drawn within this (if big enough on screen)
    constexpr float kMinAngularSize = 0.0006f; // ~1 pixel at 1080p
    constexpr int kMaxPromoted = 16;
    constexpr float kMeanSurface = 0.72f;    // the asteroid field's mean surface, as a fraction of its radius

    struct RockMesh { std::vector<float> pos, nrm; int vertexCount = 0; unsigned int vboPos = 0, vboNrm = 0; };
    RockMesh gMeshes[kLevels][kVariants];
    bool gMeshesBuilt = false;
    HotShader gRockShader{ "resources/ring_rocks.vs", "resources/ring_rocks.fs" };
    HotShader gBandShader{ "resources/ring_band.vs", "resources/ring_band.fs" };
    double gShaderCheck = -1e9;

    // An icosphere, subdivided, displaced by noise: a lumpy unit-ish rock (mean radius ~1).
    void BuildRock(RockMesh& m, int subdiv, int variant) {
        const float t = (1.0f + sqrtf(5.0f)) / 2.0f;
        std::vector<Vector3> v = { { -1, t, 0 }, { 1, t, 0 }, { -1, -t, 0 }, { 1, -t, 0 }, { 0, -1, t }, { 0, 1, t }, { 0, -1, -t }, { 0, 1, -t }, { t, 0, -1 }, { t, 0, 1 }, { -t, 0, -1 }, { -t, 0, 1 } };
        std::vector<int> f = { 0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11, 1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
                               3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9, 4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1 };
        for (auto& p : v) p = Vector3Normalize(p);
        for (int s = 0; s < subdiv; ++s) {
            std::vector<int> nf;
            std::unordered_map<uint64_t, int> mids;
            auto mid = [&](int a, int b) {
                const uint64_t key = a < b ? ((uint64_t)a << 32) | (uint32_t)b : ((uint64_t)b << 32) | (uint32_t)a;
                auto it = mids.find(key);
                if (it != mids.end()) return it->second;
                v.push_back(Vector3Normalize(Vector3Scale(Vector3Add(v[a], v[b]), 0.5f)));
                return mids[key] = (int)v.size() - 1;
            };
            for (size_t i = 0; i < f.size(); i += 3) {
                const int a = f[i], b = f[i + 1], c = f[i + 2], ab = mid(a, b), bc = mid(b, c), ca = mid(c, a);
                nf.insert(nf.end(), { a, ab, ca, b, bc, ab, c, ca, bc, ab, bc, ca });
            }
            f = nf;
        }
        // Displace (same field for both detail levels of a variant).
        const float o = variant * 17.3f;
        for (auto& p : v) {
            const float n = SimplexNoise::noise(p.x * 1.3f + o, p.y * 1.3f, p.z * 1.3f) * 0.22f + SimplexNoise::noise(p.x * 3.1f, p.y * 3.1f + o, p.z * 3.1f) * 0.08f;
            const float squash = 1.0f - 0.25f * (variant % 3) / 2.0f * p.y * p.y; // some flatter rocks
            p = Vector3Scale(p, (1.0f + n) * squash);
        }
        // Flat-shaded triangles (faceted rock).
        for (size_t i = 0; i < f.size(); i += 3) {
            const Vector3 a = v[f[i]], b = v[f[i + 1]], c = v[f[i + 2]];
            const Vector3 n = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
            for (const Vector3& q : { a, b, c }) {
                m.pos.insert(m.pos.end(), { q.x, q.y, q.z });
                m.nrm.insert(m.nrm.end(), { n.x, n.y, n.z });
            }
        }
        m.vertexCount = (int)(m.pos.size() / 3);
    }

    void EnsureMeshes() {
        if (gMeshesBuilt) return;
        gMeshesBuilt = true;
        for (int l = 0; l < kLevels; ++l)
            for (int k = 0; k < kVariants; ++k) {
                RockMesh& m = gMeshes[l][k];
                BuildRock(m, l == 0 ? 2 : 1, k);
                m.vboPos = rlLoadVertexBuffer(m.pos.data(), (int)(m.pos.size() * sizeof(float)), false);
                m.vboNrm = rlLoadVertexBuffer(m.nrm.data(), (int)(m.nrm.size() * sizeof(float)), false);
            }
    }

    // The band: an annulus in the parent's equatorial plane, `segments` around, two radii.
    struct Band { unsigned int vao = 0, vbo = 0; int vertexCount = 0; } gBand;
    void EnsureBand() {
        if (gBand.vao) return;
        // Unit annulus: x = cos, z = sin, y = 0 (inner) / 1 (outer) radial coordinate in the shader.
        std::vector<float> d;
        const int seg = 1024;
        for (int i = 0; i < seg; ++i) {
            const float a0 = 2.0f * PI * i / seg, a1 = 2.0f * PI * (i + 1) / seg;
            const float q[4][3] = { { cosf(a0), 0, sinf(a0) }, { cosf(a1), 0, sinf(a1) }, { cosf(a1), 1, sinf(a1) }, { cosf(a0), 1, sinf(a0) } };
            const int tri[6] = { 0, 1, 2, 0, 2, 3 };
            for (int t : tri) d.insert(d.end(), { q[t][0], q[t][1], q[t][2] });
        }
        gBand.vertexCount = (int)(d.size() / 3);
        gBand.vao = rlLoadVertexArray();
        rlEnableVertexArray(gBand.vao);
        gBand.vbo = rlLoadVertexBuffer(d.data(), (int)(d.size() * sizeof(float)), false);
        rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0);
        rlEnableVertexAttribute(0);
        rlDisableVertexArray();
    }

    void CheckShaders() {
        if (GetTime() - gShaderCheck < 1.0) return;
        gShaderCheck = GetTime();
        gRockShader.reloadIfChanged();
        gBandShader.reloadIfChanged();
    }

    uint32_t Hash(uint32_t h) { h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16; return h; }
}

RingSystem::RingSystem(const BodyDef& d, Planetoid& parentBody, uint32_t seed, int ringNo_)
    : ringNo(ringNo_), name(d.name), parentId(parentBody.id), parent(&parentBody) {
    const double pr = parentBody.radius * parentBody.scale;
    const double gr = parentBody.gravityRadius > 0.0f ? parentBody.gravityRadius : pr;
    const double gm = parentBody.surfaceGravity * gr * gr;
    const double base = pr + d.orbitAltitude;
    innerRadius = base - d.extent * 0.5;
    outerRadius = base + d.extent * 0.5;
    thickness = d.thickness;
    std::mt19937 rng(seed ^ 0x52494e47u ^ (uint32_t)ringNo * 7919u);
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    rocks.resize((size_t)std::max(0, d.count));
    for (size_t k = 0; k < rocks.size(); ++k) {
        RingRock& r = rocks[k];
        r.orbitRadius = base + (unit(rng) - 0.5) * d.extent;
        r.height = (unit(rng) + unit(rng) - 1.0) * 0.5 * d.thickness;  // thicker in the middle
        r.phase = unit(rng) * 2.0 * PI;
        const double dist = sqrt(r.orbitRadius * r.orbitRadius + r.height * r.height);
        r.rate = sqrt(gm / (dist * dist * dist));
        // Sizes: many small (d.radius x size range, skewed small), ~1% large (5x).
        const bool big = unit(rng) < 0.01;
        const double v = unit(rng);
        float visual = (float)(d.radius * d.sizeMin * pow(d.sizeMax / std::max(1e-3f, d.sizeMin), v * v));
        if (big) { visual *= 5.0f; ++largeCount; }
        r.radius = visual / kMeanSurface;
        r.spin = (float)(unit(rng) - 0.5) * 2.0f * d.maxSpinDegPerSec;
        r.seed = (uint32_t)rng();
        r.variant = (uint8_t)(Hash(r.seed) % kVariants);
        const uint32_t h = Hash(r.seed ^ 0xabcu);
        const unsigned char g = (unsigned char)(95 + h % 60), warm = (unsigned char)((h >> 8) % 25);
        r.tint = { (unsigned char)(g + warm), g, (unsigned char)(g - std::min<int>(g, warm / 2)), 255 };
        if (big) large.push_back((int)k);
    }
}

RingSystem::~RingSystem() {
    for (size_t i = 0; i < vao.size(); ++i) {
        if (vao[i]) rlUnloadVertexArray(vao[i]);
        if (vbo[i]) rlUnloadVertexBuffer(vbo[i]);
    }
}

double RingSystem::angleOf(const RingRock& r) const { return r.phase + r.rate * WorldTime(); }

DVec3 RingSystem::rockPosition(int i) const {
    const RingRock& r = rocks[i];
    const double a = angleOf(r);
    return parent->worldPos + DVec3(cos(a) * r.orbitRadius, r.height, sin(a) * r.orbitRadius);
}

void RingSystem::rebin() {
    if ((int)bins.size() != kBins) bins.assign(kBins, {});
    for (auto& b : bins) b.clear(); // keep their capacity: this runs every few seconds
    for (size_t k = 0; k < rocks.size(); ++k) {
        double a = fmod(angleOf(rocks[k]), 2.0 * PI);
        if (a < 0) a += 2.0 * PI;
        if (rocks[k].radius > 0.0f) bins[std::min(kBins - 1, (int)(a / (2.0 * PI) * kBins))].push_back((int)k);
    }
    binnedAt = WorldTime();
}

// Every rock in the bins covering [angle - halfWidth, angle + halfWidth] (with a bin of slack for
// the drift since the bins were built).
template <class F> static void ForBinsAround(const std::vector<std::vector<int>>& bins, double angle, double halfWidth, F&& f) {
    const double binW = 2.0 * PI / kBins;
    const int span = std::min(kBins / 2, (int)ceil(halfWidth / binW) + 1);
    double a = fmod(angle, 2.0 * PI);
    if (a < 0) a += 2.0 * PI;
    const int c = (int)(a / binW);
    for (int i = -span; i <= span; ++i)
        for (int rock : bins[((c + i) % kBins + kBins) % kBins]) f(rock);
}

Planetoid* RingSystem::promotedBody(int rock) const {
    auto it = promoted.find(rock);
    return it == promoted.end() ? nullptr : Planetoid::find(it->second);
}

int RingSystem::nearestRock(DVec3 w, double* distance) const {
    if (bins.empty()) return -1;
    const DVec3 rel = w - parent->worldPos;
    const double ang = atan2(rel.z, rel.x);
    int best = -1;
    double bestD = 1e300;
    ForBinsAround(bins, ang, 3000.0 / innerRadius, [&](int k) {
        const double d = (rockPosition(k) - w).length() - rocks[k].radius * kMeanSurface;
        if (d < bestD) { bestD = d; best = k; }
    });
    if (distance) *distance = bestD;
    return best;
}

static double PromoteReach(const RingRock& r) { return std::max(250.0, 14.0 * r.radius) + r.radius; }

void RingSystem::update(Scene& world, DVec3 focus) {
    if (Planetoid::find(parentId) != parent) return;
    if (bins.empty() || fabs(WorldTime() - binnedAt) > kRebinSeconds) rebin();
    // Demote rocks left behind (not those with a build plot: its machines need the body).
    for (auto it = promoted.begin(); it != promoted.end();) {
        Planetoid* body = Planetoid::find(it->second);
        const int rock = it->first;
        if (!body) { it = promoted.erase(it); continue; }
        const double d = (body->worldPos - focus).length();
        if (d > 1.6 * PromoteReach(rocks[rock]) && !(body->plots && body->plots->partCount() > 0)) {
            it = promoted.erase(it);
            demote(world, rock, body);
        } else ++it;
    }
    // Promote nearby rocks, nearest first, up to the cap.
    const DVec3 rel = focus - parent->worldPos;
    const double ang = atan2(rel.z, rel.x);
    const double rho = sqrt(rel.x * rel.x + rel.z * rel.z);
    if (rho < innerRadius - 3000.0 || rho > outerRadius + 3000.0 || fabs(rel.y) > thickness + 3000.0) return; // far from the band
    std::vector<std::pair<double, int>> want;
    ForBinsAround(bins, ang, 3000.0 / innerRadius, [&](int k) {
        if (promoted.count(k)) return;
        const double d = (rockPosition(k) - focus).length();
        if (d < PromoteReach(rocks[k])) want.push_back({ d, k });
    });
    std::sort(want.begin(), want.end());
    for (const auto& w : want) {
        if ((int)promoted.size() >= kMaxPromoted) break;
        promote(world, w.second);
    }
}

void RingSystem::promote(Scene& world, int rock) {
    const RingRock& r = rocks[rock];
    auto p = std::make_unique<Planetoid>(TextFormat("%s_%d", name.c_str(), rock), rockPosition(rock), Vector3{ 0, 0, 0 }, 1.0f, r.seed, r.radius);
    p->kind = BodyKind::Asteroid;
    p->spin = { 0.0f, r.spin, 0.0f };
    p->rotation = { 0.0f, (float)fmod(r.spin * WorldTime(), 360.0), 0.0f };
    p->surfaceGravity = RingRockGravity(r.radius * kMeanSurface);
    p->orbitParent = parent;
    p->orbitParentId = parentId;
    p->orbitRadius = r.orbitRadius;
    p->orbitHeight = r.height;
    p->orbitRate = r.rate;
    p->orbitPhase = r.phase;
    p->ringNo = ringNo;
    p->ringRock = rock;
    p->parent = &world.rootObject;
    auto st = dormant.find(rock);
    if (st != dormant.end()) {
        RingRockState& s = st->second;
        if (s.edits) p->edits = s.edits;
        p->dirtyRegions = std::move(s.dirtyRegions);
        p->allRegionsDirty = s.allRegionsDirty;
        p->regionsSavedTo = s.regionsSavedTo;
        p->plots = std::move(s.plots);
        dormant.erase(st);
    }
    PlaceBodyOnOrbit(*p);
    promoted[rock] = p->id;
    world.objects.push_back(std::move(p));
    world.syncRenderPositions();
}

void RingSystem::demote(Scene& world, int rock, Planetoid* body) {
    // Keep what the visit changed.
    if ((body->edits && body->edits->strokes > 0) || (body->plots && body->plots->partCount() > 0) || !body->dirtyRegions.empty()) {
        RingRockState& s = dormant[rock];
        s.edits = body->edits;
        s.dirtyRegions = std::move(body->dirtyRegions);
        s.allRegionsDirty = body->allRegionsDirty;
        s.regionsSavedTo = body->regionsSavedTo;
        s.plots = std::move(body->plots);
    }
    if (player.nearest == body) player.nearest = nullptr; // the rest validate their body by id
    for (auto it = world.objects.begin(); it != world.objects.end(); ++it)
        if (it->get() == body) { world.objects.erase(it); break; }
}

bool RingSystem::bandRange(Scene& world, double& dmin, double& dmax) const {
    if (Planetoid::find(parentId) != parent) return false;
    const DVec3 cam = world.renderOrigin + world.camera.position;
    const DVec3 rel = cam - parent->worldPos;
    const double rho = sqrt(rel.x * rel.x + rel.z * rel.z);
    const double nearRho = std::clamp(rho, innerRadius - 60000.0, outerRadius + 60000.0);
    dmin = sqrt((rho - nearRho) * (rho - nearRho) + rel.y * rel.y);
    dmax = sqrt((rho + outerRadius + 60000.0) * (rho + outerRadius + 60000.0) + rel.y * rel.y);
    dmin = std::max(dmin, 20000.0); // nearer than this it has faded out (the rocks take over)
    return dmax > dmin;
}

void RingSystem::drawRocks(Scene& world, const RenderContext& ctx, bool rebuild) {
    (void)ctx.frustum;
    if (Planetoid::find(parentId) != parent) return;
    EnsureMeshes();
    CheckShaders();
    Shader& sh = gRockShader.shader;
    if (sh.id <= 0) return;
    const int slots = kVariants * kLevels;
    if (vao.empty()) { vao.assign(slots, 0); vbo.assign(slots, 0); instCount.assign(slots, 0); }
    if (rebuild) {
        if (bins.empty()) rebin();
        std::vector<std::vector<float>> buf(slots);
        const DVec3 cam = world.renderOrigin + world.camera.position;
        const DVec3 rel = cam - parent->worldPos;
        const double ang = atan2(rel.z, rel.x);
        const double t = WorldTime();
        const Vector3 fwd = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        drawnInstances = 0;
        auto add = [&](int k) {
            if (promoted.count(k)) return;
            const RingRock& r = rocks[k];
            const DVec3 w = rockPosition(k);
            const Vector3 pos = (w - world.renderOrigin).toFloat();
            const float d = Vector3Distance(pos, world.camera.position);
            const float vr = r.radius * kMeanSurface;
            if (vr / std::max(d, 1.0f) < kMinAngularSize) return;
            // A view cone, not the pass frustum: the mid pass builds the list and its frustum starts at 1.8 km.
            const Vector3 to = Vector3Subtract(pos, world.camera.position);
            if (d > vr * 2.0f && Vector3DotProduct(to, fwd) < 0.45f * d - vr * 1.4f) return;
            const float a = (float)fmod(r.spin * t, 360.0) * DEG2RAD;
            const int level = d < 4000.0f ? 0 : 1;
            std::vector<float>& b = buf[level * kVariants + r.variant];
            b.insert(b.end(), { pos.x, pos.y, pos.z, vr, 0.0f, sinf(a * 0.5f), 0.0f, cosf(a * 0.5f), r.tint.r / 255.0f, r.tint.g / 255.0f, r.tint.b / 255.0f, 1.0f });
        };
        // Small rocks near the camera (when it is near the band), large ones from further.
        const double rho = sqrt(rel.x * rel.x + rel.z * rel.z);
        if (rho > innerRadius - kSmallReach && rho < outerRadius + kSmallReach && fabs(rel.y) < kSmallReach)
            ForBinsAround(bins, ang, kSmallReach / innerRadius, add);
        for (int k : large) add(k); // (the bin pass may have added some already: skip duplicates below)
        for (int s = 0; s < slots; ++s) {
            // Large rocks near the camera were visited twice: drop exact duplicate positions.
            std::vector<float>& b = buf[s];
            if (b.size() > 12) {
                std::vector<float> u;
                u.reserve(b.size());
                std::unordered_set<uint64_t> seen;
                for (size_t i = 0; i + 11 < b.size(); i += 12) {
                    uint64_t key; float k3[2] = { b[i], b[i + 2] }; memcpy(&key, k3, 8);
                    if (!seen.insert(key).second) continue;
                    u.insert(u.end(), b.begin() + i, b.begin() + i + 12);
                }
                b.swap(u);
            }
            instCount[s] = (int)(b.size() / 12);
            drawnInstances += instCount[s];
            if (!instCount[s]) continue;
            const RockMesh& m = gMeshes[s / kVariants][s % kVariants];
            if (!vao[s]) {
                vao[s] = rlLoadVertexArray();
                rlEnableVertexArray(vao[s]);
                rlEnableVertexBuffer(m.vboPos); rlSetVertexAttribute(0, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(0);
                rlEnableVertexBuffer(m.vboNrm); rlSetVertexAttribute(1, 3, RL_FLOAT, false, 0, 0); rlEnableVertexAttribute(1);
            } else {
                rlEnableVertexArray(vao[s]);
            }
            if (vbo[s]) rlUnloadVertexBuffer(vbo[s]);
            vbo[s] = rlLoadVertexBuffer(b.data(), (int)(b.size() * sizeof(float)), true);
            const int stride = 12 * (int)sizeof(float);
            rlSetVertexAttribute(10, 4, RL_FLOAT, false, stride, 0); rlEnableVertexAttribute(10); rlSetVertexAttributeDivisor(10, 1);
            rlSetVertexAttribute(11, 4, RL_FLOAT, false, stride, 4 * sizeof(float)); rlEnableVertexAttribute(11); rlSetVertexAttributeDivisor(11, 1);
            rlSetVertexAttribute(12, 4, RL_FLOAT, false, stride, 8 * sizeof(float)); rlEnableVertexAttribute(12); rlSetVertexAttributeDivisor(12, 1);
            rlDisableVertexArray();
            rlDisableVertexBuffer();
        }
    }
    if (!drawnInstances) return;
    rlDrawRenderBatchActive();
    rlEnableShader(sh.id);
    const Matrix vp = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
    rlSetUniformMatrix(GetShaderLocation(sh, "matViewProj"), vp);
    rlSetUniform(GetShaderLocation(sh, "sunDir"), &ctx.sunDir, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(sh, "ambient"), &ctx.ambient, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(sh, "viewPos"), &ctx.cameraPos, RL_SHADER_UNIFORM_VEC3, 1);
    const Vector3 pc = (parent->worldPos - ctx.renderOrigin).toFloat();
    const float occ[4] = { pc.x, pc.y, pc.z, parent->radius * parent->scale };
    rlSetUniform(GetShaderLocation(sh, "sunOccluder"), occ, RL_SHADER_UNIFORM_VEC4, 1);
    for (size_t s = 0; s < vao.size(); ++s) {
        if (!instCount[s]) continue;
        rlEnableVertexArray(vao[s]);
        rlDrawVertexArrayInstanced(0, gMeshes[s / kVariants][s % kVariants].vertexCount, instCount[s]);
        rlDisableVertexArray();
        if (ctx.drawCalls) ++*ctx.drawCalls;
    }
    rlDisableShader();
}

void RingSystem::drawBand(Scene& world, const RenderContext& ctx) {
    if (Planetoid::find(parentId) != parent) return;
    EnsureBand();
    CheckShaders();
    Shader& sh = gBandShader.shader;
    if (sh.id <= 0) return;
    rlDrawRenderBatchActive();
    rlEnableShader(sh.id);
    const Matrix vp = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
    rlSetUniformMatrix(GetShaderLocation(sh, "matViewProj"), vp);
    const Vector3 pc = (parent->worldPos - world.renderOrigin).toFloat();
    // The drawn band reaches a dust halo either side of the rocks (the rock band alone is a thread
    // from the ground: 3 km wide at ~270 km).
    const double halo = 60000.0;
    const float ring[4] = { pc.x, pc.y, pc.z, (float)(innerRadius - halo) };
    const float core[2] = { (float)innerRadius, (float)outerRadius };
    rlSetUniform(GetShaderLocation(sh, "coreRadii"), core, RL_SHADER_UNIFORM_VEC2, 1);
    const float outer = (float)(outerRadius + halo);
    const float occ[4] = { pc.x, pc.y, pc.z, parent->radius * parent->scale };
    rlSetUniform(GetShaderLocation(sh, "ringCentre"), ring, RL_SHADER_UNIFORM_VEC4, 1);
    rlSetUniform(GetShaderLocation(sh, "outerRadius"), &outer, RL_SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(GetShaderLocation(sh, "sunDir"), &ctx.sunDir, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(sh, "ambient"), &ctx.ambient, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(sh, "viewPos"), &ctx.cameraPos, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(sh, "sunOccluder"), occ, RL_SHADER_UNIFORM_VEC4, 1);
    const float seedF = (float)(ringNo * 13 + 7);
    rlSetUniform(GetShaderLocation(sh, "seed"), &seedF, RL_SHADER_UNIFORM_FLOAT, 1);
    rlDisableBackfaceCulling();
    rlDisableDepthMask();
    rlEnableColorBlend();
    rlSetBlendMode(BLEND_ALPHA);
    rlEnableVertexArray(gBand.vao);
    rlDrawVertexArray(0, gBand.vertexCount);
    rlDisableVertexArray();
    rlEnableDepthMask();
    rlEnableBackfaceCulling();
    rlDisableShader();
    if (ctx.drawCalls) ++*ctx.drawCalls;
}

void RingSystem::save(FILE* f, const std::string& dir) {
    // Rocks to save: promoted ones with changes, and dormant ones.
    struct Item { int rock; const TerrainDeltas* edits; PlotSet* plots; Planetoid* body; RingRockState* state; };
    std::vector<Item> items;
    for (auto& kv : promoted) {
        Planetoid* b = Planetoid::find(kv.second);
        if (b && ((b->edits && b->edits->strokes > 0) || (b->plots && b->plots->partCount() > 0))) items.push_back({ kv.first, b->edits.get(), b->plots.get(), b, nullptr });
    }
    for (auto& kv : dormant) items.push_back({ kv.first, kv.second.edits.get(), kv.second.plots.get(), nullptr, &kv.second });
    for (const Item& it : items) {
        fprintf(f, "ringrock %d %d strokes %zu\n", ringNo, it.rock, it.edits ? it.edits->strokes : (size_t)0);
        if (it.plots) it.plots->save(f); else fprintf(f, "plots2 0\n");
        if (it.body) {
            SaveRegionsFor(*it.body->edits, it.body->dirtyRegions, it.body->allRegionsDirty, it.body->regionsSavedTo, regionBodyId(it.rock), dir);
        } else if (it.state->edits) {
            SaveRegionsFor(*it.state->edits, it.state->dirtyRegions, it.state->allRegionsDirty, it.state->regionsSavedTo, regionBodyId(it.rock), dir);
        }
    }
}

bool RingSystem::loadRock(FILE* f, int rock, size_t strokes, const std::string& dir) {
    auto plots = std::make_unique<PlotSet>();
    const bool ok = plots->load(f);
    if (rock < 0 || rock >= (int)rocks.size()) return ok;
    RingRockState& s = dormant[rock];
    s.edits = LoadRegionsFor(regionBodyId(rock), dir, strokes);
    s.allRegionsDirty = false;
    s.regionsSavedTo = dir;
    s.plots = plots->plots.empty() ? nullptr : std::move(plots);
    return ok;
}

int RingSystem::clearLane(double orbitR, double h, double bodyRadius) {
    int n = 0;
    for (RingRock& r : rocks) {
        if (r.radius <= 0.0f) continue;
        const double reach = 0.72 * (bodyRadius + r.radius) + 5.0; // mean surfaces, a little clearance
        if (fabs(r.orbitRadius - orbitR) < reach && fabs(r.height - h) < reach) { r.radius = 0.0f; ++n; }
    }
    for (auto it = large.begin(); it != large.end();) it = rocks[*it].radius <= 0.0f ? large.erase(it) : it + 1;
    clearedCount += n;
    bins.clear(); // rebuilt without them
    return n;
}
