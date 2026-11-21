#include "audioscene.h"

#include "plots.h"       // kCellSize
#include "scene.h"
#include "world.h"
#include "worldtypes.h"

#include <algorithm>
#include <cmath>

AudioMat AudioMatOfVoxel(VoxelMaterial m) {
    switch (m) {
        case VoxelMaterial::Regolith: case VoxelMaterial::Gravel: case VoxelMaterial::Pebbles: return AudioMat::Gravel;
        case VoxelMaterial::Sand: case VoxelMaterial::Grass: case VoxelMaterial::Clay: case VoxelMaterial::Silt:
        case VoxelMaterial::Mud: case VoxelMaterial::Farmland: return AudioMat::Soil;
        case VoxelMaterial::Snow: return AudioMat::Snow;
        default: return AudioMat::Rock; // rock, its strata, ore, ice
    }
}

namespace {
    Vector3 Local(const AudioScene::Body& b, Vector3 p) {
        const double x = p.x, y = p.y, z = p.z;
        return { (float)(b.b.x + b.A[0][0] * x + b.A[0][1] * y + b.A[0][2] * z),
                 (float)(b.b.y + b.A[1][0] * x + b.A[1][1] * y + b.A[1][2] * z),
                 (float)(b.b.z + b.A[2][0] * x + b.A[2][1] * y + b.A[2][2] * z) };
    }
    Vector3 LocalDir(const AudioScene::Body& b, Vector3 d) {
        return { (float)(b.A[0][0] * d.x + b.A[0][1] * d.y + b.A[0][2] * d.z),
                 (float)(b.A[1][0] * d.x + b.A[1][1] * d.y + b.A[1][2] * d.z),
                 (float)(b.A[2][0] * d.x + b.A[2][1] * d.y + b.A[2][2] * d.z) };
    }
    Vector3 At(Vector3 o, Vector3 d, float t) { return { o.x + d.x * t, o.y + d.y * t, o.z + d.z * t }; }
    float Len(Vector3 v) { return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); }
    float Dot(Vector3 a, Vector3 b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
}

float AudioScene::terrainHit(const Body& b, Vector3 o, Vector3 dir, float lo, float hi, AudioMat* mat) const {
    const DensityField& f = *b.field;
    const TerrainDeltas* e = b.edits.get();
    auto dens = [&](float t) { return f.density(Local(b, At(o, dir, t)), e); };
    auto hitAt = [&](float t) {
        if (mat) {
            const Vector3 L = Local(b, At(o, dir, t));
            const int m = e ? e->materialAt(L) : -1;
            *mat = AudioMatOfVoxel(m >= 0 ? (VoxelMaterial)m : f.material(L, f.baseDensity(L)));
        }
        return t;
    };
    float t = lo, prev = lo;
    float d = dens(t);
    if (d >= ISO_LEVEL) return hitAt(t);
    while (t < hi) {
        // The density says roughly how far the surface is (DENSITY_GRADIENT per local unit across it): step half that,
        // 0.25-8 m, so open air is crossed in a few samples and a wall thicker than a quarter metre is never stepped over.
        const float gapM = (ISO_LEVEL - d) / DENSITY_GRADIENT / b.localPerMetre;
        prev = t;
        t = std::min(hi, t + std::clamp(gapM * 0.5f, 0.25f, 8.0f));
        d = dens(t);
        if (d >= ISO_LEVEL) {
            float a = prev, c = t;
            for (int i = 0; i < 6; ++i) { const float m = 0.5f * (a + c); if (dens(m) >= ISO_LEVEL) c = m; else a = m; }
            return hitAt(c);
        }
        if (t >= hi) break;
    }
    return -1.0f;
}

float AudioScene::plotHit(const Body& b, const GrassBlockPlot& p, Vector3 o, Vector3 dir, float lo, float hi, AudioMat* mat) const {
    // Into plot space (metres along east, up, north from the plot's origin, body-local units).
    const Vector3 L0 = Local(b, At(o, dir, lo));
    const DVec3 rel = DVec3(L0) - p.origin;
    const Vector3 r = { (float)rel.x, (float)rel.y, (float)rel.z };
    const Vector3 P0 = { Dot(r, p.east), Dot(r, p.up), Dot(r, p.north) };
    const Vector3 Ld = LocalDir(b, dir);
    const Vector3 D = { Dot(Ld, p.east), Dot(Ld, p.up), Dot(Ld, p.north) }; // per audio metre
    const float len = hi - lo;
    // Skip plots the segment passes far from (their cells lie within `radius` of the plot origin).
    {
        const float dd = Dot(D, D);
        const float s = dd > 0.0f ? std::clamp(-Dot(P0, D) / dd, 0.0f, len) : 0.0f;
        const Vector3 c = At(P0, D, s);
        if (Len(c) > p.radius + 2.0f * kCellSize) return -1.0f;
    }
    // Amanatides-Woo through the 0.5 m cells.
    int cell[3] = { (int)floorf(P0.x / kCellSize), (int)floorf(P0.y / kCellSize), (int)floorf(P0.z / kCellSize) };
    const float pos[3] = { P0.x, P0.y, P0.z }, dv[3] = { D.x, D.y, D.z };
    int step[3];
    float tMax[3], tDelta[3];
    for (int k = 0; k < 3; ++k) {
        if (dv[k] > 1e-9f) { step[k] = 1; tMax[k] = ((cell[k] + 1) * kCellSize - pos[k]) / dv[k]; tDelta[k] = kCellSize / dv[k]; }
        else if (dv[k] < -1e-9f) { step[k] = -1; tMax[k] = (cell[k] * kCellSize - pos[k]) / dv[k]; tDelta[k] = -kCellSize / dv[k]; }
        else { step[k] = 0; tMax[k] = 1e30f; tDelta[k] = 1e30f; }
    }
    float t = 0.0f;
    for (int it = 0; it < 4096 && t <= len; ++it) {
        const Int3 c{ cell[0], cell[1], cell[2] };
        if (p.cells.count(c)) {
            if (mat) *mat = p.solid.count(c) ? AudioMat::Brick : AudioMat::Glass; // panes, small parts: thin
            return lo + t;
        }
        const int k = tMax[0] < tMax[1] ? (tMax[0] < tMax[2] ? 0 : 2) : (tMax[1] < tMax[2] ? 1 : 2);
        t = tMax[k];
        tMax[k] += tDelta[k];
        cell[k] += step[k];
    }
    return -1.0f;
}

float AudioScene::closestHit(Vector3 o, Vector3 dir, float tMin, float tMax, AudioMat* mat, const Vector3* listenerAt) const {
    const Vector3 lis = listenerAt ? *listenerAt : listener;
    // The ends of the ray at a sound (not at the listener) get a margin.
    auto margin = [&](Vector3 end) { const Vector3 d = { end.x - lis.x, end.y - lis.y, end.z - lis.z }; return Len(d) < 0.5f ? 0.05f : kSourceMargin; };
    const float lo = tMin + margin(At(o, dir, tMin));
    const float hi = tMax - margin(At(o, dir, tMax));
    if (hi <= lo) return -1.0f;
    float best = -1.0f;
    AudioMat bestMat = AudioMat::Generic, m = AudioMat::Generic;
    for (const Body& b : bodies) {
        if (!b.field) continue;
        const float top = best >= 0.0f ? best : hi;
        const float t = terrainHit(b, o, dir, lo, top, mat ? &m : nullptr);
        if (t >= 0.0f && (best < 0.0f || t < best)) { best = t; bestMat = m; }
        if (b.plots)
            for (const GrassBlockPlot& p : *b.plots) {
                const float t2 = plotHit(b, p, o, dir, lo, best >= 0.0f ? best : hi, mat ? &m : nullptr);
                if (t2 >= 0.0f && (best < 0.0f || t2 < best)) { best = t2; bestMat = m; }
            }
    }
    if (mat && best >= 0.0f) *mat = bestMat;
    return best;
}

std::shared_ptr<const AudioScene> AudioScene_Snapshot(Scene& world, const DVec3& listenerGlobal, float range) {
    auto s = std::make_shared<AudioScene>();
    s->origin = listenerGlobal;
    for (auto& obj : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p || !p->field || p->vehicle) continue;
        const double reach = (double)p->field->maxRadius() * p->scale + range;
        if (!p->flat && (listenerGlobal - p->worldPos).length() > reach) continue;
        AudioScene::Body b;
        b.field = p->field;
        b.edits = p->edits;
        b.plots = Grass_BlockersOf(p->field.get());
        b.b = p->toLocalD(listenerGlobal);
        const DVec3 axes[3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
        for (int k = 0; k < 3; ++k) {
            const DVec3 c = p->toLocalD(listenerGlobal + axes[k]) - b.b; // column k: one global metre along axis k
            b.A[0][k] = c.x; b.A[1][k] = c.y; b.A[2][k] = c.z;
        }
        b.localPerMetre = p->scale > 0.0f ? 1.0f / p->scale : 1.0f;
        s->bodies.push_back(std::move(b));
    }
    return s;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1293): rays against a ground plane, a wall of part cells, the margins at a sound's end ----
#include "unit_helpers.h"

namespace {
// Solid below y = 0 (body-local), density changing DENSITY_GRADIENT per unit across it.
class PlaneField : public DensityField {
public:
    float baseDensity(Vector3 p) const override { return ISO_LEVEL - p.y * DENSITY_GRADIENT; }
    VoxelMaterial material(Vector3, float) const override { return VoxelMaterial::Grass; }
    Vector3 fallbackUp(Vector3) const override { return { 0, 1, 0 }; }
    float maxRadius() const override { return 1e6f; }
};

AudioScene PlaneScene(float listenerHeight) {
    AudioScene s;
    AudioScene::Body b;
    b.field = std::make_shared<PlaneField>();
    b.b = { 0, listenerHeight, 0 }; // the audio origin (the listener) stands this high over the ground
    s.bodies.push_back(b);
    return s;
}
} // namespace

TEST(AudioScene, GroundBlocksRaysThroughItOnly) {
    const AudioScene s = PlaneScene(1.7f);
    AudioMat m = AudioMat::Generic;
    // Straight down from the listener: ground 1.7 m below.
    const float t = s.closestHit({ 0, 0, 0 }, { 0, -1, 0 }, 0.0f, 10.0f, &m);
    EXPECT_NEAR(t, 1.7f, 0.05f);
    EXPECT_EQ(m, AudioMat::Soil);
    // Level, over the ground: nothing.
    EXPECT_FALSE(s.anyHit({ 0, 0, 0 }, { 1, 0, 0 }, 0.0f, 50.0f));
    // To a source 6 m down under the ground (a cave below): blocked.
    EXPECT_TRUE(s.anyHit({ 0, 0, 0 }, { 0, -1, 0 }, 0.0f, 7.7f));
}

TEST(AudioScene, SourceHalfInTheGroundIsNotOccludedByItself) {
    const AudioScene s = PlaneScene(1.7f);
    // A sound 10 m away, 0.1 m into the ground (the ray's last ~0.5 m is underground): its last 0.75 m is ignored.
    const Vector3 src = { 10.0f, -1.8f, 0.0f };
    const float d = Len(src);
    const Vector3 dir = { src.x / d, src.y / d, src.z / d };
    EXPECT_FALSE(s.anyHit({ 0, 0, 0 }, dir, 0.0f, d));
    // Traced from the sound's end instead (Steam Audio may go either way): the same.
    const Vector3 back = { -dir.x, -dir.y, -dir.z };
    EXPECT_FALSE(s.anyHit(src, back, 0.0f, d));
}

TEST(AudioScene, PartCellsBlockAndKnowGlass) {
    AudioScene s;
    AudioScene::Body b;
    b.field = std::make_shared<PlaneField>();
    b.b = { 0, 50, 0 }; // well over the ground
    auto plots = std::make_shared<GrassBlockers>();
    GrassBlockPlot p;
    p.origin = { 0, 50, 0 };
    p.east = { 1, 0, 0 }; p.up = { 0, 1, 0 }; p.north = { 0, 0, 1 };
    for (int y = -4; y < 4; ++y)
        for (int z = -4; z < 4; ++z) {
            p.cells.insert({ 10, y, z });                         // a wall at x = 5..5.5 m
            p.solid.insert({ 10, y, z });
            p.cells.insert({ -10, y, z });                        // a glass pane at x = -5..-4.5 m (not solid)
        }
    p.radius = 8.0f;
    plots->push_back(p);
    b.plots = plots;
    s.bodies.push_back(b);
    AudioMat m = AudioMat::Generic;
    EXPECT_NEAR(s.closestHit({ 0, 0, 0 }, { 1, 0, 0 }, 0.0f, 20.0f, &m), 5.0f, 1e-3f);
    EXPECT_EQ(m, AudioMat::Brick);
    EXPECT_NEAR(s.closestHit({ 0, 0, 0 }, { -1, 0, 0 }, 0.0f, 20.0f, &m), 4.5f, 1e-3f);
    EXPECT_EQ(m, AudioMat::Glass);
    EXPECT_FALSE(s.anyHit({ 0, 0, 0 }, { 0, 0, 1 }, 0.0f, 20.0f));    // along the walls
    EXPECT_FALSE(s.anyHit({ 0, 0, 0 }, { 1, 0, 0 }, 0.0f, 5.6f));     // a sound against the wall: its margin
}

TEST(AudioScene, MaterialsByClass) {
    EXPECT_EQ(AudioMatOfVoxel(VoxelMaterial::Granite), AudioMat::Rock);
    EXPECT_EQ(AudioMatOfVoxel(VoxelMaterial::Pebbles), AudioMat::Gravel);
    EXPECT_EQ(AudioMatOfVoxel(VoxelMaterial::Mud), AudioMat::Soil);
    EXPECT_EQ(AudioMatOfVoxel(VoxelMaterial::Snow), AudioMat::Snow);
}
#endif
