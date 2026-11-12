#include "world.h"
#include "ring.h"
#include "player.h"
#include "jobsystem.h"
#include "utils/logger.h"
#include "worldtypes.h"
#include "factory.h"
#include "fluid.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <deque>
#include <mutex>
#include <random>
#include <unordered_set>

const int CHUNK_SIZE = 32;

WorldGenSettings worldGen;
std::function<void(FILE*)> saveExtrasHook;
std::function<void(FILE*)> loadExtrasHook;
LodSettings lodSettings;

namespace {
    std::random_device rd;
    std::mt19937 gen(rd()); // reseeded by generateWorld()

    // ---- Async chunk pipeline -------------------------------------------------------------
    // Declared before the job system so it is destroyed after the workers are joined.
    std::mutex resultsMutex;
    std::vector<ChunkMeshData> results;        // filled by workers
    std::deque<ChunkMeshData> pendingUploads;  // main thread only
    std::unique_ptr<JobSystem> jobs;

    std::unordered_map<uint64_t, Planetoid*> planetoidRegistry; // main thread only
    uint64_t nextPlanetoidId = 1;

    WorldStats stats;
    double generationStart = -1.0; // GetTime() when a (re)generation started; <0 once reported
    double buildMsTotal = 0.0, sampleMsTotal = 0.0;
    long buildCount = 0;

    JobSystem& Jobs() {
        if (!jobs) jobs = std::make_unique<JobSystem>();
        return *jobs;
    }
}

JobSystem& WorldJobs() { return Jobs(); }

// ---------------------------------------------------------------------------------------------
// Planetoid
// ---------------------------------------------------------------------------------------------

Planetoid::Planetoid(std::string name, DVec3 wp, Vector3 rotation, float scale, unsigned int seed_, float radius_)
    : Object("planetoid", name, wp.toFloat(), rotation, WHITE, scale),
      id(nextPlanetoidId++), seed(seed_), radius(radius_),
      field(std::make_shared<PlanetoidField>(seed_, radius_)),
      edits(std::make_shared<TerrainDeltas>()) {
    worldPos = wp;
    hasWorldPos = true;
    surfaceGravity = AutoSurfaceGravity(BodyKind::Asteroid, radius_ * scale); // generateWorld may override
    gravityRadius = radius_ * scale * PlanetoidField::MEAN_SURFACE;         // the terrain's mean surface
    // Smallest level whose 2x2x2 root nodes (coords -1..0) cover the body plus room for added material.
    while (kVoxelSize * (float)(CHUNK_SIZE << rootLevel) < radius + 24.0f) rootLevel++;
    planetoidRegistry[id] = this;
}

Planetoid::~Planetoid() { planetoidRegistry.erase(id); }

Planetoid::Planetoid(std::string name, DVec3 wp, std::shared_ptr<const EmptyField> empty, float radius_)
    : Object("vehicle", name, wp.toFloat(), { 0, 0, 0 }, WHITE, 1.0f),
      id(nextPlanetoidId++), seed(0), radius(radius_), field(empty), edits(std::make_shared<TerrainDeltas>()) {
    worldPos = wp;
    hasWorldPos = true;
    surfaceGravity = 0.0f;
    gravityRadius = 1.0f;
    vehicle = true;
    while (kVoxelSize * (float)(CHUNK_SIZE << rootLevel) < radius + 4.0f) rootLevel++;
    planetoidRegistry[id] = this;
}

Planetoid::Planetoid(std::string name, DVec3 wp, std::shared_ptr<const TerrainField> terrain)
    : Object("terrain", name, wp.toFloat(), { 0, 0, 0 }, WHITE, 1.0f),
      id(nextPlanetoidId++), seed(0), radius(1e9f), field(terrain), edits(std::make_shared<TerrainDeltas>()) {
    worldPos = wp;
    hasWorldPos = true;
    flat = true;
    minHeight = terrain->minHeight;
    maxHeight = terrain->maxHeight;
    rootLevel = 4 + kVoxelShift; // 512 m root tiles, a 5x5 grid around the camera (~1.3 km view range)
    planetoidRegistry[id] = this;
}

Planetoid* Planetoid::find(uint64_t id) {
    auto it = planetoidRegistry.find(id);
    return it == planetoidRegistry.end() ? nullptr : it->second;
}

Vector3 Planetoid::toLocal(Vector3 world) const { return Vector3Transform(world, MatrixInvert(localTransform())); }
Vector3 Planetoid::toWorld(Vector3 local) const { return Vector3Transform(local, localTransform()); }

// The rotation of localTransform() (raymath's MatrixRotateXYZ) in double: world = R * (scale * local).
static void BodyRotationD(const Vector3& rotationDeg, double r[3][3]) {
    const double ax = -(double)rotationDeg.x * DEG2RAD, ay = -(double)rotationDeg.y * DEG2RAD, az = -(double)rotationDeg.z * DEG2RAD;
    const double cx = cos(ax), sx = sin(ax), cy = cos(ay), sy = sin(ay), cz = cos(az), sz = sin(az);
    // Rows of R are (m0, m4, m8), (m1, m5, m9), (m2, m6, m10) of MatrixRotateXYZ.
    r[0][0] = cz * cy;                 r[0][1] = sz * cy;                 r[0][2] = -sy;
    r[1][0] = cz * sy * sx - sz * cx;  r[1][1] = sz * sy * sx + cz * cx;  r[1][2] = cy * sx;
    r[2][0] = cz * sy * cx + sz * sx;  r[2][1] = sz * sy * cx - cz * sx;  r[2][2] = cy * cx;
}

DVec3 Planetoid::toWorldD(const DVec3& l) const {
    double r[3][3];
    BodyRotationD(rotation, r);
    const double s = scale;
    return worldPos + DVec3((r[0][0] * l.x + r[0][1] * l.y + r[0][2] * l.z) * s,
                            (r[1][0] * l.x + r[1][1] * l.y + r[1][2] * l.z) * s,
                            (r[2][0] * l.x + r[2][1] * l.y + r[2][2] * l.z) * s);
}

DVec3 Planetoid::toLocalD(const DVec3& w) const {
    double r[3][3];
    BodyRotationD(rotation, r);
    const DVec3 d = w - worldPos;
    const double inv = 1.0 / scale;
    return DVec3((r[0][0] * d.x + r[1][0] * d.y + r[2][0] * d.z) * inv,  // R is orthonormal: R^-1 = R^T
                 (r[0][1] * d.x + r[1][1] * d.y + r[2][1] * d.z) * inv,
                 (r[0][2] * d.x + r[1][2] * d.y + r[2][2] * d.z) * inv);
}

// Places a body's chunks in render space. Rotation and scale act on chunk-local vertices (small
// numbers), so float is fine there; the translation is the chunk origin's render-space position,
// computed in double so chunks far from the body centre neither jitter nor open seams as the
// camera moves. Built once per body per pass (the trigonometry is shared by every chunk).
namespace {
    struct ChunkPlacer {
        Matrix rs;
        double r[3][3];
        DVec3 offset;  // body position - render origin
        double s;
        ChunkPlacer(const Planetoid& p, const DVec3& renderOrigin)
            : rs(MatrixMultiply(MatrixScale(p.scale, p.scale, p.scale), MatrixRotateXYZ({ DEG2RAD * p.rotation.x, DEG2RAD * p.rotation.y, DEG2RAD * p.rotation.z }))),
              offset(p.worldPos - renderOrigin), s(p.scale) {
            BodyRotationD(p.rotation, r);
        }
        Matrix model(const NodeKey& key) const { const Vector3 o = key.origin(); return modelAt(DVec3(o)); }
        Matrix modelAt(const DVec3& o) const {
            Matrix m = rs;
            m.m12 = (float)(offset.x + (r[0][0] * o.x + r[0][1] * o.y + r[0][2] * o.z) * s);
            m.m13 = (float)(offset.y + (r[1][0] * o.x + r[1][1] * o.y + r[1][2] * o.z) * s);
            m.m14 = (float)(offset.z + (r[2][0] * o.x + r[2][1] * o.y + r[2][2] * o.z) * s);
            return m;
        }
    };
}

Matrix Planetoid::localModel(const DVec3& anchorLocal, const DVec3& renderOrigin) const {
    return ChunkPlacer(*this, renderOrigin).modelAt(anchorLocal);
}

Matrix Planetoid::chunkModel(const Chunk& chunk, const DVec3& renderOrigin) const {
    return ChunkPlacer(*this, renderOrigin).model(chunk.key);
}

VoxelMaterial Planetoid::materialAt(Vector3 local) const {
    if (edits) { const int m = edits->materialAt(local); if (m >= 0) return (VoxelMaterial)m; }
    return field->material(local, field->baseDensity(local));
}

float Planetoid::densityAt(Vector3 local) const {
    return field->density(local, edits.get());
}

Vector3 Planetoid::surfaceVelocity(Vector3 renderPoint) const {
    // The body's own motion (orbit) plus its spin: omega x (point - centre). Omega comes from the
    // rotation change over one second (rotation is Euler angles, so any spin axis works); taking it
    // as a vector keeps this exact far from the centre (differentiating positions in float lost
    // ~0.4 m/s at 60 km).
    if (spin.x == 0.0f && spin.y == 0.0f && spin.z == 0.0f) return linearVelocity;
    const Vector3 r1 = Vector3Add(rotation, spin);
    const Quaternion q0 = rotationQuat();
    const Quaternion q1 = QuaternionFromMatrix(MatrixRotateXYZ({ DEG2RAD * r1.x, DEG2RAD * r1.y, DEG2RAD * r1.z }));
    Vector3 axis; float angle;
    QuaternionToAxisAngle(QuaternionNormalize(QuaternionMultiply(q1, QuaternionInvert(q0))), &axis, &angle);
    if (angle > PI) angle -= 2.0f * PI;
    const Vector3 omega = Vector3Scale(axis, angle); // rad/s
    return Vector3Add(linearVelocity, Vector3CrossProduct(omega, Vector3Subtract(renderPoint, position)));
}

// Directions only rotate (no translation): transforming two points and subtracting lost ~0.2 deg in
// float with the body centre tens of km from the render origin (hundreds of metres on the surface).
Vector3 Planetoid::dirToLocal(Vector3 worldDir) const {
    double r[3][3];
    BodyRotationD(rotation, r);
    const Vector3 d = worldDir;
    return Vector3Normalize({ (float)(r[0][0] * d.x + r[1][0] * d.y + r[2][0] * d.z),
                              (float)(r[0][1] * d.x + r[1][1] * d.y + r[2][1] * d.z),
                              (float)(r[0][2] * d.x + r[1][2] * d.y + r[2][2] * d.z) });
}

Vector3 Planetoid::dirToWorld(Vector3 localDir) const {
    double r[3][3];
    BodyRotationD(rotation, r);
    const Vector3 d = localDir;
    return Vector3Normalize({ (float)(r[0][0] * d.x + r[0][1] * d.y + r[0][2] * d.z),
                              (float)(r[1][0] * d.x + r[1][1] * d.y + r[1][2] * d.z),
                              (float)(r[2][0] * d.x + r[2][1] * d.y + r[2][2] * d.z) });
}

Quaternion Planetoid::rotationQuat() const {
    return QuaternionFromMatrix(MatrixRotateXYZ({ DEG2RAD * rotation.x, DEG2RAD * rotation.y, DEG2RAD * rotation.z }));
}

bool Planetoid::useFarSphere(const RenderContext& ctx) const {
    if (!farSphere || !ctx.farShader || !farSphere->isReady()) return false;
    return Vector3Distance(position, ctx.cameraPos) / std::max(1.0f, lodSettings.zoom) > kFarSphereRadii * radius * scale;
}

void Planetoid::draw(const RenderContext& ctx) {
    if (!isActive) return;
    for (const auto& objPtr : children) objPtr->draw(ctx);
    if (farSphere && ctx.farShader) farSphere->ready(); // uploads the bake once it is done
    if (useFarSphere(ctx)) { farSphere->draw(*this, ctx, *ctx.farShader); return; }
    if (plots && ctx.partShader) {
        const Vector3 c = (worldPos - ctx.renderOrigin).toFloat();
        const float occluder[4] = { c.x, c.y, c.z, (kind == BodyKind::Planet || kind == BodyKind::Moon) ? std::max(0.0f, radius - 60.0f) * scale : 0.0f };
        SetShaderValue(*ctx.partShader, ctx.partOccluderLoc, occluder, SHADER_UNIFORM_VEC4);
        if (ctx.bindSunLighting && !plots->plots.empty()) {
            // Shadow cascades on texture slots 6-7 (DrawMesh's default material only binds slot 0).
            rlDrawRenderBatchActive();
            rlEnableShader(ctx.partShader->id);
            ctx.bindSunLighting(*ctx.partShader);
        }
        for (auto& pl : plots->plots) pl->draw(*this, ctx, *ctx.partShader);
    }
    if (drawList.empty() || !ctx.terrainMaterial) return;
    const ChunkPlacer placer(*this, ctx.renderOrigin);
    ChunkBatch batch(*ctx.terrainMaterial, ctx.terrainLayers, ctx.terrainLayersLoc); // bound once for all chunks
    // Planets use the texture-layer palette (dirt, sand, grass, ...); other bodies the scene's.
    const int palette = kind == BodyKind::Planet ? 2 : ctx.scenePalette;
    batch.setUniform(ctx.paletteLoc, &palette, SHADER_UNIFORM_INT);
    // Planets and moons shadow their own night side at any distance (the shadow map is local): a
    // sphere 60 m below sea level blocks the sun. Other bodies clear it (the uniform persists).
    float occluder[4] = { 0, 0, 0, 0 };
    if (kind == BodyKind::Planet || kind == BodyKind::Moon) {
        const Vector3 c = (worldPos - ctx.renderOrigin).toFloat();
        occluder[0] = c.x; occluder[1] = c.y; occluder[2] = c.z;
        occluder[3] = std::max(0.0f, radius - 60.0f) * scale;
    }
    batch.setUniform(ctx.sunOccluderLoc, occluder, SHADER_UNIFORM_VEC4);
    RenderContext c2 = ctx;
    c2.batch = &batch;
    for (Chunk* c : drawList) c->drawChunk(c2, placer.model(c->key));
}

void Planetoid::drawTransparent(const RenderContext& ctx) {
    if (!isActive || drawList.empty() || !ctx.waterMaterial || field->seaLevelRadius() <= 0.0f || useFarSphere(ctx)) return;
    const ChunkPlacer placer(*this, ctx.renderOrigin);
    ChunkBatch batch(*ctx.waterMaterial);
    const Vector3 c = (worldPos - ctx.renderOrigin).toFloat();
    const float occluder[4] = { c.x, c.y, c.z, std::max(0.0f, radius - 60.0f) * scale }; // night side, as for terrain
    batch.setUniform(ctx.waterOccluderLoc, occluder, SHADER_UNIFORM_VEC4);
    if (ctx.bindSunLighting) ctx.bindSunLighting(ctx.waterMaterial->shader); // cliffs and builds shade the water (slots 6-7)
    RenderContext c2 = ctx;
    c2.batch = &batch;
    // Nearest first: with the stencil set up by the scene (one water layer per pixel), the closest
    // surface claims the pixel and overlapping LOD copies of the sea behind it are skipped.
    std::vector<std::pair<float, Chunk*>> order;
    for (Chunk* ch : drawList) {
        if (ch->waterMesh.vaoId == 0) continue;
        const Matrix m = placer.model(ch->key);
        const float h = ch->key.size() * 0.5f;
        order.push_back({ Vector3DistanceSqr(Vector3Transform({ h, h, h }, m), ctx.cameraPos), ch });
    }
    std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    for (auto& o : order) o.second->drawWater(c2, placer.model(o.second->key));
    Fluid::Draw(*this, ctx, batch); // flowing water (after the static surfaces it continues)
}

void Planetoid::drawFoliage(const RenderContext& ctx) {
    if (!isActive || drawList.empty() || !ctx.grassShader) return;
    const ChunkPlacer placer(*this, ctx.renderOrigin);
    const Vector3 c = (worldPos - ctx.renderOrigin).toFloat();
    const float occluder[4] = { c.x, c.y, c.z, std::max(0.0f, radius - 60.0f) * scale }; // night side, as for terrain
    rlSetUniform(ctx.grassOccluderLoc, occluder, RL_SHADER_UNIFORM_VEC4, 1);
    for (Chunk* ch : drawList) if (ch->key.stride() <= 1.0f) ch->drawGrass(ctx, placer.model(ch->key)); // the levels that carry grass (foliage.cpp)
}

void Planetoid::drawDepthOnly(const RenderContext& ctx) {
    if (!isActive) return;
    for (const auto& objPtr : children) objPtr->drawDepthOnly(ctx);
    if (plots && ctx.depthMaterial) for (auto& pl : plots->plots) pl->drawDepth(*this, ctx, *ctx.depthMaterial); // parts cast shadows
    if (drawList.empty() || !ctx.depthMaterial) return;
    const ChunkPlacer placer(*this, ctx.renderOrigin);
    ChunkBatch batch(*ctx.depthMaterial);
    RenderContext c2 = ctx;
    c2.batch = &batch;
    // Pre-cull against the light box's column in body-local space: most chunks of the draw list lie far
    // outside it, and working out their model matrix (double precision) was most of the pass's CPU time.
    const bool column = ctx.cullAxisRadius > 0.0f;
    const Vector3 pL = column ? toLocal(ctx.cullAxisPoint) : Vector3{}, dL = column ? dirToLocal(ctx.cullAxisDir) : Vector3{};
    const float rL = column ? ctx.cullAxisRadius / std::max(scale, 1e-6f) : 0.0f;
    for (Chunk* c : drawList) {
        if (column) {
            const float h = c->key.size() * 0.5f;
            const Vector3 v = Vector3Subtract(Vector3AddValue(c->key.origin(), h), pL);
            const float along = Vector3DotProduct(v, dL), reach = rL + h * 1.75f;
            if (Vector3LengthSqr(v) - along * along > reach * reach) continue;
        }
        c->drawChunkDepth(c2, placer.model(c->key));
    }
}

// ---------------------------------------------------------------------------------------------
// Node requests / results
// ---------------------------------------------------------------------------------------------

static void RequestNode(Planetoid& p, const NodeKey& key, float priority) {
    const uint32_t version = ++p.requestedVersion[key];
    if (p.nodes.find(key) == p.nodes.end()) {
        auto chunk = std::make_unique<Chunk>();
        chunk->key = key;
        chunk->parent = &p;
        p.nodes.emplace(key, std::move(chunk));
    }
    std::shared_ptr<const DensityField> field = p.field;
    std::shared_ptr<const TerrainDeltas> edits = p.edits;
    const uint64_t id = p.id;
    std::shared_ptr<std::atomic<uint32_t>> wanted = p.nodes.at(key)->wanted;
    wanted->store(version);
    Jobs().submit(priority, [field, edits, id, key, version, wanted] {
        if (wanted->load() != version) return; // superseded or freed while queued
        PROFILE_SCOPE("Build chunk");
        ChunkMeshData data;
        data.planetoidId = id;
        data.key = key;
        data.version = version;
        BuildChunkMesh(*field, *edits, data);
        std::lock_guard<std::mutex> lock(resultsMutex);
        results.push_back(std::move(data));
    });
}

static float DistanceToBox(Vector3 p, Vector3 lo, Vector3 hi) {
    const float dx = std::max({ lo.x - p.x, 0.0f, p.x - hi.x });
    const float dy = std::max({ lo.y - p.y, 0.0f, p.y - hi.y });
    const float dz = std::max({ lo.z - p.z, 0.0f, p.z - hi.z });
    return sqrtf(dx * dx + dy * dy + dz * dz);
}

// True if the whole box lies behind the opaque sphere of radius ro (centred on the body origin) as
// seen from cam: every corner inside the sphere's shadow cone and beyond its tangent plane (that
// region is convex, so the box between the corners is hidden too).
static bool BehindHorizon(Vector3 cam, Vector3 lo, Vector3 hi, float ro) {
    if (ro <= 0.0f) return false;
    const Vector3 cv = Vector3Scale(cam, 1.0f / ro);
    const float vhMagSq = Vector3DotProduct(cv, cv) - 1.0f;
    if (vhMagSq <= 0.0f) return false; // camera inside the occluder
    for (int i = 0; i < 8; ++i) {
        const Vector3 c = { (i & 1) ? hi.x : lo.x, (i & 2) ? hi.y : lo.y, (i & 4) ? hi.z : lo.z };
        const Vector3 vt = Vector3Subtract(Vector3Scale(c, 1.0f / ro), cv);
        const float vtDotVc = -Vector3DotProduct(vt, cv);
        if (!(vtDotVc > vhMagSq && vtDotVc * vtDotVc / Vector3DotProduct(vt, vt) > vhMagSq)) return false;
    }
    return true;
}

// Octree descent from the 8 root nodes: split while the camera is closer than lodFactor * node size.
// The camera as the LOD sees it (body-local). With a zoomed view (telescope), nodes inside the view
// cone are refined as if they were `zoom` times closer.
struct LodView {
    Vector3 cam;
    Vector3 dir = { 0, 0, 1 };  // view direction (body-local)
    float zoom = 1.0f;
    float halfAngle = 0.0f;     // view cone half-angle (rad), to the screen corners
};

static float LodDistance(const LodView& v, Vector3 lo, Vector3 hi) {
    const float d = DistanceToBox(v.cam, lo, hi);
    if (v.zoom <= 1.0f || d <= 0.0f) return d;
    const Vector3 to = Vector3Subtract(Vector3Scale(Vector3Add(lo, hi), 0.5f), v.cam);
    const float len = Vector3Length(to), r = 0.87f * (hi.x - lo.x);
    if (len <= r) return d / v.zoom;
    const float ang = acosf(std::clamp(Vector3DotProduct(to, v.dir) / len, -1.0f, 1.0f));
    return ang <= v.halfAngle + asinf(std::min(1.0f, r / len)) ? d / v.zoom : d;
}

static void CollectDesired(const Planetoid& p, const NodeKey& key, const LodView& v, std::vector<NodeKey>& out) {
    const Vector3 lo = key.origin();
    const Vector3 hi = Vector3AddValue(lo, key.size());
    if (p.flat) {
        if (hi.y < p.minHeight - 8.0f || lo.y > p.maxHeight + 8.0f) return;      // outside the terrain band
    } else {
        // No surface possible: wholly above the highest terrain (+ room for built material), or
        // wholly inside the lowest possible surface (a planet's mantle).
        if (DistanceToBox({ 0, 0, 0 }, lo, hi) > p.field->maxRadius() + 24.0f) return;
        const float far = Vector3Length({ std::max(fabsf(lo.x), fabsf(hi.x)), std::max(fabsf(lo.y), fabsf(hi.y)), std::max(fabsf(lo.z), fabsf(hi.z)) });
        if (far < p.field->minRadius()) return;
        if (BehindHorizon(v.cam, lo, hi, p.field->occluderRadius() * p.scale)) return; // over the horizon
    }
    if (key.level > 0 && LodDistance(v, lo, hi) < lodSettings.lodFactor * key.size()) {
        for (int i = 0; i < 8; ++i)
            CollectDesired(p, { key.level - 1, key.x * 2 + (i & 1), key.y * 2 + ((i >> 1) & 1), key.z * 2 + ((i >> 2) & 1) }, v, out);
        return;
    }
    out.push_back(key);
}

// Decide which nodes to build and draw. Old nodes stay on screen until every desired node that
// overlaps them is ready, so LOD changes never open holes.
static void UpdateLod(Planetoid& p, Vector3 cameraRender) {
    PROFILE_FUNCTION();
    const Vector3 camLocal = p.toLocal(cameraRender);
    LodView view{ camLocal };
    if (lodSettings.zoom > 1.0f) {
        view.dir = p.dirToLocal(lodSettings.zoomDir);
        view.zoom = lodSettings.zoom;
        view.halfAngle = lodSettings.zoomHalfAngle;
    }
    // Nothing to do while the camera (in body space) stays within 1 m of the last selection and no
    // results or edits arrived: the smallest node splits at lodFactor * 32 m, so 1 m is invisible.
    if (!p.lodDirty && Vector3DistanceSqr(camLocal, p.lodCamLocal) < 1.0f) return;
    p.lodDirty = false;
    p.lodCamLocal = camLocal;
    std::vector<NodeKey> desired;
    if (p.flat) {
        // Open world: a grid of root tiles around the camera, spanning the terrain height band.
        const float rs = kVoxelSize * (float)(CHUNK_SIZE << p.rootLevel);
        const int cx = (int)floorf(camLocal.x / rs), cz = (int)floorf(camLocal.z / rs);
        const int y0 = (int)floorf(p.minHeight / rs), y1 = (int)floorf(p.maxHeight / rs);
        for (int z = cz - 2; z <= cz + 2; ++z)
            for (int x = cx - 2; x <= cx + 2; ++x)
                for (int y = y0; y <= y1; ++y) CollectDesired(p, { p.rootLevel, x, y, z }, view, desired);
    } else {
        for (int i = 0; i < 8; ++i)
            CollectDesired(p, { p.rootLevel, -1 + (i & 1), -1 + ((i >> 1) & 1), -1 + ((i >> 2) & 1) }, view, desired);
    }
    std::unordered_set<NodeKey> desiredSet(desired.begin(), desired.end());

    for (const NodeKey& k : desired) {
        if (p.requestedVersion.find(k) != p.requestedVersion.end()) continue;
        const Vector3 lo = k.origin();
        // Nearest and coarsest first: fill the view quickly, then refine.
        const float dist = LodDistance(view, lo, Vector3AddValue(lo, k.size()));
        RequestNode(p, k, dist / (1.0f + k.level));
    }

    // Previously drawn nodes that are no longer desired, kept while a desired overlap is pending.
    // Overlap tests walk ancestors through hash sets (comparing every drawn node with every desired
    // one cost ~3 ms a frame with thousands of nodes while flying).
    auto parentOf = [](const NodeKey& k) { return NodeKey{ k.level + 1, k.x >> 1, k.y >> 1, k.z >> 1 }; };
    std::unordered_set<NodeKey> pendingDesired, pendingWithin; // pending nodes, and every node containing one
    for (const NodeKey& d : desired) {
        if (p.nodes.at(d)->ready) continue;
        pendingDesired.insert(d);
        for (NodeKey k = d; k.level <= p.rootLevel; k = parentOf(k))
            if (!pendingWithin.insert(k).second) break; // its ancestors are in already
    }
    std::vector<Chunk*> retained;
    std::unordered_set<NodeKey> retainedSet;
    for (Chunk* c : p.drawList) {
        if (desiredSet.count(c->key)) continue;
        bool keep = pendingWithin.count(c->key) > 0;                // a pending desired node lies inside it
        for (NodeKey k = parentOf(c->key); !keep && k.level <= p.rootLevel; k = parentOf(k))
            keep = pendingDesired.count(k) > 0;                      // inside a pending desired node
        if (keep) { retained.push_back(c); retainedSet.insert(c->key); }
    }
    std::vector<Chunk*> draw = retained;
    for (const NodeKey& d : desired) {
        Chunk* c = p.nodes.at(d).get();
        if (!c->ready) continue;
        bool coveredByCoarser = false;
        if (!retainedSet.empty())
            for (NodeKey k = parentOf(d); !coveredByCoarser && k.level <= p.rootLevel; k = parentOf(k)) coveredByCoarser = retainedSet.count(k) > 0;
        if (!coveredByCoarser) draw.push_back(c);
    }
    // Nearest first (front to back): the depth test then rejects most of the terrain hidden behind
    // nearer hills before the (costly) terrain shader runs on it. Unsorted, much of the screen was
    // shaded two or three times over: about a quarter of the frame.
    {
        std::vector<std::pair<float, Chunk*>> byDist;
        byDist.reserve(draw.size());
        for (Chunk* c : draw) { const Vector3 lo = c->key.origin(); byDist.push_back({ DistanceToBox(camLocal, lo, Vector3AddValue(lo, c->key.size())), c }); }
        std::sort(byDist.begin(), byDist.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        for (size_t i = 0; i < draw.size(); ++i) draw[i] = byDist[i].second;
    }
    p.drawList = std::move(draw);

    // Clip the one-cell overlap on axes where the +neighbour is drawn at the same level (it would
    // duplicate that neighbour's surface and z-fight); keep it toward other levels to hide cracks.
    std::unordered_set<NodeKey> drawn;
    for (Chunk* c : p.drawList) drawn.insert(c->key);
    for (Chunk* c : p.drawList) {
        const NodeKey k = c->key;
        const Vector3 o = k.origin();
        // Keep a quarter cell past the shared face: each chunk has its own model matrix, so vertices
        // on the face can land a hair apart; both build this strip from the same samples, so the
        // duplicated sliver is identical surface and does not z-fight.
        const float s = k.size() + 0.25f * (float)k.stride();
        c->clipMax.x = drawn.count({ k.level, k.x + 1, k.y, k.z }) ? o.x + s : 1e30f;
        c->clipMax.y = drawn.count({ k.level, k.x, k.y + 1, k.z }) ? o.y + s : 1e30f;
        c->clipMax.z = drawn.count({ k.level, k.x, k.y, k.z + 1 }) ? o.z + s : 1e30f;
    }

    // Free nodes that are neither desired nor on screen (late results for them are dropped).
    std::unordered_set<const Chunk*> onScreen(p.drawList.begin(), p.drawList.end());
    for (auto it = p.nodes.begin(); it != p.nodes.end();) {
        if (!desiredSet.count(it->first) && !onScreen.count(it->second.get())) {
            p.requestedVersion.erase(it->first);
            it = p.nodes.erase(it);
        } else {
            ++it;
        }
    }
}

static void IntegrateResult(const ChunkMeshData& r) {
    Planetoid* p = Planetoid::find(r.planetoidId);
    if (!p) return;                                                        // planetoid was destroyed
    auto v = p->requestedVersion.find(r.key);
    if (v == p->requestedVersion.end() || v->second != r.version) return;  // superseded or freed
    auto it = p->nodes.find(r.key);
    if (it == p->nodes.end()) return;

    buildMsTotal += r.buildMs;
    sampleMsTotal += r.sampleMs;
    buildCount++;
    if (r.overflow) logger.logf("[ERROR] Node L%d (%d,%d,%d) exceeds 65535 vertices; skipped\n", r.key.level, r.key.x, r.key.y, r.key.z);
    it->second->setMesh(r);
    p->lodDirty = true; // readiness changed: the draw list may change
}

static bool PipelineIdle() {
    std::lock_guard<std::mutex> lock(resultsMutex);
    return Jobs().idle() && results.empty() && pendingUploads.empty();
}

static bool AnyLodPending(Scene& world) {
    for (auto& obj : world.objects)
        if (auto* p = dynamic_cast<Planetoid*>(obj.get()))
            for (auto& n : p->nodes) if (!n.second->ready) return true;
    return false;
}

void SetLodZoom(Scene& world, float zoom, Vector3 dir, float halfAngle) {
    zoom = std::max(1.0f, zoom);
    LodSettings& s = lodSettings;
    // Re-select when the zoom changes, or the aim moves a quarter of the view (finer is wasted work).
    const bool changed = fabsf(zoom - s.zoom) > 0.01f * zoom ||
                         (zoom > 1.0f && Vector3DotProduct(dir, s.zoomDir) < cosf(std::max(halfAngle, 1e-3f) * 0.25f));
    s.zoom = zoom;
    s.zoomHalfAngle = halfAngle;
    if (!changed) return;
    s.zoomDir = dir;
    for (auto& obj : world.objects)
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) p->lodDirty = true;
}

namespace { double worldTime = 0.0; }
double WorldTime() { return worldTime; }
void SetWorldTime(double t) { worldTime = t; }

// Place an orbiting body for the current world time (its parent must already be placed).
void PlaceBodyOnOrbit(Planetoid& p) {
    if (!p.orbitParent || Planetoid::find(p.orbitParentId) != p.orbitParent) { p.orbitParent = nullptr; return; }
    const Planetoid& parent = *p.orbitParent;
    const double a = p.orbitPhase + p.orbitRate * worldTime;
    const double c = cos(a), s = sin(a);
    p.worldPos = parent.worldPos + DVec3(c * p.orbitRadius, p.orbitHeight, s * p.orbitRadius);
    const double v = p.orbitRate * p.orbitRadius;
    p.linearVelocity = Vector3Add(parent.linearVelocity, { (float)(-s * v), 0.0f, (float)(c * v) });
}

void advanceBodies(Scene& world, float dt) {
    // Runs at the start of the frame, before the player/ship update, so riders follow this
    // frame's transform (spinning here after they moved left them one frame behind the surface).
    // Orbits are evaluated from the world clock (not integrated), so they never drift.
    worldTime += dt;
    bool moved = false;
    for (auto& obj : world.objects) {
        auto* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p) continue;
        if (p->spin.x != 0.0f || p->spin.y != 0.0f || p->spin.z != 0.0f) {
            float* r[3] = { &p->rotation.x, &p->rotation.y, &p->rotation.z };
            const float s[3] = { p->spin.x, p->spin.y, p->spin.z };
            for (int k = 0; k < 3; ++k) {
                // Something else (a load, a script) set the angle: take it over.
                if (!p->spinAngleSet || (float)p->spinAngle[k] != *r[k]) p->spinAngle[k] = *r[k];
                double a = fmod(p->spinAngle[k] + (double)s[k] * dt, 360.0);
                if (a < 0.0) a += 360.0;
                p->spinAngle[k] = a;
                *r[k] = (float)a;
            }
            p->spinAngleSet = true;
        }
        if (p->orbitParent) { PlaceBodyOnOrbit(*p); moved = true; } // parents come first in the list
    }
    if (moved) world.syncRenderPositions();
}

void worldUpdate(Scene& world, float uploadBudgetMs) {
    PROFILE_FUNCTION();
    // Ring rocks near the camera become voxel bodies (and back) before this frame's LOD selection.
    for (auto& ring : Rings()) ring->update(world, world.renderOrigin + world.camera.position);

    {
        std::lock_guard<std::mutex> lock(resultsMutex);
        for (ChunkMeshData& r : results) pendingUploads.push_back(std::move(r));
        results.clear();
    }

    // GPU uploads are the only main-thread cost; cap them per frame so streaming never hitches.
    const auto t0 = std::chrono::steady_clock::now();
    int uploads = 0;
    while (!pendingUploads.empty()) {
        IntegrateResult(pendingUploads.front());
        pendingUploads.pop_front();
        uploads++;
        if (std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count() > uploadBudgetMs) break;
    }

    // LOD after uploads, so this frame's draw list already includes nodes that just finished.
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            UpdateLod(*p, world.camera.position);
            if (p->trees) p->trees->update(*p, p->toLocal(world.camera.position));
        }
    }

    stats.queuedJobs = (int)Jobs().queued();
    stats.runningJobs = Jobs().running();
    stats.pendingUploads = (int)pendingUploads.size();
    stats.uploadsLastFrame = uploads;
    stats.workerThreads = (int)Jobs().threadCount();
    stats.avgBuildMs = buildCount ? (float)(buildMsTotal / buildCount) : 0.0f;
    stats.avgSampleMs = buildCount ? (float)(sampleMsTotal / buildCount) : 0.0f;
    stats.chunks = stats.nodes = stats.edits = 0;
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            stats.nodes += (int)p->nodes.size();
            for (Chunk* c : p->drawList) stats.chunks += c->hasMesh() ? 1 : 0;
            stats.edits += (int)p->edits->strokes;
        }
    }

    if (generationStart >= 0.0 && PipelineIdle() && !AnyLodPending(world)) {
        stats.lastGenerationSeconds = (float)(GetTime() - generationStart);
        logger.logf("[World] ready in %.2f s (%d meshes drawn, %d nodes, avg build %.2f ms on %d threads)\n",
                    stats.lastGenerationSeconds, stats.chunks, stats.nodes, stats.avgBuildMs, stats.workerThreads);
        TraceLog(LOG_WARNING, "World ready in %.2f s (%d meshes)", stats.lastGenerationSeconds, stats.chunks);
        generationStart = -1.0;
    }
}

// False until a (re)generation has fully streamed in, and while edits/LOD changes are building.
bool WorldIdle() { return generationStart < 0.0 && PipelineIdle(); }

const WorldStats& GetWorldStats() { return stats; }

// ---------------------------------------------------------------------------------------------
// Raycast + edits
// ---------------------------------------------------------------------------------------------

bool RaycastTerrain(Scene& world, Ray ray, float maxDistance, TerrainHit& hit) {
    PROFILE_FUNCTION();
    bool found = false;
    float best = maxDistance;
    for (auto& obj : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p) continue;
        // Work in planetoid-local space; step lengths are in local voxels.
        const Vector3 o = p->toLocal(ray.position);
        const Vector3 farPoint = p->toLocal(Vector3Add(ray.position, Vector3Scale(ray.direction, best)));
        const float localLen = Vector3Distance(o, farPoint);
        if (localLen <= 0.0f) continue;
        const Vector3 d = Vector3Scale(Vector3Subtract(farPoint, o), 1.0f / localLen);

        // Bounding sphere (highest terrain + room for added material).
        float t = 0.0f, tEnd = localLen;
        if (!p->flat) {
            const float R = p->field->maxRadius() + 16.0f;
            const float b = Vector3DotProduct(o, d);
            const float c = Vector3DotProduct(o, o) - R * R;
            const float disc = b * b - c;
            if (disc < 0.0f) continue;
            t = std::max(0.0f, -b - sqrtf(disc));
            tEnd = std::min(localLen, -b + sqrtf(disc));
        }

        // Every field has the same density gradient (DENSITY_GRADIENT per metre across a surface), so
        // the density says roughly how far the surface is: step half that gap (0.5..64 m), which
        // crosses kilometres of air above a planet in a few dozen samples instead of thousands.
        float prevT = t;
        float dens = p->densityAt(Vector3Add(o, Vector3Scale(d, t)));
        if (dens >= ISO_LEVEL) continue; // starting inside rock
        for (; t <= tEnd; prevT = t, t += std::clamp((ISO_LEVEL - dens) / DENSITY_GRADIENT * 0.5f, 0.5f, 64.0f)) {
            dens = p->densityAt(Vector3Add(o, Vector3Scale(d, t)));
            if (dens >= ISO_LEVEL) {
                float lo = prevT, hi = t; // refine the crossing
                for (int i = 0; i < 8; ++i) {
                    const float mid = 0.5f * (lo + hi);
                    if (p->densityAt(Vector3Add(o, Vector3Scale(d, mid))) >= ISO_LEVEL) hi = mid; else lo = mid;
                }
                const Vector3 local = Vector3Add(o, Vector3Scale(d, hi));
                const Vector3 worldPoint = p->toWorld(local);
                const float dist = Vector3Distance(ray.position, worldPoint);
                if (dist < best) {
                    best = dist;
                    found = true;
                    hit.planetoid = p;
                    hit.local = local;
                    hit.point = worldPoint;
                    hit.distance = dist;
                    hit.normal = Vector3Normalize(Vector3Subtract(p->toWorld(Vector3Add(local, p->field->normal(local, p->edits.get()))), worldPoint));
                    hit.material = p->materialAt(local);
                }
                break;
            }
        }
    }
    return found;
}

void ApplyTerrainEdit(Planetoid& p, const TerrainEdit& edit) { ApplyTerrainEdits(p, &edit, 1); }

void ApplyTerrainEdits(Planetoid& p, const TerrainEdit* edits, size_t count) {
    if (count == 0) return;
    // Fold the strokes into the delta lattice: a new snapshot sharing every untouched chunk (the
    // old one may still be in use by worker builds).
    TerrainDeltasBuilder b(*p.edits);
    std::unordered_set<TerrainDeltas::Key> touchedChunks;
    for (size_t i = 0; i < count; ++i) b.apply(edits[i], *p.field, &touchedChunks);
    p.edits = b.finish();
    for (TerrainDeltas::Key k : touchedChunks) p.dirtyRegions.insert(TerrainDeltas::RegionOf(k));
    p.lodDirty = true;

    // Rebuild every existing node (any level) whose samples or normals the strokes reach. Nodes
    // that do not exist yet are built later with the new deltas anyway.
    Vector3 lo = edits[0].center, hi = edits[0].center;
    float rmax = 0.0f;
    for (size_t i = 0; i < count; ++i) { lo = Vector3Min(lo, edits[i].center); hi = Vector3Max(hi, edits[i].center); rmax = std::max(rmax, edits[i].radius); }
    std::vector<NodeKey> touched;
    for (auto& entry : p.nodes) {
        const NodeKey& k = entry.first;
        const Vector3 nlo = k.origin();
        const Vector3 nhi = Vector3AddValue(nlo, k.size() * (CHUNK_SIZE + 1.0f) / CHUNK_SIZE);
        const float margin = rmax + k.stride() + 1.0f;
        if (hi.x + margin < nlo.x || lo.x - margin > nhi.x || hi.y + margin < nlo.y || lo.y - margin > nhi.y || hi.z + margin < nlo.z || lo.z - margin > nhi.z) continue;
        if (count == 1 && DistanceToBox(edits[0].center, nlo, nhi) > margin) continue;
        touched.push_back(k);
    }
    for (const NodeKey& k : touched) RequestNode(p, k, -1.0f - (float)k.level * 0.01f); // ahead of background work
    if (p.trees)
        for (size_t i = 0; i < count; i += std::max<size_t>(1, count / 64)) p.trees->onEdit(edits[i].center, edits[i].radius, p.edits);
    // Flowing water: wake it, or start it where the edit opens space next to water.
    for (size_t i = 0; i < count; i += std::max<size_t>(1, count / 64)) Fluid::OnTerrainEdit(p, edits[i].center, edits[i].radius);
}

void RequestRemesh(Planetoid& p, Vector3 c, float r) {
    std::vector<NodeKey> touched;
    for (auto& entry : p.nodes) {
        const NodeKey& k = entry.first;
        const Vector3 nlo = k.origin();
        const Vector3 nhi = Vector3AddValue(nlo, k.size() * (CHUNK_SIZE + 1.0f) / CHUNK_SIZE);
        if (DistanceToBox(c, nlo, nhi) > r + k.stride() + 1.0f) continue;
        touched.push_back(k);
    }
    for (const NodeKey& k : touched) RequestNode(p, k, -1.0f - (float)k.level * 0.01f);
}

// ---------------------------------------------------------------------------------------------
// Shared physics helpers
// ---------------------------------------------------------------------------------------------

// Approximate distance from world point c to the iso-surface of p (world units, + outside), with the
// outward surface normal. Treats density as locally linear: (ISO - density) / |gradient|.
bool TerrainSurfaceDistance(Planetoid& p, Vector3 c, float& dist, Vector3& normal) {
    const Vector3 l = p.toLocal(c);
    const float d = p.densityAt(l);
    const float e = 0.25f;
    const Vector3 g = {
        (p.densityAt({ l.x + e, l.y, l.z }) - p.densityAt({ l.x - e, l.y, l.z })) / (2 * e),
        (p.densityAt({ l.x, l.y + e, l.z }) - p.densityAt({ l.x, l.y - e, l.z })) / (2 * e),
        (p.densityAt({ l.x, l.y, l.z + e }) - p.densityAt({ l.x, l.y, l.z - e })) / (2 * e),
    };
    const float glen = Vector3Length(g);
    if (glen < 1e-5f) return false;
    dist = (ISO_LEVEL - d) / glen * p.scale;
    normal = Vector3Normalize(Vector3Subtract(p.toWorld(Vector3Add(l, Vector3Scale(g, -1.0f / glen))), c));
    return true;
}

Planetoid* NearestPlanetoid(Scene& world, Vector3 p) {
    Planetoid* best = nullptr;
    float bestDist = 1e30f;
    for (auto& obj : world.objects) {
        Planetoid* b = dynamic_cast<Planetoid*>(obj.get());
        if (!b || b->vehicle) continue; // a rocket in flight has no ground
        // Distance to its mean surface (gravityRadius: an asteroid's is ~0.72 of its bounding radius): by the
        // bounding radius a big rock 200 m off "won" over the small one underfoot, and the player fell
        // through it at the difference of their orbits.
        const float surface = b->gravityRadius > 0.0f ? b->gravityRadius : b->radius * b->scale;
        const float d = b->flat ? -1e30f : Vector3Distance(p, b->position) - surface; // terrain is always "under you"
        if (d < bestDist) { bestDist = d; best = b; }
    }
    return best;
}

Vector3 BodyGravityAt(const Planetoid& b, Vector3 p) {
    if (b.flat) return { 0.0f, -b.surfaceGravity, 0.0f }; // open-world terrain: uniform, straight down
    const Vector3 toC = Vector3Subtract(b.position, p);
    const float d = Vector3Length(toC);
    if (d < 1e-4f) return { 0, 0, 0 };
    const float R = b.gravityRadius > 0.0f ? b.gravityRadius : b.radius * b.scale;
    // Uniform-density sphere of the body's mean surface radius: inverse-square outside, falling
    // linearly to zero at the centre inside (continuous at R, where it equals surfaceGravity).
    const float g = d >= R ? b.surfaceGravity * (R * R) / (d * d) : b.surfaceGravity * d / R;
    return Vector3Scale(toC, g / d);
}

Vector3 GravityAt(Scene& world, Vector3 p) {
    Vector3 g = { 0, 0, 0 };
    for (auto& obj : world.objects)
        if (const Planetoid* b = dynamic_cast<const Planetoid*>(obj.get())) g = Vector3Add(g, BodyGravityAt(*b, p));
    return g;
}


// ---------------------------------------------------------------------------------------------
// Generation
// ---------------------------------------------------------------------------------------------

static std::vector<std::unique_ptr<RingSystem>> gRings;
std::vector<std::unique_ptr<RingSystem>>& Rings() { return gRings; }

void worldHandler(Scene& world) {
    // Rings draw through the scene's extra-geometry hooks (band in the far pass, rocks in the others).
    world.farExtraRange = [&world](double& dmin, double& dmax) {
        bool any = false;
        dmin = 1e300; dmax = 0.0;
        for (auto& ring : gRings) {
            double a, b;
            if (ring->bandRange(world, a, b)) { any = true; dmin = std::min(dmin, a); dmax = std::max(dmax, b); }
        }
        return any;
    };
    world.drawFarExtra = [&world](const RenderContext& ctx) { for (auto& ring : gRings) ring->drawBand(world, ctx); };
    world.drawPassExtra = [&world](const RenderContext& ctx, bool rebuild) { for (auto& ring : gRings) ring->drawRocks(world, ctx, rebuild); };
    world.camera.position = { 0.0f, 10.0f, 10.0f };
    world.camera.target = { 0.0f, 0.0f, 0.0f };
    world.camera.up = { 0.0f, 1.0f, 0.0f };
    world.camera.fovy = 45.0f;
    world.camera.projection = CAMERA_PERSPECTIVE;

    world.loadTerrainTextures();
    generateWorld(world);
    loadWorld(world);
    // Start out looking at the first planetoid.
    if (!world.objects.empty()) world.camera.target = world.objects[0]->position;
}

bool ReloadSystemDef() {
    SystemDef sys;
    const bool ok = LoadSystemDef(worldGen.bodiesPath, sys);
    if (!ok) {
        if (worldGen.systemLoaded) {
            logger.logf("[Bodies] keeping the current system\n");
            return false;
        }
        logger.logf("[Bodies] falling back to the built-in start system\n");
        sys = LegacySystemDef(8, 50, 1500.0f, 3.0f, { 64.0f, 0.0f, 0.0f });
    }
    // Command-line overrides (bench): first single body and first ring.
    BodyDef* start = nullptr;
    BodyDef* ring = nullptr;
    int singles = 0;
    for (BodyDef& b : sys.bodies) {
        if (b.kind == BodyKind::Ring) { if (!ring) ring = &b; }
        else { ++singles; if (!start) start = &b; }
    }
    if (start && worldGen.size > 0) start->radius = (float)worldGen.size;
    if (start && worldGen.hasFirst) { start->hasPosition = true; start->position = worldGen.first; start->orbitParent.clear(); }
    if (worldGen.count >= 0) {
        const int want = std::max(0, worldGen.count - singles);
        if (ring) ring->count = want;
        else if (want > 0) {
            BodyDef r = LegacySystemDef(2, worldGen.size > 0 ? worldGen.size : 50, 1500.0f, 3.0f, {}).bodies[1];
            r.count = want;
            sys.bodies.push_back(r);
            ring = &sys.bodies.back();
        }
    }
    if (ring && worldGen.size > 0) ring->radius = (float)worldGen.size;
    if (ring && worldGen.systemRadius > 0.0f) ring->extent = worldGen.systemRadius;
    worldGen.system = sys;
    worldGen.systemLoaded = true;
    return ok;
}

void generateWorld(Scene& world) {
    PROFILE_FUNCTION();
    if (!worldGen.systemLoaded) ReloadSystemDef();
    gRings.clear();
    worldGen.lastSeed = worldGen.seed ? worldGen.seed : rd();
    gen.seed(worldGen.lastSeed);
    logger.logf("Generating world (seed %u, %zu definitions from %s)...\n", worldGen.lastSeed, worldGen.system.bodies.size(), worldGen.bodiesPath.c_str());

    // Every body draws the same values in the same order (3 position values, its terrain seed,
    // a size value for ring members, a spin value), so a seed reproduces the system exactly.
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::unordered_map<std::string, Planetoid*> byName;
    int index = 0;
    auto add = [&](const BodyDef& d, BodyKind kind, const std::string& name, DVec3 pos, float radius, unsigned int seed, float spin) {
        auto p = std::make_unique<Planetoid>(name, pos, Vector3{ 0, 0, 0 }, 1.0f, seed, radius);
        p->kind = kind;
        if (kind == BodyKind::Planet || kind == BodyKind::Moon) {
            // Planet terrain: radius is sea level, where gravity's reference surface sits too; the
            // LOD root must also cover the highest mountains.
            auto field = std::make_shared<PlanetField>(seed, radius, kind == BodyKind::Moon);
            p->gravityRadius = radius;
            p->rootLevel = 0;
            while (kVoxelSize * (float)(CHUNK_SIZE << p->rootLevel) < field->maxRadius() + 24.0f) p->rootLevel++;
            p->field = field;
            if (kind == BodyKind::Planet) {
                field->foliage = true;
                field->foliageSeed = seed ^ 0x7aee5u;
                p->trees = std::make_unique<TreeSystem>(field, field->foliageSeed);
            }
            p->farSphere = std::make_unique<FarSphere>(field, field->foliage, field->foliageSeed); // bakes on the workers
        }
        p->spin = { 0.0f, spin, 0.0f };
        p->surfaceGravity = d.gravity >= 0.0f ? d.gravity : AutoSurfaceGravity(kind, radius);
        p->parent = &world.rootObject;
        byName[name] = p.get();
        world.objects.push_back(std::move(p));
        ++index;
    };
    // Circular orbit in the parent's equatorial plane at the speed its gravity gives:
    // omega = sqrt(GM / r^3), GM = surface gravity * gravity radius^2.
    auto setOrbit = [](Planetoid* p, Planetoid* parent, double radius, double phase, double height) {
        const double gr = parent->gravityRadius > 0.0f ? parent->gravityRadius : parent->radius * parent->scale;
        const double gm = parent->surfaceGravity * gr * gr;
        p->orbitParent = parent;
        p->orbitParentId = parent->id;
        p->orbitRadius = radius;
        p->orbitHeight = height;
        const double r = sqrt(radius * radius + height * height);
        p->orbitRate = sqrt(gm / (r * r * r));
        p->orbitPhase = phase;
        PlaceBodyOnOrbit(*p);
    };
    for (const BodyDef& d : worldGen.system.bodies) {
        Planetoid* parent = nullptr;
        if (!d.orbitParent.empty()) {
            auto it = byName.find(d.orbitParent);
            if (it != byName.end()) parent = it->second;
        }
        if (d.kind == BodyKind::Ring && parent) {
            // An orbiting ring (ring.h): its rocks are instances, promoted to voxel bodies up close.
            // Seeded from the world seed and the ring's position in the file (no draws from `gen`,
            // so the bodies after it keep their seeds).
            gRings.push_back(std::make_unique<RingSystem>(d, *parent, worldGen.lastSeed ^ (uint32_t)(gRings.size() * 2654435761u), (int)gRings.size()));
            continue;
        }
        if (d.kind == BodyKind::Ring) {
            // A loose "belt": random directions, distances up to `extent`. Doubles keep these exact.
            for (int k = 0; k < d.count; ++k) {
                const double x = (unit(gen) - 0.5) * 2.0 * d.extent;
                const double y = (unit(gen) - 0.5) * (2.0 * d.flatten) * d.extent;
                const double z = (unit(gen) - 0.5) * 2.0 * d.extent;
                const unsigned int seed = (unsigned int)gen();
                const float radius = d.radius * (float)(d.sizeMin + (d.sizeMax - d.sizeMin) * unit(gen));
                const float spin = (float)(unit(gen) - 0.5) * 2.0f * d.maxSpinDegPerSec;
                add(d, BodyKind::Asteroid, d.name + "_" + std::to_string(index), { x, y, z }, radius, seed, spin);
            }
            continue;
        }
        for (int i = 0; i < 3; ++i) (void)unit(gen); // position draws (kept so seeds reproduce)
        const unsigned int drawnSeed = (unsigned int)gen();
        (void)unit(gen);                               // spin draw
        DVec3 pos = DVec3(d.position);
        add(d, d.kind, d.name, pos, d.radius, d.hasSeed ? d.seed : drawnSeed, d.spinDegPerSec);
        if (parent) setOrbit(byName[d.name], parent, parent->radius * parent->scale + d.orbitAltitude, d.orbitPhaseDeg * DEG2RAD, 0.0);
    }
    // Big bodies orbiting inside a ring clear their lanes in it.
    for (auto& ring : gRings)
        for (auto& obj : world.objects)
            if (auto* p = dynamic_cast<Planetoid*>(obj.get()); p && p->orbitParent && p->ringRock < 0) {
                const int n = ring->clearLane(p->orbitRadius, p->orbitHeight, p->radius * p->scale);
                if (n > 0) logger.logf("[World] %s clears %d rocks from ring %s\n", p->name.c_str(), n, ring->name.c_str());
            }
    world.syncRenderPositions();
    generationStart = GetTime();
    buildMsTotal = sampleMsTotal = 0.0;
    buildCount = 0;
}

void BeginWorldStreaming() { generationStart = GetTime(); buildMsTotal = sampleMsTotal = 0.0; buildCount = 0; }

void regenerateWorld(Scene& world) {
    Jobs().clearQueued();  // in-flight jobs finish; their results are dropped (planetoid ids are gone)
    Factory::Reset();      // machines belong to the old bodies
    Fluid::Reset();        // so does flowing water
    player.nearest = nullptr; // every body goes (the player refinds its ground next update)
    gRings.clear();        // before their parent bodies go
    world.objects.clear(); // releases chunk meshes
    generateWorld(world);
}

// ---------------------------------------------------------------------------------------------
// Save / load: the world is fully determined by the generation settings + each planetoid's
// edit list, so a save is tiny and replays exactly.
// ---------------------------------------------------------------------------------------------

static std::string SAVE_PATH_STR = "saves/world.sav";
#define SAVE_PATH SAVE_PATH_STR.c_str()
void SetWorldSavePath(const char* path) { SAVE_PATH_STR = path; }
const char* WorldSavePath() { return SAVE_PATH_STR.c_str(); }

// ---------------------------------------------------------------------------------------------
// Region files: a body's edited delta chunks, 16^3 chunks per file, each chunk compressed.
//   "RGN3", int32 count, then per chunk: int32 cx, cy, cz, int32 payload bytes, run-length payload
//   (RGN2, the same in 1/256 steps, and RGN1, DEFLATE-compressed, are still read).
// ---------------------------------------------------------------------------------------------

static std::string RegionDir() { return SAVE_PATH_STR + ".regions"; }

static std::string RegionFile(const std::string& dir, int body, TerrainDeltas::Key region) {
    int rx, ry, rz;
    TerrainDeltas::Unpack(region, rx, ry, rz);
    return TextFormat("%s/b%d_%d_%d_%d.rgn", dir.c_str(), body, rx, ry, rz);
}

// RGN2/3 chunk payload: run-length coded, repeating [uint16 zero run][uint16 count][count x int16].
// Edited chunks are mostly zeros around a band of changed samples; this is compact and ~100x
// faster than DEFLATE (which made a 100k-stroke save take seconds).
static void EncodeChunk(const DeltaChunk& c, std::vector<uint16_t>& out) {
    out.clear();
    const int n = kDeltaChunk * kDeltaChunk * kDeltaChunk;
    thread_local std::vector<int16_t> dense(DeltaChunk::kN);
    c.toDense(dense.data());
    const int16_t* d = dense.data();
    int i = 0;
    while (i < n) {
        int zeros = 0;
        while (i < n && d[i] == 0 && zeros < 65535) { ++i; ++zeros; }
        const int start = i;
        while (i < n && d[i] != 0 && i - start < 65535) ++i;
        out.push_back((uint16_t)zeros);
        out.push_back((uint16_t)(i - start));
        for (int k = start; k < i; ++k) out.push_back((uint16_t)d[k]);
    }
}

static bool DecodeChunk(const uint16_t* in, size_t words, DeltaChunk& c) {
    const int n = kDeltaChunk * kDeltaChunk * kDeltaChunk;
    thread_local std::vector<int16_t> dense(DeltaChunk::kN);
    std::fill(dense.begin(), dense.end(), (int16_t)0);
    int i = 0;
    size_t w = 0;
    while (w + 1 < words) {
        const int zeros = in[w++], count = in[w++];
        if (i + zeros + count > n || w + count > words) return false;
        i += zeros;
        for (int k = 0; k < count; ++k) dense[i++] = (int16_t)in[w++];
    }
    c.fromDense(dense.data());
    return true;
}

// The placed-material channel: the same run-length layout, one uint16 per material byte.
static void EncodeMaterials(const DeltaChunk& c, std::vector<uint16_t>& out) {
    out.clear();
    const int n = kDeltaChunk * kDeltaChunk * kDeltaChunk;
    int i = 0;
    bool any = false;
    for (int k = 0; k < n && !any; ++k) any = c.mat(k) != 0;
    if (!any) return; // nothing placed: an empty channel
    while (i < n) {
        int zeros = 0;
        while (i < n && c.mat(i) == 0 && zeros < 65535) { ++i; ++zeros; }
        const int start = i;
        while (i < n && c.mat(i) != 0 && i - start < 65535) ++i;
        out.push_back((uint16_t)zeros);
        out.push_back((uint16_t)(i - start));
        for (int k = start; k < i; ++k) out.push_back((uint16_t)c.mat(k));
    }
}

static bool DecodeMaterials(const uint16_t* in, size_t words, DeltaChunk& c) {
    const int n = kDeltaChunk * kDeltaChunk * kDeltaChunk;
    int i = 0;
    size_t w = 0;
    while (w + 1 < words) {
        const int zeros = in[w++], count = in[w++];
        if (i + zeros + count > n || w + count > words) return false;
        i += zeros;
        for (int k = 0; k < count; ++k) c.setMat(i++, (uint8_t)std::min<int>(in[w++], (int)VoxelMaterial::Count));
    }
    return true;
}

static bool WriteRegion(const std::string& path, const std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>>& chunks) {
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    // RGN5: RGN4 (RGN3's run-length densities, then the placed-material channel coded the same way)
    // on the kVoxelSize lattice.
    fwrite("RGN5", 1, 4, f);
    const int32_t n = (int32_t)chunks.size();
    fwrite(&n, 4, 1, f);
    std::vector<uint16_t> enc, encM;
    for (const auto& c : chunks) {
        int cx, cy, cz;
        TerrainDeltas::Unpack(c.first, cx, cy, cz);
        EncodeChunk(*c.second, enc);
        EncodeMaterials(*c.second, encM);
        const int32_t hdr[5] = { cx, cy, cz, (int32_t)(enc.size() * 2), (int32_t)(encM.size() * 2) };
        fwrite(hdr, 4, 5, f);
        fwrite(enc.data(), 2, enc.size(), f);
        fwrite(encM.data(), 2, encM.size(), f);
    }
    fclose(f);
    return true;
}

static bool ReadRegion(const std::string& path, TerrainDeltasBuilder& b) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    char magic[4];
    int32_t n = 0;
    bool ok = fread(magic, 1, 4, f) == 4 && fread(&n, 4, 1, f) == 1;
    // RGN3: run-length, deltas in 1/64 steps. RGN2 (run-length) and RGN1 (DEFLATE) stored 1/256
    // steps: scaled down on load.
    // RGN1-4 sit on the old 1 m lattice: their edits would land in the wrong place on the kVoxelSize one
    // (a save-breaking change, 2026-11-12), so they are skipped.
    if (ok && memcmp(magic, "RGN5", 4) != 0) { logger.logf("[World] %s: terrain edits on the old 1 m lattice, skipped\n", path.c_str()); fclose(f); return true; }
    const bool v4 = ok && memcmp(magic, "RGN5", 4) == 0; // + the placed-material channel
    const bool rle = ok && (v4 || memcmp(magic, "RGN3", 4) == 0 || memcmp(magic, "RGN2", 4) == 0);
    const int legacy = ok && !v4 && memcmp(magic, "RGN3", 4) != 0 ? 4 : 1;
    auto rescale = [legacy](DeltaChunk& c) { if (legacy > 1) for (int i = 0; i < DeltaChunk::kN; ++i) c.set(i, (int16_t)lroundf(c.get(i) / (float)legacy)); };
    ok = ok && (rle || memcmp(magic, "RGN1", 4) == 0);
    std::vector<unsigned char> comp;
    for (int32_t i = 0; ok && i < n; ++i) {
        int32_t hdr[5] = {};
        if (fread(hdr, 4, v4 ? 5 : 4, f) != (size_t)(v4 ? 5 : 4) || hdr[3] <= 0 || hdr[4] < 0) { ok = false; break; }
        comp.resize((size_t)hdr[3]);
        if (fread(comp.data(), 1, comp.size(), f) != comp.size()) { ok = false; break; }
        if (rle) {
            auto chunk = std::make_shared<DeltaChunk>();
            if (!DecodeChunk((const uint16_t*)comp.data(), comp.size() / 2, *chunk)) { ok = false; break; }
            if (hdr[4] > 0) {
                std::vector<unsigned char> mats((size_t)hdr[4]);
                if (fread(mats.data(), 1, mats.size(), f) != mats.size() || !DecodeMaterials((const uint16_t*)mats.data(), mats.size() / 2, *chunk)) { ok = false; break; }
            }
            rescale(*chunk);
            b.setChunk(TerrainDeltas::KeyOf(hdr[0], hdr[1], hdr[2]), chunk);
            continue;
        }
        int size = 0;
        unsigned char* raw = DecompressData(comp.data(), (int)comp.size(), &size);
        if (raw && size == (int)(DeltaChunk::kN * sizeof(int16_t))) {
            auto chunk = std::make_shared<DeltaChunk>();
            for (int i = 0; i < DeltaChunk::kN; ++i) chunk->set(i, ((const int16_t*)raw)[i]);
            rescale(*chunk);
            b.setChunk(TerrainDeltas::KeyOf(hdr[0], hdr[1], hdr[2]), chunk);
        } else ok = false;
        if (raw) MemFree(raw);
    }
    fclose(f);
    return ok;
}

// Writes a body's changed regions (all of them, after removing old files, the first time it is
// saved to this directory). Returns the number of region files written.
int SaveRegionsFor(const TerrainDeltas& edits, std::unordered_set<TerrainDeltas::Key>& dirty, bool& allDirty, std::string& savedTo, int body, const std::string& dir) {
    const bool full = allDirty || savedTo != dir;
    if (full) {
        const std::string prefix = TextFormat("b%d_", body);
        FilePathList files = LoadDirectoryFiles(dir.c_str());
        for (unsigned int i = 0; i < files.count; ++i)
            if (std::string(GetFileName(files.paths[i])).rfind(prefix, 0) == 0) std::remove(files.paths[i]);
        UnloadDirectoryFiles(files);
    }
    std::unordered_map<TerrainDeltas::Key, std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>>> byRegion;
    for (const auto& kv : edits.chunks) {
        const TerrainDeltas::Key r = TerrainDeltas::RegionOf(kv.first);
        if (full || dirty.count(r)) byRegion[r].push_back({ kv.first, kv.second.get() });
    }
    int written = 0;
    for (const auto& kv : byRegion) written += WriteRegion(RegionFile(dir, body, kv.first), kv.second) ? 1 : 0;
    dirty.clear();
    allDirty = false;
    savedTo = dir;
    return written;
}
static int SaveBodyRegions(Planetoid& p, int body, const std::string& dir) {
    return SaveRegionsFor(*p.edits, p.dirtyRegions, p.allRegionsDirty, p.regionsSavedTo, body, dir);
}

// Reads every region file of a body into a new delta snapshot.
std::shared_ptr<const TerrainDeltas> LoadRegionsFor(int body, const std::string& dir, size_t strokes) {
    TerrainDeltasBuilder b(TerrainDeltas{});
    if (DirectoryExists(dir.c_str())) {
        const std::string prefix = TextFormat("b%d_", body);
        FilePathList files = LoadDirectoryFiles(dir.c_str());
        for (unsigned int i = 0; i < files.count; ++i)
            if (std::string(GetFileName(files.paths[i])).rfind(prefix, 0) == 0 && !ReadRegion(files.paths[i], b))
                logger.logf("[World] %s: damaged region file\n", files.paths[i]);
        UnloadDirectoryFiles(files);
    }
    b.state().strokes = strokes;
    return b.finish();
}

bool SaveWorld(Scene& world) {
    if (!DirectoryExists("saves")) MakeDirectory("saves");
    FILE* f = fopen(SAVE_PATH, "w");
    if (!f) return false;
    const DVec3 camWorld = world.renderOrigin + world.camera.position;
    const DVec3 tgtWorld = world.renderOrigin + world.camera.target;
    // v4: the system definition is stored in the save, so a save reloads the system it was made in
    // even if bodies.txt has changed since.
    fprintf(f, "RAYLIBGAME_SAVE 11\n");
    fprintf(f, "seed %u\nbodies %zu\n", worldGen.lastSeed, worldGen.system.bodies.size());
    for (const BodyDef& b : worldGen.system.bodies) fprintf(f, "%s\n", FormatBodyDef(b).c_str());
    fprintf(f, "time %.6f\n", worldTime); // v5: the world clock (orbits)
    fprintf(f, "camera %.6f %.6f %.6f %.6f %.6f %.6f\n", camWorld.x, camWorld.y, camWorld.z, tgtWorld.x, tgtWorld.y, tgtWorld.z);
    // Game-specific state (player, inventory, ...) written by whoever registered the hook.
    fprintf(f, "extras\n");
    if (saveExtrasHook) saveExtrasHook(f);
    fprintf(f, "end_extras\n");
    // v7: terrain edits live in region files next to the save (RegionDir), written incrementally:
    // only regions changed since the last save, unless the body has never been saved there.
    const std::string dir = RegionDir();
    if (!DirectoryExists(dir.c_str())) MakeDirectory(dir.c_str());
    int index = 0, regionsWritten = 0;
    size_t totalStrokes = 0;
    for (auto& obj : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p || p->ringRock >= 0 || p->vehicle) continue; // promoted ring rocks: saved with their ring (below); rockets: by the game
        const int bodyIndex = index++;
        fprintf(f, "planetoid %d rot %f %f %f strokes %zu\n", bodyIndex, p->rotation.x, p->rotation.y, p->rotation.z, p->edits->strokes);
        totalStrokes += p->edits->strokes;
        // v6: chopped / harvested plants (id, world time removed); they regrow from that time.
        const PlantRemovals none;
        const PlantRemovals& rm = p->trees ? p->trees->removals() : none;
        fprintf(f, "plants %zu\n", rm.size());
        for (const auto& kv : rm) fprintf(f, "%llu %.3f\n", (unsigned long long)kv.first, kv.second);
        // v8: build plots and their parts.
        if (p->plots) p->plots->save(f); else fprintf(f, "plots 0\n");
        Factory::Save(f, p->id); // v9: machine contents
        Fluid::Save(f, p->id);   // v11: flowing water
        regionsWritten += SaveBodyRegions(*p, bodyIndex, dir);
    }
    // v10: ring rocks a visit changed (edits, plots), whether they are voxel bodies now or not.
    for (auto& ring : gRings) ring->save(f, dir);
    fclose(f);
    logger.logf("[World] saved %zu strokes to %s (%d region files written)\n", totalStrokes, SAVE_PATH, regionsWritten);
    return true;
}

bool LoadWorldSave(Scene& world) {
    FILE* f = fopen(SAVE_PATH, "r");
    if (!f) return false;
    int version = 0;
    if (fscanf(f, "RAYLIBGAME_SAVE %d", &version) != 1 || version < 2 || version > 11) {
        fclose(f);
        logger.logf("[World] %s: unsupported save version %d\n", SAVE_PATH, version);
        return false;
    }
    Fluid::Reset(); // the save carries its own water
    unsigned int seed = 0;
    SystemDef sys;
    double c[6];
    bool ok = false;
    if (version >= 4) {
        // seed, then "bodies N" and N definition lines (same syntax as bodies.txt).
        size_t n = 0;
        ok = fscanf(f, " seed %u bodies %zu", &seed, &n) == 2;
        int ch;
        while (ok && (ch = fgetc(f)) != '\n' && ch != EOF) {} // rest of the "bodies" line
        char line[1024];
        for (size_t i = 0; ok && i < n; ++i) {
            BodyDef d;
            std::string error;
            ok = fgets(line, sizeof(line), f) != nullptr;
            if (ok) {
                std::string s(line);
                while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
                ok = ParseBodyDef(s, d, error);
                if (!ok) logger.logf("[World] %s: body %zu: %s\n", SAVE_PATH, i, error.c_str());
                else sys.bodies.push_back(d);
            }
        }
    } else {
        // v2/v3 stored the old generator's settings; rebuild the system they described.
        int count = 1, size = 50;
        float systemRadius = 0, spin = 0;
        Vector3 first{};
        ok = fscanf(f, " seed %u count %d size %d system %f %f first %f %f %f", &seed, &count, &size, &systemRadius, &spin, &first.x, &first.y, &first.z) == 8;
        if (ok) sys = LegacySystemDef(count, size, systemRadius, spin, first);
    }
    double savedTime = 0.0;
    if (ok && version >= 5) ok = fscanf(f, " time %lf", &savedTime) == 1;
    ok = ok && fscanf(f, " camera %lf %lf %lf %lf %lf %lf", &c[0], &c[1], &c[2], &c[3], &c[4], &c[5]) == 6;
    if (!ok) { fclose(f); logger.logf("[World] %s: could not read the header\n", SAVE_PATH); return false; }

    worldTime = savedTime; // before regenerating: orbiting bodies are placed for this time
    worldGen.seed = seed;
    worldGen.system = sys;
    worldGen.systemLoaded = true;
    regenerateWorld(world);
    world.setCameraWorld(DVec3(c[0], c[1], c[2]), DVec3(c[3], c[4], c[5]));

    if (version >= 3) {
        char tok[64] = {};
        if (fscanf(f, " %63s", tok) == 1 && std::string(tok) == "extras") {
            if (loadExtrasHook) loadExtrasHook(f); // must consume through "end_extras"
            else while (fscanf(f, " %63s", tok) == 1 && std::string(tok) != "end_extras") {}
        }
    }

    // Edits are attached before the next worldUpdate() requests nodes, so each node is built
    // once with its edits already applied.
    std::vector<Planetoid*> planetoids;
    for (auto& obj : world.objects) if (auto* p = dynamic_cast<Planetoid*>(obj.get())) planetoids.push_back(p);
    int index = 0;
    size_t editCount = 0;
    Vector3 rot{};
    char kind[16] = {};
    const std::string dir = RegionDir();
    while (fscanf(f, " planetoid %d rot %f %f %f %15s %zu", &index, &rot.x, &rot.y, &rot.z, kind, &editCount) == 6) {
        Planetoid* body = index >= 0 && index < (int)planetoids.size() ? planetoids[index] : nullptr;
        std::shared_ptr<const TerrainDeltas> deltas;
        if (std::string(kind) == "strokes") {
            // v7: the deltas are in the region files.
            deltas = LoadRegionsFor(index, dir, editCount);
        } else {
            // v2-v6: a brush list; fold it into deltas (written out as regions on the next save).
            TerrainDeltasBuilder b(TerrainDeltas{});
            for (size_t i = 0; i < editCount; ++i) {
                TerrainEdit e{};
                if (fscanf(f, " %f %f %f %f %f", &e.center.x, &e.center.y, &e.center.z, &e.radius, &e.strength) != 5) break;
                if (body) b.apply(e, *body->field);
            }
            deltas = b.finish();
        }
        PlantRemovals plants;
        size_t plantCount = 0;
        if (version >= 6 && fscanf(f, " plants %zu", &plantCount) == 1) {
            for (size_t i = 0; i < plantCount; ++i) {
                unsigned long long id = 0;
                double t = 0.0;
                if (fscanf(f, " %llu %lf", &id, &t) != 2) break;
                plants[(uint64_t)id] = t;
            }
        }
        // v8: build plots.
        auto plots = std::make_unique<PlotSet>();
        if (version >= 8 && !plots->load(f)) logger.logf("[World] %s: damaged plot data\n", SAVE_PATH);
        if (body) body->plots = plots->plots.empty() ? nullptr : std::move(plots);
        if (version >= 9) Factory::Load(f, body ? body->id : 0); // machine contents, applied when the layouts reach the simulation
        if (version >= 11) Fluid::Load(f, body ? body->id : 0); // flowing water (classified again as it loads)
        if (body) {
            body->edits = deltas;
            body->rotation = rot;
            body->dirtyRegions.clear();
            body->allRegionsDirty = std::string(kind) != "strokes"; // an old brush-list save: write regions next time
            body->regionsSavedTo = body->allRegionsDirty ? std::string() : dir;
            if (body->trees) body->trees->setRemovals(plants);
        }
    }
    // v10: ring rocks with edits or plots (kept as the ring's dormant state until promoted).
    int ringNo = 0, rock = 0;
    size_t rockStrokes = 0;
    while (version >= 10 && fscanf(f, " ringrock %d %d strokes %zu", &ringNo, &rock, &rockStrokes) == 3) {
        if (ringNo >= 0 && ringNo < (int)gRings.size()) gRings[ringNo]->loadRock(f, rock, rockStrokes, dir);
        else { PlotSet skip; skip.load(f); Factory::ReadSaved(f); }
    }
    fclose(f);
    logger.logf("[World] loaded %s (seed %u)\n", SAVE_PATH, seed);
    return true;
}

void loadWorld(Scene& world) {
    // Load the world scene from a file or initialize it
    logger.logf("Loading world scene...\n");
    std::ifstream file("assets/scenes/world");
    if (!file.is_open()) {
        std::cerr << "Failed to open world scene file." << std::endl;
        return;
    }

    std::string line;
    while (getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines

        logger.logf("Processing line: %s\n", line.c_str());
        std::stringstream ss(line);
        std::string cell;
        std::vector<std::string> result;
        while(std::getline(ss, cell, ',')) {
            result.push_back(cell);
        }
        if (result.size() == 13) {
            std::string type = result[0];
            std::string name = result[1];
            Vector3 position = { std::stof(result[2]), std::stof(result[3]), std::stof(result[4]) };
            Vector3 rotation = { std::stof(result[5]), std::stof(result[6]), std::stof(result[7]) };
            Color color = { static_cast<unsigned char>(std::stoi(result[8])),
                            static_cast<unsigned char>(std::stoi(result[9])),
                            static_cast<unsigned char>(std::stoi(result[10])),
                            static_cast<unsigned char>(std::stoi(result[11])) };
            float scale = std::stof(result[12]);

            if (type == "gameobject") {
                logger.logf("Creating GameObject: %s\n", name.c_str());
                std::string modelPath = "assets/models/" + name + "/" + name + ".obj";
                std::string texturePath = "assets/models/" + name + "/diffuse.png";
                Model model;
                if (!FileExists(modelPath.c_str())) {
                    std::cerr << "Model file does not exist: " << modelPath << std::endl;
                    continue;
                } else {
                    model = LoadModel(modelPath.c_str());
                    logger.logf("meshCount: %d, materialCount: %d, boneCount: %d\n", model.meshCount, model.materialCount, model.boneCount);
                    logger.logf("bones: %p, bindPose: %p\n", model.bones, model.bindPose);
                    if (model.boneCount > 0 && (model.bones == nullptr || model.bindPose == nullptr)) {
                        model.boneCount = 0;
                        model.bones = nullptr;
                        model.bindPose = nullptr;
                    }
                }

                if (model.meshCount == 0 || model.materialCount == 0) {
                    std::cerr << "Failed to load model or materials: " << modelPath << std::endl;
                    continue;
                }

                if (FileExists(texturePath.c_str()) && model.materials != nullptr && model.materials[0].maps != nullptr) {
                    model.materials[0].maps[MATERIAL_MAP_DIFFUSE].texture = LoadTexture(texturePath.c_str());
                } else {
                    logger.logf("Texture file does not exist or materials/maps not allocated: %s\n", texturePath.c_str());
                }
                world.objects.push_back(std::make_unique<GameObject>(type, name, position, rotation, color, scale, model));
                logger.logf("GameObject created: %s\n", name.c_str());
            } else if (type == "text") {
                world.uiObjects.push_back(std::make_unique<TextObject>(name, position, rotation, color, scale, static_cast<int>(scale * 20))); // Assuming scale is used for font size
            } else if (type == "cube" || type == "sphere") {
                world.objects.push_back(std::make_unique<GameObject>(type, name, position, rotation, color, scale));
            } else {
                std::cerr << "Unknown object type: " << type << std::endl;
                continue; // Skip unknown types
            }
        } else {
            std::cerr << "Invalid line format: " << line << std::endl;
            continue; // Skip invalid lines
        }
    }
    file.close();
}


bool IsUnderwater(Scene& world, Vector3 renderPoint) {
    Planetoid* p = NearestPlanetoid(world, renderPoint);
    if (!p || !p->field) return false;
    const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
    if (!pf) return false;
    const Vector3 local = p->toLocal(renderPoint);
    const float d = Vector3Length(local);
    const float sea = pf->seaLevelRadius();
    if (sea <= 0.0f || d < 1.0f) return false;
    if (Fluid::Wet(*p, local)) return true; // flowing water
    const Vector3 dir = Vector3Scale(local, 1.0f / d);
    // The water surface here: the sea, or a river or pond above it.
    float surface = sea;
    if (d >= sea && pf->inlandWaterNear(dir, 0.0f)) {
        const float w = pf->waterElevation(dir);
        if (w != PlanetField::kNoWater) surface = std::max(surface, sea + w);
    }
    if (d >= surface) return false;
    return d > pf->radius + pf->elevation(dir); // above the bed
}
