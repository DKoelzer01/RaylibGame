#include "world.h"
#include "terrainarena.h"
#include "ring.h"
#include "player.h"
#include "jobsystem.h"
#include "utils/logger.h"
#include "worldtypes.h"
#include "factory.h"
#include "fluid.h"
#include "terraingpu.h"
#include "gputimer.h"
#include <chrono>

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <climits>
#include <cstdint>
#include <cstdlib>
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
    std::atomic<long long> gRequests{ 0 }, gBuilt{ 0 }, gSkipped{ 0 }, gDropped{ 0 }; // streaming counters (#258)
    // RAYLIBGAME_STREAM_OLD=1: the streaming before #258/#261 (no request cap, no retention, no telescope keep, no free
    // cap, FIFO uploads), an A/B switch for --bench-fly.
    const bool gStreamOld = [] { const char* e = getenv("RAYLIBGAME_STREAM_OLD"); return e && e[0] == '1'; }();
    // Preload-ring builds (#296) are queued at this priority plus their distance rank: after every view build.
    constexpr float kBackgroundPriority = 1e7f;
    long long gEvicted = 0;                    // cache evictions (#295)
    bool gPreloadUrgent = false;               // spawn preload (#294): ring builds without the spare-time gate
    std::atomic<int> gForegroundBuilds{ 0 };   // view/edit builds whose closures are alive (#296)
    struct ForegroundToken {
        ForegroundToken() { ++gForegroundBuilds; }
        ~ForegroundToken() { --gForegroundBuilds; }
    };
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
        // The same placement for the terrain arena's GPU culling pass (chunk_cull.comp).
        ArenaBody arena(float texScale) const {
            ArenaBody b;
            for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) b.rot[i * 3 + j] = r[i][j];
            b.scale = s;
            b.offset[0] = offset.x; b.offset[1] = offset.y; b.offset[2] = offset.z;
            b.rs = rs;
            b.texPeriod = 1.0 / std::max(texScale, 1e-6f);
            return b;
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
    // Arena chunks are only listed: placement and culling run on the GPU (chunk_cull.comp), one multi-draw.
    for (Chunk* c : drawList) {
        if (c->arena.valid()) { if (c->isActive) TerrainArena_Queue(c->arena, c->clipLocal()); continue; }
        c->drawChunk(c2, placer.model(c->key));
    }
    ArenaCull cull;
    cull.frustum = ctx.frustum;
    const int n = TerrainArena_Flush(ctx.terrainMaterial->shader.id, false, placer.arena(ctx.texScale), cull);
    if (ctx.drawCalls) *ctx.drawCalls += n; // listed (the GPU culls some of them)
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
            // Cheap and drops most of the list: done here, so the GPU pass gets no zero-instance commands for them.
            const float h = c->key.size() * 0.5f;
            const Vector3 v = Vector3Subtract(Vector3AddValue(c->key.origin(), h), pL);
            const float along = Vector3DotProduct(v, dL), reach = rL + h * 1.75f;
            if (Vector3LengthSqr(v) - along * along > reach * reach) continue;
        }
        if (c->arena.valid()) { if (c->isActive) TerrainArena_Queue(c->arena, c->clipLocal()); continue; } // frustum-culled on the GPU (depth ignores the clip)
        c->drawChunkDepth(c2, placer.model(c->key));
    }
    ArenaCull cull;
    cull.frustum = ctx.frustum;
    TerrainArena_Flush(ctx.depthMaterial->shader.id, true, placer.arena(ctx.texScale), cull);
}

// ---------------------------------------------------------------------------------------------
// Node requests / results
// ---------------------------------------------------------------------------------------------

static void RequestNode(Planetoid& p, const NodeKey& key, float priority) {
    ++gRequests;
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
    Chunk& node = *p.nodes.at(key);
    node.background = priority >= kBackgroundPriority;
    std::shared_ptr<std::atomic<uint32_t>> wanted = node.wanted;
    wanted->store(version);
    // View and edit builds hold a token until their closures are gone (built, skipped, or dropped inside the GPU
    // pump): WorldIdle waits for those only, not for the preload ring's background builds (#296).
    std::shared_ptr<ForegroundToken> token = node.background ? nullptr : std::make_shared<ForegroundToken>();
    // Plan 6 step 4: planets sample the density block on the GPU (terraingpu.h). The worker prepares the
    // height grid and submits; the main thread's pump dispatches batches and queues the second half
    // (edits, marching cubes) back on the workers. Identical results to the CPU path.
    const bool gpu = TerrainGpu_Available();
    Jobs().submit(priority, [field, edits, id, key, version, wanted, gpu, priority, token] {
        if (wanted->load() != version) { ++gSkipped; return; } // superseded or freed while queued
        PROFILE_SCOPE("Build chunk");
        ChunkMeshData data;
        data.planetoidId = id;
        data.key = key;
        data.version = version;
        data.priority = priority;
        const auto* pf = gpu ? dynamic_cast<const PlanetField*>(field.get()) : nullptr;
        if (pf) {
            const auto t0 = std::chrono::steady_clock::now();
            GpuBlockJob job;
            job.block = ChunkBlockOf(key);
            if (pf->prepareBlock(job.block.origin, job.block.stride, job.block.E, job.block.footprint, job.grid) &&
                pf->blockIslands(job.block.origin, job.block.stride, job.block.E, job.islands, 64)) {
                const float prepMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
                job.keep = field;
                job.field = pf;
                job.wanted = [wanted, version] { return wanted->load() == version; };
                job.priority = priority; // edits (negative) go ahead of background LOD builds
                auto shared = std::make_shared<ChunkMeshData>(std::move(data));
                job.done = [field, edits, shared, wanted, version, priority, prepMs, token](std::shared_ptr<const float> batch, size_t offset) {
                    Jobs().submit(priority, [field, edits, shared, wanted, version, prepMs, batch, offset, token] {
                        if (wanted->load() != version) { ++gSkipped; return; }
                        BuildChunkMeshFromBase(*field, *edits, *shared, batch.get() + offset, prepMs);
                        std::lock_guard<std::mutex> lock(resultsMutex);
                        results.push_back(std::move(*shared));
                        ++gBuilt;
                    });
                };
                TerrainGpu_Submit(std::move(job));
                return;
            }
        }
        BuildChunkMesh(*field, *edits, data);
        std::lock_guard<std::mutex> lock(resultsMutex);
        results.push_back(std::move(data));
        ++gBuilt;
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

// The root nodes a selection descends from: the body's 2x2x2 roots, or (open world) a grid of root tiles around
// the camera spanning the terrain height band.
template <class F>
static void ForEachRoot(const Planetoid& p, Vector3 camLocal, F&& f) {
    if (p.flat) {
        const float rs = kVoxelSize * (float)(CHUNK_SIZE << p.rootLevel);
        const int cx = (int)floorf(camLocal.x / rs), cz = (int)floorf(camLocal.z / rs);
        const int y0 = (int)floorf(p.minHeight / rs), y1 = (int)floorf(p.maxHeight / rs);
        for (int z = cz - 2; z <= cz + 2; ++z)
            for (int x = cx - 2; x <= cx + 2; ++x)
                for (int y = y0; y <= y1; ++y) f(NodeKey{ p.rootLevel, x, y, z });
    } else {
        for (int i = 0; i < 8; ++i) f(NodeKey{ p.rootLevel, -1 + (i & 1), -1 + ((i >> 1) & 1), -1 + ((i >> 2) & 1) });
    }
}

// Preload ring (#294): every node CollectDesired would pick from SOME camera point of a region around the player -
// all level-0 chunks near it plus the coarser shells around them - so walking anywhere in it needs no new builds. A
// node is picked from q when q is at least lodFactor * size away and its parent is split from q; over the region that
// is (nearest point) dmin < lodFactor * parent size and (farthest) dmax >= lodFactor * size.
// The region follows the ground on planets: points within r of the centre along the surface and from kRingBelow under
// to kRingAbove over the ground (walking and jumping; the margin below covers overhangs and the height map's missing
// fine detail). A ball of radius r instead reached 300 m down into the caves: 415k nodes, 80 M vertices (#294).
// Elsewhere (asteroids, open world) it is the ball. Runs on a worker: it samples the ground once per node visited.
struct RingParams {
    std::shared_ptr<const DensityField> field;
    const PlanetField* planet = nullptr; // the ground-following region; null: the ball
    bool flat = false;
    float minHeight = 0, maxHeight = 0, occluder = 0, lodFactor = 2.0f;
    int rootLevel = 0;
    Vector3 c{}, cHigh{}, up{};  // centre, the region's highest point (horizon test), the local vertical
    float r = 0;
};
static constexpr float kRingBelow = 8.0f, kRingAbove = 8.0f, kRingBand = 12.0f;
static constexpr float kRingMaxSpeed = 15.0f; // m/s: faster than running, the ring pauses

// Nearest and farthest distance from the region to the box (both approximate, erring towards more nodes).
// Returns whether the box reaches the ground band (kRingBand around the height map): only those are preloaded, the
// rest (caves under the walking region, which the view also builds within 32 m but nobody sees) loads on demand.
static bool RingDistances(const RingParams& rp, const NodeKey& key, Vector3 lo, Vector3 hi, float& dmin, float& dmax) {
    if (!rp.planet) {
        const float d = DistanceToBox(rp.c, lo, hi);
        dmin = d - rp.r; dmax = d + rp.r;
        return true;
    }
    const float half = 0.5f * key.size(), hd = 1.7320508f * half;
    const Vector3 b = Vector3AddValue(lo, half);
    const float rb = Vector3Length(b);
    const float ground = rp.planet->radius + rp.planet->elevation(Vector3Scale(b, 1.0f / std::max(rb, 1e-3f)), key.size());
    const Vector3 w = Vector3Subtract(b, rp.c);
    const float wv = Vector3DotProduct(w, rp.up);
    const float wh = Vector3Length(Vector3Subtract(w, Vector3Scale(rp.up, wv)));
    const float dh = std::max(0.0f, wh - rp.r);
    const float dv = std::max({ 0.0f, ground - kRingBelow - rb, rb - (ground + kRingAbove) });
    dmin = sqrtf(dh * dh + dv * dv) - hd;
    const float fv = std::max(fabsf(rb - (ground - kRingBelow)), fabsf(rb - (ground + kRingAbove)));
    dmax = sqrtf((wh + rp.r) * (wh + rp.r) + fv * fv);
    return fabsf(rb - ground) < hd + kRingBand;
}

static void CollectRing(const RingParams& rp, const NodeKey& key, std::vector<NodeKey>& out) {
    const Vector3 lo = key.origin();
    const Vector3 hi = Vector3AddValue(lo, key.size());
    if (rp.flat) {
        if (hi.y < rp.minHeight - 8.0f || lo.y > rp.maxHeight + 8.0f) return;
    } else {
        if (DistanceToBox({ 0, 0, 0 }, lo, hi) > rp.field->maxRadius() + 24.0f) return;
        const float far = Vector3Length({ std::max(fabsf(lo.x), fabsf(hi.x)), std::max(fabsf(lo.y), fabsf(hi.y)), std::max(fabsf(lo.z), fabsf(hi.z)) });
        if (far < rp.field->minRadius()) return;
        if (BehindHorizon(rp.cHigh, lo, hi, rp.occluder)) return;
    }
    float dmin, dmax;
    if (!RingDistances(rp, key, lo, hi, dmin, dmax)) return;
    const float split = rp.lodFactor * key.size();
    if (key.level == 0 || dmax >= split) out.push_back(key);
    if (key.level > 0 && dmin < split)
        for (int i = 0; i < 8; ++i)
            CollectRing(rp, { key.level - 1, key.x * 2 + (i & 1), key.y * 2 + ((i >> 1) & 1), key.z * 2 + ((i >> 2) & 1) }, out);
}

// A set of node keys for the LOD selection: open addressing, buffers kept between calls (#258/#204: the selection
// rebuilt six std::unordered_sets of thousands of keys every time the camera moved, a heap allocation per insert;
// it was ~2 ms a frame flying).
struct FlatKeySet {
    std::vector<NodeKey> keys;
    std::vector<uint8_t> used;
    size_t mask = 0, n = 0;
    static size_t Hash(const NodeKey& k) {
        uint64_t h = (uint64_t)(uint32_t)k.x * 0x9E3779B97F4A7C15ull ^ (uint64_t)(uint32_t)k.y * 0xC2B2AE3D27D4EB4Full ^
                     (uint64_t)(uint32_t)k.z * 0x165667B19E3779F9ull ^ (uint64_t)(uint32_t)k.level * 0x27D4EB2F165667C5ull;
        h ^= h >> 31; h *= 0xBF58476D1CE4E5B9ull; h ^= h >> 29;
        return (size_t)h;
    }
    void reset(size_t expected) {
        size_t cap = 16;
        while (cap < expected * 2 + 16) cap <<= 1;
        if (keys.size() != cap) { keys.resize(cap); used.assign(cap, 0); }
        else std::fill(used.begin(), used.end(), 0);
        mask = cap - 1; n = 0;
    }
    bool insert(const NodeKey& k) { // true if it was not in yet
        if ((n + 1) * 2 > keys.size()) grow();
        for (size_t i = Hash(k) & mask;; i = (i + 1) & mask) {
            if (!used[i]) { used[i] = 1; keys[i] = k; ++n; return true; }
            if (keys[i] == k) return false;
        }
    }
    size_t count(const NodeKey& k) const {
        if (keys.empty()) return 0;
        for (size_t i = Hash(k) & mask;; i = (i + 1) & mask) {
            if (!used[i]) return 0;
            if (keys[i] == k) return 1;
        }
    }
    bool empty() const { return n == 0; }
private:
    void grow() {
        std::vector<NodeKey> old;
        std::vector<uint8_t> oldUsed;
        old.swap(keys); oldUsed.swap(used);
        keys.resize(old.size() * 2); used.assign(keys.size(), 0);
        mask = keys.size() - 1; n = 0;
        for (size_t i = 0; i < old.size(); ++i) if (oldUsed[i]) insert(old[i]);
    }
};

// A body's preload ring (#294/#296): its node keys (nearest/coarsest first), requested in spare streaming time,
// kept resident but not drawn. Main thread only.
struct PreloadRing {
    Vector3 centre = { 1e30f, 1e30f, 1e30f }; // body-local
    std::vector<NodeKey> keys;
    size_t cursor = 0;   // keys before it have been requested (or were already there)
    int ready = 0;       // keys flagged inRing that have a build (counted as they are flagged or built)
    bool complete = true;
    struct Job { Vector3 centre{}; std::vector<NodeKey> keys; std::atomic<bool> done{ false }; };
    std::shared_ptr<Job> job; // the next ring, being collected on a worker
    Vector3 lastCam = { 1e30f, 1e30f, 1e30f }; double lastTime = 0.0; float speed = 0.0f; // camera speed against the body, smoothed
};
static std::unordered_map<uint64_t, PreloadRing> gPreload;
static std::unordered_map<uint64_t, int> gCachedCount; // per body, from its last LOD update (#295)

// Decide which nodes to build and draw. Old nodes stay on screen until every desired node that
// overlaps them is ready, so LOD changes never open holes. Returns how many builds it requested.
static int UpdateLod(Planetoid& p, Vector3 cameraRender) {
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
    if (!p.lodDirty && Vector3DistanceSqr(camLocal, p.lodCamLocal) < 1.0f) return 0;
    p.lodDirty = false;
    p.lodCamLocal = camLocal;
    std::vector<NodeKey> desired;
    ForEachRoot(p, camLocal, [&](const NodeKey& r) { CollectDesired(p, r, view, desired); });
    static FlatKeySet desiredSet;
    desiredSet.reset(desired.size());
    for (const NodeKey& k : desired) desiredSet.insert(k);
    const double now = GetTime();

    // New builds, nearest and coarsest first, at most kMaxNewRequests a frame (#258): flying fast asked for every
    // newly desired node at once and filled the queue with builds that were stale before they started.
    const size_t kMaxNewRequests = gStreamOld ? SIZE_MAX : 96;
    std::vector<std::pair<float, NodeKey>> toRequest;
    for (const NodeKey& k : desired) {
        if (p.requestedVersion.find(k) != p.requestedVersion.end()) {
            // Queued by the preload ring (background priority) and not built yet: the view wants it now (#296).
            auto n = p.nodes.find(k);
            if (n == p.nodes.end() || n->second->ready || !n->second->background) continue;
        }
        const Vector3 lo = k.origin();
        const float dist = LodDistance(view, lo, Vector3AddValue(lo, k.size()));
        toRequest.push_back({ dist / (1.0f + k.level), k });
    }
    if (toRequest.size() > kMaxNewRequests) {
        std::partial_sort(toRequest.begin(), toRequest.begin() + kMaxNewRequests, toRequest.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        toRequest.resize(kMaxNewRequests);
        p.lodDirty = true; // the rest next frame
    }
    for (const auto& r : toRequest) RequestNode(p, r.second, r.first);
    for (const NodeKey& k : desired) if (auto it = p.nodes.find(k); it != p.nodes.end()) it->second->lastWanted = now;
    // Telescope (#261): while zoomed, the normal selection stays resident (not drawn), so closing it has nothing
    // to regenerate - closing froze the game while the coarse nodes in the view cone were built again.
    static FlatKeySet keepSet;
    keepSet.reset(64);
    if (view.zoom > 1.0f && !gStreamOld) {
        std::vector<NodeKey> unzoomed;
        LodView plain{ camLocal };
        ForEachRoot(p, camLocal, [&](const NodeKey& r) { CollectDesired(p, r, plain, unzoomed); });
        for (const NodeKey& k : unzoomed) if (auto it = p.nodes.find(k); it != p.nodes.end()) { keepSet.insert(k); it->second->lastWanted = now; }
    }

    // Previously drawn nodes that are no longer desired, kept while a desired overlap is pending.
    // Overlap tests walk ancestors through hash sets (comparing every drawn node with every desired
    // one cost ~3 ms a frame with thousands of nodes while flying).
    auto parentOf = [](const NodeKey& k) { return NodeKey{ k.level + 1, k.x >> 1, k.y >> 1, k.z >> 1 }; };
    static FlatKeySet pendingDesired, pendingWithin; // pending nodes, and every node containing one
    pendingDesired.reset(desired.size()); pendingWithin.reset(desired.size());
    for (const NodeKey& d : desired) {
        if (auto n = p.nodes.find(d); n != p.nodes.end() && n->second->ready) continue; // not requested yet (the per-frame cap) counts as pending
        pendingDesired.insert(d);
        for (NodeKey k = d; k.level <= p.rootLevel; k = parentOf(k))
            if (!pendingWithin.insert(k)) break; // its ancestors are in already
    }
    std::vector<Chunk*> retained;
    static FlatKeySet retainedSet;
    retainedSet.reset(256);
    for (Chunk* c : p.drawList) {
        if (desiredSet.count(c->key)) continue;
        bool keep = pendingWithin.count(c->key) > 0;                // a pending desired node lies inside it
        for (NodeKey k = parentOf(c->key); !keep && k.level <= p.rootLevel; k = parentOf(k))
            keep = pendingDesired.count(k) > 0;                      // inside a pending desired node
        if (keep) { retained.push_back(c); retainedSet.insert(c->key); }
    }
    std::vector<Chunk*> draw = retained;
    for (const NodeKey& d : desired) {
        auto n = p.nodes.find(d);
        if (n == p.nodes.end()) continue; // requested on a later frame (#258 cap)
        Chunk* c = n->second.get();
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
    static FlatKeySet drawn; // also "on screen" for the frees below
    drawn.reset(p.drawList.size());
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

    // Free nodes that are neither desired, on screen, kept for the telescope nor in the preload ring (late results for
    // them are dropped). A node still building is freed at once, which cancels its queued work. Built ones go to the
    // chunk cache (#295): they stay until the cached nodes hold more than cacheVertices or the body has more than
    // cacheNodes nodes, then the least recently wanted go first - turning around or doubling back reuses them
    // instead of building them again (it was a fixed 4 s hold, #258). At most kMaxFrees a frame: freeing thousands at
    // once (closing the telescope) hitched.
    const int kMaxFrees = gStreamOld ? INT_MAX : 256;
    // The full pass costs ~1-3 ms with the preload ring's ~40k nodes: at most 4 times a second (or straight on when it
    // stopped at kMaxFrees). Undesired builds are cancelled up to 0.25 s later.
    if (!gStreamOld && now - p.lastFreeSweep < 0.25) return (int)toRequest.size();
    p.lastFreeSweep = now;
    std::vector<std::pair<double, NodeKey>> cached;
    int frees = 0;

    auto vertsOf = [](const Chunk& c) -> long long { return c.arena.valid() ? c.arena.vCount : c.mesh.vertexCount; };
    long long cachedVerts = 0; // the cache's own vertices: the preload ring and the view are not counted against it
    for (auto it = p.nodes.begin(); it != p.nodes.end();) {
        const Chunk* c = it->second.get();
        if (c->inRing) { ++it; continue; } // most nodes: no lookups
        const bool keep = desiredSet.count(it->first) || drawn.count(it->first) || keepSet.count(it->first);
        if (keep) { ++it; continue; }
        if (c->ready && !gStreamOld) { cached.push_back({ c->lastWanted, it->first }); cachedVerts += vertsOf(*c); ++it; continue; }
        if (frees >= kMaxFrees) { p.lodDirty = true; p.lastFreeSweep = -1e9; ++it; continue; } // the rest next frame
        p.requestedVersion.erase(it->first);
        it = p.nodes.erase(it);
        ++frees;
    }

    size_t cachedCount = cached.size();
    auto over = [&] { return cachedVerts > lodSettings.cacheVertices || cachedCount > (size_t)lodSettings.cacheNodes; };
    if (!cached.empty() && over()) {
        std::sort(cached.begin(), cached.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        size_t i = 0;
        for (; i < cached.size() && over(); ++i) {
            if (frees >= kMaxFrees) { p.lodDirty = true; p.lastFreeSweep = -1e9; break; }
            auto it = p.nodes.find(cached[i].second);
            cachedVerts -= vertsOf(*it->second);
            p.requestedVersion.erase(it->first);
            p.nodes.erase(it);
            --cachedCount;
            ++frees;
            ++gEvicted;
        }
        cached.erase(cached.begin(), cached.begin() + (std::ptrdiff_t)i);
    }
    gCachedCount[p.id] = (int)cached.size();
    return (int)toRequest.size();
}

// Preload ring step (#294/#296) for one body: re-centre on the camera when it moved a quarter of the radius, then
// request up to `spare` ring nodes nobody asked for yet (background priority, so the view's builds go first).
static void UpdatePreload(Planetoid& p, Vector3 cameraRender, int& spare) {
    PreloadRing& ring = gPreload[p.id];
    const float r = lodSettings.preloadRadius / std::max(p.scale, 1e-6f);
    const Vector3 c = p.toLocal(cameraRender);
    bool near = r > 0.0f && !gStreamOld && !p.vehicle && p.field;
    if (near && !p.flat) {
        const float dc = Vector3Length(c); // the ball must reach the terrain band
        near = dc - r < p.field->maxRadius() + 24.0f && dc + r > p.field->minRadius();
    }
    if (!near) {
        if (!ring.keys.empty()) { // its nodes go to the cache
            for (const NodeKey& k : ring.keys) if (auto it = p.nodes.find(k); it != p.nodes.end()) it->second->inRing = false;
            ring = PreloadRing{}; p.lodDirty = true;
        }
        return;
    }
    const double now = GetTime();
    // Paused while the camera moves faster than running pace against the body (flying, falling): the ring there would
    // be stale before it was built, and its builds cost ~1.7 ms of GPU time a frame flying at 55 m/s (#296). The spawn
    // preload ignores this.
    {
        const double dt = now - ring.lastTime;
        if (ring.lastCam.x < 1e29f && dt > 0.0 && dt < 1.0) {
            const float v = Vector3Distance(c, ring.lastCam) / (float)dt, a = (float)std::min(1.0, dt * 4.0);
            ring.speed += (v - ring.speed) * a;
        }
        ring.lastCam = c; ring.lastTime = now;
    }
    const bool paused = !gPreloadUrgent && ring.speed > kRingMaxSpeed;
    // Moved a quarter of the radius from the ring (or the one being collected): collect the next one on a worker.
    const Vector3 aimed = ring.job ? ring.job->centre : ring.centre;
    if (!ring.job && !paused && Vector3DistanceSqr(c, aimed) > r * r / 16.0f) {
        RingParams rp;
        rp.field = p.field;
        rp.planet = p.flat ? nullptr : dynamic_cast<const PlanetField*>(p.field.get());
        rp.flat = p.flat; rp.minHeight = p.minHeight; rp.maxHeight = p.maxHeight;
        rp.occluder = p.field->occluderRadius() * p.scale;
        rp.lodFactor = lodSettings.lodFactor;
        rp.rootLevel = p.rootLevel;
        rp.c = c; rp.r = r;
        const float lc = Vector3Length(c);
        rp.up = p.flat ? Vector3{ 0, 1, 0 } : lc > 1e-3f ? Vector3Scale(c, 1.0f / lc) : Vector3{ 0, 1, 0 };
        rp.cHigh = p.flat ? c : Vector3Add(c, Vector3Scale(rp.up, rp.planet ? kRingAbove + 64.0f : r));
        auto job = std::make_shared<PreloadRing::Job>();
        job->centre = c;
        ring.job = job;
        ring.complete = false;
        std::vector<NodeKey> roots;
        ForEachRoot(p, c, [&](const NodeKey& k) { roots.push_back(k); });
        Jobs().submit(0.0f, [rp, job, roots] {
            PROFILE_SCOPE("Preload ring");
            for (const NodeKey& k : roots) CollectRing(rp, k, job->keys);
            std::vector<std::pair<float, NodeKey>> order;
            order.reserve(job->keys.size());
            for (const NodeKey& k : job->keys) {
                const Vector3 lo = k.origin();
                order.push_back({ DistanceToBox(rp.c, lo, Vector3AddValue(lo, k.size())) / (1.0f + k.level), k });
            }
            std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
            for (size_t i = 0; i < order.size(); ++i) job->keys[i] = order[i].second;
            job->done.store(true);
        });
    }
    if (ring.job && ring.job->done.load()) {
        for (const NodeKey& k : ring.keys) // the old ring becomes the freshest part of the cache
            if (auto it = p.nodes.find(k); it != p.nodes.end()) { it->second->lastWanted = now; it->second->inRing = false; }
        ring.centre = ring.job->centre;
        ring.keys = std::move(ring.job->keys);
        ring.job.reset();
        // Flag the new ring's existing nodes (kept from now on) and count the built ones; nodes flagged later are
        // counted when flagged (the cursor below) or when their first build arrives (IntegrateResult).
        ring.ready = 0;
        for (const NodeKey& k : ring.keys)
            if (auto it = p.nodes.find(k); it != p.nodes.end()) { it->second->inRing = true; ring.ready += it->second->ready ? 1 : 0; }
        ring.cursor = 0; ring.complete = false;
        p.lodDirty = true; // the free loop now protects the new ring
    }
    while (!paused && spare > 0 && ring.cursor < ring.keys.size()) {
        const NodeKey& k = ring.keys[ring.cursor++];
        if (p.requestedVersion.find(k) != p.requestedVersion.end()) { // built, building or cached already
            if (auto it = p.nodes.find(k); it != p.nodes.end() && !it->second->inRing) {
                it->second->inRing = true; // the view made it after the ring was installed
                ring.ready += it->second->ready ? 1 : 0;
            }
            continue;
        }
        const Vector3 lo = k.origin();
        RequestNode(p, k, kBackgroundPriority + DistanceToBox(c, lo, Vector3AddValue(lo, k.size())) / (1.0f + k.level));
        p.nodes.at(k)->inRing = true;
        --spare;
    }
    if (!ring.job) ring.complete = ring.cursor >= ring.keys.size() && ring.ready >= (int)ring.keys.size();
}

static void IntegrateResult(const ChunkMeshData& r) {
    Planetoid* p = Planetoid::find(r.planetoidId);
    if (!p) return;                                                        // planetoid was destroyed
    auto v = p->requestedVersion.find(r.key);
    if (v == p->requestedVersion.end() || v->second != r.version) { ++gDropped; return; } // superseded or freed
    auto it = p->nodes.find(r.key);
    if (it == p->nodes.end()) return;

    buildMsTotal += r.buildMs;
    sampleMsTotal += r.sampleMs;
    buildCount++;
    if (r.overflow) logger.logf("[ERROR] Node L%d (%d,%d,%d) exceeds 65535 vertices; skipped\n", r.key.level, r.key.x, r.key.y, r.key.z);
    const bool wasReady = it->second->ready;
    it->second->setMesh(r);
    if (!wasReady && it->second->inRing)
        if (auto ring = gPreload.find(p->id); ring != gPreload.end()) ++ring->second.ready; // (#294 progress)
    // Readiness changed: the draw list may change. Not for a preload-ring build (#296): the view did not ask for
    // it (it re-requests a ring node it wants), and re-selecting after every background upload cost frame time.
    if (!it->second->background) p->lodDirty = true;
}

// Idle for the view: no view or edit build alive, and none of their results waiting for upload. The preload ring's
// background builds do not count (#296): loading screens and scripted tests would wait seconds for the ring.
static bool PipelineIdle() {
    if (gForegroundBuilds.load() != 0) return false; // a result is pushed before its build's token is released
    std::lock_guard<std::mutex> lock(resultsMutex);
    for (const ChunkMeshData& r : results) if (r.priority < kBackgroundPriority) return false;
    for (const ChunkMeshData& r : pendingUploads) if (r.priority < kBackgroundPriority) return false;
    return true;
}

static bool AnyLodPending(Scene& world) {
    for (auto& obj : world.objects)
        if (auto* p = dynamic_cast<Planetoid*>(obj.get()))
            for (auto& n : p->nodes) if (!n.second->ready && !n.second->background) return true;
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
    GpuTimer_CpuMark("w-before");
    TerrainGpu_SetBackground(kBackgroundPriority, !gPreloadUrgent); // preload-ring blocks: gently, except during the spawn preload
    TerrainGpu_Pump(); // GPU density blocks (Plan 6 step 4): finished batches back to the workers, the next batch out
    GpuTimer_CpuMark("w-gpu-pump");
    // Ring rocks near the camera become voxel bodies (and back) before this frame's LOD selection.
    for (auto& ring : Rings()) ring->update(world, world.renderOrigin + world.camera.position);

    {
        std::lock_guard<std::mutex> lock(resultsMutex);
        for (ChunkMeshData& r : results) pendingUploads.push_back(std::move(r));
        results.clear();
    }

    // Uploads nearest first (#258): edits, then by the distance from the camera now (the queue was FIFO, so near
    // ground waited behind far nodes requested earlier). Nodes nobody wants any more drop out in IntegrateResult.
    if (pendingUploads.size() > 1 && !gStreamOld) {
        std::unordered_map<uint64_t, Vector3> camOf;
        for (auto& obj : world.objects)
            if (auto* p = dynamic_cast<Planetoid*>(obj.get())) camOf[p->id] = p->toLocal(world.camera.position);
        auto rank = [&](const ChunkMeshData& d) {
            if (d.priority < 0.0f) return d.priority - 1e9f;
            if (d.priority >= kBackgroundPriority) return d.priority + 1e9f; // preload ring (#296): after the view's
            auto it = camOf.find(d.planetoidId);
            if (it == camOf.end()) return 1e9f;
            const Vector3 lo = d.key.origin();
            return DistanceToBox(it->second, lo, Vector3AddValue(lo, d.key.size())) / (1.0f + d.key.level);
        };
        std::vector<std::pair<float, size_t>> order;
        order.reserve(pendingUploads.size());
        for (size_t i = 0; i < pendingUploads.size(); ++i) order.push_back({ rank(pendingUploads[i]), i });
        std::stable_sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        std::deque<ChunkMeshData> sorted;
        for (const auto& o : order) sorted.push_back(std::move(pendingUploads[o.second]));
        pendingUploads.swap(sorted);
    }
    // GPU uploads are the only main-thread cost; cap them per frame so streaming never hitches.
    const auto t0 = std::chrono::steady_clock::now();
    int uploads = 0;
    // Preload-ring results (sorted last) get at most 1 ms of it outside the spawn preload: spare time only (#296).
    const float backgroundBudgetMs = gPreloadUrgent ? uploadBudgetMs : std::min(uploadBudgetMs, 1.0f);
    while (!pendingUploads.empty()) {
        const float elapsed = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
        if (uploads > 0 && pendingUploads.front().priority >= kBackgroundPriority && elapsed > backgroundBudgetMs) break;
        IntegrateResult(pendingUploads.front());
        pendingUploads.pop_front();
        uploads++;
        if (std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count() > uploadBudgetMs) break;
    }

    GpuTimer_CpuMark("w-uploads");
    // LOD after uploads, so this frame's draw list already includes nodes that just finished.
    int viewRequests = 0;
    for (auto& obj : world.objects)
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) viewRequests += UpdateLod(*p, world.camera.position);
    GpuTimer_CpuMark("w-lod");
    // Preload rings (#294/#296) get the request budget the view left unused this frame, and only while the build
    // queue is short, so they never hold up what is on screen. During the spawn preload they run at full rate.
    {
        const int kViewCap = 96, kUrgentCap = 256;
        const int backlog = (int)Jobs().queued() + TerrainGpu_Pending();
        int spare = gPreloadUrgent ? kUrgentCap - viewRequests
                                   : (backlog < 2 * (int)Jobs().threadCount() ? kViewCap - viewRequests : 0);
        spare = std::max(spare, 0);
        std::unordered_set<uint64_t> alive;
        for (auto& obj : world.objects)
            if (auto* p = dynamic_cast<Planetoid*>(obj.get())) { UpdatePreload(*p, world.camera.position, spare); alive.insert(p->id); }
        for (auto it = gPreload.begin(); it != gPreload.end();) it = alive.count(it->first) ? std::next(it) : gPreload.erase(it);
        for (auto it = gCachedCount.begin(); it != gCachedCount.end();) it = alive.count(it->first) ? std::next(it) : gCachedCount.erase(it);
        stats.preloadNodes = stats.preloadReady = stats.cachedNodes = 0;
        bool complete = true;
        for (auto& [id, ring] : gPreload) { stats.preloadNodes += (int)ring.keys.size(); stats.preloadReady += ring.ready; complete &= ring.complete; }
        for (auto& [id, n] : gCachedCount) stats.cachedNodes += n;
        stats.evicted = gEvicted;
        if (gPreloadUrgent && complete) {
            gPreloadUrgent = false;
            // Still behind the loading screen: grow the terrain arena once for what play adds - a re-centred ring's
            // new part (about a third of it) and a full cache - so it never regrows (a GB-sized copy) mid-walk.
            const ArenaStats as = TerrainArena_Stats();
            if (as.vertices > 0) {
                const double ipv = (double)as.indices / as.vertices;
                const long long v = (long long)as.vertices + as.vertices / 3 + lodSettings.cacheVertices;
                TerrainArena_Reserve((int)std::min<long long>(v, INT_MAX / 2), (int)std::min<long long>((long long)(v * ipv), INT_MAX / 2));
            }
        }
    }
    GpuTimer_CpuMark("w-preload");
    for (auto& obj : world.objects)
        if (auto* p = dynamic_cast<Planetoid*>(obj.get()); p && p->trees) p->trees->update(*p, p->toLocal(world.camera.position));

    GpuTimer_CpuMark("w-trees");
    stats.queuedJobs = (int)Jobs().queued();
    stats.runningJobs = Jobs().running();
    stats.pendingUploads = (int)pendingUploads.size();
    stats.uploadsLastFrame = uploads;
    stats.workerThreads = (int)Jobs().threadCount();
    stats.requests = gRequests; stats.built = gBuilt; stats.skippedStale = gSkipped; stats.droppedResults = gDropped;
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

void BeginSpawnPreload() {
    gPreloadUrgent = lodSettings.preloadRadius > 0.0f && !gStreamOld;
    for (auto& [id, ring] : gPreload) { ring.centre = { 1e30f, 1e30f, 1e30f }; ring.job.reset(); } // re-centre on the camera as it is now
}

bool SpawnPreloadPending(float* progress) {
    if (progress) *progress = stats.preloadNodes ? (float)stats.preloadReady / (float)stats.preloadNodes : 1.0f;
    return gPreloadUrgent;
}

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

void unloadWorld(Scene& world) {
    Jobs().clearQueued();  // in-flight jobs finish; their results are dropped (planetoid ids are gone)
    Factory::Reset();      // machines belong to the old bodies
    Fluid::Reset();        // so does flowing water
    player.nearest = nullptr; // every body goes (the player refinds its ground next update)
    gRings.clear();        // before their parent bodies go
    world.objects.clear(); // releases chunk meshes
}

void regenerateWorld(Scene& world) {
    unloadWorld(world);
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
