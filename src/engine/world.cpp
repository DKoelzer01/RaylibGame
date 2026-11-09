#include "world.h"
#include "jobsystem.h"
#include "utils/logger.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
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
    double buildMsTotal = 0.0;
    long buildCount = 0;

    JobSystem& Jobs() {
        if (!jobs) jobs = std::make_unique<JobSystem>();
        return *jobs;
    }
}

// ---------------------------------------------------------------------------------------------
// Planetoid
// ---------------------------------------------------------------------------------------------

Planetoid::Planetoid(std::string name, DVec3 wp, Vector3 rotation, float scale, unsigned int seed_, float radius_)
    : Object("planetoid", name, wp.toFloat(), rotation, WHITE, scale),
      id(nextPlanetoidId++), seed(seed_), radius(radius_),
      field(std::make_shared<PlanetoidField>(seed_, radius_)),
      edits(std::make_shared<EditList>()) {
    worldPos = wp;
    hasWorldPos = true;
    // Smallest level whose 2x2x2 root nodes (coords -1..0) cover the body plus room for added material.
    while ((float)(CHUNK_SIZE << rootLevel) < radius + 24.0f) rootLevel++;
    planetoidRegistry[id] = this;
}

Planetoid::~Planetoid() { planetoidRegistry.erase(id); }

Planetoid::Planetoid(std::string name, DVec3 wp, std::shared_ptr<const TerrainField> terrain)
    : Object("terrain", name, wp.toFloat(), { 0, 0, 0 }, WHITE, 1.0f),
      id(nextPlanetoidId++), seed(0), radius(1e9f), field(terrain), edits(std::make_shared<EditList>()) {
    worldPos = wp;
    hasWorldPos = true;
    flat = true;
    minHeight = terrain->minHeight;
    maxHeight = terrain->maxHeight;
    rootLevel = 4; // 512-voxel root tiles, a 5x5 grid around the camera (~1.3 km view range)
    planetoidRegistry[id] = this;
}

Planetoid* Planetoid::find(uint64_t id) {
    auto it = planetoidRegistry.find(id);
    return it == planetoidRegistry.end() ? nullptr : it->second;
}

Vector3 Planetoid::toLocal(Vector3 world) const { return Vector3Transform(world, MatrixInvert(localTransform())); }
Vector3 Planetoid::toWorld(Vector3 local) const { return Vector3Transform(local, localTransform()); }

float Planetoid::densityAt(Vector3 local) const {
    return field->density(local, edits->data(), edits->size());
}

Vector3 Planetoid::surfaceVelocity(Vector3 renderPoint) const {
    if (spin.x == 0.0f && spin.y == 0.0f && spin.z == 0.0f) return { 0, 0, 0 };
    // Differentiate the body transform numerically: where will this body-fixed point be in h seconds?
    // Works for any spin axis (rotation is stored as Euler angles, so no closed-form omega).
    const float h = 0.01f;
    const Vector3 local = toLocal(renderPoint);
    const Vector3 r1 = Vector3Add(rotation, Vector3Scale(spin, h));
    Matrix m = MatrixScale(scale, scale, scale);
    m = MatrixMultiply(m, MatrixRotateXYZ({ DEG2RAD * r1.x, DEG2RAD * r1.y, DEG2RAD * r1.z }));
    m = MatrixMultiply(m, MatrixTranslate(position.x, position.y, position.z));
    return Vector3Scale(Vector3Subtract(Vector3Transform(local, m), renderPoint), 1.0f / h);
}

Vector3 Planetoid::dirToLocal(Vector3 worldDir) const {
    return Vector3Normalize(Vector3Subtract(toLocal(Vector3Add(position, worldDir)), toLocal(position)));
}

Vector3 Planetoid::dirToWorld(Vector3 localDir) const {
    return Vector3Normalize(Vector3Subtract(toWorld(localDir), toWorld({ 0, 0, 0 })));
}

Quaternion Planetoid::rotationQuat() const {
    return QuaternionFromMatrix(MatrixRotateXYZ({ DEG2RAD * rotation.x, DEG2RAD * rotation.y, DEG2RAD * rotation.z }));
}

void Planetoid::draw(const RenderContext& ctx) {
    if (!isActive) return;
    const Matrix t = localTransform();
    for (const auto& objPtr : children) objPtr->draw(ctx);
    for (Chunk* c : drawList) c->drawChunk(ctx, t);
}

void Planetoid::drawDepthOnly(const RenderContext& ctx) {
    if (!isActive) return;
    const Matrix t = localTransform();
    for (const auto& objPtr : children) objPtr->drawDepthOnly(ctx);
    for (Chunk* c : drawList) c->drawChunkDepth(ctx, t);
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
    std::shared_ptr<const EditList> edits = p.edits;
    const uint64_t id = p.id;
    Jobs().submit(priority, [field, edits, id, key, version] {
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

// Octree descent from the 8 root nodes: split while the camera is closer than lodFactor * node size.
static void CollectDesired(const Planetoid& p, const NodeKey& key, Vector3 camLocal, std::vector<NodeKey>& out) {
    const Vector3 lo = key.origin();
    const Vector3 hi = Vector3AddValue(lo, key.size());
    if (p.flat) {
        if (hi.y < p.minHeight - 8.0f || lo.y > p.maxHeight + 8.0f) return;      // outside the terrain band
    } else if (DistanceToBox({ 0, 0, 0 }, lo, hi) > p.radius + 24.0f) return;  // no surface possible here
    if (key.level > 0 && DistanceToBox(camLocal, lo, hi) < lodSettings.lodFactor * key.size()) {
        for (int i = 0; i < 8; ++i)
            CollectDesired(p, { key.level - 1, key.x * 2 + (i & 1), key.y * 2 + ((i >> 1) & 1), key.z * 2 + ((i >> 2) & 1) }, camLocal, out);
        return;
    }
    out.push_back(key);
}

// Decide which nodes to build and draw. Old nodes stay on screen until every desired node that
// overlaps them is ready, so LOD changes never open holes.
static void UpdateLod(Planetoid& p, Vector3 cameraRender) {
    PROFILE_FUNCTION();
    const Vector3 camLocal = p.toLocal(cameraRender);
    std::vector<NodeKey> desired;
    if (p.flat) {
        // Open world: a grid of root tiles around the camera, spanning the terrain height band.
        const float rs = (float)(CHUNK_SIZE << p.rootLevel);
        const int cx = (int)floorf(camLocal.x / rs), cz = (int)floorf(camLocal.z / rs);
        const int y0 = (int)floorf(p.minHeight / rs), y1 = (int)floorf(p.maxHeight / rs);
        for (int z = cz - 2; z <= cz + 2; ++z)
            for (int x = cx - 2; x <= cx + 2; ++x)
                for (int y = y0; y <= y1; ++y) CollectDesired(p, { p.rootLevel, x, y, z }, camLocal, desired);
    } else {
        for (int i = 0; i < 8; ++i)
            CollectDesired(p, { p.rootLevel, -1 + (i & 1), -1 + ((i >> 1) & 1), -1 + ((i >> 2) & 1) }, camLocal, desired);
    }
    std::unordered_set<NodeKey> desiredSet(desired.begin(), desired.end());

    for (const NodeKey& k : desired) {
        if (p.requestedVersion.find(k) != p.requestedVersion.end()) continue;
        const Vector3 lo = k.origin();
        // Nearest and coarsest first: fill the view quickly, then refine.
        const float dist = DistanceToBox(camLocal, lo, Vector3AddValue(lo, k.size()));
        RequestNode(p, k, dist / (1.0f + k.level));
    }

    // Previously drawn nodes that are no longer desired, kept while a desired overlap is pending.
    std::vector<Chunk*> retained;
    for (Chunk* c : p.drawList) {
        if (desiredSet.count(c->key)) continue;
        for (const NodeKey& d : desired) {
            if (d.overlaps(c->key) && !p.nodes.at(d)->ready) { retained.push_back(c); break; }
        }
    }
    std::vector<Chunk*> draw = retained;
    for (const NodeKey& d : desired) {
        Chunk* c = p.nodes.at(d).get();
        if (!c->ready) continue;
        bool coveredByCoarser = false;
        for (Chunk* r : retained) if (r->key.level > d.level && r->key.contains(d)) { coveredByCoarser = true; break; }
        if (!coveredByCoarser) draw.push_back(c);
    }
    p.drawList = std::move(draw);

    // Clip the one-cell overlap on axes where the +neighbour is drawn at the same level (it would
    // duplicate that neighbour's surface and z-fight); keep it toward other levels to hide cracks.
    std::unordered_set<NodeKey> drawn;
    for (Chunk* c : p.drawList) drawn.insert(c->key);
    for (Chunk* c : p.drawList) {
        const NodeKey k = c->key;
        const Vector3 o = k.origin();
        const float s = k.size() + 0.02f; // tiny overlap so the seam itself never cracks
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
    buildCount++;
    if (r.overflow) logger.logf("[ERROR] Node L%d (%d,%d,%d) exceeds 65535 vertices; skipped\n", r.key.level, r.key.x, r.key.y, r.key.z);
    it->second->setMesh(r);
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

void advanceBodies(Scene& world, float dt) {
    // Runs at the start of the frame, before the player/ship update, so riders follow this
    // frame's transform (spinning here after they moved left them one frame behind the surface).
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) p->rotation = Vector3Add(p->rotation, Vector3Scale(p->spin, dt));
    }
}

void worldUpdate(Scene& world, float uploadBudgetMs) {
    PROFILE_FUNCTION();

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
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) UpdateLod(*p, world.camera.position);
    }

    stats.queuedJobs = (int)Jobs().queued();
    stats.runningJobs = Jobs().running();
    stats.pendingUploads = (int)pendingUploads.size();
    stats.uploadsLastFrame = uploads;
    stats.workerThreads = (int)Jobs().threadCount();
    stats.avgBuildMs = buildCount ? (float)(buildMsTotal / buildCount) : 0.0f;
    stats.chunks = stats.nodes = stats.edits = 0;
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            stats.nodes += (int)p->nodes.size();
            for (Chunk* c : p->drawList) stats.chunks += c->hasMesh() ? 1 : 0;
            stats.edits += (int)p->edits->size();
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

        // Bounding sphere (radius + room for added material).
        float t = 0.0f, tEnd = localLen;
        if (!p->flat) {
            const float R = p->radius + 16.0f;
            const float b = Vector3DotProduct(o, d);
            const float c = Vector3DotProduct(o, o) - R * R;
            const float disc = b * b - c;
            if (disc < 0.0f) continue;
            t = std::max(0.0f, -b - sqrtf(disc));
            tEnd = std::min(localLen, -b + sqrtf(disc));
        }

        const float step = 0.5f;
        float prevT = t;
        if (p->densityAt(Vector3Add(o, Vector3Scale(d, t))) >= ISO_LEVEL) continue; // starting inside rock
        for (; t <= tEnd; t += step) {
            if (p->densityAt(Vector3Add(o, Vector3Scale(d, t))) >= ISO_LEVEL) {
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
                    hit.normal = Vector3Normalize(Vector3Subtract(p->toWorld(Vector3Add(local, p->field->normal(local, p->edits->data(), p->edits->size()))), worldPoint));
                    hit.material = p->field->material(local, p->field->baseDensity(local));
                }
                break;
            }
            prevT = t;
        }
    }
    return found;
}

void ApplyTerrainEdit(Planetoid& p, const TerrainEdit& edit) {
    auto next = std::make_shared<EditList>(*p.edits);
    next->push_back(edit);
    p.edits = next;

    // Rebuild every existing node (any level) whose samples or normals the brush reaches. Nodes
    // that do not exist yet are built later with the new edit list anyway.
    std::vector<NodeKey> touched;
    for (auto& entry : p.nodes) {
        const NodeKey& k = entry.first;
        const Vector3 lo = k.origin();
        const Vector3 hi = Vector3AddValue(lo, k.size() * (CHUNK_SIZE + 1.0f) / CHUNK_SIZE);
        if (DistanceToBox(edit.center, lo, hi) <= edit.radius + k.stride()) touched.push_back(k);
    }
    for (const NodeKey& k : touched) RequestNode(p, k, -1.0f - (float)k.level * 0.01f); // ahead of background work
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
        if (!b) continue;
        const float d = b->flat ? -1e30f : Vector3Distance(p, b->position) - b->radius * b->scale; // terrain is always "under you"
        if (d < bestDist) { bestDist = d; best = b; }
    }
    return best;
}


// ---------------------------------------------------------------------------------------------
// Generation
// ---------------------------------------------------------------------------------------------

void worldHandler(Scene& world) {
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

void generateWorld(Scene& world) {
    PROFILE_FUNCTION();
    worldGen.lastSeed = worldGen.seed ? worldGen.seed : rd();
    gen.seed(worldGen.lastSeed);
    logger.logf("Generating world (seed %u)...\n", worldGen.lastSeed);

    std::uniform_real_distribution<double> unit(0.0, 1.0);
    for (int i = 0; i < worldGen.planetoidCount; ++i) {
        // A loose "belt": random directions, distances up to systemRadius. Doubles keep these exact.
        DVec3 pos = { (unit(gen) - 0.5) * 2.0 * worldGen.systemRadius, (unit(gen) - 0.5) * 0.4 * worldGen.systemRadius,
                      (unit(gen) - 0.5) * 2.0 * worldGen.systemRadius };
        if (i == 0) pos = DVec3(worldGen.firstPosition); // keep the first one somewhere predictable
        const unsigned int seed = (unsigned int)gen();
        // Sizes vary per body around the configured size (the first keeps it exactly).
        const float size = (i == 0) ? (float)worldGen.planetoidSize : (float)worldGen.planetoidSize * (float)(0.5 + 2.5 * unit(gen));
        auto p = std::make_unique<Planetoid>("planetoid_" + std::to_string(i), pos, Vector3{ 0, 0, 0 }, 1.0f, seed, size);
        p->spin = { 0.0f, (float)(unit(gen) - 0.5) * 2.0f * worldGen.maxSpinDegPerSec, 0.0f };
        if (i == 0) p->spin = { 0, 0, 0 }; // the start body holds still (easier to test edits)
        p->parent = &world.rootObject;
        world.objects.push_back(std::move(p));
    }
    world.syncRenderPositions();
    generationStart = GetTime();
    buildMsTotal = 0.0;
    buildCount = 0;
}

void BeginWorldStreaming() { generationStart = GetTime(); buildMsTotal = 0.0; buildCount = 0; }

void regenerateWorld(Scene& world) {
    Jobs().clearQueued();  // in-flight jobs finish; their results are dropped (planetoid ids are gone)
    world.objects.clear(); // releases chunk meshes
    generateWorld(world);
}

// ---------------------------------------------------------------------------------------------
// Save / load: the world is fully determined by the generation settings + each planetoid's
// edit list, so a save is tiny and replays exactly.
// ---------------------------------------------------------------------------------------------

static const char* SAVE_PATH = "saves/world.sav";

bool SaveWorld(Scene& world) {
    if (!DirectoryExists("saves")) MakeDirectory("saves");
    FILE* f = fopen(SAVE_PATH, "w");
    if (!f) return false;
    const DVec3 camWorld = world.renderOrigin + world.camera.position;
    const DVec3 tgtWorld = world.renderOrigin + world.camera.target;
    fprintf(f, "RAYLIBGAME_SAVE 3\n");
    fprintf(f, "seed %u\ncount %d\nsize %d\nsystem %f %f\nfirst %f %f %f\n", worldGen.lastSeed, worldGen.planetoidCount,
            worldGen.planetoidSize, worldGen.systemRadius, worldGen.maxSpinDegPerSec,
            worldGen.firstPosition.x, worldGen.firstPosition.y, worldGen.firstPosition.z);
    fprintf(f, "camera %.6f %.6f %.6f %.6f %.6f %.6f\n", camWorld.x, camWorld.y, camWorld.z, tgtWorld.x, tgtWorld.y, tgtWorld.z);
    // Game-specific state (player, inventory, ...) written by whoever registered the hook.
    fprintf(f, "extras\n");
    if (saveExtrasHook) saveExtrasHook(f);
    fprintf(f, "end_extras\n");
    int index = 0, total = 0;
    for (auto& obj : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p) continue;
        fprintf(f, "planetoid %d rot %f %f %f edits %zu\n", index++, p->rotation.x, p->rotation.y, p->rotation.z, p->edits->size());
        for (const TerrainEdit& e : *p->edits) {
            fprintf(f, "%f %f %f %f %f\n", e.center.x, e.center.y, e.center.z, e.radius, e.strength);
            total++;
        }
    }
    fclose(f);
    logger.logf("[World] saved %d edits to %s\n", total, SAVE_PATH);
    return true;
}

bool LoadWorldSave(Scene& world) {
    FILE* f = fopen(SAVE_PATH, "r");
    if (!f) return false;
    int version = 0;
    if (fscanf(f, "RAYLIBGAME_SAVE %d", &version) != 1 || version < 2 || version > 3) {
        fclose(f);
        logger.logf("[World] %s: unsupported save version %d\n", SAVE_PATH, version);
        return false;
    }
    unsigned int seed = 0;
    int count = 1, size = 50;
    float systemRadius = 0, spin = 0;
    Vector3 first{};
    double c[6];
    bool ok = fscanf(f, " seed %u count %d size %d system %f %f first %f %f %f", &seed, &count, &size, &systemRadius, &spin, &first.x, &first.y, &first.z) == 8;
    ok = ok && fscanf(f, " camera %lf %lf %lf %lf %lf %lf", &c[0], &c[1], &c[2], &c[3], &c[4], &c[5]) == 6;
    if (!ok) { fclose(f); return false; }

    worldGen.seed = seed;
    worldGen.planetoidCount = count;
    worldGen.planetoidSize = size;
    worldGen.systemRadius = systemRadius;
    worldGen.maxSpinDegPerSec = spin;
    worldGen.firstPosition = first;
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
    while (fscanf(f, " planetoid %d rot %f %f %f edits %zu", &index, &rot.x, &rot.y, &rot.z, &editCount) == 5) {
        auto list = std::make_shared<EditList>();
        for (size_t i = 0; i < editCount; ++i) {
            TerrainEdit e{};
            if (fscanf(f, " %f %f %f %f %f", &e.center.x, &e.center.y, &e.center.z, &e.radius, &e.strength) != 5) break;
            list->push_back(e);
        }
        if (index >= 0 && index < (int)planetoids.size()) {
            planetoids[index]->edits = list;
            planetoids[index]->rotation = rot;
        }
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

