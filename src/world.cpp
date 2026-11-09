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

const int CHUNK_SIZE = 32;

WorldGenSettings worldGen;

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

Planetoid::Planetoid(std::string name, Vector3 position, Vector3 rotation, float scale, unsigned int seed_, float radius_)
    : Object("planetoid", name, position, rotation, WHITE, scale),
      id(nextPlanetoidId++), seed(seed_), radius(radius_),
      field(std::make_shared<PlanetoidField>(seed_, radius_)),
      edits(std::make_shared<EditList>()) {
    planetoidRegistry[id] = this;
}

Planetoid::~Planetoid() { planetoidRegistry.erase(id); }

Planetoid* Planetoid::find(uint64_t id) {
    auto it = planetoidRegistry.find(id);
    return it == planetoidRegistry.end() ? nullptr : it->second;
}

Vector3 Planetoid::toLocal(Vector3 world) const { return Vector3Transform(world, MatrixInvert(localTransform())); }
Vector3 Planetoid::toWorld(Vector3 local) const { return Vector3Transform(local, localTransform()); }

float Planetoid::densityAt(Vector3 local) const {
    return field->density(local, edits->data(), edits->size());
}

void Planetoid::draw(const RenderContext& ctx) {
    if (!isActive) return;
    const Matrix t = localTransform();
    for (const auto& objPtr : children) objPtr->draw(ctx);
    for (const auto& entry : chunks) entry.second->drawChunk(ctx, t);
}

void Planetoid::drawDepthOnly(const RenderContext& ctx) {
    if (!isActive) return;
    const Matrix t = localTransform();
    for (const auto& objPtr : children) objPtr->drawDepthOnly(ctx);
    for (const auto& entry : chunks) entry.second->drawChunkDepth(ctx, t);
}

// ---------------------------------------------------------------------------------------------
// Chunk requests / results
// ---------------------------------------------------------------------------------------------

static void RequestChunk(Planetoid& p, Int3 coord, float priority) {
    const uint32_t version = ++p.requestedVersion[coord];
    std::shared_ptr<const PlanetoidField> field = p.field;
    std::shared_ptr<const EditList> edits = p.edits;
    const uint64_t id = p.id;
    Jobs().submit(priority, [field, edits, id, coord, version] {
        PROFILE_SCOPE("Build chunk");
        ChunkMeshData data;
        data.planetoidId = id;
        data.coord = coord;
        data.version = version;
        BuildChunkMesh(*field, *edits, data);
        std::lock_guard<std::mutex> lock(resultsMutex);
        results.push_back(std::move(data));
    });
}

// Every chunk whose box can contain surface: the density is 0 beyond `radius`.
static void RequestInitialChunks(Planetoid& p, Vector3 cameraPos) {
    const int r = (int)ceilf(p.radius / CHUNK_SIZE) + 1;
    int requested = 0;
    for (int z = -r; z < r; ++z)
        for (int y = -r; y < r; ++y)
            for (int x = -r; x < r; ++x) {
                const Vector3 lo = { (float)(x * CHUNK_SIZE), (float)(y * CHUNK_SIZE), (float)(z * CHUNK_SIZE) };
                const Vector3 hi = Vector3AddValue(lo, (float)CHUNK_SIZE);
                const Vector3 closest = { std::max(lo.x, std::min(0.0f, hi.x)), std::max(lo.y, std::min(0.0f, hi.y)), std::max(lo.z, std::min(0.0f, hi.z)) };
                if (Vector3Length(closest) > p.radius) continue;
                // Nearest chunks first so the area around the player appears first.
                const Vector3 centreWorld = p.toWorld(Vector3AddValue(lo, CHUNK_SIZE * 0.5f));
                RequestChunk(p, { x, y, z }, Vector3Distance(centreWorld, cameraPos));
                requested++;
            }
    p.initialRequestDone = true;
    logger.logf("[World] %s: requested %d chunks\n", p.name.c_str(), requested);
}

static bool PipelineIdle();

static void IntegrateResult(const ChunkMeshData& r) {
    Planetoid* p = Planetoid::find(r.planetoidId);
    if (!p) return;                                                   // planetoid was destroyed
    auto v = p->requestedVersion.find(r.coord);
    if (v == p->requestedVersion.end() || v->second != r.version) return; // superseded by a newer request

    buildMsTotal += r.buildMs;
    buildCount++;
    if (r.overflow) logger.logf("[ERROR] Chunk (%d,%d,%d) exceeds 65535 vertices; skipped\n", r.coord.x, r.coord.y, r.coord.z);

    auto it = p->chunks.find(r.coord);
    if (r.positions.empty()) {
        if (it != p->chunks.end()) p->chunks.erase(it); // dug away entirely
        return;
    }
    if (it == p->chunks.end()) {
        auto chunk = std::make_unique<Chunk>();
        chunk->coord = r.coord;
        chunk->parent = p;
        it = p->chunks.emplace(r.coord, std::move(chunk)).first;
    }
    it->second->setMesh(r);
}

void worldUpdate(Scene& world, float uploadBudgetMs) {
    PROFILE_FUNCTION();
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            if (!p->initialRequestDone) RequestInitialChunks(*p, world.camera.position);
        }
    }

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

    stats.queuedJobs = (int)Jobs().queued();
    stats.runningJobs = Jobs().running();
    stats.pendingUploads = (int)pendingUploads.size();
    stats.uploadsLastFrame = uploads;
    stats.workerThreads = (int)Jobs().threadCount();
    stats.avgBuildMs = buildCount ? (float)(buildMsTotal / buildCount) : 0.0f;
    stats.chunks = 0;
    stats.edits = 0;
    for (auto& obj : world.objects) {
        if (auto* p = dynamic_cast<Planetoid*>(obj.get())) {
            stats.chunks += (int)p->chunks.size();
            stats.edits += (int)p->edits->size();
        }
    }

    // Initial requests were issued above, so "pipeline idle" here really means "everything built".
    if (generationStart >= 0.0 && PipelineIdle()) {
        stats.lastGenerationSeconds = (float)(GetTime() - generationStart);
        logger.logf("[World] ready in %.2f s (%d chunks, avg build %.2f ms on %d threads)\n",
                    stats.lastGenerationSeconds, stats.chunks, stats.avgBuildMs, stats.workerThreads);
        TraceLog(LOG_WARNING, "World ready in %.2f s (%d chunks)", stats.lastGenerationSeconds, stats.chunks);
        generationStart = -1.0;
    }
}

static bool PipelineIdle() {
    std::lock_guard<std::mutex> lock(resultsMutex);
    return Jobs().idle() && results.empty() && pendingUploads.empty();
}

// False until a (re)generation has fully streamed in, and while edits are being remeshed.
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
        const float R = p->radius + 16.0f;
        const float b = Vector3DotProduct(o, d);
        const float c = Vector3DotProduct(o, o) - R * R;
        const float disc = b * b - c;
        if (disc < 0.0f) continue;
        float t = std::max(0.0f, -b - sqrtf(disc));
        const float tEnd = std::min(localLen, -b + sqrtf(disc));

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

    // Remesh every chunk whose samples (or normals, +1) the brush touches, including chunks that
    // do not exist yet (adding material into empty space).
    const float reach = edit.radius + 1.0f;
    const int x0 = (int)floorf((edit.center.x - reach) / CHUNK_SIZE), x1 = (int)floorf((edit.center.x + reach) / CHUNK_SIZE);
    const int y0 = (int)floorf((edit.center.y - reach) / CHUNK_SIZE), y1 = (int)floorf((edit.center.y + reach) / CHUNK_SIZE);
    const int z0 = (int)floorf((edit.center.z - reach) / CHUNK_SIZE), z1 = (int)floorf((edit.center.z + reach) / CHUNK_SIZE);
    for (int z = z0; z <= z1; ++z)
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x)
                RequestChunk(p, { x, y, z }, -1.0f); // ahead of any background generation
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

    std::uniform_real_distribution<float> unit(0.0f, 1.0f);
    const float genRange = 500.0f;
    for (int i = 0; i < worldGen.planetoidCount; ++i) {
        Vector3 position = { (unit(gen) - 0.5f) * genRange, (unit(gen) - 0.5f) * genRange, (unit(gen) - 0.5f) * genRange };
        if (i == 0) position = worldGen.firstPosition; // keep the first one somewhere predictable
        const unsigned int seed = (unsigned int)gen();
        auto p = std::make_unique<Planetoid>("planetoid_" + std::to_string(i), position, Vector3{ 0, 0, 0 }, 1.0f,
                                             seed, (float)worldGen.planetoidSize);
        p->parent = &world.rootObject;
        world.objects.push_back(std::move(p));
    }
    // Chunks are requested on the next worldUpdate(), nearest to the camera first.
    generationStart = GetTime();
    buildMsTotal = 0.0;
    buildCount = 0;
}

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
    fprintf(f, "RAYLIBGAME_SAVE 1\n");
    fprintf(f, "seed %u\ncount %d\nsize %d\nfirst %f %f %f\n", worldGen.lastSeed, worldGen.planetoidCount,
            worldGen.planetoidSize, worldGen.firstPosition.x, worldGen.firstPosition.y, worldGen.firstPosition.z);
    fprintf(f, "camera %f %f %f %f %f %f\n", world.camera.position.x, world.camera.position.y, world.camera.position.z,
            world.camera.target.x, world.camera.target.y, world.camera.target.z);
    int index = 0, total = 0;
    for (auto& obj : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p) continue;
        fprintf(f, "planetoid %d edits %zu\n", index++, p->edits->size());
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
    if (fscanf(f, "RAYLIBGAME_SAVE %d", &version) != 1 || version != 1) { fclose(f); return false; }
    unsigned int seed = 0;
    int count = 1, size = 50;
    Vector3 first{}, camPos{}, camTarget{};
    bool ok = fscanf(f, " seed %u count %d size %d first %f %f %f", &seed, &count, &size, &first.x, &first.y, &first.z) == 6;
    ok = ok && fscanf(f, " camera %f %f %f %f %f %f", &camPos.x, &camPos.y, &camPos.z, &camTarget.x, &camTarget.y, &camTarget.z) == 6;
    if (!ok) { fclose(f); return false; }

    worldGen.seed = seed;
    worldGen.planetoidCount = count;
    worldGen.planetoidSize = size;
    worldGen.firstPosition = first;
    regenerateWorld(world);
    world.camera.position = camPos;
    world.camera.target = camTarget;

    // Edits are attached before the first worldUpdate() requests chunks, so each chunk is built
    // once with its edits already applied.
    std::vector<Planetoid*> planetoids;
    for (auto& obj : world.objects) if (auto* p = dynamic_cast<Planetoid*>(obj.get())) planetoids.push_back(p);
    int index = 0;
    size_t editCount = 0;
    while (fscanf(f, " planetoid %d edits %zu", &index, &editCount) == 2) {
        auto list = std::make_shared<EditList>();
        for (size_t i = 0; i < editCount; ++i) {
            TerrainEdit e{};
            if (fscanf(f, " %f %f %f %f %f", &e.center.x, &e.center.y, &e.center.z, &e.radius, &e.strength) != 5) break;
            list->push_back(e);
        }
        if (index >= 0 && index < (int)planetoids.size()) planetoids[index]->edits = list;
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

