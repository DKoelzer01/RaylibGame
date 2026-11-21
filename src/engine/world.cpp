#include "world.h"
#include "netid.h"
#include "gamerng.h"
#include "terrainarena.h"
#include "ring.h"
#include "player.h"
#include "jobsystem.h"
#include "utils/logger.h"
#include "worldtypes.h"
#include "factory.h"
#include "fluid.h"
#include "rivers.h"
#include "terraingpu.h"
#include "gputimer.h"
#include "scenetarget.h"
#include "waterfield.h"
#include <chrono>
#include <filesystem>

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
#include <atomic>
#include <chrono>
#include <thread>
#include "erosion.h"

const int CHUNK_SIZE = 32;

WorldGenSettings worldGen;
std::function<void(FILE*)> saveExtrasHook;
std::function<void(FILE*)> loadExtrasHook;
std::function<void()> loadDoneHook;
static std::vector<std::function<void()>>& UnloadHooks() { static std::vector<std::function<void()>> h; return h; }
void World_AddUnloadHook(std::function<void()> hook) { UnloadHooks().push_back(std::move(hook)); }
LodSettings lodSettings;
std::unordered_set<const Chunk*> gOcclusionSkip; // #134 experiment (main.cpp RAYLIBGAME_OCC_SKIP)
std::unordered_set<uint32_t> gOcclusionSkipPlots; // #134 experiment: plot ids left out (parts and machines)

namespace {
    std::random_device rd;
    std::mt19937 gen(rd()); // reseeded by generateWorld()

    // ---- Async chunk pipeline -------------------------------------------------------------
    // Declared before the job system so it is destroyed after the workers are joined.
    std::mutex resultsMutex;
    std::atomic<long long> gRequests{ 0 }, gBuilt{ 0 }, gSkipped{ 0 }, gDropped{ 0 }; // streaming counters (#258)
    long long gLodFades = 0; // nodes faded out into their replacements (#243)
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
    std::vector<ChunkMeshData> pendingUploads; // main thread only; unordered (#387: worldUpdate picks the nearest each frame)
    constexpr size_t kUploadBacklogCap = 4096; // results waiting for upload past which no new GPU batches go out (#387)
    double gUploadMaxMs = 0.0, gUploadLoopMs = 0.0, gWaterUpdMs = 0.0; // the [jobs] line: the slowest upload, the loop total, the slowest UpdateWater (#509)
    std::vector<ChunkMeshData> waterResults;   // water-only rebuilds (#360), filled by workers (resultsMutex)
    std::atomic<long long> gWaterBuilds{ 0 };  // water-only rebuilds done, and their worker time (ms x 1000)
    std::atomic<long long> gWaterUs{ 0 };
    std::atomic<int> gWaterInFlight{ 0 };      // water-only rebuilds queued or running (capped: they go ahead of view builds)
    std::atomic<long long> gWaterFullBuilds{ 0 }; // full chunk builds asked for by water (a node without its water cache)
    struct WaterToken { WaterToken() { ++gWaterInFlight; } ~WaterToken() { --gWaterInFlight; } }; // held by a water job's closure (run or dropped)
    std::unique_ptr<JobSystem> jobs;

    std::unordered_map<uint64_t, Planetoid*> planetoidRegistry; // main thread only
    uint64_t nextPlanetoidId = 1;

    WorldStats stats;
    double generationStart = -1.0; // GetTime() when a (re)generation started; <0 once reported
    double buildMsTotal = 0.0, sampleMsTotal = 0.0;
    long buildCount = 0;

    JobSystem& Jobs() {
        if (!jobs) { jobs = std::make_unique<JobSystem>(); jobs->setBackgroundFrom(kBackgroundPriority); } // ring builds are background (#296)
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

uint64_t WorldFingerprint(const Scene& world) {
    // FNV-1a over the seed and, per generated body, its id, kind and size plus the CPU density field sampled at fixed
    // directions and heights (no edits: that is what every machine regenerates from the seed). Float bits are hashed as is.
    uint64_t h = 1469598103934665603ull;
    auto mix = [&](const void* data, size_t n) { const uint8_t* b = (const uint8_t*)data; for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; } };
    mix(&worldGen.lastSeed, sizeof(worldGen.lastSeed));
    for (const auto& obj : world.objects) {
        const auto* p = dynamic_cast<const Planetoid*>(obj.get());
        if (!p || NetId_Kind(p->netId) != NetKind::Body || !p->field) continue;
        mix(&p->netId, sizeof(p->netId));
        const int kind = (int)p->kind;
        mix(&kind, sizeof(kind));
        mix(&p->radius, sizeof(p->radius));
        // Heights where the surface is: around sea level (planets, moons), else through the outer half (asteroids).
        const float sea = p->field->seaLevelRadius(), top = p->field->maxRadius();
        const float heights[3] = { sea > 0.0f ? sea - 40.0f : top * 0.5f, sea > 0.0f ? sea + 15.0f : top * 0.75f, sea > 0.0f ? sea + 150.0f : top * 0.95f };
        for (int i = 0; i < 64; ++i) {
            // A Fibonacci sphere of directions, three heights each.
            const float y = 1.0f - 2.0f * (i + 0.5f) / 64.0f, rr = sqrtf(std::max(0.0f, 1.0f - y * y)), a = 2.39996323f * (float)i;
            const Vector3 dir = { rr * cosf(a), y, rr * sinf(a) };
            for (float hgt : heights) {
                const float r = p->flat ? 0.0f : hgt;
                const float d = p->field->baseDensity(Vector3Scale(dir, r));
                mix(&d, sizeof(d));
            }
        }
    }
    return h;
}

Planetoid* Planetoid::findNet(uint64_t netId) {
    if (netId == 0) return nullptr;
    for (const auto& kv : planetoidRegistry) if (kv.second->netId == netId) return kv.second;
    return nullptr;
}

std::vector<Planetoid*> Planetoid::All() {
    std::vector<Planetoid*> out;
    for (const auto& kv : planetoidRegistry) out.push_back(kv.second);
    return out;
}

// The rotation of localTransform() (raymath's MatrixRotateXYZ) in double: world = R * (scale * local).
static void BodyRotationD(const Vector3& rotationDeg, double r[3][3]) {
    const double ax = -(double)rotationDeg.x * DEG2RAD, ay = -(double)rotationDeg.y * DEG2RAD, az = -(double)rotationDeg.z * DEG2RAD;
    const double cx = cos(ax), sx = sin(ax), cy = cos(ay), sy = sin(ay), cz = cos(az), sz = sin(az);
    // Rows of R are (m0, m4, m8), (m1, m5, m9), (m2, m6, m10) of MatrixRotateXYZ.
    r[0][0] = cz * cy;                 r[0][1] = sz * cy;                 r[0][2] = -sy;
    r[1][0] = cz * sy * sx - sz * cx;  r[1][1] = sz * sy * sx + cz * cx;  r[1][2] = cy * sx;
    r[2][0] = cz * sy * cx + sz * sx;  r[2][1] = sz * sy * cx - cz * sx;  r[2][2] = cy * cx;
}

static void BodyRotationDd(double rx, double ry, double rz, double r[3][3]);

// Plan 11 F1 (#667): RAYLIBGAME_FRAMES=1 builds a spinning body's rotation from its double spin angle (advanceBodies keeps
// phase + spin * time in double) instead of the float Euler `rotation` stored from it: float degrees near 360 step by
// 3.05e-5 deg, up to ~3 cm at Terra's surface, and frames/still saw the standing player jump 13-31 mm in single steps.
// The double angle is used while it is current - advanceBodies wrote `rotation` from it - and the float angles otherwise
// (a load, a script or a test set `rotation`).
static bool kFramesPose = [] { const char* e = getenv("RAYLIBGAME_FRAMES"); return e && e[0] == '1'; }(); // (unit tests set it)
struct BodyAngles { double x, y, z; };
static BodyAngles AnglesOf(const Planetoid& p) {
    if (kFramesPose && p.spinAngleSet && (float)p.spinAngle[0] == p.rotation.x && (float)p.spinAngle[1] == p.rotation.y && (float)p.spinAngle[2] == p.rotation.z)
        return { p.spinAngle[0], p.spinAngle[1], p.spinAngle[2] };
    return { p.rotation.x, p.rotation.y, p.rotation.z };
}
static void BodyRotationOf(const Planetoid& p, double r[3][3]) { const BodyAngles a = AnglesOf(p); BodyRotationDd(a.x, a.y, a.z, r); }
// A rotation matrix (body-local -> world, rows of R) as raymath's Matrix, times a scale: MatrixRotateXYZ's layout.
static Matrix RotationMatrix(const double r[3][3], double s = 1.0) {
    return { (float)(r[0][0] * s), (float)(r[0][1] * s), (float)(r[0][2] * s), 0.0f,
             (float)(r[1][0] * s), (float)(r[1][1] * s), (float)(r[1][2] * s), 0.0f,
             (float)(r[2][0] * s), (float)(r[2][1] * s), (float)(r[2][2] * s), 0.0f,
             0.0f, 0.0f, 0.0f, 1.0f };
}

// The same rotation, cached per thread for the last body and angles: toLocal / toWorld run in hot loops (collision, tools,
// chunk jobs), and six double cos/sin a call would add up.
static const double (*BodyRotationCached(const Planetoid& p))[3] {
    thread_local const Planetoid* body = nullptr;
    thread_local BodyAngles rot{ NAN, NAN, NAN };
    thread_local double r[3][3];
    const BodyAngles a = AnglesOf(p);
    if (body != &p || rot.x != a.x || rot.y != a.y || rot.z != a.z) { body = &p; rot = a; BodyRotationDd(a.x, a.y, a.z, r); }
    return r;
}

// Render space <-> body-local in double (#1072): the float path (MatrixInvert(localTransform()), a float rotation applied to
// a translation ~60 km from the render origin) was off by mm to ~2 cm, and differently every step as the planet turned, so
// the player's collision - its solid test, push-out and ground probe, all sampling the density through toLocal - see-sawed
// while standing still and the camera bobbed. Now R^T (c - position) / scale with the rotation in double: the body's render
// position (a float) only changes its rounding when the render origin re-centres, so the samples hold still with the ground.
Vector3 Planetoid::toLocal(Vector3 world) const {
    const double (*r)[3] = BodyRotationCached(*this);
    const double dx = (double)world.x - position.x, dy = (double)world.y - position.y, dz = (double)world.z - position.z;
    const double inv = 1.0 / scale;
    return { (float)((r[0][0] * dx + r[1][0] * dy + r[2][0] * dz) * inv), (float)((r[0][1] * dx + r[1][1] * dy + r[2][1] * dz) * inv),
             (float)((r[0][2] * dx + r[1][2] * dy + r[2][2] * dz) * inv) };
}
Vector3 Planetoid::toWorld(Vector3 local) const {
    const double (*r)[3] = BodyRotationCached(*this);
    const double s = scale, lx = local.x, ly = local.y, lz = local.z;
    return { (float)(position.x + (r[0][0] * lx + r[0][1] * ly + r[0][2] * lz) * s), (float)(position.y + (r[1][0] * lx + r[1][1] * ly + r[1][2] * lz) * s),
             (float)(position.z + (r[2][0] * lx + r[2][1] * ly + r[2][2] * lz) * s) };
}

double SecondsToSunElevation(const Planetoid& body, Vector3 localPoint, Vector3 sunDir, float elevDeg, bool rising) {
    const Vector3 s = body.spin;
    const double rate = std::max({ fabs((double)s.x), fabs((double)s.y), fabs((double)s.z) });
    const double len = sqrt((double)localPoint.x * localPoint.x + (double)localPoint.y * localPoint.y + (double)localPoint.z * localPoint.z);
    if (rate < 1e-9 || len < 1e-6) return -1.0;
    const double period = 360.0 / rate;
    const double ux = localPoint.x / len, uy = localPoint.y / len, uz = localPoint.z / len;
    const double sl = sqrt((double)sunDir.x * sunDir.x + (double)sunDir.y * sunDir.y + (double)sunDir.z * sunDir.z);
    auto elevation = [&](double t) { // degrees, after t more seconds of spin
        const Vector3 rot = { (float)fmod(body.rotation.x + s.x * t, 360.0), (float)fmod(body.rotation.y + s.y * t, 360.0), (float)fmod(body.rotation.z + s.z * t, 360.0) };
        double r[3][3];
        BodyRotationD(rot, r);
        const double wx = r[0][0] * ux + r[0][1] * uy + r[0][2] * uz, wy = r[1][0] * ux + r[1][1] * uy + r[1][2] * uz, wz = r[2][0] * ux + r[2][1] * uy + r[2][2] * uz;
        const double d = (wx * sunDir.x + wy * sunDir.y + wz * sunDir.z) / std::max(sl, 1e-9);
        return asin(std::clamp(d, -1.0, 1.0)) * RAD2DEG;
    };
    const int n = 2880; // 30 s steps on a one-hour day
    double t0 = 0.0, e0 = elevation(0.0);
    for (int i = 1; i <= n; ++i) {
        const double t1 = period * i / n, e1 = elevation(t1);
        const bool crosses = (e0 - elevDeg) * (e1 - elevDeg) <= 0.0 && e0 != e1 && (rising ? e1 > e0 : e1 < e0);
        if (crosses && i > 1) { // (not the sample we start on: "dawn" pressed at dawn goes to the next one)
            double a = t0, b = t1;
            for (int k = 0; k < 40; ++k) { const double m = 0.5 * (a + b); if ((elevation(a) - elevDeg) * (elevation(m) - elevDeg) <= 0.0) b = m; else a = m; }
            return 0.5 * (a + b);
        }
        t0 = t1; e0 = e1;
    }
    return -1.0;
}

namespace {
    const Planetoid* NearestBody(const Scene& world, Vector3& local) {
        const Planetoid* best = nullptr;
        float bestGap = 1e30f;
        for (const auto& obj : world.objects) {
            const auto* p = dynamic_cast<const Planetoid*>(obj.get());
            if (!p || p->vehicle) continue;
            const Vector3 l = p->toLocal(world.camera.position);
            const float gap = Vector3Length(l) - p->radius;
            if (gap < bestGap) { bestGap = gap; best = p; local = l; }
        }
        return best;
    }
}

double SecondsToSunElevationAtCamera(const Scene& world, float elevDeg, bool rising) {
    Vector3 local{};
    const Planetoid* p = NearestBody(world, local);
    return p ? SecondsToSunElevation(*p, local, world.sunDirection, elevDeg, rising) : -1.0;
}

double BodyDaySeconds(const Scene& world) {
    Vector3 local{};
    const Planetoid* p = NearestBody(world, local);
    if (!p) return 0.0;
    const double rate = std::max({ fabs((double)p->spin.x), fabs((double)p->spin.y), fabs((double)p->spin.z) });
    return rate > 1e-9 ? 360.0 / rate : 0.0;
}

double BodyHourSeconds(const Scene& world) { return BodyDaySeconds(world) / 24.0; }

float SunElevationAtCamera(const Scene& world) {
    Vector3 local{};
    const Planetoid* p = NearestBody(world, local);
    if (!p || Vector3Length(local) < 1e-3f) return 0.0f;
    const Vector3 up = p->dirToWorld(Vector3Normalize(local));
    return asinf(std::clamp(Vector3DotProduct(up, Vector3Normalize(world.sunDirection)), -1.0f, 1.0f)) * RAD2DEG;
}

DVec3 Planetoid::toWorldD(const DVec3& l) const {
    double r[3][3];
    BodyRotationOf(*this, r);
    const double s = scale;
    return worldPos + DVec3((r[0][0] * l.x + r[0][1] * l.y + r[0][2] * l.z) * s,
                            (r[1][0] * l.x + r[1][1] * l.y + r[1][2] * l.z) * s,
                            (r[2][0] * l.x + r[2][1] * l.y + r[2][2] * l.z) * s);
}

DVec3 Planetoid::toLocalD(const DVec3& w) const {
    double r[3][3];
    BodyRotationOf(*this, r);
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
            BodyRotationOf(p, r);
            if (kFramesPose) rs = RotationMatrix(r, p.scale); // (#667) the vertices turned by the same double pose as the translation
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

VoxelMaterial Planetoid::materialAt(Vector3 local, Vector3 normal) const {
    if (edits) { const int m = edits->materialAt(local); if (m >= 0) return (VoxelMaterial)m; }
    if (Vector3LengthSqr(normal) < 1e-6f) return field->material(local, field->baseDensity(local));
    // #636: material() assumes level ground (fallbackUp), so on a cliff the crosshair (and the dig yield) said grass or
    // gravel where the mesh, picked with its real normal, draws rock.
    return field->materialLod(local, field->baseDensity(local), 0.0f, Vector3Normalize(normal));
}

float Planetoid::densityAt(Vector3 local) const {
    ++densityCalls; // (#768)
    return field->density(local, edits.get());
}

// BodyRotationD with double angles (deg): r maps body-local into world (dirToWorld's matrix). (#872)
static void BodyRotationDd(double rx, double ry, double rz, double r[3][3]) {
    const double ax = -rx * DEG2RAD, ay = -ry * DEG2RAD, az = -rz * DEG2RAD;
    const double cx = cos(ax), sx = sin(ax), cy = cos(ay), sy = sin(ay), cz = cos(az), sz = sin(az);
    r[0][0] = cz * cy;                 r[0][1] = sz * cy;                 r[0][2] = -sy;
    r[1][0] = cz * sy * sx - sz * cx;  r[1][1] = sz * sy * sx + cz * cx;  r[1][2] = cy * sx;
    r[2][0] = cz * sy * cx + sz * sx;  r[2][1] = sz * sy * cx - cz * sx;  r[2][2] = cy * cx;
}

// Omega (rad/s, world) of a body with these Euler angles turning at `spin` deg/s: the skew part of dR/dt * R^T, dR/dt
// by a central difference over +-1 ms of the angles, all in double (#872: the float quaternion route, 2 acos(w) with w
// within 4e-7 of 1, lost up to 1.5% for bodies spinning about several axes).
static DVec3 AngularVelocityD(Vector3 rotation, Vector3 spin) {
    constexpr double h = 1e-3;
    double r0[3][3], rp[3][3], rm[3][3];
    BodyRotationDd(rotation.x, rotation.y, rotation.z, r0);
    BodyRotationDd(rotation.x + (double)spin.x * h, rotation.y + (double)spin.y * h, rotation.z + (double)spin.z * h, rp);
    BodyRotationDd(rotation.x - (double)spin.x * h, rotation.y - (double)spin.y * h, rotation.z - (double)spin.z * h, rm);
    double w[3][3] = {};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            for (int k = 0; k < 3; ++k) w[i][j] += (rp[i][k] - rm[i][k]) / (2.0 * h) * r0[j][k];
    return DVec3{ 0.5 * (w[2][1] - w[1][2]), 0.5 * (w[0][2] - w[2][0]), 0.5 * (w[1][0] - w[0][1]) };
}

// Plan 11 (frames.h, #667): a body's orientation (rows: body-local -> world) and angular velocity (world, rad/s) in double,
// from the double spin angle while it is current (whatever RAYLIBGAME_FRAMES says: only frame code asks for this).
void BodyRotationAndSpin(const Planetoid& p, double r[3][3], DVec3& omega) {
    const bool current = p.spinAngleSet && (float)p.spinAngle[0] == p.rotation.x && (float)p.spinAngle[1] == p.rotation.y && (float)p.spinAngle[2] == p.rotation.z;
    const double ax = current ? p.spinAngle[0] : p.rotation.x, ay = current ? p.spinAngle[1] : p.rotation.y, az = current ? p.spinAngle[2] : p.rotation.z;
    BodyRotationDd(ax, ay, az, r);
    if (p.spin.x == 0.0f && p.spin.y == 0.0f && p.spin.z == 0.0f) { omega = DVec3{ 0, 0, 0 }; return; }
    omega = AngularVelocityD(p.rotation, p.spin);
}

Vector3 Planetoid::surfaceVelocity(Vector3 renderPoint) const {
    // The body's own motion (orbit) plus its spin: omega x (point - centre), omega taken as a vector (differentiating
    // positions in float lost ~0.4 m/s at 60 km), in double (#872).
    if (spin.x == 0.0f && spin.y == 0.0f && spin.z == 0.0f) return linearVelocity;
    const DVec3 w = AngularVelocityD(rotation, spin);
    const double dx = (double)renderPoint.x - position.x, dy = (double)renderPoint.y - position.y, dz = (double)renderPoint.z - position.z;
    return Vector3Add(linearVelocity, Vector3{ (float)(w.y * dz - w.z * dy), (float)(w.z * dx - w.x * dz), (float)(w.x * dy - w.y * dx) });
}

bool SurfaceVelocity_SelfTest() { // (#872) every orientation and spin: omega x r against the motion of a surface point
    uint32_t rng = 777; auto rnd = [&rng] { rng = rng * 1664525u + 1013904223u; return (rng >> 8) / 16777216.0f; };
    double worst = 0.0, worstY = 0.0, worstFrame = 0.0;
    for (int i = 0; i < 1000; ++i) {
        const Vector3 rot = { rnd() * 360.0f, rnd() * 360.0f, rnd() * 360.0f };
        const bool yOnly = i < 500;
        const Vector3 spin = yOnly ? Vector3{ 0, 0.004f + 0.02f * rnd(), 0 } : Vector3{ (rnd() - 0.5f) * 2.0f, (rnd() - 0.5f) * 2.0f, (rnd() - 0.5f) * 2.0f };
        const double L[3] = { (rnd() - 0.5) * 2.0, (rnd() - 0.5) * 2.0, (rnd() - 0.5) * 2.0 };
        const double len = sqrt(L[0] * L[0] + L[1] * L[1] + L[2] * L[2]);
        auto worldAt = [&](double t, double out[3]) { // the surface point (59 km out) at time t, the angles advanced by spin * t
            double m[3][3];
            BodyRotationDd(rot.x + spin.x * t, rot.y + spin.y * t, rot.z + spin.z * t, m);
            for (int a = 0; a < 3; ++a) out[a] = (m[a][0] * L[0] + m[a][1] * L[1] + m[a][2] * L[2]) / len * 59000.0;
        };
        double a[3], b[3], c[3];
        worldAt(-0.01, a); worldAt(0.01, b); worldAt(0.0, c);
        const DVec3 w = AngularVelocityD(rot, spin);
        const double v[3] = { w.y * c[2] - w.z * c[1], w.z * c[0] - w.x * c[2], w.x * c[1] - w.y * c[0] };
        const double ref[3] = { (b[0] - a[0]) / 0.02, (b[1] - a[1]) / 0.02, (b[2] - a[2]) / 0.02 };
        const double rl = sqrt(ref[0] * ref[0] + ref[1] * ref[1] + ref[2] * ref[2]);
        const double err = sqrt((v[0] - ref[0]) * (v[0] - ref[0]) + (v[1] - ref[1]) * (v[1] - ref[1]) + (v[2] - ref[2]) * (v[2] - ref[2])) / std::max(rl, 1e-9);
        if (yOnly) worstY = std::max(worstY, err); else worst = std::max(worst, err);
        // the frame: BodyRotationDd is dirToWorld's matrix (the body's real local -> world)
        double r[3][3]; BodyRotationD(rot, r); double m[3][3]; BodyRotationDd(rot.x, rot.y, rot.z, m);
        for (int x = 0; x < 3; ++x) for (int y = 0; y < 3; ++y) worstFrame = std::max(worstFrame, fabs(r[x][y] - m[x][y]));
    }
    const bool ok = worst < 1e-3 && worstY < 1e-3 && worstFrame < 1e-9;
    printf("bench: surface velocity: worst error %.5f%% (spin about Y), %.5f%% (about all three axes) over 500 orientations each; frame %.1e: %s\n",
           worstY * 100.0, worst * 100.0, worstFrame, ok ? "PASS" : "FAIL");
    return ok;
}

// Directions only rotate (no translation): transforming two points and subtracting lost ~0.2 deg in
// float with the body centre tens of km from the render origin (hundreds of metres on the surface).
Vector3 Planetoid::dirToLocal(Vector3 worldDir) const {
    double r[3][3];
    BodyRotationOf(*this, r);
    const Vector3 d = worldDir;
    return Vector3Normalize({ (float)(r[0][0] * d.x + r[1][0] * d.y + r[2][0] * d.z),
                              (float)(r[0][1] * d.x + r[1][1] * d.y + r[2][1] * d.z),
                              (float)(r[0][2] * d.x + r[1][2] * d.y + r[2][2] * d.z) });
}

Vector3 Planetoid::dirToWorld(Vector3 localDir) const {
    double r[3][3];
    BodyRotationOf(*this, r);
    const Vector3 d = localDir;
    return Vector3Normalize({ (float)(r[0][0] * d.x + r[0][1] * d.y + r[0][2] * d.z),
                              (float)(r[1][0] * d.x + r[1][1] * d.y + r[1][2] * d.z),
                              (float)(r[2][0] * d.x + r[2][1] * d.y + r[2][2] * d.z) });
}

Quaternion Planetoid::rotationQuat() const {
    if (!kFramesPose) return QuaternionFromMatrix(MatrixRotateXYZ({ DEG2RAD * rotation.x, DEG2RAD * rotation.y, DEG2RAD * rotation.z }));
    double r[3][3];
    BodyRotationOf(*this, r); // (#667) the double pose: a grounded player's view turns with it step by step, unquantised
    return QuaternionFromMatrix(RotationMatrix(r));
}

bool Planetoid::useFarSphere(const RenderContext& ctx) const {
    if (!farSphere || !ctx.farShader || !farSphere->isReady()) return false;
    return Vector3Distance(position, ctx.cameraPos) / std::max(1.0f, lodSettings.zoom) > kFarSphereRadii * radius * scale;
}

// #830: when parts were placed or removed, the grass scatter gets the body's new part cells and the chunks under the
// changes are rebuilt (their grass redone). One remesh per changed plot, over the changed cells' extent.
static void SyncGrassBlockers(Planetoid& p) {
    bool changed = false;
    for (const auto& pl : p.plots->plots) changed = changed || !pl->changedCells.empty();
    if (!changed) return;
    auto snap = std::make_shared<GrassBlockers>();
    for (const auto& pl : p.plots->plots) {
        if (pl->cells.empty()) continue;
        GrassBlockPlot g;
        g.origin = pl->origin; g.east = pl->east; g.up = pl->up; g.north = pl->north;
        for (const auto& kv : pl->cells) {
            g.cells.insert(kv.first);
            if (PartBlocksSky(pl->parts[kv.second.part])) g.solid.insert(kv.first); // (#908) roofs shade the ground
            const Int3 c = kv.first;
            g.radius = std::max(g.radius, sqrtf((float)(c.x * c.x + c.y * c.y + c.z * c.z)) * kCellSize + kCellSize);
        }
        snap->push_back(std::move(g));
    }
    Grass_SetBlockers(static_cast<const DensityField*>(p.field.get()), snap);
    for (const auto& pl : p.plots->plots) {
        if (pl->changedCells.empty()) continue;
        Vector3 lo{ 1e30f, 1e30f, 1e30f }, hi{ -1e30f, -1e30f, -1e30f };
        for (const Int3& c : pl->changedCells) {
            const Vector3 at = pl->toBody({ (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize });
            lo = Vector3Min(lo, at); hi = Vector3Max(hi, at);
        }
        pl->changedCells.clear();
        RequestRemesh(p, Vector3Scale(Vector3Add(lo, hi), 0.5f), 0.5f * Vector3Distance(lo, hi) + kCellSize + 8.0f, 1.0f); // (#908) +8 m: the ground under a roof (its sky)
    }
}

void Planetoid::draw(const RenderContext& ctx) {
    if (!isActive) return;
    if (plots && field && profile && profile->foliage) SyncGrassBlockers(*this); // (#830) grass only grows on worlds with foliage (#1274)
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
        for (auto& pl : plots->plots) if (gOcclusionSkipPlots.empty() || !gOcclusionSkipPlots.count(pl->id)) pl->draw(*this, ctx, *ctx.partShader); // (#134 experiment skip)
        GpuTimer_Mark("parts"); // #134: parts timed apart from the terrain after them
    }
    if (drawList.empty() || !ctx.terrainMaterial) return;
    const ChunkPlacer placer(*this, ctx.renderOrigin);
    ChunkBatch batch(*ctx.terrainMaterial, ctx.terrainLayers, ctx.terrainLayersLoc); // bound once for all chunks
    // Earthlike worlds use the texture-layer palette (dirt, sand, grass, ...); other bodies the scene's (world file `palette`, #1274).
    const int palette = profile && profile->planetPalette ? 2 : ctx.scenePalette;
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
    const double now = GetTime();
    // #243: nodes done fading out leave the list here, each frame, cheaply (re-running the LOD selection every frame while
    // a fade lasted cost walking ~2 ms of CPU); the next LOD pass frees them as undesired.
    if (std::any_of(drawList.begin(), drawList.end(), [&](const Chunk* c) { return c->fadeOut && now - c->fadeT0 >= Chunk::kLodFadeSeconds; }))
        drawList.erase(std::remove_if(drawList.begin(), drawList.end(), [&](Chunk* c) {
            if (!c->fadeOut || now - c->fadeT0 < Chunk::kLodFadeSeconds) return false;
            c->fadeOut = false; c->fadeT0 = -1.0;
            return true;
        }), drawList.end());
    for (Chunk* c : drawList) {
        if (!gOcclusionSkip.empty() && gOcclusionSkip.count(c)) continue; // #134 experiment
        if (c->arena.valid()) { if (c->isActive) TerrainArena_Queue(c->arena, c->clipLocal(), c->lodFade(now)); continue; }
        if (c->fadeOut) continue; // per-chunk draws (RAYLIBGAME_ARENA=0) have no cross-fade: the new node only (#243)
        c->drawChunk(c2, placer.model(c->key));
    }
    ArenaCull cull;
    cull.frustum = ctx.frustum;
    // Hi-Z occlusion (#134): last frame's depth pyramid, tested with the placement this body had in that same frame.
    const ArenaBody bodyNow = placer.arena(ctx.texScale);
    {
        unsigned tex = 0; int w = 0, h = 0, levels = 0; Matrix vp{};
        const long long pyramidFrame = SceneTarget_HiZFrame();
        if (pyramidFrame >= 0 && hizBodyFrame == pyramidFrame && SceneTarget_HiZ(tex, w, h, levels, vp)) {
            cull.hizTex = tex; cull.hizW = w; cull.hizH = h; cull.hizLevels = levels; cull.hizViewProj = vp; cull.hizBody = hizBody;
            // (#910) Not the nodes near the camera: walking moves the eye between frames, and ground hidden from last frame's eye
            // was culled for a frame while visible from this one (grass over nothing). RAYLIBGAME_HIZ_NEAR (m, default 64).
            static const float hizNear = getenv("RAYLIBGAME_HIZ_NEAR") ? (float)atof(getenv("RAYLIBGAME_HIZ_NEAR")) : 64.0f;
            cull.hizCam = toLocal(ctx.cameraPos); cull.hizNear = hizNear;
        }
    }
    const int n = TerrainArena_Flush(ctx.terrainMaterial->shader.id, false, bodyNow, cull);
    if (hizBodyFrame != SceneTarget_Frame()) { hizBody = bodyNow; hizBodyFrame = SceneTarget_Frame(); } // the first draw this frame
    if (ctx.drawCalls) *ctx.drawCalls += n; // listed (the GPU culls some of them)
}

void Planetoid::drawTransparent(const RenderContext& ctx) {
    // Every body with water meshes, not only those with a sea (#670: simulated water on the sandbox plane, or poured on a moon,
    // was simulated but never drawn). Chunks without water add nothing.
    if (!isActive || drawList.empty() || !ctx.waterMaterial || useFarSphere(ctx)) return;
    const ChunkPlacer placer(*this, ctx.renderOrigin);
    ChunkBatch batch(*ctx.waterMaterial);
    const Vector3 c = (worldPos - ctx.renderOrigin).toFloat();
    const float occluder[4] = { c.x, c.y, c.z, std::max(0.0f, radius - 60.0f) * scale }; // night side, as for terrain
    batch.setUniform(ctx.waterOccluderLoc, occluder, SHADER_UNIFORM_VEC4);
    if (ctx.bindSunLighting) ctx.bindSunLighting(ctx.waterMaterial->shader); // cliffs and builds shade the water (slots 6-7)
    Fluid::BindKindUniforms(ctx.waterMaterial->shader); // per-vertex fluid kinds (oil, ...): their colours
    RenderContext c2 = ctx;
    c2.batch = &batch;
    // Nearest first: with the stencil set up by the scene (one water layer per pixel), the closest
    // surface claims the pixel and overlapping LOD copies of the sea behind it are skipped.
    std::vector<std::pair<float, Chunk*>> order;
    for (Chunk* ch : drawList) {
        if (ch->waterMesh.vaoId == 0 || ch->fadeOut) continue; // a node fading out (#243): its replacement draws the water
        const Matrix m = placer.model(ch->key);
        const float h = ch->key.size() * 0.5f;
        order.push_back({ Vector3DistanceSqr(Vector3Transform({ h, h, h }, m), ctx.cameraPos), ch });
    }
    std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    for (auto& o : order) o.second->drawWater(c2, placer.model(o.second->key));
}

void Planetoid::drawFoliage(const RenderContext& ctx) {
    if (!isActive || drawList.empty() || !ctx.grassShader) return;
    const ChunkPlacer placer(*this, ctx.renderOrigin);
    const Vector3 c = (worldPos - ctx.renderOrigin).toFloat();
    const float occluder[4] = { c.x, c.y, c.z, std::max(0.0f, radius - 60.0f) * scale }; // night side, as for terrain
    rlSetUniform(ctx.grassOccluderLoc, occluder, RL_SHADER_UNIFORM_VEC4, 1);
    for (Chunk* ch : drawList) if (ch->key.stride() <= 1.0f && !ch->fadeOut) ch->drawGrass(ctx, placer.model(ch->key)); // the levels that carry grass (foliage.cpp); not a node fading out (#243)
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
    JobSystem::tagNext(JobTag::ChunkPrep);
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
            // (#1058) Counted for cost checks: blocks on the GPU, and those the cave-segment limit (32) sends to the CPU path.
            static std::atomic<int> gCaveGpuBlocks{ 0 }, gCaveCpuBlocks{ 0 };
            static const bool caveCountReport = (atexit([] { printf("bench: cave blocks: %d on the GPU, %d over the 32-segment limit (CPU path)\n", gCaveGpuBlocks.load(), gCaveCpuBlocks.load()); }), true);
            (void)caveCountReport;
            const bool prepOk = pf->prepareBlock(job.block.origin, job.block.stride, job.block.E, job.block.footprint, job.grid) &&
                                pf->blockIslands(job.block.origin, job.block.stride, job.block.E, job.islands, 64);
            const bool caveOk = prepOk && pf->blockCaveSegs(job.block.origin, job.block.stride, job.block.E, job.block.footprint, job.caveSegs, 32);
            if (caveOk) ++gCaveGpuBlocks; else if (prepOk) ++gCaveCpuBlocks;
            if (caveOk) { // (#615) else the CPU path
                const float prepMs = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
                job.keep = field;
                job.field = pf;
                job.wanted = [wanted, version] { return wanted->load() == version; };
                job.priority = priority; // edits (negative) go ahead of background LOD builds
                auto shared = std::make_shared<ChunkMeshData>(std::move(data));
                auto grid = std::make_shared<const PlanetField::BlockGrid>(job.grid); // the materials read the surface from it (#148)
                job.done = [field, edits, shared, wanted, version, priority, prepMs, token, grid](std::shared_ptr<const float> batch, size_t offset) {
                    JobSystem::tagNext(JobTag::ChunkMesh);
                    Jobs().submit(priority, [field, edits, shared, wanted, version, prepMs, batch, offset, token, grid] {
                        if (wanted->load() != version) { ++gSkipped; return; }
                        BuildChunkMeshFromBase(*field, *edits, *shared, batch.get() + offset, prepMs, grid.get(), grid);
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
    float lodDist = LodDistance(v, lo, hi);
    static const bool altOn = !(getenv("RAYLIBGAME_LOD_ALT") && getenv("RAYLIBGAME_LOD_ALT")[0] == '0'); // A/B
    if (altOn && !p.flat && v.zoom <= 1.0f) {
        // Altitude counts (#257): the box reaches far above the thin surface shell inside it, so rising barely moved the
        // box distance. No surface in this node is farther out than the highest terrain (or the node's farthest
        // corner): the camera is at least |cam| - that from it. Near the ground this changes nothing; above the
        // terrain band (flying, rockets) the LOD coarsens with height.
        const float farBox = Vector3Length({ std::max(fabsf(lo.x), fabsf(hi.x)), std::max(fabsf(lo.y), fabsf(hi.y)), std::max(fabsf(lo.z), fabsf(hi.z)) });
        lodDist = std::max(lodDist, Vector3Length(v.cam) - std::min(p.field->maxRadius() + 24.0f, farBox));
    }
    if (key.level > 0 && lodDist < lodSettings.lodFactor * key.size()) {
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
    // LOD cross-fade (#243): a node leaving the list because a finer or coarser one now covers its place stays drawn for
    // Chunk::kLodFadeSeconds, dissolving into its replacement (the pair split the pixels by one screen noise, lighting.fs), so the
    // switch is no longer a pop. Nodes leaving for other reasons (out of range) just go. RAYLIBGAME_LOD_FADE=0: off.
    static const bool lodFadeOn = !(getenv("RAYLIBGAME_LOD_FADE") && getenv("RAYLIBGAME_LOD_FADE")[0] == '0');
    // Not while the camera moves fast (lodSettings.fast, > 15 m/s): flying swaps ~20 nodes a frame, and keeping each for the
    // fade cost ~0.5 ms of CPU a frame (quiet-window A/B) for pops nobody sees at that speed. Nodes dropped then lose any
    // fade-out mark, so one coming back later isn't drawn half faded.
    if (lodFadeOn && lodSettings.fast) {
        for (Chunk* c : p.drawList) if (c->fadeOut) { c->fadeOut = false; c->fadeT0 = -1.0; }
    } else if (lodFadeOn && !p.drawList.empty()) {
        const double now = GetTime();
        static FlatKeySet newSet, newWithin, oldSet, goneSet, goneWithin; // new list (+ ancestors), old list, replaced (+ ancestors)
        newSet.reset(draw.size()); newWithin.reset(draw.size());
        for (Chunk* c : draw) { newSet.insert(c->key); for (NodeKey k = parentOf(c->key); k.level <= p.rootLevel; k = parentOf(k)) if (!newWithin.insert(k)) break; }
        oldSet.reset(p.drawList.size()); goneSet.reset(64); goneWithin.reset(64);
        std::vector<Chunk*> keep;
        for (Chunk* c : p.drawList) {
            oldSet.insert(c->key);
            if (newSet.count(c->key)) { if (c->fadeOut) { c->fadeOut = false; c->fadeT0 = -1.0; } continue; } // still (or again) in the list
            if (c->fadeOut) { // fading out already
                if (now - c->fadeT0 < Chunk::kLodFadeSeconds) keep.push_back(c); else { c->fadeOut = false; c->fadeT0 = -1.0; }
                continue;
            }
            bool replaced = newWithin.count(c->key) > 0; // finer nodes now inside it
            for (NodeKey k = parentOf(c->key); !replaced && k.level <= p.rootLevel; k = parentOf(k)) replaced = newSet.count(k) > 0; // a coarser one over it
            if (!replaced) continue;
            c->fadeOut = true; c->fadeT0 = now; keep.push_back(c); ++gLodFades;
            goneSet.insert(c->key);
            for (NodeKey k = parentOf(c->key); k.level <= p.rootLevel; k = parentOf(k)) if (!goneWithin.insert(k)) break;
        }
        if (!goneSet.empty())
            for (Chunk* c : draw) { // the replacements: new in the list where a node went away
                if (oldSet.count(c->key)) continue;
                bool overlaps = goneWithin.count(c->key) > 0;
                for (NodeKey k = parentOf(c->key); !overlaps && k.level <= p.rootLevel; k = parentOf(k)) overlaps = goneSet.count(k) > 0;
                if (overlaps) { c->fadeT0 = now; c->fadeOut = false; }
            }
        for (Chunk* c : keep) draw.push_back(c);
        // (Faded-out nodes leave the list in Planetoid::draw, not by re-running this every frame: #243 walk cost.)
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
    // Also while the terrain arena is over 85% full: the cache must never make it regrow (a GB-sized copy: a 30 ms frame
    // closing the telescope, whose zoomed nodes all landed in the cache at once). The old 4 s hold had this check.
    const ArenaStats as = TerrainArena_Stats();
    long long arenaUsed = as.vertices;
    const long long arenaLimit = (long long)as.vertexCapacity * 85 / 100;
    auto over = [&] { return cachedVerts > lodSettings.cacheVertices || cachedCount > (size_t)lodSettings.cacheNodes || (as.vertexCapacity > 0 && arenaUsed > arenaLimit); };
    if (!cached.empty() && over()) {
        std::sort(cached.begin(), cached.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        size_t i = 0;
        for (; i < cached.size() && over(); ++i) {
            if (frees >= kMaxFrees) { p.lodDirty = true; p.lastFreeSweep = -1e9; break; }
            auto it = p.nodes.find(cached[i].second);
            cachedVerts -= vertsOf(*it->second);
            arenaUsed -= vertsOf(*it->second);
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
        if (ring.lastCam.x < 1e29f && dt > 0.0 && dt < 1.0 && Vector3DistanceSqr(c, ring.lastCam) < 50.0f * 50.0f) { // a bigger step is a teleport
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
        rp.lodFactor = lodSettings.detailFactor; // the walking detail: the ring is paused when moving fast anyway (#297)
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
        JobSystem::tagNext(JobTag::RingCollect);
        Jobs().submit(0.0f, /* background (#296): not counted by WorldIdle */ [rp, job, roots] {
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
        }, true);
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
        // Room for it in the terrain arena now, while the ring is still mostly unbuilt: growing later, full, is a
        // GB-sized copy (a 30 ms frame closing the telescope at Medium detail). ~770 vertices a ring node (seed 42,
        // Medium), plus a full cache and the view: a re-centred ring's old part goes to the cache.
        const ArenaStats as = TerrainArena_Stats();
        // The estimate is generous (26 M of 44.6 M used on seed 42): regrow only when it passes the capacity by more than
        // 10%, with 10% headroom. A re-centred ring 2 nodes larger than the spawn ring regrew it whole: a 57 ms frame (#364).
        long long need = (long long)ring.keys.size() * 800 + lodSettings.cacheVertices + (4 << 20);
        if (need > as.vertexCapacity + as.vertexCapacity / 10) { // (Reserve creates the arena if it is not up yet)
            need += need / 10;
            const double ipv = std::max(6.0, as.vertices > 0 ? (double)as.indices / as.vertices : 6.0);
            TerrainArena_Reserve((int)std::min<long long>(need, INT_MAX / 2), (int)std::min<long long>((long long)(need * ipv), INT_MAX / 2));
        }
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
    if (r.overflow) LOG_ERROR("World", "Node L%d (%d,%d,%d) exceeds 65535 vertices; skipped\n", r.key.level, r.key.x, r.key.y, r.key.z);
    const bool wasReady = it->second->ready;
    it->second->setMesh(r);
    if (!wasReady && it->second->inRing)
        if (auto ring = gPreload.find(p->id); ring != gPreload.end()) ++ring->second.ready; // (#294 progress)
    // Readiness changed: the draw list may change. Not for a preload-ring build (#296): the view did not ask for
    // it (it re-requests a ring node it wants), and re-selecting after every background upload cost frame time.
    if (!it->second->background) p->lodDirty = true;
}

// Idle for the view: no view or edit build alive, no other foreground job (tree tiles, water classification, the
// far-sphere bake, minimap rows: #296 review), and none of their results waiting for upload. The preload ring's
// background builds do not count (#296): loading screens and scripted tests would wait seconds for the ring.
static bool PipelineIdle() {
    if (gForegroundBuilds.load() != 0) return false; // a result is pushed before its build's token is released
    if (!Jobs().foregroundIdle()) return false;
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

namespace { double worldTime = 0.0, followOffset = 0.0, followLag = 0.0; }
double WorldTime() { return worldTime; }
void SetWorldTime(double t) { worldTime = t; followOffset = 0.0; }
// #339: a client follows the server's world clock. Offsets under 2 s are slewed in at 10% of the frame time (spin, orbits,
// the sun and weather never jump); larger ones (a join, a long stall) are taken at once.
void WorldTime_Follow(double serverNow) { followLag = serverNow - worldTime; followOffset = followLag; if (fabs(followOffset) > 2.0) { worldTime = serverNow; followOffset = 0.0; } }
double WorldTime_Lag() { return followLag; }
double WorldTime_Pending() { return followOffset; }

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
    if (followOffset != 0.0 && dt > 0.0f) { // #339: slewing toward the server's clock
        const double step = std::clamp(followOffset, -0.1 * dt, 0.1 * dt);
        worldTime += step; followOffset -= step;
    }
    bool moved = false;
    for (auto& obj : world.objects) {
        auto* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p) continue;
        if (p->spin.x != 0.0f || p->spin.y != 0.0f || p->spin.z != 0.0f) {
            float* r[3] = { &p->rotation.x, &p->rotation.y, &p->rotation.z };
            const float s[3] = { p->spin.x, p->spin.y, p->spin.z };
            const bool spinChanged = p->spin.x != p->spinPhaseFor.x || p->spin.y != p->spinPhaseFor.y || p->spin.z != p->spinPhaseFor.z;
            for (int k = 0; k < 3; ++k) {
                // The angle is a function of the world clock (#315): phase + spin * time. Something else (a load, a script)
                // set the angle, or the spin changed: take the angle over as the phase at the time before this step.
                if (!p->spinAngleSet || spinChanged || (float)p->spinAngle[k] != *r[k])
                    p->spinPhase[k] = (double)*r[k] - (double)s[k] * (worldTime - (double)dt);
                double a = fmod(p->spinPhase[k] + (double)s[k] * worldTime, 360.0);
                if (a < 0.0) a += 360.0;
                p->spinAngle[k] = a;
                *r[k] = (float)a;
            }
            p->spinPhaseFor = p->spin;
            p->spinAngleSet = true;
        }
        if (p->orbitParent) { PlaceBodyOnOrbit(*p); moved = true; } // parents come first in the list
    }
    if (moved) world.syncRenderPositions();
}

// ---------------------------------------------------------------------------------------------
// Water redraws (#360): the chunk mesher draws every fluid surface from the one water field (waterfield.h). The
// simulation publishes the boxes whose cells changed; the nodes there (every level, found by key arithmetic) are
// marked, and the drawn ones redraw their water alone - at most every 0.1 s near (0.2 / 0.5 / 1 s coarser) - from the
// terrain density they kept (ChunkWaterCache). A node without that cache gets one full rebuild instead.
// ---------------------------------------------------------------------------------------------

static void MarkWaterDirty(Planetoid& p, Vector3 lo, Vector3 hi, int mag, double now) {
    for (int level = 0; level <= p.rootLevel; ++level) {
        const float s = kVoxelSize * (float)(CHUNK_SIZE << level);
        // Nodes whose samples ([origin, origin + (CHUNK_SIZE + 1) * stride]) reach the box: exactly what the chunk's
        // ChunkSampler tests, so a marked node finds the simulated water (#382 review: a wider margin marked coarse nodes
        // that found none, and those got a full rebuild every time).
        const float margin = s / CHUNK_SIZE;
        const int x0 = (int)floorf((lo.x - margin) / s), x1 = (int)floorf(hi.x / s);
        const int y0 = (int)floorf((lo.y - margin) / s), y1 = (int)floorf(hi.y / s);
        const int z0 = (int)floorf((lo.z - margin) / s), z1 = (int)floorf(hi.z / s);
        if ((long long)(x1 - x0 + 1) * (y1 - y0 + 1) * (z1 - z0 + 1) > 4096) continue;
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) {
                    auto it = p.nodes.find(NodeKey{ level, x, y, z });
                    if (it == p.nodes.end()) continue;
                    Chunk& c = *it->second;
                    if (!c.waterDirty) { c.waterDirty = true; c.waterDirtySince = now; c.waterMag = 0; }
                    // Summed, not the largest (#469): each publish is a cell's change of 6+ units since it was last published,
                    // so a slow rise came as 6-7 at a time and never reached a coarse node's gate (8+): it waited the 3 s
                    // fallback. Jitter under 6 never publishes, so the sum doesn't add it up.
                    c.waterMag = std::min(255, c.waterMag + mag);
                }
    }
}

// Water jobs queue by distance like view builds (#382 review: at -0.5 they went ahead of every view build).
static void RequestWater(Planetoid& p, Chunk& c, float priority) {
    c.waterBusy = true;
    c.waterDirty = false;
    c.waterMag = 0;
    std::shared_ptr<const DensityField> field = p.field;
    std::shared_ptr<const ChunkWaterCache> cache = c.waterCache;
    const uint64_t id = p.id;
    const NodeKey key = c.key;
    auto token = std::make_shared<WaterToken>();
    JobSystem::tagNext(JobTag::WaterRebuild);
    Jobs().submit(priority, [field, cache, id, key, token] {
        ChunkMeshData d;
        d.planetoidId = id;
        d.key = key;
        d.waterCache = cache; // which terrain it was cut by: dropped if the chunk has newer terrain by now
        BuildChunkWater(*field, *cache, d);
        ++gWaterBuilds;
        gWaterUs += (long long)(d.buildMs * 1000.0f);
        std::lock_guard<std::mutex> lock(resultsMutex);
        waterResults.push_back(std::move(d));
    });
}

static void UpdateWater(Scene& world) {
    const double now = GetTime();
    std::vector<WaterField::Dirty> dirty;
    WaterField::TakeDirty(dirty);
    for (const WaterField::Dirty& d : dirty)
        if (Planetoid* p = Planetoid::find(d.body)) MarkWaterDirty(*p, d.lo, d.hi, d.mag, now);
    std::vector<ChunkMeshData> done;
    { std::lock_guard<std::mutex> lock(resultsMutex); done.swap(waterResults); }
    for (ChunkMeshData& r : done) {
        Planetoid* p = Planetoid::find(r.planetoidId);
        if (!p) continue;
        auto it = p->nodes.find(r.key);
        if (it == p->nodes.end()) continue;
        Chunk& c = *it->second;
        c.waterBusy = false;
        c.waterBuiltAt = now;
        // Cut by terrain the chunk no longer has (a full build after an edit landed meanwhile): drawn again (#382 review).
        if (r.waterCache != c.waterCache) { c.waterDirty = true; c.waterMag = 255; c.waterDirtySince = now; continue; }
        if (r.waterSeq >= c.waterSeq) c.uploadWater(r);
    }
    // Water caches only for nodes on screen (#382 review: ring and cached nodes kept ~160 KB each, past 128 MB): a node
    // out of the draw list for 5 s drops it; drawn again, its next water change costs one full build.
    static double lastSweep = 0.0;
    if (now - lastSweep > 2.0) {
        lastSweep = now;
        for (auto& obj : world.objects)
            if (auto* p = dynamic_cast<Planetoid*>(obj.get()))
                for (auto& n : p->nodes)
                    if (n.second->waterCache && !n.second->waterBusy && now - n.second->waterDrawnAt > 5.0) n.second->waterCache.reset();
    }
    int budget = 64; // water jobs a frame
    // In flight at once: half the workers, so streaming always has the rest (a draining pond kept hundreds a second going).
    const int maxInFlight = std::max(1, (int)Jobs().threadCount() / 2);
    for (auto& obj : world.objects) {
        auto* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p || !p->isActive) continue;
        const Vector3 cam = p->toLocal(world.camera.position);
        for (Chunk* c : p->drawList) {
            c->waterDrawnAt = now;
            if (!c->waterDirty || c->waterBusy || !c->ready || budget <= 0 || gWaterInFlight.load() >= maxInFlight) continue;
            const Vector3 lo = c->key.origin();
            const float dist = DistanceToBox(cam, lo, Vector3AddValue(lo, c->key.size()));
            // Redraw rate by level and by distance (#381): 10 a second within 10 m, down to once a second from 100 m.
            const double byLevel = c->key.level == 0 ? 0.1 : c->key.level <= 2 ? 0.2 : c->key.level <= 4 ? 0.5 : 1.0;
            const double every = std::max(byLevel, std::clamp((double)dist / 100.0, 0.1, 1.0));
            if (now - c->waterBuiltAt < every) continue;
            // Coarser nodes wait for a change they can show (their samples are 0.5-8 m apart: a 2.5 cm step is invisible
            // there), or 3 s (#381).
            static const int kShowMag[5] = { 0, 8, 13, 26, 51 };
            if (c->waterMag < kShowMag[std::min(c->key.level, 4)] && now - c->waterDirtySince < 3.0) continue;
            --budget;
            const float priority = dist / (1.0f + c->key.level);
            if (c->waterCache) { RequestWater(*p, *c, priority); continue; }
            c->waterDirty = false; c->waterBuiltAt = now;
            // No cache: a full build keeps one - unless no simulated water reaches the node (then nothing to redraw).
            const WaterField::SnapshotPtr snap = WaterField::Current();
            const WaterField::ChunkSampler probe(snap.get(), p->id, lo, Vector3AddValue(lo, (CHUNK_SIZE + 1) * c->key.stride()), true);
            if (!probe.any()) continue;
            ++gWaterFullBuilds;
            RequestNode(*p, c->key, priority);
        }
    }
}

void worldUpdate(Scene& world, float uploadBudgetMs) {
    PROFILE_FUNCTION();
    GpuTimer_CpuMark("w-before");
    TerrainGpu_SetBackground(kBackgroundPriority, !gPreloadUrgent); // preload-ring blocks: gently, except during the spawn preload
    {   // (#1122) New GPU density batches by the frame time: both slots while frames are quick, one past 7 ms, one every other
        // frame past 12 ms. Fast flight high up (300 m, 250 m/s) asked for a batch in every free slot every frame and the
        // density compute plus its uploads ate the frame (66 FPS, still falling); the terrain's detail now lags a little
        // instead. Not during the spawn preload or a load. RAYLIBGAME_STREAM_BUDGET=0: as before (A/B).
        static const bool budgetOn = !(getenv("RAYLIBGAME_STREAM_BUDGET") && getenv("RAYLIBGAME_STREAM_BUDGET")[0] == '0');
        static long frameNo = 0;
        ++frameNo;
        const float frameMs = GetFrameTime() * 1000.0f;
        int limit = 2;
        if (budgetOn && !gPreloadUrgent && generationStart < 0.0) limit = frameMs > 12.0f ? (frameNo % 2 == 0 ? 1 : 0) : frameMs > 7.0f ? 1 : 2;
        TerrainGpu_SetLaunchLimit(limit);
    }
    TerrainGpu_Pump(); // GPU density blocks (Plan 6 step 4): finished batches back to the workers, the next batch out
    GpuTimer_CpuMark("w-gpu-pump");
    // Ring rocks near the camera and near every interest centre (every player, #318) become voxel bodies (and back)
    // before this frame's LOD selection.
    if (!Rings().empty()) {
        std::vector<DVec3> foci{ world.renderOrigin + world.camera.position };
        for (const InterestCentre& c : gInterestCentres)
            if (const Planetoid* b = Planetoid::find(c.body)) foci.push_back(b->toWorldD(DVec3(c.local)));
        for (auto& ring : Rings()) ring->update(world, foci);
    }

    {
        std::lock_guard<std::mutex> lock(resultsMutex);
        for (ChunkMeshData& r : results) pendingUploads.push_back(std::move(r));
        results.clear();
    }

    // Uploads nearest first (#258): edits, then by the distance from the camera now (the queue was FIFO, so near
    // ground waited behind far nodes requested earlier). Nodes nobody wants any more drop out in IntegrateResult.
    // #387: only the few that fit this frame are picked (an O(n) rank + nth_element), not the whole backlog sorted and
    // moved every frame (15-26 ms a frame at 12-18k pending, which slowed the frames that drain it: a spiral).
    static const float kBudgetOverride = getenv("RAYLIBGAME_UPLOAD_BUDGET_MS") ? (float)atof(getenv("RAYLIBGAME_UPLOAD_BUDGET_MS")) : -1.0f; // tests: force a backlog
    if (kBudgetOverride >= 0.0f) uploadBudgetMs = kBudgetOverride;
    const size_t kPick = gPreloadUrgent ? 1024 : 64; // at most this many uploads a frame (the spawn preload takes more)
    std::vector<size_t> pick;
    if (!pendingUploads.empty()) {
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
        for (size_t i = 0; i < pendingUploads.size(); ++i) order.push_back({ gStreamOld ? (float)i : rank(pendingUploads[i]), i });
        const size_t k = std::min(kPick, order.size());
        auto byRank = [](const auto& a, const auto& b) { return a.first < b.first || (a.first == b.first && a.second < b.second); };
        if (k < order.size()) std::nth_element(order.begin(), order.begin() + (k - 1), order.end(), byRank);
        std::sort(order.begin(), order.begin() + k, byRank);
        for (size_t i = 0; i < k; ++i) pick.push_back(order[i].second);
    }
    // GPU uploads are the only main-thread cost; cap them per frame so streaming never hitches.
    const auto t0 = std::chrono::steady_clock::now();
    int uploads = 0;
    // Preload-ring results (ranked last) get at most 1 ms of it outside the spawn preload: spare time only (#296).
    const float backgroundBudgetMs = gPreloadUrgent ? uploadBudgetMs : std::min(uploadBudgetMs, 1.0f);
    std::vector<size_t> done;
    for (size_t idx : pick) {
        const float elapsed = std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count();
        if (uploads > 0 && pendingUploads[idx].priority >= kBackgroundPriority && elapsed > backgroundBudgetMs) break;
        const auto u0 = std::chrono::steady_clock::now();
        IntegrateResult(pendingUploads[idx]);
        gUploadMaxMs = std::max(gUploadMaxMs, std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - u0).count());
        done.push_back(idx);
        uploads++;
        if (std::chrono::duration<float, std::milli>(std::chrono::steady_clock::now() - t0).count() > uploadBudgetMs) break;
    }
    // Out of the backlog: swap-remove, highest index first so the indices still to remove stay valid.
    std::sort(done.begin(), done.end(), std::greater<size_t>());
    for (size_t idx : done) {
        if (idx + 1 != pendingUploads.size()) pendingUploads[idx] = std::move(pendingUploads.back());
        pendingUploads.pop_back();
    }
    // Back-pressure (#387): while the upload backlog is long, no new GPU density batches (finished ones are still
    // collected), so results can't pile up faster than frames take them.
    TerrainGpu_SetHold(pendingUploads.size() > kUploadBacklogCap);

    const auto w0 = std::chrono::steady_clock::now();
    UpdateWater(world); // simulated water: chunks redraw their water surface (#360)
    gUploadLoopMs += std::chrono::duration<double, std::milli>(w0 - t0).count();
    gWaterUpdMs = std::max(gWaterUpdMs, std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - w0).count());
    GpuTimer_CpuMark("w-uploads");
    // LOD after uploads, so this frame's draw list already includes nodes that just finished.
    // Terrain detail (#297): the player's setting, at most fastFactor while the camera moves fast against the nearest
    // body (hysteresis between slowSpeed and fastSpeed). A change re-selects every body; a changed setting also
    // re-collects the preload rings (they are built at the setting's detail).
    {
        LodSettings& s = lodSettings;
        static Vector3 lastLocal{}; static uint64_t lastBody = 0; static double lastT = 0.0; static float speed = 0.0f;
        static float lastDetail = -1.0f;
        const double now = GetTime();
        if (Planetoid* nb = NearestPlanetoid(world, world.camera.position)) {
            const Vector3 l = nb->toLocal(world.camera.position);
            const double dt = now - lastT;
            if (nb->id == lastBody && dt > 0.0 && dt < 1.0 && Vector3DistanceSqr(l, lastLocal) < 50.0f * 50.0f) // a bigger step is a teleport
                speed += (Vector3Distance(l, lastLocal) / (float)dt - speed) * (float)std::min(1.0, dt * 4.0);
            lastLocal = l; lastBody = nb->id;
        }
        lastT = now;
        if (s.fast ? speed < s.slowSpeed : speed > s.fastSpeed) s.fast = !s.fast;
        const float factor = gStreamOld ? s.detailFactor : s.fast ? std::min(s.detailFactor, s.fastFactor) : s.detailFactor;
        if (factor != s.lodFactor) {
            s.lodFactor = factor;
            for (auto& obj : world.objects) if (auto* p = dynamic_cast<Planetoid*>(obj.get())) p->lodDirty = true;
        }
        if (s.detailFactor != lastDetail) {
            if (lastDetail >= 0.0f) for (auto& [id, ring] : gPreload) { ring.centre = { 1e30f, 1e30f, 1e30f }; ring.job.reset(); }
            lastDetail = s.detailFactor;
        }
    }
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
        stats.waterBuilds = gWaterBuilds.load(); stats.waterMs = gWaterUs.load() / 1000.0;
        stats.waterFullBuilds = gWaterFullBuilds.load(); stats.waterCacheBytes = ChunkWaterCacheBytes();
        if (gPreloadUrgent && complete) {
            gPreloadUrgent = false;
            // Still behind the loading screen: grow the terrain arena once to the real ring plus a full cache and the
            // view (the estimate at the ring install may be short), so it never regrows (a GB-sized copy) mid-walk.
            const ArenaStats as = TerrainArena_Stats();
            if (as.vertices > 0) {
                const double ipv = (double)as.indices / as.vertices;
                const long long v = (long long)as.vertices + lodSettings.cacheVertices + (4 << 20);
                TerrainArena_Reserve((int)std::min<long long>(v, INT_MAX / 2), (int)std::min<long long>((long long)(v * ipv), INT_MAX / 2));
            }
        }
    }
    GpuTimer_CpuMark("w-preload");
    for (auto& obj : world.objects)
        if (auto* p = dynamic_cast<Planetoid*>(obj.get()); p && p->trees) {
            std::vector<Vector3> extra; // other interest centres on this body (#318)
            for (const InterestCentre& c : gInterestCentres) if (c.body == p->id) extra.push_back(c.local);
            p->trees->setInterest(std::move(extra));
            p->trees->update(*p, p->toLocal(world.camera.position));
        }

    GpuTimer_CpuMark("w-trees");
    stats.queuedJobs = (int)Jobs().queued();
    stats.runningJobs = Jobs().running();
    stats.pendingUploads = (int)pendingUploads.size();
    stats.uploadsLastFrame = uploads;
    stats.workerThreads = (int)Jobs().threadCount();
    stats.requests = gRequests; stats.built = gBuilt; stats.skippedStale = gSkipped; stats.droppedResults = gDropped; stats.lodFades = gLodFades;
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
    if (gJobLog) {
        // #387: every second, the job queue and what was submitted / finished per kind, to log.txt and stdout.
        static double t0 = -1.0;
        static long long sub0[(int)JobTag::Count] = {}, done0[(int)JobTag::Count] = {};
        static long long req0 = 0, built0 = 0;
        static long long staged0 = 0, direct0 = 0;
        static double ftSum = 0.0, ftMax = 0.0; static int ftN = 0;
        ftSum += GetFrameTime(); ftMax = std::max(ftMax, (double)GetFrameTime()); ++ftN;
        const double now = GetTime();
        if (t0 < 0.0) t0 = now;
        if (now - t0 >= 1.0) {
            std::string mix;
            for (int k = 0; k < (int)JobTag::Count; ++k) {
                const JobTag tg = (JobTag)k;
                const long long s = Jobs().submittedOf(tg), d = Jobs().doneOf(tg), q = Jobs().queuedOf(tg);
                if (s - sub0[k] || d - done0[k] || q)
                    mix += TextFormat(" %s q%lld +%lld -%lld", JobTagName(tg), q, s - sub0[k], d - done0[k]);
                sub0[k] = s; done0[k] = d;
            }
            const ArenaStats as = TerrainArena_Stats();
            const char* line = TextFormat("[jobs] %.0f s: frame avg %.1f max %.1f ms; queued %d running %d fg %s | gpu pending %d, uploads pending %d%s (staged +%lld direct +%lld, slowest %.1f ms, loop %.0f ms, water max %.1f ms) | requests +%lld built +%lld, nodes %d |%s\n",
                                          now, 1000.0 * ftSum / std::max(1, ftN), 1000.0 * ftMax, (int)Jobs().queued(), Jobs().running(), Jobs().foregroundIdle() ? "idle" : "busy", TerrainGpu_Pending(),
                                          (int)pendingUploads.size(), pendingUploads.size() > kUploadBacklogCap ? " (hold)" : "", as.staged - staged0, as.direct - direct0, gUploadMaxMs, gUploadLoopMs, gWaterUpdMs, (long long)gRequests - req0, (long long)gBuilt - built0, stats.nodes, mix.c_str());
            printf("%s", line);
            LOG_INFO("Perf", "%s", line);
            req0 = gRequests; built0 = gBuilt; staged0 = as.staged; direct0 = as.direct; gUploadMaxMs = gUploadLoopMs = gWaterUpdMs = 0.0; ftSum = ftMax = 0.0; ftN = 0;
            t0 = now;
        }
    }

    if (generationStart >= 0.0 && PipelineIdle() && !AnyLodPending(world)) {
        stats.lastGenerationSeconds = (float)(GetTime() - generationStart);
        LOG_INFO("World", "ready in %.2f s (%d meshes drawn, %d nodes, avg build %.2f ms on %d threads)\n",
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
                    const Vector3 nLocal = p->field->normal(local, p->edits.get());
                    hit.normal = Vector3Normalize(Vector3Subtract(p->toWorld(Vector3Add(local, nLocal)), worldPoint));
                    hit.material = p->materialAt(local, nLocal); // as drawn (#636)
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
        const float margin = rmax + (k.level > 0 ? 3.0f : 1.0f) * k.stride() + 1.0f; // coarse nodes read edits 2 strides along their columns (#679)
        if (hi.x + margin < nlo.x || lo.x - margin > nhi.x || hi.y + margin < nlo.y || lo.y - margin > nhi.y || hi.z + margin < nlo.z || lo.z - margin > nhi.z) continue;
        if (count == 1 && DistanceToBox(edits[0].center, nlo, nhi) > margin) continue;
        touched.push_back(k);
    }
    for (const NodeKey& k : touched) RequestNode(p, k, -1.0f - (float)k.level * 0.01f); // ahead of background work
    if (p.trees)
        for (size_t i = 0; i < count; i += std::max<size_t>(1, count / 64)) p.trees->onEdit(edits[i].center, edits[i].radius, p.edits);
    // Flowing water: wake it, or start it where the edit opens space next to water.
    for (size_t i = 0; i < count; i += std::max<size_t>(1, count / 64))
        Fluid::OnTerrainEdit(p, edits[i].center, edits[i].radius, !edits[i].carve && edits[i].strength > 0.0f); // adding ground: may dam a river (#361)
}

void RequestRemesh(Planetoid& p, Vector3 c, float r, float strideReach) {
    std::vector<NodeKey> touched;
    for (auto& entry : p.nodes) {
        const NodeKey& k = entry.first;
        const Vector3 nlo = k.origin();
        const Vector3 nhi = Vector3AddValue(nlo, k.size() * (CHUNK_SIZE + 1.0f) / CHUNK_SIZE);
        if (DistanceToBox(c, nlo, nhi) > r + strideReach * k.stride() + 1.0f) continue;
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
    loadWorld(world);
    // No bodies yet (#503): New Game generates them (the erosion bake behind its screen), Load Game rebuilds them from the
    // save, a multiplayer client on its first join. Benches generate at start (main.cpp).
}

bool ReloadSystemDef() {
    WorldProfile_ClearCache(); // the world files are read again with the bodies file (#1274)
    SystemDef sys;
    const bool ok = LoadSystemDef(worldGen.bodiesPath, sys);
    if (!ok) {
        if (worldGen.systemLoaded) {
            LOG_INFO("Bodies", "keeping the current system\n");
            return false;
        }
        LOG_WARN("Bodies", "falling back to the built-in start system\n");
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

std::function<void(float, const char*)> gWorldGenWait;

void generateWorld(Scene& world) {
    PROFILE_FUNCTION();
    if (!worldGen.systemLoaded) ReloadSystemDef();
    gRings.clear();
    NetId_Reset(); // a new world: allocated ids start again (a load sets the counter from its save, #317)
    worldGen.lastSeed = worldGen.seed ? worldGen.seed : rd();
    gen.seed(worldGen.lastSeed);
    GameRng_Seed(worldGen.lastSeed); // gameplay rolls follow the seed (a load restores the saved stream, #316)
    LOG_INFO("World", "Generating world (seed %u, %zu definitions from %s)...\n", worldGen.lastSeed, worldGen.system.bodies.size(), worldGen.bodiesPath.c_str());

    // Every body draws the same values in the same order (3 position values, its terrain seed,
    // a size value for ring members, a spin value), so a seed reproduces the system exactly.
    std::uniform_real_distribution<double> unit(0.0, 1.0);
    std::unordered_map<std::string, Planetoid*> byName;
    int index = 0;
    auto add = [&](const BodyDef& d, BodyKind kind, const std::string& name, DVec3 pos, float radius, unsigned int seed, float spin) {
        auto p = std::make_unique<Planetoid>(name, pos, Vector3{ 0, 0, 0 }, 1.0f, seed, radius);
        p->kind = kind;
        p->profile = WorldProfile_For(name, kind); // (#1274) planets and moons: resources/data/worlds/<name>.txt over the kind's default
        if (p->profile) {
            const WorldProfile& wp = *p->profile;
            // Planet terrain: radius is sea level, where gravity's reference surface sits too; the
            // LOD root must also cover the highest mountains.
            auto field = std::make_shared<PlanetField>(seed, radius, wp.airless);
            const bool sandbox = gSandboxWorld && kind == BodyKind::Planet; // #662: a flat test plane (no bake, no trees)
            if (wp.newWorld) field->setGenParams(WorldGenParams_ForSave(WorldSavePath())); // #504: the New World screen's choice (defaults without one)
            if (sandbox) field->setSandbox();
            if (gErosionPrototype && wp.erosion && !wp.airless && field->gen.erosion) {
                // Plan 7 (#355): the save's bake, or a new one (~60 s). On a thread (#362, the #355 review): the main
                // thread keeps drawing the loading screen through gWorldGenWait instead of freezing the window. The field
                // is only shared (p->field below) once the bake is done.
                const std::string cache = std::string(WorldSavePath()) + ".erosion";
                std::atomic<bool> baked{ false };
                std::thread bake([&] { field->bakeErosion(cache); baked.store(true); });
                while (!baked.load()) {
                    // A cache read (a save's own bake) is "loading"; only a real bake says it shapes the world (#505).
                    if (gWorldGenWait) gWorldGenWait(ErosionBaking() ? ErosionBakeProgress() : 0.0f, ErosionBaking() ? "Shaping the world (erosion)" : "Loading the world's terrain");
                    else std::this_thread::sleep_for(std::chrono::milliseconds(20));
                }
                bake.join();
            }
            p->gravityRadius = radius;
            p->rootLevel = 0;
            while (kVoxelSize * (float)(CHUNK_SIZE << p->rootLevel) < field->maxRadius() + 24.0f) p->rootLevel++;
            p->field = field;
            if (wp.foliage && !sandbox) {
                field->foliage = true;
                field->foliageSeed = seed ^ 0x7aee5u;
                p->trees = std::make_unique<TreeSystem>(field, field->foliageSeed);
            }
            p->farSphere = std::make_unique<FarSphere>(field, field->foliage, field->foliageSeed); // bakes on the workers
        }
        p->spin = { 0.0f, spin, 0.0f };
        p->surfaceGravity = d.gravity >= 0.0f ? d.gravity : AutoSurfaceGravity(kind, radius);
        p->parent = &world.rootObject;
        p->netId = NetId_Make(NetKind::Body, (uint64_t)index); // the creation order: the same for every run of the system (#317)
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
                if (n > 0) LOG_INFO("World", "%s clears %d rocks from ring %s\n", p->name.c_str(), n, ring->name.c_str());
            }
    world.syncRenderPositions();
    generationStart = GetTime();
    buildMsTotal = sampleMsTotal = 0.0;
    buildCount = 0;
}

void BeginWorldStreaming() { generationStart = GetTime(); buildMsTotal = sampleMsTotal = 0.0; buildCount = 0; }

// #387: everything queued for the old world's streaming: the GPU density queue (its blocks stayed wanted and kept the old
// world and its bake alive), finished results and the upload backlog (18k stayed pending at the menu after a quit).
void World_ClearStreamingQueues() {
    TerrainGpu_ClearQueued();
    { std::lock_guard<std::mutex> lock(resultsMutex); results.clear(); waterResults.clear(); }
    pendingUploads.clear();
    TerrainGpu_SetHold(false);
}

void unloadWorld(Scene& world) {
    for (auto& h : UnloadHooks()) h(); // (#902) pointers into the old world dropped before it goes
    Jobs().clearQueued();  // in-flight jobs finish; their results are dropped (planetoid ids are gone)
    World_ClearStreamingQueues(); // #387 (first: nothing queued to wait for; again below for what lands meanwhile)
    // The builds already running on the workers (CPU only: a GPU density block is a callback on the main thread, never a
    // blocked worker) take a few ms; give them up to 300 ms of wall clock, then go on without them (their results are
    // dropped). Never a count of short sleeps: a 1 ms sleep lasts ~15 ms on Windows, and the old 3000-step cap became a
    // 23.6 s frame (#270).
    const auto waitUntil = std::chrono::steady_clock::now() + std::chrono::milliseconds(300);
    while (Jobs().running() > 0 && std::chrono::steady_clock::now() < waitUntil) std::this_thread::yield();
    Factory::Reset();      // machines belong to the old bodies
    Fluid::Reset();        // so does flowing water
    player.nearest = nullptr; // every body goes (the player refinds its ground next update)
    gRings.clear();        // before their parent bodies go
    world.objects.clear(); // releases chunk meshes (and their arena ranges)
    World_ClearStreamingQueues(); // #387: queued GPU blocks, results and uploads of the old world
    gPreload.clear();      // the preload rings and cache counts of the old bodies
    gCachedCount.clear();
    // The GPU side of the old world (#270): the terrain arena (grown to the last world's size: hundreds of MB) and the GPU
    // density batches, abandoned, not waited for (they point at the old field; their builds' bodies are gone). Both come
    // back on the next world's first use.
    TerrainArena_Unload();
    TerrainGpu_Unload();
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
bool gSandboxWorld = false;

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
    thread_local std::vector<int16_t>* denseP = nullptr; if (!denseP) denseP = new std::vector<int16_t>(DeltaChunk::kN); auto& dense = *denseP; // (#1133, #485: MinGW frees thread_local storage before running its destructor at thread exit: never destroyed)
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
    thread_local std::vector<int16_t>* denseP = nullptr; if (!denseP) denseP = new std::vector<int16_t>(DeltaChunk::kN); auto& dense = *denseP; // (#1133, #485: MinGW frees thread_local storage before running its destructor at thread exit: never destroyed)
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
        for (int k = 0; k < count; ++k, ++w) c.setMat(i++, (uint8_t)(std::min<int>(in[w] & 0x7F, (int)VoxelMaterial::Count) | (in[w] & kMatState))); // (#1097) bit 7 kept
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
    if (ok && memcmp(magic, "RGN5", 4) != 0) { LOG_WARN("World", "%s: terrain edits on the old 1 m lattice, skipped\n", path.c_str()); fclose(f); return true; }
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
                LOG_ERROR("World", "%s: damaged region file\n", files.paths[i]);
        UnloadDirectoryFiles(files);
    }
    b.state().strokes = strokes;
    return b.finish();
}

// ---- Late-join sync (#776): regions in memory, as the save writes them ----

std::vector<TerrainDeltas::Key> Regions_List(const TerrainDeltas& d) {
    std::unordered_set<TerrainDeltas::Key> set;
    for (const auto& kv : d.chunks) set.insert(TerrainDeltas::RegionOf(kv.first));
    std::vector<TerrainDeltas::Key> out(set.begin(), set.end());
    std::sort(out.begin(), out.end());
    return out;
}

// RGN5 bytes of one region's chunks (the save's layout), chunks in key order: the same content always encodes the same,
// so its hash compares across machines.
std::vector<uint8_t> Region_Encode(const TerrainDeltas& d, TerrainDeltas::Key region) {
    std::vector<std::pair<TerrainDeltas::Key, const DeltaChunk*>> chunks;
    for (const auto& kv : d.chunks) if (TerrainDeltas::RegionOf(kv.first) == region) chunks.push_back({ kv.first, kv.second.get() });
    std::sort(chunks.begin(), chunks.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
    std::vector<uint8_t> out;
    auto put = [&](const void* p, size_t n) { const uint8_t* b = (const uint8_t*)p; out.insert(out.end(), b, b + n); };
    put("RGN5", 4);
    const int32_t n = (int32_t)chunks.size();
    put(&n, 4);
    std::vector<uint16_t> enc, encM;
    for (const auto& c : chunks) {
        int cx, cy, cz;
        TerrainDeltas::Unpack(c.first, cx, cy, cz);
        EncodeChunk(*c.second, enc);
        EncodeMaterials(*c.second, encM);
        const int32_t hdr[5] = { cx, cy, cz, (int32_t)(enc.size() * 2), (int32_t)(encM.size() * 2) };
        put(hdr, sizeof hdr);
        put(enc.data(), enc.size() * 2);
        put(encM.data(), encM.size() * 2);
    }
    return out;
}

uint64_t Region_Hash(const std::vector<uint8_t>& bytes) {
    uint64_t h = 14695981039346656037ull;
    for (uint8_t b : bytes) { h ^= b; h *= 1099511628211ull; }
    return h;
}

// Replace one region's chunks on a body with these RGN5 bytes (empty: clear it) and rebuild what shows it.
bool Region_Install(Planetoid& p, TerrainDeltas::Key region, const uint8_t* data, size_t size) {
    std::vector<std::pair<TerrainDeltas::Key, std::shared_ptr<DeltaChunk>>> got;
    if (size > 0) {
        size_t at = 0;
        auto take = [&](void* dst, size_t n) { if (at + n > size) return false; memcpy(dst, data + at, n); at += n; return true; };
        char magic[4]; int32_t n = 0;
        if (!take(magic, 4) || memcmp(magic, "RGN5", 4) != 0 || !take(&n, 4) || n < 0 || n > kRegionChunks * kRegionChunks * kRegionChunks) return false;
        for (int32_t i = 0; i < n; ++i) {
            int32_t hdr[5];
            if (!take(hdr, sizeof hdr) || hdr[3] <= 0 || hdr[4] < 0 || at + (size_t)hdr[3] + (size_t)hdr[4] > size || (hdr[3] & 1) || (hdr[4] & 1)) return false;
            auto chunk = std::make_shared<DeltaChunk>();
            std::vector<uint16_t> words((size_t)hdr[3] / 2), mats((size_t)hdr[4] / 2);
            take(words.data(), (size_t)hdr[3]);
            take(mats.data(), (size_t)hdr[4]);
            if (!DecodeChunk(words.data(), words.size(), *chunk)) return false;
            if (!mats.empty() && !DecodeMaterials(mats.data(), mats.size(), *chunk)) return false;
            const TerrainDeltas::Key k = TerrainDeltas::KeyOf(hdr[0], hdr[1], hdr[2]);
            if (TerrainDeltas::RegionOf(k) != region) return false; // (a chunk outside the region it claims)
            got.push_back({ k, chunk });
        }
        if (at != size) return false;
    }
    TerrainDeltasBuilder b(*p.edits);
    for (auto it = b.state().chunks.begin(); it != b.state().chunks.end();)
        it = TerrainDeltas::RegionOf(it->first) == region ? b.state().chunks.erase(it) : std::next(it);
    for (auto& [k, c] : got) b.setChunk(k, c);
    p.edits = b.finish();
    p.dirtyRegions.insert(region);
    p.lodDirty = true;
    // Rebuild the region's nodes (and its trees' cover) with the new deltas.
    int rx, ry, rz;
    TerrainDeltas::Unpack(region, rx, ry, rz);
    const float side = kRegionChunks * kDeltaChunk * kVoxelSize; // metres
    const Vector3 centre{ (rx + 0.5f) * side, (ry + 0.5f) * side, (rz + 0.5f) * side };
    const float reach = side * 0.8661f + 2.0f;
    RequestRemesh(p, centre, reach);
    if (p.trees) p.trees->onEdit(centre, reach, p.edits);
    return true;
}

bool SaveWorld(Scene& world) {
    if (!DirectoryExists("saves")) MakeDirectory("saves");
    // Written to <save>.tmp and renamed over the save when complete (#346): a crash mid-save leaves the last save whole.
    const std::string tmpPath = SAVE_PATH_STR + ".tmp";
    FILE* f = fopen(tmpPath.c_str(), "w");
    if (!f) return false;
    const DVec3 camWorld = world.renderOrigin + world.camera.position;
    const DVec3 tgtWorld = world.renderOrigin + world.camera.target;
    // v4: the system definition is stored in the save, so a save reloads the system it was made in
    // even if bodies.txt has changed since.
    fprintf(f, "RAYLIBGAME_SAVE 13\n"); // 13 (#346): the player's state in their record (<save>.players/), teams; 12 (#360): ponds' water bricks saved even drained
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
        // (#346) the body's NetId too: a load finds the body by it, not by its place in the list
        fprintf(f, "planetoid %d netid %llu rot %f %f %f strokes %zu\n", bodyIndex, (unsigned long long)p->netId, p->rotation.x, p->rotation.y, p->rotation.z, p->edits->strokes);
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
    const bool written = fflush(f) == 0 && !ferror(f);
    if (fclose(f) != 0 || !written) { std::error_code ec; std::filesystem::remove(tmpPath, ec); LOG_ERROR("World", "%s: could not write the save\n", SAVE_PATH); return false; }
    std::error_code ec;
    std::filesystem::rename(tmpPath, SAVE_PATH_STR, ec); // replaces the old save in one step
    if (ec) { LOG_ERROR("World", "%s: could not replace the save: %s\n", SAVE_PATH, ec.message().c_str()); return false; }
    LOG_INFO("World", "saved %zu strokes to %s (%d region files written)\n", totalStrokes, SAVE_PATH, regionsWritten);
    BackupSave();
    return true;
}

// Backups (#346, user: the newest 5): after each save, the world file and the player records (<save>.players/) are copied
// into <save>.backups/<date-time>/, and only the newest kBackups are kept. The regions are not copied (large, and written
// in place chunk by chunk). RAYLIBGAME_BACKUPS=N keeps N (0: none).
int BackupsKept() { static const int n = getenv("RAYLIBGAME_BACKUPS") ? std::max(0, atoi(getenv("RAYLIBGAME_BACKUPS"))) : 5; return n; }
std::string BackupDir() { return SAVE_PATH_STR + ".backups"; }
void BackupSave() {
    namespace fs = std::filesystem;
    const int keep = BackupsKept();
    if (keep <= 0) return;
    std::error_code ec;
    const fs::path root = BackupDir();
    // The name sorts by time: date, time, milliseconds (two saves in one second still get their own).
    const auto now = std::chrono::system_clock::now();
    const std::time_t t = std::chrono::system_clock::to_time_t(now);
    const int ms = (int)(std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000);
    std::tm tm{};
#ifdef _WIN32
    localtime_s(&tm, &t);
#else
    localtime_r(&t, &tm);
#endif
    char stamp[40];
    snprintf(stamp, sizeof(stamp), "%04d%02d%02d-%02d%02d%02d-%03d", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec, ms);
    const fs::path to = root / stamp;
    fs::create_directories(to, ec);
    const fs::path save = SAVE_PATH_STR;
    fs::copy_file(save, to / save.filename(), fs::copy_options::overwrite_existing, ec);
    if (ec) { LOG_WARN("World", "backup of %s: %s\n", SAVE_PATH, ec.message().c_str()); return; }
    const fs::path players = SAVE_PATH_STR + ".players";
    if (fs::exists(players, ec)) fs::copy(players, to / players.filename(), fs::copy_options::recursive | fs::copy_options::overwrite_existing, ec);
    std::vector<fs::path> all;
    for (const auto& e : fs::directory_iterator(root, ec)) if (e.is_directory()) all.push_back(e.path());
    std::sort(all.begin(), all.end());
    for (size_t i = 0; i + keep < all.size(); ++i) fs::remove_all(all[i], ec); // the oldest beyond the newest `keep`
}

static LoadPhases gLastLoad;
const LoadPhases& LastLoadPhases() { return gLastLoad; }

bool LoadWorldSave(Scene& world) {
    FILE* f = fopen(SAVE_PATH, "r");
    if (!f) return false;
    int version = 0;
    if (fscanf(f, "RAYLIBGAME_SAVE %d", &version) != 1 || version < 2 || version > 13) {
        fclose(f);
        LOG_ERROR("World", "%s: unsupported save version %d\n", SAVE_PATH, version);
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
                if (!ok) LOG_WARN("World", "%s: body %zu: %s\n", SAVE_PATH, i, error.c_str());
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
    if (!ok) { fclose(f); LOG_ERROR("World", "%s: could not read the header\n", SAVE_PATH); return false; }

    worldTime = savedTime; // before regenerating: orbiting bodies are placed for this time
    worldGen.seed = seed;
    worldGen.system = sys;
    worldGen.systemLoaded = true;
    // RAYLIBGAME_LOAD_PROFILE=1: where the load's time goes (--test-persist's load gate).
    using lclk = std::chrono::steady_clock;
    static const bool loadProfile = getenv("RAYLIBGAME_LOAD_PROFILE") && getenv("RAYLIBGAME_LOAD_PROFILE")[0] == '1';
    auto lms = [](lclk::time_point a) { return std::chrono::duration<double, std::milli>(lclk::now() - a).count(); };
    const auto lt0 = lclk::now();
    double regionsMs = 0.0, fluidMs = 0.0, plotsMs = 0.0;
    regenerateWorld(world);
    const double regenMs = lms(lt0);
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
    while (fscanf(f, " planetoid %d", &index) == 1) {
        // (#346) "netid <id>" since format 13's step 3: the body is found by it; older saves by their place in the list
        unsigned long long netid = 0;
        const long at = ftell(f);
        if (fscanf(f, " netid %llu", &netid) != 1) { fseek(f, at, SEEK_SET); netid = 0; }
        if (fscanf(f, " rot %f %f %f %15s %zu", &rot.x, &rot.y, &rot.z, kind, &editCount) != 5) break;
        Planetoid* body = netid ? Planetoid::findNet((NetId)netid) : nullptr;
        if (!body) body = index >= 0 && index < (int)planetoids.size() ? planetoids[index] : nullptr;
        std::shared_ptr<const TerrainDeltas> deltas;
        if (std::string(kind) == "strokes") {
            // v7: the deltas are in the region files.
            const auto rt = lclk::now();
            deltas = LoadRegionsFor(index, dir, editCount);
            regionsMs += lms(rt);
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
        const auto pt = lclk::now();
        auto plots = std::make_unique<PlotSet>();
        if (version >= 8 && !plots->load(f)) LOG_ERROR("World", "%s: damaged plot data\n", SAVE_PATH);
        if (body) body->plots = plots->plots.empty() ? nullptr : std::move(plots);
        if (version >= 9) Factory::Load(f, body ? body->id : 0); // machine contents, applied when the layouts reach the simulation
        plotsMs += lms(pt);
        const auto ft = lclk::now();
        if (version >= 11) Fluid::Load(f, body ? body->id : 0); // flowing water (classified again as it loads)
        fluidMs += lms(ft);
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
    gLastLoad = { regenMs, regionsMs, plotsMs, fluidMs };
    if (loadProfile) printf("load: %.0f ms: regenerate the world %.0f, region deltas %.0f, plots + machines %.0f, fluid %.0f, rest %.0f\n",
                            lms(lt0), regenMs, regionsMs, plotsMs, fluidMs, lms(lt0) - regenMs - regionsMs - plotsMs - fluidMs);
    // Every body is in its saved state now (rotation above): positions kept in a body's frame go back to world space here
    // (#901: the player's record was read in the extras, before the rotation came back, and stood on the wrong ground).
    if (loadDoneHook) loadDoneHook();
    LOG_INFO("World", "loaded %s (seed %u)\n", SAVE_PATH, seed);
    return true;
}

void loadWorld(Scene& world) {
    // Load the world scene from a file or initialize it
    LOG_INFO("Scene", "Loading world scene...\n");
    std::ifstream file("assets/scenes/world");
    if (!file.is_open()) {
        std::cerr << "Failed to open world scene file." << std::endl;
        return;
    }

    std::string line;
    while (getline(file, line)) {
        if (line.empty()) continue; // Skip empty lines

        LOG_DEBUG("Scene", "Processing line: %s\n", line.c_str());
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
                LOG_DEBUG("Scene", "Creating GameObject: %s\n", name.c_str());
                std::string modelPath = "assets/models/" + name + "/" + name + ".obj";
                std::string texturePath = "assets/models/" + name + "/diffuse.png";
                Model model;
                if (!FileExists(modelPath.c_str())) {
                    std::cerr << "Model file does not exist: " << modelPath << std::endl;
                    continue;
                } else {
                    model = LoadModel(modelPath.c_str());
                    LOG_DEBUG("Scene", "meshCount: %d, materialCount: %d, boneCount: %d\n", model.meshCount, model.materialCount, model.boneCount);
                    LOG_DEBUG("Scene", "bones: %p, bindPose: %p\n", model.bones, model.bindPose);
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
                    LOG_WARN("Scene", "Texture file does not exist or materials/maps not allocated: %s\n", texturePath.c_str());
                }
                world.objects.push_back(std::make_unique<GameObject>(type, name, position, rotation, color, scale, model));
                LOG_DEBUG("Scene", "GameObject created: %s\n", name.c_str());
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


Vector3 RiverCurrent(Planetoid* p, Vector3 renderPoint) {
    if (!p || !p->field) return { 0, 0, 0 };
    const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
    const RiverNetwork* net = pf ? pf->rivers() : nullptr;
    if (!net) return { 0, 0, 0 };
    const Vector3 local = p->toLocal(renderPoint);
    const float d = Vector3Length(local);
    if (d < 1.0f || Fluid::Simulated(*p, local)) return { 0, 0, 0 };
    const Vector3 dir = Vector3Scale(local, 1.0f / d);
    if (pf->waterSourceAt(dir).kind != 1) return { 0, 0, 0 };   // rivers only
    const RiverNetwork::Hit h = net->nearest(dir, 60.0f);
    if (!h.found || h.dist >= h.halfWidth + 3.0f) return { 0, 0, 0 };
    Vector3 flow = Vector3Subtract(h.down, Vector3Scale(dir, Vector3DotProduct(h.down, dir)));
    const float fl = Vector3Length(flow);
    if (fl < 1e-4f) return { 0, 0, 0 };
    const float speed = sqrtf(std::max(0.0f, 1.0f - h.dist / (h.halfWidth + 3.0f))); // as water.fs shows it (voxel.cpp)
    return p->dirToWorld(Vector3Scale(flow, speed / fl));
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
    if (Fluid::Wet(*p, local)) return true; // simulated water
    if (Fluid::Simulated(*p, local)) return false; // a simulated cell decides (#360: a drained pond is dry)
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

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): body-local <-> world frames under spin ----
#include "unit_helpers.h"
#include <random>

namespace {
// A body at `pos` (no GL: chunks are never built), turned by `rot` degrees about X, Y, Z.
std::unique_ptr<Planetoid> MakeBody(DVec3 pos, Vector3 rot, float radius = 60000.0f) {
    auto b = std::make_unique<Planetoid>("unit", pos, rot, 1.0f, 42u, radius);
    b->rotation = rot;
    return b;
}
double Dist(const DVec3& a, const DVec3& b) { return (a - b).length(); }
}

// local -> world -> local is exact to micrometres 60 km from the centre, the body 100 km from the origin, at any angle.
TEST(Frames, LocalWorldRoundTripAt60km) {
    std::mt19937 rng(1);
    std::uniform_real_distribution<float> ang(0.0f, 360.0f);
    std::uniform_real_distribution<double> dir(-1.0, 1.0);
    for (int i = 0; i < 500; ++i) {
        auto b = MakeBody({ 100000.0, -20000.0, 5000.0 }, { ang(rng), ang(rng), ang(rng) });
        DVec3 l{ dir(rng), dir(rng), dir(rng) };
        l = l * (60000.0 / std::max(1e-9, l.length()));
        const DVec3 w = b->toWorldD(l);
        EXPECT_NEAR(Dist(w, b->worldPos), 60000.0, 1e-6) << "a rotation keeps distances";
        EXPECT_LT(Dist(b->toLocalD(w), l), 1e-6);
    }
}

// The double transform is the float model matrix's (raymath MatrixRotateXYZ), so chunks and physics agree.
TEST(Frames, DoubleTransformMatchesTheModelMatrix) {
    std::mt19937 rng(2);
    std::uniform_real_distribution<float> ang(0.0f, 360.0f), c(-200.0f, 200.0f);
    for (int i = 0; i < 200; ++i) {
        auto b = MakeBody({ 0.0, 0.0, 0.0 }, { ang(rng), ang(rng), ang(rng) }, 500.0f);
        const Vector3 l = { c(rng), c(rng), c(rng) };
        const Vector3 wf = b->toWorld(l);
        const DVec3 wd = b->toWorldD(DVec3(l));
        EXPECT_LT(Dist(DVec3(wf), wd), 2e-4) << "rotation " << b->rotation.x << " " << b->rotation.y << " " << b->rotation.z;
        const Vector3 dw = b->dirToWorld(Vector3Normalize(l));
        const Vector3 dq = Vector3RotateByQuaternion(Vector3Normalize(l), b->rotationQuat());
        EXPECT_LT(Vector3Distance(dw, dq), 1e-5f);
        EXPECT_LT(Vector3Distance(b->dirToLocal(dw), Vector3Normalize(l)), 1e-5f);
    }
}

// The float toLocal (render space -> body-local) is what collision samples the density through (#1072). A point near the
// render origin, 60 km out from a turning body's centre, maps to body-local within a body-local float's rounding of the
// double path, at every angle - the old float matrix inverse was off by mm to ~2 cm and differently from angle to angle.
TEST(Frames, FloatToLocalMatchesTheDoublePathAt60km) {
    auto b = MakeBody(DVec3(0.0, -60010.0, 0.0), { 0.0f, 0.0f, 0.0f });
    b->position = { 0.0f, -60010.0f, 0.0f }; // render origin at the world origin: the body's centre 60 km below it
    double worst = 0.0;
    for (int k = 0; k < 400; ++k) {
        b->rotation = { 0.0f, 200.0f + 0.003f * (float)k, 0.0f }; // a turning planet, ~1 m of surface travel a step
        const Vector3 c = { 0.37f, 1.20f, -0.81f };              // a capsule sphere near the player (render space)
        const DVec3 ref = b->toLocalD(DVec3(c.x, c.y, c.z));   // render origin = world origin here
        const Vector3 l = b->toLocal(c);
        worst = std::max(worst, (DVec3(l.x, l.y, l.z) - ref).length());
        EXPECT_LT(Vector3Distance(b->toWorld(l), c), 0.01f);   // and back
    }
    EXPECT_LT(worst, 0.004); // body-local float rounding at 60 km: 0.8 mm measured (the old float matrix: 8.6 mm)
}

// A body-fixed point's velocity from the spin (surfaceVelocity) is the rate its world position changes: the player rides
// on it, and a wrong one slides everything standing on a turning planet (#749, #856).
namespace {
// The reference: the body-fixed point's world velocity by a central difference of the rotation in double (BodyRotationD's
// formula with double angles: the stored float angles step 1.5e-5 deg at 200 deg, too coarse to difference).
DVec3 SpinVelocityD(Vector3 rotDeg, Vector3 spinDeg, const DVec3& l) {
    auto world = [&](double t) {
        const double ax = -((double)rotDeg.x + spinDeg.x * t) * (PI / 180.0), ay = -((double)rotDeg.y + spinDeg.y * t) * (PI / 180.0),
                     az = -((double)rotDeg.z + spinDeg.z * t) * (PI / 180.0);
        const double cx = cos(ax), sx = sin(ax), cy = cos(ay), sy = sin(ay), cz = cos(az), sz = sin(az);
        const double r[3][3] = { { cz * cy, sz * cy, -sy }, { cz * sy * sx - sz * cx, sz * sy * sx + cz * cx, cy * sx }, { cz * sy * cx + sz * sx, sz * sy * cx - cz * sx, cy * cx } };
        return DVec3(r[0][0] * l.x + r[0][1] * l.y + r[0][2] * l.z, r[1][0] * l.x + r[1][1] * l.y + r[1][2] * l.z, r[2][0] * l.x + r[2][1] * l.y + r[2][2] * l.z);
    };
    const double h = 1e-3;
    return (world(h) - world(-h)) * (1.0 / (2.0 * h));
}
// Worst |surfaceVelocity - reference| / speed over `n` random orientations (0: the identity only).
double WorstSurfaceVelocityError(int n, Vector3 spin) {
    std::mt19937 rng(749);
    std::uniform_real_distribution<float> ang(0.0f, 360.0f);
    double worst = 0.0;
    for (int i = 0; i <= n; ++i) {
        const Vector3 rot = i == 0 ? Vector3{ 0, 0, 0 } : Vector3{ ang(rng), ang(rng), ang(rng) };
        auto b = MakeBody({ 0.0, 0.0, 0.0 }, rot);
        b->spin = spin;
        const DVec3 l{ 31000.0, 42000.0, -28000.0 }; // ~59 km out
        const DVec3 ref = SpinVelocityD(rot, spin, l);
        const Vector3 v = b->surfaceVelocity(b->toWorldD(l).toFloat()); // render origin = the world origin; the body sits there
        worst = std::max(worst, Dist(DVec3(v), ref) / std::max(1e-9, ref.length()));
    }
    return worst;
}
}

// A planet's spin (about its Y axis, as bodies.txt gives it) at any orientation.
TEST(Frames, SurfaceVelocityIsTheSpinOfABodyFixedPoint) {
    EXPECT_LT(WorstSurfaceVelocityError(500, { 0, 0.1f, 0 }), 0.002) << "worst relative error over 500 orientations";
}

// Disabled until #872: a spin about several Euler axes at once (ring rocks, asteroids) is off by up to 1.5%: the angular
// velocity comes from a float quaternion's acosf(w) near w = 1. Run with --gtest_also_run_disabled_tests.
TEST(Frames, DISABLED_SurfaceVelocityOfAMultiAxisSpin) {
    EXPECT_LT(WorstSurfaceVelocityError(500, { 0.02f, -0.1f, 0.05f }), 0.002) << "worst relative error over 500 orientations";
}

// A whole day of turning: the same local point comes back to the same place, and its local coordinates never drift.
TEST(Frames, ADayOfSpinKeepsLocalPointsFixed) {
    auto b = MakeBody({ 0.0, 0.0, 0.0 }, { 0, 0, 0 });
    const DVec3 l{ 1234.5, 59900.25, -321.125 };
    const DVec3 w0 = b->toWorldD(l);
    for (int deg = 0; deg <= 360; deg += 15) {
        b->rotation.y = (float)deg;
        EXPECT_LT(Dist(b->toLocalD(b->toWorldD(l)), l), 1e-6) << deg << " deg";
    }
    EXPECT_LT(Dist(b->toWorldD(l), w0), 0.1) << "360 deg is a full turn (raylib's float DEG2RAD: within 10 cm at 60 km)";
}

// Plan 11 F1 (#667): a point on a turning body moves evenly, step by step. Near 360 deg the float Euler angle steps by
// 3.05e-5 deg, so a surface point 60 km out jumped by up to ~3 cm between steps (frames/still: 13-31 mm). With the double
// pose (RAYLIBGAME_FRAMES=1) every 1/60 s step moves it the same arc, to well under a millimetre.
TEST(Frames, DoublePoseMovesASurfacePointEvenly) {
    auto stepSpread = [](bool doublePose) {
        kFramesPose = doublePose;
        auto b = MakeBody({ 0.0, 0.0, 0.0 }, { 0, 0, 0 });
        b->spin = { 0.0f, 0.1f, 0.0f }; // Terra's day: 360 deg an hour
        const DVec3 l{ 1234.5, 31000.0, 51300.0 }; // a surface point ~60 km out, off the axis
        double lo = 1e9, hi = 0.0;
        DVec3 last{};
        for (int k = 0; k <= 600; ++k) {
            const double a = 359.0 + 0.1 * k / 60.0; // 10 s of spin up to 360 deg, 1/60 s steps
            b->spinAngle[1] = a; b->spinAngleSet = true;
            b->rotation = { 0.0f, (float)a, 0.0f }; // advanceBodies stores the float from the double
            const DVec3 w = b->toWorldD(l);
            if (k > 0) { const double d = Dist(w, last); lo = std::min(lo, d); hi = std::max(hi, d); }
            last = w;
        }
        kFramesPose = false;
        return hi - lo;
    };
    const double floatSpread = stepSpread(false), doubleSpread = stepSpread(true);
    EXPECT_GT(floatSpread, 0.005) << "the float angle's steps: uneven by cm (the test still sees the old path)";
    EXPECT_LT(doubleSpread, 1e-4) << "the double pose: even steps";
}
#endif
