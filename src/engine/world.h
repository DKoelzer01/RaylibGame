#ifndef world_h
#define world_h

#include "scene.h"
#include "voxel.h"
#include "worldtypes.h"
#include "cubemarch.h"
#include "utils/logger.h"
#include "profiler.h"

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <cstdio>
#include <functional>

// What generateWorld() builds: the star system from bodiesPath (see bodies.h), plus a seed.
// seed == 0 means "pick a random seed"; lastSeed records what was actually used so an
// interesting world can be reproduced. Tweakable from the debug UI.
struct WorldGenSettings {
    unsigned int seed = 0;
    unsigned int lastSeed = 0;
    std::string bodiesPath = "resources/data/bodies.txt";
    SystemDef system;          // loaded from bodiesPath on first use (or from a save)
    bool systemLoaded = false;
    // Command-line overrides of the classic start-body-plus-ring layout, applied whenever
    // bodiesPath is (re)loaded. Negative / false = keep the file's value.
    int count = -1;            // total bodies: the first ring gets count minus the single bodies
    int size = -1;             // radius of the start body and the first ring's base size
    float systemRadius = -1.0f;// first ring's extent
    bool hasFirst = false;     // start body position
    Vector3 first = { 0, 0, 0 };
};
extern WorldGenSettings worldGen;
bool ReloadSystemDef();        // worldGen.bodiesPath -> worldGen.system (+ overrides); false on errors

struct LodSettings {
    float lodFactor = 2.0f; // split a node while the camera is closer than lodFactor * node size
    // Zoomed view (telescope): nodes inside the view cone refine as if zoom times closer. Set with
    // SetLodZoom, which re-selects LOD when the zoom or the aim changes.
    float zoom = 1.0f;
    Vector3 zoomDir = { 0, 0, 1 };  // render-space view direction
    float zoomHalfAngle = 0.0f;     // rad, to the screen corners
};
extern LodSettings lodSettings;
void SetLodZoom(Scene& world, float zoom, Vector3 dir, float halfAngle); // zoom 1: off

struct WorldStats {
    int queuedJobs = 0, runningJobs = 0, pendingUploads = 0, uploadsLastFrame = 0, workerThreads = 0;
    int chunks = 0;  // meshes drawn this frame
    int nodes = 0;   // LOD nodes alive (built, building or retained)
    int edits = 0;
    float avgBuildMs = 0.0f;            // worker time per chunk (sampling + meshing + normals)
    float avgSampleMs = 0.0f;           // of which sampling the density field
    float lastGenerationSeconds = 0.0f; // from generateWorld() until every chunk was uploaded
};

struct TerrainHit {
    Planetoid* planetoid = nullptr;
    Vector3 point{};      // world space
    Vector3 local{};      // planetoid-local voxel space
    Vector3 normal{};     // world space
    float distance = 0.0f;
    VoxelMaterial material = VoxelMaterial::Rock;
};

void worldHandler(Scene& world);
void generateWorld(Scene& world);      // creates planetoids; chunks stream in via worldUpdate()
void regenerateWorld(Scene& world);
void BeginWorldStreaming();            // for games that add bodies themselves: WorldIdle() waits for them
void loadWorld(Scene& world);          // extra objects from assets/scenes/world

// Per frame: requests chunks for new planetoids and uploads finished meshes (time-boxed).
// advanceTime = false while paused: streaming/uploads continue but bodies stop spinning.
void advanceBodies(Scene& world, float dt);  // spin and orbit bodies; call BEFORE updating anything that rides them
double WorldTime();                          // seconds of simulated time (advances only while playing; saved)
void PlaceBodyOnOrbit(Planetoid& p);             // position (and velocity) on its analytic orbit now
// Region files of a body's deltas (id = body index, or a ring rock's id): save the changed regions, load all.
int SaveRegionsFor(const TerrainDeltas& edits, std::unordered_set<TerrainDeltas::Key>& dirty, bool& allDirty, std::string& savedTo, int body, const std::string& dir);
std::shared_ptr<const TerrainDeltas> LoadRegionsFor(int body, const std::string& dir, size_t strokes);
class RingSystem;
std::vector<std::unique_ptr<RingSystem>>& Rings();
void SetWorldTime(double t);
void worldUpdate(Scene& world, float uploadBudgetMs = 3.0f); // LOD, streaming, uploads
bool WorldIdle();                      // no queued/running jobs and nothing left to upload
const WorldStats& GetWorldStats();

bool RaycastTerrain(Scene& world, Ray ray, float maxDistance, TerrainHit& hit);
void ApplyTerrainEdit(Planetoid& planetoid, const TerrainEdit& edit); // local-space brush
void ApplyTerrainEdits(Planetoid& planetoid, const TerrainEdit* edits, size_t count); // many strokes, one snapshot
// Rebuild the chunks (every level) around a body-local point: water claimed by the simulation (fluid.h).
void RequestRemesh(Planetoid& planetoid, Vector3 localCentre, float radius);

// Physics helpers (render-space points). Distance is approximate: (ISO - density) / |gradient|.
bool TerrainSurfaceDistance(Planetoid& p, Vector3 renderPoint, float& dist, Vector3& outwardNormal);
Planetoid* NearestPlanetoid(Scene& world, Vector3 renderPoint); // by distance to the surface radius
// In the sea: below sea level over open ocean and above its floor (a dry cave below sea level is not).
bool IsUnderwater(Scene& world, Vector3 renderPoint);
void SetWorldSavePath(const char* path); // bench runs save to their own file, never the player's
const char* WorldSavePath();
class JobSystem;
JobSystem& WorldJobs(); // the shared worker pool (chunk builds, foliage tiles)
// Gravity (m/s^2, render space) summed over every body, and one body's share of it: inverse-square
// outside a body, linear to zero at its centre inside, surfaceGravity at its radius.
Vector3 GravityAt(Scene& world, Vector3 renderPoint);
Vector3 BodyGravityAt(const Planetoid& body, Vector3 renderPoint);

bool SaveWorld(Scene& world);          // saves/world.sav
// Game-specific save data. The load hook must read everything up to and including "end_extras".
extern std::function<void(FILE*)> saveExtrasHook;
extern std::function<void(FILE*)> loadExtrasHook;
bool LoadWorldSave(Scene& world);

#endif
