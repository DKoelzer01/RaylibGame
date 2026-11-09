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
};
extern LodSettings lodSettings;

struct WorldStats {
    int queuedJobs = 0, runningJobs = 0, pendingUploads = 0, uploadsLastFrame = 0, workerThreads = 0;
    int chunks = 0;  // meshes drawn this frame
    int nodes = 0;   // LOD nodes alive (built, building or retained)
    int edits = 0;
    float avgBuildMs = 0.0f;            // worker time per chunk (sampling + meshing + normals)
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
void advanceBodies(Scene& world, float dt);  // spin bodies; call BEFORE updating anything that rides them
void worldUpdate(Scene& world, float uploadBudgetMs = 3.0f); // LOD, streaming, uploads
bool WorldIdle();                      // no queued/running jobs and nothing left to upload
const WorldStats& GetWorldStats();

bool RaycastTerrain(Scene& world, Ray ray, float maxDistance, TerrainHit& hit);
void ApplyTerrainEdit(Planetoid& planetoid, const TerrainEdit& edit); // local-space brush

// Physics helpers (render-space points). Distance is approximate: (ISO - density) / |gradient|.
bool TerrainSurfaceDistance(Planetoid& p, Vector3 renderPoint, float& dist, Vector3& outwardNormal);
Planetoid* NearestPlanetoid(Scene& world, Vector3 renderPoint); // by distance to the surface radius
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
