#include <raylib.h>
#include <string>
#include <vector>
#include <random>
#include "mainmenu.h"
#include "world.h"
#include "scene.h"
#include "object.h"
#include "utils/logger.h"
#include "rlights.h"
#include "raymath.h"
#include "debugui.h"
#include "profiler.h"
#include "terraintool.h"
#include "player.h"
#include "inventory.h"
#include "gameui.h"
#include "gamesave.h"
#include "audio.h"
#include "rlImGui.h"
#include "ship.h"
#include "survival.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif

using namespace std;

Logger logger("log.txt");

int gameState = 0; // 0: Main Menu, 1: Game, 2: Pause Menu

// --bench <seconds>: windowed 1920x1080, uncapped FPS, skips the menu, holds a fixed
// camera on the planetoid, prints frame stats to stdout and bench.txt, then exits.


int main(int argc, char** argv)
{
    double benchSeconds = 0.0;
    bool startWithOverlay = false;
    int benchDigStrokes = 0;     // --test-dig N: dig N strokes at the crosshair once the world is loaded
    bool benchSaveLoad = false;  // --test-saveload: then save, reload from disk, and bench the reloaded world
    double benchWalkSeconds = 0.0; // --test-walk S: spawn on the surface and auto-walk forward for S seconds
    bool benchWalkMode = false;    // --bench-walk: the player walks forward for the whole bench (frame time while streaming)
    double walkStart = -1.0;
    bool benchStandOnly = false; // --test-stand S: like --test-walk but standing still
    int benchMineStrokes = 0;     // --test-mine N: dig N strokes at the feet, craft, build (after any walk test)
    bool mineStarted = false, mineStoppedEarly = false;
    int minePhase = 0, minedStrokes = 0, mineCrafted = 0;
    double benchFlySeconds = 0.0, flyStart = -1.0; // --test-fly S: autopilot the ship to the nearest other asteroid
    bool benchSuffocate = false;   // --test-suffocate: run out of oxygen away from the ship, die, respawn
    double benchOrbitSeconds = 0.0, orbitStart = -1.0; // --test-orbit S: coast on a circular orbit for S seconds
    float orbitMin = 0, orbitMax = 0;
    bool orbitCrashed = false;
    bool benchRide = false;        // --test-ride: stand on a spinning body, then jump
    bool benchPark = false;        // --test-park: park the ship on a spinning body
    bool sunOverride = false;      // --sun x y z: sun direction (e.g. to test night-side lighting)
    int benchFps = 0;              // --fps N: cap the bench frame rate
    bool lookDown = false;         // --look-down: tilt the player view toward the ground (lighting tests)
    Vector3 sunArg = { 1, 0, 0 };
    double parkStart = -1.0;
    int parkPhase = 0;
    bool benchGravity = false;     // --test-gravity: gravity field checks, ship drop and hover
    bool benchPrecision = false;   // --test-precision: float vs double placement error on body 0 (use a large body)
    bool benchMantle = false;      // --test-mantle: planet terrain checks on body 0 (height band, mantle, ore, caves)
    bool crackCheck = false;       // --crack-check: no skybox, magenta clear; reports the share of magenta pixels
    double benchUpdateMs = 0.0, benchDrawMs = 0.0; long benchTimedFrames = 0; // CPU time in worldUpdate / drawScene
    bool benchMeasuring = false;   // this frame is inside the measured bench window
    int precisionPhase = 0;
    double precisionStart = 0.0;
    double gravStart = -1.0;
    int gravPhase = 0;
    float gravAlt0 = 0.0f, gravEnergy0 = 0.0f, gravMaxImpact = 0.0f;
    bool parkAnchored = false;
    Vector3 parkLocal0{}, parkFwd0{};
    double rideStart = -1.0;
    int ridePhase = 0;
    Vector3 rideFwd0{}, ridePos0{}, rideTakeoff{};
    float rideViewDrift = 0, ridePosDrift = 0, rideSurfaceSpeed = 0;
    int rideFrames = 0, rideGroundedFrames = 0; // grounded share while standing (ride test)
    float rideSpin = 10.0f;        // --ride-spin D: spin rate for --test-ride (deg/s)
    bool dumpSystem = false;       // --dump-system: print every body after generation
    bool loadSaveAtStart = false;  // --load: load saves/world.sav before starting
    bool spinTest = false;         // --spin-test: textures must stay attached to a spinning body
    float spinRate = 30.0f;        // --spin-rate DEG_PER_SEC
    int spinPhase = 0;
    Vector3 spinCamLocal{}, spinTgtLocal{};
    double suffocateStart = -1.0;
    DVec3 flyFrom{};
    float flyEnergy = 0.0f, flyMaxSpeed = 0.0f;
    DVec3 walkFrom{};
    float walkMinDensity = 0, walkMaxDensity = 0;
    int walkGroundedFrames = 0, walkFrames = 0;
    Vector3 benchCamOffset = { -70.0f, 35.0f, 70.0f }; // --cam dx dy dz: bench camera offset from the planetoid
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--bench") benchSeconds = (i + 1 < argc) ? atof(argv[i + 1]) : 10.0;
        if (std::string(argv[i]) == "--ui") DebugUI_Toggle(); // open the debug panel at start
        if (std::string(argv[i]) == "--seed" && i + 1 < argc) worldGen.seed = (unsigned int)strtoul(argv[i + 1], nullptr, 10);
        if (std::string(argv[i]) == "--nocull") frustumCullingEnabled = false;
        if (std::string(argv[i]) == "--overlay") startWithOverlay = true;
        if (std::string(argv[i]) == "--cam" && i + 3 < argc) { benchCamOffset = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
        if (std::string(argv[i]) == "--test-dig" && i + 1 < argc) benchDigStrokes = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-saveload") benchSaveLoad = true;
        if (std::string(argv[i]) == "--test-walk" && i + 1 < argc) benchWalkSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--bench-walk") benchWalkMode = true;
        if (std::string(argv[i]) == "--test-suffocate") benchSuffocate = true;
        if (std::string(argv[i]) == "--test-orbit" && i + 1 < argc) benchOrbitSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-ride") benchRide = true;
        if (std::string(argv[i]) == "--dump-system") dumpSystem = true;
        if (std::string(argv[i]) == "--load") loadSaveAtStart = true;
        if (std::string(argv[i]) == "--ride-spin" && i + 1 < argc) rideSpin = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-park") benchPark = true;
        if (std::string(argv[i]) == "--test-gravity") benchGravity = true;
        if (std::string(argv[i]) == "--test-precision") benchPrecision = true;
        if (std::string(argv[i]) == "--test-mantle") benchMantle = true;
        if (std::string(argv[i]) == "--crack-check") crackCheck = true;
        if (std::string(argv[i]) == "--sun" && i + 3 < argc) { sunOverride = true; sunArg = Vector3Normalize({ (float)atof(argv[i + 1]), (float)atof(argv[i + 2]), (float)atof(argv[i + 3]) }); }
        if (std::string(argv[i]) == "--look-down") lookDown = true;
        if (std::string(argv[i]) == "--fps" && i + 1 < argc) benchFps = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--headlamp") survival.headlamp = true;      // lamp on from the start
        if (std::string(argv[i]) == "--inventory") GameUI_ToggleInventory();     // inventory open (player glow)
        if (std::string(argv[i]) == "--spin-test") spinTest = true;
        if (std::string(argv[i]) == "--test-flycam") {
            // Pitch a free camera through two full loops (past straight up and down) in small steps,
            // the way customUpdateCamera does, and check it never flips: each step turns the view by
            // the same small angle and up stays perpendicular to the view.
            Camera3D c{ { 0, 0, 0 }, { 0, 0, 1 }, { 0, 1, 0 }, 60.0f, CAMERA_PERSPECTIVE };
            float worstStep = 0.0f, worstOrtho = 0.0f;
            Vector3 prev = { 0, 0, 1 };
            for (int s = 0; s < 400; ++s) {
                CameraPitch(&c, 0.02f, false, false, true);
                CameraYaw(&c, 0.003f, false);
                const Vector3 fwd = Vector3Normalize(Vector3Subtract(c.target, c.position));
                worstStep = std::max(worstStep, Vector3Angle(prev, fwd));
                worstOrtho = std::max(worstOrtho, fabsf(Vector3DotProduct(fwd, Vector3Normalize(c.up))));
                prev = fwd;
            }
            printf("bench: flycam pitched 8 rad (past both poles): largest per-step turn %.4f rad (expected ~0.02), max |dot(view, up)| %.4f\n",
                   worstStep, worstOrtho);
        }
        if (std::string(argv[i]) == "--spin-rate" && i + 1 < argc) spinRate = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-fly" && i + 1 < argc) benchFlySeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-mine" && i + 1 < argc) benchMineStrokes = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-stand" && i + 1 < argc) { benchWalkSeconds = atof(argv[i + 1]); benchStandOnly = true; }
        if (std::string(argv[i]) == "--size" && i + 1 < argc) worldGen.size = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--bodies" && i + 1 < argc) worldGen.bodiesPath = argv[i + 1];   // system definition file
        if (std::string(argv[i]) == "--count" && i + 1 < argc) worldGen.count = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--system" && i + 1 < argc) worldGen.systemRadius = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--lod" && i + 1 < argc) lodSettings.lodFactor = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--first" && i + 3 < argc) { worldGen.hasFirst = true; worldGen.first = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
    }
    const bool bench = benchSeconds > 0.0;
    benchMode = bench;

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    if (!bench) {
        SetConfigFlags(FLAG_VSYNC_HINT); // Enable VSync for smoother rendering
        SetConfigFlags(FLAG_WINDOW_RESIZABLE); // Allow window resizing
        SetConfigFlags(FLAG_WINDOW_UNDECORATED); // Remove window decorations for a cleaner look
        SetConfigFlags(FLAG_WINDOW_HIGHDPI); // Enable High DPI support
        InitWindow(GetMonitorWidth(0), GetMonitorHeight(0), "My first RAYLIB program!");
        SetTargetFPS(60);
    } else {
        InitWindow(1920, 1080, "RaylibGame bench");
        SetTargetFPS(benchFps); // --fps N caps it (frame-rate-dependence checks); 0 = uncapped
    }
    SetExitKey(KEY_NULL);
    SetTraceLogLevel(LOG_WARNING); // Set log level to warning to reduce output noise
    rlSetClipPlanes(0.5, 60000.0); // space scale: see other bodies kilometres away
    
    DebugUI_Init();
    if (!bench) Audio_Init();
    itemDB.load("resources/data/items.txt", "resources/data/recipes.txt");
    GameSave_Register();
    {   // Scenes own GPU resources: keep them in a scope that ends before CloseWindow().
    Scene mainMenu("mainMenu", true); // Main menu scene;
    Scene world("world", false);

    bool windowKeepAlive = true;

    mainMenuHandler(mainMenu);
    worldHandler(world);
    if (loadSaveAtStart && !LoadWorldSave(world)) printf("bench: --load: could not load saves/world.sav\n");
    if (dumpSystem) {
        for (auto& obj : world.objects)
            if (auto* p = dynamic_cast<Planetoid*>(obj.get()))
                printf("body %s %s pos %.3f %.3f %.3f radius %.3f seed %u spin %.4f g %.4f\n", BodyKindName(p->kind), p->name.c_str(), p->worldPos.x, p->worldPos.y, p->worldPos.z,
                       p->radius, p->seed, p->spin.y, p->surfaceGravity);
    }
    if (Planetoid* start = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get())) player.spawnOn(world, *start);
    if (!world.objects.empty()) {
        // Park the ship a few metres ahead of the player.
        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
        ship.parkNear(world, Vector3Add(feet, Vector3Scale(player.forward, 8.0f)), player.up);
    }
    world.onDraw3D = [&world] { if (!ship.piloted) TerrainTool_Draw3D(); ship.draw(world); };

    const double benchWarmup = 2.0;
    double benchStart = 0.0, benchMaxFrame = 0.0;
    long benchFrames = 0;
    if (bench) {
        gameState = 1;
        mainMenu.isActive = false;
        world.isActive = true;
        world.camera.position = Vector3Add(world.camera.target, benchCamOffset);
        world.debugDraw = startWithOverlay;
        player.flyMode = benchWalkSeconds <= 0.0 && !benchWalkMode; // the bench drives the camera directly unless walking
        player.autoWalk = benchWalkMode;                              // --bench-walk: walk through the measured window
        benchStart = -1.0; // starts once the world has finished streaming in
    }

    while (windowKeepAlive && !WindowShouldClose())
    {
        if (bench) {
            if (benchStart < 0.0) {
                if (WorldIdle() && benchDigStrokes > 0) {
                    // Scripted edit test: carve a tunnel straight in at the crosshair.
                    TerrainHit hit;
                    const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
                    if (RaycastTerrain(world, Ray{ world.camera.position, dir }, 300.0f, hit)) {
                        for (int s = 0; s < benchDigStrokes; ++s) {
                            const Vector3 local = hit.planetoid->toLocal(Vector3Add(hit.point, Vector3Scale(dir, 2.0f * s)));
                            ApplyTerrainEdit(*hit.planetoid, TerrainEdit{ local, 8.0f, -1.0f });
                        }
                        printf("bench: dug %d strokes at %.1f %.1f %.1f (%s)\n", benchDigStrokes, hit.point.x, hit.point.y, hit.point.z, VoxelMaterialName(hit.material));
                    } else {
                        printf("bench: --test-dig found no terrain at the crosshair\n");
                    }
                    benchDigStrokes = 0;
                } else if (WorldIdle() && benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && benchFlySeconds <= 0.0) {
                    const int editsBefore = GetWorldStats().edits;
                    const int itemsBefore = inventory.count("stone") + inventory.count("regolith");
                    const bool ok = SaveWorld(world) && LoadWorldSave(world); // restores camera + player from the file
                    printf("bench: save/load %s (%d edits, %d stone+regolith before / %d after)\n", ok ? "ok" : "FAILED",
                           editsBefore, itemsBefore, inventory.count("stone") + inventory.count("regolith"));
                    benchSaveLoad = false;
                } else if (benchSuffocate && benchWalkSeconds <= 0.0 && benchFlySeconds <= 0.0 && benchMineStrokes == 0) {
                    // Scripted survival: out of oxygen, far from the ship -> die -> respawn at the ship.
                    if (suffocateStart < 0.0) {
                        player.flyMode = false;
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        ship.parkNear(world, Vector3Add(feet, Vector3Scale(player.forward, 30.0f)), player.up); // out of life-support range
                        survival.oxygen = 0.0f;
                        survival.health = 15.0f;
                        suffocateStart = GetTime();
                    } else if (survival.deaths > 0 && !survival.dead) {
                        const float d = Vector3Distance((player.position - world.renderOrigin).toFloat(), ship.renderPos(world));
                        printf("bench: suffocated and respawned after %.1f s; O2 %.0f HP %.0f, %.1f m from ship, deaths %d\n",
                               GetTime() - suffocateStart, survival.oxygen, survival.health, d, survival.deaths);
                        benchSuffocate = false;
                    } else if (GetTime() - suffocateStart > 15.0) {
                        printf("bench: --test-suffocate FAILED (dead %d, HP %.1f, in life support %d)\n", survival.dead, survival.health, survival.inLifeSupport);
                        benchSuffocate = false;
                    }
                } else if (benchOrbitSeconds > 0.0 && WorldIdle()) {
                    // Scripted orbit: place the player at 1.35 radii with circular-orbit speed, no input,
                    // and check the distance from the centre stays put (no crash, no escape).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const float R = p->radius * p->scale;
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    if (orbitStart < 0.0) {
                        const float d = 1.35f * R;
                        const float g = Vector3Length(BodyGravityAt(*p, Vector3Add(p->position, { 0.0f, d, 0.0f })));
                        player.flyMode = false;
                        player.grounded = false;
                        player.position = world.renderOrigin + Vector3Add(p->position, { 0.0f, d, 0.0f });
                        player.velocity = { sqrtf(g * d), 0.0f, 0.0f };
                        orbitStart = GetTime();
                        orbitMin = 1e9f; orbitMax = 0.0f; orbitCrashed = false;
                        printf("bench: orbit test: radius %.1f, start %.1f m from centre at %.2f m/s (period %.0f s)\n", R, d, sqrtf(g * d), 2.0 * PI * d / sqrtf(g * d));
                    } else {
                        const float d = Vector3Distance(feet, p->position);
                        orbitMin = std::min(orbitMin, d);
                        orbitMax = std::max(orbitMax, d);
                        if (player.grounded) orbitCrashed = true;
                        if (GetTime() - orbitStart > benchOrbitSeconds || orbitCrashed) {
                            printf("bench: orbit after %.0f s: distance from centre %.1f..%.1f m, %s\n", GetTime() - orbitStart, orbitMin, orbitMax,
                                   orbitCrashed ? "CRASHED into the surface" : (orbitMax > 3.0f * R ? "ESCAPED" : "stable"));
                            benchOrbitSeconds = 0.0;
                        }
                    }
                } else if (benchMantle && WorldIdle()) {
                    // Planet terrain checks on body 0 (run with --bodies resources/data/bodies_planet.txt).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* f = dynamic_cast<const PlanetField*>(p->field.get());
                    if (!f) {
                        printf("bench: --test-mantle needs a planet as the first body (e.g. --bodies resources/data/bodies_planet.txt)\n");
                    } else {
                        std::mt19937 rng(1234);
                        std::normal_distribution<float> gauss(0.0f, 1.0f);
                        std::uniform_real_distribution<float> unit(0.0f, 1.0f);
                        auto randomDir = [&]() { return Vector3Normalize({ gauss(rng), gauss(rng), gauss(rng) }); };
                        // 1. Height band: march down each direction to the first solid point (gap-sized steps,
                        // like the raycast), bisect, and compare with the field's stated bounds.
                        const int N = 4000;
                        float lo = 1e9f, hi = -1e9f;
                        int ocean = 0, outside = 0;
                        for (int i = 0; i < N; ++i) {
                            const Vector3 dir = randomDir();
                            float r = f->maxRadius() + 10.0f, dens = f->baseDensity(Vector3Scale(dir, r)), prev = r;
                            while (dens < ISO_LEVEL && r > f->mantleRadius) {
                                prev = r;
                                r -= std::clamp((ISO_LEVEL - dens) / DENSITY_GRADIENT * 0.5f, 0.25f, 64.0f);
                                dens = f->baseDensity(Vector3Scale(dir, r));
                            }
                            float a = r, b = prev; // a solid, b air
                            for (int k = 0; k < 20; ++k) {
                                const float m = 0.5f * (a + b);
                                if (f->baseDensity(Vector3Scale(dir, m)) >= ISO_LEVEL) a = m; else b = m;
                            }
                            const float e = a - f->radius;
                            lo = std::min(lo, e); hi = std::max(hi, e);
                            if (e < 0.0f) ++ocean;
                            if (e < f->minElevation - 3.0f || e > f->maxElevation + 3.0f) ++outside;
                        }
                        printf("bench: planet %s: %d surface samples, elevation %.0f..%.0f m (band %.0f..%.0f), below sea level %.0f%%, outside band %d\n",
                               p->name.c_str(), N, lo, hi, f->minElevation, f->maxElevation, 100.0f * ocean / N, outside);
                        // 2. Mantle: dig just above it; the crust must change and the mantle must not.
                        const Vector3 dir = randomDir();
                        const Vector3 inMantle = Vector3Scale(dir, f->mantleRadius - 4.0f), inCrust = Vector3Scale(dir, f->mantleRadius + 12.0f);
                        const float mantleBefore = p->densityAt(inMantle), crustBefore = p->densityAt(inCrust);
                        const auto savedEdits = p->edits;
                        ApplyTerrainEdit(*p, TerrainEdit{ Vector3Scale(dir, f->mantleRadius + 4.0f), 12.0f, -20.0f });
                        const float mantleAfter = p->densityAt(inMantle), crustAfter = p->densityAt(inCrust);
                        p->edits = savedEdits; // leave the world as it was
                        const bool mantleOk = mantleAfter == mantleBefore && crustAfter < crustBefore - 1.0f && !f->editable(inMantle) && f->editable(inCrust);
                        printf("bench: mantle %.0f m below sea level: dig there changed the crust %.2f -> %.2f and the mantle %.2f -> %.2f: %s\n",
                               f->radius - f->mantleRadius, crustBefore, crustAfter, mantleBefore, mantleAfter, mantleOk ? "refused" : "NOT REFUSED");
                        // 3. Ore by depth band, and caves.
                        const float bands[4] = { 8.0f, 150.0f, 500.0f, f->radius - f->mantleRadius - 20.0f };
                        float orePct[3];
                        int caveSamples = 0, solidSamples = 0;
                        for (int band = 0; band < 3; ++band) {
                            int ore = 0, solid = 0;
                            for (int i = 0; i < 20000; ++i) {
                                const Vector3 d = randomDir();
                                const float depth = bands[band] + unit(rng) * (bands[band + 1] - bands[band]);
                                const Vector3 pt = Vector3Scale(d, f->radius + f->elevation(d) - depth);
                                if (!f->editable(pt)) continue; // below the mantle (deep ocean floors)
                                const float base = f->baseDensity(pt);
                                if (base < ISO_LEVEL) { ++caveSamples; continue; }
                                ++solid; ++solidSamples;
                                if (f->material(pt, base) == VoxelMaterial::Ore) ++ore;
                            }
                            orePct[band] = 100.0f * ore / std::max(1, solid);
                        }
                        const float cavePct = 100.0f * caveSamples / std::max(1, caveSamples + solidSamples);
                        const bool oreOk = orePct[0] > 0.0f && orePct[0] < orePct[1] && orePct[1] < orePct[2];
                        printf("bench: ore by depth: 8-150 m %.1f%%, 150-500 m %.1f%%, 500 m-mantle %.1f%%; caves %.1f%% of underground: %s\n",
                               orePct[0], orePct[1], orePct[2], cavePct, oreOk && cavePct > 0.0f ? "ok" : "BAD");
                        printf("bench: planet terrain %s\n", outside == 0 && mantleOk && oreOk && cavePct > 0.0f ? "PASS" : "FAIL");
                    }
                    benchMantle = false;
                } else if (benchPrecision && WorldIdle()) {
                    // Float vs double placement on body 0, standing on its surface (run it with a large
                    // body, e.g. --bodies with a 60 km planet): chunk translations, texture coordinates
                    // near the camera, and a body-anchor round trip, each against an exact double result.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (precisionPhase == 0) {
                        // Stand on the surface (the bench camera starts in fly mode off to the side), then
                        // wait for the terrain around the new viewpoint to stream in.
                        player.spawnOn(world, *p, Vector3Normalize({ 0.8f, 0.35f, -0.5f })); // off the spin axis and the grid axes
                        precisionStart = GetTime();
                        precisionPhase = 1;
                    } else if (GetTime() - precisionStart >= 1.5) {
                    auto ulpMm = [](float x) { x = fabsf(x); return (nextafterf(x, INFINITY) - x) * 1000.0f; };
                    const DVec3 cam = world.cameraWorld();
                    const double P = 1.0 / world.texScale;
                    double chunkOld = 0, chunkNew = 0;
                    float texOld = 0, texNew = 0;
                    int nearChunks = 0;
                    for (Chunk* c : p->drawList) {
                        const Vector3 o = c->key.origin();
                        const DVec3 exact = p->toWorldD(DVec3(o)) - world.renderOrigin;
                        const Matrix mo = MatrixMultiply(MatrixTranslate(o.x, o.y, o.z), p->localTransform());
                        const Matrix mn = p->chunkModel(*c, world.renderOrigin);
                        const float size = (float)c->key.size() * (CHUNK_SIZE + 1) / CHUNK_SIZE; // metres, incl. overlap cell
                        // Only chunks up close count: far away, a few mm of float rounding is invisible.
                        if (c->key.level > 1 || (exact + world.renderOrigin - cam).length() > 200.0 + size) continue;
                        chunkOld = std::max(chunkOld, (DVec3(mo.m12, mo.m13, mo.m14) - exact).length());
                        chunkNew = std::max(chunkNew, (DVec3(mn.m12, mn.m13, mn.m14) - exact).length());
                        ++nearChunks;
                        for (int a = 0; a < 3; ++a) {
                            const float oa = a == 0 ? o.x : a == 1 ? o.y : o.z;
                            const double wrapped = fmod(fmod((double)oa, P) + P, P);
                            texOld = std::max(texOld, std::max(ulpMm(oa), ulpMm(oa + size)));
                            texNew = std::max(texNew, ulpMm((float)wrapped + size));
                        }
                    }
                    const DVec3 abs = player.position;
                    const Vector3 feet = (abs - world.renderOrigin).toFloat();
                    const double anchorOld = ((world.renderOrigin + p->toWorld(p->toLocal(feet))) - abs).length();
                    const double anchorNew = (p->toWorldD(p->toLocalD(abs)) - abs).length();
                    printf("bench: precision on %s, %.0f m from its centre (%zu chunks drawn, %d fine chunks within 200 m):\n", p->name.c_str(),
                           (abs - p->worldPos).length(), p->drawList.size(), nearChunks);
                    printf("bench:   chunk placement error   float %.3f mm   double %.4f mm\n", chunkOld * 1000.0, chunkNew * 1000.0);
                    printf("bench:   texture coordinate step float %.3f mm   wrapped %.4f mm\n", texOld, texNew);
                    printf("bench:   anchor round trip error float %.3f mm   double %.6f mm\n", anchorOld * 1000.0, anchorNew * 1000.0);
                    benchPrecision = false;
                    }
                } else if (benchGravity && WorldIdle()) {
                    // Scripted gravity checks: the field itself (surface value, half at half depth, a
                    // quarter at 2R, no jumps between bodies), then an unpiloted ship dropped from 25 m
                    // must land and park, and a piloted ship must hover in place.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = gravStart < 0.0 ? 0.0 : GetTime() - gravStart;
                    // Height of the hull centre above the ground straight "down" (along gravity).
                    auto altitude = [&]() {
                        const Vector3 sp = ship.renderPos(world);
                        TerrainHit hit;
                        const Vector3 down = Vector3Normalize(GravityAt(world, sp));
                        return RaycastTerrain(world, Ray{ sp, down }, 500.0f, hit) ? hit.distance : -1.0f;
                    };
                    auto placeAbove = [&](float height) {
                        TerrainHit hit;
                        const float R = p->radius * p->scale;
                        if (RaycastTerrain(world, Ray{ Vector3Add(p->position, { 0, R + 80.0f, 0 }), { 0, -1, 0 } }, R + 100.0f, hit)) {
                            ship.parkNear(world, hit.point, hit.normal);
                            ship.position = world.renderOrigin + Vector3Add(hit.point, Vector3Scale(hit.normal, ship.hullRadius + height));
                        }
                        ship.velocity = { 0, 0, 0 };
                    };
                    if (gravStart < 0.0) {
                        const float R = p->gravityRadius, g0 = p->surfaceGravity; // gravity's reference radius (mean surface)
                        const Vector3 c = p->position;
                        const float gs = Vector3Length(BodyGravityAt(*p, Vector3Add(c, { R, 0, 0 })));
                        const float gh = Vector3Length(BodyGravityAt(*p, Vector3Add(c, { 0, 0.5f * R, 0 }))) / g0;
                        const float g2 = Vector3Length(BodyGravityAt(*p, Vector3Add(c, { 0, 0, -2.0f * R }))) / g0;
                        bool pass = fabsf(gs - g0) < 1e-3f && fabsf(gh - 0.5f) < 1e-3f && fabsf(g2 - 0.25f) < 1e-3f;
                        printf("bench: gravity field: surface %.3f m/s^2 (body %.3f), half depth %.3f g0, 2R %.3f g0\n", gs, g0, gh, g2);
                        // Continuity: walk the line to the nearest other body; the summed field may change
                        // no faster than its steepest gradient (2 g0 / R, just outside a surface) allows.
                        Planetoid* other = nullptr;
                        for (size_t i = 1; i < world.objects.size(); ++i)
                            if (auto* q = dynamic_cast<Planetoid*>(world.objects[i].get()))
                                if (!other || Vector3Distance(q->position, c) < Vector3Distance(other->position, c)) other = q;
                        if (other) {
                            const int N = 20000;
                            const float len = Vector3Distance(c, other->position), step = len / N;
                            float maxGrad = 0.0f;
                            for (auto& obj : world.objects)
                                if (auto* q = dynamic_cast<Planetoid*>(obj.get())) maxGrad = std::max(maxGrad, 2.0f * q->surfaceGravity / q->gravityRadius);
                            float maxJump = 0.0f, weakest = 1e9f, weakestAt = 0.0f;
                            Vector3 prev = GravityAt(world, c);
                            for (int i = 1; i <= N; ++i) {
                                const Vector3 g = GravityAt(world, Vector3Lerp(c, other->position, (float)i / N));
                                maxJump = std::max(maxJump, Vector3Distance(g, prev));
                                if (Vector3Length(g) < weakest) { weakest = Vector3Length(g); weakestAt = (float)i / N; }
                                prev = g;
                            }
                            const float bound = 1.1f * maxGrad * step;
                            pass = pass && maxJump <= bound;
                            printf("bench: gravity continuity %s -> %s (%.0f m): max change %.4f m/s^2 per %.3f m step (bound %.4f), weakest %.4f m/s^2 at %.2f of the way\n",
                                   p->name.c_str(), other->name.c_str(), len, maxJump, step, bound, weakest, weakestAt);
                        }
                        printf("bench: gravity field %s\n", pass ? "PASS" : "FAIL");
                        // Phase 1: drop an unpiloted ship from 25 m.
                        ship.piloted = false;
                        player.flyMode = true; // keep the player out of the way (fly mode: no physics)
                        placeAbove(25.0f);
                        gravAlt0 = altitude();
                        gravMaxImpact = 0.0f;
                        gravPhase = 1;
                        gravStart = GetTime();
                    } else if (gravPhase == 1) {
                        gravMaxImpact = std::max(gravMaxImpact, ship.lastImpact);
                        if (t > 10.0) {
                            printf("bench: ship drop from %.1f m (g %.2f): hull centre now %.2f m above ground, parked %s, impact %.1f m/s\n",
                                   gravAlt0, Vector3Length(GravityAt(world, ship.renderPos(world))), altitude(), ship.anchorBody ? "yes" : "NO", gravMaxImpact);
                            // Phase 2: hover a piloted ship 25 m up with no input.
                            placeAbove(25.0f);
                            ship.piloted = true;
                            gravAlt0 = altitude();
                            gravEnergy0 = ship.energy;
                            gravPhase = 2;
                            gravStart = GetTime();
                        }
                    } else if (gravPhase == 2 && t > 5.0) {
                        printf("bench: ship hover 5 s at %.1f m: altitude change %+.2f m, energy %.1f -> %.1f (hovering %s)\n",
                               gravAlt0, altitude() - gravAlt0, gravEnergy0, ship.energy, ship.hovering ? "yes" : "NO");
                        ship.piloted = false;
                        benchGravity = false;
                    }
                } else if (benchPark && WorldIdle()) {
                    // Scripted park: land the ship on the equator of a body spinning 10 deg/s; it must stay put
                    // relative to the surface (position and heading) while parked.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = parkStart < 0.0 ? 0.0 : GetTime() - parkStart;
                    if (parkStart < 0.0) {
                        p->spin = { 0.0f, 10.0f, 0.0f };
                        TerrainHit hit;
                        const float R = p->radius * p->scale;
                        if (RaycastTerrain(world, Ray{ Vector3Add(p->position, { R + 60.0f, 0, 0 }), { -1, 0, 0 } }, R + 80.0f, hit))
                            ship.parkNear(world, hit.point, hit.normal);
                        ship.piloted = false;
                        parkStart = GetTime();
                        parkPhase = 0;
                    } else if (parkPhase == 0 && t > 1.0) {
                        parkLocal0 = p->toLocal(ship.renderPos(world));
                        parkFwd0 = p->dirToLocal(ship.forward());
                        parkAnchored = ship.anchorBody != nullptr;
                        parkPhase = 1;
                    } else if (parkPhase == 1 && t > 5.0) {
                        const float moved = Vector3Distance(parkLocal0, p->toLocal(ship.renderPos(world)));
                        const float turned = Vector3Angle(parkFwd0, p->dirToLocal(ship.forward())) * RAD2DEG;
                        printf("bench: park test (10 deg/s, body turned %.0f deg in 4 s): anchored %s, ship moved %.3f m and turned %.2f deg relative to the surface\n",
                               4.0 * 10.0, parkAnchored ? "yes" : "NO", moved, turned);
                        benchPark = false;
                    }
                } else if (benchRide && WorldIdle()) {
                    // Scripted ride: stand on the equator of a body spinning 10 deg/s. The view must not drift
                    // against the surface, and a jump must land where it took off (momentum is kept).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const double t = rideStart < 0.0 ? 0.0 : GetTime() - rideStart;
                    if (rideStart < 0.0) {
                        p->spin = { 0.0f, rideSpin, 0.0f };
                        player.spawnOn(world, *p, { 1.0f, 0.0f, 0.0f });
                        rideStart = GetTime();
                        ridePhase = 0;
                    } else if (ridePhase == 0 && t > 1.5) {
                        rideFwd0 = p->dirToLocal(player.forward);
                        ridePos0 = p->toLocal(feet);
                        ridePhase = 1;
                    } else if (ridePhase == 1 && t <= 4.5) {
                        ++rideFrames; if (player.grounded) ++rideGroundedFrames;
                    } else if (ridePhase == 1 && t > 4.5) {
                        rideViewDrift = Vector3Angle(rideFwd0, p->dirToLocal(player.forward)) * RAD2DEG;
                        ridePosDrift = Vector3Distance(ridePos0, p->toLocal(feet));
                        rideTakeoff = p->toLocal(feet);
                        rideSurfaceSpeed = Vector3Length(p->surfaceVelocity(feet));
                        player.autoJump = true;
                        ridePhase = 2;
                    } else if (ridePhase == 2 && !player.grounded) {
                        ridePhase = 3; // airborne
                    } else if ((ridePhase == 3 && player.grounded) || t > 12.0) {
                        printf("bench: ride test (%.0f deg/s, surface %.1f m/s): standing 3 s (grounded %.0f%%) -> view drift %.2f deg, position drift %.2f m; jump landed %.2f m from takeoff (body frame)%s\n",
                               rideSpin, rideSurfaceSpeed, 100.0f * rideGroundedFrames / std::max(1, rideFrames), rideViewDrift, ridePosDrift, Vector3Distance(rideTakeoff, p->toLocal(feet)), ridePhase == 3 ? "" : " [never landed]");
                        benchRide = false;
                    }
                } else if (benchFlySeconds > 0.0 && benchMineStrokes == 0 && benchWalkSeconds <= 0.0) {
                    // Scripted flight: board the ship and autopilot to the nearest other asteroid.
                    if (flyStart < 0.0) {
                        Planetoid* target = nullptr;
                        float best = 1e30f;
                        for (size_t i = 1; i < world.objects.size(); ++i)
                            if (auto* p = dynamic_cast<Planetoid*>(world.objects[i].get())) {
                                const float d = Vector3Distance(p->position, ship.renderPos(world));
                                if (d < best) { best = d; target = p; }
                            }
                        if (!target) { printf("bench: --test-fly: no second asteroid\n"); benchFlySeconds = 0.0; }
                        else {
                            player.flyMode = false;
                            ship.piloted = true;
                            ship.autopilot = true;
                            // Aim for a point 1.6 radii out on the near side of the target.
                            const Vector3 toShip = Vector3Normalize(Vector3Subtract(ship.renderPos(world), target->position));
                            ship.autopilotTarget = world.renderOrigin + Vector3Add(target->position, Vector3Scale(toShip, target->radius * 1.6f));
                            flyStart = GetTime();
                            flyFrom = ship.position;
                            flyEnergy = ship.energy;
                            printf("bench: flying to %s, %.0f m away\n", target->name.c_str(), best);
                        }
                    } else {
                        const double remaining = (ship.autopilotTarget - ship.position).length();
                        if (remaining < 40.0 || GetTime() - flyStart > benchFlySeconds) {
                            ship.autopilot = false;
                            printf("bench: flew %.0f m in %.1f s, %.0f m from target, energy %.0f -> %.0f, max speed %.0f m/s, visited another asteroid: %s\n",
                                   (ship.position - flyFrom).length(), GetTime() - flyStart, remaining, flyEnergy, ship.energy, flyMaxSpeed,
                                   remaining < 40.0 ? "yes" : "no");
                            benchFlySeconds = 0.0;
                        }
                        flyMaxSpeed = std::max(flyMaxSpeed, Vector3Length(ship.velocity));
                    }
                } else if ((benchMineStrokes > 0 || minePhase == 2) && benchWalkSeconds <= 0.0) {
                    // Scripted mining: look at the ground, dig N strokes through the real tool path,
                    // craft everything possible, then build a few strokes back.
                    if (!mineStarted) {
                        player.flyMode = false;
                        player.forward = Vector3Normalize(Vector3Add(Vector3Scale(player.up, -0.8f), Vector3Scale(player.forward, 0.6f)));
                        mineStarted = true;
                        minePhase = 1;
                    } else if (minePhase == 1) {
                        const bool dug = TerrainTool_Stroke(true);
                        if (dug) { minedStrokes++; benchMineStrokes--; }
                        if (!dug || benchMineStrokes == 0) {
                            mineStoppedEarly = !dug;
                            int crafted = 0;
                            for (bool any = true; any && crafted < 100;) {
                                any = false;
                                for (const Recipe& r : itemDB.recipes()) if (inventory.craft(r)) { any = true; crafted++; }
                            }
                            mineCrafted = crafted;
                            // Look at nearby ground (not into the hole) for the build strokes next frame.
                            const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                            player.forward = Vector3Normalize(Vector3Add(Vector3Scale(flat, -0.6f), Vector3Scale(player.up, -0.8f)));
                            minePhase = 2;
                        }
                    } else {
                        const int fillersBefore = inventory.count("stone") + inventory.count("regolith") + inventory.count("concrete");
                        int built = 0;
                        for (int b = 0; b < 3; ++b) built += TerrainTool_Stroke(false) ? 1 : 0;
                        std::string inv;
                        for (const ItemDef& d : itemDB.items()) if (inventory.count(d.id)) inv += " " + d.id + "=" + std::to_string(inventory.count(d.id));
                        printf("bench: mined %d strokes%s, crafted %d times, built %d strokes (fillers %d -> %d); inventory:%s\n",
                               minedStrokes, mineStoppedEarly ? " (hole left arm's reach)" : "", mineCrafted, built, fillersBefore,
                               inventory.count("stone") + inventory.count("regolith") + inventory.count("concrete"), inv.c_str());
                        benchMineStrokes = 0;
                        minePhase = 3;
                    }
                } else if (WorldIdle() && benchWalkSeconds > 0.0) {
                    // Scripted walk: spawn on the surface and walk forward; checks gravity + collision.
                    if (walkStart < 0.0) {
                        player.autoWalk = !benchStandOnly;
                        Planetoid* start = dynamic_cast<Planetoid*>(world.objects[0].get());
                        player.spawnOn(world, *start);

                        walkStart = GetTime();
                        walkFrom = player.position;
                        walkMinDensity = 1e9f; walkMaxDensity = -1e9f; walkGroundedFrames = walkFrames = 0;
                    }
                    Planetoid* start = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const float d = start->densityAt(start->toLocal(Vector3Add(feet, Vector3Scale(player.up, 0.9f)))); // body centre
                    walkMinDensity = std::min(walkMinDensity, d);
                    walkMaxDensity = std::max(walkMaxDensity, d);
                    walkFrames++;
                    walkGroundedFrames += player.grounded ? 1 : 0;
                    if (GetTime() - walkStart > benchWalkSeconds) {
                        player.autoWalk = false;
                        const double dist = (player.position - walkFrom).length();
                        const float alt = Vector3Distance(feet, start->position);
                        printf("bench: walked %.1f m in %.1f s, grounded %.0f%% of frames, density at body %.2f..%.2f (solid >= %.2f), %.1f m from centre\n",
                               dist, benchWalkSeconds, 100.0 * walkGroundedFrames / std::max(1, walkFrames), walkMinDensity, walkMaxDensity, ISO_LEVEL, alt);
                        benchWalkSeconds = 0.0;
                    }
                } else if ((WorldIdle() && benchDigStrokes == 0 && !benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && minePhase != 2 && benchFlySeconds <= 0.0 && !benchSuffocate && benchOrbitSeconds <= 0.0 && !benchRide && !benchPark && !benchGravity && !benchPrecision && !benchMantle) || GetTime() > 180.0) {
                    benchStart = GetTime();
                    printf("bench: world ready after %.2f s (%d chunks, avg build %.2f ms of which sampling %.2f ms, %d worker threads)\n",
                        GetWorldStats().lastGenerationSeconds, GetWorldStats().chunks, GetWorldStats().avgBuildMs, GetWorldStats().avgSampleMs, GetWorldStats().workerThreads);
                }
            }
            double t = (benchStart < 0.0) ? 0.0 : GetTime() - benchStart;
            if (benchStart >= 0.0 && t > benchWarmup) {
                benchFrames++;
                ++benchTimedFrames;
                benchMeasuring = true;
                if (GetFrameTime() > benchMaxFrame) benchMaxFrame = GetFrameTime();
            }
            if (t > benchWarmup + benchSeconds) {
                double avgMs = 1000.0 * benchSeconds / (double)benchFrames;
                const char* line = TextFormat("bench: %ld frames in %.1fs | avg %.1f FPS | avg %.2f ms | worst %.2f ms | world update %.2f ms, draw %.2f ms (CPU)\n",
                    benchFrames, benchSeconds, 1000.0 / avgMs, avgMs, benchMaxFrame * 1000.0, benchUpdateMs / std::max(1L, benchTimedFrames), benchDrawMs / std::max(1L, benchTimedFrames));
                printf("%s", line);
                SaveFileText("bench.txt", (char*)line);
                if (crackCheck) {
                    // Holes between chunks show the magenta clear colour (look down at terrain that fills the view).
                    Image shot = LoadImageFromScreen();
                    Color* px = LoadImageColors(shot);
                    long holes = 0;
                    const long total = (long)shot.width * shot.height;
                    for (long k = 0; k < total; ++k) if (px[k].r > 240 && px[k].g < 16 && px[k].b > 240) ++holes;
                    printf("bench: crack check: %ld of %ld pixels show through the terrain (%.4f%%)\n", holes, total, 100.0 * holes / total);
                    UnloadImageColors(px);
                    UnloadImage(shot);
                }
                TakeScreenshot("bench.png");
                break;
            }
        } else if (!IsWindowFocused()) {
            ShowCursor();
            BeginDrawing();
            EndDrawing();
            continue;
        } else if (gameState == 1) {
            if (DebugUI_IsOpen() || GameUI_InventoryOpen()) ShowCursor(); else HideCursor();
        }
        if (!bench && IsKeyPressed(KEY_F1)) {
            DebugUI_Toggle();
            if (!DebugUI_IsOpen()) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
        }
        if (!bench && gameState == 1 && IsKeyPressed(KEY_TAB)) {
            GameUI_ToggleInventory();
            if (!GameUI_InventoryOpen()) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
        }
        const bool uiOpen = DebugUI_IsOpen() || GameUI_InventoryOpen();
        const bool gameInput = !bench && !uiOpen && gameState == 1;
        world.cameraInputEnabled = player.flyMode && !uiOpen && !ship.piloted;  // free-fly camera only in fly mode
        const bool alive = !survival.dead;
        if (gameState == 1 && gameInput && alive && IsKeyPressed(KEY_E)) {
            const Vector3 feet = (player.position - world.renderOrigin).toFloat();
            if (ship.piloted) {
                // Step out beside the ship.
                ship.piloted = false;
                player.position = world.renderOrigin + Vector3Add(ship.renderPos(world), Vector3Add(Vector3Scale(ship.right(), 3.5f), Vector3Scale(ship.up(), 0.5f)));
                player.velocity = ship.velocity;
                player.flyMode = false;
                Audio_Play(Sfx::Build, 0.4f, 1.4f);
            } else if (ship.playerInRange(world, feet, 5.0f)) {
                ship.piloted = true;
                player.flyMode = false;
                Audio_Play(Sfx::Build, 0.4f, 1.2f);
            }
        }
        advanceBodies(world, gameState == 1 ? GetFrameTime() : 0.0f); // bodies only spin while playing; before riders update
        if (gameState == 1) {
            ship.lastImpact = 0.0f;
            if (ship.piloted) {
                ship.update(world, gameInput && alive, GetFrameTime());
                ship.applyChaseCamera(world);
                player.position = ship.position; // the pilot rides inside
                player.velocity = ship.velocity;
            } else {
                ship.update(world, false, GetFrameTime());
                if (alive) player.update(world, gameInput, GetFrameTime());
                Survival_OnLanding(player.impactSpeed);
            }
            if (ship.piloted && ship.lastImpact > 12.0f) survival.health -= (ship.lastImpact - 12.0f) * 3.0f; // crash damage
            Survival_Update(world, GetFrameTime(), !ship.piloted && !player.grounded && gameInput && player.jetpackFuel > 0.0f && (IsKeyDown(KEY_SPACE) || IsKeyDown(KEY_LEFT_CONTROL) || IsKeyDown(KEY_W) || IsKeyDown(KEY_A) || IsKeyDown(KEY_S) || IsKeyDown(KEY_D)), gameInput);
            Objectives_Update(world);
            if (bench && lookDown && !ship.piloted) {
                const Vector3 t = Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)));
                player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(t), 0.8f), player.up)); // ~50 deg down
            }
            if (!ship.piloted) player.applyToCamera(world);
        }

        world.updateFloatingOrigin(); // keep the camera near render-space 0 (float precision)
        { const double w0 = GetTime(); worldUpdate(world, 3.0f); if (benchMeasuring) benchUpdateMs += (GetTime() - w0) * 1000.0; } // LOD + streaming (bodies were spun at the start of the frame)

        // Slot 1: headlamp (suit lamp from the eye, or the ship's nose light while piloting).
        // Slot 2: soft glow around the player while the inventory is open.
        if (sunOverride) world.sunDirection = sunArg;
        if (crackCheck) world.drawSkybox = false;
        while (world.lights.size() < 3) world.lights.push_back(Light{});
        {
            Light& lamp = world.lights[1];
            lamp.type = LIGHT_SPOT;
            lamp.enabled = gameState == 1 && survival.headlamp && !survival.dead;
            if (lamp.enabled) {
                const float range = ship.piloted ? 60.0f : 35.0f;
                Vector3 from, dir;
                if (ship.piloted) { from = Vector3Add(ship.renderPos(world), Vector3Scale(ship.forward(), 2.6f)); dir = ship.forward(); }
                else {
                    dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
                    from = Vector3Add(world.camera.position, Vector3Scale(world.camera.up, -0.15f)); // brow, not the eye
                }
                lamp.position = from;
                lamp.target = Vector3Add(from, Vector3Scale(dir, range));
                lamp.color = { 255, 244, 225, 110 };  // alpha = intensity
            }
            Light& glow = world.lights[2];
            glow.type = LIGHT_POINT;
            glow.enabled = gameState == 1 && GameUI_InventoryOpen() && !ship.piloted && !survival.dead;
            if (glow.enabled) {
                glow.position = Vector3Add(world.camera.position, Vector3Scale(world.camera.up, -0.3f));
                glow.target = { 6.0f, 0.0f, 0.0f };   // point light: range in target.x
                glow.color = { 255, 205, 150, 32 };   // warm, dim
            }
        }
        // --spin-test: camera locked to the spinning first asteroid; texture must turn with the rock.
        if (bench && spinTest && benchStart >= 0.0) {
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const double t = GetTime() - benchStart;
            if (spinPhase == 0) {
                p->spin = { 0.0f, spinRate, 0.0f };
                world.sunDirection = { 0.0f, 1.0f, 0.0f };      // sun on the spin axis: shading stays constant
                spinCamLocal = p->toLocal(world.camera.position);
                spinTgtLocal = p->toLocal(world.camera.target);
                spinPhase = 1;
            }
            world.camera.position = p->toWorld(spinCamLocal);    // ride along with the body
            world.camera.target = p->toWorld(spinTgtLocal);
            if (spinPhase == 1 && t > 1.0) {
                TakeScreenshot("spin_a.png");
                spinPhase = 2;
            }
            if (spinPhase == 2 && t > 3.0) {
                TakeScreenshot("spin_b.png");
                Image a = LoadImage("spin_a.png"), b = LoadImage("spin_b.png");
                ImageFormat(&a, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
                ImageFormat(&b, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
                double diff = 0.0;
                long n = 0;
                const unsigned char* pa = (const unsigned char*)a.data;
                const unsigned char* pb = (const unsigned char*)b.data;
                for (int y = a.height * 35 / 100; y < a.height * 65 / 100; ++y)
                    for (int x = a.width * 35 / 100; x < a.width * 65 / 100; ++x, ++n)
                        for (int c = 0; c < 3; ++c) diff += abs((int)pa[(y * a.width + x) * 4 + c] - (int)pb[(y * a.width + x) * 4 + c]);
                printf("bench: spin test: body turned %.0f deg between shots; mean abs pixel diff in view centre = %.2f / 255\n",
                       2.0 * spinRate, diff / (n * 3.0));
                {   // amplified difference image for inspection
                    Image d = ImageCopy(a);
                    unsigned char* pd = (unsigned char*)d.data;
                    for (int i = 0; i < a.width * a.height; ++i)
                        for (int c = 0; c < 3; ++c) pd[i * 4 + c] = (unsigned char)std::min(255, 4 * abs((int)pa[i * 4 + c] - (int)pb[i * 4 + c]));
                    ExportImage(d, "spin_diff.png");
                    UnloadImage(d);
                }
                UnloadImage(a);
                UnloadImage(b);
                spinPhase = 3;
            }
        }
        if (gameState == 1) TerrainTool_Update(world, gameInput && alive, !ship.piloted); // F5/F9 work aboard too

        switch(gameState) {
            case 0: // Main Menu
                if (IsKeyPressed(KEY_ENTER)) {
                    printf("Start Game\n");
                    HideCursor();
                    gameState = 1; // Start game
                    mainMenu.isActive = false; // Deactivate main menu
                    world.isActive = true; // Activate world scene
                }
                if (IsKeyPressed(KEY_ESCAPE)) {
                    printf("Exit game from main menu\n");
                    windowKeepAlive = false; // Exit game
                }
                break;
            case 1: // Game
                if (IsKeyPressed(KEY_ESCAPE)) {
                    ShowCursor();
                    printf("Game paused\n");
                    gameState = 2; // Pause game
                }
                break;
            case 2: // Pause Menu
                if (IsKeyPressed(KEY_R)) {
                    printf("Game unpaused\n");
                    HideCursor();
                    gameState = 1; // Resume game
                }
                if (IsKeyPressed(KEY_ESCAPE)) {
                    printf("Exit game from pause\n");
                    windowKeepAlive = false; // Exit game
                }
                break;
        }

        

        BeginDrawing();
            ClearBackground(crackCheck ? MAGENTA : BLACK); // --crack-check: holes in the terrain show magenta
            if (gameState == 0) {
                mainMenu.drawScene(gameState);
                mainMenu.drawUI(gameState);
            } else if (gameState == 1) {
                { const double d0 = GetTime(); world.drawScene(gameState); if (benchMeasuring) benchDrawMs += (GetTime() - d0) * 1000.0; }
                world.drawUI(gameState);
                if (!ship.piloted) TerrainTool_DrawHUD();
                GameUI_DrawHUD(world);
                Survival_DrawHUD();
                Objectives_DrawHUD();
                if (!ship.piloted && ship.playerInRange(world, (player.position - world.renderOrigin).toFloat(), 5.0f))
                    DrawText("[E] Board ship", GetScreenWidth() / 2 - 60, GetScreenHeight() / 2 + 40, 20, SKYBLUE);
                if (ship.piloted) DrawText("SHIP: WASD thrust  Space/Ctrl up/down  mouse steer  Z/C roll  Shift boost  X brake  L lights  [2] power cell  E exit",
                                           10, GetScreenHeight() - 30, 18, Fade(SKYBLUE, 0.9f));
                if (!bench && !uiOpen) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
            } else if (gameState == 2) {
                world.drawScene(gameState);
                world.drawUI(gameState);
                // Blur background
                DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), Fade(BLACK, 0.7f));
                // Draw pause menu text
                DrawText("Game Paused", GetScreenWidth()/2 - MeasureText("Game Paused", 20)/2, GetScreenHeight()/2 - 10, 20, WHITE);
                DrawText("Press R to resume", GetScreenWidth()/2 - MeasureText("Press R to resume", 20)/2, GetScreenHeight()/2 + 10, 20, WHITE);
                DrawText("Press ESC to exit", GetScreenWidth()/2 - MeasureText("Press ESC to exit", 20)/2, GetScreenHeight()/2 + 30, 20, WHITE);
            }
            if (gameState != 0) {
                rlImGuiBegin();
                DebugUI_Draw(world);
                GameUI_Draw(world);
                rlImGuiEnd();
            }
        EndDrawing();
        PROFILE_FRAME();
    }

    }   // scenes destroyed here, while the GL context still exists
    DebugUI_Shutdown();
    Audio_Shutdown();
    CloseWindow();
}


