#include <raylib.h>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <array>
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
#include "camera.h"
#include "survival.h"
#include "telescope.h"
#include "buildtool.h"
#include "progression_test.h"
#include "menus.h"
#include "hotbar.h"
#include "ring_test.h"
#include "rocket.h"
#include "rocket_test.h"
#include "factory.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif

using namespace std;

Logger logger("log.txt");

int gameState = 0; // 0: Main Menu, 1: Game, 2: Pause Menu

// Where to start on the first body: planets pick temperate grassland (world direction, from the
// body's seed); asteroids keep the old "top" spot.
// --telescope Z [--telescope-target name]: bench view through the telescope at zoom Z, aimed at a
// body (default: the moon) wherever it is in the sky.
static float benchTelescope = 0.0f;
static std::string benchTelescopeTarget = "luna";
static void AimAtTelescopeTarget(Scene& world) {
    for (auto& obj : world.objects) {
        auto* p = dynamic_cast<Planetoid*>(obj.get());
        if (!p || p->name != benchTelescopeTarget) continue;
        const Vector3 eye = Vector3Add((player.position - world.renderOrigin).toFloat(), Vector3Scale(player.up, player.eyeHeight));
        player.forward = Vector3Normalize(Vector3Subtract(p->position, eye));
        static bool told = false;
        if (!told) {
            told = true;
            printf("bench: telescope on %s, %.0f km away, %.1f deg above the horizon\n", p->name.c_str(),
                   Vector3Distance(eye, p->position) / 1000.0f, asinf(Vector3DotProduct(player.forward, player.up)) * RAD2DEG);
        }
        return;
    }
}

static Vector3 StartSpawnDirLocal(Planetoid& body) {
    if (auto* field = dynamic_cast<const PlanetField*>(body.field.get())) return field->findSpawnDir(body.seed);
    return { 0, 1, 0 };
}

static Vector3 StartSpawnDir(Planetoid& body) {
    if (auto* field = dynamic_cast<const PlanetField*>(body.field.get())) return body.dirToWorld(field->findSpawnDir(body.seed));
    return { 0, 1, 0 };
}

// A new game on the world as generated: the player at the start, the ship parked ahead, a
// clearing for both, the starting kit. `reset` also clears what a previous game left behind
// (inventory, quest tier, survival state).
static void StartNewWorld(Scene& world, bool reset) {
    if (reset) {
        inventory.clear();
        QuestBook_SetTier(0);
        survival = SurvivalState{};
        player.velocity = { 0, 0, 0 };
        player.flyMode = false;
        ship.piloted = false;
        Rocket_Reset();
    }
    Planetoid* start = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    if (start) {
        player.spawnOn(world, *start, StartSpawnDir(*start));
        // Park the ship a few metres ahead of the player.
        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
        ship.parkNear(world, Vector3Add(feet, Vector3Scale(player.forward, 8.0f)), player.up);
        // A clearing at the start, so the parked ship and the player don't stand in a tree.
        if (start->trees) start->trees->addClearing(start->toLocal(feet), 14.0f);
    }
    if (inventory.count("telescope") == 0) inventory.add("telescope", 1); // starting kit (a save restores its own inventory)
}

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
    float benchFlySpeed = 0.0f;    // --bench-fly V: the camera flies over body 0 at V m/s through the bench (LOD streaming)
    float benchFlyAlt = 1e9f;      // --bench-fly-alt A: fly at A m relative to sea level instead (negative: underwater, looking up)
    bool hasSpawnDir = false;      // --spawn-dir x y z: where --test-walk / --test-stand spawn on body 0
    Vector3 spawnDir = { 0, 1, 0 };
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
    int benchFps = -1;             // --fps N: cap the frame rate (0 = uncapped); default: bench uncapped, play at the monitor refresh rate
    bool lookDown = false;         // --look-down: tilt the player view toward the ground (lighting tests)
    Vector3 sunArg = { 1, 0, 0 };
    double parkStart = -1.0;
    int parkPhase = 0;
    bool benchGravity = false;     // --test-gravity: gravity field checks, ship drop and hover
    bool benchPrecision = false;   // --test-precision: float vs double placement error on body 0 (use a large body)
    bool benchMantle = false;      // --test-mantle: planet terrain checks on body 0 (height band, mantle, ore, caves)
    bool benchAirJump = false;     // --test-airjump: air control while jumping, jetpack fuel / suit power accounting
    float timeScale = 1.0f;        // --time-scale k: bodies spin and orbit k times faster (testing)
    bool benchMoon = false;        // --test-moon: moon orbit vs analytic, standing and jumping on the moon
    bool benchHole = false;        // --test-hole: walk into a fresh pit and build a pillar underfoot on the spinning planet
    bool benchWall = false;        // --test-wall: walk into shaft and tunnel walls (no launches)
    bool benchMenus = false;       // --test-menus: main menu, new game, pause, save & quit, load
    bool benchWater = false;       // --test-water: rivers and ponds near the start, their beds and clay
    bool benchFlyFollow = false;   // --test-flyfollow: the fly camera rides its reference body
    bool benchHotbar = false;      // --test-hotbar: materials and consumables from the hotbar
    bool benchRing = false, ringStarted = false; // --test-ring: milestone 6 gate
    bool benchRocket = false, rocketStarted = false; // --test-rocket: milestone 7 gate
    double rocketWait = 0.0;
    int ffPhase = 0;
    double ffStart = 0.0, ffDist0 = 0.0;
    DVec3 ffLocal0{};
    bool ffOk = true, ffBaseline = false;
    int waterPhase = 0;
    double waterStart = 0.0;
    Vector3 waterDir{};
    float waterLevel = 0.0f;
    int menuStep = 0, menuFrames = 0;
    bool menuOk = true;
    int wallPhase = 0;
    double wallStart = 0.0, wallMaxStep = 0.0, wallMaxStepAt = 0.0;
    float wallMaxImpact = 0.0f, wallMinHealth = 100.0f, wallMaxSpeed = 0.0f;
    DVec3 wallLast{};
    bool benchTrees = false;       // --test-trees: walk into a tree trunk
    bool benchChop = false;        // --test-chop: fell a tree, save/reload, regrow
    int benchPersist = 0;          // --test-persist N: milestone 2 gate (N strokes, flat builds, fast save/load)
    int benchBuild = 0;            // --test-build N: milestone 3 gate (N parts, standing on them, save/load, FPS)
    int buildPhase = 0, buildPlaced = 0;
    double buildStart = 0.0;
    bool buildOk = false;
    Plot* buildPlot = nullptr;
    int benchFactory = 0;          // --test-factory N: milestone 4 gate (N working machines, tick under 5 ms)
    bool benchProgression = false; // --test-progression: milestone 5 gate (bare hands to a working LV machine)
    int progPhase = 0;
    double progStart = 0.0;
    int factPhase = 0, factSamples = 0;
    uint32_t factPlot = 0;
    double factStart = 0.0, factSum = 0.0, factMax = 0.0;
    long long factMade0 = 0;
    std::string factChestBefore;
    int factRefundBefore = 0;
    std::vector<std::array<int, 3>> factSeeds;
    int chopPhase = 0, chopSwings = 0, chopMisses = 0, chopWood0 = 0;
    double chopStart = 0.0;
    uint64_t chopId = 0;
    bool chopFelled = false, chopSaved = false;
    Trunk chopTarget{};
    int treePhase = 0;
    double treeStart = 0.0;
    float treeMin = 1e9f;
    Trunk treeTarget{};
    int holePhase = 0;
    double holeStart = -1.0;
    float holeMaxImpact = 0.0f, holeMinHealth = 100.0f;
    Vector3 holeRim = { 0, 0, 0 };
    double moonStart = -1.0; int moonPhase = 0; DVec3 moonLocal0{}, moonTakeoff{}; bool moonAirborne = false;
    double airStart = -1.0; int airPhase = 0; float airDist = 0.0f, airTime = 0.0f, airFuel0 = 0.0f, airPower0 = 0.0f, airPowerAir = 0.0f;
    Vector3 airPrevLocal{};
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
    double walkTargetM = 0.0, walkPath = 0.0; // --test-planet-walk M: walk M metres (path, body frame), then PASS/FAIL
    DVec3 walkLast{};
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
        if (std::string(argv[i]) == "--test-planet-walk" && i + 1 < argc) { walkTargetM = atof(argv[i + 1]); benchWalkSeconds = 600.0; } // gate: 1 km on the planet
        if (std::string(argv[i]) == "--bench-walk") benchWalkMode = true;
        if (std::string(argv[i]) == "--bench-fly" && i + 1 < argc) benchFlySpeed = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--bench-fly-alt" && i + 1 < argc) benchFlyAlt = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--spawn-dir" && i + 3 < argc) { hasSpawnDir = true; spawnDir = Vector3Normalize({ (float)atof(argv[i + 1]), (float)atof(argv[i + 2]), (float)atof(argv[i + 3]) }); }
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
        if (std::string(argv[i]) == "--test-airjump") benchAirJump = true;
        if (std::string(argv[i]) == "--test-moon") benchMoon = true;
        if (std::string(argv[i]) == "--test-hole") benchHole = true;
        if (std::string(argv[i]) == "--test-wall") benchWall = true;
        if (std::string(argv[i]) == "--test-menus") benchMenus = true;
        if (std::string(argv[i]) == "--test-water") benchWater = true;
        if (std::string(argv[i]) == "--test-flyfollow") benchFlyFollow = true;
        if (std::string(argv[i]) == "--test-hotbar") benchHotbar = true;
        if (std::string(argv[i]) == "--test-ring") benchRing = true;
        if (std::string(argv[i]) == "--test-rocket") benchRocket = true;
        if (std::string(argv[i]) == "--test-trees") benchTrees = true;
        if (std::string(argv[i]) == "--test-chop") benchChop = true;
        if (std::string(argv[i]) == "--test-persist" && i + 1 < argc) benchPersist = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-build" && i + 1 < argc) benchBuild = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-factory" && i + 1 < argc) benchFactory = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-progression") benchProgression = true;
        if (std::string(argv[i]) == "--time-scale" && i + 1 < argc) timeScale = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--crack-check") crackCheck = true;
        if (std::string(argv[i]) == "--no-foliage") foliageEnabled = false;
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
                FlyInput in;
                in.pitch = 0.02f;
                in.yaw = 0.003f;
                ApplyFlyInput(&c, in, 0.0f);
                const Vector3 fwd = Vector3Normalize(Vector3Subtract(c.target, c.position));
                worstStep = std::max(worstStep, Vector3Angle(prev, fwd));
                worstOrtho = std::max(worstOrtho, fabsf(Vector3DotProduct(fwd, Vector3Normalize(c.up))));
                prev = fwd;
            }
            printf("bench: flycam pitched 8 rad (past both poles): largest per-step turn %.4f rad (expected ~0.02), max |dot(view, up)| %.4f\n",
                   worstStep, worstOrtho);
            // Movement is camera-relative at any orientation: after random yaw / pitch / roll, each key
            // moves exactly along the camera's own forward / right / up, and roll (Q/E) turns up about
            // the view without moving the view.
            std::mt19937 rng(7);
            std::uniform_real_distribution<float> ang(-0.3f, 0.3f);
            float worstAxis = 1.0f, worstRollView = 0.0f, rolled = 0.0f;
            for (int s = 0; s < 300; ++s) {
                FlyInput turn;
                turn.yaw = ang(rng); turn.pitch = ang(rng); turn.roll = ang(rng);
                ApplyFlyInput(&c, turn, 0.0f);
                const Vector3 f = Vector3Normalize(Vector3Subtract(c.target, c.position));
                const Vector3 u = c.up, r = Vector3CrossProduct(f, u);
                const Vector3 axes[3] = { f, r, u };
                for (int a = 0; a < 3; ++a) {
                    FlyInput mv;
                    (a == 0 ? mv.forward : a == 1 ? mv.right : mv.up) = 1.0f;
                    const Vector3 before = c.position;
                    ApplyFlyInput(&c, mv, 0.1f);
                    worstAxis = std::min(worstAxis, Vector3DotProduct(Vector3Normalize(Vector3Subtract(c.position, before)), axes[a]));
                }
                FlyInput roll;
                roll.roll = 0.05f; // E
                const Vector3 upBefore = c.up;
                ApplyFlyInput(&c, roll, 0.0f);
                worstRollView = std::max(worstRollView, Vector3Angle(f, Vector3Normalize(Vector3Subtract(c.target, c.position))));
                rolled += Vector3Angle(upBefore, c.up);
            }
            printf("bench: flycam moves along its own axes after random turns: worst alignment %.6f (1 = exact); roll turns up %.3f rad/step, view moves %.6f rad\n",
                   worstAxis, rolled / 300.0f, worstRollView);
        }
        if (std::string(argv[i]) == "--spin-rate" && i + 1 < argc) spinRate = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-fly" && i + 1 < argc) benchFlySeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-mine" && i + 1 < argc) benchMineStrokes = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-stand" && i + 1 < argc) { benchWalkSeconds = atof(argv[i + 1]); benchStandOnly = true; }
        if (std::string(argv[i]) == "--telescope" && i + 1 < argc) benchTelescope = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--telescope-target" && i + 1 < argc) benchTelescopeTarget = argv[i + 1];
        if (std::string(argv[i]) == "--size" && i + 1 < argc) worldGen.size = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--bodies" && i + 1 < argc) worldGen.bodiesPath = argv[i + 1];   // system definition file
        if (std::string(argv[i]) == "--count" && i + 1 < argc) worldGen.count = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--system" && i + 1 < argc) worldGen.systemRadius = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--lod" && i + 1 < argc) lodSettings.lodFactor = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--first" && i + 3 < argc) { worldGen.hasFirst = true; worldGen.first = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
    }
    const bool bench = benchSeconds > 0.0;
    benchMode = bench;
    if (bench) SetWorldSavePath("saves/bench.sav"); // tests save and load; keep the player's F5 save out of it

    SetConfigFlags(FLAG_MSAA_4X_HINT);
    if (!bench) {
        SetConfigFlags(FLAG_WINDOW_RESIZABLE); // Allow window resizing
        SetConfigFlags(FLAG_WINDOW_UNDECORATED); // Remove window decorations for a cleaner look
        SetConfigFlags(FLAG_WINDOW_HIGHDPI); // Enable High DPI support
        InitWindow(GetMonitorWidth(0), GetMonitorHeight(0), "My first RAYLIB program!");
        // No VSync: cap at the monitor's refresh rate instead (VSync on a 60 Hz desktop locked the game to 60).
        SetTargetFPS(benchFps >= 0 ? benchFps : std::max(60, GetMonitorRefreshRate(GetCurrentMonitor())));
    } else {
        InitWindow(1920, 1080, "RaylibGame bench");
        SetTargetFPS(std::max(0, benchFps)); // --fps N caps it (frame-rate-dependence checks); 0 = uncapped
    }
    SetExitKey(KEY_NULL);
    SetTraceLogLevel(LOG_WARNING); // Set log level to warning to reduce output noise
    rlSetClipPlanes(0.5, 60000.0); // space scale: see other bodies kilometres away (0.2 m near walls: see the frame loop)
    
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
    Rocket_Init(&world);
    if (loadSaveAtStart) { // the player's save, read-only; later bench saves still go to bench.sav
        SetWorldSavePath("saves/world.sav");
        if (!LoadWorldSave(world)) printf("bench: --load: could not load saves/world.sav\n");
        if (bench) SetWorldSavePath("saves/bench.sav");
    }
    if (dumpSystem) {
        for (auto& obj : world.objects)
            if (auto* p = dynamic_cast<Planetoid*>(obj.get()))
                printf("body %s %s pos %.3f %.3f %.3f radius %.3f seed %u spin %.4f g %.4f\n", BodyKindName(p->kind), p->name.c_str(), p->worldPos.x, p->worldPos.y, p->worldPos.z,
                       p->radius, p->seed, p->spin.y, p->surfaceGravity);
    }
    StartNewWorld(world, false);
    if (!bench) Settings_Load();
    bool worldFresh = !loadSaveAtStart; // the world generated at launch: New Game uses it as it is
    MenuResult pendingMenu;
    bool justPaused = false;
    world.onDraw3D = [&world] { if (!ship.piloted && !buildTool.active && !Rocket_Boarded()) TerrainTool_Draw3D(); BuildTool_Draw3D(world); ship.draw(world); };

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
        player.flyFollow = false;                                     // the bench drives the camera itself
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
                        survival.breathableAltitude = -1e9f; // no air to breathe (the planet's would refill the suit)
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
                } else if (benchChop && (chopPhase > 0 || WorldIdle())) {
                    // --test-chop: walk up to the nearest tree, chop it down (wood), save and reload (still
                    // gone), then wind the clock past regrowth (back, full size).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    auto present = [&](uint64_t id) {
                        std::vector<Trunk> dummy;
                        PlantHit h;
                        // Re-aim at the trunk: present if the raycast still finds this id.
                        const Vector3 eye = p->toLocal(Vector3Add(feet, Vector3Scale(player.up, player.eyeHeight)));
                        const Vector3 aim = Vector3Normalize(Vector3Subtract(Vector3Add(chopTarget.base, Vector3Scale(chopTarget.up, 1.2f)), eye));
                        return p->trees->raycast(eye, aim, 8.0f, h) && h.plant.id == id;
                    };
                    if (chopPhase == 0) {
                        if (!p || !p->trees) { printf("bench: --test-chop needs a planet with trees\n"); benchChop = false; }
                        else { player.spawnOn(world, *p, StartSpawnDir(*p)); chopPhase = 1; chopStart = GetTime(); }
                    } else if (chopPhase == 1 && GetTime() - chopStart > 1.0 && player.grounded) {
                        std::vector<Trunk> near;
                        const Vector3 lf = p->toLocal(feet);
                        p->trees->trunksNear(lf, 400.0f, near);
                        float best = 1e9f;
                        for (const Trunk& t : near) { const float d = Vector3Distance(t.base, lf); if (d < best && t.radius > 0.15f) { best = d; chopTarget = t; } }
                        if (best > 1e8f) { printf("bench: --test-chop: no tree nearby\n"); benchChop = false; }
                        else {
                            const Vector3 up = Vector3Normalize(chopTarget.base);
                            Vector3 away = Vector3Subtract(lf, chopTarget.base);
                            away = Vector3Normalize(Vector3Subtract(away, Vector3Scale(up, Vector3DotProduct(away, up))));
                            player.spawnOn(world, *p, p->dirToWorld(Vector3Normalize(Vector3Add(chopTarget.base, Vector3Scale(away, 2.0f)))));
                            chopPhase = 2; chopStart = GetTime();
                        }
                    } else if (chopPhase == 2 && GetTime() - chopStart > 1.0) {
                        // Face the trunk at chest height and swing.
                        const Vector3 eye = Vector3Add(feet, Vector3Scale(player.up, player.eyeHeight));
                        const Vector3 aim = p->toWorld(Vector3Add(chopTarget.base, Vector3Scale(chopTarget.up, 1.2f)));
                        player.forward = Vector3Normalize(Vector3Subtract(aim, eye));
                        player.applyToCamera(world);
                        TerrainTool_Update(world, false, true); // aim only
                        uint64_t id = 0;
                        if (!TerrainTool_AimingPlant(&id)) {
                            if (++chopMisses > 30) { printf("bench: --test-chop: could not aim at the trunk: FAIL\n"); benchChop = false; }
                        } else {
                            if (!chopId) { chopId = id; chopWood0 = inventory.count("wood"); }
                            ++chopSwings;
                            if (TerrainTool_Chop()) {
                                const int wood = inventory.count("wood") - chopWood0;
                                printf("bench: chop: felled a %.1f m tree in %d swings, +%d wood; gone: %s\n", chopTarget.height, chopSwings, wood, present(chopId) ? "no" : "yes");
                                chopFelled = wood > 0 && !present(chopId);
                                SaveWorld(world);
                                LoadWorldSave(world);
                                chopPhase = 3; chopStart = GetTime();
                            }
                        }
                    } else if (chopPhase == 3 && GetTime() - chopStart > 3.0) {
                        p = dynamic_cast<Planetoid*>(world.objects[0].get());
                        const bool savedGone = p->trees->removals().count(chopId) && !present(chopId);
                        printf("bench: chop: after save + reload the tree is %s\n", savedGone ? "still gone" : "BACK");
                        SetWorldTime(WorldTime() + kRegrowDelay + kRegrowTime + 1.0);
                        p->trees->setRemovals(p->trees->removals()); // rebuild its tiles at the new time
                        chopPhase = 4; chopStart = GetTime();
                        chopSaved = savedGone;
                    } else if (chopPhase == 4 && GetTime() - chopStart > 3.0) {
                        const bool back = present(chopId);
                        printf("bench: chop: %.0f min later it has regrown: %s -> %s\n", (kRegrowDelay + kRegrowTime) / 60.0, back ? "yes" : "no",
                               chopFelled && chopSaved && back ? "PASS" : "FAIL");
                        benchChop = false;
                    }
                } else if (benchTrees && (treePhase > 0 || WorldIdle())) {
                    // --test-trees: walk into the nearest tree trunk; the player must stop at its bark.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (!p || !p->trees) { printf("bench: --test-trees needs a planet with trees\n"); benchTrees = false; }
                    else if (treePhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        treePhase = 1;
                        treeStart = GetTime();
                    } else if (treePhase == 1 && GetTime() - treeStart > 1.0 && player.grounded) {
                        std::vector<Trunk> near;
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        p->trees->trunksNear(p->toLocal(feet), 400.0f, near);
                        const Vector3 lf = p->toLocal(feet);
                        float best = 1e9f;
                        for (const Trunk& t : near) { const float d = Vector3Distance(t.base, lf); if (d < best && t.radius > 0.15f) { best = d; treeTarget = t; } }
                        if (best > 1e8f) { printf("bench: no tree trunk within 400 m (%d tiles, %d instances)\n", p->trees->loadedTiles(), p->trees->instanceCount()); benchTrees = false; }
                        else {
                            // Stand 6 m from it, facing it.
                            const Vector3 up = Vector3Normalize(treeTarget.base);
                            Vector3 away = Vector3Subtract(lf, treeTarget.base);
                            away = Vector3Normalize(Vector3Subtract(away, Vector3Scale(up, Vector3DotProduct(away, up))));
                            player.spawnOn(world, *p, p->dirToWorld(Vector3Normalize(Vector3Add(treeTarget.base, Vector3Scale(away, 6.0f)))));
                            {
                                const Vector3 nl = p->toLocal((player.position - world.renderOrigin).toFloat());
                                printf("bench: trees: target %.1f m from the start, spawned %.1f m from it (trunk height %.1f, base r %.1f, surface r %.1f)\n",
                                       best, Vector3Distance(nl, treeTarget.base), treeTarget.height, Vector3Length(treeTarget.base), Vector3Length(nl));
                            }
                            treePhase = 2;
                            treeStart = GetTime();
                            treeMin = 1e9f;
                        }
                    } else if (treePhase == 2) {
                        const Vector3 lf = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 d = Vector3Subtract(lf, treeTarget.base);
                        const Vector3 side = Vector3Subtract(d, Vector3Scale(treeTarget.up, Vector3DotProduct(d, treeTarget.up)));
                        if (GetTime() - treeStart > 1.0 && treeMin > 1e8f) printf("bench: trees: start %.1f m from the trunk axis, grounded %d\n", Vector3Length(side), player.grounded);
                        if (GetTime() - treeStart > 1.0) {
                            player.forward = p->dirToWorld(Vector3Normalize(Vector3Negate(side)));
                            player.autoWalk = true;
                            treeMin = std::min(treeMin, Vector3Length(side));
                        }
                        if (GetTime() - treeStart > 5.0) {
                            player.autoWalk = false;
                            const float limit = treeTarget.radius + player.radius - 0.05f;
                            printf("bench: trees: %d tiles, %d instances drawn; walked into a %.1f m trunk (r %.2f m): closest %.2f m from its axis (limit %.2f): %s\n",
                                   p->trees->loadedTiles(), p->trees->instanceCount(), treeTarget.height, treeTarget.radius, treeMin, limit, treeMin >= limit && treeMin < limit + 0.3f ? "PASS" : "FAIL");
                            benchTrees = false;
                        }
                    }
                } else if (benchPersist > 0 && WorldIdle()) {
                    // --test-persist N (milestone 2 gate): N strokes over ~600x600 m around the start;
                    // chunk builds must stay flat, a save and a load must take under a second.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { printf("bench: --test-persist needs a planet\n"); benchPersist = 0; }
                    else {
                        using clk = std::chrono::steady_clock;
                        auto ms = [](clk::time_point a) { return std::chrono::duration<double, std::milli>(clk::now() - a).count(); };
                        const Vector3 centre = Vector3Normalize(p->toLocal((player.position - world.renderOrigin).toFloat()));
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(centre, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, centre);
                        // Builds: 24 finest nodes around the area, timed before and after.
                        std::vector<NodeKey> keys;
                        for (int k = 0; k < 24; ++k) {
                            const float a = k * 0.26f, rr = 40.0f + 10.0f * k;
                            const Vector3 d = Vector3Normalize(Vector3Add(centre, Vector3Add(Vector3Scale(t1, cosf(a) * rr / pf->radius), Vector3Scale(t2, sinf(a) * rr / pf->radius))));
                            const Vector3 s = Vector3Scale(d, pf->radius + pf->elevation(d));
                            keys.push_back({ 0, (int)floorf(s.x / CHUNK_SIZE), (int)floorf(s.y / CHUNK_SIZE), (int)floorf(s.z / CHUNK_SIZE) });
                        }
                        auto timeBuilds = [&]() {
                            const auto t = clk::now();
                            for (const NodeKey& k : keys) { ChunkMeshData m; m.key = k; BuildChunkMesh(*p->field, *p->edits, m); }
                            return ms(t) / keys.size();
                        };
                        const double buildBefore = timeBuilds();
                        std::mt19937 rng(7);
                        std::uniform_real_distribution<float> u(0.0f, 1.0f);
                        std::vector<TerrainEdit> strokes((size_t)benchPersist);
                        for (TerrainEdit& e : strokes) {
                            const Vector3 d = Vector3Normalize(Vector3Add(centre, Vector3Add(Vector3Scale(t1, (u(rng) - 0.5f) * 600.0f / pf->radius), Vector3Scale(t2, (u(rng) - 0.5f) * 600.0f / pf->radius))));
                            const float depth = -8.0f * u(rng) + 1.0f;
                            e.center = Vector3Scale(d, pf->radius + pf->elevation(d) + depth);
                            e.radius = 1.5f + 2.5f * u(rng);
                            e.strength = u(rng) < 0.8f ? -0.12f - 0.18f * u(rng) : 0.12f + 0.1f * u(rng);
                        }
                        auto t = clk::now();
                        ApplyTerrainEdits(*p, strokes.data(), strokes.size());
                        const double foldMs = ms(t);
                        const double buildAfter = timeBuilds();
                        t = clk::now();
                        SaveWorld(world);
                        const double saveFull = ms(t);
                        std::vector<TerrainEdit> more(strokes.begin(), strokes.begin() + std::min<size_t>(100, strokes.size()));
                        ApplyTerrainEdits(*p, more.data(), more.size());
                        t = clk::now();
                        SaveWorld(world);
                        const double saveIncr = ms(t);
                        const size_t chunks = p->edits->chunks.size();
                        const float probe = p->densityAt(strokes[0].center);
                        t = clk::now();
                        const bool loaded = LoadWorldSave(world);
                        const double loadMs = ms(t);
                        p = dynamic_cast<Planetoid*>(world.objects[0].get());
                        const bool same = loaded && p->edits->chunks.size() == chunks && fabsf(p->densityAt(strokes[0].center) - probe) < 1e-4f;
                        const bool pass = same && saveFull < 1000.0 && loadMs < 1000.0 && buildAfter < buildBefore * 1.5 + 0.5;
                        printf("bench: persist: %d strokes folded in %.0f ms into %zu delta chunks (%zu strokes recorded)\n", benchPersist, foldMs, chunks, p->edits->strokes);
                        printf("bench: persist: chunk build %.2f ms before, %.2f ms after; save %.0f ms (full), %.0f ms (100 more strokes); load %.0f ms; reloaded identical: %s -> %s\n",
                               buildBefore, buildAfter, saveFull, saveIncr, loadMs, same ? "yes" : "NO", pass ? "PASS" : "FAIL");
                        benchPersist = 0;
                    }
                } else if (benchBuild > 0 && (buildPhase > 0 || WorldIdle())) {
                    // --test-build N (milestone 3 gate): N parts on a plot by the start, the player standing
                    // on them, a save/load round trip, then the frame rate with the factory in view.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (buildPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        buildPhase = 1;
                        buildStart = GetTime();
                    } else if (buildPhase == 1 && GetTime() - buildStart > 1.0 && player.grounded) {
                        using clk = std::chrono::steady_clock;
                        const auto t0 = clk::now();
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        TerrainHit th;
                        const Vector3 above = Vector3Add(Vector3Add(feet, Vector3Scale(flat, 12.0f)), Vector3Scale(player.up, 20.0f));
                        RaycastTerrain(world, Ray{ above, Vector3Negate(player.up) }, 60.0f, th);
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(th.local, 0.0f, true);
                        buildPlot = plot;
                        std::mt19937 rng(11);
                        int placed = 0, tries = 0;
                        auto put = [&](PartType t, Int3 c, int o) { ++tries; if (BuildTool_Place(*p, *plot, { t, c, (uint8_t)o }, false)) ++placed; };
                        // A 48 x 48 floor slab (cubes), walls 12 high round it, then a mix inside.
                        for (int z = 0; z < 48 && placed < benchBuild; ++z) for (int x = 0; x < 48 && placed < benchBuild; ++x) put(PartType::Cube, { x, 0, z }, 0);
                        for (int y = 1; y <= 12 && placed < benchBuild; ++y)
                            for (int k = 0; k < 48 && placed < benchBuild; ++k) {
                                put(PartType::Cube, { k, y, 0 }, 0); put(PartType::Cube, { k, y, 47 }, 0);
                                put(PartType::Cube, { 0, y, k }, 0); put(PartType::Cube, { 47, y, k }, 0);
                            }
                        const PartType mix[6] = { PartType::Block, PartType::Beam4, PartType::Slope, PartType::Corner, PartType::Beam8, PartType::Cube };
                        while (placed < benchBuild && tries < benchBuild * 20) {
                            const Int3 c = { 2 + (int)(rng() % 44), 1 + (int)(rng() % 20), 2 + (int)(rng() % 44) };
                            put(mix[rng() % 6], c, (int)(rng() % 24));
                        }
                        const double ms = std::chrono::duration<double, std::milli>(clk::now() - t0).count();
                        printf("bench: build: %d parts placed in %.0f ms (%d tries), %zu cells on plot %u\n", placed, ms, tries, plot->cells.size(), plot->id);
                        // A clear 6 x 6 pad outside the walls to stand on.
                        for (int z = 0; z < 6; ++z) for (int x = -8; x < -2; ++x) put(PartType::Cube, { x, 0, z }, 0);
                        buildPlaced = plot->partCount();
                        // Stand on the slab: drop the player 2 m above a free corner of it.
                        const Vector3 top = p->toWorld(plot->toBody({ -2.75f, 0.5f + 2.0f, 1.5f })); // over the pad
                        player.placeAt(world, *p, top);
                        player.velocity = p->surfaceVelocity(top);
                        buildPhase = 2;
                        buildStart = GetTime();
                    } else if (buildPhase == 2 && GetTime() - buildStart > 2.0) {
                        Plot* plot = buildPlot;
                        const Vector3 lf = plot->toPlot(p->toLocal((player.position - world.renderOrigin).toFloat()));
                        const bool standing = player.grounded && fabsf(lf.y - 0.5f) < 0.15f;
                        printf("bench: build: the player stands on the slab: %s (feet %.2f m above the plot floor, grounded %d)\n", standing ? "yes" : "NO", lf.y, player.grounded);
                        const bool saved = SaveWorld(world), loaded = LoadWorldSave(world);
                        p = dynamic_cast<Planetoid*>(world.objects[0].get());
                        const int after = p->plots ? p->plots->partCount() : 0;
                        printf("bench: build: save %s, load %s: %d -> %d parts\n", saved ? "ok" : "FAILED", loaded ? "ok" : "FAILED", buildPlaced, after);
                        buildOk = standing && saved && loaded && after == buildPlaced;
                        // Look at the factory from 30 m for the frame-time measurement.
                        Plot* reloaded = p->plots ? p->plots->plots.front().get() : nullptr;
                        if (reloaded) {
                            const Vector3 centre = p->toWorld(reloaded->toBody({ 12.0f, 3.0f, 12.0f }));
                            player.spawnOn(world, *p, p->dirToWorld(Vector3Normalize(reloaded->toBody({ 12.0f, 0.0f, -18.0f }))));
                            const Vector3 eye = Vector3Add((player.position - world.renderOrigin).toFloat(), Vector3Scale(player.up, player.eyeHeight));
                            player.forward = Vector3Normalize(Vector3Subtract(centre, eye));
                        }
                        printf("bench: build gate (parts, standing, save/load): %s; frame time with the factory in view follows\n", buildOk ? "PASS" : "FAIL");
                        benchBuild = 0;
                    }
                } else if (benchFactory > 0 && (factPhase > 0 || WorldIdle())) {
                    // --test-factory N (milestone 4 gate): N machines in groups of generator + 2 furnaces +
                    // chest, joined by a cable and a conveyor, all working: the 20 Hz tick must stay
                    // under 5 ms, and plates must reach the chests.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (factPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        factPhase = 1;
                        factStart = GetTime();
                    } else if (factPhase == 1 && GetTime() - factStart > 1.0 && player.grounded) {
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(p->toLocal(Vector3Add(Vector3Add(feet, Vector3Scale(flat, 20.0f)), Vector3Scale(player.up, 3.0f))), 0.0f, true);
                        factPlot = plot->id;
                        const int groups = benchFactory / 4;
                        const int perRow = (int)ceilf(sqrtf((float)groups));
                        for (int g = 0; g < groups; ++g) {
                            const int bx = (g % perRow) * 12, bz = (g / perRow) * 5, by = 0;
                            const int gen = plot->place({ PartType::Generator, { bx, by, bz }, 0 });
                            const int f1 = plot->place({ PartType::Furnace, { bx + 3, by, bz }, 0 });
                            const int f2 = plot->place({ PartType::Furnace, { bx + 6, by, bz }, 0 });
                            plot->place({ PartType::Chest, { bx + 9, by, bz }, 0 });
                            for (int x = 0; x < 8; ++x) plot->place({ PartType::Cable, { bx + x, by, bz - 1 }, 0 });
                            for (int x = 3; x < 11; ++x) plot->place({ PartType::Conveyor, { bx + x, by, bz + 2 }, 0 });
                            factSeeds.push_back({ gen, f1, f2 });
                        }
                        printf("bench: factory: %d groups, %d parts on plot %u\n", groups, plot->partCount(), plot->id);
                        factPhase = 2; factStart = GetTime();
                    } else if (factPhase == 2 && GetTime() - factStart > 1.0) {
                        for (const auto& s : factSeeds) {
                            Factory::Seed(p->id, factPlot, s[0], "wood", 64);
                            Factory::Seed(p->id, factPlot, s[1], "iron_ore", 32);
                            Factory::Seed(p->id, factPlot, s[2], "iron_ore", 32);
                        }
                        factPhase = 3; factStart = GetTime();
                    } else if (factPhase == 3 && GetTime() - factStart > 2.0) {
                        factMade0 = Factory::Stats().itemsMade;
                        Factory::View(p->id, factPlot, factSeeds[0][2] + 1); // watch the first chest from now on
                        factPhase = 4; factStart = GetTime(); factSamples = 0; factSum = 0.0; factMax = 0.0;
                    } else if (factPhase == 4) {
                        const FactoryStats s = Factory::Stats();
                        factSum += s.avgTickMs; ++factSamples;
                        if (GetTime() - factStart > 12.0) {
                            const double avg = factSum / std::max(1, factSamples);
                            const MachineView chest = Factory::View(p->id, factPlot, factSeeds[0][2] + 1);
                            const bool pass = avg < 5.0 && s.activeMachines >= benchFactory / 2 && s.itemsMade > 0;
                            printf("bench: factory: %d machines (%d active), %d item / %d power networks; tick %.2f ms avg, %.2f ms worst; %lld items made in the last 12 s (%lld total); first chest: %s -> %s\n",
                                   s.machines, s.activeMachines, s.itemNetworks, s.powerNetworks, avg, s.maxTickMs, s.itemsMade - factMade0, s.itemsMade, chest.input.c_str(), pass ? "PASS" : "FAIL");
                            factChestBefore = chest.input;
                            SaveWorld(world);
                            const bool loaded = LoadWorldSave(world);
                            if (!loaded) printf("bench: factory: load FAILED\n");
                            factPhase = 5; factStart = GetTime();
                        }
                    } else if (factPhase == 5) {
                        Planetoid* p2 = dynamic_cast<Planetoid*>(world.objects[0].get());
                        const MachineView chest = Factory::View(p2->id, factPlot, factSeeds[0][2] + 1);
                        if (GetTime() - factStart > 2.0) {
                            printf("bench: factory: first chest before save %s, after reload %s -> %s\n", factChestBefore.c_str(), chest.input.c_str(),
                                   chest.valid && chest.input != "-" && !chest.input.empty() ? "PASS" : "FAIL");
                            // Then remove that chest: what it held must come back to the inventory.
                            factRefundBefore = inventory.count("iron_plate");
                            if (Plot* pl = p2->plots ? p2->plots->byId(factPlot) : nullptr) pl->remove(factSeeds[0][2] + 1);
                            factPhase = 6; factStart = GetTime();
                        }
                    } else if (factPhase == 6 && GetTime() - factStart > 2.0) {
                        const int got = inventory.count("iron_plate") - factRefundBefore;
                        printf("bench: factory: removed the chest: %d iron plates back in the inventory -> %s\n", got, got > 0 ? "PASS" : "FAIL");
                        benchFactory = 0;
                    }
                } else if (benchProgression && (progPhase > 0 || WorldIdle())) {
                    // --test-progression (milestone 5 gate): bare hands to a working LV machine.
                    if (progPhase == 0) {
                        Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        progPhase = 1; progStart = GetTime();
                    } else if (progPhase == 1 && GetTime() - progStart > 2.0) {
                        progPhase = 2;
                    } else if (progPhase == 2) {
                        if (Progression_Update(world) != 0) benchProgression = false;
                    }
                } else if (benchWall && (wallPhase > 0 || WorldIdle())) {
                    // --test-wall: a 7 m shaft with a tunnel off its bottom; walk into the walls at eight
                    // headings, jumping now and then. Nothing may launch you: the largest per-frame move
                    // relative to the body, impacts, health and deaths are reported.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - wallStart;
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const DVec3 local = p->toLocalD(player.position);
                    if (wallPhase >= 2) {
                        const double step = (local - wallLast).length();
                        if (step > wallMaxStep) { wallMaxStep = step; wallMaxStepAt = t; }
                        if (step > 0.2) printf("bench: wall: %.2f m move at %.2f s (grounded %d, dt %.4f, v %.1f, contact dot %.2f, body density at feet %.2f)\n", step, t, player.grounded ? 1 : 0, GetFrameTime(), Vector3Length(player.velocity), player.lastContactDot, p->densityAt(p->toLocal(Vector3Add(feet, Vector3Scale(player.up, player.radius)))));
                        wallMaxImpact = std::max(wallMaxImpact, player.impactSpeed);
                        wallMinHealth = std::min(wallMinHealth, survival.health);
                        wallMaxSpeed = std::max(wallMaxSpeed, Vector3Length(Vector3Subtract(player.velocity, player.grounded ? Vector3{ 0, 0, 0 } : p->surfaceVelocity(feet))));
                    }
                    wallLast = local;
                    if (wallPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        wallStart = GetTime();
                        wallPhase = 1;
                    } else if (wallPhase == 1 && t > 1.5 && player.grounded) {
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        const Vector3 side = Vector3CrossProduct(flat, player.up);
                        for (int k = 0; k < 6; ++k)      // the shaft, ~3.5 m wide
                            ApplyTerrainEdit(*p, TerrainEdit{ p->toLocal(Vector3Add(feet, Vector3Add(Vector3Scale(flat, 4.0f), Vector3Scale(player.up, -1.2f * k)))), 1.8f, -10.0f });
                        for (int k = 1; k < 5; ++k)      // a tunnel off its bottom, 2.4 m high
                            ApplyTerrainEdit(*p, TerrainEdit{ p->toLocal(Vector3Add(feet, Vector3Add(Vector3Add(Vector3Scale(flat, 4.0f), Vector3Scale(side, 1.6f * k)), Vector3Scale(player.up, -5.6f)))), 1.25f, -10.0f });
                        player.autoWalk = true;
                        wallPhase = 2;
                        wallStart = GetTime();
                    } else if (wallPhase == 2 && t < 12.0) {
                        // Turn 45 degrees a second (and jump every other second) while walking.
                        const int k = (int)t;
                        static int lastK = -1;
                        if (k != lastK && t > 3.0) {
                            lastK = k;
                            const Matrix r = MatrixRotate(player.up, 45.0f * DEG2RAD);
                            player.forward = Vector3Normalize(Vector3Transform(player.forward, r));
                            if (k % 2 == 0) player.autoJump = true;
                        }
                    } else if (wallPhase == 2) {
                        player.autoWalk = false;
                        const bool pass = survival.deaths == 0 && wallMinHealth > 90.0f && wallMaxStep < 1.0;
                        printf("bench: wall shove 9 s: largest frame move %.2f m (at %.1f s), top speed vs ground %.1f m/s, worst impact %.1f m/s, health %.0f, deaths %d: %s\n",
                               wallMaxStep, wallMaxStepAt, wallMaxSpeed, wallMaxImpact, wallMinHealth, survival.deaths, pass ? "PASS" : "FAIL");
                        benchWall = false;
                    }
                } else if (benchWater && (waterPhase > 0 || WorldIdle())) {
                    // --test-water: inland water near the start (rivers, ponds): find some, look at it
                    // (bench_water.png), dig its bed for clay, and check the underwater test inside it.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { printf("bench: --test-water needs the planet\n"); benchWater = false; }
                    else if (waterPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed); // around the start (the bench camera may be over the sea)
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        int wet = 0, samples = 0;
                        float best = 1e30f;
                        for (float e = -4000.0f; e <= 4000.0f; e += 20.0f)
                            for (float n = -4000.0f; n <= 4000.0f; n += 20.0f) {
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                ++samples;
                                if (!pf->inlandWaterNear(d, 0.0f)) continue;
                                const float w = pf->waterElevation(d);
                                if (w == PlanetField::kNoWater || w <= 0.0f || pf->elevation(d) >= w - 0.2f) continue;
                                ++wet;
                                if (pf->elevation(d) > w - 1.0f) continue; // look at (and dig) a spot at least 1 m deep
                                const float dist = sqrtf(e * e + n * n);
                                if (dist < best) { best = dist; waterDir = d; waterLevel = w; }
                            }
                        printf("bench: water: %.2f%% of 64 km^2 around the start is inland water (%d of %d samples); nearest %.0f m away\n", 100.0f * wet / samples, wet, samples, best);
                        if (best > 1e29f) { printf("bench: water: FAIL (none)\n"); benchWater = false; }
                        else { waterPhase = 1; waterStart = GetTime(); player.flyMode = true; sunOverride = true; } // daylight for the picture
                    } else {
                        // Camera 22 m above the water, 30 m off to one side, looking at it.
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(waterDir, { 0, 1, 0 }));
                        const Vector3 at = Vector3Scale(waterDir, pf->radius + waterLevel);
                        const Vector3 eye = Vector3Add(Vector3Scale(waterDir, pf->radius + waterLevel + 22.0f), Vector3Scale(t1, 30.0f));
                        world.camera.position = p->toWorld(eye);
                        world.camera.target = p->toWorld(at);
                        world.camera.up = p->dirToWorld(waterDir);
                        sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(waterDir, Vector3Scale(t1, -0.6f))));
                        player.position = world.renderOrigin + world.camera.position;
                        if (waterPhase == 1 && GetTime() - waterStart > 3.0 && WorldIdle()) {
                            TakeScreenshot("bench_water.png");
                            const Vector3 inside = Vector3Scale(waterDir, pf->radius + waterLevel - 0.3f);
                            const bool under = IsUnderwater(world, p->toWorld(inside));
                            const bool above = IsUnderwater(world, p->toWorld(Vector3Scale(waterDir, pf->radius + waterLevel + 1.0f)));
                            const float bedE = pf->elevation(waterDir);
                            const Vector3 bed = Vector3Scale(waterDir, pf->radius + bedE - 0.5f);
                            const float base = pf->baseDensity(bed);
                            const VoxelMaterial m = pf->material(bed, base);
                            const int clay0 = inventory.count("clay");
                            for (int k = 0; k < 12; ++k) TerrainTool_DigAt(*p, Vector3Subtract(bed, Vector3Scale(waterDir, 0.4f * k)));
                            const int clay = inventory.count("clay") - clay0;
                            const bool pass = under && !above && clay > 0 && (m == VoxelMaterial::Sand || m == VoxelMaterial::Gravel);
                            printf("bench: water: surface %.1f m, bed %.1f m (%s); underwater 0.3 m down: %s, 1 m up: %s; 12 digs in the bed: %d clay -> %s\n",
                                   waterLevel, bedE, VoxelMaterialName(m), under ? "yes" : "no", above ? "yes" : "no", clay, pass ? "PASS" : "FAIL");
                            benchWater = false;
                        }
                    }
                } else if (benchFlyFollow && (ffPhase > 0 || WorldIdle())) {
                    // --test-flyfollow: the fly camera with no input stays put relative to its reference
                    // body: over the spinning planet (100 m up), near the orbiting moon (200 m up), and
                    // at 3 moon radii (inside its sphere of influence: its motion, not its spin).
                    Planetoid* planet = dynamic_cast<Planetoid*>(world.objects[0].get());
                    Planetoid* moon = nullptr;
                    for (auto& obj : world.objects) if (auto* q = dynamic_cast<Planetoid*>(obj.get())) if (q->kind == BodyKind::Moon && q->orbitParent) { moon = q; break; }
                    struct Case { Planetoid* body; Vector3 local; const char* name; };
                    const Case cases[3] = {
                        { planet, Vector3Scale(StartSpawnDirLocal(*planet), planet->radius + 400.0f), "planet, 100 m+ up" },
                        { moon, moon ? Vector3Scale(Vector3Normalize({ 0.3f, 1.0f, 0.2f }), moon->radius + 600.0f) : Vector3{}, "moon, low" },
                        { moon, moon ? Vector3Scale(Vector3Normalize({ 0.3f, 1.0f, 0.2f }), moon->radius * 3.0f) : Vector3{}, "moon, 3 radii" } };
                    if (!moon) { printf("bench: --test-flyfollow needs the moon\n"); benchFlyFollow = false; }
                    else if (ffPhase % 2 == 0) {
                        if (ffPhase / 2 >= 3) {
                            printf("bench: fly follow: %s\n", ffOk ? "PASS" : "FAIL");
                            benchFlyFollow = false;
                        } else {
                            const Case& c = cases[ffPhase / 2];
                            player.flyMode = true;
                            player.flyFollow = true;
                            world.camera.position = c.body->toWorld(c.local);
                            world.camera.target = Vector3Add(world.camera.position, c.body->dirToWorld({ 1, 0, 0 }));
                            world.camera.up = c.body->dirToWorld(Vector3Normalize(c.local));
                            ffStart = GetTime();
                            ffPhase++;
                            ffBaseline = true; // taken next frame, once the carry has a recorded pose
                        }
                    } else if (ffBaseline) {
                        const Case& c = cases[ffPhase / 2];
                        ffLocal0 = c.body->toLocalD(world.renderOrigin + world.camera.position);
                        ffDist0 = (world.renderOrigin + world.camera.position - c.body->worldPos).length();
                        ffBaseline = false;
                    } else if (GetTime() - ffStart > 3.0) {
                        const Case& c = cases[ffPhase / 2];
                        const DVec3 cam = world.renderOrigin + world.camera.position;
                        const double localDrift = (c.body->toLocalD(cam) - ffLocal0).length();
                        const double distDrift = fabs((cam - c.body->worldPos).length() - ffDist0);
                        const bool rotating = ffPhase / 2 < 2;
                        const bool ok = rotating ? localDrift < 0.1 : distDrift < 0.1;
                        ffOk &= ok;
                        printf("bench: fly follow (%s): 3 s later %.3f m from its spot on the body, distance to the centre changed %.3f m (body moved %.0f m) -> %s\n",
                               c.name, localDrift, distDrift, 3.0 * Vector3Length(c.body->linearVelocity), ok ? "ok" : "BAD");
                        ffPhase++;
                    }
                } else if (benchHotbar && WorldIdle()) {
                    // --test-hotbar: materials from the hotbar (any shape costs its material), a consumable,
                    // the foundation not digging outside the part, and materials kept through a save.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    bool ok = pf != nullptr;
                    if (ok) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        inventory.clear();
                        Hotbar_Reset();
                        inventory.add("concrete", 20); inventory.add("wood", 10); inventory.add("oxygen_canister", 1);
                        const Vector3 home = StartSpawnDirLocal(*p);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 }));
                        const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Scale(t1, 10.0f / pf->radius)));
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(Vector3Scale(d, pf->radius + pf->elevation(d)), 0.0f, true);
                        Hotbar_Press(1, false); // key 2: concrete
                        const bool concreteMode = buildTool.active && GetBuildMaterial(buildTool.material).item == std::string("concrete");
                        const int c0 = inventory.count("concrete");
                        const bool slope = BuildTool_Place(*p, *plot, { PartType::Slope, { 0, 0, 0 }, 0, (uint8_t)buildTool.material }, true);
                        const int c1 = inventory.count("concrete");
                        Hotbar_Press(3, false); // key 4: wood
                        const int w0 = inventory.count("wood");
                        const bool beam = BuildTool_Place(*p, *plot, { PartType::Beam4, { 2, 0, 0 }, 0, (uint8_t)buildTool.material }, true);
                        const int w1 = inventory.count("wood");
                        // Foundation: the ground a little outside a 1 m block keeps its density.
                        auto outside = [&](float off) { return plot->toBody({ -0.2f + off, 0.3f, 2.5f }); };
                        const float before = p->densityAt(outside(0.0f));
                        Hotbar_Press(1, false);
                        const bool block = BuildTool_Place(*p, *plot, { PartType::Block, { 0, -1, 4 }, 0, (uint8_t)buildTool.material }, true);
                        const float after = p->densityAt(outside(0.0f));
                        survival.oxygen = 50.0f;
                        Hotbar_Press(8, false); // key 9: oxygen canister
                        const bool used = inventory.count("oxygen_canister") == 0 && survival.oxygen > 90.0f;
                        Hotbar_Press(0, false); // key 1: terrain tool
                        const bool toolMode = !buildTool.active;
                        const bool saved = SaveWorld(world), loaded = saved && LoadWorldSave(world);
                        Planetoid* q = dynamic_cast<Planetoid*>(world.objects[0].get());
                        int woodParts = 0, concreteParts = 0;
                        if (q && q->plots)
                            for (auto& pl : q->plots->plots)
                                for (const Part& part : pl->parts) {
                                    if (part.type == PartType::Count) continue;
                                    if (PartItem(part) == std::string("wood")) ++woodParts;
                                    if (PartItem(part) == std::string("concrete")) ++concreteParts;
                                }
                        ok = concreteMode && slope && c0 - c1 == 1 && beam && w0 - w1 == 2 && block && fabsf(after - before) < 0.02f && used && toolMode && loaded && woodParts == 1 && concreteParts == 2;
                        printf("bench: hotbar: concrete mode %d, slope -%d concrete, wood beam4 -%d wood, block placed %d (ground outside %.3f -> %.3f), O2 used %d, tool mode %d, save/load %d/%d: %d wood + %d concrete parts -> %s\n",
                               concreteMode, c0 - c1, w0 - w1, block, before, after, used, toolMode, saved, loaded, woodParts, concreteParts, ok ? "PASS" : "FAIL");
                    }
                    benchHotbar = false;
                } else if (benchRing && (ringStarted || WorldIdle())) {
                    if (!ringStarted) { ringStarted = true; Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get()); player.spawnOn(world, *p, StartSpawnDir(*p)); }
                    else if (RingTest_Update(world) != 0) benchRing = false;
                } else if (benchRocket && (rocketStarted || WorldIdle())) {
                    if (!rocketStarted) { rocketStarted = true; Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get()); player.spawnOn(world, *p, StartSpawnDir(*p)); rocketWait = GetTime(); }
                    else if (GetTime() - rocketWait > 1.5 && RocketTest_Update(world) != 0) benchRocket = false;
                } else if (benchHole && (holePhase > 0 || WorldIdle())) { // edits keep the world busy: run on
                    // --test-hole: on the spinning planet, dig a pit ahead and walk into it, then build a
                    // pillar under your feet while jumping. Neither may hurt you (impacts are measured
                    // against the moving terrain, not the world frame).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = holeStart < 0.0 ? 0.0 : GetTime() - holeStart;
                    holeMaxImpact = std::max(holeMaxImpact, player.impactSpeed);
                    holeMinHealth = std::min(holeMinHealth, survival.health);
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    if (holePhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        holeStart = GetTime();
                        holePhase = 1;
                    } else if (holePhase == 1 && t > 1.5 && player.grounded) {
                        // A 3 m wide, ~5 m deep pit starting 3 m ahead.
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        for (int k = 0; k < 4; ++k)
                            for (int s = 0; s < 3; ++s)
                                ApplyTerrainEdit(*p, TerrainEdit{ p->toLocal(Vector3Add(feet, Vector3Add(Vector3Scale(flat, 3.5f + 1.5f * s), Vector3Scale(player.up, -1.2f * k)))), 1.8f, -10.0f });
                        player.autoWalk = true;
                        holePhase = 2;
                        holeStart = GetTime();
                    } else if (holePhase == 2 && t > 4.0) {
                        player.autoWalk = false;
                        printf("bench: walked into a fresh pit: worst impact %.1f m/s, health %.0f, deaths %d, %.1f m below the rim\n",
                               holeMaxImpact, holeMinHealth, survival.deaths, Vector3DotProduct(Vector3Subtract(feet, holeRim), player.up) * -1.0f);
                        holeMaxImpact = 0.0f;
                        holePhase = 3;
                        holeStart = GetTime();
                    } else if (holePhase == 3 && t < 4.0) {
                        // Build a pillar under the feet a stroke every 0.25 s, jumping now and then.
                        static double nextStroke = 0.0;
                        if (GetTime() >= nextStroke) {
                            ApplyTerrainEdit(*p, TerrainEdit{ p->toLocal(Vector3Subtract(feet, Vector3Scale(player.up, 0.6f))), 1.6f, 6.0f });
                            nextStroke = GetTime() + 0.25;
                            if (((int)(t * 4.0)) % 4 == 0) player.autoJump = true;
                        }
                    } else if (holePhase == 3) {
                        const bool pass = survival.deaths == 0 && holeMinHealth > 90.0f;
                        printf("bench: built a pillar underfoot for 4 s: worst impact %.1f m/s, health %.0f, deaths %d, speed vs ground %.1f m/s: %s\n",
                               holeMaxImpact, holeMinHealth, survival.deaths,
                               Vector3Length(Vector3Subtract(player.velocity, player.grounded ? Vector3{ 0, 0, 0 } : p->surfaceVelocity(feet))),
                               pass ? "PASS" : "FAIL");
                        benchHole = false;
                    }
                    if (holePhase == 1) holeRim = feet;
                } else if (benchMoon && WorldIdle()) {
                    // The moon (first body of kind moon, with an orbit): its position must match the
                    // analytic orbit; a player standing on it must stay put relative to it, and a jump must
                    // land near takeoff in its frame (the moon's orbital speed carries over on takeoff).
                    Planetoid* moon = nullptr;
                    for (auto& obj : world.objects)
                        if (auto* q = dynamic_cast<Planetoid*>(obj.get())) if (q->kind == BodyKind::Moon && q->orbitParent) { moon = q; break; }
                    if (!moon) { printf("bench: --test-moon needs an orbiting moon (the default resources/data/bodies.txt has one)\n"); benchMoon = false; }
                    else {
                        const double t = moonStart < 0.0 ? 0.0 : GetTime() - moonStart;
                        const DVec3 local = moon->toLocalD(player.position);
                        if (moonPhase == 0) {
                            player.spawnOn(world, *moon, Vector3Normalize({ 0.3f, 1.0f, 0.2f }));
                            moonStart = GetTime();
                            moonPhase = 1;
                        } else if (moonPhase == 1 && t > 2.5 && player.grounded) {
                            moonLocal0 = local;
                            moonPhase = 2;
                            moonStart = GetTime();
                        } else if (moonPhase == 2 && t > 3.0) {
                            const Planetoid* parent = moon->orbitParent;
                            const double a = moon->orbitPhase + moon->orbitRate * WorldTime();
                            const DVec3 expected = parent->worldPos + DVec3(cos(a) * moon->orbitRadius, 0.0, sin(a) * moon->orbitRadius);
                            printf("bench: moon orbit: radius %.0f km, period %.2f h, speed %.1f m/s (expected %.1f); off the analytic orbit by %.3f m at t = %.1f s\n",
                                   moon->orbitRadius / 1000.0, 2.0 * PI / moon->orbitRate / 3600.0, Vector3Length(moon->linearVelocity),
                                   moon->orbitRate * moon->orbitRadius, (moon->worldPos - expected).length(), WorldTime());
                            printf("bench: standing on the moon 3 s while it moved %.0f m: drift %.3f m relative to its surface (g %.2f m/s^2)\n",
                                   3.0 * Vector3Length(moon->linearVelocity), (local - moonLocal0).length(), player.gravity);
                            moonTakeoff = local;
                            player.autoJump = true;
                            moonAirborne = false;
                            moonPhase = 3;
                            moonStart = GetTime();
                        } else if (moonPhase == 3) {
                            if (!player.grounded) moonAirborne = true;
                            else if (moonAirborne || t > 20.0) {
                                printf("bench: jump on the moon: %.1f s in the air, landed %.2f m from takeoff in the moon's frame%s\n",
                                       t, (local - moonTakeoff).length(), moonAirborne ? "" : " [never took off]");
                                benchMoon = false;
                            }
                        }
                    }
                } else if (benchAirJump && WorldIdle()) {
                    // Walk forward and jump while still holding forward: air control must keep you moving at
                    // about walking speed relative to the ground, burning no fuel and no suit power. Then
                    // land with half a tank: recharging it is what costs suit power.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = airStart < 0.0 ? 0.0 : GetTime() - airStart;
                    const DVec3 local = p->toLocalD(player.position);
                    if (airPhase == 0) {
                        player.spawnOn(world, *p, Vector3Normalize({ 0.6f, 0.5f, 0.6f }));
                        player.autoWalk = true;
                        survival.lastJetpackFuel = player.jetpackFuel = 1.0f;
                        airStart = GetTime();
                        airPhase = 1;
                    } else if (airPhase == 1 && t > 2.0 && player.grounded) {
                        player.autoJump = true;
                        airFuel0 = player.jetpackFuel;
                        airPower0 = survival.power;
                        airDist = airTime = 0.0f;
                        airPrevLocal = local.toFloat();
                        airPhase = 2;
                    } else if (airPhase == 2) {
                        if (!player.grounded) {
                            const Vector3 l = local.toFloat();
                            const Vector3 step = Vector3Subtract(l, airPrevLocal);
                            const Vector3 n = Vector3Normalize(l);
                            airDist += Vector3Length(Vector3Subtract(step, Vector3Scale(n, Vector3DotProduct(step, n)))); // level distance, body frame
                            airTime += GetFrameTime();
                        } else if (airTime > 0.2f) {
                            airPowerAir = survival.power;
                            printf("bench: air jump: %.2f s airborne, %.1f m level relative to the ground (%.2f m/s, walk %.1f m/s); fuel %.2f -> %.2f, suit power %.2f -> %.2f\n",
                                   airTime, airDist, airDist / airTime, player.walkSpeed, airFuel0, player.jetpackFuel, airPower0, airPowerAir);
                            player.autoWalk = false;
                            player.jetpackFuel = 0.5f; // half a tank used; recharging it on the ground costs power
                            survival.lastJetpackFuel = 0.5f;
                            airPower0 = survival.power;
                            airStart = GetTime();
                            airPhase = 3;
                        }
                        airPrevLocal = local.toFloat();
                    } else if (airPhase == 3 && t > 3.0) {
                        printf("bench: jetpack recharge 0.50 -> %.2f cost %.2f suit power (expected ~%.1f plus %.2f baseline)\n",
                               player.jetpackFuel, airPower0 - survival.power, 0.5f * (player.jetpackFuel - 0.5f) / 0.5f * survival.jetpackPowerPerTank,
                               3.0f * survival.powerDrainPerSec);
                        benchAirJump = false;
                    }
                } else if (benchMantle && WorldIdle()) {
                    // Planet terrain checks on body 0 (the default system).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* f = dynamic_cast<const PlanetField*>(p->field.get());
                    if (!f) {
                        printf("bench: --test-mantle needs a planet as the first body (the default resources/data/bodies.txt has one)\n");
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
                        // Surface materials: the share of each over random surface points (with real slopes).
                        int counts[(int)VoxelMaterial::Count] = {};
                        for (int i = 0; i < 2000; ++i) {
                            const Vector3 dir = randomDir();
                            const Vector3 pt = Vector3Scale(dir, f->radius + f->elevation(dir));
                            const Vector3 n = f->normal(pt, nullptr, 0, 0.5f);
                            ++counts[(int)f->materialLod(pt, f->baseDensity(pt), 0.0f, n)];
                        }
                        std::string mix;
                        for (int m = 0; m < (int)VoxelMaterial::Count; ++m)
                            if (counts[m]) mix += TextFormat("%s %.1f%%  ", m == (int)VoxelMaterial::Regolith ? "dirt" : VoxelMaterialName((VoxelMaterial)m), 100.0f * counts[m] / 2000.0f);
                        printf("bench: surface materials: %s\n", mix.c_str());
                        // 2. Mantle: dig just above it; the crust must change and the mantle must not.
                        const Vector3 dir = randomDir();
                        const Vector3 inMantle = Vector3Scale(dir, f->mantleRadius - 4.0f), inCrust = Vector3Scale(dir, f->mantleRadius + 12.0f);
                        const float mantleBefore = p->densityAt(inMantle), crustBefore = p->densityAt(inCrust);
                        const auto savedEdits = p->edits;
                        ApplyTerrainEdit(*p, TerrainEdit{ Vector3Scale(dir, f->mantleRadius + 4.0f), 12.0f, -20.0f });
                        const float mantleAfter = p->densityAt(inMantle), crustAfter = p->densityAt(inCrust);
                        p->edits = savedEdits; // leave the world as it was
                        // A hole 700 m down (above the mantle) must open: deltas saturated at ~510 m of rock.
                        bool deepOk = true;
                        {
                            const Vector3 dd = randomDir();
                            const float depth = std::min(700.0f, f->radius + f->elevation(dd) - f->mantleRadius - 30.0f);
                            const Vector3 deep = Vector3Scale(dd, f->radius + f->elevation(dd) - depth);
                            ApplyTerrainEdit(*p, TerrainEdit{ deep, 6.0f, -(DENSITY_GRADIENT * depth + 20.0f) });
                            const float after = p->densityAt(deep);
                            p->edits = savedEdits;
                            deepOk = after < ISO_LEVEL;
                            printf("bench: a hole %.0f m down: density %.1f (air below %.2f): %s\n", depth, after, ISO_LEVEL, deepOk ? "open" : "STILL SOLID");
                        }
                        const bool mantleOk = deepOk && mantleAfter == mantleBefore && crustAfter < crustBefore - 1.0f && !f->editable(inMantle) && f->editable(inCrust);
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
                        // Which kind of cave: share of underground (30 m .. mantle) inside a tunnel / a cavern core.
                        int inTunnel = 0, inCavern = 0, probes = 0;
                        for (int i = 0; i < 40000; ++i) {
                            const Vector3 d = randomDir();
                            const float depth = 30.0f + unit(rng) * (f->radius - f->mantleRadius - 80.0f);
                            const Vector3 pt = Vector3Scale(d, f->radius + f->elevation(d) - depth);
                            if (!f->editable(pt)) continue;
                            ++probes;
                            if (std::max(fabsf(f->cavesA.fractal(2, pt.x, pt.y, pt.z)), fabsf(f->cavesB.fractal(2, pt.x, pt.y, pt.z))) < PlanetField::kTunnelLevel) {
                                static bool shown = false;
                                ++inTunnel;
                                if (!shown && depth < 80.0f && f->baseDensity(pt) < ISO_LEVEL && (shown = true))
                                    printf("bench: a tunnel %.0f m down at local %.1f %.1f %.1f (use --cam to look inside)\n", depth, pt.x, pt.y, pt.z);
                            }
                            if (f->caverns.fractal(2, pt.x, pt.y * 1.4f, pt.z) > PlanetField::kCavernLevel) ++inCavern;
                        }
                        printf("bench: caves 30 m .. mantle: tunnels %.2f%% of the volume, caverns %.2f%%\n",
                               100.0f * inTunnel / std::max(1, probes), 100.0f * inCavern / std::max(1, probes));
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
                        if (hasSpawnDir) player.spawnOn(world, *start, spawnDir); else player.spawnOn(world, *start, StartSpawnDir(*start));

                        walkStart = GetTime();
                        walkFrom = start->toLocalD(player.position); // body frame: the spin carries the walker too
                        walkMinDensity = 1e9f; walkMaxDensity = -1e9f; walkGroundedFrames = walkFrames = 0;
                        walkLast = walkFrom; walkPath = 0.0;
                    }
                    Planetoid* start = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const float d = start->densityAt(start->toLocal(Vector3Add(feet, Vector3Scale(player.up, 0.9f)))); // body centre
                    walkMinDensity = std::min(walkMinDensity, d);
                    walkMaxDensity = std::max(walkMaxDensity, d);
                    walkFrames++;
                    walkGroundedFrames += player.grounded ? 1 : 0;
                    { const DVec3 here = start->toLocalD(player.position); walkPath += (here - walkLast).length(); walkLast = here; }
                    if (GetTime() - walkStart > benchWalkSeconds || (walkTargetM > 0.0 && walkPath >= walkTargetM)) {
                        player.autoWalk = false;
                        const double dist = (start->toLocalD(player.position) - walkFrom).length();
                        const float alt = Vector3Distance(feet, start->position);
                        printf("bench: walked %.1f m in %.1f s, grounded %.0f%% of frames, density at body %.2f..%.2f (solid >= %.2f), %.1f m from centre\n",
                               dist, GetTime() - walkStart, 100.0 * walkGroundedFrames / std::max(1, walkFrames), walkMinDensity, walkMaxDensity, ISO_LEVEL, alt);
                        if (walkTargetM > 0.0) {
                            const double pct = 100.0 * walkGroundedFrames / std::max(1, walkFrames);
                            const bool pass = walkPath >= walkTargetM && pct >= 90.0 && walkMaxDensity < ISO_LEVEL;
                            printf("bench: planet walk %.0f m of %.0f m (path) in %.1f s, grounded %.0f%%, deepest body density %.2f (solid >= %.2f): %s\n",
                                   walkPath, walkTargetM, GetTime() - walkStart, pct, walkMaxDensity, ISO_LEVEL, pass ? "PASS" : "FAIL");
                            walkTargetM = 0.0;
                        }
                        benchWalkSeconds = 0.0;
                    }
                } else if ((WorldIdle() && benchDigStrokes == 0 && !benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && minePhase != 2 && benchFlySeconds <= 0.0 && !benchSuffocate && benchOrbitSeconds <= 0.0 && !benchRide && !benchPark && !benchGravity && !benchPrecision && !benchMantle && !benchAirJump && !benchMoon && !benchHole && !benchWall && !benchMenus && !benchWater && !benchFlyFollow && !benchHotbar && !benchRing && !benchRocket && !benchTrees && !benchChop && benchPersist == 0 && benchBuild == 0 && benchFactory == 0 && !benchProgression) || GetTime() > 600.0) {
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
            if ((benchFlySpeed > 0.0f || benchFlyAlt < 1e8f) && benchStart >= 0.0) {
                // --bench-fly V: the camera flies a great circle over body 0 at V m/s, 150 m above the
                // terrain, looking ahead and down; node counts each second show whether LOD keeps up.
                Planetoid* b = dynamic_cast<Planetoid*>(world.objects[0].get());
                const auto* pf = b ? dynamic_cast<const PlanetField*>(b->field.get()) : nullptr;
                if (pf) {
                    const double ang = benchFlySpeed * t / (double)pf->radius;
                    const Vector3 dir = { (float)cos(ang) * 0.8f, 0.6f, (float)sin(ang) * 0.8f };       // body-local, 37 deg north
                    const Vector3 ahead = { -(float)sin(ang) * 0.8f, 0.0f, (float)cos(ang) * 0.8f };
                    const float alt = benchFlyAlt < 1e8f ? pf->radius + benchFlyAlt : pf->radius + std::max(0.0f, pf->elevation(dir)) + 150.0f;
                    const DVec3 pos = b->toWorldD(DVec3(Vector3Scale(dir, alt)));
                    const float look = benchFlyAlt < 0.0f ? alt + 80.0f : alt - 60.0f; // underwater: look up at the surface
                    const DVec3 tgt = b->toWorldD(DVec3(Vector3Add(Vector3Scale(dir, look), Vector3Scale(Vector3Normalize(ahead), 200.0f))));
                    world.setCameraWorld(pos, tgt);
                    world.camera.up = b->dirToWorld(dir);
                    if (benchFlyAlt > 5000.0f && benchFlyAlt < 1e8f) { // high up: look down at the planet, a little off-centre
                        world.setCameraWorld(pos, b->toWorldD(DVec3(Vector3Add(Vector3Scale(dir, pf->radius * 0.5f), Vector3Scale(Vector3Normalize(ahead), pf->radius * 0.3f)))));
                        world.camera.up = b->dirToWorld(Vector3Normalize(ahead));
                    }
                    static int lastSec = -1;
                    if ((int)t != lastSec) {
                        lastSec = (int)t;
                        const WorldStats& s = GetWorldStats();
                        printf("bench: fly t=%2d s  %5.1f km flown  nodes %5d  drawn %4d  queued %4d  %.0f FPS\n", (int)t,
                               benchFlySpeed * t / 1000.0, s.nodes, s.chunks, s.queuedJobs, (double)GetFPS());
                    }
                }
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
        if (!bench && !telescope.active && !DebugUI_IsOpen()) world.camera.fovy = settings.fov; // Settings (the telescope zooms from it)
        // Near plane: 0.5 m keeps depth precision out to the 2 km pass split, but the frustum corners
        // then sit 0.65 m from the eye, past the 0.45 m player radius, so walls clipped into view.
        // Close to terrain or a build it drops to 0.2 m (looking at a wall, far precision matters less).
        {
            float nearPlane = 0.5f;
            if (Planetoid* nb = world.objects.empty() ? nullptr : NearestPlanetoid(world, world.camera.position)) { // (player.nearest can be stale right after a load)
                float dist; Vector3 n;
                if (TerrainSurfaceDistance(*nb, world.camera.position, dist, n) && dist < 2.0f) nearPlane = 0.2f;
                if (nearPlane > 0.2f && nb->plots)
                    for (auto& pl : nb->plots->plots) {
                        Vector3 c = nb->toLocal(world.camera.position), nl;
                        if (pl->collideSphere(c, 2.0f / nb->scale, nl)) { nearPlane = 0.2f; break; }
                    }
            }
            rlSetClipPlanes(nearPlane, 60000.0);
        }
        const bool gameInput = !bench && !uiOpen && gameState == 1;
        world.cameraInputEnabled = player.flyMode && !uiOpen && !ship.piloted;  // free-fly camera only in fly mode
        const bool alive = !survival.dead;
        if (gameState == 1 && gameInput && alive && !player.flyMode && !Rocket_Boarded() && !Rocket_AimingCockpit() && IsKeyPressed(KEY_E)) { // E rolls in free-fly
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
        // One clamped frame time for everything that moves (bodies, ship, player): clamping only the
        // player let the ground spin away under it during a slow frame.
        const float frameDt = std::min(GetFrameTime(), 0.1f);
        advanceBodies(world, gameState == 1 ? frameDt * timeScale : 0.0f); // bodies only spin while playing; before riders update
        if (gameState == 1) {
            ship.lastImpact = 0.0f;
            if (ship.piloted) {
                ship.update(world, gameInput && alive, frameDt);
                ship.applyChaseCamera(world);
                player.position = ship.position; // the pilot rides inside
                player.velocity = ship.velocity;
            } else if (Rocket_Boarded()) {
                ship.update(world, false, frameDt);
                Rocket_Update(world, gameInput && alive, frameDt * timeScale); // flight (and its camera)
            } else {
                ship.update(world, false, frameDt);
                if (alive) player.update(world, gameInput, frameDt);
                Survival_OnLanding(player.impactSpeed);
                Rocket_Update(world, gameInput && alive, frameDt); // on foot: aim at a rocket (E boards, I fuels)
            }
            if (ship.piloted && ship.lastImpact > 12.0f) survival.health -= (ship.lastImpact - 12.0f) * 3.0f; // crash damage
            Survival_Update(world, frameDt, gameInput);
            // Autosave every 5 minutes of play (incremental: only regions changed since the last save).
            static double lastAutosave = GetTime();
            if (!bench && alive && GetTime() - lastAutosave > 300.0) {
                lastAutosave = GetTime();
                if (SaveWorld(world)) GameUI_Toast("Autosaved");
            }
            Factory::Sync(world);                                // machine layouts in, items out
            Machine_Update(world, gameInput && alive && !ship.piloted);
            Objectives_Update(world);
            if (bench && lookDown && !ship.piloted) {
                const Vector3 t = Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)));
                player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(t), 0.8f), player.up)); // ~50 deg down
            }
            if (bench && benchTelescope > 0.0f && !ship.piloted) AimAtTelescopeTarget(world);
            if (!ship.piloted && !Rocket_Boarded()) player.applyToCamera(world);
            Telescope_Update(world, gameInput && alive, !ship.piloted && !player.flyMode && alive && !buildTool.active);
            if (bench && benchTelescope > 0.0f) { Telescope_SetActive(world, true); telescope.zoom = benchTelescope; Telescope_Update(world, false, !ship.piloted); }
        }

        world.updateFloatingOrigin(); // keep the camera near render-space 0 (float precision)
        { const double w0 = GetTime(); worldUpdate(world, 3.0f); if (benchMeasuring) benchUpdateMs += (GetTime() - w0) * 1000.0; } // LOD + streaming (bodies were spun at the start of the frame)

        // Slot 1: headlamp (suit lamp from the eye, or the ship's nose light while piloting).
        // Slot 2: soft glow around the player while the inventory is open.
        if (sunOverride) world.sunDirection = sunArg;
        // Atmosphere: the nearest planet (moons and asteroids have none) colours the sky and haze.
        world.atmosphere.enabled = false;
        {
            float best = 1e30f;
            for (auto& obj : world.objects)
                if (auto* q = dynamic_cast<Planetoid*>(obj.get()))
                    if (q->kind == BodyKind::Planet) {
                        const float alt = Vector3Distance(world.camera.position, q->position) - q->radius * q->scale;
                        if (alt < best) {
                            best = alt;
                            world.atmosphere.enabled = true;
                            world.atmosphere.centre = q->position;
                            world.atmosphere.radius = q->radius * q->scale;
                        }
                    }
        }
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
        if (gameState == 1) Hotbar_Update(world, gameInput && alive, ship.piloted);
        if (gameState == 1) BuildTool_Update(world, gameInput && alive, !ship.piloted && !Rocket_Boarded() && !player.flyMode && !telescope.active && alive);
        if (gameState == 1) TerrainTool_Update(world, gameInput && alive, !ship.piloted && !Rocket_Boarded() && !telescope.active && !buildTool.active); // F5/F9 work aboard too

        // --test-menus: the menu flow without clicks: main menu (screenshot), New Game, pause
        // (screenshot), Save & Quit to Menu, Load Game; each step's state is checked.
        if (benchMenus) {
            ++menuFrames;
            if (menuStep == 0) { gameState = 0; mainMenu.isActive = true; world.isActive = false; menuStep = 1; menuFrames = 0; }
            else if (menuStep == 1 && menuFrames == 30) { TakeScreenshot("bench_menu_main.png"); pendingMenu.action = MenuAction::NewGame; menuStep = 2; menuFrames = 0; }
            else if (menuStep == 2 && menuFrames > 60 && WorldIdle()) {
                menuOk &= gameState == 1;
                printf("bench: menus: New Game -> in game: %s\n", gameState == 1 ? "yes" : "NO");
                Menu_ResetPause(); gameState = 2; justPaused = true; menuStep = 3; menuFrames = 0;
            } else if (menuStep == 3 && menuFrames == 30) { TakeScreenshot("bench_menu_pause.png"); pendingMenu.action = MenuAction::SaveAndQuit; menuStep = 4; menuFrames = 0; }
            else if (menuStep == 4 && menuFrames == 5) {
                menuOk &= gameState == 0;
                printf("bench: menus: Save & Quit -> main menu: %s\n", gameState == 0 ? "yes" : "NO");
                pendingMenu.action = MenuAction::LoadGame; pendingMenu.path = "saves/bench.sav"; menuStep = 5; menuFrames = 0;
            } else if (menuStep == 5 && menuFrames == 5) {
                menuOk &= gameState == 1;
                printf("bench: menus: Load Game -> in game: %s -> %s\n", gameState == 1 ? "yes" : "NO", menuOk ? "PASS" : "FAIL");
                benchMenus = false;
            }
        }
        // Menus: the result of last frame's clicks (the menus draw, and see clicks, in BeginDrawing).
        {
            auto enterGame = [&] {
                HideCursor();
                gameState = 1;
                mainMenu.isActive = false;
                world.isActive = true;
            };
            auto toMenu = [&] {
                ShowCursor();
                gameState = 0;
                mainMenu.isActive = true;
                world.isActive = false;
            };
            const MenuResult m = pendingMenu;
            pendingMenu = MenuResult{};
            switch (m.action) {
                case MenuAction::NewGame:
                    SetWorldSavePath(bench ? "saves/bench.sav" : Menu_NewSavePath().c_str()); // bench runs never touch the player's saves
                    if (!worldFresh) { worldGen.seed = 0; regenerateWorld(world); StartNewWorld(world, true); }
                    worldFresh = false;
                    enterGame();
                    break;
                case MenuAction::LoadGame:
                    SetWorldSavePath(m.path.c_str());
                    if (LoadWorldSave(world)) { worldFresh = false; enterGame(); }
                    else GameUI_Toast(TextFormat("Could not load %s", m.path.c_str()));
                    break;
                case MenuAction::Resume: enterGame(); break;
                case MenuAction::Save: GameUI_Toast(SaveWorld(world) ? "Saved" : "Save failed"); break;
                case MenuAction::SaveAndQuit:
                    if (!SaveWorld(world)) GameUI_Toast("Save failed");
                    toMenu();
                    break;
                case MenuAction::Exit: windowKeepAlive = false; break;
                case MenuAction::None: break;
            }
            // Esc in game pauses (unless it closes the inventory or build mode first); the pause menu
            // handles Esc itself (resume, or back out of Settings).
            if (gameState == 1 && m.action == MenuAction::None && IsKeyPressed(KEY_ESCAPE)) {
                if (GameUI_InventoryOpen()) { GameUI_ToggleInventory(); SetMousePosition(GetScreenWidth() / 2, GetScreenHeight() / 2); }
                else if (buildTool.active) buildTool.active = false;
                else { ShowCursor(); Menu_ResetPause(); gameState = 2; justPaused = true; }
            }
        }

        

        BeginDrawing();
            ClearBackground(crackCheck ? MAGENTA : BLACK); // --crack-check: holes in the terrain show magenta
            if (gameState == 0) {
                pendingMenu = Menu_Main();
            } else if (gameState == 1) {
                { const double d0 = GetTime(); world.drawScene(gameState); if (benchMeasuring) benchDrawMs += (GetTime() - d0) * 1000.0; }
                player.recordFly(world); // after the fly camera's input (in drawScene): the pose to carry next frame
                world.drawUI(gameState);
                Telescope_DrawHUD(world);
                BuildTool_DrawHUD();
                Hotbar_DrawHUD(ship.piloted);
                Machine_DrawHUD();
                Rocket_DrawHUD(world);
                if (!ship.piloted && !telescope.active && !buildTool.active && !Rocket_Boarded()) TerrainTool_DrawHUD();
                GameUI_DrawHUD(world);
                Survival_DrawHUD();
                Objectives_DrawHUD();
                if (!ship.piloted && ship.playerInRange(world, (player.position - world.renderOrigin).toFloat(), 5.0f))
                    DrawText("[E] Board ship", GetScreenWidth() / 2 - 60, GetScreenHeight() / 2 + 40, 20, SKYBLUE);
                if (ship.piloted) DrawText("SHIP: WASD thrust  Space/Ctrl up/down  mouse steer  Z/C roll  Shift boost  X brake  L lights  0 power cell (hotbar)  E exit",
                                           10, GetScreenHeight() - 30, 18, Fade(SKYBLUE, 0.9f));
                if (!bench && !uiOpen) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
            } else if (gameState == 2) {
                world.drawScene(gameState);
                world.drawUI(gameState);
                const MenuResult m = Menu_Pause();
                if (!justPaused) pendingMenu = m; // the Esc that opened it is not also a "resume"
            }
            justPaused = false;
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
    Factory::Stop();
    CloseWindow();
}


