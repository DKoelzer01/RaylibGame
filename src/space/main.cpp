#include <set>
#include <map>
#include <cstring>
#include <raylib.h>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <array>
#include "terrainarena.h"
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
#include "gputimer.h"
#include "ship.h"
#include "camera.h"
#include "survival.h"
#include "telescope.h"
#include "glcaps.h"
#include "buildtool.h"
#include "progression_test.h"
#include "menus.h"
#include "hotbar.h"
#include "itemicons.h"
#include "particles.h"
#include "ambience.h"
#include "machinemodels.h"
#include "gpucompute.h"
#include "terraingpu.h"
#ifdef _WIN32
extern "C" __declspec(dllimport) int __stdcall MessageBoxA(void*, const char*, const char*, unsigned int); // the GL 4.3 error box (no windows.h: it clashes with raylib)
#endif
#include "instancebatch.h"
#include "navigation.h"
#include "viewmodel.h"
#include "farming.h"
#include "animals.h"
#include "suit.h"
#include "toolwear.h"
#include "recipebook.h"
#include "station.h"
#include "weather.h"
#include "storage.h"
#include "tools.h"
#include "rivers.h"
#include "questtree.h"
#include "ring_test.h"
#include "rocket.h"
#include "jobsystem.h"
#include "rocket_test.h"
#include "blueprint.h"
#include "ring.h"
#include "factory.h"
#include "fluid.h"
#include "waterworks.h"

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

// A new game on the world as generated: the player at the start, a
// clearing (no tree on the spot), the starting kit. `reset` also clears what a previous game left behind
// (inventory, quest tier, survival state).
static void StartNewWorld(Scene& world, bool reset) {
    if (reset) {
        inventory.clear();
        QuestBook_Reset();
        survival = SurvivalState{};
        player.velocity = { 0, 0, 0 };
        player.flyMode = false;
        ship.piloted = false;
        Rocket_Reset();
        Nav_Reset();
        Farm_Reset();
        Animals_Reset();
        Particles_Clear();
        Ambience_Reset();
        Station_Reset();
        Suit_Reset();
        ToolWear_Reset();
        RecipeBook_Reset();
    }
    Planetoid* start = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    if (start) {
        player.spawnOn(world, *start, StartSpawnDir(*start));
        Nav_SetStart(start->name, StartSpawnDir(*start));
        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
        // A clearing at the start, so the player does not start in a tree.
        if (start->trees) start->trees->addClearing(start->toLocal(feet), 14.0f);
    }
    // The starting kit by difficulty (a save restores its own inventory and difficulty).
    if (reset) worldDifficulty = settings.difficulty;
    if (inventory.count("telescope") == 0) inventory.add("telescope", 1);
    if (inventory.count("map") == 0 && worldDifficulty < 2) inventory.add("map", 1); // hard: find your own way
    if (reset && worldDifficulty == 0) inventory.add("berries", 6);
}

// --bench <seconds>: windowed 1920x1080, uncapped FPS, skips the menu, holds a fixed
// camera on the planetoid, prints frame stats to stdout and bench.txt, then exits.


int main(int argc, char** argv)
{
    // Data and shaders are read from resources/ relative to the working directory. Launched from Explorer
    // (or a build dir) that is the exe's folder, which has no resources/ next to it: the game found nothing and
    // crashed (#234). Use the first of the working directory, the exe's folder and its parents that has it.
    if (!DirectoryExists("resources")) {
        std::string dir = GetApplicationDirectory();
        for (int up = 0; up < 6 && !dir.empty(); ++up) {
            if (DirectoryExists((dir + "resources").c_str())) { ChangeDirectory(dir.c_str()); break; }
            const size_t cut = dir.find_last_of("/\\", dir.size() >= 2 ? dir.size() - 2 : std::string::npos);
            if (cut == std::string::npos) break;
            dir = dir.substr(0, cut + 1);
        }
    }
    double benchSeconds = 0.0;
    bool startWithOverlay = false;
    int benchDigStrokes = 0;     // --test-dig N: dig N strokes at the crosshair once the world is loaded
    bool benchPaint = false;     // --test-paint: placed gravel stays gravel (and through a save and load)
    bool benchSaveLoad = false;  // --test-saveload: then save, reload from disk, and bench the reloaded world
    double benchWalkSeconds = 0.0; // --test-walk S: spawn on the surface and auto-walk forward for S seconds
    bool benchWalkMode = false;    // --bench-walk: the player walks forward for the whole bench (frame time while streaming)
    bool benchStandMode = false;   // --bench-stand: as --bench-walk, but standing still at the start (#204: walking vs standing, same view)
    bool cursorLocked = false;     // DisableCursor while playing (#205)
    bool perfLog = false;          // --perf-log: per-pass CPU/GPU averages to log.txt every 5 s in normal play, tagged moving/still (#204)
    float benchFlySpeed = 0.0f;    // --bench-fly V: the camera flies over body 0 at V m/s through the bench (LOD streaming)
    bool benchFlyBack = false;     // --bench-fly-back: fly out for half the bench, then back the same way (chunk cache, #295)
    float benchFlyAlt = 1e9f;      // --bench-fly-alt A: fly at A m relative to sea level instead (negative: underwater, looking up)
    bool hasSpawnDir = false;      // --spawn-dir x y z: where --test-walk / --test-stand spawn on body 0
    Vector3 spawnDir = { 0, 1, 0 };
    double walkStart = -1.0;
    bool benchStandOnly = false; // --test-stand S: like --test-walk but standing still
    int walkJolts = 0, walkJoltFrames = 0, walkSteepFrames = 0, walkSteepJolts = 0; double walkJoltMax = 0.0, walkSlopeMax = 0.0;
    int walkUnsticks = 0; // --test-planet-walk: times the walker had to jump or turn to get unstuck
    double standCamMax = 0.0; float standFwdMax = 0.0f, standUpMax = 0.0f, standTexelSlide = 0.0f; int standFrames = 0, standJitterFrames = 0;
    int benchMineStrokes = 0;     // --test-mine N: dig N strokes at the feet, craft, build (after any walk test)
    bool mineStarted = false, mineStoppedEarly = false;
    int minePhase = 0, minedStrokes = 0, mineCrafted = 0;
    bool benchWaterStats = false; // --water-stats: how much the water shaping raises/cuts the ground (#251)
    bool benchPreloadSpawn = false; double preloadT0 = -1.0; // --preload-spawn: time the spawn preload ring (#294)
    double benchFlySeconds = 0.0, flyStart = -1.0; // --test-fly S: autopilot the ship to the nearest other asteroid
    bool benchSuffocate = false;   // --test-suffocate: run out of oxygen away from the ship, die, respawn
    double benchOrbitSeconds = 0.0, orbitStart = -1.0; // --test-orbit S: coast on a circular orbit for S seconds
    float orbitMin = 0, orbitMax = 0;
    bool orbitCrashed = false;
    bool benchRide = false;        // --test-ride: stand on a spinning body, then jump
    bool benchPark = false;        // --test-park: park the ship on a spinning body
    bool sunOverride = false;      // --sun x y z: sun direction (e.g. to test night-side lighting)
    bool benchDay = false;         // --day: a mid-morning sun over the start (perf benches: shadows and all, #291)
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
    bool benchInventory = false;   // --test-inventory: grid stacking, fractions, sorting, saves
    bool benchMap = false;         // --test-map: compass, exploration, waypoints, map images, saves
    bool benchQuests = false;      // --test-quests: tiers open by their tasks; quest state through a save
    bool benchHunger = false, benchFarm = false; // --test-hunger, --test-farm
    bool benchRivers = false, riverOk = false; int riverPhase = 0; double riverStart = 0.0; float riverLongest = 0.0f, riverBuild = 0.0f, riverWet = 0.0f, riverLevel = 0.0f; Vector3 riverDir{ 0, 0, 0 }; // --test-rivers
    bool benchTools = false, toolsOk = true; int toolsPhase = 0, toolsLife = 0, toolsTilled = 0; double toolsStart = 0.0; // --test-tools
    bool benchSuit = false, suitOk = true; int suitPhase = 0; double suitStart = 0.0; float suitFoodRatio = 0.0f, suitWater = 0.0f, suitRebreather = 0.0f; // --test-suit
    int benchShotTabs = 0; // --shot-tabs
    bool benchBiomes = false; int benchShotBiomes = -1; // --test-biomes, --shot-biomes
    bool benchOil = false; int oilPhase = 0; double oilT0 = 0.0; // --test-oil
    bool benchSpring = false; int springPhase = 0; double springT0 = 0.0; // --test-spring
    bool benchShotAnimals = false; int animalShotPhase = 0; // --shot-animals
    bool benchShotMachines = false; int machineShotPhase = 0; // --shot-machines
    bool benchParticles = false; int particlePhase = 0; double particleT0 = 0.0; Vector3 particleEye{}, particleAt{}; // --shot-particles
    bool benchShadow = false; int shadowPhase = 0; double shadowT0 = 0.0; Vector3 shadowEye{}, shadowAt{}; // --shot-shadow
    bool benchFarForest = false; int farForestPhase = 0; double farForestT0 = 0.0; Vector3 farForestEye{}, farForestAt{}; // --shot-far-forest
    bool benchTerrain = false; int terrainPhase = 0; double terrainT0 = 0.0; Vector3 terrainAt{}, terrainOut{}; // --shot-terrain
    bool benchPump = false; int pumpPhase = 0, pumpCable = -1; double pumpT0 = 0.0, pumpedAtCut = 0.0; uint32_t pumpPlot = 0; Vector3 pumpBasin{}; float pumpRise = 0.0f; // --test-pump
    bool benchBucket = false; int bucketPhase = 0, bucketTrips = 0; double bucketT0 = 0.0; Vector3 bucketWater{}, bucketBasin{}, bucketField{}; // --test-bucket
    bool benchPond = false; int pondPhase = 0; double pondT0 = 0.0; Vector3 pondDirB{}, pondOutB{}, pondPit{}, pondMid{}; float pondLevelB = 0.0f, pondRB = 0.0f; // --test-pond
    bool benchFlow = false; int flowPhase = 0; double flowStart = 0.0, flowEdit = 0.0, flowReachT = -1.0, flowCaveT = -1.0, flowSleepT = -1.0; Vector3 flowDir{}, flowOut{}, flowEnd{}, flowCave{}; float flowLevel = 0.0f; // --test-flow
    bool benchForest = false; int forestPhase = 0; double forestStart = 0.0; Vector3 forestCam{}, forestTarget{}, forestUp{}; bool outcropHigh = false; // --shot-forest (--outcrop-high: the highest outcrop in range, not the nearest)
    bool benchRecipes = false; // --test-recipes
    bool benchWeather = false, weatherOk = true; int weatherPhase = 0, weatherCounts[6] = {}, weatherChanges = 0; double weatherStart = 0.0; float weatherChill = 0.0f; // --test-weather
    bool benchDecor = false, decorOk = true; int decorPhase = 0, decorGems = 0, decorGemKinds = 0, decorDigs = 0, decorPlaced = 0; double decorStart = 0.0; Vector3 decorCamPos{}, decorCamTarget{}, decorCamUp{ 0, 1, 0 }; // --test-decor
    bool benchMaterials = false, matOk = true; int matPhase = 0, matPlaced = 0; double matStart = 0.0; std::string matKinds; Vector3 matCamPos{}, matCamTarget{}, matCamUp{ 0, 1, 0 }; // --test-materials
    bool benchDamage = false, damageOk = true; int damagePhase = 0; double damageStart = 0.0; float damageThrust = 0.0f, damageRepaired = 0.0f, damageVacuum = 0.0f; // --test-damage
    bool benchCrafting = false, craftOk = true; int craftPhase = 0, craftChest = -1, craftInv0 = 0, craftWood0 = 0, craftDeposited = 0, craftInChest = 0; double craftStart = 0.0; Plot* craftPlot = nullptr; // --test-crafting
    bool benchThreats = false, threatsOk = true; int threatPhase = 0; double threatStart = 0.0; float threatBite = 0.0f, threatFire = 0.0f, threatRam = 0.0f; int threatSpawned = 0, threatLeft = 0, threatPeaceful = 0; // --test-threats
    bool benchAnimals = false; int animalPhase = 0, animalCount = 0; double animalStart = 0.0, animalMs = 0.0; int animalNatural = 0; float penWorst = 0.0f, animalLured = 0.0f, animalSpread = 0.0f, animalFled = 0.0f, animalFollow = 0.0f; // --test-animals
    int hungerPhase = 0, farmPhase = 0; double hungerStart = 0.0, farmStart = 0.0; bool hungerOk = false; float hungerTemp = 0.0f, hungerDrain = 0.0f, farmQuarter = 0.0f, farmGrowth = 0.0f;
    bool benchViewmodel = false;   // --test-viewmodel: hands still at rest, bob within limits walking, screenshots
    int vmPhase = 0, vmSigns = 0; double vmStart = 0.0; float vmStillMax = 0.0f, vmMaxX = 0.0f, vmMaxY = 0.0f, vmWalked = 0.0f; Vector3 vmLast{}; DVec3 vmFrom{};
    int mapPhase = 0, mapTiles = 0; double mapStart = 0.0; float mapHeading = -1.0f, mapBearing = -1.0f; bool mapImageOk = false; DVec3 mapFrom{};
    bool benchRing = false, ringStarted = false; // --test-ring: milestone 6 gate
    bool benchRocket = false, rocketStarted = false; // --test-rocket: milestone 7 gate
    bool benchStation = false, stationStarted = false; double stationWait = 0.0; // --test-station: plan 3 M20 gate
    double rocketWait = 0.0;
    bool benchBlueprint = false;   // --test-blueprint: place, fill in, save and load a blueprint
    bool benchSmallBody = false;   // --test-smallbody: land on a huge and a small ring rock
    int sbPhase = 0, sbRock = -1, sbGrounded = 0, sbFrames = 0, sbDeaths = 0;
    double sbStart = 0.0;
    float sbMaxImpact = 0.0f, sbHealth = 100.0f;
    bool sbDropped = false, sbOk = true;
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
    int wallPhase = 0; bool wallOk = true;
    double wallStart = 0.0, wallMaxStep = 0.0, wallMaxStepAt = 0.0;
    float wallMaxImpact = 0.0f, wallMinHealth = 100.0f, wallMaxSpeed = 0.0f; int wallCamInside = 0; // frames with the camera inside terrain (#255)
    int wallCamShake = 0; double wallCamShakeMax = 0.0; DVec3 wallCamPrev{}, wallCamStepPrev{}; int wallCamKnown = 0, wallCamShakePushed = 0, wallCamShakeGrounded = 0; // camera shake (#255)
    bool benchLedge = false; int ledgePhase = 0, ledgePass = 0; bool ledgeOk = true, ledgeAir = false; double ledgeStart = 0.0, ledgeAirStart = 0.0; float ledgeV0 = 0.0f, ledgeExcess = 0.0f; float ledgeWalk0 = 0.0f;
    bool benchTor = false; int torPhase = 0, torHeading = 0, torFlicker = 0, torPops = 0; double torStart = 0.0, torAirAt = -1.0, torMaxStep = 0.0; float torMaxUp = 0.0f; Vector3 torDir{}; DVec3 torLast{}; bool torWasGrounded = true; // --test-tor (#242)
    bool benchMachines = false; int machPhase = 0; double machStart = 0.0; // --test-machines (#193 #194 #229)
    bool benchFall = false; int fallPhase = 0, fallIdx = 0; bool fallOk = true; double fallStart = 0.0; // --test-fall (#263)
    bool benchTeleClose = false; int teleCPhase = 0; double teleCStart = 0.0, teleCWorst = 0.0; long long teleCReq0 = 0; // --test-telescope-close (#261)
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
    bool progressionFull = false;  // --test-progression full: milestone 8 gate (on to titanium alloy)
    bool progressionEv = false;    // --test-progression ev: M24 gate (tiers 9-10 from a moon and ring haul)
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
    double moonStart = -1.0; int moonPhase = 0; DVec3 moonLocal0{}, moonTakeoff{}; bool moonAirborne = false; float moonMaxUp = 0.0f; double moonMaxRise = 0.0, moonMinRise = 0.0; int moonStrokes = 0;
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
    bool benchAtSave = false;      // --at-save
    bool spinTest = false;         // --spin-test: textures must stay attached to a spinning body
    float spinRate = 30.0f;        // --spin-rate DEG_PER_SEC
    int spinPhase = 0;
    Vector3 spinCamLocal{}, spinTgtLocal{};
    double suffocateStart = -1.0;
    double suffocateDrift = 0.0;    // --test-suffocate: how far the dead camera drifted in the body frame (#231)
    double suffocateSettle = 0.0;   // ...and how far the eye smoothing moved it in the first 0.2 s
    DVec3 flyFrom{};
    float flyEnergy = 0.0f, flyMaxSpeed = 0.0f;
    DVec3 walkFrom{};
    float walkMinDensity = 0, walkMaxDensity = 0;
    int walkGroundedFrames = 0, walkFrames = 0;
    double walkTargetM = 0.0, walkPath = 0.0; // --test-planet-walk M: walk M metres (path, body frame), then PASS/FAIL
    DVec3 walkLast{};
    Vector3 benchCamOffset = { -70.0f, 35.0f, 70.0f }; // --cam dx dy dz: bench camera offset from the planetoid
    bool benchCamGiven = false;   // --cam given: use the offset (otherwise a planet bench looks along the surface from 40 m up, #291)
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--bench") benchSeconds = (i + 1 < argc) ? atof(argv[i + 1]) : 10.0;
        if (std::string(argv[i]) == "--ui") DebugUI_Toggle(); // open the debug panel at start
        if (std::string(argv[i]) == "--seed" && i + 1 < argc) worldGen.seed = (unsigned int)strtoul(argv[i + 1], nullptr, 10);
        if (std::string(argv[i]) == "--nocull") frustumCullingEnabled = false;
        if (std::string(argv[i]) == "--overlay") startWithOverlay = true;
        if (std::string(argv[i]) == "--cam" && i + 3 < argc) { benchCamGiven = true; benchCamOffset = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
        if (std::string(argv[i]) == "--test-dig" && i + 1 < argc) benchDigStrokes = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-paint") benchPaint = true;
        if (std::string(argv[i]) == "--test-saveload") benchSaveLoad = true;
        if (std::string(argv[i]) == "--test-walk" && i + 1 < argc) benchWalkSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--no-eye-smooth") player.eyeSmoothing = false;
        if (std::string(argv[i]) == "--test-planet-walk" && i + 1 < argc) { walkTargetM = atof(argv[i + 1]); benchWalkSeconds = 600.0; } // gate: 1 km on the planet
        if (std::string(argv[i]) == "--bench-walk") benchWalkMode = true;
        if (std::string(argv[i]) == "--bench-stand") { benchWalkMode = true; benchStandMode = true; }
        if (std::string(argv[i]) == "--perf-log") perfLog = true;
        if (std::string(argv[i]) == "--bench-fly" && i + 1 < argc) benchFlySpeed = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--bench-fly-back") benchFlyBack = true;
        if (std::string(argv[i]) == "--water-stats") benchWaterStats = true;
        if (std::string(argv[i]) == "--bench-pond") gSpawnPond = true; // start at the nearest pond's rim (#251 crack checks)
        if (std::string(argv[i]) == "--bench-biome" && i + 1 < argc) { // start the bench (and the player) in the nearest deep patch of a biome (#133 checks)
            for (int b = 0; b < (int)Biome::Count; ++b) if (!strcmp(BiomeName((Biome)b), argv[i + 1])) gSpawnBiome = b;
            if (gSpawnBiome < 0) printf("bench: --bench-biome: unknown biome %s\n", argv[i + 1]);
        }
        if (std::string(argv[i]) == "--bench-fly-alt" && i + 1 < argc) benchFlyAlt = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--spawn-dir" && i + 3 < argc) { hasSpawnDir = true; spawnDir = Vector3Normalize({ (float)atof(argv[i + 1]), (float)atof(argv[i + 2]), (float)atof(argv[i + 3]) }); }
        if (std::string(argv[i]) == "--test-suffocate") benchSuffocate = true;
        if (std::string(argv[i]) == "--test-orbit" && i + 1 < argc) benchOrbitSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-ride") benchRide = true;
        if (std::string(argv[i]) == "--dump-system") dumpSystem = true;
        if (std::string(argv[i]) == "--load") loadSaveAtStart = true;
        if (std::string(argv[i]) == "--at-save") { loadSaveAtStart = true; benchAtSave = true; } // --load, and bench from the save's own camera (no respawn)
        if (std::string(argv[i]) == "--ride-spin" && i + 1 < argc) rideSpin = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-park" || std::string(argv[i]) == "--test-fly") { printf("bench: %s was removed with the placeholder ship\n", argv[i]); if (i + 1 < argc && argv[i + 1][0] != '-') ++i; }
        if (std::string(argv[i]) == "--test-gravity") benchGravity = true;
        if (std::string(argv[i]) == "--test-precision") benchPrecision = true;
        if (std::string(argv[i]) == "--test-mantle") benchMantle = true;
        if (std::string(argv[i]) == "--test-airjump") benchAirJump = true;
        if (std::string(argv[i]) == "--test-moon") benchMoon = true;
        if (std::string(argv[i]) == "--test-hole") benchHole = true;
        if (std::string(argv[i]) == "--test-wall") benchWall = true;
        if (std::string(argv[i]) == "--test-ledge") benchLedge = true; // walk and sprint off a pit edge (#240)
        if (std::string(argv[i]) == "--test-tor") benchTor = true; // walk across a low (beach) tor crown, logging collision pops (#242)
        if (std::string(argv[i]) == "--test-machines") benchMachines = true; // fuel slot, crafts consumed at the end, the boiler's water and heat
        if (std::string(argv[i]) == "--test-fall") benchFall = true; // drops from 3-100 m: fall damage by height
        if (std::string(argv[i]) == "--test-telescope-close") benchTeleClose = true; // zoom 16x at the horizon, close: worst frame after (#261)
        if (std::string(argv[i]) == "--test-menus") benchMenus = true;
        if (std::string(argv[i]) == "--test-water") benchWater = true;
        if (std::string(argv[i]) == "--test-flyfollow") benchFlyFollow = true;
        if (std::string(argv[i]) == "--test-hotbar") benchHotbar = true;
        if (std::string(argv[i]) == "--test-inventory") benchInventory = true;
        if (std::string(argv[i]) == "--test-map") benchMap = true;
        if (std::string(argv[i]) == "--test-quests") benchQuests = true;
        if (std::string(argv[i]) == "--test-hunger") benchHunger = true;
        if (std::string(argv[i]) == "--test-farm") benchFarm = true;
        if (std::string(argv[i]) == "--test-animals") benchAnimals = true;
        if (std::string(argv[i]) == "--test-threats") benchThreats = true;
        if (std::string(argv[i]) == "--shot-tabs") benchShotTabs = 1;
        if (std::string(argv[i]) == "--shot-forest") benchForest = true;
        if (std::string(argv[i]) == "--outcrop-high") outcropHigh = true;
        if (std::string(argv[i]) == "--test-flow") benchFlow = true;
        if (std::string(argv[i]) == "--test-pond") benchPond = true;
        if (std::string(argv[i]) == "--test-bucket") benchBucket = true;
        if (std::string(argv[i]) == "--test-pump") benchPump = true;
        if (std::string(argv[i]) == "--shot-terrain") benchTerrain = true;
        if (std::string(argv[i]) == "--shot-shadow") benchShadow = true;
        if (std::string(argv[i]) == "--shot-particles") benchParticles = true;
        if (std::string(argv[i]) == "--test-biomes") benchBiomes = true;
        if (std::string(argv[i]) == "--test-oil") benchOil = true;
        if (std::string(argv[i]) == "--test-spring") benchSpring = true;
        if (std::string(argv[i]) == "--shot-animals") benchShotAnimals = true;
        if (std::string(argv[i]) == "--shot-machines") benchShotMachines = true;
        if (std::string(argv[i]) == "--shot-biomes") benchShotBiomes = 0;
        if (std::string(argv[i]) == "--shot-far-forest") benchFarForest = true;
        if (std::string(argv[i]) == "--test-recipes") benchRecipes = true;
        if (std::string(argv[i]) == "--test-weather") benchWeather = true;
        if (std::string(argv[i]) == "--test-decor") benchDecor = true;
        if (std::string(argv[i]) == "--test-materials") benchMaterials = true;
        if (std::string(argv[i]) == "--test-damage") benchDamage = true;
        if (std::string(argv[i]) == "--test-crafting") benchCrafting = true;
        if (std::string(argv[i]) == "--test-suit") benchSuit = true;
        if (std::string(argv[i]) == "--test-tools") benchTools = true;
        if (std::string(argv[i]) == "--test-rivers") benchRivers = true;
        if (std::string(argv[i]) == "--test-viewmodel") benchViewmodel = true;
        if (std::string(argv[i]) == "--test-ring") benchRing = true;
        if (std::string(argv[i]) == "--test-rocket") benchRocket = true;
        if (std::string(argv[i]) == "--test-station") benchStation = true;
        if (std::string(argv[i]) == "--test-blueprint") benchBlueprint = true;
        if (std::string(argv[i]) == "--test-smallbody") benchSmallBody = true;
        if (std::string(argv[i]) == "--test-trees") benchTrees = true;
        if (std::string(argv[i]) == "--test-chop") benchChop = true;
        if (std::string(argv[i]) == "--test-persist" && i + 1 < argc) benchPersist = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-build" && i + 1 < argc) benchBuild = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-factory" && i + 1 < argc) benchFactory = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-progression") { benchProgression = true; progressionFull = i + 1 < argc && std::string(argv[i + 1]) == "full"; progressionEv = i + 1 < argc && std::string(argv[i + 1]) == "ev"; }
        if (std::string(argv[i]) == "--time-scale" && i + 1 < argc) timeScale = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--crack-check") crackCheck = true;
        if (std::string(argv[i]) == "--no-foliage") foliageEnabled = false;
        if (std::string(argv[i]) == "--grass-dist" && i + 1 < argc) kGrassDrawDistance = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--tree-dist" && i + 1 < argc) kTreeDrawDistance = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--sun" && i + 3 < argc) { sunOverride = true; sunArg = Vector3Normalize({ (float)atof(argv[i + 1]), (float)atof(argv[i + 2]), (float)atof(argv[i + 3]) }); }
        if (std::string(argv[i]) == "--day") benchDay = true;
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
        if (std::string(argv[i]) == "--test-mine" && i + 1 < argc) benchMineStrokes = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-stand" && i + 1 < argc) { benchWalkSeconds = atof(argv[i + 1]); benchStandOnly = true; }
        if (std::string(argv[i]) == "--telescope" && i + 1 < argc) benchTelescope = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--telescope-target" && i + 1 < argc) benchTelescopeTarget = argv[i + 1];
        if (std::string(argv[i]) == "--size" && i + 1 < argc) worldGen.size = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--bodies" && i + 1 < argc) worldGen.bodiesPath = argv[i + 1];   // system definition file
        if (std::string(argv[i]) == "--count" && i + 1 < argc) worldGen.count = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--system" && i + 1 < argc) worldGen.systemRadius = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--lod" && i + 1 < argc) lodSettings.detailFactor = (float)atof(argv[i + 1]); // terrain detail (#297)
        if (std::string(argv[i]) == "--lod-fast" && i + 1 < argc) lodSettings.fastFactor = (float)atof(argv[i + 1]); // when flying (99: no drop)
        if (std::string(argv[i]) == "--preload" && i + 1 < argc) lodSettings.preloadRadius = (float)atof(argv[i + 1]); // m, 0 = off (#294)
        if (std::string(argv[i]) == "--cache-verts" && i + 1 < argc) lodSettings.cacheVertices = (int)(atof(argv[i + 1]) * 1048576.0); // millions (#295)
        if (std::string(argv[i]) == "--preload-spawn") benchPreloadSpawn = true; // bench: the spawn preload as a new game does it, timed
        if (std::string(argv[i]) == "--first" && i + 3 < argc) { worldGen.hasFirst = true; worldGen.first = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
    }
    const bool bench = benchSeconds > 0.0;
    benchMode = bench;
    if (bench) SetWorldSavePath("saves/bench.sav"); // tests save and load; keep the player's F5 save out of it
    if (bench && !benchThreats) Animals_SetThreats(false); // hostile animals and night spawns would upset other tests

    bool msaa = true;
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--no-msaa") msaa = false;
    if (msaa) SetConfigFlags(FLAG_MSAA_4X_HINT);
    bool nativeWindow = false; // --native-window: a bench in the game's own window (borderless, monitor-sized), uncapped
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--native-window") nativeWindow = true;
    bool exclusive = false; // --fullscreen: exclusive fullscreen (bypasses the desktop compositor)
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--fullscreen") { exclusive = true; nativeWindow = true; }
    if (!bench || nativeWindow) {
        SetConfigFlags(FLAG_WINDOW_RESIZABLE); // Allow window resizing
        SetConfigFlags(FLAG_WINDOW_UNDECORATED); // Remove window decorations for a cleaner look
        SetConfigFlags(FLAG_WINDOW_HIGHDPI); // Enable High DPI support
        if (exclusive) SetConfigFlags(FLAG_FULLSCREEN_MODE);
        InitWindow(GetMonitorWidth(0), GetMonitorHeight(0), "RaylibGame");
        // No VSync: cap at the monitor's refresh rate instead (VSync on a 60 Hz desktop locked the game to 60).
        SetTargetFPS(bench ? std::max(0, benchFps) : benchFps >= 0 ? benchFps : std::max(60, GetMonitorRefreshRate(GetCurrentMonitor())));
    } else {
        InitWindow(1920, 1080, "RaylibGame bench");
        SetTargetFPS(std::max(0, benchFps)); // --fps N caps it (frame-rate-dependence checks); 0 = uncapped
    }
    if (!IsWindowReady()) {
        // raylib is built for an OpenGL 4.3 core context (compute shaders, multi-draw indirect, SSBOs): without
        // it there is no window and no fallback renderer (#185). Say so instead of exiting silently.
        const char* msg = "RaylibGame needs a graphics card and driver with OpenGL 4.3 or newer.\n\n"
                          "The OpenGL 4.3 context could not be created. Update the graphics driver, or run the game on the "
                          "dedicated GPU if this machine has two.";
        fprintf(stderr, "%s\n", msg);
#ifdef _WIN32
        MessageBoxA(nullptr, msg, "RaylibGame: OpenGL 4.3 required", 0x10u /* MB_ICONERROR */);
#endif
        return 1;
    }
    GL_Caps(); // log the GL version and the optional 4.4-4.6 paths available (#286)
    SetExitKey(KEY_NULL);
    SetTraceLogLevel(LOG_WARNING); // Set log level to warning to reduce output noise
    rlSetClipPlanes(0.5, 60000.0); // space scale: see other bodies kilometres away (0.2 m near walls: see the frame loop)
    
    DebugUI_Init();
    if (!bench) Audio_Init();
    itemDB.load("resources/data/items.txt", "resources/data/recipes.txt");
    Fluid::LoadKinds("resources/data/fluids.txt"); // water, crude oil, ... (fluid.h)
    ItemIcons_Init(); // one painted icon per item (itemicons.h)
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--test-compute") {
            // GPU compute self-test: 1M values through a compute shader, read back and checked.
            const unsigned n = 1u << 20;
            std::vector<unsigned> v(n);
            for (unsigned k = 0; k < n; ++k) v[k] = k * 3u;
            ComputeProgram prog; ComputeBuffer buf;
            bool ok = Compute_Available() && prog.load("resources/compute/selftest.comp");
            if (ok) {
                buf.create(n * 4, v.data());
                prog.use(); prog.setUint("count", n); buf.bind(0);
                const double t0 = GetTime();
                prog.dispatch(n / 64);
                std::vector<unsigned> r(n);
                buf.download(r.data(), n * 4);
                const double ms = (GetTime() - t0) * 1000.0;
                int bad = 0;
                for (unsigned k = 0; k < n; ++k) if (r[k] != k * 3u * 2u + k) ++bad;
                ok = bad == 0;
                printf("bench: compute: GL %s, 1M values dispatched and read back in %.2f ms, %d wrong\n", Compute_Available() ? "4.3" : "3.3", ms, bad);
                buf.destroy(); prog.unload();
            }
            printf("bench: compute: %s\n", ok ? "PASS" : "FAIL (no compute)");
        }
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--shot-icons") printf("bench: icons: %d items -> bench_icons.png (%s)\n", (int)itemDB.items().size(), ItemIcons_SaveSheet("bench_icons.png") ? "saved" : "FAILED");
    GameSave_Register();
    {   // Scenes own GPU resources: keep them in a scope that ends before CloseWindow().
    Scene mainMenu("mainMenu", true); // Main menu scene;
    Scene world("world", false);

    bool windowKeepAlive = true;

    mainMenuHandler(mainMenu);
    worldHandler(world);
    Rocket_Init(&world);
    Station_Init(&world);
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
    if (!benchAtSave) StartNewWorld(world, false);
    if (!bench) Settings_Load();
    bool worldFresh = !loadSaveAtStart; // the world generated at launch: New Game uses it as it is
    MenuResult pendingMenu;
    bool justPaused = false;
    bool loading = false;          // the loading screen is up
    double loadStart = 0.0;
    size_t loadMax = 1;
    bool spawnLoad = false;        // raised by a spawn or teleport (not entering the game): no minimum time
    int loadFrames = 0;
    int seenSpawns = 0;
    world.onDraw3D = [&world] { Farm_Draw3D(world); Animals_Draw3D(world); Machines_Draw3D(world); InstanceBatch_Flush(world); Weather_Draw3D(world); Particles_Draw(world); if (!ship.piloted && !buildTool.active && !Rocket_Boarded()) TerrainTool_Draw3D(); BuildTool_Draw3D(world); };
    terrainTool.till = [&world](Planetoid& body, Vector3 ground) { return Farm_Till(world, body, ground); }; // a shovel's RMB
    terrainTool.swapHeld = [](const std::string& from, const std::string& to, const char* ev) { // a bucket filled or emptied
        const int s = Hotbar_Selected();
        if (Hotbar_Slot(s) == from) Hotbar_Bind(s, to);
        QuestBook_Count(ev);
    };

    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--test-dig-shape") {
            // Digging is what the preview shows: 40 digs (straight down, the tools' radii 0.35-0.85 m, sunk
            // half a radius as the tool does) at random spots; the ground each one removes is compared with
            // the previewed sphere (both within the ground before), sampled every 5 cm: intersection over union.
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
            if (pf) {
                const Vector3 home = pf->findSpawnDir(p->seed);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                std::mt19937 rng(11);
                std::uniform_real_distribution<float> u(-1.0f, 1.0f);
                const float radii[4] = { 0.35f, 0.45f, 0.6f, 0.85f };
                double sum = 0.0, worst = 1.0; int n = 0;
                for (int k = 0; k < 40; ++k) {
                    const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, 150.0f * u(rng) / pf->radius), Vector3Scale(t2, 150.0f * u(rng) / pf->radius))));
                    // The ground along the radial line (the real surface, 3D band and all).
                    float h = pf->radius + pf->elevation(d) + 20.0f;
                    while (h > pf->radius - 100.0f && p->densityAt(Vector3Scale(d, h)) < ISO_LEVEL) h -= 0.02f;
                    const float r = radii[k % 4];
                    const Vector3 centre = Vector3Scale(d, h - 0.5f * r);
                    const float step = 0.05f, box = r + 0.3f;
                    std::vector<char> before;
                    std::vector<Vector3> pts;
                    for (float z = -box; z <= box; z += step)
                        for (float y = -box; y <= box; y += step)
                            for (float x = -box; x <= box; x += step) {
                                const Vector3 q = Vector3Add(centre, { x, y, z });
                                pts.push_back(q);
                                before.push_back(p->densityAt(q) >= ISO_LEVEL);
                            }
                    TerrainEdit e{};
                    e.center = centre; e.radius = r; e.carve = true;
                    ApplyTerrainEdit(*p, e);
                    long inter = 0, uni = 0;
                    for (size_t j = 0; j < pts.size(); ++j) {
                        if (!before[j]) continue;
                        const bool want = Vector3Distance(pts[j], centre) < r;
                        const bool got = p->densityAt(pts[j]) < ISO_LEVEL;
                        inter += want && got; uni += want || got;
                    }
                    const double iou = uni ? (double)inter / uni : 1.0;
                    sum += iou; worst = std::min(worst, iou); ++n;
                }
                printf("bench: dig shape: %d digs (radius 0.35-0.85 m, lattice %.2f m): mean IoU with the previewed sphere %.3f, worst %.3f -> %s\n",
                       n, kVoxelSize, sum / std::max(1, n), worst, sum / std::max(1, n) >= 0.85 && worst >= 0.7 ? "PASS" : "FAIL");
            }
            return 0;
        }
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--test-terrain-gpu") {
            // Plan 6 step 3 (#146): the GPU port of the height map (resources/compute/terrain.glsl) against
            // PlanetField::elevationBase over 1M points: random directions (half near the start, half anywhere)
            // and footprints 0-1024 m. Gate: every height within 1e-4 m.
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
            ComputeProgram prog;
            if (!pf || !prog.load("resources/compute/terrain_test.comp")) { printf("bench: terrain gpu: no planet or shader -> FAIL\n"); return 1; }
            const int N = 1 << 20;
            std::vector<float> in(N * 4);
            std::mt19937 rng(11);
            std::uniform_real_distribution<float> u(-1.0f, 1.0f);
            const Vector3 home = pf->findSpawnDir(p->seed);
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
            const float fps[] = { 0.0f, 0.0f, 0.25f, 1.0f, 2.0f, 4.0f, 8.0f, 16.0f, 32.0f, 64.0f, 128.0f, 256.0f, 512.0f, 1024.0f };
            for (int k = 0; k < N; ++k) {
                Vector3 d;
                if (k & 1) d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, 5000.0f * u(rng) / pf->radius), Vector3Scale(t2, 5000.0f * u(rng) / pf->radius))));
                else { do d = { u(rng), u(rng), u(rng) }; while (Vector3LengthSqr(d) < 0.01f || Vector3LengthSqr(d) > 1.0f); d = Vector3Normalize(d); }
                in[k * 4] = d.x; in[k * 4 + 1] = d.y; in[k * 4 + 2] = d.z;
                in[k * 4 + 3] = fps[(rng() >> 4) % (sizeof(fps) / sizeof(fps[0]))];
            }
            std::vector<float> params;
            pf->gpuNoiseParams(params);
            ComputeBuffer bin, bout;
            bin.create(N * 16, in.data(), false);
            bout.create(N * 4, nullptr, true);
            prog.use();
            rlSetUniform(prog.uniform("noiseParams"), params.data(), RL_SHADER_UNIFORM_VEC4, (int)params.size() / 4);
            prog.setFloat("planetRadius", pf->radius);
            prog.setInt("planetMoon", 0);
            prog.setUint("count", (unsigned)N);
            bin.bind(0); bout.bind(1);
            const double g0 = GetTime();
            prog.dispatch((N + 63) / 64);
            std::vector<float> gpu(N);
            bout.download(gpu.data(), N * 4);
            const double gpuMs = (GetTime() - g0) * 1000.0;
            const double c0 = GetTime();
            std::vector<float> cpu(N);
            for (int k = 0; k < N; ++k) cpu[k] = pf->elevationBase({ in[k * 4], in[k * 4 + 1], in[k * 4 + 2] }, in[k * 4 + 3]);
            const double cpuMs = (GetTime() - c0) * 1000.0;
            double worst = 0.0; int over = 0, exact = 0, worstK = 0;
            for (int k = 0; k < N; ++k) {
                const double e = fabs((double)gpu[k] - cpu[k]);
                if (e > worst) { worst = e; worstK = k; }
                over += e > 1e-4; exact += e == 0.0;
            }
            printf("bench: terrain gpu: %d points, %d exact, %d over 1e-4 m, worst %.3g m (cpu %.4f gpu %.4f at fp %.2f); cpu %.0f ms (1 thread), gpu %.1f ms (incl. readback) -> %s\n",
                   N, exact, over, worst, cpu[worstK], gpu[worstK], in[worstK * 4 + 3], cpuMs, gpuMs, over == 0 ? "PASS" : "FAIL");
            bin.destroy(); bout.destroy(); prog.unload();
            // Step 4 (#147): whole density blocks (as chunk builds sample them: CHUNK_SIZE + 4 per axis), CPU
            // PlanetField::sampleBlock vs the GPU, at every level of detail, around the start and anywhere.
            int blocks = 0, skipped = 0;
            long long samples = 0, bOver = 0, bExact = 0;
            double bWorst = 0.0, cpuBlockMs = 0.0, gpuBlockMs = 0.0;
            const int E = CHUNK_SIZE + 4;
            std::vector<float> cpuB((size_t)E * E * E), gpuB((size_t)E * E * E), isl;
            for (int b = 0; b < 240; ++b) {
                const float stride = 0.25f * (float)(1 << (b % 9));  // 0.25 .. 64 m
                Vector3 d;
                if (b % 3) d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, 3000.0f * u(rng) / pf->radius), Vector3Scale(t2, 3000.0f * u(rng) / pf->radius))));
                else { do d = { u(rng), u(rng), u(rng) }; while (Vector3LengthSqr(d) < 0.01f || Vector3LengthSqr(d) > 1.0f); d = Vector3Normalize(d); }
                const float r = pf->radius + pf->elevation(d) + stride * E * (u(rng) * 0.6f - 0.2f) - 30.0f * (b % 5 == 0);
                const Vector3 o = Vector3AddValue(Vector3Scale(d, r), -0.5f * stride * (E - 1));
                const float fp = std::max(stride, 1.0f);
                PlanetField::BlockGrid g;
                if (!pf->prepareBlock(o, stride, E, fp, g) || !pf->blockIslands(o, stride, E, isl, 64)) { ++skipped; continue; }
                const double c1 = GetTime();
                pf->sampleBlock(o, stride, E, fp, g, cpuB.data());
                cpuBlockMs += (GetTime() - c1) * 1000.0;
                const double g1 = GetTime();
                if (!TerrainGpu_SampleBlock(*pf, o, stride, E, fp, g, isl, gpuB.data())) { printf("bench: terrain gpu: block dispatch failed -> FAIL\n"); return 1; }
                gpuBlockMs += (GetTime() - g1) * 1000.0;
                ++blocks;
                for (size_t k = 0; k < cpuB.size(); ++k) {
                    const double e = fabs((double)gpuB[k] - cpuB[k]);
                    bWorst = std::max(bWorst, e); bOver += e > 1e-4; bExact += e == 0.0; ++samples;
                }
            }
            printf("bench: terrain gpu: %d blocks (%d skipped: off the band), %lld samples, %lld exact, %lld over 1e-4, worst %.3g; per block cpu %.2f ms, gpu %.2f ms (sync, incl. upload + readback) -> %s\n",
                   blocks, skipped, samples, bExact, bOver, bWorst, cpuBlockMs / std::max(1, blocks), gpuBlockMs / std::max(1, blocks), bOver == 0 ? "PASS" : "FAIL");
            return over == 0 && bOver == 0 ? 0 : 1;
        }
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--profile-density") {
            // Plan item 7: what terrain sampling costs, part by part, over 200k points within 20 m of the
            // surface around the start (one thread, ns per call).
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
            if (pf) {
                const Vector3 home = pf->findSpawnDir(p->seed);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                std::vector<Vector3> pts; std::vector<float> surf;
                std::mt19937 rng(7);
                std::uniform_real_distribution<float> u(-1.0f, 1.0f);
                for (int k = 0; k < 200000; ++k) {
                    const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, 3000.0f * u(rng) / pf->radius), Vector3Scale(t2, 3000.0f * u(rng) / pf->radius))));
                    const float s = pf->radius + pf->elevation(d);
                    surf.push_back(s); pts.push_back(Vector3Scale(d, s + 20.0f * u(rng)));
                }
                auto timeIt = [&](const char* name, auto&& f) {
                    volatile float sink = 0.0f;
                    const auto t0 = std::chrono::steady_clock::now();
                    for (size_t k = 0; k < pts.size(); ++k) sink = sink + f(k);
                    const double ns = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / pts.size();
                    printf("bench: density profile: %-26s %7.0f ns\n", name, ns);
                };
                {
                    // noise4 must match noise() exactly (the terrain is sampled both ways nowhere, but
                    // determinism across builds matters).
                    std::uniform_real_distribution<float> big(-5000.0f, 5000.0f);
                    double worst = 0.0; int diff = 0;
                    for (int k = 0; k < 250000; ++k) {
                        float xs[4], ys[4], zs[4], o[4];
                        for (int l = 0; l < 4; ++l) { xs[l] = big(rng); ys[l] = big(rng); zs[l] = big(rng); }
                        SimplexNoise::noise4(xs, ys, zs, o);
                        for (int l = 0; l < 4; ++l) { const double e = fabs((double)o[l] - SimplexNoise::noise(xs[l], ys[l], zs[l])); worst = std::max(worst, e); diff += e > 0.0; }
                    }
                    printf("bench: density profile: noise4 vs noise over 1M points: %d differ, worst %.3g -> %s\n", diff, worst, diff == 0 ? "PASS" : (worst < 1e-5 ? "PASS (rounding)" : "FAIL"));
                }
                timeIt("elevation (height map)", [&](size_t k) { return pf->elevation(Vector3Normalize(pts[k])); });
                timeIt("terrain3D (3D band)", [&](size_t k) { const float d = Vector3Length(pts[k]); return pf->terrain3D(pts[k], d, surf[k], 1.0f); });
                timeIt("rockFeature (tors, boulders)", [&](size_t k) { const float d = Vector3Length(pts[k]); return pf->rockFeature(pts[k], d, surf[k], 1.0f); });
                timeIt("density (all, no edits)", [&](size_t k) { return pf->density(pts[k], nullptr, 1.0f); });
                timeIt("waterSourceAt", [&](size_t k) { return (float)pf->waterSourceAt(Vector3Normalize(pts[k])).kind; });
            }
            return 0;
        }

    const double benchWarmup = 2.0;
    double benchStart = 0.0, benchMaxFrame = 0.0;
    std::vector<float> benchFrameTimes; // every timed frame (s): percentiles at the end (#258)
    long benchFrames = 0;
    if (bench) {
        gameState = 1;
        mainMenu.isActive = false;
        world.isActive = true;
        Planetoid* benchPlanet = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
        if (!benchAtSave && !benchCamGiven && benchPlanet && dynamic_cast<const PlanetField*>(benchPlanet->field.get())) {
            // A planet (#291): 40 m above the start, looking ahead along the surface and a little down (the asteroid
            // offset below put the camera under the terrain). It rides the planet from here (#127).
            const Vector3 feet = (player.position - world.renderOrigin).toFloat();
            const Vector3 ahead = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
            world.camera.position = Vector3Add(feet, Vector3Scale(player.up, 40.0f));
            world.camera.target = Vector3Add(world.camera.position, Vector3Subtract(Vector3Scale(ahead, 100.0f), Vector3Scale(player.up, 25.0f)));
            world.camera.up = player.up;
        } else if (!benchAtSave) world.camera.position = Vector3Add(world.camera.target, benchCamOffset);
        if (benchDay) { // mid-morning over the start: 45 degrees up, from ahead-left (shadows fall toward the camera's right)
            const Vector3 ahead = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
            const Vector3 side = Vector3CrossProduct(ahead, player.up);
            sunOverride = true;
            sunArg = Vector3Normalize(Vector3Add(Vector3Scale(player.up, 0.7f), Vector3Add(Vector3Scale(ahead, 0.5f), Vector3Scale(side, -0.5f))));
        }
        world.debugDraw = startWithOverlay;
        player.flyMode = benchWalkSeconds <= 0.0 && !benchWalkMode; // the bench drives the camera directly unless walking
        player.autoWalk = benchWalkMode && !benchStandMode;           // --bench-walk: walk through the measured window (--bench-stand: stand)
        player.flyFollow = benchFlySpeed <= 0.0f && benchFlyAlt >= 1e8f; // #298: --bench-fly places the camera itself every frame (the carry undid it). #127: the bench camera rides the body under it (sphere of influence), as the game's fly camera does; it
                                                                      // stood still in absolute space while the planet spun past (benches that place the camera every frame override it)
        benchStart = -1.0; // starts once the world has finished streaming in
        if (benchPreloadSpawn) { BeginSpawnPreload(); preloadT0 = GetTime(); }
        if (benchWaterStats) {
            // --water-stats (#251 investigation): how far the water shaping (rivers, lakes, ponds) raises or cuts the
            // ground: elevation() vs elevationBase() (the land before any water) on a 12 km grid around the start.
            Planetoid* wp = dynamic_cast<Planetoid*>(world.objects[0].get());
            const PlanetField* wf = wp ? dynamic_cast<const PlanetField*>(wp->field.get()) : nullptr;
            if (wf) {
                const Vector3 home = wf->findSpawnDir(wp->seed);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                const char* zone[4] = { "river zone (within half + 12 m)", "lake", "inside a pond", "elsewhere (pond embankments, river valleys)" };
                long n[4] = {}, fill1[4] = {}, fill3[4] = {}, fill6[4] = {}, cut3[4] = {};
                float maxFill[4] = {};
                const float step = 20.0f;
                for (float e = -6000.0f; e <= 6000.0f; e += step)
                    for (float nn = -6000.0f; nn <= 6000.0f; nn += step) {
                        const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / wf->radius), Vector3Scale(t2, nn / wf->radius))));
                        const float diff = wf->elevation(d, 1.0f) - wf->elevationBase(d, 1.0f);
                        int z = 3;
                        float lv, lo;
                        const RiverNetwork* net = wf->rivers();
                        if (net && net->lakeAt(d, lv, lo) && lo < RiverNetwork::kLakeMargin * 2.0f) z = 1;
                        else if (net) { const RiverNetwork::Hit h = net->nearest(d, 60.0f); if (h.found && h.dist < h.halfWidth + 12.0f) z = 0; }
                        if (z == 3) {
                            const PlanetField::WaterSource ws = wf->waterSourceAt(d);
                            if (ws.kind == 3) z = 2;
                        }
                        ++n[z];
                        if (diff > 1.0f) ++fill1[z];
                        if (diff > 3.0f) ++fill3[z];
                        if (diff > 6.0f) ++fill6[z];
                        if (diff < -3.0f) ++cut3[z];
                        maxFill[z] = std::max(maxFill[z], diff);
                    }
                const float cell = step * step / 1e4f; // hectares per sample
                for (int z = 0; z < 4; ++z)
                    printf("bench: water stats: %s: %.0f ha; raised > 1 m %.1f ha, > 3 m %.1f ha, > 6 m %.1f ha (max %.1f m); cut > 3 m %.1f ha\n",
                           zone[z], n[z] * cell, fill1[z] * cell, fill3[z] * cell, fill6[z] * cell, maxFill[z], cut3[z] * cell);
            }
            benchWaterStats = false;
        }
    }

    while (windowKeepAlive && !WindowShouldClose())
    {
        GpuTimer_CpuFrameStart(); GpuTimer_BeginFrame();
        if (preloadT0 >= 0.0 && !SpawnPreloadPending()) {
            const WorldStats& ws = GetWorldStats();
            printf("bench: spawn preload: %d ring nodes built in %.2f s; arena %.1f M vertices, %d nodes alive\n",
                   ws.preloadNodes, GetTime() - preloadT0, TerrainArena_Stats().vertices / 1048576.0, ws.nodes);
            preloadT0 = -1.0;
        }
        if (bench && !loading) { // scripted tests wait out the loading screen (a half-built world is no game state)
            if (benchStart < 0.0 && preloadT0 >= 0.0) {
                // --preload-spawn: the timed window starts after the spawn preload, as play does
            } else if (benchStart < 0.0) {
                if (WorldIdle() && benchDigStrokes > 0) {
                    // Scripted edit test: carve a shaft straight down under the player.
                    TerrainHit hit;
                    const Vector3 eye = Vector3Add((player.position - world.renderOrigin).toFloat(), Vector3Scale(player.up, 6.0f)), dir = Vector3Negate(player.up); // straight down: the bench camera is not placed yet
                    if (RaycastTerrain(world, Ray{ eye, dir }, 300.0f, hit)) {
                        for (int s = 0; s < benchDigStrokes; ++s) {
                            const Vector3 local = hit.planetoid->toLocal(Vector3Add(hit.point, Vector3Scale(dir, 2.0f * s)));
                            ApplyTerrainEdit(*hit.planetoid, TerrainEdit{ local, 8.0f, -1.0f });
                        }
                        printf("bench: dug %d strokes at %.1f %.1f %.1f (%s)\n", benchDigStrokes, hit.point.x, hit.point.y, hit.point.z, VoxelMaterialName(hit.material));
                    } else {
                        printf("bench: --test-dig found no terrain at the crosshair\n");
                    }
                    benchDigStrokes = 0;
                } else if (WorldIdle() && benchPaint) {
                    // Placed ground keeps its material: build a gravel mound at the crosshair, look at it,
                    // then save, load and look again.
                    TerrainHit hit;
                    const Vector3 eye = Vector3Add((player.position - world.renderOrigin).toFloat(), Vector3Scale(player.up, 6.0f)), dir = Vector3Negate(player.up); // straight down from above the player
                    if (RaycastTerrain(world, Ray{ eye, dir }, 300.0f, hit)) {
                        const VoxelMaterial before = hit.material;
                        Planetoid& p = *hit.planetoid;
                        const Vector3 up = Vector3Normalize(hit.local);
                        for (int s = 0; s < 4; ++s) {
                            TerrainEdit e{};
                            e.center = Vector3Add(hit.local, Vector3Scale(up, 0.5f + 0.4f * s));
                            e.radius = 2.0f;
                            e.strength = 0.6f;
                            e.paint = (int)VoxelMaterial::Gravel;
                            ApplyTerrainEdit(p, e);
                        }
                        TerrainHit a, b;
                        const bool hitA = RaycastTerrain(world, Ray{ eye, dir }, 300.0f, a);
                        const bool saved = SaveWorld(world) && LoadWorldSave(world);
                        const Vector3 eye2 = Vector3Add((player.position - world.renderOrigin).toFloat(), Vector3Scale(player.up, 6.0f)), dir2 = Vector3Negate(player.up);
                        const bool hitB = RaycastTerrain(world, Ray{ eye2, dir2 }, 300.0f, b);
                        const bool ok = hitA && a.material == VoxelMaterial::Gravel && saved && hitB && b.material == VoxelMaterial::Gravel;
                        printf("bench: paint: ground was %s; placed gravel reads %s, after save/load %s -> %s\n", VoxelMaterialName(before),
                               hitA ? VoxelMaterialName(a.material) : "nothing", hitB ? VoxelMaterialName(b.material) : "nothing", ok ? "PASS" : "FAIL");
                    } else {
                        printf("bench: --test-paint found no terrain at the crosshair -> FAIL\n");
                    }
                    benchPaint = false;
                } else if (WorldIdle() && benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && benchFlySeconds <= 0.0) {
                    const int editsBefore = GetWorldStats().edits;
                    const int itemsBefore = inventory.count("stone") + inventory.count("regolith");
                    const bool ok = SaveWorld(world) && LoadWorldSave(world); // restores camera + player from the file
                    printf("bench: save/load %s (%d edits, %d stone+regolith before / %d after)\n", ok ? "ok" : "FAILED",
                           editsBefore, itemsBefore, inventory.count("stone") + inventory.count("regolith"));
                    benchSaveLoad = false;
                } else if (benchSuffocate && benchWalkSeconds <= 0.0 && benchFlySeconds <= 0.0 && benchMineStrokes == 0) {
                    // Scripted survival: out of oxygen -> die -> respawn at the start, at rest.
                    if (suffocateStart < 0.0) {
                        player.flyMode = false;
                        player.velocity = Vector3Scale(player.up, 30.0f); // moving fast when it dies: none of it may survive the respawn
                        survival.oxygen = 0.0f;
                        survival.health = 15.0f;
                        survival.breathableAltitude = -1e9f; // no air to breathe (the planet's would refill the suit)
                        suffocateStart = GetTime();
                    } else if (survival.dead && player.nearest) {
                        // #231: while dead the camera holds its place in the body's frame (drift over the death).
                        static DVec3 firstLocal{};
                        static bool have = false;
                        static double deadAt = -1.0, settleDrift = 0.0;
                        if (deadAt < 0.0) deadAt = GetTime();
                        const DVec3 camLocal = player.nearest->toLocalD(world.renderOrigin + world.camera.position);
                        if (!have) firstLocal = camLocal; // the first dead frame: how far the eye smoothing then settles
                        if (GetTime() - deadAt < 0.2) { settleDrift = std::max(settleDrift, (camLocal - firstLocal).length()); firstLocal = camLocal; have = true; }
                        else suffocateDrift = std::max(suffocateDrift, (camLocal - firstLocal).length()); // held after that
                        suffocateSettle = settleDrift;
                    } else if (survival.deaths > 0 && !survival.dead) {
                        printf("bench: suffocate: while dead the camera drifted %.3f m in the body's frame (#231) -> %s\n", suffocateDrift, suffocateDrift < 0.01 ? "PASS" : "FAIL");
                        const float d = Vector3Length(Vector3Subtract(player.velocity, !player.grounded && player.nearest ? player.nearest->surfaceVelocity((player.position - world.renderOrigin).toFloat()) : Vector3{ 0, 0, 0 }));
                        printf("bench: suffocated and respawned after %.1f s; O2 %.0f HP %.0f, moving %.1f m/s against the ground, deaths %d\n",
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
                            if (!chopId) { chopId = id; chopWood0 = inventory.count("any_log"); }
                            ++chopSwings;
                            if (TerrainTool_Chop()) {
                                const int wood = inventory.count("any_log") - chopWood0;
                                printf("bench: chop: felled a %.1f m tree in %d swings, +%d logs; gone: %s\n", chopTarget.height, chopSwings, wood, present(chopId) ? "no" : "yes");
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
                        size_t deltaBytes = 0;
                        for (const auto& [k, c] : p->edits->chunks) deltaBytes += c->bytes();
                        const double deltaMB = deltaBytes / 1048576.0, denseMB = chunks * (DeltaChunk::kN * 2.0) / 1048576.0;
                        const float probe = p->densityAt(strokes[0].center);
                        t = clk::now();
                        const bool loaded = LoadWorldSave(world);
                        const double loadMs = ms(t);
                        p = dynamic_cast<Planetoid*>(world.objects[0].get());
                        const bool same = loaded && p->edits->chunks.size() == chunks && fabsf(p->densityAt(strokes[0].center) - probe) < 1e-4f;
                        const bool pass = same && saveFull < 1000.0 && loadMs < 1000.0 && buildAfter < buildBefore * 1.5 + 0.5;
                        printf("bench: persist: %d strokes folded in %.0f ms into %zu delta chunks (%zu strokes recorded)\n", benchPersist, foldMs, chunks, p->edits->strokes);
                        printf("bench: persist: delta memory %.0f MB (%.0f MB if dense)\n", deltaMB, denseMB);
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
                            Factory::Seed(p->id, factPlot, s[0], "oak_log", 32);
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
                            // Remove the plot's first part before saving: the chest's index must survive the reload
                            // (saves kept only live parts, so every index after a removed one shifted down).
                            if (Plot* pl = p->plots ? p->plots->byId(factPlot) : nullptr) pl->remove(0);
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
                        if (Progression_Update(world, progressionFull, progressionEv) != 0) benchProgression = false;
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
                        if (p->densityAt(p->toLocal(world.camera.position)) >= ISO_LEVEL) ++wallCamInside;
                        {   // Shake (#255): the change of the camera's per-frame step on the body (smooth motion: ~0; a jumping
                            // keep-out correction or push-out: centimetres).
                            const DVec3 cam = p->toLocalD(world.renderOrigin + world.camera.position);
                            const DVec3 stepNow = cam - wallCamPrev;
                            if (wallCamKnown >= 2) {
                                const double s = (stepNow - wallCamStepPrev).length();
                                if (s > 0.005) { ++wallCamShake; if (player.lastEyePush > 0.0f) ++wallCamShakePushed; if (player.grounded) ++wallCamShakeGrounded; }
                                if (s > 0.02 && getenv("WALL_DEBUG")) printf("SHAKE t=%.3f phase %d s=%.1f mm push=%.3f grounded=%d feetstep=%.1f mm\n", t, wallPhase, s * 1000.0, player.lastEyePush, player.grounded ? 1 : 0, (local - wallLast).length() * 1000.0);
                                wallCamShakeMax = std::max(wallCamShakeMax, s);
                            }
                            wallCamPrev = cam; wallCamStepPrev = stepNow; ++wallCamKnown;
                        }
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
                        const bool pass = survival.deaths == 0 && wallMinHealth > 90.0f && wallMaxStep < 1.0;
                        printf("bench: wall shove 9 s: largest frame move %.2f m (at %.1f s), top speed vs ground %.1f m/s, worst impact %.1f m/s, health %.0f, deaths %d: %s\n",
                               wallMaxStep, wallMaxStepAt, wallMaxSpeed, wallMaxImpact, wallMinHealth, survival.deaths, pass ? "PASS" : "FAIL");
                        wallOk = pass;
                        // Tight corners: two crevices narrowing to under the body's width meet in a V ahead;
                        // walk and jump into its point, turning a little either way.
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        const Vector3 side = Vector3CrossProduct(flat, player.up);
                        for (int arm = -1; arm <= 1; arm += 2)
                            for (int k = 0; k < 8; ++k) {
                                const float r = std::max(0.42f, 1.2f - 0.11f * k);
                                const Vector3 along = Vector3Add(Vector3Scale(flat, 0.8f * k), Vector3Scale(side, arm * (1.6f - 0.2f * k)));
                                ApplyTerrainEdit(*p, TerrainEdit{ p->toLocal(Vector3Add(feet, Vector3Add(along, Vector3Scale(player.up, r)))), r, -10.0f });
                            }
                        player.forward = flat;
                        wallMaxStep = 0.0; wallMaxImpact = 0.0f; wallMaxSpeed = 0.0f; wallMinHealth = survival.health;
                        wallPhase = 3;
                        wallStart = GetTime();
                    } else if (wallPhase == 3 && t < 7.0) {
                        const int k = (int)(t * 2.0);
                        static int lastK = -1;
                        if (k != lastK && t > 1.0) {
                            lastK = k;
                            const Matrix r = MatrixRotate(player.up, (k % 2 ? 20.0f : -20.0f) * DEG2RAD);
                            player.forward = Vector3Normalize(Vector3Transform(player.forward, r));
                            if (k % 3 == 0) player.autoJump = true;
                        }
                    } else if (wallPhase == 3) {
                        player.autoWalk = false;
                        const bool pass = survival.deaths == 0 && wallMinHealth > 90.0f && wallMaxStep < 1.0;
                        printf("bench: tight corner 6 s: largest frame move %.2f m (at %.1f s), top speed vs ground %.1f m/s, worst impact %.1f m/s, health %.0f, deaths %d: %s\n",
                               wallMaxStep, wallMaxStepAt, wallMaxSpeed, wallMaxImpact, wallMinHealth, survival.deaths, pass ? "PASS" : "FAIL");
                        printf("bench: wall: camera inside terrain in %d frames (#255)\n", wallCamInside);
                        printf("bench: wall: camera shake: %d frames over 5 mm (change of its per-frame step; %d with the eye pushed out of the terrain, %d grounded), worst %.1f mm (#255)\n", wallCamShake, wallCamShakePushed, wallCamShakeGrounded, wallCamShakeMax * 1000.0);
                        printf("bench: wall: %s\n", wallOk && pass && wallCamInside == 0 ? "PASS" : "FAIL");
                        benchWall = false;
                    }
                } else if (benchTeleClose && (teleCPhase > 0 || WorldIdle())) {
                    // --test-telescope-close (#261): look at the horizon through the telescope at 16x until the zoomed view
                    // has loaded, close it, and report the worst frame and the builds asked for over the next 4 s.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - teleCStart;
                    if (teleCPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        inventory.add("telescope", 1);
                        teleCPhase = 1; teleCStart = GetTime();
                    } else if (teleCPhase == 1 && t > 1.5 && player.grounded) {
                        player.forward = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))); // level: the horizon
                        Telescope_SetActive(world, true);
                        telescope.zoom = 16.0f;
                        teleCPhase = 2; teleCStart = GetTime();
                    } else if (teleCPhase == 2) {
                        telescope.zoom = 16.0f;
                        if ((t > 3.0 && WorldIdle()) || t > 20.0) {
                            printf("bench: telescope close: zoomed view loaded in %.1f s (%d nodes)\n", t, GetWorldStats().nodes);
                            teleCReq0 = GetWorldStats().requests;
                            Telescope_SetActive(world, false);
                            teleCWorst = 0.0;
                            teleCPhase = 3; teleCStart = GetTime();
                        }
                    } else if (teleCPhase == 3) {
                        teleCWorst = std::max(teleCWorst, (double)GetFrameTime());
                        if (t > 4.0) {
                            const long long req = GetWorldStats().requests - teleCReq0;
                            printf("bench: telescope close: worst frame in the 4 s after closing %.1f ms, %lld builds asked for, %d nodes -> %s\n",
                                   teleCWorst * 1000.0, req, GetWorldStats().nodes, teleCWorst < 0.025 ? "PASS" : "FAIL");
                            benchTeleClose = false;
                        }
                    }
                } else if (benchFall && (fallPhase > 0 || WorldIdle())) {
                    // --test-fall (#263): drops onto flat ground near the start. Expected on Terra (g 5.9, safe 6 m/s, x8.5):
                    // 3 m none, 10 m ~41, 20 m ~79, 30 m and 100 m fatal.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const float heights[] = { 3.0f, 10.0f, 20.0f, 30.0f, 100.0f };
                    const float lo[] = { 0.0f, 25.0f, 60.0f, 100.0f, 100.0f }, hi[] = { 0.0f, 55.0f, 95.0f, 1e9f, 1e9f };
                    if (fallPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        fallPhase = 1; fallStart = GetTime();
                    } else if (fallPhase == 1 && GetTime() - fallStart > 1.5 && player.grounded) {
                        survival.health = 100.0f; survival.dead = false;
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        player.position = player.position + DVec3(Vector3Scale(player.up, heights[fallIdx]));
                        fallPhase = 2; fallStart = GetTime();
                    } else if (fallPhase == 2 && (player.grounded || survival.dead) && GetTime() - fallStart > 0.3) {
                        const float lost = 100.0f - std::max(survival.health, 0.0f);
                        const bool dead = survival.dead || survival.health <= 0.0f;
                        const bool pass = lost >= lo[fallIdx] - 0.01f && lost <= hi[fallIdx] + 0.01f;
                        printf("bench: fall: %.0f m: lost %.0f health%s (expected %.0f-%s) -> %s\n", heights[fallIdx], lost, dead ? " (dead)" : "", lo[fallIdx],
                               hi[fallIdx] > 1e8f ? "fatal" : TextFormat("%.0f", hi[fallIdx]), pass ? "ok" : "FAIL");
                        fallOk = fallOk && pass;
                        if (dead) Survival_Respawn(world);
                        if (++fallIdx < 5) { fallPhase = 1; fallStart = GetTime() - 1.0; }
                        else { printf("bench: fall: %s\n", fallOk ? "PASS" : "FAIL"); survival.health = 100.0f; benchFall = false; }
                    } else if (fallPhase == 2 && GetTime() - fallStart > 30.0) {
                        printf("bench: fall: %.0f m: never landed -> FAIL\n", heights[fallIdx]); benchFall = false;
                    }
                } else if (benchMachines && (machPhase > 0 || WorldIdle())) {
                    // --test-machines: the fuel slot (#194), inputs consumed when a craft ends and a craft cancelled when
                    // its input is taken (#193), the boiler's water, heat and banked fire without a cable (#229).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (machPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        machPhase = 1; machStart = GetTime();
                    } else if (machPhase == 1 && GetTime() - machStart > 1.5 && player.grounded) {
                        bool ok = true;
                        auto check = [&](bool c, const char* what) { printf("bench: machines: %s: %s\n", what, c ? "ok" : "FAILED"); ok = ok && c; };
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(p->toLocal(Vector3Add(Vector3Add(feet, Vector3Scale(flat, 12.0f)), Vector3Scale(player.up, 2.0f))), 0.0f, true);
                        const int fire = plot->place({ PartType::Campfire, { 0, 0, 0 }, 0 });
                        const int kiln = plot->place({ PartType::Kiln, { 4, 0, 0 }, 0 });
                        const int boiler = plot->place({ PartType::Boiler, { 8, 0, 0 }, 0 });
                        Factory::Sync(world);
                        Factory::WaitTicks(3);
                        Factory::SetSpeed(20);
                        auto view = [&](int part) { Factory::View(p->id, plot->id, part); Factory::WaitTicks(2); Factory::Sync(world); return Factory::View(p->id, plot->id, part); };
                        auto inputCount = [](const MachineView& v, const char* item) { int n = 0; for (const auto& s : v.inSlots) if (s.first == item) n += s.second; return n; };
                        auto outputCount = [](const MachineView& v, const char* item) { int n = 0; for (const auto& s : v.outSlots) if (s.first == item) n += s.second; return n; };
                        // #194: logs in a campfire (it also chars them): the first load lights the empty fuel slot, the next goes to the inputs.
                        Factory::Seed(p->id, plot->id, fire, "oak_log", 5);
                        MachineView v = view(fire);
                        check(v.hasFuelSlot && v.fuelSlot.first == "oak_log" && v.fuelSlot.second >= 4 && inputCount(v, "oak_log") == 0, "logs light an empty fuel slot");
                        Factory::Seed(p->id, plot->id, fire, "oak_log", 3);
                        v = view(fire);
                        check(inputCount(v, "oak_log") > 0 || outputCount(v, "charcoal") > 0, "more logs go to the inputs (to char), not the fuel");
                        // #193: kiln, ore in, coal as fuel; mid-craft the ore is still in the input slot; taking it cancels the craft.
                        Factory::Seed(p->id, plot->id, kiln, "coal", 4);
                        Factory::Seed(p->id, plot->id, kiln, "iron_ore", 1);
                        Factory::WaitTicks(60);
                        v = view(kiln);
                        check(v.fuelSlot.first == "coal" && v.duration > 0 && v.progress > 0 && inputCount(v, "iron_ore") == 1, "mid-craft the ore is still in the input slot");
                        Factory::TakeSlot(p->id, plot->id, kiln, 0, 0);
                        Factory::WaitTicks(400);
                        v = view(kiln);
                        check(v.duration == 0 && outputCount(v, "iron_ingot") == 0, "taking the ore out cancels the craft");
                        Factory::Seed(p->id, plot->id, kiln, "iron_ore", 1);
                        Factory::WaitTicks(400);
                        v = view(kiln);
                        check(outputCount(v, "iron_ingot") == 1 && inputCount(v, "iron_ore") == 0, "a finished craft takes its ore and makes an ingot");
                        // #229: a boiler without a cable: fuel alone does nothing (no water); a bucket of water lights it, it heats
                        // up, then banks its fire (no load) and the empty bucket comes back.
                        Factory::Seed(p->id, plot->id, boiler, "coal", 8);
                        Factory::WaitTicks(40);
                        v = view(boiler);
                        check(v.status == MachineStatus::NoWater && v.heat == 0 && v.waterMax > 0, "fuel without water: no water, cold");
                        inventory.add("water_bucket", 1);
                        const int buckets0 = inventory.count("bucket");
                        if (inventory.remove("water_bucket", 1)) Factory::Insert(p->id, plot->id, boiler, "water_bucket", 1);
                        Factory::WaitTicks(200);
                        v = view(boiler);
                        const int heatMid = v.heat;
                        check(heatMid > 0 && heatMid < v.heatMax && v.water < 1000 && inventory.count("bucket") == buckets0 + 1, "a bucket of water: it heats up; the empty bucket comes back");
                        Factory::WaitTicks(400);
                        v = view(boiler);
                        const int water1 = v.water;
                        Factory::WaitTicks(400);
                        v = view(boiler);
                        printf("bench: machines: boiler heat %d / %d, water %d -> %d L over 400 banked ticks, status %s\n", v.heat, v.heatMax, water1, v.water, MachineStatusName(v.status));
                        check(v.heat == v.heatMax && v.status == MachineStatus::NoLoad && water1 - v.water <= 60, "hot with no load: the fire is banked (about 1 tick in 8)");
                        Factory::SetSpeed(1);
                        printf("bench: machines: %s\n", ok ? "PASS" : "FAIL");
                        benchMachines = false;
                    }
                } else if (benchTor && (torPhase > 0 || WorldIdle())) {
                    // --test-tor (#242): the lowest-lying tor near the start (beach country, 4-40 m up), the player
                    // put on its crown and walked out across it at eight headings. Walking on a crown must not pop:
                    // per-frame moves stay small, no upward speed without a jump, grounding doesn't flicker.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - torStart;
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const DVec3 local = p->toLocalD(player.position);
                    if (torPhase == 0) {
                        const Vector3 up = StartSpawnDirLocal(*p);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                        const Vector3 t2 = Vector3CrossProduct(t1, up);
                        float best = 1e9f;
                        for (int i = -240; pf && i <= 240; ++i)
                            for (int j = -240; j <= 240; ++j) {
                                const float dd = 25.0f * sqrtf((float)(i * i + j * j));
                                if (dd >= best) continue;
                                const Vector3 dir = Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, 25.0f * i / pf->radius), Vector3Scale(t2, 25.0f * j / pf->radius))));
                                const float e = pf->elevation(dir);
                                if (e < 4.0f || e > 40.0f) continue;
                                const float surf = pf->radius + e;
                                if (pf->rockFeature(Vector3Scale(dir, surf + 2.0f), surf + 2.0f, surf, 1.0f) > 1.0f) { best = dd; torDir = dir; }
                            }
                        if (best > 1e8f) { printf("bench: tor: no low tor within 6 km -> FAIL\n"); benchTor = false; }
                        else {
                            printf("bench: tor: a tor %.0f m from the start at %.1f m elevation\n", best, pf->elevation(torDir));
                            if (getenv("TOR_PROFILE")) { // debug: the top solid height across the tor, two directions
                                const float base = pf->radius + pf->elevation(torDir);
                                for (int ax = 0; ax < 2; ++ax) {
                                    const Vector3 a = ax ? t2 : t1;
                                    printf("bench: tor: profile %d:", ax);
                                    for (float x = -14.0f; x <= 14.0f; x += 0.5f) {
                                        const Vector3 d = Vector3Normalize(Vector3Add(torDir, Vector3Scale(a, x / pf->radius)));
                                        float top = -99.0f;
                                        for (float hh = 16.0f; hh > -4.0f; hh -= 0.05f)
                                            if (p->densityAt(Vector3Scale(d, base + hh)) >= ISO_LEVEL) { top = hh; break; }
                                        printf(" %.2f", top);
                                    }
                                    printf("\n");
                                }
                            }
                            player.spawnOn(world, *p, p->dirToWorld(torDir));
                            torPhase = 1; torStart = GetTime();
                        }
                    } else if (torPhase == 1 && t > 3.0 && WorldIdle() && player.grounded) {
                        // Respawn on the crown for each heading; walk out from its middle.
                        player.spawnOn(world, *p, p->dirToWorld(torDir));
                        torPhase = 2; torStart = GetTime();
                    } else if (torPhase == 2 && t > 1.0 && player.grounded) {
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        player.forward = Vector3Transform(flat, MatrixRotate(player.up, 45.0f * DEG2RAD * torHeading));
                        player.autoWalk = true;
                        torPhase = 3; torStart = GetTime(); torLast = local; torWasGrounded = true; torAirAt = -1.0;
                    } else if (torPhase == 3) {
                        const DVec3 dl = local - torLast;
                        const double step = dl.length();
                        const float upMove = (float)(dl.x * player.up.x + dl.y * player.up.y + dl.z * player.up.z);
                        const float vUp = Vector3DotProduct(Vector3Subtract(player.velocity, player.grounded ? Vector3{ 0, 0, 0 } : p->surfaceVelocity(feet)), player.up);
                        if (step > torMaxStep) torMaxStep = step;
                        if (upMove > 0.06f) { ++torPops; printf("bench: tor: heading %d pop %.3f m up at %.2f s (grounded %d, contact dot %.2f, v up %.2f)\n", torHeading, upMove, t, player.grounded ? 1 : 0, player.lastContactDot, vUp); }
                        if (!player.grounded && vUp > torMaxUp) torMaxUp = vUp;
                        if (torWasGrounded && !player.grounded) torAirAt = GetTime();
                        if (!torWasGrounded && player.grounded && torAirAt > 0.0 && GetTime() - torAirAt < 0.15) {
                            ++torFlicker;
                            printf("bench: tor: heading %d ground flicker (%.0f ms in the air) at %.2f s, contact dot %.2f\n", torHeading, (GetTime() - torAirAt) * 1000.0, t, player.lastContactDot);
                        }
                        torWasGrounded = player.grounded;
                        torLast = local;
                        // Off the crown (a fall started) or 2.5 s: next heading.
                        if (t > 2.5 || (!player.grounded && vUp < -3.0f)) {
                            player.autoWalk = false;
                            if (++torHeading < 8) { torPhase = 1; torStart = GetTime() - 3.0; }
                            else {
                                const bool pass = torPops <= 1 && torMaxStep < 0.6; // the #242 defect: repeated pops (sinking into the crown, shoved out); one-frame flickers at the crown edge and ramp hops are reported only
                                printf("bench: tor: 8 headings: %d pops, %d ground flickers, top upward speed %.2f m/s, largest frame move %.2f m -> %s\n",
                                       torPops, torFlicker, torMaxUp, torMaxStep, pass ? "PASS" : "FAIL");
                                benchTor = false;
                            }
                        }
                    }
                } else if (benchLedge && (ledgePhase > 0 || WorldIdle())) {
                    // --test-ledge (#240): walk, then sprint, off the sharp lip of a dug pit. In the first
                    // 0.2 s of the fall the downward speed must come from gravity alone (g*t): the edge
                    // must not throw the player down.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - ledgeStart;
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    if (ledgePhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        ledgeWalk0 = player.walkSpeed;
                        ledgeStart = GetTime();
                        ledgePhase = 1;
                    } else if (ledgePhase == 1 && t > 1.5 && player.grounded) {
                        if (ledgePass == 0) {
                            const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                            player.forward = flat;
                            for (int k = 0; k < 5; ++k) // a pit ~5 m wide and deep, its lip 2.5 m ahead
                                ApplyTerrainEdit(*p, TerrainEdit{ p->toLocal(Vector3Add(feet, Vector3Add(Vector3Scale(flat, 5.0f), Vector3Scale(player.up, -1.5f * k)))), 2.5f, -10.0f });
                        }
                        player.walkSpeed = ledgeWalk0 * (ledgePass == 1 ? player.runMultiplier : 1.0f);
                        player.autoWalk = true;
                        ledgeAir = false; ledgeExcess = 0.0f;
                        ledgePhase = 2; ledgeStart = GetTime();
                    } else if (ledgePhase == 2) {
                        const float vUp = Vector3DotProduct(Vector3Subtract(player.velocity, player.grounded ? Vector3{ 0, 0, 0 } : p->surfaceVelocity(feet)), player.up);
                        if (!ledgeAir && !player.grounded) { ledgeAir = true; ledgeAirStart = GetTime(); ledgeV0 = vUp; }
                        if (ledgeAir) {
                            const float ta = (float)(GetTime() - ledgeAirStart);
                            if (ta <= 0.2f) ledgeExcess = std::max(ledgeExcess, -vUp - player.gravity * ta);
                        }
                        if ((ledgeAir && GetTime() - ledgeAirStart > 0.3) || t > 6.0) {
                            player.autoWalk = false;
                            player.walkSpeed = ledgeWalk0;
                            const bool pass = ledgeAir && ledgeExcess < 0.5f;
                            printf("bench: ledge %s: left the ground at %.2f m/s vertical, worst extra downward speed in 0.2 s %.2f m/s: %s\n",
                                   ledgePass == 0 ? "walk" : "sprint", ledgeV0, ledgeExcess, pass ? "PASS" : "FAIL");
                            ledgeOk = ledgeOk && pass;
                            if (++ledgePass < 2) { ledgePhase = 0; }
                            else { printf("bench: ledge: %s\n", ledgeOk ? "PASS" : "FAIL"); benchLedge = false; }
                        }
                    }
                } else if (benchRivers && (riverPhase > 0 || WorldIdle())) {
                    // --test-rivers: the drainage network. Along every river the water level never rises
                    // downstream and the width never shrinks; walking the longest one from source to sea
                    // through the terrain queries, the water surface only falls and lies in its channel.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (riverPhase == 0) {
                        bool ok = pf != nullptr;
                        auto check = [&](bool c, const char* what) { if (!c) { printf("bench: rivers: FAILED %s\n", what); ok = false; } };
                        const double t0 = GetTime();
                        const RiverNetwork* net = pf ? pf->rivers() : nullptr;
                        const double buildS = GetTime() - t0;
                        check(net != nullptr, "a network");
                        if (net) {
                            const RiverNetwork::Stats s = net->stats();
                            float rise = 0.0f, shrink = 0.0f;
                            const int links = net->checkDownstream(rise, shrink);
                            check(links > 100 && s.mouths > 5, "rivers reach the sea");
                            check(rise <= 1e-4f, "levels never rise downstream (network)");
                            check(shrink <= 1e-4f, "widths never shrink downstream");
                            std::vector<Vector3> dirs; std::vector<float> levels;
                            net->longestPath(5.0f, dirs, levels);
                            // Through the queries: waterElevation along the centreline; the surface may only
                            // fall (a confluence's neighbouring branch may differ by a little), the bed lies
                            // under it, and the path ends at the sea.
                            float worstUp = 0.0f, prev = 1e30f, dryRun = 0.0f, worstDry = 0.0f;
                            int wet = 0;
                            for (size_t i = 0; i < dirs.size(); ++i) {
                                const float w = pf->waterElevation(dirs[i]);
                                const float e = pf->elevation(dirs[i]);
                                if (w != PlanetField::kNoWater) {
                                    if (prev < 1e29f) worstUp = std::max(worstUp, w - prev);
                                    prev = std::min(prev, w);
                                }
                                if (w != PlanetField::kNoWater && e < w - 0.2f) { ++wet; dryRun = 0.0f; }
                                else { dryRun += 5.0f; worstDry = std::max(worstDry, dryRun); }
                            }
                            const float wetShare = dirs.empty() ? 0.0f : (float)wet / dirs.size();
                            check(dirs.size() > 100, "a long river");
                            check(worstUp < 0.05f, "the surface never rises walking downstream");
                            check(wetShare > 0.9f, "water in the channel along the river");
                            check(!dirs.empty() && pf->elevation(dirs.back()) < 3.0f, "it ends at the sea (or the polar ice shelf)");
                            riverLongest = s.longestKm; riverBuild = (float)buildS; riverWet = wetShare;
                            printf("bench: rivers: network built in %.2f s: %d cells of %.0f m, %d land, %d river (%.1f%%), %d lake cells, %d mouths; %d links: worst rise %.5f m, worst narrowing %.5f m\n",
                                   buildS, s.cells, s.cellSize, s.landCells, s.riverCells, 100.0f * s.riverCells / std::max(1, s.landCells), s.lakeCells, s.mouths, links, rise, shrink);
                            printf("bench: rivers: longest %.1f km, walked every 5 m: %.0f%% wet, longest dry run %.0f m, worst rise %.3f m, from %.1f m to %.1f m\n",
                                   s.longestKm, 100.0f * wetShare, worstDry, worstUp, levels.empty() ? 0.0f : levels.front(), levels.empty() ? 0.0f : levels.back());
                            // A look at a river: the nearest to the start (every 150 m within 8 km).
                            {
                                const Vector3 home = pf->findSpawnDir(p->seed);
                                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                                float bestD = 1e30f;
                                for (float e = -8000.0f; e <= 8000.0f; e += 150.0f)
                                    for (float n2 = -8000.0f; n2 <= 8000.0f; n2 += 150.0f) {
                                        const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n2 / pf->radius))));
                                        const RiverNetwork::Hit h = net->nearest(d, 75.0f);
                                        if (!h.found || h.level < 2.0f || e * e + n2 * n2 >= bestD) continue;
                                        bestD = e * e + n2 * n2; riverDir = d; riverLevel = h.level;
                                    }
                                if (bestD < 1e29f) printf("bench: rivers: nearest river to the start %.1f km away, its water at %.1f m\n", sqrtf(bestD) / 1000.0f, riverLevel);
                            }
                        }
                        riverOk = ok;
                        if (ok && Vector3Length(riverDir) > 0.5f) { riverPhase = 1; riverStart = GetTime(); player.flyMode = true; sunOverride = true; }
                        else { printf("bench: rivers: %s\n", ok ? "PASS" : "FAIL"); benchRivers = false; }
                    } else {
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(riverDir, { 0, 1, 0 }));
                        const Vector3 at = Vector3Scale(riverDir, pf->radius + riverLevel);
                        const Vector3 eye = Vector3Add(Vector3Scale(riverDir, pf->radius + riverLevel + 70.0f), Vector3Scale(t1, 110.0f));
                        world.camera.position = p->toWorld(eye);
                        world.camera.target = p->toWorld(at);
                        world.camera.up = p->dirToWorld(riverDir);
                        sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(riverDir, Vector3Scale(t1, -0.6f))));
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - riverStart > 3.0 && WorldIdle()) {
                            TakeScreenshot("bench_river.png");
                            printf("bench: rivers: %s\n", riverOk ? "PASS" : "FAIL");
                            benchRivers = false;
                        }
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
                        {
                            // Water surface smoothness: the level's change between neighbours 4 m apart, where
                            // both are water (a river's surface should be nearly flat across and gently sloped along).
                            double sum = 0.0; int pairs = 0, steep = 0; float worst = 0.0f;
                            auto lvl = [&](float e, float n, float& w) {
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                if (!pf->inlandWaterNear(d, 0.0f)) return false;
                                w = pf->waterElevation(d);
                                return w != PlanetField::kNoWater && w > 0.0f && pf->elevation(d) < w - 0.2f;
                            };
                            for (float e = -2000.0f; e <= 2000.0f; e += 4.0f)
                                for (float n = -2000.0f; n <= 2000.0f; n += 4.0f) {
                                    float a, b, c;
                                    if (!lvl(e, n, a)) continue;
                                    for (int k = 0; k < 2; ++k) {
                                        const float de = k ? 0.0f : 4.0f, dn = k ? 4.0f : 0.0f;
                                        if (!lvl(e + de, n + dn, b) || !lvl(e - de, n - dn, c)) continue;
                                        const float s = fabsf(b - 2.0f * a + c); // bumps and kinks, not an even slope
                                        sum += s; ++pairs; worst = std::max(worst, s); if (s > 0.2f) ++steep;
                                    }
                                }
                            printf("bench: water surface roughness (second difference over 4 m): %d samples, mean %.3f m, worst %.2f m, %.1f%% over 0.2 m\n",
                                   pairs, pairs ? sum / pairs : 0.0, worst, pairs ? 100.0 * steep / pairs : 0.0);
                        }
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
                } else if (benchHunger && (hungerPhase > 0 || WorldIdle())) {
                    // --test-hunger: the food bar drains at its set rate at rest, starving hurts slowly,
                    // eating restores it, the poles are cold and make you hungrier, food is saved.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (hungerPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        hungerPhase = 1; hungerStart = GetTime();
                    } else if (hungerPhase == 1 && GetTime() - hungerStart > 1.5 && player.grounded) {
                        bool ok = true;
                        auto check = [&](bool c, const char* what) { if (!c) { printf("bench: hunger: FAILED %s\n", what); ok = false; } };
                        player.velocity = { 0, 0, 0 };
                        survival.food = 100.0f; survival.health = 100.0f; survival.oxygen = 100.0f; survival.sickUntil = 0.0f;
                        Survival_Update(world, 60.0f, false);
                        const float tempHere = survival.temperature, drained = 100.0f - survival.food;
                        const float expected = 60.0f * survival.foodDrainPerSec * (tempHere < 5.0f ? 1.0f + (5.0f - tempHere) / 15.0f : 1.0f);
                        check(fabsf(drained - expected) < expected * 0.05f + 0.01f, "drain at rest");
                        survival.food = 0.0f; survival.health = 100.0f;
                        Survival_Update(world, 10.0f, false);
                        const float starve = 100.0f - survival.health;
                        check(starve > 4.0f && starve < 6.0f, "starving: 0.5 a second");
                        survival.food = 50.0f; survival.health = 100.0f;
                        inventory.add("berries", 2);
                        Survival_UseItem("berries");
                        check(fabsf(survival.food - (50.0f + Survival_FoodValue("berries"))) < 0.01f && Survival_FoodValue("berries") > 0.0f, "eating berries");
                        // Save and load the food bar.
                        survival.food = 37.0f;
                        const std::string path = std::string(GetWorkingDirectory()) + "/hunger_test.tmp";
                        if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                        survival.food = 100.0f;
                        if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                        remove(path.c_str());
                        check(fabsf(survival.food - 37.0f) < 0.01f, "food saved");
                        hungerOk = ok; hungerTemp = tempHere; hungerDrain = drained;
                        // Then to the pole: cold.
                        player.spawnOn(world, *p, p->dirToWorld(Vector3Normalize({ 0.02f, 1.0f, 0.03f })));
                        hungerPhase = 2; hungerStart = GetTime();
                    } else if (hungerPhase == 2 && GetTime() - hungerStart > 1.5 && player.grounded) {
                        player.velocity = { 0, 0, 0 };
                        survival.food = 100.0f; survival.health = 100.0f;
                        Survival_Update(world, 60.0f, false);
                        const float poleTemp = survival.temperature, poleDrain = 100.0f - survival.food;
                        const bool cold = poleTemp < 0.0f && poleDrain > hungerDrain * 1.5f;
                        printf("bench: hunger: %.1f C at the start, food -%.2f a minute at rest; at the pole %.1f C, -%.2f a minute; starving, eating, saving checked -> %s\n",
                               hungerTemp, hungerDrain, poleTemp, poleDrain, hungerOk && cold ? "PASS" : "FAIL");
                        survival.food = 100.0f; survival.health = 100.0f;
                        benchHunger = false;
                    }
                } else if (benchTools && (toolsPhase > 0 || WorldIdle())) {
                    // --test-tools: a tool wears with each stroke (counted through the swing counter) and
                    // breaks at its durability; repairing restores it for its head material; wear is
                    // saved; an iron hoe tills 3x3 m of farmland in one stroke.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: tools: FAILED %s\n", what); toolsOk = false; } };
                    if (toolsPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        toolsPhase = 1; toolsStart = GetTime(); toolsOk = true;
                    } else if (toolsPhase == 1 && GetTime() - toolsStart > 1.5 && player.grounded) {
                        ToolWear_Reset();
                        inventory.add("stone_pickaxe", 1);
                        terrainTool.heldTool = "stone_pickaxe";
                        terrainTool.swings += 7; // seven strokes, counted by the main loop's wear hook next frame
                        toolsPhase = 2;
                    } else if (toolsPhase == 2) {
                        check(ToolWear_Used("stone_pickaxe") == 7, "strokes wear the held tool");
                        const int life = ToolDurability(ToolOf("stone_pickaxe"));
                        ToolWear_OnStroke("stone_pickaxe", life - 8);
                        check(inventory.count("stone_pickaxe") == 1 && ToolWear_Left("stone_pickaxe") > 0.0f, "one stroke from breaking");
                        // Save and load the wear.
                        const std::string path = std::string(GetWorkingDirectory()) + "/tools_test.tmp";
                        if (FILE* fp = fopen(path.c_str(), "w")) { saveExtrasHook(fp); fprintf(fp, "end_extras\n"); fclose(fp); }
                        ToolWear_Reset();
                        if (FILE* fp = fopen(path.c_str(), "r")) { loadExtrasHook(fp); fclose(fp); }
                        remove(path.c_str());
                        check(ToolWear_Used("stone_pickaxe") == life - 1, "wear saved");
                        // Repair (2 stone for stone tools), then wear it out.
                        const int stone0 = inventory.count("stone");
                        check(!ToolWear_Repair("stone_pickaxe") || stone0 >= 2, "no repair without stone");
                        inventory.add("stone", 2);
                        const int stone1 = inventory.count("stone");
                        check(ToolWear_Repair("stone_pickaxe") && ToolWear_Used("stone_pickaxe") == 0 && inventory.count("stone") == stone1 - 2, "repaired for 2 stone");
                        ToolWear_OnStroke("stone_pickaxe", life);
                        check(inventory.count("stone_pickaxe") == 0, "broke at its durability");
                        toolsLife = life;
                        // An iron hoe: look at the ground just ahead and till.
                        inventory.add("iron_hoe", 1);
                        terrainTool.heldTool = "iron_hoe";
                        const Vector3 level = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        player.forward = Vector3Normalize(Vector3Subtract(level, Vector3Scale(player.up, 0.8f)));
                        player.applyToCamera(world);
                        TerrainTool_Update(world, false, true); // aim only
                        int before = 0;
                        if (p->plots) for (auto& pl : p->plots->plots) for (auto& part : pl->parts) before += part.type == PartType::Farmland ? 1 : 0;
                        const bool stroked = TerrainTool_Stroke(false);
                        int after = 0;
                        if (p->plots) for (auto& pl : p->plots->plots) for (auto& part : pl->parts) after += part.type == PartType::Farmland ? 1 : 0;
                        toolsTilled = after - before;
                        check(stroked && toolsTilled >= 7, "an iron hoe tills 3x3 m");
                        printf("bench: tools: wear counted, saved, repaired, a stone pickaxe broke after %d strokes; an iron hoe tilled %d squares in one stroke -> %s\n",
                               toolsLife, toolsTilled, toolsOk ? "PASS" : "FAIL");
                        terrainTool.heldTool.clear();
                        benchTools = false;
                    }
                } else if (benchSuit && (suitPhase > 0 || WorldIdle())) {
                    // --test-suit: each module changes its stat by the stated amount; underwater and the
                    // cold drain as set; swapping suits keeps the modules that fit and returns the rest;
                    // the setup survives a save and load.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: suit: FAILED %s\n", what); suitOk = false; } };
                    auto near = [](float a, float b) { return fabsf(a - b) < 0.02f * fabsf(b) + 1e-4f; };
                    auto fitOnly = [&](const char* id) {
                        for (int i = 0; i < (int)Suit_Modules().size(); ++i) Suit_Unfit(i);
                        if (id) { inventory.add(id, 1); Suit_Fit(id); }
                        Suit_Apply(world, 0.0f);
                    };
                    auto drainOver = [&](float seconds, float& oxygen, float& food) {
                        survival.oxygen = 100.0f; survival.food = 100.0f; survival.health = 100.0f; survival.sickUntil = 0.0f;
                        player.velocity = { 0, 0, 0 };
                        Survival_Update(world, seconds, false);
                        oxygen = 100.0f - survival.oxygen; food = 100.0f - survival.food;
                    };
                    if (suitPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        suitPhase = 1; suitStart = GetTime(); suitOk = true;
                    } else if (suitPhase == 1 && GetTime() - suitStart > 1.5 && player.grounded) {
                        Suit_Reset();
                        Suit_Apply(world, 0.0f);
                        check(Suit_Worn() == "basic_suit" && Suit_Modules().size() == 2, "the basic suit, 2 slots");
                        const float walk0 = player.walkSpeed, fuel0 = player.jetpackFuelPerSec, thrust0 = player.jetpackAccel;
                        float o, f0, f;
                        drainOver(60.0f, o, f0);
                        // Wear the explorer suit: 6 slots, 1.5x oxygen and power.
                        inventory.add("explorer_suit", 1);
                        check(Suit_Wear("explorer_suit") && Suit_Modules().size() == 6 && near(Suit_Stats().oxygen, 1.5f), "the explorer suit");
                        check(inventory.count("basic_suit") == 1, "the old suit back in the inventory");
                        fitOnly("servo_legs");     check(near(player.walkSpeed, walk0 * 1.2f), "servo legs: +20% walk");
                        fitOnly("thruster_pack");  check(near(player.jetpackFuelPerSec, fuel0 / 1.6f) && near(player.jetpackAccel, thrust0 * 1.6f), "thruster pack: +60% fuel and thrust");
                        fitOnly("oxygen_tank_module"); check(near(Suit_Stats().oxygen, 2.0f), "oxygen tank: +50%");
                        fitOnly("battery_pack");   check(near(Suit_Stats().power, 2.0f), "battery pack: +50%");
                        fitOnly("lamp_upgrade");   check(near(Suit_Stats().lamp, 2.0f), "lamp upgrade");
                        fitOnly("grip_boots");     check(player.walkableCos < 0.4f, "grip boots");
                        fitOnly("food_pouch");     drainOver(60.0f, o, f); check(near(f, f0 * 0.75f), "food pouch: hunger 25% slower");
                        suitFoodRatio = f / f0;
                        fitOnly("ore_scanner");    check(Suit_Stats().scanner, "ore scanner");
                        fitOnly(nullptr);
                        check(near(player.walkSpeed, walk0) && near(Suit_Stats().oxygen, 1.5f), "modules out: the suit's own figures");
                        // Fill all six, swap to the basic suit: two stay, four come back.
                        const char* six[6] = { "servo_legs", "thruster_pack", "oxygen_tank_module", "battery_pack", "lamp_upgrade", "food_pouch" };
                        int before = 0;
                        for (const char* id : six) before += inventory.count(id);
                        for (const char* id : six) { inventory.add(id, 1); Suit_Fit(id); }
                        int fittedN = 0; for (const auto& m : Suit_Modules()) fittedN += m.empty() ? 0 : 1;
                        check(fittedN == 6, "six modules fitted");
                        check(Suit_Wear("basic_suit"), "back to the basic suit");
                        fittedN = 0; for (const auto& m : Suit_Modules()) fittedN += m.empty() ? 0 : 1;
                        int after = 0;
                        for (const char* id : six) after += inventory.count(id);
                        check(fittedN == 2 && after == before + 4 && inventory.count("explorer_suit") == 1, "two kept, four and the explorer suit returned");
                        // Save and load.
                        const std::string w = Suit_Worn(); const std::vector<std::string> mods = Suit_Modules();
                        const std::string path = std::string(GetWorkingDirectory()) + "/suit_test.tmp";
                        if (FILE* fp = fopen(path.c_str(), "w")) { saveExtrasHook(fp); fprintf(fp, "end_extras\n"); fclose(fp); }
                        Suit_Reset();
                        if (FILE* fp = fopen(path.c_str(), "r")) { loadExtrasHook(fp); fclose(fp); }
                        remove(path.c_str());
                        check(Suit_Worn() == w && Suit_Modules() == mods, "the suit saved");
                        Suit_Reset(); Suit_Apply(world, 0.0f);
                        // Next: under the sea.
                        Vector3 sea{ 0, 0, 0 };
                        for (int i = 0; i < 4000 && sea.x == 0.0f && sea.y == 0.0f; ++i) {
                            const Vector3 d = Vector3Normalize({ sinf(i * 2.39996f) * sqrtf(1.0f - powf(1.0f - 2.0f * (i + 0.5f) / 4000.0f, 2.0f)), 1.0f - 2.0f * (i + 0.5f) / 4000.0f, cosf(i * 2.39996f) * sqrtf(1.0f - powf(1.0f - 2.0f * (i + 0.5f) / 4000.0f, 2.0f)) });
                            const float e = pf->elevation(d);
                            if (e < -8.0f && e > -30.0f && fabsf(d.y) < 0.6f) sea = d;
                        }
                        check(Vector3Length(sea) > 0.5f, "a shallow sea found");
                        player.spawnOn(world, *p, p->dirToWorld(sea));
                        suitPhase = 2; suitStart = GetTime();
                    } else if (suitPhase == 2 && GetTime() - suitStart > 2.0) {
                        float base, rebreathed, f;
                        fitOnly(nullptr);
                        drainOver(10.0f, base, f);
                        check(survival.underwater, "the head is underwater");
                        fitOnly("rebreather");
                        drainOver(10.0f, rebreathed, f);
                        suitWater = base; suitRebreather = rebreathed;
                        check(near(base, survival.oxygenDrainPerSec * 2.0f * 10.0f), "underwater: oxygen drains twice as fast");
                        check(near(rebreathed, base * 0.25f), "rebreather: 4x longer underwater");
                        fitOnly(nullptr);
                        player.spawnOn(world, *p, p->dirToWorld(Vector3Normalize({ 0.02f, 1.0f, 0.03f })));
                        suitPhase = 3; suitStart = GetTime();
                    } else if (suitPhase == 3 && GetTime() - suitStart > 1.5 && player.grounded) {
                        float o, cold, insulated, regulated;
                        fitOnly(nullptr);            drainOver(60.0f, o, cold);
                        const float T = survival.temperature;
                        fitOnly("insulation_module"); drainOver(60.0f, o, insulated);
                        fitOnly("thermal_regulator"); drainOver(60.0f, o, regulated);
                        const float base = 60.0f * survival.foodDrainPerSec;
                        check(T < 0.0f, "the pole is cold");
                        check(near(cold, base * (1.0f + (5.0f - T) / 15.0f)), "cold: hungrier");
                        check(near(insulated, base * (1.0f + 0.4f * (5.0f - T) / 15.0f)), "insulation: 40% of the cold felt");
                        check(near(regulated, base), "thermal regulator: no cold");
                        fitOnly(nullptr);
                        printf("bench: suit: modules checked (food pouch x%.2f), suits swapped, saved; underwater O2 -%.1f in 10 s, %.2f with a rebreather; at the pole (%.0f C) food -%.2f, insulated -%.2f, regulated -%.2f a minute -> %s\n",
                               suitFoodRatio, suitWater, suitRebreather, T, cold, insulated, regulated, suitOk ? "PASS" : "FAIL");
                        Suit_Reset();
                        survival.food = 100.0f; survival.oxygen = 100.0f; survival.health = 100.0f;
                        benchSuit = false;
                    }
                } else if (benchCrafting && (craftPhase > 0 || WorldIdle())) {
                    // --test-crafting: recipes are known once their ingredients are discovered (and any unlock
                    // is met); search; a chest beside the player counts for crafting and is pulled from; the
                    // inventory deposits into chests that already hold an item.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: crafting: FAILED %s\n", what); craftOk = false; } };
                    auto recipeFor = [](const char* out) -> const Recipe* { for (const Recipe& r : itemDB.recipes()) if (r.output == out) return &r; return nullptr; };
                    if (craftPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        craftPhase = 1; craftStart = GetTime(); craftOk = true;
                    } else if (craftPhase == 1 && GetTime() - craftStart > 1.5 && player.grounded) {
                        const Recipe* rod = recipeFor("steel_rod");
                        const Recipe* suit = recipeFor("explorer_suit");
                        check(rod && suit, "recipes found");
                        if (rod && suit) {
                            check(!RecipeBook_Known(*rod), "steel rods hidden before steel is found");
                            inventory.add("steel_ingot", 1);
                            check(RecipeBook_Known(*rod), "finding steel reveals steel rods");
                            check(RecipeBook_Matches(*rod, "rod") && RecipeBook_Matches(*rod, "steel ingot") && !RecipeBook_Matches(*rod, "zzz"), "search by name or ingredient");
                            inventory.add("titanium_alloy_ingot", 1); inventory.add("advanced_circuit", 1); inventory.add("cloth", 1);
                            check(!RecipeBook_Known(*suit), "the explorer suit also needs its unlock (three nights)");
                            QuestBook_Count("nights", 3.0);
                            check(RecipeBook_Known(*suit), "unlocked after three nights");
                        }
                        // A chest 2 m ahead.
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        const Vector3 d = Vector3Normalize(Vector3Add(me, Vector3Scale(Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), 2.5f)));
                        const Vector3 at = Vector3Scale(d, pf->radius + pf->elevation(d));
                        craftPlot = p->plots->plotFor(at, 0.0f, true);
                        craftChest = craftPlot->place(Part{ PartType::Chest, craftPlot->cellAt(craftPlot->toPlot(Vector3Add(at, Vector3Scale(d, 0.1f)))), 0, 0 });
                        check(craftChest >= 0 && Storage_NearbyChests(world).size() == 1, "a chest within reach");
                        craftPhase = 2; craftStart = GetTime();
                    } else if (craftPhase == 2 && GetTime() - craftStart > 1.0) {
                        Factory::Seed(p->id, craftPlot->id, craftChest, "iron_ingot", 4);
                        craftPhase = 3; craftStart = GetTime();
                    } else if (craftPhase == 3 && GetTime() - craftStart > 0.5) {
                        check(Storage_CountNearby(world, "iron_ingot") == 4, "the chest's contents counted");
                        craftInv0 = inventory.count("iron_ingot");
                        check(Storage_Pull(world, "iron_ingot", 3), "pull 3 iron from the chest");
                        craftPhase = 4; craftStart = GetTime();
                    } else if (craftPhase == 4 && GetTime() - craftStart > 0.5) {
                        check(inventory.count("iron_ingot") == craftInv0 + 3 && Storage_CountNearby(world, "iron_ingot") == 1, "the iron arrived; one left in the chest");
                        inventory.add("oak_planks", 5);
                        craftWood0 = inventory.count("oak_planks");
                        craftDeposited = Storage_Deposit(world);
                        check(craftDeposited == 1 && inventory.count("iron_ingot") == 0 && inventory.count("oak_planks") == craftWood0, "deposit: iron into the chest that holds iron; planks stay");
                        craftPhase = 5; craftStart = GetTime();
                    } else if (craftPhase == 5 && GetTime() - craftStart > 0.5) {
                        craftInChest = Storage_CountNearby(world, "iron_ingot");
                        check(craftInChest == craftInv0 + 4, "the chest holds the deposited iron");
                        if (craftPlot && craftChest >= 0) craftPlot->remove(craftChest);
                        printf("bench: crafting: discovery reveals recipes, an unlock gates one, search works; a nearby chest counted, pulled from, deposited into (%d iron there) -> %s\n",
                               craftInChest, craftOk ? "PASS" : "FAIL");
                        benchCrafting = false;
                    }
                } else if (benchDamage && (damagePhase > 0 || WorldIdle())) {
                    // --test-damage (plan 3, M18): rocket part damage lowers thrust, a broken tank leaks, broken
                    // legs stop counting; damage is saved with the plot; a repair kit mends the worst part;
                    // vacuum hurts after 60 s in an unrated suit, never in the explorer suit.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: damage: FAILED %s\n", what); damageOk = false; } };
                    if (damagePhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        damagePhase = 1; damageStart = GetTime(); damageOk = true;
                    } else if (damagePhase == 1 && GetTime() - damageStart > 1.5 && player.grounded) {
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        const Vector3 d = Vector3Normalize(Vector3Add(me, Vector3Scale(Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), 4.0f)));
                        Plot* r = p->plots->plotFor(Vector3Scale(d, pf->radius + pf->elevation(d)), 0.0f, true);
                        int legs[4] = { -1, -1, -1, -1 };
                        for (int k = 0; k < 4; ++k) legs[k] = r->place(Part{ PartType::LandingLeg, { k & 1 ? 1 : -2, 0, k & 2 ? 1 : -2 }, 0, 0 });
                        const int engine = r->place(Part{ PartType::RocketEngine, { -1, 0, -1 }, 0, 0 });
                        const int tank = r->place(Part{ PartType::FuelTank, { -1, 2, -1 }, 0, 0 });
                        r->place(Part{ PartType::Guidance, { 0, 4, 0 }, 0, 0 });
                        r->place(Part{ PartType::Cockpit, { -1, 5, -1 }, 0, 0 });
                        check(engine >= 0 && tank >= 0 && legs[0] >= 0, "a rocket built");
                        const RocketStats s0 = Rocket_Stats(*r, p);
                        Rocket_DamagePart(*r, engine, 0.5f);
                        Rocket_DamagePart(*r, tank, 1.0f);
                        Rocket_DamagePart(*r, legs[0], 1.0f);
                        const RocketStats s1 = Rocket_Stats(*r, p);
                        check(fabsf(s1.thrust - s0.thrust * 0.5f) < 1.0f, "half-damaged engine: half thrust");
                        check(s1.leak > 1.9f && s1.leak < 2.1f, "a broken tank leaks 2 kg/s");
                        check(s1.legs == s0.legs - 1, "a broken leg stops counting");
                        check(s1.damagedParts == 3 && s1.damage > 0.0f, "damage counted");
                        damageThrust = s1.thrust / std::max(1.0f, s0.thrust);
                        // Saved with the plot.
                        const std::string path = std::string(GetWorkingDirectory()) + "/damage_test.tmp";
                        if (FILE* f = fopen(path.c_str(), "w")) { p->plots->save(f); fclose(f); }
                        PlotSet back;
                        bool loaded = false;
                        if (FILE* f = fopen(path.c_str(), "r")) { loaded = back.load(f); fclose(f); }
                        remove(path.c_str());
                        const Plot* rb = loaded ? back.byId(r->id) : nullptr;
                        check(rb && rb->damage.count(engine) && fabsf(rb->damage.at(engine) - 0.5f) < 0.01f && rb->damage.count(tank), "damage saved with the plot");
                        // A repair kit: the worst part (the tank or the leg, both broken) mended by half.
                        inventory.add("repair_kit", 1);
                        const bool used = Survival_UseItem("repair_kit");
                        const RocketStats s2 = Rocket_Stats(*r, p);
                        check(used && inventory.count("repair_kit") == 0 && s2.damage < s1.damage - 0.05f, "a repair kit mends the worst part");
                        damageRepaired = s1.damage - s2.damage;
                        for (int k = (int)r->parts.size() - 1; k >= 0; --k) r->remove(k);
                        // Vacuum: high above the air in the basic suit, then the explorer suit.
                        const DVec3 keep = player.position;
                        player.position = p->worldPos + DVec3(Vector3Scale(p->dirToWorld(up), pf->radius + 40000.0f));
                        Suit_Reset(); Suit_Apply(world, 0.0f);
                        survival.health = 100.0f; survival.oxygen = 100.0f; survival.vacuumTime = 0.0f;
                        Survival_Update(world, 30.0f, false);
                        const float after30 = survival.health;
                        Survival_Update(world, 40.0f, false);
                        damageVacuum = 100.0f - survival.health;
                        check(survival.inVacuum && after30 >= 99.9f && damageVacuum > 10.0f, "vacuum: 60 s, then it hurts (basic suit)");
                        inventory.add("explorer_suit", 1);
                        Suit_Wear("explorer_suit"); Suit_Apply(world, 0.0f);
                        survival.health = 100.0f; survival.oxygen = 100.0f; survival.vacuumTime = 0.0f;
                        Survival_Update(world, 30.0f, false); Survival_Update(world, 40.0f, false);
                        check(survival.health >= 99.9f, "the explorer suit is vacuum-rated");
                        player.position = keep;
                        Suit_Reset(); Suit_Apply(world, 0.0f);
                        survival.health = 100.0f; survival.oxygen = 100.0f; survival.vacuumTime = 0.0f;
                        printf("bench: damage: engine at half thrust (%.2f), tank leaking, leg out, saved; a repair kit mended %.2f; vacuum cost %.0f health past 60 s in the basic suit, none in the explorer suit -> %s\n",
                               damageThrust, damageRepaired, damageVacuum, damageOk ? "PASS" : "FAIL");
                        benchDamage = false;
                    }
                } else if (benchMaterials && (matPhase > 0 || WorldIdle())) {
                    // --test-materials (plan 3, M21): trees give logs of their species (planks), rock gives the
                    // stone of its region and depth, the detail shapes (slab, panel, pane, post) build in the new
                    // materials and keep them through a save; a screenshot of a little hut built from them.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: materials: FAILED %s\n", what); matOk = false; } };
                    if (matPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        matPhase = 1; matStart = GetTime(); matOk = true;
                    } else if (matPhase == 1 && GetTime() - matStart > 1.5 && player.grounded) {
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        // Logs from a tree, sawn into planks.
                        int logs = 0;
                        for (int k = 0; k < 4 && logs == 0; ++k) {
                            TerrainTool_HarvestNearest(*p, me, true);
                            logs = inventory.count("oak_log") + inventory.count("birch_log") + inventory.count("pine_log");
                        }
                        check(logs > 0, "a tree gives logs of its species");
                        int planks = 0;
                        for (const char* l : { "oak", "birch", "pine" }) {
                            const std::string log = std::string(l) + "_log", plank = std::string(l) + "_planks";
                            for (const Recipe& r : itemDB.recipes()) if (r.output == plank && inventory.count(log) > 0) inventory.craft(r);
                            planks += inventory.count(plank);
                        }
                        check(planks >= 4, "logs saw into planks");
                        // Stone kinds over the region and with depth.
                        std::set<std::string> kinds;
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
                        for (int i = 0; i < 400; ++i) {
                            const float a = ((i * 37) % 200 - 100) * 200.0f, b = ((i * 91) % 200 - 100) * 200.0f; // +-20 km
                            const Vector3 d = Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, a / pf->radius), Vector3Scale(t2, b / pf->radius))));
                            const float surf = pf->radius + pf->elevation(d);
                            for (float depth : { 5.0f, 30.0f, 80.0f }) kinds.insert(pf->stoneKind(Vector3Scale(d, surf - depth)));
                        }
                        std::string list;
                        for (const auto& k : kinds) list += (list.empty() ? "" : ", ") + k;
                        check(kinds.size() >= 4 && kinds.count("basalt"), "several stone kinds (basalt deep down)");
                        matKinds = list;
                        // Digging rock gives its kind too.
                        double kindsBefore = 0.0;
                        for (const auto& k : kinds) kindsBefore += inventory.obtainedTotal(k);
                        inventory.add("stone_pickaxe", 1);
                        terrainTool.heldTool = "stone_pickaxe";
                        for (int k = 0; k < 6; ++k) TerrainTool_DigAt(*p, Vector3Subtract(Vector3Scale(up, pf->radius + pf->elevation(up)), Vector3Scale(up, 6.0f + k * 0.6f)));
                        double kindsAfter = 0.0;
                        for (const auto& k : kinds) kindsAfter += inventory.obtainedTotal(k);
                        check(kindsAfter > kindsBefore, "digging rock gives its stone kind");
                        terrainTool.heldTool.clear();
                        // A little hut 6 m ahead: a marble slab floor, oak panel walls, glass panes, granite posts.
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        const Vector3 hd = Vector3Normalize(Vector3Add(me, Vector3Scale(t1, 6.0f)));
                        Plot* hut = p->plots->plotFor(Vector3Scale(hd, pf->radius + pf->elevation(hd)), 0.0f, true);
                        auto mat = [](const char* item) { return (uint8_t)std::max(0, BuildMaterialForItem(item)); };
                        int placed = 0, tried = 0;
                        auto put = [&](PartType t, Int3 c, int orient, const char* m) { ++tried; placed += hut->place(Part{ t, c, (uint8_t)orient, mat(m) }) >= 0 ? 1 : 0; };
                        for (int x = 0; x < 6; ++x) for (int z = 0; z < 6; ++z) put(PartType::Slab, { x, 0, z }, 0, "marble");
                        for (int y = 1; y <= 4; ++y) {
                            for (int x = 1; x < 5; ++x) { put(y == 2 || y == 3 ? PartType::Pane : PartType::Panel, { x, y, 0 }, RotateOrientation(0, 1), y == 2 || y == 3 ? "glass" : "oak_planks"); put(PartType::Panel, { x, y, 5 }, RotateOrientation(0, 1), "birch_planks"); }
                            for (int z = 1; z < 5; ++z) { put(PartType::Panel, { 0, y, z }, 0, "pine_planks"); put(PartType::Panel, { 5, y, z }, 0, "oak_planks"); }
                            put(PartType::Post, { 0, y, 0 }, 0, "granite"); put(PartType::Post, { 5, y, 0 }, 0, "granite"); put(PartType::Post, { 0, y, 5 }, 0, "granite"); put(PartType::Post, { 5, y, 5 }, 0, "granite");
                        }
                        for (int x = 0; x < 6; ++x) for (int z = 0; z < 6; ++z) put(PartType::Slab, { x, 5, z }, 0, (x + z) % 2 ? "brick" : "basalt");
                        check(placed == tried, "the hut's parts placed");
                        bool partial = true;
                        for (const auto& kv : hut->cells) partial = partial && !kv.second.shape.full;
                        check(partial && !hut->cells.empty(), "detail shapes fill part of a cell");
                        // Save and load: the materials come back.
                        const std::string path = std::string(GetWorkingDirectory()) + "/materials_test.tmp";
                        if (FILE* f = fopen(path.c_str(), "w")) { p->plots->save(f); fclose(f); }
                        PlotSet back; bool loaded = false;
                        if (FILE* f = fopen(path.c_str(), "r")) { loaded = back.load(f); fclose(f); }
                        remove(path.c_str());
                        const Plot* hb = loaded ? back.byId(hut->id) : nullptr;
                        bool same = hb && hb->parts.size() == hut->parts.size();
                        for (size_t i = 0; same && i < hut->parts.size(); ++i) same = hb->parts[i].type == hut->parts[i].type && hb->parts[i].material == hut->parts[i].material;
                        check(same, "the hut saved with its materials");
                        matPlaced = placed;
                        // Look at the hut.
                        const Vector3 hutC = hut->toBody({ 1.5f, 1.2f, 1.5f });
                        world.camera.position = p->toWorld(Vector3Add(hutC, Vector3Add(Vector3Scale(up, 3.5f), Vector3Add(Vector3Scale(t1, -4.0f), Vector3Scale(t2, -5.0f)))));
                        world.camera.target = p->toWorld(hutC);
                        world.camera.up = p->dirToWorld(up);
                        player.flyMode = true;
                        player.position = world.renderOrigin + world.camera.position;
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 0.5f))));
                        matCamPos = Vector3Add(hutC, Vector3Add(Vector3Scale(up, 3.5f), Vector3Add(Vector3Scale(t1, -4.0f), Vector3Scale(t2, -5.0f)))); matCamTarget = hutC; matCamUp = up; // body-local
                        matPhase = 2; matStart = GetTime();
                    } else if (matPhase == 2 && GetTime() - matStart <= 2.0) {
                        world.camera.position = p->toWorld(matCamPos); world.camera.target = p->toWorld(matCamTarget); world.camera.up = p->dirToWorld(matCamUp); // hold the view on the hut
                        player.position = world.renderOrigin + world.camera.position;
                    } else if (matPhase == 2 && GetTime() - matStart > 2.0) {
                        world.camera.position = p->toWorld(matCamPos); world.camera.target = p->toWorld(matCamTarget); world.camera.up = p->dirToWorld(matCamUp);
                        TakeScreenshot("bench_materials.png");
                        printf("bench: materials: logs and planks of three woods; stone kinds %s; a hut of %d detail parts in 8 materials, saved -> %s\n", matKinds.c_str(), matPlaced, matOk ? "PASS" : "FAIL");
                        benchMaterials = false;
                    }
                } else if (benchDecor && (decorPhase > 0 || WorldIdle())) {
                    // --test-decor (plan 3, M22): rock gives gems by its stone kind; dyes, plaster in eight colours,
                    // stained glass, lamps and gem blocks craft and build; a screenshot of a decorated wall.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: decor: FAILED %s\n", what); decorOk = false; } };
                    if (decorPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        decorPhase = 1; decorStart = GetTime(); decorOk = true;
                    } else if (decorPhase == 1 && GetTime() - decorStart > 1.5 && player.grounded) {
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
                        // A spot of each stone kind within 20 km, then dig it until its gem turns up.
                        std::map<std::string, Vector3> spot;
                        for (int i = 0; i < 2000 && spot.size() < 5; ++i) {
                            const float a = ((i * 37) % 200 - 100) * 200.0f, b = ((i * 91) % 200 - 100) * 200.0f;
                            const Vector3 d = Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, a / pf->radius), Vector3Scale(t2, b / pf->radius))));
                            const float surf = pf->radius + pf->elevation(d);
                            for (float depth : { 8.0f, 30.0f, 80.0f }) { const Vector3 q = Vector3Scale(d, surf - depth); spot.emplace(pf->stoneKind(q), q); }
                        }
                        inventory.add("stone_pickaxe", 1);
                        terrainTool.heldTool = "stone_pickaxe";
                        const char* gems[5] = { "quartz", "amethyst", "ruby", "sapphire", "emerald" };
                        terrainTool.gemBoost = 75; // gems are rare in play (about 1 in 3000 strokes); the test checks each stone holds one
                        int g0[5]; for (int k = 0; k < 5; ++k) g0[k] = inventory.count(gems[k]);
                        int digs = 0;
                        for (const auto& kv : spot) {
                            const Vector3 dir = Vector3Normalize(kv.second);
                            for (int k = 0; k < 160; ++k, ++digs) TerrainTool_DigAt(*p, Vector3Subtract(kv.second, Vector3Scale(dir, 0.5f * (k % 40))));
                        }
                        terrainTool.heldTool.clear();
                        terrainTool.gemBoost = 1;
                        int kinds = 0;
                        for (int k = 0; k < 5; ++k) { decorGems += inventory.count(gems[k]) - g0[k]; kinds += inventory.count(gems[k]) > g0[k] ? 1 : 0; }
                        decorGemKinds = kinds; decorDigs = digs;
                        check(kinds >= 3, "rock gives gems of several kinds");
                        // Dyes, plaster, stained glass, gem blocks, lamps.
                        inventory.add("limestone", 40); inventory.add("sand", 30); inventory.add("charcoal", 6); inventory.add("berries", 24); inventory.add("grain", 16);
                        inventory.add("fibre", 18); inventory.add("flax_seeds", 24); inventory.add("roots", 8); inventory.add("glass", 24); inventory.add("power_cell", 2);
                        for (const char* g : gems) inventory.add(g, 4);
                        auto make = [&](const char* item, int n) {
                            for (int tries = 0; tries < 24 && inventory.count(item) < n; ++tries) {
                                bool any = false;
                                for (const Recipe& r : itemDB.recipes()) if (r.output == item && inventory.canCraft(r)) { inventory.craft(r); any = true; break; }
                                if (!any) break;
                            }
                            return inventory.count(item) >= n;
                        };
                        bool blocks = true; for (const char* g : gems) { const bool ok1 = make(TextFormat("%s_block", g), 1); if (!ok1) printf("bench: decor: no %s block (%d %s)\n", g, inventory.count(g), g); blocks = ok1 && blocks; }
                        check(blocks, "gem blocks");
                        const char* colours[8] = { "white", "black", "red", "orange", "yellow", "green", "blue", "purple" };
                        bool allDyes = true, allPlaster = true;
                        for (const char* c : colours) allDyes = make(TextFormat("%s_dye", c), c[0] == 'w' ? 10 : 6) && allDyes;
                        allPlaster = make("white_plaster", 72);
                        for (int k = 1; k < 8; ++k) allPlaster = make(TextFormat("%s_plaster", colours[k]), 4) && allPlaster;
                        check(allDyes, "every dye crafts");
                        check(allPlaster, "plaster in eight colours");
                        inventory.clear(); // room for the next crafts (a full 45-slot grid refuses outputs)
                        inventory.add("glass", 24); inventory.add("quartz", 2); inventory.add("power_cell", 2);
                        for (const char* c : { "red", "yellow", "green", "blue" }) inventory.add(TextFormat("%s_dye", c), 2);
                        bool glass = true; for (const char* c : { "red", "yellow", "green", "blue" }) glass = make(TextFormat("%s_glass", c), 2) && glass;
                        check(glass, "stained glass");
                        check(make("lamp", 2), "lamps");
                        // A wall 6 m ahead: plaster stripes, a stained-glass window, lamps and gem blocks.
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        const Vector3 wd = Vector3Normalize(Vector3Add(me, Vector3Scale(t1, 6.0f)));
                        Plot* wall = p->plots->plotFor(Vector3Scale(wd, pf->radius + pf->elevation(wd)), 0.0f, true);
                        auto mat = [](const char* item) { return (uint8_t)std::max(0, BuildMaterialForItem(item)); };
                        int placed = 0, tried = 0;
                        auto put = [&](PartType t, Int3 c, const char* m) { ++tried; placed += wall->place(Part{ t, c, (uint8_t)RotateOrientation(0, 1), mat(m) }) >= 0 ? 1 : 0; };
                        for (int x = 0; x < 8; ++x) for (int y = 0; y < 2; ++y) put(PartType::Panel, { x, y, 0 }, TextFormat("%s_plaster", colours[x]));
                        const char* gl[4] = { "red_glass", "yellow_glass", "green_glass", "blue_glass" };
                        for (int x = 0; x < 8; ++x) put(PartType::Pane, { x, 2, 0 }, gl[x % 4]);
                        for (int x = 0; x < 8; ++x) put(PartType::Panel, { x, 3, 0 }, x % 2 ? "lamp" : "white_plaster");
                        const char* gb[5] = { "quartz_block", "amethyst_block", "ruby_block", "sapphire_block", "emerald_block" };
                        for (int x = 0; x < 5; ++x) put(PartType::Cube, { 1 + x, 4, 0 }, gb[x]);
                        check(placed == tried, "the wall built");
                        decorPlaced = placed;
                        const Vector3 wc = wall->toBody({ 2.0f, 1.2f, 0.0f });
                        decorCamPos = Vector3Add(wc, Vector3Add(Vector3Scale(up, 0.8f), Vector3Scale(Vector3Normalize(Vector3Subtract(wall->toBody({ 2.0f, 1.2f, -1.0f }), wc)), 5.5f)));
                        decorCamTarget = wc; decorCamUp = up;
                        player.flyMode = true;
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(Vector3Normalize(Vector3Subtract(decorCamPos, wc)), 0.6f))));
                        decorPhase = 2; decorStart = GetTime();
                    } else if (decorPhase == 2) {
                        world.camera.position = p->toWorld(decorCamPos); world.camera.target = p->toWorld(decorCamTarget); world.camera.up = p->dirToWorld(decorCamUp);
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - decorStart > 2.0) {
                            TakeScreenshot("bench_decor.png");
                            printf("bench: decor: %d gems of %d kinds in %d digs; 8 dyes, 8 plasters, 4 stained glasses, lamps, 5 gem blocks; a wall of %d parts -> %s\n",
                                   decorGems, decorGemKinds, decorDigs, decorPlaced, decorOk ? "PASS" : "FAIL");
                            benchDecor = false;
                        }
                    }
                } else if (benchWeather && (weatherPhase > 0 || WorldIdle())) {
                    // --test-weather (plan 3, M23): the forecast over the planet gives every kind of weather (snow
                    // toward the poles), is the same for the same place and time, and moves on over a day; rain
                    // greys the sky and chills the air; screenshots of rain, snow, a storm and fog.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: weather: FAILED %s\n", what); weatherOk = false; } };
                    if (weatherPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        weatherPhase = 1; weatherStart = GetTime(); weatherOk = true;
                    } else if (weatherPhase == 1 && GetTime() - weatherStart > 1.5 && player.grounded) {
                        int count[6] = {}, polarSnow = 0, polar = 0;
                        const double t0 = WorldTime();
                        for (int i = 0; i < 3000; ++i) {
                            const float z = 1.0f - 2.0f * (i + 0.5f) / 3000.0f, r = sqrtf(1.0f - z * z), a = i * 2.39996f;
                            const Vector3 d = { r * cosf(a), z, r * sinf(a) };
                            const WeatherState w = Weather_At(world, d, t0 + (i % 7) * 3600.0);
                            ++count[(int)w.kind];
                            if (fabsf(z) > 0.8f) { ++polar; polarSnow += w.kind == WeatherKind::Snow ? 1 : 0; }
                        }
                        for (int k = 0; k < 6; ++k) weatherCounts[k] = count[k];
                        check(count[0] > 0 && count[1] > 0 && count[2] > 0, "clear, overcast and rain all occur");
                        check(count[4] > 0 && polarSnow > 0, "snow (toward the poles)");
                        check(count[3] > 0, "storms occur");
                        check(count[5] > 0, "fog occurs");
                        const Vector3 me = Vector3Normalize(p->toLocal((player.position - world.renderOrigin).toFloat()));
                        const WeatherState a1 = Weather_At(world, me, t0), a2 = Weather_At(world, me, t0);
                        check(a1.cloud == a2.cloud && a1.kind == a2.kind, "the same place and time, the same weather");
                        std::set<int> overDay;
                        for (int h = 0; h < 24; ++h) overDay.insert((int)Weather_At(world, me, t0 + h * 600.0).kind);
                        weatherChanges = (int)overDay.size();
                        check(weatherChanges >= 2, "the weather moves on over a day");
                        // Forced rain: the sky greys and the air chills.
                        Weather_Force((int)WeatherKind::Clear); Weather_Update(world, 1.0f); Survival_Update(world, 0.01f, false);
                        const float tClear = survival.temperature;
                        Weather_Force((int)WeatherKind::Rain); Weather_Update(world, 1.0f); Survival_Update(world, 0.01f, false);
                        const float tRain = survival.temperature;
                        check(world.weatherCloud > 0.6f && tRain < tClear - 2.0f, "rain greys the sky and chills the air");
                        weatherChill = tClear - tRain;
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(me, { 0.0f, 0.15f, 0.0f }))); // mid-morning for the pictures
                        weatherPhase = 2; weatherStart = GetTime();
                    } else if (weatherPhase >= 2 && weatherPhase <= 5) {
                        // A look up the slope in each weather (the screenshot after 2 s of each).
                        const WeatherKind shots[4] = { WeatherKind::Rain, WeatherKind::Snow, WeatherKind::Storm, WeatherKind::Fog };
                        const char* names[4] = { "bench_weather_rain.png", "bench_weather_snow.png", "bench_weather_storm.png", "bench_weather_fog.png" };
                        const int k = weatherPhase - 2;
                        Weather_Force((int)shots[k]);
                        if (GetTime() - weatherStart > 2.0) {
                            TakeScreenshot(names[k]);
                            ++weatherPhase; weatherStart = GetTime();
                        }
                    } else if (weatherPhase == 6) {
                        Weather_Force(-1);
                        printf("bench: weather: over 3000 places clear %d, overcast %d, rain %d, storm %d, snow %d, fog %d; %d kinds at the start over a day; rain %.1f C colder -> %s\n",
                               weatherCounts[0], weatherCounts[1], weatherCounts[2], weatherCounts[3], weatherCounts[4], weatherCounts[5], weatherChanges, weatherChill, weatherOk ? "PASS" : "FAIL");
                        benchWeather = false;
                    }
                } else if (benchRecipes) {
                    // --test-recipes: every item can be had from what the world gives (digging, chopping,
                    // bushes, grass, animals, crops) through hand recipes, machines (once their own item
                    // can be made) and the starting kit. Lists what cannot, and items nothing uses.
                    auto split = [](const std::string& s, char sep) {
                        std::vector<std::string> out; size_t a = 0, b;
                        while ((b = s.find(sep, a)) != std::string::npos) { out.push_back(s.substr(a, b - a)); a = b + 1; }
                        out.push_back(s.substr(a));
                        for (auto& c : out) { const size_t x = c.find_first_not_of(" \t"), y = c.find_last_not_of(" \t\r\n"); c = x == std::string::npos ? "" : c.substr(x, y - x + 1); }
                        return out;
                    };
                    auto ids = [&](const std::string& list) { std::vector<std::string> out; for (const std::string& e : split(list, ',')) { char id[64] = {}; if (e != "-" && sscanf(e.c_str(), "%63s", id) == 1) out.push_back(id); } return out; };
                    std::set<std::string> have = { "anorthite", "nickel_iron_ore", "platinum_ore", "iridium_ore", "maple_log", "spruce_log", "fibre", "berries", "grain_seeds", "root_seeds", "flax_seeds", "dirt", "gravel", "stone", "sand", "clay",
                                                   "ice", "regolith", "iron_ore", "copper_ore", "tin_ore", "coal", "bauxite_ore", "ilmenite_ore", "telescope", "map", "basic_suit",
                                                   "oak_log", "birch_log", "pine_log", "granite", "basalt", "sandstone", "limestone", "marble",
                                                   "quartz", "amethyst", "ruby", "sapphire", "emerald" };
                    std::set<std::string> used;
                    struct Rule { std::string machine; std::vector<std::string> in, out; };
                    std::vector<Rule> rules;
                    for (const Recipe& r : itemDB.recipes()) { Rule u; u.out = { r.output }; for (const auto& in : r.inputs) u.in.push_back(in.first); rules.push_back(u); }
                    for (const ItemDef& d : itemDB.items()) // wildcards: any planks (any log) is had once one kind is
                        for (const char* w : { "any_planks", "any_log" }) if (d.id != w && ItemMatches(w, d.id)) { Rule u; u.in = { d.id }; u.out = { w }; rules.push_back(u); }
                    { Rule u; u.in = { "bucket" }; u.out = { "water_bucket" }; rules.push_back(u); } // scooping water
                    // Plan 5 biome resources: oil seeps (a bucket), desert scrub, jungle trees, swamp digging.
                    { Rule u; u.in = { "bucket" }; u.out = { "crude_oil" }; rules.push_back(u); }
                    { Rule u; u.in = {}; u.out = { "cactus", "rubber", "peat" }; rules.push_back(u); }
                    char buf[1024];
                    if (FILE* f = fopen("resources/data/machine_recipes.txt", "r")) {
                        while (fgets(buf, sizeof buf, f)) {
                            if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
                            const auto c = split(buf, '|');
                            if (c.size() < 3) continue;
                            Rule u; u.machine = c[0] == "furnace" ? "electric_furnace" : c[0]; u.out = ids(c[1]); u.in = ids(c[2]);
                            rules.push_back(u);
                        }
                        fclose(f);
                    }
                    if (FILE* f = fopen("resources/data/animals.txt", "r")) { // drops; products need their tool
                        while (fgets(buf, sizeof buf, f)) {
                            if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
                            const auto c = split(buf, '|');
                            if (c.size() < 13) continue;
                            for (const std::string& d : ids(c[8])) have.insert(d);
                            char prod[64] = {}, tool[64] = {};
                            if (sscanf(c[9].c_str(), "%63s %63s", prod, tool) >= 1 && strcmp(prod, "-") != 0) { Rule u; u.out = { prod }; if (tool[0] && strcmp(tool, "-") != 0) u.in = { tool }; u.machine = "#tool"; rules.push_back(u); }
                            for (const std::string& d : ids(c[7])) used.insert(d);
                        }
                        fclose(f);
                    }
                    if (FILE* f = fopen("resources/data/crops.txt", "r")) {
                        while (fgets(buf, sizeof buf, f)) {
                            if (buf[0] == '#' || buf[0] == '\n' || buf[0] == '\r') continue;
                            const auto c = split(buf, '|');
                            if (c.size() < 3) continue;
                            Rule u; u.in = { c[1] }; u.out = ids(c[2]); u.machine = "#farmland"; rules.push_back(u);
                        }
                        fclose(f);
                    }
                    for (bool grew = true; grew;) {
                        grew = false;
                        for (const Rule& u : rules) {
                            if (!u.machine.empty() && u.machine[0] != '#' && !have.count(u.machine)) continue;
                            if (u.machine == "#farmland" && !have.count("farmland")) continue;
                            bool ok = true;
                            for (const std::string& in : u.in) ok = ok && (u.machine == "#tool" ? have.count(in) > 0 : have.count(in) > 0);
                            if (!ok) continue;
                            for (const std::string& o : u.out) if (have.insert(o).second) grew = true;
                        }
                    }
                    for (const Rule& u : rules) { for (const std::string& in : u.in) used.insert(in); if (!u.machine.empty() && u.machine[0] != '#') used.insert(u.machine); }
                    std::vector<std::string> missing, unused;
                    for (const ItemDef& d : itemDB.items()) {
                        if (!have.count(d.id)) missing.push_back(d.id);
                        const bool consumed = used.count(d.id) || BuildMaterialForItem(d.id.c_str()) >= 0 || Survival_IsFood(d.id) || Survival_IsUsable(d.id) || Suit_IsSuit(d.id) || Suit_IsModule(d.id) ||
                                              ToolOf(d.id).type != ToolType::Hands || KindOf(d.id) == ItemKind::Block || KindOf(d.id) == ItemKind::Tool;
                        if (!consumed) unused.push_back(d.id);
                    }
                    printf("bench: recipes: %zu items, %zu reachable from the world;", itemDB.items().size(), have.size());
                    if (!missing.empty()) { printf(" unreachable:"); for (const auto& m : missing) printf(" %s", m.c_str()); }
                    printf("\nbench: recipes: %zu items nothing uses (end products or loose ends):", unused.size());
                    for (const auto& m : unused) printf(" %s", m.c_str());
                    printf("\nbench: recipes: %s\n", missing.empty() ? "PASS" : "FAIL");
                    benchRecipes = false;
                } else if (benchShotTabs > 0 && (benchShotTabs > 1 || WorldIdle())) {
                    // --shot-tabs: screenshots of the Quest tree (some progress) and Suit tabs.
                    static int frames = 0;
                    if (benchShotTabs == 1) {
                        if (Planetoid* p0 = dynamic_cast<Planetoid*>(world.objects[0].get())) { player.flyMode = false; player.spawnOn(world, *p0, StartSpawnDir(*p0)); } // stands by the containers placed below
                        QuestBook_SetTier(3);
                        QuestBook_AtLeast("explored", 30.0); QuestBook_Count("till"); QuestBook_Count("feed");
                        inventory.add("food_pouch", 1); inventory.add("oxygen_tank_module", 1); inventory.add("grip_boots", 1); inventory.add("field_suit", 1);
                        Suit_Reset(); Suit_Fit("food_pouch");
                        QuestBook_Update();
                        GameUI_ShowTab("Quest tree");
                        benchShotTabs = 2; frames = 0;
                    } else if (benchShotTabs == 2 && ++frames == 40) {
                        TakeScreenshot("bench_questtree.png");
                        GameUI_ShowTab("Suit");
                        benchShotTabs = 3; frames = 0;
                    } else if (benchShotTabs == 3 && ++frames == 40) {
                        TakeScreenshot("bench_suit.png");
                        GameUI_ShowTab("Crafting");
                        benchShotTabs = 4; frames = 0;
                    } else if (benchShotTabs == 4 && ++frames == 40) {
                        TakeScreenshot("bench_crafting.png");
                        // Container screens: a chest and a furnace next to the player, filled, opened.
                        Planetoid* p = player.nearest;
                        static uint32_t plotId = 0; static int chestPart = -1, furnacePart = -1; static int stoneBefore = 0;
                        if (p) {
                            const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                            const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                            if (!p->plots) p->plots = std::make_unique<PlotSet>();
                            Plot* plot = p->plots->plotFor(p->toLocal(Vector3Add(Vector3Add(feet, Vector3Scale(flat, 5.0f)), Vector3Scale(player.up, 1.0f))), 0.0f, true);
                            plotId = plot->id;
                            chestPart = plot->place({ PartType::Chest, { 0, 0, 0 }, 0 });
                            furnacePart = plot->place({ PartType::Furnace, { 3, 0, 0 }, 0 });
                        }
                        benchShotTabs = 5; frames = 0;
                        if (!p || chestPart < 0 || furnacePart < 0) { printf("bench: container: could not place a chest and furnace -> FAIL\n"); benchShotTabs = 0; GameUI_ToggleInventory(); }
                        (void)stoneBefore;
                    } else if (benchShotTabs == 5 && ++frames == 30) {
                        // The layout reached the simulation: fill both, open the chest.
                        Planetoid* p = player.nearest;
                        static uint32_t plotId2 = 0;
                        (void)plotId2;
                        Plot* plot = p && p->plots ? p->plots->plots.back().get() : nullptr;
                        if (plot) {
                            int chest = -1, furnace = -1;
                            for (int i = 0; i < (int)plot->parts.size(); ++i) { if (plot->parts[i].type == PartType::Chest) chest = i; if (plot->parts[i].type == PartType::Furnace) furnace = i; }
                            Factory::Seed(p->id, plot->id, chest, "stone", 40); Factory::Seed(p->id, plot->id, chest, "oak_log", 20); Factory::Seed(p->id, plot->id, chest, "iron_ore", 7);
                            Factory::Seed(p->id, plot->id, furnace, "iron_ore", 5); Factory::Seed(p->id, plot->id, furnace, "oak_log", 8);
                            GameUI_OpenContainer(p->id, plot->id, chest, (int)PartType::Chest);
                        }
                        benchShotTabs = 6; frames = 0;
                    } else if (benchShotTabs == 6 && ++frames == 40) {
                        TakeScreenshot("bench_chest.png");
                        // Take the chest's first slot (stone) through the screen's own command.
                        Planetoid* p = player.nearest;
                        Plot* plot = p->plots->plots.back().get();
                        int chest = -1;
                        for (int i = 0; i < (int)plot->parts.size(); ++i) if (plot->parts[i].type == PartType::Chest) chest = i;
                        const int before = inventory.count("stone");
                        Factory::TakeSlot(p->id, plot->id, chest, 0, false);
                        Factory::WaitTicks(3);
                        Factory::Sync(world);
                        const MachineView v = Factory::View(p->id, plot->id, chest);
                        printf("bench: container: took slot 0 of the chest: stone %d -> %d, chest slots %d (slot 0 now %s x%d): %s\n", before, inventory.count("stone"),
                               (int)v.inSlots.size(), v.inSlots.empty() ? "?" : v.inSlots[0].first.c_str(), v.inSlots.empty() ? 0 : v.inSlots[0].second,
                               inventory.count("stone") == before + 40 && v.inSlots.size() == 16 ? "PASS" : "FAIL");
                        for (int i = 0; i < (int)plot->parts.size(); ++i) if (plot->parts[i].type == PartType::Furnace) GameUI_OpenContainer(p->id, plot->id, i, (int)PartType::Furnace);
                        benchShotTabs = 7; frames = 0;
                    } else if (benchShotTabs == 7 && ++frames == 40) {
                        TakeScreenshot("bench_machine.png");
                        printf("bench: shot-tabs: bench_questtree.png, bench_suit.png, bench_crafting.png, bench_chest.png, bench_machine.png\n");
                        GameUI_ToggleInventory();
                        benchShotTabs = 0;
                    }
                } else if (benchShotBiomes >= 0 && (benchShotBiomes > 0 || WorldIdle())) {
                    // --shot-biomes (Plan 5 B): the nearest deep patch of each biome, seen from 30 m up looking
                    // out across it (bench_biome_<name>.png).
                    static const Biome kShots[] = { Biome::Desert, Biome::Badlands, Biome::Savanna, Biome::Jungle, Biome::Swamp, Biome::Taiga, Biome::Meadow, Biome::Tundra, Biome::Forest, Biome::Alpine };
                    constexpr int kShotCount = (int)(sizeof(kShots) / sizeof(kShots[0]));
                    static Vector3 dirs[kShotCount];
                    static double t0 = 0.0; static bool waiting = false;
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { benchShotBiomes = -1; }
                    else if (benchShotBiomes == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        std::mt19937 rng(11);
                        std::uniform_real_distribution<float> u(-1.0f, 1.0f);
                        float best[kShotCount];
                        for (int k = 0; k < kShotCount; ++k) { best[k] = -2.0f; dirs[k] = home; }
                        for (int i = 0; i < 60000; ++i) {
                            Vector3 d = { u(rng), u(rng), u(rng) };
                            if (Vector3LengthSqr(d) > 1.0f || Vector3LengthSqr(d) < 1e-4f) continue;
                            d = Vector3Normalize(d);
                            const float close = Vector3DotProduct(d, home);
                            const float e = pf->elevation(d, 32.0f);
                            if (e < 6.0f) continue;
                            const Biome b = pf->biomeAt(d, e);
                            for (int k = 0; k < kShotCount; ++k) {
                                if (b != kShots[k] || close <= best[k]) continue;
                                // Deep in it: the same biome 200 m out on four sides.
                                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(d, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, d);
                                bool deep = true;
                                for (int s = 0; s < 4 && deep; ++s) {
                                    const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(s < 2 ? t1 : t2, (s & 1 ? 200.0f : -200.0f) / pf->radius)));
                                    deep = pf->biomeAt(q, pf->elevation(q, 32.0f)) == b;
                                }
                                if (deep) { best[k] = close; dirs[k] = d; }
                            }
                        }
                        for (int k = 0; k < kShotCount; ++k) printf("bench: biomes: nearest %s %.0f km away%s\n", BiomeName(kShots[k]), acosf(std::clamp(best[k], -1.0f, 1.0f)) * pf->radius / 1000.0f, best[k] < -1.5f ? " (none found)" : "");
                        benchShotBiomes = 1; waiting = false;
                    } else {
                        const int k = benchShotBiomes - 1;
                        if (k == kShotCount) {
                            // One more: the nearest floating island (Plan 5 C), from 70 m off at its height.
                            static Vector3 isl{}; static int islPhase = 0; static double islT0 = 0.0;
                            if (islPhase == 0) {
                                std::vector<Vector3> found;
                                const Vector3 home = pf->findSpawnDir(p->seed);
                                pf->floatingIslandsNear(home, 40000.0f, found);
                                float best = 1e30f;
                                for (const Vector3& c : found) { const float dd = Vector3Distance(Vector3Normalize(c), home); if (dd < best) { best = dd; isl = c; } }
                                if (found.empty()) { printf("bench: biomes: no floating island within 40 km\n"); ++benchShotBiomes; }
                                else { player.spawnOn(world, *p, Vector3Normalize(isl)); player.flyMode = true; islPhase = 1; islT0 = GetTime(); }
                            } else {
                                const Vector3 up = Vector3Normalize(isl);
                                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
                                const Vector3 eye = Vector3Add(Vector3Add(isl, Vector3Scale(t1, 330.0f)), Vector3Scale(up, -15.0f)); // from the side, a little below the top: the whole island and its underside
                                world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(isl); world.camera.up = p->dirToWorld(up);
                                player.position = world.renderOrigin + world.camera.position;
                                if (GetTime() - islT0 > 6.0 && WorldIdle()) { TakeScreenshot("bench_biome_island.png"); printf("bench: biomes: bench_biome_island.png\n"); ++benchShotBiomes; islPhase = 0; }
                            }
                        } else
                        if (k >= kShotCount) { printf("bench: biomes: shots done\n"); benchShotBiomes = -1; }
                        else if (!waiting) {
                            player.spawnOn(world, *p, dirs[k]);
                            player.flyMode = true;
                            waiting = true; t0 = GetTime();
                        } else {
                            const Vector3 up = dirs[k];
                            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
                            const float e = std::max(pf->elevation(up), 0.0f);
                            const Vector3 eye = Vector3Scale(up, pf->radius + e + 30.0f);
                            const Vector3 look = Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 160.0f / pf->radius)));
                            const Vector3 at = Vector3Scale(look, pf->radius + std::max(pf->elevation(look), 0.0f));
                            world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(at); world.camera.up = p->dirToWorld(up);
                            player.position = world.renderOrigin + world.camera.position;
                            sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(Vector3CrossProduct(up, t1), 0.7f))));
                            if (GetTime() - t0 > 4.0 && WorldIdle()) {
                                if (GetTime() - t0 > 6.0) {
                                    TakeScreenshot(TextFormat("bench_biome_%s.png", BiomeName(kShots[k])));
                                    printf("bench: biomes: bench_biome_%s.png\n", BiomeName(kShots[k]));
                                    ++benchShotBiomes; waiting = false;
                                }
                            } else if (GetTime() - t0 > 60.0) { printf("bench: biomes: %s never loaded\n", BiomeName(kShots[k])); ++benchShotBiomes; waiting = false; }
                        }
                    }
                } else if (benchOil && (oilPhase > 0 || WorldIdle())) {
                    // --test-oil (Plan 5 B): the nearest oil seep to the start wells up crude oil (fluid sim,
                    // oil cells), a bucket takes some, it survives save + load; bench_oil.png.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    static Vector3 seep{}; static float before = 0.0f;
                    const double t = GetTime() - oilT0;
                    if (!pf) { printf("bench: oil: needs the planet -> FAIL\n"); benchOil = false; }
                    else if (oilPhase == 0) {
                        const Vector3 home = Vector3Scale(pf->findSpawnDir(p->seed), pf->radius);
                        std::vector<Vector3> found;
                        float best = 1e30f;
                        for (float range : { 4000.0f, 12000.0f, 30000.0f }) {
                            Ambience_SeepsNear(*p, home, range, found);
                            for (const Vector3& s : found) { const float d = Vector3Distance(s, home); if (d < best) { best = d; seep = s; } }
                            if (best < 1e29f) break;
                        }
                        if (best > 1e29f) { printf("bench: oil: no seep within 30 km -> FAIL\n"); benchOil = false; }
                        else {
                            const Vector3 up = Vector3Normalize(seep);
                            printf("bench: oil: nearest seep %.1f km from the start (%s)\n", best / 1000.0f, BiomeName(pf->biomeAt(up, pf->elevation(up))));
                            player.spawnOn(world, *p, up);
                            player.flyMode = true;
                            oilPhase = 1; oilT0 = GetTime();
                        }
                    } else {
                        const Vector3 up = Vector3Normalize(seep);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
                        const Vector3 eye = Vector3Add(Vector3Add(seep, Vector3Scale(t1, 7.0f)), Vector3Scale(up, 4.0f));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(seep); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 0.6f))));
                        const Vector3 pool = Vector3Add(seep, Vector3Scale(up, 0.3f));
                        const float oil = Fluid::FluidAt(*p, pool, 4.0f, Fluid::KindByName("crude_oil"));
                        if (oilPhase == 1 && t > 25.0) {

                            before = oil;
                            const float took = Fluid::Take(*p, pool, 1.0f, 2.5f, Fluid::KindByName("crude_oil"));
            { const Fluid::Stats fs = Fluid::GetStats(); printf("bench: oil: %.1f m^3 in the pool after 25 s (water near: %.1f m^3); a bucket takes %.1f; sources %d, sites %d, bricks %d, other fluid %.1f m^3, kinds %d\n", oil, Fluid::WaterAt(*p, pool, 2.0f), took, Fluid::SourceCount(), fs.sites, fs.bricks, fs.otherVolume, Fluid::KindCount()); }
                            oilPhase = 2; oilT0 = GetTime();
                        } else if (oilPhase == 2 && t > 3.0) {
                            const bool saved = SaveWorld(world), loaded = saved && LoadWorldSave(world);
                            printf("bench: oil: after the bucket %.1f m^3; saved %d, loaded %d\n", oil, (int)saved, (int)loaded);
                            oilPhase = 3; oilT0 = GetTime(); Fluid::SetSpeed(8); // fast-forward: the seep fills its basin to its head and rests
                        } else if (oilPhase == 3 && t > 6.0 && Fluid::Asleep()) {
                            Fluid::SetSpeed(1);
                            TakeScreenshot("bench_oil.png");
                            const bool ok = before >= 1.0f && oil >= 3.0f && Fluid::SourceCount() >= 1;
                            printf("bench: oil: after the load the seep (%d source) filled its basin and rested at %.1f m^3 (equilibrium) -> %s\n", Fluid::SourceCount(), oil, ok ? "PASS" : "FAIL");
                            benchOil = false;
                        } else if (oilPhase == 3 && t > 90.0) { Fluid::SetSpeed(1); printf("bench: oil: no equilibrium after the load (%.1f m^3) -> FAIL\n", oil); benchOil = false; }
                    }
                } else if (benchSpring && (springPhase > 0 || WorldIdle())) {
                    // --test-spring (Plan 5): the nearest water spring to the start (Fluid::AddSource, head above the
                    // ground) runs off as a stream: after a while most of its water is well away from it; bench_spring.png.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    static Vector3 spring{};
                    const double t = GetTime() - springT0;
                    if (!pf) { printf("bench: spring: needs the planet -> FAIL\n"); benchSpring = false; }
                    else if (springPhase == 0) {
                        const Vector3 home = Vector3Scale(pf->findSpawnDir(p->seed), pf->radius);
                        std::vector<Vector3> found;
                        float best = 1e30f;
                        for (float range : { 3000.0f, 10000.0f }) {
                            Ambience_SpringsNear(*p, home, range, found);
                            for (const Vector3& s : found) { const float d = Vector3Distance(s, home); if (d < best) { best = d; spring = s; } }
                            if (best < 1e29f) break;
                        }
                        if (best > 1e29f) { printf("bench: spring: none within 10 km -> FAIL\n"); benchSpring = false; }
                        else {
                            const Vector3 up = Vector3Normalize(spring);
                            printf("bench: spring: nearest %.1f km from the start (%s, %.0f m up)\n", best / 1000.0f, BiomeName(pf->biomeAt(up, pf->elevation(up))), pf->elevation(up));
                            player.spawnOn(world, *p, up);
                            player.flyMode = true;
                            springPhase = 1; springT0 = GetTime();
                        }
                    } else {
                        const Vector3 up = Vector3Normalize(spring);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
                        const Vector3 eye = Vector3Add(Vector3Add(spring, Vector3Scale(t1, 18.0f)), Vector3Scale(up, 26.0f)); // well above: the ground can stand above the height map
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(spring); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 0.6f))));
                        if (springPhase == 1 && t > 5.0) { Fluid::SetSpeed(6); springPhase = 2; springT0 = GetTime(); }
                        else if (springPhase == 2 && t > 30.0) {
                            Fluid::SetSpeed(1);
                            TakeScreenshot("bench_spring.png");
                            const float total = (float)Fluid::GetStats().volume, nearSrc = Fluid::FluidAt(*p, spring, 3.0f, 0);
                            const bool ok = Fluid::SourceCount() >= 1 && total > 5.0f && nearSrc < 0.6f * total;
                            printf("bench: spring: %.1f m^3 of moving water, %.1f of it within 3 m of the spring (a stream: most runs off) -> %s\n", total, nearSrc, ok ? "PASS" : "FAIL");
                            printf("bench: spring: bench_spring.png\n");
                            benchSpring = false;
                        }
                    }
                } else if (benchShotAnimals && (animalShotPhase > 0 || WorldIdle())) {
                    // --shot-animals (Plan 5 A): one of each kind in a row, side-on, 9 m in front (bench_animals_models.png).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    static double t0 = 0.0; static Vector3 eye{}, at{};
                    if (!pf) benchShotAnimals = false;
                    else if (animalShotPhase == 0) {
                        Animals_Reset(); Animals_SetThreats(false);
                        player.flyMode = true;
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
                        auto ground = [&](Vector3 q) { const Vector3 d = Vector3Normalize(q); return Vector3Scale(d, pf->radius + pf->elevation(d)); };
                        const char* kinds[] = { "woolback", "plains_runner", "clucker", "hornhide", "thornback", "dusk_stalker" };
                        for (int k = 0; k < 6; ++k) Animals_Spawn(world, *p, kinds[k], ground(Vector3Add(me, Vector3Add(Vector3Scale(t1, -6.5f + 2.6f * k), Vector3Scale(t2, 9.0f)))), 1, false, 0.0f);
                        at = ground(Vector3Add(me, Vector3Scale(t2, 9.0f)));
                        eye = Vector3Add(ground(me), Vector3Scale(up, 2.2f));
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(Vector3Subtract(t1, t2), 0.5f)))); // daylight, from behind the camera
                        animalShotPhase = 1; t0 = GetTime();
                    } else {
                        const Vector3 up = Vector3Normalize(at);
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(Vector3Add(at, Vector3Scale(up, 0.6f))); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - t0 > 4.0) {
                            TakeScreenshot("bench_animals_models.png");
                            printf("bench: animals: bench_animals_models.png (%d animals)\n", Animals_Count());
                            Animals_SetThreats(true);
                            benchShotAnimals = false;
                        }
                    }
                } else if (benchShotMachines && (machineShotPhase > 0 || WorldIdle())) {
                    // --shot-machines (Plan 5): every machine type in two rows on a plot ahead (bench_machines.png).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    static double t0 = 0.0; static Vector3 eye{}, at{};
                    if (!p) benchShotMachines = false;
                    else if (machineShotPhase == 0) {
                        player.flyMode = true;
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 fwd = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(fwd, Vector3Scale(player.up, Vector3DotProduct(fwd, player.up))));
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(p->toLocal(Vector3Add(Vector3Add(feet, Vector3Scale(flat, 10.0f)), Vector3Scale(player.up, 1.0f))), 0.0f, true);
                        const PartType types[] = { PartType::Chest, PartType::Generator, PartType::Furnace, PartType::Macerator, PartType::Battery, PartType::Campfire, PartType::Kiln,
                                                   PartType::CokeOven, PartType::BlastFurnace, PartType::Boiler, PartType::ForgeHammer, PartType::AlloySmelter, PartType::SteamMacerator,
                                                   PartType::OreWasher, PartType::Centrifuge, PartType::WaterPump, PartType::MvGenerator, PartType::Electrolyzer, PartType::ChemicalReactor,
                                                   PartType::Miner, PartType::HvGenerator, PartType::ArcFurnace, PartType::Refinery, PartType::PrecisionAssembler, PartType::FusionReactor, PartType::MassDriver };
                        int k = 0, placed = 0;
                        for (PartType t : types) { placed += plot->place({ t, { (k % 13) * 4 - 24, 2, (k / 13) * 6 }, 0 }) >= 0 ? 1 : 0; ++k; }
                        printf("bench: machines: placed %d of %d types\n", placed, k);
                        at = p->toLocal(Vector3Add(Vector3Add(feet, Vector3Scale(flat, 12.0f)), Vector3Scale(player.up, 1.5f)));
                        eye = p->toLocal(Vector3Add(feet, Vector3Scale(player.up, 6.0f)));
                        sunOverride = true; sunArg = Vector3Normalize(Vector3Add(player.up, Vector3Scale(Vector3Negate(flat), 0.7f))); // daylight from behind
                        machineShotPhase = 1; t0 = GetTime();
                    } else {
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(at); world.camera.up = p->dirToWorld(Vector3Normalize(eye));
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - t0 > 3.0) {
                            TakeScreenshot("bench_machines.png");
                            printf("bench: machines: bench_machines.png\n");
                            benchShotMachines = false;
                        }
                    }
                } else if (benchBiomes && WorldIdle()) {
                    // --test-biomes (Plan 5 B): biome shares over the land, the biomes within 5 km of the
                    // start (regional: several in reach), and the start's own (a temperate one).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    benchBiomes = false;
                    if (!pf) printf("bench: biomes: needs the planet -> FAIL\n");
                    else {
                        int count[(int)Biome::Count] = {}, land = 0;
                        std::mt19937 rng(7);
                        std::uniform_real_distribution<float> u(-1.0f, 1.0f);
                        for (int i = 0; i < 40000; ++i) {
                            Vector3 d = { u(rng), u(rng), u(rng) };
                            if (Vector3LengthSqr(d) > 1.0f || Vector3LengthSqr(d) < 1e-4f) { --i; continue; }
                            d = Vector3Normalize(d);
                            const float e = pf->elevation(d, 64.0f);
                            if (e < 0.0f) continue;
                            ++land; ++count[(int)pf->biomeAt(d, e)];
                        }
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Biome start = pf->biomeAt(home, pf->elevation(home));
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        int nearCount[(int)Biome::Count] = {};
                        for (float a = -5000.0f; a <= 5000.0f; a += 100.0f)
                            for (float b = -5000.0f; b <= 5000.0f; b += 100.0f) {
                                if (a * a + b * b > 25e6f) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, a / pf->radius), Vector3Scale(t2, b / pf->radius))));
                                const float e = pf->elevation(d, 64.0f);
                                ++nearCount[(int)pf->biomeAt(d, e)];
                            }
                        {   // Plan 5 C: floating islands near the start
                            std::vector<Vector3> isl;
                            pf->floatingIslandsNear(home, 25000.0f, isl);
                            float bestI = 1e30f;
                            for (const Vector3& c : isl) bestI = std::min(bestI, Vector3Distance(Vector3Scale(Vector3Normalize(c), pf->radius), Vector3Scale(home, pf->radius)));
                            printf("bench: biomes: %d floating islands within 25 km, the nearest %.1f km away\n", (int)isl.size(), bestI < 1e29f ? bestI / 1000.0f : -1.0f);
                        }
                        int mouths = 0, probes = 0; // Plan 5 C: open ground 2 m under the surface (cave mouths), on a 25 m grid within 2 km
                        for (float a = -2000.0f; a <= 2000.0f; a += 25.0f)
                            for (float b = -2000.0f; b <= 2000.0f; b += 25.0f) {
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, a / pf->radius), Vector3Scale(t2, b / pf->radius))));
                                const float e = pf->elevation(d);
                                if (e < 2.0f) continue;
                                ++probes;
                                if (pf->baseDensity(Vector3Scale(d, pf->radius + e - 2.0f)) < ISO_LEVEL && pf->baseDensity(Vector3Scale(d, pf->radius + e - 6.0f)) < ISO_LEVEL) ++mouths;
                            }
                        printf("bench: biomes: cave openings: %d of %d land columns within 2 km open 2-6 m down\n", mouths, probes);
                        int kinds = 0, nearKinds = 0;
                        printf("bench: biomes: land share:");
                        for (int b = 1; b < (int)Biome::Count; ++b) {
                            printf(" %s %.1f%%", BiomeName((Biome)b), 100.0f * count[b] / std::max(land, 1));
                            if (b != (int)Biome::Polar && count[b] > 0) ++kinds;
                            if (b != (int)Biome::Ocean && nearCount[b] > 20) ++nearKinds;
                        }
                        printf("\nbench: biomes: within 5 km of the start:");
                        for (int b = 0; b < (int)Biome::Count; ++b) if (nearCount[b]) printf(" %s %d", BiomeName((Biome)b), nearCount[b]);
                        const bool ok = kinds == 10 && nearKinds >= 3 && (start == Biome::Forest || start == Biome::Meadow);
                        printf("\nbench: biomes: %d of 10 land biomes on the planet, %d near the start, start in %s -> %s\n", kinds, nearKinds, BiomeName(start), ok ? "PASS" : "FAIL");
                    }
                } else if (benchParticles && (particlePhase > 0 || WorldIdle())) {
                    // --shot-particles: every particle kind in front of a standing camera: a smoke column,
                    // a flame with sparks, a dig puff with clods, water splashes and fireflies (bench_particles.png).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - particleT0;
                    if (particlePhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        particlePhase = 1; particleT0 = GetTime();
                    } else if (particlePhase == 1 && t > 1.5) {
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        const auto* pfp = dynamic_cast<const PlanetField*>(p->field.get());
                        const Vector3 ground = pfp ? Vector3Scale(up, pfp->radius + pfp->elevation(up)) : me;
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 0.8f))));
                        player.flyMode = true;
                        particleEye = Vector3Add(ground, Vector3Scale(up, 1.7f));
                        particleAt = Vector3Add(Vector3Add(ground, Vector3Scale(t1, 6.0f)), Vector3Scale(up, 1.0f)); // above the grass
                        particlePhase = 2; particleT0 = GetTime();
                    } else if (particlePhase == 2) {
                        const Vector3 up = Vector3Normalize(particleEye);
                        world.camera.position = p->toWorld(particleEye);
                        world.camera.target = p->toWorld(Vector3Add(particleAt, Vector3Scale(up, 0.6f)));
                        world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        const Vector3 fwd = Vector3Normalize(Vector3Subtract(particleAt, particleEye));
                        const Vector3 side = Vector3Normalize(Vector3CrossProduct(fwd, up));
                        auto at = [&](float s, float lift) { const Vector3 q = Vector3Add(particleAt, Vector3Scale(side, s)); const Vector3 d = Vector3Normalize(q);
                                                             const auto* pf = dynamic_cast<const PlanetField*>(p->field.get()); return Vector3Scale(d, (pf ? pf->radius + pf->elevation(d) : Vector3Length(q)) + lift + 1.0f); };
                        static double lastEmit = 0.0;
                        if (GetTime() - lastEmit > 0.05 && t < 2.6) {
                            lastEmit = GetTime();
                            ParticleBurst b;
                            b = {}; b.kind = ParticleKind::Smoke; b.color = { 190, 188, 185, 255 }; b.count = 2; b.size = 0.5f; b.life = 4.0f; b.speed = 0.9f; b.spread = 0.25f;
                            Particles_Emit(*p, at(-4.0f, 0.3f), b);
                            b = {}; b.kind = ParticleKind::Fire; b.color = { 255, 120, 30, 255 }; b.count = 4; b.size = 0.45f; b.life = 0.6f; b.speed = 1.0f; b.spread = 0.2f;
                            Particles_Emit(*p, at(-1.3f, 0.1f), b);
                            b = {}; b.kind = ParticleKind::Spark; b.color = { 255, 170, 60, 255 }; b.count = 1; b.size = 0.05f; b.life = 1.4f; b.speed = 2.5f; b.spread = 0.3f;
                            Particles_Emit(*p, at(-1.3f, 0.3f), b);
                            b = {}; b.kind = ParticleKind::Splash; b.color = { 175, 195, 225, 200 }; b.count = 3; b.size = 0.06f; b.life = 0.4f; b.speed = 1.4f; b.spread = 0.45f;
                            Particles_Emit(*p, at(4.0f, 0.0f), b);
                            if (GetRandomValue(0, 4) == 0) {
                                b = {}; b.kind = ParticleKind::Dust; b.color = { 122, 92, 62, 200 }; b.count = 10; b.size = 0.35f; b.life = 1.1f; b.speed = 1.2f; b.spread = 0.9f;
                                Particles_Emit(*p, at(1.4f, 0.0f), b);
                                b = {}; b.kind = ParticleKind::Debris; b.color = { 122, 92, 62, 255 }; b.count = 7; b.size = 0.07f; b.life = 1.6f; b.speed = 3.2f; b.spread = 0.5f;
                                Particles_Emit(*p, at(1.4f, 0.0f), b);
                                b = {}; b.kind = ParticleKind::Firefly; b.color = { 200, 255, 110, 255 }; b.count = 3; b.size = 0.1f; b.life = 6.0f; b.speed = 0.3f; b.spread = 1.0f;
                                Particles_Emit(*p, at(GetRandomValue(-50, 50) / 10.0f, 1.2f), b);
                            }
                        }
                        if (t > 2.4) {
                            TakeScreenshot("bench_particles.png");
                            printf("bench: particles: %d live -> bench_particles.png\n", Particles_Count());
                            benchParticles = false;
                        }
                    }
                } else if (benchShadow && (shadowPhase > 0 || WorldIdle())) {
                    // --shot-shadow: a tree's shadow on the ground from 3 m away (bench_shadow.png), and a
                    // close look straight at flat sunlit ground (bench_shadow_ground.png: no acne pattern).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - shadowT0;
                    if (shadowPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        shadowPhase = 1; shadowT0 = GetTime();
                    } else if (shadowPhase == 1 && t > 1.5) {
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        std::vector<Trunk> near;
                        if (p->trees) p->trees->trunksNear(me, 120.0f, near);
                        float best = 1e9f; Trunk tr{};
                        for (const Trunk& k : near) { const float d = Vector3Distance(k.base, me); if (d < best && k.height > 3.0f) { best = d; tr = k; } }
                        // Morning sun from one side, low enough for long shadows.
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 }));
                        const Vector3 sunL = Vector3Normalize(Vector3Add(up, Vector3Scale(t1, 1.2f)));
                        sunOverride = true; sunArg = p->dirToWorld(sunL);
                        player.flyMode = true;
                        if (best > 1e8f) { printf("bench: shadow: no tree nearby\n"); benchShadow = false; }
                        else {
                            const Vector3 away = Vector3Normalize(Vector3Subtract(Vector3Negate(t1), Vector3Scale(up, Vector3DotProduct(Vector3Negate(t1), up))));
                            shadowAt = Vector3Add(tr.base, Vector3Scale(away, 6.0f)); // past the crown: the trunk's and crown's shadow edges
                            const Vector3 side = Vector3CrossProduct(up, away);
                            shadowEye = Vector3Add(Vector3Add(shadowAt, Vector3Scale(side, 4.0f)), Vector3Scale(up, 2.2f));
                            shadowPhase = 2; shadowT0 = GetTime();
                        }
                    } else if (shadowPhase >= 2) {
                        const Vector3 up = Vector3Normalize(shadowAt);
                        const Vector3 eye = shadowPhase == 2 ? shadowEye : Vector3Add(Vector3Add(shadowEye, Vector3Scale(Vector3Normalize(Vector3Subtract(shadowEye, shadowAt)), 6.0f)), Vector3Scale(up, 0.5f));
                        const Vector3 at = shadowPhase == 2 ? shadowAt : Vector3Add(shadowEye, Vector3Scale(Vector3Normalize(Vector3Subtract(shadowEye, shadowAt)), 9.0f)); // bare ground ahead, looked down at
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(Vector3Subtract(at, Vector3Scale(up, shadowPhase == 2 ? 0.0f : 3.0f))); world.camera.up = p->dirToWorld(up);
                        if (shadowPhase == 4) {
                            // Standing height (1.7 m over the ground), looking out across the field, a little down.
                            const Vector3 out = Vector3Normalize(Vector3Subtract(shadowEye, shadowAt));
                            const Vector3 flat = Vector3Normalize(Vector3Subtract(out, Vector3Scale(up, Vector3DotProduct(out, up))));
                            const auto* pf4 = dynamic_cast<const PlanetField*>(p->field.get());
                            const Vector3 ground = pf4 ? Vector3Scale(up, pf4->radius + pf4->elevation(up)) : shadowAt;
                            const Vector3 e4 = Vector3Add(ground, Vector3Scale(up, 1.7f));
                            world.camera.position = p->toWorld(e4);
                            world.camera.target = p->toWorld(Vector3Add(e4, Vector3Subtract(Vector3Scale(flat, 20.0f), Vector3Scale(up, 2.5f))));
                        }
                        player.position = world.renderOrigin + world.camera.position;
                        if (t > 3.0) {
                            TakeScreenshot(shadowPhase == 2 ? "bench_shadow.png" : shadowPhase == 3 ? "bench_shadow_ground.png" : "bench_grass.png");
                            if (shadowPhase == 4) { printf("bench: shadow: bench_shadow.png, bench_shadow_ground.png, bench_grass.png\n"); benchShadow = false; }
                            else { ++shadowPhase; shadowT0 = GetTime(); }
                        }
                    }
                } else if (benchFarForest && (farForestPhase > 0 || WorldIdle())) {
                    // --shot-far-forest: from 30 m over the start, look toward the most wooded direction
                    // 1.5 km out (bench_far_forest.png): real trees out to the tree distance, canopy tint past it.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { printf("bench: far forest: needs the planet\n"); benchFarForest = false; }
                    else if (farForestPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = -1.0f; Vector3 bestDir = t1;
                        for (int k = 0; k < 72; ++k) {
                            const float a = k * 2.0f * PI / 72.0f;
                            const Vector3 flat = Vector3Add(Vector3Scale(t1, cosf(a)), Vector3Scale(t2, sinf(a)));
                            float cover = 0.0f;
                            for (float r = 1000.0f; r <= 2200.0f; r += 200.0f) {
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Scale(flat, r / pf->radius)));
                                cover += TreeCanopy(d, pf->radius, pf->elevation(d), pf->foliageSeed, pf);
                            }
                            if (cover > best) { best = cover; bestDir = flat; }
                        }
                        const Vector3 ground = Vector3Scale(home, pf->radius + pf->elevation(home));
                        farForestEye = Vector3Add(ground, Vector3Scale(home, 30.0f));
                        farForestAt = Vector3Add(farForestEye, Vector3Subtract(Vector3Scale(bestDir, 1500.0f), Vector3Scale(home, 60.0f)));
                        printf("bench: far forest: most wooded direction has canopy %.2f (of 7) between 1 and 2.2 km\n", best);
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(home, Vector3Scale(bestDir, -0.6f)))); // sun behind the camera
                        farForestPhase = 1; farForestT0 = GetTime();
                    } else {
                        world.camera.position = p->toWorld(farForestEye); world.camera.target = p->toWorld(farForestAt);
                        world.camera.up = p->dirToWorld(Vector3Normalize(farForestEye));
                        player.position = world.renderOrigin + world.camera.position;
                        player.flyMode = true;
                        if (GetTime() - farForestT0 > 10.0) {
                            TakeScreenshot("bench_far_forest.png");
                            printf("bench: far forest: bench_far_forest.png (trees drawn to %.0f m)\n", kTreeDrawDistance);
                            benchFarForest = false;
                        }
                    }
                } else if (benchTerrain && (terrainPhase > 0 || WorldIdle())) {
                    // --shot-terrain (plan item 6): the nearest strongly rugged land (the 3D band at full
                    // strength), seen from the air (bench_terrain.png) and from the ground (bench_terrain_low.png).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - terrainT0;
                    if (!pf) { printf("bench: terrain: needs the planet\n"); benchTerrain = false; }
                    else if (terrainPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f; int rugged = 0, samples = 0;
                        for (float e = -5000.0f; e <= 5000.0f; e += 50.0f)
                            for (float n = -5000.0f; n <= 5000.0f; n += 50.0f) {
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const float s = pf->radius + pf->elevation(d);
                                float amp = 0.0f;
                                for (float h : { -8.0f, 0.0f, 8.0f }) amp = std::max(amp, fabsf(pf->terrain3D(Vector3Scale(d, s + h), s + h, s, 1.0f)));
                                ++samples;
                                if (amp < 6.0f) continue;
                                ++rugged;
                                const float dist = sqrtf(e * e + n * n);
                                if (dist < best) { best = dist; terrainAt = Vector3Scale(d, s); terrainOut = t1; }
                            }
                        printf("bench: terrain: %.1f%% of 100 km^2 around the start is rugged (3D band over 6 m); nearest %.0f m away\n", 100.0f * rugged / std::max(1, samples), best);
                        if (best > 1e29f) { printf("bench: terrain: no rugged land -> FAIL\n"); benchTerrain = false; }
                        else {
                            player.flyMode = true; sunOverride = true;
                            const Vector3 up = Vector3Normalize(terrainAt);
                            sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(terrainOut, 0.6f))));
                            terrainPhase = 1; terrainT0 = GetTime();
                        }
                    } else {
                        const Vector3 up = Vector3Normalize(terrainAt);
                        const Vector3 side = Vector3CrossProduct(up, terrainOut);
                        const bool low = terrainPhase == 2;
                        Vector3 eye = Vector3Add(terrainAt, Vector3Add(Vector3Scale(up, low ? 2.0f : 45.0f), Vector3Add(Vector3Scale(terrainOut, low ? -30.0f : -60.0f), Vector3Scale(side, low ? 8.0f : 20.0f))));
                        if (low) { // down from 40 m to the first ground, then 2 m over it
                            eye = Vector3Add(eye, Vector3Scale(up, 38.0f));
                            for (int k = 0; k < 80 && p->densityAt(eye) < ISO_LEVEL; ++k) eye = Vector3Subtract(eye, Vector3Scale(up, 0.5f));
                            eye = Vector3Add(eye, Vector3Scale(up, 2.0f));
                        }
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(low ? Vector3Add(eye, Vector3Add(Vector3Scale(terrainOut, 50.0f), Vector3Scale(up, -4.0f))) : terrainAt); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        if (t > 4.0 && WorldIdle()) { // the chunks there built (a 4 s timer alone shot a half-loaded world)
                            TakeScreenshot(low ? "bench_terrain_low.png" : "bench_terrain.png");
                            if (low) { printf("bench: terrain: bench_terrain.png, bench_terrain_low.png\n"); benchTerrain = false; }
                            else { terrainPhase = 2; terrainT0 = GetTime(); }
                        }
                    }
                } else if (benchPump && (pumpPhase > 0 || WorldIdle())) {
                    // --test-pump (Plan 4, F6): a powered water pump at a river's edge sends water through pipes
                    // to an outlet over a bowl dug 30 m inland; the bowl fills; cutting the cable stops it.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - pumpT0;
                    if (!pf) { printf("bench: pump: needs the planet\n"); benchPump = false; }
                    else if (pumpPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f, level = 0.0f, shoreS = 0.0f; Vector3 wdir{}, out{};
                        for (float e = -3000.0f; e <= 3000.0f; e += 16.0f)
                            for (float n = -3000.0f; n <= 3000.0f; n += 16.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if (!ws.kind || ws.kind == 3 || ws.level <= 2.0f || pf->elevation(d) > ws.level - 1.5f) continue;
                                const Vector3 outs[4] = { t1, Vector3Negate(t1), t2, Vector3Negate(t2) };
                                for (const Vector3& o : outs)
                                    for (float s = 2.0f; s <= 20.0f; s += 1.0f) {
                                        const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(o, s / pf->radius)));
                                        if (pf->waterSourceAt(q).kind) continue;
                                        const float eq = pf->elevation(q);
                                        if (eq > ws.level + 0.2f && eq < ws.level + 6.0f) { best = dist; wdir = d; out = o; level = ws.level; shoreS = s; }
                                        break;
                                    }
                            }
                        if (best > 1e29f) { printf("bench: pump: no water with a shore -> FAIL\n"); benchPump = false; }
                        else {
                            auto along = [&](float s) { return Vector3Normalize(Vector3Add(wdir, Vector3Scale(out, s / pf->radius))); };
                            const Vector3 intakeAt = Vector3Scale(along(shoreS - 2.5f), pf->radius + level + 0.15f); // the pump's floor, just over the water
                            const Vector3 bd = along(shoreS + 28.0f);
                            const Vector3 b = Vector3Scale(bd, pf->radius + pf->elevation(bd));
                            const TerrainEdit bowl{ Vector3Subtract(b, Vector3Scale(bd, 0.8f)), 2.6f, -1.0f, true };
                            ApplyTerrainEdit(*p, bowl);
                            pumpBasin = Vector3Subtract(b, Vector3Scale(bd, 2.4f));
                            pumpRise = pf->elevation(bd) - level;
                            if (!p->plots) p->plots = std::make_unique<PlotSet>();
                            Plot* plot = p->plots->plotFor(intakeAt, 0.0f, true);
                            pumpPlot = plot->id;
                            const Int3 pc = plot->cellAt(plot->toPlot(intakeAt));
                            // The pipe runs along the plot axis closest to inland.
                            const Vector3 o = Vector3Subtract(plot->toPlot(Vector3Add(intakeAt, out)), plot->toPlot(intakeAt));
                            const bool alongX = fabsf(o.x) >= fabsf(o.z);
                            const int dir = (alongX ? o.x : o.z) > 0.0f ? 1 : -1;
                            plot->place({ PartType::WaterPump, pc, 0 });
                            // Power: a generator behind the pump, joined by a cable cell.
                            auto off = [&](Int3 c, int a) { return alongX ? Int3{ c.x + a, c.y, c.z } : Int3{ c.x, c.y, c.z + a }; };
                            const int gen = plot->place({ PartType::Generator, off(pc, dir > 0 ? -3 : 3), 0 });
                            pumpCable = plot->place({ PartType::Cable, off(pc, dir > 0 ? -1 : 2), 0 });
                            // Pipes: out of the pump's top on the inland side, along the axis, across, then down to
                            // the outlet 1 m over the bowl.
                            const Int3 tc = plot->cellAt(plot->toPlot(Vector3Add(b, Vector3Scale(bd, 1.0f))));
                            Int3 c = off(Int3{ pc.x, pc.y + 1, pc.z }, dir > 0 ? 2 : -1);
                            int pipes = 0;
                            auto pipe = [&](Int3 at) { if (plot->place({ PartType::Pipe, at, 0 }) >= 0) ++pipes; };
                            pipe(c);
                            const int ty = std::max(tc.y + 1, c.y);
                            while (c.y < ty) { ++c.y; pipe(c); }
                            while (c.x != tc.x) { c.x += c.x < tc.x ? 1 : -1; pipe(c); }
                            while (c.z != tc.z) { c.z += c.z < tc.z ? 1 : -1; pipe(c); }
                            while (c.y > tc.y + 1) { --c.y; pipe(c); }
                            const int outlet = plot->place({ PartType::Outlet, { c.x, c.y - 1, c.z }, 0 });
                            Factory::Seed(p->id, plot->id, gen, "oak_log", 64);
                            printf("bench: pump: river %.0f m from the start; pump at its edge, %d pipes, outlet %s over a bowl 30 m inland and %.1f m higher; generator %d, cable %d\n",
                                   best, pipes, outlet >= 0 ? "placed" : "BLOCKED", pumpRise, gen, pumpCable);
                            player.flyMode = true; sunOverride = true;
                            sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(bd, Vector3Scale(out, -0.5f))));
                            pumpPhase = 1; pumpT0 = GetTime();
                        }
                    } else if (pumpPhase >= 1) {
                        const Vector3 up = Vector3Normalize(pumpBasin);
                        const Vector3 eye = Vector3Add(pumpBasin, Vector3Add(Vector3Scale(up, 12.0f), Vector3Scale(Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), 14.0f)));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(pumpBasin); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        static double lastPrint = 0.0;
                        const WaterworksStats ws = Waterworks_Stats();
                        if (t - lastPrint >= 5.0) {
                            lastPrint = t;
                            Plot* plot = p->plots ? p->plots->byId(pumpPlot) : nullptr;
                            int pumpPart = -1;
                            if (plot) for (int i = 0; i < (int)plot->parts.size(); ++i) if (plot->parts[i].type == PartType::WaterPump) pumpPart = i;
                            const MachineView v = pumpPart >= 0 ? Factory::View(p->id, pumpPlot, pumpPart) : MachineView{};
                            printf("bench: pump t=%3.0f s: %s, %d networks, pumped %.1f m^3 (%d dry), bowl holds %.1f m^3\n", t, v.valid ? MachineStatusName(v.status) : "-", ws.networks, ws.pumped, ws.dry,
                                   Fluid::WaterAt(*p, pumpBasin, 2.5f));
                        }
                        if (pumpPhase == 1 && t > 60.0) {
                            TakeScreenshot("bench_pump.png");
                            // Cut the power.
                            if (Plot* plot = p->plots ? p->plots->byId(pumpPlot) : nullptr) plot->remove(pumpCable);
                            pumpedAtCut = ws.pumped;
                            pumpPhase = 2;
                        } else if (pumpPhase == 2 && t > 75.0) {
                            const double after = ws.pumped - pumpedAtCut;
                            const float bowl = Fluid::WaterAt(*p, pumpBasin, 2.5f);
                            const bool pass = pumpedAtCut >= 10.0 && bowl >= 5.0f && after <= 1.0;
                            printf("bench: pump: pumped %.1f m^3 in 60 s (%.1f m up), the bowl holds %.1f m^3; after cutting the cable %.1f m^3 more -> %s\n",
                                   pumpedAtCut, pumpRise, bowl, after, pass ? "PASS" : "FAIL");
                            benchPump = false;
                        }
                    }
                } else if (benchBucket && (bucketPhase > 0 || WorldIdle())) {
                    // --test-bucket (Plan 4, F5): scoop river water with a bucket, carry it 50 m inland to a dug
                    // basin and pour it (twice); the basin holds the water and a field beside it counts as watered.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - bucketT0;
                    auto lookFrom = [&](Vector3 standLocal, Vector3 targetLocal) {
                        player.flyMode = false;
                        player.placeAt(world, *p, p->toWorld(standLocal));
                        const Vector3 eye = Vector3Add(standLocal, Vector3Scale(Vector3Normalize(standLocal), player.eyeHeight));
                        player.forward = p->dirToWorld(Vector3Normalize(Vector3Subtract(targetLocal, eye)));
                        player.applyToCamera(world);
                        TerrainTool_Update(world, false, true); // aim only
                    };
                    if (!pf) { printf("bench: bucket: needs the planet\n"); benchBucket = false; }
                    else if (bucketPhase == 0) {
                        // The nearest river or lake with a low shore (as --test-flow).
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f, level = 0.0f, shoreS = 0.0f; Vector3 wdir{}, out{};
                        for (float e = -3000.0f; e <= 3000.0f; e += 16.0f)
                            for (float n = -3000.0f; n <= 3000.0f; n += 16.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if (!ws.kind || ws.kind == 3 || ws.level <= 2.0f || pf->elevation(d) > ws.level - 1.0f) continue;
                                const Vector3 outs[4] = { t1, Vector3Negate(t1), t2, Vector3Negate(t2) };
                                for (const Vector3& o : outs)
                                    for (float s = 2.0f; s <= 20.0f; s += 1.0f) {
                                        const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(o, s / pf->radius)));
                                        if (pf->waterSourceAt(q).kind) continue;
                                        const float eq = pf->elevation(q);
                                        if (eq > ws.level + 0.2f && eq < ws.level + 2.5f) { best = dist; wdir = d; out = o; level = ws.level; shoreS = s; }
                                        break;
                                    }
                            }
                        if (best > 1e29f) { printf("bench: bucket: no water with a low shore -> FAIL\n"); benchBucket = false; }
                        else {
                            auto along = [&](float s) { return Vector3Normalize(Vector3Add(wdir, Vector3Scale(out, s / pf->radius))); };
                            auto ground = [&](float s) { const Vector3 d = along(s); return Vector3Scale(d, pf->radius + pf->elevation(d)); };
                            bucketWater = Vector3Scale(along(shoreS - 3.0f), pf->radius + level - 0.4f);
                            // The basin: a bowl 50 m inland; the field 3 m beyond it.
                            const Vector3 b = ground(shoreS + 50.0f);
                            const TerrainEdit bowl{ Vector3Subtract(b, Vector3Scale(Vector3Normalize(b), 0.8f)), 2.2f, -1.0f, true };
                            ApplyTerrainEdit(*p, bowl);
                            bucketBasin = Vector3Subtract(b, Vector3Scale(Vector3Normalize(b), 2.4f)); // 0.6 m over the bowl's floor
                            bucketField = ground(shoreS + 53.0f);
                            inventory.add("bucket", 1);
                            terrainTool.heldTool = "bucket";
                            printf("bench: bucket: water %.0f m from the start (level %.1f m); basin 50 m inland; field watered before: %s\n", best, level,
                                   Fluid::StandingWaterNear(*p, bucketField, 4.0f) ? "yes" : "no");
                            bucketPhase = 1; bucketT0 = GetTime();
                        }
                    } else if (bucketPhase == 1 && t > 1.0) {
                        // At the shore, looking at the water: scoop.
                        auto along = [&](Vector3 from, float s) { (void)from; return s; };
                        (void)along;
                        const Vector3 up = Vector3Normalize(bucketWater);
                        const Vector3 outDir = Vector3Normalize(Vector3Subtract(bucketBasin, bucketWater));
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(outDir, Vector3Scale(up, Vector3DotProduct(outDir, up))));
                        Vector3 stand = Vector3Add(bucketWater, Vector3Scale(flat, 2.5f)); // at the water's edge
                        stand = Vector3Scale(Vector3Normalize(stand), pf->radius + pf->elevation(Vector3Normalize(stand)) + 0.3f);
                        lookFrom(stand, bucketWater);
                        {
                            const Vector3 dir = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
                            float firstWet = -1.0f;
                            for (float s = 0.4f; s <= 8.0f && firstWet < 0.0f; s += 0.25f) if (IsUnderwater(world, Vector3Add(world.camera.position, Vector3Scale(dir, s)))) firstWet = s;
                            const float aimDist = Vector3Distance(p->toLocal(world.camera.position), bucketWater);
                            printf("bench: bucket: aiming %.1f m at the water, first water along the view at %.2f m, %.1f m of water available\n", aimDist, firstWet,
                                   Fluid::WaterAt(*p, bucketWater, 1.5f));
                        }
                        TerrainTool_Bucket(world);
                        printf("bench: bucket: trip %d: scooped -> water bucket %d\n", bucketTrips + 1, inventory.count("water_bucket"));
                        bucketPhase = 2; bucketT0 = GetTime();
                    } else if (bucketPhase == 2 && t > 1.0) {
                        // At the basin: pour into it.
                        const Vector3 up = Vector3Normalize(bucketBasin);
                        const Vector3 side = Vector3Normalize(Vector3Subtract(bucketField, bucketBasin));
                        Vector3 stand = Vector3Add(bucketBasin, Vector3Scale(side, 3.0f));
                        stand = Vector3Scale(Vector3Normalize(stand), pf->radius + pf->elevation(Vector3Normalize(stand)) + 0.3f);
                        (void)up;
                        lookFrom(stand, bucketBasin);
                        TerrainTool_Bucket(world);
                        ++bucketTrips;
                        printf("bench: bucket: trip %d: poured -> empty bucket %d\n", bucketTrips, inventory.count("bucket"));
                        bucketPhase = bucketTrips < 2 ? 1 : 3; bucketT0 = GetTime();
                    } else if (bucketPhase == 3 && t > 4.0) {
                        const bool basinWet = Fluid::WaterAt(*p, bucketBasin, 1.6f) >= 1.0f;
                        const bool field = Fluid::StandingWaterNear(*p, bucketField, 4.0f);
                        const float vol = (float)Fluid::GetStats().volume;
                        const bool pass = basinWet && field && QuestBook_CountOf("pour") >= 2.0 && inventory.count("bucket") == 1;
                        printf("bench: bucket: basin %s, field 3 m away watered %s, pour events %.0f, %.1f m^3 of moving water -> %s\n", basinWet ? "wet" : "dry", field ? "yes" : "no",
                               QuestBook_CountOf("pour"), vol, pass ? "PASS" : "FAIL");
                        benchBucket = false;
                    }
                } else if (benchPond && (pondPhase > 0 || WorldIdle())) {
                    // --test-pond (Plan 4, F4): drain the nearest pond into a pit dug beside it. The pond is
                    // taken over (no static surface), its level drops, the pit fills; then save, go 300 m away,
                    // load and come back: the pond stays drained and the pit full.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    auto at = [&](Vector3 dir, float h) { return Vector3Scale(dir, pf->radius + h); };
                    auto along = [&](float s) { return Vector3Normalize(Vector3Add(pondDirB, Vector3Scale(pondOutB, s / pf->radius))); };
                    const double t = GetTime() - pondT0;
                    if (!pf) { printf("bench: pond: needs the planet\n"); benchPond = false; }
                    else if (pondPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f;
                        for (float e = -6000.0f; e <= 6000.0f; e += 24.0f)
                            for (float n = -6000.0f; n <= 6000.0f; n += 24.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if (ws.kind != 3) continue;
                                best = dist; pondDirB = ws.pondDir; pondLevelB = ws.level; pondRB = ws.pondR;
                            }
                        if (best > 1e29f) { printf("bench: pond: no pond within 6 km -> FAIL\n"); benchPond = false; }
                        else {
                            pondOutB = Vector3Normalize(Vector3CrossProduct(pondDirB, { 0, 1, 0 }));
                            const float bed = pf->elevation(pondDirB);
                            printf("bench: pond: %.0f m from the start, radius %.0f m, level %.1f m, %.1f m deep at the middle\n", best, pondRB, pondLevelB, pondLevelB - bed);
                            // A trench from the middle out past the rim, its floor below the pond's bed, then a pit.
                            std::vector<TerrainEdit> cut;
                            const float floorH = bed - 1.5f;
                            for (float s = 0.0f; s <= pondRB + 12.0f; s += 0.8f)
                                for (float h = floorH + 1.2f; h <= pondLevelB + 3.0f; h += 1.2f) cut.push_back({ at(along(s), h), 1.4f, -1.0f, true });
                            for (int k = 0; k < 3; ++k) cut.push_back({ at(along(pondRB + 14.0f), floorH - 1.5f - 1.5f * k), 3.5f, -1.0f, true });
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            pondPit = at(along(pondRB + 14.0f), floorH - 4.5f);
                            pondMid = at(along(pondRB * 0.3f), pondLevelB - 0.3f);
                            player.flyMode = true; sunOverride = true;
                            sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(pondDirB, Vector3Scale(pondOutB, 0.5f))));
                            pondPhase = 1; pondT0 = GetTime();
                        }
                    } else if (pondPhase == 1) {
                        // Watch from above the rim.
                        const Vector3 up = pondDirB;
                        const Vector3 eye = Vector3Add(at(along(pondRB + 6.0f), pondLevelB + 18.0f), Vector3Scale(Vector3CrossProduct(up, pondOutB), 22.0f));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(at(along(pondRB * 0.5f), pondLevelB)); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        static double lastPrint = 0.0;
                        if (t - lastPrint >= 5.0) {
                            lastPrint = t;
                            const Fluid::Stats s = Fluid::GetStats();
                            printf("bench: pond t=%3.0f s: %d bricks (%d awake), %.0f wet cells, %.1f m^3, %d ponds taken over; middle %s, pit %s\n", t, s.bricks, s.activeBricks, s.wetCells, s.volume, s.claimedPonds,
                                   Fluid::Wet(*p, pondMid) ? "wet" : "dry", Fluid::Wet(*p, pondPit) ? "wet" : "dry");
                        }
                        if (t > 8.0 && t < 8.1) TakeScreenshot("bench_pond.png");
                        if ((t > 15.0 && Fluid::Asleep()) || t > 150.0) {
                            TakeScreenshot("bench_pond_drained.png");
                            const bool claimed = Fluid::GetStats().claimedPonds > 0;
                            const bool midDry = !Fluid::Wet(*p, pondMid), pitWet = Fluid::Wet(*p, pondPit);
                            const float vol = (float)Fluid::GetStats().volume;
                            printf("bench: pond: settled after %.0f s: taken over %d, the top 0.3 m of the pond %s, the pit %s, %.1f m^3 of water\n", t, (int)claimed, midDry ? "drained" : "still full", pitWet ? "filled" : "dry", vol);
                            const bool saved = SaveWorld(world);
                            // Away 300 m, load, back.
                            player.placeAt(world, *p, p->toWorld(at(along(pondRB + 300.0f), pf->elevation(along(pondRB + 300.0f)) + 2.0f)));
                            const bool loaded = saved && LoadWorldSave(world);
                            printf("bench: pond: saved %d, went 300 m away, loaded %d\n", (int)saved, (int)loaded);
                            pondPhase = 2; pondT0 = GetTime();
                            pondRB = claimed && midDry && pitWet ? pondRB : -pondRB; // remember the first result
                        }
                    } else if (pondPhase == 2) {
                        const Vector3 up = pondDirB;
                        const float r = fabsf(pondRB);
                        const Vector3 eye = Vector3Add(at(along(r + 6.0f), pondLevelB + 18.0f), Vector3Scale(Vector3CrossProduct(up, pondOutB), 22.0f));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(at(along(r * 0.5f), pondLevelB)); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        if (t > 12.0 && Fluid::Asleep()) {
                            TakeScreenshot("bench_pond_loaded.png");
                            const PlanetField::WaterSource ws = pf->waterSourceAt(pondDirB);
                            const bool stillClaimed = ws.kind != 3;
                            const bool midDry = !Fluid::Wet(*p, pondMid), pitWet = Fluid::Wet(*p, pondPit);
                            const bool pass = pondRB > 0.0f && stillClaimed && midDry && pitWet;
                            printf("bench: pond: after the load the pond is %s, its top %s, the pit %s, %.1f m^3 -> %s\n", stillClaimed ? "still taken over" : "back to static water",
                                   midDry ? "drained" : "full", pitWet ? "full" : "dry", Fluid::GetStats().volume, pass ? "PASS" : "FAIL");
                            printf("bench: pond: bench_pond.png, bench_pond_drained.png, bench_pond_loaded.png\n");
                            // Swimming: into the flooded pit, 1.5 m under its surface.
                            player.flyMode = false;
                            player.placeAt(world, *p, p->toWorld(Vector3Add(pondPit, Vector3Scale(Vector3Normalize(pondPit), 1.0f))));
                            pondPhase = 3; pondT0 = GetTime();
                        } else if (t > 120.0) { printf("bench: pond: did not settle after the load -> FAIL\n"); benchPond = false; }
                    } else if (pondPhase == 3) {
                        static float h0 = 0.0f; static int inWaterFrames = 0, frames = 0;
                        const float h = Vector3Length(p->toLocal((player.position - world.renderOrigin).toFloat()));
                        if (t < 0.5 || h0 == 0.0f) { h0 = h; inWaterFrames = frames = 0; pondT0 = GetTime() - std::min(t, 0.5); } // the first frame here can come late (a hitch)
                        else { ++frames; inWaterFrames += player.inWater ? 1 : 0; }
                        if (t > 3.0) {
                            const float sink = (h0 - h) / 2.5f;
                            const bool pass = frames > 0 && inWaterFrames * 10 >= frames * 9 && sink < 1.0f;
                            printf("bench: pond: swimming in the pit: in water %d%% of frames, sank %.2f m/s (falling in air would be ~%.0f m/s) -> %s\n",
                                   frames ? 100 * inWaterFrames / frames : 0, sink, player.gravity * 2.5f * 0.5f, pass ? "PASS" : "FAIL");
                            benchPond = false;
                            h0 = 0.0f;
                        }
                    }
                } else if (benchFlow && (flowPhase > 0 || WorldIdle())) {
                    // --test-flow (Plan 4, F2/F3): cut a channel out of the nearest inland water, ending in a
                    // shaft down to a sealed chamber. Water must reach the channel's end, flood the chamber up
                    // the shaft toward the feeding level, keep the step under budget, and sleep once settled.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    auto at = [&](Vector3 dir, float h) { return Vector3Scale(dir, pf->radius + h); }; // body-local
                    if (!pf) { printf("bench: flow: needs the planet\n"); benchFlow = false; }
                    else if (flowPhase == 0) {
                        // Nearest inland water at least 1.5 m deep with a shore within 25 m.
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f;
                        for (float e = -3000.0f; e <= 3000.0f; e += 16.0f)
                            for (float n = -3000.0f; n <= 3000.0f; n += 16.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if (!ws.kind || ws.kind == 3 || ws.level <= 2.0f || pf->elevation(d) > ws.level - 1.5f) continue;
                                // A shore within 25 m in one of four directions, low enough to cut through.
                                const Vector3 outs[4] = { t1, Vector3Negate(t1), t2, Vector3Negate(t2) };
                                for (const Vector3& o : outs) {
                                    for (float s = 2.0f; s <= 25.0f; s += 1.0f) {
                                        const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(o, s / pf->radius)));
                                        const float eq = pf->elevation(q);
                                        if (pf->waterSourceAt(q).kind) continue;
                                        if (eq > ws.level + 0.3f && eq < ws.level + 4.0f) { best = dist; flowDir = d; flowOut = o; flowLevel = ws.level; flowStart = s; }
                                        break;
                                    }
                                }
                            }
                        if (best > 1e29f) { printf("bench: flow: no inland water with a low shore near the start -> FAIL\n"); benchFlow = false; }
                        else {
                            const float shore = (float)flowStart;
                            printf("bench: flow: water %.0f m from the start, level %.1f m, shore %.0f m out\n", best, flowLevel, shore);
                            // The channel: from 2 m inside the water to 14 m past the shore, its floor 1.5 m under the
                            // water; then a 10 m shaft and a chamber of radius 3.
                            std::vector<TerrainEdit> cut;
                            for (float s = shore - 2.0f; s <= shore + 14.0f; s += 0.8f) {
                                const Vector3 q = Vector3Normalize(Vector3Add(flowDir, Vector3Scale(flowOut, s / pf->radius)));
                                for (float h = flowLevel - 0.3f; h <= flowLevel + 4.0f; h += 1.2f) cut.push_back({ at(q, h), 1.6f, -1.0f, true });
                            }
                            const Vector3 endDir = Vector3Normalize(Vector3Add(flowDir, Vector3Scale(flowOut, (shore + 14.0f) / pf->radius)));
                            for (float h = flowLevel - 1.0f; h >= flowLevel - 12.0f; h -= 0.8f) cut.push_back({ at(endDir, h), 1.3f, -1.0f, true });
                            cut.push_back({ at(endDir, flowLevel - 13.0f), 3.0f, -1.0f, true });
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            flowEnd = at(Vector3Normalize(Vector3Add(flowDir, Vector3Scale(flowOut, (shore + 10.0f) / pf->radius))), flowLevel - 0.6f);
                            flowCave = at(endDir, flowLevel - 14.5f);
                            flowEdit = GetTime();
                            // Camera over the channel's end, looking back at the water.
                            player.flyMode = true; sunOverride = true;
                            sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(flowDir, Vector3Scale(flowOut, -0.5f))));
                            flowPhase = 1;
                        }
                    } else if (flowPhase == 1) {
                        const Vector3 up = Vector3Normalize(flowEnd);
                        const Vector3 eye = Vector3Add(Vector3Add(flowEnd, Vector3Scale(up, 9.0f)), Vector3Scale(flowOut, 9.0f));
                        world.camera.position = p->toWorld(eye);
                        world.camera.target = p->toWorld(Vector3Subtract(flowEnd, Vector3Scale(flowOut, 6.0f)));
                        world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        const double t = GetTime() - flowEdit;
                        if (flowReachT < 0.0 && Fluid::Wet(*p, flowEnd)) flowReachT = t;
                        if (flowCaveT < 0.0 && Fluid::Wet(*p, flowCave)) flowCaveT = t;
                        static double lastPrint = 0.0;
                        if (t - lastPrint >= 5.0) {
                            lastPrint = t;
                            const Fluid::Stats s = Fluid::GetStats();
                            printf("bench: flow t=%3.0f s: %d sites, %d bricks (%d awake), %.0f wet cells, %.1f m^3, step %.2f ms (max %.2f), mesh %.2f ms, %d meshes; changes d%lld s%lld r%lld l%lld\n",
                                   t, s.sites, s.bricks, s.activeBricks, s.wetCells, s.volume, s.stepMs, s.maxStepMs, s.meshMs, s.meshes, s.dbg[0], s.dbg[1], s.dbg[2], s.dbg[3]);
                        }
                        if (flowReachT >= 0.0 && flowCaveT >= 0.0 && flowSleepT < 0.0 && Fluid::Asleep()) flowSleepT = t;
                        if (t > 6.0 && t < 7.0) TakeScreenshot("bench_flow.png");
                        if (getenv("FLOWDBG") && t > 8.0 && t < 8.1) {
                            // Cells along the channel at three heights, and down the shaft.
                            const Vector3 up = Vector3Normalize(flowEnd);
                            for (float s = -6.0f; s <= 6.0f; s += 2.0f)
                                for (float h : { -0.5f, 0.5f, 1.5f }) {
                                    const Vector3 q = Vector3Add(Vector3Add(flowEnd, Vector3Scale(flowOut, s - 10.0f)), Vector3Scale(up, h));
                                    printf("bench: flow dbg along %+.0f m, %+.1f m: %s\n", s - 10.0f, h, Fluid::DebugCell(*p, q).c_str());
                                }
                            for (float h = 1.0f; h >= -15.0f; h -= 1.0f) {
                                const Vector3 q = Vector3Add(Vector3Add(flowEnd, Vector3Scale(flowOut, 4.0f)), Vector3Scale(up, h));
                                printf("bench: flow dbg shaft %+.0f m: %s  density %.2f\n", h, Fluid::DebugCell(*p, q).c_str(), p->densityAt(q));
                            }
                        }
                        if ((flowSleepT >= 0.0 && t > flowSleepT + 1.0) || t > 120.0) {
                            TakeScreenshot("bench_flow_end.png");
                            const Fluid::Stats s = Fluid::GetStats();
                            const float caveTop = Fluid::SurfaceAbove(*p, flowCave, 16.0f);
                            const bool pass = flowReachT >= 0.0 && flowCaveT >= 0.0 && flowSleepT >= 0.0 && s.maxStepMs < 30.0;
                            printf("bench: flow: reached the channel end after %.1f s, the chamber after %.1f s; water stands %.1f m over the chamber floor point (feed level %.1f m above it); asleep after %.1f s; step max %.2f ms -> %s\n",
                                   flowReachT, flowCaveT, caveTop, 14.5f, flowSleepT, s.maxStepMs, pass ? "PASS" : "FAIL");
                            printf("bench: flow: bench_flow.png, bench_flow_end.png\n");
                            benchFlow = false;
                        }
                    }
                } else if (benchForest && (forestPhase > 0 || WorldIdle())) {
                    // --shot-forest: the species mix around the start (60/20/20 per area) and a daylight view
                    // from 30 m up (bench_forest.png).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (forestPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        forestPhase = 1; forestStart = GetTime();
                    } else if (forestPhase == 1 && GetTime() - forestStart > 2.0) {
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                        int n[(int)WoodSpecies::Count];
                        if (p->trees) {
                            for (const float range : { 150.0f, 600.0f }) {
                                p->trees->speciesCounts(me, range, n);
                                int total = 0; for (int s = 0; s < (int)WoodSpecies::Count; ++s) total += n[s];
                                printf("bench: forest: %d trees within %.0f m:", total, range);
                                for (int s = 0; s < (int)WoodSpecies::Count; ++s) printf(" %s %d%%", WoodSpeciesInfo((WoodSpecies)s).name, total ? 100 * n[s] / total : 0);
                                printf("\n");
                            }
                        }
                        forestCam = Vector3Add(me, Vector3Add(Vector3Scale(up, 30.0f), Vector3Scale(t1, -15.0f)));
                        forestTarget = Vector3Add(me, Vector3Scale(t1, 70.0f));
                        forestUp = up;
                        player.flyMode = true;
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(t1, -0.4f))));
                        forestPhase = 2; forestStart = GetTime();
                    } else if (forestPhase == 2) {
                        world.camera.position = p->toWorld(forestCam); world.camera.target = p->toWorld(forestTarget); world.camera.up = p->dirToWorld(forestUp);
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - forestStart > 3.0) {
                            TakeScreenshot("bench_forest.png");
                            // The nearest rock outcrop (within 3 km, on a 25 m grid): look at it next.
                            const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                            const Vector3 up = Vector3Normalize(forestTarget);
                            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                            const Vector3 t2 = Vector3CrossProduct(t1, up);
                            float bestD = 1e9f, bestE = -1e9f; Vector3 bestP{};
                            const int span = outcropHigh ? 240 : 120;
                            for (int i = -span; pf && i <= span; ++i)
                                for (int j = -span; j <= span; ++j) {
                                    const float dd = 25.0f * sqrtf((float)(i * i + j * j));
                                    if (!outcropHigh && dd >= bestD) continue;
                                    const Vector3 dir = Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, 25.0f * i / pf->radius), Vector3Scale(t2, 25.0f * j / pf->radius))));
                                    const float surf = pf->radius + pf->elevation(dir);
                                    if (outcropHigh && (surf <= bestE || surf - pf->radius > 500.0f)) continue;  // grassland, below the snow
                                    const Vector3 q = Vector3Scale(dir, surf + 2.0f);
                                    if (pf->rockFeature(q, surf + 2.0f, surf, 1.0f) > 0.5f) { bestD = dd; bestE = surf; bestP = q; }
                                }
                            if (outcropHigh && bestE > 0.0f) printf("bench: forest: highest outcrop at %.0f m elevation\n", bestE - pf->radius);
                            printf("bench: forest: nearest outcrop %.0f m away\n", bestD);
                            if (bestD < 1e8f) {
                                const Vector3 bu = Vector3Normalize(bestP);
                                const Vector3 bt = Vector3Normalize(Vector3CrossProduct(bu, fabsf(bu.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                                forestCam = Vector3Add(bestP, Vector3Add(Vector3Scale(bu, 14.0f), Vector3Scale(bt, -30.0f)));
                                forestTarget = bestP; forestUp = bu;
                                sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(bu, Vector3Scale(bt, -0.5f))));
                                forestPhase = 3; forestStart = GetTime();
                            } else {
                                printf("bench: forest: bench_forest.png\n");
                                benchForest = false;
                            }
                        }
                    } else if (forestPhase == 3) {
                        world.camera.position = p->toWorld(forestCam); world.camera.target = p->toWorld(forestTarget); world.camera.up = p->dirToWorld(forestUp);
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - forestStart > 5.0) {
                            TakeScreenshot("bench_outcrop.png");
                            printf("bench: forest: bench_forest.png, bench_outcrop.png\n");
                            benchForest = false;
                        }
                    }
                } else if (benchThreats && (threatPhase > 0 || WorldIdle())) {
                    // --test-threats: after dark a dusk stalker comes 45-60 m out and hunts (a bite hurts); a
                    // campfire keeps it off; daylight sends it away; peaceful worlds get none; a thornback
                    // charges when you come close and can be fought off.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: threats: FAILED %s\n", what); threatsOk = false; } };
                    if (threatPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        threatPhase = 1; threatStart = GetTime();
                    } else if (threatPhase == 1 && GetTime() - threatStart > 1.5 && player.grounded) {
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        auto step = [&](float seconds) { for (float t = 0.0f; t < seconds; t += 0.1f) Animals_Update(world, 0.1f, 0.1f, false); };
                        auto nearestNight = [&]() { float best = 1e9f; for (int i = 0; i < Animals_Count(); ++i) if (Animals_Hostile(i)) best = std::min(best, Vector3Distance(Animals_Position(i), me)); return best; };
                        worldDifficulty = 1;
                        Animals_Reset();
                        survival.health = 100.0f;
                        // Night: a stalker comes.
                        world.sunDirection = p->dirToWorld(Vector3Negate(up));
                        // A fixed spawn direction: animals walk the real terrain now (#223) and have no pathfinding yet
                        // (#253), so one that comes out beyond a cliff can't reach the player.
                        SetRandomSeed(3);
                        Animals_RefreshNow();
                        Animals_Update(world, 0.1f, 0.1f, false);
                        threatSpawned = Animals_NightCount();
                        const float d0 = nearestNight();
                        check(threatSpawned == 1 && d0 > 40.0f && d0 < 65.0f, "a stalker comes out of the dark 45-60 m away");
                        step(15.0f);
                        threatBite = 100.0f - survival.health;
                        if (threatBite <= 0.0f) printf("bench: threats: the stalker ended %.1f m away after 15 s\n", nearestNight());
                        check(threatBite > 0.0f, "it hunts the player down and bites");
                        // A campfire beside the player: no more bites; the stalker keeps its distance.
                        survival.health = 100.0f;
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        const Vector3 fireDir = Vector3Normalize(Vector3Add(me, Vector3Scale(Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), 3.0f)));
                        const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                        const Vector3 fireAt = Vector3Scale(fireDir, pf->radius + pf->elevation(fireDir));
                        Plot* fp = p->plots->plotFor(fireAt, 0.0f, true);
                        const int firePart = fp->place(Part{ PartType::Campfire, fp->cellAt(fp->toPlot(Vector3Add(fireAt, Vector3Scale(fireDir, 0.1f)))), 0, 0 });
                        check(firePart >= 0, "a campfire placed");
                        Animals_RefreshNow();
                        step(15.0f);
                        threatFire = nearestNight();
                        check(survival.health >= 100.0f && threatFire > 12.0f, "a campfire keeps it off");
                        if (firePart >= 0) fp->remove(firePart);
                        Animals_RefreshNow();
                        // Daylight: it runs off and is gone.
                        world.sunDirection = p->dirToWorld(up);
                        step(45.0f);
                        threatLeft = Animals_NightCount();
                        check(threatLeft == 0, "gone by day");
                        // Peaceful: nothing comes.
                        worldDifficulty = 0;
                        world.sunDirection = p->dirToWorld(Vector3Negate(up));
                        Animals_RefreshNow();
                        step(2.0f);
                        threatPeaceful = Animals_NightCount();
                        check(threatPeaceful == 0, "no threats on peaceful");
                        // A thornback 8 m away by day: it charges; fought off with an axe.
                        worldDifficulty = 1;
                        world.sunDirection = p->dirToWorld(up);
                        survival.health = 100.0f;
                        Animals_Reset();
                        Animals_Spawn(world, *p, "thornback", Vector3Add(me, Vector3Scale(Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), -8.0f)), 1, false, 0.5f);
                        step(6.0f);
                        threatRam = 100.0f - survival.health;
                        check(threatRam > 0.0f, "a thornback charges when you come close");
                        bool fell = false;
                        for (int k = 0; k < 12 && !fell && Animals_Count() > 0; ++k) fell = Animals_Hit(0, 6.0f);
                        check(fell || Animals_Count() == 0, "fought off");
                        printf("bench: threats: night stalker at %.0f m, bit for %.0f; by a campfire it kept %.0f m off, no bites; gone by day (%d left); peaceful %d; thornback rammed for %.0f -> %s\n",
                               d0, threatBite, threatFire, threatLeft, threatPeaceful, threatRam, threatsOk ? "PASS" : "FAIL");
                        Animals_Reset();
                        survival.health = 100.0f;
                        benchThreats = false;
                    }
                } else if (benchAnimals && (animalPhase > 0 || WorldIdle())) {
                    // --test-animals: a herd spawns on the ground and stays there, a hit one flees, a fed one
                    // follows, two fed adults breed, 200 animals update in under 1 ms, animals are saved.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (animalPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        animalPhase = 1; animalStart = GetTime();
                    } else if (animalPhase == 1 && GetTime() - animalStart > 1.5 && player.grounded) {
                        bool ok = true;
                        auto check = [&](bool c, const char* what) { if (!c) { printf("bench: animals: FAILED %s\n", what); ok = false; } };
                        const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, up));
                        const Vector3 north = Vector3CrossProduct(up, east);
                        auto onGround = [&](int i) { const Vector3 a = Animals_Position(i), u = Vector3Normalize(a); return p->densityAt(Vector3Add(a, Vector3Scale(u, 0.15f))) < ISO_LEVEL && p->densityAt(Vector3Subtract(a, Vector3Scale(u, 0.15f))) >= ISO_LEVEL; }; // on the density surface (#223), not the height map
                        auto step = [&](float seconds) { for (float t = 0.0f; t < seconds; t += 0.1f) Animals_Update(world, 0.1f, 0.1f, false); };
                        animalNatural = Animals_Count(); // herds placed around the player as it arrived
                        check(animalNatural > 0, "herds spawn around the player");
                        Animals_Reset();
                        world.sunDirection = p->dirToWorld(up); // noon: no resting
                        const Vector3 herdAt = Vector3Add(me, Vector3Scale(east, 30.0f));
                        const int made = Animals_Spawn(world, *p, "woolback", herdAt, 6);
                        check(made >= 4, "a herd spawned");
                        bool grounded = true;
                        for (int i = 0; i < Animals_Count(); ++i) grounded = grounded && onGround(i);
                        check(grounded, "spawned on the ground");
                        step(60.0f);
                        float spread = 0.0f; grounded = true;
                        for (int i = 0; i < Animals_Count(); ++i) { spread = std::max(spread, Vector3Distance(Animals_Position(i), herdAt)); grounded = grounded && onGround(i); }
                        check(grounded && spread < 80.0f, "a minute of grazing: still on the ground and near");
                        animalSpread = spread;
                        // Hit one: it runs from the player.
                        const int h = Animals_Nearest(me, 200.0f);
                        const float d0 = Vector3Distance(Animals_Position(h), me);
                        check(!Animals_Hit(h, 1.0f), "a light hit does not kill");
                        step(3.0f);
                        animalFled = Vector3Distance(Animals_Position(h), me) - d0;
                        check(animalFled > 5.0f, "the hit one fled");
                        // Feed one: it follows the player.
                        inventory.add("grain", 4);
                        const int fdx = Animals_Nearest(me, 200.0f);
                        const float f0 = Vector3Distance(Animals_Position(fdx), me);
                        check(Animals_Feed(fdx), "fed");
                        step(12.0f);
                        animalFollow = Vector3Distance(Animals_Position(fdx), me);
                        check(animalFollow < std::max(4.0f, f0 - 8.0f), "the fed one came to the player");
                        // Holding its food (a selected food slot) brings one in without feeding it; the
                        // second press eats.
                        {
                            Animals_Reset();
                            // Behind the player (a ~1.2 m natural ledge lies east of the start, seed 42): animals now go round ledges
                            // and build parts instead of walking through them (#223).
                            const Vector3 ahead = Vector3Normalize(p->dirToLocal(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))));
                            Animals_Spawn(world, *p, "woolback", Vector3Subtract(me, Vector3Scale(ahead, 12.0f)), 1, false, 0.5f);
                            const float h0 = Animals_Count() ? Vector3Distance(Animals_Position(0), me) : 0.0f;
                            inventory.add("berries", 2);
                            Hotbar_Bind(8, "berries");
                            Hotbar_Press(8, false);
                            check(Hotbar_HeldFood() == "berries", "holding berries");
                            step(8.0f);
                            animalLured = Animals_Count() ? Vector3Distance(Animals_Position(0), me) : 99.0f;
                            check(animalLured < 4.0f && h0 > 8.0f, "held food lures an animal");
                            const int g0 = inventory.count("berries");
                            survival.food = 50.0f;
                            Hotbar_Press(8, false);
                            check(inventory.count("berries") == g0 - 1, "the second press eats it");
                            Hotbar_Reset();
                        }
                        // A pen: fences (and a gate) 24 m out around a herd. Scared, they run, but stay inside;
                        // the player walks through the gate, not through a fence.
                        {
                            Animals_Reset();
                            if (!p->plots) p->plots = std::make_unique<PlotSet>();
                            const Vector3 penDir = Vector3Normalize(Vector3Add(me, Vector3Scale(north, 40.0f)));
                            const Vector3 penAt = Vector3Scale(penDir, pf->radius + pf->elevation(penDir));
                            Plot* pen = p->plots->plotFor(penAt, 0.0f, true);
                            int walls = 0, wallsTried = 0;
                            const uint8_t across = (uint8_t)RotateOrientation(0, 1); // a fence turned to run along x
                            Int3 gateCell{ 0, 0, 0 }, fenceCell{ 0, 0, 0 };
                            for (int a = -48; a < 48; a += 2)
                                for (int side = 0; side < 4; ++side) {
                                    const bool alongX = side < 2;
                                    const int cx = alongX ? a : side == 2 ? -48 : 48, cz = alongX ? (side == 0 ? -48 : 49) : a; // turned fences sit a cell toward -z
                                    const Vector3 col = Vector3Normalize(pen->toBody({ (cx + 0.5f) * kCellSize, 0.0f, (cz + 0.5f) * kCellSize }));
                                    const int g = pen->cellAt(pen->toPlot(Vector3Scale(col, pf->radius + pf->elevation(col)))).y;
                                    const bool gate = side == 0 && a == 0;
                                    for (int y = g - 2; y <= g; y += 2) {
                                        Part fp{ gate ? PartType::Gate : PartType::Fence, { cx, y, cz }, alongX ? across : (uint8_t)0, 0 };
                                        ++wallsTried;
                                        walls += pen->place(fp) >= 0 ? 1 : 0;
                                    }
                                    if (gate) gateCell = { cx, g, cz };
                                    if (side == 2 && a == 0) fenceCell = { cx, g, cz };
                                }
                            // The player: a 0.45 m sphere at a gate cell passes, at a fence cell is pushed out.
                            auto pushed = [&](Int3 c) {
                                Vector3 at = pen->toBody({ (c.x + 0.5f) * kCellSize, (c.y + 0.5f) * kCellSize, (c.z + 0.5f) * kCellSize }), n;
                                return pen->collideSphere(at, 0.45f, n);
                            };
                            check(walls == wallsTried, "every fence placed");
                            check(!pushed(gateCell) && pushed(fenceCell), "the player passes the gate, not a fence");
                            Animals_Spawn(world, *p, "woolback", penAt, 6);
                            for (int i = 0; i < Animals_Count(); ++i) Animals_Hit(i, 0.0f);
                            step(30.0f);
                            bool inside = Animals_Count() >= 4;
                            for (int i = 0; i < Animals_Count(); ++i) {
                                const Vector3 q = pen->toPlot(Animals_Position(i));
                                penWorst = std::max(penWorst, std::max(fabsf(q.x), fabsf(q.z)));
                                inside = inside && fabsf(q.x) < 24.0f && fabsf(q.z) < 24.0f;
                            }
                            check(inside, "a fenced pen holds a scared herd");
                            for (int i = 0; i < (int)pen->parts.size(); ++i) pen->remove(i);
                        }
                        // Two fed adults side by side breed.
                        Animals_Reset();
                        Animals_Spawn(world, *p, "woolback", Vector3Add(me, Vector3Scale(north, 10.0f)), 2);
                        check(Animals_Count() == 2, "a pair");
                        if (Animals_Count() == 2) {
                            Animals_Feed(0); Animals_Feed(1);
                            step(15.0f); // they follow the player and meet
                            check(Animals_YoungCount() == 1, "a young one born");
                        }
                        // Kill one: drops.
                        const int meat0 = inventory.count("raw_meat");
                        check(Animals_Hit(0, 999.0f) && inventory.count("raw_meat") > meat0, "meat when it falls");
                        // Save and load.
                        const int before = Animals_Count();
                        const Vector3 pos0 = Animals_Position(0);
                        const std::string path = std::string(GetWorkingDirectory()) + "/animals_test.tmp";
                        if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                        Animals_Reset();
                        if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                        remove(path.c_str());
                        check(Animals_Count() == before && Vector3Distance(Animals_Position(0), pos0) < 0.01f, "saved and loaded");
                        // Aim (#222): an animal 4 m ahead in the crosshair is aimed (the terrain tool holds off); with
                        // ground raised between, the ground is what's aimed at.
                        {
                            Animals_Reset();
                            const Vector3 fwdL = Vector3Normalize(p->dirToLocal(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))));
                            Animals_Spawn(world, *p, "woolback", Vector3Add(me, Vector3Scale(fwdL, 4.0f)), 1, false, 0.0f);
                            const Vector3 eye = Vector3Add(me, Vector3Scale(up, 1.6f));
                            const Vector3 body = Vector3Add(Animals_Position(0), Vector3Scale(Vector3Normalize(Animals_Position(0)), 0.4f));
                            const Camera3D keep = world.camera;
                            world.camera.position = p->toWorld(eye);
                            world.camera.target = p->toWorld(body);
                            Animals_Aim(world, false);
                            const bool aimedOpen = Animals_Aiming();
                            ApplyTerrainEdit(*p, TerrainEdit{ Vector3Lerp(eye, body, 0.5f), 0.7f, 10.0f });
                            Animals_Aim(world, false);
                            const bool aimedBehind = Animals_Aiming();
                            world.camera = keep;
                            check(aimedOpen && !aimedBehind, "aim: an animal in view is aimed, one behind a mound is not");
                            Animals_Reset();
                        }
                        // Terrain (#223): a pit dug under a standing animal drops it onto the pit floor; animals
                        // grazing around a 4 m pit for a minute never end up inside the ground or in the air.
                        {
                            Animals_Reset();
                            const Vector3 pitAt = Vector3Add(me, Vector3Scale(north, 25.0f));
                            Animals_Spawn(world, *p, "woolback", pitAt, 1, false, 0.0f);
                            const float r0 = Vector3Length(Animals_Position(0));
                            ApplyTerrainEdit(*p, TerrainEdit{ Animals_Position(0), 2.0f, -10.0f });
                            step(2.5f);
                            const float dropped = r0 - Vector3Length(Animals_Position(0));
                            check(Animals_Count() == 1 && dropped > 0.8f && onGround(0), "an animal over a dug hole drops onto its floor");
                            Animals_Reset();
                            const Vector3 u = Vector3Normalize(pitAt);
                            for (int k = 1; k < 4; ++k) ApplyTerrainEdit(*p, TerrainEdit{ Vector3Subtract(pitAt, Vector3Scale(u, 1.2f * k)), 2.4f, -10.0f });
                            Animals_Spawn(world, *p, "woolback", Vector3Add(pitAt, Vector3Scale(east, 4.0f)), 6, false, 3.0f);
                            int bad = 0, samples = 0;
                            for (float t = 0.0f; t < 60.0f; t += 0.1f) {
                                Animals_Update(world, 0.1f, 0.1f, false);
                                for (int i = 0; i < Animals_Count(); ++i) { ++samples; bad += onGround(i) ? 0 : 1; }
                            }
                            printf("bench: animals: around a dug pit for a minute: %d of %d samples off the ground, hole drop %.2f m\n", bad, samples, dropped);
                            check(samples > 0 && bad <= samples / 100, "animals stay on the ground around a dug pit");
                        }
                        // 200 within 300 m: the update stays under 1 ms.
                        Animals_Reset();
                        for (int k = 0; k < 20; ++k) Animals_Spawn(world, *p, "woolback", Vector3Add(me, Vector3Add(Vector3Scale(east, -150.0f + 15.0f * k), Vector3Scale(north, 40.0f))), 10);
                        double worst = 0.0;
                        for (int k = 0; k < 60; ++k) { Animals_Update(world, 1.0f / 60.0f, 1.0f / 60.0f, false); worst = std::max(worst, Animals_LastUpdateMs()); }
                        animalCount = Animals_Count(); animalMs = worst;
                        check(animalCount >= 150 && worst < 1.0, "200 animals under 1 ms");
                        Animals_Reset();
                        printf("bench: animals: %d placed around the start; herd grazed within %.0f m for a minute, hit one fled %.1f m, fed one came to %.1f m, held berries lured one to %.1f m, in a pen at most %.1f m from its centre, a pair bred, meat, saved; %d animals in %.3f ms -> %s\n",
                               animalNatural, animalSpread, animalFled, animalFollow, animalLured, penWorst, animalCount, animalMs, ok ? "PASS" : "FAIL");
                        // A look at them: one of each kind ahead of the player, then a screenshot.
                        const Vector3 fwd = Vector3Normalize(p->dirToLocal(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))));
                        const char* kinds[] = { "woolback", "plains_runner", "clucker", "hornhide" };
                        for (int k = 0; k < 4; ++k) Animals_Spawn(world, *p, kinds[k], Vector3Add(me, Vector3Scale(fwd, 14.0f + 4.0f * k)), 2);
                        {   // and a fence with a gate, 6 m ahead
                            const Vector3 fd = Vector3Normalize(Vector3Add(me, Vector3Scale(fwd, 6.0f)));
                            Plot* fp = p->plots->plotFor(Vector3Scale(fd, pf->radius + pf->elevation(fd)), 0.0f, true);
                            for (int k = -3; k <= 3; ++k) fp->place(Part{ k == 0 ? PartType::Gate : PartType::Fence, { 2 * k, 0, 0 }, (uint8_t)RotateOrientation(0, 1), 0 });
                        }
                        world.sunDirection = Vector3Normalize(Vector3Add(p->dirToWorld(up), p->dirToWorld(east)));
                        animalPhase = 2; animalStart = GetTime();
                    } else if (animalPhase == 2 && GetTime() - animalStart > 0.6) {
                        TakeScreenshot("bench_animals.png");
                        Animals_Reset();
                        benchAnimals = false;
                    }
                } else if (benchFarm && (farmPhase > 0 || WorldIdle())) {
                    // --test-farm: till with a shovel's hook, plant grain, no growth at night, ripe in
                    // daylight, harvest gives grain and seeds; a growing crop survives a save.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (farmPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        farmPhase = 1; farmStart = GetTime();
                    } else if (farmPhase == 1 && GetTime() - farmStart > 1.5 && player.grounded) {
                        bool ok = true;
                        auto check = [&](bool c, const char* what) { if (!c) { printf("bench: farm: FAILED %s\n", what); ok = false; } };
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 ahead = Vector3Add(feet, Vector3Scale(Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))), 3.0f));
                        TerrainHit hit;
                        const bool found = RaycastTerrain(world, Ray{ Vector3Add(ahead, Vector3Scale(player.up, 3.0f)), Vector3Negate(player.up) }, 10.0f, hit);
                        check(found, "ground ahead");
                        const bool tilled = found && Farm_Till(world, *p, hit.local);
                        check(tilled, "tilled");
                        Plot* plot = nullptr; int part = -1;
                        if (tilled && p->plots)
                            for (auto& pl : p->plots->plots)
                                for (int i = (int)pl->parts.size() - 1; i >= 0 && part < 0; --i)
                                    if (pl->parts[i].type == PartType::Farmland) { plot = pl.get(); part = i; }
                        check(part >= 0, "farmland part");
                        if (plot) {
                            inventory.add("grain_seeds", 2);
                            check(Farm_Plant(*p, *plot, part, "grain_seeds"), "planted");
                            const Vector3 up = p->dirToWorld(plot->up);
                            std::string crop; float g = 0.0f;
                            world.sunDirection = Vector3Negate(up); // night
                            Farm_Update(world, 2000.0f);
                            Farm_Crop(*p, *plot, part, crop, g);
                            check(crop == "grain" && g == 0.0f, "no growth at night");
                            world.sunDirection = up; // noon
                            Farm_Update(world, 900.0f);
                            Farm_Crop(*p, *plot, part, crop, g);
                            farmQuarter = g;
                            check(g > 0.2f && g < 0.6f, "a quarter (or half, by water) after 900 s of day");
                            Farm_Update(world, 4000.0f);
                            Farm_Crop(*p, *plot, part, crop, g);
                            check(g >= 1.0f, "ripe");
                            const int grain0 = inventory.count("grain"), seeds0 = inventory.count("grain_seeds");
                            check(Farm_Harvest(*p, *plot, part) > 0 && inventory.count("grain") == grain0 + 3 && inventory.count("grain_seeds") > seeds0, "harvest: grain and seeds");
                            // Plant again, half grow, save and load.
                            check(Farm_Plant(*p, *plot, part, "grain_seeds"), "replanted");
                            Farm_Update(world, 600.0f);
                            float before = 0.0f; Farm_Crop(*p, *plot, part, crop, before);
                            const std::string path = std::string(GetWorkingDirectory()) + "/farm_test.tmp";
                            if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                            Farm_Reset();
                            if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                            remove(path.c_str());
                            float after = -1.0f; Farm_Crop(*p, *plot, part, crop, after);
                            check(crop == "grain" && fabsf(after - before) < 0.001f, "a growing crop saved");
                            farmGrowth = before;
                        }
                        printf("bench: farm: tilled, planted; night 0, 900 s of day %.2f, ripe, harvested; replanted at %.2f survived a save -> %s\n", farmQuarter, farmGrowth, ok ? "PASS" : "FAIL");
                        benchFarm = false;
                    }
                } else if (benchQuests && WorldIdle()) {
                    // --test-quests: tier 0's tasks done open tier 1; part of tier 1 done survives a save
                    // and load (tier, done tasks, obtained totals, placed machines); the rest opens tier 2.
                    bool ok = true;
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: quests: FAILED %s\n", what); ok = false; } };
                    inventory.clear();
                    QuestBook_Reset();
                    check(QuestBook_Tier() == 0 && QuestBook_TaskCount(0) >= 4, "tier 0 has tasks");
                    inventory.add("oak_log", 2); inventory.add("fibre", 2); inventory.addAmount("gravel", 4.0); inventory.add("stone_pickaxe", 1);
                    inventory.add("stone", 8); inventory.add("fired_crucible", 1);
                    inventory.remove("oak_log", 2); // used up: obtained totals still count
                    QuestBook_NotePlaced("campfire");
                    QuestBook_Update();
                    check(QuestBook_Tier() == 1, "tier 0 done opens tier 1");
                    QuestBook_Count("eat"); // an event task ("do eat 1")
                    QuestBook_NotePlaced("kiln");
                    inventory.add("copper_ore", 4);
                    QuestBook_Update();
                    check(QuestBook_Tier() == 1 && QuestBook_TaskDone(1, 0) && QuestBook_TaskDone(1, 1) && QuestBook_TaskDone(1, 2) && !QuestBook_TaskDone(1, 3), "part of tier 1");
                    const std::string path = std::string(GetWorkingDirectory()) + "/quest_test.tmp";
                    if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                    QuestBook_Reset(); inventory.clear();
                    if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                    remove(path.c_str());
                    check(QuestBook_Tier() == 1 && QuestBook_TaskDone(1, 0) && QuestBook_TaskDone(1, 1) && QuestBook_TaskDone(1, 2) && !QuestBook_TaskDone(1, 3) && inventory.obtainedTotal("any_log") >= 2.0 && QuestBook_CountOf("eat") >= 1.0, "saved and loaded");
                    inventory.add("tin_ore", 2); inventory.add("copper_ingot", 3); inventory.add("bronze_ingot", 1);
                    QuestBook_Update();
                    check(QuestBook_Tier() == 2, "tier 1 done opens tier 2");
                    // The quest tree: side nodes open by their requirements, complete by their task, pay a
                    // reward, and are saved.
                    check(QuestTree_NodeCount() >= 12, "side lines loaded");
                    check(QuestTree_IsOpen("f1") && QuestTree_IsOpen("h1") && QuestTree_IsOpen("s1") && !QuestTree_IsOpen("s2") && QuestTree_IsOpen("e1"), "nodes open by tier and by their requirements");
                    const int berries0 = inventory.count("berries");
                    QuestBook_AtLeast("explored", 25.0);
                    QuestBook_Update();
                    check(QuestTree_IsDone("e1") && inventory.count("berries") == berries0 + 4 && QuestTree_IsOpen("e2") && !QuestTree_IsOpen("e3"), "a side node done pays its reward and opens the next");
                    if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                    QuestBook_Reset();
                    check(!QuestTree_IsDone("e1"), "reset clears the tree");
                    if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                    remove(path.c_str());
                    check(QuestTree_IsDone("e1"), "the tree saved");
                    printf("bench: quests: %d tiers, tier 0 -> 1 -> 2 by their tasks, save/load kept the tier and tasks -> %s\n", QuestBook_TierCount(), ok ? "PASS" : "FAIL");
                    inventory.clear(); QuestBook_Reset();
                    benchQuests = false;
                } else if (benchViewmodel && (vmPhase > 0 || WorldIdle())) {
                    // --test-viewmodel: standing still the hands don't move; walking they bob within the
                    // set amplitude, once a step; screenshots of the hands and each tool type.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - vmStart;
                    const Vector3 b = Viewmodel_BobOffset();
                    if (vmPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        terrainTool.heldTool.clear();
                        vmStart = GetTime(); vmPhase = 1;
                    } else if (vmPhase == 1 && t > 1.5 && player.grounded) {
                        vmLast = b; vmStillMax = 0.0f; vmPhase = 2; vmStart = GetTime();
                    } else if (vmPhase == 2) {
                        vmStillMax = std::max(vmStillMax, Vector3Distance(b, vmLast)); vmLast = b;
                        if (t > 2.0) { vmFrom = p->toLocalD(player.position); player.autoWalk = true; vmMaxX = vmMaxY = 0.0f; vmSigns = 0; vmPhase = 3; vmStart = GetTime(); }
                    } else if (vmPhase == 3) {
                        vmMaxX = std::max(vmMaxX, fabsf(b.x)); vmMaxY = std::max(vmMaxY, fabsf(b.y));
                        if ((b.x > 0.0f) != (vmLast.x > 0.0f) && fabsf(b.x) > 1e-5f) ++vmSigns;
                        vmLast = b;
                        if ((p->toLocalD(player.position) - vmFrom).length() > 10.0 || t > 6.0) {
                            player.autoWalk = false;
                            vmWalked = (float)(p->toLocalD(player.position) - vmFrom).length();
                            player.forward = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, 0.35f)));
                            vmPhase = 4; vmStart = GetTime();
                        }
                    } else if (vmPhase >= 4 && vmPhase <= 7 && t > 1.2) {
                        static const char* shots[4] = { "bench_hands.png", "bench_hands_pickaxe.png", "bench_hands_axe.png", "bench_hands_shovel.png" };
                        static const char* next[4] = { "stone_pickaxe", "iron_axe", "steel_shovel", "" };
                        TakeScreenshot(shots[vmPhase - 4]);
                        if (next[vmPhase - 4][0]) { // into the hand as a player does: on a hotbar key, pressed (setting heldTool alone drew bare hands)
                            inventory.add(next[vmPhase - 4], 1);
                            Hotbar_Bind(vmPhase - 3, next[vmPhase - 4]);
                            Hotbar_Press(vmPhase - 3, false);
                        }
                        vmPhase++; vmStart = GetTime();
                    } else if (vmPhase == 8) {
                        const bool pass = vmStillMax < 1e-5f && vmMaxX <= 0.0101f && vmMaxY <= 0.0151f && vmMaxX > 0.003f && vmSigns >= 4; // a side swing every 1.8 m walked: about 5 in 10 m
                        printf("bench: viewmodel: standing still moved %.5f m; walking %.1f m bobbed up to %.4f m across, %.4f m down, %d step swings; drawn %d -> %s\n",
                               vmStillMax, vmWalked, vmMaxX, vmMaxY, vmSigns, Viewmodel_Visible() ? 1 : 0, pass ? "PASS" : "FAIL");
                        benchViewmodel = false;
                    }
                } else if (benchMap && (mapPhase > 0 || WorldIdle())) {
                    // --test-map: the compass reads east when facing east; walking marks the ground;
                    // a waypoint north bears 0; the full map and the minimap render; exploration and
                    // waypoints survive a save and load.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - mapStart;
                    const Vector3 up = player.up;
                    const Vector3 axis = p->dirToWorld({ 0, 1, 0 });
                    const Vector3 east = Vector3Normalize(Vector3CrossProduct(axis, up));
                    const Vector3 north = Vector3CrossProduct(up, east);
                    if (mapPhase == 0) {
                        if (inventory.count("map") == 0) inventory.add("map", 1);
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        player.forward = east;
                        mapStart = GetTime();
                        mapPhase = 1;
                    } else if (mapPhase == 1 && t > 1.0 && player.grounded && player.nearest) {
                        player.forward = east;
                        mapHeading = Nav_Heading(world);
                        mapFrom = p->toLocalD(player.position);
                        player.autoWalk = true;
                        mapPhase = 2; mapStart = GetTime();
                    } else if (mapPhase == 2 && t > 40.0) {
                        player.autoWalk = false;
                        const double walked = (p->toLocalD(player.position) - mapFrom).length();
                        mapTiles = Nav_ExploredTiles(p->name);
                        const Vector3 local = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 lup = Vector3Normalize(local), lnorth = p->dirToLocal(north);
                        const Vector3 wp = Vector3Normalize(Vector3Add(Vector3Scale(lup, p->radius * p->scale), Vector3Scale(lnorth, 200.0f)));
                        Nav_AddWaypoint(p->name, wp, "test");
                        mapBearing = Nav_Bearing(world, wp);
                        printf("bench: map: facing east the compass reads %.1f deg; walked %.0f m, %d tiles explored; a waypoint 200 m north bears %.1f deg\n",
                               mapHeading, walked, mapTiles, mapBearing);
                        Nav_ToggleMap();
                        mapPhase = 3; mapStart = GetTime();
                    } else if (mapPhase == 3 && (Nav_MapImageReady() || t > 20.0)) {
                        mapImageOk = Nav_MapImageReady();
                        if (mapImageOk) TakeScreenshot("bench_map.png");
                        mapPhase = 4; mapStart = GetTime();
                    } else if (mapPhase == 4 && t > 0.5) {
                        Nav_CloseMap();
                        mapPhase = 5; mapStart = GetTime();
                    } else if (mapPhase == 5 && (Nav_MinimapReady() && t > 1.0 || t > 15.0)) {
                        const bool miniOk = Nav_MinimapReady();
                        if (miniOk) TakeScreenshot("bench_minimap.png");
                        // Save and load the extras: exploration and waypoints must come back.
                        const int tiles = Nav_ExploredTiles(p->name), wps = Nav_WaypointCount();
                        const std::string path = std::string(GetWorkingDirectory()) + "/nav_test.tmp";
                        if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                        Nav_Reset();
                        if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                        remove(path.c_str());
                        const bool saved = Nav_ExploredTiles(p->name) == tiles && Nav_WaypointCount() == wps;
                        const bool pass = fabsf(mapHeading - 90.0f) < 2.0f && mapTiles >= 2 && (mapBearing < 3.0f || mapBearing > 357.0f) && mapImageOk && miniOk && saved;
                        printf("bench: map: full map rendered %s, minimap %s, save/load %d tiles %d waypoints -> %s: %s\n", mapImageOk ? "yes" : "NO", miniOk ? "yes" : "NO",
                               Nav_ExploredTiles(p->name), Nav_WaypointCount(), saved ? "kept" : "LOST", pass ? "PASS" : "FAIL");
                        benchMap = false;
                    }
                } else if (benchInventory && WorldIdle()) {
                    // --test-inventory: the grid's stacking rules, fractions, crafting room, sorting and
                    // the save round trip.
                    bool ok = true;
                    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: inventory: FAILED %s\n", what); ok = false; } };
                    inventory.clear();
                    inventory.add("telescope", 3);
                    int tools = 0;
                    for (const InvSlot& s : inventory.grid()) if (s.id == "telescope") { ++tools; check(s.centi == 100, "tools stack to 1"); }
                    check(tools == 3, "three telescopes in three slots");
                    inventory.add("conveyor", 300);
                    int conv = 0;
                    for (const InvSlot& s : inventory.grid()) if (s.id == "conveyor") ++conv;
                    check(conv == 2 && inventory.count("conveyor") == 300, "blocks stack to 256");
                    inventory.add("stone", 9000);
                    check(inventory.count("stone") == 9000 && inventory.grid()[5].centi == 819200, "resources stack to 8192");
                    for (int i = 0; i < 1000; ++i) inventory.addFractional("iron_ore", 0.0123f);
                    check(fabs(inventory.amount("iron_ore") - 12.30) < 0.02 && inventory.count("iron_ore") == 12, "fractional amounts to 0.01");
                    check(inventory.removeAmount("iron_ore", 0.25) && fabs(inventory.amount("iron_ore") - 12.05) < 0.02, "remove a fraction");
                    check(!inventory.remove("stone", 9001), "can't remove more than held");
                    inventory.sort();
                    check(KindOf(inventory.grid()[0].id) == ItemKind::Tool && inventory.count("stone") == 9000, "sort: tools first, nothing lost");
                    // Full: 45 slots of 8192 regolith; the rest is refused and flagged.
                    inventory.clear();
                    const double left = inventory.addAmount("regolith", 8192.0 * Inventory::kSlots + 10.0);
                    check(fabs(left - 10.0) < 0.01 && inventory.overflowed, "a full grid refuses the rest");
                    inventory.overflowed = false;
                    // Crafting must not run when the output has no room.
                    const Recipe* rr = nullptr;
                    for (const Recipe& r : itemDB.recipes()) if (r.output == "rope") { rr = &r; break; }
                    if (rr) {
                        // Inputs that leave a remainder free no cell: with the rest full, no room for rope.
                        inventory.clear();
                        for (auto& in : rr->inputs) inventory.add(in.first, in.second * 2);
                        inventory.addAmount("regolith", 1e9); inventory.overflowed = false;
                        check(!inventory.canCraft(*rr), "no craft without room for the output");
                        // Exactly the inputs: using them up frees their cells for the output.
                        inventory.clear();
                        for (auto& in : rr->inputs) inventory.add(in.first, in.second);
                        inventory.addAmount("regolith", 1e9); inventory.overflowed = false;
                        check(inventory.canCraft(*rr) && inventory.craft(*rr) && inventory.count("rope") == rr->outputCount, "craft into the inputs' freed cells");
                    }
                    // Save round trip through the extras hooks.
                    inventory.clear();
                    inventory.add("telescope", 1); inventory.addAmount("copper_ore", 12.34); inventory.add("kiln", 2);
                    inventory.swapSlots(0, 7);
                    terrainTool.creative = true; // saved with the world
                    const std::string path = std::string(GetWorkingDirectory()) + "/inv_test.tmp";
                    if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                    inventory.clear();
                    if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                    remove(path.c_str());
                    check(inventory.grid()[7].id == "telescope" && fabs(inventory.amount("copper_ore") - 12.34) < 0.001 && inventory.count("kiln") == 2, "save and load keep slots and fractions");
                    check(terrainTool.creative, "creative mode is saved");
                    terrainTool.creative = false;
                    printf("bench: inventory: %s\n", ok ? "PASS" : "FAIL");
                    inventory.clear();
                    benchInventory = false;
                } else if (benchHotbar && WorldIdle()) {
                    // --test-hotbar: materials from the hotbar (any shape costs its material), a consumable,
                    // the foundation not digging outside the part, and materials kept through a save.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    bool ok = pf != nullptr;
                    if (ok) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        inventory.clear();
                        Hotbar_Reset(); // empty: bind what the test uses (slot 1 stays empty: bare hands)
                        Hotbar_Bind(1, "concrete"); Hotbar_Bind(3, "oak_planks"); Hotbar_Bind(8, "oxygen_canister");
                        inventory.add("concrete", 20); inventory.add("oak_planks", 10); inventory.add("oxygen_canister", 1);
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
                        const int w0 = inventory.count("oak_planks");
                        const bool beam = BuildTool_Place(*p, *plot, { PartType::Beam4, { 2, 0, 0 }, 0, (uint8_t)buildTool.material }, true);
                        const int w1 = inventory.count("oak_planks");
                        // Foundation: the ground a little outside a 1 m block keeps its density.
                        auto outside = [&](float off) { return plot->toBody({ -0.2f + off, 0.3f, 2.5f }); };
                        const float before = p->densityAt(outside(0.0f));
                        Hotbar_Press(1, false);
                        const bool block = BuildTool_Place(*p, *plot, { PartType::Block, { 0, -1, 4 }, 0, (uint8_t)buildTool.material }, true);
                        const float after = p->densityAt(outside(0.0f));
                        survival.oxygen = 50.0f;
                        Hotbar_Press(8, false); // key 9: oxygen canister
                        const bool used = inventory.count("oxygen_canister") == 0 && survival.oxygen > 90.0f;
                        // #272: in build mode, a key bound to a pickaxe not carried changes nothing (the selection, build
                        // mode); carried, it goes into the hand and build mode ends.
                        Hotbar_Press(1, false);
                        Hotbar_Bind(5, "titanium_pickaxe");
                        const int sel0 = Hotbar_Selected();
                        Hotbar_Press(5, false);
                        const bool unownedKept = Hotbar_Selected() == sel0 && buildTool.active && terrainTool.heldTool != "titanium_pickaxe";
                        inventory.add("titanium_pickaxe", 1);
                        Hotbar_Press(5, false);
                        const bool ownedTaken = Hotbar_Selected() == 5 && !buildTool.active && terrainTool.heldTool == "titanium_pickaxe";
                        printf("bench: hotbar: unowned pickaxe key in build mode: selection and build mode kept %d; owned: in hand, build mode off %d (#272)\n", unownedKept, ownedTaken);
                        Hotbar_Press(0, false); // key 1: an empty slot, bare hands
                        const bool toolMode = !buildTool.active && unownedKept && ownedTaken;
                        const bool saved = SaveWorld(world), loaded = saved && LoadWorldSave(world);
                        Planetoid* q = dynamic_cast<Planetoid*>(world.objects[0].get());
                        int woodParts = 0, concreteParts = 0;
                        if (q && q->plots)
                            for (auto& pl : q->plots->plots)
                                for (const Part& part : pl->parts) {
                                    if (part.type == PartType::Count) continue;
                                    if (PartItem(part) == std::string("oak_planks")) ++woodParts;
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
                } else if (benchStation && (stationStarted || WorldIdle())) {
                    if (!stationStarted) { stationStarted = true; stationWait = GetTime(); }
                    else if (GetTime() - stationWait > 1.5 && StationTest_Update(world) != 0) benchStation = false;
                } else if (benchRocket && (rocketStarted || WorldIdle())) {
                    if (!rocketStarted) { rocketStarted = true; Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get()); player.spawnOn(world, *p, StartSpawnDir(*p)); rocketWait = GetTime(); }
                    else if (GetTime() - rocketWait > 1.5 && RocketTest_Update(world) != 0) benchRocket = false;
                } else if (benchBlueprint && WorldIdle()) {
                    // --test-blueprint: the starter rocket preset, placed from the inventory (one landing leg
                    // short: it places the rest and says what is missing; placing again fills it in), then
                    // saved as a blueprint and loaded back.
                    benchBlueprint = false;
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    int idx = -1;
                    for (size_t i = 0; i < Blueprint_List().size(); ++i) if (Blueprint_List()[i].name == "Starter rocket") idx = (int)i;
                    if (!pf || idx < 0) { printf("bench: blueprint: FAIL (no planet or no Starter rocket preset)\n"); }
                    else {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        const Vector3 home = StartSpawnDirLocal(*p);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 }));
                        const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Scale(t1, 16.0f / pf->radius)));
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(Vector3Scale(d, pf->radius + pf->elevation(d)), 0.0f, true);
                        inventory.clear();
                        inventory.add("rocket_engine", 1); inventory.add("landing_leg", 3); inventory.add("fuel_tank", 4);
                        inventory.add("guidance_computer", 1); inventory.add("cockpit", 1);
                        const Blueprint& bp = Blueprint_List()[idx];
                        std::string missing;
                        const int first = Blueprint_Place(*p, *plot, bp, { 0, 0, 0 }, true, missing);
                        inventory.add("landing_leg", 1);
                        std::string missing2;
                        const int second = Blueprint_Place(*p, *plot, bp, { 0, 0, 0 }, true, missing2);
                        const bool saved = Blueprint_SaveFile(*plot, "saves/bench_blueprint.bp", "bench copy");
                        Blueprint back;
                        const bool loaded = Blueprint_LoadFile("saves/bench_blueprint.bp", back);
                        const bool ok = first == 10 && missing == "1 landing_leg" && second == 1 && plot->partCount() == (int)bp.parts.size() && saved && loaded && back.parts.size() == bp.parts.size();
                        printf("bench: blueprint: %s (%zu parts): placed %d (missing: %s), then %d more; saved and loaded %zu parts -> %s\n",
                               bp.name.c_str(), bp.parts.size(), first, missing.c_str(), second, back.parts.size(), ok ? "PASS" : "FAIL");
                    }
                } else if (benchSmallBody && (sbPhase > 0 || WorldIdle())) {
                    // --test-smallbody: drop out of fly mode 30 m over a huge ring rock (pebble), then over a
                    // small ring rock: land, stand, and no death; then dig at the feet on the moon (no launch).
                    const double t = GetTime() - sbStart;
                    auto findBody = [&](const char* name) -> Planetoid* { for (auto& o : world.objects) if (auto* q = dynamic_cast<Planetoid*>(o.get()); q && q->name == name) return q; return nullptr; };
                    auto dropOver = [&](Planetoid& b, float height) {
                        const DVec3 c = b.worldPos;
                        const Vector3 dir = Vector3Normalize({ 0.3f, 1.0f, 0.2f });
                        float r = b.radius * b.scale;
                        while (r > 1.0f && b.densityAt(Vector3Scale(b.dirToLocal(dir), r)) < ISO_LEVEL) r -= 0.5f;
                        player.flyMode = false;
                        player.placeAt(world, b, (c + DVec3(Vector3Scale(dir, r + height)) - world.renderOrigin).toFloat());
                        sbMaxImpact = 0.0f; sbHealth = survival.health; sbDeaths = survival.deaths; sbGrounded = 0; sbFrames = 0;
                    };
                    if (sbPhase > 0 && sbPhase % 2 == 1) { sbMaxImpact = std::max(sbMaxImpact, player.impactSpeed); ++sbFrames; if (player.grounded) ++sbGrounded; }
                    if (sbPhase == 0) {
                        Planetoid* pebble = findBody("pebble");
                        if (!pebble) { printf("bench: smallbody: no pebble\n"); benchSmallBody = false; }
                        else { player.flyMode = true; world.setCameraWorld(pebble->worldPos + DVec3(0, pebble->radius * 2.0, 0), pebble->worldPos); sbPhase = 1; sbStart = GetTime(); }
                    } else if (sbPhase == 1 && t > 2.0 && !sbDropped) {
                        sbDropped = true; dropOver(*findBody("pebble"), 30.0f); sbStart = GetTime();
                    } else if (sbPhase == 1 && sbDropped && t > 12.0) {
                        printf("bench: smallbody: pebble: worst impact %.1f m/s, grounded %d%% of the last frames, deaths %d, health %.0f, gravity %.2f\n",
                               sbMaxImpact, 100 * sbGrounded / std::max(1, sbFrames), survival.deaths - sbDeaths, survival.health, player.gravity);
                        sbOk &= survival.deaths == sbDeaths && player.grounded;
                        // A small ring rock next: the nearest one to pebble's lane, promoted when the camera is near.
                        RingSystem& ring = *Rings()[0];
                        int k = -1; float best = 0.0f;
                        const DVec3 here = player.position;
                        for (size_t i = 0; i < ring.rocks.size(); ++i) {
                            if (ring.rocks[i].radius <= 0.0f || ring.rocks[i].radius * 0.72f > 30.0f) continue;
                            const double d = (ring.rockPosition((int)i) - here).length();
                            if (k < 0 || d < best) { best = (float)d; k = (int)i; }
                        }
                        sbRock = k;
                        player.flyMode = true;
                        world.setCameraWorld(ring.rockPosition(k) + DVec3(0, 60, 0), ring.rockPosition(k));
                        sbPhase = 2; sbStart = GetTime(); sbDropped = false;
                    } else if (sbPhase == 2) {
                        RingSystem& ring = *Rings()[0];
                        world.setCameraWorld(ring.rockPosition(sbRock) + DVec3(0, 60, 0), ring.rockPosition(sbRock));
                        if (Planetoid* b = ring.promotedBody(sbRock); b && t > 1.0) { dropOver(*b, 3.0f); sbPhase = 3; sbStart = GetTime(); }
                        else if (t > 10.0) { printf("bench: smallbody: rock %d never promoted\n", sbRock); sbOk = false; sbPhase = 4; }
                    } else if (sbPhase == 3 && t < 12.0) {
                    } else if (sbPhase == 3 && t >= 12.0) {
                        printf("bench: smallbody: ring rock %d (%.0f m across): worst impact %.1f m/s, grounded %d%%, deaths %d, health %.0f, gravity %.2f\n", sbRock,
                               Rings()[0]->rocks[sbRock].radius * 1.44f, sbMaxImpact, 100 * sbGrounded / std::max(1, sbFrames), survival.deaths - sbDeaths, survival.health, player.gravity);
                        sbOk &= survival.deaths == sbDeaths && player.grounded;
                        sbPhase = 4;
                    } else if (sbPhase == 4) {
                        printf("bench: smallbody: %s\n", sbOk ? "PASS" : "FAIL");
                        benchSmallBody = false;
                    }
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
                                // Then dig at the feet for 4 s: the ground going away must not launch the player.
                                inventory.add("steel_shovel", 1); terrainTool.heldTool = "steel_shovel"; // a wide bite: the ground goes
                                player.forward = Vector3Normalize(Vector3Add(Vector3Scale(player.up, -0.9f), Vector3Scale(player.forward, 0.3f)));
                                moonTakeoff = local; moonMaxUp = 0.0f; moonMaxRise = 0.0; moonMinRise = 0.0; moonStrokes = 0;
                                moonPhase = 4;
                                moonStart = GetTime();
                            }
                        } else if (moonPhase == 4) {
                            const Vector3 fe = (player.position - world.renderOrigin).toFloat();
                            const Vector3 rel = player.grounded ? player.velocity : Vector3Subtract(player.velocity, moon->surfaceVelocity(fe));
                            moonMaxUp = std::max(moonMaxUp, Vector3DotProduct(rel, player.up));
                            moonMaxRise = std::max(moonMaxRise, local.length() - moonTakeoff.length()); moonMinRise = std::min(moonMinRise, local.length() - moonTakeoff.length());
                            if (TerrainTool_Stroke(true)) ++moonStrokes;
                            if (t > 4.0) {
                                const bool ok = moonMaxUp < 3.0f && moonMaxRise < 1.5;
                                printf("bench: dig at the feet on the moon: %d strokes, worst upward speed %.2f m/s, rose at most %.2f m, sank %.2f m: %s\n",
                                       moonStrokes, moonMaxUp, moonMaxRise, -moonMinRise, ok ? "PASS" : "FAIL");
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
                        // Ore: found at every depth, most of it near the surface (deposits plus pockets of the common ores).
                        const bool oreOk = orePct[2] > 0.3f && orePct[0] >= orePct[1] && orePct[1] >= orePct[2] * 0.8f;
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
                    // quarter at 2R, no jumps between bodies).
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
                        benchGravity = false; // (the placeholder ship and its drop / hover checks are gone)
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
                        for (const ItemDef& d : itemDB.items()) if (inventory.amount(d.id) > 0.0) inv += " " + d.id + "=" + TextFormat("%.2f", inventory.amount(d.id));
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
                    if (!benchStandOnly && GetTime() - walkStart > 1.0) {
                        // Walking: the eye height (ground frame) should change smoothly. A jolt is a change of
                        // more than 3 mm in its per-frame step while grounded (bumps felt as camera shake).
                        // The step comes from the player's own ground-frame bookkeeping: converting the camera
                        // into the spinning planet's frame here added cm of float rounding at 60 km.
                        static double lvz = 0.0; static int have = 0;
                        if (player.eyeStepKnown) {
                            const double vz = player.eyeStep * 1000.0; // mm this frame
                            if (have >= 1 && player.grounded) {
                                const double j = fabs(vz - lvz);
                                walkJoltMax = std::max(walkJoltMax, j);
                                if (j > 3.0) ++walkJolts;
                                ++walkJoltFrames;
                                walkSlopeMax = std::max(walkSlopeMax, (double)acosf(std::clamp(Vector3DotProduct(player.groundNormal, player.up), -1.0f, 1.0f)) * RAD2DEG);
                                if (acosf(std::clamp(Vector3DotProduct(player.groundNormal, player.up), -1.0f, 1.0f)) * RAD2DEG > 20.0f) { ++walkSteepFrames; if (j > 3.0) ++walkSteepJolts; }
                            }
                            lvz = vz; have = 1;
                        } else {
                            have = 0;
                        }
                    }
                    if (benchStandOnly &&GetTime() - walkStart > 1.0) {
                        // Standing still: the camera must hold still in the ground's frame (no jitter).
                        static DVec3 lastCam{}; static Vector3 lastFwd{}, lastUp{}; static bool have = false;
                        const DVec3 cam = start->toLocalD(world.renderOrigin + world.camera.position);
                        const Vector3 fwd = start->dirToLocal(Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position)));
                        const Vector3 up = start->dirToLocal(world.camera.up);
                        if (have) {
                            standCamMax = std::max(standCamMax, (cam - lastCam).length());
                            if ((cam - lastCam).length() > 0.002) ++standJitterFrames;
                            standFwdMax = std::max(standFwdMax, acosf(std::clamp(Vector3DotProduct(fwd, lastFwd), -1.0f, 1.0f)) * RAD2DEG);
                            standUpMax = std::max(standUpMax, acosf(std::clamp(Vector3DotProduct(up, lastUp), -1.0f, 1.0f)) * RAD2DEG);
                            standFrames++;
                        }
                        {
                            // A fixed point on the ground (the feet at 1 s) in shadow-map texels: it must not
                            // slide across the texel grid while standing (only the sun's slow motion).
                            static bool havePt = false;
                            if (!havePt) {
                                havePt = true;
                                const DVec3 groundPt = start->toLocalD(player.position);
                                world.shadowProbe = [start, groundPt] { return start->toWorldD(groundPt); };
                                world.shadowProbeSlide = 0.0f;
                            }
                            standTexelSlide = world.shadowProbeSlide;
                        }
                        lastCam = cam; lastFwd = fwd; lastUp = up; have = true;
                    }
                    { const DVec3 here = start->toLocalD(player.position); walkPath += (here - walkLast).length(); walkLast = here; }
                    if (!benchStandOnly && player.autoWalk) {
                        // Getting unstuck, as a player would (overhangs, cliffs and boulders are part of the world, the
                        // walker adapts): every 1.5 s, under 1 m of progress means stuck - jump first; still stuck, turn
                        // aside (alternating sides, wider each time) and walk on.
                        static double checkAt = 0.0; static DVec3 checkPos{}; static int stuck = 0; static int unsticks = 0;
                        if (GetTime() - walkStart < 0.1) { checkAt = GetTime() + 1.5; checkPos = start->toLocalD(player.position); stuck = 0; unsticks = 0; }
                        if (GetTime() > checkAt) {
                            const DVec3 here = start->toLocalD(player.position);
                            if ((here - checkPos).length() < 1.0) {
                                ++stuck; ++unsticks;
                                if (stuck == 1) player.autoJump = true;
                                else player.forward = Vector3Normalize(Vector3RotateByQuaternion(player.forward, QuaternionFromAxisAngle(player.up, (stuck % 2 ? 1.0f : -1.0f) * DEG2RAD * std::min(70.0f * (stuck - 1), 160.0f))));
                            } else stuck = 0;
                            checkPos = here; checkAt = GetTime() + 1.5;
                        }
                        walkUnsticks = unsticks;
                    }
                    if (GetTime() - walkStart > benchWalkSeconds || (walkTargetM > 0.0 && walkPath >= walkTargetM)) {
                        player.autoWalk = false;
                        const double dist = (start->toLocalD(player.position) - walkFrom).length();
                        const float alt = Vector3Distance(feet, start->position);
                        printf("bench: walked %.1f m in %.1f s, grounded %.0f%% of frames, density at body %.2f..%.2f (solid >= %.2f), %.1f m from centre\n",
                               dist, GetTime() - walkStart, 100.0 * walkGroundedFrames / std::max(1, walkFrames), walkMinDensity, walkMaxDensity, ISO_LEVEL, alt);
                        if (!benchStandOnly) printf("bench: walk: got unstuck %d times (jumped, then turned aside)\n", walkUnsticks);
                        if (!benchStandOnly) printf("bench: walk camera: %d jolts (eye height step changes > 3 mm between frames) in %d grounded frames (%d of %d on ground over 20 deg), worst %.1f mm, steepest ground %.0f deg\n",
                                                    walkJolts, walkJoltFrames, walkSteepJolts, walkSteepFrames, walkJoltMax, walkSlopeMax);
                        if (benchStandOnly) printf("bench: standing still %d frames: camera moved at most %.4f m per frame in the ground's frame, view turned %.4f deg, up %.4f deg; the ground slid %.2f shadow texels/s: %s\n", standFrames, standCamMax, standFwdMax, standUpMax, standTexelSlide / std::max(0.1, GetTime() - walkStart - 1.0),
                                                   standJitterFrames * 100 <= standFrames && standTexelSlide / std::max(0.1, GetTime() - walkStart - 1.0) < 10.0 ? "PASS" : "FAIL");
                        if (walkTargetM > 0.0) {
                            const double pct = 100.0 * walkGroundedFrames / std::max(1, walkFrames);
                            const bool pass = walkPath >= walkTargetM && pct >= 90.0 && walkMaxDensity < ISO_LEVEL;
                            printf("bench: planet walk %.0f m of %.0f m (path) in %.1f s, grounded %.0f%%, deepest body density %.2f (solid >= %.2f): %s\n",
                                   walkPath, walkTargetM, GetTime() - walkStart, pct, walkMaxDensity, ISO_LEVEL, pass ? "PASS" : "FAIL");
                            walkTargetM = 0.0;
                        }
                        benchWalkSeconds = 0.0;
                    }
                } else if ((WorldIdle() && benchDigStrokes == 0 && !benchPaint && !benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && minePhase != 2 && benchFlySeconds <= 0.0 && !benchSuffocate && benchOrbitSeconds <= 0.0 && !benchRide && !benchPark && !benchGravity && !benchPrecision && !benchMantle && !benchAirJump && !benchMoon && !benchHole && !benchWall && !benchLedge && !benchTor && !benchMachines && !benchFall && !benchTeleClose && !benchMenus && !benchWater && !benchFlyFollow && !benchHotbar && !benchInventory && !benchMap && !benchViewmodel && !benchQuests && !benchHunger && !benchFarm && !benchAnimals && !benchThreats && benchShotTabs == 0 && !benchForest && !benchFlow && !benchPond && !benchBucket && !benchPump && !benchTerrain && !benchShadow && !benchParticles && !benchBiomes && !benchOil && !benchSpring && !benchShotAnimals && !benchShotMachines && benchShotBiomes < 0 && !benchFarForest && !benchRecipes && !benchWeather && !benchDecor && !benchMaterials && !benchDamage && !benchCrafting && !benchSuit && !benchTools && !benchRivers && !benchRing && !benchRocket && !benchStation && !benchBlueprint && !benchSmallBody && !benchTrees && !benchChop && benchPersist == 0 && benchBuild == 0 && benchFactory == 0 && !benchProgression) || GetTime() > 600.0) {
                    benchStart = GetTime();
                    printf("bench: world ready after %.2f s (%d chunks, avg build %.2f ms of which sampling %.2f ms, %d worker threads)\n",
                        GetWorldStats().lastGenerationSeconds, GetWorldStats().chunks, GetWorldStats().avgBuildMs, GetWorldStats().avgSampleMs, GetWorldStats().workerThreads);
                }
            }
            double t = (benchStart < 0.0) ? 0.0 : GetTime() - benchStart;
            // #127: how far the camera moves against the planet over the timed window (a still bench view: ~0).
            static DVec3 benchCamLocal0{};
            static bool benchCamLocalSet = false;
            Planetoid* benchBody0 = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
            if (benchStart >= 0.0 && t > benchWarmup) {
                benchFrames++;
                if (!benchCamLocalSet && benchBody0) { benchCamLocal0 = benchBody0->toLocalD(world.renderOrigin + world.camera.position); benchCamLocalSet = true; }
                ++benchTimedFrames;
                benchMeasuring = true;
                if (GetFrameTime() > benchMaxFrame) benchMaxFrame = GetFrameTime();
                benchFrameTimes.push_back(GetFrameTime());
            }
            if ((benchFlySpeed > 0.0f || benchFlyAlt < 1e8f) && benchStart >= 0.0) {
                // --bench-fly V: the camera flies a great circle over body 0 at V m/s, 150 m above the
                // terrain, looking ahead and down; node counts each second show whether LOD keeps up.
                Planetoid* b = dynamic_cast<Planetoid*>(world.objects[0].get());
                const auto* pf = b ? dynamic_cast<const PlanetField*>(b->field.get()) : nullptr;
                if (pf) {
                    const double half = benchWarmup + benchSeconds * 0.5;
                    const bool back = benchFlyBack && t > half;
                    const double ang = benchFlySpeed * (back ? 2.0 * half - t : t) / (double)pf->radius;
                    const Vector3 dir = { (float)cos(ang) * 0.8f, 0.6f, (float)sin(ang) * 0.8f };       // body-local, 37 deg north
                    const Vector3 ahead = Vector3Scale(Vector3{ -(float)sin(ang) * 0.8f, 0.0f, (float)cos(ang) * 0.8f }, back ? -1.0f : 1.0f);
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
                        if (getenv("FLY_SHOTS") ? lastSec >= 3 : (lastSec == 8 || lastSec == 10 || lastSec == 12)) TakeScreenshot(TextFormat("bench_fly_%02d.png", lastSec)); // mid-flight views
                        const WorldStats& s = GetWorldStats();
                        printf("bench: fly t=%2d s  %5.1f km flown  nodes %5d  drawn %4d  queued %4d  %.0f FPS\n", (int)t,
                               benchFlySpeed * t / 1000.0, s.nodes, s.chunks, s.queuedJobs, (double)GetFPS());
                    }
                    // --crack-check while flying: every 10th frame count the magenta (cleared, nothing drawn) pixels below the
                    // horizon band; a frame with holes is logged and the worst saved (streaming holes are brief).
                    static long flyFrame = 0;
                    static int worstHoles = 0;
                    if (crackCheck && (++flyFrame % 10) == 0) {
                        Image shot = LoadImageFromScreen();
                        Color* px = LoadImageColors(shot);
                        int holes = 0;
                        for (int y = shot.height / 3; y < shot.height; y += 2)
                            for (int x = 0; x < shot.width; x += 2) { const Color c = px[y * shot.width + x]; holes += (c.r > 200 && c.g < 60 && c.b > 200) ? 1 : 0; }
                        UnloadImageColors(px);
                        if (holes > 50) {
                            printf("bench: fly holes: t=%.2f s, %d magenta samples\n", t, holes);
                            if (holes > worstHoles) { worstHoles = holes; ExportImage(shot, "bench_fly_holes.png"); }
                        }
                        UnloadImage(shot);
                    }
                }
            }
            if (t > benchWarmup + benchSeconds) {
                double avgMs = 1000.0 * benchSeconds / (double)benchFrames;
                const char* line = TextFormat("bench: %ld frames in %.1fs | avg %.1f FPS | avg %.2f ms | worst %.2f ms | world update %.2f ms, draw %.2f ms (CPU)\n",
                    benchFrames, benchSeconds, 1000.0 / avgMs, avgMs, benchMaxFrame * 1000.0, benchUpdateMs / std::max(1L, benchTimedFrames), benchDrawMs / std::max(1L, benchTimedFrames));
                printf("%s", line);
                SaveFileText("bench.txt", (char*)line);
                if (benchCamLocalSet && benchBody0) printf("bench: camera drift against the planet over the timed window: %.2f m (#127)\n", (benchBody0->toLocalD(world.renderOrigin + world.camera.position) - benchCamLocal0).length());
                GpuTimer_Print();
                GpuTimer_PrintCpu();
                if (TerrainArena_Available()) { const ArenaStats as = TerrainArena_Stats(); printf("bench: terrain arena: %d / %d vertices, %d / %d indices, %.1f chunk draws per multi-draw; uploads %lld staged, %lld direct\n", as.vertices, as.vertexCapacity, as.indices, as.indexCapacity, as.draws / (double)std::max(1, as.flushes), as.staged, as.direct); }
                { const TerrainGpuStats ts = TerrainGpu_Stats(); const WorldStats& ws = GetWorldStats(); printf("bench: terrain gpu: %lld blocks in %lld batches, avg batch latency %.2f ms; chunk build avg %.2f ms (sampling part %.2f ms)\n", ts.blocks, ts.batches, ts.avgLatencyMs, ws.avgBuildMs, ws.avgSampleMs); }
                {   // Streaming (#258): how much work flying asked for and threw away; frame time spread.
                    const WorldStats& ws = GetWorldStats();
                    std::vector<float> ft = benchFrameTimes;
                    std::sort(ft.begin(), ft.end());
                    auto pct = [&](double q) { return ft.empty() ? 0.0f : ft[std::min(ft.size() - 1, (size_t)(q * ft.size()))] * 1000.0f; };
                    int over20 = 0; for (float f : ft) over20 += f > 0.020f ? 1 : 0;
                    printf("bench: streaming: %lld requests, %lld built, %lld skipped stale, %lld results dropped; %d nodes alive; frame ms p50 %.2f p95 %.2f p99 %.2f, %d frames over 20 ms\n",
                           ws.requests, ws.built, ws.skippedStale, ws.droppedResults, ws.nodes, pct(0.5), pct(0.95), pct(0.99), over20);
                    printf("bench: preload: %d ring nodes, %d built; cache: %d nodes, %lld evicted; arena %.1f M vertices\n",
                           ws.preloadNodes, ws.preloadReady, ws.cachedNodes, ws.evicted, TerrainArena_Stats().vertices / 1048576.0);
                }
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
            if (cursorLocked) { EnableCursor(); cursorLocked = false; }
            BeginDrawing();
            EndDrawing();
            continue;
        }
        {
            // The cursor (#204, #205): locked while playing (DisableCursor: hidden, raw relative motion, no
            // warping) and free in menus and UIs, switched only when that changes. It used to be hidden and warped
            // back to the centre with SetMousePosition every frame: a blocking OS call per frame (benches skip it,
            // so they never showed the frame-rate drop while moving and looking) that also lost the motion between
            // warps (choppy look at low frame rates).
            const bool wantLock = !bench && gameState == 1 && !(DebugUI_IsOpen() || GameUI_InventoryOpen() || Nav_MapOpen());
            if (wantLock != cursorLocked) {
                if (wantLock) DisableCursor(); else EnableCursor();
                cursorLocked = wantLock;
            }
        }
        if (!bench && IsKeyPressed(KEY_F1)) {
            DebugUI_Toggle();
            // (relocking the cursor centres it)
        }
        if (!bench && gameState == 1 && IsKeyPressed(KEY_TAB)) {
            GameUI_ToggleInventory();
            // (relocking the cursor centres it)
        }
        // M: the map (aboard a rocket M is its orbit map).
        if (!bench && gameState == 1 && IsKeyPressed(KEY_M) && !Rocket_Boarded() && !GameUI_InventoryOpen() && !DebugUI_IsOpen()) {
            Nav_ToggleMap();
            // (relocking the cursor centres it)
        }
        const bool uiOpen = DebugUI_IsOpen() || GameUI_InventoryOpen() || Nav_MapOpen();
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
        const bool gameInput = !bench && !uiOpen && gameState == 1 && !loading;
        world.cameraInputEnabled = player.flyMode && !uiOpen && !ship.piloted;  // free-fly camera only in fly mode
        const bool alive = !survival.dead;
        // One clamped frame time for everything that moves (bodies, ship, player): clamping only the
        // player let the ground spin away under it during a slow frame.
        const float frameDt = std::min(GetFrameTime(), 0.1f);
        // A spawn or teleport: the loading screen until the terrain there is built (play and tests pause).
        if (player.spawns != seenSpawns) { seenSpawns = player.spawns; if (gameState == 1) { loading = true; spawnLoad = true; loadStart = GetTime(); loadMax = 1; loadFrames = 0; } }
        GpuTimer_CpuMark("bench-logic");
        // Benches that place the camera themselves (shots, tours): take that as the pose the fly-follow carries, or the
        // carry in player.update put the camera back where it was last frame (#298: since 4720ba4 the biome and machine
        // shots all showed the start, and --bench-fly never moved).
        if (bench && player.flyMode && player.flyFollow) player.recordFly(world);
        advanceBodies(world, gameState == 1 && !loading ? frameDt * timeScale * Rocket_Warp() : 0.0f); // bodies only spin while playing (warped with a coasting rocket); before riders update
        if (gameState == 1) {
            ship.lastImpact = 0.0f;
            Station_Update(world); // stations on their orbits, before anything rides on them
            if (ship.piloted) {
                ship.update(world, gameInput && alive, frameDt);
                ship.applyChaseCamera(world);
                player.position = ship.position; // the pilot rides inside
                player.velocity = ship.velocity;
            } else if (Rocket_Boarded()) {
                ship.update(world, false, frameDt);
                Rocket_Update(world, gameInput && alive, frameDt * timeScale * Rocket_Warp()); // flight (and its camera)
            } else {
                ship.update(world, false, frameDt);
                Suit_Apply(world, frameDt); // the worn suit and modules: walk, jetpack, grip, the scanner
                if (alive && !loading) player.update(world, gameInput, frameDt);
                {
                    // Dead (until the respawn): the view stays exactly where the player died, held in the frame of the
                    // body they were on - player.update is skipped, so nothing else carried it with the spinning,
                    // orbiting planet (it drifted off at ~100 m/s, or kept the death velocity; #231).
                    static uint64_t deathId = 0; // re-resolved each frame: the body may be destroyed while dead (ring rocks, regen)
                    static DVec3 deathLocal{};
                    static Vector3 deathForward{}, deathUp{};
                    static bool held = false;
                    if (!alive && !loading) {
                        if (!held) {
                            Planetoid* deathBody = player.nearest;
                            deathId = deathBody ? deathBody->id : 0;
                            if (deathBody) {
                                deathLocal = deathBody->toLocalD(player.position);
                                deathForward = deathBody->dirToLocal(player.forward);
                                deathUp = deathBody->dirToLocal(player.up);
                            }
                            held = true;
                        }
                        if (Planetoid* deathBody = deathId ? Planetoid::find(deathId) : nullptr) {
                            player.position = deathBody->toWorldD(deathLocal);
                            player.forward = deathBody->dirToWorld(deathForward);
                            player.up = deathBody->dirToWorld(deathUp);
                        }
                        player.velocity = { 0, 0, 0 };
                    } else {
                        held = false;
                        deathId = 0;
                    }
                }
                Survival_OnLanding(player.landingSpeed); // floor touchdowns only (#263, #255)
                Rocket_Update(world, gameInput && alive, frameDt); // on foot: aim at a rocket (E boards, I fuels)
            }
            if (ship.piloted && ship.lastImpact > 12.0f) survival.health -= (ship.lastImpact - 12.0f) * 3.0f; // crash damage
            if (!loading) Survival_Update(world, frameDt, gameInput);
            if (!loading) Nav_Update(world, !bench && Nav_MapOpen()); // exploration, the map
            if (!loading) Viewmodel_Update(world, frameDt);
            if (!loading) Farm_Update(world, frameDt * timeScale * Rocket_Warp()); // crops grow on world time
            if (!loading) Animals_Update(world, frameDt, frameDt * timeScale * Rocket_Warp(), gameInput && alive && !ship.piloted);
            if (!loading) Weather_Update(world, frameDt * timeScale); // the sky, rain and snow, lightning
            if (!loading) Ambience_Update(world, frameDt); // rain splashes, fireflies, machine smoke
            Particles_Update(world, frameDt); // smoke, dust, exhaust, splashes (particles.h)
            // Autosave every 5 minutes of play (incremental: only regions changed since the last save).
            static double lastAutosave = GetTime();
            if (!bench && alive && GetTime() - lastAutosave > 300.0) {
                lastAutosave = GetTime();
                if (SaveWorld(world)) GameUI_Toast("Autosaved");
            }
            Factory::Sync(world);                                // machine layouts in, items out
            Fluid::Update(world, frameDt);                       // flowing water: commands in, a step job, meshes up
            Waterworks_Update(world, frameDt);                   // pumps send water down pipes to outlets
            Machine_Update(world, gameInput && alive && !ship.piloted);
            Objectives_Update(world);
            if (bench && lookDown && !ship.piloted) {
                const Vector3 t = Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)));
                player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(t), 0.8f), player.up)); // ~50 deg down
            }
            if (bench && benchTelescope > 0.0f && !ship.piloted) AimAtTelescopeTarget(world);
            if (!ship.piloted && !Rocket_Boarded()) player.applyToCamera(world);
            Telescope_Update(world, gameInput && alive && !Rocket_Boarded(), !ship.piloted && !player.flyMode && alive && !buildTool.active); // aboard, T moves fuel
            if (bench && benchTelescope > 0.0f) { Telescope_SetActive(world, true); telescope.zoom = benchTelescope; Telescope_Update(world, false, !ship.piloted); }
        }

        GpuTimer_CpuMark("game-update");
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
                            world.atmosphere.toLocalX = q->dirToLocal({ 1, 0, 0 });
                            world.atmosphere.toLocalY = q->dirToLocal({ 0, 1, 0 });
                            world.atmosphere.toLocalZ = q->dirToLocal({ 0, 0, 1 });
                        }
                    }
        }
        if (crackCheck) world.drawSkybox = false;
        if (crackCheck) world.clearColor = MAGENTA;
        while (world.lights.size() < 3) world.lights.push_back(Light{});
        {
            Light& lamp = world.lights[1];
            lamp.type = LIGHT_SPOT;
            lamp.enabled = gameState == 1 && survival.headlamp && !survival.dead;
            if (lamp.enabled) {
                const float range = ship.piloted ? 60.0f : 35.0f * Suit_Stats().lamp; // the lamp upgrade reaches further
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
        GpuTimer_CpuMark("world+lights");
        if (gameState == 1) Hotbar_Update(world, gameInput && alive, ship.piloted);
        if (gameState == 1) BuildTool_Update(world, gameInput && alive, !ship.piloted && !Rocket_Boarded() && !player.flyMode && !telescope.active && alive);
        if (gameState == 1 && !loading) Animals_Aim(world, gameInput && alive && !ship.piloted); // this frame's camera (#222)
        if (gameState == 1) TerrainTool_Update(world, gameInput && alive && !Animals_Aiming(), !ship.piloted && !Rocket_Boarded() && !telescope.active && !buildTool.active); // F5/F9 work aboard too
        {   // Tool wear: each stroke, chop or hit with the held tool (the swing counter covers them all).
            static int lastSwings = 0;
            if (terrainTool.swings > lastSwings && !terrainTool.creative && !terrainTool.terraform && !terrainTool.heldTool.empty())
                ToolWear_OnStroke(terrainTool.heldTool, terrainTool.swings - lastSwings);
            lastSwings = terrainTool.swings;
        }

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
                printf("bench: menus: Save & Quit -> main menu: %s; old world unloaded (%zu bodies left, #270): %s\n", gameState == 0 ? "yes" : "NO", world.objects.size(), world.objects.empty() ? "yes" : "NO");
                menuOk &= world.objects.empty();
                pendingMenu.action = MenuAction::LoadGame; pendingMenu.path = "saves/bench.sav"; menuStep = 5; menuFrames = 0;
            } else if (menuStep == 5 && menuFrames == 5) {
                menuOk &= gameState == 1;
                printf("bench: menus: Load Game -> in game: %s -> %s\n", gameState == 1 ? "yes" : "NO", menuOk ? "PASS" : "FAIL");
                benchMenus = false;
            }
        }
        // Menus: the result of last frame's clicks (the menus draw, and see clicks, in BeginDrawing).
        {
            auto enterGame = [&](bool fresh = true) {
                HideCursor();
                gameState = 1;
                mainMenu.isActive = false;
                world.isActive = true;
                loading = !bench; loadStart = GetTime(); loadMax = 1; // the loading screen until the terrain around you is in
                if (fresh && !bench) BeginSpawnPreload(); // a new or loaded world: the ring around the spawn first (#294)
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
                case MenuAction::Resume: enterGame(false); loading = false; break;
                case MenuAction::Save: GameUI_Toast(SaveWorld(world) ? "Saved" : "Save failed"); break;
                case MenuAction::SaveAndQuit:
                    if (!SaveWorld(world)) GameUI_Toast("Save failed");
                    // Nothing of the old world keeps running behind the menu (#270): jobs, the factory, water,
                    // streaming, bodies and their GPU meshes, animals, particles. New Game and Load build a world again.
                    if (telescope.active) Telescope_SetActive(world, false);
                    player.flyMode = false;
                    unloadWorld(world);
                    Animals_Reset();
                    Particles_Clear();
                    worldFresh = false;
                    toMenu();
                    break;
                case MenuAction::Exit: windowKeepAlive = false; break;
                case MenuAction::None: break;
            }
            // Esc in game pauses (unless it closes the inventory or build mode first); the pause menu
            // handles Esc itself (resume, or back out of Settings).
            if (gameState == 1 && m.action == MenuAction::None && IsKeyPressed(KEY_ESCAPE)) {
                if (Nav_MapOpen()) Nav_CloseMap();
                else if (GameUI_InventoryOpen()) GameUI_ToggleInventory();
                else if (buildTool.active) buildTool.active = false;
                else { ShowCursor(); Menu_ResetPause(); gameState = 2; justPaused = true; }
            }
        }

        

        GpuTimer_Mark("pre-draw (uploads, updates)");
        BeginDrawing();
            ClearBackground(crackCheck ? MAGENTA : BLACK); // --crack-check: holes in the terrain show magenta
            GpuTimer_Enable(benchMeasuring || perfLog); GpuTimer_Mark("clear");
            GpuTimer_CpuMark("begin-draw");
            if (gameState == 0) {
                pendingMenu = Menu_Main();
            } else if (gameState == 1) {
                { const double d0 = GetTime(); world.drawScene(gameState); if (benchMeasuring) benchDrawMs += (GetTime() - d0) * 1000.0; }
                player.recordFly(world); // after the fly camera's input (in drawScene): the pose to carry next frame
                Viewmodel_Draw(world);   // the hands, over the world
                world.drawUI(gameState);
                Telescope_DrawHUD(world);
                BuildTool_DrawHUD();
                Hotbar_DrawHUD(ship.piloted);
                Machine_DrawHUD();
                if (!ship.piloted && !buildTool.active) Animals_DrawHUD();
                Suit_DrawHUD();
                Weather_DrawHUD();
                if (!ship.piloted && !telescope.active && !buildTool.active && !Rocket_Boarded()) TerrainTool_DrawHUD();
                GameUI_DrawHUD(world);
                Survival_DrawHUD();
                Objectives_DrawHUD();
                if (!Rocket_Boarded()) Nav_DrawHUD(world);
                Rocket_DrawHUD(world); // the orbit map covers the rest
                Nav_DrawMap(world);    // last: the full map covers everything
                if (loading) {
                    // Loading screen: black, with the streaming progress, until the world around you is built
                    // (it used to show you standing in a void among trees and water while the ground popped in).
                    const size_t pending = WorldJobs().queued() + WorldJobs().running();
                    loadMax = std::max(loadMax, pending);
                    float ringProgress = 1.0f;
                    const bool preparing = SpawnPreloadPending(&ringProgress); // first load: the 300 m ring around the spawn (#294)
                    const float progress = preparing ? ringProgress : WorldIdle() ? 1.0f : 1.0f - (float)pending / (float)loadMax;
                    DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), BLACK);
                    const int bw = GetScreenWidth() / 3, bx = (GetScreenWidth() - bw) / 2, by = GetScreenHeight() / 2;
                    DrawText(preparing ? "Preparing spawn" : "Building the world", bx, by - 40, 28, RAYWHITE);
                    DrawRectangle(bx, by, bw, 14, Fade(RAYWHITE, 0.15f));
                    DrawRectangle(bx, by, (int)(bw * std::clamp(progress, 0.02f, 1.0f)), 14, SKYBLUE);
                    // Entering the game holds it half a second; a spawn only until the streamer has had a few
                    // frames to queue the new place and has built it (where it is built already: at once).
                    ++loadFrames;
                    const bool built = !preparing && (spawnLoad ? WorldIdle() && loadFrames >= 3 : WorldIdle() && GetTime() - loadStart > 0.5);
                    if (built || GetTime() - loadStart > (preparing ? 90.0 : 20.0)) {
                        if (bench) printf("bench: loading screen after a spawn: %d frames, %.2f s\n", loadFrames, GetTime() - loadStart);
                        loading = false; spawnLoad = false;
                    }
                }
                // (No per-frame cursor warp: the cursor is locked with DisableCursor while playing, #204 / #205.)
                if (perfLog) {
                    // --perf-log (#204): every 5 s, the per-pass averages to log.txt, tagged by how the player moved (average
                    // speed, and mouse-look turning) so "standing" and "moving" can be compared from a real session.
                    static double perfT0 = GetTime(), perfSpeed = 0.0, perfLook = 0.0; static long perfN = 0;
                    perfSpeed += Vector3Length(Vector3Subtract(player.velocity, !player.grounded && player.nearest ? player.nearest->surfaceVelocity((player.position - world.renderOrigin).toFloat()) : Vector3{ 0, 0, 0 })); /* against the ground: airborne velocity is world-frame and includes the planet's spin (grounded it is surface-relative already) */ perfLook += Vector2Length(GetMouseDelta()); ++perfN;
                    if (GetTime() - perfT0 >= 5.0) {
                        const double sp = perfSpeed / std::max(1L, perfN), lk = perfLook / std::max(1L, perfN);
                        GpuTimer_LogAndReset(TextFormat("%s (speed %.1f m/s, look %.1f px/frame, %d FPS)", sp > 0.5 || lk > 0.5 ? "MOVING" : "still", sp, lk, GetFPS()));
                        perfT0 = GetTime(); perfSpeed = perfLook = 0.0; perfN = 0;
                    }
                }
            } else if (gameState == 2) {
                world.drawScene(gameState);
                world.drawUI(gameState);
                const MenuResult m = Menu_Pause();
                if (!justPaused) pendingMenu = m; // the Esc that opened it is not also a "resume"
            }
            justPaused = false;
            if (gameState != 0) {
                GpuTimer_Mark("hud");
                rlImGuiBegin();
                DebugUI_Draw(world);
                GameUI_Draw(world);
                rlImGuiEnd();
                GpuTimer_Mark("imgui");
            }
            rlDrawRenderBatchActive(); // the batched 2D (HUD text, bars, minimap) flushes here, not inside EndDrawing
            GpuTimer_Mark("2d-flush");
        EndDrawing();
        GpuTimer_Mark("swap+poll");
        PROFILE_FRAME();
    }

    }   // scenes destroyed here, while the GL context still exists
    DebugUI_Shutdown();
    Audio_Shutdown();
    Factory::Stop();
    CloseWindow();
}


