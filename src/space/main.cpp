#include <set>
#include <map>
#include <cstring>
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
#include "gputimer.h"
#include "ship.h"
#include "camera.h"
#include "survival.h"
#include "telescope.h"
#include "buildtool.h"
#include "progression_test.h"
#include "menus.h"
#include "hotbar.h"
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

// A new game on the world as generated: the player at the start, the ship parked ahead, a
// clearing for both, the starting kit. `reset` also clears what a previous game left behind
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
    double benchSeconds = 0.0;
    bool startWithOverlay = false;
    int benchDigStrokes = 0;     // --test-dig N: dig N strokes at the crosshair once the world is loaded
    bool benchPaint = false;     // --test-paint: placed gravel stays gravel (and through a save and load)
    bool benchSaveLoad = false;  // --test-saveload: then save, reload from disk, and bench the reloaded world
    double benchWalkSeconds = 0.0; // --test-walk S: spawn on the surface and auto-walk forward for S seconds
    bool benchWalkMode = false;    // --bench-walk: the player walks forward for the whole bench (frame time while streaming)
    float benchFlySpeed = 0.0f;    // --bench-fly V: the camera flies over body 0 at V m/s through the bench (LOD streaming)
    float benchFlyAlt = 1e9f;      // --bench-fly-alt A: fly at A m relative to sea level instead (negative: underwater, looking up)
    bool hasSpawnDir = false;      // --spawn-dir x y z: where --test-walk / --test-stand spawn on body 0
    Vector3 spawnDir = { 0, 1, 0 };
    double walkStart = -1.0;
    bool benchStandOnly = false; // --test-stand S: like --test-walk but standing still
    int walkJolts = 0, walkJoltFrames = 0, walkSteepFrames = 0, walkSteepJolts = 0; double walkJoltMax = 0.0, walkSlopeMax = 0.0;
    double standCamMax = 0.0; float standFwdMax = 0.0f, standUpMax = 0.0f, standTexelSlide = 0.0f; int standFrames = 0, standJitterFrames = 0;
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
    bool benchInventory = false;   // --test-inventory: grid stacking, fractions, sorting, saves
    bool benchMap = false;         // --test-map: compass, exploration, waypoints, map images, saves
    bool benchQuests = false;      // --test-quests: tiers open by their tasks; quest state through a save
    bool benchHunger = false, benchFarm = false; // --test-hunger, --test-farm
    bool benchRivers = false, riverOk = false; int riverPhase = 0; double riverStart = 0.0; float riverLongest = 0.0f, riverBuild = 0.0f, riverWet = 0.0f, riverLevel = 0.0f; Vector3 riverDir{ 0, 0, 0 }; // --test-rivers
    bool benchTools = false, toolsOk = true; int toolsPhase = 0, toolsLife = 0, toolsTilled = 0; double toolsStart = 0.0; // --test-tools
    bool benchSuit = false, suitOk = true; int suitPhase = 0; double suitStart = 0.0; float suitFoodRatio = 0.0f, suitWater = 0.0f, suitRebreather = 0.0f; // --test-suit
    int benchShotTabs = 0; // --shot-tabs
    bool benchTerrain = false; int terrainPhase = 0; double terrainT0 = 0.0; Vector3 terrainAt{}, terrainOut{}; // --shot-terrain
    bool benchPump = false; int pumpPhase = 0, pumpCable = -1; double pumpT0 = 0.0, pumpedAtCut = 0.0; uint32_t pumpPlot = 0; Vector3 pumpBasin{}; float pumpRise = 0.0f; // --test-pump
    bool benchBucket = false; int bucketPhase = 0, bucketTrips = 0; double bucketT0 = 0.0; Vector3 bucketWater{}, bucketBasin{}, bucketField{}; // --test-bucket
    bool benchPond = false; int pondPhase = 0; double pondT0 = 0.0; Vector3 pondDirB{}, pondOutB{}, pondPit{}, pondMid{}; float pondLevelB = 0.0f, pondRB = 0.0f; // --test-pond
    bool benchFlow = false; int flowPhase = 0; double flowStart = 0.0, flowEdit = 0.0, flowReachT = -1.0, flowCaveT = -1.0, flowSleepT = -1.0; Vector3 flowDir{}, flowOut{}, flowEnd{}, flowCave{}; float flowLevel = 0.0f; // --test-flow
    bool benchForest = false; int forestPhase = 0; double forestStart = 0.0; Vector3 forestCam{}, forestTarget{}, forestUp{}; // --shot-forest
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
        if (std::string(argv[i]) == "--test-paint") benchPaint = true;
        if (std::string(argv[i]) == "--test-saveload") benchSaveLoad = true;
        if (std::string(argv[i]) == "--test-walk" && i + 1 < argc) benchWalkSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--no-eye-smooth") player.eyeSmoothing = false;
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
        if (std::string(argv[i]) == "--test-flow") benchFlow = true;
        if (std::string(argv[i]) == "--test-pond") benchPond = true;
        if (std::string(argv[i]) == "--test-bucket") benchBucket = true;
        if (std::string(argv[i]) == "--test-pump") benchPump = true;
        if (std::string(argv[i]) == "--shot-terrain") benchTerrain = true;
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
    if (bench && !benchThreats) Animals_SetThreats(false); // hostile animals and night spawns would upset other tests

    bool msaa = true;
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--no-msaa") msaa = false;
    if (msaa) SetConfigFlags(FLAG_MSAA_4X_HINT);
    bool nativeWindow = false; // --native-window: a bench in the game's own window (borderless, monitor-sized), uncapped
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--native-window") nativeWindow = true;
    if (!bench || nativeWindow) {
        SetConfigFlags(FLAG_WINDOW_RESIZABLE); // Allow window resizing
        SetConfigFlags(FLAG_WINDOW_UNDECORATED); // Remove window decorations for a cleaner look
        SetConfigFlags(FLAG_WINDOW_HIGHDPI); // Enable High DPI support
        InitWindow(GetMonitorWidth(0), GetMonitorHeight(0), "My first RAYLIB program!");
        // No VSync: cap at the monitor's refresh rate instead (VSync on a 60 Hz desktop locked the game to 60).
        SetTargetFPS(bench ? std::max(0, benchFps) : benchFps >= 0 ? benchFps : std::max(60, GetMonitorRefreshRate(GetCurrentMonitor())));
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
    world.onDraw3D = [&world] { Farm_Draw3D(world); Animals_Draw3D(world); Weather_Draw3D(world); if (!ship.piloted && !buildTool.active && !Rocket_Boarded()) TerrainTool_Draw3D(); BuildTool_Draw3D(world); };
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
    long benchFrames = 0;
    if (bench) {
        gameState = 1;
        mainMenu.isActive = false;
        world.isActive = true;
        if (!benchAtSave) world.camera.position = Vector3Add(world.camera.target, benchCamOffset);
        world.debugDraw = startWithOverlay;
        player.flyMode = benchWalkSeconds <= 0.0 && !benchWalkMode; // the bench drives the camera directly unless walking
        player.autoWalk = benchWalkMode;                              // --bench-walk: walk through the measured window
        player.flyFollow = false;                                     // the bench drives the camera itself
        benchStart = -1.0; // starts once the world has finished streaming in
    }

    while (windowKeepAlive && !WindowShouldClose())
    {
        GpuTimer_CpuFrameStart(); GpuTimer_BeginFrame();
        if (bench && !loading) { // scripted tests wait out the loading screen (a half-built world is no game state)
            if (benchStart < 0.0) {
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
                    } else if (survival.deaths > 0 && !survival.dead) {
                        const float d = Vector3Length(Vector3Subtract(player.velocity, player.nearest ? player.nearest->surfaceVelocity((player.position - world.renderOrigin).toFloat()) : Vector3{ 0, 0, 0 }));
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
                        printf("bench: wall: %s\n", wallOk && pass ? "PASS" : "FAIL");
                        benchWall = false;
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
                        if (t > 4.0) {
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
                        if (t < 0.5) { h0 = h; inWaterFrames = frames = 0; }
                        else { ++frames; inWaterFrames += player.inWater ? 1 : 0; }
                        if (t > 3.0) {
                            const float sink = (h0 - h) / 2.5f;
                            const bool pass = frames > 0 && inWaterFrames * 10 >= frames * 9 && sink < 1.0f;
                            printf("bench: pond: swimming in the pit: in water %d%% of frames, sank %.2f m/s (falling in air would be ~%.0f m/s) -> %s\n",
                                   frames ? 100 * inWaterFrames / frames : 0, sink, player.gravity * 2.5f * 0.5f, pass ? "PASS" : "FAIL");
                            benchPond = false;
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
                            float bestD = 1e9f; Vector3 bestP{};
                            for (int i = -120; pf && i <= 120; ++i)
                                for (int j = -120; j <= 120; ++j) {
                                    const float dd = 25.0f * sqrtf((float)(i * i + j * j));
                                    if (dd >= bestD) continue;
                                    const Vector3 dir = Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, 25.0f * i / pf->radius), Vector3Scale(t2, 25.0f * j / pf->radius))));
                                    const float surf = pf->radius + pf->elevation(dir);
                                    const Vector3 q = Vector3Scale(dir, surf + 2.0f);
                                    if (pf->rockFeature(q, surf + 2.0f, surf, 1.0f) > 0.5f) { bestD = dd; bestP = q; }
                                }
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
                        Animals_RefreshNow();
                        Animals_Update(world, 0.1f, 0.1f, false);
                        threatSpawned = Animals_NightCount();
                        const float d0 = nearestNight();
                        check(threatSpawned == 1 && d0 > 40.0f && d0 < 65.0f, "a stalker comes out of the dark 45-60 m away");
                        step(15.0f);
                        threatBite = 100.0f - survival.health;
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
                        auto onGround = [&](int i) { const Vector3 a = Animals_Position(i); return fabsf(Vector3Length(a) - (pf->radius + pf->elevation(Vector3Normalize(a)))) < 0.05f; };
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
                            Animals_Spawn(world, *p, "woolback", Vector3Add(me, Vector3Scale(east, 12.0f)), 1, false, 0.5f);
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
                        if (next[vmPhase - 4][0]) { inventory.add(next[vmPhase - 4], 1); terrainTool.heldTool = next[vmPhase - 4]; }
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
                        Hotbar_Press(0, false); // key 1: an empty slot, bare hands
                        const bool toolMode = !buildTool.active;
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
                    if (GetTime() - walkStart > benchWalkSeconds || (walkTargetM > 0.0 && walkPath >= walkTargetM)) {
                        player.autoWalk = false;
                        const double dist = (start->toLocalD(player.position) - walkFrom).length();
                        const float alt = Vector3Distance(feet, start->position);
                        printf("bench: walked %.1f m in %.1f s, grounded %.0f%% of frames, density at body %.2f..%.2f (solid >= %.2f), %.1f m from centre\n",
                               dist, GetTime() - walkStart, 100.0 * walkGroundedFrames / std::max(1, walkFrames), walkMinDensity, walkMaxDensity, ISO_LEVEL, alt);
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
                } else if ((WorldIdle() && benchDigStrokes == 0 && !benchPaint && !benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && minePhase != 2 && benchFlySeconds <= 0.0 && !benchSuffocate && benchOrbitSeconds <= 0.0 && !benchRide && !benchPark && !benchGravity && !benchPrecision && !benchMantle && !benchAirJump && !benchMoon && !benchHole && !benchWall && !benchMenus && !benchWater && !benchFlyFollow && !benchHotbar && !benchInventory && !benchMap && !benchViewmodel && !benchQuests && !benchHunger && !benchFarm && !benchAnimals && !benchThreats && benchShotTabs == 0 && !benchForest && !benchFlow && !benchPond && !benchBucket && !benchPump && !benchTerrain && !benchRecipes && !benchWeather && !benchDecor && !benchMaterials && !benchDamage && !benchCrafting && !benchSuit && !benchTools && !benchRivers && !benchRing && !benchRocket && !benchStation && !benchBlueprint && !benchSmallBody && !benchTrees && !benchChop && benchPersist == 0 && benchBuild == 0 && benchFactory == 0 && !benchProgression) || GetTime() > 600.0) {
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
                GpuTimer_Print();
                GpuTimer_PrintCpu();
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
            if (DebugUI_IsOpen() || GameUI_InventoryOpen() || Nav_MapOpen()) ShowCursor(); else HideCursor();
        }
        if (!bench && IsKeyPressed(KEY_F1)) {
            DebugUI_Toggle();
            if (!DebugUI_IsOpen()) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
        }
        if (!bench && gameState == 1 && IsKeyPressed(KEY_TAB)) {
            GameUI_ToggleInventory();
            if (!GameUI_InventoryOpen()) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
        }
        // M: the map (aboard a rocket M is its orbit map).
        if (!bench && gameState == 1 && IsKeyPressed(KEY_M) && !Rocket_Boarded() && !GameUI_InventoryOpen() && !DebugUI_IsOpen()) {
            Nav_ToggleMap();
            if (!Nav_MapOpen()) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
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
                Survival_OnLanding(player.impactSpeed);
                Rocket_Update(world, gameInput && alive, frameDt); // on foot: aim at a rocket (E boards, I fuels)
            }
            if (ship.piloted && ship.lastImpact > 12.0f) survival.health -= (ship.lastImpact - 12.0f) * 3.0f; // crash damage
            if (!loading) Survival_Update(world, frameDt, gameInput);
            if (!loading) Nav_Update(world, !bench && Nav_MapOpen()); // exploration, the map
            if (!loading) Viewmodel_Update(world, frameDt);
            if (!loading) Farm_Update(world, frameDt * timeScale * Rocket_Warp()); // crops grow on world time
            if (!loading) Animals_Update(world, frameDt, frameDt * timeScale * Rocket_Warp(), gameInput && alive && !ship.piloted);
            if (!loading) Weather_Update(world, frameDt * timeScale); // the sky, rain and snow, lightning
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
                loading = !bench; loadStart = GetTime(); loadMax = 1; // the loading screen until the terrain around you is in
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
                case MenuAction::Resume: enterGame(); loading = false; break;
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
                if (Nav_MapOpen()) { Nav_CloseMap(); SetMousePosition(GetScreenWidth() / 2, GetScreenHeight() / 2); }
                else if (GameUI_InventoryOpen()) { GameUI_ToggleInventory(); SetMousePosition(GetScreenWidth() / 2, GetScreenHeight() / 2); }
                else if (buildTool.active) buildTool.active = false;
                else { ShowCursor(); Menu_ResetPause(); gameState = 2; justPaused = true; }
            }
        }

        

        GpuTimer_Mark("pre-draw (uploads, updates)");
        BeginDrawing();
            ClearBackground(crackCheck ? MAGENTA : BLACK); // --crack-check: holes in the terrain show magenta
            GpuTimer_Enable(benchMeasuring); GpuTimer_Mark("clear");
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
                    const float progress = WorldIdle() ? 1.0f : 1.0f - (float)pending / (float)loadMax;
                    DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), BLACK);
                    const int bw = GetScreenWidth() / 3, bx = (GetScreenWidth() - bw) / 2, by = GetScreenHeight() / 2;
                    DrawText("Building the world", bx, by - 40, 28, RAYWHITE);
                    DrawRectangle(bx, by, bw, 14, Fade(RAYWHITE, 0.15f));
                    DrawRectangle(bx, by, (int)(bw * std::clamp(progress, 0.02f, 1.0f)), 14, SKYBLUE);
                    // Entering the game holds it half a second; a spawn only until the streamer has had a few
                    // frames to queue the new place and has built it (where it is built already: at once).
                    ++loadFrames;
                    const bool built = spawnLoad ? WorldIdle() && loadFrames >= 3 : WorldIdle() && GetTime() - loadStart > 0.5;
                    if (built || GetTime() - loadStart > 20.0) {
                        if (bench) printf("bench: loading screen after a spawn: %d frames, %.2f s\n", loadFrames, GetTime() - loadStart);
                        loading = false; spawnLoad = false;
                    }
                }
                if (!bench && !uiOpen) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
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


