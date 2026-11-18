#include <set>
#include <map>
#include <cstring>
#include <raylib.h>
#include <string>
#include <vector>
#include <random>
#include <chrono>
#include <thread>
#include <array>
#include "terrainarena.h"
#include "mainmenu.h"
#include "world.h"
#include "playerinput.h"
#include "playersession.h"
#include "netid.h"
#include "net.h"
#include "gamerng.h"
#include "headless.h"
#include "hostserver.h"
#include "joininfo.h"
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
#include "helpmanual.h"
#include "gamesave.h"
#include "audio.h"
#include "rlImGui.h"
#include "gputimer.h"
#include "ship.h"
#include "camera.h"
#include "survival.h"
#include "combat.h"
#include "telescope.h"
#include "glcaps.h"
#include "buildtool.h"
#include "recorder.h"
#include "slideshow.h"
#include "bugreport.h"
#include "progression_test.h"
#include "menus.h"
#include "newworld.h"
#include "seedsurvey.h"
#include "hotbar.h"
#include "itemicons.h"
#include "particles.h"
#include "ambience.h"
#include "waterfalls.h"
#include "gpumem.h"
#include "structures.h"
#include "fishing.h"
#include "trapping.h"
#include "events.h"
#include "machinemodels.h"
#include "gpucompute.h"
#include "terraingpu.h"
#include "erosion.h"
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
#include "wind.h"
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
#include "waterfield.h"
#include "waterworks.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif

using namespace std;
static bool gLogLevelFromArgs = false; // --log-level given: settings.txt's log_level does not apply (#514)


int gameState = 0; // 0: Main Menu, 1: Game, 2: Pause Menu, 3: New World screen (#504)

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

// RAYLIBGAME_SPAWN_DIR="x,y,z" (body-local) moves the start there, so a bench photographs a chosen place (the erosion
// compare dump prints its window's direction, #358).
static bool SpawnDirOverride(Vector3& d) {
    const char* e = getenv("RAYLIBGAME_SPAWN_DIR");
    if (!e || sscanf(e, "%f,%f,%f", &d.x, &d.y, &d.z) != 3 || Vector3Length(d) < 1e-6f) return false;
    d = Vector3Normalize(d);
    return true;
}

static bool gInFrame = false; // between the main loop's BeginDrawing and EndDrawing (#362: gWorldGenWait must not nest a frame)
// --perf-log (#387): once a second, the frames per screen state (menu, newworld, bake, loading, play, paused): count,
// average and worst frame, to stdout and log.txt. Called once per drawn frame.
static bool gFpsLog = false;
static void FpsLog_Frame(const char* state) {
    if (!gFpsLog) return;
    static const char* cur = nullptr;
    static double t0 = 0.0, last = 0.0, sum = 0.0, worst = 0.0;
    static int n = 0;
    const double now = GetTime();
    auto flush = [&] {
        if (cur && n > 0) {
            const char* line = TextFormat("[fps] %.0f s %s: %d frames in %.1f s (%.0f FPS), avg %.1f ms, worst %.1f ms\n", now, cur, n, now - t0, n / std::max(1e-3, now - t0), 1000.0 * sum / n, 1000.0 * worst);
            printf("%s", line);
            LOG_INFO("%s", line);
        }
        cur = state; t0 = now; sum = worst = 0.0; n = 0;
    };
    if (!cur || strcmp(cur, state) != 0 || now - t0 >= 1.0) flush();
    if (last > 0.0) { const double dt = now - last; sum += dt; worst = std::max(worst, dt); ++n; }
    last = now;
}


static Vector3 StartSpawnDirLocal(Planetoid& body) {
    if (Vector3 d; SpawnDirOverride(d)) return d;
    if (auto* field = dynamic_cast<const PlanetField*>(body.field.get())) return field->findSpawnDir(body.seed);
    return { 0, 1, 0 };
}

static Vector3 StartSpawnDir(Planetoid& body) {
    if (Vector3 d; SpawnDirOverride(d)) return body.dirToWorld(d);
    if (auto* field = dynamic_cast<const PlanetField*>(body.field.get())) return body.dirToWorld(field->findSpawnDir(body.seed));
    return { 0, 1, 0 };
}

// Tests with a second player or centre (#314, #318): a world direction about `dist` m from `from` over dry land (the first
// of 16 headings whose ground stands 5 m over the sea; else the first heading).
static Vector3 LandDirAt(Planetoid& body, Vector3 from, float dist) {
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    Vector3 e = Vector3CrossProduct(from, { 0, 1, 0 });
    if (Vector3Length(e) < 0.1f) e = { 1, 0, 0 };
    e = Vector3Normalize(e);
    Vector3 firstDir{};
    for (int k = 0; k < 16; ++k) {
        const Vector3 side = Vector3RotateByQuaternion(e, QuaternionFromAxisAngle(from, k * (2.0f * PI / 16.0f)));
        const Vector3 d = Vector3Normalize(Vector3Add(from, Vector3Scale(side, dist / (body.radius * body.scale))));
        if (k == 0) firstDir = d;
        if (!pf) return d;
        const Vector3 l = body.dirToLocal(d);
        if (pf->radius + pf->elevation(l) > pf->seaLevelRadius() + 5.0f) return d;
    }
    return firstDir;
}

// A new game on the world as generated: the player at the start, a
// clearing (no tree on the spot), the starting kit. `reset` also clears what a previous game left behind
// (inventory, quest tier, survival state).
// dawn (#259): a new world starts just before sunrise at the spawn (not for benches: they keep their fixed start).
static void StartNewWorld(Scene& world, bool reset, bool dawn) {
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
        Structures_Reset();
        Fishing_Reset();
        Trapping_Reset();
        Combat_Reset();
        Events_Reset();
        Particles_Clear();
        Ambience_Reset();
        Waterfalls_Reset();
        Station_Reset();
        Suit_Reset();
        ToolWear_Reset();
        RecipeBook_Reset();
    }
    Planetoid* start = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
    if (start) {
        // Dawn at the spawn (#259): the planet's spin phase set so the sun is 4 degrees under the spawn's horizon and
        // rising (~16 game minutes before sunrise), before the player is placed (the spawn direction is body-local).
        // The phase is solved numerically from the spin, so its sign and the rotation conventions don't matter.
        if (dawn) {
            const double t = SecondsToSunElevation(*start, StartSpawnDirLocal(*start), world.sunDirection, -4.0f, true);
            if (t >= 0.0) {
                const Vector3 s = start->spin;
                auto wrap = [](double a) { a = fmod(a, 360.0); return a < 0.0 ? a + 360.0 : a; };
                start->rotation = { (float)wrap(start->rotation.x + s.x * t), (float)wrap(start->rotation.y + s.y * t), (float)wrap(start->rotation.z + s.z * t) };
                world.syncRenderPositions();
            }
        }
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


extern const char* kGameVersion; // "<branch> <commit>", generated at build time by cmake/game_version.cmake (#501)

// --test-build (#602): every live part (plot, type, cell, orientation), to name what appears, or what a save/load adds or loses.
using PartCensusMap = std::map<std::tuple<uint32_t, int, int, int, int, int>, int>;
static PartCensusMap PartCensus(const Planetoid* b) {
    PartCensusMap m;
    if (b && b->plots) for (const auto& pl : b->plots->plots)
        for (const Part& q : pl->parts) if (q.type != PartType::Count) ++m[{ pl->id, (int)q.type, q.origin.x, q.origin.y, q.origin.z, q.orient }];
    return m;
}
static int PartCensusDiff(const PartCensusMap& was, const PartCensusMap& now, const char* what, int limit) {
    int n = 0;
    for (const auto& kv : now) {
        const auto it = was.find(kv.first);
        const int w = it == was.end() ? 0 : it->second;
        if (kv.second > w && n++ < limit)
            printf("bench: build: %s: plot %u, part type %d at %d,%d,%d orient %d\n", what, std::get<0>(kv.first), std::get<1>(kv.first),
                   std::get<2>(kv.first), std::get<3>(kv.first), std::get<4>(kv.first), std::get<5>(kv.first));
    }
    return n;
}
static PartCensusMap gBuildCensus;

int main(int argc, char** argv)
{
    // Scripted runs (--bench, --test-*, --shot-*): stdout unbuffered, before any output, so a run killed by a timeout
    // still shows every result line it reached (#556: the last lines were lost and a finished test looked hung).
    // (Line buffering is full buffering in the Windows CRT.)
    for (int i = 1; i < argc; ++i)
        if (!strncmp(argv[i], "--bench", 7) || !strncmp(argv[i], "--test-", 7) || !strncmp(argv[i], "--shot-", 7)) { setvbuf(stdout, nullptr, _IONBF, 0); break; }
    Log_Open("log.txt");
    Log_SetThreadName("main");
    // --log-level LEVEL (#514): before anything logs. It wins over settings.txt's log_level.
    for (int i = 1; i + 1 < argc; ++i)
        if (std::string(argv[i]) == "--log-level") {
            LogLevel l;
            if (Log_ParseLevel(argv[i + 1], l)) { Log_SetLevel(l); gLogLevelFromArgs = true; }
            else printf("--log-level: unknown level '%s' (trace, debug, info, warn, error, fatal, off)\n", argv[i + 1]);
        }
    Log_StartupChecks();
    // --to-slide <in.png> <out.jpg> [width]: a screenshot made into a menu slide (#386): scaled to the width (1600 px), JPEG.
    if (argc >= 4 && std::string(argv[1]) == "--to-slide") {
        Image img = LoadImage(argv[2]);
        if (!img.data) { printf("to-slide: cannot read %s\n", argv[2]); return 1; }
        const int w = argc >= 5 ? atoi(argv[4]) : 1600;
        if (img.width > w) ImageResize(&img, w, img.height * w / img.width);
        ImageFormat(&img, PIXELFORMAT_UNCOMPRESSED_R8G8B8);
        const bool ok = ExportImage(img, argv[3]);
        printf("to-slide: %s -> %s (%dx%d) %s\n", argv[2], argv[3], img.width, img.height, ok ? "ok" : "FAILED");
        UnloadImage(img);
        return ok ? 0 : 1;
    }
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
    // --seed-survey COUNT FIRSTSEED OUTDIR [THREADS] [SUBSET] (#532): how alike different seeds' start planets are. No window.
    if (argc >= 5 && std::string(argv[1]) == "--seed-survey") return SeedSurvey_Run(argc, argv);
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
    double playNewSeconds = 0.0;   // --play-new S (#387): start a new game as the menu does (settings, spawn preload), quit S seconds after entering
    double playLoadSeconds = 0.0;  // --play-load S SAVE (#509): load SAVE as the menu does, then walk and sweep the view for S seconds
    std::string playLoadPath;
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
    float benchUnfocus = 0.0f;    // --bench-unfocus S (#510): unfocused S s during the spawn preload, and while playing (FPS before/after)
    int unfocusPhase = 0; double unfocusT0 = 0.0; long unfocusFrames = 0; double unfocusFpsBefore = 0.0;
    double unfocusSim0 = 0.0; bool unfocusSimOk = true; // #525: sim time while unfocused (s)
    bool benchUnfocused = false;  // the game sees its window unfocused (scripted: the window itself is never minimized or raised)
    double benchFlySeconds = 0.0, flyStart = -1.0; // --test-fly S: autopilot the ship to the nearest other asteroid
    bool benchSuffocate = false;   // --test-suffocate: run out of oxygen away from the ship, die, respawn
    double benchOrbitSeconds = 0.0, orbitStart = -1.0; // --test-orbit S: coast on a circular orbit for S seconds
    float orbitMin = 0, orbitMax = 0;
    bool orbitCrashed = false;
    bool benchRide = false;        // --test-ride: stand on a spinning body, then jump
    bool benchPark = false;        // --test-park: park the ship on a spinning body
    bool sunOverride = false;      // --sun x y z: sun direction (e.g. to test night-side lighting)
    bool benchDay = false;         // --day: a mid-morning sun over the start (perf benches: shadows and all, #291)
    float benchSunElev = -999.0f;  // --sun-elev DEG: the sun straight ahead of the bench camera at this elevation (sunrise / sunset shots, #276)
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
    bool benchDoorway = false; int doorPhase = 0, doorLane = 0; double doorStart = 0.0; uint32_t doorPlot = 0; bool doorPass[2] = { false, false }; // --test-doorway (#255): 2 m openings pass, 1.5 m ones block
    bool benchMenus = false;       // --test-menus: main menu, new game, pause, save & quit, load
    bool benchMeteor = false;      // --test-meteor (#177)
    bool benchShower = false;      // --test-shower (#316): one world-wide meteor shower, rolled by the authority
    std::string hostSave;          // --host <save> (#311): a dedicated server: load (or create) the save, run headless, listen
    uint16_t hostPort = kDefaultHostPort; // --port N
    bool benchHost = false;        // --test-host (#311): with --host, a loopback client joins
    bool seedGiven = false;        // --seed N was passed (a --host new world: else a random seed)
    double autosaveMinutes = kDefaultAutosaveMinutes; // --autosave-minutes N (--host; 0 = off)
    unsigned hostSaveSeed = 0;     // the seed an existing --host save holds
    bool hostNewWorld = false;     // --host: the save did not exist yet
    bool benchHeadless = false;    // --test-headless (#311): the sim runs with --headless, nothing drawn
    bool benchTrapping = false;    // --test-trapping (#176)
    bool benchCombat = false;      // --test-combat: knockback, shield blocks and stamina, deaths by difficulty, the grave (#165)
    bool benchFishing = false;     // --test-fishing (#176)
    bool benchStructures = false;  // --test-structures: ruins, wrecks and wonders (#164)
    bool vramReport = false;       // --vram: GPU memory per subsystem at the end of a bench (#302)
    bool benchFill = false;        // --test-fill: an exact sphere of dirt (#228)
    bool benchLamps = false;       // --test-lamps: 100 lamp cubes light the ground (#139)
    bool benchGiant = false;       // --shot-giant: one tree of each size class (#274)
    bool benchTorch = false;       // --test-torch: a held torch lights the ground (#218)
    bool benchPosture = false, postureStarted = false; // --test-posture (#578)
    bool benchPause = false;       // --test-pause (#511): behind the pause menu nothing moves (the game clock and the sim stop)
    bool benchBite = false, biteStarted = false; // --test-bite (#165)
    bool benchShotDig = false;     // --shot-dig (#533): the dig preview ring, looking down at the ground
    bool benchF2 = false;          // --test-f2 (#579): take an F2 screenshot
    bool benchTime = false;        // --test-time: dawn start (#259), a +6 h jump carries the player (#275)
    bool benchPickups = false;     // --test-pickups: the pickup feed (#208): repeats merge, a screenshot while it shows
    bool shotInventory = false;    // --shot-inventory: long item names in their cells (#219)
    bool shotGlow = false;         // --shot-glow: hot items glow in the inventory and on the hotbar (#580)
    bool shotChest = false;        // --shot-chest: a filled chest open beside the inventory (#495: 50 cells, 10 wide)
    bool benchFoundation = false;  // --test-foundation: a part half in the ground clears its cell and nothing else (#210)
    bool shotPreview = false;      // --shot-preview: the build preview in the part's real shape (#197)
    bool testInvCursor = false;    // --test-inv-cursor (#192)
    bool noHud = false;            // --no-hud: no HUD, hands or overlay in play (clean screenshots: the menu slides, #386)
    bool testVein = false;         // --test-vein (#207)
    bool testContainerCursor = false; // --test-container-cursor (#192): a chest with the same slot behaviour as the inventory
    int testBugReport = 0;         // --test-bug-report (#169): frames since the scripted report was sent (0: not asked)
    int pickupFrames = 0;
    int timeStep = 0, timeFrames = 0; float timeElev0 = 0.0f; Vector3 timeLocal0 = { 0, 0, 0 }; bool timeOk = true;
    bool benchWater = false;       // --test-water: rivers and ponds near the start, their beds and clay
    bool benchFlyFollow = false;   // --test-flyfollow: the fly camera rides its reference body
    bool benchHotbar = false;      // --test-hotbar: materials and consumables from the hotbar
    bool benchInventory = false;   // --test-inventory: grid stacking, fractions, sorting, saves
    bool benchMap = false;         // --test-map: compass, exploration, waypoints, map images, saves
    bool benchQuests = false;      // --test-quests: tiers open by their tasks; quest state through a save
    bool benchHunger = false, benchFarm = false; // --test-hunger, --test-farm
    bool benchViewFall = false; Vector3 benchFallLip{ 0, 1, 0 }, benchFallDown{ 1, 0, 0 }; float benchFallDist = 80.0f; // --view-fall (#540), body-local
    bool benchRivers = false, riverOk = false; int riverPhase = 0; double riverStart = 0.0; float riverLongest = 0.0f, riverBuild = 0.0f, riverWet = 0.0f, riverLevel = 0.0f; Vector3 riverDir{ 0, 0, 0 }; // --test-rivers
    bool benchTools = false, toolsOk = true; int toolsPhase = 0, toolsLife = 0, toolsTilled = 0; double toolsStart = 0.0; // --test-tools
    bool benchSuit = false, suitOk = true; int suitPhase = 0; double suitStart = 0.0; float suitFoodRatio = 0.0f, suitWater = 0.0f, suitRebreather = 0.0f; // --test-suit
    int benchShotTabs = 0; // --shot-tabs
    bool benchBiomes = false; int benchShotBiomes = -1; // --test-biomes, --shot-biomes
    bool benchOil = false, oilRateOk = false; int oilPhase = 0; double oilT0 = 0.0; // --test-oil
    bool benchSpring = false; int springPhase = 0; double springT0 = 0.0; // --test-spring
    bool benchShotAnimals = false; int animalShotPhase = 0; // --shot-animals
    bool benchShotMachines = false; int machineShotPhase = 0; // --shot-machines
    bool benchParticles = false; int particlePhase = 0; double particleT0 = 0.0; Vector3 particleEye{}, particleAt{}; // --shot-particles
    bool benchShadow = false; int shadowPhase = 0; double shadowT0 = 0.0; Vector3 shadowEye{}, shadowAt{}; // --shot-shadow
    bool benchFarForest = false; int farForestPhase = 0; double farForestT0 = 0.0; Vector3 farForestEye{}, farForestAt{}; // --shot-far-forest
    bool benchTerrain = false; int terrainPhase = 0; double terrainT0 = 0.0; Vector3 terrainAt{}, terrainOut{}; // --shot-terrain
    bool benchWind = false; // --test-wind (#650)
    bool benchCampfire = false; int campfirePhase = 0; double campfireT0 = 0.0; uint32_t campfirePlot = 0; int campfirePart = -1; Vector3 campfireAt{}, campfireOut{}; // --shot-campfire (#552)
    bool benchCaveMouth = false; int mouthPhase = 0; double mouthT0 = 0.0; Vector3 mouthAt{}, mouthOut{}; // --shot-cave-mouth (#549)
    bool benchMountain = false; int mountainPhase = 0; double mountainT0 = 0.0; Vector3 mountainAt{}, mountainOut{}; // --shot-mountain (#545)
    bool benchShelter = false; int shelterPhase = 0; double shelterT0 = 0.0; Vector3 shelterEye{}, shelterLook{}; // --shot-shelter (#369)
    bool benchRoom = false; int roomPhase = 0; double roomT0 = 0.0; uint32_t roomPlot = 0; Vector3 roomEye{}, roomLook{}; // --shot-room (#369)
    bool benchFacing = false; // --test-facing (#196)
    bool benchPlacement = false; // --test-placement (#550, #551): parts centred on the aim, no flicker between two cells
    bool benchPumpHand = false; // --test-pump-hand (#280): the pump test with a hand pump worked for 60 s
    bool benchPumpBoiler = false; int pumpBoilerPart = -1; // --test-pump-boiler (#229 review): the pipes end at a full boiler; the world loses no water
    bool benchCaveWater = false; int caveWaterPhase = 0; double caveWaterT0 = 0.0; Vector3 caveAt{}, caveOpenAt{}, caveRiverDir{}; float caveLevel = 0.0f, caveDist = 0.0f; // --test-cave-water (#559)
    bool benchPump = false; int pumpPhase = 0, pumpCable = -1; double pumpT0 = 0.0, pumpedAtCut = 0.0; uint32_t pumpPlot = 0; Vector3 pumpBasin{}; float pumpRise = 0.0f; // --test-pump
    bool benchBucket = false; int bucketPhase = 0, bucketTrips = 0; double bucketT0 = 0.0; Vector3 bucketWater{}, bucketBasin{}, bucketField{}; // --test-bucket
    float benchFluidLoad = 0.0f; int flPhase = 0; double flT0 = 0.0; Vector3 flMid{}, flE{}, flN{}; // --test-fluid-load S (#638/#642 perf)
    bool benchDamBreak = false; int dbPhase = 0; double dbT0 = 0.0; float dbPoured = 0.0f, dbFloor = -2.5f; // --test-dam-break (#680: Ritter dam break, physics reference)
    bool benchLeak = false; int lkPhase = 0; double lkT0 = 0.0; float lkPoured = 0.0f; // --test-leak (#680: a hole low in a wall)
    bool benchLeakParts = false; Plot* lkPlot = nullptr; std::vector<int> lkHole; // --test-leak-parts (#692): the wall is placed stone blocks
    bool benchFlowShade = false; int fsPhase = 0; double fsT0 = 0.0; Vector3 fsMid{}, fsE{}, fsN{}; // --test-flow-shade (#677, #696)
    bool benchWaterSpeed = false; int wsPhase = 0; double wsT0 = 0.0, wsFront[3] = { -1, -1, -1 }, wsSleep = -1.0; Vector3 wsSite{}, wsEye{}; float wsMaxFront = 0.0f, wsVol0 = 0.0f; // --test-water-speed (#638)
    bool benchPond = false; int pondPhase = 0; double pondT0 = 0.0; Vector3 pondDirB{}, pondOutB{}, pondPit{}, pondMid{}; float pondLevelB = 0.0f, pondRB = 0.0f; // --test-pond
    bool benchFlowSea = false; // --test-flow-sea (#190): the flow test from the sea instead of inland water
    bool benchFlow = false; int flowPhase = 0; double flowStart = 0.0, flowEdit = 0.0, flowReachT = -1.0, flowCaveT = -1.0, flowSleepT = -1.0; Vector3 flowDir{}, flowOut{}, flowEnd{}, flowCave{}; float flowLevel = 0.0f; // --test-flow
    float damDrift = 0.0f; // --test-dam: the current mid-channel (0: the check failed)
    bool benchPressure = false; int prPhase = 0; double prT0 = 0.0; Vector3 prDir{}; float prGround = 0.0f; // --test-pressure (#393)
    bool benchFlowMeter = false; int fmPhase = 0; double fmT0 = 0.0; Plot* fmPlot = nullptr; // --test-flow-meter (#392)
    bool benchPartWall = false; int pwPhase = 0; double pwT0 = 0.0; Plot* pwPlot = nullptr; std::vector<int> pwParts; bool pwHeld = false, pwWestWet = false; // --test-part-wall (#391)
    bool benchDam = false; int damPhase = 0; double damT0 = 0.0; Vector3 damDir{}, damDown{}, damSide{}, damUp{}, damDn{}; float damLevel = 0.0f, damHW = 0.0f, damDepth = 0.0f, damCrest = 0.0f, damLvU = 0.0f, damLvD = 0.0f, damMaxU = -1e9f, damMinD = 1e9f, damLateD = -1e9f; // --test-dam (#361)
    bool benchDamDivert = false, benchDamDivertFull = false; int divPhase = 0; Vector3 divDir{}; float divLevel0 = 0.0f, divDepth = 0.0f; int divPiece = -1; double divDryT = -1.0, divOpenT = -1.0; std::vector<TerrainEdit> divWall; // --test-dam-divert (#656)
    bool benchForest = false; int forestPhase = 0; double forestStart = 0.0; Vector3 forestCam{}, forestTarget{}, forestUp{}; bool outcropHigh = false; // --shot-forest (--outcrop-high: the highest outcrop in range, not the nearest)
    bool benchRecipes = false; // --test-recipes
    bool benchWeather = false, weatherOk = true; int weatherPhase = 0, weatherCounts[6] = {}, weatherChanges = 0; double weatherStart = 0.0; float weatherChill = 0.0f; // --test-weather
    bool benchDecor = false, decorOk = true; int decorPhase = 0, decorGems = 0, decorGemKinds = 0, decorDigs = 0, decorPlaced = 0; double decorStart = 0.0; Vector3 decorCamPos{}, decorCamTarget{}, decorCamUp{ 0, 1, 0 }; // --test-decor
    bool benchMaterials = false, matOk = true; int matPhase = 0, matPlaced = 0; double matStart = 0.0; std::string matKinds; Vector3 matCamPos{}, matCamTarget{}, matCamUp{ 0, 1, 0 }; // --test-materials
    bool benchDamage = false, damageOk = true; int damagePhase = 0; double damageStart = 0.0; float damageThrust = 0.0f, damageRepaired = 0.0f, damageVacuum = 0.0f; // --test-damage
    bool benchCrafting = false, craftOk = true; int craftPhase = 0, craftChest = -1, craftInv0 = 0, craftWood0 = 0, craftDeposited = 0, craftInChest = 0; double craftStart = 0.0; Plot* craftPlot = nullptr; // --test-crafting
    int threatSeed = 3; // --threat-seed N (#253)
    bool benchAnimalPath = false; int animalPathPhase = 0; double animalPathStart = 0.0; // --test-animal-path (#253)
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
    bool menuEsc = false; // --test-menus: an Esc press for the in-game Esc check, as if the key went down this frame
    int menuBakesBeforeLoad = 0; // --test-menus: erosion bakes before Load Game (#505: the load must not bake)
    double lastQuitMs = 0.0;     // how long the last quit to the menu took to unload the world (#270)
    int wallPhase = 0; bool wallOk = true;
    double wallStart = 0.0, wallMaxStep = 0.0, wallMaxStepAt = 0.0;
    float wallMaxImpact = 0.0f, wallMinHealth = 100.0f, wallMaxSpeed = 0.0f; int wallCamInside = 0; // frames with the camera inside terrain (#255)
    int wallCamShake = 0; double wallCamShakeMax = 0.0; DVec3 wallCamPrev{}, wallCamStepPrev{}; int wallCamKnown = 0, wallCamShakePushed = 0, wallCamShakeGrounded = 0; // camera shake (#255)
    bool benchLedge = false; int ledgePhase = 0, ledgePass = 0; bool ledgeOk = true, ledgeAir = false; double ledgeStart = 0.0, ledgeAirStart = 0.0; float ledgeV0 = 0.0f, ledgeExcess = 0.0f; float ledgeWalk0 = 0.0f;
    bool benchTor = false; int torPhase = 0, torHeading = 0, torFlicker = 0, torPops = 0; double torStart = 0.0, torAirAt = -1.0, torMaxStep = 0.0; float torMaxUp = 0.0f; Vector3 torDir{}; DVec3 torLast{}; bool torWasGrounded = true; // --test-tor (#242)
    bool benchMachines = false; int machPhase = 0; double machStart = 0.0; // --test-machines (#193 #194 #229)
    bool benchHover = false; int hoverPhase = 0; double hoverStart = 0.0; // --test-hover (#163)
    bool benchNetIds = false;      // --test-netids: stable network ids over a save and load (#317)
    bool benchNet = false;         // --test-net: connect by IP over loopback (GameNetworkingSockets, #310)
    bool benchFingerprint = false; // --test-fingerprint: world fingerprint and the gameplay RNG (#316)
    bool benchInterest = false; int interestPhase = 0; // --test-interest: a second interest centre (#318)
    bool benchTwoPlayers = false; int twoPhase = 0; // --test-two-players: a second Player in its own frame (#314)
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
    // Planet walk legs (#363): M metres as 100 m legs at places drawn from the seed, teleporting between them.
    constexpr double kWalkLegM = 100.0, kWalkLegMaxS = 45.0;
    int walkLegs = 0, walkLeg = 0, walkLegsDone = 0, walkSettleFrames = 0, walkLegGrounded = 0, walkLegFrames = 0;
    bool walkLegReady = true;
    double walkLegStart = 0.0, walkLegPath = 0.0, walkLegSlope = 0.0;
    int walkLegJumps = 0, walkFellThrough = 0; bool walkLegDied = false; Vector3 walkPickDir{ 0, 1, 0 };
    DVec3 walkLast{};
    Vector3 benchCamOffset = { -70.0f, 35.0f, 70.0f }; // --cam dx dy dz: bench camera offset from the planetoid
    float benchViewHeight = 40.0f; // --bench-height H: a planet bench's view height over the start (m)
    float benchPitch = -999.0f; bool benchUphill = false; // --bench-pitch DEG (up > 0), --bench-uphill: face up the eroded slope (#357 cliff shots)
    bool benchCamGiven = false;   // --cam given: use the offset (otherwise a planet bench looks along the surface from 40 m up, #291)
    if (const char* e = getenv("RAYLIBGAME_EROSION"); e && e[0]) gErosionPrototype = e[0] != '0'; // on by default (#362, user 2026-11-15); 0: the noise terrain (A/B)
    if (const char* e3 = getenv("RAYLIBGAME_3D")) gTerrain3D = e3[0] == '1'; else if (gErosionPrototype) gTerrain3D = false; // no 3D band in erosion tests (RAYLIBGAME_3D=0 for baselines)
    if (const char* g = getenv("RAYLIBGAME_GULLIES")) gGullies = g[0] == '1';
    if (const char* g = getenv("RAYLIBGAME_CLIFF_FACES")) gCliffFaces = g[0] == '1'; // #588 A/B
    if (const char* g = getenv("RAYLIBGAME_FOLDS")) gStrataFolds = g[0] == '1'; // #591 A/B
    if (const char* c = getenv("RAYLIBGAME_SKY_PROBES")) gSkyProbes = c[0] == '1'; // #369 A/B
    if (const char* c = getenv("RAYLIBGAME_AO")) gTerrainAO = c[0] != '0'; // #132 A/B
    for (int i = 1; i < argc; ++i) {
        if (Recorder_ParseArg(argc, argv, i)) continue; // --record and friends (#375)
        if (std::string(argv[i]) == "--bench") benchSeconds = (i + 1 < argc) ? atof(argv[i + 1]) : 10.0;
        if (std::string(argv[i]) == "--ui") DebugUI_Toggle(); // open the debug panel at start
        if (std::string(argv[i]) == "--seed" && i + 1 < argc) worldGen.seed = (unsigned int)strtoul(argv[i + 1], nullptr, 10);
        if (std::string(argv[i]) == "--host" && i + 1 < argc) hostSave = argv[i + 1];
        if (std::string(argv[i]) == "--port" && i + 1 < argc) hostPort = (uint16_t)atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-host") benchHost = true;
        if (std::string(argv[i]) == "--seed") seedGiven = true;
        if (std::string(argv[i]) == "--autosave-minutes" && i + 1 < argc) autosaveMinutes = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--nocull") frustumCullingEnabled = false;
        if (std::string(argv[i]) == "--overlay") startWithOverlay = true;
        if (std::string(argv[i]) == "--cam" && i + 3 < argc) { benchCamGiven = true; benchCamOffset = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
        if (std::string(argv[i]) == "--test-dig" && i + 1 < argc) benchDigStrokes = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-paint") benchPaint = true;
        if (std::string(argv[i]) == "--test-saveload") benchSaveLoad = true;
        if (std::string(argv[i]) == "--test-walk" && i + 1 < argc) benchWalkSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--no-eye-smooth") player.eyeSmoothing = false;
        if (std::string(argv[i]) == "--test-planet-walk" && i + 1 < argc) { walkTargetM = atof(argv[i + 1]); benchWalkSeconds = 600.0; } // gate: 1 km on the planet, in 100 m legs (#363)
        if (std::string(argv[i]) == "--bench-walk") benchWalkMode = true;
        if (std::string(argv[i]) == "--bench-stand") { benchWalkMode = true; benchStandMode = true; }
        if (std::string(argv[i]) == "--perf-log") { perfLog = true; gJobLog = true; gFpsLog = true; }
        if (std::string(argv[i]) == "--play-new" && i + 1 < argc) playNewSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--play-load" && i + 2 < argc) { playLoadSeconds = atof(argv[i + 1]); playLoadPath = argv[i + 2]; }
        if (std::string(argv[i]) == "--bench-fly" && i + 1 < argc) benchFlySpeed = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--bench-fly-back") benchFlyBack = true;
        if (std::string(argv[i]) == "--water-stats") benchWaterStats = true;
        if (std::string(argv[i]) == "--bench-height" && i + 1 < argc) benchViewHeight = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--bench-pitch" && i + 1 < argc) benchPitch = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--bench-uphill") benchUphill = true;
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
        if (std::string(argv[i]) == "--test-doorway") benchDoorway = true;
        if (std::string(argv[i]) == "--test-ledge") benchLedge = true; // walk and sprint off a pit edge (#240)
        if (std::string(argv[i]) == "--test-tor") benchTor = true; // walk across a low (beach) tor crown, logging collision pops (#242)
        if (std::string(argv[i]) == "--test-machines") benchMachines = true; // fuel slot, crafts consumed at the end, the boiler's water and heat
        if (std::string(argv[i]) == "--test-hover") benchHover = true; // the flight pack (#163)
        if (std::string(argv[i]) == "--test-netids") benchNetIds = true; // stable network ids (#317)
        if (std::string(argv[i]) == "--test-net") benchNet = true; // loopback connect (#310)
        if (std::string(argv[i]) == "--test-fingerprint") benchFingerprint = true; // world fingerprint, game RNG (#316)
        if (std::string(argv[i]) == "--test-interest") benchInterest = true; // a second interest centre (#318)
        if (std::string(argv[i]) == "--test-two-players") benchTwoPlayers = true; // two players side by side (#314)
        if (std::string(argv[i]) == "--test-fall") benchFall = true; // drops from 3-100 m: fall damage by height
        if (std::string(argv[i]) == "--test-telescope-close") benchTeleClose = true; // zoom 16x at the horizon, close: worst frame after (#261)
        if (std::string(argv[i]) == "--test-menus") benchMenus = true;
        if (std::string(argv[i]) == "--test-time") benchTime = true;
        if (std::string(argv[i]) == "--test-torch") benchTorch = true;
        if (std::string(argv[i]) == "--test-posture") benchPosture = true;
        if (std::string(argv[i]) == "--test-pause") benchPause = true;
        if (std::string(argv[i]) == "--test-bite") benchBite = true;
        if (std::string(argv[i]) == "--shot-dig") benchShotDig = true;
        if (std::string(argv[i]) == "--test-f2") benchF2 = true;
        if (std::string(argv[i]) == "--shot-giant") benchGiant = true;
        if (std::string(argv[i]) == "--test-lamps") benchLamps = true;
        if (std::string(argv[i]) == "--test-fill") benchFill = true;
        if (std::string(argv[i]) == "--vram") vramReport = true;
        if (std::string(argv[i]) == "--test-structures") benchStructures = true;
        if (std::string(argv[i]) == "--test-fishing") benchFishing = true;
        if (std::string(argv[i]) == "--test-trapping") benchTrapping = true;
        if (std::string(argv[i]) == "--test-combat") benchCombat = true;
        if (std::string(argv[i]) == "--test-meteor") benchMeteor = true;
        if (std::string(argv[i]) == "--test-shower") benchShower = true;
        if (std::string(argv[i]) == "--test-headless") benchHeadless = true;
        if (std::string(argv[i]) == "--test-pickups") benchPickups = true;
        if (std::string(argv[i]) == "--shot-inventory") shotInventory = true;
        if (std::string(argv[i]) == "--shot-glow") shotGlow = true;
        if (std::string(argv[i]) == "--shot-chest") shotChest = true;
        if (std::string(argv[i]) == "--test-foundation") benchFoundation = true;
        if (std::string(argv[i]) == "--shot-preview") shotPreview = true;
        if (std::string(argv[i]) == "--test-inv-cursor") testInvCursor = true;
        if (std::string(argv[i]) == "--no-hud") noHud = true;
        if (std::string(argv[i]) == "--test-vein") testVein = true;
        if (std::string(argv[i]) == "--test-container-cursor") testContainerCursor = true;
        if (std::string(argv[i]) == "--test-bug-report") testBugReport = 1;
        if (std::string(argv[i]) == "--test-saves") Menu_TestSaves();
        if (std::string(argv[i]) == "--test-machine-stack") { const bool ok = Machines_TestStack(); fflush(stdout); return ok ? 0 : 1; }
        if (std::string(argv[i]) == "--test-water") benchWater = true;
        if (std::string(argv[i]) == "--test-cave-water") benchCaveWater = true;
        if (std::string(argv[i]) == "--test-flyfollow") benchFlyFollow = true;
        if (std::string(argv[i]) == "--test-hotbar") benchHotbar = true;
        if (std::string(argv[i]) == "--test-inventory") benchInventory = true;
        if (std::string(argv[i]) == "--test-map") benchMap = true;
        if (std::string(argv[i]) == "--test-quests") benchQuests = true;
        if (std::string(argv[i]) == "--test-hunger") benchHunger = true;
        if (std::string(argv[i]) == "--test-farm") benchFarm = true;
        if (std::string(argv[i]) == "--test-animals") benchAnimals = true;
        if (std::string(argv[i]) == "--test-threats") benchThreats = true;
        if (std::string(argv[i]) == "--threat-seed" && i + 1 < argc) threatSeed = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-animal-path") benchAnimalPath = true;
        if (std::string(argv[i]) == "--shot-tabs") benchShotTabs = 1;
        if (std::string(argv[i]) == "--shot-forest") benchForest = true;
        if (std::string(argv[i]) == "--outcrop-high") outcropHigh = true;
        if (std::string(argv[i]) == "--test-flow") benchFlow = true;
        if (std::string(argv[i]) == "--test-dam") benchDam = true; // #361: a dam backs a river up and it spills over
        if (std::string(argv[i]) == "--test-dam-divert") benchDam = benchDamDivert = true; // #656: a dam that holds the river dries it downstream; removed, the river comes back
        if (std::string(argv[i]) == "--dam-divert-full") benchDamDivertFull = true; // #656: wait for the full return (flow 0.9, up to 480 s) instead of the 90 s rate check
        if (std::string(argv[i]) == "--test-part-wall") benchPartWall = true; // #391: a wall of build blocks holds water
        if (std::string(argv[i]) == "--test-flow-meter") benchFlowMeter = true; // #392: the flow readout reads a known stream
        if (std::string(argv[i]) == "--test-pressure") benchPressure = true; // #393: hydrostatic pressure in a shaft and a pond
        if (std::string(argv[i]) == "--test-flow-sea") { benchFlow = true; benchFlowSea = true; } // #190: the channel from the sea
        if (std::string(argv[i]) == "--test-pond") benchPond = true;
        if (std::string(argv[i]) == "--test-bucket") benchBucket = true;
        if (std::string(argv[i]) == "--test-water-speed") benchWaterSpeed = true;
        if (std::string(argv[i]) == "--sandbox") gSandboxWorld = true; // #662: a flat test plane
        if (std::string(argv[i]) == "--test-fluid-load" && i + 1 < argc) benchFluidLoad = (float)atof(argv[i + 1]); // basin side (m)
        if (std::string(argv[i]) == "--test-dam-break") benchDamBreak = true; // best with --sandbox
        if (std::string(argv[i]) == "--test-leak") benchLeak = true; // best with --sandbox
        if (std::string(argv[i]) == "--test-leak-parts") benchLeak = benchLeakParts = true; // #692: the wall is placed blocks
        if (std::string(argv[i]) == "--test-flow-shade") benchFlowShade = true; // #677: a pool, a chute and a fall, for the water's look
        if (std::string(argv[i]) == "--test-pump") benchPump = true;
        if (std::string(argv[i]) == "--test-pump-hand") { benchPump = true; benchPumpHand = true; }
        if (std::string(argv[i]) == "--test-pump-boiler") { benchPump = true; benchPumpBoiler = true; }
        if (std::string(argv[i]) == "--shot-terrain") benchTerrain = true;
        if (std::string(argv[i]) == "--test-wind") benchWind = true;
        if (std::string(argv[i]) == "--shot-campfire") benchCampfire = true;
        if (std::string(argv[i]) == "--shot-cave-mouth") benchCaveMouth = true;
        if (std::string(argv[i]) == "--shot-mountain") benchMountain = true;
        if (std::string(argv[i]) == "--shot-shelter") benchShelter = true;
        if (std::string(argv[i]) == "--shot-room") benchRoom = true;
        if (std::string(argv[i]) == "--test-facing") benchFacing = true;
        if (std::string(argv[i]) == "--test-placement") benchPlacement = true;
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
        if (std::string(argv[i]) == "--view-fall-dist" && i + 1 < argc) benchFallDist = (float)atof(argv[i + 1]); // --view-fall's distance (m)
        if (std::string(argv[i]) == "--view-fall" && i + 6 < argc) { // --view-fall LIP DOWN (#540): the bench camera looks at a waterfall
            benchViewFall = true;
            benchFallLip = Vector3Normalize({ (float)atof(argv[i + 1]), (float)atof(argv[i + 2]), (float)atof(argv[i + 3]) });
            benchFallDown = Vector3Normalize({ (float)atof(argv[i + 4]), (float)atof(argv[i + 5]), (float)atof(argv[i + 6]) });
        }
        if (std::string(argv[i]) == "--test-viewmodel") benchViewmodel = true;
        if (std::string(argv[i]) == "--test-ring") benchRing = true;
        if (std::string(argv[i]) == "--test-rocket") benchRocket = true;
        if (std::string(argv[i]) == "--test-station") benchStation = true;
        if (std::string(argv[i]) == "--test-blueprint") benchBlueprint = true;
        if (std::string(argv[i]) == "--test-smallbody") benchSmallBody = true;
        if (std::string(argv[i]) == "--test-trees") benchTrees = true;
        if (std::string(argv[i]) == "--test-chop") benchChop = true;
        if (std::string(argv[i]) == "--test-persist" && i + 1 < argc) benchPersist = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-build" && i + 1 < argc) { benchBuild = atoi(argv[i + 1]); Structures_SetEnabled(false); } // no world structures: it counts every part on the planet
        if (std::string(argv[i]) == "--test-factory" && i + 1 < argc) benchFactory = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--test-progression") { benchProgression = true; progressionFull = i + 1 < argc && std::string(argv[i + 1]) == "full"; progressionEv = i + 1 < argc && std::string(argv[i + 1]) == "ev"; }
        if (std::string(argv[i]) == "--time-scale" && i + 1 < argc) timeScale = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--crack-check") crackCheck = true;
        if (std::string(argv[i]) == "--no-foliage") foliageEnabled = false;
        if (std::string(argv[i]) == "--grass-dist" && i + 1 < argc) kGrassDrawDistance = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--tree-dist" && i + 1 < argc) kTreeDrawDistance = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--sun" && i + 3 < argc) { sunOverride = true; sunArg = Vector3Normalize({ (float)atof(argv[i + 1]), (float)atof(argv[i + 2]), (float)atof(argv[i + 3]) }); }
        if (std::string(argv[i]) == "--day") benchDay = true;
        if (std::string(argv[i]) == "--sun-elev" && i + 1 < argc) benchSunElev = (float)atof(argv[i + 1]);
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
            // Near a body (#232): with a tilted local up and a camera started rolled, random turns (roll pressed too) keep the
            // horizon level once eased in, the pitch inside +-87 deg, W/S/A/D level and Space along up.
            {
                const Vector3 g = Vector3Normalize({ 0.3f, 1.0f, -0.2f });
                Camera3D gc{ { 0, 0, 0 }, { 0, 0.2f, 1 }, Vector3Normalize({ 0.5f, 1, 0 }), 60.0f, CAMERA_PERSPECTIVE };
                for (int s = 0; s < 60; ++s) ApplyFlyInput(&gc, FlyInput{}, 1.0f / 60.0f, &g); // ease the roll out
                float worstRoll = 0.0f, maxPitch = 0.0f, worstLevel = 0.0f, worstUp = 1.0f;
                for (int s = 0; s < 300; ++s) {
                    FlyInput turn;
                    turn.yaw = ang(rng); turn.pitch = 2.0f * ang(rng); turn.roll = ang(rng);
                    ApplyFlyInput(&gc, turn, 1.0f, &g); // dt 1: the up settles at once
                    const Vector3 f = Vector3Normalize(Vector3Subtract(gc.target, gc.position));
                    const Vector3 r = Vector3Normalize(Vector3CrossProduct(f, gc.up));
                    worstRoll = std::max(worstRoll, fabsf(Vector3DotProduct(r, g)));            // right stays level
                    maxPitch = std::max(maxPitch, fabsf(asinf(std::clamp(Vector3DotProduct(f, g), -1.0f, 1.0f))) * RAD2DEG);
                    for (int a = 0; a < 3; ++a) {
                        FlyInput mv;
                        (a == 0 ? mv.forward : a == 1 ? mv.right : mv.up) = 1.0f;
                        const Vector3 before = gc.position;
                        ApplyFlyInput(&gc, mv, 0.1f, &g);
                        const Vector3 d = Vector3Normalize(Vector3Subtract(gc.position, before));
                        if (a < 2) worstLevel = std::max(worstLevel, fabsf(Vector3DotProduct(d, g)));
                        else worstUp = std::min(worstUp, Vector3DotProduct(d, g));
                    }
                }
                const bool ok = worstRoll < 1e-3f && maxPitch <= 87.01f && worstLevel < 1e-3f && worstUp > 0.999f;
                printf("bench: flycam near a body (#232): horizon tilt %.5f, max pitch %.2f deg, WASD off level %.5f, Space along up %.5f -> %s\n",
                       worstRoll, maxPitch, worstLevel, worstUp, ok ? "PASS" : "FAIL");
                // #538: quick turns at 60 FPS (no settling time), looking 40 deg up or down: the horizon never tilts, not even
                // for the few frames the up used to take to catch up.
                float worstQuick = 0.0f;
                for (int s = 0; s < 240; ++s) {
                    FlyInput turn;
                    turn.yaw = (s / 20) % 2 ? -0.15f : 0.15f;                   // ~9 rad/s flicks, reversing
                    turn.pitch = s == 0 ? 0.7f : s == 120 ? -1.4f : 0.0f;      // 40 deg up, then down
                    ApplyFlyInput(&gc, turn, 1.0f / 60.0f, &g);
                    const Vector3 f = Vector3Normalize(Vector3Subtract(gc.target, gc.position));
                    const Vector3 r = Vector3Normalize(Vector3CrossProduct(f, gc.up));
                    worstQuick = std::max(worstQuick, fabsf(Vector3DotProduct(r, g)));
                }
                const bool quickOk = worstQuick < 1e-3f;
                printf("bench: flycam quick turns at 60 FPS (#538): worst horizon tilt %.5f (sin of the roll) -> %s\n", worstQuick, quickOk ? "PASS" : "FAIL");
            }
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
        if (std::string(argv[i]) == "--bench-unfocus" && i + 1 < argc) benchUnfocus = (float)atof(argv[i + 1]); // #510
        if (std::string(argv[i]) == "--first" && i + 3 < argc) { worldGen.hasFirst = true; worldGen.first = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
    }
    const bool bench = benchSeconds > 0.0;
    benchMode = bench;
    if (bench) SetWorldSavePath(gSandboxWorld ? "saves/sandbox.sav" : "saves/bench.sav"); // tests save and load; keep the player's F5 save out of it (#662: the sandbox its own)
    if (bench && !benchThreats && !benchAnimalPath) Animals_SetThreats(false); // hostile animals and night spawns would upset other tests

    bool msaa = true;
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--no-msaa") msaa = false;
    if (!bench) { Settings_Load(); if (settings.terrainDetail == 0) msaa = false; } // #302 (user 2026-11-15): Low runs without MSAA (VRAM); benches use --no-msaa
    if (!bench && !gLogLevelFromArgs && !settings.logLevel.empty()) {
        LogLevel l;
        if (Log_ParseLevel(settings.logLevel.c_str(), l)) { Log_SetLevel(l); LOG_INFO("Log", "log level %s (settings.txt)", Log_LevelName(l)); }
        else LOG_WARN("Log", "settings.txt: unknown log_level '%s'", settings.logLevel.c_str());
    }
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--headless") gHeadless = true; // #311: simulate, draw nothing
    if (benchHeadless || !hostSave.empty()) gHeadless = true; // a dedicated server (--host) draws nothing either
    if (!hostSave.empty()) { // --host: the world's seed (user, 2026-11-15) - before the world is generated
        if (FILE* fs = fopen(hostSave.c_str(), "rb")) { // an existing save keeps its own seed (its "seed N" line)
            char line[512];
            while (fgets(line, sizeof line, fs)) if (sscanf(line, "seed %u", &hostSaveSeed) == 1) break;
            fclose(fs);
            if (seedGiven) printf("host: %s keeps its own seed %u (--seed ignored)\n", hostSave.c_str(), hostSaveSeed);
        } else {
            hostNewWorld = true;
            if (!seedGiven) { // a random seed (never 0: that is the generator's own "random" value, not printable up front)
                std::random_device rd;
                do worldGen.seed = (unsigned int)rd(); while (worldGen.seed == 0);
                printf("host: no --seed for the new world: random seed %u\n", worldGen.seed);
            }
        }
    }
    if (gHeadless) msaa = false;
    if (msaa) SetConfigFlags(FLAG_MSAA_4X_HINT);
    bool nativeWindow = false; // --native-window: a bench in the game's own window (borderless, monitor-sized), uncapped
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--native-window") nativeWindow = true;
    // --bench-window: the 1920x1080 bench window for runs that are not benches (scripted play, --play-new): the
    // borderless monitor-sized one covers the screen and takes the input of a user at the machine.
    bool benchWindow = playNewSeconds > 0.0 || playLoadSeconds > 0.0;
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--bench-window") benchWindow = true;
    bool exclusive = false; // --fullscreen: exclusive fullscreen (bypasses the desktop compositor)
    for (int i = 1; i < argc; ++i) if (std::string(argv[i]) == "--fullscreen") { exclusive = true; nativeWindow = true; }
    // #520: who is driving this exe, for the window title: --session <name>, else RAYLIBGAME_SESSION, else the checkout
    // the exe was built in (an agent's worktree: <checkout>/build/<dir>/bin/game.exe). The user's root game.exe sits in
    // the main checkout itself (no build/ in its path) and shows no driver.
    std::string driver;
    for (int i = 1; i + 1 < argc; ++i) if (std::string(argv[i]) == "--session") driver = argv[i + 1];
    if (driver.empty() && getenv("RAYLIBGAME_SESSION")) driver = getenv("RAYLIBGAME_SESSION");
    if (driver.empty()) {
        std::string dir = GetApplicationDirectory();
        for (char& c : dir) if (c == '\\') c = '/';
        const size_t b = dir.rfind("/build/");
        if (b != std::string::npos) {
            const std::string checkout = dir.substr(0, b);
            driver = checkout.substr(checkout.find_last_of('/') + 1);
        }
    }
    // "RaylibGame[ <mode>] - <branch> <commit>[ - <driver>]" (#501 version, #520 mode + driver)
    auto windowTitle = [&](const char* mode) {
        return std::string(TextFormat("RaylibGame%s%s - %s", mode[0] ? " " : "", mode, kGameVersion)) + (driver.empty() ? "" : " - " + driver);
    };
    if (gHeadless) { // #311: a hidden window for the GL context the loaders still need; nothing is drawn into it
        SetConfigFlags(FLAG_WINDOW_HIDDEN);
        InitWindow(320, 180, windowTitle("headless").c_str());
        SetTargetFPS(bench ? std::max(0, benchFps) : kSimHz); // a server ticks at the sim rate
    } else if ((!bench && !benchWindow && playNewSeconds <= 0.0) || nativeWindow) {
        SetConfigFlags(FLAG_WINDOW_RESIZABLE); // Allow window resizing
        SetConfigFlags(FLAG_WINDOW_UNDECORATED); // Remove window decorations for a cleaner look
        SetConfigFlags(FLAG_WINDOW_HIGHDPI); // Enable High DPI support
        if (exclusive) SetConfigFlags(FLAG_FULLSCREEN_MODE);
        InitWindow(GetMonitorWidth(0), GetMonitorHeight(0), windowTitle("").c_str()); // #501: branch + commit (#520: + driver)
        // No VSync: cap at the monitor's refresh rate instead (VSync on a 60 Hz desktop locked the game to 60).
        SetTargetFPS(bench ? std::max(0, benchFps) : benchFps >= 0 ? benchFps : std::max(60, GetMonitorRefreshRate(GetCurrentMonitor())));
    } else {
        // Benches and --play-new (#387): a plain window that opens without taking focus, so an agent run never grabs the
        // keyboard and mouse of whoever is at the machine.
        SetConfigFlags(FLAG_WINDOW_UNFOCUSED);
        int winW = playNewSeconds > 0.0 ? 1280 : 1920, winH = playNewSeconds > 0.0 ? 720 : 1080;
        for (int i = 1; i + 2 < argc; ++i) if (std::string(argv[i]) == "--window") { winW = std::max(320, atoi(argv[i + 1])); winH = std::max(240, atoi(argv[i + 2])); } // --window W H (#512: menus at 1280x720)
        InitWindow(winW, winH, windowTitle(playNewSeconds > 0.0 ? "play-new" : playLoadSeconds > 0.0 ? "play-load" : "bench").c_str()); // #520
        SetTargetFPS(std::max(0, benchFps)); // --fps N caps it (frame-rate-dependence checks); 0 = uncapped
    }
    if (vramReport && IsWindowReady()) GpuMem_Install(); // --vram (#302): count GPU allocations per subsystem from here on
    // A long world-generation step (the erosion bake, #362) runs on a thread; meanwhile draw the loading screen with its
    // progress, so the window stays responsive. From inside a frame (the debug menu's Regenerate) just wait.
    gWorldGenWait = [bench](float progress, const char* what) {
        if (gInFrame || !IsWindowReady()) { std::this_thread::sleep_for(std::chrono::milliseconds(20)); return; }
        FpsLog_Frame("bake");
        static double lastPrint = -10.0;
        if (bench && GetTime() - lastPrint > 5.0) { lastPrint = GetTime(); printf("bench: world generation: %s %.0f%%\n", what, progress * 100.0f); }
        BeginDrawing();
        ClearBackground(BLACK);
        const int bw = GetScreenWidth() / 3, bx = (GetScreenWidth() - bw) / 2, by = GetScreenHeight() / 2;
        DrawText(what, bx, by - 40, 28, RAYWHITE);
        DrawRectangle(bx, by, bw, 14, Fade(RAYWHITE, 0.15f));
        DrawRectangle(bx, by, (int)(bw * std::clamp(progress, 0.02f, 1.0f)), 14, SKYBLUE);
        if (strncmp(what, "Shaping", 7) == 0) // a real bake (a cache read is quick and says "Loading", #505)
            DrawText("A new world's terrain is shaped once (about a minute); loading it later is quick.", bx, by + 26, 16, Fade(RAYWHITE, 0.6f));
        if (settings.showFps) DrawFPS(10, 10); // every screen when the setting is on (#572)
        EndDrawing();
    };
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
    GameUI_InitPickupFeed(); // "+N item" in the bottom right (#208)
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
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) != "--test-heat") continue;
        // --test-heat (#580): hot stacks merge by amount, glow by heat, cool to nothing; cold adds stay cold.
        Inventory inv;
        inv.addAmount("iron_ingot", 2, 1.0f);
        const float fresh = inv.grid()[0].heat;
        const ItemGlow g0 = Item_Glow(inv.grid()[0]);
        inv.addAmount("iron_ingot", 2, 0.0f);
        const float mixed = inv.grid()[0].heat;
        inv.addAmount("stone", 5);
        const float stone = inv.grid()[1].heat;
        inv.cool(Inventory::kCoolSeconds * 0.25f);
        const float cooled = inv.grid()[0].heat;
        inv.cool(Inventory::kCoolSeconds);
        const ItemGlow g1 = Item_Glow(inv.grid()[0]);
        const bool ok = fresh == 1.0f && fabsf(mixed - 0.5f) < 1e-4f && stone == 0.0f && fabsf(cooled - 0.25f) < 1e-3f && inv.grid()[0].heat == 0.0f &&
                        g0.strength > 0.9f && g0.color.r == 255 && g1.strength == 0.0f;
        printf("bench: heat: fresh %.2f (glow %.2f, rgb %d,%d,%d), merged with as many cold %.2f, stone %.2f, after a quarter of the cooling %.2f, cooled %.2f (glow %.2f) -> %s\n",
               fresh, g0.strength, g0.color.r, g0.color.g, g0.color.b, mixed, stone, cooled, inv.grid()[0].heat, g1.strength, ok ? "PASS" : "FAIL");
    }
    GameSave_Register();
    {   // Scenes own GPU resources: keep them in a scope that ends before CloseWindow().
    Scene mainMenu("mainMenu", true); // Main menu scene;
    Scene world("world", false);

    bool windowKeepAlive = true;

    mainMenuHandler(mainMenu);
    worldHandler(world);
    Rocket_Init(&world);
    Station_Init(&world);
    // The main menu has no world (#503: the erosion bake ran before the title screen). Benches skip the menu: they
    // generate here (or load, below).
    if (bench && !loadSaveAtStart && !benchMenus) { // --test-menus starts at the title screen with no world, as a launch does
        generateWorld(world);
        if (!world.objects.empty()) world.camera.target = world.objects[0]->position;
    }
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
    bool hostLoaded = false;
    if (!hostSave.empty()) { // --host (#311): the server's world, saved back to the same file
        SetWorldSavePath(hostSave.c_str());
        if (FILE* probe = fopen(hostSave.c_str(), "rb")) {
            fclose(probe);
            hostLoaded = LoadWorldSave(world);
            printf("host: %s %s\n", hostLoaded ? "loaded" : "could NOT load", hostSave.c_str());
        } else printf("host: %s does not exist: a new world from seed %u, saved there\n", hostSave.c_str(), worldGen.seed);
        // No world exists at launch any more (#503: only benches generate one above), so a new host world is made here.
        if (!hostLoaded && world.objects.empty()) generateWorld(world);
    }
    if (!benchAtSave && !hostLoaded && !world.objects.empty()) StartNewWorld(world, false, !bench);
    if (!hostSave.empty() && !hostLoaded) SaveWorld(world);
    if (!bench) Settings_Load();
    bool worldFresh = !world.objects.empty() && !loadSaveAtStart; // a bench's world generated at launch: New Game uses it as it is
    MenuResult pendingMenu;
    NewWorldChoice pendingNewWorld; // the New World screen's choice, handled at the start of the next frame (#504)
    int newWorldFrames = 0;
    bool justPaused = false;
    bool loading = false;          // the loading screen is up
    double loadStart = 0.0;
    bool spawnLoad = false;        // raised by a spawn or teleport (not entering the game): no minimum time
    int loadFrames = 0;
    int seenSpawns = 0;
    world.onDraw3D = [&world] { Farm_Draw3D(world); Animals_Draw3D(world); Machines_Draw3D(world); InstanceBatch_Flush(world); Weather_Draw3D(world); Waterfalls_Draw3D(world); Particles_Draw(world); if (!ship.piloted && !buildTool.active && !Rocket_Boarded()) TerrainTool_Draw3D(); Fishing_Draw3D(world); Events_Draw3D(world); Combat_DrawGraves(world); BuildTool_Draw3D(world); };
    terrainTool.till = [&world](Planetoid& body, Vector3 ground) { return Farm_Till(world, body, ground); }; // a shovel's RMB
    terrainTool.swapHeld = [](const std::string& from, const std::string& to, const char* ev) { // a bucket filled or emptied
        const int s = Hotbar_Selected();
        if (Hotbar_Slot(s) == from) Hotbar_Bind(s, to);
        QuestBook_Count(ev);
    };

    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--probe-films") {
            // #249: thin near-iso films over the ground. An 81 x 81 m square at the start (1 m apart): each column from 3 m
            // over the surface down to 1 m under it in 5 cm steps; a solid layer above an air gap counts as a film.
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
            if (pf) {
                const Vector3 home = pf->findSpawnDir(p->seed);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                int films = 0, cols = 0; float thickest = 0.0f; Vector3 worstDir = home;
                for (int a = -40; a <= 40; ++a)
                    for (int b = -40; b <= 40; ++b) {
                        const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, a / pf->radius), Vector3Scale(t2, b / pf->radius))));
                        const float surf = pf->radius + pf->elevation(d);
                        ++cols;
                        bool seenSolid = false, seenAirBelowSolid = false; float layer = 0.0f;
                        for (float h = surf + 3.0f; h > surf - 1.0f; h -= 0.05f) { // top down: solid, then air, then solid = a film over a gap
                            const bool solid = p->densityAt(Vector3Scale(d, h)) >= ISO_LEVEL;
                            if (solid && !seenSolid) { seenSolid = true; layer = 0.0f; }
                            else if (solid && seenSolid && !seenAirBelowSolid) layer += 0.05f;
                            else if (!solid && seenSolid) seenAirBelowSolid = true;
                            else if (solid && seenAirBelowSolid) { ++films; if (layer + 0.05f > thickest) { thickest = layer + 0.05f; worstDir = d; } break; }
                        }
                    }
                printf("bench: films (#249): %d of %d columns have a solid layer over an air gap within 3 m of the ground (thickest %.2f m) at dir %.6f,%.6f,%.6f -> %s\n",
                       films, cols, thickest, worstDir.x, worstDir.y, worstDir.z, films == 0 ? "none" : "FOUND");
            }
            return 0;
        }
        else if (std::string(argv[i]) == "--test-dig-shape") {
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
                       n, kVoxelSize, sum / std::max(1, n), worst, sum / std::max(1, n) >= 0.83 && worst >= 0.65 ? "PASS" : "FAIL"); // worst 0.65 (#532: the floor depends on the ground the 40 digs land on: 0.664-0.715 over seeds 1, 7, 42, 100, 2024 with the carve unchanged)
                // (#425: the bar was the first measurement itself, 0.85 / 0.73 on 2026-11-12. The carve hasn't changed since, but
                // the ground at the 40 random spots has (erosion, lips, rock features), and on a 0.25 m lattice the small spheres
                // score ~0.85 +- 0.01 wherever they land: 0.847 on 2026-11-15. A 1 m lattice scored 0.28, so 0.83 still catches
                // any real regression.)
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
            { const float gp[4] = { pf->gen.ocean, pf->gen.mountains, pf->gen.roughness, 1.0f / std::max(pf->gen.biomeSize, 0.1f) }; rlSetUniform(prog.uniform("genParams"), gp, RL_SHADER_UNIFORM_VEC4, 1); } // #504
            prog.setUint("count", (unsigned)N);
            bin.bind(0); bout.bind(1);
            const double g0 = GetTime();
            prog.dispatch((N + 63) / 64);
            std::vector<float> gpu(N);
            bout.download(gpu.data(), N * 4);
            const double gpuMs = (GetTime() - g0) * 1000.0;
            const double c0 = GetTime();
            std::vector<float> cpu(N);
            // The point kernel ports the noise terrain only (no bake): compare against the CPU without it (#561: with erosion
            // on, eroded CPU heights against the GPU's noise read as a 605 m mismatch).
            PlanetField::tIgnoreBake = true;
            for (int k = 0; k < N; ++k) cpu[k] = pf->elevationBase({ in[k * 4], in[k * 4 + 1], in[k * 4 + 2] }, in[k * 4 + 3]);
            PlanetField::tIgnoreBake = false;
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
                timeIt("land (before rivers, ponds)", [&](size_t k) { return pf->profileLandHeight(Vector3Normalize(pts[k]), 1.0f); });
                for (float fpq : { 1.0f, 30.0f, 100.0f }) { double sum = 0, mx = 0; int nland = 0; for (size_t k = 0; k < pts.size(); ++k) { const Vector3 dq = Vector3Normalize(pts[k]); gGullies = false; const float a = pf->profileLandHeight(dq, fpq); gGullies = true; const float b = pf->profileLandHeight(dq, fpq); if (a > 20.0f) { sum += fabs(b - a); mx = std::max(mx, (double)fabs(b - a)); ++nland; } } printf("bench: gullies at fp %.0f m: mean |dh| %.2f m, max %.1f m over %d land points\n", fpq, nland ? sum / nland : 0.0, mx, nland); }
                timeIt("terrain3D (3D band)", [&](size_t k) { const float d = Vector3Length(pts[k]); return pf->terrain3D(pts[k], d, surf[k], 1.0f); });
                timeIt("rockFeature (tors, boulders)", [&](size_t k) { const float d = Vector3Length(pts[k]); return pf->rockFeature(pts[k], d, surf[k], 1.0f); });
                timeIt("density (all, no edits)", [&](size_t k) { return pf->density(pts[k], nullptr, 1.0f); });
                timeIt("waterSourceAt", [&](size_t k) { return (float)pf->waterSourceAt(Vector3Normalize(pts[k])).kind; });
                if (const RiverNetwork* net = pf->rivers()) {
                    // The river current (#361) asks this per river water vertex in the chunk mesher: its cost by a river.
                    timeIt("rivers nearest (60 m)", [&](size_t k) { return net->nearest(Vector3Normalize(pts[k]), 60.0f).dist; });
                    std::vector<Vector3> onRiver;
                    for (size_t k = 0; k < pts.size() && onRiver.size() < 20000; ++k) if (pf->waterSourceAt(Vector3Normalize(pts[k])).kind == 1) onRiver.push_back(pts[k]);
                    if (!onRiver.empty()) {
                        const auto t0 = std::chrono::steady_clock::now();
                        float sink2 = 0.0f;
                        for (int rep = 0; rep < 10; ++rep) for (const Vector3& q : onRiver) sink2 += net->nearest(Vector3Normalize(q), 60.0f).dist;
                        printf("bench: density profile: %-26s %7.0f ns (%zu river points)%s\n", "rivers nearest on rivers", std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / (10.0 * onRiver.size()), onRiver.size(), sink2 < -1.0f ? " " : "");
                    }
                }
            }
            return 0;
        }
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--test-dig-lod") {
            // #679: a 64 x 64 m box dug 3 m deep at the start (best with --sandbox: a flat plane), then the node over the
            // floor's middle built at every level 0-6 with the edit, and the height of its floor vertices (the ones well
            // inside the walls) against the dug floor: mean, rms and worst. Coarse nodes point-sampled the edit lattice,
            // so their floors came out lumpy by up to half a stride, or not dug at all. Also (#693): no lattice point inside
            // the cut may stay solid (specks over a trench); --dig-lod-flush cuts the box the way --test-fluid-load does.
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
            if (!pf || !p->edits) { printf("bench: dig-lod: no planet -> FAIL\n"); return 1; }
            const Vector3 home = pf->findSpawnDir(p->seed);
            const Vector3 eE = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), eN = Vector3CrossProduct(eE, home);
            const Vector3 mid = Vector3Scale(home, pf->radius + pf->elevation(home));
            const float half = 32.0f, depth = 3.0f;
            bool flush = false;
            for (int k = 1; k < argc; ++k) if (std::string(argv[k]) == "--dig-lod-flush") flush = true;
            TerrainDeltasBuilder bld(*p->edits);
            for (int a = 0; a < 4; ++a)
                for (int c = 0; c < 4; ++c) {
                    TerrainEdit e{};
                    // 16 m tiles overlapping by 1 m, their tops 1 m over the ground; --dig-lod-flush: as --test-fluid-load cuts
                    // its basin, tiles sharing their faces and the top flush with the ground. Such faces left films of ground
                    // where they passed within a rounding step of lattice points (#693): they count as floor errors of ~3 m.
                    if (flush) {
                        e.center = Vector3Add(mid, Vector3Add(Vector3Add(Vector3Scale(eE, 16.0f * (a - 1.5f)), Vector3Scale(eN, 16.0f * (c - 1.5f))), Vector3Scale(home, -0.5f * depth)));
                        e.half = { 8.0f, 8.0f, 0.5f * depth };
                    } else {
                        e.center = Vector3Add(mid, Vector3Add(Vector3Add(Vector3Scale(eE, 16.0f * (a - 1.5f)), Vector3Scale(eN, 16.0f * (c - 1.5f))), Vector3Scale(home, 1.0f - (0.5f * depth + 0.5f))));
                        e.half = { 8.5f, 8.5f, 0.5f * depth + 0.5f };
                    }
                    e.carve = true; e.box = true;
                    e.ax = eE; e.ay = eN; e.az = home;
                    e.radius = Vector3Length(e.half);
                    bld.apply(e, *p->field);
                }
            const std::shared_ptr<const TerrainDeltas> dug = bld.finish();
            bool ok = true;
            {   // #693: lattice points left solid inside the cut, from 0.5 m over the floor to 0.5 m over the old ground,
                // over a 20 m square in the middle (specks of ground floating over the trench)
                int specks = 0, checked = 0;
                const int c0x = (int)lroundf(mid.x / kVoxelSize), c0y = (int)lroundf(mid.y / kVoxelSize), c0z = (int)lroundf(mid.z / kVoxelSize);
                for (int dz = -40; dz <= 40; ++dz) for (int dy = -40; dy <= 40; ++dy) for (int dx = -40; dx <= 40; ++dx) {
                    const Vector3 q = { (c0x + dx) * kVoxelSize, (c0y + dy) * kVoxelSize, (c0z + dz) * kVoxelSize };
                    const float h = Vector3DotProduct(Vector3Subtract(q, mid), home);
                    if (h < 0.5f - depth || h > 0.5f) continue;
                    ++checked;
                    if (p->field->baseDensity(q) + dug->at(c0x + dx, c0y + dy, c0z + dz) >= ISO_LEVEL) ++specks;
                }
                ok = specks == 0;
                printf("bench: dig-lod: %d of %d lattice points in the cut still solid -> %s\n", specks, checked, ok ? "PASS" : "FAIL");
            }
            for (int level = 0; level <= 6; ++level) {
                const NodeKey probe{ level, 0, 0, 0 };
                const float size = probe.size();
                const Vector3 floorMid = Vector3Add(mid, Vector3Scale(home, -depth));
                ChunkMeshData out;
                out.key = { level, (int)floorf(floorMid.x / size), (int)floorf(floorMid.y / size), (int)floorf(floorMid.z / size) };
                BuildChunkMesh(*p->field, *dug, out);
                const Vector3 org = out.key.origin();
                const float stride = out.key.stride(), inner = half - 2.0f * stride - 1.0f;
                double sum = 0.0, sum2 = 0.0, worst = 0.0; int n = 0, raw = 0;
                for (size_t v = 0; v + 2 < out.positions.size(); v += 3) {
                    const Vector3 q = Vector3Subtract(Vector3Add(org, { out.positions[v], out.positions[v + 1], out.positions[v + 2] }), mid);
                    const float x = Vector3DotProduct(q, eE), y = Vector3DotProduct(q, eN), h = Vector3DotProduct(q, home);
                    if (fabsf(x) > inner || fabsf(y) > inner) continue;
                    ++raw;
                    const double dev = h + depth;
                    if (fabs(dev) > depth + 2.0 * stride) continue; // another surface (none expected)
                    sum += dev; sum2 += dev * dev; worst = std::max(worst, fabs(dev)); ++n;
                }
                // The node's own share of the basin's interior: none when the interior is narrower than 2 strides of margin.
                const bool covers = inner > stride;
                const double mean = n ? sum / n : 0.0, rms = n ? sqrt(sum2 / n) : 0.0;
                const double bar = 0.15 + 0.05 * stride; // point sampling: 0.30 m at 2 m, 0.91 at 4 m, 2.47 at 8 m; filtered 0.04-0.08, 0.11-0.15, 0.20-0.24
                const bool pass = !covers || (n > 0 && worst <= bar);
                ok = ok && pass;
                printf("bench: dig-lod: level %d (stride %.2f m): %d floor vertices (%d inside), mean %+.3f m, rms %.3f m, worst %.3f m (bar %.2f)%s -> %s\n",
                       level, stride, n, raw, mean, rms, worst, bar, covers ? "" : " (interior narrower than the margin: not judged)", pass ? "PASS" : "FAIL");
            }
            printf("bench: dig-lod: %s\n", ok ? "PASS" : "FAIL");
            return ok ? 0 : 1;
        }
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--shot-hillshade") {
            // #678 evidence: a hillshade of the height map (PlanetField::elevation, what the minimap and the terrain show)
            // around the start, 1024 px at --hillshade-mpp m per pixel (default 1), lit from a low sun so the slope steps
            // of the bake-grid interpolation show. Writes bench_hillshade.png; prints the mean |second difference| of the
            // height (m per px^2), which grows with every slope jump.
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
            if (!pf) { printf("bench: hillshade: no planet -> FAIL\n"); return 1; }
            float mpp = 1.0f;
            for (int k = 1; k + 1 < argc; ++k) if (std::string(argv[k]) == "--hillshade-mpp") mpp = (float)atof(argv[k + 1]);
            const int N = 1024;
            const Vector3 home = pf->findSpawnDir(p->seed);
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
            std::vector<float> h((size_t)N * N);
            std::vector<std::thread> team;
            for (int t = 0; t < 8; ++t)
                team.emplace_back([&, t] {
                    for (int j = t; j < N; j += 8)
                        for (int k = 0; k < N; ++k) {
                            const float x = (k - 0.5f * N) * mpp, y = (j - 0.5f * N) * mpp;
                            h[(size_t)j * N + k] = pf->elevation(Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, x / pf->radius), Vector3Scale(t2, y / pf->radius)))), mpp);
                        }
                });
            for (auto& th : team) th.join();
            Image img = GenImageColor(N, N, BLACK);
            const Vector3 sun = Vector3Normalize({ -0.8f, 0.45f, -0.4f });
            double curv = 0.0; long long nc = 0;
            for (int j = 1; j < N - 1; ++j)
                for (int k = 1; k < N - 1; ++k) {
                    auto H = [&](int a, int b) { return h[(size_t)b * N + a]; };
                    const float gx = (H(k + 1, j) - H(k - 1, j)) / (2 * mpp), gy = (H(k, j + 1) - H(k, j - 1)) / (2 * mpp);
                    const float s = std::clamp(0.15f + 0.85f * Vector3DotProduct(Vector3Normalize({ -gx, 1.0f, -gy }), sun), 0.0f, 1.0f);
                    const float e = H(k, j);
                    const Vector3 c = e <= 0.0f ? Vector3{ 0.2f, 0.35f, 0.7f } : Vector3{ 0.85f, 0.75f, 0.55f };
                    ImageDrawPixel(&img, k, j, { (unsigned char)(255 * c.x * s), (unsigned char)(255 * c.y * s), (unsigned char)(255 * c.z * s), 255 });
                    curv += fabs(H(k + 1, j) - 2 * H(k, j) + H(k - 1, j)) + fabs(H(k, j + 1) - 2 * H(k, j) + H(k, j - 1)); nc += 2;
                }
            ExportImage(img, "bench_hillshade.png");
            UnloadImage(img);
            printf("bench: hillshade: %d px at %.2f m/px, mean |d2h| %.5f m -> bench_hillshade.png\n", N, mpp, curv / std::max<long long>(1, nc));
            return 0;
        }

    const double benchWarmup = 2.0;
    double benchStart = 0.0, benchMaxFrame = 0.0;
    std::vector<float> benchFrameTimes; // every timed frame (s): percentiles at the end (#258)
    long benchFrames = 0;
    if (!hostSave.empty()) { // --host (#311): straight into the world (no menu), listening
        gameState = 1;
        mainMenu.isActive = false;
        world.isActive = true;
        std::string status;
        Host_SetAutosaveMinutes(autosaveMinutes);
        printf("host: autosave %s\n", autosaveMinutes > 0.0 ? TextFormat("every %g minutes, and on shutdown", autosaveMinutes) : "off (still on shutdown)");
        Host_Start(hostPort, status);
        printf("host: %s\n", status.c_str());
    }
    if (bench) {
        gameState = benchMenus ? 0 : 1; // --test-menus: from the title screen
        mainMenu.isActive = benchMenus;
        world.isActive = !benchMenus;
        Planetoid* benchPlanet = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
        if (!benchAtSave && !benchCamGiven && benchPlanet && dynamic_cast<const PlanetField*>(benchPlanet->field.get())) {
            // A planet (#291): 40 m above the start, looking ahead along the surface and a little down (the asteroid
            // offset below put the camera under the terrain). It rides the planet from here (#127).
            const Vector3 feet = (player.position - world.renderOrigin).toFloat();
            Vector3 ahead = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
            if (benchUphill) // --bench-uphill (#357): face up the eroded slope, toward the cliff above
                if (const auto* pfb = dynamic_cast<const PlanetField*>(benchPlanet->field.get())) {
                    const Vector3 dl = Vector3Normalize(benchPlanet->toLocal(feet));
                    const Vector3 g = benchPlanet->dirToWorld(pfb->uphillAt(dl));
                    const Vector3 gt = Vector3Subtract(g, Vector3Scale(player.up, Vector3DotProduct(g, player.up)));
                    if (Vector3Length(gt) > 1e-4f) ahead = Vector3Normalize(gt);
                }
            world.camera.position = Vector3Add(feet, Vector3Scale(player.up, benchViewHeight)); // --bench-height (overviews, #352)
            world.camera.target = Vector3Add(world.camera.position, Vector3Subtract(Vector3Scale(ahead, std::max(100.0f, 2.5f * benchViewHeight)), Vector3Scale(player.up, std::max(25.0f, 0.8f * benchViewHeight))));
            if (benchPitch > -90.0f) { // --bench-pitch DEG: level view tilted up (> 0) or down
                const float a = benchPitch * DEG2RAD;
                world.camera.target = Vector3Add(world.camera.position, Vector3Add(Vector3Scale(ahead, 100.0f * cosf(a)), Vector3Scale(player.up, 100.0f * sinf(a))));
            }
            world.camera.up = player.up;
            if (benchViewFall)
                if (const auto* pfb = dynamic_cast<const PlanetField*>(benchPlanet->field.get())) {
                    // --view-fall (#540): from --view-fall-dist (80) m downstream of the lip and 0.3 of that up, looking back at it.
                    const DVec3 lip = DVec3(benchFallLip) * (double)(pfb->radius + pfb->elevation(benchFallLip));
                    const Vector3 at = (benchPlanet->toWorldD(lip) - world.renderOrigin).toFloat();
                    const Vector3 upW = benchPlanet->dirToWorld(benchFallLip);
                    Vector3 downW = benchPlanet->dirToWorld(benchFallDown);
                    downW = Vector3Normalize(Vector3Subtract(downW, Vector3Scale(upW, Vector3DotProduct(downW, upW))));
                    world.camera.position = Vector3Add(at, Vector3Add(Vector3Scale(downW, benchFallDist), Vector3Scale(upW, benchFallDist * 0.3f)));
                    world.camera.target = Vector3Subtract(at, Vector3Scale(upW, 10.0f));
                    world.camera.up = upW;
                }
        } else if (!benchAtSave) world.camera.position = Vector3Add(world.camera.target, benchCamOffset);
        if (benchDay) { // mid-morning over the start: 45 degrees up, from ahead-left (shadows fall toward the camera's right)
            const Vector3 ahead = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
            const Vector3 side = Vector3CrossProduct(ahead, player.up);
            sunOverride = true;
            sunArg = Vector3Normalize(Vector3Add(Vector3Scale(player.up, 0.7f), Vector3Add(Vector3Scale(ahead, 0.5f), Vector3Scale(side, -0.5f))));
        }
        if (benchSunElev > -90.0f) { // --sun-elev: the sun straight ahead at that elevation (#276)
            const Vector3 ahead = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
            const float a = benchSunElev * DEG2RAD;
            sunOverride = true;
            sunArg = Vector3Normalize(Vector3Add(Vector3Scale(player.up, sinf(a)), Vector3Scale(ahead, cosf(a))));
        }
        world.debugDraw = startWithOverlay;
        player.flyMode = benchWalkSeconds <= 0.0 && !benchWalkMode; // the bench drives the camera directly unless walking
        player.autoWalk = benchWalkMode && !benchStandMode;           // --bench-walk: walk through the measured window (--bench-stand: stand)
        player.flyFollow = benchFlySpeed <= 0.0f && benchFlyAlt >= 1e8f; // #298: --bench-fly places the camera itself every frame (the carry undid it). #127: the bench camera rides the body under it (sphere of influence), as the game's fly camera does; it
                                                                      // stood still in absolute space while the planet spun past (benches that place the camera every frame override it)
        benchStart = -1.0; // starts once the world has finished streaming in
        if (benchPreloadSpawn) { BeginSpawnPreload(); preloadT0 = GetTime(); }
        if (benchUnfocus > 0.0f && benchPreloadSpawn) { benchUnfocused = true; unfocusPhase = 1; unfocusT0 = GetTime(); } // #510: alt-tabbed during the preload
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

    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--test-newworld") {
            // The New World screen (#504) on its own: the preview of a random seed builds; then other parameters (more
            // ocean, taller mountains) rebuild it. Screenshots of both; --record captures the turning planet.
            Recorder_Start();
            NewWorld_Open();
            const double t0 = GetTime();
            double firstBuild = -1.0, secondBuild = -1.0;
            bool changed = false, shot1 = false, shot2 = false;
            while (!WindowShouldClose()) {
                const double t = GetTime() - t0;
                Recorder_Update(true);
                BeginDrawing();
                NewWorld_Draw();
                if (Recorder_WantsFrame(true)) { rlDrawRenderBatchActive(); Recorder_Capture(); } // the batch drawn first (the planet, the text)
                EndDrawing();
                Recorder_FrameDone();
                if (firstBuild < 0.0 && NewWorld_PreviewBuildSeconds() >= 0.0) firstBuild = NewWorld_PreviewBuildSeconds();
                if (firstBuild >= 0.0 && !shot1 && t > 3.0) { TakeScreenshot("bench_newworld.png"); shot1 = true; }
                if (shot1 && !changed) {
                    WorldGenParams p; p.ocean = 0.6f; p.mountains = 1.8f; p.roughness = 1.5f;
                    NewWorld_TestSetParams(p);
                    changed = true;
                }
                if (changed && secondBuild < 0.0 && NewWorld_PreviewBuildSeconds() != firstBuild && NewWorld_PreviewBuildSeconds() >= 0.0) secondBuild = NewWorld_PreviewBuildSeconds();
                if (secondBuild >= 0.0 && !shot2 && t > 7.0) { TakeScreenshot("bench_newworld_changed.png"); shot2 = true; }
                if ((shot2 && (!Recorder_Active() || Recorder_Finished())) || t > 40.0) break;
            }
            Recorder_Shutdown();
            NewWorld_Shutdown();
            // The .worldgen sidecar round-trips every field (the bake key hashes them: a lost field rebakes the wrong world).
            WorldGenParams wp; wp.ocean = 0.3f; wp.mountains = 1.7f; wp.roughness = 0.6f; wp.biomeSize = 1.3f; wp.resources = 0.8f;
            wp.erosion = true; wp.erosionStrength = 1.4f; wp.cliffs = 0.35f; wp.karst = false; wp.volcanoes = false;
            const std::string sidecar = std::string(GetApplicationDirectory()) + "newworld_test";
            const bool saved = WorldGenParams_Save(sidecar, wp);
            const WorldGenParams back = WorldGenParams_ForSave(sidecar);
            remove((sidecar + ".worldgen").c_str());
            const bool roundTrip = saved && back.key() == wp.key() && back.cliffs == wp.cliffs && !back.karst && !back.volcanoes;
            printf("bench: newworld: world settings round-trip through the .worldgen sidecar: %s\n", roundTrip ? "yes" : "NO");
            const bool pass = firstBuild >= 0.0 && secondBuild >= 0.0 && shot1 && shot2 && roundTrip;
            printf("bench: newworld: preview built in %.2f s, rebuilt after a change in %.2f s; bench_newworld.png, bench_newworld_changed.png -> %s\n",
                   firstBuild, secondBuild, pass ? "PASS" : "FAIL");
            CloseWindow();
            return pass ? 0 : 1;
        }
    for (int i = 1; i < argc; ++i)
        if (std::string(argv[i]) == "--shot-settings") {
            // The settings page (#512) in both menus, every tab: a screenshot each, and every row and the Back button on
            // screen (rows above Back). Run it at more than one size: --window 1280 720.
            const char* tabs[] = { "graphics", "controls", "audio", "gameplay" };
            bool pass = true;
            for (int menu = 0; menu < 2; ++menu)
                for (int tab = 0; tab < 4; ++tab) {
                    Menu_TestSettings(tab, menu == 1);
                    for (int f = 0; f < 3 && !WindowShouldClose(); ++f) { // a few frames: the layout settles, the shot is of a drawn one
                        BeginDrawing();
                        ClearBackground(BLACK);
                        if (menu == 0) Menu_Main(); else Menu_Pause();
                        EndDrawing();
                    }
                    float bottom = 0.0f; Rectangle back{};
                    Menu_SettingsLayout(bottom, back);
                    const bool ok = bottom <= back.y && back.y + back.height <= GetScreenHeight() && back.y >= 0.0f;
                    pass &= ok;
                    const char* file = TextFormat("bench_settings_%s_%s.png", menu == 0 ? "main" : "pause", tabs[tab]);
                    TakeScreenshot(file);
                    printf("bench: settings: %s %s (%dx%d): rows end at %.0f, Back %.0f-%.0f -> %s\n", menu == 0 ? "main" : "pause", tabs[tab],
                           GetScreenWidth(), GetScreenHeight(), bottom, back.y, back.y + back.height, ok ? "ok" : "OFF SCREEN");
                }
            printf("bench: settings (#512): every tab fits with Back on screen -> %s\n", pass ? "PASS" : "FAIL");
            CloseWindow();
            return pass ? 0 : 1;
        }
    Recorder_Start(); // --record (#375): fixed frame time, seeded random values, the capture folder
    // --play-new S (#387): the menus as a player uses them, for --perf-log. 10 s on the title screen; New Game (the New World
    // screen creates after a moment) -> bake -> loading -> S s of play; Save & Quit -> 10 s on the title screen; Load that
    // save -> loading -> 30 s of play; exit.
    int playStep = 0;
    double playStepAt = 0.0;
    std::string playSave;
    double playLoadStart = -2.0; // --play-load (#509): -2 load not issued yet, -1 loading, else the walk's start
    while (windowKeepAlive && !WindowShouldClose())
    {
        if (Recorder_Finished()) break; // a --record capture is complete
        if (playNewSeconds > 0.0) {
            const double now = GetTime();
            if (playStep == 0 && gameState == 0 && now > 10.0) { pendingMenu.action = MenuAction::NewGame; playStep = 1; }
            else if (playStep == 1 && gameState == 1 && !loading) { playStep = 2; playStepAt = now; }
            else if (playStep == 2 && now - playStepAt > playNewSeconds) { playSave = WorldSavePath(); pendingMenu.action = MenuAction::SaveAndQuit; playStep = 3; }
            else if (playStep == 3 && gameState == 0) { playStep = 4; playStepAt = now; }
            else if (playStep == 4 && now - playStepAt > 10.0) { pendingMenu.action = MenuAction::LoadGame; pendingMenu.path = playSave; playStep = 5; }
            else if (playStep == 5 && gameState == 1 && !loading) { playStep = 6; playStepAt = now; }
            else if (playStep == 6 && now - playStepAt > 30.0) break;
        }
        // --play-load (#509): the playtest's freezes came walking ~6 m/s through a loaded world while looking around fast
        // (18-37 px/frame). Walk forward; the view sweeps +-90 degrees every 4 s so the walk still goes somewhere.
        if (playLoadSeconds > 0.0 && gameState == 0 && pendingMenu.action == MenuAction::None && playLoadStart < -1.5) {
            pendingMenu.action = MenuAction::LoadGame; pendingMenu.path = playLoadPath; playLoadStart = -1.0;
        }
        if (playLoadSeconds > 0.0 && gameState == 1) {
            if (playLoadStart < 0.0 && !loading) playLoadStart = GetTime();
            if (playLoadStart >= 0.0 && GetTime() - playLoadStart > playLoadSeconds) break;
            const double t = playLoadStart >= 0.0 ? GetTime() - playLoadStart : 0.0;
            gFakeHeld = loading ? 0 : (1ull << (int)Act::Forward) | (fmod(t, 1.5) < 0.1 ? 1ull << (int)Act::Jump : 0ull); // a hop every 1.5 s: a scripted walk gets stuck on rough ground
            gFakeLook = loading ? Vector2{ 0.0f, 0.0f } : Vector2{ (float)(60.0 * sin(t * PI * 0.5)), 0.0f };
        }
        Recorder_Update(WorldIdle() && !loading); // its start delay (after the loading screen), F10 in normal play
        GpuTimer_CpuFrameStart(!loading); GpuTimer_BeginFrame();
        FpsLog_Frame(gameState == 0 ? "menu" : gameState == 3 ? "newworld" : gameState == 2 ? "paused" : loading ? "loading" : "play");
        // --bench-unfocus S (#510): 1 unfocused during the preload for S s; then, playing: 2 measure 5 s, 3 unfocused S s,
        // 4 settle 2 s, 5 measure 5 s, then FPS before/after.
        if (benchUnfocus > 0.0f) {
            const double ut = GetTime() - unfocusT0;
            if (unfocusPhase == 1 && ut > benchUnfocus) { benchUnfocused = false; printf("bench: unfocus: refocused after %.1f s unfocused during the preload\n", ut); unfocusPhase = 2; unfocusT0 = -1.0; }
            else if (unfocusPhase == 0 && !benchPreloadSpawn) unfocusPhase = 2, unfocusT0 = -1.0;
            else if (unfocusPhase == 2 && preloadT0 < 0.0 && WorldIdle()) {
                if (unfocusT0 < 0.0) { unfocusT0 = GetTime(); unfocusFrames = 0; } // the measurement starts
                else {
                    ++unfocusFrames;
                    if (ut > 5.0) { unfocusFpsBefore = unfocusFrames / ut; benchUnfocused = true; unfocusPhase = 3; unfocusT0 = GetTime(); unfocusSim0 = SimTime(); }
                }
            }
            else if (unfocusPhase == 3 && ut > benchUnfocus) {
                // #525: the game keeps running behind another window (time of day, physics and the factory all step on it).
                const double simRan = SimTime() - unfocusSim0;
                unfocusSimOk = simRan > 0.9 * ut;
                printf("bench: unfocus: the simulation ran %.2f s of %.2f s unfocused (#525): %s\n", simRan, ut, unfocusSimOk ? "yes" : "NO");
                benchUnfocused = false; unfocusPhase = 4; unfocusT0 = GetTime();
            }
            else if (unfocusPhase == 4 && ut > 2.0) { unfocusPhase = 5; unfocusT0 = GetTime(); unfocusFrames = 0; }
            else if (unfocusPhase == 5) { ++unfocusFrames; if (ut > 5.0) { printf("bench: unfocus: %.1f FPS before %.0f s unfocused while playing, %.1f FPS after -> %s\n", unfocusFpsBefore, benchUnfocus, unfocusFrames / ut, unfocusSimOk ? "PASS" : "FAIL"); unfocusPhase = 6; } }
        }
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
                        {   // #260: bushes and stones are picked on their model: from 2 m off, a ray at the middle of the
                            // nearest of each hits it, one 0.3 m over its top misses it.
                            std::vector<Plant> plants;
                            p->trees->plantsNear(lf, 300.0f, plants);
                            int checked = 0, ok = 0;
                            for (int want = 0; want < 2; ++want) { // 0 bush, 1 stone
                                const Plant* best = nullptr; float bd = 1e9f;
                                for (const Plant& q : plants) {
                                    if (q.tree() || q.stone() != (want == 1)) continue;
                                    const float d = Vector3Distance(q.base, lf);
                                    if (d < bd) { bd = d; best = &q; }
                                }
                                if (!best) continue;
                                const Vector3 side = Vector3Normalize(Vector3CrossProduct(best->up, Vector3{ 0.3f, 1.0f, 0.2f }));
                                const float top = best->height; // bushes ~1.2 x scale, stones ~0.5 x scale
                                for (int k = 0; k < 2; ++k) {
                                    const Vector3 target = Vector3Add(best->base, Vector3Scale(best->up, k == 0 ? top * 0.3f : top * 1.4f + 0.3f));
                                    const Vector3 from = Vector3Add(target, Vector3Scale(side, 2.0f));
                                    PlantHit h;
                                    const bool hit = p->trees->raycast(from, Vector3Normalize(Vector3Subtract(target, from)), 4.0f, h) && h.plant.id == best->id;
                                    ++checked; ok += (hit == (k == 0)) ? 1 : 0;
                                }
                            }
                            printf("bench: chop: bush and stone picking (#260): %d of %d rays as expected -> %s\n", ok, checked, ok == checked && checked > 0 ? "PASS" : "FAIL");
                        }
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
                        TerrainTool_Update(world, PlayerInput::None(), true); // aim only
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
                    } else if (chopPhase == 4 && GetTime() - chopStart > 1.0) {
                        // The tiles rebuild on worker threads: poll up to 20 s (a fixed 3 s missed it under load, #465).
                        const bool back = present(chopId);
                        if (back || GetTime() - chopStart >= 20.0) {
                            printf("bench: chop: %.0f min later it has regrown: %s (after %.1f s) -> %s\n", (kRegrowDelay + kRegrowTime) / 60.0, back ? "yes" : "no",
                                   GetTime() - chopStart, chopFelled && chopSaved && back ? "PASS" : "FAIL");
                            benchChop = false;
                        }
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
                        gBuildCensus = PartCensus(p);
                        // Stand on the slab: drop the player 2 m above a free corner of it.
                        const Vector3 top = p->toWorld(plot->toBody({ -2.75f, 0.5f + 2.0f, 1.5f })); // over the pad
                        player.placeAt(world, *p, top);
                        player.velocity = p->surfaceVelocity(top);
                        buildPhase = 2;
                        buildStart = -1.0; // the 2 s count starts when play resumes: placeAt raises the loading screen (it took 2.4 s once)
                    } else if (buildPhase == 2 && buildStart < 0.0) {
                        buildStart = GetTime();
                    } else if (buildPhase == 2 && GetTime() - buildStart > 2.0) {
                        Plot* plot = buildPlot;
                        const Vector3 lf = plot->toPlot(p->toLocal((player.position - world.renderOrigin).toFloat()));
                        const bool standing = player.grounded && fabsf(lf.y - 0.5f) < 0.15f;
                        printf("bench: build: the player stands on the slab: %s (feet %.2f m above the plot floor, grounded %d)\n", standing ? "yes" : "NO", lf.y, player.grounded);
                        // Every live part (plot, type, cell, orientation), to name what a save/load adds or loses (#602).
                        // Parts that appeared during the 2 s of play (#602), then the save/load round trip itself.
                        const PartCensusMap before = PartCensus(p);
                        const int appeared = PartCensusDiff(gBuildCensus, before, "appeared during play", 5);
                        const bool saved = SaveWorld(world), loaded = LoadWorldSave(world);
                        p = dynamic_cast<Planetoid*>(world.objects[0].get());
                        const PartCensusMap now = PartCensus(p);
                        const int added = PartCensusDiff(before, now, "added by the load", 5), lost = PartCensusDiff(now, before, "lost by the load", 5);
                        const int after = p->plots ? p->plots->partCount() : 0;
                        printf("bench: build: save %s, load %s: %d placed, %d appeared during play, %d added and %d lost by the load (%d -> %d parts)\n",
                               saved ? "ok" : "FAILED", loaded ? "ok" : "FAILED", buildPlaced, appeared, added, lost, buildPlaced, after);
                        buildOk = standing && saved && loaded && added == 0 && lost == 0 && appeared == 0;
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
                        Factory::Cable(p->id, factPlot, Int3{ 0, 0, -1 }); // and the first group's cable (#570)
                        factPhase = 4; factStart = GetTime(); factSamples = 0; factSum = 0.0; factMax = 0.0;
                    } else if (factPhase == 4) {
                        const FactoryStats s = Factory::Stats();
                        factSum += s.avgTickMs; ++factSamples;
                        if (GetTime() - factStart > 12.0) {
                            const double avg = factSum / std::max(1, factSamples);
                            const MachineView chest = Factory::View(p->id, factPlot, factSeeds[0][2] + 1);
                            const CableView cable = Factory::Cable(p->id, factPlot, Int3{ 0, 0, -1 });
                            const bool cableOk = cable.valid && cable.voltage > 0 && cable.generators == 1 && cable.consumers == 2 && cable.supply > 0.0f;
                            printf("bench: factory: the first cable's network (#570): %d V, %d generator, %d machines, %.1f of %.1f EU/t -> %s\n", cable.voltage, cable.generators, cable.consumers, cable.supply, cable.demand, cableOk ? "ok" : "WRONG");
                            const bool pass = avg < 5.0 && s.activeMachines >= benchFactory / 2 && s.itemsMade > 0 && cableOk;
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
                        // Wait for the reloaded layout to reach the factory (10,000 machines: no fixed time; #495), then 1 s more.
                        static double validAt = -1.0;
                        if (chest.valid && validAt < 0.0) { validAt = GetTime(); printf("bench: factory: the chest is back in the factory %.2f s after the reload\n", validAt - factStart); }
                        if ((validAt >= 0.0 && GetTime() - validAt > 1.0) || GetTime() - factStart > 30.0) {
                            printf("bench: factory: first chest before save %s, after reload %s -> %s\n", factChestBefore.c_str(), chest.input.c_str(),
                                   chest.valid && chest.input != "-" && !chest.input.empty() ? "PASS" : "FAIL");
                            // Then remove that chest: what it held must come back to the inventory.
                            factRefundBefore = inventory.count("iron_ingot") + inventory.count("iron_plate"); // the furnaces make ingots since #212 (#434)
                            if (Plot* pl = p2->plots ? p2->plots->byId(factPlot) : nullptr) pl->remove(factSeeds[0][2] + 1);
                            factPhase = 6; factStart = GetTime();
                        }
                    } else if (factPhase == 6 && GetTime() - factStart > 2.0) {
                        const int got = inventory.count("iron_ingot") + inventory.count("iron_plate") - factRefundBefore;
                        printf("bench: factory: removed the chest: %d iron ingots or plates back in the inventory -> %s\n", got, got > 0 ? "PASS" : "FAIL");
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
                } else if (benchDoorway && (doorPhase > 0 || WorldIdle())) {
                    // --test-doorway (#255, user 2026-11-14): the player is 1.7 m tall. Two walls on a floor slab, each with a
                    // 1.5 m wide opening: 2.0 m high (must walk through) and 1.5 m high (must stop). Walk straight at each.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const double t = GetTime() - doorStart;
                    Plot* plot = p->plots ? p->plots->byId(doorPlot) : nullptr;
                    const float laneZ[2] = { 2.25f, 8.25f }; // plot metres: the openings' centres (cells 3-5 and 15-17)
                    constexpr int kDoorY = 60; constexpr float kDoorFloor = (kDoorY + 1) * 0.5f; // the slab 30 m up: no ground or trees in the lanes
                    auto placeInLane = [&](int lane) {
                        const Vector3 at = p->toWorld(plot->toBody({ 1.5f, kDoorFloor + 0.3f, laneZ[lane] }));
                        player.placeAt(world, *p, at); // sets the velocity to the surface's there (render coordinates)
                        player.forward = Vector3Normalize(Vector3Subtract(p->toWorld(plot->toBody({ 2.0f, kDoorFloor, laneZ[lane] })), p->toWorld(plot->toBody({ 1.0f, kDoorFloor, laneZ[lane] }))));
                    };
                    if (doorPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        doorPhase = 1; doorStart = GetTime();
                    } else if (doorPhase == 1 && t > 1.0 && player.grounded) {
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        TerrainHit th;
                        RaycastTerrain(world, Ray{ Vector3Add(Vector3Add(feet, Vector3Scale(flat, 10.0f)), Vector3Scale(player.up, 20.0f)), Vector3Negate(player.up) }, 60.0f, th);
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        plot = p->plots->plotFor(th.local, 0.0f, true);
                        doorPlot = plot->id;
                        int placed = 0;
                        auto put = [&](Int3 c) { if (BuildTool_Place(*p, *plot, { PartType::Cube, { c.x, c.y + kDoorY, c.z }, 0 }, false)) ++placed; };
                        for (int z = 0; z < 21; ++z) for (int x = 0; x < 16; ++x) put({ x, 0, z });     // floor slab, 8 x 10.5 m (cells relative to the slab row)
                        // Walls at x cell 8 (4.0-4.5 m), 4 m high; the openings: z cells 3-5 up to y 4 (2.0 m) and 15-17 up to y 3 (1.5 m).
                        for (int y = 1; y <= 8; ++y)
                            for (int z = 0; z < 21; ++z) {
                                if (z >= 3 && z <= 5 && y <= 4) continue;
                                if (z >= 15 && z <= 17 && y <= 3) continue;
                                put({ 8, y, z });
                            }
                        for (int y = 1; y <= 8; ++y) for (int x = 0; x < 16; ++x) put({ x, y, 10 }); // a divider between the lanes
                        for (int y = 1; y <= 4; ++y) {                                                   // rails: nobody walks off 30 m up
                            for (int x = 0; x < 16; ++x) { put({ x, y, 0 }); put({ x, y, 20 }); }
                            for (int z = 0; z < 21; ++z) { put({ 0, y, z }); put({ 15, y, z }); }
                        }
                        printf("bench: doorway: %d parts placed on plot %u\n", placed, plot->id);
                        doorLane = 0; placeInLane(0);
                        doorPhase = 2; doorStart = -1.0; // the clock starts when play resumes: placeAt raises the loading screen
                    } else if (doorPhase == 2 && doorStart < 0.0) {
                        doorStart = GetTime();
                    } else if (doorPhase == 2 && t > 1.0) {
                        const Vector3 lf = plot->toPlot(p->toLocal((player.position - world.renderOrigin).toFloat()));
                        // Walk the lane; stop once past the wall (pressing the far rail would only bank impact speed).
                        player.autoWalk = lf.x < 6.0f && t < 7.0;
                        player.forward = Vector3Normalize(Vector3Subtract(p->toWorld(plot->toBody({ 2.0f, kDoorFloor, laneZ[doorLane] })), p->toWorld(plot->toBody({ 1.0f, kDoorFloor, laneZ[doorLane] }))));
                        if (getenv("DOOR_DEBUG")) {
                            static float lastHealth = 100.0f;
                            if (survival.health < lastHealth - 0.01f)
                                printf("DOOR hurt t=%.2f: health %.0f -> %.0f, landing %.2f m/s, impact %.2f m/s, grounded %d, y %.2f, v %.2f along up\n", t, lastHealth, survival.health,
                                       player.landingSpeed, player.impactSpeed, player.grounded ? 1 : 0, lf.y - kDoorFloor, Vector3DotProduct(player.velocity, player.up));
                            lastHealth = survival.health;
                            static int lastTick = -1;
                            if ((int)(t * 2.0) != lastTick) {
                                lastTick = (int)(t * 2.0);
                                printf("DOOR lane %d t=%.1f: x %.2f y %.2f z %.2f grounded %d health %.0f deaths %d\n", doorLane, t, lf.x, lf.y - kDoorFloor, lf.z, player.grounded ? 1 : 0, survival.health, survival.deaths);
                            }
                        }
                        if (t > 7.5) { // half a second standing still before the verdict and the next placement
                            const bool through = lf.x > 5.5f;
                            doorPass[doorLane] = survival.deaths == 0 && (doorLane == 0 ? through : (!through && lf.x < 4.0f));
                            printf("bench: doorway: %s opening: the player ended %.2f m along the lane (the wall is at 4.0-4.5 m), feet %.2f m over the slab top, %.2f m across (opening centre %.2f), deaths %d: %s -> %s\n",
                                   doorLane == 0 ? "2.0 m" : "1.5 m", lf.x, lf.y - kDoorFloor, lf.z, laneZ[doorLane], survival.deaths, through ? "walked through" : "stopped", doorPass[doorLane] ? "PASS" : "FAIL");
                            if (doorLane == 0) { doorLane = 1; placeInLane(1); doorStart = -1.0; }
                            else {
                                printf("bench: doorway: player %.2f m tall: %s\n", player.height, doorPass[0] && doorPass[1] ? "PASS" : "FAIL");
                                benchDoorway = false;
                            }
                        }
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
                        wallMinHealth = std::min(wallMinHealth, survival.health + survival.fallDamageTaken); // hard landings in the shaft are fall damage, not the wall (#439)
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
                        printf("bench: wall shove 9 s: largest frame move %.2f m (at %.1f s), top speed vs ground %.1f m/s, worst impact %.1f m/s, health %.0f not counting fall damage (%.0f so far), deaths %d: %s\n",
                               wallMaxStep, wallMaxStepAt, wallMaxSpeed, wallMaxImpact, wallMinHealth, survival.fallDamageTaken, survival.deaths, pass ? "PASS" : "FAIL");
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
                        wallMaxStep = 0.0; wallMaxImpact = 0.0f; wallMaxSpeed = 0.0f; wallMinHealth = survival.health + survival.fallDamageTaken;
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
                } else if (benchHover && (hoverPhase > 0 || WorldIdle())) {
                    // --test-hover (#163): with the flight pack fitted, a drop from 100 m is caught by the fall safety (no
                    // damage), then hovering holds the position (no drift) and the height.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    static DVec3 held{}; static double heldT = 0.0; static float ground = 0.0f;
                    const double t = GetTime() - hoverStart;
                    if (hoverPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        inventory.add("flight_pack", 1); inventory.add("battery_pack", 1); // flight needs suit power (#217)
                        Suit_Reset(); Suit_Fit("flight_pack"); Suit_Fit("battery_pack");
                        hoverPhase = 1; hoverStart = GetTime();
                    } else if (hoverPhase == 1 && t > 1.5 && player.grounded) {
                        survival.health = 100.0f; survival.dead = false; survival.power = 100.0f;
                        ground = (float)p->toLocalD(player.position).length();
                        player.position = player.position + DVec3(Vector3Scale(player.up, 100.0f));
                        hoverPhase = 2; hoverStart = GetTime();
                    } else if (hoverPhase == 2) {
                        // Measured from 1.5 s after the catch: the teleport up leaves a grounded (surface-relative) velocity read as
                        // world-frame (~84 m/s of spin), which the hover brakes first.
                        static double caughtT = -1.0;
                        if (player.hovering && caughtT < 0.0) caughtT = GetTime();
                        if (caughtT > 0.0 && GetTime() - caughtT > 1.5) { held = p->toLocalD(player.position); heldT = GetTime(); hoverPhase = 3; }
                        else if ((player.grounded && t > 0.1) || survival.dead || t > 20.0) { // (t: a frame may run no sim step, #315)
                            printf("bench: hover: the fall was not caught (grounded %d, dead %d, health %.0f) -> FAIL\n", player.grounded, survival.dead, survival.health);
                            benchHover = false;
                        }
                    } else if (hoverPhase == 3 && GetTime() - heldT <= 3.0) {
                        static double lastP = 0.0;
                        if (GetTime() - lastP > 0.5) {
                            lastP = GetTime();
                            const Vector3 at = (player.position - world.renderOrigin).toFloat();
                            const Vector3 rel = Vector3Subtract(player.velocity, p->surfaceVelocity(at));
                            const float vUp = Vector3DotProduct(rel, player.up);
                            printf("bench: hover dbg: rel speed %.2f m/s (up %.2f), drift so far %.2f m, grounded %d\n", Vector3Length(rel), vUp, (float)(p->toLocalD(player.position) - held).length(), player.grounded);
                        }
                    } else if (hoverPhase == 3 && GetTime() - heldT > 3.0) {
                        const DVec3 now = p->toLocalD(player.position);
                        const float drift = (float)(now - held).length(), height = (float)now.length() - ground;
                        const bool pass = player.hovering && !player.grounded && survival.health > 99.0f && drift < 0.5f && height > 5.0f;
                        printf("bench: hover: caught the 100 m fall, health %.0f; hovering %.0f m up, drifted %.2f m in 3 s, suit power %.1f -> %s\n",
                               survival.health, height, drift, survival.power, pass ? "PASS" : "FAIL");
                        player.hovering = false; benchHover = false;
                    }
                } else if (benchTwoPlayers && (twoPhase > 0 || WorldIdle())) {
                    // --test-two-players (#314): a second Player 2 km away walks at the same time as the local one, each stepped
                    // in its own frame (RenderOriginScope) on the sim clock; the camera stays with the local player. Both must
                    // walk at about walking speed, mostly grounded, never inside solid ground.
                    static Player other;
                    static long long otherTick = 0;
                    static DVec3 start0{}, start1{};
                    static double twoStart = 0.0;
                    static int frames = 0, grounded0 = 0, grounded1 = 0;
                    static float worst0 = 0.0f, worst1 = 0.0f, twoDistUsed = 2000.0f;
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    auto inside = [&](const Player& pl) { // body density at the capsule's middle
                        const DVec3 mid = pl.position + DVec3(Vector3Scale(pl.up, 0.9f));
                        return p->densityAt(p->toLocalD(mid).toFloat());
                    };
                    if (twoPhase == 0) {
                        const Vector3 d0 = StartSpawnDir(*p);
                        const char* env = getenv("RAYLIBGAME_TWO_DIST"); // bench knob: the second player's distance (m)
                        twoDistUsed = env ? (float)atof(env) : 2000.0f;
                        const float twoDist = twoDistUsed;
                        const Vector3 d1 = LandDirAt(*p, d0, twoDist); // over dry land
                        player.spawnOn(world, *p, d0);
                        other.spawnOn(world, *p, d1);
                        otherTick = SimTick();
                        twoPhase = 1; twoStart = -1.0;
                    } else {
                        // The second player: one update per simulation step since the last frame, in its own frame.
                        PlayerInput walk;
                        if (twoPhase == 2) walk.hold(Act::Forward);
                        for (; otherTick < SimTick(); ++otherTick) {
                            RenderOriginScope frame(world, other.position);
                            other.update(world, walk, SimDt());
                        }
                        if (twoPhase == 1) {
                            if (twoStart < 0.0) twoStart = GetTime();
                            if (GetTime() - twoStart > 2.0 && player.grounded && other.grounded) {
                                start0 = p->toLocalD(player.position); start1 = p->toLocalD(other.position);
                                player.autoWalk = true;
                                twoPhase = 2; twoStart = GetTime();
                            }
                        } else if (twoPhase == 2) {
                            ++frames; grounded0 += player.grounded ? 1 : 0; grounded1 += other.grounded ? 1 : 0;
                            worst0 = std::max(worst0, inside(player)); worst1 = std::max(worst1, inside(other));
                            if (GetTime() - twoStart > 5.0) {
                                player.autoWalk = false;
                                const double w0 = (p->toLocalD(player.position) - start0).length() * p->scale;
                                const double w1 = (p->toLocalD(other.position) - start1).length() * p->scale;
                                const double apart = (other.position - player.position).length();
                                const float g0 = 100.0f * grounded0 / std::max(1, frames), g1 = 100.0f * grounded1 / std::max(1, frames);
                                const bool pass = w0 > 15.0 && w1 > 15.0 && g0 > 50.0f && g1 > 50.0f && worst0 < ISO_LEVEL && worst1 < ISO_LEVEL && apart > 0.75 * fabs(twoDistUsed);
                                printf("bench: two players %.0f m apart: local walked %.1f m (grounded %.0f%%, densest %.2f), second walked %.1f m (grounded %.0f%%, densest %.2f) in 5 s -> %s\n",
                                       apart, w0, g0, worst0, w1, g1, worst1, pass ? "PASS" : "FAIL");
                                benchTwoPlayers = false;
                            }
                        }
                    }
                } else if (benchInterest && (interestPhase > 0 || WorldIdle())) {
                    // --test-interest (#318): a second interest centre 2 km from the player (a remote player on a server). Tree
                    // tiles load around it, and animals placed there are simulated (they move); herds may spawn there too.
                    static Vector3 centre{};
                    static int first = 0, placed = 0, before = 0, maxSim = 0;
                    static std::vector<Vector3> start;
                    static double t0 = 0.0;
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (interestPhase == 0) {
                        const Vector3 d0 = StartSpawnDir(*p);
                        const Vector3 d1 = LandDirAt(*p, d0, 3000.0f); // over dry land, past the camera's 2 km tree range
                        player.spawnOn(world, *p, d0);
                        Player probe; // finds the ground there
                        probe.spawnOn(world, *p, d1);
                        centre = p->toLocalD(probe.position).toFloat();
                        if (!getenv("RAYLIBGAME_NO_INTEREST")) gInterestCentres = { InterestCentre{ p->id, centre } }; // (the A/B: without it)
                        before = Animals_Count();
                        first = Animals_Count();
                        placed = Animals_Spawn(world, *p, "woolback", centre, 3, false, 6.0f);
                        start.clear();
                        for (int i = first; i < first + placed; ++i) start.push_back(Animals_Position(i));
                        interestPhase = 1; t0 = -1.0;
                    } else if (interestPhase == 1) {
                        if (t0 < 0.0) t0 = GetTime();
                        maxSim = std::max(maxSim, Animals_SimulatedNear(*p, centre, 50.0f)); // the ones placed there, simulated
                        if (GetTime() - t0 > 8.0) {
                            float moved = 0.0f;
                            for (int k = 0; k < placed && first + k < Animals_Count(); ++k) moved = std::max(moved, Vector3Distance(Animals_Position(first + k), start[k]));
                            int nearCentre = 0;
                            for (int i = 0; i < Animals_Count(); ++i) nearCentre += Vector3Distance(Animals_Position(i), centre) < 300.0f ? 1 : 0;
                            const int tiles = p->trees ? p->trees->readyTilesNear(centre, 128.0f) : 0;
                            const bool pass = placed > 0 && maxSim >= placed && tiles > 0;
                            printf("bench: interest: a centre 3 km away: %d tree tiles built around it, %d animals placed there (%d of them simulated, moved up to %.1f m), %d animals within 300 m of it (%d before) -> %s\n",
                                   tiles, placed, maxSim, moved, nearCentre, before, pass ? "PASS" : "FAIL");
                            gInterestCentres.clear();
                            benchInterest = false;
                        }
                    }
                } else if (benchFingerprint && WorldIdle()) {
                    // --test-fingerprint (#316): the world fingerprint is the same for the same seed (regenerated) and differs for
                    // another seed; the gameplay RNG replays from a seed and continues over a save and load.
                    benchFingerprint = false;
                    const unsigned seedWas = worldGen.seed;
                    const uint64_t fp1 = WorldFingerprint(world), fp1b = WorldFingerprint(world);
                    regenerateWorld(world);
                    const uint64_t fp2 = WorldFingerprint(world);
                    worldGen.seed = worldGen.lastSeed + 1;
                    regenerateWorld(world);
                    const uint64_t fp3 = WorldFingerprint(world);
                    worldGen.seed = seedWas;
                    regenerateWorld(world);
                    const uint64_t fp4 = WorldFingerprint(world);
                    // The RNG: a seed replays; a save and load continue the stream where it was.
                    GameRng_Seed(1234); int a[8]; for (int& v : a) v = GameRng_Int(0, 1000000);
                    GameRng_Seed(1234); bool replay = true; for (int v : a) replay &= GameRng_Int(0, 1000000) == v;
                    for (int k = 0; k < 5; ++k) GameRng_Int(0, 9);
                    const uint64_t before = GameRng_State();
                    const bool saved = SaveWorld(world) && LoadWorldSave(world);
                    const bool kept = GameRng_State() == before;
                    // The join handshake's world description round-trips and carries this world's fingerprint.
                    const JoinInfo sent = JoinInfo_FromWorld(world);
                    JoinInfo got; std::string why;
                    const bool join = JoinInfo_Read(JoinInfo_Write(sent), got, why) && got.seed == sent.seed && got.bodies == sent.bodies &&
                                      got.fingerprint == sent.fingerprint && got.fingerprint == WorldFingerprint(world) && !got.bodies.empty();
                    if (!join) printf("bench: fingerprint: join info: %s\n", why.c_str());
                    printf("bench: fingerprint: join info %zu bytes, %zu bodies, round trip %s\n", JoinInfo_Write(sent).size(), sent.bodies.size(), join ? "yes" : "no");
                    const bool pass = fp1 == fp1b && fp1 == fp2 && fp3 != fp1 && fp4 == fp1 && replay && saved && kept && join;
                    printf("bench: fingerprint: %016llx, same after regenerating %s, another seed differs %s (%016llx), back %s; game RNG replays %s, kept over save/load %s -> %s\n",
                           (unsigned long long)fp1, fp1 == fp2 && fp1 == fp1b ? "yes" : "no", fp3 != fp1 ? "yes" : "no", (unsigned long long)fp3,
                           fp4 == fp1 ? "yes" : "no", replay ? "yes" : "no", saved && kept ? "yes" : "no", pass ? "PASS" : "FAIL");
                } else if (benchNet && WorldIdle()) {
                    // --test-net (#310): host on a loopback port and connect to it from the same process; a reliable message
                    // goes out and its echo comes back. SKIP when built without ENABLE_GNS.
                    benchNet = false;
                    std::string err;
                    if (!Net_Available()) printf("bench: net: built without ENABLE_GNS -> SKIP\n");
                    else if (!Net_Init(err)) printf("bench: net: init failed (%s) -> FAIL\n", err.c_str());
                    else {
                        const uint16_t port = 27315;
                        const bool hosting = Net_Host(port);
                        const uint32_t client = hosting ? Net_Connect(TextFormat("127.0.0.1:%d", port)) : 0;
                        bool sent = false, echoed = false, served = false;
                        const double t0 = GetTime();
                        while (client && GetTime() - t0 < 5.0 && !echoed) {
                            Net_Poll();
                            if (!sent && Net_IsConnected(client)) { const char hi[] = "hello"; sent = Net_Send(client, hi, sizeof(hi), true); }
                            std::vector<NetMessage> in;
                            Net_Receive(in);
                            for (const NetMessage& m : in) {
                                const std::string text(m.data.begin(), m.data.end());
                                if (m.conn != client && text.rfind("hello", 0) == 0) { const char back[] = "echo"; served = Net_Send(m.conn, back, sizeof(back), true); }
                                else if (m.conn == client && text.rfind("echo", 0) == 0) echoed = true;
                            }
                            std::this_thread::sleep_for(std::chrono::milliseconds(5));
                        }
                        const int conns = Net_ConnectionCount();
                        printf("bench: net: listen on %d %s, client connected %s, server echoed %s, echo back %s in %.0f ms (%d connections) -> %s\n", port,
                               hosting ? "yes" : "no", sent ? "yes" : "no", served ? "yes" : "no", echoed ? "yes" : "no", (GetTime() - t0) * 1000.0, conns,
                               hosting && sent && served && echoed ? "PASS" : "FAIL");
                        Net_Shutdown();
                    }
                } else if (benchNetIds && WorldIdle()) {
                    // --test-netids (#317): bodies carry their creation-order ids; animals get allocated ids that survive a save
                    // and load unchanged; allocation after the load never reuses one; a body is found by its id again.
                    benchNetIds = false;
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    int bodies = 0, bodyOk = 0;
                    std::set<uint64_t> bodyIds;
                    for (auto& obj : world.objects)
                        if (auto* q = dynamic_cast<Planetoid*>(obj.get()); q && q->ringRock < 0 && !q->vehicle) {
                            const bool ok = NetId_Kind(q->netId) == NetKind::Body && NetId_Value(q->netId) == (uint64_t)bodies && Planetoid::findNet(q->netId) == q;
                            bodyOk += ok ? 1 : 0; ++bodies; bodyIds.insert(q->netId);
                        }
                    player.spawnOn(world, *p, StartSpawnDir(*p));
                    const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                    Animals_Spawn(world, *p, "woolback", me, 4);
                    std::map<uint64_t, std::string> before; // id -> species/body (the order may change on a load)
                    for (int i = 0; i < Animals_Count(); ++i) before[Animals_NetId(i)] = Animals_Body(i);
                    const bool unique = (int)before.size() == Animals_Count() && !before.count(0);
                    const bool saved = SaveWorld(world) && LoadWorldSave(world);
                    std::map<uint64_t, std::string> after;
                    for (int i = 0; i < Animals_Count(); ++i) after[Animals_NetId(i)] = Animals_Body(i);
                    std::set<uint64_t> bodyIdsAfter;
                    for (auto& obj : world.objects)
                        if (auto* q = dynamic_cast<Planetoid*>(obj.get()); q && q->ringRock < 0 && !q->vehicle) bodyIdsAfter.insert(q->netId);
                    const NetId fresh = NetId_Allocate(NetKind::Animal);
                    const bool freshNew = !after.count(fresh) && NetId_Kind(fresh) == NetKind::Animal;
                    const bool pass = bodies > 0 && bodyOk == bodies && bodyIds == bodyIdsAfter && unique && saved && Animals_Count() > 0 && before == after && freshNew;
                    printf("bench: netids: %d bodies (%d with their creation-order id, same after the load: %s); %d animals, unique ids %s, "
                           "kept over save/load %s; a new id after the load is fresh: %s -> %s\n",
                           bodies, bodyOk, bodyIds == bodyIdsAfter ? "yes" : "no", Animals_Count(), unique ? "yes" : "no",
                           before == after ? "yes" : "no", freshNew ? "yes" : "no", pass ? "PASS" : "FAIL");
                } else if (benchFall && (fallPhase > 0 || WorldIdle())) {
                    // --test-fall (#263): drops onto flat ground near the start. Expected on Terra (g 5.9, safe 6 m/s, x8.5):
                    // 3 m none, 5 m none (14 would be under the 15-point grace, #536), 10 m ~41, 20 m ~79, 30 m and 100 m fatal.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const float heights[] = { 3.0f, 5.0f, 10.0f, 20.0f, 30.0f, 100.0f };
                    const float lo[] = { 0.0f, 0.0f, 25.0f, 60.0f, 100.0f, 100.0f }, hi[] = { 0.0f, 0.0f, 55.0f, 95.0f, 1e9f, 1e9f };
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
                        if (++fallIdx < 6) { fallPhase = 1; fallStart = GetTime() - 1.0; }
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
                                else {
                                    dryRun += 5.0f; worstDry = std::max(worstDry, dryRun);
                                    if (dryRun == 500.0f) { // where a long dry run is (diagnosis)
                                        float lv = 0.0f, lo = 0.0f;
                                        const bool lk = net->lakeAt(dirs[i], lv, lo);
                                        printf("bench: rivers: dry 500 m at dir %.6f,%.6f,%.6f: path level %.2f, water %.2f, ground %.2f, lake %d (level %.2f, outside %.0f m)\n",
                                               dirs[i].x, dirs[i].y, dirs[i].z, levels[i], w == PlanetField::kNoWater ? -999.0f : w, e, lk ? 1 : 0, lv, lo);
                                    }
                                }
                            }
                            const float wetShare = dirs.empty() ? 0.0f : (float)wet / dirs.size();
                            check(dirs.size() > 100, "a long river");
                            check(worstUp < 0.05f, "the surface never rises walking downstream");
                            check(wetShare > 0.9f, "water in the channel along the river");
                            check(!dirs.empty() && pf->elevation(dirs.back()) < 3.0f, "it ends at the sea (or the polar ice shelf)");
                            riverLongest = s.longestKm; riverBuild = (float)buildS; riverWet = wetShare;
                            printf("bench: rivers: network built in %.2f s: %d cells of %.0f m, %d land, %d river (%.1f%%), %d lake cells, %d mouths; %d links: worst rise %.5f m, worst narrowing %.5f m\n",
                                   buildS, s.cells, s.cellSize, s.landCells, s.riverCells, 100.0f * s.riverCells / std::max(1, s.landCells), s.lakeCells, s.mouths, links, rise, shrink);
                            printf("bench: rivers: %d river cells above 700 m, highest source at %.0f m (#359)\n", s.riverCellsAbove700, s.highestSource);
                            if (net->bakeNative()) {
                                // Heads (#625): a river starts as a seep and a trickle, never at a waterfall (unless the ground drops there).
                                printf("bench: rivers: heads (#625): %d sources, %d with a waterfall in their first 3 links\n", s.heads, s.headFalls);
                                check(s.headFalls == 0, "no river starts at a waterfall (#625)");
                                Vector3 hd, dn;
                                if (net->nearestHead(pf->findSpawnDir(p->seed), 2.0f, hd, dn))
                                    printf("bench: rivers: a river head %.1f km from the start; look downstream from behind it: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f --view-fall-dist 40 (up it from below: the same with %.5f %.5f %.5f --view-fall-dist 30)\n",
                                           Vector3Distance(hd, pf->findSpawnDir(p->seed)) * pf->radius / 1000.0f, hd.x, hd.y, hd.z, -dn.x, -dn.y, -dn.z, dn.x, dn.y, dn.z);
                            }
                            printf("bench: rivers: longest %.1f km, walked every 5 m: %.0f%% wet, longest dry run %.0f m, worst rise %.3f m, from %.1f m to %.1f m\n",
                                   s.longestKm, 100.0f * wetShare, worstDry, worstUp, levels.empty() ? 0.0f : levels.front(), levels.empty() ? 0.0f : levels.back());
                            {   // Walls of water (#540): along every river, water standing more than 3 m deeper than its channel - over
                                // ground that fell away under a level drawn as a straight line between the nodes (a cliff).
                                long n = 0, walls = 0, walls10 = 0, walls50 = 0, lakeWalls = 0; float worst = 0.0f; Vector3 worstDir{};
                                net->forEachCentrePoint(10.0f, [&](Vector3 d, float level, float depth) {
                                    const float excess = level - pf->elevation(d) - depth;
                                    ++n;
                                    if (excess > 3.0f) {
                                        ++walls; walls10 += excess > 10.0f; walls50 += excess > 50.0f;
                                        float lv, lo;
                                        if (net->lakeAt(d, lv, lo) && lo < 2.0f * RiverNetwork::kLakeMargin) ++lakeWalls; // a lake's own water, not a river's
                                    }
                                    if (excess > worst) { worst = excess; worstDir = d; }
                                });
                                printf("bench: rivers: walls of water (#540): %ld of %ld centreline points every 10 m (%.2f%%) stand > 3 m deeper than the channel "
                                       "(> 10 m: %ld, > 50 m: %ld; in or by a lake: %ld), worst %.1f m at dir %.5f,%.5f,%.5f\n",
                                       walls, n, 100.0 * walls / std::max(1L, n), walls10, walls50, lakeWalls, worst, worstDir.x, worstDir.y, worstDir.z);
                                {   // Plan 10 baseline (#638, #637, #631): river sizes, length-weighted, outside lakes; at each link's middle
                                    // whether the generated water is there and deep enough to draw (> 0.15 m over the ground), and whether
                                    // the ground 2 m beyond either edge lies under the level (water that would stand as a wall or spill).
                                    constexpr int kWB = 8;
                                    const float wEdge[kWB] = { 2, 5, 10, 20, 30, 50, 90, 1e9f }, dEdge[kWB] = { 0.5f, 1, 2, 3, 4, 5, 6, 1e9f };
                                    double wLen[kWB] = {}, dLen[kWB] = {}, total = 0.0, visible = 0.0, spill = 0.0, spillDeep[5] = {}, spillWide[5] = {};
                                    float widest = 0.0f;
                                    net->forEachLink([&](const RiverNetwork::Link& L) {
                                        if (L.inLake) return;
                                        const float len = Vector3Distance(L.a, L.b), w = 0.5f * (L.wa + L.wb), dp = 0.5f * (L.da + L.db);
                                        int i = 0; while (w >= wEdge[i]) ++i;
                                        int j = 0; while (dp >= dEdge[j]) ++j;
                                        wLen[i] += len; dLen[j] += len; total += len; widest = std::max(widest, std::max(L.wa, L.wb));
                                        const Vector3 mid = Vector3Scale(Vector3Add(L.a, L.b), 0.5f), up = Vector3Normalize(mid);
                                        const Vector3 side = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(L.b, L.a), up));
                                        const float lvl = 0.5f * (L.la + L.lb), wl = pf->waterElevation(up);
                                        if (wl != PlanetField::kNoWater && wl - pf->elevation(up) > 0.15f) visible += len;
                                        for (int s = -1; s <= 1; s += 2) {
                                            const Vector3 q = Vector3Normalize(Vector3Add(mid, Vector3Scale(side, s * (0.5f * w + 2.0f))));
                                            const float deficit = lvl - pf->elevation(q);
                                            if (deficit > 0.3f) {
                                                spill += 0.5 * len;
                                                // How far under the water, and how far out the ground stays under it (#640: lower the level, or
                                                // fill the low ground beside the river?).
                                                static const float kDeep[4] = { 1.0f, 2.0f, 5.0f, 10.0f }, kWide[4] = { 6.0f, 14.0f, 30.0f, 60.0f };
                                                int b = 0; while (b < 4 && deficit >= kDeep[b]) ++b;
                                                spillDeep[b] += 0.5 * len;
                                                float out = 2.0f;
                                                while (out < 60.0f && pf->elevation(Vector3Normalize(Vector3Add(mid, Vector3Scale(side, s * (0.5f * w + out))))) < lvl) out += 4.0f;
                                                int c = 0; while (c < 4 && out >= kWide[c]) ++c;
                                                spillWide[c] += 0.5 * len;
                                            }
                                        }
                                    });
                                    printf("bench: rivers: ground under the level 2 m past an edge (km): by depth under it <1 m %.0f, 1-2 m %.0f, 2-5 m %.0f, 5-10 m %.0f, >10 m %.0f; "
                                           "low out to <6 m %.0f, 6-14 m %.0f, 14-30 m %.0f, 30-60 m %.0f, >60 m %.0f (#640)\n",
                                           spillDeep[0] / 1000.0, spillDeep[1] / 1000.0, spillDeep[2] / 1000.0, spillDeep[3] / 1000.0, spillDeep[4] / 1000.0,
                                           spillWide[0] / 1000.0, spillWide[1] / 1000.0, spillWide[2] / 1000.0, spillWide[3] / 1000.0, spillWide[4] / 1000.0);
                                    {   // Cameras for captures (#638): the nearest link to the start in each width class, seen from 30 m
                                        // downstream and 9 m up, looking back up the river.
                                        const Vector3 home = pf->findSpawnDir(p->seed);
                                        const float lo[3] = { 2.0f, 10.0f, 20.0f }, hi[3] = { 5.0f, 20.0f, 50.0f };
                                        float best[3] = { 1e30f, 1e30f, 1e30f }; Vector3 mid[3]{}, down[3]{};
                                        net->forEachLink([&](const RiverNetwork::Link& L) {
                                            if (L.inLake) return;
                                            const float w = 0.5f * (L.wa + L.wb);
                                            const Vector3 m = Vector3Normalize(Vector3Add(L.a, L.b));
                                            const float d = Vector3Distance(m, home);
                                            for (int c = 0; c < 3; ++c)
                                                if (w >= lo[c] && w < hi[c] && d < best[c]) { best[c] = d; mid[c] = m; down[c] = Vector3Normalize(Vector3Subtract(L.b, L.a)); }
                                        });
                                        for (int c = 0; c < 3; ++c)
                                            if (best[c] < 1e29f)
                                                printf("bench: rivers: a %g-%g m river %.1f km from the start: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f --view-fall-dist 30\n",
                                                       lo[c], hi[c], best[c] * pf->radius / 1000.0f, mid[c].x, mid[c].y, mid[c].z, down[c].x, down[c].y, down[c].z);
                                    }
                                    auto row = [&](const char* name, const float* edge, const double* lenB) {
                                        printf("bench: rivers: %s (km of river, outside lakes):", name);
                                        for (int i = 0; i < kWB; ++i) printf(" <%g: %.0f", edge[i] > 1e8f ? 999.0f : edge[i], lenB[i] / 1000.0);
                                        printf("\n");
                                    };
                                    row("width m", wEdge, wLen);
                                    row("depth m", dEdge, dLen);
                                    printf("bench: rivers: %.0f km of river, widest %.1f m; drawn water at link middles %.1f%%; ground under the level 2 m past an edge %.1f%% (#631, #610)\n",
                                           total / 1000.0, widest, 100.0 * visible / std::max(1.0, total), 100.0 * spill / std::max(1.0, total));
                                    // Junctions (#632): where two or more links end at one node, the generated water's surface over a
                                    // disc of the outflow's width around it: the spread of levels (a step in the surface) and the water
                                    // standing over the ground there more than the channel depth (raised water).
                                    std::unordered_map<int64_t, std::pair<int, RiverNetwork::Link>> ends;
                                    auto keyOf = [](Vector3 v) { return ((int64_t)lroundf(v.x) * 73856093LL) ^ ((int64_t)lroundf(v.y) * 19349663LL) ^ ((int64_t)lroundf(v.z) * 83492791LL); };
                                    net->forEachLink([&](const RiverNetwork::Link& L) { if (!L.inLake) { auto& e = ends[keyOf(L.b)]; ++e.first; if (e.first == 1 || L.wb > e.second.wb) e.second = L; } });
                                    int junctions = 0, stepped = 0, raised = 0; float worstStep = 0.0f, worstRaise = 0.0f; Vector3 worstJ{};
                                    for (const auto& kv : ends) {
                                        if (kv.second.first < 2 || junctions >= 4000) continue;
                                        const RiverNetwork::Link& L = kv.second.second;
                                        ++junctions;
                                        const Vector3 up = Vector3Normalize(L.b);
                                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
                                        float lo = 1e30f, hi = -1e30f, raise = 0.0f;
                                        const float r = std::max(2.0f, 0.75f * L.wb);
                                        for (float e = -r; e <= r; e += r / 3.0f)
                                            for (float n = -r; n <= r; n += r / 3.0f) {
                                                const Vector3 q = Vector3Normalize(Vector3Add(L.b, Vector3Add(Vector3Scale(t1, e), Vector3Scale(t2, n))));
                                                const float w = pf->waterElevation(q);
                                                if (w == PlanetField::kNoWater) continue;
                                                lo = std::min(lo, w); hi = std::max(hi, w);
                                                raise = std::max(raise, w - pf->elevation(q) - std::max(L.db, 0.5f));
                                            }
                                        if (hi < lo) continue;
                                        if (hi - lo > 0.3f) ++stepped;
                                        if (raise > 0.5f) ++raised;
                                        if (hi - lo > worstStep) { worstStep = hi - lo; worstJ = up; }
                                        worstRaise = std::max(worstRaise, raise);
                                    }
                                    printf("bench: rivers: junctions (#632): %d checked, %d with a surface step > 0.3 m (worst %.2f m at dir %.5f,%.5f,%.5f), %d with water > 0.5 m over its channel depth (worst %.2f m)\n",
                                           junctions, stepped, worstStep, worstJ.x, worstJ.y, worstJ.z, raised, worstRaise);
                                }
                                // The waterfalls, and the nearest big ones to the start (their --view-fall, to look).
                                const std::vector<RiverNetwork::Fall> falls = net->waterfalls();
                                int f10 = 0, f30 = 0; float tallest = 0.0f;
                                for (const auto& f : falls) { f10 += f.drop > 10.0f; f30 += f.drop > 30.0f; tallest = std::max(tallest, f.drop); }
                                printf("bench: rivers: %zu waterfalls (> 10 m: %d, > 30 m: %d, tallest %.0f m)\n", falls.size(), f10, f30, tallest);
                                {   // Beds and floodplains (#595), every 500 m along every river: the bed's material by how fast the
                                    // reach runs, and the ground 25 m to each side where it is a floodplain (flat, < 2.5 m over the water).
                                    // Banks (#609): the ground beside the water at 1-25 m out from its edge, where it lies near the
                                    // water's height (-0.3..+2 m, where the old silt band was), by the river's size (stream: half width
                                    // < 1.5 m, mid, river: >= 6 m): the share taking a wet material (silt, mud, river stones, clay). A
                                    // road-like band shows as a high share far out; a trickle should wet almost nothing.
                                    int cls[3] = {}, own[3] = {};
                                    constexpr int kD = 5;
                                    const float kDist[kD] = { 1.0f, 3.0f, 6.0f, 12.0f, 25.0f };
                                    int bankN[3][kD] = {}, bankWet[3][kD] = {}, bankSilt[3][kD] = {};
                                    const Vector3 home609 = pf->findSpawnDir(p->seed);
                                    float bestBank[2] = { 1e30f, 1e30f }; Vector3 bankDir[2]{}, bankDown[2]{};
                                    net->forEachCentrePoint(500.0f, [&](Vector3 d, float level, float) {
                                        const RiverNetwork::Hit h = net->nearest(d, 5.0f);
                                        const float e = pf->elevation(d);
                                        if (!h.found || level < 1.0f) return; // inland
                                        float lv, lo;
                                        if (net->lakeAt(d, lv, lo) && lo < 2.0f * RiverNetwork::kLakeMargin) return; // a lake's bed or shore
                                        if (e <= level - 0.5f) { // a bed under water
                                            const Vector3 bed = Vector3Scale(d, pf->radius + e - 0.2f);
                                            const VoxelMaterial m = pf->material(bed, pf->baseDensity(bed));
                                            if (m != VoxelMaterial::Ice && m != VoxelMaterial::Snow && m != VoxelMaterial::Rock) { // not frozen, or a rock outcrop
                                                const int k = h.slope > 0.015f ? 0 : h.slope > 0.003f ? 1 : 2;
                                                ++cls[k];
                                                own[k] += k == 0 ? (m == VoxelMaterial::Pebbles || m == VoxelMaterial::Gravel)
                                                        : k == 1 ? (m == VoxelMaterial::Sand || m == VoxelMaterial::Gravel)
                                                                 : (m == VoxelMaterial::Silt || m == VoxelMaterial::Sand || m == VoxelMaterial::Clay || m == VoxelMaterial::Gravel);
                                            }
                                        }
                                        const int sz = h.halfWidth < 1.5f ? 0 : h.halfWidth < 6.0f ? 1 : 2;
                                        if (sz > 0 && h.slope < 0.015f && level > 25.0f) { // the nearest calm mid-size and big river inland of the start (to look at its banks)
                                            const float ds = Vector3Distance(d, home609);
                                            if (ds < bestBank[sz - 1]) { bestBank[sz - 1] = ds; bankDir[sz - 1] = d; bankDown[sz - 1] = h.down; }
                                        }
                                        const Vector3 side = Vector3Normalize(Vector3CrossProduct(d, h.down));
                                        for (int sgn = -1; sgn <= 1; sgn += 2)
                                            for (int i = 0; i < kD; ++i) {
                                                const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(side, sgn * (h.halfWidth + kDist[i]) / pf->radius)));
                                                const float eq = pf->elevation(q);
                                                if (eq < h.level - 0.3f || eq > h.level + 2.0f || pf->waterElevation(q) != PlanetField::kNoWater) continue;
                                                const Vector3 g = Vector3Scale(q, pf->radius + eq - 0.2f);
                                                const VoxelMaterial m = pf->material(g, pf->baseDensity(g));
                                                ++bankN[sz][i];
                                                bankWet[sz][i] += m == VoxelMaterial::Silt || m == VoxelMaterial::Mud || m == VoxelMaterial::Pebbles || m == VoxelMaterial::Clay;
                                                bankSilt[sz][i] += m == VoxelMaterial::Silt;
                                            }
                                    });
                                    auto share = [](int a, int b) { return b ? 100.0f * a / b : 0.0f; };
                                    printf("bench: rivers: beds (#595): fast %d (%.0f%% river stones or gravel), moderate %d (%.0f%% sand or gravel), slow %d (%.0f%% silt, sand, clay or gravel)\n",
                                           cls[0], share(own[0], cls[0]), cls[1], share(own[1], cls[1]), cls[2], share(own[2], cls[2]));
                                    const char* szName[3] = { "streams (half width < 1.5 m)", "mid (1.5-6 m)", "rivers (>= 6 m)" };
                                    for (int s2 = 0; s2 < 3; ++s2) {
                                        printf("bench: rivers: banks (#609), %s: wet ground (silt) at", szName[s2]);
                                        for (int i = 0; i < kD; ++i) printf(" %.0f m %.0f%% (%.0f%%) n %d%s", kDist[i], share(bankWet[s2][i], bankN[s2][i]), share(bankSilt[s2][i], bankN[s2][i]), bankN[s2][i], i + 1 < kD ? "," : "\n");
                                    }
                                    for (int s2 = 0; s2 < 2; ++s2)
                                        if (bestBank[s2] < 1e29f)
                                            printf("bench: rivers: nearest calm %s river %.1f km from the start; look upstream along it: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f --view-fall-dist 30\n",
                                                   s2 ? "big" : "mid-size", bestBank[s2] * pf->radius / 1000.0f, bankDir[s2].x, bankDir[s2].y, bankDir[s2].z, bankDown[s2].x, bankDown[s2].y, bankDown[s2].z);
                                    check(share(own[0], cls[0]) > 60.0f && share(own[1], cls[1]) > 60.0f && share(own[2], cls[2]) > 60.0f, "beds by the river's speed (#595)");
                                    // A trickle wets almost nothing beyond its edge, and no river lays a road: little wet ground 25 m out.
                                    check(share(bankWet[0][1] + bankWet[0][2] + bankWet[0][3] + bankWet[0][4], bankN[0][1] + bankN[0][2] + bankN[0][3] + bankN[0][4]) < 3.0f,
                                          "streams wet almost nothing 3 m+ from their edge (#609)");
                                    check(share(bankWet[2][4], bankN[2][4]) < 20.0f && share(bankWet[1][3] + bankWet[1][4], bankN[1][3] + bankN[1][4]) < 10.0f,
                                          "no band of wet ground far from a river (#609)");
                                    check(share(bankSilt[2][1] + bankSilt[2][2], bankN[2][1] + bankN[2][2]) > 5.0f, "silt in patches on big rivers' floodplains (#595)");
                                }
                                {   // Perched rivers (#547; the measure from the water session's m102): a river running beside lower
                                    // ground, as along a canyon's edge instead of in it. Every 100 m along every river (outside lakes):
                                    // the ground 30 m out from either bank more than 5 m under the river's bed. Reported.
                                    long n = 0, perched = 0;
                                    net->forEachCentrePoint(100.0f, [&](Vector3 d, float level, float) {
                                        const RiverNetwork::Hit h = net->nearest(d, 5.0f);
                                        float lv, lo;
                                        if (!h.found || level < 1.0f || (net->lakeAt(d, lv, lo) && lo < RiverNetwork::kLakeMargin)) return;
                                        const float bed = pf->elevation(d);
                                        const Vector3 side = Vector3Normalize(Vector3CrossProduct(d, h.down));
                                        bool low = false;
                                        for (int sgn = -1; sgn <= 1 && !low; sgn += 2)
                                            low = pf->elevation(Vector3Normalize(Vector3Add(d, Vector3Scale(side, sgn * (h.halfWidth + 30.0f) / pf->radius)))) < bed - 5.0f;
                                        ++n;
                                        perched += low;
                                    });
                                    printf("bench: rivers: perched (#547): %ld of %ld points every 100 m have ground 5 m+ under the bed 30 m out from a bank (%.2f%%)\n",
                                           perched, n, 100.0 * perched / std::max(1L, n));
                                }
                                {   // Slots (#547: "water carves large chunks"): river points (every 100 m, outside lakes) in a trench - on both
                                    // sides the ground steps 8 m+ up within 20 m of the bank and then stays near flat to 60 m (a natural V valley
                                    // keeps rising).
                                    long n = 0, slots = 0;
                                    net->forEachCentrePoint(100.0f, [&](Vector3 d, float level, float) {
                                        const RiverNetwork::Hit h = net->nearest(d, 5.0f);
                                        float lv, lo;
                                        if (!h.found || level < 1.0f || (net->lakeAt(d, lv, lo) && lo < RiverNetwork::kLakeMargin)) return;
                                        const Vector3 side = Vector3Normalize(Vector3CrossProduct(d, h.down));
                                        bool both = true;
                                        for (int sgn = -1; sgn <= 1 && both; sgn += 2) {
                                            auto at = [&](float off) { return pf->elevation(Vector3Normalize(Vector3Add(d, Vector3Scale(side, sgn * (h.halfWidth + off) / pf->radius)))); };
                                            const float e20 = at(20.0f), e60 = at(60.0f);
                                            both = e20 > level + 8.0f && e60 - e20 < 0.25f * (e20 - level);
                                        }
                                        ++n; slots += both;
                                    });
                                    printf("bench: rivers: slots (#547): %ld of %ld points every 100 m sit in a trench (8 m+ up within 20 m of both banks, near flat to 60 m) (%.2f%%)\n",
                                           slots, n, 100.0 * slots / std::max(1L, n));
                                }
                                if (gErosionPrototype) {
                                    // Canyons carry rivers (#544, #547), over the whole planet: points ~100 m apart on the canyon lines
                                    // (notch >= 0.8, on land), how many have a river within 150 m - over the whole line, and past
                                    // its first 3 km from the head (a stream needs ~2 km^2 of drainage to start). And rim rivers:
                                    // river points (every 50 m) with ground 15 m+ under their bed 60-250 m to the side, in a canyon.
                                    const int M = 1 << 22;
                                    long cn = 0, cwet = 0;
                                    float canyonLookD = 1e30f; Vector3 canyonLook{}, canyonLookDown{};
                                    for (int i = 0; i < M; ++i) {
                                        const float y = 1.0f - 2.0f * (i + 0.5f) / M, rr = sqrtf(std::max(0.0f, 1.0f - y * y)), ph = i * 2.39996323f;
                                        const Vector3 d = { rr * cosf(ph), y, rr * sinf(ph) };
                                        const Vector3 s = Vector3Scale(d, pf->radius);
                                        if (pf->canyonFault(s.x, s.y, s.z) < 0.8f) continue;
                                        if (pf->elevationBase(d, 50.0f) < 5.0f) continue;
                                        ++cn;
                                        const RiverNetwork::Hit ch = net->nearest(d, 150.0f);
                                        cwet += ch.found;
                                        // The nearest to the start, for a look (a river in a canyon, from downstream).
                                        const float dist = Vector3Distance(d, pf->findSpawnDir(p->seed));
                                        if (ch.found && ch.level > 2.0f && dist < canyonLookD) { canyonLookD = dist; canyonLook = d; canyonLookDown = ch.down; }
                                    }
                                    if (canyonLookD < 1e29f)
                                        printf("bench: canyons: a river in a canyon %.1f km from the start; look: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f --view-fall-dist 250\n",
                                               canyonLookD * pf->radius / 1000.0f, canyonLook.x, canyonLook.y, canyonLook.z, canyonLookDown.x, canyonLookDown.y, canyonLookDown.z);
                                    long rn = 0, rim = 0;
                                    net->forEachCentrePoint(50.0f, [&](Vector3 d, float level, float) {
                                        const Vector3 s = Vector3Scale(d, pf->radius);
                                        if (level < 1.0f || pf->canyonLine(s.x, s.y, s.z) > 0.06f) return; // only near the canyon lines
                                        const RiverNetwork::Hit h = net->nearest(d, 5.0f);
                                        if (!h.found) return;
                                        float lv, lo;
                                        if (net->lakeAt(d, lv, lo) && lo < RiverNetwork::kLakeMargin) return;
                                        const float bed = pf->elevation(d);
                                        const Vector3 side = Vector3Normalize(Vector3CrossProduct(d, h.down));
                                        bool low = false;
                                        for (int sgn = -1; sgn <= 1 && !low; sgn += 2)
                                            for (float off = 60.0f; off <= 250.0f && !low; off += 38.0f) {
                                                const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(side, sgn * (h.halfWidth + off) / pf->radius)));
                                                const Vector3 qs = Vector3Scale(q, pf->radius);
                                                low = pf->canyonFault(qs.x, qs.y, qs.z) > 0.5f && pf->elevationBase(q, 8.0f) < bed - 15.0f;
                                            }
                                        ++rn; rim += low;
                                    });
                                    const float wetPct = 100.0f * cwet / std::max(1L, cn), rimPct = 100.0f * rim / std::max(1L, rn);
                                    printf("bench: canyons: %ld points on canyon lines (planet), %.1f%% with a river within 150 m; %ld river points near canyon lines, %ld (%.2f%%) on a rim above a lower canyon floor\n",
                                           cn, wetPct, rn, rim, rimPct);
                                    static const bool graded = !(getenv("RAYLIBGAME_CANYON_GRADE") && getenv("RAYLIBGAME_CANYON_GRADE")[0] == '0');
                                    if (graded) {
                                        check(wetPct > 20.0f, "rivers run in the canyons (#544; 18.6% before the grading on seed 42, 25% after; bake-native rivers raise it)");
                                        check(rimPct < 5.0f, "no rivers on a canyon's rim above its floor (#547)");
                                    }
                                }
                                const Vector3 home = pf->findSpawnDir(p->seed);
                                std::vector<std::pair<float, int>> near;
                                for (int i = 0; i < (int)falls.size(); ++i)
                                    if (falls[i].drop >= 8.0f && falls[i].top > 2.0f) near.push_back({ Vector3Distance(home, falls[i].lip) * pf->radius, i });
                                std::sort(near.begin(), near.end());
                                for (size_t k = 0; k < std::min<size_t>(3, near.size()); ++k) {
                                    const auto& f = falls[near[k].second];
                                    printf("bench: rivers: waterfall %.1f km from the start: drops %.1f m from %.1f m, %.0f m wide; look: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f\n",
                                           near[k].first / 1000.0f, f.drop, f.top, 2.0f * f.halfWidth, f.lip.x, f.lip.y, f.lip.z, f.down.x, f.down.y, f.down.z);
                                }
                                // And a small one near the start (#614: its effects should be small too): the least drop x width within 6 km.
                                int small = -1; float smallD = 1e30f;
                                for (int i = 0; i < (int)falls.size(); ++i)
                                    if (falls[i].top > 2.0f && Vector3Distance(home, falls[i].lip) * pf->radius < 6000.0f) { const float sz = falls[i].drop * falls[i].halfWidth; if (sz < smallD) { smallD = sz; small = i; } }
                                if (small >= 0) {
                                    const auto& f = falls[small];
                                    printf("bench: rivers: small waterfall %.1f km from the start: drops %.1f m from %.1f m, %.0f m wide; look: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f --view-fall-dist 25\n",
                                           Vector3Distance(home, f.lip) * pf->radius / 1000.0f, f.drop, f.top, 2.0f * f.halfWidth, f.lip.x, f.lip.y, f.lip.z, f.down.x, f.down.y, f.down.z);
                                }
                            }
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
                            // Canyons (#544): the bake's canyon lines within 30 km of the start - how deep the baked ground
                            // lies under the land 400 m to either side, and how many carry a river. Information only.
                            if (gErosionPrototype) {
                                const Vector3 home = pf->findSpawnDir(p->seed);
                                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                                auto at = [&](Vector3 d, float e, float n) { return Vector3Normalize(Vector3Add(d, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius)))); };
                                int n = 0, deep = 0, wet = 0; double sum = 0.0; float bestD = 1e30f; Vector3 bestDir{}, bestSide{};
                                for (float e = -30000.0f; e <= 30000.0f; e += 200.0f)
                                    for (float n2 = -30000.0f; n2 <= 30000.0f; n2 += 200.0f) {
                                        const Vector3 d = at(home, e, n2);
                                        const Vector3 s = Vector3Scale(d, pf->radius);
                                        if (pf->canyonFault(s.x, s.y, s.z) < 0.8f) continue;
                                        const float h = pf->elevationBase(d, 8.0f);
                                        if (h < 5.0f) continue;
                                        float rim = 0.0f; Vector3 side{};
                                        for (int k = 0; k < 8; ++k) {
                                            const float a = k * 0.785398f;
                                            const Vector3 q = at(d, 400.0f * cosf(a), 400.0f * sinf(a));
                                            const float hq = pf->elevationBase(q, 8.0f);
                                            if (hq - h > rim) { rim = hq - h; side = Vector3Subtract(q, d); }
                                        }
                                        ++n; sum += rim; deep += rim > 25.0f;
                                        const RiverNetwork::Hit hit = net->nearest(d, 150.0f);
                                        wet += hit.found;
                                        if (rim > 40.0f && e * e + n2 * n2 < bestD) { bestD = e * e + n2 * n2; bestDir = d; bestSide = Vector3Normalize(side); }
                                    }
                                printf("bench: canyons: %d samples on canyon lines within 30 km, mean depth under the higher side %.1f m, %.0f%% over 25 m, %.0f%% with a river within 150 m\n",
                                       n, n ? sum / n : 0.0, 100.0 * deep / std::max(1, n), 100.0 * wet / std::max(1, n));
                                if (bestD < 1e29f) printf("bench: canyons: nearest deep one %.1f km from the start; look: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f\n",
                                                          sqrtf(bestD) / 1000.0f, bestDir.x, bestDir.y, bestDir.z, -bestSide.x, -bestSide.y, -bestSide.z);
                                // Fault scarps (#590): over the whole planet (~1.5 km apart), places where the designed uplift
                                // steps by > 40 m within 800 m, and the relief the bake left across the same span.
                                const int N = 1 << 20;
                                int nf = 0; double designed = 0.0, baked = 0.0; float fD = 1e30f; Vector3 fDir{}, fSide{};
                                for (int i = 0; i < N; ++i) {
                                    const float y = 1.0f - 2.0f * (i + 0.5f) / N, rr = sqrtf(std::max(0.0f, 1.0f - y * y)), ph = i * 2.39996323f;
                                    const Vector3 d = { rr * cosf(ph), y, rr * sinf(ph) };
                                    const Vector3 s = Vector3Scale(d, pf->radius);
                                    if (pf->faultUplift(s) <= 0.0f) continue;
                                    const Vector3 u1 = Vector3Normalize(Vector3CrossProduct(d, { 0, 1, 0 })), u2 = Vector3CrossProduct(u1, d);
                                    float lo = 1e9f, hi = -1e9f; Vector3 dlo{}, dhi{};
                                    for (int k = 0; k < 4; ++k) {
                                        const float a = k * 1.5707963f;
                                        const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(Vector3Add(Vector3Scale(u1, cosf(a)), Vector3Scale(u2, sinf(a))), 800.0f / pf->radius)));
                                        const float up = pf->faultUplift(Vector3Scale(q, pf->radius));
                                        if (up < lo) { lo = up; dlo = q; }
                                        if (up > hi) { hi = up; dhi = q; }
                                    }
                                    if (hi - lo < 40.0f) continue;
                                    const float hl = pf->elevationBase(dlo, 50.0f), hh = pf->elevationBase(dhi, 50.0f);
                                    if (hl < 5.0f || hh < 5.0f) continue;

                                    ++nf; designed += hi - lo; baked += hh - hl;
                                    const float dist = Vector3Distance(d, pf->findSpawnDir(p->seed)) * pf->radius;
                                    if (dist < fD) { fD = dist; fDir = dhi; fSide = Vector3Normalize(Vector3Subtract(dlo, dhi)); } // look from the top, down the scarp
                                }
                                printf("bench: faults: %d land samples on fault scarps, designed step %.0f m over 1.6 km, baked relief across it %.0f m\n",
                                       nf, nf ? designed / nf : 0.0, nf ? baked / nf : 0.0);
                                if (fD < 1e29f) printf("bench: faults: nearest %.1f km from the start; look: --view-fall %.6f %.6f %.6f %.5f %.5f %.5f\n",
                                                       fD / 1000.0f, fDir.x, fDir.y, fDir.z, fSide.x, fSide.y, fSide.z);
                                // Folds (#591): local relief (highest - lowest of 5 points 600 m apart) on land inside fold
                                // regions against land outside them, below 600 m (no mountain ranges).
                                {
                                    const int M = 1 << 18;
                                    int nIn = 0, nOut = 0; double rIn = 0.0, rOut = 0.0, iIn = 0.0;
                                    for (int i = 0; i < M; ++i) {
                                        const float y = 1.0f - 2.0f * (i + 0.5f) / M, rr = sqrtf(std::max(0.0f, 1.0f - y * y)), ph = i * 2.39996323f;
                                        const Vector3 d = { rr * cosf(ph), y, rr * sinf(ph) };
                                        const float h0 = pf->elevationBase(d, 50.0f);
                                        if (h0 < 5.0f || h0 > 600.0f) continue;
                                        const bool foldsOn = gStrataFolds; gStrataFolds = true; // the regions even in the A/B (bake done)
                                        const float w = pf->foldWeight(Vector3Scale(d, pf->radius));
                                        gStrataFolds = foldsOn;
                                        if (w > 0.0f && w < 0.9f) continue;
                                        const Vector3 u1 = Vector3Normalize(Vector3CrossProduct(d, { 0, 1, 0 })), u2 = Vector3CrossProduct(u1, d);
                                        // Ridge-scale roughness: |centre - mean of 4 points 150 m away| (base relief cancels).
                                        auto rough = [&](bool input) {
                                            PlanetField::tIgnoreBake = PlanetField::tBakeInput = input;
                                            float sum = 0.0f;
                                            for (int k = 0; k < 4; ++k)
                                                sum += pf->elevationBase(Vector3Normalize(Vector3Add(d, Vector3Scale(k & 1 ? u2 : u1, (k < 2 ? 150.0f : -150.0f) / pf->radius))), 50.0f);
                                            const float r = fabsf(pf->elevationBase(d, 50.0f) - 0.25f * sum);
                                            PlanetField::tIgnoreBake = PlanetField::tBakeInput = false;
                                            return r;
                                        };
                                        if (w >= 0.9f) { ++nIn; rIn += rough(false); iIn += rough(true); }
                                        else { ++nOut; rOut += rough(false); }
                                    }
                                    printf("bench: folds: ridge roughness (150 m) on low land: %.2f m in fold regions (%d samples; bake input %.2f m), %.2f m elsewhere (%d)\n",
                                           nIn ? rIn / nIn : 0.0, nIn, nIn ? iIn / nIn : 0.0, nOut ? rOut / nOut : 0.0, nOut);
                                    // Karst (#592): the same roughness inside karst regions, and how many points have a river near.
                                    int nk = 0, wetK = 0, nl = 0, wetL = 0; double rk = 0.0, ik = 0.0; float kD = 1e30f; Vector3 kDir{};
                                    for (int i = 0; i < M; ++i) {
                                        const float y = 1.0f - 2.0f * (i + 0.5f) / M, rr = sqrtf(std::max(0.0f, 1.0f - y * y)), ph = i * 2.39996323f;
                                        const Vector3 d = { rr * cosf(ph), y, rr * sinf(ph) };
                                        const float h0 = pf->elevationBase(d, 50.0f);
                                        if (h0 < 5.0f) continue;
                                        const bool wet = net->nearest(d, 150.0f).found;
                                        if (pf->karstWeight(Vector3Scale(d, pf->radius)) < 0.9f) { ++nl; wetL += wet; continue; }
                                        const Vector3 u1 = Vector3Normalize(Vector3CrossProduct(d, { 0, 1, 0 })), u2 = Vector3CrossProduct(u1, d);
                                        auto rough = [&](bool input) {
                                            PlanetField::tIgnoreBake = PlanetField::tBakeInput = input;
                                            float sum = 0.0f;
                                            for (int k = 0; k < 4; ++k)
                                                sum += pf->elevationBase(Vector3Normalize(Vector3Add(d, Vector3Scale(k & 1 ? u2 : u1, (k < 2 ? 150.0f : -150.0f) / pf->radius))), 50.0f);
                                            const float r = fabsf(pf->elevationBase(d, 50.0f) - 0.25f * sum);
                                            PlanetField::tIgnoreBake = PlanetField::tBakeInput = false;
                                            return r;
                                        };
                                        ++nk; wetK += wet; rk += rough(false); ik += rough(true);
                                        const float dist = Vector3Distance(d, pf->findSpawnDir(p->seed)) * pf->radius;
                                        if (dist < kD) { kD = dist; kDir = d; }
                                    }
                                    printf("bench: karst: %d land samples in karst regions: roughness (150 m) %.2f m (bake input %.2f m); a river within 150 m: %.1f%% (elsewhere %.1f%%)\n",
                                           nk, nk ? rk / nk : 0.0, nk ? ik / nk : 0.0, 100.0 * wetK / std::max(1, nk), 100.0 * wetL / std::max(1, nl));
                                    if (kD < 1e29f) printf("bench: karst: nearest %.1f km from the start, dir %.6f %.6f %.6f\n", kD / 1000.0f, kDir.x, kDir.y, kDir.z);
                                    // Volcanoes (#593): the summits (local maxima of the designed relief on the same sphere of
                                    // samples, > 900 m), the baked height there against the design, and gullies on the flanks.
                                    {
                                        int nv = 0; double designed = 0.0, baked = 0.0, flankRough = 0.0; int nfl = 0, lakes = 0;
                                        float vD = 1e30f; Vector3 vDir{};
                                        std::vector<Vector3> peaks;
                                        for (int i = 0; i < M; ++i) {
                                            const float y = 1.0f - 2.0f * (i + 0.5f) / M, rr = sqrtf(std::max(0.0f, 1.0f - y * y)), ph = i * 2.39996323f;
                                            const Vector3 d = { rr * cosf(ph), y, rr * sinf(ph) };
                                            const float v = pf->volcano(Vector3Scale(d, pf->radius));
                                            if (v <= 0.0f) continue;
                                            const Vector3 u1 = Vector3Normalize(Vector3CrossProduct(d, { 0, 1, 0 })), u2 = Vector3CrossProduct(u1, d);
                                            if (v > 300.0f && v < 0.6f * 3000.0f) { // flanks: gullies = roughness at 150 m
                                                float sum = 0.0f;
                                                for (int k = 0; k < 4; ++k)
                                                    sum += pf->elevationBase(Vector3Normalize(Vector3Add(d, Vector3Scale(k & 1 ? u2 : u1, (k < 2 ? 150.0f : -150.0f) / pf->radius))), 50.0f);
                                                flankRough += fabsf(pf->elevationBase(d, 50.0f) - 0.25f * sum); ++nfl;
                                            }
                                            if (v < 900.0f) continue;
                                            bool top = true;
                                            for (const Vector3& q : peaks) if (Vector3Distance(q, d) * pf->radius < 3000.0f) { top = false; break; }
                                            if (!top) continue;
                                            // The summit: the highest designed point within 3 km (a coarse climb).
                                            Vector3 b = d; float bv = v;
                                            for (int step = 0; step < 60; ++step) {
                                                bool moved = false;
                                                for (int k = 0; k < 8; ++k) {
                                                    const float a = k * 0.785398f;
                                                    const Vector3 q = Vector3Normalize(Vector3Add(b, Vector3Scale(Vector3Add(Vector3Scale(u1, cosf(a)), Vector3Scale(u2, sinf(a))), 100.0f / pf->radius)));
                                                    const float qv = pf->volcano(Vector3Scale(q, pf->radius));
                                                    if (qv > bv) { bv = qv; b = q; moved = true; }
                                                }
                                                if (!moved) break;
                                            }
                                            bool dup = false;
                                            for (const Vector3& q : peaks) if (Vector3Distance(q, b) * pf->radius < 3000.0f) { dup = true; break; }
                                            if (dup) continue;
                                            peaks.push_back(b);
                                            ++nv; designed += bv; baked += pf->elevationBase(b, 50.0f);
                                            int lid = -1; float lake = 0.0f, lo = 0.0f;
                                            for (int k = 0; k < 9 && lid < 0; ++k) {
                                                const Vector3 q = k == 0 ? b : Vector3Normalize(Vector3Add(b, Vector3Scale(k & 1 ? u2 : u1, ((k - 1) / 2 + 1) * (k & 2 ? 150.0f : -150.0f) / pf->radius)));
                                                if (net->lakeAt(q, lake, lo, lid) && lo < RiverNetwork::kLakeMargin) break;
                                                lid = -1;
                                            }
                                            lakes += lid >= 0;
                                            const float dist = Vector3Distance(b, pf->findSpawnDir(p->seed)) * pf->radius;
                                            if (dist < vD) { vD = dist; vDir = b; }
                                        }
                                        printf("bench: volcanoes: %d summits, designed %.0f m above the base on average, baked ground there %.0f m; %d with a crater lake; flank roughness (150 m) %.2f m\n",
                                               nv, nv ? designed / nv : 0.0, nv ? baked / nv : 0.0, lakes, nfl ? flankRough / nfl : 0.0);
                                        if (vD < 1e29f) printf("bench: volcanoes: nearest %.1f km from the start, dir %.6f %.6f %.6f\n", vD / 1000.0f, vDir.x, vDir.y, vDir.z);
                                    }
                                }
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
                            Vector3 bed = Vector3Scale(waterDir, pf->radius + bedE - 0.5f);
                            const VoxelMaterial m = pf->material(bed, pf->baseDensity(bed));
                            // #206: beds are clay, with gravel bars in the deeper middle and sand on the shores; on a bar,
                            // dig the nearest clay bed within 30 m instead.
                            if (m != VoxelMaterial::Clay) {
                                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(waterDir, { 0, 1, 0 })), t2 = Vector3CrossProduct(waterDir, t1);
                                bool found = false;
                                for (float r = 2.0f; r <= 30.0f && !found; r += 2.0f)
                                    for (int a = 0; a < 24 && !found; ++a) {
                                        const float ang = a * 2.0f * PI / 24.0f;
                                        const Vector3 dir = Vector3Normalize(Vector3Add(Vector3Scale(waterDir, pf->radius), Vector3Add(Vector3Scale(t1, cosf(ang) * r), Vector3Scale(t2, sinf(ang) * r))));
                                        const float e = pf->elevation(dir);
                                        const Vector3 q = Vector3Scale(dir, pf->radius + e - 0.5f);
                                        if (pf->material(q, pf->baseDensity(q)) == VoxelMaterial::Clay) { bed = q; found = true; }
                                    }
                            }
                            const int clay0 = inventory.count("clay");
                            for (int k = 0; k < 12; ++k) TerrainTool_DigAt(*p, Vector3Subtract(bed, Vector3Scale(waterDir, 0.4f * k)));
                            const int clay = inventory.count("clay") - clay0;
                            const bool pass = under && !above && clay > 0 && (m == VoxelMaterial::Clay || m == VoxelMaterial::Sand || m == VoxelMaterial::Gravel || m == VoxelMaterial::Pebbles); // #206: clay beds (under river stones and sand too, #595)
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
                        Survival_Update(world, 60.0f, PlayerInput::None());
                        const float tempHere = survival.temperature, drained = 100.0f - survival.food;
                        const float expected = 60.0f * survival.foodDrainPerSec * (tempHere < 5.0f ? 1.0f + (5.0f - tempHere) / 15.0f : 1.0f);
                        check(fabsf(drained - expected) < expected * 0.05f + 0.01f, "drain at rest");
                        survival.food = 0.0f; survival.health = 100.0f;
                        Survival_Update(world, 10.0f, PlayerInput::None());
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
                        Survival_Update(world, 60.0f, PlayerInput::None());
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
                        Hotbar_Bind(3, "stone_pickaxe");
                        ToolWear_OnStroke("stone_pickaxe", life);
                        check(inventory.count("stone_pickaxe") == 0, "broke at its durability");
                        check(Hotbar_Slot(3).empty(), "the broken tool left the hotbar (#209)");
                        toolsLife = life;
                        // An iron hoe: look at the ground just ahead and till.
                        inventory.add("iron_hoe", 1);
                        terrainTool.heldTool = "iron_hoe";
                        const Vector3 level = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        player.forward = Vector3Normalize(Vector3Subtract(level, Vector3Scale(player.up, 0.8f)));
                        player.applyToCamera(world);
                        TerrainTool_Update(world, PlayerInput::None(), true); // aim only
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
                        Survival_Update(world, seconds, PlayerInput::None());
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
                        // #217: the basic suit has no jetpack and no power; a thruster pack needs a battery pack's power.
                        check(thrust0 == 0.0f && !Suit_Stats().powered && survival.power == 0.0f, "the basic suit: no jetpack, no suit power");
                        fitOnly("thruster_pack");  check(player.jetpackAccel == 0.0f, "a thruster pack alone: no jetpack (no suit power)");
                        inventory.add("battery_pack", 1); Suit_Fit("battery_pack"); Suit_Apply(world, 0.0f);
                        check(player.jetpackAccel > 0.0f && near(player.jetpackFuelPerSec, fuel0 / 1.6f), "thruster + battery pack: a jetpack with +60% fuel");
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
                        // #283: a chest takes anything, also items no recipe or fuel list names (they were refused).
                        check(Factory::Accepts(PartType::Chest, "cooked_meat") && Factory::Accepts(PartType::Chest, "ruby") && Factory::Accepts(PartType::Chest, "bench_never_seen_item_283"), "a chest accepts food, gems and unknown items (#283)");
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
                        Survival_Update(world, 30.0f, PlayerInput::None());
                        const float after30 = survival.health;
                        Survival_Update(world, 40.0f, PlayerInput::None());
                        damageVacuum = 100.0f - survival.health;
                        check(survival.inVacuum && after30 >= 99.9f && damageVacuum > 10.0f, "vacuum: 60 s, then it hurts (basic suit)");
                        inventory.add("explorer_suit", 1);
                        Suit_Wear("explorer_suit"); Suit_Apply(world, 0.0f);
                        survival.health = 100.0f; survival.oxygen = 100.0f; survival.vacuumTime = 0.0f;
                        Survival_Update(world, 30.0f, PlayerInput::None()); Survival_Update(world, 40.0f, PlayerInput::None());
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
                        inventory.clear(); // room for the next crafts (a full grid refuses outputs)
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
                        Weather_Force((int)WeatherKind::Clear); Weather_Update(world, 1.0f); Survival_Update(world, 0.01f, PlayerInput::None());
                        const float tClear = survival.temperature;
                        Weather_Force((int)WeatherKind::Rain); Weather_Update(world, 1.0f); Survival_Update(world, 0.01f, PlayerInput::None());
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
                        for (const char* w : { "any_planks", "any_log", "any_stone" }) if (d.id != w && ItemMatches(w, d.id)) { Rule u; u.in = { d.id }; u.out = { w }; rules.push_back(u); }
                    { Rule u; u.in = { "bucket" }; u.out = { "water_bucket" }; rules.push_back(u); } // scooping water
                    // Plan 5 biome resources: oil seeps (a bucket), desert scrub, jungle trees, swamp digging.
                    { Rule u; u.in = { "bucket" }; u.out = { "crude_oil" }; rules.push_back(u); }
                    { Rule u; u.in = {}; u.out = { "cactus", "rubber", "peat" }; rules.push_back(u); }
                    { Rule u; u.in = { "fishing_rod" }; u.out = { "raw_cod", "raw_trout", "raw_catfish" }; rules.push_back(u); } // fishing (#176)
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
                    // #214: ores only go into smelting and ore processing; every product is made from ingots (plates, wire...).
                    static const char* kOreMachines[] = { "kiln", "electric_furnace", "steam_macerator", "macerator", "arc_furnace", "blast_furnace",
                                                          "electrolyzer", "chemical_reactor", "ore_washer", "thermal_centrifuge" };
                    std::vector<std::string> oreUses;
                    for (const Rule& u : rules)
                        for (const std::string& in : u.in) {
                            if (in.size() < 4 || in.compare(in.size() - 4, 4, "_ore") != 0) continue;
                            bool ok = false;
                            for (const char* m : kOreMachines) ok = ok || u.machine == m;
                            if (!ok) oreUses.push_back((u.machine.empty() ? std::string("hand") : u.machine) + ": " + (u.out.empty() ? std::string("?") : u.out[0]) + " from " + in);
                        }
                    printf("\nbench: recipes: ore outside smelting/processing (#214): %zu", oreUses.size());
                    for (const auto& s : oreUses) printf("  [%s]", s.c_str());
                    printf("\nbench: recipes: %s\n", missing.empty() && oreUses.empty() ? "PASS" : "FAIL");
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
                        HelpManual_OpenItem("iron_ingot"); // a page with every section: made by hand and in machines, used in, ...
                        GameUI_ShowTab("Help"); // the manual (#354)
                        benchShotTabs = 30; frames = 0;
                    } else if (benchShotTabs == 30 && ++frames == 40) {
                        TakeScreenshot("bench_help.png");
                        printf("bench: help manual: %d topics (bench_help.png)\n", HelpManual_TopicCount());
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
                        static float oilAt10 = -1.0f, rateBefore = 0.0f, oilAfter5 = -1.0f; // #651: the fill rate before the save and after the load
                        if (oilPhase == 1 && t > 10.0 && oilAt10 < 0.0f) oilAt10 = oil;
                        if (oilPhase == 1 && t > 25.0) {
                            rateBefore = (oil - oilAt10) / 15.0f;
                            before = oil;
                            const float took = Fluid::Take(*p, pool, 1.0f, 2.5f, Fluid::KindByName("crude_oil"));
            { const Fluid::Stats fs = Fluid::GetStats(); printf("bench: oil: %.1f m^3 in the pool after 25 s (water near: %.1f m^3); a bucket takes %.1f; sources %d, sites %d, bricks %d, other fluid %.1f m^3, kinds %d\n", oil, Fluid::WaterAt(*p, pool, 2.0f), took, Fluid::SourceCount(), fs.sites, fs.bricks, fs.otherVolume, Fluid::KindCount()); }
                            oilPhase = 2; oilT0 = GetTime();
                        } else if (oilPhase == 2 && t > 3.0) {
                            const bool saved = SaveWorld(world), loaded = saved && LoadWorldSave(world);
                            printf("bench: oil: after the bucket %.1f m^3; saved %d, loaded %d\n", oil, (int)saved, (int)loaded);
                            oilPhase = 4; oilT0 = GetTime();
                        } else if (oilPhase == 4 && t > 5.0 && oilAfter5 < 0.0f) {
                            oilAfter5 = oil;
                        } else if (oilPhase == 4 && t > 20.0) {
                            // The same seep fills at the same rate after the load (#651: user saw it speed up).
                            const float rateAfter = (oil - oilAfter5) / 15.0f;
                            oilRateOk = Fluid::SourceCount() == 1 && rateAfter <= 1.25f * rateBefore + 0.003f;
                            printf("bench: oil: fill rate %.4f m^3/s before the save (10-25 s), %.4f after the load (5-20 s); %d source(s) (#651) -> %s\n",
                                   rateBefore, rateAfter, Fluid::SourceCount(), oilRateOk ? "ok" : "FASTER");
                            oilPhase = 3; oilT0 = GetTime(); Fluid::SetSpeed(8); // fast-forward: the seep fills its basin to its head and rests
                        } else if (oilPhase == 3 && t > 6.0 && Fluid::Asleep()) {
                            Fluid::SetSpeed(1);
                            TakeScreenshot("bench_oil.png");
                            const bool ok = oilRateOk && before >= 1.0f && oil >= 0.9f * before && Fluid::SourceCount() >= 1; // refilled to its level before the bucket (#603: the 3 m^3 fitted the old seep; seed 42's nearest basin now rests at 1.7)
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
                            const float want = std::max(3.0f, 5.0f * WaterField::CellSize()); // films hold less water with smaller cells (W4c)
                            const bool ok = Fluid::SourceCount() >= 1 && total > want && nearSrc < 0.6f * total;
                            printf("bench: spring: %.1f m^3 of moving water, %.1f of it within 3 m of the spring (a stream: most runs off); films shown/hidden %lld times (#660) -> %s\n",
                                   total, nearSrc, Fluid::GetStats().filmFlips, ok ? "PASS" : "FAIL");
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
                } else if (benchFacing && WorldIdle()) {
                    // --test-facing (#196): parts face the player until R/V. Ground aims from four headings put the front
                    // (local +z) toward the camera; a wall aim puts it out of the wall; a hand rotation sticks while the
                    // camera turns; another part faces the player again.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    player.spawnOn(world, *p, StartSpawnDir(*p));
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                    TerrainHit th;
                    RaycastTerrain(world, Ray{ Vector3Add(Vector3Add(feet, Vector3Scale(flat, 6.0f)), Vector3Scale(player.up, 20.0f)), Vector3Negate(player.up) }, 60.0f, th);
                    if (!p->plots) p->plots = std::make_unique<PlotSet>();
                    Plot* plot = p->plots->plotFor(th.local, 0.0f, true);
                    BuildTool_Place(*p, *plot, { PartType::Cube, { 0, 1, 0 }, 0 }, false); // a block standing on the ground
                    const Camera3D keep = world.camera;
                    const bool flyKeep = player.flyMode;
                    buildTool.active = true;
                    buildTool.part = (int)PartType::Furnace;
                    player.nearest = p; // the build tool aims on the nearest body (set by the player's own update)
                    auto lookFrom = [&](Vector3 eyePlot, Vector3 atPlot) {
                        world.camera.position = p->toWorld(plot->toBody(eyePlot));
                        world.camera.target = p->toWorld(plot->toBody(atPlot));
                        world.camera.up = p->dirToWorld(plot->up);
                        BuildTool_Update(world, PlayerInput::None(), true);
                        return BuildTool_AimOrient();
                    };
                    auto front = [](int o) { const int* m = Orientation(std::max(o, 0)); return Int3{ m[2], m[5], m[8] }; };
                    int good = 0;
                    const Int3 heads[4] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
                    for (const Int3& h : heads) { // the camera 4 m out along h, aiming at the ground 2.5 m out: the front should point back at it
                        const Vector3 at = { 0.25f + 2.5f * h.x, 0.5f, 0.25f + 2.5f * h.z };
                        const int o = lookFrom({ at.x + 2.0f * h.x, 2.2f, at.z + 2.0f * h.z }, at);
                        const Int3 f = front(o);
                        const bool ok = o >= 0 && f.x == h.x && f.y == 0 && f.z == h.z;
                        good += ok;
                        printf("bench: facing: ground, camera toward (%d, 0, %d): orient %d, front (%d, %d, %d) %s\n", h.x, h.z, o, f.x, f.y, f.z, ok ? "ok" : "WRONG");
                    }
                    const int ow = lookFrom({ 3.0f, 0.75f, 0.25f }, { 0.25f, 0.75f, 0.25f }); // the block's +x face
                    const bool wallOk = front(ow).x == 1 && front(ow).y == 0 && front(ow).z == 0;
                    printf("bench: facing: on the block's +x face: front (%d, %d, %d) %s\n", front(ow).x, front(ow).y, front(ow).z, wallOk ? "ok" : "WRONG");
                    buildTool.orientManual = true; // as after R
                    const int om = lookFrom({ 0.25f, 2.2f, -4.0f }, { 0.25f, 0.5f, -2.0f });
                    const bool stuck = om == ow;
                    printf("bench: facing: rotated by hand, then turned: orient %d -> %d %s\n", ow, om, stuck ? "(kept) ok" : "CHANGED");
                    buildTool.part = (int)PartType::Generator; // another part: automatic again
                    const int oa = lookFrom({ 0.25f, 2.2f, -4.0f }, { 0.25f, 0.5f, -2.0f });
                    const bool reset = front(oa).z == -1 && !buildTool.orientManual;
                    printf("bench: facing: another part: front (%d, %d, %d), automatic %s\n", front(oa).x, front(oa).y, front(oa).z, reset ? "ok" : "WRONG");
                    world.camera = keep; player.flyMode = flyKeep; buildTool.active = false;
                    printf("bench: facing: %s\n", good == 4 && wallOk && stuck && reset ? "PASS" : "FAIL");
                    benchFacing = false;
                } else if (benchPlacement && WorldIdle()) {
                    // --test-placement (#550, #551): a furnace (1 m cube) aimed at bare ground 20 m out goes down centred on
                    // the aim point (it grew from its corner there); aimed at the ground beside it, the box's cell does not
                    // flip back and forth while the aim wobbles 3 cm either side of the point where it moves on.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    player.spawnOn(world, *p, StartSpawnDir(*p));
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                    const Vector3 fl = p->dirToLocal(flat);
                    TerrainHit th;
                    RaycastTerrain(world, Ray{ Vector3Add(Vector3Add(feet, Vector3Scale(flat, 20.0f)), Vector3Scale(player.up, 20.0f)), Vector3Negate(player.up) }, 60.0f, th);
                    const Vector3 g = th.local, up = Vector3Normalize(g);
                    const Camera3D keep = world.camera;
                    buildTool.active = true; buildTool.orientManual = false;
                    buildTool.part = (int)PartType::Furnace;
                    player.nearest = p;
                    auto aimAt = [&](Vector3 at) { // from 3 m back along the view and 2 m up
                        const Vector3 eye = Vector3Add(Vector3Subtract(at, Vector3Scale(fl, 3.0f)), Vector3Scale(up, 2.0f));
                        world.camera.position = p->toWorld(eye);
                        world.camera.target = p->toWorld(at);
                        world.camera.up = p->dirToWorld(up);
                        BuildTool_Update(world, PlayerInput::None(), true);
                    };
                    bool ok = true;
                    aimAt(g);
                    Vector3 newOrigin{}; Int3 o{};
                    const bool fresh = BuildTool_AimPlacement(&newOrigin, &o) == nullptr;
                    float off = 1e9f;
                    if (fresh) {
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(newOrigin, 0.0f, true);
                        const Part part{ PartType::Furnace, o, (uint8_t)BuildTool_AimOrient(), 0 };
                        std::vector<std::pair<Int3, CellShape>> fp;
                        Plot::Footprint(part, fp);
                        Vector3 mid{}; int lowY = INT32_MAX;
                        for (const auto& c : fp) { mid = Vector3Add(mid, Vector3{ c.first.x + 0.5f, 0.0f, c.first.z + 0.5f }); lowY = std::min(lowY, c.first.y); }
                        mid = Vector3Scale(mid, (float)kCellSize / (float)fp.size());
                        mid.y = (float)lowY * kCellSize;
                        off = Vector3Distance(plot->toBody(mid), g); // the bottom face's middle from the aim
                        BuildTool_Place(*p, *plot, part, false);
                    }
                    printf("bench: placement: on bare ground, a new plot: %s; the furnace's bottom middle is %.3f m from the aim (#550)\n", fresh ? "yes" : "NO", off);
                    ok = ok && fresh && off < 0.05f;
                    // Beside it (3 m on), step the aim along the view in 1 cm steps until the box moves a cell, then wobble.
                    const Vector3 side = Vector3Normalize(Vector3CrossProduct(fl, up));
                    const Vector3 base = Vector3Add(g, Vector3Scale(side, 3.0f));
                    auto originAt = [&](float s) { aimAt(Vector3Add(base, Vector3Scale(fl, s))); Int3 r{}; BuildTool_AimPlacement(nullptr, &r); return r; };
                    const Int3 o0 = originAt(0.0f);
                    float flip = -1.0f;
                    for (float s = 0.01f; s < 1.0f; s += 0.01f) { const Int3 r = originAt(s); if (!(r == o0)) { flip = s; break; } }
                    int changes = 0;
                    if (flip > 0.0f) {
                        Int3 prev = originAt(flip);
                        for (int k = 0; k < 20; ++k) {
                            const Int3 r = originAt(flip + (k % 2 ? 0.03f : -0.03f));
                            if (!(r == prev)) ++changes;
                            prev = r;
                        }
                    }
                    printf("bench: placement: beside it the box moves a cell %.2f m on; wobbling 3 cm either side of that: %d of 20 aims moved it (#551)\n", flip, changes);
                    ok = ok && flip > 0.0f && changes <= 1; // the first step back may move it once, never back and forth
                    world.camera = keep; buildTool.active = false;
                    printf("bench: placement: %s\n", ok ? "PASS" : "FAIL");
                    benchPlacement = false;
                } else if (benchRoom && (roomPhase > 0 || WorldIdle())) {
                    // --shot-room (#369): a stone room 5 x 5 m, 3 m high, roofed, a doorway 1.5 m wide in one wall; at night, from
                    // inside by the far wall, looking back at the doorway. Unlit, the room should be near black.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (roomPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        roomPhase = 1; roomT0 = GetTime();
                    } else if (roomPhase == 1 && GetTime() - roomT0 > 1.0 && player.grounded) {
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        TerrainHit th;
                        RaycastTerrain(world, Ray{ Vector3Add(Vector3Add(feet, Vector3Scale(flat, 8.0f)), Vector3Scale(player.up, 20.0f)), Vector3Negate(player.up) }, 60.0f, th);
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(th.local, 0.0f, true);
                        roomPlot = plot->id;
                        constexpr int kY = 4; // two metres up: no ground inside the room
                        int placed = 0;
                        auto put = [&](int x, int y, int z) { if (BuildTool_Place(*p, *plot, { PartType::Cube, { x, y + kY, z }, 0 }, false)) ++placed; };
                        for (int z = 0; z < 10; ++z) for (int x = 0; x < 10; ++x) { put(x, 0, z); put(x, 7, z); } // floor and roof
                        for (int y = 1; y <= 6; ++y)
                            for (int k = 0; k < 10; ++k) {
                                put(k, y, 9); put(0, y, k); put(9, y, k);
                                if (!(k >= 4 && k <= 6 && y <= 4)) put(k, y, 0); // the doorway: 1.5 m wide, 2 m high
                            }
                        printf("bench: room: %d parts\n", placed);
                        roomEye = plot->toBody({ 2.5f, (kY + 1) * kCellSize + 1.6f, 4.25f });
                        roomLook = plot->toBody({ 2.5f, (kY + 1) * kCellSize + 1.0f, 0.0f });
                        player.flyMode = true;
                        roomPhase = 2; roomT0 = GetTime();
                    } else if (roomPhase == 2) {
                        world.camera.position = p->toWorld(roomEye); world.camera.target = p->toWorld(roomLook); world.camera.up = p->dirToWorld(Vector3Normalize(roomEye));
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - roomT0 > 3.0) {
                            TakeScreenshot("bench_room.png");
                            printf("bench: room: bench_room.png (sky probes %s)\n", gSkyProbes ? "on" : "off");
                            benchRoom = false;
                        }
                    }
                } else if (benchShelter && (shelterPhase > 0 || WorldIdle())) {
                    // --shot-shelter (#369): a pit 4 m deep with a tunnel off its side under ~2.5 m of rock; at night, from the
                    // pit, looking into the tunnel. Depth alone keeps the tunnel ~75% sky-lit; the overhead rays darken it.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (shelterPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        shelterPhase = 1; shelterT0 = GetTime();
                    } else if (shelterPhase == 1 && GetTime() - shelterT0 > 1.0) {
                        const Vector3 feet = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(feet);
                        const Vector3 fl = p->dirToLocal(player.forward);
                        const Vector3 fw = Vector3Normalize(Vector3Subtract(fl, Vector3Scale(up, Vector3DotProduct(fl, up))));
                        const Vector3 pit = Vector3Add(feet, Vector3Scale(fw, 4.0f));
                        for (int k = 0; k < 3; ++k) ApplyTerrainEdit(*p, TerrainEdit{ Vector3Subtract(pit, Vector3Scale(up, 1.2f * k)), 2.2f, -10.0f });
                        for (int k = 1; k <= 6; ++k) ApplyTerrainEdit(*p, TerrainEdit{ Vector3Add(Vector3Subtract(pit, Vector3Scale(up, 3.0f)), Vector3Scale(fw, 1.6f * k)), 1.4f, -10.0f });
                        shelterEye = Vector3Subtract(Vector3Subtract(pit, Vector3Scale(up, 1.8f)), Vector3Scale(fw, 1.0f));
                        shelterLook = Vector3Add(Vector3Subtract(pit, Vector3Scale(up, 3.2f)), Vector3Scale(fw, 6.0f));
                        player.flyMode = true;
                        shelterPhase = 2; shelterT0 = GetTime();
                    } else if (shelterPhase == 2) {
                        world.camera.position = p->toWorld(shelterEye); world.camera.target = p->toWorld(shelterLook); world.camera.up = p->dirToWorld(Vector3Normalize(shelterEye));
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - shelterT0 > 4.0 && WorldIdle()) {
                            TakeScreenshot("bench_shelter.png");
                            printf("bench: shelter: bench_shelter.png (sky probes %s)\n", gSkyProbes ? "on" : "off");
                            benchShelter = false;
                        }
                    }
                } else if (benchWind && WorldIdle()) {
                    // --test-wind (#650): the field is deterministic, its gusts travel with the mean wind, bursts come more
                    // often in storms, the day cycle runs calm dawn < afternoon, and the wind grows with height.
                    benchWind = false;
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf || !Wind_Enabled()) printf("bench: wind: %s -> SKIP\n", pf ? "RAYLIBGAME_WIND=0" : "needs the planet");
                    else {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 at = Vector3Scale(home, pf->radius + pf->elevation(home));
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        const Vector3 mean = Vector3Scale(t1, 8.0f); // 8 m/s along t1
                        Wind_SetFrame(p, at, mean, 0.0f, 1.0f, 0.0f);
                        const Vector3 q = Vector3Add(at, Vector3Scale(t2, 37.0f));
                        const Vector3 w0 = Wind_At(*p, q, 10.0f, false), w0b = Wind_At(*p, q, 10.0f, false);
                        const bool same = Vector3Distance(w0, w0b) < 1e-6f;
                        // Advection: after 3 s the pattern seen at q is found 24 m downwind.
                        Wind_SetFrame(p, at, mean, 0.0f, 1.0f, 3.0f);
                        const Vector3 w1 = Wind_At(*p, Vector3Add(q, Vector3Scale(t1, 24.0f)), 10.0f, false);
                        const float drift = fabsf(Vector3Length(w1) - Vector3Length(w0));
                        // Bursts: the share of a 400 m square in a burst, calm vs storm; speed range.
                        auto burstShare = [&](float storm, float& lo, float& hi) {
                            Wind_SetFrame(p, at, mean, storm, 1.0f, 0.0f);
                            int n = 0, b = 0; lo = 1e9f; hi = 0.0f;
                            for (float x = -200.0f; x <= 200.0f; x += 4.0f)
                                for (float y = -200.0f; y <= 200.0f; y += 4.0f) {
                                    const Vector3 s = Vector3Add(at, Vector3Add(Vector3Scale(t1, x), Vector3Scale(t2, y)));
                                    ++n; b += Wind_Gust01(*p, s) > 0.5f;
                                    const float v = Vector3Length(Wind_At(*p, s, 10.0f, false)); lo = std::min(lo, v); hi = std::max(hi, v);
                                }
                            return (float)b / (float)n;
                        };
                        float cLo, cHi, sLo, sHi;
                        const float calmShare = burstShare(0.0f, cLo, cHi), stormShare = burstShare(1.0f, sLo, sHi);
                        // The day: the sun rising in the east, mid-afternoon, at night (the spin axis is +y).
                        const Vector3 up = home, east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, up));
                        auto sunAt = [&](float th) { return Vector3Add(Vector3Scale(east, cosf(th)), Vector3Scale(up, sinf(th))); };
                        const float dawn = Wind_DailyFactor(up, sunAt(0.02f), { 0, 1, 0 }), afternoon = Wind_DailyFactor(up, sunAt(0.7f * PI), { 0, 1, 0 }), night = Wind_DailyFactor(up, sunAt(1.5f * PI), { 0, 1, 0 });
                        Wind_SetFrame(p, at, mean, 0.0f, 1.0f, 0.0f);
                        const float low = Vector3Length(Wind_At(*p, q, 1.0f, false)), high = Vector3Length(Wind_At(*p, q, 20.0f, false));
                        const float terrain = Vector3Length(Wind_At(*p, q, 10.0f, true));
                        const bool ok = same && drift < 0.15f && stormShare > calmShare && calmShare > 0.0f && calmShare < 0.2f &&
                                        dawn < night && night < afternoon && high > low && cLo > 0.5f && cHi > 12.0f && terrain > 0.0f;
                        printf("bench: wind: deterministic %d; gust 24 m downwind 3 s later differs %.3f m/s; bursts calm %.1f%% storm %.1f%% of the ground; "
                               "speed calm %.1f..%.1f storm %.1f..%.1f m/s; day dawn %.2f afternoon %.2f night %.2f; at 1 m %.1f, 20 m %.1f m/s; with terrain %.1f -> %s\n",
                               (int)same, drift, 100.0f * calmShare, 100.0f * stormShare, cLo, cHi, sLo, sHi, dawn, afternoon, night, low, high, terrain, ok ? "PASS" : "FAIL");
                    }
                } else if (benchCampfire && (campfirePhase > 0 || WorldIdle())) {
                    // --shot-campfire (#552): a lit campfire 2.5 m from the camera, at dusk (bench_campfire.png).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (!p) benchCampfire = false;
                    else if (campfirePhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        campfirePhase = 1; campfireT0 = GetTime();
                    } else if (campfirePhase == 1 && GetTime() - campfireT0 > 1.5 && player.grounded) {
                        const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                        const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* plot = p->plots->plotFor(p->toLocal(Vector3Add(Vector3Add(feet, Vector3Scale(flat, 4.0f)), Vector3Scale(player.up, 1.0f))), 0.0f, true);
                        const int fire = plot->place({ PartType::Campfire, { 0, 0, 0 }, 0 });
                        Factory::Sync(world);
                        Factory::WaitTicks(3);
                        if (fire >= 0) { Factory::Seed(p->id, plot->id, fire, "oak_log", 5); Factory::Seed(p->id, plot->id, fire, "oak_log", 3); }
                        campfirePlot = plot->id; campfirePart = fire;
                        campfireAt = plot->toBody({ 0.5f, 0.0f, 0.5f });
                        campfireOut = p->dirToLocal(flat);
                        const Vector3 up = Vector3Normalize(campfireAt);
                        sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(campfireOut, -3.0f)))); // low sun behind the camera
                        printf("bench: campfire: placed part %d\n", fire);
                        campfirePhase = 2; campfireT0 = GetTime();
                    } else if (campfirePhase == 2) {
                        const Vector3 up = Vector3Normalize(campfireAt);
                        const Vector3 eye = Vector3Add(campfireAt, Vector3Add(Vector3Scale(campfireOut, -2.3f), Vector3Scale(up, 1.3f)));
                        player.flyMode = true;
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(Vector3Add(campfireAt, Vector3Scale(up, 0.25f))); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - campfireT0 > 3.0) {
                            TakeScreenshot("bench_campfire.png");
                            printf("bench: campfire: working %d -> bench_campfire.png\n", (int)Factory::Working(p->id, campfirePlot, campfirePart));
                            benchCampfire = false;
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
                } else if (benchMountain && (mountainPhase > 0 || WorldIdle())) {
                    // --shot-mountain (#545): the highest ground within 8 km of the start, seen from 1.5 km away and 250 m
                    // over it, sun from the side (bench_mountain.png): the mid/far LOD transitions lie across its slopes.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { printf("bench: mountain: needs the planet\n"); benchMountain = false; }
                    else if (mountainPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float bestE = -1e30f; Vector3 bestD = home;
                        for (float a = -8000.0f; a <= 8000.0f; a += 100.0f)
                            for (float b = -8000.0f; b <= 8000.0f; b += 100.0f) {
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, a / pf->radius), Vector3Scale(t2, b / pf->radius))));
                                const float e = pf->elevation(d);
                                if (e > bestE) { bestE = e; bestD = d; }
                            }
                        mountainAt = Vector3Scale(bestD, pf->radius + bestE);
                        mountainOut = Vector3Normalize(Vector3Subtract(Vector3Scale(home, pf->radius), Vector3Scale(bestD, pf->radius))); // from the peak toward the start
                        if (Vector3Length(mountainOut) < 0.5f) mountainOut = t1;
                        printf("bench: mountain: highest ground within 8 km: %.0f m\n", bestE);
                        player.flyMode = true; sunOverride = true;
                        const Vector3 up = Vector3Normalize(mountainAt);
                        sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(Vector3CrossProduct(up, mountainOut), 1.2f))));
                        mountainPhase = 1; mountainT0 = GetTime();
                    } else {
                        const Vector3 up = Vector3Normalize(mountainAt);
                        const Vector3 eye = Vector3Add(mountainAt, Vector3Add(Vector3Scale(mountainOut, 1500.0f), Vector3Scale(up, 250.0f)));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(Vector3Subtract(mountainAt, Vector3Scale(up, 150.0f))); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        if (GetTime() - mountainT0 > 5.0 && WorldIdle()) {
                            TakeScreenshot("bench_mountain.png");
                            printf("bench: mountain: bench_mountain.png\n");
                            benchMountain = false;
                        }
                    }
                } else if (benchCaveMouth && (mouthPhase > 0 || WorldIdle())) {
                    // --shot-cave-mouth (#549): the nearest cave mouth to the start (open ground 2 m and 6 m under the
                    // surface), how open it is per LOD footprint, and shots from 20, 80 and 250 m (bench_mouth_<m>.png).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - mouthT0;
                    static const float kMouthDist[3] = { 20.0f, 80.0f, 250.0f };
                    if (!pf) { printf("bench: cave mouth: needs the planet\n"); benchCaveMouth = false; }
                    else if (mouthPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f; int mouths = 0, probes = 0;
                        for (float a = -3000.0f; a <= 3000.0f; a += 20.0f)
                            for (float b = -3000.0f; b <= 3000.0f; b += 20.0f) {
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, a / pf->radius), Vector3Scale(t2, b / pf->radius))));
                                const float e = pf->elevation(d);
                                if (e < 2.0f) continue;
                                ++probes;
                                const float s = pf->radius + e;
                                if (pf->baseDensity(Vector3Scale(d, s - 2.0f)) >= ISO_LEVEL || pf->baseDensity(Vector3Scale(d, s - 6.0f)) >= ISO_LEVEL) continue;
                                ++mouths;
                                const float dist = sqrtf(a * a + b * b);
                                if (dist < best) { best = dist; mouthAt = Vector3Scale(d, s); }
                            }
                        printf("bench: cave mouth: %d of %d land columns within 3 km open 2-6 m down; nearest %.0f m from the start\n", mouths, probes, best);
                        if (best > 1e29f) { printf("bench: cave mouth: none found -> FAIL\n"); benchCaveMouth = false; }
                        else {
                            const Vector3 up = Vector3Normalize(mouthAt);
                            // Look in from the side the mouth opens to: the horizontal direction with the most open air 3 m down.
                            const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), e2 = Vector3CrossProduct(up, e1);
                            float bestOpen = -1.0f;
                            for (int k = 0; k < 16; ++k) {
                                const float ang = 2.0f * PI * k / 16.0f;
                                const Vector3 h = Vector3Add(Vector3Scale(e1, cosf(ang)), Vector3Scale(e2, sinf(ang)));
                                int open = 0;
                                for (float r = 2.0f; r <= 30.0f; r += 2.0f) open += pf->baseDensity(Vector3Add(mouthAt, Vector3Add(Vector3Scale(h, r), Vector3Scale(up, -3.0f)))) < ISO_LEVEL;
                                if (open > bestOpen) { bestOpen = (float)open; mouthOut = h; }
                            }
                            // How open the mouth is to the meshes of each LOD footprint: samples in a 12 m disc, 2 m down.
                            printf("bench: cave mouth: open samples 2 m down by footprint:");
                            for (float fp : { 0.25f, 1.0f, 2.0f, 4.0f, 8.0f, 12.0f, 16.0f }) {
                                int open = 0, n = 0;
                                for (float x = -12.0f; x <= 12.0f; x += 2.0f)
                                    for (float y = -12.0f; y <= 12.0f; y += 2.0f) {
                                        if (x * x + y * y > 144.0f) continue;
                                        const Vector3 q = Vector3Add(mouthAt, Vector3Add(Vector3Add(Vector3Scale(e1, x), Vector3Scale(e2, y)), Vector3Scale(up, -2.0f)));
                                        const Vector3 qd = Vector3Normalize(q);
                                        const float surf = pf->radius + pf->elevation(qd, fp);
                                        ++n; open += pf->baseDensityLod(Vector3Scale(qd, surf - 2.0f), fp) < ISO_LEVEL;
                                    }
                                printf(" fp%g %d/%d", fp, open, n);
                            }
                            printf("\n");
                            player.flyMode = true; sunOverride = true;
                            sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(up, Vector3Scale(mouthOut, 0.5f)))); // light into the mouth
                            mouthPhase = 1; mouthT0 = GetTime();
                        }
                    } else {
                        const Vector3 up = Vector3Normalize(mouthAt);
                        const float r = kMouthDist[mouthPhase - 1];
                        const Vector3 eye = Vector3Add(mouthAt, Vector3Add(Vector3Scale(mouthOut, r), Vector3Scale(up, 4.0f + 0.35f * r)));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(Vector3Add(mouthAt, Vector3Scale(up, -3.0f))); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        if (t > 4.0 && WorldIdle()) {
                            TakeScreenshot(TextFormat("bench_mouth_%d.png", (int)r));
                            printf("bench: cave mouth: bench_mouth_%d.png\n", (int)r);
                            if (++mouthPhase > 3) benchCaveMouth = false;
                            mouthT0 = GetTime();
                        }
                    }
                } else if (benchCaveWater && (caveWaterPhase > 0 || WorldIdle())) {
                    // --test-cave-water (#559): a cave hollowed in rock on a river's bed (roofed, under the water's level) with a torch in it
                    // stays dry; a hole dug at the bed, open to the river, fills from it.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { printf("bench: cave water: needs the planet -> FAIL\n"); benchCaveWater = false; }
                    else if (caveWaterPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f; Vector3 wdir{}; float level = 0.0f;
                        for (float span = 3000.0f; span <= 12000.0f && best > 1e29f; span *= 2.0f) // out to 12 km (coarser) if need be
                        for (float e = -span; e <= span; e += 16.0f * span / 3000.0f)
                            for (float n = -span; n <= span; n += 16.0f * span / 3000.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if ((ws.kind != 1 && ws.kind != 2) || ws.level <= 2.0f || pf->elevation(d) > ws.level - 7.0f) continue; // a river or a lake 7 m deep or more: over the 6.5 m rock below (#623: with bake-native lakes the nearest 4 m site was 4.7 m deep, the rock stood out of the water and its tunnel stayed dry; 9 m deep in the same lake it fills)
                                // ...with an open channel: no ground standing in the water over the height map's bed.
                                const float eb = pf->elevation(d);
                                bool open = true;
                                for (float hgt = eb + 0.5f; hgt < ws.level - 0.3f && open; hgt += 0.75f)
                                    open = pf->density(Vector3Scale(d, pf->radius + hgt), p->edits.get()) < ISO_LEVEL;
                                if (!open) continue;
                                // ...and calm and wide: the same water 15 m away all round (no fall, #540), 1.5 m deep or more
                                // there (a narrow mountain stream left the test's tunnel in the air; any generated water had the bug).
                                bool calm = true;
                                const Vector3 probes[4] = { Vector3Scale(t1, 15.0f), Vector3Scale(t1, -15.0f), Vector3Scale(t2, 15.0f), Vector3Scale(t2, -15.0f) };
                                for (const Vector3& o : probes) {
                                    const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(o, 1.0f / pf->radius)));
                                    const PlanetField::WaterSource wq = pf->waterSourceAt(q);
                                    calm = calm && wq.kind == ws.kind && fabsf(wq.level - ws.level) < 0.3f && pf->elevation(q) < ws.level - 1.5f;
                                }
                                if (!calm) continue;
                                best = dist; wdir = d; level = ws.level;
                            }
                        if (best > 1e29f) { printf("bench: cave water: no calm, deep river or lake within 12 km of the start -> FAIL\n"); benchCaveWater = false; }
                        else {
                            // The water is simulated around the player: go there first (the river can be km from the start).
                            caveRiverDir = wdir; caveLevel = level; caveDist = best;
                            player.spawnOn(world, *p, Vector3Normalize(Vector3Add(wdir, Vector3Scale(t1, 30.0f / pf->radius))));
                            caveWaterPhase = 1; caveWaterT0 = GetTime();
                        }
                    } else if (caveWaterPhase == 1) {
                        if (GetTime() - caveWaterT0 > 4.0 && WorldIdle()) {
                            const Vector3 wdir = caveRiverDir;
                            const float level = caveLevel, best = caveDist;
                            const Vector3 home = pf->findSpawnDir(p->seed);
                            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 }));
                            const float bed = pf->elevation(wdir);
                            // Rock on the bed (a boulder, its top under the surface: the site is 7 m deep), hollowed: a cave inside it whose
                            // cells lie between the bed and the river's level - the generated water's band - roofed by rock.
                            TerrainEdit rock{ Vector3Scale(wdir, pf->radius + bed + 2.5f), 4.0f, 0.0f };
                            rock.fill = true;
                            ApplyTerrainEdit(*p, rock);
                            // Two caves in the rock, across the river from each other: one roofed and walled (the torch in it),
                            // the control tunnelled out to the river under its surface - the river pours in through the fluid
                            // simulation. Same rock, same water (a control elsewhere read the released reach's drainage, #361).
                            const RiverNetwork::Hit rh = pf->rivers() ? pf->rivers()->nearest(wdir, 20.0f) : RiverNetwork::Hit{};
                            const Vector3 along = rh.found ? rh.down : t1;
                            const Vector3 side = Vector3Normalize(Vector3CrossProduct(wdir, along));
                            auto offset = [&](float s, float h) { return Vector3Scale(Vector3Normalize(Vector3Add(wdir, Vector3Scale(side, s / pf->radius))), pf->radius + bed + h); };
                            caveAt = offset(-1.6f, 1.2f);
                            ApplyTerrainEdit(*p, TerrainEdit{ caveAt, 1.0f, -1.0f, true });
                            caveOpenAt = offset(1.6f, 1.2f);
                            ApplyTerrainEdit(*p, TerrainEdit{ caveOpenAt, 1.0f, -1.0f, true });
                            for (float s = 2.2f; s <= 6.0f; s += 0.5f) ApplyTerrainEdit(*p, TerrainEdit{ offset(s, 1.2f), 1.0f, -1.0f, true }); // the tunnel out: 2 m across (water cells are 1 m)
                            if (!p->plots) p->plots = std::make_unique<PlotSet>();
                            const Vector3 torchAt = Vector3Subtract(caveAt, Vector3Scale(wdir, 0.6f)); // on the cave's floor
                            Plot* plot = p->plots->plotFor(torchAt, 0.0f, true);
                            plot->place({ PartType::Torch, plot->cellAt(plot->toPlot(torchAt)), 0 });
                            Fluid::OnTerrainEdit(*p, caveAt, 6.0f); // as placing a part does: the bricks around it are reclassified
                            printf("bench: cave water: river or lake %.1f km from the start, water at %.1f m over a bed at %.1f m; a rock on the bed with two caves under the water's level, one roofed (a torch in it), one tunnelled out to the river\n",
                                   best / 1000.0f, level, bed);
                            caveWaterPhase = 2; caveWaterT0 = GetTime();
                        }
                    } else if (GetTime() - caveWaterT0 > 8.0) {
                        const float inCave = Fluid::WaterAt(*p, caveAt, 0.9f), inHole = Fluid::WaterAt(*p, caveOpenAt, 0.9f);
                        const bool ok = inCave < 0.05f && inHole > 1.0f;
                        {   // what the generated water is at the tunnel's mouth and beyond it (#623 diagnosis)
                            const Vector3 up = Vector3Normalize(caveOpenAt), t = Vector3Normalize(Vector3Subtract(caveOpenAt, caveAt));
                            for (float s = 3.0f; s <= 12.0f; s += 3.0f) {
                                const Vector3 q = Vector3Normalize(Vector3Add(up, Vector3Scale(t, s / pf->radius)));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(q);
                                printf("bench: cave water: %.0f m past the cave: source kind %d level %.2f, ground %.2f, simulated water %.2f m^3\n",
                                       s, ws.kind, ws.level, pf->elevation(q), Fluid::WaterAt(*p, Vector3Scale(q, Vector3Length(caveOpenAt)), 0.9f));
                            }
                        }
                        printf("bench: cave water: %.2f m^3 of water in the roofed cave (want ~0), %.2f m^3 in the one tunnelled out to the river (want > 1) -> %s\n",
                               inCave, inHole, ok ? "PASS" : "FAIL");
                        benchCaveWater = false;
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
                            plot->place({ benchPumpHand ? PartType::HandPump : PartType::WaterPump, pc, 0 }); // --test-pump-hand: the hand pump (#280)
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
                            const int outlet = benchPumpBoiler ? -1 : plot->place({ PartType::Outlet, { c.x, c.y - 1, c.z }, 0 });
                            if (benchPumpBoiler) { // a boiler under the pipe's end instead of the outlet, filled to the brim (four buckets)
                                const Int3 nbrs[6] = { { c.x, c.y - 1, c.z }, { c.x + 1, c.y, c.z }, { c.x - 1, c.y, c.z }, { c.x, c.y, c.z + 1 }, { c.x, c.y, c.z - 1 }, { c.x, c.y + 1, c.z } };
                                for (const Int3& at : nbrs) if (pumpBoilerPart < 0) pumpBoilerPart = plot->place({ PartType::Boiler, at, 0 }); // its origin cell touches the pipe's end
                                if (pumpBoilerPart >= 0) Factory::Seed(p->id, plot->id, pumpBoilerPart, "water_bucket", 4);
                                printf("bench: pump: boiler %s under the pipe's end\n", pumpBoilerPart >= 0 ? "placed" : "BLOCKED");
                            }
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
                        if (benchPumpHand && pumpPhase == 1) // working the hand pump (RMB held), until the "power" is cut
                            if (Plot* hp = p->plots ? p->plots->byId(pumpPlot) : nullptr)
                                for (int i = 0; i < (int)hp->parts.size(); ++i)
                                    if (hp->parts[i].type == PartType::HandPump) Waterworks_HandPump(p->id, pumpPlot, i, GameFrameTime());
                        if (t - lastPrint >= 5.0) {
                            lastPrint = t;
                            Plot* plot = p->plots ? p->plots->byId(pumpPlot) : nullptr;
                            int pumpPart = -1;
                            if (plot) for (int i = 0; i < (int)plot->parts.size(); ++i) if (plot->parts[i].type == PartType::WaterPump || plot->parts[i].type == PartType::HandPump) pumpPart = i;
                            const MachineView v = pumpPart >= 0 && !benchPumpBoiler ? Factory::View(p->id, pumpPlot, pumpPart) : MachineView{};
                            printf("bench: pump t=%3.0f s: %s, %d networks, pumped %.1f m^3 (%d dry), bowl holds %.1f m^3\n", t, v.valid ? MachineStatusName(v.status) : "-", ws.networks, ws.pumped, ws.dry,
                                   Fluid::WaterAt(*p, pumpBasin, 2.5f));
                        }
                        static MachineView bv; // the boiler, watched every frame (Factory::View fills from the tick after it is asked)
                        if (benchPumpBoiler && pumpBoilerPart >= 0) { const MachineView v2 = Factory::View(p->id, pumpPlot, pumpBoilerPart); if (v2.valid) bv = v2; }
                        if (benchPumpBoiler && t > 20.0) {
                            const bool pass = pumpBoilerPart >= 0 && ws.pumped <= 0.05 && bv.water == bv.waterMax && bv.waterMax > 0;
                            printf("bench: pump-boiler: a full boiler on the network; the pump took %.2f m^3 from the world in 20 s, the boiler holds %d of %d L -> %s\n",
                                   ws.pumped, bv.water, bv.waterMax, pass ? "PASS" : "FAIL");
                            benchPump = false;
                        } else if (pumpPhase == 1 && t > 60.0) {
                            TakeScreenshot("bench_pump.png");
                            // Cut the power.
                            if (Plot* plot = p->plots ? p->plots->byId(pumpPlot) : nullptr) plot->remove(pumpCable);
                            pumpedAtCut = ws.pumped;
                            pumpPhase = 2;
                        } else if (pumpPhase == 2 && t > 75.0) {
                            const double after = ws.pumped - pumpedAtCut;
                            const float bowl = Fluid::WaterAt(*p, pumpBasin, 2.5f);
                            const bool pass = benchPumpHand ? (pumpedAtCut >= 2.0 && bowl >= 1.0f && after <= 0.5) : (pumpedAtCut >= 10.0 && bowl >= 5.0f && after <= 1.0); // by hand: 50 L/s
                            printf("bench: pump: pumped %.1f m^3 in 60 s (%.1f m up), the bowl holds %.1f m^3; after cutting the cable %.1f m^3 more -> %s\n",
                                   pumpedAtCut, pumpRise, bowl, after, pass ? "PASS" : "FAIL");
                            benchPump = false;
                        }
                    }
                } else if (benchLeak && (lkPhase > 0 || WorldIdle())) {
                    // --test-leak (#680, best with --sandbox; user: 'would the column change still allow a hole in a wall to leak?'):
                    // a 20 m reservoir 4 m wide, 2 m deep behind a 1 m wall of ground, and an empty channel 6 m deep beyond; a 1 x 1 m hole is
                    // cut through the wall 0.5 m over the floor. The early outflow (the volume past the wall, 1-3 s) against an
                    // orifice's sqrt(2 g h) x area (a real one passes ~0.6 of it), and after 60 s the pool against an orifice-then-weir model of the drain to the hole's
                    // bottom (0.5 m deep).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    static Vector3 home{}, ax{}, ay{};
                    static float r0 = 0.0f, h0 = 0.0f, vol1 = -1.0f, vol3 = -1.0f, pastBefore = 0.0f;
                    const double t = SimTime() - lkT0;
                    auto pt = [&](float x, float y, float h) {
                        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
                        return Vector3Scale(d, r0 + h);
                    };
                    auto box = [&](float x0, float x1, float y0, float y1, float h0b, float h1b) {
                        TerrainEdit e{};
                        e.center = pt(0.5f * (x0 + x1), 0.5f * (y0 + y1), 0.5f * (h0b + h1b));
                        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 0.5f * (y1 - y0), 0.5f * (h1b - h0b) };
                        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
                        e.radius = Vector3Length(e.half);
                        return e;
                    };
                    auto depthAt = [&](float x) { const float s = Fluid::SurfaceAbove(*p, pt(x, 0.0f, -2.5f + 0.55f), 4.0f); return s > 0.0f ? s + 0.55f : 0.0f; };
                    auto past = [&]() { float v = 0.0f; for (float x = 4.0f; x < 40.0f; x += 4.0f) for (float h : { -5.0f, -1.0f }) v += Fluid::WaterAt(*p, pt(x, 0.0f, h), 2.0f); return v; };
                    if (!pf) { printf("bench: leak: needs the planet -> FAIL\n"); benchLeak = false; }
                    else {
                        if (lkPhase == 0) {
                            home = pf->findSpawnDir(p->seed);
                            ax = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })); ay = Vector3CrossProduct(ax, home);
                            r0 = pf->radius + pf->elevation(home);
                            if (benchLeakParts) {
                                // #692: the frame is a build plot's (its east and north along the ground, its origin on it), so the
                                // blocks sit exactly where the trench is dug.
                                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                                lkPlot = p->plots->plotFor(Vector3Scale(home, r0), 0.0f, true);
                                const Vector3 o = lkPlot->toBody({ 0.0f, 0.0f, 0.0f });
                                home = Vector3Normalize(o); r0 = Vector3Length(o); ax = lkPlot->east; ay = lkPlot->north;
                            }
                            std::vector<TerrainEdit> cut = benchLeakParts
                                ? std::vector<TerrainEdit>{ box(-20.0f, 2.0f, -2.0f, 2.0f, -2.5f, 3.0f), box(2.0f, 40.0f, -2.0f, 2.0f, -6.0f, 3.0f) } // the dam's slot dug too: blocks fill it
                                : std::vector<TerrainEdit>{ box(-20.0f, 0.0f, -2.0f, 2.0f, -2.5f, 3.0f), box(2.0f, 40.0f, -2.0f, 2.0f, -6.0f, 3.0f) }; // up through the air (#693: a top at the ground's height left slivers where it stood higher); a 2 m wall (#691: 1 m - one cell at 1 m cells - leaked); the channel 6 m deep: it never backs up to the hole
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            if (benchLeakParts) {
                                // The dam: 1 m stone blocks, 2 thick (x 0-2 m), from 1 m under the floor to 0.5 m over the ground
                                // (plot cells are 0.5 m: y -7..-1 up), 1 m into each bank (z -3..3 m). The hole's two blocks (1-2 m over
                                // the floor, z 0-1 m) come out later.
                                const uint8_t stone = (uint8_t)std::max(0, BuildMaterialForItem("stone"));
                                int placed = 0;
                                lkHole.clear();
                                for (int cx : { 0, 2 })
                                    for (int cy = -7; cy <= -1; cy += 2)
                                        for (int cz = -6; cz <= 4; cz += 2) {
                                            const int idx = lkPlot->place({ PartType::Block, { cx, cy, cz }, 0, stone });
                                            if (idx < 0) continue;
                                            ++placed;
                                            if (cy == -3 && cz == 0) lkHole.push_back(idx);
                                        }
                                printf("bench: leak: a 20 m reservoir behind a 2 m dam of %d placed stone blocks, a channel beyond (cells %.2f m)\n", placed, WaterField::CellSize());
                            } else printf("bench: leak: a 20 m reservoir behind a 2 m wall, a channel beyond (cells %.2f m)\n", WaterField::CellSize());
                            lkPhase = 1; lkT0 = SimTime();
                        }
                        // The camera 45 degrees off the channel on its downstream side, 45 degrees up (#691, user: 'see both the
                        // trench before and after the dam, the hole in the dam, and the edges of the dam'): the hole's outlet faces it.
                        world.camera.position = p->toWorld(pt(6.0f, -6.0f, 6.0f)); world.camera.target = p->toWorld(pt(0.5f, 0.0f, -1.5f));
                        world.camera.up = p->dirToWorld(home);
                        player.flyMode = true; player.position = world.renderOrigin + world.camera.position;
                        // Water on the grass beside the reservoir and the channel (#691): a bad fill, or a leak around the edges.
                        // (3.5 m out, 0.5 m round: clear of the dug edge, which the lattice rounds out to ~2.5 m at 1 m cells)
                        auto onGrass = [&]() { float v = 0.0f; for (float x = -19.0f; x < 40.0f; x += 3.0f) for (float y : { -3.5f, 3.5f }) v += Fluid::WaterAt(*p, pt(x, y, 0.25f), 0.5f); return v; };
                        if (lkPhase == 1 && t > 1.0 && WorldIdle()) {
                            // Fill at 5 m^3/s, 1 m^3 at a time moving over the whole floor (#691: a pour every frame - 60 m^3/s -
                            // piled mounds over the rim faster than water spreads, and they spilled onto the grass).
                            static int k = 0;
                            static double nextPour = 0.0;
                            if (t >= nextPour) {
                                nextPour = t + 0.2;
                                const float x = -19.0f + 2.0f * (float)(k % 10), y = (k / 10) % 2 ? 1.0f : -1.0f;
                                ++k;
                                lkPoured += Fluid::Pour(*p, pt(x, y, -2.3f), 1.0f);
                            }
                            if (lkPoured >= 160.0f || depthAt(-10.0f) >= 2.0f) { lkPhase = 2; lkT0 = SimTime(); }
                        } else if (lkPhase == 2 && ((t > 8.0 && Fluid::Asleep()) || t > 40.0)) {
                            const float grass = onGrass();
                            printf("bench: leak: after the fill %.2f m^3 on the grass beside the trench%s\n", grass, grass > 0.05f ? " -> FAIL (a bad fill)" : "");
                            h0 = depthAt(-10.0f);
                            pastBefore = past(); // the wall (or its seams with the ground) must hold until the hole opens
                            if (benchLeakParts) {
                                for (int idx : lkHole) lkPlot->remove(idx); // the hole: two blocks out, 1 x 1 m through the dam, 1-2 m over the floor
                            } else {
                                std::vector<TerrainEdit> cut = { box(-0.2f, 2.2f, -0.5f, 0.5f, -2.0f, -1.0f) }; // the hole: 1 x 1 m through the wall, 0.5-1.5 m over the floor
                                ApplyTerrainEdits(*p, cut.data(), cut.size());
                            }
                            printf("bench: leak: pool %.2f m deep (%.1f m^3 already past the wall); a 1 x 1 m hole opened %.1f m over the floor%s\n", h0, past(), benchLeakParts ? 1.0f : 0.5f,
                                   benchLeakParts ? TextFormat(" (%zu blocks taken out)", lkHole.size()) : "");
                            lkPhase = 3; lkT0 = SimTime();
                        } else if (lkPhase == 3) {
                            if (t > 1.0 && vol1 < 0.0f) vol1 = past();
                            if (t > 3.0 && vol3 < 0.0f) vol3 = past();
                            static double lastPrint = 0.0;
                            if (t - lastPrint >= 20.0) {
                                lastPrint = t;
                                const Fluid::Stats fs = Fluid::GetStats();
                                printf("bench: leak t=%3.0f s: pool %.2f m deep, %.1f m^3 past the wall; %d of %d bricks awake, %d steady, %d settled\n",
                                       t, depthAt(-10.0f), past(), fs.activeBricks, fs.bricks, fs.steadyBricks, fs.settled);
                            }
                            constexpr float kDrain = 60.0f; // s judged after the hole opens (was 120: the test ran over 3 min at 0.25 m cells)
                            if (t > kDrain) {
                                const float lip = benchLeakParts ? 1.0f : 0.5f; // the hole's bottom over the floor (blocks: block-aligned)
                                const float g = 9.81f, head = std::max(0.0f, h0 - (lip + 0.5f)); // over the hole's centre
                                const float qIdeal = sqrtf(2.0f * g * head) * 1.0f, q = (vol3 - vol1) / 2.0f, final = depthAt(-10.0f);
                                // The expected depth after kDrain: the pool (its area as dug: the water poured over its depth) through the hole as an orifice while
                                // it is under water (0.6 sqrt(2 g h) x 1 m^2, h over its centre), then a weir over its lip (1.7 x 1 m x h^1.5) - it slows as it
                                // drains: some water over the lip is left at the end, not none.
                                const float poolArea = h0 > 0.1f ? lkPoured / h0 : 80.0f; // ~80 m^2; the lattice rounds the dug box (#691)
                                float d = h0;
                                for (int k = 0; k < (int)(kDrain / 0.01f); ++k) {
                                    const float qd = d >= lip + 1.0f ? 0.6f * sqrtf(2.0f * g * (d - (lip + 0.5f))) : d > lip ? 1.7f * powf(d - lip, 1.5f) : 0.0f;
                                    d -= qd * 0.01f / poolArea;
                                }
                                const float grass = onGrass();
                                const bool ok = h0 > 1.5f && q > 0.4f * qIdeal && q < 1.6f * qIdeal && fabsf(final - d) < 0.12f && grass < 0.05f && pastBefore < 0.2f;
                                printf("bench: leak: outflow %.2f m^3/s early (orifice ideal %.2f: %.2f of it; a real one ~0.6); after %.0f s the pool is %.2f m deep (an orifice then a weir over the lip: %.2f); %.2f m^3 on the grass -> %s\n",
                                       q, qIdeal, q / std::max(qIdeal, 0.01f), kDrain, final, d, grass, ok ? "PASS" : "FAIL");
                                benchLeak = false;
                            }
                        }
                    }
                } else if (benchDamBreak && (dbPhase > 0 || WorldIdle())) {
                    // --test-dam-break (#680, best with --sandbox): a physics reference. A channel 4 m wide and 2.5 m deep; its first
                    // 20 m hold water 2 m deep behind a 0.5 m wall of ground, which is cut away. Ritter's solution for a dam break
                    // on a dry, frictionless bed: the front runs at 2 sqrt(g h0), the depth at the dam is 4/9 h0 and the
                    // discharge (8/27) sqrt(g) h0^1.5 per metre of width. Measured in the first 3 s (the wave reflected off the
                    // reservoir's back reaches the dam at ~4.5 s). Friction makes a real front slower; the ratios are what to
                    // compare across cell sizes.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    static Vector3 home{}, ax{}, ay{};
                    static float r0 = 0.0f, h0 = 0.0f, depthDam = -1.0f, qDam = 0.0f;
                    static double front[3] = { -1.0, -1.0, -1.0 };
                    const double t = SimTime() - dbT0;
                    auto pt = [&](float x, float y, float h) { // x along the channel, y across, h over the ground (m)
                        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
                        return Vector3Scale(d, r0 + h);
                    };
                    auto box = [&](float x0, float x1) { // the channel's cut from x0 to x1, 4 m wide, 2.5 m deep
                        TerrainEdit e{};
                        e.center = pt(0.5f * (x0 + x1), 0.0f, 0.25f); // up through the air to 3 m (#693: slivers where the ground stood higher)
                        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 2.0f, 2.75f };
                        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
                        e.radius = Vector3Length(e.half);
                        return e;
                    };
                    if (!pf) { printf("bench: dam-break: needs the planet -> FAIL\n"); benchDamBreak = false; }
                    else {
                        if (dbPhase == 0) {
                            home = pf->findSpawnDir(p->seed);
                            ax = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })); ay = Vector3CrossProduct(ax, home);
                            r0 = pf->radius + pf->elevation(home);
                            std::vector<TerrainEdit> cut = { box(-20.0f, 0.0f), box(2.0f, 50.0f) }; // a 2 m wall (#691: 1 m - one cell at 1 m cells - leaked)
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            printf("bench: dam-break: a 20 m reservoir behind a 1 m wall, a 50 m channel, 4 m wide, 2.5 m deep (cells %.2f m)\n", WaterField::CellSize());
                            dbPhase = 1; dbT0 = SimTime();
                        }
                        // The camera over the wall, looking down at the reservoir and the channel's start (the subject, CLAUDE.md).
                        world.camera.position = p->toWorld(pt(-4.0f, -9.0f, 9.0f)); world.camera.target = p->toWorld(pt(-4.0f, 0.0f, -2.0f));
                        world.camera.up = p->dirToWorld(home);
                        player.flyMode = true; player.position = world.renderOrigin + world.camera.position;
                        if (dbPhase == 1 && t > 1.0 && WorldIdle()) {
                            // Fill the reservoir to 2 m: 160 m^3 at 5 m^3/s, 1 m^3 at a time moving over the floor (#691: a pour every
                            // frame - 60 m^3/s - piled mounds over the rim faster than water spreads, and they spilled onto the grass).
                            static int k = 0;
                            static double nextPour = 0.0;
                            if (t >= nextPour) {
                                nextPour = t + 0.2;
                                const float x = -19.0f + 2.0f * (float)(k % 10), y = (k / 10) % 2 ? 1.0f : -1.0f;
                                ++k;
                                dbPoured += Fluid::Pour(*p, pt(x, y, -2.3f), 1.0f);
                            }
                            if (dbPoured >= 160.0f || [&] { const float s = Fluid::SurfaceAbove(*p, pt(-10.0f, 0.0f, -2.5f + 0.55f), 4.0f); return s > 0.0f && s + 0.55f >= 2.0f; }()) { dbPhase = 2; dbT0 = SimTime(); }
                        } else if (dbPhase == 2 && ((t > 8.0 && Fluid::Asleep()) || t > 40.0)) {
                            // The floor as the ground has it (the box's -2.5 m is rounded by the lattice), the surface over it.
                            for (float h = -1.0f; h > -3.0f; h -= 0.02f) if (p->densityAt(pt(-10.0f, 0.0f, h)) >= ISO_LEVEL) { dbFloor = h; break; }
                            const float s = Fluid::SurfaceAbove(*p, pt(-10.0f, 0.0f, dbFloor + 0.55f), 4.0f); // over a cell that is water at either size
                            h0 = s > 0.0f ? s + 0.55f : 0.0f;
                            const float beyond = Fluid::WaterAt(*p, pt(10.0f, 0.0f, dbFloor + 0.5f), 4.0f);
                            std::vector<TerrainEdit> cut = { box(-0.2f, 2.2f) }; // the wall goes
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            printf("bench: dam-break: reservoir %.2f m deep over its floor at %.2f m (%.0f m^3 poured, %.0f in the simulation, %.1f m^3 already past the wall); wall removed\n",
                                   h0, dbFloor, dbPoured, Fluid::GetStats().volume, beyond);
                            printf("bench: dam-break: depth along the reservoir (x = -18 -14 -10 -6 -2 m):");
                            for (float x : { -18.0f, -14.0f, -10.0f, -6.0f, -2.0f }) { const float d = Fluid::SurfaceAbove(*p, pt(x, 0.0f, dbFloor + 0.55f), 4.0f); printf(" %.2f", d > 0.0f ? d + 0.55f : 0.0f); }
                            printf("\n");
                            dbPhase = 3; dbT0 = SimTime();
                        } else if (dbPhase == 3) {
                            const float dist[3] = { 5.0f, 10.0f, 20.0f };
                            for (int i = 0; i < 3; ++i)
                                if (front[i] < 0.0 && Fluid::WaterAt(*p, pt(2.0f + dist[i], 0.0f, dbFloor + 0.5f), 0.5f) > 0.02f) front[i] = t; // over the floor cell at either size
                            if (t > 1.4 && depthDam < 0.0f) { const float sd = Fluid::SurfaceAbove(*p, pt(1.0f, 0.0f, dbFloor + 0.55f), 4.0f); depthDam = sd > 0.0f ? sd + 0.55f : 0.0f; }
                            // The discharge just below the dam (body-local, along the channel): asked every frame (the meter smooths
                            // over ~0.5 s), its peak between 1 and 2 s.
                            const float qNow = Vector3DotProduct(Fluid::FlowAt(*p, pt(2.5f, 0.0f, -1.8f), 2.0f), ax);
                            if (t > 1.0 && t < 2.0) qDam = std::max(qDam, qNow);
                            // The same discharge from the volume past x = 2.5 m (cubes +-2 m along the channel), between 1 and 2 s:
                            // independent of the flow meter.
                            static float vol1 = -1.0f, vol2 = -1.0f;
                            auto downstream = [&]() { float v = 0.0f; for (float x = 4.5f; x < 50.0f; x += 4.0f) v += Fluid::WaterAt(*p, pt(x, 0.0f, dbFloor + 2.0f), 2.0f); return v; };
                            if (t > 1.0 && vol1 < 0.0f) vol1 = downstream();
                            if (t > 2.0 && vol2 < 0.0f) vol2 = downstream();
                            if (t > 3.0) {
                                const float g = 9.81f, c = 2.0f * sqrtf(g * std::max(h0, 0.01f));
                                const float dTheory = 4.0f / 9.0f * h0, qTheory = 8.0f / 27.0f * sqrtf(g) * powf(h0, 1.5f) * 4.0f;
                                const float speed = front[1] > 0.0 ? 10.0f / (float)front[1] : 0.0f;
                                printf("bench: dam-break: front at 5 m after %.2f s, 10 m after %.2f s, 20 m after %.2f s (theory %.2f / %.2f / %.2f s: %.1f m/s); mean speed to 10 m %.1f m/s = %.2f of theory\n",
                                       front[0], front[1], front[2], 5.0f / c, 10.0f / c, 20.0f / c, c, speed, speed / c);
                                printf("bench: dam-break: depth at the dam %.2f m (theory %.2f: %.2f); discharge %.1f m^3/s (theory %.1f: %.2f); by volume past 2.5 m (1-2 s) %.1f m^3/s (%.2f)\n",
                                       depthDam, dTheory, depthDam / std::max(dTheory, 0.01f), qDam, qTheory, qDam / std::max(qTheory, 0.01f), vol2 - vol1, (vol2 - vol1) / std::max(qTheory, 0.01f));
                                const bool ok = h0 > 1.5f && front[1] > 0.0 && speed / c > 0.4f && speed / c < 1.3f;
                                printf("bench: dam-break: -> %s\n", ok ? "PASS" : "FAIL");
                                benchDamBreak = false;
                            }
                        }
                    }
                } else if (benchFlowShade && (fsPhase > 0 || WorldIdle())) {
                    // --test-flow-shade (#677, #696; best with --sandbox): a pool fed by a source spills down an 18 m chute (2 m drop)
                    // and falls 4 m into a pit - still, running and falling water side by side, watched from the pit's far side against
                    // the sun so the reflection, the streaks and the foam all show. PASS: the fall's steep water is drawn with
                    // its current running down its surface at 1 m/s or more (it pointed into the surface: nothing was drawn).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - fsT0;
                    auto at = [&](float x, float y, float h) { return Vector3Add(fsMid, Vector3Add(Vector3Add(Vector3Scale(fsE, x), Vector3Scale(fsN, y)), Vector3Scale(Vector3Normalize(fsMid), h))); };
                    if (!pf) { printf("bench: flow-shade: needs the planet -> FAIL\n"); benchFlowShade = false; }
                    else {
                        if (fsPhase == 0) {
                            const Vector3 home = pf->findSpawnDir(p->seed);
                            fsE = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })); fsN = Vector3CrossProduct(fsE, home);
                            fsMid = Vector3Scale(home, pf->radius + pf->elevation(home));
                            const Vector3 up = Vector3Normalize(fsMid);
                            auto box = [&](Vector3 c, Vector3 half, Vector3 ax, Vector3 ay, Vector3 az) {
                                TerrainEdit e{}; e.center = c; e.carve = true; e.box = true; e.half = half; e.ax = ax; e.ay = ay; e.az = az; e.radius = Vector3Length(half);
                                return e;
                            };
                            std::vector<TerrainEdit> cut;
                            cut.push_back(box(at(-22.0f, 0.0f, -1.0f), { 4.5f, 4.0f, 1.0f }, fsE, fsN, up)); // the pool: 9 x 8 m, 2 m deep
                            // The chute: its bed from -1 m at the pool down to -3 m at the pit, 3 m wide (a tilted box, 2 m tall over the bed).
                            const Vector3 ax = Vector3Normalize(Vector3Subtract(Vector3Scale(fsE, 18.0f), Vector3Scale(up, 2.0f)));
                            const Vector3 az = Vector3Normalize(Vector3Add(Vector3Scale(up, 18.0f), Vector3Scale(fsE, 2.0f)));
                            cut.push_back(box(Vector3Add(at(-9.0f, 0.0f, -2.0f), az), { 9.6f, 1.5f, 1.0f }, ax, fsN, az));
                            cut.push_back(box(at(8.0f, 0.0f, -3.5f), { 8.0f, 8.0f, 3.5f }, fsE, fsN, up)); // the pit: 16 x 16 m, 7 m deep
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            sunOverride = true; // low, behind the water seen from the camera: a glint and the sky in it
                            sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(Vector3Scale(up, 0.45f), Vector3Add(Vector3Scale(fsE, -1.0f), Vector3Scale(fsN, 0.25f)))));
                            player.flyMode = true;
                            printf("bench: flow-shade: a 9 m pool, an 18 m chute (2 m drop) and a 4 m fall into a 16 m pit\n");
                            fsPhase = 1; fsT0 = GetTime();
                        }
                        // Over the pit's far corner, looking back up the chute to the pool: the fall below, the running chute in the middle.
                        world.camera.position = p->toWorld(at(10.0f, -9.0f, 11.0f)); world.camera.target = p->toWorld(at(-9.0f, 0.0f, -2.5f));
                        world.camera.up = p->dirToWorld(Vector3Normalize(fsMid));
                        player.position = world.renderOrigin + world.camera.position;
                        if (fsPhase == 1 && t > 1.0 && WorldIdle()) {
                            Fluid::AddSource(*p, at(-22.0f, 0.0f, -1.9f), 0, 10.0f, 1.4f); // the pool fills to 0.5 m under the rim and spills down the chute
                            fsPhase = 2; fsT0 = GetTime();
                        } else if (fsPhase == 2 && t > 60.0) {
                            long long fv[3]; ChunkWaterFlowStats(fv);
                            const double ok = fv[1] > 0 ? (double)fv[2] / (double)fv[1] : 0.0;
                            const Fluid::Stats s = Fluid::GetStats();
                            printf("bench: flow-shade: %.0f m^3 of water; moving water vertices meshed %lld, steep %lld, steep drawn running down their surface at >= 1 m/s %lld (%.0f%%) -> %s\n",
                                   s.volume, fv[0], fv[1], fv[2], ok * 100.0, fv[1] >= 30 && ok >= 0.6 ? "PASS" : "FAIL");
                            TakeScreenshot("bench_flow_shade.png");
                            benchFlowShade = false;
                        }
                    }
                } else if (benchFluidLoad > 0.0f && (flPhase > 0 || WorldIdle())) {
                    // --test-fluid-load S (Plan 10, #642 perf; best with --sandbox): an S x S m basin 3 m deep cut into the ground
                    // by the start and filled from a source in every 16 m tile (8 m^3/s each, head 2.5 m), so every brick of it is
                    // awake while it fills. Between 10 and 30 s: the awake bricks, the simulation job's time, the frame times.
                    // The fluid's GPU and main-thread share prints with the bench's per-pass timings ("fluid").
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - flT0;
                    static std::vector<float> frameMs;
                    static double jobMs = 0.0, meshMs = 0.0, jobMax = 0.0; static long long samples = 0; static int bricksMax = 0, activeMax = 0;
                    auto tile = [&](float x, float y, float depthBelow) { return Vector3Add(flMid, Vector3Add(Vector3Add(Vector3Scale(flE, x), Vector3Scale(flN, y)), Vector3Scale(Vector3Normalize(flMid), -depthBelow))); };
                    if (!pf) { printf("bench: fluid-load: needs the planet -> FAIL\n"); benchFluidLoad = 0.0f; }
                    else if (flPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        flE = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })); flN = Vector3CrossProduct(flE, home);
                        flMid = Vector3Scale(home, pf->radius + pf->elevation(home));
                        const int tiles = std::max(1, (int)lroundf(benchFluidLoad / 16.0f));
                        std::vector<TerrainEdit> cut;
                        for (int i = 0; i < tiles; ++i)
                            for (int j = 0; j < tiles; ++j) {
                                TerrainEdit e{};
                                e.center = tile(16.0f * (i - 0.5f * (tiles - 1)), 16.0f * (j - 0.5f * (tiles - 1)), 1.5f);
                                e.carve = true; e.box = true; e.half = { 8.0f, 8.0f, 1.5f };
                                e.ax = flE; e.ay = flN; e.az = Vector3Normalize(flMid);
                                e.radius = Vector3Length(e.half);
                                cut.push_back(e);
                            }
                        ApplyTerrainEdits(*p, cut.data(), cut.size());
                        const Vector3 eye = Vector3Add(tile(0.0f, -0.3f * tiles * 16.0f, 0.0f), Vector3Scale(Vector3Normalize(flMid), 12.0f + 0.9f * tiles * 16.0f));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(flMid); world.camera.up = p->dirToWorld(Vector3Normalize(flMid));
                        player.flyMode = true; player.position = world.renderOrigin + world.camera.position;
                        printf("bench: fluid-load: a %d x %d m basin 3 m deep (%d tiles of 16 m)\n", tiles * 16, tiles * 16, tiles * tiles);
                        flPhase = 1; flT0 = GetTime();
                    } else {
                        world.camera.position = p->toWorld(Vector3Add(tile(0.0f, -0.3f * benchFluidLoad, 0.0f), Vector3Scale(Vector3Normalize(flMid), 12.0f + 0.9f * benchFluidLoad)));
                        world.camera.target = p->toWorld(flMid); world.camera.up = p->dirToWorld(Vector3Normalize(flMid));
                        player.position = world.renderOrigin + world.camera.position;
                        if (flPhase == 1 && t > 1.0 && WorldIdle()) {
                            const int tiles = std::max(1, (int)lroundf(benchFluidLoad / 16.0f));
                            for (int i = 0; i < tiles; ++i)
                                for (int j = 0; j < tiles; ++j)
                                    Fluid::AddSource(*p, tile(16.0f * (i - 0.5f * (tiles - 1)), 16.0f * (j - 0.5f * (tiles - 1)), 2.7f), 0, 8.0f, 2.5f);
                            flPhase = 2; flT0 = GetTime();
                        } else if (flPhase == 2) {
                            if (t > 10.0 && t <= 30.0) {
                                const Fluid::Stats s = Fluid::GetStats();
                                frameMs.push_back(GetFrameTime() * 1000.0f);
                                jobMs += s.stepMs; meshMs += s.meshMs; jobMax = std::max(jobMax, s.maxStepMs); ++samples;
                                bricksMax = std::max(bricksMax, s.bricks); activeMax = std::max(activeMax, s.activeBricks);
                            } else if (t > 30.0) {
                                std::sort(frameMs.begin(), frameMs.end());
                                auto pct = [&](float q) { return frameMs.empty() ? 0.0f : frameMs[std::min(frameMs.size() - 1, (size_t)(q * frameMs.size()))]; };
                                const Fluid::Stats s = Fluid::GetStats();
                                printf("bench: fluid-load: %d bricks (%d awake at most), %.0f m^3; step job %.2f ms mean (max %.2f), publish %.2f ms mean; frame ms p50 %.2f p95 %.2f p99 %.2f over %zu frames; GPU slots %d of %d -> PASS\n",
                                       bricksMax, activeMax, s.volume, samples ? jobMs / samples : 0.0, jobMax, samples ? meshMs / samples : 0.0, pct(0.5f), pct(0.95f), pct(0.99f), frameMs.size(), s.gpuSlots, s.gpuCapacity);
                                benchFluidLoad = 0.0f;
                            }
                        }
                    }
                } else if (benchWaterSpeed && (wsPhase > 0 || WorldIdle())) {
                    // --test-water-speed (Plan 10 W0, #638): a dam break. A 3 x 3 m column of water 6 m tall poured on flat dry
                    // ground near the start; how fast its front spreads (5, 10, 15 m from the middle), how far it gets, how long
                    // until the water sleeps, and the volume kept. A shallow-water front from a 6 m column runs ~2 sqrt(g h)
                    // (~15 m/s); the cellular sim today spreads by diffusion. Measurement only: PASS = it ran (W4 sets the bar).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const double t = GetTime() - wsT0;
                    if (!pf) { printf("bench: water-speed: needs the planet\n"); benchWaterSpeed = false; }
                    else if (wsPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        auto at = [&](Vector3 d, float e, float n) { return Vector3Normalize(Vector3Add(d, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius)))); };
                        Vector3 site{}; float bestRange = 1e30f;
                        for (float r = 0.0f; r <= 600.0f && bestRange > 0.6f; r += 20.0f)
                            for (float a = 0.0f; a < 2.0f * PI; a += r > 0.0f ? 20.0f / r : 7.0f) {
                                const Vector3 d = at(home, r * cosf(a), r * sinf(a));
                                const float e0 = pf->elevation(d);
                                if (e0 < 5.0f || pf->inlandWaterNear(d, 40.0f)) continue;
                                float lo = e0, hi = e0;
                                for (float rr = 5.0f; rr <= 25.0f; rr += 5.0f)
                                    for (int k = 0; k < 8; ++k) { const float e = pf->elevation(at(d, rr * cosf(k * PI / 4), rr * sinf(k * PI / 4))); lo = std::min(lo, e); hi = std::max(hi, e); }
                                if (hi - lo < bestRange) { bestRange = hi - lo; site = d; }
                            }
                        if (bestRange > 2.0f) { printf("bench: water-speed: no flat dry ground within 600 m -> FAIL\n"); benchWaterSpeed = false; }
                        else {
                            wsSite = Vector3Scale(site, pf->radius + pf->elevation(site));
                            // Watch from nearly overhead, 28 m up (over the trees), looking down at the column; pour once the view has
                            // loaded (phase 3), so a capture sees the whole spread.
                            const Vector3 eyeDir = at(site, -6.0f, 0.0f);
                            player.flyMode = true;
                            const Vector3 eye = Vector3Scale(eyeDir, pf->radius + pf->elevation(site) + 28.0f);
                            wsEye = eye;
                            player.placeAt(world, *p, p->toWorld(eye));
                            player.forward = p->dirToWorld(Vector3Normalize(Vector3Subtract(wsSite, eye)));
                            player.applyToCamera(world);
                            printf("bench: water-speed: site %.0f m from the start, ground within %.1f m over 25 m\n", Vector3Distance(site, home) * pf->radius, bestRange);
                            wsPhase = 3; wsT0 = GetTime();
                        }
                    } else if (wsPhase == 3) {
                        world.camera.position = p->toWorld(wsEye); world.camera.target = p->toWorld(wsSite);
                        world.camera.up = p->dirToWorld(Vector3Normalize(wsSite));
                        player.position = world.renderOrigin + world.camera.position;
                        if (t > 1.0 && WorldIdle()) {
                            const Vector3 site = Vector3Normalize(wsSite);
                            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(site, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, site);
                            for (int x = -1; x <= 1; ++x)
                                for (int y = -1; y <= 1; ++y) { // 6 m^3 per pour: a cell column filled 6 cells up
                                    const Vector3 d = Vector3Normalize(Vector3Add(site, Vector3Add(Vector3Scale(t1, x / pf->radius), Vector3Scale(t2, y / pf->radius))));
                                    Fluid::Pour(*p, Vector3Add(Vector3Scale(d, pf->radius + pf->elevation(d)), Vector3Scale(site, 0.3f)), 6.0f);
                                }
                            printf("bench: water-speed: poured 54 m^3 as a 3 x 3 x 6 m column after %.1f s of loading\n", t);
                            wsPhase = 1; wsT0 = GetTime();
                        }
                    } else if (wsPhase == 1) {
                        // Hold the camera (the player drifts down otherwise; captures watch the spread). Set directly, as --test-dam
                        // does: placing the player every frame kept the loading screen up.
                        world.camera.position = p->toWorld(wsEye); world.camera.target = p->toWorld(wsSite);
                        world.camera.up = p->dirToWorld(Vector3Normalize(wsSite));
                        player.position = world.renderOrigin + world.camera.position;
                        // The front: the farthest of 8 rays where water stands within 0.6 m of the ground (sampled every 4 frames).
                        static int frame = 0;
                        if ((++frame & 3) == 0) {
                            const Vector3 up = Vector3Normalize(wsSite);
                            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, up);
                            float front = 0.0f;
                            for (int k = 0; k < 8; ++k)
                                for (float s = 1.0f; s <= 30.0f; s += 1.0f) {
                                    const Vector3 q = Vector3Normalize(Vector3Add(up, Vector3Add(Vector3Scale(t1, s * cosf(k * PI / 4) / pf->radius), Vector3Scale(t2, s * sinf(k * PI / 4) / pf->radius))));
                                    if (Fluid::WaterAt(*p, Vector3Scale(q, pf->radius + pf->elevation(q) + 0.2f), 0.6f) > 0.005f) front = std::max(front, s);
                                }
                            wsMaxFront = std::max(wsMaxFront, front);
                            const float edges[3] = { 5.0f, 10.0f, 15.0f };
                            for (int i = 0; i < 3; ++i) if (wsFront[i] < 0.0 && front >= edges[i]) wsFront[i] = t;
                            const Fluid::Stats fs = Fluid::GetStats();
                            if (wsVol0 <= 0.0f && fs.volume > 1.0) wsVol0 = (float)fs.volume;
                            if (t > 1.0 && wsSleep < 0.0 && Fluid::Asleep()) wsSleep = t;
                            if (wsSleep >= 0.0 || t > 120.0) {
                                printf("bench: water-speed: front at 5 m after %.2f s, 10 m after %.2f s, 15 m after %.2f s (-1: never); farthest %.0f m; asleep after %.1f s (-1: not in 120 s)\n",
                                       wsFront[0], wsFront[1], wsFront[2], wsMaxFront, wsSleep);
                                printf("bench: water-speed: volume %.1f m^3 at the start, %.1f m^3 at the end; step %.2f ms (max %.2f), %d bricks (%d active)\n",
                                       wsVol0, fs.volume, fs.stepMs, fs.maxStepMs, fs.bricks, fs.activeBricks);
                                printf("bench: water-speed: %s\n", wsFront[0] >= 0.0 ? "PASS" : "FAIL");
                                benchWaterSpeed = false;
                            }
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
                        TerrainTool_Update(world, PlayerInput::None(), true); // aim only
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
                        for (int pass = 0; pass < 2 && best > 1e29f; ++pass) // 6 km at 24 m; then 20 km at 60 m (the bake's small lakes are 200 m+, rarer, #360)
                            for (float e = pass ? -20000.0f : -6000.0f; e <= (pass ? 20000.0f : 6000.0f); e += pass ? 60.0f : 24.0f)
                            for (float n = pass ? -20000.0f : -6000.0f; n <= (pass ? 20000.0f : 6000.0f); n += pass ? 60.0f : 24.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if (ws.kind != 3) continue;
                                best = dist; pondDirB = ws.pondDir; pondLevelB = ws.level; pondRB = ws.pondR;
                            }
                        if (best > 1e29f) {
                            // Eroded worlds (#487): the bake's ponds are few (309 on seed 42) and may all be far from the
                            // start; the nearest of them, wherever it is (the test travels there).
                            // (One inside a lake reads as the lake: the next one then.)
                            int skipped = 0;
                            for (const Vector3& bd : pf->bakePondsByDistance(home)) {
                                const PlanetField::WaterSource ws = pf->waterSourceAt(bd);
                                if (ws.kind != 3) { ++skipped; continue; }
                                best = acosf(std::clamp(Vector3DotProduct(home, bd), -1.0f, 1.0f)) * pf->radius;
                                pondDirB = ws.pondDir; pondLevelB = ws.level; pondRB = ws.pondR;
                                printf("bench: pond: none within 20 km of the start; the bake's nearest pond is %.1f km away (%d nearer ones lie in lakes)\n", best / 1000.0f, skipped);
                                break;
                            }
                        }
                        if (best > 1e29f) { printf("bench: pond: no pond within 6 km -> FAIL\n"); benchPond = false; }
                        else {
                            pondOutB = Vector3Normalize(Vector3CrossProduct(pondDirB, { 0, 1, 0 }));
                            const float bed = pf->elevation(pondDirB);
                            printf("bench: pond: %.0f m from the start, radius %.0f m, level %.1f m, %.1f m deep at the middle\n", best, pondRB, pondLevelB, pondLevelB - bed);
                            printf("bench: pond: middle column: %s\n", pf->debugColumn(pondDirB).c_str());
                            // A trench from the middle out past the rim, its floor below the pond's bed, then a pit.
                            std::vector<TerrainEdit> cut;
                            const float floorH = bed - 1.5f;
                            for (float s = 0.0f; s <= pondRB + 12.0f; s += 0.8f)
                                for (float h = floorH + 1.2f; h <= pondLevelB + 3.0f; h += 1.2f) cut.push_back({ at(along(s), h), 1.4f, -1.0f, true });
                            for (int k = 0; k < 3; ++k) cut.push_back({ at(along(pondRB + 14.0f), floorH - 1.5f - 1.5f * k), 3.5f, -1.0f, true });
                            {   // The pit holds at least the pond's top 0.5 m (#360: the bake's ponds reach 60 m across; three
                                // spheres held ~500 m^3, a 60 m pond's top 0.3 m is ~3400): more spheres around it.
                                const float need = PI * pondRB * pondRB * 0.5f, each = 4.0f / 3.0f * PI * 3.5f * 3.5f * 3.5f * 0.6f;
                                const int extra = std::max(0, (int)ceilf(need / each) - 3);
                                const Vector3 side = Vector3Normalize(Vector3CrossProduct(pondDirB, pondOutB));
                                for (int k = 0; k < extra; ++k) {
                                    const float ring = 5.0f * (1 + k / 8), ang = (k % 8) * (PI / 4.0f);
                                    const Vector3 d = Vector3Normalize(Vector3Add(along(pondRB + 14.0f + ring * cosf(ang)), Vector3Scale(side, ring * sinf(ang) / pf->radius)));
                                    cut.push_back({ at(d, floorH - 0.5f), 3.5f, -1.0f, true }); // shallower than the core: it fills first (swimming)
                                }
                            }
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
                            printf("bench: pond t=%3.0f s: %d bricks (%d awake), %.0f wet cells, %.1f m^3, %d bricks of pond water; middle %s, pit %s; %d settled, woken gpu %lld cpu %lld; changes (total) d%lld s%lld r%lld l%lld; GPU slots %d of %d (%lld released)\n", t, s.bricks, s.activeBricks, s.wetCells, s.volume, s.basinBricks,
                                   Fluid::Wet(*p, pondMid) ? "wet" : "dry", Fluid::Wet(*p, pondPit) ? "wet" : "dry", s.settled, s.wakes[0], s.wakes[1], s.dbgTotal[0], s.dbgTotal[1], s.dbgTotal[2], s.dbgTotal[3],
                                   s.gpuSlots, s.gpuCapacity, s.slotReleases);
                        }
                        if (t > 8.0 && t < 8.1) TakeScreenshot("bench_pond.png");
                        // The drain into the pit runs ~3 min on seed 42's 35 m pond (levels meet ~180 s); the water sleeps
                        // within ~10 s of that (#366). A 150 s cap judged it mid-drain, which read as "never settles".
                        if ((t > 15.0 && Fluid::Asleep()) || t > 240.0) {
                            if (!Fluid::Asleep()) printf("bench: pond: not asleep after 240 s (%d of %d bricks awake)\n", Fluid::GetStats().activeBricks, Fluid::GetStats().bricks);
                            TakeScreenshot("bench_pond_drained.png");
                            const bool claimed = Fluid::GetStats().basinBricks > 0; // the pond is simulated water (#360: no take-over)
                            const bool midDry = !Fluid::Wet(*p, pondMid), pitWet = Fluid::Wet(*p, pondPit);
                            const float vol = (float)Fluid::GetStats().volume;
                            printf("bench: pond: settled after %.0f s: taken over %d, the top 0.3 m of the pond %s, the pit %s, %.1f m^3 of water, %lld pond faces pinned at the site edge\n", t, (int)claimed, midDry ? "drained" : "still full", pitWet ? "filled" : "dry", vol, Fluid::GetStats().pinnedFaces);
                            // #385: farming, animals and the map ask GeneratedWaterValid 0.2 m under the generated level. The
                            // pond only drops ~0.35 m (its middle still holds water, and they rightly say so): judge a shore
                            // column the drain left dry, its bed 0.05-0.25 m under the old level.
                            bool readersDry = true;
                            {
                                float shore = -1.0f;
                                for (float s = 0.0f; s <= pondRB + 4.0f && shore < 0.0f; s += 0.5f) {
                                    const float bed = pf->elevation(along(-s));
                                    if (bed > pondLevelB - 0.25f && bed < pondLevelB - 0.05f && !Fluid::Wet(*p, at(along(-s), bed + 0.05f))) shore = s;
                                }
                                if (!Fluid::Asleep()) printf("bench: pond: still draining at the cap: the readers' check is skipped (a draining shore still holds films)\n");
                                else if (shore >= 0.0f) {
                                    readersDry = !Fluid::GeneratedWaterValid(p->id, at(along(-shore), pondLevelB - 0.2f));
                                    printf("bench: pond: farming, animals and the map read the drained shore (%.1f m out) as %s -> %s\n", shore, readersDry ? "dry" : "water", readersDry ? "PASS" : "FAIL");
                                } else printf("bench: pond: no drained shore column found to check the readers on\n");
                            }
                            {   // #189: the taken-over pond keeps its bed: no trees grow in it.
                                std::vector<Trunk> trunks;
                                if (p->trees) p->trees->trunksNear(at(pondDirB, pondLevelB), pondRB * 0.8f, trunks);
                                // In the bed: within 0.8 of the radius across, standing no higher than the water (trunksNear
                                // also returned trees on the bank 6 m above it, 35-38 m from the middle of a 34 m pond).
                                int inBed = 0;
                                for (const Trunk& tr : trunks) {
                                    const Vector3 db = Vector3Normalize(tr.base);
                                    inBed += Vector3Distance(db, pondDirB) * pf->radius < pondRB * 0.8f && Vector3Length(tr.base) - pf->radius < pondLevelB + 0.3f;
                                }
                                printf("bench: pond: %d trees inside the taken-over pond (%zu near it, on its banks) -> %s\n", inBed, trunks.size() - inBed, inBed == 0 ? "PASS" : "FAIL");
                            }
                            const bool saved = SaveWorld(world);
                            // Away 300 m, load, back.
                            player.placeAt(world, *p, p->toWorld(at(along(pondRB + 300.0f), pf->elevation(along(pondRB + 300.0f)) + 2.0f)));
                            const bool loaded = saved && LoadWorldSave(world);
                            printf("bench: pond: saved %d, went 300 m away, loaded %d\n", (int)saved, (int)loaded);
                            pondPhase = 2; pondT0 = GetTime();
                            pondRB = claimed && midDry && pitWet && readersDry ? pondRB : -pondRB; // remember the first result
                        }
                    } else if (pondPhase == 2) {
                        const Vector3 up = pondDirB;
                        const float r = fabsf(pondRB);
                        const Vector3 eye = Vector3Add(at(along(r + 6.0f), pondLevelB + 18.0f), Vector3Scale(Vector3CrossProduct(up, pondOutB), 22.0f));
                        world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(at(along(r * 0.5f), pondLevelB)); world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        // Judged once asleep, or after 120 s regardless (#360: the bake's 60 m ponds keep a few bricks trickling;
                        // that is reported, #366): the state is what this checks.
                        if ((t > 12.0 && Fluid::Asleep()) || t > 120.0) {
                            if (!Fluid::Asleep()) printf("bench: pond: not asleep after 120 s (%d of %d bricks awake): judged as it stands\n", Fluid::GetStats().activeBricks, Fluid::GetStats().bricks);
                            TakeScreenshot("bench_pond_loaded.png");
                            const PlanetField::WaterSource ws = pf->waterSourceAt(pondDirB);
                            const bool stillClaimed = Fluid::Simulated(*p, pondMid); (void)ws; // its cells still decide (drained), not the generated pond
                            const bool midDry = !Fluid::Wet(*p, pondMid), pitWet = Fluid::Wet(*p, pondPit);
                            const bool pass = pondRB > 0.0f && stillClaimed && midDry && pitWet;
                            printf("bench: pond: after the load the pond is %s, its top %s, the pit %s, %.1f m^3 -> %s\n", stillClaimed ? "still taken over" : "back to static water",
                                   midDry ? "drained" : "full", pitWet ? "full" : "dry", Fluid::GetStats().volume, pass ? "PASS" : "FAIL");
                            printf("bench: pond: bench_pond.png, bench_pond_drained.png, bench_pond_loaded.png\n");
                            pondPhase = 4; pondT0 = GetTime(); // swimming next (topped up first when shallow)
                        } else if (t > 120.0) { printf("bench: pond: did not settle after the load -> FAIL\n"); benchPond = false; }
                    } else if (pondPhase == 4) {
                        // Swimming: into the flooded pit, 1.5 m under its surface. Water shallower than a player (the drained
                        // volume depends on the cell size: quarter-metre cells follow the ground and hold less, W4c) is topped up
                        // by pouring into the pit first (under 1.2 m: up to 1.5 m), so the swim is judged in water deep enough to swim in.
                        const float surf = Fluid::SurfaceAbove(*p, pondPit, 12.0f);
                        static int topping = -1; // decided on the first frame here: under 1.2 m it is poured up to 1.5 m
                        if (topping < 0) topping = surf < 1.2f ? 1 : 0;
                        static double lastPour = -10.0; // poured water mounds where it lands: judged 2 s after a pour, once spread
                        if (topping && t < 40.0 && (t - lastPour < 2.0 || surf < 1.6f)) {
                            if (t - lastPour >= 2.0) { lastPour = t; Fluid::Pour(*p, Vector3Add(pondPit, Vector3Scale(Vector3Normalize(pondPit), 0.3f)), 8.0f); }
                        } else {
                            player.flyMode = false;
                            const float lift = surf > 0.0f ? std::clamp(surf - 1.5f, 0.3f, 1.0f) : 1.0f;
                            printf("bench: pond: the pit's water stands %.1f m over its floor point%s; swimming from %.1f m up\n", surf, topping ? " (topped up)" : "", lift);
                            player.placeAt(world, *p, p->toWorld(Vector3Add(pondPit, Vector3Scale(Vector3Normalize(pondPit), lift))));
                            pondPhase = 3; pondT0 = GetTime();
                        }
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
                } else if (benchPressure && (prPhase > 0 || WorldIdle())) {
                    // --test-pressure (#393): a 12 m shaft by the start filled by a spring: 5 m under its surface reads
                    // ~49 kPa (rho g h); 2 m under the nearest pond's generated surface reads ~19.6 kPa.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    auto at = [&](Vector3 dir, float h) { return Vector3Scale(dir, pf->radius + h); };
                    if (!pf) { printf("bench: pressure: needs the planet -> FAIL\n"); benchPressure = false; }
                    else if (prPhase == 0) {
                        prDir = StartSpawnDir(*p);
                        prGround = pf->elevation(prDir) + 8.0f;
                        while (prGround > pf->elevation(prDir) - 8.0f && p->densityAt(at(prDir, prGround)) < ISO_LEVEL) prGround -= 0.1f;
                        std::vector<TerrainEdit> cut;
                        for (float h = prGround - 0.5f; h >= prGround - 12.0f; h -= 0.8f) cut.push_back({ at(prDir, h), 1.6f, -1.0f, true });
                        ApplyTerrainEdits(*p, cut.data(), cut.size());
                        Fluid::AddSource(*p, at(prDir, prGround - 12.5f), 0, 3.0f, 11.0f); // fills to ~1.5 m under the rim
                        printf("bench: pressure: a 12 m shaft dug, a spring at its bottom\n");
                        prPhase = 1; prT0 = GetTime();
                    } else if (GetTime() - prT0 > 50.0) {
                        // The shaft's surface: the top of the wet run up the middle.
                        float surf = -1e9f;
                        for (float h = prGround - 12.0f; h <= prGround + 1.0f; h += 0.05f) { if (Fluid::Wet(*p, at(prDir, h))) surf = h; else if (surf > -1e8f) break; }
                        const float shaftP = surf > -1e8f ? Fluid::PressureAt(*p, at(prDir, surf - 5.0f)) : 0.0f;
                        // The nearest pond (generated water): 2 m under its level at its middle.
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(prDir, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, prDir);
                        float best = 1e30f, pondP = -1.0f, pondDepth = 0.0f;
                        for (float e = -6000.0f; e <= 6000.0f; e += 48.0f)
                            for (float n = -6000.0f; n <= 6000.0f; n += 48.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(prDir, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if (ws.kind != 3 || ws.level - pf->elevation(ws.pondDir) < 2.5f) continue;
                                best = dist; pondDepth = ws.level - pf->elevation(ws.pondDir);
                                pondP = Fluid::PressureAt(*p, at(ws.pondDir, ws.level - 2.0f));
                            }
                        const float wantShaft = 5.0f * 9810.0f, wantPond = 2.0f * 9810.0f;
                        const bool okShaft = fabsf(shaftP - wantShaft) < 0.1f * wantShaft, okPond = pondP < 0.0f || fabsf(pondP - wantPond) < 0.1f * wantPond;
                        printf("bench: pressure: the shaft's water stands %.1f m over its bottom; 5 m under its surface %.1f kPa (want %.1f)\n", surf - (prGround - 12.0f), shaftP / 1000.0f, wantShaft / 1000.0f);
                        if (pondP >= 0.0f) printf("bench: pressure: a pond %.0f m away, %.1f m deep: 2 m under its surface %.1f kPa (want %.1f)\n", best, pondDepth, pondP / 1000.0f, wantPond / 1000.0f);
                        else printf("bench: pressure: no pond 2.5 m+ deep within 6 km (generated water not checked)\n");
                        printf("bench: pressure: -> %s\n", okShaft && okPond ? "PASS" : "FAIL");
                        benchPressure = false;
                    }
                } else if (benchFlowMeter && (fmPhase > 0 || WorldIdle())) {
                    // --test-flow-meter (#392): a spring of 1 m^3/s at the head of a 20 m trench that ends in a deep pit; once
                    // the stream runs steady, Fluid::FlowAt across the middle of the trench reads ~1 m^3/s along it.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { printf("bench: flow-meter: needs the planet -> FAIL\n"); benchFlowMeter = false; }
                    else {
                        auto at = [&](float x, float y) { return fmPlot->toBody({ x, y, 0.0f }); };
                        const double t = GetTime() - fmT0;
                        if (fmPhase == 0) {
                            const Vector3 d = StartSpawnDir(*p);
                            float gh = pf->elevation(d) + 8.0f;
                            while (gh > pf->elevation(d) - 8.0f && p->densityAt(Vector3Scale(d, pf->radius + gh)) < ISO_LEVEL) gh -= 0.1f;
                            if (!p->plots) p->plots = std::make_unique<PlotSet>();
                            fmPlot = p->plots->plotFor(Vector3Scale(d, pf->radius + gh), 0.0f, true); // a frame on the real ground
                            std::vector<TerrainEdit> cut;
                            for (float x = -10.0f; x <= 10.0f; x += 0.8f) { cut.push_back({ at(x, -1.0f), 1.4f, -1.0f, true }); cut.push_back({ at(x, 0.2f), 1.4f, -1.0f, true }); }
                            for (int k = 0; k < 4; ++k) cut.push_back({ at(14.0f, -3.0f - 2.5f * k), 4.0f, -1.0f, true }); // the pit: ~600 m^3
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            Fluid::AddSource(*p, at(-9.5f, -2.2f), 0, 1.0f, 1.5f);
                            sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(at(0.0f, 50.0f)));
                            player.flyMode = true;
                            printf("bench: flow-meter: trench and pit dug, a 1 m^3/s spring at the trench's head\n");
                            fmPhase = 1; fmT0 = GetTime();
                        } else {
                            const Vector3 up = Vector3Normalize(at(0.0f, 0.0f));
                            world.camera.position = p->toWorld(Vector3Add(at(0.0f, 8.0f), Vector3Scale(fmPlot->north, 8.0f))); world.camera.target = p->toWorld(at(0.0f, -1.5f));
                            world.camera.up = p->dirToWorld(up);
                            player.position = world.renderOrigin + world.camera.position;
                            const Vector3 q = Fluid::FlowAt(*p, at(0.0f, -1.8f), 2.0f); // asked every frame: the readout smooths
                            static double lastPrint = 0.0;
                            const Fluid::Stats fs = Fluid::GetStats();
                            if (t - lastPrint >= 10.0) { lastPrint = t; printf("bench: flow-meter t=%3.0f s: %.2f m^3/s along the trench (%.2f across); %.1f m^3 of water; %d bricks awake, %d steady, publish %.2f ms\n", t, Vector3DotProduct(q, fmPlot->east), Vector3Length(Vector3Subtract(q, Vector3Scale(fmPlot->east, Vector3DotProduct(q, fmPlot->east)))), fs.volume, fs.activeBricks, fs.steadyBricks, fs.meshMs); }
                            // Steady flow (#674): the trench's bricks turn steady (running, read back once a second); a pour into the
                            // trench at 45 s wakes them (fewer steady), and the meter still reads the stream at 60 s.
                            static int steadyBefore = 0, steadyAfter = 1 << 30;
                            static bool poured = false;
                            if (t > 30.0 && t < 45.0) steadyBefore = std::max(steadyBefore, fs.steadyBricks);
                            if (t > 45.0 && !poured) { poured = true; Fluid::Pour(*p, at(0.0f, -1.5f), 3.0f); }
                            if (t > 45.5 && t < 48.0) steadyAfter = std::min(steadyAfter, fs.steadyBricks);
                            if (t > 60.0) {
                                const float along = Vector3DotProduct(q, fmPlot->east);
                                TakeScreenshot("bench_flow_meter.png");
                                const bool gpu = Fluid::GpuMode();
                                const bool steadyOk = !gpu || (steadyBefore > 0 && steadyAfter < steadyBefore);
                                printf("bench: flow-meter: steady bricks %d before a pour, %d after (#674)%s\n", steadyBefore, steadyAfter, gpu ? "" : " (CPU mode: none)");
                                printf("bench: flow-meter: %.2f m^3/s along the trench from a 1 m^3/s spring -> %s\n", along, along > 0.6f && along < 1.4f && steadyOk ? "PASS" : "FAIL");
                                benchFlowMeter = false;
                            }
                        }
                    }
                } else if (benchPartWall && (pwPhase > 0 || WorldIdle())) {
                    // --test-part-wall (#391): a trench on dry ground with a wall of stone blocks across its middle; water
                    // poured into the west half stays there while the wall stands, and runs east once it is taken down.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    if (!pf) { printf("bench: part-wall: needs the planet -> FAIL\n"); benchPartWall = false; }
                    else {
                        auto at = [&](float x, float y) { return pwPlot->toBody({ x, y, 0.0f }); }; // along the trench (plot east), metres
                        const double t = GetTime() - pwT0;
                        if (pwPhase == 0) {
                            const Vector3 d = StartSpawnDir(*p);
                            // The real ground (the 3D terrain stands up to metres over the height map): first solid down the column.
                            float gh = pf->elevation(d) + 8.0f;
                            while (gh > pf->elevation(d) - 8.0f && p->densityAt(Vector3Scale(d, pf->radius + gh)) < ISO_LEVEL) gh -= 0.1f;
                            const Vector3 here = Vector3Scale(d, pf->radius + gh);
                            if (!p->plots) p->plots = std::make_unique<PlotSet>();
                            pwPlot = p->plots->plotFor(here, 0.0f, true);
                            // The trench: 16 m along the plot's east axis, ~2.8 m wide, its floor ~2.4 m down.
                            std::vector<TerrainEdit> cut;
                            for (float x = -8.0f; x <= 8.0f; x += 0.8f) { cut.push_back({ at(x, -1.0f), 1.4f, -1.0f, true }); cut.push_back({ at(x, 0.2f), 1.4f, -1.0f, true }); } // open to the sky
                            ApplyTerrainEdits(*p, cut.data(), cut.size());
                            // The wall: 1 m stone blocks across it at x = 0, from 4 m down to 1 m up, 4 m to either side.
                            for (int y = -8; y <= 2; y += 2)
                                for (int z = -8; z <= 6; z += 2) {
                                    const int idx = pwPlot->place({ PartType::Block, { -1, y, z }, 0, (uint8_t)std::max(0, BuildMaterialForItem("stone")) });
                                    if (idx >= 0) pwParts.push_back(idx);
                                }
                            for (float x : { -6.5f, -4.5f, -2.5f }) Fluid::Pour(*p, at(x, -1.5f), 6.0f); // into the west half (a pour fills one column, <= 6 m^3)
                            sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(Vector3Normalize(here), Vector3Scale(pwPlot->north, 0.4f)))); // midday
                            printf("bench: part-wall: trench dug, %zu blocks in the wall, 18 m^3 poured west of it\n", pwParts.size());
                            player.flyMode = true;
                            pwPhase = 1; pwT0 = GetTime();
                        } else {
                            // Watch from above the wall, looking along the trench.
                            const Vector3 up = Vector3Normalize(at(0.0f, 0.0f));
                            world.camera.position = p->toWorld(Vector3Add(at(5.0f, 6.0f), Vector3Scale(pwPlot->north, 6.0f))); world.camera.target = p->toWorld(at(-2.0f, -1.5f));
                            world.camera.up = p->dirToWorld(up);
                            player.position = world.renderOrigin + world.camera.position;
                            auto wetAt = [&](float x) { for (float y = -2.6f; y <= -0.8f; y += 0.2f) if (Fluid::Wet(*p, at(x, y))) return true; return false; };
                            const bool west = wetAt(-4.0f), east = wetAt(4.0f);
                            if (pwPhase == 1 && t > 25.0) {
                                pwWestWet = west; pwHeld = west && !east;
                                TakeScreenshot("bench_part_wall.png");
                                printf("bench: part-wall: after 25 s the west half is %s, the east half %s -> %s\n", west ? "wet" : "dry", east ? "wet" : "dry", pwHeld ? "held" : "NOT held");
                                {   // Water perched up the trench's walls (#454): wet points 0.3-2 m over the pool's surface in the west half.
                                    auto at3 = [&](float x, float y, float z) { return pwPlot->toBody({ x, y, z }); };
                                    float surf = -99.0f;
                                    for (float y = -2.6f; y <= 0.5f; y += 0.05f) if (Fluid::Wet(*p, at3(-4.0f, y, 0.0f))) surf = y;
                                    // Films a cell or more over the pool: water in cells whose bottom is 0.5 m+ over its surface.
                                    const float perched = surf > -99.0f ? Fluid::VolumeAbove(*p, at3(-4.0f, surf, 0.0f), 5.0f, 0.5f) : -1.0f;
                                    printf("bench: part-wall: pool surface %.2f m; %.2f m^3 of water perched 0.5 m+ over it (#454)\n", surf, perched);
                                }
                                for (int idx : pwParts) pwPlot->remove(idx);
                                pwPhase = 2; pwT0 = GetTime();
                            } else if (pwPhase == 2 && (east || t > 30.0)) {
                                TakeScreenshot("bench_part_wall_open.png");
                                printf("bench: part-wall: wall removed: the east half is %s after %.1f s -> %s\n", east ? "wet" : "dry", t,
                                       pwHeld && east ? "PASS" : "FAIL");
                                printf("bench: part-wall: bench_part_wall.png, bench_part_wall_open.png\n");
                                benchPartWall = false;
                            }
                        }
                    }
                } else if (benchDam && (damPhase > 0 || WorldIdle())) {
                    // --test-dam (#361, Plan 7 step 7 gate): build a dam across the nearest small river. Adding ground into the
                    // river releases its reach: the water backs up behind the dam, the reach below it drains, and once the pool
                    // reaches the crest it spills over and the water below comes back.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
                    const RiverNetwork* net = pf ? pf->rivers() : nullptr;
                    auto at = [&](Vector3 dir, float h) { return Vector3Scale(dir, pf->radius + h); };
                    auto off = [&](Vector3 dir, Vector3 t, float s) { return Vector3Normalize(Vector3Add(dir, Vector3Scale(t, s / pf->radius))); };
                    // The water surface over a column (m above the radius): the highest wet point from just over the bed up to
                    // 5 m over the generated level (10 cm steps); dry: -1e9.
                    auto surfaceAt = [&](Vector3 dir, float level) {
                        float top = -1e9f;
                        for (float h = pf->elevation(dir) + 0.15f; h <= level + 5.0f; h += 0.1f)
                            if (Fluid::Wet(*p, at(dir, h))) top = h; else if (top > -1e8f) break;
                        return top;
                    };
                    if (!pf || !net) { printf("bench: dam: needs a planet with rivers -> FAIL\n"); benchDam = false; }
                    else if (damPhase == 0) {
                        const Vector3 home = pf->findSpawnDir(p->seed);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
                        float best = 1e30f;
                        for (float e = -8000.0f; e <= 8000.0f; e += 32.0f)
                            for (float n = -8000.0f; n <= 8000.0f; n += 32.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                                if (ws.kind != 1 || ws.level < 2.0f) continue;
                                const RiverNetwork::Hit h = net->nearest(d, 1e4f);
                                if (!h.found || h.dist > 1.0f + h.halfWidth * 0.3f || h.halfWidth < 3.0f || h.halfWidth > 12.0f || h.depth < 1.5f) continue;
                                // In a valley: the valley sides 15 m past the channel stand 1.2 m+ over the water (eroded terrain; noise-terrain rivers are perched), so the pool can rise
                                // to the crest instead of spreading over a flat floodplain around the dam's ends.
                                const Vector3 dn = Vector3Normalize(Vector3Subtract(h.down, Vector3Scale(d, Vector3DotProduct(h.down, d))));
                                const Vector3 sd = Vector3Normalize(Vector3CrossProduct(d, dn));
                                const float bank = std::min(pf->elevation(off(d, sd, h.halfWidth + 15.0f)), pf->elevation(off(d, sd, -h.halfWidth - 15.0f)));
                                if (bank < ws.level + (benchDamDivert ? 2.5f : 1.2f)) continue; // #656: a dam that holds the river a while
                                // A gentle reach (#729): the generated level 12 m up and down the river within 0.6 m (a 0.8 m crest then
                                // stands over the pool's start), and for --test-dam-divert under 6 m of drop to the reach it watches 200 m
                                // downstream, clear of the sea (a reach the sea backs up never dries). After bake v26 seed 42's pick fell
                                // 4 m over 24 m: the pool stood over the crest from the start.
                                if (pf->waterSourceAt(off(d, dn, -12.0f)).level - pf->waterSourceAt(off(d, dn, 12.0f)).level > 0.6f) continue;
                                if (benchDamDivert) {
                                    Vector3 w = d;
                                    for (int s = 0; s < 20; ++s) { const RiverNetwork::Hit hw = net->nearest(w, 1e4f); w = off(w, hw.down, 10.0f); }
                                    const float lw = net->nearest(w, 1e4f).level;
                                    if (ws.level - lw > 6.0f || lw < 3.0f) continue;
                                }
                                best = dist; damDir = d; damLevel = ws.level; damHW = h.halfWidth; damDepth = h.depth; damCrest = benchDamDivert ? bank - 0.3f : std::min(ws.level + 0.8f, bank - 0.3f); // a low dam: a released reach fills at what its channel carries (#432)
                                damDown = dn; // along the ground
                            }
                        if (best > 1e29f) {
                            const char* e = getenv("RAYLIBGAME_EROSION");
                            if (!e || atoi(e) == 0) printf("bench: dam: no river in a valley within 8 km (noise-terrain rivers are perched over their banks; run with RAYLIBGAME_EROSION=1) -> SKIP\n");
                            else printf("bench: dam: no small river in a valley within 8 km of the start -> FAIL\n");
                            benchDam = false;
                        }
                        else {
                            damSide = Vector3Normalize(Vector3CrossProduct(damDir, damDown));
                            printf("bench: dam: river %.0f m from the start: level %.2f m, half-width %.1f m, depth %.1f m\n", best, damLevel, damHW, damDepth);
                            // First the current (#361) at the dam site, before the dam (nothing released yet): what
                            // carries a swimmer and pushes a wader (player.cpp) and what water.fs draws. Mid-channel ~1 m/s
                            // downstream, slower near the bank, none on dry land.
                            const Vector3 dd = damDir; // mid-channel by the search (a point 60 m along a straight line left the meandering river)
                            const float lv = pf->waterSourceAt(dd).level;
                            const Vector3 mid = p->dirToLocal(RiverCurrent(p, p->toWorld(at(dd, lv - 0.5f))));
                            const Vector3 bank = p->dirToLocal(RiverCurrent(p, p->toWorld(at(off(dd, damSide, damHW + 1.5f), lv - 0.3f))));
                            const Vector3 dry = p->dirToLocal(RiverCurrent(p, p->toWorld(at(off(dd, damSide, damHW + 20.0f), lv + 3.0f))));
                            const float ms = Vector3Length(mid), bs = Vector3Length(bank), ds = Vector3Length(dry);
                            const float along = ms > 0.0f ? Vector3DotProduct(mid, damDown) / ms : 0.0f;
                            damDrift = (ms > 0.6f && ms < 1.5f && along > 0.9f && bs < ms && ds == 0.0f) ? ms : 0.0f;
                            printf("bench: dam: the current at the dam site: %.2f m/s mid-channel (%.2f along the river), %.2f by the bank, %.2f on dry land -> %s\n",
                                   ms, along, bs, ds, damDrift > 0.0f ? "PASS" : "FAIL");
                            damPhase = 10;
                        }
                    } else if (damPhase == 10) {
                        {
                            // The dam: a wall across the river and 15 m into each bank, from under the bed up to the crest.
                            std::vector<TerrainEdit> wall;
                            for (float w = -(damHW + 15.0f); w <= damHW + 15.0f; w += 0.8f) {
                                const Vector3 dw = off(damDir, damSide, w);
                                const float bed = std::min(pf->elevation(dw), damLevel - damDepth);
                                for (float h = bed - 0.5f; h <= damCrest - 0.6f; h += 0.8f) wall.push_back({ at(dw, h), 1.3f, 5.0f, false });
                            }
                            ApplyTerrainEdits(*p, wall.data(), wall.size());
                            divWall = wall;
                            damUp = off(damDir, damDown, -12.0f); damDn = off(damDir, damDown, 12.0f);
                            damLvU = pf->waterSourceAt(damUp).level; damLvD = pf->waterSourceAt(damDn).level;
                            printf("bench: dam: dam crest %.2f m (%zu edits); generated level 12 m upstream %.2f, downstream %.2f\n", damCrest, wall.size(), damLvU, damLvD);
                            player.flyMode = true; sunOverride = true;
                            sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(damDir, Vector3Scale(damSide, 0.5f))));
                            damPhase = 1; damT0 = GetTime();
                        }
                    } else if (benchDamDivert) {
                        // --test-dam-divert (#656): the dam holds the river (the pool fills behind it). The river 200 m downstream
                        // (followed along the network) must dry to its bed; then the dam is dug out and the river must come back.
                        const double t = GetTime() - damT0;
                        if (divPhase == 0) {
                            Vector3 d = damDir;
                            for (int s = 0; s < 20; ++s) { const RiverNetwork::Hit h = net->nearest(d, 1e4f); d = off(d, h.down, 10.0f); }
                            const RiverNetwork::Hit h = net->nearest(d, 1e4f);
                            divDir = d; divPiece = h.piece; divLevel0 = h.level; divDepth = h.depth;
                            printf("bench: dam-divert: watching the river 200 m downstream (piece %d, %.1f m off its centreline): level %.2f m, depth %.2f m, half-width %.1f m\n", h.piece, h.dist, h.level, h.depth, h.halfWidth);
                            divPhase = 1;
                        }
                        // Watch the downstream river from above its bank.
                        const Vector3 dn = Vector3Normalize(Vector3CrossProduct(divDir, damSide));
                        world.camera.position = p->toWorld(at(off(off(divDir, damSide, damHW + 25.0f), dn, -15.0f), divLevel0 + 25.0f));
                        world.camera.target = p->toWorld(at(divDir, divLevel0 - divDepth));
                        world.camera.up = p->dirToWorld(divDir);
                        if (getenv("RAYLIBGAME_DIVERT_CAM")) { // diagnostics/captures: over the first dam, looking down the reach
                            world.camera.position = p->toWorld(at(off(off(damDir, damSide, damHW + 20.0f), damDown, -25.0f), damLevel + 22.0f));
                            world.camera.target = p->toWorld(at(off(damDir, damDown, 20.0f), damLevel - damDepth));
                            world.camera.up = p->dirToWorld(damDir);
                        }
                        player.position = world.renderOrigin + world.camera.position;
                        const RiverNetwork::Hit h = net->nearest(divDir, 1e4f);
                        const float flow = net->flowOf(divPiece);
                        static double lastPrint = -10.0;
                        if (t - lastPrint >= 5.0) {
                            lastPrint = t;
                            if (getenv("RAYLIBGAME_RIVERFLOW_LOG")) {
                                // Diagnostics: the water surface along the river through the first reach, 50 m up to 60 m down
                                // (a straight line along the dam site's downstream; "-": no simulated water), the natural level and the ground.
                                const RiverNetwork::NaturalScope natural;
                                std::string prof;
                                for (int s = -5; s <= 6; ++s) {
                                    const Vector3 d = off(damDir, damDown, s * 10.0f);
                                    const RiverNetwork::Hit hh = net->nearest(d, 1e4f);
                                    const float su = surfaceAt(d, hh.level);
                                    const std::string sv = su < -1e8f ? std::string("-") : std::string(TextFormat("%.2f", su));
                                    prof += TextFormat(" %+d:%s/%.2f/%.2f", s * 10, sv.c_str(), hh.level, pf->elevation(d));
                                }
                                printf("bench: dam-divert profile t=%3.0f:%s\n", t, prof.c_str());
                            }
                            printf("bench: dam-divert t=%3.0f s: downstream flow %.2f, level %.2f m (was %.2f), half-width %.1f; upstream surface %s; %d reaches released\n",
                                   t, flow, h.level, divLevel0, h.halfWidth, [&] { const float su = surfaceAt(damUp, damLvU); return su < -1e8f ? std::string("dry") : std::string(TextFormat("%.2f", su)); }().c_str(), Fluid::ReleasedReaches());
                        }
                        // Dig a dam out down to about the old ground (centres no lower than 1 m over it: a sphere 1.4 m across left
                        // a 2 m pit under the bed that had to fill before any water went on). elevation() is the carved channel's
                        // bed (7.59 m at seed 42's dam, level 8.56 m): the cut leaves no weir (#656 follow-up checked it).
                        auto digOut = [&](std::vector<TerrainEdit>& wall) {
                            for (TerrainEdit& e : wall) {
                                const Vector3 d = Vector3Normalize(e.center);
                                e.center = at(d, std::max(Vector3Length(e.center) - pf->radius, pf->elevation(d) + 1.0f));
                                e.carve = true; e.radius = 1.4f;
                            }
                            ApplyTerrainEdits(*p, wall.data(), wall.size());
                        };
                        static std::vector<TerrainEdit> divWall2;
                        static double phaseT = 0.0, upCutT = -1.0;
                        static bool loadKept = false;
                        if (divPhase == 1 && flow < 0.05f) {
                            divDryT = t;
                            TakeScreenshot("bench_dam_dry.png");
                            // A second dam 300 m upstream (followed along the network): the first reach's springs must follow
                            // the cut (they bring what arrives from upstream).
                            Vector3 d = damDir;
                            for (int s = 0; s < 30; ++s) { const RiverNetwork::Hit hu = net->nearest(d, 1e4f); d = off(d, hu.down, -10.0f); }
                            const RiverNetwork::Hit hu = net->nearest(d, 1e4f);
                            const Vector3 dn = Vector3Normalize(Vector3Subtract(hu.down, Vector3Scale(d, Vector3DotProduct(hu.down, d))));
                            const Vector3 sd = Vector3Normalize(Vector3CrossProduct(d, dn));
                            for (float w = -(hu.halfWidth + 15.0f); w <= hu.halfWidth + 15.0f; w += 0.8f) {
                                const Vector3 dw = off(d, sd, w);
                                const float bed = std::min(pf->elevation(dw), hu.level - hu.depth);
                                for (float hh = bed - 0.5f; hh <= hu.level + 2.0f; hh += 0.8f) divWall2.push_back({ at(dw, hh), 1.3f, 5.0f, false });
                            }
                            ApplyTerrainEdits(*p, divWall2.data(), divWall2.size());
                            printf("bench: dam-divert: dry 200 m downstream after %.0f s (bench_dam_dry.png); a second dam 300 m upstream (level %.2f m, %.1f m off the centreline)\n", t, hu.level, hu.dist);
                            divPhase = 3; phaseT = t;
                        } else if (divPhase == 3) {
                            const float sc = Fluid::ReleaseInflowScale();
                            if (t - phaseT >= 5.0 && fmod(t - phaseT, 5.0) < 0.05) printf("bench: dam-divert t=%3.0f s: the first reach's inflow %.2f of natural (%d reaches released)\n", t, sc, Fluid::ReleasedReaches());
                            if (sc >= 0.0f && sc < 0.1f) {
                                upCutT = t - phaseT;
                                // Save and load: the river's flow state must come back with the world.
                                const bool ok = SaveWorld(world) && LoadWorldSave(world);
                                printf("bench: dam-divert: the first reach's springs followed the cut upstream after %.0f s; save/load %s\n", upCutT, ok ? "ok" : "FAILED");
                                divPhase = 4; phaseT = t;
                            } else if (t - phaseT > 120.0) { printf("bench: dam-divert: the first reach's springs never followed the cut upstream (inflow %.2f)\n", sc); divPhase = 4; phaseT = t; }
                        } else if (divPhase == 4 && t - phaseT > 2.0) {
                            loadKept = flow < 0.1f;
                            printf("bench: dam-divert: after the load the river 200 m below reads flow %.2f (%s)\n", flow, loadKept ? "kept" : "lost");
                            digOut(divWall); digOut(divWall2);
                            printf("bench: dam-divert: both dams dug out\n");
                            divPhase = 2; phaseT = t;
                        } else if (divPhase == 2 && flow > 0.9f && divOpenT < 0.0) {
                            divOpenT = t;
                            TakeScreenshot("bench_dam_back.png");
                        }
                        // The return (#656 follow-up): by default a rate check, so the test fits ~3 min - the flow 200 m below
                        // passes 0.5 within 60 s of the dig-out and is still rising (or at 0.9) at 90 s, the level with it.
                        // --dam-divert-full waits for the flow to reach 0.9 and the level its natural height (up to 480 s).
                        static float flow60 = -1.0f, level0Dig = 0.0f, maxFlow60 = 0.0f;
                        static bool rateDone = false, rateOk = false;
                        if (divPhase == 2) {
                            const double since = t - phaseT;
                            if (since < 0.1) level0Dig = h.level;
                            if (since <= 60.0) maxFlow60 = std::max(maxFlow60, flow);
                            if (since >= 60.0 && flow60 < 0.0) { flow60 = flow; printf("bench: dam-divert: 60 s after the dig-out: flow %.2f (best %.2f), level %.2f m\n", flow, maxFlow60, h.level); }
                            if (!benchDamDivertFull && !rateDone && (since >= 90.0 || divOpenT >= 0.0)) {
                                rateDone = true;
                                const bool rising = divOpenT >= 0.0 || (flow60 >= 0.0f && flow >= flow60 - 0.02f);
                                rateOk = (divOpenT >= 0.0 || (maxFlow60 >= 0.5f && flow >= 0.5f)) && rising && h.level > level0Dig;
                                printf("bench: dam-divert: %.0f s after the dig-out: flow %.2f (best by 60 s %.2f, at 60 s %.2f), level %.2f m (at the dig-out %.2f) -> returning %s\n",
                                       since, flow, maxFlow60, flow60, h.level, level0Dig, rateOk ? "yes" : "no");
                            }
                        }
                        const bool full = benchDamDivertFull;
                        if ((full && (divOpenT >= 0.0 || t > 480.0)) || (!full && (rateDone || t > 300.0))) {
                            const bool back = full ? divOpenT >= 0.0 && fabsf(h.level - divLevel0) < 0.3f : rateOk;
                            printf("bench: dam-divert: dried after %s s, upstream cut reached the reach after %s s, state kept over save/load: %s, back %s s after the removal; level %.2f m (natural %.2f)%s -> %s\n",
                                   divDryT < 0.0 ? "never" : TextFormat("%.0f", divDryT), upCutT < 0.0 ? "never" : TextFormat("%.0f", upCutT), loadKept ? "yes" : "no",
                                   divOpenT < 0.0 ? (full ? "never" : "(not yet)") : TextFormat("%.0f", divOpenT - phaseT), h.level, divLevel0, full ? "" : TextFormat(", returning: %s", rateOk ? "yes" : "no"),
                                   divDryT >= 0.0 && upCutT >= 0.0 && loadKept && back ? "PASS" : "FAIL");
                            printf("bench: dam-divert: bench_dam_dry.png, bench_dam_back.png\n");
                            benchDam = false;
                        }
                    } else {
                        // Watch from low over the channel downstream, looking up the river over the dam into the pool behind it: the
                        // generated river in front, the simulated pool beyond at a grazing angle (#696: their reflections must match).
                        // The bank trees stand beside the channel, not in it (from 24 m over the bank they hid most of the water).
                        const Vector3 up = damDir;
                        const Vector3 eyeDir = off(damDir, damDown, 14.0f);
                        world.camera.position = p->toWorld(at(eyeDir, damCrest + 3.5f)); world.camera.target = p->toWorld(at(off(damDir, damDown, -20.0f), damLevel));
                        world.camera.up = p->dirToWorld(up);
                        player.position = world.renderOrigin + world.camera.position;
                        const double t = GetTime() - damT0;
                        const float su = surfaceAt(damUp, damLvU), sd = surfaceAt(damDn, damLvD);
                        if (t > 3.0) { damMaxU = std::max(damMaxU, su); damMinD = std::min(damMinD, sd < -1e8f ? damLvD - damDepth : sd); }
                        static double lastPrint = 0.0, spillT = -1.0, stallT = -1.0, bypassT = -1.0;
                        static float stallLevel = -1e9f;
                        if (spillT < 0.0 && damMaxU >= damCrest - 0.2f) spillT = t;
                        // Around the dam (#666, user: 'until it overtops or goes around the dam'): the pool rose, then held under the
                        // crest for 45 s while the springs kept feeding it - the water leaves another way. That passes too, and ends
                        // the test (it ran its full 420 s).
                        if (t > 20.0 && su > -1e8f) {
                            if (stallT < 0.0 || su > stallLevel + 0.05f) { stallT = t; stallLevel = su; }
                            else if (bypassT < 0.0 && spillT < 0.0 && t - stallT > 45.0 && damMaxU >= damLvU + 0.2f) bypassT = t;
                        }
                        if (t - lastPrint >= 5.0) {
                            lastPrint = t;
                            const Fluid::Stats s = Fluid::GetStats();
                            printf("bench: dam t=%3.0f s: upstream surface %s, downstream %s; %d reaches released, %d springs, %d bricks (%d awake), %.0f m^3, step %.2f ms (max %.2f); upstream cell %s\n",
                                   t, su < -1e8f ? "dry" : TextFormat("%.2f", su), sd < -1e8f ? "dry" : TextFormat("%.2f", sd), Fluid::ReleasedReaches(), Fluid::SourceCount(), s.bricks, s.activeBricks, s.volume, s.stepMs, s.maxStepMs,
                                   Fluid::DebugCell(*p, at(damUp, damLvU - 0.3f)).c_str());
                        }
                        if (t > 10.0 && t < 10.1) TakeScreenshot("bench_dam.png");
                        if ((spillT >= 0.0 && t > spillT + 25.0) || bypassT >= 0.0 || t > 420.0) {
                            damLateD = sd < -1e8f ? damLvD - damDepth : sd;
                            TakeScreenshot("bench_dam_end.png");
                            const int released = Fluid::ReleasedReaches();
                            // Backs up to the crest; spills: the pool then holds there (water leaves over the crest) instead of rising
                            // on. The reach below keeps the downstream river's level (its sink holds it: a flat river backs up from
                            // below too), so it is reported, not judged.
                            const bool backs = damMaxU >= damLvU + 0.2f && damMaxU >= damCrest - 0.25f, spills = spillT >= 0.0 && su > -1e8f && su <= damCrest + 0.35f;
                            const bool around = bypassT >= 0.0; // rose, then held under the crest while fed: it went around (#666)
                            printf("bench: dam: %d reaches released; upstream rose %.2f -> %.2f m (crest %.2f), spilled after %.0f s and holds at %.2f m; %s; downstream %.2f -> lowest %.2f, now %.2f m -> %s\n",
                                   released, damLvU, damMaxU, damCrest, spillT, su, around ? TextFormat("went around the dam (held at %.2f m from %.0f s)", stallLevel, stallT) : "did not go around",
                                   damLvD, damMinD, damLateD, released >= 1 && ((backs && spills) || around) && damDrift > 0.0f ? "PASS" : "FAIL");
                            printf("bench: dam: bench_dam.png, bench_dam_end.png\n");
                            benchDam = false;
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
                        // Out to 3 km, then 6 and 12 km (coarser) if nothing qualifies: a seed's geography is not fixed (#532
                        // moved every seed's land; seed 42's sea went out of the old 3 km square).
                        for (float span = 3000.0f; span <= 12000.0f && best > 1e29f; span *= 2.0f)
                        for (float e = -span; e <= span; e += 16.0f * span / 3000.0f)
                            for (float n = -span; n <= span; n += 16.0f * span / 3000.0f) {
                                const float dist = sqrtf(e * e + n * n);
                                if (dist >= best) continue;
                                const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                                PlanetField::WaterSource ws;
                                if (benchFlowSea) { ws.kind = 4; ws.level = 0.0f; if (pf->elevation(d) > -1.5f) continue; } // the sea, 1.5 m+ deep (#190)
                                else {
                                    ws = pf->waterSourceAt(d);
                                    if (!ws.kind || ws.kind == 3 || ws.level <= 2.0f || pf->elevation(d) > ws.level - 1.5f) continue;
                                    // Flat, wide water (#723): 6 m out, at least two of the four directions are water and none of it
                                    // stands more than 0.3 m off this level. After the 2026-11-17 terrain merges seed 42's nearest pick
                                    // was a 2 m mountain stream climbing 0.4 m a metre, which a level channel cut off its side never drains.
                                    int wide = 0; bool flat = true;
                                    for (const Vector3& o : { t1, Vector3Negate(t1), t2, Vector3Negate(t2) }) {
                                        const PlanetField::WaterSource wn = pf->waterSourceAt(Vector3Normalize(Vector3Add(d, Vector3Scale(o, 6.0f / pf->radius))));
                                        if (!wn.kind) continue;
                                        ++wide; flat = flat && fabsf(wn.level - ws.level) < 0.3f;
                                    }
                                    if (wide < 2 || !flat) continue;
                                }
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
                        if (best > 1e29f) { printf("bench: flow: no water with a low shore within 12 km of the start -> FAIL\n"); benchFlow = false; }
                        else {
                            const float shore = (float)flowStart;
                            printf("bench: flow: water %.0f m from the start, level %.1f m, shore %.0f m out\n", best, flowLevel, shore);
                            if (getenv("FLOWDBG"))
                                for (float s = shore - 8.0f; s <= shore + 2.0f; s += 2.0f) {
                                    const Vector3 q = Vector3Normalize(Vector3Add(flowDir, Vector3Scale(flowOut, s / pf->radius)));
                                    const PlanetField::WaterSource wq = pf->waterSourceAt(q);
                                    printf("bench: flow dbg column %+.0f m from the shore: source kind %d level %.2f; %s\n", s - shore, wq.kind, wq.level, pf->debugColumn(q).c_str());
                                }
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
                        // Every step's time seen until the water sleeps (#627): a single worst step is one-off spikes (a brick wake, the OS), so the gate is a rolling mean.
                        static std::vector<float> flowSteps; static std::vector<double> flowStepT; static long long flowLastStep = -1;
                        { const Fluid::Stats fs = Fluid::GetStats(); if (fs.steps != flowLastStep) { flowLastStep = fs.steps; if (fs.steps > 0 && flowSleepT < 0.0) { flowSteps.push_back((float)fs.stepMs); flowStepT.push_back(t); } } }
                        static double lastPrint = 0.0;
                        if (t - lastPrint >= 5.0) {
                            lastPrint = t;
                            const Fluid::Stats s = Fluid::GetStats();
                            printf("bench: flow t=%3.0f s: %d sites, %d bricks (%d awake), %.0f wet cells, %.1f m^3, step %.2f ms (max %.2f), publish %.2f ms; changes d%lld s%lld r%lld l%lld\n",
                                   t, s.sites, s.bricks, s.activeBricks, s.wetCells, s.volume, s.stepMs, s.maxStepMs, s.meshMs, s.dbg[0], s.dbg[1], s.dbg[2], s.dbg[3]);
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
                            // The gate (#627): the worst 1 s rolling mean of the step time (a sustained cost, not one spike); the
                            // working steps' percentiles and the single max are printed beside it.
                            float roll = 0.0f;
                            for (size_t a = 0, b = 0; b < flowStepT.size(); ++b) {
                                while (flowStepT[b] - flowStepT[a] > 1.0) ++a;
                                if (flowStepT[b] - flowStepT[0] < 1.0) continue; // whole windows only (the first step alone is no mean)
                                float sum = 0.0f; for (size_t k = a; k <= b; ++k) sum += flowSteps[k];
                                roll = std::max(roll, sum / (float)(b - a + 1));
                            }
                            std::vector<float> sorted;
                            for (float v : flowSteps) if (v > 0.01f) sorted.push_back(v);
                            std::sort(sorted.begin(), sorted.end());
                            auto pct = [&](double q) { return sorted.empty() ? 0.0f : sorted[std::min(sorted.size() - 1, (size_t)(q * sorted.size()))]; };
                            const int over30 = (int)(sorted.end() - std::upper_bound(sorted.begin(), sorted.end(), 30.0f));
                            const bool pass = flowReachT >= 0.0 && flowCaveT >= 0.0 && flowSleepT >= 0.0 && roll < 30.0f;
                            printf("bench: flow: %zu steps (%zu working) until asleep: worst 1 s mean %.2f ms (the gate, < 30 ms); working p50 %.2f, p95 %.2f, p99 %.2f ms, %d over 30 ms, max %.2f ms (#627)\n",
                                   flowSteps.size(), sorted.size(), roll, pct(0.5), pct(0.95), pct(0.99), over30, s.maxStepMs);
                            printf("bench: flow: reached the channel end after %.1f s, the chamber after %.1f s; water stands %.1f m over the chamber floor point (feed level %.1f m above it); asleep after %.1f s; step 1 s mean %.2f ms (max %.2f) -> %s\n",
                                   flowReachT, flowCaveT, caveTop, 14.5f, flowSleepT, roll, s.maxStepMs, pass ? "PASS" : "FAIL");
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
                        int n[(int)WoodSpecies::Count], cls[kTreeClasses];
                        if (p->trees) {
                            for (const float range : { 150.0f, 600.0f }) {
                                p->trees->speciesCounts(me, range, n, cls);
                                int total = 0; for (int s = 0; s < (int)WoodSpecies::Count; ++s) total += n[s];
                                printf("bench: forest: %d trees within %.0f m:", total, range);
                                for (int s = 0; s < (int)WoodSpecies::Count; ++s) printf(" %s %d%%", WoodSpeciesInfo((WoodSpecies)s).name, total ? 100 * n[s] / total : 0);
                                printf("; size classes S %d M %d L %d XL %d\n", cls[0], cls[1], cls[2], cls[3]);
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
                } else if (benchAnimalPath && (animalPathPhase > 0 || WorldIdle())) {
                    // --test-animal-path (#253): at night a dusk stalker inside a U-shaped pocket of build cubes (a 3.5 m
                    // wall 16 m wide, 10 m off, wings 10 m long reaching away; open on the far side) must find its way
                    // out and round to bite the player within 40 s. Without route planning it pushes into the wall.
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    if (animalPathPhase == 0) {
                        player.spawnOn(world, *p, StartSpawnDir(*p));
                        animalPathPhase = 10; animalPathStart = GetTime();
                    } else if (animalPathPhase == 10 && GetTime() - animalPathStart > 1.0) {
                        // A flat, dry spot near the start for the pocket (the eroded terrain is rough round the start): the
                        // test is about getting out of the pocket, not about the hills (#253).
                        const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                        const Vector3 up0 = Vector3Normalize(p->toLocal((player.position - world.renderOrigin).toFloat()));
                        const Vector3 u1 = Vector3Normalize(Vector3CrossProduct(up0, fabsf(up0.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                        const Vector3 u2 = Vector3CrossProduct(u1, up0);
                        auto dirAt = [&](Vector3 c, float a, float b) { return Vector3Normalize(Vector3Add(c, Vector3Scale(Vector3Add(Vector3Scale(u1, a), Vector3Scale(u2, b)), 1.0f / pf->radius))); };
                        bool found = false;
                        for (float r = 0.0f; r <= 600.0f && !found; r += 25.0f)
                            for (int k = 0; k < (r > 0.0f ? 12 : 1) && !found; ++k) {
                                const Vector3 c = dirAt(up0, r * cosf(k * 0.5236f), r * sinf(k * 0.5236f));
                                float lo = 1e9f, hi = -1e9f; bool wet = false;
                                for (int i = -2; i <= 2; ++i) for (int j = -2; j <= 2; ++j) {
                                    const Vector3 d = dirAt(c, i * 12.0f, j * 12.0f);
                                    const float e = pf->elevation(d);
                                    lo = std::min(lo, e); hi = std::max(hi, e);
                                    wet = wet || e < 1.5f || pf->inlandWaterNear(d, 0.0f);
                                }
                                if (!wet && hi - lo < 3.0f) {
                                    found = true;
                                    player.placeAt(world, *p, p->toWorld(Vector3Scale(c, pf->radius + pf->elevation(c) + 1.5f)));
                                    printf("bench: animal-path: a flat spot %.0f m from the start (%.1f m of relief)\n", r, hi - lo);
                                }
                            }
                        animalPathPhase = 1; animalPathStart = GetTime();
                    } else if (animalPathPhase == 1 && GetTime() - animalPathStart > 1.5 && player.grounded) {
                        const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                        const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                        const Vector3 up = Vector3Normalize(me);
                        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                        const Vector3 t2 = Vector3CrossProduct(t1, up);
                        auto onGround = [&](float a, float b) { const Vector3 d = Vector3Normalize(Vector3Add(me, Vector3Add(Vector3Scale(t1, a), Vector3Scale(t2, b)))); return Vector3Scale(d, pf->radius + pf->elevation(d)); };
                        if (!p->plots) p->plots = std::make_unique<PlotSet>();
                        Plot* pl = p->plots->plotFor(onGround(0.0f, 15.0f), 0.0f, true);
                        int cubes = 0;
                        auto wall = [&](float a0, float b0, float a1, float b1) {
                            const float len = sqrtf((a1 - a0) * (a1 - a0) + (b1 - b0) * (b1 - b0));
                            for (float s = 0.0f; s <= len; s += 0.25f) {
                                const Vector3 g = onGround(a0 + (a1 - a0) * s / len, b0 + (b1 - b0) * s / len);
                                const Vector3 gu = Vector3Normalize(g);
                                for (float h = -1.0f; h <= 3.5f; h += 0.25f) {
                                    const Int3 c = pl->cellAt(pl->toPlot(Vector3Add(g, Vector3Scale(gu, h))));
                                    if (pl->partAt(c) < 0 && pl->place(Part{ PartType::Cube, c, 0, 0 }) >= 0) ++cubes;
                                }
                            }
                        };
                        wall(-8.0f, 10.0f, 8.0f, 10.0f);  // the wall facing the player
                        wall(-8.0f, 10.0f, -8.0f, 20.0f); // the wings, away from the player
                        wall(8.0f, 10.0f, 8.0f, 20.0f);
                        worldDifficulty = 1;
                        Animals_SetThreats(true);
                        Animals_Reset();
                        survival.health = 100.0f;
                        world.sunDirection = p->dirToWorld(Vector3Negate(up));
                        const int n = Animals_Spawn(world, *p, "dusk_stalker", onGround(0.0f, 15.0f), 1, false, 0.1f);
                        float t = 0.0f, closest = 1e9f;
                        while (t < 40.0f && survival.health >= 100.0f) {
                            Animals_Update(world, 0.1f, 0.1f, false);
                            t += 0.1f;
                            if (Animals_Count() > 0) closest = std::min(closest, Vector3Distance(Animals_Position(0), me));
                        }
                        int plans, found; double avgMs, worstMs; Animals_PathStats(plans, found, avgMs, worstMs);
                        const bool ok = n == 1 && cubes > 100 && survival.health < 100.0f;
                        printf("bench: animal-path: %d cubes; the stalker %s after %.1f s (closest %.1f m); route plans %d (%d reached), %.2f ms each, worst %.2f ms: %s\n",
                               cubes, survival.health < 100.0f ? "came round and bit" : "never reached the player", t, closest, plans, found, avgMs, worstMs, ok ? "PASS" : "FAIL");
                        Animals_Reset();
                        survival.health = 100.0f;
                        benchAnimalPath = false;
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
                        // A fixed spawn direction by default: animals walk the real terrain (#223); with route planning (#253)
                        // any direction should work (loop --threat-seed to check).
                        SetRandomSeed((unsigned)threatSeed); // --threat-seed N: other spawn directions (#253; default 3)
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
                        const float fireStart = nearestNight(); // it was biting: right beside the player when the fire is lit
                        step(15.0f);
                        threatFire = nearestNight();
                        // No bites is the rule. It also backs off, but without pathfinding (#253) it can catch on the terrain
                        // on the way (#181: once 3 m off after 15 s, no bites): away past 12 m, or at least no closer.
                        check(survival.health >= 100.0f && (threatFire > 12.0f || threatFire >= fireStart - 0.5f), "a campfire keeps it off");
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
                        { int plans, found; double avgMs, worstMs; Animals_PathStats(plans, found, avgMs, worstMs);
                          printf("bench: threats: route plans %d (%d reached the goal), %.2f ms each, worst %.2f ms (#253)\n", plans, found, avgMs, worstMs); }
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
                        auto onGround = [&](int i, float tol = 0.15f) { const Vector3 a = Animals_Position(i), u = Vector3Normalize(a); return p->densityAt(Vector3Add(a, Vector3Scale(u, tol))) < ISO_LEVEL && p->densityAt(Vector3Subtract(a, Vector3Scale(u, tol))) >= ISO_LEVEL; }; // on the density surface (#223), not the height map
                        auto step = [&](float seconds) { for (float t = 0.0f; t < seconds; t += 0.1f) Animals_Update(world, 0.1f, 0.1f, false); };
                        animalNatural = Animals_Count(); // herds placed around the player as it arrived
                        check(animalNatural > 0, "herds spawn around the player");
                        Animals_Reset(); Animals_SetNaturalSpawns(false); // #373: on the sim clock the spawn check runs every second of the steps below: new herds would join the measures
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
                            check(inventory.count("berries") == g0, "the second press doesn't eat (#202)");
                            check(!Hotbar_TestEatHold(0.5f) && inventory.count("berries") == g0, "a short right-click doesn't eat");
                            check(Hotbar_TestEatHold(1.3f) && inventory.count("berries") == g0 - 1, "holding right-click 1.2 s eats one");
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
                                    for (int y = g - 2; y <= g + 2; y += 2) { // to g+3: the density surface can stand over the height map, and fleeing animals with routes (#253) find any gap
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
                            Animals_Aim(world, PlayerInput::None());
                            const bool aimedOpen = Animals_Aiming();
                            ApplyTerrainEdit(*p, TerrainEdit{ Vector3Lerp(eye, body, 0.5f), 0.7f, 10.0f });
                            Animals_Aim(world, PlayerInput::None());
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
                            // Animals find the ground in 0.25 m steps with a linear guess between them (animals.cpp Ground), so on
                            // the pit's dug edges they stand up to ~0.25 m off the true surface (#373: every "inside" sample was
                            // within 0.25 m of air). Count only what is off by more than that step.
                            constexpr float kPitTol = 0.3f;
                            int bad = 0, samples = 0, inAir = 0, buried = 0, worstRun = 0;
                            std::vector<int> run(64, 0); // consecutive off-ground samples per animal
                            for (float t = 0.0f; t < 60.0f; t += 0.1f) {
                                Animals_Update(world, 0.1f, 0.1f, false);
                                for (int i = 0; i < Animals_Count() && i < 64; ++i) {
                                    ++samples;
                                    if (onGround(i, kPitTol)) { run[i] = 0; continue; }
                                    ++bad; worstRun = std::max(worstRun, ++run[i]);
                                    const Vector3 a = Animals_Position(i), u = Vector3Normalize(a);
                                    if (p->densityAt(Vector3Subtract(a, Vector3Scale(u, kPitTol))) < ISO_LEVEL) ++inAir; else ++buried;
                                }
                            }
                            printf("bench: animals: around a dug pit for a minute: %d of %d samples more than %.2f m off the ground (%d in the air, %d inside it; longest run %d), hole drop %.2f m\n",
                                   bad, samples, kPitTol, inAir, buried, worstRun, dropped);
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
                    check(QuestBook_Tier() == 1 && QuestBook_TaskDone(1, 0) && QuestBook_TaskDone(1, 1) && !QuestBook_TaskDone(1, 3), "part of tier 1"); // kiln, copper ore done; copper ingots not (#230: eating moved to the side quests)
                    const std::string path = std::string(GetWorkingDirectory()) + "/quest_test.tmp";
                    if (FILE* f = fopen(path.c_str(), "w")) { saveExtrasHook(f); fprintf(f, "end_extras\n"); fclose(f); }
                    QuestBook_Reset(); inventory.clear();
                    if (FILE* f = fopen(path.c_str(), "r")) { loadExtrasHook(f); fclose(f); }
                    remove(path.c_str());
                    check(QuestBook_Tier() == 1 && QuestBook_TaskDone(1, 0) && QuestBook_TaskDone(1, 1) && !QuestBook_TaskDone(1, 3) && inventory.obtainedTotal("any_log") >= 2.0 && QuestBook_CountOf("eat") >= 1.0, "saved and loaded");
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
                    check(inventory.grid()[Inventory::kGrid].id == "telescope" && inventory.grid()[Inventory::kGrid + 2].id == "telescope", "new tools go to empty hotbar slots (#209)");
                    inventory.add("conveyor", 300);
                    int conv = 0;
                    for (const InvSlot& s : inventory.grid()) if (s.id == "conveyor") ++conv;
                    check(conv == 2 && inventory.count("conveyor") == 300, "blocks stack to 256");
                    inventory.add("stone", 9000);
                    check(inventory.count("stone") == 9000 && inventory.grid()[2].id == "stone" && inventory.grid()[2].centi == 819200, "resources stack to 8192 (main grid first; the hotbar holds the tools)");
                    for (int i = 0; i < 1000; ++i) inventory.addFractional("iron_ore", 0.0123f);
                    check(fabs(inventory.amount("iron_ore") - 12.30) < 0.02 && inventory.count("iron_ore") == 12, "fractional amounts to 0.01");
                    check(inventory.removeAmount("iron_ore", 0.25) && fabs(inventory.amount("iron_ore") - 12.05) < 0.02, "remove a fraction");
                    check(!inventory.remove("stone", 9001), "can't remove more than held");
                    inventory.sort();
                    check(KindOf(inventory.grid()[0].id) == ItemKind::Block && inventory.grid()[Inventory::kGrid].id == "telescope" && inventory.count("stone") == 9000, "sort: the main grid by kind, the hotbar kept, nothing lost");
                    // Full: every grid and hotbar cell holds 8192 regolith; the rest is refused and flagged.
                    inventory.clear();
                    const double left = inventory.addAmount("regolith", 8192.0 * Inventory::kSlots + 10.0);
                    check(fabs(left - 10.0) < 0.01 && inventory.overflowed, "a full grid refuses the rest");
                    check(inventory.offHand().empty(), "nothing overflows into the off hand (#165)");
                    // #495: 10 wide like the hotbar; a 9-wide save's hotbar (45-54) and off hand (55) load at 50-59 and 60.
                    check(Inventory::kGrid == 50 && Inventory::kOffHand == 60 && Inventory::SavedCell(44, 9) == 44 && Inventory::SavedCell(45, 9) == 50 &&
                          Inventory::SavedCell(55, 9) == 60 && Inventory::SavedCell(55, 10) == 55, "10 x 5 grid; old saves' hotbar and off hand remapped (#495)");
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
                    check(inventory.grid()[7].id == "copper_ore" && inventory.grid()[Inventory::kGrid].id == "telescope" && fabs(inventory.amount("copper_ore") - 12.34) < 0.001 && inventory.count("kiln") == 2, "save and load keep slots and fractions");
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
                        Hotbar_Reset(); // empty: put what the test uses on the hotbar (slot 1 stays empty: bare hands)
                        inventory.add("concrete", 20); inventory.add("oak_planks", 10); inventory.add("oxygen_canister", 1);
                        Hotbar_Bind(1, "concrete"); Hotbar_Bind(3, "oak_planks"); Hotbar_Bind(8, "oxygen_canister"); // the stacks move there (#209)
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
                        // #272 / #209: slots hold real stacks, so an item not carried can't be in one: an empty slot is bare
                        // hands (build mode ends). A new pickaxe lands in a free hotbar slot; moved to slot 6 and pressed,
                        // it goes into the hand and build mode ends.
                        Hotbar_Press(1, false);
                        Hotbar_Press(5, false);
                        const bool unownedKept = Hotbar_Selected() == 5 && !buildTool.active && terrainTool.heldTool.empty();
                        Hotbar_Press(1, false);
                        inventory.add("titanium_pickaxe", 1);
                        const bool landed = inventory.grid()[Inventory::kGrid].id == "titanium_pickaxe"; // the first free slot (key 1)
                        Hotbar_Bind(5, "titanium_pickaxe");
                        Hotbar_Press(5, false);
                        const bool ownedTaken = landed && Hotbar_Selected() == 5 && !buildTool.active && terrainTool.heldTool == "titanium_pickaxe" && Hotbar_Slot(0).empty();
                        printf("bench: hotbar: an empty slot in build mode is bare hands %d; a new pickaxe lands in a free slot, moved to key 6 and taken in hand %d (#272, #209)\n", unownedKept, ownedTaken);
                        Hotbar_Press(0, false); // key 1: an empty slot, bare hands
                        const bool toolMode = !buildTool.active && unownedKept && ownedTaken;
                        // #554: placing the last chest of the slot in hand ends build mode (the hand is empty).
                        inventory.add("chest", 1);
                        Hotbar_Bind(6, "chest");
                        Hotbar_Press(6, false);
                        Hotbar_Update(world, PlayerInput::None(), false);
                        const bool chestMode = buildTool.active && buildTool.part == (int)PartType::Chest;
                        const bool chestPlaced = BuildTool_Place(*p, *plot, { PartType::Chest, { 6, 0, 0 }, 0, 0 }, true);
                        Hotbar_Update(world, PlayerInput::None(), false);
                        const bool lastEnds = chestMode && chestPlaced && inventory.count("chest") == 0 && !buildTool.active;
                        printf("bench: hotbar: the last chest in hand placed: build mode %d -> %d (#554)\n", (int)chestMode, (int)buildTool.active);
                        // #531: F flies only in creative mode.
                        const bool creative0 = terrainTool.creative, fly0 = player.flyMode;
                        PlayerInput f; f.press(Act::Fly);
                        player.flyMode = false;
                        terrainTool.creative = false; player.flyAllowed = false; player.update(world, f, 1.0f / 60.0f);
                        const bool flyLocked = !player.flyMode;
                        terrainTool.creative = true; player.flyAllowed = true; player.update(world, f, 1.0f / 60.0f);
                        const bool flyCreative = player.flyMode;
                        player.flyAllowed = false; player.update(world, PlayerInput::None(), 1.0f / 60.0f); // creative off mid-flight: down
                        const bool flyDropped = !player.flyMode;
                        terrainTool.creative = creative0; player.flyMode = fly0; player.flyAllowed = true;
                        printf("bench: hotbar: F out of creative flies %d, in creative %d, creative off mid-flight lands %d (#531)\n", (int)!flyLocked, (int)flyCreative, (int)flyDropped);
                        Hotbar_Press(0, false);
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
                        ok = lastEnds && flyLocked && flyCreative && flyDropped && concreteMode && slope && c0 - c1 == 1 && beam && w0 - w1 == 4 /* a beam is 4 single cubes now, one plank each (#200) */ && block && (before >= ISO_LEVEL ? after >= ISO_LEVEL : fabsf(after - before) < 0.02f) /* the ground beside the block stays ground: #210's exact box carve pulls its density toward the surface (CSG), it doesn't dig */ && used && toolMode && loaded && woodParts == 4 && concreteParts == 2;
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
                            airTime += GameFrameTime();
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
                } else if (WorldIdle() && benchWalkSeconds > 0.0 && walkLegs > 0 && !walkLegReady) {
                    // Planet walk, after a teleport (#363): the leg starts once the new place has loaded (idle) and a
                    // few frames have settled the player on the ground.
                    Planetoid* start = dynamic_cast<Planetoid*>(world.objects[0].get());
                    // First idle frame: the place has loaded; set the player on it and let gravity back on.
                    if (++walkSettleFrames == 1) { player.spawnOn(world, *start, walkPickDir); player.flyMode = false; }
                    if (walkSettleFrames >= 30) {
                        walkLegReady = true; player.autoWalk = true;
                        walkLegStart = GetTime(); walkLegPath = walkPath; walkLast = start->toLocalD(player.position);
                        walkLegGrounded = walkLegFrames = 0; walkLegSlope = 0.0;
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
                        walkLegStart = walkStart; walkLegPath = 0.0; walkLeg = 0; walkLegReady = true;
                        walkLegs = walkTargetM > 0.0 ? std::max(1, (int)lround(walkTargetM / kWalkLegM)) : 0;
                    }
                    Planetoid* start = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const float d = start->densityAt(start->toLocal(Vector3Add(feet, Vector3Scale(player.up, 0.9f)))); // body centre
                    walkMinDensity = std::min(walkMinDensity, d);
                    walkMaxDensity = std::max(walkMaxDensity, d);
                    walkFrames++;
                    walkGroundedFrames += player.grounded ? 1 : 0;
                    walkLegFrames++; walkLegGrounded += player.grounded ? 1 : 0;
                    if (walkLegs > 0) {
                        // Deaths (falls off what the walker strides into) lose the leg, not the test (user, 2026-11-14).
                        if (survival.dead) walkLegDied = true;
                        // Falling through the world is the failure: the body 30 m under the height-map surface.
                        if (const PlanetField* pfw = dynamic_cast<const PlanetField*>(start->field.get())) {
                            const Vector3 lp = start->toLocal(feet);
                            const float alt = Vector3Length(lp);
                            if (alt > 1.0f && alt < pfw->radius + pfw->elevation(Vector3Scale(lp, 1.0f / alt)) - 30.0f) ++walkFellThrough;
                        }
                    }
                    if (player.grounded) walkLegSlope = std::max(walkLegSlope, (double)acosf(std::clamp(Vector3DotProduct(player.groundNormal, player.up), -1.0f, 1.0f)) * RAD2DEG);
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
                    {   // A jump of more than 20 m in a frame is a teleport (a respawn), not walking: not counted, reported.
                        const DVec3 here = start->toLocalD(player.position);
                        const double step = (here - walkLast).length();
                        if (step < 20.0) walkPath += step; else ++walkLegJumps;
                        walkLast = here;
                    }
                    if (!benchStandOnly && player.autoWalk) {
                        // Getting unstuck, as a player would (overhangs, cliffs and boulders are part of the world, the
                        // walker adapts): every 1.5 s, under 1 m of progress means stuck - jump first; still stuck, turn
                        // aside (alternating sides, wider each time) and walk on.
                        static double checkAt = 0.0; static DVec3 checkPos{}; static int stuck = 0; static int unsticks = 0;
                        if (GetTime() - walkLegStart < 0.1) { checkAt = GetTime() + 1.5; checkPos = start->toLocalD(player.position); stuck = 0; unsticks = 0; }
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
                    if (walkLegs > 0 && (walkPath - walkLegPath >= kWalkLegM || GetTime() - walkLegStart > kWalkLegMaxS || walkLegDied)) {
                        // A leg is done (#363): note it, then teleport to the next place, or finish.
                        const PlanetField* pfw = dynamic_cast<const PlanetField*>(start->field.get());
                        const Vector3 here = Vector3Normalize(start->toLocal((player.position - world.renderOrigin).toFloat()));
                        const double legM = walkPath - walkLegPath, legS = GetTime() - walkLegStart;
                        const double legGround = 100.0 * walkLegGrounded / std::max(1, walkLegFrames);
                        walkLegsDone += legM >= kWalkLegM && walkLegJumps == 0 && !walkLegDied ? 1 : 0;
                        printf("bench: planet walk leg %d/%d: %.0f m in %.1f s at %.0f m elevation, grounded %.0f%%, unstuck %d times, steepest ground %.0f deg%s%s\n",
                               walkLeg + 1, walkLegs, legM, legS, pfw ? pfw->elevation(here) : 0.0f, legGround, walkUnsticks, walkLegSlope,
                               walkLegDied ? " (died)" : legM >= kWalkLegM ? "" : " (timed out)", walkLegJumps ? " (displaced mid-leg)" : "");
                        if (++walkLeg < walkLegs) {
                            // The next place: land drawn from the seed; every other leg prefers ground above 300 m (slopes).
                            std::mt19937 rng(9173u * (unsigned)(walkLeg + 1) + worldGen.seed);
                            std::uniform_real_distribution<float> U(-1.0f, 1.0f);
                            Vector3 pick = here;
                            for (int tries = 0; tries < 400; ++tries) {
                                const Vector3 dq = Vector3Normalize({ U(rng), U(rng) * 0.8f, U(rng) });
                                const float e = pfw ? pfw->elevation(dq, 4.0f) : 10.0f;
                                if (e < 3.0f || e > 2500.0f || (pfw && pfw->waterElevation(dq, 4.0f) > e - 0.5f)) continue;
                                pick = dq;
                                if (walkLeg % 2 == 0 || e > 300.0f || tries > 300) break;
                            }
                            // Hover (fly mode) while the place loads: walking on half-loaded ground fell through it.
                            walkPickDir = pick;
                            player.spawnOn(world, *start, pick);
                            player.flyMode = true; walkLegJumps = 0; walkLegDied = false;
                            survival.health = 100.0f; survival.dead = false; // a death on the last leg: back to life here
                            player.autoWalk = false; walkLegReady = false; walkSettleFrames = 0;
                            walkLast = start->toLocalD(player.position);
                        } else {
                            player.autoWalk = false;
                            const double pct = 100.0 * walkGroundedFrames / std::max(1, walkFrames);
                            // Pass (user, 2026-11-14): most legs walked, and never through the world (body under the surface
                            // or inside solid ground). Steep or deadly legs are part of the terrain, not failures.
                            const bool pass = walkLegsDone * 2 > walkLegs && walkFellThrough == 0 && walkMaxDensity < ISO_LEVEL;
                            printf("bench: planet walk %d of %d legs of %.0f m completed in %.1f s, grounded %.0f%%, %d frames under the surface, deepest body density %.2f (solid >= %.2f): %s\n",
                                   walkLegsDone, walkLegs, kWalkLegM, GetTime() - walkStart, pct, walkFellThrough, walkMaxDensity, ISO_LEVEL, pass ? "PASS" : "FAIL");
                            walkTargetM = 0.0; walkLegs = 0; benchWalkSeconds = 0.0;
                        }
                    } else if (walkLegs == 0 && (GetTime() - walkStart > benchWalkSeconds || (walkTargetM > 0.0 && walkPath >= walkTargetM))) {
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
                } else if ((WorldIdle() && benchDigStrokes == 0 && !benchPaint && !benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && minePhase != 2 && benchFlySeconds <= 0.0 && !benchSuffocate && benchOrbitSeconds <= 0.0 && !benchRide && !benchPark && !benchGravity && !benchPrecision && !benchMantle && !benchAirJump && !benchMoon && !benchHole && !benchWall && !benchDoorway && !benchLedge && !benchTor && !benchMachines && !benchFall && !benchHover && !benchNetIds && !benchNet && !benchTwoPlayers && !benchFingerprint && !benchInterest && !benchTeleClose && !benchMenus && !benchWater && !benchFlyFollow && !benchHotbar && !benchInventory && !benchMap && !benchViewmodel && !benchQuests && !benchHunger && !benchFarm && !benchAnimals && !benchThreats && !benchAnimalPath && !benchDam && !benchFlowMeter && !benchPartWall && !benchPressure && benchShotTabs == 0 && !benchForest && !benchFlow && !benchPond && !benchBucket && !benchWaterSpeed && benchFluidLoad <= 0.0f && !benchFlowShade && !benchPump && !benchTerrain && !benchWind && !benchCampfire && !benchCaveMouth && !benchMountain && !benchShelter && !benchRoom && !benchFacing && !benchPosture && !benchPause && !benchBite && !benchPlacement && !benchShotDig && !benchShadow && !benchParticles && !benchBiomes && !benchOil && !benchSpring && !benchShotAnimals && !benchShotMachines && benchShotBiomes < 0 && !benchFarForest && !benchRecipes && !benchWeather && !benchDecor && !benchMaterials && !benchDamage && !benchCrafting && !benchSuit && !benchTools && !benchRivers && !benchRing && !benchRocket && !benchStation && !benchBlueprint && !benchSmallBody && !benchTrees && !benchChop && benchPersist == 0 && benchBuild == 0 && benchFactory == 0 && !benchCombat && !benchHost && !benchShower && !benchProgression && (benchUnfocus <= 0.0f || unfocusPhase >= 6)) || GetTime() > 600.0) {
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
                GpuMem_Report("end of bench"); // --vram (#302)
                if (TerrainArena_Available()) { const ArenaStats as = TerrainArena_Stats(); printf("bench: terrain arena: %d / %d vertices, %d / %d indices, %.1f chunk draws per multi-draw; uploads %lld staged, %lld direct\n", as.vertices, as.vertexCapacity, as.indices, as.indexCapacity, as.draws / (double)std::max(1, as.flushes), as.staged, as.direct); }
                if (TerrainArena_Available()) { long long pu[17]; TerrainArena_PaletteUse(pu); printf("bench: material palette (#480): chunks by distinct materials:"); for (int k = 1; k <= 16; ++k) if (pu[k]) printf(" %d:%lld", k, pu[k]); printf("\n"); }
                { const TerrainGpuStats ts = TerrainGpu_Stats(); const WorldStats& ws = GetWorldStats(); printf("bench: terrain gpu: %lld blocks in %lld batches, avg batch latency %.2f ms; chunk build avg %.2f ms (sampling part %.2f ms)\n", ts.blocks, ts.batches, ts.avgLatencyMs, ws.avgBuildMs, ws.avgSampleMs); }
                {   // Streaming (#258): how much work flying asked for and threw away; frame time spread.
                    const WorldStats& ws = GetWorldStats();
                    std::vector<float> ft = benchFrameTimes;
                    std::sort(ft.begin(), ft.end());
                    auto pct = [&](double q) { return ft.empty() ? 0.0f : ft[std::min(ft.size() - 1, (size_t)(q * ft.size()))] * 1000.0f; };
                    int over20 = 0; for (float f : ft) over20 += f > 0.020f ? 1 : 0;
                    printf("bench: streaming: %lld requests, %lld built, %lld skipped stale, %lld results dropped; %d nodes alive; frame ms p50 %.2f p95 %.2f p99 %.2f, %d frames over 20 ms\n",
                           ws.requests, ws.built, ws.skippedStale, ws.droppedResults, ws.nodes, pct(0.5), pct(0.95), pct(0.99), over20);
                    printf("bench: preload: %d ring nodes, %d built; cache: %d nodes, %lld evicted; arena %.1f M vertices; LOD cross-fades %lld\n",
                           ws.preloadNodes, ws.preloadReady, ws.cachedNodes, ws.evicted, TerrainArena_Stats().vertices / 1048576.0, ws.lodFades);
                    if (ws.waterBuilds || ws.waterFullBuilds)
                        printf("bench: water redraws: %lld chunk water rebuilds, avg %.2f ms on a worker, total %.0f ms; %lld full builds for water; caches %.1f MB (#360)\n",
                               ws.waterBuilds, ws.waterBuilds ? ws.waterMs / ws.waterBuilds : 0.0, ws.waterMs, ws.waterFullBuilds, ws.waterCacheBytes / 1048576.0);
                    if (ws.waterBuilds || ws.waterFullBuilds) {
                        long long ph[6]; ChunkWaterPhases(ph);
                        printf("bench: water mesh phases: sampling %.0f ms, marching %.0f ms, shore clip %.0f ms, vertices %.0f ms; %.1f%% of %lld density samples exact (#381)\n",
                               ph[0] / 1e6, ph[1] / 1e6, ph[2] / 1e6, ph[3] / 1e6, ph[5] ? 100.0 * ph[4] / ph[5] : 0.0, ph[5]);
                    }
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
        }
        // Unfocused (alt-tab, #510): the frame still runs, so streaming, uploads, the spawn preload and the loading screen
        // keep going (skipping whole frames stalled the loading, and the backlog cost more on return), and so does the
        // whole game (#525: time of day, physics, the factory; only the pause menu pauses); input and the cursor lock are off. Benches, --bench-window and --play-new runs open
        // unfocused (#387) and ignore focus; --bench-unfocus scripts it (its own window is never minimized or raised).
        const bool focused = benchUnfocus > 0.0f ? !benchUnfocused : (bench || benchWindow || playNewSeconds > 0.0 || IsWindowFocused());
        {
            // The cursor (#204, #205): locked while playing (DisableCursor: hidden, raw relative motion, no
            // warping) and free in menus and UIs, switched only when that changes. It used to be hidden and warped
            // back to the centre with SetMousePosition every frame: a blocking OS call per frame (benches skip it,
            // so they never showed the frame-rate drop while moving and looking) that also lost the motion between
            // warps (choppy look at low frame rates).
            const bool wantLock = !bench && focused && gameState == 1 && !(DebugUI_IsOpen() || GameUI_InventoryOpen() || Nav_MapOpen());
            if (wantLock != cursorLocked) {
                if (wantLock) DisableCursor(); else EnableCursor();
                cursorLocked = wantLock;
            }
        }
        if (!bench && IsKeyPressed(KEY_F7)) Fluid::SetDebugView(!Fluid::DebugView()); // #695: water coloured by its speed
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
        const bool gameInput = !bench && focused && !uiOpen && gameState == 1 && !loading;
        world.cameraInputEnabled = player.flyMode && !uiOpen && !ship.piloted && focused;  // free-fly camera only in fly mode; not behind another window (#539)
        const bool alive = !survival.dead;
        // One clamped frame time for everything that moves (bodies, ship, player): clamping only the
        // player let the ground spin away under it during a slow frame.
        const float frameDt = std::min(GameFrameTime(), 0.1f); // fixed while recording (#375)
        // The animation clock (GameClock) and the simulation steps stop behind the single-player pause menu (#511: the world,
        // the factory, fluids and time of day all wait); a host keeps running for its clients, and so does an unfocused window (#525).
        const bool paused = gameState == 2 && gNetRole == NetRole::Offline;
        if (!paused) gGameClock += frameDt;
        // A spawn or teleport: the loading screen until the terrain there is built (play and tests pause).
        if (player.spawns != seenSpawns) { seenSpawns = player.spawns; if (gameState == 1) { loading = true; spawnLoad = true; loadStart = GetTime(); loadFrames = 0; } }
        GpuTimer_CpuMark("bench-logic");
        // Benches that place the camera themselves (shots, tours): take that as the pose the fly-follow carries, or the
        // carry in player.update put the camera back where it was last frame (#298: since 4720ba4 the biome and machine
        // shots all showed the start, and --bench-fly never moved).
        if (bench && player.flyMode && player.flyFollow) player.recordFly(world);
        // The simulation steps (#315): the gameplay below runs in fixed steps of SimDt() (60 Hz), as many as this frame's
        // time owes (none on most frames above 60 FPS, two at 30 FPS). Key and button edges are latched for them in
        // PlayerInput_Step() (playerinput.h); the mouse look turns the view once per frame, the camera interpolates.
        PlayerInput_BeginFrame();
        const int simSteps = SimClock_BeginFrame(paused ? 0.0f : frameDt); // unfocused too (#525): the world keeps running behind another window
        if (gameState == 1 && gameInput && alive && !loading && !ship.piloted && !Rocket_Boarded() && !player.flyMode) player.look(PlayerInput_Frame().look);
        for (int simStep = 0; simStep < simSteps; ++simStep) {
        const float simDt = SimDt();
        PlayerSession& session = LocalSession(); // single-player: one session (#313); a server steps each one
        session.input = PlayerInput_Step();
        const PlayerInput& stepIn = session.input; // the local player's input for this step (#312)
        const bool alive = !survival.dead; // (may change within the frame's steps)
        // Time of day (#275): the debug UI's time scale multiplies --time-scale; its jumps (+1 h, dawn, ...) are applied here
        // in one step, so orbits and world time move too, and the riders below are carried along.
        const bool playing = gameState == 1 && !loading;
        const double timeSkip = playing ? gTimeSkipSeconds : 0.0;
        if (playing) gTimeSkipSeconds = 0.0;
        advanceBodies(world, playing ? (float)(simDt * timeScale * gDebugTimeScale * Rocket_Warp() + timeSkip) : 0.0f); // bodies only spin while playing (warped with a coasting rocket); before riders update
        if (gameState == 1) {
            ship.lastImpact = 0.0f;
            Station_Update(world); // stations on their orbits, before anything rides on them
            if (ship.piloted) {
                ship.update(world, PlayerInput_If(gameInput && alive, stepIn), simDt);
                player.position = ship.position; // the pilot rides inside
                player.velocity = ship.velocity;
            } else if (Rocket_Boarded()) {
                ship.update(world, PlayerInput::None(), simDt);
                Rocket_Update(world, PlayerInput_If(gameInput && alive, stepIn), simDt * timeScale * Rocket_Warp()); // flight (and its camera)
            } else {
                ship.update(world, PlayerInput::None(), simDt);
                Suit_Apply(world, simDt); // the worn suit and modules: walk, jetpack, grip, the scanner
                if (alive && !loading) {
                    // Benches drive the player through the same input struct as the keyboard (#312): autoWalk holds W,
                    // autoJump presses Space once.
                    PlayerInput scripted;
                    if (player.autoWalk) scripted.hold(Act::Forward);
                    if (player.autoJump) { scripted.press(Act::Jump); player.autoJump = false; }
                    if (player.autoCrouch) scripted.hold(Act::Down);
                    const PlayerInput& moveIn = gameInput ? stepIn : scripted.enabled ? scripted : PlayerInput::None();
                    Combat_Step(session, simDt, moveIn); // stamina and the block before this step's movement (#165)
                    // Physics in the player's own frame (#314): render space centred on the feet for the step, whatever the
                    // camera and the floating origin do (a server runs every player this way).
                    player.flyAllowed = terrainTool.creative || bench; // F flies in creative mode only (#531; benches fly their cameras)
                    RenderOriginScope frame(world, player.position);
                    player.update(world, moveIn, simDt);
                }
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
                Rocket_Update(world, PlayerInput_If(gameInput && alive, stepIn), simDt); // on foot: aim at a rocket (E boards, I fuels)
            }
            if (ship.piloted && ship.lastImpact > 12.0f) survival.health -= (ship.lastImpact - 12.0f) * 3.0f; // crash damage
            if (!loading) Survival_Update(world, simDt, PlayerInput_If(gameInput, stepIn));
            if (!loading) inventory.cool(simDt); // hot items cool (#580)
            if (!loading) Farm_Update(world, simDt * timeScale * gDebugTimeScale * Rocket_Warp()); // crops grow on world time
            if (!loading) Animals_Update(world, simDt, simDt * timeScale * gDebugTimeScale * Rocket_Warp(), gameInput && alive && !ship.piloted);
            if (!loading) Weather_Update(world, simDt * timeScale * gDebugTimeScale); // the sky, rain and snow, lightning
            Factory::Sync(world);                                // machine layouts in, items out
            Factory::Advance(simDt);                             // the factory thread runs its ticks on the sim clock (#315)
            GpuTimer_Mark("w-sim-pre");                          // (benches: the fluid's own share below)
            Fluid::Update(world, simDt);                         // flowing water: commands in, a step job, meshes up
            GpuTimer_Mark("fluid");                              // its GPU step dispatches and readback (#642 perf)
            Waterworks_Update(world, simDt);                     // pumps send water down pipes to outlets
            Objectives_Update(world);
        }
        PlayerInput_EndStep();
        SimClock_EndStep();
        }
        if (gameState == 1) {
            // Per frame: the view, the HUD's helpers and client-only effects (not simulation, #315).
            if (ship.piloted) ship.applyChaseCamera(world);
            if (!loading) Nav_Update(world, !bench && Nav_MapOpen(), !bench && gameInput && !Nav_MapOpen()); // exploration, the map, the minimap zoom (#278)
            if (!loading) Viewmodel_Update(world, Combat_HitStop(LocalSession(), settings.reduceMotion) ? 0.0f : frameDt); // hit-stop (#165)
            if (!loading) Structures_Update(world); // ruins, wrecks and wonders near the player (#164)
            if (!loading) Fishing_Update(world, frameDt, gameInput && alive && !ship.piloted); // (#176)
            if (!loading) Trapping_Update(world); // snares (#176)
            if (!loading) Events_Update(world, frameDt); // meteor showers (#177)
            if (!hostSave.empty()) Host_Update(world); // --host (#311): answer joins
            if (!loading) Ambience_Update(world, frameDt); // rain splashes, fireflies, machine smoke
            if (!loading) Waterfalls_Update(world, frameDt); // falls near the camera: sheets, mist and spray (#540)
            Particles_Update(world, frameDt); // smoke, dust, exhaust, splashes (particles.h)
            // Autosave every 5 minutes of play (incremental: only regions changed since the last save).
            static double lastAutosave = GetTime();
            if (!bench && alive && GetTime() - lastAutosave > 300.0) {
                lastAutosave = GetTime();
                if (SaveWorld(world)) GameUI_Toast("Autosaved");
            }
            if (bench && lookDown && !ship.piloted) {
                const Vector3 t = Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)));
                player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(t), 0.8f), player.up)); // ~50 deg down
            }
            if (bench && benchTelescope > 0.0f && !ship.piloted) AimAtTelescopeTarget(world);
            if (!ship.piloted && !Rocket_Boarded()) { player.applyToCamera(world); Combat_ShakeCamera(world, LocalSession(), settings.reduceMotion); } // hit shake (#165)
            else if (Rocket_Boarded()) Rocket_View(world, PlayerInput_If(gameInput && alive, PlayerInput_Frame()), SimAlpha());
            Machine_Update(world, PlayerInput_If(gameInput && alive && !ship.piloted, PlayerInput_Frame())); // after the camera: aims with this frame's view (#285, as #222)
            Telescope_Update(world, gameInput && alive && !Rocket_Boarded(), !ship.piloted && !player.flyMode && alive && !buildTool.active); // aboard, T moves fuel
            if (bench && benchTelescope > 0.0f) { Telescope_SetActive(world, true); telescope.zoom = benchTelescope; Telescope_Update(world, false, !ship.piloted); }
        }

        GpuTimer_CpuMark("game-update");
        world.updateFloatingOrigin(); // keep the camera near render-space 0 (float precision)
        // Upload time a frame: 3 ms in play. On the loading screen (#387) most of the frame: 70% of a frame of at least
        // 16.7 ms (a 240 FPS cap would leave 2.9 ms; there is nothing to show but progress), so uploads keep up with the
        // workers instead of following the frame rate.
        // Unfocused (#510) takes the loading budget too (its frames are capped at ~30 a second).
        const float uploadBudgetMs = (loading || !focused) ? std::min(30.0f, 0.7f * std::max(1000.0f / 60.0f, GetFrameTime() * 1000.0f)) : 3.0f; // at most 30 ms: Esc stays responsive
        { const double w0 = GetTime(); worldUpdate(world, uploadBudgetMs); if (benchMeasuring) benchUpdateMs += (GetTime() - w0) * 1000.0; } // LOD + streaming (bodies were spun at the start of the frame)

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
        while (world.lights.size() < 4) world.lights.push_back(Light{});
        if (gameState == 1) { Ambience_GatherLights(world); Events_AddLights(world); } else world.localLights.clear(); // lamps, campfires, fire machines (#139); meteor impacts (#177)
        {   // Slot 3: a torch in hand (#218), a flickering point light at the hand; it burns down while held in the game.
            static float burnt = 0.0f;
            Light& torch = world.lights[3];
            torch.type = LIGHT_POINT;
            const std::string& held = Hotbar_Selected() >= 0 ? Hotbar_Slot(Hotbar_Selected()) : std::string();
            const bool fibre = held == "fibre_torch", lit = (fibre || held == "torch") && inventory.count(held) > 0;
            torch.enabled = gameState == 1 && lit && !ship.piloted && !survival.dead;
            if (torch.enabled) {
                const float t = (float)GetTime();
                const float flicker = settings.fireFlicker ? 0.85f + 0.1f * sinf(t * 13.0f) + 0.05f * sinf(t * 31.7f + 1.3f) : 0.85f; // #558: steady when off
                const Vector3 fwd = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
                const Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, world.camera.up));
                torch.position = Vector3Add(Vector3Add(world.camera.position, Vector3Scale(right, 0.3f)), Vector3Add(Vector3Scale(world.camera.up, -0.15f), Vector3Scale(fwd, 0.35f)));
                torch.target = { fibre ? 7.0f : 12.0f, 0.0f, 0.0f };   // range (m)
                torch.color = fibre ? Color{ 255, 165, 85, (unsigned char)(110 * flicker) } : Color{ 255, 180, 105, (unsigned char)(150 * flicker) };
                burnt += GetFrameTime();
                if (burnt >= (fibre ? 90.0f : 300.0f)) { // one burns out; the next in the stack is lit
                    burnt = 0.0f;
                    inventory.remove(held, 1);
                    if (inventory.count(held) == 0) GameUI_Toast(fibre ? "Your fibre torch burnt out" : "Your torch burnt out");
                }
            }
        }
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
        if (gameState == 1) Hotbar_Update(world, PlayerInput_If(gameInput && alive, PlayerInput_Frame()), ship.piloted);
        if (gameState == 1) BuildTool_Update(world, PlayerInput_If(gameInput && alive, PlayerInput_Frame()), !ship.piloted && !Rocket_Boarded() && !player.flyMode && !telescope.active && alive);
        if (gameState == 1 && !loading) Animals_Aim(world, PlayerInput_If(gameInput && alive && !ship.piloted, PlayerInput_Frame())); // this frame's camera (#222)
        if (gameState == 1) TerrainTool_Update(world, PlayerInput_If(gameInput && alive && !Animals_Aiming(), PlayerInput_Frame()), !ship.piloted && !Rocket_Boarded() && !telescope.active && !buildTool.active); // F5/F9 work aboard too
        {   // Tool wear: each stroke, chop or hit with the held tool (the swing counter covers them all).
            static int lastSwings = 0;
            if (terrainTool.swings > lastSwings && !terrainTool.creative && !terrainTool.terraform && !terrainTool.heldTool.empty())
                ToolWear_OnStroke(terrainTool.heldTool, terrainTool.swings - lastSwings);
            lastSwings = terrainTool.swings;
        }

        // --test-foundation (#210): a block placed half into the ground 8 m ahead clears its own cell, and the ground just
        // outside it (0.3 m beside and below) is untouched.
        if (benchFoundation && gameState == 1 && !loading && WorldIdle()) {
            benchFoundation = false;
            Planetoid* p = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
            const Vector3 feet = (player.position - world.renderOrigin).toFloat();
            const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
            TerrainHit th;
            const Vector3 above = Vector3Add(Vector3Add(feet, Vector3Scale(flat, 8.0f)), Vector3Scale(player.up, 6.0f));
            if (!p || !RaycastTerrain(world, Ray{ above, Vector3Negate(player.up) }, 200.0f, th)) { printf("bench: foundation: no ground -> FAIL\n"); }
            else {
                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                Plot* plot = p->plots->plotFor(th.local, 0.0f, true);
                const Int3 cell = plot->cellAt(plot->toPlot(th.local));
                auto at = [&](float fx, float fy, float fz) { return plot->toBody({ (cell.x + fx) * kCellSize, (cell.y + fy) * kCellSize, (cell.z + fz) * kCellSize }); };
                // Probes: the cell's centre and lower middle (inside); 0.3 m outside each side face and under the bottom.
                const float o = 0.3f / kCellSize;
                const Vector3 inside[2] = { at(0.5f, 0.5f, 0.5f), at(0.5f, 0.15f, 0.5f) };
                const Vector3 outside[5] = { at(-o, 0.2f, 0.5f), at(1.0f + o, 0.2f, 0.5f), at(0.5f, 0.2f, -o), at(0.5f, 0.2f, 1.0f + o), at(0.5f, -o, 0.5f) };
                float before[5];
                for (int k = 0; k < 5; ++k) before[k] = p->densityAt(outside[k]);
                const bool placed = BuildTool_Place(*p, *plot, { PartType::Cube, cell, 0 }, false);
                bool cleared = true, kept = true;
                for (const Vector3& v : inside) cleared &= p->densityAt(v) < ISO_LEVEL;
                float worst = 0.0f;
                for (int k = 0; k < 5; ++k) { const float d = fabsf(p->densityAt(outside[k]) - before[k]); worst = std::max(worst, d); kept &= d < 0.02f; }
                printf("bench: foundation: placed %s, its cell cleared %s, ground 0.3 m outside unchanged %s (largest change %.3f) -> %s\n",
                       placed ? "yes" : "NO", cleared ? "yes" : "NO", kept ? "yes" : "NO", worst, placed && cleared && kept ? "PASS" : "FAIL");
            }
        }
        if (testInvCursor && gameState == 1 && !loading && WorldIdle()) {
            testInvCursor = false;
            if (Planetoid* p0 = dynamic_cast<Planetoid*>(world.objects[0].get())) { player.flyMode = false; player.spawnOn(world, *p0, StartSpawnDir(*p0)); player.nearest = p0; } // on the ground: a drop lands at the feet (#564)
            GameUI_TestCursorStack();
        }
        // --shot-preview (#197): a slope previewed on the ground next to a block (green), then into the block (red).
        if (shotPreview && gameState == 1 && !loading && WorldIdle()) {
            static int f = 0;
            static Plot* plot = nullptr;
            static Planetoid* p = nullptr;
            ++f;
            if (f == 1) {
                p = dynamic_cast<Planetoid*>(world.objects[0].get());
                player.spawnOn(world, *p, StartSpawnDir(*p));
                const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                TerrainHit th;
                RaycastTerrain(world, Ray{ Vector3Add(Vector3Add(feet, Vector3Scale(flat, 6.0f)), Vector3Scale(player.up, 20.0f)), Vector3Negate(player.up) }, 60.0f, th);
                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                plot = p->plots->plotFor(th.local, 0.0f, true);
                BuildTool_Place(*p, *plot, { PartType::Cube, { 0, 1, 0 }, 0 }, false);
                buildTool.active = true;
                buildTool.part = (int)PartType::Slope;
                player.nearest = p;
            }
            if (plot && f >= 2) {
                const bool red = f > 40;
                const Vector3 eye = red ? Vector3{ 2.6f, 1.6f, 1.9f } : Vector3{ 2.2f, 2.0f, -1.4f }, at = red ? Vector3{ 0.25f, 0.7f, 0.25f } : Vector3{ 1.2f, 0.5f, 0.25f };
                world.camera.position = p->toWorld(plot->toBody(eye));
                world.camera.target = p->toWorld(plot->toBody(at));
                world.camera.up = p->dirToWorld(plot->up);
                BuildTool_Update(world, PlayerInput::None(), true);
                if (f == 38) { TakeScreenshot("bench_preview.png"); printf("bench: shot: bench_preview.png (a slope on the ground)\n"); }
                if (f == 78) { TakeScreenshot("bench_preview_red.png"); printf("bench: shot: bench_preview_red.png (into the block)\n"); shotPreview = false; }
            }
        }
        if (testVein && gameState == 1 && !loading && WorldIdle()) { // #207
            testVein = false;
            if (auto* vp = dynamic_cast<Planetoid*>(world.objects[0].get())) TerrainTool_TestVein(*vp, StartSpawnDirLocal(*vp));
        }
        // --test-bug-report (#169): a scripted report to RAYLIBGAME_BUG_RELAY; BugReport_Update prints the relay's answer.
        if (testBugReport > 0 && gameState == 1 && !loading && WorldIdle()) {
            if (testBugReport == 1) BugReport_TestSend(world);
            ++testBugReport;
            if (testBugReport > 3 && !BugReport_Sending()) { printf("bench: bug report: sent (see the 'bug report:' line)\n"); testBugReport = 0; }
        }
        // --test-container-cursor (#192): a chest placed 6 m ahead; its slots behave as the inventory's (cursor, halves, ones, Shift).
        if (testContainerCursor && gameState == 1 && !loading && WorldIdle()) {
            testContainerCursor = false;
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const Vector3 feet = (player.position - world.renderOrigin).toFloat();
            const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
            TerrainHit th;
            if (p && RaycastTerrain(world, Ray{ Vector3Add(Vector3Add(feet, Vector3Scale(flat, 6.0f)), Vector3Scale(player.up, 6.0f)), Vector3Negate(player.up) }, 200.0f, th)) {
                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                Plot* plot = p->plots->plotFor(th.local, 0.0f, true);
                const Int3 cell = plot->cellAt(plot->toPlot(Vector3Add(th.local, Vector3Scale(plot->up, 0.1f))));
                if (BuildTool_Place(*p, *plot, { PartType::Chest, cell, 0 }, false)) {
                    Factory::Sync(world); Factory::WaitTicks(2); Factory::Sync(world);
                    GameUI_TestContainerCursor(world, p->id, plot->id, plot->partAt(cell), (int)PartType::Chest);
                } else printf("bench: container cursor: the chest would not go down -> FAIL\n");
            } else printf("bench: container cursor: no ground -> FAIL\n");
        }
        // --shot-chest (#495): a chest 6 m ahead holding a few stacks, open beside the inventory.
        if (shotChest && gameState == 1 && !loading && WorldIdle()) {
            static int f = 0;
            static uint64_t chestBody = 0; static uint32_t chestPlot = 0; static int chestPart = -1;
            if (++f == 2) { // out of the bench's free-fly camera, onto the ground (a chest screen closes out of reach)
                if (Planetoid* p0 = dynamic_cast<Planetoid*>(world.objects[0].get())) { player.flyMode = false; player.spawnOn(world, *p0, StartSpawnDir(*p0)); }
            }
            if (f == 40) {
                Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                TerrainHit th;
                if (p && RaycastTerrain(world, Ray{ Vector3Add(Vector3Add(feet, Vector3Scale(flat, 6.0f)), Vector3Scale(player.up, 6.0f)), Vector3Negate(player.up) }, 200.0f, th)) {
                    if (!p->plots) p->plots = std::make_unique<PlotSet>();
                    Plot* plot = p->plots->plotFor(th.local, 0.0f, true);
                    const Int3 cell = plot->cellAt(plot->toPlot(Vector3Add(th.local, Vector3Scale(plot->up, 0.1f))));
                    if (BuildTool_Place(*p, *plot, { PartType::Chest, cell, 0 }, false)) { chestBody = p->id; chestPlot = plot->id; chestPart = plot->partAt(cell); }
                }
                if (chestPart < 0) { printf("bench: shot: no chest placed -> FAIL\n"); shotChest = false; }
            }
            if (f == 45 && chestPart >= 0) {
                Factory::Sync(world);
                Factory::Insert(chestBody, chestPlot, chestPart, "oak_log", 40);
                Factory::Insert(chestBody, chestPlot, chestPart, "stone", 64);
                Factory::Insert(chestBody, chestPlot, chestPart, "iron_ore", 12);
                inventory.add("stone_pickaxe", 1); inventory.add("oak_log", 12);
                GameUI_OpenContainer(chestBody, chestPlot, chestPart, (int)PartType::Chest);
            }
            if (f == 90) { TakeScreenshot("bench_chest.png"); printf("bench: shot: bench_chest.png\n"); shotChest = false; }
        }
        // --shot-inventory (#219): the inventory with long item names (cut to fit; they scroll on hover).
        // --shot-glow (#580): hot items glow in the inventory and on the hotbar: a fresh iron ingot stack (yellow-white),
        // a half-cooled copper one (orange) and a nearly cool one (dull red), the fresh one also on the hotbar; the
        // inventory open (bench_glow.png).
        if (shotGlow && gameState == 1 && !loading && WorldIdle()) {
            static int f = 0;
            if (++f == 5) {
                inventory.addAmount("iron_ingot", 20, 1.0f); inventory.addAmount("copper_ingot", 12, 0.5f); inventory.addAmount("bronze_ingot", 6, 0.12f);
                Hotbar_Bind(2, "iron_ingot");
                inventory.addAmount("iron_ingot", 5, 1.0f); // a second stack in the grid (the first went to the hotbar)
                GameUI_ToggleInventory();
            }
            if (f == 30) {
                TakeScreenshot("bench_glow.png");
                const ItemGlow hot = Item_Glow(inventory.grid()[Inventory::kGrid + 2]);
                printf("bench: glow: hotbar iron ingot heat %.2f, glow %.2f (r %d g %d b %d); screenshot bench_glow.png\n",
                       inventory.grid()[Inventory::kGrid + 2].heat, hot.strength, hot.color.r, hot.color.g, hot.color.b);
                GameUI_ToggleInventory();
                shotGlow = false;
            }
        }
        if (shotInventory && gameState == 1 && !loading && WorldIdle()) {
            static int f = 0;
            if (++f == 5) { inventory.add("oxygen_tank_module", 1); inventory.add("stone_pickaxe", 1); inventory.add("steam_macerator", 1); inventory.add("oak_log", 12); inventory.grid()[Inventory::kOffHand] = InvSlot{ "wooden_shield", 100 }; GameUI_ToggleInventory(); }
            if (f == 30) { TakeScreenshot("bench_inventory.png"); printf("bench: shot: bench_inventory.png\n"); shotInventory = false; }
        }
        // --test-pickups (#208): dug stone in three bits (one merged line), a log and a tool; a screenshot while they show.
        if (benchPickups && gameState == 1 && !loading && WorldIdle()) {
            ++pickupFrames;
            if (pickupFrames == 5) { inventory.addFractional("stone", 1.25f); inventory.add("oak_log", 3); }
            if (pickupFrames == 15) inventory.addFractional("stone", 1.25f);
            if (pickupFrames == 25) { inventory.addFractional("stone", 0.5f); inventory.add("stone_pickaxe", 1); }
            if (pickupFrames == 40) {
                TakeScreenshot("bench_pickups.png");
                printf("bench: pickups: screenshot bench_pickups.png (expect +3 Stone, +3 Oak Log, +1 Stone Pickaxe)\n");
                benchPickups = false;
            }
        }
        // --shot-giant (#274): the nearest tree of each size class, side by side in turn (bench_tree_<class>.png), with
        // its trunk, chop hits (stone axe) and logs.
        if (benchGiant && gameState == 1 && !loading && WorldIdle()) {
            static int gf = 0, cls = 0;
            static Plant tree{};
            static bool have = false;
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            if (gf == 0) {
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                have = p->trees && p->trees->nearestPlant(me, 1500.0f, true, tree, false, cls);
                if (!have) { printf("bench: giant: no class %d tree within 1.5 km\n", cls); }
                else {
                    const float v = tree.chopSize();
                    printf("bench: giant: class %d (%s) %.0f m away: trunk %.1f m x r %.2f m, chop size x%.2f -> %d stone-axe hits, %d logs\n", cls,
                           WoodSpeciesInfo(SpeciesOfVariant(tree.variant)).name, Vector3Distance(me, tree.base), tree.height, tree.radius, v,
                           (int)lroundf(4.0f + 9.0f * tree.trunkVolume()), std::clamp((int)lroundf(3.0f * v), 1, 48));
                    player.flyMode = true;
                }
            }
            if (have) { // stand off by twice its height, at a third of it up, looking at its middle
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(tree.up, fabsf(tree.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const float d = std::max(8.0f, tree.height * 1.6f);
                const Vector3 eye = Vector3Add(tree.base, Vector3Add(Vector3Scale(t1, d), Vector3Scale(tree.up, tree.height * 0.35f)));
                const Vector3 at = Vector3Add(tree.base, Vector3Scale(tree.up, tree.height * 0.55f));
                world.camera.position = p->toWorld(eye); world.camera.target = p->toWorld(at); world.camera.up = p->dirToWorld(tree.up);
                player.position = world.renderOrigin + world.camera.position;
                player.forward = Vector3Normalize(p->dirToWorld(Vector3Subtract(at, eye)));
                sunOverride = true; sunArg = p->dirToWorld(Vector3Normalize(Vector3Add(tree.up, Vector3Scale(t1, 0.6f))));
            }
            if (++gf == 90 || !have) {
                if (have) { TakeScreenshot(TextFormat("bench_tree_%d.png", cls)); }
                gf = 0; have = false;
                if (++cls == kTreeClasses) { printf("bench: giant: done (bench_tree_0..3.png)\n"); benchGiant = false; }
            }
        }
        // --test-structures (#164): every set built in a row beside the start, its chests filled; then the nearest natural
        // site within 2 km, and a look at the first one (bench_structures.png).
        if (benchStructures && gameState == 1 && !loading && WorldIdle()) {
            static int sf = 0;
            static bool built = false;
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            ++sf;
            if (sf == 5) { player.flyMode = false; player.spawnOn(world, *p, StartSpawnDir(*p)); }
            if (sf == 20) {
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                built = Structures_Test(world, *p, me);
                uint64_t layout = 0;
                built = Structures_TestPlacement(*p, 2000.0f, layout) && built; // the world's placement round the start (#164)
                Vector3 at; std::string set;
                if (Structures_NearestSite(*p, me, 2000.0f, at, set)) printf("bench: structures: the nearest natural site is a %s, %.0f m from the start\n", set.c_str(), Vector3Distance(at, me));
                else printf("bench: structures: no natural site within 2 km\n");
            }
            if (sf == 80) Structures_TestFill(*p);
            static Vector3 shotEye{}, shotAt{}, shotUp{};
            // A tour of the sets built round the start (#164: each new set looked at): set i stands on side i % 4 (t1, -t1, t2,
            // -t2), 30 + 22 * (i / 4) m out (Structures_Test); 30 frames on each from 12 m back toward the start and 10 m up,
            // bench_structure_<i>.png on the last.
            static Vector3 home{};
            if (sf == 20) home = p->toLocal((player.position - world.renderOrigin).toFloat());
            const int nSets = Structures_SetCount(), tourEnd = 81 + 30 * nSets;
            if (sf >= 81 && sf < tourEnd) {
                const int i = (sf - 81) / 30;
                if ((sf - 81) % 30 == 0) {
                    shotUp = Vector3Normalize(home);
                    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(shotUp, fabsf(shotUp.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                    const Vector3 t2 = Vector3CrossProduct(t1, shotUp);
                    const Vector3 side = i % 4 == 0 ? t1 : i % 4 == 1 ? Vector3Negate(t1) : i % 4 == 2 ? t2 : Vector3Negate(t2);
                    shotAt = Vector3Add(home, Vector3Scale(side, 30.0f + 22.0f * (i / 4)));
                    shotEye = Vector3Add(shotAt, Vector3Add(Vector3Scale(side, -12.0f), Vector3Scale(shotUp, 10.0f)));
                    player.flyMode = true;
                }
                world.camera.position = p->toWorld(shotEye); world.camera.target = p->toWorld(shotAt); world.camera.up = p->dirToWorld(shotUp);
                player.position = world.renderOrigin + world.camera.position;
                if ((sf - 81) % 30 == 29) TakeScreenshot(TextFormat("bench_structure_%d.png", i));
            }
            if (sf == tourEnd) {
                TakeScreenshot("bench_structures.png");
                int items = 0;
                const int full = Structures_TestCheck(*p, items);
                const bool ok = built && full == Structures_SetCount() && items > 0;
                printf("bench: structures: %d of %d chests hold loot (%d items): %s\n", full, Structures_SetCount(), items, ok ? "PASS" : "FAIL");
                // Then the world's own placement: walk up to the nearest natural site (100 m short of it) and let it build.
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                static Vector3 siteAt{};
                std::string set;
                if (Structures_NearestSite(*p, me, 3000.0f, siteAt, set)) {
                    const Vector3 toward = Vector3Normalize(Vector3Subtract(me, siteAt));
                    player.flyMode = false;
                    { const auto* spf = dynamic_cast<const PlanetField*>(p->field.get()); const Vector3 d = Vector3Normalize(Vector3Add(siteAt, Vector3Scale(toward, 100.0f)));
                      player.placeAt(world, *p, p->toWorld(Vector3Scale(d, spf->radius + spf->elevation(d) + 1.5f))); } // exactly there (spawnOn looks for a safe spot)
                    printf("bench: structures: going to a %s 100 m away\n", set.c_str());
                    sf = 1000;
                } else benchStructures = false;
            }
            if (sf == 1200) { // ~3 s later: built (or skipped as unsuitable: water, a slope, a build nearby)?
                static Vector3 dummy{};
                Vector3 at; std::string set;
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                const bool found = Structures_NearestSite(*p, me, 400.0f, at, set);
                const bool built = found && p->plots && p->plots->plotFor(at, 25.0f, false);
                printf("bench: structures: the %s %.0f m from the player is %s (sites built so far %d; last left out: %s; fly %d)\n", found ? set.c_str() : "site",
                       found ? Vector3Distance(at, me) : -1.0f, built ? "built" : "not built", Structures_SitesBuilt(), Structures_LastSkip(), (int)player.flyMode);
                (void)dummy;
                benchStructures = false;
            }
        }
        // --test-meteor (#177): a forced shower at the (dark) dawn start: streaks in the sky, and a meteor brought down
        // 40 m off: a crater (the ground's middle now air) with a meteorite at its bottom. bench_meteor.png looks at it.
        if (benchMeteor && gameState == 1 && !loading && WorldIdle()) {
            static int mf = 0;
            static bool landed = false, crater = false, meteorRolled = false;
            static Vector3 ground{};
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
            ++mf;
            if (mf == 5) { player.flyMode = false; player.spawnOn(world, *p, StartSpawnDir(*p)); }
            if (mf == 20) Events_ForceShower();
            if (mf == 30) {
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                const Vector3 up = Vector3Normalize(me);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const Vector3 t2 = Vector3CrossProduct(t1, up);
                for (int k = 0; k < 8 && !landed; ++k) {
                    const float b = k * 0.785f;
                    const uint64_t rngBefore = GameRng_State();
                    if (Events_TestImpact(world, *p, 40.0f, b)) {
                        landed = true;
                        meteorRolled = GameRng_State() != rngBefore; // the meteorite's roll came from the authority's stream (#316)
                        const Vector3 d = Vector3Normalize(Vector3Add(up, Vector3Scale(Vector3Add(Vector3Scale(t1, cosf(b)), Vector3Scale(t2, sinf(b))), 40.0f / pf->radius)));
                        ground = Vector3Scale(d, pf->radius + pf->elevation(d));
                        crater = p->densityAt(Vector3Subtract(ground, Vector3Scale(d, 1.5f))) < ISO_LEVEL; // 1.5 m under the old ground: now air
                    }
                }
            }
            static Vector3 craterEye{}, craterAt{}, craterUp{};
            if (mf == 32 && landed) { // look at the crater from 16 m back and 9 m up
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                craterUp = Vector3Normalize(ground);
                const Vector3 back = Vector3Normalize(Vector3Subtract(me, ground));
                craterEye = Vector3Add(ground, Vector3Add(Vector3Scale(back, 16.0f), Vector3Scale(craterUp, 9.0f)));
                craterAt = Vector3Subtract(ground, Vector3Scale(craterUp, 2.0f));
                player.flyMode = true;
            }
            if (mf >= 32 && landed) { // the camera itself, each frame (the player's isn't applied in this bench)
                world.camera.position = p->toWorld(craterEye); world.camera.target = p->toWorld(craterAt); world.camera.up = p->dirToWorld(craterUp);
                player.position = world.renderOrigin + world.camera.position;
            }
            if (mf == 75) {
                TakeScreenshot("bench_meteor.png");
                int meteorites = 0;
                if (p->plots) for (auto& pl : p->plots->plots) for (const Part& part : pl->parts) meteorites += part.type == PartType::Meteorite;
                const bool ok = Events_ShowerOn() && Events_Streaks() > 0 && landed && crater && meteorites >= 1 && meteorRolled;
                printf("bench: meteor: shower on %d, %d streaks drawn, landed %d, crater %d, meteorites %d, game rng rolled %d (bench_meteor.png) -> %s\n",
                       (int)Events_ShowerOn(), Events_Streaks(), (int)landed, (int)crater, meteorites, (int)meteorRolled, ok ? "PASS" : "FAIL");
                benchMeteor = false;
            }
        }
        // --test-shower (#316): meteor showers are one world-wide event rolled by the authority. The local player flies 3 km
        // up (no longer a target, no sky to see); another player (an interest centre) stands 1.5 km away on the ground.
        //   1. crossing a 10-minute bucket rolls the gameplay RNG (no darkness or position needed);
        //   2. a client crossing a bucket rolls nothing and starts nothing; it applies the host's shower as sent;
        //   3. a forced shower starts although the local player is in the air, and both meteorites land at the centre;
        //   4. the shower is saved and loaded with the world.
        if (benchShower && gameState == 1 && !loading && WorldIdle()) {
            static int sf = 0, phase = 0;
            static bool rolled = false, clientQuiet = false, clientApplied = false, started = false, landed = false, kept = false;
            static uint32_t count0 = 0;
            static uint64_t rng0 = 0;
            static Vector3 centre{};
            static EventsShower host{};
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
            ++sf;
            auto fly = [&] { // the local player 3 km over the spawn, flying
                const Vector3 d = p->dirToLocal(StartSpawnDir(*p));
                player.flyMode = true;
                player.position = world.renderOrigin + p->toWorld(Vector3Scale(d, pf->radius + 3000.0f));
                player.velocity = { 0, 0, 0 };
            };
            if (phase == 0 && sf == 5) {
                fly();
                const Vector3 d0 = p->dirToLocal(StartSpawnDir(*p));
                const Vector3 side = Vector3Normalize(Vector3CrossProduct(d0, fabsf(d0.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                Vector3 d = d0;
                for (int k = 0; k < 16; ++k) { // 1.5 km away on dry land
                    d = Vector3Normalize(Vector3Add(d0, Vector3Scale(Vector3RotateByAxisAngle(side, d0, k * 0.39f), 1500.0f / pf->radius)));
                    if (pf->elevation(d) > 5.0f) break;
                }
                centre = Vector3Scale(d, pf->radius + pf->elevation(d) + 1.0f);
                gInterestCentres = { InterestCentre{ p->id, centre } };
                rng0 = GameRng_State(); count0 = Events_ShowerCount();
                SetWorldTime((floor(WorldTime() / 600.0) + 1.0) * 600.0 + 1.0); // into the next bucket
                phase = 1; sf = 0;
            } else if (phase == 1 && sf == 3) {
                fly();
                rolled = GameRng_State() != rng0;
                // As a client: another bucket rolls nothing; the host's shower applies as sent.
                gNetRole = NetRole::Client;
                rng0 = GameRng_State(); count0 = Events_ShowerCount();
                SetWorldTime(WorldTime() + 600.0);
                phase = 2; sf = 0;
            } else if (phase == 2 && sf == 3) {
                clientQuiet = GameRng_State() == rng0 && Events_ShowerCount() == count0;
                Events_ApplyShower({ WorldTime() - 1.0, WorldTime() + 119.0, count0 + 7 });
                clientApplied = Events_ShowerOn() && Events_ShowerCount() == count0 + 7;
                Events_ApplyShower({ -1.0, -1.0, count0 + 7 }); // over
                gNetRole = NetRole::Offline;
                count0 = Events_ShowerCount();
                Events_ForceShower();
                phase = 3; sf = 0;
            } else if (phase == 3) {
                fly();
                if (sf == 3) started = Events_ShowerOn() && Events_ShowerCount() == count0 + 1;
                const EventsShower s = Events_GetShower();
                if (sf == 4) SetWorldTime(s.start + 0.40 * (s.end - s.start));
                if (sf == 8) SetWorldTime(s.start + 0.80 * (s.end - s.start));
                if (sf == 12) {
                    int atPlayer = 0, atCentres = 0;
                    Events_ImpactTargets(atPlayer, atCentres);
                    int near = 0;
                    if (p->plots) for (auto& pl : p->plots->plots) for (const Part& part : pl->parts)
                        if (part.type == PartType::Meteorite && Vector3Distance(pl->toBody({ (part.origin.x + 0.5f) * kCellSize, 0.0f, (part.origin.z + 0.5f) * kCellSize }), centre) < 260.0f) ++near;
                    landed = atCentres == 2 && atPlayer == 0 && near >= 1;
                    printf("bench: shower: meteorites at the centre %d, at the local player %d, meteorite parts near the centre %d (shower %.0f..%.0f s, now %.0f s, started %d)\n", atCentres, atPlayer, near, s.start, s.end, WorldTime(), (int)started);
                    host = s;
                    const bool saved = SaveWorld(world);
                    Events_Reset(); // forget it: the load must bring it back
                    const bool loaded = saved && LoadWorldSave(world);
                    kept = loaded;
                    phase = 4; sf = 0;
                }
            } else if (phase == 4 && sf == 3) {
                const EventsShower s = Events_GetShower();
                kept = kept && s.id == host.id && fabs(s.start - host.start) < 0.01 && fabs(s.end - host.end) < 0.01;
                gInterestCentres.clear();
                count0 = Events_ShowerCount();
                phase = 5; sf = 0;
            } else if (phase == 5) {
                // 5. the rate: one bucket per frame for kRateBuckets buckets; the share that started a shower should be near
                //    kMeteorShowerOdds (5%, user 2026-11-15). Each jump also ends the previous shower before its meteorites.
                constexpr int kRateBuckets = 600;
                SetWorldTime((floor(WorldTime() / kMeteorShowerBucketSeconds) + 1.0) * kMeteorShowerBucketSeconds + 1.0);
                if (sf > kRateBuckets) {
                    const int showers = (int)(Events_ShowerCount() - count0);
                    const double expect = kMeteorShowerOdds * kRateBuckets, sigma = sqrt(expect * (1.0 - kMeteorShowerOdds));
                    const bool rateOk = fabs(showers - expect) < 4.0 * sigma;
                    printf("bench: shower: rate %d showers in %d buckets (%.1f%%, expected %.0f%% = %.0f +- %.0f)\n", showers, kRateBuckets,
                           100.0 * showers / kRateBuckets, 100.0 * kMeteorShowerOdds, expect, 4.0 * sigma);
                    const bool ok = rolled && clientQuiet && clientApplied && started && landed && kept && rateOk;
                    printf("bench: shower: bucket rolled the game RNG %d, client rolled nothing %d, client applied the host's shower %d, forced shower with the player in the air %d, both meteorites at the other player %d, kept over save/load %d, rate %d -> %s\n",
                           (int)rolled, (int)clientQuiet, (int)clientApplied, (int)started, (int)landed, (int)kept, (int)rateOk, ok ? "PASS" : "FAIL");
                    benchShower = false;
                }
            }
        }
        // --test-host (#311, with --host <save>): the dedicated server loaded or created its save and runs headless.
        //   - seed (user, 2026-11-15): a new world takes --seed N, else a random seed (printed); an existing save keeps its own;
        //   - autosave: with a short --autosave-minutes (<= 1) the test waits for one autosave;
        //   - a client opened in this process over loopback sends "join" and gets the world description (fingerprint and
        //     seed must match the host's world) and the current meteor shower. That part SKIPs without ENABLE_GNS.
        if (benchHost && gameState == 1 && !loading && WorldIdle()) {
            static uint32_t client = 0;
            static double t0 = -1.0;
            static bool sent = false, gotJoin = false, gotShower = false, sameWorld = false, sameShower = false;
            if (t0 < 0.0) t0 = GetTime();
            bool saveOnDisk = false;
            if (FILE* fsv = hostSave.empty() ? nullptr : fopen(hostSave.c_str(), "rb")) { fclose(fsv); saveOnDisk = true; }
            const unsigned expectSeed = hostNewWorld ? worldGen.seed : hostSaveSeed; // given or rolled up front / the save's
            const bool seedOk = expectSeed != 0 && worldGen.lastSeed == expectSeed;
            const char* seedHow = !hostNewWorld ? "the save's own" : seedGiven ? "--seed" : "random";
            const bool wantAutosave = autosaveMinutes > 0.0 && autosaveMinutes <= 1.0;
            const bool autosaveDone = !wantAutosave || Host_Autosaves() >= 1;
            if (hostSave.empty()) { printf("bench: host: --test-host needs --host <save> -> FAIL\n"); benchHost = false; }
            else if (!Host_Listening()) {
                if (autosaveDone || GetTime() - t0 > 90.0) {
                    const bool ok = saveOnDisk && gHeadless && Net_IsAuthority() && seedOk && autosaveDone;
                    printf("bench: host: save %s on disk %d, headless %d, authority %d, seed %u (%s) ok %d, autosaves %d; not listening (%s) -> %s\n",
                           hostSave.c_str(), (int)saveOnDisk, (int)gHeadless, (int)Net_IsAuthority(), worldGen.lastSeed, seedHow, (int)seedOk, Host_Autosaves(),
                           Net_Available() ? "network failed" : "built without ENABLE_GNS", ok ? "SKIP (network part)" : "FAIL");
                    benchHost = false;
                }
            } else {
                static unsigned joinSeed = 0;
                if (!client) { client = Net_Connect(TextFormat("127.0.0.1:%d", hostPort)); Host_SetLocalClient(client); }
                if (client && !sent && Net_IsConnected(client)) sent = Net_Send(client, "join", 4, true);
                for (const NetMessage& m : Host_TakeLocalClientMessages()) {
                    const std::string text(m.data.begin(), m.data.end());
                    if (!text.empty() && text[0] == 'J') {
                        JoinInfo info; std::string why;
                        gotJoin = JoinInfo_Read(text.substr(1), info, why);
                        sameWorld = gotJoin && info.fingerprint == WorldFingerprint(world);
                        joinSeed = gotJoin ? info.seed : 0;
                        if (!gotJoin) printf("bench: host: join info unreadable: %s\n", why.c_str());
                    } else if (!text.empty() && text[0] == 'S') {
                        double s = 0.0, e = 0.0; unsigned id = 0;
                        gotShower = sscanf(text.c_str() + 1, " %lf %lf %u", &s, &e, &id) == 3;
                        sameShower = gotShower && id == Events_GetShower().id;
                    }
                }
                if ((gotJoin && gotShower && autosaveDone) || GetTime() - t0 > 90.0) {
                    const bool ok = saveOnDisk && gHeadless && sent && sameWorld && sameShower && Host_Joins() == 1 && seedOk && joinSeed == worldGen.lastSeed && autosaveDone;
                    printf("bench: host: save %s on disk %d, headless %d, listening on %d; seed %u (%s) ok %d, sent to the client %d; autosaves %d; loopback client joined %d, world description %d (fingerprint matches %d), shower %d (same %d), joins answered %d after %.1f s -> %s\n",
                           hostSave.c_str(), (int)saveOnDisk, (int)gHeadless, hostPort, worldGen.lastSeed, seedHow, (int)seedOk, (int)(joinSeed == worldGen.lastSeed),
                           Host_Autosaves(), (int)sent, (int)gotJoin, (int)sameWorld, (int)gotShower, (int)sameShower, Host_Joins(), GetTime() - t0, ok ? "PASS" : "FAIL");
                    benchHost = false;
                }
            }
        }
        // --test-headless (#311, with --bench 12): the world runs with nothing drawn - the sim steps at its rate, the player
        // stands on streamed ground (collision uses the CPU field), animals live; drawScene is skipped every frame.
        if (benchHeadless && gameState == 1 && !loading && WorldIdle()) {
            static double t0 = -1.0;
            static long long tick0 = 0, draws0 = 0;
            static DVec3 pos0{}; // body-local (the planet spins and orbits)
            Planetoid* hb = dynamic_cast<Planetoid*>(world.objects[0].get());
            const DVec3 local = hb ? hb->toLocalD(player.position) : player.position;
            static double settle = -1.0; // stand the player on the spawn first (the bench starts it 40 m up, flying)
            if (settle < 0.0) { settle = GetTime(); player.flyMode = false; if (hb) player.spawnOn(world, *hb, StartSpawnDir(*hb)); }
            else if (GetTime() - settle < 2.0) {}
            else if (t0 < 0.0) { t0 = GetTime(); tick0 = SimTick(); draws0 = gHeadlessSkippedDraws; pos0 = local; }
            else if (GetTime() - t0 > 8.0) {
                const double secs = GetTime() - t0, simSecs = (double)(SimTick() - tick0) / kSimHz;
                const long long skipped = gHeadlessSkippedDraws - draws0;
                const float moved = (float)(local - pos0).length();
                const bool ok = gHeadless && simSecs > 0.8 * secs && skipped > 0 && survival.deaths == 0 && !survival.dead && moved < 2.0f;
                printf("bench: headless: %.1f s of sim in %.1f s, %lld scene draws skipped, grounded %d, moved %.2f m, %d animals, deaths %d -> %s\n",
                       simSecs, secs, skipped, (int)player.grounded, moved, Animals_Count(), survival.deaths, ok ? "PASS" : "FAIL");
                benchHeadless = false;
            }
        }
        // --test-trapping (#176): a snare set down from the hand at the ground ahead; aged 2 h of world time it holds a
        // catch, which the player (standing beside it) collects; then save/load keeps a fresh snare's state.
        if (benchCombat && gameState == 1 && !loading && WorldIdle() && Combat_TestFrame(world)) benchCombat = false;
        if (benchTrapping && gameState == 1 && !loading && WorldIdle()) {
            static int tf = 0;
            static int meat0 = 0, placed = 0;
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            ++tf;
            if (tf == 5) { player.flyMode = false; player.spawnOn(world, *p, StartSpawnDir(*p)); }
            if (tf >= 8 && tf < 30) player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))), 0.8f), Vector3Scale(player.up, 0.6f)));
            if (tf == 30) { inventory.add("snare", 2); Hotbar_Bind(6, "snare"); Hotbar_Press(6, false); }
            if (tf == 35) { const int before = inventory.count("snare"); placed = Hotbar_PlaceTorch() && inventory.count("snare") == before - 1; meat0 = inventory.count("raw_meat"); }
            if (tf == 40) Trapping_TestAge(world, 7200.0);
            if (tf == 50) Trapping_Update(world);
            if (tf == 60) {
                const int got = inventory.count("raw_meat") - meat0;
                const bool ok = placed && got >= 1;
                printf("bench: trapping: snare placed %d; aged 2 h, it caught and the player collected %d raw meat (hide %d) -> %s\n", placed, got, inventory.count("hide"), ok ? "PASS" : "FAIL");
                benchTrapping = false;
            }
        }
        // --test-fishing (#176): a cast into the sea near the start lands a cod; reeling in early lands nothing.
        if (benchFishing && gameState == 1 && !loading && WorldIdle()) {
            static int ff = 0;
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
            if (++ff == 5) player.spawnOn(world, *p, StartSpawnDir(*p));
            if (ff == 20) {
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                const Vector3 up = Vector3Normalize(me);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const Vector3 t2 = Vector3CrossProduct(t1, up);
                Vector3 sea{}; bool found = false;
                for (float rr = 20.0f; rr < 1500.0f && !found; rr += 20.0f)
                    for (int k = 0; k < 16 && !found; ++k) {
                        const float a = k * 0.3927f;
                        const Vector3 d = Vector3Normalize(Vector3Add(up, Vector3Scale(Vector3Add(Vector3Scale(t1, cosf(a)), Vector3Scale(t2, sinf(a))), rr / pf->radius)));
                        if (pf->elevation(d) < -2.0f) { sea = Vector3Scale(d, pf->radius - 0.3f); found = true; }
                    }
                const std::string caught = found ? Fishing_TestCast(world, *p, sea, false) : std::string();
                const std::string early = found ? Fishing_TestCast(world, *p, sea, true) : std::string("?");
                const std::string land = Fishing_CatchAt(*p, me);
                const bool ok = found && caught == "raw_cod" && early.empty() && inventory.count("raw_cod") >= 1 && (land == "raw_trout" || land == "raw_catfish");
                printf("bench: fishing: sea %s; a full cast caught \"%s\", reeling in early \"%s\"; inland here it would be %s -> %s\n",
                       found ? "found" : "NOT found", caught.c_str(), early.c_str(), land.c_str(), ok ? "PASS" : "FAIL");
                benchFishing = false;
            }
        }
        // --test-fill (#228): standing at the start looking at the ground 3 m ahead, an exact sphere of dirt.
        if (benchFill && gameState == 1 && !loading && WorldIdle()) {
            static int ff = 0;
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            ++ff;
            if (ff == 5) { player.flyMode = false; player.spawnOn(world, *p, StartSpawnDir(*p)); }
            if (ff >= 8 && ff < 30) player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))), 0.8f), Vector3Scale(player.up, 0.5f)));
            if (ff == 30) TerrainTool_TestFill();
            if (ff == 70) { TakeScreenshot("bench_fill.png"); printf("bench: fill: bench_fill.png\n"); benchFill = false; }
        }
        // --test-lamps (#139): at the dawn start (dark), 100 lamp cubes on a 2 m grid round the player light the ground.
        if (benchLamps && gameState == 1 && !loading && WorldIdle()) {
            static int lf = 0;
            static double lumOff = 0.0;
            static int placed = 0;
            auto groundLum = [](const char* file) {
                Image im = LoadImage(file);
                ImageFormat(&im, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
                const Color* px = (const Color*)im.data;
                double sum = 0.0; int n = 0;
                for (int y = im.height * 5 / 10; y < im.height * 9 / 10; y += 2)
                    for (int x = im.width * 2 / 10; x < im.width * 8 / 10; x += 2) { const Color c = px[y * im.width + x]; sum += 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b; ++n; }
                UnloadImage(im);
                return n ? sum / n : 0.0;
            };
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            ++lf;
            if (lf == 5) { player.flyMode = false; player.spawnOn(world, *p, StartSpawnDir(*p)); }
            if (lf >= 8 && lf < 70) player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))), 0.8f), Vector3Scale(player.up, 0.6f)));
            if (lf == 25) {
                TakeScreenshot("bench_lamps_off.png"); lumOff = groundLum("bench_lamps_off.png");
                const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
                const Vector3 me = p->toLocal((player.position - world.renderOrigin).toFloat());
                const Vector3 up = Vector3Normalize(me);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const Vector3 t2 = Vector3CrossProduct(t1, up);
                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                const int lampMat = BuildMaterialForItem("lamp");
                Plot* pl = p->plots->plotFor(me, 0.0f, true);
                for (int i = -5; i < 5; ++i)
                    for (int j = -5; j < 5; ++j) {
                        if (abs(2 * i + 1) <= 1 && abs(2 * j + 1) <= 1) continue; // not on the player
                        const Vector3 d = Vector3Normalize(Vector3Add(me, Vector3Add(Vector3Scale(t1, 2.0f * i + 1.0f), Vector3Scale(t2, 2.0f * j + 1.0f))));
                        const Vector3 at = Vector3Scale(d, pf->radius + pf->elevation(d) + 0.3f);
                        placed += pl->place(Part{ PartType::Cube, pl->cellAt(pl->toPlot(at)), 0, (uint8_t)std::max(0, lampMat) }) >= 0 ? 1 : 0;
                    }
            }
            if (lf == 60) {
                TakeScreenshot("bench_lamps.png");
                const double lumOn = groundLum("bench_lamps.png");
                const bool ok = placed >= 90 && LocalLights_Count() >= 90 && lumOn > lumOff + 6.0;
                printf("bench: lamps: %d lamp cubes placed, %d local lights uploaded (busiest tile %d), ground luminance %.1f -> %.1f: %s\n",
                       placed, LocalLights_Count(), LocalLights_MaxPerTile(), lumOff, lumOn, ok ? "PASS" : "FAIL");
                benchLamps = false;
            }
        }
        // --test-torch (#218): at the dawn start (still dark), the ground in front gets brighter with a torch in hand.
        // --test-posture (#578): crouch and prone shorten the capsule and slow the walk; getting up needs room (a roof 1 m
        // up keeps you prone until it is taken away). Distances walked in 1 s against standing: crouched ~1/2, prone ~1/4.
        if (benchPosture && gameState == 1 && !loading && (postureStarted || WorldIdle())) {
            postureStarted = true;
            static int phase = 0; static double t0 = 0.0; static DVec3 from{}; /* body-local */ static float dStand = 0.0f, dCrouch = 0.0f, dProne = 0.0f;
            static bool ok = true; static std::vector<int> roof; static uint32_t roofPlot = 0;
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            auto check = [&](bool c, const char* what) { printf("bench: posture: %s: %s\n", what, c ? "ok" : "NO"); ok = ok && c; };
            auto walked = [&]() { return (float)(p->toLocalD(player.position) - from).length() * p->scale; };
            const double t = GetTime() - t0;
            if (phase == 0) { player.flyMode = false; player.spawnOn(world, *p, StartSpawnDir(*p)); phase = 1; t0 = GetTime(); }
            else if (phase == 1 && t > 1.5 && player.grounded) { from = p->toLocalD(player.position); player.autoWalk = true; phase = 2; t0 = GetTime(); }
            else if (phase == 2 && t > 1.0) { dStand = walked(); player.autoWalk = false; player.autoCrouch = true; phase = 3; t0 = GetTime(); }
            else if (phase == 3 && t > 0.5) {
                check(player.posture == Player::Posture::Crouch && fabsf(player.height - Player::kCrouchHeight) < 1e-3f && fabsf(player.eyeHeight - Player::kCrouchEye) < 0.05f, "Ctrl held: crouched, 1.25 m, eye 1.15 m");
                from = p->toLocalD(player.position); player.autoWalk = true; phase = 4; t0 = GetTime();
            } else if (phase == 4 && t > 1.0) {
                dCrouch = walked(); player.autoWalk = false; player.autoCrouch = false; player.proneWanted = true; phase = 5; t0 = GetTime();
            } else if (phase == 5 && t > 0.5) {
                check(player.posture == Player::Posture::Prone && fabsf(player.height - Player::kProneHeight) < 1e-3f, "Z: prone, 0.9 m");
                from = p->toLocalD(player.position); player.autoWalk = true; phase = 6; t0 = GetTime();
            } else if (phase == 6 && t > 1.0) {
                dProne = walked(); player.autoWalk = false;
                printf("bench: posture: 1 s walked: standing %.2f m, crouched %.2f m, prone %.2f m\n", dStand, dCrouch, dProne);
                check(dStand > 3.0f && fabsf(dCrouch / dStand - 0.5f) < 0.15f && fabsf(dProne / dStand - 0.25f) < 0.1f, "crouched ~1/2 and prone ~1/4 of the walk");
                // A roof 1.0-1.5 m over the feet, 2.5 m square.
                if (!p->plots) p->plots = std::make_unique<PlotSet>();
                const Vector3 feet = p->toLocal((player.position - world.renderOrigin).toFloat());
                Plot* pl = p->plots->plotFor(feet, 0.0f, true);
                roofPlot = pl->id;
                const Int3 c = pl->cellAt(pl->toPlot(Vector3Add(feet, Vector3Scale(Vector3Normalize(feet), 0.1f / p->scale))));
                for (int z = -2; z <= 2; ++z) for (int x = -2; x <= 2; ++x) { const int i = pl->place({ PartType::Cube, { c.x + x, c.y + 2, c.z + z }, 0 }); if (i >= 0) roof.push_back(i); }
                player.proneWanted = false; // Z again: get up
                phase = 7; t0 = GetTime();
            } else if (phase == 7 && t > 0.6) {
                check(!roof.empty() && player.posture == Player::Posture::Prone, TextFormat("under a roof 1 m up (%d cubes): still prone", (int)roof.size()));
                if (Plot* pl = p->plots->byId(roofPlot)) for (int i : roof) pl->remove(i);
                phase = 8; t0 = GetTime();
            } else if (phase == 8 && t > 0.6) {
                check(player.posture == Player::Posture::Stand && fabsf(player.height - player.standHeight) < 1e-3f, "roof gone: standing again");
                printf("bench: posture: %s\n", ok ? "PASS" : "FAIL");
                benchPosture = false;
            }
        }
        // --test-f2 (#579): the F2 screenshot path once the world is up (the recorder prints "screenshot: <dir> saved").
        if (benchF2 && gameState == 1 && !loading && WorldIdle()) { benchF2 = false; Recorder_Screenshot(); }
        // --shot-dig (#533): standing on the ground, looking 40 deg down: the dig preview's ring (bench_dig.png).
        if (benchShotDig && gameState == 1 && !loading && WorldIdle()) {
            static int df = 0;
            ++df;
            if (df == 5) if (Planetoid* p0 = dynamic_cast<Planetoid*>(world.objects[0].get())) { player.flyMode = false; player.spawnOn(world, *p0, StartSpawnDir(*p0)); }
            if (df >= 8) {
                const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
                player.forward = Vector3Normalize(Vector3Subtract(flat, Vector3Scale(player.up, 0.84f))); // tan 40 deg
            }
            if (df == 40) { TakeScreenshot("bench_dig.png"); printf("bench: shot-dig: bench_dig.png\n"); benchShotDig = false; }
        }
        // --test-bite (#165): a thornback set down 2 m ahead charges, warns, and strikes no sooner than its wind-up; the hit
        // goes through Combat_Damage (8 health on normal) and shoves the player away from it.
        if (benchBite && gameState == 1 && !loading && (biteStarted || WorldIdle())) {
            biteStarted = true;
            static int phase = 0; static double t0 = 0.0, tHit = -1.0; static DVec3 atHit{}; // body-local
            static float hp0 = 100.0f; static double sim0 = 0.0; // the spawn, in sim time (the wind-up runs on it)
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const double t = GetTime() - t0;
            if (phase == 0) { player.flyMode = false; player.spawnOn(world, *p, StartSpawnDir(*p)); phase = 1; t0 = GetTime(); }
            else if (phase == 1 && t > 1.5 && player.grounded) {
                Animals_SetThreats(true);
                survival.health = 100.0f; hp0 = survival.health;
                const Vector3 feet = p->toLocal((player.position - world.renderOrigin).toFloat());
                const Vector3 flat = Vector3Normalize(p->dirToLocal(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))));
                const int n = Animals_Spawn(world, *p, "thornback", Vector3Add(feet, Vector3Scale(flat, 2.0f / p->scale)), 1, false, 0.0f);
                printf("bench: bite: %d thornback set down 2 m ahead\n", n);
                phase = 2; t0 = GetTime();
            } else if (phase == 2) {
                if (tHit < 0.0 && survival.health < hp0 - 0.01f) { tHit = SimTime() - sim0; atHit = p->toLocalD(player.position); }
                if (tHit >= 0.0 && SimTime() - sim0 > tHit + 0.4) {
                    const float lost = hp0 - survival.health;
                    const float shoved = (float)(p->toLocalD(player.position) - atHit).length() * p->scale;
                    const bool ok = tHit >= Animals_BiteWindUp() && fabsf(lost - 8.0f * Difficulty_ThreatScale()) < 0.5f && shoved > 0.3f;
                    printf("bench: bite: first hit %.2f s after it was set down (wind-up %.1f s), %.1f health lost, shoved %.2f m -> %s\n",
                           tHit, Animals_BiteWindUp(), lost, shoved, ok ? "PASS" : "FAIL");
                    Animals_SetThreats(false); benchBite = false;
                } else if (t > 10.0) { printf("bench: bite: no hit within 10 s -> FAIL\n"); Animals_SetThreats(false); benchBite = false; }
            }
        }
        if (benchTorch && gameState == 1 && !loading && WorldIdle()) {
            static int tf = 0;
            static double lumOff = 0.0;
            static bool torchPlacedOk = false;
            auto groundLum = [](const char* file) {
                Image im = LoadImage(file);
                ImageFormat(&im, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
                const Color* px = (const Color*)im.data;
                double sum = 0.0; int n = 0;
                for (int y = im.height * 6 / 10; y < im.height * 9 / 10; y += 2)
                    for (int x = im.width * 3 / 10; x < im.width * 7 / 10; x += 2) { const Color c = px[y * im.width + x]; sum += 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b; ++n; }
                UnloadImage(im);
                return n ? sum / n : 0.0;
            };
            ++tf;
            if (tf == 5) { if (Planetoid* p0 = dynamic_cast<Planetoid*>(world.objects[0].get())) { player.flyMode = false; player.spawnOn(world, *p0, StartSpawnDir(*p0)); } inventory.add("torch", 3); Hotbar_Bind(7, "torch"); Hotbar_Press(8, false); }
            if (tf >= 8 && tf < 40) player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up)))), 0.8f), Vector3Scale(player.up, 0.6f))); // look down at the ground in front
            if (tf == 20) { TakeScreenshot("bench_torch_off.png"); lumOff = groundLum("bench_torch_off.png"); Hotbar_Press(7, false); }
            if (tf == 40) {
                TakeScreenshot("bench_torch.png");
                const double lumOn = groundLum("bench_torch.png");
                const bool lit = world.lights.size() > 3 && world.lights[3].enabled;
                const bool ok = lit && lumOn > lumOff + 4.0;
                printf("bench: torch: selected %d \"%s\", torches %d, fly %d, piloted %d\n", Hotbar_Selected(), Hotbar_Selected() >= 0 ? Hotbar_Slot(Hotbar_Selected()).c_str() : "", inventory.count("torch"), (int)player.flyMode, (int)ship.piloted);
                printf("bench: torch: light %d, ground luminance %.1f -> %.1f with a torch in hand: %s\n", (int)lit, lumOff, lumOn, ok ? "PASS" : "FAIL");
                // Then set one down where the crosshair meets the ground (#218), and put the hand away: it lights the ground.
                const int before = inventory.count("torch");
                const bool placedOk = Hotbar_PlaceTorch();
                torchPlacedOk = placedOk && inventory.count("torch") == before - 1;
                Hotbar_Press(8, false);
            }
            if (tf == 70) {
                printf("bench: torch: %d particles alive (the placed torch's flame, #558)\n", Particles_Count());
                TakeScreenshot("bench_torch_placed.png");
                const double lumPlaced = groundLum("bench_torch_placed.png");
                const bool handOff = !(world.lights.size() > 3 && world.lights[3].enabled);
                const bool ok = torchPlacedOk && handOff && LocalLights_Count() >= 1 && lumPlaced > lumOff + 4.0;
                printf("bench: torch: placed one %d (light in hand now %d), %d local light(s), ground luminance %.1f -> %.1f from the placed torch: %s\n",
                       (int)torchPlacedOk, (int)!handOff, LocalLights_Count(), lumOff, lumPlaced, ok ? "PASS" : "FAIL");
                benchTorch = false;
            }
        }
        // --test-pause (#511): Esc into the pause menu, then two screen grabs 3 s apart: no pixel should move (grass, trees,
        // water, clouds, particles), and neither the animation clock nor the simulation may advance.
        if (benchPause && !loading && WorldIdle() && (gameState == 1 || gameState == 2)) {
            static int pf = 0;
            static double t0 = -1.0, clock0 = 0.0, sim0 = 0.0;
            static Image first{};
            ++pf;
            if (pf == 30) menuEsc = true; // as the key going down
            if (pf == 60) { first = LoadImageFromScreen(); clock0 = GameClock(); sim0 = SimTime(); t0 = GetTime(); }
            if (t0 >= 0.0 && GetTime() - t0 > 3.0) {
                Image now = LoadImageFromScreen();
                long moved = 0, total = 0;
                if (now.width == first.width && now.height == first.height) {
                    const Color* a = (const Color*)first.data; const Color* b = (const Color*)now.data;
                    total = (long)now.width * now.height;
                    for (long i = 0; i < total; ++i)
                        if (abs(a[i].r - b[i].r) > 8 || abs(a[i].g - b[i].g) > 8 || abs(a[i].b - b[i].b) > 8) ++moved;
                }
                UnloadImage(first); UnloadImage(now);
                const double pct = total ? 100.0 * moved / total : 100.0;
                const bool ok = gameState == 2 && GameClock() == clock0 && SimTime() == sim0 && pct < 0.02; // (grass swaying in the bench view: ~0.07%)
                printf("bench: pause: paused %d; over %.1f s the animation clock moved %.3f s, the simulation %.3f s, %.3f%% of pixels changed: %s\n",
                       (int)(gameState == 2), GetTime() - t0, GameClock() - clock0, SimTime() - sim0, pct, ok ? "PASS" : "FAIL");
                benchPause = false;
            }
        }
        // --test-time (#259, #275): a new game's dawn start (the sun ~4 deg under the spawn's horizon, rising), then a +6 h
        // jump through the debug time skip: the sun higher, and the player carried along (same body-local spot, grounded).
        if (benchTime && gameState == 1 && !loading && WorldIdle()) {
            ++timeFrames;
            Planetoid* body = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get());
            auto feetLocal = [&] { return body->toLocal((player.position - world.renderOrigin).toFloat()); };
            if (!body) { printf("bench: time: no body -> FAIL\n"); benchTime = false; }
            else if (timeStep == 0 && timeFrames == 5) { StartNewWorld(world, false, true); timeStep = 1; timeFrames = 0; }
            else if (timeStep == 1 && timeFrames == 60) {
                timeElev0 = SunElevationAtCamera(world);
                const double rise = SecondsToSunElevationAtCamera(world, 0.0f, true), hour = BodyHourSeconds(world);
                const bool ok = timeElev0 > -6.0f && timeElev0 < -2.0f && rise > 0.0 && rise < 2.0 * hour;
                timeOk &= ok;
                printf("bench: time: dawn start: sun %.2f deg, sunrise in %.0f s (%.2f h): %s\n", timeElev0, rise, rise / std::max(hour, 1e-9), ok ? "yes" : "NO");
                timeLocal0 = feetLocal();
                gTimeSkipSeconds = 6.0 * hour;
                timeStep = 2; timeFrames = 0;
            } else if (timeStep == 2 && timeFrames == 90) {
                const float e = SunElevationAtCamera(world);
                const float moved = Vector3Distance(feetLocal(), timeLocal0);
                const bool ok = e > timeElev0 + 20.0f && moved < 1.0f && player.grounded;
                timeOk &= ok;
                printf("bench: time: +6 h: sun %.2f -> %.2f deg, player moved %.2f m on the body, grounded %s: %s -> %s\n", timeElev0, e, moved,
                       player.grounded ? "yes" : "no", ok ? "yes" : "NO", timeOk ? "PASS" : "FAIL");
                TakeScreenshot("bench_time.png");
                benchTime = false;
            }
        }
        // --test-menus: the menu flow without clicks: main menu (screenshot), New Game, pause
        // (screenshot), Save & Quit to Menu, Load Game; each step's state is checked.
        if (benchMenus) {
            ++menuFrames;
            if (menuStep == 0) { gameState = 0; mainMenu.isActive = true; world.isActive = false; menuStep = 1; menuFrames = 0; }
            else if (menuStep == 1 && menuFrames == 30) {
                // #503: the title screen comes first, with no world: no bodies (no erosion bake ran), no world jobs.
                const bool noWorld = world.objects.empty() && WorldJobs().idle() && TerrainArena_Stats().vertexCapacity == 0;
                menuOk &= noWorld;
                printf("bench: menus: title screen with no world (no bodies, no jobs, no terrain arena, #503): %s\n", noWorld ? "yes" : "NO");
                TakeScreenshot("bench_menu_main.png"); pendingMenu.action = MenuAction::NewGame; menuStep = 2; menuFrames = 0;
            }
            else if (menuStep == 2 && menuFrames > 60 && WorldIdle() && !loading) { // (New Game spawns the player: wait out its loading screen)
                menuOk &= gameState == 1;
                printf("bench: menus: New Game -> in game: %s\n", gameState == 1 ? "yes" : "NO");
                menuEsc = true; menuStep = 3; menuFrames = 0; // pause through the in-game Esc check, as a key press does
            } else if (menuStep == 3 && menuFrames == 30) {
                menuOk &= gameState == 2;
                printf("bench: menus: Esc in play -> pause menu: %s\n", gameState == 2 ? "yes" : "NO");
                TakeScreenshot("bench_menu_pause.png"); pendingMenu.action = MenuAction::SaveAndQuit; menuStep = 4; menuFrames = 0;
            }
            else if (menuStep == 4 && menuFrames == 5) {
                menuOk &= gameState == 0;
                printf("bench: menus: Save & Quit -> main menu: %s; old world unloaded (%zu bodies left, #270): %s\n", gameState == 0 ? "yes" : "NO", world.objects.size(), world.objects.empty() ? "yes" : "NO");
                menuOk &= world.objects.empty();
                // #270: nothing of it runs or stays allocated behind the menu.
                const bool down = WorldJobs().idle() && TerrainGpu_Pending() == 0 && TerrainArena_Stats().vertexCapacity == 0;
                menuOk &= down;
                printf("bench: menus: behind the menu: %zu jobs queued, %d running, %d GPU batches, arena %d vertices (#270): %s\n",
                       WorldJobs().queued(), WorldJobs().running(), TerrainGpu_Pending(), TerrainArena_Stats().vertexCapacity, down ? "torn down" : "NOT torn down");
                menuBakesBeforeLoad = ErosionBakeCount();
                pendingMenu.action = MenuAction::LoadGame; pendingMenu.path = "saves/bench.sav"; menuStep = 5; menuFrames = 0;
            } else if (menuStep == 5 && menuFrames == 5) {
                menuOk &= gameState == 1;
                printf("bench: menus: Load Game -> in game: %s\n", gameState == 1 ? "yes" : "NO");
                // #505: a save with its own erosion cache loads it; it never bakes again.
                const bool cached = !gErosionPrototype || ErosionBakeCount() == menuBakesBeforeLoad;
                menuOk &= cached;
                printf("bench: menus: the load read the save's erosion cache, no bake (#505): %s\n", cached ? "yes" : "NO");
                // #388: the loading screen as entering a world raises it (benches load without one), then Esc on it.
                loading = true; spawnLoad = false; loadStart = GetTime(); loadFrames = 0;
                menuEsc = true; menuStep = 6; menuFrames = 0;
            } else if (menuStep == 6 && menuFrames == 2) {
                const bool ignored = gameState == 1 && loading;
                menuOk &= ignored;
                printf("bench: menus: Esc on the loading screen does not pause (#388): %s\n", ignored ? "yes" : "NO");
                TakeScreenshot("bench_menu_loading_esc.png");
                // Resume from a pause taken on the loading screen (older builds allowed one) keeps the loading screen.
                pendingMenu.action = MenuAction::Resume; menuStep = 7; menuFrames = 0;
            } else if (menuStep == 7 && menuFrames == 2) {
                const bool kept = loading && gameState == 1;
                menuOk &= kept;
                printf("bench: menus: Resume during loading keeps the loading screen (#388): %s\n", kept ? "yes" : "NO");
                menuStep = 8; menuFrames = 0;
            } else if (menuStep == 8 && (!loading || menuFrames > 1800)) {
                // The loading screen still ends on its own once the world is built (0.5 s minimum on entering).
                const bool ended = !loading && gameState == 1;
                menuOk &= ended;
                printf("bench: menus: the loading screen then ends by itself after %.2f s: %s\n", GetTime() - loadStart, ended ? "yes" : "NO");
                loading = false;
                // #270: quit while the world streams hard (a far teleport: thousands of builds queued and running).
                if (Planetoid* p0 = world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get())) {
                    player.flyMode = false;
                    player.spawnOn(world, *p0, p0->dirToWorld(Vector3Negate(StartSpawnDirLocal(*p0)))); // the far side
                }
                menuStep = 9; menuFrames = 0;
            } else if (menuStep == 9 && menuFrames == 20) {
                printf("bench: menus: quitting while streaming: %zu jobs queued, %d running, %d GPU batches\n", WorldJobs().queued(), WorldJobs().running(), TerrainGpu_Pending());
                pendingMenu.action = MenuAction::SaveAndQuit; menuStep = 10; menuFrames = 0;
            } else if (menuStep == 10 && menuFrames == 3) {
                const bool fast = lastQuitMs < 1000.0;
                // (#556: the teleport's loading screen was still up at the quit; it must not outlive the world)
                const bool down = gameState == 0 && !loading && WorldJobs().queued() == 0 && TerrainGpu_Pending() == 0 && TerrainArena_Stats().vertexCapacity == 0;
                menuOk &= fast && down;
                printf("bench: menus: quit to the menu while streaming took %.0f ms (under 1 s), torn down (#270): %s -> %s\n", lastQuitMs,
                       fast && down ? "yes" : "NO", menuOk ? "PASS" : "FAIL");
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
                loading = !bench; loadStart = GetTime(); // the loading screen until the terrain around you is in
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
                    // The New World screen first (#504): seed, terrain parameters, a preview. Its choice is handled below.
                    NewWorld_Open();
                    gameState = 3; newWorldFrames = 0;
                    break;
                case MenuAction::LoadGame:
                    SetWorldSavePath(m.path.c_str());
                    if (LoadWorldSave(world)) { worldFresh = false; enterGame(); }
                    else GameUI_Toast(TextFormat("Could not load %s", m.path.c_str()));
                    break;
                case MenuAction::Resume: { const bool wasLoading = loading; enterGame(false); loading = wasLoading; break; } // a pause during loading (#388) keeps the loading screen
                case MenuAction::Save: GameUI_Toast(SaveWorld(world) ? "Saved" : "Save failed"); break;
                case MenuAction::SaveAndQuit:
                    if (!SaveWorld(world)) GameUI_Toast("Save failed");
                    // Nothing of the old world keeps running behind the menu (#270): jobs, the factory, water,
                    // streaming, bodies and their GPU meshes, animals, particles. New Game and Load build a world again.
                    if (telescope.active) Telescope_SetActive(world, false);
                    player.flyMode = false;
                    { const double q0 = GetTime(); unloadWorld(world); lastQuitMs = (GetTime() - q0) * 1000.0; } // (#270: must never block for seconds)
                    if (bench) printf("bench: quit to the menu: the world unloaded in %.0f ms\n", lastQuitMs);
                    // A quit from the loading screen (a spawn or teleport in progress, #556) takes it down too: only the in-game
                    // draw ends it, so it stayed up behind the menu (benches wait out a loading screen and never finished).
                    loading = false; spawnLoad = false;
                    Animals_Reset();
                    Particles_Clear();
                    worldFresh = false;
                    toMenu();
                    break;
                case MenuAction::Exit: windowKeepAlive = false; break;
                case MenuAction::ReportBug: BugReport_Open(world); break; // stays paused (#169)
                case MenuAction::None: break;
            }
            // The New World screen's choice from last frame (#504). --test-menus creates with its own seed after a moment.
            NewWorldChoice nw = pendingNewWorld;
            pendingNewWorld = NewWorldChoice{};
            if (gameState == 3 && (benchMenus || playNewSeconds > 0.0) && nw.kind == NewWorldChoice::None) { // --play-new (#387) too
                ++newWorldFrames;
                if (newWorldFrames == 15 && benchMenus) { TakeScreenshot("bench_menu_newworld.png"); printf("bench: menus: New Game -> the New World screen (#504): yes\n"); }
                if (newWorldFrames == 20) { nw.kind = NewWorldChoice::Create; nw.seed = worldGen.seed; }
            }
            if (gameState == 3 && nw.kind == NewWorldChoice::Create) {
                SetWorldSavePath(bench ? "saves/bench.sav" : Menu_NewSavePath().c_str()); // bench runs never touch the player's saves
                if (!worldFresh) {
                    if (bench) worldGen.seed = nw.seed; // benches write no parameters beside the save (default terrain)
                    else NewWorld_Apply(nw);           // the seed, and <save>.worldgen: generateWorld reads it
                    if (!bench) { settings.difficulty = nw.difficulty; Settings_Save(); } // the world keeps it (#526); the next New World screen starts from it
                    regenerateWorld(world);            // the erosion bake (or its cache) behind its progress screen
                    StartNewWorld(world, true, !bench);
                }
                worldFresh = false;
                enterGame();
            } else if (gameState == 3 && nw.kind == NewWorldChoice::Cancel) {
                gameState = 0; // back to the title screen; no world was made
            }
            // Esc in game pauses (unless it closes the inventory or build mode first); the pause menu
            // handles Esc itself (resume, or back out of Settings). Not on the loading screen (#388: Resume skipped it into a
            // half-built world).
            const bool escPressed = IsKeyPressed(KEY_ESCAPE) || menuEsc;
            menuEsc = false;
            if (gameState == 1 && !loading && m.action == MenuAction::None && escPressed) {
                if (Nav_MapOpen()) Nav_CloseMap();
                else if (GameUI_InventoryOpen()) GameUI_ToggleInventory();
                else if (buildTool.active) buildTool.active = false;
                else { ShowCursor(); Menu_ResetPause(); gameState = 2; justPaused = true; }
            }
        }

        

        GpuTimer_Mark("pre-draw (uploads, updates)");
        BeginDrawing();
        gInFrame = true; // (#362: world generation started from inside a frame - the debug menu - must not draw its own)
            ClearBackground(crackCheck ? MAGENTA : BLACK); // --crack-check: holes in the terrain show magenta
            GpuTimer_Enable(benchMeasuring || perfLog); GpuTimer_Mark("clear");
            GpuTimer_CpuMark("begin-draw");
            if (gameState == 0) {
                pendingMenu = Menu_Main();
            } else if (gameState == 3) {
                pendingNewWorld = NewWorld_Draw(); // the New World screen (#504)
            } else if (gameState == 1) {
                { const double d0 = GetTime(); world.drawScene(gameState); if (benchMeasuring) benchDrawMs += (GetTime() - d0) * 1000.0; }
                player.recordFly(world); // after the fly camera's input (in drawScene): the pose to carry next frame
                if (!noHud && !gHeadless) Viewmodel_Draw(world);   // the hands, over the world
                if (Recorder_WantsFrame(false)) { rlDrawRenderBatchActive(); Recorder_Capture(); } // --record-hud 0: before the UI
                if (!noHud && !gHeadless) { // --no-hud: clean shots (the menu slides, #386); --headless: nothing (#311)
                world.drawUI(gameState);
                Telescope_DrawHUD(world);
                BuildTool_DrawHUD();
                Hotbar_DrawHUD(ship.piloted);
                Machine_DrawHUD();
                if (!ship.piloted && !buildTool.active) Animals_DrawHUD();
                Suit_DrawHUD();
                Weather_DrawHUD();
                if (!ship.piloted && !telescope.active && !buildTool.active && !Rocket_Boarded()) { TerrainTool_DrawHUD(); Fishing_DrawHUD(); }
                GameUI_DrawHUD(world);
                Survival_DrawHUD();
                Objectives_DrawHUD();
                if (!Rocket_Boarded()) Nav_DrawHUD(world);
                Rocket_DrawHUD(world); // the orbit map covers the rest
                Nav_DrawMap(world);    // last: the full map covers everything
                }
                if (loading) {
                    // Loading screen: black, with the streaming progress, until the world around you is built
                    // (it used to show you standing in a void among trees and water while the ground popped in).
                    // The bar (#494) is work done: the world jobs finished since this screen opened over those plus the ones
                    // still pending (the old 1 - pending / peak went back and forth as the streamer queued more, and rescaled
                    // at each new peak). The spawn ring ("Preparing spawn") fills the first 80%. The target never goes back;
                    // the drawn bar eases toward it and fills to the end before the screen closes.
                    static double barFor = -1.0;    // the loadStart this bar belongs to
                    static uint64_t barDone0 = 0;
                    static float barTarget = 0.0f, barShown = 0.0f, barLowest = 0.0f, barJump = 0.0f; // the worst step back and forward (tests)
                    static bool barSawRing = false, barClosing = false;
                    if (barFor != loadStart) { barFor = loadStart; barDone0 = WorldJobs().done(); barTarget = barShown = 0.0f; barLowest = barJump = 0.0f; barSawRing = barClosing = false; }
                    const size_t pending = WorldJobs().queued() + WorldJobs().running();
                    float ringProgress = 1.0f;
                    const bool preparing = SpawnPreloadPending(&ringProgress); // first load: the 300 m ring around the spawn (#294)
                    barSawRing |= preparing;
                    const double doneJobs = (double)(WorldJobs().done() - barDone0);
                    const float built01 = WorldIdle() ? 1.0f : (float)(doneJobs / std::max(1.0, doneJobs + (double)pending));
                    float want = preparing ? 0.8f * std::clamp(ringProgress, 0.0f, 1.0f) : barSawRing ? 0.8f + 0.2f * built01 : built01;
                    if (barClosing) want = 1.0f;
                    barTarget = std::max(barTarget, want);
                    const float before = barShown;
                    barShown += (barTarget - barShown) * (1.0f - expf(-(barClosing ? 25.0f : 8.0f) * std::min(GameFrameTime(), 0.1f))); // the last fill is quick (~0.2 s)
                    if (barTarget - barShown < 0.002f) barShown = barTarget;
                    barLowest = std::min(barLowest, barShown - before);
                    if (!barClosing) barJump = std::max(barJump, barShown - before); // the biggest one-frame step while loading (a jump)
                    DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), BLACK);
                    Slideshow_Draw(0.5f); // screenshots behind the loading bar (#386), carrying on from the menu
                    const int bw = GetScreenWidth() / 3, bx = (GetScreenWidth() - bw) / 2, by = GetScreenHeight() / 2;
                    DrawText(preparing ? "Preparing spawn" : "Building the world", bx, by - 40, 28, RAYWHITE);
                    DrawRectangle(bx, by, bw, 14, Fade(RAYWHITE, 0.15f));
                    DrawRectangle(bx, by, (int)(bw * std::clamp(barShown, 0.02f, 1.0f)), 14, SKYBLUE);
                    // Entering the game holds it half a second; a spawn only until the streamer has had a few
                    // frames to queue the new place and has built it (where it is built already: at once).
                    ++loadFrames;
                    const bool built = !preparing && (spawnLoad ? WorldIdle() && loadFrames >= 3 : WorldIdle() && GetTime() - loadStart > 0.5);
                    if (built) barClosing = true; // fill to the end first (a few frames), then close
                    const bool timedOut = GetTime() - loadStart > (preparing ? 90.0 : 20.0);
                    if ((barClosing && barShown >= 0.999f) || timedOut) {
                        if (bench) printf("bench: loading screen after a spawn: %d frames, %.2f s; the bar never went back: %s (worst step %.4f), biggest one-frame step %.3f\n",
                                          loadFrames, GetTime() - loadStart, barLowest >= 0.0f ? "yes" : "NO", barLowest, barJump);
                        loading = false; spawnLoad = false;
                    }
                }
                // (No per-frame cursor warp: the cursor is locked with DisableCursor while playing, #204 / #205.)
                if (perfLog) {
                    // --perf-log (#204): every 5 s, the per-pass averages to log.txt, tagged by how the player moved (average
                    // speed, and mouse-look turning) so "standing" and "moving" can be compared from a real session.
                    static double perfT0 = GetTime(), perfSpeed = 0.0, perfLook = 0.0; static long perfN = 0, perfLoading = 0;
                    // On the loading screen the player is held (no nearest body yet: the velocity is the planet's spin, which
                    // read as "MOVING 89 m/s", #387): those frames count as loading, not movement.
                    if (loading) ++perfLoading;
                    else perfSpeed +=Vector3Length(Vector3Subtract(player.velocity, !player.grounded && player.nearest ? player.nearest->surfaceVelocity((player.position - world.renderOrigin).toFloat()) : Vector3{ 0, 0, 0 })); /* against the ground: airborne velocity is world-frame and includes the planet's spin (grounded it is surface-relative already) */ perfLook += Vector2Length(PlayerInput_Frame().look); ++perfN;
                    if (GetTime() - perfT0 >= 5.0) {
                        const long played = std::max(1L, perfN - perfLoading);
                        const double sp = perfSpeed / played, lk = perfLook / std::max(1L, perfN);
                        const char* tag = perfLoading * 2 > perfN ? "LOADING" : sp > 0.5 || lk > 0.5 ? "MOVING" : "still";
                        GpuTimer_LogAndReset(TextFormat("%s (speed %.1f m/s, look %.1f px/frame, %d FPS)", tag, sp, lk, GetFPS()));
                        perfT0 = GetTime(); perfSpeed = perfLook = 0.0; perfN = 0; perfLoading = 0;
                    }
                }
            } else if (gameState == 2) {
                world.drawScene(gameState);
                if (justPaused) BugReport_CaptureScene(); // the world as it paused, for a bug report's screenshot (#169)
                world.drawUI(gameState);
                const MenuResult m = BugReport_IsOpen() ? MenuResult{} : Menu_Pause(); // the bug report window takes the input (#169)
                if (!justPaused) pendingMenu = m; // the Esc that opened it is not also a "resume"
            }
            justPaused = false;
            if (gameState != 0 && !gHeadless) {
                GpuTimer_Mark("hud");
                rlImGuiBegin();
                DebugUI_Draw(world);
                GameUI_Draw(world);
                BugReport_Draw();
                BugReport_Update();
                rlImGuiEnd();
                GpuTimer_Mark("imgui");
            }
            if (settings.showFps) DrawFPS(10, 10); // every screen when the setting is on (#572: menus, loading, play), over everything
            rlDrawRenderBatchActive(); // the batched 2D (HUD text, bars, minimap) flushes here, not inside EndDrawing
            GpuTimer_Mark("2d-flush");
            if (Recorder_WantsFrame(true)) Recorder_Capture(); // the finished frame, HUD included (#375)
            Recorder_FrameDone();
        EndDrawing();
        gInFrame = false;
        {   // Unfocused (#510): at most ~30 frames a second (the streaming still gets every frame's budget), so a background
            // window does not burn the GPU. A wait, not SetTargetFPS: the FPS cap the player set is never touched.
            static double lastFrameEnd = 0.0;
            if (!focused) { const double el = GetTime() - lastFrameEnd; if (el < 1.0 / 30.0) WaitTime(1.0 / 30.0 - el); }
            lastFrameEnd = GetTime();
        }
        GpuTimer_Mark("swap+poll");
        PROFILE_FRAME();
    }
    Recorder_Shutdown(); // the frames still queued, the GIF and the contact sheet (#375)
    if (!hostSave.empty()) Host_SaveOnShutdown(world); // --host: the world on disk before it goes away

    }   // scenes destroyed here, while the GL context still exists
    NewWorld_Shutdown(); // a preview build in flight (#504)
    DebugUI_Shutdown();
    Audio_Shutdown();
    Factory::Stop();
    Host_Stop();
    CloseWindow();
}


