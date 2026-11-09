#include <raylib.h>
#include <string>
#include <vector>
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
    double rideStart = -1.0;
    int ridePhase = 0;
    Vector3 rideFwd0{}, ridePos0{}, rideTakeoff{};
    float rideViewDrift = 0, ridePosDrift = 0, rideSurfaceSpeed = 0;
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
        if (std::string(argv[i]) == "--test-suffocate") benchSuffocate = true;
        if (std::string(argv[i]) == "--test-orbit" && i + 1 < argc) benchOrbitSeconds = atof(argv[i + 1]);
        if (std::string(argv[i]) == "--test-ride") benchRide = true;
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
        if (std::string(argv[i]) == "--size" && i + 1 < argc) worldGen.planetoidSize = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--count" && i + 1 < argc) worldGen.planetoidCount = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--system" && i + 1 < argc) worldGen.systemRadius = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--lod" && i + 1 < argc) lodSettings.lodFactor = (float)atof(argv[i + 1]);
        if (std::string(argv[i]) == "--first" && i + 3 < argc) worldGen.firstPosition = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) };
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
        SetTargetFPS(0);
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
        player.flyMode = benchWalkSeconds <= 0.0; // the bench drives the camera directly unless walking
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
                        const float g = Player::SurfaceGravity(R) * (R * R) / (d * d);
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
                } else if (benchRide && WorldIdle()) {
                    // Scripted ride: stand on the equator of a body spinning 10 deg/s. The view must not drift
                    // against the surface, and a jump must land where it took off (momentum is kept).
                    Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
                    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
                    const double t = rideStart < 0.0 ? 0.0 : GetTime() - rideStart;
                    if (rideStart < 0.0) {
                        p->spin = { 0.0f, 10.0f, 0.0f };
                        player.spawnOn(world, *p, { 1.0f, 0.0f, 0.0f });
                        rideStart = GetTime();
                        ridePhase = 0;
                    } else if (ridePhase == 0 && t > 1.5) {
                        rideFwd0 = p->dirToLocal(player.forward);
                        ridePos0 = p->toLocal(feet);
                        ridePhase = 1;
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
                        printf("bench: ride test (10 deg/s, surface %.1f m/s): standing 3 s -> view drift %.2f deg, position drift %.2f m; jump landed %.2f m from takeoff (body frame)%s\n",
                               rideSurfaceSpeed, rideViewDrift, ridePosDrift, Vector3Distance(rideTakeoff, p->toLocal(feet)), ridePhase == 3 ? "" : " [never landed]");
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
                } else if ((WorldIdle() && benchDigStrokes == 0 && !benchSaveLoad && benchWalkSeconds <= 0.0 && benchMineStrokes == 0 && minePhase != 2 && benchFlySeconds <= 0.0 && !benchSuffocate && benchOrbitSeconds <= 0.0 && !benchRide) || GetTime() > 180.0) {
                    benchStart = GetTime();
                    printf("bench: world ready after %.2f s (%d chunks, avg build %.2f ms, %d worker threads)\n",
                        GetWorldStats().lastGenerationSeconds, GetWorldStats().chunks, GetWorldStats().avgBuildMs, GetWorldStats().workerThreads);
                }
            }
            double t = (benchStart < 0.0) ? 0.0 : GetTime() - benchStart;
            if (benchStart >= 0.0 && t > benchWarmup) {
                benchFrames++;
                if (GetFrameTime() > benchMaxFrame) benchMaxFrame = GetFrameTime();
            }
            if (t > benchWarmup + benchSeconds) {
                double avgMs = 1000.0 * benchSeconds / (double)benchFrames;
                const char* line = TextFormat("bench: %ld frames in %.1fs | avg %.1f FPS | avg %.2f ms | worst %.2f ms\n",
                    benchFrames, benchSeconds, 1000.0 / avgMs, avgMs, benchMaxFrame * 1000.0);
                printf("%s", line);
                SaveFileText("bench.txt", (char*)line);
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
            if (!ship.piloted) player.applyToCamera(world);
        }

        world.updateFloatingOrigin(); // keep the camera near render-space 0 (float precision)
        worldUpdate(world, 3.0f, gameState == 1); // LOD + streaming always; bodies only spin while playing
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
        if (gameState == 1) TerrainTool_Update(world, gameInput && !ship.piloted && alive);

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
            ClearBackground(BLACK);
            if (gameState == 0) {
                mainMenu.drawScene(gameState);
                mainMenu.drawUI(gameState);
            } else if (gameState == 1) {
                world.drawScene(gameState);
                world.drawUI(gameState);
                if (!ship.piloted) TerrainTool_DrawHUD();
                GameUI_DrawHUD(world);
                Survival_DrawHUD();
                Objectives_DrawHUD();
                if (!ship.piloted && ship.playerInRange(world, (player.position - world.renderOrigin).toFloat(), 5.0f))
                    DrawText("[E] Board ship", GetScreenWidth() / 2 - 60, GetScreenHeight() / 2 + 40, 20, SKYBLUE);
                if (ship.piloted) DrawText("SHIP: WASD thrust  Space/Ctrl up/down  mouse steer  Z/C roll  Shift boost  X brake  [2] power cell  E exit",
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


