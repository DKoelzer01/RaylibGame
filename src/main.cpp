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

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif

using namespace std;

Logger logger("log.txt");

int gameState = 0; // 0: Main Menu, 1: Game, 2: Pause Menu

// --bench <seconds>: windowed 1920x1080, uncapped FPS, skips the menu, holds a fixed
// camera on the planetoid, prints frame stats to stdout and bench.txt, then exits.
bool benchMode = false;

int main(int argc, char** argv)
{
    double benchSeconds = 0.0;
    bool startWithOverlay = false;
    int benchDigStrokes = 0;     // --test-dig N: dig N strokes at the crosshair once the world is loaded
    bool benchSaveLoad = false;  // --test-saveload: then save, reload from disk, and bench the reloaded world
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
        if (std::string(argv[i]) == "--size" && i + 1 < argc) worldGen.planetoidSize = atoi(argv[i + 1]);
        if (std::string(argv[i]) == "--count" && i + 1 < argc) worldGen.planetoidCount = atoi(argv[i + 1]);
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
    
    DebugUI_Init();
    {   // Scenes own GPU resources: keep them in a scope that ends before CloseWindow().
    Scene mainMenu("mainMenu", true); // Main menu scene;
    Scene world("world", false);

    bool windowKeepAlive = true;

    mainMenuHandler(mainMenu);
    worldHandler(world);
    world.onDraw3D = TerrainTool_Draw3D;

    const double benchWarmup = 2.0;
    double benchStart = 0.0, benchMaxFrame = 0.0;
    long benchFrames = 0;
    if (bench) {
        gameState = 1;
        mainMenu.isActive = false;
        world.isActive = true;
        world.camera.position = Vector3Add(world.camera.target, benchCamOffset);
        world.debugDraw = startWithOverlay;
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
                } else if (WorldIdle() && benchSaveLoad) {
                    const Vector3 camPos = world.camera.position, camTarget = world.camera.target;
                    const int editsBefore = GetWorldStats().edits;
                    const bool ok = SaveWorld(world) && LoadWorldSave(world);
                    world.camera.position = camPos;
                    world.camera.target = camTarget;
                    printf("bench: save/load %s (%d edits saved)\n", ok ? "ok" : "FAILED", editsBefore);
                    benchSaveLoad = false;
                } else if ((WorldIdle() && benchDigStrokes == 0 && !benchSaveLoad) || GetTime() > 120.0) {
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
            if (DebugUI_IsOpen()) ShowCursor(); else HideCursor();
        }
        if (!bench && IsKeyPressed(KEY_F1)) {
            DebugUI_Toggle();
            if (!DebugUI_IsOpen()) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
        }
        world.cameraInputEnabled = !DebugUI_IsOpen();

        worldUpdate(world); // streams chunks in even while the menu is showing
        if (gameState == 1) TerrainTool_Update(world, !bench && !DebugUI_IsOpen());

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
                TerrainTool_DrawHUD();
                if (!bench && !DebugUI_IsOpen()) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
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
            if (gameState != 0) DebugUI_Draw(world);
        EndDrawing();
        PROFILE_FRAME();
    }

    }   // scenes destroyed here, while the GL context still exists
    DebugUI_Shutdown();
    CloseWindow();
}


