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
    Vector3 benchCamOffset = { -70.0f, 35.0f, 70.0f }; // --cam dx dy dz: bench camera offset from the planetoid
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--bench") benchSeconds = (i + 1 < argc) ? atof(argv[i + 1]) : 10.0;
        if (std::string(argv[i]) == "--cam" && i + 3 < argc) { benchCamOffset = { (float)atof(argv[i+1]), (float)atof(argv[i+2]), (float)atof(argv[i+3]) }; }
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
    
    {   // Scenes own GPU resources: keep them in a scope that ends before CloseWindow().
    Scene mainMenu("mainMenu", true); // Main menu scene;
    Scene world("world", false);

    bool windowKeepAlive = true;

    mainMenuHandler(mainMenu);
    worldHandler(world);

    const double benchWarmup = 2.0;
    double benchStart = 0.0, benchMaxFrame = 0.0;
    long benchFrames = 0;
    if (bench) {
        gameState = 1;
        mainMenu.isActive = false;
        world.isActive = true;
        world.camera.position = Vector3Add(world.camera.target, benchCamOffset);
        benchStart = GetTime();
    }

    while (windowKeepAlive && !WindowShouldClose())
    {
        if (bench) {
            double t = GetTime() - benchStart;
            if (t > benchWarmup) {
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
            HideCursor();
        } 
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
                if (!bench) SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
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
        EndDrawing();
    }
    
    }   // scenes destroyed here, while the GL context still exists
    CloseWindow();
}


