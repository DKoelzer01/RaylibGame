#include <raylib.h>
#include <string>
#include <vector>
#include "mainmenu.h"
#include "world.h"
#include "scene.h"
#include "object.h"
#include "utils/logger.h"
#include "rlights.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#endif

using namespace std;

Logger logger("log.txt");

int gameState = 0; // 0: Main Menu, 1: Game, 2: Pause Menu

int main() 
{
    SetConfigFlags(FLAG_MSAA_4X_HINT);
    SetConfigFlags(FLAG_VSYNC_HINT); // Enable VSync for smoother rendering
    SetConfigFlags(FLAG_WINDOW_RESIZABLE); // Allow window resizing
    SetConfigFlags(FLAG_WINDOW_UNDECORATED); // Remove window decorations for a cleaner look
    SetConfigFlags(FLAG_WINDOW_HIGHDPI); // Enable High DPI support
    InitWindow(GetMonitorWidth(0), GetMonitorHeight(0), "My first RAYLIB program!");
    SetTargetFPS(60);
    SetExitKey(KEY_NULL);
    SetTraceLogLevel(LOG_WARNING); // Set log level to warning to reduce output noise
    
    Scene mainMenu("mainMenu", true); // Main menu scene;
    Scene world("world", false);;

    bool windowKeepAlive = true;

    mainMenuHandler(mainMenu);
    worldHandler(world);

    while (windowKeepAlive)
    {
        if (!IsWindowFocused()) {
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
                SetMousePosition(GetScreenWidth()/2, GetScreenHeight()/2);
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
    
    CloseWindow();
}


