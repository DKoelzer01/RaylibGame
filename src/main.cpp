#include <raylib.h>
#include <string>
#include <vector>
using namespace std;

struct GameState {
    int state = 0; // 0 = main menu, 1 = game, 2 = pause menu
};

struct textObject {
    string text = ""; // Text to display
    string type = "UI"; // Type of the object (e.g., "text" or "UI")
    Vector3 position = {0.0,0.0,0.0}; // Position of the text relative to the object
    Color color = WHITE; // Color of the text
    float size = 10; // Size of the text
};

struct Object {
    string type; // Type of the object (e.g., "cube", "sphere")
    Vector3 position = {0.0,0.0,0.0}; // Position of the object
    Vector3 rotation = {0.0,0.0,0.0}; // Rotation of the object
    Color color; // Color of the object
    float scale = 1.0; // Scale of the object
    textObject text = {}; // Text object for displaying text on the object
};

struct Scene {
    vector<Object> objects; // Array of objects in the scene
    vector<Object> uiObjects; // Array of UI objects in the scene
    Camera3D camera; // Camera for the scene
};

void drawScene(Scene *scene);

void drawUI(vector<Object> *uiObjects) {
    for (int i = 0; i < uiObjects->size(); i++) {
        Object obj = (*uiObjects)[i];
        DrawText(obj.text.text.c_str(), obj.text.position.x, obj.text.position.y, obj.text.size, obj.text.color);
    }
}

Scene mainMenu;
Scene world;




int main() 
{
    GameState gameState = {0};

    bool windowKeepAlive = true;

    mainMenu.camera = {0};
    mainMenu.camera.position = (Vector3){ 0.0f, 10.0f, 10.0f };
    mainMenu.camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    mainMenu.camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    mainMenu.camera.fovy = 45.0f;
    mainMenu.camera.projection = CAMERA_PERSPECTIVE;
    mainMenu.uiObjects.push_back({ "text", { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, WHITE, 1.0f,  {"Welcome to the Game!","UI"} });
    mainMenu.uiObjects.push_back({ "text", { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, WHITE, 1.0f, {"Press ENTER to start","UI",{0,20,0}} });
    mainMenu.uiObjects.push_back({ "text", { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, WHITE, 1.0f, {"Press ESC to exit","UI",{0,40,0}} });

    world.objects.push_back({ "cube", { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, RED, 1.0f });
    world.objects.push_back({ "cube", { 0.0f, 4.0f, 0.0f }, { 0.0f, 2.0f, 2.0f }, GREEN, 1.0f });
    world.objects.push_back({ "cube", { 5.0f, 0.0f, 0.0f }, { 2.0f, 2.0f, 0.0f }, BLUE, 1.0f });
    world.camera = {0};
    world.camera.position = (Vector3){ 0.0f, 10.0f, 10.0f };
    world.camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    world.camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    world.camera.fovy = 45.0f;
    world.camera.projection = CAMERA_PERSPECTIVE;

    InitWindow(GetMonitorWidth(0), GetMonitorHeight(0), "My first RAYLIB program!");
    SetTargetFPS(60);
    
    SetExitKey(KEY_NULL);

    if (!IsWindowFullscreen()) {
        ToggleFullscreen();
    }

    while (windowKeepAlive)
    {
        switch(gameState.state) {
            case 0: // Main Menu
                if (IsKeyPressed(KEY_ENTER)) {
                    printf("Start Game\n");
                    gameState.state = 1; // Start game
                }
                if (IsKeyPressed(KEY_ESCAPE)) {
                    printf("Exit game from main menu\n");
                    windowKeepAlive = false; // Exit game
                }
                break;
            case 1: // Game
                if (IsKeyPressed(KEY_ESCAPE)) {
                    printf("Game paused\n");
                    gameState.state = 2; // Pause game
                }
                break;
            case 2: // Pause Menu
                if (IsKeyPressed(KEY_R)) {
                    printf("Game unpaused\n");
                    gameState.state = 1; // Resume game
                }
                if (IsKeyPressed(KEY_ESCAPE)) {
                    printf("Exit game from pause\n");
                    windowKeepAlive = false; // Exit game
                }
                break;
        }

        BeginDrawing();
            ClearBackground(BLACK);
            DrawFPS(100,100);
            if (gameState.state == 0) {
                drawScene(&mainMenu);
                drawUI(&mainMenu.uiObjects);
            } else if (gameState.state == 1) {
                drawScene(&world);
                drawUI(&world.uiObjects);
            } else if (gameState.state == 2) {
                drawScene(&world);
                drawUI(&world.uiObjects);
                EnableCursor();
                DrawText("Game Paused", GetScreenWidth()/2 - MeasureText("Game Paused", 20)/2, GetScreenHeight()/2 - 10, 20, WHITE);
                DrawText("Press R to resume", GetScreenWidth()/2 - MeasureText("Press R to resume", 20)/2, GetScreenHeight()/2 + 10, 20, WHITE);
                DrawText("Press ESC to exit", GetScreenWidth()/2 - MeasureText("Press ESC to exit", 20)/2, GetScreenHeight()/2 + 30, 20, WHITE);
            }
        EndDrawing();
    }
    
    CloseWindow();
}

void drawScene(Scene *scene) {
    BeginMode3D(scene->camera);
    UpdateCamera(&scene->camera, CAMERA_FREE);
    for (int i = 0; i < scene->objects.size(); i++) {
        Object obj = scene->objects[i];
        if( obj.type == "text") {
            // TODO
            // Draw 3D text
            continue;
        } else if (obj.type == "cube") {
            DrawCube(obj.position, obj.scale, obj.scale, obj.scale, obj.color);
            DrawCubeWires(obj.position, obj.scale, obj.scale, obj.scale, BLACK);
        }
    }
    
    EndMode3D();
}