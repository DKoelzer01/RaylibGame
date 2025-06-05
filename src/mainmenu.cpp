#include "mainmenu.h"

void mainMenuHandler(Scene& mainMenu) {
    mainMenu.camera.position = { 0.0f, 10.0f, 10.0f };
    mainMenu.camera.target = { 0.0f, 0.0f, 0.0f };
    mainMenu.camera.up = { 0.0f, 1.0f, 0.0f };
    mainMenu.camera.fovy = 45.0f;
    mainMenu.camera.projection = CAMERA_PERSPECTIVE;

    
    mainMenu.uiObjects.push_back(std::make_unique<TextObject>(
        "Test Game", 
        Vector3{ static_cast<float>(GetScreenWidth()/2 - MeasureText("Test Game", 40)/2), static_cast<float>(GetScreenHeight()/2 - 30), 0.0f },
        Vector3{ 0.0f, 0.0f, 0.0f },
        WHITE,
        1.0f,
        40));

    mainMenu.uiObjects.push_back(std::make_unique<TextObject>(
        "Press Enter to Start", 
        Vector3{ static_cast<float>(GetScreenWidth()/2 - MeasureText("Press Enter to Start", 40)/2), static_cast<float>(GetScreenHeight()/2 + 10), 0.0f },
        Vector3{ 0.0f, 0.0f, 0.0f },
        WHITE,
        1.0f,
        40));
    
    mainMenu.uiObjects.push_back(std::make_unique<TextObject>(
        "Press ESC to exit", 
        Vector3{ static_cast<float>(GetScreenWidth()/2 - MeasureText("Press ESC to exit", 40)/2), static_cast<float>(GetScreenHeight()/2 + 50), 0.0f },
        Vector3{ 0.0f, 0.0f, 0.0f },
        WHITE,
        1.0f,
        40));
}