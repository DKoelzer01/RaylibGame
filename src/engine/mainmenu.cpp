#include "mainmenu.h"

// The main menu scene only holds a camera; its buttons are drawn by src/space/menus.cpp.
void mainMenuHandler(Scene& mainMenu) {
    mainMenu.camera.position = { 0.0f, 10.0f, 10.0f };
    mainMenu.camera.target = { 0.0f, 0.0f, 0.0f };
    mainMenu.camera.up = { 0.0f, 1.0f, 0.0f };
    mainMenu.camera.fovy = 45.0f;
    mainMenu.camera.projection = CAMERA_PERSPECTIVE;
}
