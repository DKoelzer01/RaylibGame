#include "scene.h"

Scene::Scene(std::string name, bool isActive) : name(name), isActive(isActive) { 
    std::cout << "Scene created: " << name << std::endl;
}

void Scene::drawScene(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active
    BeginMode3D(camera);
    if(gamestate != 0 && gamestate != 2) { // If not in main menu or pause menu
        UpdateCamera(&camera, CAMERA_FREE);
    }
    
    for (const auto& objPtr : objects) { objPtr->draw(); }
    EndMode3D();
}

void Scene::drawUI(int gamestate) {
    if (!isActive) return; // Skip drawing if the scene is not active
    for (const auto& objPtr : uiObjects) { objPtr->draw(); }
}
