#ifndef SCENE_H
#define SCENE_H

#include "object.h"
#include "textobject.h"
#include "gameobject.h"
#include "camera.h"

#include <raylib.h>
#include "rlgl.h"
#include "raymath.h"

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <random>
#include <array>
#include <unordered_map>


class Scene {
    public:
        std::string name; // Name of the scene
        Model skybox;
        bool isActive; // Current game state
        std::vector<std::unique_ptr<Object>> objects; // Array of objects in the scene
        std::vector<std::unique_ptr<Object>> uiObjects; // Array of UI objects in the scene
        Camera3D camera; // Camera for the scene

        Object rootObject;

        Scene(std::string name, bool isActive = false);
        void drawScene(int gamestate);
        void drawUI(int gamestate);
};

#endif