#ifndef SCENE_H
#define SCENE_H

#include "object.h"
#include "textobject.h"
#include "gameobject.h"
#include "camera.h"
#include "frustum.h"

#include <raylib.h>
#include "rlgl.h"
#include "raymath.h"
#include "rlights.h"

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <random>
#include <array>
#include <unordered_map>

// A shader loaded from a vs/fs pair that reloads itself when either file changes on disk.
// A failed compile keeps the previous working shader.
struct HotShader {
    std::string vsPath, fsPath;
    Shader shader{};
    long vsTime = 0, fsTime = 0;

    // Returns true when a new shader was swapped in (callers must re-apply locations/uniforms).
    bool reloadIfChanged(bool force = false);
    void unload();
};

class Scene {
    public:
        std::string name; // Name of the scene
        Model skybox{};
        bool isActive; // Current game state
        std::vector<std::unique_ptr<Object>> objects; // Array of objects in the scene
        std::vector<std::unique_ptr<Object>> uiObjects; // Array of UI objects in the scene
        Camera3D camera{}; // Camera for the scene

        std::vector<Light> lights; // Array of lights in the scene
        HotShader lighting{ "resources/lighting.vs", "resources/lighting.fs" };
        HotShader depth{ "resources/depth.vs", "resources/depth.fs" };

        Material terrainMaterial{};   // lighting shader + albedo/normal/rough-ao maps + shadow map
        Material depthMaterial{};     // depth shader only
        RenderTexture2D shadowMap{};  // depth-only FBO; shadowMap.depth is the sampled texture
        Matrix lightSpaceMatrix = MatrixIdentity();

        // Sun: direction *towards* the light, and shadow-volume parameters.
        Vector3 sunDirection = { 1.0f, 0.0f, 0.0f };
        float shadowHalfExtent = 200.0f;   // half-width of the orthographic shadow box
        float sunDistance = 500.0f;

        bool debugDraw = false;  // F3: light marker, shadow box, shadow-map preview
        int lastDrawCalls = 0;

        Object rootObject;

        Scene(std::string name, bool isActive = false);
        virtual ~Scene();
        // Loads the planetoid textures (only the world scene needs them).
        void loadTerrainTextures();
        void drawScene(int gamestate);
        void drawUI(int gamestate);

    private:
        void applyLightingShader();   // (re)binds material + static uniforms after a (re)load
        int locViewPos = -1, locLightSpace = -1, locAmbient = -1, locShadowRes = -1, locTexScale = -1;
        bool terrainTexturesLoaded = false;
};

#endif
