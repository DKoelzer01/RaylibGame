#ifndef SCENE_H
#define SCENE_H

#include "object.h"
#include "textobject.h"
#include "gameobject.h"
#include "camera.h"
#include "frustum.h"
#include "terrainlayers.h"

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
#include <functional>

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

extern bool benchMode; // true in --bench runs: no input-driven camera, scripted tests

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
        HotShader water{ "resources/water.vs", "resources/water.fs" };
        HotShader grass{ "resources/grass.vs", "resources/grass.fs" };
        HotShader trees{ "resources/trees.vs", "resources/trees.fs" };
        HotShader farSphere{ "resources/farsphere.vs", "resources/farsphere.fs" };
        HotShader parts{ "resources/parts.vs", "resources/parts.fs" };
        HotShader atmosphereShader{ "resources/skybox.vs", "resources/atmosphere.fs" };
        Material waterMaterial{};
        Model atmosphereDome{};       // a cube like the skybox, drawn over it with the atmosphere shader

        // Atmosphere of the planet the camera is near (set by the game each frame). Tints the sky over
        // the star skybox by sun angle and altitude, and drives day/night haze (fog) and sky ambient.
        struct Atmosphere {
            bool enabled = false;
            Vector3 centre = { 0, 0, 0 }; // render space
            float radius = 0.0f;          // sea level
            float height = 25000.0f;      // fades out from sea level to this altitude
            float hazeHeight = 8000.0f;   // haze (fog) lies in this layer: looking up through it, distant bodies stay clear
        } atmosphere;

        Material terrainMaterial{};   // lighting shader + albedo/normal/rough-ao maps + shadow map
        TextureArray terrainLayers{};     // dirt/sand/gravel/grass/snow/ice/stone layers (planet palette)
        Material depthMaterial{};     // depth shader only
        RenderTexture2D shadowMap{};  // depth-only FBO; shadowMap.depth is the sampled texture
        Matrix lightSpaceMatrix = MatrixIdentity();
        // Tests: a world point whose sub-texel drift across the shadow map is summed (shadows crawling).
        std::function<DVec3()> shadowProbe;
        float shadowProbeSlide = 0.0f;
        Vector2 shadowProbeLast{};
        bool shadowProbeHave = false;

        // Sun: direction *towards* the light, and shadow-volume parameters.
        Vector3 sunDirection = { 1.0f, 0.0f, 0.0f };
        float shadowHalfExtent = 200.0f;   // half-width of the orthographic shadow box
        float sunDistance = 500.0f;

        bool debugDraw = false;  // F3: light marker, shadow box, shadow-map preview
        bool cameraInputEnabled = true;  // false while the debug UI owns the mouse
        float texScale = 0.05f;          // terrain texture repeats per world unit
        Vector3 ambientColor = { 1.0f / 255.0f, 1.0f / 255.0f, 1.0f / 255.0f }; // linear; night side near-black (use the headlamp)
        void reloadShaders();            // force reload (debug UI button)
        int lastDrawCalls = 0;
        int materialPalette = 0;         // lighting.fs: 0 = space materials, 1 = fantasy
        bool drawSkybox = true;          // false: the game clears to a sky colour instead
        Vector3 fogColor = { 0, 0, 0 };  // linear-space fog colour (fantasy: sky haze)
        float fogDensity = 0.0f;         // 0 = no fog
        float weatherCloud = 0.0f;       // weather (set by the game each frame): cloud cover 0..1 greys the sky, dims the sun
        float weatherFog = 0.0f;         // ground fog 0..1 (thick at 1: ~60 m visibility)
        std::function<void()> onDraw3D;  // extra world-space drawing at the end of the main pass
        // Extra geometry in the depth passes (rings, ring.h): range and draw in the far pass, draw in
        // the mid and near passes (rebuild = the first of the frame).
        std::function<bool(double& dmin, double& dmax)> farExtraRange;
        std::function<void(const RenderContext&)> drawFarExtra;
        std::function<void(const RenderContext&, bool rebuild)> drawPassExtra;

        // Floating origin: render-space (float) = world (double) - renderOrigin. The camera is kept
        // near the render-space origin so vertex positions stay precise anywhere in the system.
        DVec3 renderOrigin{};
        float originRebaseDistance = 1024.0f;
        void syncRenderPositions();               // position = worldPos - renderOrigin for top-level objects
        void updateFloatingOrigin();              // re-centre on the camera when it strays too far
        void setCameraWorld(DVec3 position, DVec3 target);
        DVec3 cameraWorld() const { return renderOrigin + camera.position; }

        Object rootObject;

        Scene(std::string name, bool isActive = false);
        virtual ~Scene();
        // Loads the planetoid textures (only the world scene needs them).
        void loadTerrainTextures();
        void drawScene(int gamestate);
        void drawUI(int gamestate);

    private:
        void applyLightingShader();   // (re)binds material + static uniforms after a (re)load
        int locTerrainLayers = -1, locSunOccluder = -1;
        int locViewPos = -1, locLightSpace = -1, locAmbient = -1, locShadowRes = -1, locTexScale = -1, locPalette = -1, locFogColor = -1, locFogDensity = -1, locFogShell = -1, locChunkOrigin = -1, locChunkClip = -1;
        bool terrainTexturesLoaded = false;
};

#endif
