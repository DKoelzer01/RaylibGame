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
#include "locallights.h"

#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <sstream>
#include <random>
#include <array>
#include <unordered_map>

// (#1089) The terrain tool's dig / fill sphere for this frame: the terrain shader draws where it cuts the ground (a line on
// the real surface, a faint tint inside). Body-local, so a render-origin shift between the tool's update and the draw is
// harmless; body = nullptr: none. Set by TerrainTool_Update, read by Scene when it binds the lighting shader.
class Planetoid;
struct ToolSphereHint { const Planetoid* body = nullptr; Vector3 local{ 0, 0, 0 }; float radius = 0.0f; Color color{ 0, 0, 0, 0 }; double time = 0.0; }; // time: when set (stale after 0.1 s: modes without the tool)
inline ToolSphereHint gToolSphere;
#include <functional>

// A shader loaded from a vs/fs pair that reloads itself when either file changes on disk.
// A failed compile keeps the previous working shader.
struct HotShader {
    std::string vsPath, fsPath;
    Shader shader{};
    long vsTime = 0, fsTime = 0;
    double lastCheck = -1e9;     // reloadIfChanged looks at the files at most twice a second
    std::vector<std::string> includes; // files pulled in by `#include "..."` (watched with the fs)

    // Returns true when a new shader was swapped in (callers must re-apply locations/uniforms).
    bool reloadIfChanged(bool force = false);
    void unload();
};

extern bool benchMode; // true in --bench runs: no input-driven camera, scripted tests
void SetShadowCompare(unsigned int depthTexture, bool on); // shadow maps: hardware compare on/off (debug previews read raw depth)

class Scene {
    public:
        std::string name; // Name of the scene
        Model skybox{};
        bool isActive; // Current game state
        std::vector<std::unique_ptr<Object>> objects; // Array of objects in the scene
        std::vector<std::unique_ptr<Object>> uiObjects; // Array of UI objects in the scene
        Camera3D camera{}; // Camera for the scene

        std::vector<Light> lights; // Array of lights in the scene
        std::vector<LocalLight> localLights; // lamps, campfires, glowing machines, render space: filled by the game each frame (#139)
        HotShader lighting{ "resources/lighting.vs", "resources/lighting.fs" };
        HotShader depth{ "resources/depth.vs", "resources/depth.fs" };
        HotShader water{ "resources/water.vs", "resources/water.fs" };
        HotShader grass{ "resources/grass.vs", "resources/grass.fs" };
        HotShader trees{ "resources/trees.vs", "resources/trees.fs" };
        HotShader treeDepth{ "resources/trees.vs", "resources/tree_depth.fs" }; // trees and bushes in the shadow pass (fading with distance, #541)
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
            Vector3 toLocalX = { 1, 0, 0 }, toLocalY = { 0, 1, 0 }, toLocalZ = { 0, 0, 1 }; // planet-local images of the world axes (clouds)
            float cloudHeight = 1500.0f;  // the cloud layer, above sea level (m)
        } atmosphere;

        Material terrainMaterial{};   // lighting shader + albedo/normal/rough-ao maps + shadow map
        TextureArray terrainLayers{};     // dirt/sand/gravel/grass/snow/ice/stone layers (planet palette)
        Material depthMaterial{};     // depth shader only
        RenderTexture2D shadowMap{};  // depth-only FBO; shadowMap.depth is the sampled texture
        RenderTexture2D shadowNear{}; // the near cascade (sharp shadows close by)
        // #138: a wide terrain-only cascade (+-2 km, ~2 m texels) past the far one, so mountains shadow valleys at dawn and
        // dusk. Read by the terrain shader only (lighting.fs, material slot EMISSION). RAYLIBGAME_WIDE_SHADOWS=0: off.
        RenderTexture2D shadowWide{};
        Matrix lightSpaceWide{};
        float shadowWideHalfExtent = 2000.0f;
        float shadowWideDistance = 4000.0f; // the light camera back along the sun (depth range twice this)
        bool wideShadows = true;
        int wideEvery = 4;          // render it every Nth frame (world-anchored: 2 m texels hardly notice a few frames of spin)
        int wideFrame = 0;
        DVec3 wideOrigin{};         // the render origin it was drawn with (the matrix is moved along if the origin shifts)
        Matrix lightSpaceWideDrawn{};
        Matrix lightSpaceNear = MatrixIdentity();
        float shadowNearHalfExtent = 32.0f; // half-width of the near cascade's box (2048 texels: ~3 cm)
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
        Vector3 shadowSunLocal{ 0, 1, 0 }; uint64_t shadowSunBody = 0; // the stepped shadow light direction (ground frame)

        bool debugDraw = false;  // F3: light marker, shadow box, shadow-map preview
        bool cameraInputEnabled = true;  // false while the debug UI owns the mouse
        float texScale = 0.05f;          // terrain texture repeats per world unit
        Vector3 ambientColor = { 1.0f / 255.0f, 1.0f / 255.0f, 1.0f / 255.0f }; // linear; night side near-black (use the headlamp)
        void reloadShaders();            // force reload (debug UI button)
        int lastDrawCalls = 0;
        int materialPalette = 0;         // lighting.fs: 0 = space materials, 1 = fantasy
        bool drawSkybox = true;          // false: the game clears to a sky colour instead
        Color clearColor = BLACK;        // the scene target's clear (scenetarget.h; --crack-check: magenta)
        Vector3 fogColor = { 0, 0, 0 };  // linear-space fog colour (fantasy: sky haze)
        float fogDensity = 0.0f;         // 0 = no fog
        float weatherCloud = 0.0f;       // weather (set by the game each frame): cloud cover 0..1 greys the sky, dims the sun
        float weatherFog = 0.0f;         // ground fog 0..1 (thick at 1: ~60 m visibility)
        float weatherRain = 0.0f;        // rain falling 0..1 (set by the game): rings on water surfaces (#407)
        Vector3 cloudDrift = { 0, 0, 0 }; // weather: how far the wind has carried the clouds (planet-local m)
        std::function<void()> onDraw3D;  // extra world-space drawing at the end of the main pass
        std::function<void(Shader&)> frameLighting;  // set during onDraw3D: binds shadow cascades + headlamp to a shader
        Vector3 frameSun{ 0, 1, 0 }, frameAmbient{}, frameFog{}, frameAmbientUp{}; /* frameAmbientUp: hemisphere ambient up, #130 */ float frameFogDensity = 0.0f;
        Vector3 sunFade{ 0, 0, 0 };            // 1 - the sun light's colour this frame (#276: warm and dim when low), for shadow.glsl's sunTint()
        Vector3 lightDirection{ 1, 0, 0 };     // the directional light this frame: the sun, or the moon at night (#224); the sky keeps sunDirection
        float moonLight = 0.0f, moonPhase = 0.0f; // the moon's share as the light (0..1) and its lit fraction (#224)
        Vector4 frameFogShell{};
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

// Re-centres render space on `centre` for the scope's lifetime (#314, multiplayer): simulation code inside it (player physics)
// sees render-space floats relative to its own player, so it is precise and independent of where the camera and the client's
// floating origin are; a server runs each player in its own scope. The camera moves along; the origin is restored after.
class RenderOriginScope {
public:
    RenderOriginScope(Scene& scene, const DVec3& centre);
    ~RenderOriginScope();
    RenderOriginScope(const RenderOriginScope&) = delete;
    RenderOriginScope& operator=(const RenderOriginScope&) = delete;
private:
    Scene& scene;
    DVec3 saved;
};

#endif
