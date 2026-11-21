#ifndef LOCALLIGHTS_H
#define LOCALLIGHTS_H

// Local point lights (#139): lamps, campfires, glowing machines, torches. The game fills a list each frame (render
// space); LocalLights_Upload keeps the nearest kMaxLocalLights within reach of the camera, bins them into 16 px screen
// tiles on the CPU (each light's projected sphere) and uploads two small textures: the lights (2 texels each) and the
// tiles (a count, then up to 15 light indices). shadow.glsl's localLights() loops over its tile's lights, so the cost per
// pixel is the lights that actually reach it. Textures, not SSBOs: the world shaders stay #version 330.

#include <raylib.h>
#include <vector>

struct LocalLight {
    Vector3 position;   // render space
    float range;        // m: the light reaches zero here
    Vector3 color;      // linear 0..1
    float intensity;
};

constexpr int kMaxLocalLights = 128;
constexpr int kLightTileSize = 16;      // px
constexpr int kLightsPerTile = 15;

// After the 3D camera and projection are set (the view-projection and viewport are read from rlgl/GL).
void LocalLights_Upload(const std::vector<LocalLight>& lights, Vector3 cameraPos);
// The textures at units 13 and 14 and the uniforms (localLightTex, localTileTex, localTileInfo) for one shader.
void LocalLights_Bind(Shader& shader);
int LocalLights_Count();        // lights uploaded this frame
int LocalLights_MaxPerTile();   // the busiest tile's count this frame (tests)

#endif
