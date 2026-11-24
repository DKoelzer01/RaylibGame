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
// (#1739) The static light grid (lightgrid.h) around the camera: a 128 m cube of it in a wrap-around 3D texture (unit 18,
// RGBA8: colour, level; 16^3 bricks updated as they change or as the camera moves into new ones), sampled by shadow.glsl's gridLight().
// `anchor` is a body-local cell whose coordinates are multiples of 128 near the camera; `gridFromRender` maps render space
// to body-local metres minus the anchor (inverse of the body's localModel at the anchor). `cam` is the camera's cell.
// grid null or gain <= 0: off (the shaders skip it). LocalLights_Bind sets its uniforms with the rest.
class LightGrid;
void LocalLights_UploadGrid(LightGrid* grid, const int cam[3], const int anchor[3], const Matrix& gridFromRender, float gain);
int LocalLights_GridBricksUploaded(); // bricks sent this frame (tests, profiling)

int LocalLights_Count();        // lights uploaded this frame
int LocalLights_MaxPerTile();   // the busiest tile's count this frame (tests)

#endif
