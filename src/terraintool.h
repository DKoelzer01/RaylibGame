#ifndef TERRAINTOOL_H
#define TERRAINTOOL_H

class Scene;

// First-person terrain editing: raycast from the crosshair, hold LMB to dig / RMB to add,
// mouse wheel resizes the brush, F5 saves, F9 loads.
struct TerrainToolSettings {
    float brushRadius = 4.0f;    // world units
    float digStrength = 0.6f;    // density change per stroke at the brush centre
    float strokesPerSecond = 10.0f;
    float reach = 150.0f;        // max raycast distance
};
extern TerrainToolSettings terrainTool;

void TerrainTool_Update(Scene& world, bool inputEnabled); // call once per frame before drawing
void TerrainTool_Draw3D();                                // inside BeginMode3D (Scene::onDraw3D)
void TerrainTool_DrawHUD();                               // 2D, after the scene

#endif
