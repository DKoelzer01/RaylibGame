#ifndef TERRAINTOOL_H
#define TERRAINTOOL_H

class Scene;

// First-person terrain editing: raycast from the crosshair, hold LMB to dig / RMB to add,
// mouse wheel resizes the brush, F5 saves, F9 loads.
struct TerrainToolSettings {
    float brushRadius = 4.0f;    // world units
    float digStrength = 0.12f;   // density change per stroke at the brush centre (~0.5-1 m of rock)
    float strokesPerSecond = 10.0f;
    float reach = 12.0f;         // max raycast distance (arm's reach; fly mode uses 150)
    bool creative = false;       // building costs nothing
};
extern TerrainToolSettings terrainTool;

void TerrainTool_Update(Scene& world, bool inputEnabled); // call once per frame before drawing
bool TerrainTool_Stroke(bool dig);                         // one dig/build stroke at the last crosshair hit
void TerrainTool_Draw3D();                                // inside BeginMode3D (Scene::onDraw3D)
void TerrainTool_DrawHUD();                               // 2D, after the scene

#endif
