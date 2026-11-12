#ifndef TERRAINTOOL_H
#define TERRAINTOOL_H

#include <cstdint>
#include <raylib.h>
#include <functional>
#include <string>

class Scene;

// First-person terrain editing: raycast from the crosshair, hold LMB to dig / RMB to add,
// mouse wheel resizes the brush, F5 saves, F9 loads.
struct TerrainToolSettings {
    float brushRadius = 4.0f;    // world units
    float digStrength = 0.12f;   // density change per stroke at the brush centre (~0.5-1 m of rock)
    float strokesPerSecond = 10.0f;
    float reach = 12.0f;         // max raycast distance (arm's reach; fly mode uses 150)
    bool creative = false;       // building costs nothing
    bool terraform = false;      // debug (F1): the old brush tool - dig and fill any material, wheel resizes
    std::string heldTool;        // the tool item in hand (hotbar); empty or not carried: bare hands
    int swings = 0;              // counts mining strokes and chops (the hands animate a swing on each)
    // A shovel's RMB on soft ground tills it into farmland (set by the game; returns true if it did).
    std::function<bool(class Planetoid&, Vector3 localGround)> till;
};
extern TerrainToolSettings terrainTool;

void TerrainTool_Update(Scene& world, bool inputEnabled, bool editEnabled = true); // once per frame before drawing; editEnabled=false keeps save/load (F5/F9) only
bool TerrainTool_Stroke(bool dig);                         // one dig/build stroke at the last crosshair hit
bool TerrainTool_IsFiller(const std::string& item);        // bulk ground placed with RMB (dirt, gravel, sand, stone...)
bool TerrainTool_Chop();     // one swing at the aimed plant (scripted tests); true when it came down
bool TerrainTool_AimingPlant(uint64_t* id = nullptr);
class Planetoid;
bool TerrainTool_DigAt(Planetoid& p, Vector3 local);              // tests: dig a point, gain what it held
bool TerrainTool_HarvestNearest(Planetoid& p, Vector3 local, bool tree); // tests: fell / harvest the nearest plant
void TerrainTool_Draw3D();                                // inside BeginMode3D (Scene::onDraw3D)
void TerrainTool_DrawHUD();                               // 2D, after the scene

#endif
