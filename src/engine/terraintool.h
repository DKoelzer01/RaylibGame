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
    bool holdingItem = false;    // food or feed in hand (hotbar): RMB never places terrain with it (#201)
    bool terraform = false;      // debug (F1): the old brush tool - dig and fill any material, wheel resizes
    std::string heldTool;        // the tool item in hand (hotbar); empty or not carried: bare hands
    int swings = 0;              // counts mining strokes and chops (the hands animate a swing on each)
    float swingTime = 0.3f;      // how long the current swing takes (its cooldown): the hands animate over it
    int gemBoost = 1;            // tests: multiplies the (rare) odds of a gem or quartz per rock stroke
    // A shovel's RMB on soft ground tills it into farmland (set by the game; returns true if it did).
    std::function<bool(class Planetoid&, Vector3 localGround)> till;
    // The item in hand became another (an empty bucket filled, a full one emptied): the game rebinds
    // the hotbar slot and counts quest events ("scoop", "pour").
    std::function<void(const std::string& from, const std::string& to, const char* event)> swapHeld;
};
extern TerrainToolSettings terrainTool;

void TerrainTool_Update(Scene& world, bool inputEnabled, bool editEnabled = true); // once per frame before drawing; editEnabled=false keeps save/load (F5/F9) only
bool TerrainTool_Stroke(bool dig);                         // one dig/build stroke at the last crosshair hit
bool TerrainTool_IsFiller(const std::string& item);        // bulk ground placed with RMB (dirt, gravel, sand, stone...)
float TerrainTool_Cooldown();                              // 1 just after a swing .. 0 ready for the next
bool TerrainTool_Chop();     // one swing at the aimed plant (scripted tests); true when it came down
bool TerrainTool_Bucket(class Scene& world); // the bucket in hand: scoop or pour (RMB; scripted tests call it after aiming)
bool TerrainTool_AimingPlant(uint64_t* id = nullptr);
class Planetoid;
bool TerrainTool_DigAt(Planetoid& p, Vector3 local);              // tests: dig a point, gain what it held
bool TerrainTool_HarvestNearest(Planetoid& p, Vector3 local, bool tree, bool stone = false); // tests: fell / harvest the nearest plant (or pick up stones)
void TerrainTool_Draw3D();                                // inside BeginMode3D (Scene::onDraw3D)
void TerrainTool_DrawHUD();                               // 2D, after the scene

#endif
