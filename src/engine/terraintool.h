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
    bool sphereFill = false;     // #228: the fill shape is an exact sphere (else the soft blob); Ctrl+wheel with a filler in hand picks it, paid by volume
    float fillRadius = 1.0f;     // its radius (m)
    std::string heldTool;        // the tool item in hand (hotbar); empty or not carried: bare hands
    bool rmbBlocks = false;      // a shield in the off hand and nothing RMB uses in hand: RMB blocks, never builds (#165)
    int swings = 0;              // counts mining strokes and chops (the hands animate a swing on each)
    float swingTime = 0.3f;      // how long the current swing takes (its cooldown): the hands animate over it
    int gemBoost = 1;            // tests: multiplies the (rare) odds of a gem or quartz per rock stroke
    // A shovel's RMB on soft ground tills it into farmland (set by the game; returns true if it did).
    std::function<bool(class Planetoid&, Vector3 localGround)> till;
    // The item in hand became another (an empty bucket filled, a full one emptied): the game rebinds
    // the hotbar slot and counts quest events ("scoop", "pour").
    std::function<void(const std::string& from, const std::string& to, const char* event)> swapHeld;
};
extern TerrainToolSettings& terrainTool; // the local PlayerSession's (#313)

struct PlayerInput;
void TerrainTool_Update(Scene& world, const PlayerInput& in, bool editEnabled = true); // once per frame before drawing (input: #312); editEnabled=false keeps save/load (F5/F9) only
bool TerrainTool_Stroke(bool dig);                         // one dig/build stroke at the last crosshair hit
bool TerrainTool_IsFiller(const std::string& item);        // bulk ground placed with RMB (dirt, gravel, sand, stone...)
float TerrainTool_Cooldown();                              // 1 just after a swing .. 0 ready for the next
bool TerrainTool_Chop();     // one swing at the aimed plant (scripted tests); true when it came down
bool TerrainTool_Bucket(class Scene& world); // the bucket in hand: scoop or pour (RMB; scripted tests call it after aiming)
bool TerrainTool_AimingPlant(uint64_t* id = nullptr);
bool TerrainTool_LastHit(class Planetoid*& body, Vector3& local, Vector3& normalWorld); // the crosshair's terrain hit (placing torches, #218)
bool TerrainTool_TestFill(); // --test-fill (#228): an exact sphere of dirt, paid by volume
class Planetoid;
bool TerrainTool_DigAt(Planetoid& p, Vector3 local);              // tests: dig a point, gain what it held
bool TerrainTool_HarvestNearest(Planetoid& p, Vector3 local, bool tree, bool stone = false); // tests: fell / harvest the nearest plant (or pick up stones)
void TerrainTool_Draw3D();                                // inside BeginMode3D (Scene::onDraw3D)
bool TerrainTool_TestVein(Planetoid& body, Vector3 homeDir); // --test-vein (#207): a deposit mined 93% by hand comes loose
void TerrainTool_DrawHUD();                               // 2D, after the scene

// Multiplayer (#329): what a dig stroke yields, with no side effects (the local stroke and the server's validated one
// both use it, so a client's inventory gets exactly what the server credits). `centre`/`r` the carve sphere, `hitLocal`
// the aimed surface point (materials as drawn there), all body-local. Gems and seeds roll on GameRng (the server's).
#include <map>
struct DigYield {
    std::map<std::string, float> items;
    const char* gem = nullptr;   // a rare gem or quartz (+1)
    const char* seed = nullptr;  // grass turned up a seed (+1)
    bool ore = false;            // some ore came out (a vein may be nearly done)
};
DigYield TerrainTool_DigYield(Planetoid& p, Vector3 centre, float r, Vector3 hitLocal, bool softOnly, int gemBoost);
int TerrainTool_FillerMaterial(const std::string& item); // the VoxelMaterial a filler places (-1: not a filler)
const char* TerrainTool_ItemAt(Planetoid& p, Vector3 local, int material); // the item a material there digs out (ore: its kind)
// Set while joined to a server (#329): a dig or place stroke is sent as a request instead of applied and paid here;
// the server validates it, applies it in its order and credits or charges the inventory (deltas come back).
// dig: the carve sphere; place: the soft blob, `item` the filler that pays for it.
using TerrainToolNetEdit = std::function<void(Planetoid& p, bool dig, Vector3 hitLocal, Vector3 centre, float radius, const std::string& item)>;
extern TerrainToolNetEdit gTerrainToolNetEdit;

#endif
