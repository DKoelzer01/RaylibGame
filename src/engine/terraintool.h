#ifndef TERRAINTOOL_H
#define TERRAINTOOL_H

#include <cstdint>
#include <raylib.h>
#include <functional>

// (#827) Draws the key hints line at (x, y) (the game: its bindings, size and fade); unset, the tool draws a fixed line.
inline std::function<void(int x, int y)> gKeyHintsHook;
#include <string>

class Scene;

// First-person terrain editing: raycast from the crosshair, hold LMB to dig / RMB to add,
// mouse wheel resizes the brush, F5 saves, F9 loads.
struct TerrainToolSettings {
    float brushRadius = 4.0f;    // world units
    float digStrength = 0.12f;   // density change per stroke at the brush centre (~0.5-1 m of rock)
    float strokesPerSecond = 10.0f;
    float reach = 12.0f;         // max raycast distance (arm's reach; fly mode uses 150)
    bool noAim = false;          // (#1071) no aim at all (no hit, no dig preview): scripted pictures that use no tool
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
    std::function<bool(class Planetoid&, Vector3 localGround, float radius)> till; // (#1085) a free brush: body-local radius
    // The item in hand became another (an empty bucket filled, a full one emptied): the game rebinds
    // the hotbar slot and counts quest events ("scoop", "pour").
    std::function<void(const std::string& from, const std::string& to, const char* event)> swapHeld;
};
extern TerrainToolSettings& terrainTool; // the local PlayerSession's (#313)
inline float gHudScale = 1.0f; // (#891) the interface scale for the crosshair's text (the game sets it each frame)
// Skills (#883): the game sets these (unset: no XP, every multiplier 1, as a server and tests run). A stroke's events go
// out (rock_stroke, ore_stroke, soil_stroke, gem, chop, tree_felled; `local`: the body-local point, for the anti-farming
// cell); the perks' multipliers come in (rock_speed, ore_speed, rock_radius, gem_chance, ore_yield, stone_yield; #919:
// soil_speed, soil_radius, chop_speed, log_yield, fibre_yield, stone_find).
inline std::function<void(const char* event, double n, Vector3 local)> gSkillEvent;
inline std::function<float(const char* stat)> gSkillStat;
inline float SkillStat(const char* stat) { return gSkillStat ? gSkillStat(stat) : 1.0f; }
inline std::function<float()> gTerrainStrokeScale; // (#1254) a potion of mining haste: strokes a second times this (unset: 1)
// (#919) n items times a yield stat: the whole part, and the fraction left by chance (the gameplay RNG): x1.1 on 3 logs
// gives 3 most times and 4 one time in three.
// by: whose stats (#920: a server's remote player); empty = the local player's (SkillStat).
int SkillScaled(int n, const char* stat, const std::function<float(const char*)>& by = nullptr);

struct PlayerInput;
void TerrainTool_Update(Scene& world, const PlayerInput& in, bool editEnabled = true); // once per frame before drawing (input: #312); editEnabled=false keeps save/load (F5/F9) only
bool TerrainTool_Stroke(bool dig);                         // one dig/build stroke at the last crosshair hit
bool TerrainTool_IsFiller(const std::string& item);        // bulk ground placed with RMB (dirt, gravel, sand, stone...)
float TerrainTool_Cooldown();                              // 1 just after a swing .. 0 ready for the next
float TerrainTool_StrokeRate();                            // strokes a second at the aimed ground with the tool in hand (#883: with the perks)
bool TerrainTool_Chop();     // one swing at the aimed plant (scripted tests); true when it came down
bool TerrainTool_Bucket(class Scene& world); // the bucket in hand: scoop or pour (RMB; scripted tests call it after aiming)
bool TerrainTool_BucketFromOffHand(class Scene& world); // (#836) the off hand's bucket where the player aims (the off-hand key)
bool TerrainTool_AimingPlant(uint64_t* id = nullptr);
bool TerrainTool_AimedGround(class Planetoid** body, Vector3* local); // the ground the tool aims at (density-marched: no chunks needed), body-local
float TerrainTool_ChopProgress(); // (#1094) the share of the aimed plant cut so far (0..1; 0: none being cut)
bool TerrainTool_LastHit(class Planetoid*& body, Vector3& local, Vector3& normalWorld); // the crosshair's terrain hit (placing torches, #218)
bool TerrainTool_TestFill(); // --test-fill (#228): an exact sphere of dirt, paid by volume
bool TerrainTool_TestFillOffHand(); // (#836) a tool in hand: RMB places the off hand's block stack
class Planetoid;
bool TerrainTool_DigAt(Planetoid& p, Vector3 local);              // tests: dig a point, gain what it held
bool TerrainTool_HarvestNearest(Planetoid& p, Vector3 local, bool tree, bool stone = false); // tests: fell / harvest the nearest plant (or pick up stones)
void TerrainTool_Draw3D();                                // inside BeginMode3D (Scene::onDraw3D)
void TerrainTool_ForgetWorld();                           // drop the hit and plant pointers into a world being unloaded (#902)
bool TerrainTool_TestVein(Planetoid& body, Vector3 homeDir); // --test-vein (#207): a deposit mined 93% by hand comes loose
void TerrainTool_DrawHUD();                               // 2D, after the scene

// Multiplayer (#329): what a dig stroke yields, with no side effects (the local stroke and the server's validated one
// both use it, so a client's inventory gets exactly what the server credits). `centre`/`r` the carve sphere, `hitLocal`
// the aimed surface point (materials as drawn there), all body-local. Gems and seeds roll on GameRng (the server's).
#include <map>
struct DigYield {
    std::map<std::string, float> items;
    // (#1097) Placed ground's share: 1.8 of its filler a m^3 (it cost 2), no extras and no perks. TerrainTool_ScaleYield
    // adds it to `items` after scaling them.
    std::map<std::string, float> placed;
    const char* gem = nullptr;   // a rare gem or quartz (+1)
    const char* seed = nullptr;  // grass turned up a seed (+1)
    bool ore = false;            // some ore came out (a vein may be nearly done)
};
DigYield TerrainTool_DigYield(Planetoid& p, Vector3 centre, float r, Vector3 hitLocal, bool softOnly, float gemBoost); // gemBoost: x the gem odds (tests, perks)
// (#1250) A rare find in rock other than a gem (magic: a resonant crystal on a ley line), asked after the gem roll when
// none came up: an item id or nullptr. Unset: nothing.
inline std::function<const char*(class Planetoid& body, Vector3 hitLocal)> gRockFind;
float TerrainTool_PlaceCost(const Planetoid& p, const struct TerrainEdit& e); // (#1097) a soft placing blob's filler: 2 a m^3 of ground it adds
// (#883, #1043) The perks on a stroke, from one player's stats (`stat`: a multiplier by name, 1 = none). The local stroke and
// the server's check of a joined player's stroke both use these, so a perked stroke is never refused and both credit alike.
// material: the VoxelMaterial hit; pickaxe: the tool in hand is a pick.
struct DigPerkScale { float radius = 1.0f, rate = 1.0f; };
DigPerkScale TerrainTool_DigPerkScale(bool pickaxe, int material, const std::function<float(const char*)>& stat);
void TerrainTool_ScaleYield(DigYield& got, float oreMul, float stoneMul); // Mining's ore_yield and stone_yield on a stroke's items
int TerrainTool_FillerMaterial(const std::string& item); // the VoxelMaterial a filler places (-1: not a filler)
const char* TerrainTool_ItemAt(Planetoid& p, Vector3 local, int material); // the item a material there digs out (ore: its kind)
// Set while joined to a server (#329): a dig or place stroke is sent as a request instead of applied and paid here;
// the server validates it, applies it in its order and credits or charges the inventory (deltas come back).
// dig: the carve sphere; place: the soft blob, `item` the filler that pays for it.
using TerrainToolNetEdit = std::function<void(Planetoid& p, bool dig, Vector3 hitLocal, Vector3 centre, float radius, const std::string& item)>;
extern TerrainToolNetEdit gTerrainToolNetEdit;

// Multiplayer (#335): what a felled plant (tree, bush, loose stone) gives, from its id and the biome it grew in, with no
// side effects (the local chop and the server's both use it).
#include <vector>
struct Plant;
struct ToolSpec;
struct PlantDrops {
    std::vector<std::pair<std::string, int>> items;
    std::string toast;  // "+4 oak logs", ...
    int sound = 0;      // 0 pickup, 1 a tree's (build) thud, 2 a stone
};
PlantDrops TerrainTool_PlantDrops(const Planetoid& p, const Plant& plant);
// (#919, #920) Woodcutting's and Foraging's perks on a felled plant's drops (log_yield, fibre_yield, stone_find), from
// `stat` (one player's stats; empty = the local player's).
void TerrainTool_ScalePlantDrops(PlantDrops& d, const Plant& plant, const std::function<float(const char*)>& stat = nullptr);
int TerrainTool_HitsToFell(const Plant& plant, const ToolSpec& tool);
// Set while joined (#335): a plant that comes down is a request (body, plant id); the server removes it for everyone and
// credits the drops (deltas), and its removal comes back (TerrainTool_PlantRemoved).
inline std::function<void(Planetoid& p, uint64_t plantId)> gPlantRequest;
// The server's removal of a plant (#335): gone here too, at its time (regrowth runs on the shared clock); `mine`: this
// client felled it (the toast and the sound of its drops).
void TerrainTool_PlantRemoved(Planetoid& p, uint64_t plantId, double when, bool mine);

#endif
