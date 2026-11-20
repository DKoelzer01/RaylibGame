#ifndef FARMING_H
#define FARMING_H

// Farming (plan 2, milestone 11). Farmland is a build part (1 m square of tilled soil: craft it
// from dirt and place it, or till soft ground with a shovel's RMB). Look at farmland and press I to
// plant seeds you carry, O to harvest a ripe crop. Crops (resources/data/crops.txt) grow through four
// stages in daylight only, twice as fast with water within 4 m; harvesting gives the crop and seeds
// back. Seeds come from harvesting bushes. Crop state is saved per plot part.

#include <cstdio>
#include <cstdint>
#include <string>
#include <raylib.h>

class Scene;
class Planetoid;
class Plot;

void Farm_Update(Scene& world, float worldDt);           // growth (world seconds, warped with time)
void Farm_Draw3D(Scene& world);                           // the crops (inside BeginMode3D)
// Aim helpers for the machine/farmland aim (main thread): plant the first seed carried / harvest.
bool Farm_Plant(Planetoid& body, Plot& plot, int part, const std::string& seed = "");
int Farm_Harvest(Planetoid& body, Plot& plot, int part);  // items gained (0: nothing ripe)
// What a farmland part holds: crop name ("" none) and growth 0..1.
bool Farm_Crop(const Planetoid& body, const Plot& plot, int part, std::string& crop, float& growth);
void Farm_DrawAimHUD(const Planetoid& body, const Plot& plot, int part);
// A shovel's RMB on soft ground: farmland there (no cost). Returns true if it tilled.
bool Farm_Till(Scene& world, Planetoid& body, Vector3 localGround, float radius = 0.6f); // (#1085) paint the soil there farmland
void Farm_Reset();

// Multiplayer (#335): the server owns the fields. A joined client's I / O / till are requests (set while joined); the
// server sows and harvests from the player's inventory as it keeps it (no toasts, sounds or quest counts there: the
// caller does those) and sends each client the fields near it (Farm_FieldsNear -> MsgFarmFields -> Farm_SetRemoteFields).
#include <functional>
#include <vector>
class Inventory;
struct NetFarmField;
inline std::function<void(Planetoid& body, const Plot& plot, int part, bool harvest, const std::string& seed)> gFarmRequest;
inline std::function<bool(Planetoid& body, Vector3 localGround, float radius)> gFarmTillRequest; // joined (#335, #1085): the server tills
bool Farm_SowFor(Inventory& inv, Planetoid& body, Plot& plot, int part, const std::string& seed, std::string* why = nullptr);
// stat: the harvester's stats for Farming's perks (#919; #920: a server's remote player); empty = no perks.
int Farm_HarvestFor(Inventory& inv, Planetoid& body, Plot& plot, int part, std::string* why = nullptr, const std::function<float(const char*)>& stat = nullptr);
void Farm_FieldsNear(const Planetoid& body, Vector3 local, float range, std::vector<NetFarmField>& out);
void Farm_SetRemoteFields(const Planetoid& body, const std::vector<NetFarmField>& fields); // replaces that body's fields

void Farm_Save(FILE* f);
bool Farm_LoadLine(const char* tok, FILE* f);

#endif
