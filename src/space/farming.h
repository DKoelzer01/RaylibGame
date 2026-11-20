#ifndef FARMING_H
#define FARMING_H

// Farming (plan 2, milestone 11; #1085). Farmland is a state of dirt: a shovel's or hoe's RMB tills the soil where it aims
// (a free brush that paints the ground Farmland, grass and all). Look at farmland and press I to plant seeds you carry, O
// to harvest a ripe crop. A crop stands where it was planted (body-local, on a 5 cm lattice), at least its spacing from
// the next (crops.txt); it grows while the ground under it is still farmland. Crops (resources/data/crops.txt) grow
// through four stages in daylight only, twice as fast with water within 4 m; harvesting gives the crop and seeds back.
// Seeds come from harvesting bushes.

#include <cstdio>
#include <cstdint>
#include <string>
#include <raylib.h>

class Scene;
class Planetoid;

void Farm_Update(Scene& world, float worldDt);           // growth (world seconds, warped with time)
void Farm_Draw3D(Scene& world);                           // the crops (inside BeginMode3D)
// The aim (main thread), at a body-local ground point: plant the first seed carried there / harvest the crop there.
bool Farm_Plant(Planetoid& body, Vector3 localGround, const std::string& seed = "");
int Farm_Harvest(Planetoid& body, Vector3 localGround);  // items gained (0: nothing ripe)
// The crop at the aimed ground (within 0.4 m): its name ("" none) and growth 0..1.
bool Farm_Crop(const Planetoid& body, Vector3 localGround, std::string& crop, float& growth);
bool Farm_Aims(const Planetoid& body, Vector3 localGround); // farmland or a crop there: the farm aim's
void Farm_DrawAimHUD(const Planetoid& body, Vector3 localGround);
int Farm_CropCount();                                     // crops growing (all bodies; tests)
// A shovel's or hoe's RMB on soil: farmland there (no cost). Returns true if it tilled.
bool Farm_Till(Scene& world, Planetoid& body, Vector3 localGround, float radius = 0.6f);
void Farm_Reset();

// Multiplayer (#335): the server owns the fields. A joined client's I / O / till are requests (set while joined); the
// server sows and harvests from the player's inventory as it keeps it (no toasts, sounds or quest counts there: the
// caller does those) and sends each client the fields near it (Farm_FieldsNear -> MsgFarmFields -> Farm_SetRemoteFields).
#include <functional>
#include <vector>
class Inventory;
struct NetFarmField;
inline std::function<void(Planetoid& body, Vector3 localGround, bool harvest, const std::string& seed)> gFarmRequest;
inline std::function<bool(Planetoid& body, Vector3 localGround, float radius)> gFarmTillRequest; // joined (#335, #1085): the server tills
bool Farm_SowFor(Inventory& inv, Planetoid& body, Vector3 localGround, const std::string& seed, std::string* why = nullptr);
// stat: the harvester's stats for Farming's perks (#919; #920: a server's remote player); empty = no perks.
int Farm_HarvestFor(Inventory& inv, Planetoid& body, Vector3 localGround, std::string* why = nullptr, const std::function<float(const char*)>& stat = nullptr);
void Farm_FieldsNear(const Planetoid& body, Vector3 local, float range, std::vector<NetFarmField>& out);
void Farm_SetRemoteFields(const Planetoid& body, const std::vector<NetFarmField>& fields); // replaces that body's fields

void Farm_Save(FILE* f);
bool Farm_LoadLine(const char* tok, FILE* f);
// After a load (bodies, terrain edits and plots in): an old save's farmland parts become tilled ground (#1085).
void Farm_LoadDone();

#endif
