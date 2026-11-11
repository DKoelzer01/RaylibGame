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
bool Farm_Till(Scene& world, Planetoid& body, Vector3 localGround);
void Farm_Reset();
void Farm_Save(FILE* f);
bool Farm_LoadLine(const char* tok, FILE* f);

#endif
