#ifndef FISHING_H
#define FISHING_H

// Fishing (#176): with a fishing rod in hand, RMB on water within reach casts a bobber; after 4-12 s a fish bites
// (a splash, a sound, "Bite!") and RMB within 1.5 s lands it. Reeling in early catches nothing; walking more than 4 m
// off or putting the rod away loses the line. The catch depends on the water and the biome: cod in the sea, trout in
// rivers and lakes, catfish in warm waters (jungle, swamp, savanna). All cook into cooked fish at a campfire.

#include <raylib.h>
#include <string>

class Scene;
class Planetoid;

void Fishing_Update(Scene& world, float dt, bool inputEnabled);
void Fishing_Draw3D(Scene& world);
void Fishing_DrawHUD();
void Fishing_Reset();
// What a catch at this water point would be (body-local point at the water): "raw_cod", "raw_trout" or "raw_catfish".
std::string Fishing_CatchAt(const Planetoid& body, Vector3 local);
// Tests (--test-fishing): cast at a body-local water point, wait out the bite with simulated time, land it. Returns
// the item caught ("" for none); `early` reels in before the bite (nothing caught).
std::string Fishing_TestCast(Scene& world, Planetoid& body, Vector3 waterLocal, bool early);

#endif
