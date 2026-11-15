#ifndef EVENTS_H
#define EVENTS_H

// World events (#177), first: meteor showers. At night, each 10 minutes of world time has a 15% chance (from a hash of
// the world seed and the time) of a shower: two minutes of streaks across the sky, and two meteors that come down
// 80-180 m from the player on dry land. Each impact flashes, throws up dust, carves a crater and leaves a hot,
// glowing meteorite (PartType::Meteorite) at its bottom: broken apart it gives nickel-iron and now and then platinum
// or iridium ore, otherwise found only in space.

#include <cstdio>

class Scene;
class Planetoid;

void Events_Update(Scene& world, float dt);
void Events_Draw3D(Scene& world);          // the streaks
void Events_AddLights(Scene& world);       // impact flashes and hot meteorites into world.localLights (#139)
void Events_Reset();
void Events_Save(FILE* f);
bool Events_LoadLine(const char* tok, FILE* f);
bool Events_ShowerOn();
void Events_ForceShower();                 // tests and debugging: a shower now
// Tests: a meteor comes down `dist` metres from the player on a bearing (radians): true if it landed (dry ground).
bool Events_TestImpact(Scene& world, Planetoid& body, float dist, float bearing);
int Events_Streaks();                      // streaks drawn so far (tests)

#endif
