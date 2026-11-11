#ifndef ANIMALS_H
#define ANIMALS_H

// Animals (plan 2, milestone 12): passive placeholder animals - coloured boxes with a head - from
// resources/data/animals.txt. Herds are placed from the world seed in 256 m tiles near the player
// (once seen, an animal is saved); only animals within 300 m are simulated.
//
// Behaviour: wander and graze, rest at night, flee when hit (wild ones also when the player comes
// close). They walk on the terrain surface, avoid water and steep ground, and treat build parts as
// walls, so a ring of blocks or beams is a pen. Look at one: I feeds it its food (it then follows
// you for 30 s - lead it into a pen; two fed adults near each other breed a young one that grows up
// over a day), O collects its product (wool with shears, milk with a bucket, eggs), LMB attacks
// (meat and hide when it falls).

#include <cstdio>
#include <string>
#include <raylib.h>

class Scene;
class Planetoid;

void Animals_Update(Scene& world, float dt, float worldDt, bool inputEnabled);
void Animals_Draw3D(Scene& world);
void Animals_DrawHUD();
bool Animals_Aiming();                    // the crosshair is on an animal (the terrain tool holds off)
void Animals_Reset();
void Animals_Save(FILE* f);
bool Animals_LoadLine(const char* tok, FILE* f);

// Tests.
int Animals_Count();
int Animals_Spawn(Scene& world, Planetoid& body, const std::string& species, Vector3 bodyLocal, int count, bool young = false, float spread = 18.0f);
int Animals_Nearest(Vector3 bodyLocal, float range);          // index or -1
Vector3 Animals_Position(int index);                           // body-local
bool Animals_Hit(int index, float damage);                     // true if it fell
bool Animals_Feed(int index);                                  // with the player's food for it
int Animals_YoungCount();
double Animals_LastUpdateMs();
const std::string& Animals_Body(int index);                   // the body's name
bool Animals_Wild(int index);
bool Animals_Hostile(int index);
bool Animals_IsFeed(const std::string& item);                   // some animal eats it
void Animals_SetThreats(bool on);
void Animals_RefreshNow();                                     // tests: rescan campfires, allow a night spawn, recheck pens
int Animals_NightCount();                              // hostile behaviour and night spawns (bench runs: off)
int Animals_PennedCount();                                   // near the player, last checked (every 2 s)

#endif
