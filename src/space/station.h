#ifndef STATION_H
#define STATION_H

// Space stations (plan 3, milestone 20): bodies that carry build plots on a fixed orbit. A station
// is founded when a rocket in orbit deploys its station core (rocket.cpp); from then on its body moves
// on an exact two-body (Kepler) orbit around the planet it was founded over, from the state it had
// then. Rockets dock to it (their plot joins the station's body) and undock (they launch from it).
// Stations are saved with the world.

#include <cstdint>
#include <cstdio>
#include <raylib.h>
#include "dvec3.h"

class Scene;
class Planetoid;

// A new station body at world position `pos`, moving at `vel` (world), orbiting `parent`. Returns
// the body (owned by the scene); it carries no plots yet.
Planetoid* Station_Create(Scene& world, Planetoid& parent, DVec3 pos, Vector3 vel, Quaternion rot, float extent);
void Station_Init(Scene* world);                // once at start (loads need the scene)
void Station_Update(Scene& world);             // before the player and the rocket each frame
bool Station_Is(const Planetoid& body);
int Station_Count();
Planetoid* Station_Nearest(Scene& world, DVec3 pos, float range); // a station body within range
void Station_Reset();                          // a new world / a load: remove every station
void Station_Save(FILE* f);
bool Station_LoadLine(const char* tok, FILE* f);
// Tests: the orbit's apoapsis/periapsis altitudes now (m above the parent's radius).
bool Station_Orbit(int index, double& apoapsis, double& periapsis);

#endif
