#ifndef AMBIENCE_H
#define AMBIENCE_H

// Environmental particles around the player: rain splashes, fireflies on clear nights, smoke from
// working fire machines and campfire flames (ambience.cpp; particles.h does the drawing).

#include <raylib.h>
class Scene;
void Ambience_Update(Scene& world, float dt);
// Oil seeps (Plan 5 B) within range of a body-local point: their surface points (body-local).
#include <vector>
class Planetoid;
int Ambience_SeepsNear(const Planetoid& body, Vector3 local, float range, std::vector<Vector3>& out);
void Ambience_Reset();
int Ambience_SpringsNear(const Planetoid& body, Vector3 local, float range, std::vector<Vector3>& out); // water springs (Plan 5)

#endif
