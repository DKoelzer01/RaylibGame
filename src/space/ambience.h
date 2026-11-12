#ifndef AMBIENCE_H
#define AMBIENCE_H

// Environmental particles around the player: rain splashes, fireflies on clear nights, smoke from
// working fire machines and campfire flames (ambience.cpp; particles.h does the drawing).

class Scene;
void Ambience_Update(Scene& world, float dt);

#endif
