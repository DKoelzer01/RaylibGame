#ifndef world_h
#define world_h

#include "scene.h"
#include <sstream>
#include <random>
#include "SimplexNoise.h"
#include "cubemarch.h"

void worldHandler(Scene& world);
void generateWorld(Scene& world);
void loadWorld(Scene& world);

#endif