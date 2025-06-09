#ifndef world_h
#define world_h

#include "scene.h"
#include "SimplexNoise.h"
#include "cubemarch.h"
#include "utils/logger.h"
#include "worldtypes.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <memory>
#include <unordered_map>
#include <tuple>
#include <functional>
#include <algorithm>


void worldHandler(Scene& world);
int generatePlanetoids(float randScale, Scene& world, float genRange, float minSize, float maxSize);
void generateWorld(Scene& world);
void loadWorld(Scene& world);

#endif