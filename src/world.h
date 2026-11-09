#ifndef world_h
#define world_h

#include "scene.h"
#include "SimplexNoise.h"
#include "cubemarch.h"
#include "utils/logger.h"
#include "worldtypes.h"
#include "profiler.h"

#include <vector>
#include <string>
#include <queue>
#include <unordered_set>
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

// Tweakable from the debug UI. seed == 0 means "pick a random seed"; lastSeed records what was
// actually used so an interesting world can be reproduced.
struct WorldGenSettings {
    unsigned int seed = 0;
    unsigned int lastSeed = 0;
    int planetoidCount = 1;
    int planetoidSize = 50;
    Vector3 firstPosition = { 64.0f, 0.0f, 0.0f };
};
extern WorldGenSettings worldGen;

void worldHandler(Scene& world);
int generatePlanetoids(float randScale, Scene& world, float genRange, float minSize, float maxSize);
void generateWorld(Scene& world);
void regenerateWorld(Scene& world);
void loadWorld(Scene& world);

#endif
