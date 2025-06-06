#ifndef world_h
#define world_h

#include "scene.h"
#include "SimplexNoise.h"
#include "cubemarch.h"
#include "utils/logger.h"

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
#include <unordered_map>
#include <functional>

struct Int3 {
    int x, y, z;
    bool operator==(const Int3& other) const { return x == other.x && y == other.y && z == other.z; }
};

struct Chunk {
    std::vector<float> noiseValues;
    std::vector<Vector3> vertices;
    std::vector<size_t> indices;
    Int3 position; // Position of the chunk in the world
    Mesh mesh;     // Store mesh for chunk lifetime
    Model model;   // Store model for chunk lifetime
    // Add more as needed (e.g., mesh, cache)
};

void worldHandler(Scene& world);
int generatePlanetoids(float randScale, Scene& world, float genRange, float minSize, float maxSize);
void generateWorld(Scene& world);
void loadWorld(Scene& world);

#endif