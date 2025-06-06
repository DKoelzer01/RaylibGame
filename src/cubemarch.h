#ifndef cubemarch_h
#define cubemarch_h

#include "cubemarchconsts.h"
#include "utils/logger.h"

#include <raylib.h>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <cassert>

struct EdgeKey {
    int x, y, z, edge;
    bool operator==(const EdgeKey& other) const {
        return x == other.x && y == other.y && z == other.z && edge == other.edge;
    }
};

namespace std {
    template <>
    struct hash<EdgeKey> {
        std::size_t operator()(const EdgeKey& k) const {
            return ((std::hash<int>()(k.x) ^ (std::hash<int>()(k.y) << 1)) >> 1)
                 ^ (std::hash<int>()(k.z) << 1) ^ (std::hash<int>()(k.edge) << 2);
        }
    };
}

std::array<int,15> get_triangulation(size_t x, size_t y, size_t z, std::vector<float>& noiseValues,size_t size, float threshold);

void marchCube(
    size_t x, size_t y, size_t z,
    std::vector<float>& noiseValues, size_t size,
    std::vector<Vector3>& vertices,
    std::vector<size_t>& indices,
    std::vector<int>& edgeCacheFlat,
    float threshold,
    int cx, int cy, int cz, int chunkGrid,
    int cacheType
);

extern std::vector<int> edgeCacheFlat;

#endif