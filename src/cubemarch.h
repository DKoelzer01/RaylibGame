#ifndef cubemarch_h
#define cubemarch_h

#include "cubemarchconsts.h"
#include "utils/logger.h"
#include "SimplexNoise.h"

#include <raylib.h>
#include <vector>
#include <array>
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <tuple>
#include <cassert>

struct EdgeCacheEntry {
    int vertexIdx;
    Vector3 normal;
    EdgeCacheEntry() : vertexIdx(-1), normal({0,0,0}) {}
};

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

std::array<int,15> get_triangulation(int x, int y, int z, std::vector<float>& noiseValues,int size, float threshold);

void marchCube(
    int x, int y, int z,
    std::vector<float>& noiseValues, int size,
    std::vector<Vector3>& vertices,
    std::vector<int>& indices,
    std::vector<Vector3>& normals,
    std::vector<EdgeCacheEntry>* edgeCacheLocal,
    std::vector<EdgeCacheEntry>* edgeCacheX,
    std::vector<EdgeCacheEntry>* edgeCacheY,
    std::vector<EdgeCacheEntry>* edgeCacheZ,
    float threshold,
    int cx, int cy, int cz,
    const Vector3& chunkOrigin,
    SimplexNoise* simplexNoise
);

extern std::vector<int> edgeCacheFlat;

// Canonical mapping for shared face edge cache
// faceDir: 0=X, 1=Y, 2=Z
// side: +1 for positive face (+X, +Y, +Z), -1 for negative face (-X, -Y, -Z)
void getCanonicalFaceEdgeIndex(
    int faceDir, int side, // 0=X, 1=Y, 2=Z; side=+1 or -1
    int x, int y, int z, int edge,
    int& a, int& b, int& canonicalEdge
);
// Returns flatIdx for a face, given (a, b, edge) and stride
int getFaceFlatIdx(int a, int b, int edge, int stride);

// Add this struct for edge cache entries with normals


#endif