#include "edgecache.h"
#include "cubemarch.h"
#include "world.h" // Include world.h to bring in CHUNK_SIZE definition
#include "worldtypes.h"

extern const int CHUNK_SIZE; // Declare CHUNK_SIZE as extern

// Canonical edge remapping tables for each face direction
// For each face, maps local edge index to canonical edge index (for both + and - sides)
// edgeRemap[faceDir][side][edge]
// faceDir: 0=X, 1=Y, 2=Z; side: 0=+1, 1=-1
static const int edgeRemap[3][2][12] = {
    // X face
    {   // +X (side=+1)
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
        // -X (side=-1)
        {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10}
    },
    // Y face
    {   // +Y
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
        // -Y
        {3, 2, 1, 0, 7, 6, 5, 4, 11, 10, 9, 8}
    },
    // Z face
    {   // +Z
        {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
        // -Z
        {4, 5, 6, 7, 0, 1, 2, 3, 8, 9, 10, 11}
    }
};

void getCanonicalFaceEdgeIndex(
    int faceDir, int side, // 0=X, 1=Y, 2=Z; side=+1 or -1
    int x, int y, int z, int edge,
    int& a, int& b, int& canonicalEdge
) {
    // For each face, a/b are the two axes orthogonal to the face
    // X face: a=y, b=z
    // Y face: a=x, b=z
    // Z face: a=x, b=y
    if (faceDir == 0) { // X
        a = y;
        b = z;
    } else if (faceDir == 1) { // Y
        a = x;
        b = z;
    } else { // Z
        a = x;
        b = y;
    }
    int sideIdx = (side == +1) ? 0 : 1;
    canonicalEdge = edgeRemap[faceDir][sideIdx][edge];
}

int getFaceFlatIdx(int a, int b, int edge, int stride) {
    return (a + b * stride) * 12 + edge;
}

std::array<int,15> get_triangulation(int x, int y, int z, std::vector<float>& noiseValues, int size, float threshold) {
    int cubeIndex = 0;
    // Bounds check for all 8 cube corners
    assert(x + 1 < size && y + 1 < size && z + 1 < size);

    auto idx = [size](int x, int y, int z) -> int {
        return x + y * size + z * size * size;
    };

    if (noiseValues[idx(x, y, z)] < threshold) cubeIndex |= 1 << 0; // (0,0,0)
    if (noiseValues[idx(x + 1, y, z)] < threshold) cubeIndex |= 1 << 1; // (1,0,0)
    if (noiseValues[idx(x + 1, y + 1, z)] < threshold) cubeIndex |= 1 << 2; // (1,1,0)
    if (noiseValues[idx(x, y + 1, z)] < threshold) cubeIndex |= 1 << 3; // (0,1,0)
    if (noiseValues[idx(x, y, z + 1)] < threshold) cubeIndex |= 1 << 4; // (0,0,1)
    if (noiseValues[idx(x + 1, y, z + 1)] < threshold) cubeIndex |= 1 << 5; // (1,0,1)
    if (noiseValues[idx(x + 1, y + 1, z + 1)] < threshold) cubeIndex |= 1 << 6; // (1,1,1)
    if (noiseValues[idx(x, y + 1, z + 1)] < threshold) cubeIndex |= 1 << 7; // (0,1,1)

    if (cubeIndex == 0 || cubeIndex == 255) {
        // No triangles to create
        return std::array<int, 15>{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    }

    std::array<int, 15> triangulation;
    std::copy(std::begin(TRIANGULATIONS[cubeIndex]), std::end(TRIANGULATIONS[cubeIndex]), triangulation.begin());
    return triangulation;
}

// Helper to sample noise at floating-point position
// Manual clamp implementation for C++11/14 compatibility
template<typename T>
T clamp(T v, T lo, T hi) {
    if (v < lo) return lo;
    if (v > hi) return hi;
    return v;
}

// Call this for each cube
void marchCube(
    int x, int y, int z,
    std::vector<float>& noiseValues, int size,
    std::vector<Vector3>& vertices,
    std::vector<int>& indices,
    std::vector<EdgeCacheEntry>* edgeCacheLocal,
    std::vector<EdgeCacheEntry>* edgeCacheX,
    std::vector<EdgeCacheEntry>* edgeCacheY,
    std::vector<EdgeCacheEntry>* edgeCacheZ,
    float threshold,
    int cx, int cy, int cz,
    const Vector3& chunkOrigin,
    SimplexNoise* simplexNoise,
    Planetoid* planetoid
) {
    assert(x + 1 < size && y + 1 < size && z + 1 < size);
    std::array<int, 15> triangulation = get_triangulation(x, y, z, noiseValues, size, threshold);
    for (int i = 0; i < 15 && triangulation[i] != -1; i += 3) {
        size_t triIdx[3];
        for (int v = 0; v < 3; ++v) {
            int edgeIndex = triangulation[i + v];
            if (edgeIndex < 0 || edgeIndex >= 12) continue;
            int v0i = CUBE_EDGES[edgeIndex][0];
            int v1i = CUBE_EDGES[edgeIndex][1];
            int gx0 = (int)x + CUBE_VERTICES[v0i][0] + cx * (int)(size-1);
            int gy0 = (int)y + CUBE_VERTICES[v0i][1] + cy * (int)(size-1);
            int gz0 = (int)z + CUBE_VERTICES[v0i][2] + cz * (int)(size-1);
            bool onX = (gx0 % (int)(size-1) == 0) && (cx != 0) && edgeCacheX;
            bool onY = (gy0 % (int)(size-1) == 0) && (cy != 0) && edgeCacheY;
            bool onZ = (gz0 % (int)(size-1) == 0) && (cz != 0) && edgeCacheZ;
            std::vector<EdgeCacheEntry>* cache = edgeCacheLocal;
            size_t flatIdx = 0;
            int a=0, b=0, canonicalEdge=0;
            if (onX) {
                // For +X face, use y and z as local face coordinates
                int faceY = y;
                int faceZ = z;
                getCanonicalFaceEdgeIndex(0, +1, (int)(size-2), faceY, faceZ, edgeIndex, a, b, canonicalEdge);
                flatIdx = getFaceFlatIdx(a, b, canonicalEdge, (int)(size-1));
                cache = edgeCacheX;
            } else if (onY) {
                // For +Y face, use x and z as local face coordinates
                int faceX = x;
                int faceZ = z;
                getCanonicalFaceEdgeIndex(1, +1, faceX, (int)(size-2), faceZ, edgeIndex, a, b, canonicalEdge);
                flatIdx = getFaceFlatIdx(a, b, canonicalEdge, (int)(size-1));
                cache = edgeCacheY;
            } else if (onZ) {
                // For +Z face, use x and y as local face coordinates
                int faceX = x;
                int faceY = y;
                getCanonicalFaceEdgeIndex(2, +1, faceX, faceY, (int)(size-2), edgeIndex, a, b, canonicalEdge);
                flatIdx = getFaceFlatIdx(a, b, canonicalEdge, (int)(size-1));
                cache = edgeCacheZ;
            } else {
                flatIdx = ((z * (size-1) + y) * (size-1) + x) * 12 + edgeIndex;
                cache = edgeCacheLocal;
            }
            int cachedIdx = -1;
            if (flatIdx < cache->size() && (*cache)[flatIdx].vertexIdx != -1) {
                cachedIdx = (*cache)[flatIdx].vertexIdx;
            }
            if (cachedIdx != -1) {
                triIdx[v] = static_cast<size_t>(cachedIdx);
                // No normals.push_back or normal computation here!
            } else {
                Vector3 pos_a = { static_cast<float>(x + CUBE_VERTICES[v0i][0]), static_cast<float>(y + CUBE_VERTICES[v0i][1]), static_cast<float>(z + CUBE_VERTICES[v0i][2]) };
                Vector3 pos_b = { static_cast<float>(x + CUBE_VERTICES[v1i][0]), static_cast<float>(y + CUBE_VERTICES[v1i][1]), static_cast<float>(z + CUBE_VERTICES[v1i][2]) };
                auto idx_local = [size](int x, int y, int z) -> int {
                    return x + y * size + z * size * size;
                };
                float val_a = noiseValues[idx_local(static_cast<int>(pos_a.x), static_cast<int>(pos_a.y), static_cast<int>(pos_a.z))];
                float val_b = noiseValues[idx_local(static_cast<int>(pos_b.x), static_cast<int>(pos_b.y), static_cast<int>(pos_b.z))];
                float t = (threshold - val_a) / (val_b - val_a + 1e-8f);
                t = clamp(t, 0.0f, 1.0f);
                Vector3 position = {
                    pos_a.x + t * (pos_b.x - pos_a.x),
                    pos_a.y + t * (pos_b.y - pos_a.y),
                    pos_a.z + t * (pos_b.z - pos_a.z)
                };
                // Store mesh vertices in chunk-local space (no chunkOrigin offset)
                vertices.push_back(position);
                size_t idx = vertices.size() - 1;
                if (flatIdx < cache->size()) {
                    (*cache)[flatIdx].vertexIdx = static_cast<int>(idx);
                }
                triIdx[v] = idx;
            }
        }
        indices.push_back(triIdx[0]);
        indices.push_back(triIdx[1]);
        indices.push_back(triIdx[2]);
    }
    // No normals generated here; handled by Chunk::calculateNormals
}