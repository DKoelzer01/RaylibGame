#include "cubemarch.h"

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

std::array<__int8,15> get_triangulation(unsigned __int8 x, unsigned __int8 y, unsigned __int8 z, std::vector<float>& noiseValues, unsigned __int8 size, float threshold) {
    unsigned __int8 cubeIndex = 0;
    // Bounds check for all 8 cube corners
    assert(x + 1 < size && y + 1 < size && z + 1 < size);

    auto idx = [size](unsigned __int8 x, unsigned __int8 y, unsigned __int8 z) -> unsigned __int16 {
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
        return std::array<__int8, 15>{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    }

    std::array<__int8, 15> triangulation;
    std::copy(std::begin(TRIANGULATIONS[cubeIndex]), std::end(TRIANGULATIONS[cubeIndex]), triangulation.begin());
    return triangulation;
}

// Call this for each cube
void marchCube(
    unsigned __int8 x, unsigned __int8 y, unsigned __int8 z,
    std::vector<float>& noiseValues, int size,
    std::vector<Vector3>& vertices,
    std::vector<size_t>& indices,
    std::vector<int>* edgeCacheLocal,
    std::vector<int>* edgeCacheX,
    std::vector<int>* edgeCacheY,
    std::vector<int>* edgeCacheZ,
    float threshold,
    int cx, int cy, int cz, int chunkGrid
) {
    assert(x + 1 < size && y + 1 < size && z + 1 < size);
    std::array<__int8, 15> triangulation = get_triangulation(x, y, z, noiseValues, size, threshold);
    for (int i = 0; i < 15 && triangulation[i] != -1; i += 3) {
        size_t triIdx[3];
        for (int v = 0; v < 3; ++v) {
            int edgeIndex = triangulation[i + v];
            if (edgeIndex < 0 || edgeIndex >= 12) continue;
            // Compute global (chunk-space) coordinates for the edge
            int v0i = CUBE_EDGES[edgeIndex][0];
            int v1i = CUBE_EDGES[edgeIndex][1];
            int gx0 = (int)x + CUBE_VERTICES[v0i][0] + cx * (int)(size-1);
            int gy0 = (int)y + CUBE_VERTICES[v0i][1] + cy * (int)(size-1);
            int gz0 = (int)z + CUBE_VERTICES[v0i][2] + cz * (int)(size-1);
            int gx1 = (int)x + CUBE_VERTICES[v1i][0] + cx * (int)(size-1);
            int gy1 = (int)y + CUBE_VERTICES[v1i][1] + cy * (int)(size-1);
            int gz1 = (int)z + CUBE_VERTICES[v1i][2] + cz * (int)(size-1);
            // The edge is between (gx0,gy0,gz0) and (gx1,gy1,gz1)
            // Canonical owner: the chunk with the minimum (cx,cy,cz) that touches the edge
            // For each axis, if the edge lies on a chunk border, use the cache for the lowest chunk index
            bool onX = (gx0 % (int)(size-1) == 0) && edgeCacheX;
            bool onY = (gy0 % (int)(size-1) == 0) && edgeCacheY;
            bool onZ = (gz0 % (int)(size-1) == 0) && edgeCacheZ;
            std::vector<int>* cache = edgeCacheLocal;
            size_t flatIdx = 0;
            int a=0, b=0, canonicalEdge=0;
            // Priority: X > Y > Z for corners/edges
            if (onX) {
                // Use X face cache, canonicalize to face at min X
                int faceX = gx0 / (int)(size-1);
                int localY = gy0 % (int)(size-1);
                int localZ = gz0 % (int)(size-1);
                getCanonicalFaceEdgeIndex(0, +1, (int)(size-2), localY, localZ, edgeIndex, a, b, canonicalEdge);
                flatIdx = getFaceFlatIdx(a, b, canonicalEdge, (int)(size-1));
                cache = edgeCacheX;
            } else if (onY) {
                int faceY = gy0 / (int)(size-1);
                int localX = gx0 % (int)(size-1);
                int localZ = gz0 % (int)(size-1);
                getCanonicalFaceEdgeIndex(1, +1, localX, (int)(size-2), localZ, edgeIndex, a, b, canonicalEdge);
                flatIdx = getFaceFlatIdx(a, b, canonicalEdge, (int)(size-1));
                cache = edgeCacheY;
            } else if (onZ) {
                int faceZ = gz0 / (int)(size-1);
                int localX = gx0 % (int)(size-1);
                int localY = gy0 % (int)(size-1);
                getCanonicalFaceEdgeIndex(2, +1, localX, localY, (int)(size-2), edgeIndex, a, b, canonicalEdge);
                flatIdx = getFaceFlatIdx(a, b, canonicalEdge, (int)(size-1));
                cache = edgeCacheZ;
            } else {
                // Local cache for interior edges
                flatIdx = ((z * (size-1) + y) * (size-1) + x) * 12 + edgeIndex;
                cache = edgeCacheLocal;
            }
            int cachedIdx = -1;
            if (flatIdx < cache->size()) cachedIdx = (*cache)[flatIdx];
            if (cachedIdx != -1) {
                triIdx[v] = static_cast<size_t>(cachedIdx);
            } else {
                // Compute vertex position as before
                size_t idx0 = (x + CUBE_VERTICES[v0i][0]) + (y + CUBE_VERTICES[v0i][1]) * size + (z + CUBE_VERTICES[v0i][2]) * size * size;
                size_t idx1 = (x + CUBE_VERTICES[v1i][0]) + (y + CUBE_VERTICES[v1i][1]) * size + (z + CUBE_VERTICES[v1i][2]) * size * size;
                assert(idx0 < noiseValues.size());
                assert(idx1 < noiseValues.size());
                float v0 = noiseValues[idx0];
                float v1 = noiseValues[idx1];
                float denom = v1 - v0;
                if (fabs(denom) < 1e-6f) denom = 1e-6f;
                float t = (threshold - v0) / denom;
                t = fmaxf(0.0f, fminf(1.0f, t));
                Vector3 pos_a = { static_cast<float>(x + CUBE_VERTICES[v0i][0]), static_cast<float>(y + CUBE_VERTICES[v0i][1]), static_cast<float>(z + CUBE_VERTICES[v0i][2]) };
                Vector3 pos_b = { static_cast<float>(x + CUBE_VERTICES[v1i][0]), static_cast<float>(y + CUBE_VERTICES[v1i][1]), static_cast<float>(z + CUBE_VERTICES[v1i][2]) };
                Vector3 position = {
                    pos_a.x + t * (pos_b.x - pos_a.x),
                    pos_a.y + t * (pos_b.y - pos_a.y),
                    pos_a.z + t * (pos_b.z - pos_a.z)
                };
                vertices.push_back(position);
                size_t idx = vertices.size() - 1;
                if (flatIdx < cache->size()) (*cache)[flatIdx] = static_cast<int>(idx);
                triIdx[v] = idx;
            }
        }
        indices.push_back(triIdx[0]);
        indices.push_back(triIdx[1]);
        indices.push_back(triIdx[2]);
    }
}