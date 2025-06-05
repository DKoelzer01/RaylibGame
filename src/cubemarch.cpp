#include "cubemarch.h"

std::array<int,15> get_triangulation(int x, int y, int z, std::vector<float>& noiseValues,int size, float threshold) {
    int cubeIndex = 0;
    if (noiseValues[(x    ) + (y    ) * size + (z    ) * size * size] < threshold) cubeIndex |= 1 << 0; // (0,0,0)
    if (noiseValues[(x + 1) + (y    ) * size + (z    ) * size * size] < threshold) cubeIndex |= 1 << 1; // (1,0,0)
    if (noiseValues[(x + 1) + (y + 1) * size + (z    ) * size * size] < threshold) cubeIndex |= 1 << 2; // (1,1,0)
    if (noiseValues[(x    ) + (y + 1) * size + (z    ) * size * size] < threshold) cubeIndex |= 1 << 3; // (0,1,0)
    if (noiseValues[(x    ) + (y    ) * size + (z + 1) * size * size] < threshold) cubeIndex |= 1 << 4; // (0,0,1)
    if (noiseValues[(x + 1) + (y    ) * size + (z + 1) * size * size] < threshold) cubeIndex |= 1 << 5; // (1,0,1)
    if (noiseValues[(x + 1) + (y + 1) * size + (z + 1) * size * size] < threshold) cubeIndex |= 1 << 6; // (1,1,1)
    if (noiseValues[(x    ) + (y + 1) * size + (z + 1) * size * size] < threshold) cubeIndex |= 1 << 7; // (0,1,1)

    if (cubeIndex == 0 || cubeIndex == 255) {
        // No triangles to create
        return std::array<int, 15>{-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1};
    }

    std::array<int, 15> triangulation;
    std::copy(std::begin(TRIANGULATIONS[cubeIndex]), std::end(TRIANGULATIONS[cubeIndex]), triangulation.begin());
    return triangulation;
}

// Call this for each cube
void marchCube(
    int x, int y, int z,
    std::vector<float>& noiseValues, int size,
    std::vector<Vector3>& vertices,
    std::vector<unsigned int>& indices,
    std::unordered_map<EdgeKey, unsigned int>& edgeCache,
    float threshold
) {
    std::array<int, 15> triangulation = get_triangulation(x, y, z, noiseValues, size, threshold);

    for (int i = 0; i < 15 && triangulation[i] != -1; i += 3) {
        unsigned int triIdx[3];
        for (int v = 0; v < 3; ++v) { // <--- loop forward!
            int edgeIndex = triangulation[i + v];
            if (edgeIndex < 0 || edgeIndex >= 12) continue;

            int nx = x, ny = y, nz = z, nedge = edgeIndex;
            const int* neighbor = EDGE_NEIGHBORS[edgeIndex];
            if (neighbor[0] < 0 && x > 0) {
                nx += neighbor[0];
                nedge = neighbor[3];
            } else if (neighbor[1] < 0 && y > 0) {
                ny += neighbor[1];
                nedge = neighbor[3];
            } else if (neighbor[2] < 0 && z > 0) {
                nz += neighbor[2];
                nedge = neighbor[3];
            }

            EdgeKey key{nx, ny, nz, nedge};
            auto it = edgeCache.find(key);
            if (it != edgeCache.end()) {
                triIdx[v] = it->second;
            } else {
                int v0i = CUBE_EDGES[edgeIndex][0];
                int v1i = CUBE_EDGES[edgeIndex][1];
                int x0 = CUBE_VERTICES[v0i][0], y0 = CUBE_VERTICES[v0i][1], z0 = CUBE_VERTICES[v0i][2];
                int x1 = CUBE_VERTICES[v1i][0], y1 = CUBE_VERTICES[v1i][1], z1 = CUBE_VERTICES[v1i][2];

                float v0 = noiseValues[(x + x0) + (y + y0) * size + (z + z0) * size * size];
                float v1 = noiseValues[(x + x1) + (y + y1) * size + (z + z1) * size * size];
                float denom = v1 - v0;
                if (fabs(denom) < 1e-6f) denom = 1e-6f;
                float t = (threshold - v0) / denom;
                t = fmaxf(0.0f, fminf(1.0f, t));

                Vector3 pos_a = { static_cast<float>(x + x0), static_cast<float>(y + y0), static_cast<float>(z + z0) };
                Vector3 pos_b = { static_cast<float>(x + x1), static_cast<float>(y + y1), static_cast<float>(z + z1) };
                Vector3 position = {
                    pos_a.x + t * (pos_b.x - pos_a.x),
                    pos_a.y + t * (pos_b.y - pos_a.y),
                    pos_a.z + t * (pos_b.z - pos_a.z)
                };
                vertices.push_back(position);
                unsigned int idx = static_cast<unsigned int>(vertices.size() - 1);
                edgeCache[key] = idx;
                triIdx[v] = idx;
            }
        }
        // Reverse winding order here:
        indices.push_back(triIdx[0]);
        indices.push_back(triIdx[1]);
        indices.push_back(triIdx[2]);
    }
}