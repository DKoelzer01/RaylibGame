#include "cubemarch.h"

// Definition for the global edgeCacheFlat
std::vector<int> edgeCacheFlat;

std::array<int,15> get_triangulation(size_t x, size_t y, size_t z, std::vector<float>& noiseValues, size_t size, float threshold) {
    int cubeIndex = 0;
    // Bounds check for all 8 cube corners
    assert(x + 1 < size && y + 1 < size && z + 1 < size);

    auto idx = [size](size_t x, size_t y, size_t z) -> size_t {
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

// Call this for each cube
void marchCube(
    size_t x, size_t y, size_t z,
    std::vector<float>& noiseValues, size_t size,
    std::vector<Vector3>& vertices,
    std::vector<size_t>& indices,
    std::unordered_map<EdgeKey, size_t>& edgeCache,
    float threshold
) {
    // Bounds check for cube
    assert(x + 1 < size && y + 1 < size && z + 1 < size);

    std::array<int, 15> triangulation = get_triangulation(x, y, z, noiseValues, size, threshold);

    // Flat array cache: size*size*size*12, initialized to -1 elsewhere
    extern std::vector<int> edgeCacheFlat;

    auto edgeCacheIndex = [size](size_t x, size_t y, size_t z, size_t edge) -> size_t {
        size_t dim = size - 1;
        return (((z * dim + y) * dim + x) * 12) + edge;
    };

    for (int i = 0; i < 15 && triangulation[i] != -1; i += 3) {
        size_t triIdx[3];
        for (int v = 0; v < 3; ++v) {
            int edgeIndex = triangulation[i + v];
            if (edgeIndex < 0 || edgeIndex >= 12) continue;

            size_t nx = x, ny = y, nz = z, nedge = static_cast<size_t>(edgeIndex);
            const int* neighbor = EDGE_NEIGHBORS[edgeIndex];
            if (neighbor[0] < 0 && x > 0) {
                nx += neighbor[0];
                nedge = static_cast<size_t>(neighbor[3]);
            } else if (neighbor[1] < 0 && y > 0) {
                ny += neighbor[1];
                nedge = static_cast<size_t>(neighbor[3]);
            } else if (neighbor[2] < 0 && z > 0) {
                nz += neighbor[2];
                nedge = static_cast<size_t>(neighbor[3]);
            }

            // Debug: Log info for z > 213
            // if (z > 213) {
            //     printf("[DEBUG] z=%zu, y=%zu, x=%zu, nx=%zu, ny=%zu, nz=%zu, nedge=%zu\n", z, y, x, nx, ny, nz, nedge);
            // }

            bool inBounds = (nx < size - 1 && ny < size - 1 && nz < size - 1 && nedge < 12);

            size_t flatIdx = 0;
            int cachedIdx = -1;
            if (inBounds) {
                flatIdx = edgeCacheIndex(nx, ny, nz, nedge);
                assert(flatIdx < edgeCacheFlat.size());
                cachedIdx = edgeCacheFlat[flatIdx];
                // if (z > 213) {
                //     printf("[DEBUG] edgeCache access: flatIdx=%zu, cachedIdx=%d\n", flatIdx, cachedIdx);
                // }
                if (cachedIdx != -1) {
                    triIdx[v] = static_cast<size_t>(cachedIdx);
                    continue;
                }
            }

            int v0i = CUBE_EDGES[edgeIndex][0];
            int v1i = CUBE_EDGES[edgeIndex][1];
            int x0 = CUBE_VERTICES[v0i][0], y0 = CUBE_VERTICES[v0i][1], z0 = CUBE_VERTICES[v0i][2];
            int x1 = CUBE_VERTICES[v1i][0], y1 = CUBE_VERTICES[v1i][1], z1 = CUBE_VERTICES[v1i][2];
            size_t idx0 = (x + x0) + (y + y0) * size + (z + z0) * size * size;
            size_t idx1 = (x + x1) + (y + y1) * size + (z + z1) * size * size;
            // if (z > 213) {
            //     printf("[DEBUG] idx0=%zu, idx1=%zu, size=%zu\n", idx0, idx1, size);
            // }
            assert(idx0 < noiseValues.size());
            assert(idx1 < noiseValues.size());

            float v0 = noiseValues[idx0];
            float v1 = noiseValues[idx1];
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
            size_t idx = vertices.size() - 1;
            if (inBounds) {
                edgeCacheFlat[flatIdx] = static_cast<int>(idx);
            }
            triIdx[v] = idx;
        }
        indices.push_back(triIdx[0]);
        indices.push_back(triIdx[1]);
        indices.push_back(triIdx[2]);
    }
}