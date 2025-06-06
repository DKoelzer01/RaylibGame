#include "cubemarch.h"

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
    std::vector<int>& edgeCacheFlat,
    float threshold,
    int cx, int cy, int cz, int chunkGrid,
    int cacheType // <-- new parameter
) {
    // Bounds check for cube
    assert(x + 1 < size && y + 1 < size && z + 1 < size);

    std::array<int, 15> triangulation = get_triangulation(x, y, z, noiseValues, size, threshold);

    // Flat array cache: size*size*size*12, initialized to -1 elsewhere
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

            // Only skip the cache for the outermost planetoid boundary
            bool isOutermostX = (cx == chunkGrid - 1) && (nx == size - 1);
            bool isOutermostY = (cy == chunkGrid - 1) && (ny == size - 1);
            bool isOutermostZ = (cz == chunkGrid - 1) && (nz == size - 1);
            bool inBounds = (nx < size - 1 && ny < size - 1 && nz < size - 1 && nedge < 12) &&
                !(isOutermostX || isOutermostY || isOutermostZ);

            size_t flatIdx = 0;
            if (cacheType == 1) { // X face
                if (ny >= size - 1 || nz >= size - 1) {
                    logger.logf("Skipping out-of-bounds X face: nx=%zu ny=%zu nz=%zu\n", nx, ny, nz);
                    continue; // skip out-of-bounds
                }
                flatIdx = (ny + nz * (size - 1)) * 12 + nedge;
            } else if (cacheType == 2) { // Y face
                if (nx >= size - 1 || nz >= size - 1) {
                    logger.logf("Skipping out-of-bounds Y face: nx=%zu ny=%zu nz=%zu\n", nx, ny, nz);
                    continue; // skip out-of-bounds
                }
                flatIdx = (nx + nz * (size - 1)) * 12 + nedge;
            } else if (cacheType == 3) { // Z face
                if (nx >= size - 1 || ny >= size - 1) {
                    logger.logf("Skipping out-of-bounds Z face: nx=%zu ny=%zu nz=%zu\n", nx, ny, nz);
                    continue; // skip out-of-bounds
                }
                flatIdx = (nx + ny * (size - 1)) * 12 + nedge;
            } else { // local
                flatIdx = edgeCacheIndex(nx, ny, nz, nedge);
            }
            int cachedIdx = -1;
            if (inBounds) {
                assert(flatIdx < edgeCacheFlat.size());
                cachedIdx = edgeCacheFlat[flatIdx];
            }

            // if (cacheType > 0 && inBounds && cachedIdx != -1) {
            //     const Vector3& v = vertices[cachedIdx];
            //     logger.logf("READ: Chunk (%d,%d,%d) cacheType=%d flatIdx=%zu -> vertexIdx=%d pos=(%.6f,%.6f,%.6f)\n",
            //         cx, cy, cz, cacheType, flatIdx, cachedIdx, v.x, v.y, v.z);
            // }


            if (cachedIdx != -1) {
                triIdx[v] = static_cast<size_t>(cachedIdx);
            } else {
                int v0i = CUBE_EDGES[edgeIndex][0];
                int v1i = CUBE_EDGES[edgeIndex][1];
                int x0 = CUBE_VERTICES[v0i][0], y0 = CUBE_VERTICES[v0i][1], z0 = CUBE_VERTICES[v0i][2];
                int x1 = CUBE_VERTICES[v1i][0], y1 = CUBE_VERTICES[v1i][1], z1 = CUBE_VERTICES[v1i][2];
                size_t idx0 = (x + x0) + (y + y0) * size + (z + z0) * size * size;
                size_t idx1 = (x + x1) + (y + y1) * size + (z + z1) * size * size;
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
                
                // if (cacheType > 0 && inBounds) {
                //     const Vector3& v = vertices[idx];
                //     logger.logf("WRITE: Chunk (%d,%d,%d) cacheType=%d flatIdx=%zu -> vertexIdx=%zu pos=(%.6f,%.6f,%.6f)\n",
                //         cx, cy, cz, cacheType, flatIdx, idx, v.x, v.y, v.z);
                // }

                triIdx[v] = idx;
            }

            for (int v = 0; v < 3; ++v) {
                if (triIdx[v] == -1) {
                    logger.logf("WARNING: Using -1 vertex index at chunk (%d,%d,%d), cube (%d,%d,%d), edge %d, cacheType %d\n",
                        cx, cy, cz, x, y, z, nedge, cacheType);
                }
            }
        }
        indices.push_back(triIdx[0]);
        indices.push_back(triIdx[1]);
        indices.push_back(triIdx[2]);

        if (cacheType > 0) {
            logger.logf("BORDER TRI: chunk (%d,%d,%d) face %d cube (%zu,%zu,%zu) indices [%zu,%zu,%zu]\n",
                cx, cy, cz, cacheType, x, y, z, triIdx[0], triIdx[1], triIdx[2]);
            logger.logf("  v0: (%.6f, %.6f, %.6f)\n", vertices[triIdx[0]].x, vertices[triIdx[0]].y, vertices[triIdx[0]].z);
            logger.logf("  v1: (%.6f, %.6f, %.6f)\n", vertices[triIdx[1]].x, vertices[triIdx[1]].y, vertices[triIdx[1]].z);
            logger.logf("  v2: (%.6f, %.6f, %.6f)\n", vertices[triIdx[2]].x, vertices[triIdx[2]].y, vertices[triIdx[2]].z);
        }
    }

    
    if ((cacheType == 1 && x == 0) ||
        (cacheType == 2 && y == 0) ||
        (cacheType == 3 && z == 0)) {
        logger.logf("ERROR: Using shared cache on 0 face! x=%zu y=%zu z=%zu cacheType=%d\n", x, y, z, cacheType);
    }
}