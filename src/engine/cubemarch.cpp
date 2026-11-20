#include "cubemarch.h"

#include <algorithm>
#include <cassert>

void marchCube(
    int x, int y, int z,
    const std::vector<float>& density, int size,
    std::vector<Vector3>& vertices,
    std::vector<int>& indices,
    std::vector<int>& edgeCache,
    float threshold
) {
    assert(x + 1 < size && y + 1 < size && z + 1 < size);
    auto idx = [size](int ix, int iy, int iz) { return ix + iy * size + iz * size * size; };

    int cubeIndex = 0;
    for (int c = 0; c < 8; ++c) {
        if (density[idx(x + CUBE_VERTICES[c][0], y + CUBE_VERTICES[c][1], z + CUBE_VERTICES[c][2])] < threshold) {
            cubeIndex |= 1 << c;
        }
    }
    if (cubeIndex == 0 || cubeIndex == 255) return;

    const signed char* tri = TRIANGULATIONS[cubeIndex];
    for (int i = 0; i < 15 && tri[i] != -1; i += 3) {
        int triIdx[3];
        for (int v = 0; v < 3; ++v) {
            const int edge = tri[i + v];
            const int* c0 = CUBE_VERTICES[CUBE_EDGES[edge][0]];
            const int* c1 = CUBE_VERTICES[CUBE_EDGES[edge][1]];

            // Identify the edge by its lower grid corner + axis so neighbouring cells share it.
            const int lx = x + std::min(c0[0], c1[0]);
            const int ly = y + std::min(c0[1], c1[1]);
            const int lz = z + std::min(c0[2], c1[2]);
            const int axis = (c0[0] != c1[0]) ? 0 : (c0[1] != c1[1]) ? 1 : 2;
            const int key = idx(lx, ly, lz) * 3 + axis;

            int vi = edgeCache[key];
            if (vi < 0) {
                // Always interpolate from the lower corner to the upper one so the result is
                // bit-identical no matter which cell (or which chunk) computes it.
                const int hx = lx + (axis == 0), hy = ly + (axis == 1), hz = lz + (axis == 2);
                const float da = density[idx(lx, ly, lz)];
                const float db = density[idx(hx, hy, hz)];
                float t = (threshold - da) / (db - da + 1e-8f);
                t = std::max(0.0f, std::min(1.0f, t));
                vertices.push_back({ lx + t * (hx - lx), ly + t * (hy - ly), lz + t * (hz - lz) });
                vi = static_cast<int>(vertices.size()) - 1;
                edgeCache[key] = vi;
            }
            triIdx[v] = vi;
        }
        indices.push_back(triIdx[0]);
        indices.push_back(triIdx[1]);
        indices.push_back(triIdx[2]);
    }
}
