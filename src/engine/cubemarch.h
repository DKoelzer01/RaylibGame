#ifndef cubemarch_h
#define cubemarch_h

#include "cubemarchconsts.h"

#include <raylib.h>
#include <vector>

// Marching cubes over one chunk's density grid.
//
// `density` is a (size x size x size) grid indexed x + y*size + z*size*size, where size = CHUNK_SIZE + 1
// (the +1 is the shared boundary sample with the neighbouring chunk).
//
// `edgeCache` must hold size*size*size*3 ints initialised to -1. Each grid edge is identified by its
// lower corner and its axis, so a vertex is created once and shared by every cell touching that edge.
// The cache is per chunk on purpose: vertex indices refer to this chunk's `vertices` array only.
// Seams still line up because neighbouring chunks sample identical density values on the shared face.
void marchCube(
    int x, int y, int z,
    const std::vector<float>& density, int size,
    std::vector<Vector3>& vertices,
    std::vector<int>& indices,
    std::vector<int>& edgeCache,
    float threshold
);

#endif
