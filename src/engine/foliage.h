#ifndef FOLIAGE_H
#define FOLIAGE_H

#include <vector>
#include <raylib.h>

class DensityField;
struct ChunkMeshData;

// Foliage (design doc section 11). Step 1: grass tufts.
//
// Placement is a pure function of the chunk's mesh: after a finest-level chunk is meshed (on the
// worker), tufts are scattered over its triangles from a hash of (node, triangle, k) and kept only
// where grass is the dominant material, the ground is open to the sky, the slope is gentle and the
// point is above sea level. Same seed, same grass; nothing is stored. Digging re-meshes the chunk,
// and the new (dug) surface is dirt, so the grass there is gone.
//
// Each tuft instance is 8 floats: position (chunk-local, mesh frame), scale (m), body-local up
// (xyz) and a hash in [0, 1) for rotation, lean and colour variation.
constexpr int kGrassFloats = 8;
constexpr float kGrassPerSquareMetre = 4.0f;
constexpr float kGrassDrawDistance = 70.0f; // fully faded out by here (m)
extern bool foliageEnabled;  // draw foliage (bench: --no-foliage)

void ScatterGrass(const DensityField& field, ChunkMeshData& out); // worker thread

// GPU side, per chunk: a VAO pairing the shared tuft mesh with this chunk's instance buffer.
struct FoliageGpu {
    unsigned int vao = 0, vbo = 0;
    int count = 0;
};
void Foliage_Upload(FoliageGpu& g, const std::vector<float>& instances); // main thread
void Foliage_Release(FoliageGpu& g);
// Draws one chunk's tufts; the grass shader must be bound and its view-projection set.
void Foliage_Draw(const FoliageGpu& g, Shader shader, int modelLoc, const Matrix& model);

#endif
