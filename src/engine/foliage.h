#ifndef FOLIAGE_H
#define FOLIAGE_H

#include <vector>
#include <memory>
#include <unordered_set>
#include "dvec3.h"
#include "voxel.h"
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
// Each tuft instance is 8 floats: position (chunk-local, mesh frame), scale (m) + 10 * variant, body-local up
// (xyz) and a hash in [0, 1) for rotation, lean and colour variation.
constexpr int kGrassFloats = 8;
constexpr float kGrassPerSquareMetre = 7.0f;
constexpr int kGrassVariants = 3;            // sparse tuft, clump, dense clump (foliage.cpp)
inline float kGrassDrawDistance = 100.0f;    // fully faded out by here (m); --grass-dist
extern bool foliageEnabled;  // draw foliage (bench: --no-foliage)

void ScatterGrass(const DensityField& field, ChunkMeshData& out); // worker thread

// #830: no grass in the cells of placed parts (a bed, a floor, a machine). A body's snapshot of its plots (frame and occupied
// cells, plot space), swapped whole on the main thread when parts change and read by ScatterGrass on the workers.
struct GrassBlockPlot { DVec3 origin{}; Vector3 east{}, up{}, north{}; float radius = 0.0f; std::unordered_set<Int3> cells;
                        std::unordered_set<Int3> solid; }; // (#908) solid: building blocks (not glass, not small parts): they shade the ground under them
using GrassBlockers = std::vector<GrassBlockPlot>;
void Grass_SetBlockers(const DensityField* field, std::shared_ptr<const GrassBlockers> blockers);
std::shared_ptr<const GrassBlockers> Grass_BlockersOf(const DensityField* field); // (#908) the chunk builder's sky probe reads it too

// GPU side, per chunk: VAOs pairing the shared tuft meshes with this chunk's instance buffer.
struct FoliageGpu {
    unsigned int vao[kGrassVariants] = {}, vbo = 0;  // a VAO per variant over one instance buffer
    int count[kGrassVariants] = {};
};
void Foliage_Upload(FoliageGpu& g, const std::vector<float>& instances); // main thread
void Foliage_Release(FoliageGpu& g);
// Draws one chunk's tufts; the grass shader must be bound and its view-projection set.
void Foliage_Draw(const FoliageGpu& g, Shader shader, int modelLoc, const Matrix& model);

#endif
