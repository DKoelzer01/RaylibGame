#ifndef TERRAINGPU_H
#define TERRAINGPU_H

// Planet density blocks on the GPU (Plan 6 step 4, #147): the CPU prepares a block's gnomonic height grid
// (PlanetField::prepareBlock: rivers, lakes and ponds stay on the CPU) and the islands near it; the GPU
// computes every sample (resources/compute/terrain_density.comp), bit-identical to PlanetField::sampleBlock,
// so CPU queries (collision, digging, fluids) still match the drawn surface exactly.
//
// Chunk builds use it asynchronously: a worker prepares the grid and submits; the main thread's pump batches
// the queued blocks into one dispatch and, once its fence has passed (never stalling), hands each block's
// samples to `done` (which queues the build's second half on the workers). Off with RAYLIBGAME_GPU_TERRAIN=0.

#include "voxel.h"
#include <functional>
#include <memory>
#include <vector>

struct GpuBlockJob {
    std::shared_ptr<const DensityField> keep;       // keeps the field alive while queued
    const PlanetField* field = nullptr;
    ChunkBlock block;
    PlanetField::BlockGrid grid;
    std::vector<float> islands;                     // PlanetField::blockIslands
    std::function<bool()> wanted;                   // false: dropped (superseded)
    float priority = 0.0f;                          // the build job's (lower first: edits are negative)
    // Main thread: the batch's samples (shared by its blocks) and this block's first one (E^3 floats).
    std::function<void(std::shared_ptr<const float>, size_t)> done; // samples (a batch, offset in floats): mapped readback memory or a CPU buffer
};

bool TerrainGpu_Available();
void TerrainGpu_Submit(GpuBlockJob&& job);          // any thread
void TerrainGpu_Pump();                             // main thread, once a frame
int TerrainGpu_Pending();                           // queued + in flight
// Synchronous (tests): dispatch, wait, read back. False if unavailable.
bool TerrainGpu_SampleBlock(const PlanetField& pf, Vector3 origin, float stride, int S, float footprint,
                            const PlanetField::BlockGrid& grid, const std::vector<float>& islands, float* out);
struct TerrainGpuStats { long long blocks = 0, batches = 0; double avgLatencyMs = 0.0; };
TerrainGpuStats TerrainGpu_Stats();
void TerrainGpu_Unload();

#endif
