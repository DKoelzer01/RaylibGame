#ifndef TERRAINARENA_H
#define TERRAINARENA_H

// Terrain mesh arena (Plan 6 steps 1-2, #144 / #145): every chunk's terrain mesh lives in shared buffers
// (ranges allocated per chunk, freed on rebuild), and a body's chunks draw with ONE
// glMultiDrawElementsIndirect per pass.
//
// Positions sit in their own buffer (the shadow passes read only them), the other attributes interleaved.
// Each chunk also has a slot in a GPU-resident chunk table (origin, size, overlap clip, its index range).
// Per pass the CPU only lists the slots of the LOD-selected chunks; a compute pass (chunk_cull.comp) places
// them (double precision), culls them against the view frustum or the shadow light's column and writes the
// indirect commands (culled: instanceCount 0, so the near-first order holds) and per-draw records (model rows,
// texturing origin, clip). The records are instanced attributes (locations 9-13) picked by each command's
// baseInstance; lighting.vs / depth.vs take this path when the uniform `arenaDraw` is 1, and the fragment
// shader's matModel is the body's rotation (shared by its chunks).
//
// The game requires GL 4.3 (#185). TerrainArena_Available() is false only with RAYLIBGAME_ARENA=0 or if the cull shader
// fails to compile (e.g. no fp64): then chunks keep their own raylib meshes and per-chunk draws.

#include <raylib.h>
#include <cstdint>
#include <vector>

struct ChunkMeshData;
struct Frustum;

struct ArenaAlloc {
    int vFirst = -1, vCount = 0;   // vertices
    int iFirst = -1, iCount = 0;   // indices (16-bit, relative to vFirst: the command's baseVertex)
    int slot = -1;                 // chunk table entry
    Vector3 clip = { -1.0f, -1.0f, -1.0f }; // the clip last written to the table (chunk-local)
    bool valid() const { return slot >= 0 && iCount > 0; }
};

// How a body places its chunks (ChunkPlacer): model = [rot * scale | offset + rot * origin * scale].
struct ArenaBody {
    double rot[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };  // row-major
    double scale = 1.0;
    double offset[3] = { 0, 0, 0 };                 // body position - render origin
    Matrix rs{};                                    // rot * scale as a raylib matrix (the fragment shader's matModel)
    double texPeriod = 20.0;                        // texture repeat (m)
};
// What a pass culls against: the view frustum (render space; null = none) and, for shadow passes, the light
// box's column in body-local space (radius <= 0: none).
struct ArenaCull {
    const Frustum* frustum = nullptr;
    Vector3 colPoint{}, colDir{};
    float colRadius = -1.0f;
};

bool TerrainArena_Available();
ArenaAlloc TerrainArena_Upload(const ChunkMeshData& data);  // main thread; fills the table entry from data.key
void TerrainArena_Free(ArenaAlloc& a);

// Drawing: queue the LOD-selected chunks of one body (with their current chunk-local clip), then flush with
// the pass's shader already bound (ChunkBatch did that: textures, palette, occluder).
void TerrainArena_Queue(ArenaAlloc& a, Vector3 clip);
// Culls and draws what was queued (one compute dispatch + one MDI), then clears the queue. Returns the
// number queued (drawn or culled).
int TerrainArena_Flush(unsigned int shaderId, bool positionsOnly, const ArenaBody& body, const ArenaCull& cull);

struct ArenaStats { int vertices = 0, vertexCapacity = 0, indices = 0, indexCapacity = 0, draws = 0, flushes = 0; long long staged = 0, direct = 0; }; // uploads through the staging ring / glBufferSubData (#286)
ArenaStats TerrainArena_Stats();
void TerrainArena_Unload();

#endif
