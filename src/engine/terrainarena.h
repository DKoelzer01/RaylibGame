#ifndef TERRAINARENA_H
#define TERRAINARENA_H

// Terrain mesh arena (Plan 6 step 1, #144): every chunk's terrain mesh lives in one shared vertex buffer and
// one shared index buffer (ranges allocated per chunk, freed on rebuild), so a body's visible chunks draw
// with ONE glMultiDrawElementsIndirect per pass instead of a draw call (plus two uniform uploads and a
// matrix upload) per chunk.
//
// Positions sit in their own buffer (the shadow passes read only them), the other attributes interleaved.
// Per-draw data (model rows, texturing origin, overlap clip) are instanced attributes (locations 9-13) that
// each command's baseInstance selects (the GL 4.3 base-instance trick: one path on 4.3 and 4.6, and the
// data is fetched once per draw rather than per vertex). lighting.vs / depth.vs take this path when the
// uniform `arenaDraw` is 1; the fragment shader's matModel is the body's rotation (shared by its chunks).
//
// Needs GL 4.3 (multi-draw indirect). Without it TerrainArena_Available() is false and chunks keep their
// own raylib meshes and per-chunk draws (the old path, also used with RAYLIBGAME_ARENA=0).

#include <raylib.h>
#include <cstdint>
#include <vector>

struct ChunkMeshData;

struct ArenaAlloc {
    int vFirst = -1, vCount = 0;   // vertices
    int iFirst = -1, iCount = 0;   // indices (16-bit, relative to vFirst: the command's baseVertex)
    bool valid() const { return vFirst >= 0 && iCount > 0; }
};

bool TerrainArena_Available();
ArenaAlloc TerrainArena_Upload(const ChunkMeshData& data);  // main thread
void TerrainArena_Free(ArenaAlloc& a);

// Drawing: queue the visible chunks of one body, then flush with the pass's shader already bound
// (ChunkBatch did that: textures, palette, occluder). Per-draw records: model (render space), the chunk's
// wrapped texturing origin and its chunk-local clip bound.
void TerrainArena_Add(const ArenaAlloc& a, const Matrix& model, Vector3 texOrigin, Vector3 clip);
void TerrainArena_Flush(unsigned int shaderId, bool positionsOnly = false); // one MDI of what was queued, then clears (positionsOnly: shadow passes)

struct ArenaStats { int vertices = 0, vertexCapacity = 0, indices = 0, indexCapacity = 0, draws = 0, flushes = 0; };
ArenaStats TerrainArena_Stats();
void TerrainArena_Unload();

#endif
