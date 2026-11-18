#ifndef WATERSURFACE_H
#define WATERSURFACE_H

// One water surface (Plan 10 W1, #639): the generated inland water (rivers, lakes, ponds) as a 2D map over a planet,
// read by the chunk water mesher at every level of detail.
//
// Levels (m over sea level) sit on the nodes of a cube-face grid whose spacing matches the samples asking: 0.5 m x 2^k at
// a face's centre. Nodes are computed lazily in tiles of 64 x 64 cells and cached (LRU). A node's value depends only on
// the field and its own position, never on which chunk asked first, so chunks of one level of detail meet exactly and the
// map is deterministic.
//
// The old mesher interpolated levels from a lattice every 4 samples per chunk: 16 m apart or more in distant chunks, so a
// narrow stream fell between nodes (#631) and the slope between a wet node and a dry one crossed the banks as teeth (#635).
// Here the map's spacing is the samples' own.

#include <raylib.h>
#include <cstdint>
#include <memory>

class PlanetField;

namespace WaterSurface {

constexpr int kTile = 64; // cells per tile edge (tiles hold kTile + 1 nodes a side)

struct Tile;

// Reads one field's map at one level of detail; keeps the last tiles it used (no lock on repeat hits). One per mesher
// call (not shared between threads).
class Reader {
public:
    // fp: the sample spacing (m) the map serves; it picks the grid level.
    Reader(const PlanetField& pf, float fp);
    ~Reader();
    // The drawn level at a unit direction (m over sea level): a wet node's level; a dry node beside wet ones takes the
    // lower of their highest level and just under its own ground (the surface stays level out to the shore instead of
    // sloping across it); a dry node with no wet neighbour reads 0 (sea level: under the ground inland). Bilinear between
    // nodes.
    float drawLevel(Vector3 dir);
    float cellMetres() const { return cell; } // the grid's spacing at a face's centre (m)
private:
    const PlanetField& pf;
    int level = 0, n = 1;
    float cell = 0.5f, fpEval = 0.5f;
    struct Slot { uint64_t key = ~0ull; std::shared_ptr<const Tile> tile; };
    Slot slots[4];
    int next = 0;
    const Tile* tileFor(int face, int ti, int tj);
};

void Clear();               // a new world or a load
// Drops the cached tiles (every level of detail) within radius m of dir (a unit direction): the river flow changed there (#656).
void Invalidate(const PlanetField& pf, Vector3 dir, float radius);
long long CacheBytes();     // benches
long long TilesBuilt();     // benches: tiles computed since the start (empty ones too)

} // namespace WaterSurface

#endif
