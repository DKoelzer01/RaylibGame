#ifndef EROSION_H
#define EROSION_H

// Erosion bake (Plan 7, #351): stream-power erosion with uplift over the whole planet, run once per world at creation
// and stored beside the save.
//
// Step 1 (#355), the planet at ~200 m on a cube-sphere grid: the input height (the noise land) is the uplift map; the
// land rises where it is high, rivers carve it until uplift and erosion balance, and the result is rescaled to the
// input's height range. Per step: implicit stream power (Braun & Willett 2013, n = 1), receivers first; a talus limit
// (slopes steeper than the rock holds slump). Every `floodEvery` steps: a priority flood from the sea (filled heights,
// receivers, an order) and the drainage area. Rock: hard regions hold steeper faces; strata regions alternate hard
// caprock and soft layers (mesas, plateaus, cliff bands).
//
// Step 2 (#356), each cube face again at 4x the resolution (~50 m), from the coarse result, with a margin past the
// face edges so rivers cross them; sampling blends the faces over the margins. land() reads the height: the fine grid
// for fine footprints, the coarse one for coarse footprints (no aliasing far away).
//
// Deterministic: the same input gives the same grids, hashed (the multiplayer world fingerprint, #316).

#include <raylib.h>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>
#include "strata.h"

using StrataFn = std::function<StrataColumn(Vector3)>; // #587: the rock column at a direction (PlanetField::strataColumn)

// Bump with any change to the bake's code (Simulate, SetRock, the refinement, the pond extraction): a cache from
// another version or other parameters is baked again rather than loaded (#355 review).
constexpr uint32_t kErosionBakeVersion = 27; // 27: hillslope diffusion, random receivers (#1884); 26: varied talus limits, lateral erosion on the fine grid (#633); 25: final drainage, diagonal channel corners worn, channels in the bake (#547); 22: graded canyon floors (#544), divides past 400 m of cut; 19: formations kept inside the elevation bound; 18: // 18: glacial 2.0 (#548); 17: // 16: volcanoes (#593); 15: karst (#592); // 13: fold outcrops raised in the bake input; 12: folded strata (#591); 11: narrower fault steps; 10: fault scarps in the bake input (#590); 9: canyons in the bake input (#544); 8: rock strata (#587); 7: the noise seed offsets changed (#532); 6: wind (aeolian) erosion, #370
struct ErosionParams {
    float cell = 200.0f;      // m, coarse cell at a face centre (gnomonic cells shrink toward the face edges)
    int steps = 400;
    int floodEvery = 10;      // steps between drainage recomputes
    float k = 0.01f;          // stream-power coefficient (per step, area in m^2, m = 0.5)
    float m = 0.5f;
    float uplift = 2.0f;      // m per step at the highest input (6 with start 0.1 flattened narrow coastal ranges)
    float talusSoft = 0.7f;   // max slope (rise / run) in soft rock (~35 deg)
    float talusHard = 1.7f;   // in hard rock (~60 deg): rock faces
    int refine = 4;           // fine cells per coarse cell (0: no refinement)
    int refineSteps = 100;
    int margin = 64;          // fine cells past each face edge
    float glacial = 2.0f;     // (#548: 1.5 -> 2, slightly stronger carving) glacial erosion rate (x the fluvial rate at 1 km^2 of drainage); 0: none (#358)
    float glacialLine = 900.0f; // m: the ice-age snowline at the equator (falls to sea level by latitude ~0.8)
    float wind = 1.0f;        // aeolian erosion strength on dry land (#370; 0: none, e.g. airless bodies). Scales the
                              // per-step deflation on windward faces and the sand carried downwind.
    uint64_t landKey = 0;     // the input land's generation parameters (WorldGenParams::key, #504; 0: the defaults)
    float start = 1.0f;       // the land starts at this fraction of the input: 1 carves into the input relief (0.1 grew
                              // it from low land: heights then follow river length, so narrow coastal ranges came out low)
    float hillslope = 40.0f;  // #1884: hillslope diffusion (soil creep), m^2 per step: rounds crests and the shoulders of
                              // V-notches; small next to a coarse cell, about a cell over the fine pass (GUESS; 0: off)
    bool randomReceivers = true; // #1884: each flood gives a cell a receiver among its lower neighbours, odds by slope,
                              // instead of the steepest one: channels stop following the 8 grid directions (false: steepest)
};
uint64_t ErosionParamsKey(const ErosionParams& p); // FNV-1a over every field and kErosionBakeVersion
bool BakeRiversOld(); // #547 A/B: RAYLIBGAME_RIVERS_OLD=1 = the 200 m river network and its carving (no final drainage in the bake)
float NoiseSeedSpanUnits(); // voxel.cpp (#532): the noise seed offset span, in the cache key when not the default

class PlanetErosion {
public:
    // height(dir, footprint): the input height (m above sea level), thread-safe. `threads` for sampling and faces.
    // strata (#587): the rock column per cell (layer hardness and erodibility); empty: the old hardness noise and bands.
    void build(const std::function<float(Vector3, float)>& height, float planetRadius, unsigned int seed, const ErosionParams& p, int threads, const StrataFn& strata = StrataFn());
    // The baked height (m above sea level) at dir for a sample footprint (m): fine grid under 2 fine cells, else coarse.
    float heightAt(Vector3 dir, float footprint) const;
    // The height and its uphill gradient (m per m, in the tangent plane at dir): central differences one cell apart, so
    // the gradient is continuous (the gullies' direction must not jump at cell borders).
    struct Slope { float h; Vector3 grad; };
    Slope slopeAt(Vector3 dir, float footprint) const;
    bool ready() const { return n > 0; }
    bool save(const std::string& path) const;          // ERO6: the params key, deltas, split bytes, deflated in 16 pieces; the ponds; the hash; written to path.tmp, then renamed
    bool load(const std::string& path, float planetRadius, unsigned int seed, uint64_t paramsKey); // false: missing, other seed, radius, parameters or code version
    uint64_t hash() const { return hashValue; }
    float buildSeconds = 0.0f, coarseSeconds = 0.0f, refineSeconds = 0.0f, maxCut = 0.0f, maxRaise = 0.0f;
    std::string profile;      // #605: seconds per pass of the last build()
    float heightScale = 1.0f; // #587: the bake's result = its simulation heights x this (strata are laid in simulation heights)
    int cells() const { return 6 * n * n; }
    int fineCells() const { return 6 * fe * fe; }
    float fineCell() const { return fe > 0 ? fineMetres : 0.0f; }
    // Ponds (#360): the small closed depressions the 50 m refinement left (at most 12 cells): centre, a radius, the
    // spill level (m above sea level) and depth. Water levels come from these; PlanetField places its ponds here.
    struct Pond { Vector3 dir; float r, level, depth; };
    const std::vector<Pond>& ponds() const { return pondList; }

    // Drainage on the final terrain (#547, bake-native rivers): after the faces, one priority flood over every face's own
    // 50 m cells (a single grid across the face edges) gives each cell its receiver (steepest descent on the filled
    // surface), drainage area and filled height. Channels are worn into the stored heights right there (a light incision
    // growing with the drainage), so the terrain and the rivers agree by construction. Cells: g = (face * nf + j) * nf + i.
    static constexpr float kFlowArea0 = 1.0e5f; // m^2: area class 1 starts here; each class doubles
    bool hasFlow() const { return !flow.empty(); }
    int flowCells() const { return 6 * nf * nf; }
    int flowCellOf(Vector3 dir) const;
    Vector3 flowDir(int g) const;
    void flowNeighbours(int g, int out[8]) const;     // -1 where none (8 directions: E W N S NE SE NW SW in i, j)
    float flowHeight(int g) const;                    // the stored (incised) height, m
    int flowReceiver(int g) const;                    // the cell it drains to, or -1 (the sea, an outlet)
    float flowArea(int g) const;                      // m^2 drained through it (its class's middle), 0 under kFlowArea0
    static float IncisionDepth(float area);           // m the bake wore a channel cell down by
    // Filled cells (lakes and ponds of the flood): cell -> the water level there (m).
    const std::vector<std::pair<int, float>>& filledCells() const { return filled; }

private:
    std::vector<uint8_t> flow;                       // per cell: receiver direction (low 4 bits, 15: none) | area class << 4
    std::vector<std::pair<int, float>> filled;       // sorted by cell
    void routeFinal(int threads);                    // the drainage and incision above (build's last pass)
public:
    float routeSeconds = 0.0f;

private:
    int n = 0;                     // coarse cells per face side
    float R = 0.0f;
    unsigned int seed = 0;
    std::vector<int16_t> coarse;   // baked height, 0.1 m, 6 n^2
    int nf = 0, fm = 0, fe = 0;    // fine cells per face side, margin, extended side (nf + 2 fm)
    float fineMetres = 0.0f;
    std::vector<int16_t> fine;     // baked height, 0.1 m, 6 fe^2 (each face with its margin)
    std::vector<Pond> pondList;
    uint64_t hashValue = 0;
    uint64_t paramsKey = 0;  // ErosionParamsKey of the build (saved with the cache)
    void computeHash();
    float coarseAt(Vector3 dir) const;
    float fineAt(Vector3 dir) const;
};

// Gullies (Plan 7 step 3, #357): an erosion-filter noise. Each octave (900, 420, 190, 85, 38, 17 m) is a sum of cosine stripes
// running down the slope, one kernel per corner of a 3D lattice (no cube-face seams), so the ribs and gullies follow
// the baked drainage; each octave's slope bends the next one's direction (small gullies branch into larger ones).
// Strength follows the slope (flat ground stays flat); octaves finer than the sample footprint fade out.
// p = surface point (m), dir = its unit direction, grad = slopeAt's gradient. Returns a height offset (m).
float ErosionGullies(Vector3 p, Vector3 dir, Vector3 grad, float footprint, unsigned int seed, float strength);
float ErosionBakeProgress(); // 0..1 while a bake runs (the loading screen, #362)
bool ErosionBaking();        // a bake is running (not a cache read: #505)
int ErosionBakeCount();      // bakes run so far (tests: a load must read the cache, #505)

#endif
