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
    float start = 1.0f;       // the land starts at this fraction of the input: 1 carves into the input relief (0.1 grew
                              // it from low land: heights then follow river length, so narrow coastal ranges came out low)
};

class PlanetErosion {
public:
    // height(dir, footprint): the input height (m above sea level), thread-safe. `threads` for sampling and faces.
    void build(const std::function<float(Vector3, float)>& height, float planetRadius, unsigned int seed, const ErosionParams& p, int threads);
    // The baked height (m above sea level) at dir for a sample footprint (m): fine grid under 2 fine cells, else coarse.
    float heightAt(Vector3 dir, float footprint) const;
    bool ready() const { return n > 0; }
    bool save(const std::string& path) const;          // ERO3: deltas, split bytes, deflated, with the hash
    bool load(const std::string& path, float planetRadius, unsigned int seed); // false: missing, other seed or radius
    uint64_t hash() const { return hashValue; }
    float buildSeconds = 0.0f, coarseSeconds = 0.0f, refineSeconds = 0.0f, maxCut = 0.0f, maxRaise = 0.0f;
    int cells() const { return 6 * n * n; }
    int fineCells() const { return 6 * fe * fe; }
    float fineCell() const { return fe > 0 ? fineMetres : 0.0f; }

private:
    int n = 0;                     // coarse cells per face side
    float R = 0.0f;
    unsigned int seed = 0;
    std::vector<int16_t> coarse;   // baked height, 0.1 m, 6 n^2
    int nf = 0, fm = 0, fe = 0;    // fine cells per face side, margin, extended side (nf + 2 fm)
    float fineMetres = 0.0f;
    std::vector<int16_t> fine;     // baked height, 0.1 m, 6 fe^2 (each face with its margin)
    uint64_t hashValue = 0;
    void computeHash();
    float coarseAt(Vector3 dir) const;
    float fineAt(Vector3 dir) const;
};

#endif
