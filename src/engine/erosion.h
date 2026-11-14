#ifndef EROSION_H
#define EROSION_H

// Erosion bake (Plan 7, #351; step 1 = #355): stream-power erosion with uplift over the whole planet on a cube-sphere
// grid, run once per world at creation, stored beside the save. The input height (today's noise) is the uplift map:
// the land rises where it is high, rivers carve it until uplift and erosion balance (steady state), and the result is
// rescaled to the input's height range. land() adds the change (baked - input) to the noise, so the fine detail stays.
//
// Per step: implicit stream power (Braun & Willett 2013, n = 1): h = (h + dt U + f h_r) / (1 + f), f = K A^m / L,
// receivers first; then a talus limit (slopes steeper than the rock allows slump to it). Every `floodEvery` steps:
// a priority flood from the sea (filled heights, receivers, a processing order) and the drainage area.
// Rock: erodibility varies by region (hard regions keep steep faces), and strata regions alternate hard caprock with
// soft layers in height, which erode into mesas, plateaus and cliff bands. The sea (input below sea level) is the base
// level. Deterministic: the same input gives the same grid, hashed (the world fingerprint for multiplayer, #316).

#include <raylib.h>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

struct ErosionParams {
    float cell = 200.0f;      // m, at a face centre (the gnomonic cells shrink toward the face edges)
    int steps = 400;
    int floodEvery = 10;      // steps between drainage recomputes
    float k = 0.01f;          // stream-power coefficient (per step, area in m^2, m = 0.5)
    float m = 0.5f;
    float uplift = 6.0f;      // m per step at the highest input
    float talusSoft = 0.7f;   // max slope (rise / run) in soft rock (~35 deg)
    float talusHard = 1.7f;   // in hard rock (~60 deg): rock faces
};

class PlanetErosion {
public:
    // height(dir): the input height (m above sea level) at a unit direction, smooth at about `cell`. Threads: the
    // sampling runs on `threads` threads (height must be thread-safe); the simulation is single-threaded.
    void build(const std::function<float(Vector3)>& height, float planetRadius, unsigned int seed, const ErosionParams& p, int threads);
    // The height change at dir (baked - input, m), bilinear. 0 before build.
    float sample(Vector3 dir) const;
    bool ready() const { return n > 0; }
    bool save(const std::string& path) const;          // binary, with the hash
    bool load(const std::string& path, float planetRadius, unsigned int seed); // false: missing, other seed or radius
    uint64_t hash() const { return hashValue; }
    float buildSeconds = 0.0f, sampleSeconds = 0.0f, maxCut = 0.0f, maxRaise = 0.0f;
    int cells() const { return 6 * n * n; }

    // Debug: a face's hillshade (with drainage) before / after (RAYLIBGAME_EROSION_DUMP=1).
    std::vector<float> debugBefore, debugAfter, debugArea;

private:
    int n = 0;                 // cells per face side
    float R = 0.0f;
    unsigned int seed = 0;
    std::vector<float> delta;  // baked - input, per cell (6 n^2)
    uint64_t hashValue = 0;
    void computeHash();
};

#endif
