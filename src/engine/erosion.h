#ifndef EROSION_H
#define EROSION_H

// Erosion bake (#352 prototype, milestone #351): stream-power fluvial erosion with hillslope diffusion on a regional
// grid, run once. The result is used as a height *change* (eroded minus the smooth input) added to the noise
// terrain, so the fine detail stays and the valleys, ridges and drainage come from the simulation.
//
// Grid: gnomonic (tangent plane at `centre`), N x N cells of dx metres. Each step: priority-flood fill (every cell
// drains to the sea or the grid edge), steepest-descent receivers, drainage area, the implicit stream-power update
// (Braun & Willett 2013: h = (h + dt (U + K A^m h_r / L)) / (1 + dt K A^m / L), n = 1), then diffusion. U pulls
// toward the input height (ridges keep their height while valleys incise). Deterministic, single-threaded.

#include <raylib.h>
#include <functional>
#include <vector>

struct ErosionParams {
    float extent = 25000.0f;  // m, square side
    float dx = 100.0f;        // m, cell
    int steps = 120;
    float k = 0.004f;         // stream-power coefficient (per step, area in m^2, m = 0.5)
    float m = 0.5f;
    float relax = 0.01f;      // pull toward the input height per step
    float diffusion = 0.08f;  // hillslope diffusion (fraction of the Laplacian per step, < 0.25)
    float fade = 2000.0f;     // m: the bake fades out toward the grid edge
    bool upliftMode = false;  // the input as an uplift map from low ground (steady-state relief), rescaled to its height range
};

class ErosionBake {
public:
    // height(dir): the input height (m above sea level) at a unit direction, smooth at about dx.
    void build(const std::function<float(Vector3)>& height, Vector3 centre, float planetRadius, const ErosionParams& p);
    // The height change at dir (eroded - input, m) and its weight (1 inside, 0 outside the grid). False outside.
    bool sample(Vector3 dir, float& delta, float& weight) const;
    bool ready() const { return n > 0; }
    float buildSeconds = 0.0f;
    float maxCut = 0.0f, maxFill = 0.0f; // stats

private:
    Vector3 nrm{}, t1{}, t2{};
    float R = 0.0f, half = 0.0f, dx = 1.0f, fade = 1.0f;
    int n = 0;
    std::vector<float> delta; // eroded - input, per cell
};

#endif
