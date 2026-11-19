#pragma once
// Plan 6 steps 3-4 (#146 #147): the GPU height map (1M points, footprints 0-1024 m) and whole density blocks (240, every
// level of detail) against the CPU's, on body 0. Gate: every value within 1e-4 m. Prints the "bench: terrain gpu:" lines.
// Used by --test-terrain-gpu and the visual/terrain-gpu scenario (#870); needs the GL context (compute shaders).
#include <string>
class Scene;
struct TerrainGpuCheck {
    bool ok = false; std::string error;
    int points = 0, exact = 0, over = 0; double worst = 0.0, cpuMs = 0.0, gpuMs = 0.0;
    int blocks = 0, skipped = 0; long long samples = 0, bExact = 0, bOver = 0; double bWorst = 0.0, cpuBlockMs = 0.0, gpuBlockMs = 0.0;
};
TerrainGpuCheck TerrainGpu_CheckAgainstCpu(Scene& world);
