// visual/terrain-gpu (#870; moved from --test-terrain-gpu, Plan 6 steps 3-4, #146 #147): the GPU's height map (1M points,
// footprints 0-1024 m) and whole density blocks (240, every level of detail) against the CPU's on the planet at seed 42.
// PASS as the flag: no value off by more than 1e-4 m (the GPU terrain is bit-identical by design).
#include "scenario.h"
#include "terraingpucheck.h"
#include "world.h"

SCENARIO("visual/terrain-gpu") {
    s.tags = { "visual", "terrain", "gpu", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 180.0;
    s.legacyFlag = "--test-terrain-gpu";
    s.record = false; // (the check takes seconds of CPU on one frame: nothing to watch)
    s.step("compare GPU and CPU terrain", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const TerrainGpuCheck r = TerrainGpu_CheckAgainstCpu(c.world);
        if (!r.error.empty()) { c.fail("terrain gpu: " + r.error); return true; }
        c.gt("points compared", r.points, 0.5);
        c.eq("points over 1e-4 m", r.over, 0.0);
        c.metric("points exact", r.exact);
        c.metric("worst point error", r.worst, "m");
        c.gt("blocks compared", r.blocks, 0.5);
        c.eq("block samples over 1e-4 m", (double)r.bOver, 0.0);
        c.metric("worst block error", r.bWorst, "m");
        c.metric("block cpu", r.cpuBlockMs, "ms");
        c.metric("block gpu", r.gpuBlockMs, "ms");
        return true;
    });
}
