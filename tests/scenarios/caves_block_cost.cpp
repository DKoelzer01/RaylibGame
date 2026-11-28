// caves/block-cost (#2070 step 0): what caves cost a chunk build today, by depth, before caves reach the mantle. Blocks
// as chunk builds sample them (CHUNK_SIZE + 4 per axis) are centred 0-600 m under the ground within 3 km of seed 42's
// start, at two strides. Per block: its cave segments (blockCaveSegs, no cap: a histogram against the GPU's
// kGpuCaveSegsMax, over which the block takes the CPU path) and the density's cost on the CPU (sampleBlock, one thread)
// and the GPU (synchronous, incl. upload and readback), split into blocks with and without cave segments.
// The budget #2070 sets: deep blocks no slower than today's shallow blocks with caves, and few over the cap.
// Reported only (timing is meaningful on a quiet -O3 build: tools/slot.sh quiet, -DRAYLIBGAME_FAST_BUILD=OFF).
#include "scenario.h"
#include "terraingpu.h"
#include "voxel.h"
#include "world.h"
#include "raymath.h"
#include <random>
#include <string>
#include <vector>

namespace {
constexpr float kDepths[] = { 0.0f, 30.0f, 60.0f, 100.0f, 150.0f, 250.0f, 400.0f, 600.0f }; // m under the ground (block centre)
constexpr float kStrides[] = { 0.25f, 1.0f };     // m: the finest chunk level and a mid one
constexpr int kBlocksPer = 24;                    // blocks per depth and stride
constexpr int kSegBins[] = { 0, 1, 9, 17, 33, 65 }; // histogram lower edges: 0, 1-8, 9-16, 17-32, 33-64, 65+
constexpr int kFloatsPerSeg = 24;                 // 6 vec4 a segment (voxel.cpp GatherCaveSegs; was 20, counting 1.2x high)
}

SCENARIO("caves/block-cost") {
    s.tags = { "caves", "perf" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 300.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("blocks by depth", 200.0, [](ScenarioCtx& c) {
        if (!c.field || !c.body) { c.fail("no planet"); return true; }
        const PlanetField& pf = *c.field;
        const bool gpuOn = TerrainGpu_Available();
        const Vector3 home = pf.findSpawnDir(c.body->seed);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
        const int E = CHUNK_SIZE + 4;
        std::vector<float> out((size_t)E * E * E), isl, cav;
        int measured = 0;
        for (float stride : kStrides) {
            std::mt19937 rng(2070); // the same places at every depth and stride
            std::uniform_real_distribution<float> u(-1.0f, 1.0f);
            for (float depth : kDepths) {
                int n[2] = {}, hist[6] = {}, over = 0, skipped = 0;
                double cpuMs[2] = {}, gpuMs[2] = {}; int gpuN[2] = {};
                long long segSum = 0; int segMax = 0;
                for (int b = 0; b < kBlocksPer; ++b) {
                    const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, 3000.0f * u(rng) / pf.radius), Vector3Scale(t2, 3000.0f * u(rng) / pf.radius))));
                    const float r = pf.radius + pf.elevation(d) - depth;
                    const Vector3 o = Vector3AddValue(Vector3Scale(d, r), -0.5f * stride * (E - 1));
                    const float fp = std::max(stride, 1.0f);
                    PlanetField::BlockGrid g;
                    if (!pf.prepareBlock(o, stride, E, fp, g) || !pf.blockIslands(o, stride, E, isl, 64)) { ++skipped; continue; }
                    pf.blockCaveSegs(o, stride, E, fp, cav, INT32_MAX);
                    const int segs = (int)cav.size() / kFloatsPerSeg;
                    int bin = 0;
                    while (bin + 1 < 6 && segs >= kSegBins[bin + 1]) ++bin;
                    ++hist[bin]; segSum += segs; segMax = std::max(segMax, segs);
                    over += segs > PlanetField::kGpuCaveSegsMax;
                    const int k = segs > 0; // 0: no cave segments, 1: with
                    const double c0 = GetTime();
                    pf.sampleBlock(o, stride, E, fp, g, out.data());
                    cpuMs[k] += (GetTime() - c0) * 1000.0;
                    ++n[k];
                    if (gpuOn && segs <= PlanetField::kGpuCaveSegsMax) {
                        const double g0 = GetTime();
                        if (TerrainGpu_SampleBlock(pf, o, stride, E, fp, g, isl, out.data(), &cav)) { gpuMs[k] += (GetTime() - g0) * 1000.0; ++gpuN[k]; }
                    }
                }
                measured += n[0] + n[1];
                const std::string at = TextFormat("%.0f m s%.2f", depth, stride);
                auto mean = [](double sum, int cnt) { return cnt > 0 ? sum / cnt : 0.0; };
                c.metric(at + " cpu ms (no caves)", mean(cpuMs[0], n[0]), "ms");
                c.metric(at + " cpu ms (caves)", mean(cpuMs[1], n[1]), "ms");
                if (gpuOn) {
                    c.metric(at + " gpu ms (no caves)", mean(gpuMs[0], gpuN[0]), "ms");
                    c.metric(at + " gpu ms (caves)", mean(gpuMs[1], gpuN[1]), "ms");
                }
                c.metric(at + " blocks with caves", n[1], "");
                c.metric(at + " blocks over the cap", over, "");
                c.note(TextFormat("%s: %d blocks (%d skipped), %d with caves; segments mean %.1f, max %d; histogram 0 / 1-8 / 9-16 / 17-32 / 33-64 / 65+: %d %d %d %d %d %d; %d over the cap (%d)",
                                  at.c_str(), n[0] + n[1], skipped, n[1], mean((double)segSum, n[0] + n[1]), segMax,
                                  hist[0], hist[1], hist[2], hist[3], hist[4], hist[5], over, PlanetField::kGpuCaveSegsMax));
            }
        }
        c.check("blocks measured (>= half)", measured * 2 >= kBlocksPer * (int)(sizeof(kDepths) / sizeof(kDepths[0])) * (int)(sizeof(kStrides) / sizeof(kStrides[0])),
                std::to_string(measured));
        c.note(gpuOn ? "GPU density: on" : "GPU density: unavailable (CPU times only)");
        return true;
    });
}
