// --test-terrain-gpu and the visual/terrain-gpu scenario (Plan 6 steps 3-4, #146 #147; shared since #870): the GPU ports of
// the height map and of whole density blocks against the CPU. See terraingpucheck.h.
#include "terraingpucheck.h"
#include "gpucompute.h"
#include "terraingpu.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

TerrainGpuCheck TerrainGpu_CheckAgainstCpu(Scene& world) {
    TerrainGpuCheck R;
    auto fail = [&R](const char* why) { R.error = why; return R; };
            Planetoid* p = dynamic_cast<Planetoid*>(world.objects[0].get());
            const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
            ComputeProgram prog;
            if (!pf || !prog.load("resources/compute/terrain_test.comp")) { printf("bench: terrain gpu: no planet or shader -> FAIL\n"); return fail("no planet or shader"); }
            const int N = 1 << 20;
            std::vector<float> in(N * 4);
            std::mt19937 rng(11);
            std::uniform_real_distribution<float> u(-1.0f, 1.0f);
            const Vector3 home = pf->findSpawnDir(p->seed);
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
            const float fps[] = { 0.0f, 0.0f, 0.25f, 1.0f, 2.0f, 4.0f, 8.0f, 16.0f, 32.0f, 64.0f, 128.0f, 256.0f, 512.0f, 1024.0f };
            for (int k = 0; k < N; ++k) {
                Vector3 d;
                if (k & 1) d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, 5000.0f * u(rng) / pf->radius), Vector3Scale(t2, 5000.0f * u(rng) / pf->radius))));
                else { do d = { u(rng), u(rng), u(rng) }; while (Vector3LengthSqr(d) < 0.01f || Vector3LengthSqr(d) > 1.0f); d = Vector3Normalize(d); }
                in[k * 4] = d.x; in[k * 4 + 1] = d.y; in[k * 4 + 2] = d.z;
                in[k * 4 + 3] = fps[(rng() >> 4) % (sizeof(fps) / sizeof(fps[0]))];
            }
            std::vector<float> params;
            pf->gpuNoiseParams(params);
            ComputeBuffer bin, bout;
            bin.create(N * 16, in.data(), false);
            bout.create(N * 4, nullptr, true);
            prog.use();
            rlSetUniform(prog.uniform("noiseParams"), params.data(), RL_SHADER_UNIFORM_VEC4, (int)params.size() / 4);
            prog.setFloat("planetRadius", pf->radius);
            prog.setInt("planetMoon", 0);
            { const float gp[4] = { pf->gen.ocean, pf->gen.mountains, pf->gen.roughness, 1.0f / std::max(pf->gen.biomeSize, 0.1f) }; rlSetUniform(prog.uniform("genParams"), gp, RL_SHADER_UNIFORM_VEC4, 1); } // #504
            prog.setUint("count", (unsigned)N);
            bin.bind(0); bout.bind(1);
            const double g0 = GetTime();
            prog.dispatch((N + 63) / 64);
            std::vector<float> gpu(N);
            bout.download(gpu.data(), N * 4);
            const double gpuMs = (GetTime() - g0) * 1000.0;
            const double c0 = GetTime();
            std::vector<float> cpu(N);
            // The point kernel ports the noise terrain only (no bake): compare against the CPU without it (#561: with erosion
            // on, eroded CPU heights against the GPU's noise read as a 605 m mismatch).
            PlanetField::tIgnoreBake = true;
            for (int k = 0; k < N; ++k) cpu[k] = pf->elevationBase({ in[k * 4], in[k * 4 + 1], in[k * 4 + 2] }, in[k * 4 + 3]);
            PlanetField::tIgnoreBake = false;
            const double cpuMs = (GetTime() - c0) * 1000.0;
            double worst = 0.0; int over = 0, exact = 0, worstK = 0;
            for (int k = 0; k < N; ++k) {
                const double e = fabs((double)gpu[k] - cpu[k]);
                if (e > worst) { worst = e; worstK = k; }
                over += e > 1e-4; exact += e == 0.0;
            }
            printf("bench: terrain gpu: %d points, %d exact, %d over 1e-4 m, worst %.3g m (cpu %.4f gpu %.4f at fp %.2f); cpu %.0f ms (1 thread), gpu %.1f ms (incl. readback) -> %s\n",
                   N, exact, over, worst, cpu[worstK], gpu[worstK], in[worstK * 4 + 3], cpuMs, gpuMs, over == 0 ? "PASS" : "FAIL");
            bin.destroy(); bout.destroy(); prog.unload();
            R.points = N; R.exact = exact; R.over = over; R.worst = worst; R.cpuMs = cpuMs; R.gpuMs = gpuMs;
            // Step 4 (#147): whole density blocks (as chunk builds sample them: CHUNK_SIZE + 4 per axis), CPU
            // PlanetField::sampleBlock vs the GPU, at every level of detail, around the start and anywhere.
            int blocks = 0, skipped = 0;
            long long samples = 0, bOver = 0, bExact = 0;
            double bWorst = 0.0, cpuBlockMs = 0.0, gpuBlockMs = 0.0;
            const int E = CHUNK_SIZE + 4;
            std::vector<float> cpuB((size_t)E * E * E), gpuB((size_t)E * E * E), isl;
            for (int b = 0; b < 240; ++b) {
                const float stride = 0.25f * (float)(1 << (b % 9));  // 0.25 .. 64 m
                Vector3 d;
                if (b % 3) d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, 3000.0f * u(rng) / pf->radius), Vector3Scale(t2, 3000.0f * u(rng) / pf->radius))));
                else { do d = { u(rng), u(rng), u(rng) }; while (Vector3LengthSqr(d) < 0.01f || Vector3LengthSqr(d) > 1.0f); d = Vector3Normalize(d); }
                const float r = pf->radius + pf->elevation(d) + stride * E * (u(rng) * 0.6f - 0.2f) - 30.0f * (b % 5 == 0);
                const Vector3 o = Vector3AddValue(Vector3Scale(d, r), -0.5f * stride * (E - 1));
                const float fp = std::max(stride, 1.0f);
                PlanetField::BlockGrid g;
                if (!pf->prepareBlock(o, stride, E, fp, g) || !pf->blockIslands(o, stride, E, isl, 64)) { ++skipped; continue; }
                const double c1 = GetTime();
                pf->sampleBlock(o, stride, E, fp, g, cpuB.data());
                cpuBlockMs += (GetTime() - c1) * 1000.0;
                const double g1 = GetTime();
                if (!TerrainGpu_SampleBlock(*pf, o, stride, E, fp, g, isl, gpuB.data())) { printf("bench: terrain gpu: block dispatch failed -> FAIL\n"); return fail("block dispatch failed"); }
                gpuBlockMs += (GetTime() - g1) * 1000.0;
                ++blocks;
                for (size_t k = 0; k < cpuB.size(); ++k) {
                    const double e = fabs((double)gpuB[k] - cpuB[k]);
                    bWorst = std::max(bWorst, e); bOver += e > 1e-4; bExact += e == 0.0; ++samples;
                }
            }
            printf("bench: terrain gpu: %d blocks (%d skipped: off the band), %lld samples, %lld exact, %lld over 1e-4, worst %.3g; per block cpu %.2f ms, gpu %.2f ms (sync, incl. upload + readback) -> %s\n",
                   blocks, skipped, samples, bExact, bOver, bWorst, cpuBlockMs / std::max(1, blocks), gpuBlockMs / std::max(1, blocks), bOver == 0 ? "PASS" : "FAIL");
            R.blocks = blocks; R.skipped = skipped; R.samples = samples; R.bExact = bExact; R.bOver = bOver; R.bWorst = bWorst;
            R.cpuBlockMs = cpuBlockMs / std::max(1, blocks); R.gpuBlockMs = gpuBlockMs / std::max(1, blocks);
    R.ok = R.over == 0 && R.bOver == 0;
    return R;
}
