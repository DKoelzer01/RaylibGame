// Shared by worldgen/metrics-1 (1 seed, the PR path) and worldgen/metrics-8 (8 seeds, nightly) (#1882): the measuring stick for
// terrain generation changes (#1883 spawn scoring, #1884 erosion bake, #351). Report-only, apart from one check: every land
// biome must appear (on the one seed, or across the 8).
//
// Each seed's planet is built in this process (PlanetField, default generation parameters, the whole-planet erosion bake
// with no cache: bakeErosion("")) on a worker thread while the frames keep running; no world is loaded for it. The first
// seed is the scenario world's own planet (body 0's terrain seed: seed 42's Terra), the rest are fixed terrain seeds 1..7.
// Over a Fibonacci sphere of directions (one per ~0.75 km^2) it measures:
//   land %             elevation above sea level (0 m); lakes count as land
//   biome shares       biomeAt over the land samples, % of land
//   slope histogram    over land, from the heights 8 m east, west, north and south (a 16 m baseline: a building plot's size)
//   flat buildable %   land whose four 8 m neighbours lie within kFlatTolM of the centre (GUESS: a plot levels with <= 1 m of cut)
//   bake hash          PlanetErosion::hash (low 32 bits as a metric, the full hash as a note): stable unless the bake changes
//   spawn              findSpawnDir's elevation and biome, and its spawnScore (#1883): flat share, nearest fresh water,
//                      ore kinds within 300 m; and how long the search took
#pragma once
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

namespace worldgen_metrics {

constexpr int kSamples = 60000;          // directions over the sphere (60 km radius: ~0.75 km^2 each)
constexpr float kStepM = 8.0f;           // the neighbours' distance (m)
constexpr float kFlatTolM = 1.0f;        // GUESS: flat buildable when the neighbours are within this of the centre (m)
constexpr int kSlopeBins = 6;
constexpr float kSlopeEdges[kSlopeBins - 1] = { 5.0f, 10.0f, 20.0f, 30.0f, 45.0f }; // degrees
inline const char* SlopeBinName(int i) {
    static const char* n[kSlopeBins] = { "0-5", "5-10", "10-20", "20-30", "30-45", "45+" };
    return n[i];
}

struct SeedStats {
    unsigned seed = 0;
    uint64_t hash = 0;
    double bakeSec = 0, sampleSec = 0;
    int samples = 0, land = 0, flat = 0;
    int biome[(int)Biome::Count] = {};
    int slope[kSlopeBins] = {};
    float spawnElev = 0;
    Biome spawnBiome = Biome::Ocean;
    PlanetField::SpawnScore spawn; // #1883: what the start has within 300 m
    double spawnSec = 0, riverSec = 0;
};

inline SeedStats Measure(unsigned seed, float radius) {
    SeedStats st;
    st.seed = seed;
    st.samples = kSamples;
    PlanetField field(seed, radius, false);
    const auto t0 = std::chrono::steady_clock::now();
    field.bakeErosion(""); // no cache: always the bake this code makes
    const auto t1 = std::chrono::steady_clock::now();
    st.bakeSec = std::chrono::duration<double>(t1 - t0).count();
    st.hash = field.erosionHash();

    const int threads = (int)std::max(1u, std::thread::hardware_concurrency());
    std::vector<SeedStats> part(threads);
    std::vector<std::thread> pool;
    const double golden = 3.14159265358979 * (3.0 - std::sqrt(5.0));
    const float stepAngle = kStepM / radius;
    for (int t = 0; t < threads; ++t) pool.emplace_back([&, t] {
        SeedStats& p = part[t];
        for (int i = t; i < kSamples; i += threads) {
            const double y = 1.0 - 2.0 * (i + 0.5) / kSamples, r = std::sqrt(std::max(0.0, 1.0 - y * y)), phi = golden * i;
            const Vector3 dir = { (float)(r * std::cos(phi)), (float)y, (float)(r * std::sin(phi)) };
            const float h = field.elevation(dir);
            if (h <= 0.0f) continue;
            ++p.land;
            ++p.biome[(int)field.biomeAt(dir, h)];
            // A tangent basis at dir (any: the histogram ignores direction).
            const Vector3 a = std::fabs(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 };
            Vector3 e1 = { a.y * dir.z - a.z * dir.y, a.z * dir.x - a.x * dir.z, a.x * dir.y - a.y * dir.x };
            const float l1 = std::sqrt(e1.x * e1.x + e1.y * e1.y + e1.z * e1.z);
            e1 = { e1.x / l1, e1.y / l1, e1.z / l1 };
            const Vector3 e2 = { dir.y * e1.z - dir.z * e1.y, dir.z * e1.x - dir.x * e1.z, dir.x * e1.y - dir.y * e1.x };
            auto at = [&](Vector3 t, float s) {
                Vector3 d = { dir.x + t.x * s, dir.y + t.y * s, dir.z + t.z * s };
                const float l = std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
                return field.elevation({ d.x / l, d.y / l, d.z / l });
            };
            const float hE = at(e1, stepAngle), hW = at(e1, -stepAngle), hN = at(e2, stepAngle), hS = at(e2, -stepAngle);
            const float gx = (hE - hW) / (2.0f * kStepM), gy = (hN - hS) / (2.0f * kStepM);
            const float deg = std::atan(std::sqrt(gx * gx + gy * gy)) * 57.2957795f;
            int bin = 0;
            while (bin < kSlopeBins - 1 && deg >= kSlopeEdges[bin]) ++bin;
            ++p.slope[bin];
            const float dev = std::max(std::max(std::fabs(hE - h), std::fabs(hW - h)), std::max(std::fabs(hN - h), std::fabs(hS - h)));
            if (dev <= kFlatTolM) ++p.flat;
        }
    });
    for (auto& th : pool) th.join();
    for (const SeedStats& p : part) {
        st.land += p.land; st.flat += p.flat;
        for (int b = 0; b < (int)Biome::Count; ++b) st.biome[b] += p.biome[b];
        for (int b = 0; b < kSlopeBins; ++b) st.slope[b] += p.slope[b];
    }
    st.sampleSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t1).count();
    // The river network is built on first use (any water query; a world's first chunks make it): timed apart from the
    // spawn search, which would otherwise carry it.
    const auto tr = std::chrono::steady_clock::now();
    (void)field.waterSourceAt({ 0, 1, 0 });
    const auto t2 = std::chrono::steady_clock::now();
    st.riverSec = std::chrono::duration<double>(t2 - tr).count();
    const Vector3 spawn = field.findSpawnDir(seed);
    st.spawnSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t2).count();
    st.spawnElev = field.elevation(spawn);
    st.spawnBiome = field.biomeAt(spawn, st.spawnElev);
    st.spawn = field.spawnScore(spawn);
    return st;
}

struct Run {
    std::mutex mx;
    std::vector<SeedStats> done;  // in seed order
    std::atomic<bool> started{ false }, finished{ false };
    size_t reported = 0;
};

inline void Report(ScenarioCtx& c, const SeedStats& st) {
    char hash[32];
    snprintf(hash, sizeof hash, "%016llx", (unsigned long long)st.hash);
    const std::string p = "seed " + std::to_string(st.seed) + " ";
    c.note(p + "bake hash " + hash + ", baked in " + std::to_string((int)st.bakeSec) + " s, sampled in " + std::to_string((int)st.sampleSec) + " s");
    c.metric(p + "bake hash low32", (double)(uint32_t)st.hash, "", false);
    c.metric(p + "bake", st.bakeSec, "s");
    const double land = std::max(st.land, 1);
    c.metric(p + "land", 100.0 * st.land / st.samples, "%", false);
    c.metric(p + "flat buildable", 100.0 * st.flat / land, "% of land", false);
    for (int b = 0; b < kSlopeBins; ++b) c.metric(p + "slope " + SlopeBinName(b) + " deg", 100.0 * st.slope[b] / land, "% of land", false);
    for (int b = (int)Biome::Jungle; b < (int)Biome::Count; ++b) c.metric(p + BiomeName((Biome)b), 100.0 * st.biome[b] / land, "% of land", false);
    c.metric(p + "spawn elevation", st.spawnElev, "m", false);
    c.note(p + "spawn biome " + BiomeName(st.spawnBiome));
    c.metric(p + "spawn flat within 300 m", 100.0 * st.spawn.flat, "%", false);
    c.metric(p + "spawn water distance", std::min(st.spawn.waterM, 1000.0f), "m (1000: none within 300 m)");
    c.metric(p + "spawn ore kinds within 300 m", st.spawn.ores, "of 4", false);
    c.metric(p + "spawn score", st.spawn.score, "", false);
    c.metric(p + "spawn search", st.spawnSec * 1000.0, "ms");
    c.metric(p + "river network build", st.riverSec * 1000.0, "ms");
}

// The scenario: seeds = how many planets (the world's own first, then terrain seeds 1, 2, ...).
inline void Define(ScenarioDef& s, int seeds) {
    s.world = ScenarioWorld::Sandbox; // (the flat test plane: nothing to stream; the planets are built in-process)
    s.record = false;                 // (nothing to watch)
    // Sim time runs at <= 0.83 x wall while waiting (a 20 ms sleep per 1/60 s frame): a bake is ~60-130 s.
    s.timeout = 120.0 + 280.0 * seeds;
    auto run = std::make_shared<Run>();
    s.step("bake and measure " + std::to_string(seeds) + (seeds == 1 ? " planet" : " planets"), s.timeout, [run, seeds](ScenarioCtx& c) {
        if (!run->started.exchange(true)) {
            std::vector<unsigned> list = { c.body->seed };
            for (unsigned k = 1; (int)list.size() < seeds; ++k) if (k != c.body->seed) list.push_back(k);
            const float radius = c.field->radius;
            std::thread([run, list, radius] {
                for (unsigned seed : list) {
                    SeedStats st = Measure(seed, radius);
                    std::lock_guard<std::mutex> lk(run->mx);
                    run->done.push_back(st);
                }
                run->finished.store(true);
            }).detach();
        }
        {
            std::lock_guard<std::mutex> lk(run->mx);
            while (run->reported < run->done.size()) Report(c, run->done[run->reported++]);
            c.waitingFor("planet " + std::to_string(run->done.size() + 1) + " of " + std::to_string(seeds) + " (bake, then samples)");
        }
        if (!run->finished.load()) { std::this_thread::sleep_for(std::chrono::milliseconds(20)); return false; }
        // The one check: no land biome vanished (across all the seeds).
        int total[(int)Biome::Count] = {};
        for (const SeedStats& st : run->done) for (int b = 0; b < (int)Biome::Count; ++b) total[b] += st.biome[b];
        for (int b = (int)Biome::Jungle; b < (int)Biome::Count; ++b) c.gt(std::string(BiomeName((Biome)b)) + " land samples", total[b], 0.5);
        return true;
    });
}

} // namespace worldgen_metrics
