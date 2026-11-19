// perf/profiler-overhead (#871): what a profiling zone costs. The zone microbenchmark (200k zones on this thread, on and
// off) and the zones the run actually records a frame, which together bound the profiler's share of a frame; the frame-level
// check (profiling off vs on, < 1% of the mean frame) is the A/B runner's:
//   game.exe --scenario-ab perf/planet-view --ab-a RAYLIBGAME_PROFILE=0 --ab-b RAYLIBGAME_PROFILE=1 --ab-runs 3 --ab-max-delta 1
// (in a quiet window: CLAUDE.md).
#include "scenario.h"
#include "profiler.h"
#include "world.h"
#include <chrono>
#include <cstdio>

namespace {
struct State { long frames = 0; double frameMsSum = 0.0; std::chrono::steady_clock::time_point last{}; } S;
}

SCENARIO("perf/profiler-overhead") {
    s.tags = { "perf", "profiler" };
    s.record = false;
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 2.0 && WorldIdle(); });
    s.step("micro", 5.0, [](ScenarioCtx& c) {
        double onNs = 0.0, offNs = 0.0;
        Profiler_MeasureZoneCost(onNs, offNs);
        c.metric("zone cost, profiling on", onNs, "ns");
        c.metric("zone cost, profiling off", offNs, "ns");
        c.lt("zone cost on (ns; plan: 20-40)", onNs, 60.0);
        c.lt("zone cost off (ns)", offNs, 5.0);
        return true;
    });
    // The zones recorded per frame in a running world, times their cost: the profiler's share of a frame.
    s.step("frames", 15.0, [](ScenarioCtx& c) {
        const auto now = std::chrono::steady_clock::now();
        if (c.t < 0.02) { Profiler_Reset(); S.last = now; return false; }
        S.frameMsSum += std::chrono::duration<double, std::milli>(now - S.last).count(); ++S.frames; S.last = now;
        c.waitingFor("10 s of frames");
        if (c.t <= 10.0) return false;
        const std::string j = Profiler_SummaryJson();
        // Zones recorded: the sum of the zones' counts (the "count" fields of the summary's zone list, cpu:/gpu: samples excluded).
        double zones = 0.0;
        for (size_t at = 0; (at = j.find("{\"name\": \"", at)) != std::string::npos; ) {
            at += 10;
            const bool sample = j.compare(at, 4, "cpu:") == 0 || j.compare(at, 4, "gpu:") == 0;
            const size_t k = j.find("\"count\": ", at);
            if (k == std::string::npos) break;
            if (!sample) zones += atof(j.c_str() + k + 9);
        }
        double onNs = 0.0, offNs = 0.0;
        Profiler_MeasureZoneCost(onNs, offNs);
        const double perFrame = S.frames ? zones / S.frames : 0.0, meanMs = S.frames ? S.frameMsSum / S.frames : 0.0;
        const double sharePct = meanMs > 0.0 ? 100.0 * perFrame * onNs / 1e6 / meanMs : 0.0;
        c.metric("zones per frame", perFrame, "");
        c.metric("frame mean", meanMs, "ms");
        c.metric("profiler share of a frame", sharePct, "%");
        c.lt("profiler share of a frame (%, from zone counts)", sharePct, 1.0);
        return true;
    });
}
