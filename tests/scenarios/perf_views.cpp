// Perf scenarios (#871): a fixed view held for a timed window; the frame times are measured by the scenario itself (so a
// run with the profiler off has them too: the A/B overhead check compares exactly that), and the profiler's summary goes
// into results.json. Perf policy (user): always reported; the report fails a run only on a large regression against the
// recent passing runs (+15% mean frame time, or a new hitch over 30 ms).
#include "scenario.h"
#include "fluid.h"
#include "profiler.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <vector>

namespace {
struct Timing {
    std::vector<double> ms;
    std::chrono::steady_clock::time_point last{};
    bool started = false;
    void frame() {
        const auto now = std::chrono::steady_clock::now();
        if (started) ms.push_back(std::chrono::duration<double, std::milli>(now - last).count());
        last = now; started = true;
    }
    void report(ScenarioCtx& c) const {
        if (ms.empty()) { c.check("frames timed", false); return; }
        std::vector<double> s = ms;
        std::sort(s.begin(), s.end());
        double sum = 0.0; for (double v : s) sum += v;
        auto pct = [&](double q) { return s[std::min(s.size() - 1, (size_t)(q * s.size()))]; };
        c.metric("frame mean", sum / s.size(), "ms");
        c.metric("frame p50", pct(0.5), "ms");
        c.metric("frame p95", pct(0.95), "ms");
        c.metric("frame p99", pct(0.99), "ms");
        c.metric("frame max", s.back(), "ms");
        c.metric("frames over 30 ms", (double)(s.end() - std::upper_bound(s.begin(), s.end(), 30.0)), "");
        c.gt("frames timed", (double)s.size(), 60.0);
    }
};

// Shared shape: settle (the world idle), then hold a view for `seconds` of sim time, timing every frame.
void PerfView(ScenarioDef& s, double seconds, float height, float lookAhead) {
    struct V { Timing timing; Vector3 eye{}, at{}, up{}; };
    static V v;
    s.step("settle", 60.0, [height, lookAhead](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        v.up = c.spawnDirLocal();
        const Vector3 ax = Vector3Normalize(Vector3CrossProduct(v.up, { 0, 1, 0 }));
        const float r = c.field->radius + std::max(0.0f, c.field->elevation(v.up));
        v.eye = Vector3Scale(v.up, r + height);
        v.at = Vector3Add(Vector3Scale(v.up, r), Vector3Scale(ax, lookAhead));
        c.camera(v.eye, v.at, v.up);
        return true;
    });
    s.step("warm", 10.0, [](ScenarioCtx& c) { c.camera(v.eye, v.at, v.up); c.waitingFor("streaming to settle"); return c.t > 3.0 && WorldIdle(); });
    s.step("measure", seconds + 5.0, [seconds](ScenarioCtx& c) {
        if (c.t < 0.02) Profiler_Reset(); // the summary covers the timed window only
        c.camera(v.eye, v.at, v.up);
        v.timing.frame();
        c.waitingFor("the timed window");
        if (c.t <= seconds) return false;
        v.timing.report(c);
        return true;
    });
}
}

SCENARIO("perf/sandbox-view") {
    s.tags = { "perf" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.legacyFlag = "--sandbox --bench 10";
    s.record = false;
    PerfView(s, 15.0, 40.0f, 200.0f);
}

SCENARIO("perf/planet-view") {
    s.tags = { "perf", "planet", "budget" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.legacyFlag = "--seed 42 --bench 10";
    s.record = false;
    PerfView(s, 20.0, 40.0f, 300.0f); // as a plain planet --bench: 40 m over the start, looking along the ground
}
