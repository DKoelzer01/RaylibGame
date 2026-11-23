// Shared by the budget scenes (#1690): frame times of a timed window, reported as "frame mean" and "frame p95" (the names
// tools/ci/perf_budgets.txt keys on), with the scene's own cost reported by the scene.
#pragma once
#include "scenario.h"
#include <algorithm>
#include <chrono>
#include <vector>

struct BudgetFrames {
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
        c.metric("frame mean", sum / s.size(), "ms");
        c.metric("frame p95", s[std::min(s.size() - 1, (size_t)(0.95 * s.size()))], "ms");
        c.metric("frame max", s.back(), "ms");
        c.gt("frames timed", (double)s.size(), 60.0);
    }
};

constexpr double kBudgetWarmup = 5.0, kBudgetMeasure = 20.0; // s: not measured, then the camera path
