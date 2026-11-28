// Perf scenarios (#871): a fixed view held for a timed window; the frame times are measured by the scenario itself (so a
// run with the profiler off has them too: the A/B overhead check compares exactly that), and the profiler's summary goes
// into results.json. Perf policy (user): always reported; the report fails a run only on a large regression against the
// recent passing runs (+15% mean frame time, or a new hitch over 30 ms).
#include "scenario.h"
#include "biomevariants.h"
#include "fluid.h"
#include "profiler.h"
#include "trees.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <string>
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

// (#1995 phase 3) Dense flora at eye height: the start's bushes, stones, grass and trees filling the view (the small
// flora's size classes and their instancing budget, before/after).
SCENARIO("perf/flora-view") {
    s.tags = { "perf", "planet", "foliage" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.env = { { "RAYLIBGAME_WEATHER_FORCE", "0" } };
    s.timeout = 120.0;
    s.record = false;
    PerfView(s, 20.0, 1.8f, 60.0f);
}

// (#1995 phase 4) A view over a forest: the densest old-growth patch near the start, and the heading (of 16) and lowest
// eye height (70-370 m) whose sight line runs clear of the terrain to a point 1.5 km out (a patch in hills looked
// straight into a slope 50 m off) over the most trees (the clearest line otherwise ran out over the sea).
namespace {
struct ForestView { Vector3 up{}, ax{}; float ground = 0.0f, height = 70.0f; };
bool FindForestView(ScenarioCtx& c, ForestView& fv) {
    const int v = BiomeVariants().find("forest.old_growth");
    const Vector3 home = c.field->findSpawnDir(c.body->seed);
    Vector3 best = home;
    int most = -1;
    if (v >= 0) {
        const Vector3 base = c.field->nearestVariantDir(home, v);
        const Vector3 b1 = Vector3Normalize(Vector3CrossProduct(base, fabsf(base.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })), b2 = Vector3CrossProduct(base, b1);
        for (int j = 0; j < 17; ++j) { // the patch and points 300 / 600 m around it: the one growing the most trees
            const float a = (j - 1) * 0.785398f, r = j == 0 ? 0.0f : (j <= 8 ? 300.0f : 600.0f) / c.field->radius;
            const Vector3 cand = Vector3Normalize(Vector3Add(base, Vector3Add(Vector3Scale(b1, cosf(a) * r), Vector3Scale(b2, sinf(a) * r))));
            int t = 0, b = 0, sn = 0;
            Trees_CountScatter(*c.field, c.field->foliageSeed, cand, 96.0f, t, b, sn);
            if (t > most) { most = t; best = cand; }
        }
    }
    c.note(TextFormat("forest patch: %d trees generated within 96 m%s", most, v < 0 ? " (no forest.old_growth variant: the start)" : ""));
    fv.up = best;
    const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(fv.up, fabsf(fv.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })), e2 = Vector3CrossProduct(fv.up, e1);
    // The ground by density, down from 400 m over the elevation field (cliff tops stand well above it).
    float r = c.field->radius + std::max(0.0f, c.field->elevation(fv.up));
    for (float y = r + 400.0f; y > r - 100.0f; y -= 2.0f)
        if (c.field->density(Vector3Scale(fv.up, y), nullptr) >= ISO_LEVEL) { r = std::max(y, c.field->radius); break; }
    fv.ground = r;
    int along[16] = {};
    for (int d = 0; d < 16; ++d)
        for (float m : { 300.0f, 700.0f, 1100.0f, 1500.0f }) {
            const float a = d * 0.3926991f;
            const Vector3 p = Vector3Normalize(Vector3Add(fv.up, Vector3Scale(Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a))), m / c.field->radius)));
            int t = 0, b = 0, sn = 0;
            Trees_CountScatter(*c.field, c.field->foliageSeed, p, 64.0f, t, b, sn);
            along[d] += t;
        }
    float bestClear = -1.0f;
    int bestTrees = -1;
    for (float h = 70.0f; h <= 370.0f && bestClear < 1200.0f; h += 50.0f)
        for (int d = 0; d < 16; ++d) {
            const float a = d * 0.3926991f;
            const Vector3 ax = Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a)));
            const Vector3 eye = Vector3Scale(fv.up, r + h), at = Vector3Add(Vector3Scale(fv.up, r), Vector3Scale(ax, 1500.0f));
            float clear = 1500.0f;
            for (int i = 1; i <= 60; ++i) // rock by density (cliffs and tors are not in the elevation field)
                if (c.field->density(Vector3Lerp(eye, at, i / 60.0f), nullptr) >= ISO_LEVEL) { clear = 1500.0f * (i - 1) / 60.0f; break; }
            const bool better = clear >= 1200.0f ? (bestClear < 1200.0f || along[d] > bestTrees) : clear > bestClear;
            if (better) { bestClear = clear; bestTrees = along[d]; fv.ax = ax; fv.height = h; }
        }
    c.note(TextFormat("view: eye %.0f m over the patch, sight line clear for %.0f m, %d trees along it", fv.height, bestClear, bestTrees));
    return c.gt("trees in the patch", (double)most, 10.0);
}
}

// (#1995 phase 4) The long view over a forest: from FindForestView's eye, looking 1.5 km across the forest, so the 2 km
// tree distance is full of far trees (the far LOD budget, before/after).
SCENARIO("perf/forest-view") {
    s.tags = { "perf", "planet", "foliage" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--no-hud" };
    s.env = { { "RAYLIBGAME_WEATHER_FORCE", "0" } };
    s.timeout = 180.0;
    s.record = false;
    struct F { Timing timing; Vector3 eye{}, at{}, up{}; };
    static F f;
    s.step("find a forest", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        ForestView fv;
        FindForestView(c, fv);
        f.up = fv.up;
        f.eye = Vector3Scale(fv.up, fv.ground + fv.height);
        f.at = Vector3Add(Vector3Scale(fv.up, fv.ground), Vector3Scale(fv.ax, 1500.0f));
        c.camera(f.eye, f.at, f.up);
        return true;
    });
    s.step("warm", 60.0, [](ScenarioCtx& c) {
        c.camera(f.eye, f.at, f.up);
        c.waitingFor("streaming to settle");
        return c.t > 5.0 && (WorldIdle() || c.t > 45.0);
    });
    s.step("measure", 25.0, [](ScenarioCtx& c) {
        if (c.t < 0.02) Profiler_Reset();
        c.camera(f.eye, f.at, f.up);
        f.timing.frame();
        c.waitingFor("the timed window");
        if (c.t <= 20.0) return false;
        if (c.body->trees) {
            c.metric("tree instances drawn", (double)c.body->trees->instanceCount(), "");
            c.metric("tree instances visible (after the GPU cull)", (double)c.body->trees->visibleInstances(), "");
        }
        c.metric("gpu:trees", std::max(0.0, Profiler_ZoneMeanMs("gpu:trees")), "ms");
        f.timing.report(c);
        c.shot("forest view", 100.0); // after the timed window: what was measured, to judge the far trees by eye
        return true;
    });
}

// (#1995 phase 4, the user 2026-11-28: can simple far trees reach the horizon?) The same forest seen toward the geometric
// horizon from the ground (2 m), 200 m and 1500 m up. Per height: the tiles' load (time, tile builds and their cost),
// then 10 s timed: frame times, the GPU passes, tree instances and tiles, memory. Run it with --tree-dist N to cost a
// farther tree distance (2000 m today; the horizon is ~490 m from 2 m, ~4.9 km from 200 m, ~13.4 km from 1500 m).
SCENARIO("perf/horizon-view") {
    s.tags = { "perf", "planet", "foliage" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--no-hud" };
    s.env = { { "RAYLIBGAME_WEATHER_FORCE", "0" } };
    s.timeout = 420.0;
    s.record = false;
    struct H { ForestView fv; Vector3 eye{}, at{}; int tiles0 = 0; };
    static H h;
    s.step("find a forest", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        FindForestView(c, h.fv);
        return true;
    });
    for (float height : { 2.0f, 200.0f, 1500.0f }) {
        const std::string tag = TextFormat("%.0f m", height);
        s.step("to " + tag, 1.0, [height](ScenarioCtx& c) {
            const float R = h.fv.ground, dip = acosf(R / (R + height)); // the geometric horizon's angle below level
            h.eye = Vector3Scale(h.fv.up, R + height);
            h.at = Vector3Add(h.eye, Vector3Scale(Vector3Subtract(Vector3Scale(h.fv.ax, cosf(dip)), Vector3Scale(h.fv.up, sinf(dip))), 1000.0f));
            h.tiles0 = c.body->trees ? c.body->trees->loadedTiles() : 0;
            c.camera(h.eye, h.at, h.fv.up);
            Profiler_Reset();
            return true;
        });
        s.step("load " + tag, 90.0, [tag](ScenarioCtx& c) {
            c.camera(h.eye, h.at, h.fv.up);
            c.waitingFor("streaming to settle");
            if (c.t < 3.0 || (!WorldIdle() && c.t < 80.0)) return false;
            c.metric(tag + ": load time", c.t, "s");
            const int tiles = c.body->trees ? c.body->trees->loadedTiles() : 0;
            c.metric(tag + ": tree tiles added", (double)(tiles - h.tiles0), "");
            c.metric(tag + ": tree tile build (one job)", std::max(0.0, Profiler_ZoneMeanMs("trees: build tile")), "ms");
            c.metric(tag + ": tree batch rebuild (main thread)", std::max(0.0, Profiler_ZoneMeanMs("trees: batch start")) + std::max(0.0, Profiler_ZoneMeanMs("trees: batch upload")), "ms");
            Profiler_Reset();
            return true;
        });
        s.step("measure " + tag, 15.0, [tag](ScenarioCtx& c) {
            static Timing t;
            if (c.t < 0.02) { t = Timing(); Profiler_Reset(); }
            c.camera(h.eye, h.at, h.fv.up);
            t.frame();
            c.waitingFor("the timed window");
            if (c.t <= 10.0) return false;
            std::vector<double> s = t.ms;
            std::sort(s.begin(), s.end());
            double sum = 0.0; for (double v : s) sum += v;
            c.metric(tag + ": frame mean", s.empty() ? 0.0 : sum / s.size(), "ms");
            c.metric(tag + ": frame p95", s.empty() ? 0.0 : s[std::min(s.size() - 1, (size_t)(0.95 * s.size()))], "ms");
            for (const char* pass : { "gpu:trees", "gpu:shadow", "gpu:near", "gpu:far", "gpu:grass", "Main pass" })
                c.metric(tag + ": " + pass, std::max(0.0, Profiler_ZoneMeanMs(pass)), "ms");
            if (c.body->trees) {
                c.metric(tag + ": tree instances drawn", (double)c.body->trees->instanceCount(), "");
                c.metric(tag + ": tree tiles loaded", (double)c.body->trees->loadedTiles(), "");
            }
            c.metric(tag + ": commit", Profiler_CommitMB(), "MB");
            c.metric(tag + ": VRAM used", Profiler_VramUsedMB(), "MB");
            c.shot("horizon " + tag, 100.0);
            return true;
        });
    }
    // (#2133) Flying: 200 m up at 100 m/s for 10 s along the view (the tree batches rebuild every 40 m of camera movement;
    // tiles stream in ahead). Frame p95/p99/max and the batch rebuild's parts.
    s.step("fly 200 m", 20.0, [](ScenarioCtx& c) {
        static Timing t;
        static Vector3 start{};
        if (c.t < 0.02) { t = Timing(); Profiler_Reset(); start = Vector3Scale(h.fv.up, h.fv.ground + 200.0f); }
        const Vector3 eye = Vector3Add(start, Vector3Scale(h.fv.ax, 100.0f * (float)c.t));
        const Vector3 up = Vector3Normalize(eye);
        c.camera(Vector3Scale(up, h.fv.ground + 200.0f), Vector3Add(Vector3Scale(up, h.fv.ground + 150.0f), Vector3Scale(h.fv.ax, 1000.0f)), up);
        t.frame();
        c.waitingFor("the flight");
        if (c.t <= 10.0) return false;
        std::vector<double> s = t.ms;
        std::sort(s.begin(), s.end());
        auto pct = [&](double q) { return s.empty() ? 0.0 : s[std::min(s.size() - 1, (size_t)(q * s.size()))]; };
        c.metric("fly: frame p50", pct(0.5), "ms");
        c.metric("fly: frame p95", pct(0.95), "ms");
        c.metric("fly: frame p99", pct(0.99), "ms");
        c.metric("fly: frame max", s.empty() ? 0.0 : s.back(), "ms");
        for (const char* z : { "trees: batch start", "trees: batch sort", "trees: batch upload" })
            c.metric(std::string("fly: ") + z, std::max(0.0, Profiler_ZoneMeanMs(z)), "ms");
        return true;
    });
}

// (#2141) A view where a hill hides a forest: a ridge crest near FindForestView's patch with the most trees 300-1500 m
// behind it (the ground falling away there), the eye 150 m in front of the crest and 6 m below its line, looking over it
// at the forest. 10 s timed: gpu:trees and the instances the GPU cull keeps, the Hi-Z occlusion's saving
// (RAYLIBGAME_TREE_HIZ=0 for the A/B). Crests are often cliff faces: the view need not be pretty, only hide the forest.
SCENARIO("perf/hill-view") {
    s.tags = { "perf", "planet", "foliage" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--no-hud" };
    s.env = { { "RAYLIBGAME_WEATHER_FORCE", "0" } };
    s.timeout = 240.0;
    s.record = false;
    struct V { Vector3 up{}, ax{}, eye{}, at{}; float hit = 0.0f; Timing timing; };
    static V v;
    s.step("find a hidden forest", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        ForestView fv;
        FindForestView(c, fv);
        const float R = c.field->radius;
        // The ground by density (cliffs and tors are not in the elevation field), down from 300 m over the elevation.
        auto groundAt = [&](Vector3 dir) {
            const float g = R + std::max(0.0f, c.field->elevation(dir));
            for (float y = g + 300.0f; y > g - 100.0f; y -= 1.0f)
                if (c.field->density(Vector3Scale(dir, y), nullptr) >= ISO_LEVEL) return y;
            return g;
        };
        const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(fv.up, fabsf(fv.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })), e2 = Vector3CrossProduct(fv.up, e1);
        float best = 0.0f;
        for (int i = -4; i <= 4; ++i)
            for (int j = -4; j <= 4; ++j) { // crest candidates on a 400 m grid round the patch
                const Vector3 crest = Vector3Normalize(Vector3Add(fv.up, Vector3Scale(Vector3Add(Vector3Scale(e1, (float)i), Vector3Scale(e2, (float)j)), 400.0f / R)));
                if (c.field->elevation(crest) < 3.0f) continue;
                const float top = R + c.field->elevation(crest);
                const Vector3 f1 = Vector3Normalize(Vector3CrossProduct(crest, fabsf(crest.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })), f2 = Vector3CrossProduct(crest, f1);
                for (int d = 0; d < 16; ++d) {
                    const float a = d * 0.3926991f;
                    const Vector3 ax = Vector3Add(Vector3Scale(f1, cosf(a)), Vector3Scale(f2, sinf(a)));
                    auto dirAt = [&](float m) { return Vector3Normalize(Vector3Add(crest, Vector3Scale(ax, m / R))); };
                    auto elev = [&](float m) { return R + c.field->elevation(dirAt(m)); };
                    if (elev(-150.0f) > top - 4.0f) continue; // the eye's side lies below the crest
                    bool falls = true;                          // and the far side falls away: hidden from the eye
                    for (float m = 150.0f; m <= 450.0f && falls; m += 100.0f) falls = elev(m) < top - 8.0f;
                    if (!falls) continue;
                    int trees = 0;
                    for (float m : { 300.0f, 500.0f, 800.0f, 1100.0f, 1400.0f }) {
                        int t = 0, b = 0, sn = 0;
                        Trees_CountScatter(*c.field, c.field->foliageSeed, dirAt(m), 96.0f, t, b, sn);
                        trees += t;
                    }
                    if ((float)trees <= best) continue;
                    // The eye 150 m before the crest, 6 m under its line, in the open (by density) and its sight line to
                    // the crest clear of rock.
                    const Vector3 eyeDir = dirAt(-150.0f);
                    const float crestTop = groundAt(crest);
                    const Vector3 eye = Vector3Scale(eyeDir, std::max(crestTop - 6.0f, groundAt(eyeDir) + 2.0f));
                    const Vector3 toCrest = Vector3Subtract(Vector3Scale(crest, crestTop - 2.0f), eye);
                    bool open = true;
                    for (float k = 0.0f; k <= 0.9f && open; k += 0.03f) open = c.field->density(Vector3Add(eye, Vector3Scale(toCrest, k)), nullptr) < ISO_LEVEL;
                    if (!open) continue;
                    // The view's centre line (to the forest) meets rock at the crest, not a wall close by.
                    const Vector3 at = Vector3Scale(dirAt(800.0f), R + c.field->elevation(dirAt(800.0f)) + 10.0f);
                    const Vector3 look = Vector3Normalize(Vector3Subtract(at, eye));
                    float hit = 0.0f;
                    while (hit < 400.0f && c.field->density(Vector3Add(eye, Vector3Scale(look, hit)), nullptr) < ISO_LEVEL) hit += 2.0f;
                    if (hit < 120.0f) continue;
                    best = (float)trees; v.up = eyeDir; v.ax = ax; v.eye = eye; v.at = at;
                    v.hit = hit;
                }
            }
        if (best <= 0.0f) { c.note("no ridge hides a forest near the patch: the forest view's eye at 2 m"); v.up = fv.up; v.ax = fv.ax; v.eye = Vector3Scale(fv.up, fv.ground + 2.0f); v.at = Vector3Add(v.eye, Vector3Scale(fv.ax, 1000.0f)); }
        c.note(TextFormat("hill view: %.0f trees generated behind the crest, the view meets rock at %.0f m", best, v.hit));
        c.camera(v.eye, v.at, v.up);
        c.gt("trees behind the crest", best, 20.0);
        return true;
    });
    s.step("warm", 60.0, [](ScenarioCtx& c) {
        c.camera(v.eye, v.at, v.up);
        c.waitingFor("streaming to settle");
        return c.t > 5.0 && (WorldIdle() || c.t > 45.0);
    });
    s.step("measure", 15.0, [](ScenarioCtx& c) {
        if (c.t < 0.02) Profiler_Reset();
        c.camera(v.eye, v.at, v.up);
        v.timing.frame();
        c.waitingFor("the timed window");
        if (c.t <= 10.0) return false;
        if (c.body->trees) {
            c.metric("tree instances drawn", (double)c.body->trees->instanceCount(), "");
            c.metric("tree instances visible (after the GPU cull)", (double)c.body->trees->visibleInstances(), "");
        }
        c.metric("gpu:trees", std::max(0.0, Profiler_ZoneMeanMs("gpu:trees")), "ms");
        v.timing.report(c);
        c.shot("hill view", 100.0);
        return true;
    });
}

// (#2141) Hill edges sweeping across trees, recorded: from the forest view's eye (70 m over the patch, where the terrain
// hides about half the forest's trees) the camera sinks to 12 m in 6 s, looking at the same point 1.5 km out, so ridge
// lines climb across the forest behind them. Run it with RAYLIBGAME_TREE_HIZ=0 and =1 (--record-compare): the pyramid is
// a frame old, so a tree must not pop in or out at an edge. The recording starts at the sink (--record-on-cue): runs line
// up frame by frame.
SCENARIO("visual/hill-edge-trees") {
    s.tags = { "visual", "planet", "foliage" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--no-hud", "--record-on-cue", "--record-every", "1" };
    s.env = { { "RAYLIBGAME_WEATHER_FORCE", "0" } };
    s.timeout = 180.0;
    struct E { Vector3 up{}, at{}; float ground = 0.0f; };
    static E e;
    auto eyeAt = [](float h) { return Vector3Scale(e.up, e.ground + h); };
    s.step("find a forest", 60.0, [eyeAt](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        ForestView fv;
        FindForestView(c, fv);
        e.up = fv.up;
        e.ground = fv.ground;
        e.at = Vector3Add(Vector3Scale(fv.up, fv.ground), Vector3Scale(fv.ax, 1500.0f));
        c.camera(eyeAt(70.0f), e.at, e.up);
        return true;
    });
    s.step("warm", 60.0, [eyeAt](ScenarioCtx& c) {
        c.camera(eyeAt(70.0f), e.at, e.up);
        c.waitingFor("streaming to settle");
        return c.t > 5.0 && (WorldIdle() || c.t > 45.0);
    });
    s.step("sink", 8.0, [eyeAt](ScenarioCtx& c) {
        if (c.t < 0.02) c.cue();
        c.camera(eyeAt(70.0f - 58.0f * (float)std::min(c.t / 6.0, 1.0)), e.at, e.up);
        c.waitingFor("the sink");
        if (c.t <= 6.0) return false;
        c.keepCapture("hill edges sweeping across the forest: trees must not pop in or out (#2141)");
        return true;
    });
}
