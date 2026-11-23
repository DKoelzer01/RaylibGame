// perf/audio-water-cave (#1736, #1690): the Audio v2 water-and-caves update at a river bank by a cave mouth. Planet, audio v2
// on; the camera stands 1.7 m over the bank of the nearest river within 150 m of a cave mouth (a mouth and a river in one
// place, as the update scans rivers, simulated flow and the mouth's wind every 0.5-1 s), for 20 s after the warm-up.
// The stall was WaterAt(16 m) in the update's "look first": ~2 M cell reads a site. This scene measures the replacement
// (Fluid::AnyWaterNear) and the old call side by side at the player's spot, and checks they agree; the frame max shows any
// hitch the update causes (the sandbox may have no audio device: then the scene's own timings are the measurement).
// Limits live in tools/ci/perf_budgets.txt (runner only).
#include "perf_budget.h"
#include "fluid.h"
#include "profiler.h"
#include "rivers.h"
#include "voxel.h"
#include "waterfield.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <vector>

namespace {
struct St {
    BudgetFrames frames;
    bool found = false, both = false;
    Vector3 dir{}, up{};
    double nearMax = 0.0, nearSum = 0.0; int nearN = 0; // AnyWaterNear(16 m) per frame, ms
} S;

double Ms(std::chrono::steady_clock::time_point a) { return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - a).count(); }

void Stand(ScenarioCtx& c) {
    const float R = c.field->radius;
    const float e = std::max(0.0f, c.field->elevation(S.dir));
    const Vector3 eye = Vector3Scale(S.dir, R + e + 1.7f);
    const Vector3 side = Vector3Normalize(Vector3CrossProduct(S.dir, fabsf(S.dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    c.camera(eye, Vector3Add(eye, Vector3Scale(side, 20.0f)), S.dir);
}
}

SCENARIO("perf/audio-water-cave") {
    s.tags = { "perf", "budget" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.env = { { "RAYLIBGAME_AUDIO_V2", "1" } };
    s.timeout = 240.0;
    s.record = false;

    // An idle world, or 30 s at most: a slow headless box may never get there (#1690), and the search and the side-by-side
    // timing don't need it (warm-up settles again before the per-frame measure).
    s.step("settle", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("2 s of play and an idle world (30 s at most)");
        if (c.T >= 30.0 && !WorldIdle()) { c.note("settle: the world was not idle after 30 s; going on"); return true; }
        return c.T >= 2.0 && WorldIdle();
    });
    s.step("find a river bank by a cave mouth", 5.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        const RiverNetwork* net = c.field->rivers();
        std::vector<Vector3> mouths;
        float bestD = 1e30f;
        for (int kind : { 0, 2, 1 }) { // hillside mouths, resurgences, shafts
            c.field->caveMouthsList(home, 4000.0f, kind, mouths);
            for (const Vector3& m : mouths) {
                const Vector3 dir = Vector3Normalize(m);
                if (!net) break;
                const RiverNetwork::Hit h = net->nearest(dir, 150.0f);
                if (!h.found || h.piece < 0 || h.depth <= 0.05f) continue;
                const float d = Vector3Distance(Vector3Scale(dir, c.field->radius), Vector3Scale(home, c.field->radius));
                if (d < bestD) { bestD = d; S.dir = dir; S.found = S.both = true; }
            }
        }
        if (!S.found) { // no river by a mouth: the nearest mouth alone (the river scan then idles; the flow scan and wind still run)
            Vector3 at{}, out{};
            if (c.field->nearestCaveMouth(home, 5000.0f, at, out)) { S.dir = Vector3Normalize(at); S.found = true; }
        }
        if (!c.check("a cave mouth within 4 km of the start", S.found)) c.fail("none found");
        c.note(S.both ? "river and cave mouth found together" : "no river within 150 m of a mouth: mouth only");
        return true;
    });
    s.step("warm-up", kBudgetSettleCap + 5.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        Stand(c);
        c.waitingFor("the world to settle");
        static BudgetSettle settle;
        return settle.done(c);
    });
    s.step("old call against new", 30.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        Stand(c);
        const Vector3 here = Vector3Scale(S.dir, c.field->radius + std::max(0.0f, c.field->elevation(S.dir)) + 1.7f);
        const auto t0 = std::chrono::steady_clock::now();
        const float wat = Fluid::WaterAt(*c.body, here, 16.0f);
        const double oldMs = Ms(t0);
        const auto t1 = std::chrono::steady_clock::now();
        const bool any = Fluid::AnyWaterNear(*c.body, here, 16.0f, 0.2f);
        const double newMs = Ms(t1);
        c.metric("WaterAt(16 m), the old look", oldMs, "ms", false);
        c.metric("AnyWaterNear(16 m), the new look", newMs, "ms", false);
        c.metric("water within 16 m (m^3)", wat, "m^3", false);
        c.check("AnyWaterNear agrees with WaterAt > 0.2 (a budget overrun may say yes to a no)", any || !(wat > 0.2f));
        c.lt("AnyWaterNear(16 m) ms", newMs, 2.0);
        return true;
    });
    // The dense case: inside simulated water every brick within 16 m exists, so a brick walk still reads the whole cube.
    // The wet brick with the most bricks within 4 (the 16 m cube at 0.25 m cells) anywhere on the body; both calls at its centre.
    s.step("old call against new, in water", 5.0, [](ScenarioCtx& c) {
        const WaterField::SnapshotPtr snap = WaterField::Current();
        const WaterField::Site* best = nullptr; int bi = 0, bj = 0, bk = 0, bestN = -1;
        for (const WaterField::Site& s : snap->sites) {
            if (s.body != c.body->id) continue;
            for (const auto& kv : s.bricks) {
                if (!kv.second || !(kv.second->kinds & 1)) continue;
                const int i = (int)(kv.first >> 22) - 1024, j = (int)((kv.first >> 11) & 2047) - 1024, k = (int)(kv.first & 2047) - 1024;
                int n = 0;
                for (int dk = -4; dk <= 4; ++dk) for (int dj = -4; dj <= 4; ++dj) for (int di = -4; di <= 4; ++di) n += s.brick(i + di, j + dj, k + dk) != nullptr;
                if (n > bestN) { bestN = n; best = &s; bi = i; bj = j; bk = k; }
            }
        }
        if (!best) { c.note("no simulated water on the body: the in-water case is skipped"); return true; }
        Vector3 lo, hi;
        best->brickBounds(bi, bj, bk, lo, hi);
        const Vector3 here = Vector3Scale(Vector3Add(lo, hi), 0.5f);
        const auto t0 = std::chrono::steady_clock::now();
        const float wat = Fluid::WaterAt(*c.body, here, 16.0f);
        const double oldMs = Ms(t0);
        const auto t1 = std::chrono::steady_clock::now();
        const bool any = Fluid::AnyWaterNear(*c.body, here, 16.0f, 0.2f);
        const double newMs = Ms(t1);
        c.metric("in water: bricks in the 16 m cube (of 729)", bestN, "", false);
        c.metric("in water: water within 16 m (m^3)", wat, "m^3", false);
        c.metric("in water: WaterAt(16 m), the old look", oldMs, "ms", false);
        c.metric("in water: AnyWaterNear(16 m), the new look", newMs, "ms", false);
        c.check("in water: AnyWaterNear agrees with WaterAt > 0.2", any || !(wat > 0.2f));
        c.lt("in water: AnyWaterNear(16 m) ms", newMs, 2.0);
        return true;
    });
    s.step("measure", kBudgetMeasure + 10.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        if (c.t < 0.02) Profiler_Reset();
        Stand(c);
        S.frames.frame();
        const Vector3 here = Vector3Scale(S.dir, c.field->radius + std::max(0.0f, c.field->elevation(S.dir)) + 1.7f);
        const auto t0 = std::chrono::steady_clock::now();
        Fluid::AnyWaterNear(*c.body, here, 16.0f, 0.2f);
        const double ms = Ms(t0);
        S.nearMax = std::max(S.nearMax, ms); S.nearSum += ms; ++S.nearN;
        c.waitingFor("the timed window");
        if (c.t <= kBudgetMeasure) return false;
        S.frames.report(c);
        c.metric("AnyWaterNear per-frame max", S.nearMax, "ms");
        c.metric("audio:water and caves zone mean", std::max(0.0, Profiler_ZoneMeanMs("audio:water and caves")), "ms");
        c.metric("audio:flow zone mean", std::max(0.0, Profiler_ZoneMeanMs("audio:flow")), "ms", false);
        c.metric("audio:rivers zone mean", std::max(0.0, Profiler_ZoneMeanMs("audio:rivers")), "ms", false);
        c.metric("audio:caves zone mean", std::max(0.0, Profiler_ZoneMeanMs("audio:caves")), "ms", false);
        return true;
    });
}
