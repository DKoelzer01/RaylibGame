// perf/water-shore (#1690): a budget scene. Planet; the nearest sea shore (found as in water/beach), the camera 1.7 m over
// the beach walking along the waterline at 3 m/s for 20 s after 5 s of warm-up (waves, shore foam, the fluid bricks of the
// shore and any river mouth in view). Frame mean and p95, and the water render pass's CPU and GPU time (the scene's own cost).
// Limits live in tools/ci/perf_budgets.txt (runner only).
#include "perf_budget.h"
#include "profiler.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
struct St { BudgetFrames frames; bool found = false; Vector3 sea{}, up{}; } S;
void Walk(ScenarioCtx& c, double t) {
    const float R = c.field->radius;
    const Vector3 up = Vector3Normalize(S.up);
    const Vector3 out = Vector3Normalize(Vector3Subtract(S.sea, S.up));
    const Vector3 along = Vector3Normalize(Vector3CrossProduct(up, out));
    const float d = (float)(t * 3.0 - 20.0); // m along the shore, from -20 to +60
    const Vector3 beach = Vector3Normalize(Vector3Add(Vector3Scale(S.up, R), Vector3Add(Vector3Scale(along, d), Vector3Scale(out, -8.0f))));
    const float e = std::max(0.0f, c.field->elevation(beach));
    const Vector3 eye = Vector3Scale(beach, R + e + 1.7f);
    const Vector3 target = Vector3Scale(Vector3Normalize(Vector3Add(Vector3Scale(beach, R), Vector3Add(Vector3Scale(along, 25.0f), Vector3Scale(out, 15.0f)))), R);
    c.camera(eye, target, up);
}
}

SCENARIO("perf/water-shore") {
    s.tags = { "perf", "budget" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 180.0;
    s.record = false;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("find the nearest shore", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        const Vector3 ax = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), ay = Vector3CrossProduct(ax, home);
        const float R = c.field->radius;
        float best = 1e30f;
        for (int k = 0; k < 16; ++k) {
            const float a = 6.2831853f * k / 16.0f;
            const Vector3 dirT = Vector3Add(Vector3Scale(ax, cosf(a)), Vector3Scale(ay, sinf(a)));
            Vector3 prev = home;
            for (float d = 10.0f; d < 4000.0f && d < best; d += 10.0f) {
                const Vector3 dd = Vector3Normalize(Vector3Add(Vector3Scale(home, R), Vector3Scale(dirT, d)));
                if (c.field->elevation(dd) < -0.5f) { best = d; S.found = true; S.up = prev; S.sea = dd; break; }
                prev = dd;
            }
        }
        if (!c.check("a sea shore within 4 km of the start", S.found)) c.fail("no shore");
        return true;
    });
    s.step("warm-up", kBudgetSettleCap + 5.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        Walk(c, 0.0);
        c.waitingFor("the world to settle");
        static BudgetSettle settle;
        return settle.done(c);
    });
    s.step("measure", kBudgetMeasure + 10.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        if (c.t < 0.02) Profiler_Reset(); // the zone means below cover the timed window only
        Walk(c, c.t);
        S.frames.frame();
        c.waitingFor("the timed window");
        if (c.t <= kBudgetMeasure) return false;
        S.frames.report(c);
        // (Fluid::GetStats().stepMs is the local fluid sim only: 0 on a plain shore.) The sea's cost is its render pass.
        c.metric("water pass gpu mean", std::max(0.0, Profiler_ZoneMeanMs("gpu:water")), "ms");
        c.metric("water pass cpu mean", std::max(0.0, Profiler_ZoneMeanMs("cpu:water")), "ms");
        return true;
    });
}
