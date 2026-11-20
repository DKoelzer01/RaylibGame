// fluid/flow (#870, moved from --test-flow, Plan 4 F2/F3): cut a channel out of the nearest inland water, ending in a shaft
// down to a sealed chamber. The water must reach the channel's end, flood the chamber, and sleep once settled, with the step
// under budget: the gate is the worst 1 s rolling mean of the fluid step time (< 30 ms, #627). PASS as the flag.
#include "scenario.h"
#include "fluid.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
struct State {
    Vector3 dir{}, out{}, end{}, cave{};
    float level = 0.0f, shore = 0.0f;
    double reachT = -1.0, caveT = -1.0, sleepT = -1.0;
    std::vector<float> steps; std::vector<double> stepT; long long lastStep = -1;
} S;
}

SCENARIO("fluid/flow") {
    s.tags = { "fluid", "water", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 180.0;
    s.legacyFlag = "--seed 42 --test-flow";

    s.step("dig", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (!WorldIdle()) return false;
        const PlanetField* pf = c.field;
        auto at = [&](Vector3 d, float h) { return Vector3Scale(d, pf->radius + h); };
        const Vector3 home = pf->findSpawnDir(c.body->seed);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(home, { 0, 1, 0 })), t2 = Vector3CrossProduct(t1, home);
        float best = 1e30f;
        // Nearest inland water 1.5 m+ deep, flat and wide (#723), with a low shore within 25 m: out to 3 km, then 6 and 12.
        for (float span = 3000.0f; span <= 12000.0f && best > 1e29f; span *= 2.0f)
            for (float e = -span; e <= span; e += 16.0f * span / 3000.0f)
                for (float n = -span; n <= span; n += 16.0f * span / 3000.0f) {
                    const float dist = sqrtf(e * e + n * n);
                    if (dist >= best) continue;
                    const Vector3 d = Vector3Normalize(Vector3Add(home, Vector3Add(Vector3Scale(t1, e / pf->radius), Vector3Scale(t2, n / pf->radius))));
                    const PlanetField::WaterSource ws = pf->waterSourceAt(d);
                    if (!ws.kind || ws.kind == 3 || ws.level <= 2.0f || pf->elevation(d) > ws.level - 1.5f) continue;
                    int wide = 0; bool flat = true;
                    for (const Vector3& o : { t1, Vector3Negate(t1), t2, Vector3Negate(t2) }) {
                        const PlanetField::WaterSource wn = pf->waterSourceAt(Vector3Normalize(Vector3Add(d, Vector3Scale(o, 6.0f / pf->radius))));
                        if (!wn.kind) continue;
                        ++wide; flat = flat && fabsf(wn.level - ws.level) < 0.3f;
                    }
                    if (wide < 2 || !flat) continue;
                    const Vector3 outs[4] = { t1, Vector3Negate(t1), t2, Vector3Negate(t2) };
                    for (const Vector3& o : outs)
                        for (float sh = 2.0f; sh <= 25.0f; sh += 1.0f) {
                            const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(o, sh / pf->radius)));
                            const float eq = pf->elevation(q);
                            if (pf->waterSourceAt(q).kind) continue;
                            if (eq > ws.level + 0.3f && eq < ws.level + 4.0f) { best = dist; S.dir = d; S.out = o; S.level = ws.level; S.shore = sh; }
                            break;
                        }
                }
        if (best > 1e29f) { c.fail("no water with a low shore within 12 km of the start"); return false; }
        c.note("water " + std::to_string((int)best) + " m from the start, level " + std::to_string(S.level) + " m, shore " + std::to_string((int)S.shore) + " m out");
        std::vector<TerrainEdit> cut;
        for (float sh = S.shore - 2.0f; sh <= S.shore + 14.0f; sh += 0.8f) {
            const Vector3 q = Vector3Normalize(Vector3Add(S.dir, Vector3Scale(S.out, sh / pf->radius)));
            for (float h = S.level - 0.3f; h <= S.level + 4.0f; h += 1.2f) cut.push_back({ at(q, h), 1.6f, -1.0f, true });
        }
        const Vector3 endDir = Vector3Normalize(Vector3Add(S.dir, Vector3Scale(S.out, (S.shore + 14.0f) / pf->radius)));
        // A 6 m shaft to a 2 m chamber (#935 / #1032, user 2026-11-19: shrink the scene, keep ~2 min): ~200 m^3 of shaft and
        // chamber settled only after ~220-230 s. The same scene as the legacy --test-flow (main.cpp kFlowShaftM...).
        constexpr float kShaftM = 6.0f, kChamberR = 2.0f, kCaveM = 8.0f;
        for (float h = S.level - 1.0f; h >= S.level - kShaftM; h -= 0.8f) cut.push_back({ at(endDir, h), 1.3f, -1.0f, true });
        cut.push_back({ at(endDir, S.level - kShaftM - 1.0f), kChamberR, -1.0f, true });
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        S.end = at(Vector3Normalize(Vector3Add(S.dir, Vector3Scale(S.out, (S.shore + 10.0f) / pf->radius))), S.level - 0.6f);
        S.cave = at(endDir, S.level - kCaveM);
        return true;
    });
    s.step("flow", 125.0, [](ScenarioCtx& c) {
        const Vector3 up = Vector3Normalize(S.end);
        c.camera(Vector3Add(Vector3Add(S.end, Vector3Scale(up, 9.0f)), Vector3Scale(S.out, 9.0f)), Vector3Subtract(S.end, Vector3Scale(S.out, 6.0f)), up);
        const double t = c.t;
        if (S.reachT < 0.0 && Fluid::Wet(*c.body, S.end)) S.reachT = t;
        if (S.caveT < 0.0 && Fluid::Wet(*c.body, S.cave)) S.caveT = t;
        const Fluid::Stats fs = Fluid::GetStats();
        if (fs.steps != S.lastStep) { S.lastStep = fs.steps; if (fs.steps > 0 && S.sleepT < 0.0) { S.steps.push_back((float)fs.stepMs); S.stepT.push_back(t); } }
        if (S.reachT >= 0.0 && S.caveT >= 0.0 && S.sleepT < 0.0 && Fluid::Asleep()) S.sleepT = t;
        c.waitingFor(S.reachT < 0.0 ? "the channel's end wet" : S.caveT < 0.0 ? "the chamber wet" : "the water asleep");
        if (!((S.sleepT >= 0.0 && t > S.sleepT + 1.0) || t > 120.0)) return false;
        float roll = 0.0f; // the worst 1 s rolling mean of the step time (#627)
        for (size_t a = 0, b = 0; b < S.stepT.size(); ++b) {
            while (S.stepT[b] - S.stepT[a] > 1.0) ++a;
            if (S.stepT[b] - S.stepT[0] < 1.0) continue;
            float sum = 0.0f; for (size_t k = a; k <= b; ++k) sum += S.steps[k];
            roll = std::max(roll, sum / (float)(b - a + 1));
        }
        c.metric("reached the channel end", S.reachT, "s");
        c.metric("reached the chamber", S.caveT, "s");
        c.metric("asleep after", S.sleepT, "s");
        c.metric("fluid step, worst 1 s mean", roll, "ms");
        c.metric("fluid step, max", fs.maxStepMs, "ms");
        c.metric("water over the chamber floor", Fluid::SurfaceAbove(*c.body, S.cave, 16.0f), "m", false);
        c.check("water reached the channel's end", S.reachT >= 0.0);
        c.check("water reached the chamber", S.caveT >= 0.0);
        c.check("the water went to sleep", S.sleepT >= 0.0);
        c.timingLt("fluid step, worst 1 s rolling mean (ms, #627)", roll, 30.0); // timing: judged on the quiet PC only (#1034)
        return true;
    });
}
