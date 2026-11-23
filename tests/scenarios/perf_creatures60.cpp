// perf/creatures-60 (#1690): a budget scene. Sandbox; 60 mixed creatures (five species, hostile ones with threats on) in a
// ring 12-45 m around the player, the camera circling them 6 m up. 5 s of warm-up not measured, then 20 s timed: frame
// mean and p95, and the animal update's mean (the scene's own cost). Limits live in tools/ci/perf_budgets.txt (runner only).
#include "perf_budget.h"
#include "animals.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
struct St { BudgetFrames frames; double sum = 0.0, worst = 0.0; int n = 0, count0 = 0; Vector3 feet{}, up{}, e{}, nn{}; } S;
void Orbit(ScenarioCtx& c, double t) {
    const float a = (float)(t * 0.2), r = 14.0f;
    const Vector3 eye = Vector3Add(S.feet, Vector3Add(Vector3Scale(S.up, 6.0f), Vector3Add(Vector3Scale(S.e, r * cosf(a)), Vector3Scale(S.nn, r * sinf(a)))));
    c.camera(eye, Vector3Add(S.feet, Vector3Scale(S.up, 1.0f)), S.up);
}
}

SCENARIO("perf/creatures-60") {
    s.tags = { "perf", "budget" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    s.record = false;

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false);
        Animals_SetThreats(true);
        worldDifficulty = 1;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("60 mixed creatures around", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        S.feet = c.body->toLocal((player.position - c.world.renderOrigin).toFloat());
        S.up = Vector3Normalize(S.feet);
        S.e = Vector3Normalize(Vector3CrossProduct(S.up, fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        S.nn = Vector3CrossProduct(S.e, S.up);
        static const char* kSpecies[5] = { "woolback", "plains_runner", "clucker", "hornhide", "thornback" };
        S.count0 = Animals_Count();
        for (int k = 0; k < 60; ++k) {
            const float a = 2.0f * PI * k / 60.0f, r = 12.0f + 33.0f * (k % 6) / 5.0f;
            const Vector3 dir = Vector3Normalize(Vector3Add(S.feet, Vector3Add(Vector3Scale(S.e, r * cosf(a)), Vector3Scale(S.nn, r * sinf(a)))));
            Animals_Spawn(c.world, *c.body, kSpecies[k % 5], Vector3Scale(dir, c.field->radius + c.field->elevation(dir)), 1, false, 0.01f);
        }
        c.gt("creatures spawned", Animals_Count() - S.count0, 50.0);
        return true;
    });
    s.step("warm-up", kBudgetSettleCap + 5.0, [](ScenarioCtx& c) {
        survival.health = 100.0f; survival.dead = false;
        Orbit(c, c.t);
        c.waitingFor("the world to settle");
        static BudgetSettle settle;
        return settle.done(c);
    });
    s.step("measure", kBudgetMeasure + 10.0, [](ScenarioCtx& c) {
        survival.health = 100.0f; survival.dead = false;
        Orbit(c, c.t + kBudgetWarmup);
        S.frames.frame();
        const double ms = Animals_LastUpdateMs(); S.sum += ms; S.worst = std::max(S.worst, ms); ++S.n;
        c.waitingFor("the timed window");
        if (c.t <= kBudgetMeasure) return false;
        S.frames.report(c);
        c.metric("animal update mean", S.n ? S.sum / S.n : 0.0, "ms");
        c.metric("animal update worst", S.worst, "ms");
        c.metric("simulated animals", Animals_SimulatedCount(), "", false);
        return true;
    });
}
