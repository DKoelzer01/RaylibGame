// perf/creatures-20 (#1223, Plan 12 testing): 20 creatures on the brain active around a player, the animal update's CPU
// cost per step. Sandbox; 20 thornbacks in a ring 15-40 m out (threats on: they perceive, charge, wind up and ram; the
// player's health is held up). 10 s of steps; the update's mean and worst are reported, with a soft gate on the mean.
// A timing run: a quiet window on the PC (perf scenarios stay off the shared runner).
#include "scenario.h"
#include "animals.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
struct St { double sum = 0.0, worst = 0.0; int n = 0; int count0 = 0; } S;
}

SCENARIO("perf/creatures-20") {
    s.tags = { "perf", "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
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
    s.step("20 thornbacks around", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        const Vector3 feet = c.body->toLocal((player.position - c.world.renderOrigin).toFloat());
        const Vector3 up = Vector3Normalize(feet);
        const Vector3 e = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 n = Vector3CrossProduct(e, up);
        S.count0 = Animals_Count();
        for (int k = 0; k < 20; ++k) {
            const float a = 2.0f * PI * k / 20.0f, r = 15.0f + 25.0f * (k % 5) / 4.0f;
            const Vector3 dir = Vector3Normalize(Vector3Add(feet, Vector3Add(Vector3Scale(e, r * cosf(a)), Vector3Scale(n, r * sinf(a)))));
            Animals_Spawn(c.world, *c.body, "thornback", Vector3Scale(dir, c.field->radius + c.field->elevation(dir)), 1, false, 0.01f);
        }
        c.eq("20 creatures", Animals_Count() - S.count0, 20);
        return true;
    });
    s.step("10 s of steps", 20.0, [](ScenarioCtx& c) {
        survival.health = 100.0f; survival.dead = false;
        if (c.t > 1.0) { const double ms = Animals_LastUpdateMs(); S.sum += ms; S.worst = std::max(S.worst, ms); ++S.n; }
        if (c.t < 11.0) return false;
        const double mean = S.n ? S.sum / S.n : 0.0;
        c.metric("animal update, mean (ms)", mean, "ms");
        c.metric("animal update, worst (ms)", S.worst, "ms");
        c.metric("simulated animals", Animals_SimulatedCount(), "", false);
        c.timingLt("animal update mean under 1 ms with 20 creatures", mean, 1.0);
        return true;
    });
}
