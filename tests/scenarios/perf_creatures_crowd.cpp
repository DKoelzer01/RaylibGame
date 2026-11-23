// perf/creatures-crowd (#1649): a crowd chasing the player, the creature update's CPU cost per step (the user's world_30
// had dozens of cave creatures out on the meadow, #1633: the step's p99 was 19-23 ms). The planet; 64 creatures in a ring
// 8-50 m out (24 lost miners, 24 tunnel lurkers, 16 thornbacks; threats on, all provoked every 3 s so they chase; the
// player's health held up), 10 s of steps;
// mean, p95, p99 and worst reported. A timing run: a quiet window (perf scenarios stay off the shared runner).
#include "scenario.h"
#include "animals.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <vector>

namespace {
struct St { std::vector<double> ms; int count0 = 0; } S;
double Pct(std::vector<double> v, double p) { if (v.empty()) return 0.0; std::sort(v.begin(), v.end()); return v[std::min(v.size() - 1, (size_t)(p * v.size()))]; }
}

SCENARIO("perf/creatures-crowd") {
    s.tags = { "perf", "creatures" };
    s.world = ScenarioWorld::Planet; // real terrain: the ground probes cost what they do in play
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
    s.step("64 creatures around", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        const Vector3 feet = c.body->toLocal((player.position - c.world.renderOrigin).toFloat());
        const Vector3 up = Vector3Normalize(feet);
        const Vector3 e = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 n = Vector3CrossProduct(e, up);
        S.count0 = Animals_Count();
        for (int k = 0; k < 64; ++k) {
            const float a = 2.0f * PI * k / 64.0f + 0.37f * (k % 3), r = 8.0f + 42.0f * (k % 8) / 7.0f;
            const Vector3 dir = Vector3Normalize(Vector3Add(feet, Vector3Add(Vector3Scale(e, r * cosf(a)), Vector3Scale(n, r * sinf(a)))));
            const char* id = k % 8 < 3 ? "lost_miner" : k % 8 < 6 ? "tunnel_lurker" : "thornback";
            Animals_Spawn(c.world, *c.body, id, Vector3Scale(dir, c.field->radius + c.field->elevation(dir)), 1, false, 0.01f);
        }
        c.eq("64 creatures", Animals_Count() - S.count0, 64);
        return true;
    });
    s.step("10 s of steps", 25.0, [](ScenarioCtx& c) {
        survival.health = 100.0f; survival.dead = false;
        static double provokedAt = -10.0;
        if (c.t < 0.1) provokedAt = -10.0;
        if (c.t - provokedAt > 3.0) { Animals_TestProvokeAll(); provokedAt = c.t; } // all of them on the player, chasing
        if (c.t > 1.0) S.ms.push_back(Animals_LastUpdateMs());
        if (c.t < 11.0) return false;
        double sum = 0.0; for (double v : S.ms) sum += v;
        c.metric("creature update, mean (ms)", S.ms.empty() ? 0.0 : sum / S.ms.size(), "ms");
        c.metric("creature update, p95 (ms)", Pct(S.ms, 0.95), "ms");
        c.metric("creature update, p99 (ms)", Pct(S.ms, 0.99), "ms");
        c.metric("creature update, worst (ms)", Pct(S.ms, 1.0), "ms");
        c.metric("creatures alive", Animals_Count() - S.count0, "", false);
        c.metric("simulated", Animals_SimulatedCount(), "", false);
        return true;
    });
}
