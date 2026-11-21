// perf/flyby-gaps (#1418, the #545 guard): flying fast over the planet with the frame rate capped (30 FPS, as a capped or
// vsynced player has it), the view's terrain keeps up: the desired LOD nodes that are neither built nor covered by a kept
// node (WorldStats::lodUncovered: gaps on screen, sky through the mountains) stay near zero in steady flight. Before #1417
// the streaming budget took the cap's wait for work and throttled the GPU density batches: ~2000 nodes uncovered.
// The camera flies a great circle from the start 1300 m over sea level, 2 m a frame (60 m/s at 30 FPS), for 1200 frames;
// the gap count is averaged from frame 450 on (after the first 15 s, the streaming in after the jump).
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace {
struct State { Vector3 home{}, axis{}; long frames = 0; double sum = 0.0; int samples = 0, worst = 0; } S;
constexpr long kFlightFrames = 1200, kSteadyFrom = 450;
constexpr float kMetresPerFrame = 2.0f, kHeight = 1300.0f;

void Fly(ScenarioCtx& c, long f) {
    const float R = c.field->radius;
    const float a = (float)f * kMetresPerFrame / (R + kHeight);
    const Vector3 dir = Vector3RotateByAxisAngle(S.home, S.axis, a);
    const Vector3 ahead = Vector3RotateByAxisAngle(S.home, S.axis, a + 400.0f / R); // looking 400 m ahead and down at the ground
    const Vector3 eye = Vector3Scale(dir, R + kHeight);
    const Vector3 target = Vector3Scale(ahead, R + c.field->elevation(ahead));
    c.camera(eye, target, dir);
}
}

SCENARIO("perf/flyby-gaps") {
    s.tags = { "perf", "terrain" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--fps", "30", "--lod", "2.5" };
    s.record = false; // a recording fixes the frame time and its writes cost frame time: the cap must be the only limit
    s.timeout = 180.0;

    s.step("settle", 90.0, [](ScenarioCtx& c) {
        c.waitingFor("2 s of play and an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        S.home = c.field->findSpawnDir(c.body->seed);
        S.axis = Vector3Normalize(Vector3CrossProduct(S.home, fabsf(S.home.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        return true;
    });
    s.step("the first view loaded", 60.0, [](ScenarioCtx& c) {
        Fly(c, 0);
        c.waitingFor("the view 1300 m up built (an idle world)");
        return c.t > 1.0 && WorldIdle();
    });
    s.step("fly 1200 frames at 2 m a frame", 120.0, [](ScenarioCtx& c) {
        Fly(c, S.frames);
        if (S.frames >= kSteadyFrom) {
            const int u = GetWorldStats().lodUncovered;
            S.sum += u; ++S.samples; S.worst = std::max(S.worst, u);
        }
        if (++S.frames < kFlightFrames) return false;
        const double mean = S.samples ? S.sum / S.samples : 0.0;
        c.metric("uncovered desired nodes in steady flight (mean)", mean, "nodes");
        c.metric("uncovered desired nodes in steady flight (worst)", S.worst, "nodes");
        c.note(TextFormat("%d samples over %ld frames of flight, %.0f m flown", S.samples, S.frames, S.frames * kMetresPerFrame));
        c.lt("uncovered desired nodes in steady flight, mean (gaps; ~2000 before #1417)", mean, 60.0);
        return true;
    });
}
