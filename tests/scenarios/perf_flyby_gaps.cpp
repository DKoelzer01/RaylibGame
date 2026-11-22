// perf/flyby-gaps (#1418, the #545 guard): flying fast over the planet with the frame rate capped (30 FPS, as a capped or
// vsynced player has it), the view's terrain keeps up: the GPU density queue stays short in steady flight (the check), and the
// desired LOD nodes neither built nor covered (WorldStats::lodUncovered) are reported. Before #1417 the queue passed 30 000
// blocks: the streaming budget took the cap's wait for work and throttled the GPU density batches (gaps in the mountains).
// The camera flies a great circle from the start 1300 m over sea level, 2 m a frame (60 m/s at 30 FPS), for 1200 frames;
// the gap count is averaged from frame 450 on (after the first 15 s, the streaming in after the jump).
#include "scenario.h"
#include "voxel.h"
#include "terraingpu.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace {
struct State { Vector3 home{}, axis{}; long frames = 0; double sum = 0.0; int samples = 0, worst = 0, gpuPeak = 0; } S;
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
    s.args = { "--fps", "30", "--lod", "2.5", "--preload", "300", "--cache-verts", "8" }; // the preload ring and cache as in play: their builds share the GPU queue
    s.record = false; // a recording fixes the frame time and its writes cost frame time: the cap must be the only limit
    s.timeout = 180.0;

    s.step("settle", 90.0, [](ScenarioCtx& c) {
        c.waitingFor("2 s of play and an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        S.home = c.field->findSpawnDir(c.body->seed);
        // Toward the highest ground within 8 km (seed 42: a 1650 m peak), so the flight crosses mountains, as the gaps were
        // seen: a heading out over the sea (the start is near it) builds little and never backs up.
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(S.home, fabsf(S.home.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 t2 = Vector3CrossProduct(S.home, t1);
        float bestE = -1e30f; Vector3 peak = Vector3Add(S.home, Vector3Scale(t1, 1e-3f));
        for (float a = -8000.0f; a <= 8000.0f; a += 200.0f)
            for (float b = -8000.0f; b <= 8000.0f; b += 200.0f) {
                if (a * a + b * b < 2000.0f * 2000.0f) continue; // at least 2 km out: a heading, not a nearby bump
                const Vector3 d = Vector3Normalize(Vector3Add(S.home, Vector3Add(Vector3Scale(t1, a / c.field->radius), Vector3Scale(t2, b / c.field->radius))));
                const float e = c.field->elevation(d);
                if (e > bestE) { bestE = e; peak = d; }
            }
        S.axis = Vector3Normalize(Vector3CrossProduct(S.home, peak));
        c.metric("highest ground on the heading (m)", bestE, "m");
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
            S.sum += u; ++S.samples; S.worst = std::max(S.worst, u); S.gpuPeak = std::max(S.gpuPeak, TerrainGpu_Pending());
        }
        if (++S.frames < kFlightFrames) return false;
        const double mean = S.samples ? S.sum / S.samples : 0.0;
        c.metric("uncovered desired nodes in steady flight (mean)", mean, "nodes");
        c.metric("uncovered desired nodes in steady flight (worst)", S.worst, "nodes");
        c.note(TextFormat("%d samples over %ld frames of flight, %.0f m flown; GPU density queue peak %d", S.samples, S.frames, S.frames * kMetresPerFrame, S.gpuPeak));
        c.metric("GPU density queue peak in steady flight", S.gpuPeak, "blocks");
        // The throttled streaming itself: with the cap's wait taken for work (before #1417) the queue grew past 30 000 blocks
        // in this flight (the view's terrain falling behind: gaps); with the fix it peaks under 1000.
        c.lt("GPU density queue peak in steady flight (blocks; > 30 000 before #1417)", S.gpuPeak, 5000.0);
        return true;
    });
}
