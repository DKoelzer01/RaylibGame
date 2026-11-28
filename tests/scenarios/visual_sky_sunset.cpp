// visual/sky-sunset (#2143): the physical sky's LUTs (sky view + haze, ~0.5 ms of GPU a rebuild) follow the moving sun in
// steps of 4e-4 in cos(zenith) instead of 1e-4, which the turning planet's sun crossed nearly every frame. At sunset, where
// the sky's colour changes fastest, the steps must not show. Sandbox, no --day (the sun moves): the world clock is stepped
// until the sun is 1-5 degrees up, then 6 s of real-rate time looking at the horizon under the sun (recorded from the cue).
// Checks: the LUTs rebuilt on under an eighth of the frames (before: a quarter). The look: record with
// RAYLIBGAME_SKY_SUN_STEP=1e-4 (before) and compare (capture skill -Compare).
#include "scenario.h"
#include "player.h"
#include "skyphys.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
struct St { Vector3 feet{}, up{}; long rebuilds0 = 0; int frames = 0; } S;

float SunUp(ScenarioCtx& c) { return Vector3DotProduct(Vector3Normalize(c.world.sunDirection), c.body->dirToWorld(S.up)); }

void Look(ScenarioCtx& c) { // level, toward the sun's bearing, a little up: the horizon glow fills the upper half
    const Vector3 sunL = c.body->dirToLocal(Vector3Normalize(c.world.sunDirection));
    Vector3 flat = Vector3Subtract(sunL, Vector3Scale(S.up, Vector3DotProduct(sunL, S.up)));
    flat = Vector3Length(flat) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(S.up, { 0, 0, 1 })) : Vector3Normalize(flat);
    const Vector3 eye = Vector3Add(S.feet, Vector3Scale(S.up, 4.0f));
    c.camera(eye, Vector3Add(eye, Vector3Add(Vector3Scale(flat, 20.0f), Vector3Scale(S.up, 3.0f))), S.up);
}
}

SCENARIO("visual/sky-sunset") {
    s.tags = { "visual", "sky" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--no-hud" };
    s.timeout = 90.0;

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        S.feet = c.body->toLocal((player.position - c.world.renderOrigin).toFloat());
        S.up = Vector3Normalize(S.feet);
        return true;
    });
    s.step("the clock to sunset: the sun 1-5 degrees up and going down", 30.0, [](ScenarioCtx& c) {
        Look(c);
        static float last = 2.0f;
        const float up = SunUp(c);
        c.waitingFor(TextFormat("sun up %.3f", up));
        if (up > 0.017f && up < 0.087f && up < last) return true;
        last = up;
        SetWorldTime(WorldTime() + (up > 0.3f || up < 0.0f ? 60.0 : 4.0)); // 6 or 0.4 degrees a step
        return false;
    });
    s.step("a few frames for the view", 2.0, [](ScenarioCtx& c) {
        Look(c);
        if (c.t < 0.5) return false;
        S.rebuilds0 = SkyPhys_Rebuilds(); S.frames = 0;
        c.cue();
        return true;
    });
    s.step("6 s of sunset", 10.0, [](ScenarioCtx& c) {
        Look(c);
        ++S.frames;
        if (c.t < 6.0) return false;
        const long n = SkyPhys_Rebuilds() - S.rebuilds0;
        c.metric("sky LUT rebuilds", (double)n, "");
        c.metric("frames", S.frames, "", false);
        c.gt("the LUTs followed the sun (rebuilt at all)", (double)n, 0.5);
        c.lt("rebuilt on under an eighth of the frames", (double)n, S.frames / 8.0); // 23 of 361; at 1e-4 (before): 90
        c.shot("sunset", 6.0);
        return true;
    });
}
