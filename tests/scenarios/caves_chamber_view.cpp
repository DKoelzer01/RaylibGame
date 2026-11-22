// caves/chamber-view (#1116): a look at a cave chamber's walls - the eye 3 m over the nearest chamber's centre to the start,
// lit by a light at the eye (RAYLIBGAME_CAMERA_LIGHT, 60 m), turned to four sides in turn once the terrain there is built.
// Kept for a person to judge (floating rock on the upper walls before #1116; the hall's character after).
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
struct State { bool found = false; Vector3 at{}, up{}, e1{}, e2{}; double idleAt = -1.0; } S;
}

SCENARIO("caves/chamber-view") {
    s.tags = { "caves" };
    s.world = ScenarioWorld::Planet;
    s.env = { { "RAYLIBGAME_CAMERA_LIGHT", "1" } };
    s.timeout = 180.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("the nearest chamber", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        Vector3 out{};
        S.found = c.field->nearestCaveMouth(home, 3000.0f, S.at, out, 3); // kind 3: a chamber's centre
        if (!c.check("a chamber within 3 km of the start", S.found)) { c.fail("no chamber"); return true; }
        S.up = Vector3Normalize(S.at);
        S.e1 = Vector3Normalize(Vector3CrossProduct(S.up, fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        S.e2 = Vector3CrossProduct(S.up, S.e1);
        c.metric("chamber depth under the surface (m)", c.field->radius + c.field->elevation(S.up) - Vector3Length(S.at), "m");
        return true;
    });
    s.step("built", 60.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        const Vector3 eye = Vector3Add(S.at, Vector3Scale(S.up, 3.0f));
        c.camera(eye, Vector3Add(eye, Vector3Scale(S.e1, 20.0f)), S.up);
        c.waitingFor("the chamber's terrain built (an idle world)");
        return c.t > 2.0 && WorldIdle();
    });
    for (int k = 0; k < 4; ++k) {
        s.step("look " + std::to_string(k), 4.0, [k](ScenarioCtx& c) {
            if (!S.found) return true;
            if (k == 0 && c.t < 0.05) c.cue();
            const float a = 1.5707963f * k;
            const Vector3 dir = Vector3Add(Vector3Scale(S.e1, cosf(a)), Vector3Scale(S.e2, sinf(a)));
            const Vector3 eye = Vector3Add(S.at, Vector3Scale(S.up, 3.0f));
            c.camera(eye, Vector3Add(eye, Vector3Add(Vector3Scale(dir, 20.0f), Vector3Scale(S.up, 4.0f))), S.up); // a little up: the upper walls
            if (c.t > 1.5) c.grab("look" + std::to_string(k)); // (#1071) a slideshow picture, when asked for
            return c.t > 2.0;
        });
    }
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("a cave chamber's walls (#1116): for a person to judge"); return true; });
}
