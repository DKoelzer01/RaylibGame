// caves/cavern-view (#1553): a look inside the nearest large cavern (kind 4) to the start, as caves/chamber-view does for
// chambers: lit by a light at the eye (RAYLIBGAME_CAMERA_LIGHT, 60 m), four looks from 3 m over its centre (the dome, the
// pillars and the breakdown floor) and two from 0.7 r out toward its wall (the wall at the user's rules: outcrops capped,
// nothing floating). Kept for a person to judge.
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>
#include <vector>

namespace {
struct State { bool found = false; Vector3 at{}, up{}, e1{}, e2{}; float r = 0.0f; } S;
}

SCENARIO("caves/cavern-view") {
    s.tags = { "caves" };
    s.world = ScenarioWorld::Planet;
    s.env = { { "RAYLIBGAME_CAMERA_LIGHT", "1" } };
    s.timeout = 200.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("the nearest cavern", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        std::vector<Vector3> at; std::vector<float> radii;
        c.field->caveMouthsList(home, 3000.0f, 4, at, &radii); // kind 4: a cavern's centre
        float best = 1e30f;
        for (size_t i = 0; i < at.size(); ++i) {
            const float d = Vector3Distance(Vector3Scale(Vector3Normalize(at[i]), c.field->radius), Vector3Scale(home, c.field->radius));
            if (d < best) { best = d; S.at = at[i]; S.r = radii[i]; S.found = true; }
        }
        if (!c.check("a cavern within 3 km of the start", S.found)) { c.fail("no cavern"); return true; }
        S.up = Vector3Normalize(S.at);
        S.e1 = Vector3Normalize(Vector3CrossProduct(S.up, fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        S.e2 = Vector3CrossProduct(S.up, S.e1);
        c.metric("cavern radius (m)", S.r, "m");
        c.metric("cavern depth under the surface (m)", c.field->radius + c.field->elevation(S.up) - Vector3Length(S.at), "m");
        c.metric("cavern distance from the start (m)", best, "m");
        return true;
    });
    s.step("built", 90.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        const Vector3 eye = Vector3Add(S.at, Vector3Scale(S.up, 3.0f));
        c.camera(eye, Vector3Add(eye, Vector3Scale(S.e1, 20.0f)), S.up);
        c.waitingFor("the cavern's terrain built (an idle world)");
        return c.t > 2.0 && WorldIdle();
    });
    for (int k = 0; k < 6; ++k) {
        s.step("look " + std::to_string(k), 8.0, [k](ScenarioCtx& c) {
            if (!S.found) return true;
            if (k == 0 && c.t < 0.05) c.cue();
            const bool wall = k >= 4; // the last two: 0.7 r out, facing the wall
            const float a = wall ? 3.1415927f * (k - 4) : 1.5707963f * k;
            const Vector3 dir = Vector3Add(Vector3Scale(S.e1, cosf(a)), Vector3Scale(S.e2, sinf(a)));
            const Vector3 eye = Vector3Add(Vector3Add(S.at, Vector3Scale(S.up, 3.0f)), Vector3Scale(dir, wall ? 0.7f * S.r : 0.0f));
            c.camera(eye, Vector3Add(eye, Vector3Add(Vector3Scale(dir, 20.0f), Vector3Scale(S.up, wall ? 6.0f : 4.0f))), S.up);
            if (wall && c.t < 3.0) return false; // (the terrain near the new eye builds)
            if (c.t > (wall ? 4.5 : 1.5)) c.grab("look" + std::to_string(k));
            return c.t > (wall ? 5.0 : 2.0) && (!wall || WorldIdle());
        });
    }
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("inside a large cavern (#1553): for a person to judge"); return true; });
}
