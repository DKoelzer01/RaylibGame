// water/beach (#1164): the look of an undisturbed sea shore - the depth colour and the shore foam on generated water. The
// nearest shore to the start: walked out from it in 16 directions in 10 m steps until the ground drops under sea level;
// the eye 1.7 m over the beach 8 m inland of the waterline, looking out along it and a little toward the sea. Kept for a
// person to judge.
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
struct State { bool found = false; Vector3 land{}, sea{}, up{}; } S;
}

SCENARIO("water/beach") {
    s.tags = { "water", "light" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 150.0;

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
                if (c.field->elevation(dd) < -0.5f) {                       // the sea (under the sea level)
                    best = d; S.found = true;
                    S.up = prev;
                    const Vector3 out = Vector3Normalize(Vector3Subtract(dd, prev));
                    S.land = Vector3Scale(Vector3Normalize(Vector3Subtract(Vector3Scale(prev, R), Vector3Scale(out, 8.0f * R / R))), 1.0f);
                    S.sea = dd;
                    break;
                }
                prev = dd;
            }
        }
        if (!c.check("a sea shore within 4 km of the start", S.found)) { c.fail("no shore"); return true; }
        c.metric("nearest shore from the start (m)", best, "m");
        return true;
    });
    s.step("look along the waterline", 40.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        const float R = c.field->radius;
        const Vector3 up = Vector3Normalize(S.up);
        const Vector3 out = Vector3Normalize(Vector3Subtract(S.sea, S.up));            // toward the sea, along the ground
        const Vector3 along = Vector3Normalize(Vector3CrossProduct(up, out));
        const Vector3 beach = Vector3Normalize(Vector3Subtract(Vector3Scale(S.up, R), Vector3Scale(out, 8.0f)));
        const float e = std::max(0.0f, c.field->elevation(beach));
        const Vector3 eye = Vector3Scale(beach, R + e + 1.7f);
        const Vector3 target = Vector3Scale(Vector3Normalize(Vector3Add(Vector3Scale(S.up, R), Vector3Add(Vector3Scale(along, 25.0f), Vector3Scale(out, 10.0f)))), R);
        c.camera(eye, target, up);
        c.waitingFor("the shore built (an idle world), then 3 s");
        static double idleAt = -1.0;
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) { idleAt = c.t; c.cue(); }
        return idleAt >= 0.0 && c.t > idleAt + 3.0;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("a sea shore's look (#1164): for a person to judge"); return true; });
}
