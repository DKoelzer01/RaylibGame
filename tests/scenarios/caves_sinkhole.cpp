// caves/sinkhole (#1181): flat-ground cave entrances. Hillside mouths need a real face, so a swamp or meadow start had
// none for kilometres (#1115: the nearest 1.1 km from seed 42's start at the best odds). A sinkhole - a collapse pit with
// a passage from its foot - must be found near the start, open down its shaft (the ground's density there is air), and
// show from above (the camera looks down into it; the capture is kept for a person to judge).
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct State {
    bool found = false;
    Vector3 at{}, out{}, up{}, side{};
} S;
}

SCENARIO("caves/sinkhole") {
    s.tags = { "caves", "world" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 150.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("find a sinkhole near the start", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        S.found = c.field->nearestCaveMouth(home, 3000.0f, S.at, S.out, 1); // kind 1: shafts (sinkholes, karst sinks, skylights)
        if (!c.check("a sinkhole or shaft within 3 km of the start", S.found)) { c.fail("none found"); return true; }
        const float dist = Vector3Distance(Vector3Scale(Vector3Normalize(S.at), c.field->radius), Vector3Scale(home, c.field->radius));
        c.metric("nearest shaft entrance from the start (m)", dist, "m");
        c.lt("nearest shaft entrance from the start (m)", dist, 3000.0);
        S.up = Vector3Normalize(S.at);
        S.side = Vector3Normalize(Vector3CrossProduct(S.up, fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        // Open down the shaft: air 3 m under the rim and 6 m down (every sinkhole is 8 m+ deep).
        const float d3 = c.field->density(Vector3Subtract(S.at, Vector3Scale(S.up, 3.0f)), nullptr);
        const float d6 = c.field->density(Vector3Subtract(S.at, Vector3Scale(S.up, 6.0f)), nullptr);
        c.lt("density 3 m down the shaft (air < iso)", d3, ISO_LEVEL);
        c.lt("density 6 m down the shaft (air < iso)", d6, ISO_LEVEL);
        return true;
    });
    s.step("look down into it", 40.0, [](ScenarioCtx& c) {
        if (!S.found) return true;
        c.camera(Vector3Add(S.at, Vector3Add(Vector3Scale(S.up, 14.0f), Vector3Scale(S.side, 8.0f))), Vector3Subtract(S.at, Vector3Scale(S.up, 6.0f)), S.up);
        c.waitingFor("the terrain there built (an idle world), then 2 s");
        static double idleAt = -1.0;
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) { idleAt = c.t; c.cue(); }
        return idleAt >= 0.0 && c.t > idleAt + 2.0;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) {
        c.keepCapture("a sinkhole seen from above (#1181): for a person to judge");
        return true;
    });
}
