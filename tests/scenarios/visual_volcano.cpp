// visual/volcano (#1591): the nearest volcano (as /locate volcano finds it) seen from above its crater rim, then from its
// flank: an ash cone with no snow, the crater floor's lava crust glowing, basalt flanks over the climate's snowy skirt.
// For a person to judge (a look change); run again with --sun-elev -30 to see the crust light the crater at night.
#include "scenario.h"
#include "locate.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>
#include <vector>

namespace {
Vector3 vent{}, side{};
bool found = false;
}

SCENARIO("visual/volcano") {
    s.tags = { "visual", "terrain" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day", "--record-hud", "0" };
    s.timeout = 240.0;
    s.step("find the nearest volcano's vent", 20.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        // Climb the volcano weight from /locate's summit to its vent (the crater's centre)
        std::vector<LocateHit> hits; std::string why;
        const Vector3 from = Vector3Scale(c.field->findSpawnDir(c.body->seed), c.field->radius + 2.0f);
        if (!c.check("/locate volcano finds one", Locate_Find(*c.body, from, "volcano", 1, hits, why), why)) return true;
        Vector3 d = Vector3Normalize(hits[0].at);
        float w = c.field->volcanoWeight(Vector3Scale(d, c.field->radius));
        for (int it = 0; it < 400; ++it) { // a greedy climb, 20 m steps
            bool moved = false;
            for (int k = 0; k < 8; ++k) {
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(d, { 0, 1, 0 })), t2 = Vector3CrossProduct(d, t1);
                const float a = 0.785398f * k;
                const Vector3 q = Vector3Normalize(Vector3Add(d, Vector3Scale(Vector3Add(Vector3Scale(t1, cosf(a)), Vector3Scale(t2, sinf(a))), 20.0f / c.field->radius)));
                const float qw = c.field->volcanoWeight(Vector3Scale(q, c.field->radius));
                if (qw > w) { w = qw; d = q; moved = true; }
            }
            if (!moved) break;
        }
        vent = d; found = true;
        side = Vector3Normalize(Vector3CrossProduct(vent, { 0, 1, 0 }));
        c.note(TextFormat("vent weight %.3f, %.0f m up", w, c.field->elevation(vent)));
        return c.check("the vent found (weight near 1)", w > 0.9f);
    });
    s.step("over the crater, looking down into it", 30.0, [](ScenarioCtx& c) {
        if (!found) return true;
        const float R = c.field->radius, e = c.field->elevation(vent);
        const Vector3 at = Vector3Scale(vent, R + e);
        const Vector3 eye = Vector3Add(Vector3Scale(vent, R + e + 260.0f), Vector3Scale(side, 420.0f));
        c.camera(eye, at, vent);
        c.waitingFor("the crater to stream in");
        if (c.t < 6.0 || !WorldIdle()) return false;
        c.shot("crater", 100.0);
        return true;
    });
    s.step("from the flank: the bare cone over the snowy skirt", 30.0, [](ScenarioCtx& c) {
        if (!found) return true;
        const float R = c.field->radius;
        const Vector3 foot = Vector3Normalize(Vector3Add(vent, Vector3Scale(side, 9000.0f / R)));
        const Vector3 eye = Vector3Scale(foot, R + c.field->elevation(foot) + 600.0f);
        c.camera(eye, Vector3Scale(vent, R + c.field->elevation(vent) - 300.0f), foot);
        c.waitingFor("the cone to stream in");
        if (c.t < 6.0 || !WorldIdle()) return false;
        c.shot("flank", 100.0);
        return true;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) { c.keepCapture("a volcano's look (#1591): for a person to judge"); return true; });
}
