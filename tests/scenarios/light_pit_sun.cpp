// light/pit-sun (#898, #1067): a pit dug a few metres down with the sun nearly overhead must be lit like the ground
// around it - by the sun where it shines in and by the sky seen through the opening - and must not go dark with depth
// alone (user 2026-11-20: "A pit dug a few meters down quickly becomes dark even if the sun is overhead"). A 6 x 6 x 6 m box
// is carved at the start with the sun 80 degrees up; the recording holds a view from the rim and one from the floor, then
// one from inside the nearest cave mouth (a closed passage must stay dark: no light leaking in). Kept as evidence: the
// look is for a person to judge; the checks only gate that the views happened.
#include "scenario.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f;
    bool mouth = false;
    Vector3 mouthAt{}, mouthOut{};
    Vector3 pt(float x, float y, float h) const {
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    TerrainEdit box(float x0, float x1, float y0, float y1, float h0b, float h1b) const {
        TerrainEdit e{};
        e.center = pt(0.5f * (x0 + x1), 0.5f * (y0 + y1), 0.5f * (h0b + h1b));
        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 0.5f * (y1 - y0), 0.5f * (h1b - h0b) };
        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        return e;
    }
} S;
}

SCENARIO("light/pit-sun") {
    s.tags = { "light", "world" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--sun-elev", "80" };
    s.timeout = 150.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("dig a 6 m pit", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        // The big pit, and a narrow one 12 m off: 2 x 2 m, 5 m down (the hand-dug case the user saw go dark).
        std::vector<TerrainEdit> cut = { S.box(-3.0f, 3.0f, -3.0f, 3.0f, -6.0f, 3.0f), S.box(11.0f, 13.0f, -1.0f, 1.0f, -5.0f, 3.0f) };
        // The user's shaft (2026-11-20 recording 204045): dug by hand, ~1.5 m across and 10 m down - tool strokes (carved
        // spheres) stacked every 0.5 m, 12 m the other way from the big pit.
        for (float d = 0.0f; d <= 10.0f; d += 0.5f) {
            TerrainEdit e{};
            e.center = S.pt(-12.0f, 0.0f, -d);
            e.radius = 0.8f; e.carve = true;
            cut.push_back(e);
        }
        ApplyTerrainEdits(*c.body, cut.data(), cut.size());
        c.camera(S.pt(6.0f, -6.0f, 5.0f), S.pt(0.0f, 0.0f, -5.5f), S.home);
        c.note("a 6 x 6 x 6 m box carved at the start, the sun 80 degrees up");
        return true;
    });
    s.step("remesh", 30.0, [](ScenarioCtx& c) {
        c.camera(S.pt(6.0f, -6.0f, 5.0f), S.pt(0.0f, 0.0f, -5.5f), S.home);
        c.waitingFor("the pit meshed (an idle world)");
        return c.t > 1.0 && WorldIdle();
    });
    s.step("from the rim", 3.0, [](ScenarioCtx& c) {
        static bool cued = false;
        if (!cued) { cued = true; c.cue(); }
        c.camera(S.pt(6.0f, -6.0f, 5.0f), S.pt(0.0f, 0.0f, -5.5f), S.home);
        return c.t > 2.0;
    });
    s.step("from the floor", 3.0, [](ScenarioCtx& c) {
        c.camera(S.pt(-1.5f, 2.0f, -4.4f), S.pt(1.5f, -3.0f, -5.0f), S.home); // eye 1.6 m over the floor, toward a wall's foot
        return c.t > 2.0;
    });
    s.step("narrow pit from the rim", 3.0, [](ScenarioCtx& c) {
        c.camera(S.pt(12.0f, -2.2f, 2.5f), S.pt(12.0f, 0.3f, -5.0f), S.home);  // looking down into the 2 x 2 m hole
        return c.t > 2.0;
    });
    s.step("narrow pit from its floor", 3.0, [](ScenarioCtx& c) {
        c.camera(S.pt(11.6f, -0.6f, -3.6f), S.pt(12.6f, 0.9f, -4.8f), S.home); // eye 1.4 m over its floor
        return c.t > 2.0;
    });
    s.step("hand-dug shaft from its bottom, looking up", 3.0, [](ScenarioCtx& c) {
        c.camera(S.pt(-12.0f, 0.1f, -9.0f), S.pt(-11.7f, 0.4f, 0.0f), S.ax); // looking up the shaft: the screen's up across it
        return c.t > 2.0;
    });
    s.step("hand-dug shaft from the rim, looking down", 3.0, [](ScenarioCtx& c) {
        c.camera(S.pt(-12.0f, -2.0f, 1.8f), S.pt(-12.0f, 0.0f, -9.0f), S.home);
        return c.t > 2.0;
    });
    s.step("inside a cave mouth", 45.0, [](ScenarioCtx& c) {
        static bool looked = false;
        static double idleAt = -1.0;
        if (!looked) {
            looked = true;
            S.mouth = c.field->nearestCaveMouth(S.home, 4000.0f, S.mouthAt, S.mouthOut, 0);
            c.check("a hillside cave mouth within 4 km", S.mouth);
        }
        if (!S.mouth) return true;
        const Vector3 up = Vector3Normalize(S.mouthAt);
        const Vector3 in = Vector3Negate(S.mouthOut);
        c.camera(Vector3Add(S.mouthAt, Vector3Scale(in, 8.0f)), Vector3Add(S.mouthAt, Vector3Scale(in, 20.0f)), up);
        c.waitingFor("the terrain at the mouth built (an idle world), then 2 s");
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) idleAt = c.t;
        return idleAt >= 0.0 && c.t > idleAt + 2.0;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) {
        c.check("the views were held", true);
        c.keepCapture("pit and cave-mouth lighting with the sun overhead (#898/#1067): for a person to judge");
        return true;
    });
}
