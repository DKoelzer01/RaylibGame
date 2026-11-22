// slides/* (#1071): hand-placed views for the menu slideshow, on the real planet at seed 42 by day. Each scenario finds its
// subject (a river near the start, the nearest cave chamber, a lava tube's skylight on the nearest volcano flank, the planet,
// its ring and the moon from space), places the camera, waits for the terrain there to be built and takes its pictures
// with c.grab (only when RAYLIBGAME_SLIDE_SHOTS names a folder; tools/slides/shots_c.sh runs them). Tagged visual: they
// judge nothing and run on the PC only.
#include "scenario.h"
#include "terraintool.h"
#include "plots.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace {
Vector3 Tangent(Vector3 up) { return Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })); }
// The direction toward the sun in the planet's frame, through the same toLocal as everything else (two points' difference:
// dirToLocal with either sign framed the night side in the first reshoots).
Vector3 SunLocal(ScenarioCtx& c) {
    return Vector3Normalize(Vector3Subtract(c.body->toLocal(Vector3Scale(c.world.sunDirection, 1000.0f)), c.body->toLocal({ 0, 0, 0 })));
}

// One view per step: the camera is set every frame; the picture is taken once the world there is built (or at 30 s).
void View(ScenarioDef& s, const std::string& name, std::function<bool(ScenarioCtx&, Vector3&, Vector3&, Vector3&)> place) {
    s.step("view " + name, 45.0, [name, place](ScenarioCtx& c) {
        Vector3 eye{}, target{}, up{};
        if (!place(c, eye, target, up)) return true; // nothing to show (noted by the finder)
        c.camera(eye, target, up);
        c.waitingFor("the terrain at the view built (an idle world)");
        if (c.t < 4.0 || (!WorldIdle() && c.t < 30.0)) return false;
        c.grab(name);
        return c.t > 4.5;
    });
}

void Common(ScenarioDef& s) {
    s.tags = { "visual", "slides" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 300.0;
    s.record = false;
    s.step("settle", 60.0, [](ScenarioCtx& c) {
        // No tool aim in the pictures: with the eye near a wall the dig preview drew its red "can't dig" ring in the
        // middle of the cave chamber's slide (#1071). These views use no tool. (reach doesn't do it: c.camera flies, and
        // a flying aim reaches 150 m.)
        terrainTool.noAim = true;
        c.waitingFor("2 s of play and an idle world");
        return c.T >= 2.0 && WorldIdle();
    });
}

// The river: the nearest river cell to the start on a spiral of directions out to 8 km.
struct River { bool found = false; Vector3 at{}, up{}; } R;
// The chamber: its centre and floor; lamps placed around it so the hall can be seen.
struct Chamber { bool found = false; Vector3 floor{}, up{}, e1{}, e2{}; } C;
// The lava tube: a skylight (kind 1 mouth) near the nearest volcano flank.
struct Tube { bool found = false; Vector3 at{}, up{}, out{}; } T;
}

SCENARIO("slides/river-valley") {
    Common(s);
    s.step("the nearest river", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        const Vector3 k1 = Tangent(home), k2 = Vector3CrossProduct(home, k1);
        for (float r = 200.0f; r <= 8000.0f && !R.found; r += 100.0f)
            for (int a = 0; a < 64 && !R.found; ++a) {
                const float th = 6.2831853f * a / 64.0f;
                const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, c.field->radius), Vector3Add(Vector3Scale(k1, r * cosf(th)), Vector3Scale(k2, r * sinf(th)))));
                if (c.field->waterSourceAt(d).kind == 1) {
                    R.found = true; R.up = d; R.at = Vector3Scale(d, c.field->radius + c.field->elevation(d));
                    c.metric("river from the start (m)", r, "m");
                    // (#1438) where, for whoever looks at the river there: the body-local direction, lat/lon (y up) and height
                    c.note(TextFormat("river cell: dir (%.6f, %.6f, %.6f), lat %.5f lon %.5f, ground %.1f m over the radius, %.0f m from the start at %.0f deg",
                                      d.x, d.y, d.z, asinf(d.y) * RAD2DEG, atan2f(d.z, d.x) * RAD2DEG, c.field->elevation(d), r, th * RAD2DEG));
                }
            }
        c.check("a river within 8 km of the start", R.found);
        return true;
    });
    // From 60 m up and 220 m away, four sides: the valley along or across the river.
    for (int k = 0; k < 4; ++k)
        View(s, "side" + std::to_string(k), [k](ScenarioCtx& c, Vector3& eye, Vector3& target, Vector3& up) {
            if (!R.found) return false;
            const Vector3 e1 = Tangent(R.up), e2 = Vector3CrossProduct(R.up, e1);
            const float a = 1.5707963f * k;
            const Vector3 side = Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a)));
            const Vector3 ground = Vector3Add(R.at, Vector3Scale(side, 220.0f));
            const Vector3 gd = Vector3Normalize(ground);
            eye = Vector3Scale(gd, c.field->radius + std::max(c.field->elevation(gd), c.field->elevation(R.up)) + 60.0f);
            target = R.at; up = R.up;
            return true;
        });
}

SCENARIO("slides/cave-chamber") {
    Common(s);
    // The lamps alone left the hall black (batch C, first pass): a light at the eye as well (60 m).
    s.env = { { "RAYLIBGAME_CAMERA_LIGHT", "1" } };
    s.step("the nearest chamber, lit", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        Vector3 at{}, out{};
        C.found = c.field->nearestCaveMouth(home, 3000.0f, at, out, 3);
        if (!c.check("a chamber within 3 km of the start", C.found)) return true;
        C.up = Vector3Normalize(at); C.e1 = Tangent(C.up); C.e2 = Vector3CrossProduct(C.up, C.e1);
        // Its floor: down from the centre to the first solid sample (density > 0 is rock).
        C.floor = at;
        for (float d = 0.0f; d < 30.0f; d += 0.25f) {
            const Vector3 p = Vector3Subtract(at, Vector3Scale(C.up, d));
            if (c.field->density(p, nullptr) > 0.0f) { C.floor = Vector3Add(p, Vector3Scale(C.up, 0.3f)); break; }
        }
        // Lamps in a ring 6-14 m around the centre on the floor: a hall lit from within, not by a light at the eye.
        Planetoid* p = c.body;
        if (!p->plots) p->plots = std::make_unique<PlotSet>();
        const int lampMat = std::max(0, BuildMaterialForItem("lamp"));
        int placed = 0;
        for (int i = 0; i < 10; ++i) {
            const float a = 6.2831853f * i / 10.0f, r = i % 2 ? 6.0f : 14.0f;
            Vector3 q = Vector3Add(C.floor, Vector3Add(Vector3Scale(C.e1, r * cosf(a)), Vector3Scale(C.e2, r * sinf(a))));
            for (float d = -4.0f; d < 4.0f; d += 0.25f) { // settle onto the floor there
                const Vector3 t = Vector3Subtract(q, Vector3Scale(C.up, d));
                if (c.field->density(t, nullptr) <= 0.0f && c.field->density(Vector3Subtract(t, Vector3Scale(C.up, 0.3f)), nullptr) > 0.0f) { q = t; break; }
            }
            Plot* pl = p->plots->plotFor(q, 0.0f, true);
            placed += pl->place(Part{ PartType::Cube, pl->cellAt(pl->toPlot(q)), 0, (uint8_t)lampMat }) >= 0 ? 1 : 0;
        }
        c.metric("lamps placed in the chamber", placed, "", false);
        return true;
    });
    for (int k = 0; k < 3; ++k)
        View(s, "hall" + std::to_string(k), [k](ScenarioCtx& c, Vector3& eye, Vector3& target, Vector3& up) {
            if (!C.found) return false;
            const float a = 2.0943951f * k; // three sides, from 2.5 m over the floor at the hall's edge, across it
            const Vector3 dir = Vector3Add(Vector3Scale(C.e1, cosf(a)), Vector3Scale(C.e2, sinf(a)));
            eye = Vector3Add(C.floor, Vector3Add(Vector3Scale(dir, -10.0f), Vector3Scale(C.up, 2.5f)));
            target = Vector3Add(C.floor, Vector3Add(Vector3Scale(dir, 10.0f), Vector3Scale(C.up, 3.0f)));
            up = C.up;
            return true;
        });
}

SCENARIO("slides/lava-tube") {
    Common(s);
    s.step("a skylight on the nearest volcano flank", 1.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        const Vector3 k1 = Tangent(home), k2 = Vector3CrossProduct(home, k1);
        float best = 1e30f; Vector3 near = home; // (main.cpp's RAYLIBGAME_CAVES_TUBE search: a 1 km grid out to 60 km)
        for (float a = -60000.0f; a <= 60000.0f; a += 1000.0f)
            for (float b = -60000.0f; b <= 60000.0f; b += 1000.0f) {
                const float dd = a * a + b * b;
                if (dd >= best) continue;
                const Vector3 vd = Vector3Normalize(Vector3Add(Vector3Scale(home, c.field->radius), Vector3Add(Vector3Scale(k1, a), Vector3Scale(k2, b))));
                const float w = c.field->volcanoWeight(Vector3Scale(vd, c.field->radius));
                if (w > 0.4f && w < 0.8f && c.field->elevation(vd) > 10.0f) { best = dd; near = vd; }
            }
        T.found = best < 1e29f && c.field->nearestCaveMouth(near, 3000.0f, T.at, T.out, 1);
        if (!c.check("a lava tube skylight near a volcano flank", T.found)) return true;
        T.up = Vector3Normalize(T.at);
        c.metric("skylight from the start (m)", sqrtf(best), "m");
        return true;
    });
    // From inside the tube, 8 m under the skylight and 6 m along, looking back up at the light coming in; and from outside.
    View(s, "inside", [](ScenarioCtx& c, Vector3& eye, Vector3& target, Vector3& up) {
        if (!T.found) return false;
        const Vector3 along = Vector3Normalize(Vector3Subtract(T.out, Vector3Scale(T.up, Vector3DotProduct(T.out, T.up))));
        eye = Vector3Add(T.at, Vector3Add(Vector3Scale(T.up, -8.0f), Vector3Scale(along, 6.0f)));
        target = Vector3Add(T.at, Vector3Scale(T.up, -2.0f));
        up = T.up;
        return true;
    });
}

SCENARIO("slides/space") {
    Common(s);
    // The planet's frame: its axis is local y, the ring in its equatorial (xz) plane. The sun's direction in that frame
    // decides where the day side is: the first pass, over the start, looked at the night side.
    View(s, "limb", [](ScenarioCtx& c, Vector3& eye, Vector3& target, Vector3& up) {
        const Vector3 sun = SunLocal(c);
        const Vector3 side = Tangent(sun);
        const Vector3 over = Vector3Normalize(Vector3Add(Vector3Scale(sun, 0.6f), Vector3Scale(side, 0.8f))); // afternoon light
        const Vector3 away = Vector3Normalize(Vector3Subtract(side, Vector3Scale(over, Vector3DotProduct(side, over)))); // the sun behind
        const float R = c.field->radius;
        eye = Vector3Scale(over, R + 25000.0f); // 25 km up, the sun behind the camera
        // From 25 km over a 60 km planet the horizon dips ~45 degrees: aim 35 degrees down, the limb in the upper frame.
        target = Vector3Add(eye, Vector3Scale(Vector3Subtract(Vector3Scale(away, 0.819f), Vector3Scale(over, 0.574f)), 100000.0f));
        up = over;
        return true;
    });
    View(s, "planet-ring", [](ScenarioCtx& c, Vector3& eye, Vector3& target, Vector3& up) {
        const float R = c.field->radius;
        eye = { R * 3.2f, R * 0.55f, R * 1.2f }; // a little above the ring plane, the whole planet in view
        target = { 0, 0, 0 }; up = { 0, 1, 0 };
        return true;
    });
    View(s, "moon", [](ScenarioCtx& c, Vector3& eye, Vector3& target, Vector3& up) {
        const Planetoid* moon = nullptr;
        for (const auto& o : c.world.objects)
            if (const Planetoid* p = dynamic_cast<const Planetoid*>(o.get()); p && p != c.body && p->name == "luna") moon = p;
        if (!c.check("a body named luna", moon != nullptr)) return false;
        const Vector3 m = c.body->toLocalD(moon->worldPos).toFloat();
        // From the moon's sunlit side, a little off the sun line (the first pass looked at its night side).
        const Vector3 sun = SunLocal(c);
        eye = Vector3Add(m, Vector3Add(Vector3Scale(sun, moon->radius * 3.0f), Vector3Scale(Tangent(sun), moon->radius * 1.5f)));
        target = m; up = { 0, 1, 0 };
        return true;
    });
}
