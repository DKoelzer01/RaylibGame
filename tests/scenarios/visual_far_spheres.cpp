// visual/far-sphere-* and visual/dust-burst (#1859 L5, L6, L8, L9): views to record, not checks (run with the capture skill).
// - far-sphere-day / far-sphere-night: Luna from Terra's surface, where the moon stands high and the sun is up (or down).
//   Luna draws as a far sphere: its east-west relief (L5) and which side the sun lights (L6: at night it was lit from the
//   moon's own direction).
// - asteroid-night: Goliath's night side from 60 m up (L9: water and grass took a night-side shadow sphere, as on planets).
// - dust-burst: default dust bursts (dir zero: away from the body centre, L8) at seed 42's start, 25-45 degrees from the
//   equator; before L8 they went along body +Y, sideways there.
#include "scenario.h"
#include "particles.h"
#include "player.h"
#include "voxel.h"
#include "weather.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <cmath>

namespace {
constexpr double kHold = 4.0; // s of the view held for the recording

Planetoid* BodyOfKind(ScenarioCtx& c, BodyKind kind, const char* name = nullptr) {
    for (auto& o : c.world.objects)
        if (auto* p = dynamic_cast<Planetoid*>(o.get()); p && p->kind == kind && (!name || p->name == name)) return p;
    return nullptr;
}
// Any tangent of a unit direction (the camera's up when it looks along the direction).
Vector3 Tangent(Vector3 d) { return Vector3Normalize(Vector3CrossProduct(d, fabsf(d.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })); }
// The camera at a body-local eye on `b`, looking along a body-local direction.
void Look(ScenarioCtx& c, Planetoid& b, Vector3 eyeLocal, Vector3 lookLocal, Vector3 upLocal) {
    gScenarioCamera = true;
    c.world.camera.position = b.toWorld(eyeLocal);
    c.world.camera.target = b.toWorld(Vector3Add(eyeLocal, Vector3Scale(lookLocal, 100.0f)));
    c.world.camera.up = b.dirToWorld(upLocal);
    player.flyMode = true;
    player.position = c.world.renderOrigin + c.world.camera.position;
}

struct St { bool set = false; double t0 = 0.0; Vector3 eye{}, look{}, up{}; } S;

void FarSphereView(ScenarioDef& s, bool night) {
    s.tags = { "visual" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 180.0;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 1.0 && WorldIdle(); });
    s.step(night ? "Luna overhead at night" : "Luna overhead by day", 60.0, [night](ScenarioCtx& c) {
        Planetoid* moon = BodyOfKind(c, BodyKind::Moon);
        if (!moon) { c.fail("no moon"); return true; }
        if (!S.set) {
            // The eye under Luna (the moon at the zenith); the sun 90 degrees round from it in azimuth, 40 degrees up by day
            // and 30 down at night, so the moon is half lit and the terminator shows its relief. Nothing else sets the sun
            // (no --sun override), and the weather is held clear so no cloud hides the moon.
            const Vector3 up = Vector3Normalize(c.body->toLocalD(moon->worldPos).toFloat());
            const Vector3 side = Tangent(up);
            const float a = (night ? -30.0f : 40.0f) * DEG2RAD;
            c.world.sunDirection = c.body->dirToWorld(Vector3Normalize(Vector3Add(Vector3Scale(side, cosf(a)), Vector3Scale(up, sinf(a)))));
            Weather_Force((int)WeatherKind::Clear);
            const float e = c.field ? c.field->elevation(up) : 0.0f;
            S.eye = Vector3Scale(up, c.body->radius + std::max(e, 0.0f) + 150.0f);
            S.up = Vector3CrossProduct(side, up); // the sun to the frame's side
            S.set = true; S.t0 = c.T;
            c.cue();
        }
        // (the moon moves a little in 4 s: re-aimed each frame)
        S.look = Vector3Normalize(Vector3Subtract(c.body->toLocalD(moon->worldPos).toFloat(), S.eye));
        Look(c, *c.body, S.eye, S.look, S.up);
        return c.T - S.t0 >= kHold;
    });
}
}

SCENARIO("visual/far-sphere-day") { FarSphereView(s, false); }
SCENARIO("visual/far-sphere-night") { FarSphereView(s, true); }

SCENARIO("visual/asteroid-night") {
    s.tags = { "visual" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 180.0;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 1.0 && WorldIdle(); });
    s.step("Goliath's night side", 60.0, [](ScenarioCtx& c) {
        Planetoid* a = BodyOfKind(c, BodyKind::Asteroid, "goliath");
        if (!a) { c.fail("no goliath"); return true; }
        if (!S.set) {
            // Over the night side, 40 degrees past the terminator, looking along the ground toward the sunlit limb.
            const Vector3 sun = a->dirToLocal(c.world.sunDirection);
            const Vector3 t = Tangent(sun);
            const Vector3 up = Vector3Normalize(Vector3Add(Vector3Scale(sun, -sinf(40.0f * DEG2RAD)), Vector3Scale(t, cosf(40.0f * DEG2RAD))));
            S.eye = Vector3Scale(up, a->radius + 60.0f);
            S.look = Vector3Normalize(Vector3Add(Vector3Scale(up, -0.35f), Vector3Scale(sun, 1.0f)));
            S.up = up;
            S.set = true; S.t0 = c.T;
            c.cue();
        }
        Look(c, *a, S.eye, S.look, S.up);
        return c.T - S.t0 >= kHold;
    });
}

SCENARIO("visual/dust-burst") {
    s.tags = { "visual" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 180.0;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 1.0 && WorldIdle(); });
    s.step("default dust bursts at the start", 60.0, [](ScenarioCtx& c) {
        const Vector3 up = c.spawnDirLocal();
        const float ground = c.body->radius + c.field->elevation(up);
        const Vector3 at = Vector3Scale(up, ground + 0.2f);
        if (!S.set) {
            S.up = up;
            S.look = Tangent(up); // level: the bursts rise straight up the frame when they go away from the centre
            S.eye = Vector3Subtract(Vector3Scale(up, ground + 1.6f), Vector3Scale(S.look, 6.0f));
            S.set = true; S.t0 = c.T;
            c.note(TextFormat("start %.0f deg from the equator", asinf(fabsf(up.y)) * RAD2DEG));
            c.cue();
        }
        static double last = -1.0;
        if (c.T - last > 0.25) { // a default burst (only the look set: dir stays zero)
            last = c.T;
            ParticleBurst b;
            b.count = 40; b.size = 0.3f; b.life = 1.6f; b.speed = 2.0f; b.spread = 0.15f; b.color = { 235, 225, 210, 255 };
            Particles_Emit(*c.body, at, b);
        }
        Look(c, *c.body, S.eye, S.look, S.up);
        return c.T - S.t0 >= kHold;
    });
}
