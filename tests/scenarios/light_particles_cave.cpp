// light/particles-cave (#2087): smoke and dust underground take the cave's light, not the sun's angle. A sealed 16 x 16 x 4 m
// room is carved 10 m down at noon; the player stands at one end while a dust cloud is thrown up 5 m ahead. Before the fix the
// dust was lit by the noon sun through the rock (bright in a black room) and, at night, near-black even in the headlamp's
// beam. Checks: the sky a burst sees is ~0 in the room and ~1 over the ground; the lit shader drew the dust. The capture shows
// the cloud dark with the lamp off, then lit by the headlamp: for a person to judge.
#include "scenario.h"
#include "inventory.h"
#include "particles.h"
#include "player.h"
#include "suit.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>

namespace {
struct State {
    Vector3 home{}, ax{}, ay{};
    float r0 = 0.0f;
    Vector3 pt(float x, float y, float h) const {
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(home, r0), Vector3Add(Vector3Scale(ax, x), Vector3Scale(ay, y))));
        return Vector3Scale(d, r0 + h);
    }
    TerrainEdit box(float x0, float x1, float y0, float y1, float h0, float h1) const {
        TerrainEdit e{};
        e.center = pt(0.5f * (x0 + x1), 0.5f * (y0 + y1), 0.5f * (h0 + h1));
        e.carve = true; e.box = true; e.half = { 0.5f * (x1 - x0), 0.5f * (y1 - y0), 0.5f * (h1 - h0) };
        e.ax = ax; e.ay = ay; e.az = Vector3Normalize(e.center);
        e.radius = Vector3Length(e.half);
        return e;
    }
    bool litSeen = false;
} S;
constexpr float kFloor = -10.0f, kCeil = -6.0f; // m under the surface
// Stand at the room's west end facing east (the player's own eye and headlamp), and keep the dust coming 5 m ahead.
void Stand(ScenarioCtx& c) {
    if (c.t < 0.05) { player.flyMode = false; player.placeAt(c.world, *c.body, c.body->toWorld(S.pt(-6.0f, 0.0f, kFloor + 0.1f))); }
    player.forward = c.body->dirToWorld(S.ax);
}
void Dust(ScenarioCtx& c) {
    ParticleBurst b;
    b.kind = ParticleKind::Dust; b.color = { 196, 184, 160, 220 }; b.count = 3; b.size = 0.6f; b.life = 2.0f; b.speed = 0.4f; b.spread = 1.0f;
    Particles_Emit(*c.body, S.pt(-1.0f, 0.0f, kFloor + 1.4f), b);
    S.litSeen = S.litSeen || Particles_LastDrawLit();
}
}

SCENARIO("light/particles-cave") {
    s.tags = { "light", "visual" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--sun-elev", "60" }; // noon: the old shading lit cave dust fully
    s.timeout = 150.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("carve a sealed room", 1.0, [](ScenarioCtx& c) {
        S.home = c.field->findSpawnDir(c.body->seed);
        S.ax = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 })); S.ay = Vector3CrossProduct(S.ax, S.home);
        S.r0 = c.field->radius + c.field->elevation(S.home);
        TerrainEdit room = S.box(-8.0f, 8.0f, -8.0f, 8.0f, kFloor, kCeil);
        ApplyTerrainEdits(*c.body, &room, 1);
        inventory.add("battery_pack", 1); // suit power: the basic suit's headlamp needs one (#217)
        c.check("a battery pack fitted", Suit_Fit("battery_pack"));
        survival.headlamp = false;
        return true;
    });
    s.step("in the room, meshed", 30.0, [](ScenarioCtx& c) {
        Stand(c);
        c.waitingFor("an idle world round the player");
        if (!(c.t > 1.0 && WorldIdle())) return false;
        c.lt("the sky a burst sees in the room", Particles_SkyAt(*c.body, S.pt(0.0f, 0.0f, kFloor + 1.4f)), 0.2);
        c.gt("the sky a burst sees over the ground", Particles_SkyAt(*c.body, S.pt(0.0f, 0.0f, 3.0f)), 0.9);
        return true;
    });
    s.step("lamp off: the dust dark", 3.0, [](ScenarioCtx& c) {
        static bool cued = false;
        if (!cued) { cued = true; c.cue(); }
        Stand(c); Dust(c);
        survival.headlamp = false;
        return c.t > 2.5;
    });
    s.step("headlamp on: the dust lit", 3.0, [](ScenarioCtx& c) {
        Stand(c); Dust(c);
        survival.headlamp = true; survival.power = 100.0f;
        return c.t > 2.5;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) {
        c.check("the lit particle shader drew the dust", S.litSeen);
        c.keepCapture("cave dust dark, then lit by the headlamp (#2087): for a person to judge");
        return true;
    });
}
