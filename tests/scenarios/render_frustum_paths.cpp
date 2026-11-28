// render/frustum-paths (#2140, #2172, #2177): the extra 3D paths skip what is out of the main view (Scene::inView). Before,
// creatures to 160 m were drawn all around the camera, and particles to 1.5 km everywhere but behind it.
// Sandbox: 48 creatures in a ring 12-40 m around the player and 48 still smoke puffs in a ring 15 m out. The camera stands
// 3 m up and looks level along six headings 60 degrees apart with culling on, then once with the frustum toggle off
// (creatures: all drawn; particles: the front half, the old behind-the-camera test). Each heading draws well under that,
// and the six together see every one (neighbouring views overlap).
#include "scenario.h"
#include "animals.h"
#include "object.h"
#include "particles.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace {
constexpr int kCount = 48;
constexpr int kHeadings = 6;
constexpr int kFramesPerView = 6; // the counters read the frame before: let a few frames go by on each view
struct St { Vector3 feet{}, up{}, e{}, n{}; int view = 0, frames = 0; int animals[kHeadings + 1] = {}, particles[kHeadings + 1] = {}; } S;

void Look(ScenarioCtx& c, int k) {
    const float a = 2.0f * PI * k / kHeadings;
    const Vector3 eye = Vector3Add(S.feet, Vector3Scale(S.up, 3.0f));
    const Vector3 d = Vector3Add(Vector3Scale(S.e, cosf(a)), Vector3Scale(S.n, sinf(a)));
    c.camera(eye, Vector3Add(eye, Vector3Scale(d, 20.0f)), S.up);
}

void Puffs(ScenarioCtx& c) { // still smoke, long-lived: one puff each, 15 m out and 2 m up
    for (int k = 0; k < kCount; ++k) {
        const float a = 2.0f * PI * (k + 0.5f) / kCount;
        const Vector3 local = Vector3Add(S.feet, Vector3Add(Vector3Scale(S.up, 2.0f), Vector3Add(Vector3Scale(S.e, 15.0f * cosf(a)), Vector3Scale(S.n, 15.0f * sinf(a)))));
        ParticleBurst b;
        b.kind = ParticleKind::Smoke; b.color = { 200, 200, 200, 200 }; b.count = 1;
        b.size = 0.4f; b.life = 60.0f; b.speed = 0.0f; b.spread = 0.0f;
        Particles_Emit(*c.body, local, b);
    }
}
}

SCENARIO("render/frustum-paths") {
    s.tags = { "visual", "culling" }; // it reads draw counters: the PC (the runner's headless parts draw nothing)
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 60.0;

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false);
        Animals_SetThreats(false);
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("creatures and smoke in a ring", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        S.feet = c.body->toLocal((player.position - c.world.renderOrigin).toFloat());
        S.up = Vector3Normalize(S.feet);
        S.e = Vector3Normalize(Vector3CrossProduct(S.up, fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        S.n = Vector3CrossProduct(S.e, S.up);
        const int before = Animals_Count();
        for (int k = 0; k < kCount; ++k) {
            const float a = 2.0f * PI * (k + 0.5f) / kCount, r = 12.0f + 28.0f * (k % 5) / 4.0f;
            const Vector3 dir = Vector3Normalize(Vector3Add(S.feet, Vector3Add(Vector3Scale(S.e, r * cosf(a)), Vector3Scale(S.n, r * sinf(a)))));
            Animals_Spawn(c.world, *c.body, "woolback", Vector3Scale(dir, c.field->radius + c.field->elevation(dir)), 1, false, 0.01f);
        }
        Particles_Clear();
        Puffs(c);
        c.gt("creatures spawned", Animals_Count() - before, kCount * 0.9);
        return true;
    });
    s.step("six headings, culling on; then culling off", 20.0, [](ScenarioCtx& c) {
        survival.health = 100.0f; survival.dead = false;
        frustumCullingEnabled = S.view < kHeadings;
        Look(c, S.view % kHeadings);
        if (++S.frames < kFramesPerView) return false;
        S.animals[S.view] = Animals_DrawnCount();
        S.particles[S.view] = Particles_DrawnCount();
        S.frames = 0;
        if (++S.view <= kHeadings) return false;
        frustumCullingEnabled = true;
        const int allA = S.animals[kHeadings], frontP = S.particles[kHeadings], aliveP = Particles_Count();
        c.gt("creatures drawn, culling off", allA, kCount * 0.9);
        c.gt("particles alive", aliveP, kCount * 0.9);
        c.gt("particles drawn, culling off (the front half)", frontP, aliveP * 0.4);
        int sumA = 0, sumP = 0, maxP = 0;
        for (int k = 0; k < kHeadings; ++k) {
            const std::string h = std::to_string(k * 360 / kHeadings);
            c.lt("creatures drawn looking " + h + " deg", S.animals[k], allA * 0.5);
            c.lt("particles drawn looking " + h + " deg", S.particles[k], aliveP * 0.5);
            sumA += S.animals[k]; sumP += S.particles[k]; maxP = std::max(maxP, S.particles[k]);
        }
        c.gt("particles beside the view culled (the front half over the most in one view)", frontP, maxP * 1.5);
        c.gt("creatures seen over the six headings", sumA, allA - 0.5);
        c.gt("particles seen over the six headings", sumP, aliveP - 0.5);
        c.cue();
        return true;
    });
    // For the recording (pop-in at the edges): a slow pan, 45 deg/s, across the ring. With RAYLIBGAME_FRAME_CULL=0 the
    // frames must come out the same.
    s.step("a slow pan for the recording", 6.0, [](ScenarioCtx& c) {
        survival.health = 100.0f; survival.dead = false;
        const float a = (float)c.t * PI / 4.0f;
        const Vector3 eye = Vector3Add(S.feet, Vector3Scale(S.up, 3.0f));
        const Vector3 d = Vector3Add(Vector3Scale(S.e, cosf(a)), Vector3Scale(S.n, sinf(a)));
        c.camera(eye, Vector3Add(eye, Vector3Scale(d, 20.0f)), S.up);
        return c.t >= 4.0;
    });
}
