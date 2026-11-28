// fx/dust-near-eye (#1126): dig dust drifting into the eye no longer fills the screen. Mining a wall in a tight tunnel
// threw puffs of the ore's colour (copper green, rusty orange) into the camera, where a 0.35-0.8 m billboard a few cm
// from the near plane covered the view. Now particles fade out near the eye (Particles_NearFade).
// On the sandbox, looking at the ground ahead: a reference frame, then a dig stroke's dust (the same burst as
// StrokeParticles, in copper ore's green) let loose 0.3 m in front of the eye. Over the next half second the share of
// pixels that changed against the reference stays small; with RAYLIBGAME_PARTICLE_NEAR_FADE=0 it is most of the screen.
// A second burst 2 m away still shows (the fade is only near the eye).
#include "scenario.h"
#include "particles.h"
#include "playersession.h"
#include "world.h"
#include <raylib.h>
#include <raymath.h>
#include <algorithm>
#include <string>

namespace {
struct State { Image ref{}; bool haveRef = false; double worstNear = 0.0, worstFar = 0.0; int frames = 0; } S;

void Look(ScenarioCtx& c) { // the eye 1.6 m over the ground, looking ahead and a little down; set each frame
    Player& p = LocalSession().player;
    const Vector3 feet = (p.position - c.world.renderOrigin).toFloat();
    const Vector3 up = Vector3Normalize(p.up);
    Vector3 fwd = Vector3Subtract(p.forward, Vector3Scale(up, Vector3DotProduct(p.forward, up)));
    fwd = Vector3Length(fwd) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(fwd);
    c.world.camera.position = Vector3Add(feet, Vector3Scale(up, 1.6f));
    c.world.camera.target = Vector3Add(feet, Vector3Scale(fwd, 6.0f));
    c.world.camera.up = up;
}

// Percent of pixels (every 4th in each direction) differing from the reference by more than 24 in any channel.
double ChangedPct(const Image& ref) {
    Image now = LoadImageFromScreen();
    if (now.width != ref.width || now.height != ref.height) { UnloadImage(now); return 100.0; }
    Color* a = LoadImageColors(ref);
    Color* b = LoadImageColors(now);
    long changed = 0, total = 0;
    for (int y = 0; y < now.height; y += 4)
        for (int x = 0; x < now.width; x += 4) {
            const Color p = a[y * now.width + x], q = b[y * now.width + x];
            const int d = std::max({ abs(p.r - q.r), abs(p.g - q.g), abs(p.b - q.b) });
            changed += d > 24 ? 1 : 0; ++total;
        }
    UnloadImageColors(a); UnloadImageColors(b); UnloadImage(now);
    return total ? 100.0 * changed / total : 0.0;
}

void Burst(ScenarioCtx& c, float ahead) { // StrokeParticles' dust in copper ore's green, `ahead` m in front of the eye
    const Vector3 eye = c.world.camera.position;
    const Vector3 fwd = Vector3Normalize(Vector3Subtract(c.world.camera.target, eye));
    const DVec3 at = c.world.renderOrigin + Vector3Add(eye, Vector3Scale(fwd, ahead));
    ParticleBurst dust;
    dust.kind = ParticleKind::Dust; dust.color = { 70, 150, 90, 200 }; dust.count = 10;
    dust.size = 0.35f; dust.life = 1.1f; dust.speed = 0.4f; dust.spread = 0.9f; // slower than a stroke's: it stays in front of the eye
    Particles_EmitWorld(*c.body, at, dust, Vector3{ 0, 0, 0 }, Vector3{ 0, 0, 0 });
}
}

SCENARIO("fx/dust-near-eye") {
    s.tags = { "fx" };
    s.needsRender = true; // (#2130) it compares screen grabs
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" }; // mid-morning: dust lit as by day (at night it is dark on dark)
    s.timeout = 40.0;
    s.step("settle, a reference frame", 15.0, [](ScenarioCtx& c) {
        c.waitingFor("the player standing on the sandbox");
        Look(c);
        if (c.T < 2.0) return false;
        Particles_Clear();
        if (c.t < 0.5) return false; // a few frames with no particles drawn
        S.ref = LoadImageFromScreen(); S.haveRef = true;
        return true;
    });
    s.step("dust let loose 0.3 m from the eye: it does not fill the screen", 5.0, [](ScenarioCtx& c) {
        Look(c);
        if (S.frames == 0) { Burst(c, 0.3f); S.frames = 1; c.cue(); return false; }
        S.worstNear = std::max(S.worstNear, ChangedPct(S.ref));
        if (c.t < 0.5) return false;
        c.note(TextFormat("worst share of the screen changed by the near dust: %.1f %%", S.worstNear));
        c.lt("near dust: share of the screen changed (%)", S.worstNear, 15.0);
        Particles_Clear();
        S.frames = 0;
        return true;
    });
    s.step("dust 2 m away still shows", 5.0, [](ScenarioCtx& c) {
        Look(c);
        if (S.frames == 0) { if (c.t < 0.3) return false; Burst(c, 2.0f); S.frames = 1; return false; }
        S.worstFar = std::max(S.worstFar, ChangedPct(S.ref));
        if (c.t < 0.8) return false;
        c.note(TextFormat("worst share of the screen changed by the far dust: %.1f %%", S.worstFar));
        c.gt("far dust: some of the screen changed (%)", S.worstFar, 0.3);
        if (S.haveRef) { UnloadImage(S.ref); S.haveRef = false; }
        return true;
    });
}
