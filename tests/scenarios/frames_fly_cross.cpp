// frames/fly-cross (Plan 11 F4, #1412): the fly camera crosses its reference body's altitude line without a step. Following
// the moon (it orbits the planet, so a frame without the carry leaves the camera a frame of orbital motion behind), the
// camera is nudged straight down and up again, 0.2 moon radii a second, from 2.6 to 1.4 radii and back: it crosses the line at
// 2 radii twice, where the carry switches between the body's spinning frame (inside) and its centre (outside). The view looks
// at the moon. Each frame, the camera's move beyond the nudge is measured both in the moon's spinning frame and from its
// centre; the smaller of the two is the step (the right frame gives ~0 on either side of the line, a missed or wrong carry
// gives both ~the orbital speed x a frame). Gate: no frame's step over 0.05 m. frames/fly-cross-pose runs it with
// RAYLIBGAME_FRAMES=1. (A ring rock's sphere of influence, the gate's second case, is not covered: ring rocks are not
// fly-follow references.)
#include "scenario.h"
#include "player.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
constexpr double kRate = 0.2;            // moon radii a second
constexpr double kFrom = 2.6, kTo = 1.4; // moon radii from its centre
struct CrossState {
    Planetoid* moon = nullptr;
    bool placed = false, have = false;
    double t0 = 0.0, r = 0.0;            // the scenario's own radius (moon radii): what the nudges add up to
    DVec3 prevLocal{}, prevOff{};        // last frame, after its nudge: body-local, and from the moon's centre (world axes)
    double lastT = 0.0;
    int frames = 0, crossings = 0; bool wasInside = false;
    double worst = 0.0, worstAt = 0.0, sum = 0.0;
};
CrossState st;

DVec3 CamWorld(ScenarioCtx& c) { return c.world.renderOrigin + DVec3(c.world.camera.position); }

void Build(ScenarioDef& s) {
    s.tags = { "frames", "camera" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 200.0;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        return c.T > 1.0 && WorldIdle();
    });
    s.step("fly down through the moon's altitude line and back", 40.0, [](ScenarioCtx& c) {
        if (!st.placed) {
            st = CrossState{};
            for (auto& obj : c.world.objects)
                if (auto* q = dynamic_cast<Planetoid*>(obj.get())) if (q->kind == BodyKind::Moon && q->orbitParent) { st.moon = q; break; }
            if (!st.moon) { c.fail("no moon"); return true; }
            Planetoid& m = *st.moon;
            player.flyMode = true;
            player.flyFollow = true;
            st.r = kFrom;
            const Vector3 dir = Vector3Normalize({ 0.3f, 1.0f, 0.2f });
            c.world.camera.position = m.toWorld(Vector3Scale(dir, (float)(m.radius * m.scale * st.r)));
            c.world.camera.target = Vector3Add(c.world.camera.position, m.dirToWorld(Vector3Negate(dir))); // looking at the moon
            c.world.camera.up = m.dirToWorld(Vector3Normalize(Vector3CrossProduct(dir, { 1, 0, 0 })));
            st.placed = true;
            st.t0 = c.T;
            st.wasInside = false;
            return false;
        }
        Planetoid& m = *st.moon;
        const double R = m.radius * m.scale;
        const DVec3 cam = CamWorld(c);
        const DVec3 local = m.toLocalD(cam), off = cam - m.worldPos;
        const double run = c.T - st.t0;
        if (st.have && run > 0.5) { // (the first half second: the carry takes its first recorded pose)
            const double stepLocal = (local - st.prevLocal).length();
            const double stepOff = (off - st.prevOff).length();
            const double step = std::min(stepLocal, stepOff);
            ++st.frames; st.sum += step;
            if (step > st.worst) { st.worst = step; st.worstAt = off.length() / R; }
        }
        const bool inside = off.length() < 2.0 * R;
        if (st.have && inside != st.wasInside) ++st.crossings;
        st.wasInside = inside;
        // The nudge (as input moves the fly camera): straight down for the first half, then up.
        const double half = (kFrom - kTo) / kRate;
        if (run >= 2.0 * half) {
            c.note(TextFormat("%d frames, %d crossings of the 2-radius line; the camera's step beyond its nudge: worst %.4f m (at %.2f radii), mean %.5f m (%s)",
                              st.frames, st.crossings, st.worst, st.worstAt, st.frames ? st.sum / st.frames : 0.0, Frames_On() ? "frames on" : "frames off"));
            c.metric("worst step at a frame", st.worst, "m");
            c.check("crossed the altitude line both ways", st.crossings >= 2, TextFormat("%d crossings", st.crossings));
            c.lt("the fly camera's worst step beyond its nudge (m)", st.worst, 0.05);
            return true;
        }
        const double frameDt = st.have ? std::max(1e-4, c.T - st.lastT) : 1.0 / 60.0;
        const double dr = kRate * frameDt * (run < half ? -1.0 : 1.0); // moon radii this frame
        const DVec3 dirW = off * (1.0 / off.length());
        const DVec3 move = dirW * (dr * R);
        c.world.camera.position = (cam + move - c.world.renderOrigin).toFloat();
        c.world.camera.target = Vector3Add(c.world.camera.position, (m.worldPos - (cam + move)).toFloat());
        st.prevLocal = m.toLocalD(cam + move); // the pose after the nudge: what the carry must keep
        st.prevOff = cam + move - m.worldPos;
        st.lastT = c.T;
        st.have = true;
        c.waitingFor("down to 1.4 moon radii and back up");
        return false;
    });
}
}

SCENARIO("frames/fly-cross") { Build(s); }
SCENARIO("frames/fly-cross-pose") { s.env = { { "RAYLIBGAME_FRAMES", "1" } }; Build(s); }
