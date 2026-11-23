// frames/ride (Plan 11 F2, #667, #1143): a deterministic --test-ride. On the asteroid system's first body, set turning
// 10 deg/s about its Y axis, the player stands at a fixed body-local spot on its equator; once grounded and settled, the
// view and the feet are watched for 3 s of sim time in the body's own frame (the view must turn with the ground, not
// drift), then the player jumps and the landing is measured from the takeoff (body frame: momentum kept, it lands where it
// took off). --test-ride spawned at a world direction at the wall-clock moment the world went idle, so its spot (and its
// numbers) changed from run to run; slot 2 of #1144's gate saw 12.6 deg of view drift once with RAYLIBGAME_FRAMES=1. Gates:
// view drift < 0.5 deg and feet < 5 cm standing; the jump's landing is reported (spin forces move it metres here).
// frames/ride-pose runs the same with RAYLIBGAME_FRAMES=1.
#include "scenario.h"
#include "player.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
constexpr float kSpinDeg = 10.0f;
struct RideState {
    bool placed = false, settled = false, jumped = false, airborne = false;
    double t0 = 0.0;
    Vector3 fwd0{};
    DVec3 feet0{}, takeoff{};
    double maxStepTurn = 0.0, viewDrift = 0.0, posDrift = 0.0;
    Vector3 lastFwd{};
    int frames = 0, grounded = 0;
};
RideState st;

static void Build(ScenarioDef& s) {
    s.tags = { "frames", "physics" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--bodies", "resources/data/bodies_asteroids.txt" };
    s.timeout = 240.0;
    s.record = false;
    s.needsPlanetField = false; // (the asteroid system's first body)
    s.step("an idle world", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        return c.T > 1.0 && WorldIdle();
    });
    s.step("stand on the spinning body", 60.0, [](ScenarioCtx& c) {
        Planetoid& b = *c.body;
        if (!st.placed) {
            st = RideState{};
            b.spin = { 0.0f, kSpinDeg, 0.0f };
            player.spawnOn(c.world, b, b.dirToWorld({ 1.0f, 0.0f, 0.0f })); // body-local +X: the equator, always the same spot
            st.placed = true;
            st.t0 = c.T;
            return false;
        }
        if (!st.settled) {
            c.waitingFor("grounded and 2 s settled");
            if (!player.grounded) { st.t0 = c.T; return false; }
            if (c.T - st.t0 < 2.0) return false;
            st.settled = true;
            st.t0 = c.T;
            st.fwd0 = st.lastFwd = b.dirToLocal(player.forward);
            st.feet0 = b.toLocalD(player.position);
            return false;
        }
        const Vector3 fwd = b.dirToLocal(player.forward);
        st.maxStepTurn = std::max(st.maxStepTurn, (double)Vector3Angle(st.lastFwd, fwd) * RAD2DEG);
        st.lastFwd = fwd;
        ++st.frames; if (player.grounded) ++st.grounded;
        c.waitingFor("3 s standing");
        if (c.T - st.t0 < 3.0) return false;
        st.viewDrift = Vector3Angle(st.fwd0, fwd) * RAD2DEG;
        st.posDrift = (b.toLocalD(player.position) - st.feet0).length();
        const Vector3 feet = (player.position - c.world.renderOrigin).toFloat();
        c.note(TextFormat("standing 3 s at %.0f deg/s, surface %.2f m/s, grounded %.0f%%: view drift %.3f deg (largest single step %.4f deg), feet %.4f m",
                          kSpinDeg, Vector3Length(b.surfaceVelocity(feet)), 100.0 * st.grounded / std::max(1, st.frames), st.viewDrift, st.maxStepTurn, st.posDrift));
        if (Frames_On()) // (F2: is the frame path in use here? the frame must be Fixed on this body)
            c.note(TextFormat("frame: body %llx (this body %llx), %s; frame speed %.4f m/s", (unsigned long long)player.frame.frame.body,
                              (unsigned long long)b.netId, player.frame.frame.kind == FrameKind::Fixed ? "Fixed" : "Centred", player.frame.vel.length()));
        c.metric("view drift", st.viewDrift, "deg");
        c.metric("largest single-step view turn", st.maxStepTurn, "deg");
        c.metric("feet drift", st.posDrift, "m");
        c.lt("view drift", st.viewDrift, 0.5);
        c.lt("feet drift", st.posDrift, 0.05);
        st.takeoff = b.toLocalD(player.position);
        player.autoJump = true;
        st.t0 = c.T;
        return true;
    });
    s.step("jump and land", 30.0, [](ScenarioCtx& c) {
        Planetoid& b = *c.body;
        if (!st.airborne) { c.waitingFor("takeoff"); if (!player.grounded) st.airborne = true; return false; }
        c.waitingFor("landing");
        if (!player.grounded) return false;
        const double land = (b.toLocalD(player.position) - st.takeoff).length();
        c.note(TextFormat("jump: landed %.3f m from takeoff after %.2f s (body frame)", land, c.T - st.t0));
        c.metric("jump landing from takeoff", land, "m");
        // (Reported, not gated: on a 10 deg/s body of r 50 m the spin forces are not small - Coriolis deflects a vertical jump by
        // ~w g t^3 / 6 = ~3.7 m, and the centrifugal push is half of g - so ~2.7 m from takeoff is physics, not drift. The
        // scenario compares the switch off and on: frames/ride against frames/ride-pose.)
        return true;
    });
}
}

SCENARIO("frames/ride") { Build(s); }
SCENARIO("frames/ride-pose") { s.env = { { "RAYLIBGAME_FRAMES", "1" } }; Build(s); }
