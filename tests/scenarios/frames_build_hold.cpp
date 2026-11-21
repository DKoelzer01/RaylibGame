// frames/build-hold (Plan 11 F3, #1277; the #551 "placement boxes jitter" class): the player stands still on Terra with the
// build tool out, looking at the ground 3 m ahead, for 10 s of sim time while the planet turns under them (3 deg a minute,
// ~1.7 m/s of ground speed at the start). The camera follows the player the game's own way (applyToCamera, every frame) and
// the build tool aims from it every frame. Nothing moves in the body's frame, so the aimed placement must not either:
// counted are the frames whose aimed cell (or new plot origin) differs from the first, and how far the new origin moved.
// Reported, not gated yet (no baseline); frames/build-hold-pose runs the same with RAYLIBGAME_FRAMES=1. The recording
// keeps the box in view (the camera is the player's).
#include "scenario.h"
#include "buildtool.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <cmath>
#include <cstring>

namespace {
constexpr double kHoldSeconds = 10.0;
struct HoldState {
    bool placed = false, settled = false, aimed = false;
    double t0 = 0.0;
    const Plot* plot = nullptr; Int3 origin{}; Vector3 newOrigin{};
    int frames = 0, changed = 0, none = 0;
    double maxMove = 0.0; // the new plot origin's largest distance from the first (body-local m)
    // What the aim is made from, in the body's frame and double: how far each wandered from its first value while still.
    DVec3 feet0{}, cam0{}; Vector3 fwd0{};
    double feetWob = 0.0, camWob = 0.0, fwdWob = 0.0; // m, m, deg
};
HoldState st;

void Build(ScenarioDef& s) {
    s.tags = { "frames", "build" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 180.0;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        return c.T > 1.0 && WorldIdle();
    });
    s.step("hold the build tool still while the planet turns", 60.0, [](ScenarioCtx& c) {
        Planetoid& b = *c.body;
        if (!st.placed) {
            st = HoldState{};
            player.spawnOn(c.world, b, c.spawnDirWorld());
            st.placed = true;
            st.t0 = c.T;
            return false;
        }
        if (!st.settled) {
            c.waitingFor("grounded and 2 s settled");
            if (!player.grounded) { st.t0 = c.T; return false; }
            if (c.T - st.t0 < 2.0) return false;
            // Look at the ground 3 m ahead (the view is carried with the body from here on, as in play).
            const Vector3 up = Vector3Normalize(player.up);
            Vector3 flat = Vector3Subtract(player.forward, Vector3Scale(up, Vector3DotProduct(player.forward, up)));
            flat = Vector3Length(flat) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(flat);
            player.forward = Vector3Normalize(Vector3Subtract(Vector3Scale(flat, 3.0f), Vector3Scale(up, player.eyeHeight)));
            buildTool.active = true; buildTool.orientManual = false;
            buildTool.part = (int)PartType::Cube;
            player.nearest = &b;
            st.settled = true;
            st.t0 = c.T;
            return false;
        }
        player.applyToCamera(c.world); // the camera at the player's eye, as every frame of play
        BuildTool_Update(c.world, PlayerInput::None(), true);
        Vector3 newOrigin{}; Int3 origin{};
        const Plot* plot = BuildTool_AimPlacement(&newOrigin, &origin);
        const bool any = plot || newOrigin.x != 0.0f || newOrigin.y != 0.0f || newOrigin.z != 0.0f;
        const DVec3 feetL = b.toLocalD(player.position);
        const DVec3 camL = b.toLocalD(c.world.renderOrigin + DVec3(c.world.camera.position));
        const Vector3 fwdL = b.dirToLocal(Vector3Normalize(Vector3Subtract(c.world.camera.target, c.world.camera.position)));
        if (st.frames == 0) { st.feet0 = feetL; st.cam0 = camL; st.fwd0 = fwdL; }
        st.feetWob = std::max(st.feetWob, (feetL - st.feet0).length());
        st.camWob = std::max(st.camWob, (camL - st.cam0).length());
        st.fwdWob = std::max(st.fwdWob, (double)Vector3Angle(fwdL, st.fwd0) * RAD2DEG);
        ++st.frames;
        if (!any) ++st.none;
        else if (!st.aimed) { st.aimed = true; st.plot = plot; st.origin = origin; st.newOrigin = newOrigin; }
        else {
            const bool same = plot == st.plot && origin.x == st.origin.x && origin.y == st.origin.y && origin.z == st.origin.z &&
                              std::memcmp(&newOrigin, &st.newOrigin, sizeof newOrigin) == 0;
            if (!same) ++st.changed;
            if (!plot) st.maxMove = std::max(st.maxMove, (double)Vector3Distance(newOrigin, st.newOrigin));
        }
        c.waitingFor("10 s holding still");
        if (c.T - st.t0 < kHoldSeconds) return false;
        c.note(TextFormat("%d frames held still: the aim changed in %d (no aim in %d); the new plot origin moved up to %.2f mm (%s)",
                          st.frames, st.changed, st.none, st.maxMove * 1000.0, Frames_On() ? "frames on" : "frames off"));
        c.note(TextFormat("while still, in the body's frame: feet wandered %.3f mm, the camera %.3f mm, the view %.4f deg",
                          st.feetWob * 1000.0, st.camWob * 1000.0, st.fwdWob));
        c.metric("feet wobble", st.feetWob * 1000.0, "mm");
        c.metric("camera wobble", st.camWob * 1000.0, "mm");
        c.metric("view wobble", st.fwdWob, "deg");
        c.metric("aim changes while still", st.changed, "frames");
        c.metric("new plot origin movement", st.maxMove * 1000.0, "mm");
        c.check("the build tool aimed at the ground", st.aimed);
        // With the frame path (F3 plots, #1381's new plot origin in double) a still player's aim holds exactly; the old path
        // flips by float steps (hundreds of frames on the runner) and is reported only.
        if (Frames_On()) c.lt("aim changes while still (frame path)", st.changed, 1);
        buildTool.active = false;
        return true;
    });
}
}

SCENARIO("frames/build-hold") { Build(s); }
SCENARIO("frames/build-hold-pose") { s.env = { { "RAYLIBGAME_FRAMES", "1" } }; Build(s); }
