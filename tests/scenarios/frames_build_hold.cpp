// frames/build-hold (Plan 11 F3, #1277; the #551 "placement boxes jitter" class): the player stands still on Terra with the
// build tool out, looking at the ground 3 m ahead, for 10 s of sim time while the planet turns under them (3 deg a minute,
// ~1.7 m/s of ground speed at the start). The camera follows the player the game's own way (applyToCamera, every frame) and
// the build tool aims from it every frame. Nothing moves in the body's frame, so the aimed placement must not either:
// counted are the frames whose aimed cell (or new plot origin) differs from the first, and how far the new origin moved.
// The recording keeps the box in view (the camera is the player's).
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
#include <algorithm>

namespace {
constexpr double kHoldSeconds = 10.0;
struct HoldState {
    bool placed = false, settled = false, aimed = false;
    double t0 = 0.0;
    const Plot* plot = nullptr; Int3 origin{}; Vector3 newOrigin{};
    int frames = 0, changed = 0, none = 0, cellFlips = 0;
    double maxMove = 0.0; // the new plot origin's largest distance from the first (body-local m)
    // What the aim is made from, in the body's frame and double: how far each wandered from its first value while still.
    DVec3 feet0{}, cam0{}, aim0{}; Vector3 fwd0{};
    double feetWob = 0.0, camWob = 0.0, fwdWob = 0.0, aimWob = 0.0; // m, m, deg, m
    double maxGrains = 0.0; // the new plot origin's largest move on one axis, in float steps of its largest coordinate
};
HoldState st;

// The float grain of a body-local point: the spacing of floats at its largest coordinate (3.9 mm at 32-64 km, 1.95 mm at
// 16-32 km). The transforms that make the point round all three axes at about that size, so a still point flips by up to one
// grain on any axis - a small coordinate's own (much finer) ulps say nothing about it.
double Grain(Vector3 p) {
    const float m = std::max({ fabsf(p.x), fabsf(p.y), fabsf(p.z) });
    return (double)nextafterf(m, INFINITY) - (double)m;
}

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
        DVec3 aimL; Vector3 aimDirL, aimDirW; player.aimRay(c.world, b, aimL, aimDirL, aimDirW); // (F4 #1412: what the tools aim from)
        if (st.frames == 0) { st.feet0 = feetL; st.cam0 = camL; st.fwd0 = fwdL; st.aim0 = aimL; }
        st.aimWob = std::max(st.aimWob, (aimL - st.aim0).length());
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
            // A visible flip: another plot, or another cell of the same one (what #551 looked like).
            if (plot != st.plot || (plot && (origin.x != st.origin.x || origin.y != st.origin.y || origin.z != st.origin.z))) ++st.cellFlips;
            if (!plot) {
                st.maxMove = std::max(st.maxMove, (double)Vector3Distance(newOrigin, st.newOrigin));
                const Vector3 d = Vector3Subtract(newOrigin, st.newOrigin);
                st.maxGrains = std::max(st.maxGrains, std::max({ fabsf(d.x), fabsf(d.y), fabsf(d.z) }) / Grain(st.newOrigin));
            }
        }
        c.waitingFor("10 s holding still");
        if (c.T - st.t0 < kHoldSeconds) return false;
        c.note(TextFormat("%d frames held still: the aim changed in %d (no aim in %d); the new plot origin moved up to %.2f mm",
                          st.frames, st.changed, st.none, st.maxMove * 1000.0));
        c.note(TextFormat("while still, in the body's frame: feet wandered %.3f mm, the camera %.3f mm, the view %.4f deg",
                          st.feetWob * 1000.0, st.camWob * 1000.0, st.fwdWob));
        c.metric("feet wobble", st.feetWob * 1000.0, "mm");
        c.metric("camera wobble", st.camWob * 1000.0, "mm");
        c.metric("aim origin wobble", st.aimWob * 1000.0, "mm");
        c.note(TextFormat("the aim's origin (Player::aimRay) wandered %.4f mm", st.aimWob * 1000.0));
        c.metric("view wobble", st.fwdWob, "deg");
        c.metric("aim changes while still", st.changed, "frames");
        c.metric("new plot origin movement", st.maxMove * 1000.0, "mm");
        c.check("the build tool aimed at the ground", st.aimed);
        // (#1381, re-scoped) gated on what shows: no flip to another plot or cell, and a new plot's ghost box moving no more
        // than one float grain (Grain) on any axis. Bit-exact equality was too strict: a still camera wobbles ~0.06 mm (eye
        // smoothing, interpolation; F4 #1412), and rounding a wobbling value flips between neighbouring floats now and then.
        // (#1759) A 5 mm bound on the distance failed ~1 run in 10: one grain (3.9 mm) on each of two axes is 5.5 mm.
        c.lt("plot or cell flips while still", st.cellFlips, 1);
        c.metric("new plot origin movement (float grains, largest axis)", st.maxGrains, "grains");
        c.note(TextFormat("the new plot origin's grain: %.2f mm; it moved up to %.2f grains on an axis", Grain(st.newOrigin) * 1000.0, st.maxGrains));
        c.lt("a new plot's origin moved more than one float grain on an axis while still", st.maxGrains, 1.01);
        buildTool.active = false;
        return true;
    });
}
}

SCENARIO("frames/build-hold") { Build(s); } // (F8b #2113: frames/build-hold-pose folded in)
