// frames/still (Plan 11 F0, #667): a player standing still for 3 minutes of sim time, at the start's latitude and on the
// equator, while the planet turns under it (Terra: 3 deg a minute). Measured in the body's own frame: how far the feet
// moved and how far the view turned. Standing still must stay still on the ground (today's ground carry; Plan 11: v = 0
// in the body-fixed frame). F0 records today's numbers; the bound is loose (5 cm, 0.5 deg) and tightens at F2.
#include "scenario.h"
#include "player.h"
#include "world.h"
#include "worldtypes.h"
#include "voxel.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
constexpr double kStandSeconds = 180.0;
struct StillState { bool placed = false, settled = false; double t0 = 0.0; DVec3 feet{}; Vector3 fwd{}, up{}; double maxStep = 0.0; DVec3 last{}; };
StillState st;

Vector3 Equator(const Planetoid& b, Vector3 dir) {
    const Vector3 s = b.spin;
    if (Vector3Length(s) <= 0.0f) return dir;
    const Vector3 axis = b.dirToLocal(Vector3Normalize(s));
    Vector3 d = Vector3Subtract(dir, Vector3Scale(axis, Vector3DotProduct(dir, axis)));
    d = Vector3Normalize(d);
    for (int k = 0; k < 360; ++k) { // a dry spot along the equator
        const Vector3 e = Vector3RotateByAxisAngle(d, axis, k * 0.25f * DEG2RAD);
        if (const auto* pf = dynamic_cast<const PlanetField*>(b.field.get()); pf && pf->elevation(e) > 3.0f) return e;
    }
    return d;
}

void StandStep(ScenarioDef& s, const char* name, bool equator, double maxStepBound) {
    s.step(std::string("stand still ") + name, kStandSeconds + 90.0, [name, equator, maxStepBound](ScenarioCtx& c) {
        Planetoid& b = *c.body;
        if (!st.placed) {
            st = StillState{};
            const Vector3 dir = equator ? Equator(b, c.spawnDirLocal()) : c.spawnDirLocal();
            player.spawnOn(c.world, b, b.dirToWorld(dir));
            st.placed = true;
            st.t0 = c.T;
            return false;
        }
        if (!st.settled) { // grounded, and 3 s for the hold and the landing to pass
            c.waitingFor("grounded at the spot");
            if (!player.grounded) { st.t0 = c.T; return false; }
            if (c.T - st.t0 < 3.0) return false;
            st.settled = true;
            st.t0 = c.T;
            st.feet = st.last = b.toLocalD(player.position);
            st.fwd = b.dirToLocal(player.forward);
            st.up = b.dirToLocal(player.up);
            return false;
        }
        const DVec3 feet = b.toLocalD(player.position);
        st.maxStep = std::max(st.maxStep, (feet - st.last).length());
        st.last = feet;
        c.waitingFor("the stand to end");
        if (c.T - st.t0 < kStandSeconds) return false;
        const double drift = (feet - st.feet).length();
        const Vector3 fwd = b.dirToLocal(player.forward), up = b.dirToLocal(player.up);
        const double turn = acos(std::clamp((double)Vector3DotProduct(fwd, st.fwd), -1.0, 1.0)) * 180.0 / PI;
        const double tilt = acos(std::clamp((double)Vector3DotProduct(up, st.up), -1.0, 1.0)) * 180.0 / PI;
        const std::string n = name;
        c.note(n + TextFormat(": %.0f s standing: feet moved %.4f m (largest step %.4f m), view turned %.4f deg, up tilted %.4f deg, grounded %d",
                              kStandSeconds, drift, st.maxStep, turn, tilt, (int)player.grounded));
        c.metric(n + " feet drift", drift, "m");
        c.metric(n + " largest step", st.maxStep, "m");
        c.metric(n + " view turn", turn, "deg");
        c.metric(n + " up tilt", tilt, "deg");
        c.lt(n + " feet drift", drift, 0.05);
        if (maxStepBound > 0.0) c.lt(n + " largest single step", st.maxStep, maxStepBound);
        c.lt(n + " view turn", turn, 0.5);
        c.check(n + " still grounded", player.grounded);
        if (Frames_On()) { // F2a (#1146): standing still is v = 0 in the Fixed frame of the body underfoot
            const double v = player.frame.vel.length();
            c.metric(n + " frame speed", v, "m/s");
            c.check(n + " frame: Fixed on the body underfoot", player.frame.frame.kind == FrameKind::Fixed && player.frame.frame.body == b.netId);
            c.lt(n + " frame speed standing", v, 0.01);
        }
        st.placed = false;
        return true;
    });
}
}

static void Build(ScenarioDef& s, double maxStepBound) {
    s.tags = { "frames", "physics", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 900.0;
    s.record = false;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        return c.T > 1.0 && WorldIdle();
    });
    StandStep(s, "start", false, maxStepBound);
    StandStep(s, "equator", true, maxStepBound);
}

// Today's frames (the float Euler body angle): single steps are only reported (13-31 mm at F0).
SCENARIO("frames/still") { Build(s, 0.0); }

// Plan 11 F1 (#667): RAYLIBGAME_FRAMES=1 - body poses from the double spin angle; no single step over 1 mm.
SCENARIO("frames/still-pose") {
    s.env = { { "RAYLIBGAME_FRAMES", "1" } };
    Build(s, 0.001);
}
