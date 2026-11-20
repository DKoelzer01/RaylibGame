// frames/long-fall (Plan 11 F0, #667): drops from rest (co-rotating with the ground) of 100 m and 1 km, at the start's
// latitude and on the equator. Measured at touchdown: the landing speed against the inverse-square free fall, and the
// eastward drift (body-local, along the ground) against the Coriolis deflection of a dropped body, w g t^3 cos(lat) / 3.
// F0 records today's numbers; Plan 11's frames must keep them (the user's decision: spin forces for anything airborne).
#include "scenario.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include "worldtypes.h"
#include "voxel.h"
#include <raymath.h>
#include <cmath>
#include <cstdio>
#include <string>

namespace {
double Dot(const DVec3& a, const DVec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
struct Drop { const char* name; bool equator; double height; };
const Drop kDrops[] = { { "start-100m", false, 100.0 }, { "start-1km", false, 1000.0 }, { "equator-100m", true, 100.0 }, { "equator-1km", true, 1000.0 } };

struct DropState {
    bool started = false, falling = false;
    DVec3 startLocal{};
    double t0 = 0.0, groundR = 0.0, lat = 0.0, g0 = 0.0, maxLanding = 0.0, fallSpeed = 0.0;
};
DropState st;

// The body's spin axis in its own frame (unit) and its rate (rad/s): spin is degrees per second about the Euler axes.
DVec3 SpinAxisLocal(const Planetoid& b, double& omega) {
    const Vector3 s = b.spin;
    omega = std::sqrt((double)s.x * s.x + (double)s.y * s.y + (double)s.z * s.z) * PI / 180.0;
    if (omega <= 0.0) return DVec3(0, 1, 0);
    // An axis fixed in the body is the same in its own frame whatever the angle: the world axis turned back.
    const Vector3 wa = Vector3Normalize(s);
    const Vector3 la = b.dirToLocal(wa);
    return DVec3(la.x, la.y, la.z);
}

// The real ground in a body-local direction (the density surface, searched down from 60 m over the height map), radius.
double GroundRadius(const Planetoid& b, const PlanetField& f, Vector3 dir) {
    const double r0 = f.radius + f.elevation(dir);
    for (double r = r0 + 60.0; r > r0 - 60.0; r -= 0.05)
        if (b.densityAt(Vector3Scale(dir, (float)r)) >= ISO_LEVEL) return r;
    return r0;
}

// Land within a few degrees of a direction (not the sea): the first dry spot stepping east along the parallel.
Vector3 DryNear(const PlanetField& f, Vector3 dir, const DVec3& axis) {
    const Vector3 ax = { (float)axis.x, (float)axis.y, (float)axis.z };
    for (int k = 0; k < 360; ++k) {
        const Vector3 d = Vector3RotateByAxisAngle(dir, ax, k * 0.25f * DEG2RAD);
        if (f.elevation(d) > 3.0f) return d;
    }
    return dir;
}
}

static void BuildLongFall(ScenarioDef& s) {
    s.tags = { "frames", "physics", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 420.0;
    s.record = false;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        return c.T > 1.0 && WorldIdle();
    });
    for (const Drop& d : kDrops) {
        s.step(std::string("drop ") + d.name, 90.0, [d](ScenarioCtx& c) {
            Planetoid& b = *c.body;
            double omega = 0.0;
            const DVec3 axis = SpinAxisLocal(b, omega);
            if (!st.started) {
                st = DropState{};
                Vector3 dir = c.spawnDirLocal();
                if (d.equator) { // the start's meridian, down to the equator
                    const DVec3 dd(dir.x, dir.y, dir.z);
                    const DVec3 eq = dd - axis * Dot(dd, axis);
                    dir = Vector3Normalize({ (float)eq.x, (float)eq.y, (float)eq.z });
                }
                dir = DryNear(*c.field, dir, axis);
                st.groundR = GroundRadius(b, *c.field, dir);
                st.lat = std::asin(std::clamp(Dot(DVec3(dir.x, dir.y, dir.z), axis), -1.0, 1.0));
                st.startLocal = DVec3(dir.x, dir.y, dir.z) * (st.groundR + d.height);
                const DVec3 groundW = b.toWorldD(DVec3(dir.x, dir.y, dir.z) * st.groundR);
                st.g0 = Vector3Length(BodyGravityAt(b, (groundW - c.world.renderOrigin).toFloat()));
                player.placeAt(c.world, b, (b.toWorldD(st.startLocal) - c.world.renderOrigin).toFloat());
                st.t0 = c.T;
                st.started = true;
                return false;
            }
            // The spawn hold keeps the player at the drop point while the place loads: the fall starts when it moves.
            const double fallen = (b.toLocalD(player.position) - st.startLocal).length();
            if (!st.falling) {
                c.waitingFor("the hold to release the player");
                if (fallen < 0.01) { st.t0 = c.T; return false; }
                st.falling = true;
            }
            c.waitingFor("touchdown");
            if (fallen < 0.5 * d.height) { st.maxLanding = 0.0; return false; } // (a landing flag left from before the drop)
            st.maxLanding = std::max(st.maxLanding, (double)player.landingSpeed); // into the floor (its normal: less on a slope)
            if (!player.grounded) { // the fall speed against the turning ground, the last airborne frame before touchdown
                const Vector3 feet = (player.position - c.world.renderOrigin).toFloat();
                st.fallSpeed = Vector3Length(Vector3Subtract(player.velocity, b.surfaceVelocity(feet)));
            }
            if (!player.grounded) return false;
            const double t = c.T - st.t0;
            const DVec3 land = b.toLocalD(player.position);
            const DVec3 up = st.startLocal * (1.0 / st.startLocal.length());
            DVec3 east(axis.y * up.z - axis.z * up.y, axis.z * up.x - axis.x * up.z, axis.x * up.y - axis.y * up.x); // axis x up
            east = east * (1.0 / std::max(1e-9, east.length()));
            const DVec3 moved = land - st.startLocal;
            const double driftEast = Dot(moved, east), driftSide = (moved - up * Dot(moved, up) - east * driftEast).length();
            const double rTop = st.groundR + d.height;
            const double vAnalytic = std::sqrt(2.0 * st.g0 * st.groundR * st.groundR * (1.0 / st.groundR - 1.0 / rTop));
            const double driftAnalytic = omega * st.g0 * t * t * t * std::cos(st.lat) / 3.0;
            const std::string n = d.name;
            c.note(n + TextFormat(": latitude %.1f deg, fall %.2f s, speed %.2f m/s (free fall %.2f; into the floor %.2f), drift east %.3f m (Coriolis %.3f), sideways %.3f m",
                                  st.lat * 180.0 / PI, t, st.fallSpeed, vAnalytic, st.maxLanding, driftEast, driftAnalytic, driftSide));
            c.metric(n + " landing speed into the floor", st.maxLanding, "m/s", false);
            c.metric(n + " fall speed", st.fallSpeed, "m/s", false);
            c.metric(n + " fall speed / free fall", st.fallSpeed / vAnalytic, "", false);
            c.metric(n + " drift east", driftEast, "m", false);
            c.metric(n + " Coriolis drift", driftAnalytic, "m", false);
            c.metric(n + " drift sideways", driftSide, "m");
            // Baseline gates, loose (F0): the landing speed within 5% of free fall, the drift within 25% (+3 cm) of Coriolis.
            c.near(n + " fall speed vs free fall", st.fallSpeed, vAnalytic, 0.05 * vAnalytic + 0.2);
            c.near(n + " eastward drift vs Coriolis", driftEast, driftAnalytic, 0.25 * std::fabs(driftAnalytic) + 0.03);
            // A hard landing kills (fall damage): revived here, before its respawn could carry the player off mid-drop.
            survival.dead = false; survival.health = 100.0f;
            st.started = false;
            return true;
        });
    }
}

SCENARIO("frames/long-fall") { BuildLongFall(s); }
// Plan 11 F2b (#1146): the same with RAYLIBGAME_FRAMES=1 (the player's velocity in the Fixed frame, Frames_Acceleration's spin forces).
SCENARIO("frames/long-fall-pose") { s.env = { { "RAYLIBGAME_FRAMES", "1" } }; BuildLongFall(s); }
