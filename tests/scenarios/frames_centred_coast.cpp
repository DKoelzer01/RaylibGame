// frames/centred-coast (Plan 11 F8b-2, #2113): the player airborne above Terra's altitude line flies in its Centred frame
// (the velocity against the centre, carried with the centre only, felt gravity) and comes back down into the Fixed frame.
// Placed 37.5 km up (over the 36.75 km Fixed -> Centred line), pushed 300 m/s down once Centred, followed through the
// Centred -> Fixed switch at 33.25 km. Gates: the path against an RK4 reference in the Centred frame, the global velocity
// continuous from update to update (no frame-conversion kick), and both frame kinds seen.
#include "scenario.h"
#include "frames.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace {
struct CoastState {
    bool placed = false, released = false, pushed = false;
    DVec3 startLocal{};
    double groundR = 0.0, lastT = 0.0, worstPos = 0.0, worstVel = 0.0, worstKick = 0.0, worstKickAt = 0.0;
    int centredUpdates = 0, fixedUpdates = 0, switches = 0;
    FrameKind lastKind = FrameKind::Fixed;
    FramePos ref{};   // the reference, in Terra's Centred frame
    DVec3 lastVg{};   // the player's global velocity at the last callback
    bool haveLast = false;
};
CoastState st;

DVec3 AccelCentred(const FramePos& f) { return Frames_Acceleration(Planetoid::All(), f); }

// RK4 in the Centred frame (felt gravity only: no spin forces there), over dt in steps of at most 1/480 s.
void Integrate(FramePos& f, double dt) {
    const int n = std::max(1, (int)std::ceil(dt * 480.0));
    const double h = dt / n;
    for (int i = 0; i < n; ++i) {
        auto at = [&f](const DVec3& p, const DVec3& v) { FramePos q = f; q.pos = p; q.vel = v; return AccelCentred(q); };
        const DVec3 p0 = f.pos, v0 = f.vel;
        const DVec3 a1 = at(p0, v0);
        const DVec3 a2 = at(p0 + v0 * (h / 2), v0 + a1 * (h / 2));
        const DVec3 a3 = at(p0 + (v0 + a1 * (h / 2)) * (h / 2), v0 + a2 * (h / 2));
        const DVec3 a4 = at(p0 + (v0 + a2 * (h / 2)) * h, v0 + a3 * h);
        f.pos = p0 + v0 * h + (a1 + a2 + a3) * (h * h / 6); // (RK4 for x'' = a: h^2/6 of the first three stages)
        f.vel = v0 + (a1 + a2 * 2.0 + a3 * 2.0 + a4) * (h / 6);
    }
}
}

SCENARIO("frames/centred-coast") {
    s.tags = { "frames", "physics", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 300.0;
    s.record = false;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        return c.T > 1.0 && WorldIdle();
    });
    s.step("coast down through the altitude line", 90.0, [](ScenarioCtx& c) {
        Planetoid& b = *c.body;
        survival.dead = false; survival.health = 100.0f; // (no air up there: the gate is the flight, not the survival)
        if (!st.placed) {
            st = CoastState{};
            const Vector3 dir = c.spawnDirLocal();
            st.groundR = c.field->radius + c.field->elevation(dir);
            st.startLocal = DVec3(dir.x, dir.y, dir.z) * (st.groundR + 37500.0);
            player.placeAt(c.world, b, (b.toWorldD(st.startLocal) - c.world.renderOrigin).toFloat());
            st.placed = true;
            return false;
        }
        // The spawn hold keeps the player co-rotating at the drop point while the place loads: the fall starts when it moves.
        if (!st.released) {
            c.waitingFor("the hold to release the player");
            if ((b.toLocalD(player.position) - st.startLocal).length() < 0.01) return false;
            st.released = true;
        }
        const FrameRef centred{ b.netId, FrameKind::Centred };
        if (!st.pushed) {
            c.waitingFor("the Centred frame");
            if (player.frame.frame != centred) return false;
            // Straight down at 300 m/s, against the centre: set in the frame's own terms (Centred axes are the world's).
            const DVec3 down = player.frame.pos * (-1.0 / player.frame.pos.length());
            player.velocity = Vector3Add(player.velocity, (down * 300.0).toFloat());
            player.frame.vel = player.frame.vel + down * 300.0;
            st.ref = player.frame;
            st.pushed = true;
            st.lastT = c.T;
            return false;
        }
        // The player's state this callback, in the reference's frame and globally.
        FramePos now;
        DVec3 pg, vg;
        if (!Frames_Convert(player.frame, centred, now) || !Frames_ToGlobal(player.frame, pg, vg)) { c.fail("no frame state"); return true; }
        const double dt = c.T - st.lastT;
        if (dt > 0.0) {
            const DVec3 vBefore = st.ref.vel;
            Integrate(st.ref, dt);
            // A kick: the global velocity's change against the reference's (gravity over dt); a conversion done twice or
            // missed is the centre's speed or the ground's (w x r ~ 170 m/s here).
            if (st.haveLast) {
                const double kick = ((vg - st.lastVg) - (st.ref.vel - vBefore)).length();
                if (kick > st.worstKick) { st.worstKick = kick; st.worstKickAt = (now.pos.length() - st.groundR) / 1000.0; }
            }
            st.worstPos = std::max(st.worstPos, (now.pos - st.ref.pos).length());
            st.worstVel = std::max(st.worstVel, (now.vel - st.ref.vel).length());
        }
        st.lastVg = vg; st.haveLast = dt > 0.0 || st.haveLast; st.lastT = c.T;
        const FrameKind k = player.frame.frame.kind;
        if (k != st.lastKind && st.centredUpdates + st.fixedUpdates > 0) ++st.switches;
        (k == FrameKind::Centred ? st.centredUpdates : st.fixedUpdates)++;
        st.lastKind = k;
        const double alt = now.pos.length() - st.groundR;
        c.waitingFor(TextFormat("the Fixed frame below 33.25 km (at %.2f km, %s)", alt / 1000.0, k == FrameKind::Centred ? "Centred" : "Fixed"));
        if (player.grounded) { c.fail("landed before the test ended"); return true; }
        if (!(k == FrameKind::Fixed && alt < 32000.0)) return false;
        c.note(TextFormat("Centred for %d callbacks, then Fixed for %d (%d switch); worst position %.3f m and velocity %.4f m/s off the "
                          "RK4 reference; worst velocity kick %.4f m/s at %.2f km",
                          st.centredUpdates, st.fixedUpdates, st.switches, st.worstPos, st.worstVel, st.worstKick, st.worstKickAt));
        c.metric("worst position off the reference", st.worstPos, "m");
        c.metric("worst velocity kick", st.worstKick, "m/s");
        c.gt("callbacks in the Centred frame", st.centredUpdates, 10);
        c.gt("callbacks in the Fixed frame", st.fixedUpdates, 0);
        c.eq("frame switches", st.switches, 1);
        // Semi-implicit Euler at 1/120 s against RK4: ~g h t / 2 = 0.2 m over the ~12 s. A wrong conversion is ~170 m/s.
        c.lt("position off the RK4 reference (m)", st.worstPos, 1.0);
        c.lt("velocity off the RK4 reference (m/s)", st.worstVel, 0.1);
        c.lt("global velocity kick between updates (m/s)", st.worstKick, 0.1);
        st.placed = false;
        return true;
    });
    // (F8b-2) A placeAt mid-flight: the place's velocity (at rest on the turning ground), not the last flight's. placeAt
    // anchors the player, so before the spawn count was checked the Fixed path took this for the same flight and rebuilt the
    // stored 300+ m/s fall.
    s.step("teleport mid-flight", 60.0, [](ScenarioCtx& c) {
        Planetoid& b = *c.body;
        survival.dead = false; survival.health = 100.0f;
        static int spawns0 = 0, after = 0;
        static double fastest = 0.0;
        if (!st.placed) {
            const double v = Vector3Length(player.groundVelocity(c.world));
            c.note(TextFormat("before the teleport: %.1f m/s against the ground, airborne %s", v, player.grounded ? "no" : "yes"));
            c.gt("falling fast before the teleport (m/s)", v, 100.0);
            const DVec3 dir = st.startLocal * (1.0 / st.startLocal.length());
            st.startLocal = dir * (st.groundR + 2000.0);
            player.placeAt(c.world, b, (b.toWorldD(st.startLocal) - c.world.renderOrigin).toFloat());
            spawns0 = player.spawns; after = 0; fastest = 0.0;
            st.placed = true;
            return false;
        }
        if (player.spawns != spawns0) { c.fail("placed again (a respawn) during the check"); return true; }
        fastest = std::max(fastest, (double)Vector3Length(player.groundVelocity(c.world)));
        c.waitingFor("a second after the teleport");
        if (++after < 2 || c.t < 1.0) return false;
        // From rest a second's fall at ~9 m/s^2 is ~10 m/s; the old flight would be 300+.
        c.lt("fastest against the ground in the second after the teleport (m/s)", fastest, 30.0);
        return true;
    });
}
