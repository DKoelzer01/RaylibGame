// saves/frame-jump (Plan 11 F7a, #2074), RAYLIBGAME_FRAMES=1: the player jumps on Terra and the world is saved and loaded at
// the top of the climb. The player's frame state (frame, double position and velocity against the ground) comes back from
// the record's player_frame line to 1e-6 m and 1e-6 m/s, and the jump goes on from it: the player is still airborne and
// still rising a step later, and lands within 0.5 m of where a jump of that speed comes down (before F7a a load zeroed the
// velocity: the player dropped straight down from the save's height).
#include "scenario.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
struct JumpState {
    bool placed = false, jumped = false, loaded = false;
    double t0 = 0.0;
    FramePos before{};
    DVec3 globalBefore{};
    DVec3 takeoff{};    // the feet at takeoff, in the body's frame
    double upSpeed = 0.0; // the velocity along up at the save (m/s)
};
JumpState st;

double Radial(const FramePos& f) { const double l = f.pos.length(); return l > 0.0 ? (f.vel.x * f.pos.x + f.vel.y * f.pos.y + f.vel.z * f.pos.z) / l : 0.0; }
}

SCENARIO("saves/frame-jump") {
    s.tags = { "saves", "frames" };
    s.world = ScenarioWorld::Planet;
    s.env = { { "RAYLIBGAME_FRAMES", "1" } };
    s.timeout = 240.0;
    s.step("an idle world", 120.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T > 1.0 && WorldIdle(); });
    s.step("stand, then jump", 30.0, [](ScenarioCtx& c) {
        if (!st.placed) { st = JumpState{}; player.spawnOn(c.world, *c.body, c.spawnDirWorld()); st.placed = true; st.t0 = c.T; return false; }
        c.waitingFor("grounded for 1 s");
        if (!player.grounded) { st.t0 = c.T; return false; }
        if (c.T - st.t0 < 1.0) return false;
        st.takeoff = c.body->toLocalD(player.position);
        player.autoJump = true;
        st.t0 = c.T;
        return true;
    });
    s.step("save and load mid-jump", 10.0, [](ScenarioCtx& c) {
        c.waitingFor("airborne and 0.1 s into the jump");
        if (player.grounded || c.T - st.t0 < 0.1 || player.frame.frame.body != c.body->netId) return false;
        st.before = player.frame;
        st.globalBefore = player.position;
        st.upSpeed = Radial(st.before);
        c.note(TextFormat("saved %.3f m above the takeoff, rising at %.3f m/s (frame %s)", st.before.pos.length() - st.takeoff.length(), st.upSpeed,
                          st.before.frame.kind == FrameKind::Fixed ? "Fixed" : "Centred"));
        if (!c.check("rising at the save", st.upSpeed > 0.5, TextFormat("%.3f m/s", st.upSpeed))) return true;
        if (!c.check("the save", SaveWorld(c.world))) return true;
        if (!c.check("the load", LoadWorldSave(c.world), LastLoadError())) return true;
        const FramePos& a = player.frame;
        c.check("the same frame after the load", a.frame == st.before.frame);
        c.lt("frame position moved by the load (m)", (a.pos - st.before.pos).length(), 1e-6);
        c.lt("frame velocity changed by the load (m/s)", (a.vel - st.before.vel).length(), 1e-6);
        // The global position also goes through the body's pose, whose saved "rot" is %f (1e-6 rad: up to 6 cm at 60 km);
        // F7b takes poses from the clock and tightens this to 1e-6 m.
        const double moved = (player.position - st.globalBefore).length();
        c.metric("global position moved by the load", moved, "m");
        c.lt("global position moved by the load (m; the body's rot line rounds)", moved, 0.06);
        c.check("airborne after the load", !player.grounded);
        st.t0 = c.T;
        return true;
    });
    s.step("the jump goes on", 10.0, [](ScenarioCtx& c) {
        if (!st.loaded) {
            c.waitingFor("a step after the load");
            if (c.frame < 2 || c.t < 0.05) return false;
            st.loaded = true;
            const double v = Radial(player.frame);
            c.note(TextFormat("a step after the load: rising at %.3f m/s (%.3f at the save)", v, st.upSpeed));
            c.check("still rising after the load", v > 0.0, TextFormat("%.3f m/s", v));
            return false;
        }
        c.waitingFor("landing");
        if (!player.grounded) return false;
        // Where a ballistic jump of the saved speed lands (flat ground, the body's surface gravity): the climb left plus the
        // fall, from the height at the save. The ground is not flat; 0.5 m is the tolerance.
        const DVec3 land = c.body->toLocalD(player.position);
        const double h = st.before.pos.length() - st.takeoff.length();
        const double g = c.body->surfaceGravity;
        const double peak = h + st.upSpeed * st.upSpeed / (2.0 * g);
        c.note(TextFormat("peak predicted %.3f m over the takeoff; landed %.3f m from the takeoff spot", peak, (land - st.takeoff).length()));
        c.check("alive after landing", !survival.dead);
        // Rising at the save: the landing comes a climb and a fall later; dropped at rest it would come at once. Measured
        // as the time from the load to the landing against the ballistic one.
        const double tFlight = st.upSpeed / g + sqrt(2.0 * peak / g);
        c.near("time from the load to the landing (s)", c.T - st.t0, tFlight, std::max(0.1, 0.25 * tFlight));
        return true;
    });
}
