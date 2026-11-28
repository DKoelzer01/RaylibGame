// frames/net-correction (#1855 M3, M4; Plan 11 F2b): the player falls from 60 m over the start. One
// second into the fall its net state is taken (as a server sends it), given 5 m/s more upward (a server's correction) and applied.
// Checked: right after, the velocity against the ground is the old one plus the 5 m/s (M3: the frame path's airborne velocity,
// against the ground already, had the ground's motion taken off it again - ~105 m/s near a spinning planet), and one frame
// later the correction still holds, less gravity (M4: the next step rebuilt the old flight from the frame path's stored velocity).
#include "scenario.h"
#include "player.h"
#include "world.h"
#include "frames.h"
#include <raymath.h>
#include <cmath>

namespace {
struct State {
    Vector3 dir{}, before{}, applied{}, after{};
    DVec3 placed{};
    double tFall = -1.0, tApply = 0.0;
    int phase = 0; // 0 placed, 1 corrected (wait a frame), 2 done
} S;
}

SCENARIO("frames/net-correction") {
    s.tags = { "frames", "physics", "planet", "multiplayer" };
    s.world = ScenarioWorld::Planet;
    s.step("place 60 m over the ground", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle() || !c.field) return false;
        S.dir = c.spawnDirWorld();
        const float r = c.field->radius + c.field->elevation(S.dir) + 60.0f;
        S.placed = DVec3(S.dir.x * r, S.dir.y * r, S.dir.z * r);
        player.flyMode = false; // (scenarios start flying: the fly branch returns before the fall and the frame path)
        player.placeAt(c.world, *c.body, (c.body->toWorldD(S.placed) - c.world.renderOrigin).toFloat());
        return true;
    });
    s.step("fall, correct, one frame on", 60.0, [](ScenarioCtx& c) {
        const Vector3 up = c.body->dirToWorld(S.dir);
        // (No c.camera: it flies the player to the eye. The player's own view shows the fall, as in frames/long-fall.)
        if (S.phase == 0) {
            if (player.grounded) return false;
            c.waitingFor("the spawn hold to let the player fall");
            if (S.tFall < 0.0 && (c.body->toLocalD(player.position) - S.placed).length() > 0.01) S.tFall = c.t; // (the hold keeps it still while the place loads)
            if (S.tFall < 0.0) return false;
            if (c.t - S.tFall < 1.0) return false;
            S.before = player.groundVelocity(c.world);
            // The frame path needs the player in the body's Fixed frame (else nothing here is on the path).
            const bool onPath = player.frame.frame.kind == FrameKind::Fixed && player.frame.frame.body == c.body->netId;
            c.check("on the frame path (the body's Fixed frame)", onPath,
                    TextFormat("kind %d, frame body %llu, the body %llu", (int)player.frame.frame.kind,
                               (unsigned long long)player.frame.frame.body, (unsigned long long)c.body->netId));
            Player::NetState ns;
            if (!player.netState(c.world, ns)) { c.check("the player has a net state", false, "netState false"); return true; }
            ns.vel = Vector3Add(ns.vel, Vector3Scale(S.dir, 5.0f)); // body-local: 5 m/s more up (the server's correction)
            c.check("the corrected state applies", player.applyNetState(c.world, ns), "applyNetState false");
            S.applied = player.groundVelocity(c.world);
            S.tApply = c.t;
            S.phase = 1;
            return false;
        }
        if (S.phase == 1) { S.phase = 2; return false; } // (one frame steps with the corrected state)
        S.after = player.groundVelocity(c.world);
        const float upBefore = Vector3DotProduct(S.before, up), upApplied = Vector3DotProduct(S.applied, up), upAfter = Vector3DotProduct(S.after, up);
        const float sideErr = Vector3Length(Vector3Subtract(Vector3Subtract(S.applied, Vector3Scale(up, upApplied)), Vector3Subtract(S.before, Vector3Scale(up, upBefore))));
        c.note(TextFormat("against the ground: up %.2f m/s before, %.2f applied, %.2f a frame later; sideways change at the apply %.2f m/s",
                          upBefore, upApplied, upAfter, sideErr));
        c.between("the applied state is the old velocity + 5 m/s up (m/s up, M3)", upApplied - upBefore, 4.9, 5.1);
        c.lt("no sideways jump at the apply (m/s, M3: ~the ground's speed when the ground's motion was taken off twice)", sideErr, 0.1);
        // Gravity (~9.8 m/s^2 at the surface) acts over the frames since the apply; without M4 the old flight came back (~0).
        const double expect = 5.0 - 9.8 * (c.t - S.tApply);
        c.note(TextFormat("%.3f s after the apply: expected %.2f m/s over before (the correction less gravity)", c.t - S.tApply, expect));
        c.between("frames later the correction holds (m/s up against before, M4: within 0.3 of the correction less gravity)", upAfter - upBefore, expect - 0.3, expect + 0.3);
        return true;
    });
}
