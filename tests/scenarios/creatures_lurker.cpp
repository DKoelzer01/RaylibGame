// The tunnel lurker (#1224, Plan 12 first wave).
//   creatures/arena-lurker (sandbox): blind, it hunts by sound. 12 m from a player who keeps still and crouched (silent), it
//     stays put; one dig-loud noise at the player and it comes, in bursts, and strikes within a few seconds.
//   creatures/lurker-ceiling (seed 42, in the nearest chamber): idle under a roof, it climbs onto the ceiling; a player walking
//     under it is heard, it rattles and drops on them (x1.5, a heavy hit).
#include "scenario.h"
#include "animals.h"
#include "creatures.h"
#include "player.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct St { uint64_t net = 0; float d0 = 0.0f; float health0 = 100.0f; Vector3 floor{}, up{}; int drops0 = 0; double hitAt = -1.0; } S;

Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
float DistToIt(ScenarioCtx& c) { const int i = Animals_FindNet(S.net); return i < 0 ? -1.0f : Vector3Distance(Animals_Position(i), FeetLocal(c)); }

void ThreatsOn() { Animals_SetNaturalSpawns(false); Animals_SetThreats(true); worldDifficulty = 1; }
}

SCENARIO("creatures/arena-lurker") {
    s.tags = { "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        ThreatsOn();
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("a lurker 12 m off", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        const Vector3 feet = FeetLocal(c), up = Vector3Normalize(feet);
        const Vector3 f = c.body->dirToLocal(player.forward);
        const Vector3 flat = Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
        const Vector3 dir = Vector3Normalize(Vector3Add(feet, Vector3Scale(flat, 12.0f)));
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "tunnel_lurker", Vector3Scale(dir, c.field->radius + c.field->elevation(dir)), 1, false, 0.01f);
        if (!c.check("spawned", Animals_Count() > before)) { c.fail("no lurker"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        player.posture = Player::Posture::Crouch;
        survival.health = 100.0f;
        return true;
    });
    s.step("silent: it stays", 6.0, [](ScenarioCtx& c) {
        player.posture = Player::Posture::Crouch;
        if (c.t < 0.2) { S.d0 = DistToIt(c); return false; }
        if (c.t < 4.0) return false;
        const float d = DistToIt(c);
        c.gt("still about 12 m off after 4 s of silence (m)", d, S.d0 - 3.0f);
        c.near("no hit", 100.0f - survival.health, 0.0, 0.5);
        return true;
    });
    s.step("a loud noise: it comes and strikes", 10.0, [](ScenarioCtx& c) {
        player.posture = Player::Posture::Crouch;
        static bool made = false;
        if (!made) { made = true; S.health0 = survival.health; Noise_Emit(c.body->id, FeetLocal(c), NoiseTune::kDig, Animals_Clock(), kPlayerSource); }
        if (survival.health < S.health0 - 0.5f && S.hitAt < 0.0) S.hitAt = c.t;
        // (The noise rides the animals' clock: emit it again each step for the first 0.5 s, so the clock base doesn't matter.)
        if (c.t < 0.5) Noise_Emit(c.body->id, FeetLocal(c), NoiseTune::kDig, Animals_Clock(), kPlayerSource);
        if (S.hitAt < 0.0) { c.waitingFor("the lurker's strike"); return false; }
        c.lt("struck within 8 s of the noise", S.hitAt, 8.0);
        c.near("a slash for 9", S.health0 - survival.health, 9.0, 0.6);
        survival.health = 100.0f;
        return true;
    });
}

SCENARIO("creatures/lurker-ceiling") {
    s.tags = { "creatures", "caves" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;
    s.record = false;

    s.step("settle", 120.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("a lurker on the floor of the nearest chamber", 60.0, [](ScenarioCtx& c) {
        static bool placed = false;
        if (!placed) {
            ThreatsOn();
            Vector3 at, out;
            if (!c.check("a chamber within 3 km", c.field->nearestCaveMouth(c.field->findSpawnDir(c.body->seed), 3000.0f, at, out, 3))) { c.fail("no chamber"); return true; }
            S.up = Vector3Normalize(at);
            bool found = false;
            for (float d = 0.0f; d < 30.0f && !found; d += 0.25f) {
                const Vector3 q = Vector3Subtract(at, Vector3Scale(S.up, d));
                if (c.field->density(q, nullptr) >= ISO_LEVEL) { S.floor = Vector3Add(q, Vector3Scale(S.up, 0.3f)); found = true; }
            }
            if (!c.check("the chamber's floor", found)) { c.fail("no floor"); return true; }
            // The player stands off to one side on the floor (so the chamber builds); the lurker 6 m away.
            player.flyMode = false;
            player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(S.floor, Vector3Scale(S.up, 0.2f))));
            placed = true;
            return false;
        }
        c.waitingFor("the chamber built (an idle world)");
        if (c.t < 2.0 || !WorldIdle()) return false;
        const Vector3 e = Vector3Normalize(Vector3CrossProduct(S.up, fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "tunnel_lurker", Vector3Add(S.floor, Vector3Scale(e, 6.0f)), 1, false, 0.01f);
        // Animals_Spawn drops from above (the surface): put it on the chamber floor instead.
        if (!c.check("spawned", Animals_Count() > before)) { c.fail("no lurker"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        Animals_TestPin(Animals_Count() - 1, Vector3Add(S.floor, Vector3Scale(e, 6.0f)), e);
        player.posture = Player::Posture::Crouch; // silent while it settles
        return true;
    });
    s.step("it climbs onto the ceiling", 15.0, [](ScenarioCtx& c) {
        player.posture = Player::Posture::Crouch;
        const int i = Animals_FindNet(S.net);
        if (i < 0) { c.fail("the lurker is gone"); return true; }
        if (!Animals_Clinging(i)) { c.waitingFor("the lurker on the ceiling"); return false; }
        c.check("clinging", true);
        const Vector3 p = Animals_Position(i);
        c.gt("up off the floor (m)", Vector3DotProduct(Vector3Subtract(p, S.floor), S.up), 1.5);
        return true;
    });
    s.step("walk under it: the drop", 10.0, [](ScenarioCtx& c) {
        static bool moved = false;
        const int i = Animals_FindNet(S.net);
        if (i < 0) { c.fail("the lurker is gone"); return true; }
        if (!moved) {
            moved = true;
            player.posture = Player::Posture::Stand;
            // Right under it, on the floor below.
            const Vector3 p = Animals_Position(i);
            Vector3 under = p;
            for (float d = 0.0f; d < 12.0f; d += 0.25f) {
                const Vector3 q = Vector3Subtract(p, Vector3Scale(S.up, d + 0.5f));
                if (c.field->density(q, nullptr) >= ISO_LEVEL) { under = Vector3Add(q, Vector3Scale(S.up, 0.3f)); break; }
            }
            player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(under, Vector3Scale(S.up, 0.2f))));
            S.health0 = survival.health = 100.0f;
            S.drops0 = Animals_LurkerDrops();
            return false;
        }
        Noise_Emit(c.body->id, FeetLocal(c), NoiseTune::kWalk, Animals_Clock(), kPlayerSource); // footsteps (a walking player under it)
        if (Animals_LurkerDrops() == S.drops0) { c.waitingFor("the drop"); return false; }
        c.check("it dropped on the player", true);
        c.near("x1.5: 13.5", S.health0 - survival.health, 13.5, 0.8);
        c.check("off the ceiling", !Animals_Clinging(i));
        survival.health = 100.0f;
        return true;
    });
}
