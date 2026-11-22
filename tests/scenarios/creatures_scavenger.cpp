// creatures/arena-scavenger (#1225): a group of three cave scavengers (sandbox) 12 m from a player.
//   The thrower keeps 5 m+ away and throws stones (thrown counted; hits or blocks land on the player); the other two close
//   in and slash. One hurt under 40% retreats (gets further away). A pile dropped by a scavenger is taken whole by it (a
//   theft), and when the thief falls the player gets it back with its loot.
#include "scenario.h"
#include "animals.h"
#include "combat.h"
#include "inventory.h"
#include "player.h"
#include "playersession.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <string>

namespace {
struct St { Vector3 stand{}; uint64_t net[3] = {}; float throwerMin = 1e9f, bladeMin = 1e9f; int thrown0 = 0; float health0 = 100.0f; float d0 = 0.0f; int thefts0 = 0; int stone0 = 0; } S;

Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
int Idx(int k) { return Animals_FindNet(S.net[k]); }
float Dist(ScenarioCtx& c, int k) { const int i = Idx(k); return i < 0 ? -1.0f : Vector3Distance(Animals_Position(i), FeetLocal(c)); }
Vector3 Ground(ScenarioCtx& c, Vector3 near, Vector3 dirFlat, float d) {
    const Vector3 dir = Vector3Normalize(Vector3Add(near, Vector3Scale(dirFlat, d)));
    return Vector3Scale(dir, c.field->radius + c.field->elevation(dir));
}
Vector3 Flat(ScenarioCtx& c) {
    const Vector3 feet = FeetLocal(c), up = Vector3Normalize(feet), f = c.body->dirToLocal(player.forward);
    return Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
}
}

SCENARIO("creatures/arena-scavenger") {
    s.tags = { "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    s.args = { "--day" }; // daylight: their sight is their own, not the dark's (a separate matter: torches light only 9 m)

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false); Animals_SetThreats(true); worldDifficulty = 1;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("a group of three, 12 m off", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        const Vector3 feet = FeetLocal(c), flat = Flat(c), side = Vector3Normalize(Vector3CrossProduct(flat, Vector3Normalize(feet)));
        for (int k = 0; k < 3; ++k) {
            const int before = Animals_Count();
            Animals_Spawn(c.world, *c.body, "cave_scavenger", Ground(c, Ground(c, feet, flat, 12.0f), side, (k - 1) * 2.0f), 1, false, 0.01f);
            if (!c.check("spawned " + std::to_string(k), Animals_Count() > before)) { c.fail("no scavenger"); return true; }
            S.net[k] = Animals_NetId(Animals_Count() - 1);
        }
        for (int k = 0; k < 3; ++k) Animals_TestSetGroup(Idx(k), S.net[0], k == 0);
        S.thrown0 = Animals_ProjectilesThrown();
        S.stand = FeetLocal(c);
        survival.health = 100.0f; S.health0 = 100.0f;
        return true;
    });
    s.step("the thrower throws from a distance; the others close in", 15.0, [](ScenarioCtx& c) {
        survival.health = std::max(survival.health, 40.0f); survival.dead = false;
        // The player holds their ground (the rams' knockback would carry them out of the thrower's range).
        if (Vector3Distance(FeetLocal(c), S.stand) > 0.5f) player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(S.stand, Vector3Scale(Vector3Normalize(S.stand), 0.05f))));
        const float dt = Dist(c, 0);
        if (dt >= 0.0f && c.t > 1.0) S.throwerMin = std::min(S.throwerMin, dt);
        for (int k = 1; k < 3; ++k) { const float d = Dist(c, k); if (d >= 0.0f) S.bladeMin = std::min(S.bladeMin, d); }
        if (c.t < 10.0) return false;
        c.gt("stones thrown in 10 s", Animals_ProjectilesThrown() - S.thrown0, 1.5);
        c.gt("the thrower kept its distance (closest, m)", S.throwerMin, 4.0);
        c.lt("a blade-carrier came close (closest, m)", S.bladeMin, 2.5);
        c.metric("stones that hit", Animals_ProjectileHits(), "", false);
        survival.health = 100.0f;
        return true;
    });
    s.step("hurt under 40%, one retreats", 6.0, [](ScenarioCtx& c) {
        survival.health = std::max(survival.health, 40.0f);
        static bool hit = false;
        if (!hit) { hit = true; DamageEvent e; e.amount = 26.0f; Animals_HitEvent(Idx(1), e); S.d0 = Dist(c, 1); return false; }
        if (c.t < 3.0) return false;
        c.gt("it got further away (m)", Dist(c, 1) - S.d0, 3.0);
        return true;
    });
    s.step("a dropped pile: taken, and back when the thief falls", 15.0, [](ScenarioCtx& c) {
        survival.health = std::max(survival.health, 40.0f);
        static int phase = 0;
        if (phase == 0) {
            // Kill the other two (the test is about the thief), then drop 5 stone 4 m from the thief, out of the player's reach.
            for (int k = 0; k < 2; ++k) { const int i = Idx(k); if (i >= 0) { DamageEvent e; e.amount = 999.0f; Animals_HitEvent(i, e); } }
            const int t = Idx(2);
            if (t < 0) { c.fail("the thief is gone"); return true; }
            InvSlot st; st.id = "stone"; st.centi = 500;
            const Vector3 p = Animals_Position(t), up = Vector3Normalize(p);
            const Vector3 sideways = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            Combat_DropFor(c.body->name, DVec3(Ground(c, p, sideways, 4.0f)), st, "");
            S.thefts0 = Animals_Thefts();
            phase = 1;
            return false;
        }
        if (phase == 1) {
            if (Animals_Thefts() == S.thefts0) { c.waitingFor("the theft"); return false; }
            c.check("the scavenger took the pile", true);
            S.stone0 = LocalSession().inventory.count("stone");
            const int t = Idx(2);
            if (t < 0) { c.fail("the thief is gone"); return true; }
            DamageEvent e; e.amount = 999.0f; Animals_HitEvent(t, e);
            phase = 2;
            return false;
        }
        c.eq("the 5 stone back when it fell", LocalSession().inventory.count("stone") - S.stone0, 5);
        survival.health = 100.0f;
        return true;
    });
}
