// creatures/arena-thornback and creatures/arena-thornback-old (#1223, the acceptance gate): a thornback on the creature
// brain behaves as the hostile-animal code before it did. The same arena runs twice: once on the brain (creatures.txt), once
// with RAYLIBGAME_CREATURES=0 (the old code); both must pass the same checks.
//   It charges a player standing 8 m off (inside its 12 m), telegraphs, and rams for 8 (normal difficulty) within a few
//   seconds; struck at 0.5 s after the telegraph with a 1.2 s cooldown (two hits at least 1.6 s apart). The player then
//   steps 40 m away: past its 12 m (and not angry: never hit), it stops coming.
#include "scenario.h"
#include "animals.h"
#include "creatures.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct St { uint64_t net = 0; double firstHit = -1, secondHit = -1; float health0 = 100.0f; int tele0 = 0; Vector3 far{}; float dFar0 = 0.0f; } S;

Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }

Vector3 GroundAt(ScenarioCtx& c, Vector3 near, Vector3 flat, float d) {
    const Vector3 dir = Vector3Normalize(Vector3Add(near, Vector3Scale(flat, d)));
    return Vector3Scale(dir, c.field->radius + c.field->elevation(dir));
}

Vector3 Flat(ScenarioCtx& c) {
    const Vector3 feet = FeetLocal(c), up = Vector3Normalize(feet), f = c.body->dirToLocal(player.forward);
    return Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
}

void Build(ScenarioDef& s, bool old) {
    s.tags = { "creatures", "combat" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    if (old) s.env = { { "RAYLIBGAME_CREATURES", "0" } };
    s.step("settle, threats on", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false);
        Animals_SetThreats(true);
        worldDifficulty = 1; // normal: x1 damage
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("a thornback 8 m off", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        const Vector3 at = GroundAt(c, FeetLocal(c), Flat(c), 8.0f);
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "thornback", at, 1, false, 0.01f);
        if (!c.check("spawned", Animals_Count() > before)) { c.fail("no thornback"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        survival.health = 100.0f; survival.dead = false;
        S.health0 = survival.health;
        S.tele0 = Creature_EventCount(CreatureEvent::Telegraph);
        return true;
    });
    s.step("it charges and rams", 12.0, [](ScenarioCtx& c) {
        if (survival.health < S.health0 - 0.5f) {
            if (S.firstHit < 0) { S.firstHit = c.t; c.near("a ram: 8 at normal", S.health0 - survival.health, 8.0, 0.6); S.health0 = survival.health; }
            else if (S.secondHit < 0) { S.secondHit = c.t; S.health0 = survival.health; }
        }
        survival.health = std::max(survival.health, 30.0f); S.health0 = std::max(S.health0, 30.0f);
        if (S.secondHit < 0) { c.waitingFor("two rams"); return false; }
        c.lt("first ram (s after the spawn)", S.firstHit, 6.0);
        c.between("between the rams (s): wind-up 0.5 + cooldown 1.2", S.secondHit - S.firstHit, 1.6, 4.0);
        if (Creatures_On()) c.gt("telegraphs raised (creature events)", Creature_EventCount(CreatureEvent::Telegraph) - S.tele0, 1.5);
        survival.health = 100.0f;
        // Step 40 m back from it.
        const int i = Animals_FindNet(S.net);
        if (i < 0) { c.fail("the thornback is gone"); return true; }
        const Vector3 feet = FeetLocal(c), them = Animals_Position(i);
        Vector3 away = Vector3Subtract(feet, them);
        const Vector3 up = Vector3Normalize(feet);
        away = Vector3Normalize(Vector3Subtract(away, Vector3Scale(up, Vector3DotProduct(away, up))));
        S.far = GroundAt(c, feet, away, 40.0f);
        player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(S.far, Vector3Scale(Vector3Normalize(S.far), 0.2f))));
        return true;
    });
    s.step("40 m away, it stops coming", 8.0, [](ScenarioCtx& c) {
        const int i = Animals_FindNet(S.net);
        if (i < 0) { c.fail("the thornback is gone"); return true; }
        const float d = Vector3Distance(Animals_Position(i), FeetLocal(c));
        if (c.t < 1.0) { S.dFar0 = d; return false; }
        if (c.t < 5.0) return false;
        c.gt("still far after 4 s (m)", d, 25.0);
        c.gt("not closing in (m: now minus then)", d - S.dFar0, -3.0);
        c.near("no hit out there", 100.0f - survival.health, 0.0, 0.5);
        return true;
    });
}
}

SCENARIO("creatures/arena-thornback") { Build(s, false); }
SCENARIO("creatures/arena-thornback-old") { Build(s, true); }
