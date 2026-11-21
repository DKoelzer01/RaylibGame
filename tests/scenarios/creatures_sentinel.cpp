// creatures/arena-sentinel (#1228): a ruin sentinel (sandbox, by day) 15 m from a player who holds their ground.
//   It charges (the telegraph) and shoots darts that hit; it doesn't come in. Blunt hits do x1.5 (10 -> 15, less armour 2),
//   a slash is mostly turned (40% resisted, armour 2); a heavy slash doesn't make it bleed, a bite doesn't poison it.
//   Destroyed, it drops a basic circuit.
#include "scenario.h"
#include "animals.h"
#include "combat.h"
#include "creatures.h"
#include "damage.h"
#include "inventory.h"
#include "player.h"
#include "playersession.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <string>

namespace {
struct St { uint64_t net = 0; Vector3 stand{}; float d0 = 0.0f; int thrown0 = 0, hits0 = 0, tele0 = 0, circuit0 = 0; } S;
Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
int Idx() { return Animals_FindNet(S.net); }
float Hit(DamageType t, bool heavy = false, DamageKind k = DamageKind::Melee) {
    DamageEvent e; e.amount = 10.0f; e.type = t; e.heavy = heavy; e.kind = k;
    float dealt = 0.0f; Animals_HitEvent(Idx(), e, &dealt); return dealt;
}
}

SCENARIO("creatures/arena-sentinel") {
    s.tags = { "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.args = { "--day" };

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false); Animals_SetThreats(true); worldDifficulty = 1;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("a sentinel 15 m off", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        const Vector3 feet = FeetLocal(c), up = Vector3Normalize(feet);
        const Vector3 f = c.body->dirToLocal(player.forward);
        const Vector3 flat = Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
        const Vector3 dir = Vector3Normalize(Vector3Add(feet, Vector3Scale(flat, 15.0f)));
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "ruin_sentinel", Vector3Scale(dir, c.field->radius + c.field->elevation(dir)), 1, false, 0.01f);
        if (!c.check("spawned", Animals_Count() > before)) { c.fail("no sentinel"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        S.stand = feet;
        S.thrown0 = Animals_ProjectilesThrown(); S.hits0 = Animals_ProjectileHits(); S.tele0 = Creature_EventCount(CreatureEvent::Telegraph);
        survival.health = 100.0f;
        return true;
    });
    s.step("it charges and shoots, and stays put", 12.0, [](ScenarioCtx& c) {
        survival.health = std::max(survival.health, 40.0f); survival.dead = false;
        if (Vector3Distance(FeetLocal(c), S.stand) > 0.5f) player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(S.stand, Vector3Scale(Vector3Normalize(S.stand), 0.05f))));
        if (c.t < 0.2) { S.d0 = Vector3Distance(Animals_Position(Idx()), FeetLocal(c)); return false; }
        if (c.t < 9.0) return false;
        c.gt("charge-ups (telegraphs)", Creature_EventCount(CreatureEvent::Telegraph) - S.tele0, 1.5);
        c.gt("darts shot in 9 s", Animals_ProjectilesThrown() - S.thrown0, 1.5);
        c.gt("darts that hit", Animals_ProjectileHits() - S.hits0, 0.5);
        c.gt("it stayed back (m)", Vector3Distance(Animals_Position(Idx()), FeetLocal(c)), S.d0 - 3.0f);
        survival.health = 100.0f;
        return true;
    });
    s.step("weak to blunt, hard to cut, no bleed or poison", 3.0, [](ScenarioCtx& c) {
        if (Idx() < 0) { c.fail("the sentinel is gone"); return true; }
        c.near("blunt: x1.5 less armour (10 -> 13)", Hit(DamageType::Blunt), 13.0, 0.01);
        c.near("slash: 40% resisted less armour (10 -> 4)", Hit(DamageType::Slash), 4.0, 0.01);
        Hit(DamageType::Slash, true);
        const StatusState* st = Animals_Status(Idx());
        c.check("a heavy slash: no bleed", st && st->bleed <= 0.0f);
        Hit(DamageType::Pierce, false, DamageKind::Bite);
        st = Animals_Status(Idx());
        c.check("a bite: no poison", st && st->poison <= 0.0f);
        S.circuit0 = LocalSession().inventory.count("basic_circuit");
        { DamageEvent e; e.amount = 999.0f; e.type = DamageType::Blunt; Animals_HitEvent(Idx(), e); }
        c.eq("destroyed, it drops a basic circuit", LocalSession().inventory.count("basic_circuit") - S.circuit0, 1);
        return true;
    });
}
