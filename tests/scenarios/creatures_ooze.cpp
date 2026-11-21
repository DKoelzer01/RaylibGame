// creatures/arena-ooze (#1226): the mineral ooze (sandbox, by day).
//   Blind, it comes to a walking player 4 m off (their steps), and its touch corrodes the axe in their hand (15 strokes).
//   Cut (a slash) it splits in two: the new one has half its health, a generation on, the same ore. A blunt hit, or a slash
//   while it burns, doesn't split it. Killed, it drops the ore it ate.
#include "scenario.h"
#include "animals.h"
#include "combat.h"
#include "creatures.h"
#include "damage.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "playersession.h"
#include "survival.h"
#include "toolwear.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct St { uint64_t net = 0; int wear0 = 0; float health0 = 0.0f; int count0 = 0, splits0 = 0; std::string ore; int ore0 = 0; } S;
Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
int Idx() { return Animals_FindNet(S.net); }
}

SCENARIO("creatures/arena-ooze") {
    s.tags = { "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.args = { "--day" };

    s.step("settle, an axe in hand", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false); Animals_SetThreats(true); worldDifficulty = 1;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        InvSlot& cell = LocalSession().inventory.grid()[Inventory::kGrid + 0];
        cell.id = "stone_axe"; cell.centi = 100;
        LocalSession().hotbar.selected = 0;
        return true;
    });
    s.step("an ooze 4 m off", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        const Vector3 feet = FeetLocal(c), up = Vector3Normalize(feet);
        const Vector3 f = c.body->dirToLocal(player.forward);
        const Vector3 flat = Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
        const Vector3 dir = Vector3Normalize(Vector3Add(feet, Vector3Scale(flat, 4.0f)));
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "mineral_ooze", Vector3Scale(dir, c.field->radius + c.field->elevation(dir)), 1, false, 0.01f);
        if (!c.check("spawned", Animals_Count() > before)) { c.fail("no ooze"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        S.wear0 = ToolWear_Used("stone_axe");
        S.health0 = survival.health = 100.0f;
        return true;
    });
    s.step("it comes to the steps and its touch corrodes the axe", 15.0, [](ScenarioCtx& c) {
        Noise_Emit(c.body->id, FeetLocal(c), NoiseTune::kSprint * 0.4f, Animals_Clock(), kPlayerSource); // footsteps (a walk)
        survival.health = std::max(survival.health, 50.0f);
        if (ToolWear_Used("stone_axe") - S.wear0 < 15) { c.waitingFor("its touch on the axe"); return false; }
        c.eq("15 strokes of wear from one touch", ToolWear_Used("stone_axe") - S.wear0, 15);
        return true;
    });
    s.step("a cut splits it; blunt and burning don't", 3.0, [](ScenarioCtx& c) {
        const int i = Idx();
        if (i < 0) { c.fail("the ooze is gone"); return true; }
        S.count0 = Animals_Count(); S.splits0 = Animals_OozeSplits();
        S.ore = Animals_Ore(i);
        c.check("it has an ore", !S.ore.empty(), S.ore);
        { DamageEvent e; e.amount = 2.0f; e.type = DamageType::Blunt; Animals_HitEvent(i, e); }
        c.eq("a blunt hit: no split", Animals_OozeSplits() - S.splits0, 0);
        const float before = Animals_HealthOf(Idx());
        { DamageEvent e; e.amount = 2.0f; e.type = DamageType::Slash; Animals_HitEvent(Idx(), e); }
        c.eq("a slash: split", Animals_OozeSplits() - S.splits0, 1);
        c.eq("two oozes", Animals_Count() - S.count0, 1);
        const int b = Animals_Count() - 1;
        c.near("each has half", Animals_HealthOf(b), (before - 2.0f) * 0.5f, 0.6);
        c.eq("the new one a generation on", Animals_Generation(b), 1);
        c.check("the same ore", Animals_Ore(b) == S.ore);
        { DamageEvent e; e.amount = 1.0f; e.type = DamageType::Heat; Animals_HitEvent(Idx(), e); } // set it burning
        { DamageEvent e; e.amount = 2.0f; e.type = DamageType::Slash; Animals_HitEvent(Idx(), e); }
        c.eq("a slash while it burns: no split", Animals_OozeSplits() - S.splits0, 1);
        // Killed: the ore. (The split-off one: 1.)
        S.ore0 = LocalSession().inventory.count(S.ore);
        { DamageEvent e; e.amount = 999.0f; Animals_HitEvent(b, e); }
        c.eq("killed, it drops its ore", LocalSession().inventory.count(S.ore) - S.ore0, 1);
        return true;
    });
}
