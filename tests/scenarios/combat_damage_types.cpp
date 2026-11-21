// combat/damage-types (#1215, combat step 1): a creature hit with each damage type. Each hit does its amount (a hornhide has
// no resistances yet, so 4 of 4) and brings its effect: heat burns, cold chills (and slows it), energy shocks (the wind-up
// interrupt), a slash heavy bleeds, a bite poisons; blunt and pierce bring none. Then a second of standing still: the burn,
// bleed and poison take health over time.
#include "scenario.h"
#include "animals.h"
#include "combat.h"
#include "damage.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct State { int idx = -1; uint64_t net = 0; float hpAfterHits = 0.0f; } S; // (net: its stable id; an index shifts when another animal goes)

int Fresh(ScenarioCtx& c) { // a new hornhide 6 m from the start (the last one spawned)
    const Vector3 home = c.spawnDirLocal();
    const Vector3 at = Vector3Scale(home, c.field->radius + c.field->elevation(home));
    const int before = Animals_Count();
    Animals_Spawn(c.world, *c.body, "hornhide", at, 1, false, 6.0f);
    return Animals_Count() > before ? Animals_Count() - 1 : -1;
}
}

SCENARIO("combat/damage-types") {
    s.tags = { "combat" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;

    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("1 s of play and an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("each type on a fresh hornhide", 1.0, [](ScenarioCtx& c) {
        struct Case { DamageType type; bool heavy; DamageKind kind; StatusKind want; const char* name; };
        const Case cases[] = {
            { DamageType::Blunt, false, DamageKind::Melee, StatusKind::None, "blunt" },
            { DamageType::Pierce, false, DamageKind::Melee, StatusKind::None, "pierce" },
            { DamageType::Slash, true, DamageKind::Melee, StatusKind::Bleed, "slash heavy" },
            { DamageType::Heat, false, DamageKind::Melee, StatusKind::Burn, "heat" },
            { DamageType::Cold, false, DamageKind::Melee, StatusKind::Chill, "cold" },
            { DamageType::Energy, false, DamageKind::Melee, StatusKind::Shock, "energy" },
            { DamageType::Pierce, false, DamageKind::Bite, StatusKind::Poison, "bite" },
        };
        for (const Case& k : cases) {
            const int i = Fresh(c);
            if (!c.check(std::string("a hornhide for ") + k.name, i >= 0)) continue;
            DamageEvent e; e.amount = 4.0f; e.type = k.type; e.heavy = k.heavy; e.kind = k.kind;
            float dealt = 0.0f;
            Animals_HitEvent(i, e, &dealt);
            c.near(std::string(k.name) + ": damage dealt", dealt, 4.0, 1e-4);
            const StatusState* st = Animals_Status(i);
            if (!st) { c.fail("no status"); return true; }
            bool got = false;
            switch (k.want) {
                case StatusKind::None: got = !st->any() && !st->interrupted; break;
                case StatusKind::Burn: got = st->burn > 0.0f; break;
                case StatusKind::Chill: got = st->chill > 0.0f && Status_SpeedFactor(*st) < 1.0f; break;
                case StatusKind::Bleed: got = st->bleed > 0.0f && st->bleedStacks == 1; break;
                case StatusKind::Poison: got = st->poison > 0.0f; break;
                case StatusKind::Shock: got = st->interrupted; break;
            }
            c.check(std::string(k.name) + ": brings " + StatusKind_Name(k.want), got);
        }
        // One hornhide that burns, bleeds and is poisoned at once: its health falls over time.
        S.idx = Fresh(c);
        if (S.idx < 0) { c.fail("no hornhide"); return true; }
        for (DamageType t : { DamageType::Heat }) { DamageEvent e; e.amount = 1.0f; e.type = t; Animals_HitEvent(S.idx, e); }
        { DamageEvent e; e.amount = 1.0f; e.type = DamageType::Slash; e.heavy = true; Animals_HitEvent(S.idx, e); }
        { DamageEvent e; e.amount = 1.0f; e.type = DamageType::Pierce; e.kind = DamageKind::Bite; Animals_HitEvent(S.idx, e); }
        S.net = Animals_NetId(S.idx);
        S.hpAfterHits = Animals_HealthOf(S.idx);
        return true;
    });
    s.step("a second of burn, bleed and poison", 5.0, [](ScenarioCtx& c) {
        if (c.t < 1.0) return false;
        const int i = Animals_FindNet(S.net); // (by its id: the list may have changed under the index)
        if (!c.check("the burning hornhide is still there", i >= 0, "index then " + std::to_string(S.idx) + ", animals " + std::to_string(Animals_Count()))) return true;
        c.note("simulated animals " + std::to_string(Animals_SimulatedCount()) + "; its index " + std::to_string(S.idx) + " -> " + std::to_string(i));
        const float lost = S.hpAfterHits - Animals_HealthOf(i);
        // burn 3/s + bleed 2/s + poison 1.5/s = 6.5 a second (the animal near the start simulates every step)
        c.between("health lost to effects over ~1 s", lost, 5.0, 8.5);
        return true;
    });
}
