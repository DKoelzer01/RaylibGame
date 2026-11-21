// combat/poise (#1217, combat step 3): poise wears down, breaks into a stagger, and the next hit does x1.5.
//   creature: a hornhide (40 health: 21 poise) takes hits of 10 poise; the third breaks it (staggered 1 s, its attack off);
//             the next hit of 4 does 6, the one after 4 again. While it staggers it stands still (a hit makes it flee
//             otherwise); once the stagger ends it runs.
//   small:    a clucker staggers on one heavy hit, however little poise it carries.
//   player:   no armour (15 poise): a strong bite (8 damage: 20 poise) staggers; in a heavy wind-up (hyper-armour) the same
//             bite lands without a stagger.
#include "scenario.h"
#include "animals.h"
#include "combat.h"
#include "damage.h"
#include "player.h"
#include "playersession.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct St { int idx = -1; uint64_t net = 0; Vector3 still{}; float poiseBefore = 0.0f; } S; // (net: across steps an index can shift)

int Fresh(ScenarioCtx& c, const char* species, float spread) {
    const Vector3 home = c.spawnDirLocal();
    const Vector3 at = Vector3Scale(home, c.field->radius + c.field->elevation(home));
    const int before = Animals_Count();
    Animals_Spawn(c.world, *c.body, species, at, 1, false, spread);
    return Animals_Count() > before ? Animals_Count() - 1 : -1;
}

float Hit(int i, float amount, float poise, bool heavy = false) {
    DamageEvent e; e.amount = amount; e.poise = poise; e.heavy = heavy;
    float dealt = 0.0f;
    Animals_HitEvent(i, e, &dealt);
    return dealt;
}
}

SCENARIO("combat/poise") {
    s.tags = { "combat" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;

    // The player stands at the start: animals are simulated only near a player (a bench camera in the air is none).
    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("stand", 5.0, [](ScenarioCtx& c) { return c.t >= 0.5; });
    s.step("a hornhide's poise breaks on the third hit", 2.0, [](ScenarioCtx& c) {
        S.idx = Fresh(c, "hornhide", 8.0f);
        if (S.idx < 0) { c.fail("no hornhide"); return true; }
        Hit(S.idx, 1.0f, 10.0f);
        Hit(S.idx, 1.0f, 10.0f);
        const PoiseState* p = Animals_Poise(S.idx);
        c.near("its poise (5 + 0.4 x 40)", p->max, 21.0, 1e-3);
        c.check("two hits: not staggered", !Poise_Staggered(*p));
        Hit(S.idx, 1.0f, 10.0f);
        c.check("the third: staggered", Poise_Staggered(*p));
        c.near("the stagger (s)", p->stagger, 1.0, 1e-3);
        c.near("the next hit does x1.5", Hit(S.idx, 4.0f, 1.0f), 6.0, 1e-3);
        c.near("only the next", Hit(S.idx, 4.0f, 1.0f), 4.0, 1e-3);
        S.still = Animals_Position(S.idx);
        S.net = Animals_NetId(S.idx);
        return true;
    });
    s.step("it stands while it staggers, then runs", 4.0, [](ScenarioCtx& c) {
        const PoiseState* p = Animals_Poise(Animals_FindNet(S.net));
        if (!p) { c.fail("the hornhide is gone"); return true; }
        if (c.t < 0.5) return false;
        static bool checked = false;
        if (!checked) {
            checked = true;
            c.check("still staggered at 0.5 s", Poise_Staggered(*p));
            c.lt("moved while staggered (m)", Vector3Distance(Animals_Position(Animals_FindNet(S.net)), S.still), 0.05);
        }
        if (c.t < 2.5) return false;
        c.check("the stagger is over", !Poise_Staggered(*p));
        c.gt("then it runs (m from where it stood)", Vector3Distance(Animals_Position(Animals_FindNet(S.net)), S.still), 0.5);
        return true;
    });
    s.step("a small creature staggers on any heavy", 2.0, [](ScenarioCtx& c) {
        const int i = Fresh(c, "clucker", 4.0f);
        if (i < 0) { c.fail("no clucker"); return true; }
        Hit(i, 0.5f, 0.5f);
        c.check("a light: no stagger", !Poise_Staggered(*Animals_Poise(i)));
        Hit(i, 0.5f, 0.5f, true);
        c.check("a heavy: staggered", Poise_Staggered(*Animals_Poise(i)));
        return true;
    });
    s.step("a strong bite staggers an unarmoured player", 2.0, [](ScenarioCtx& c) {
        survival.hyperArmour = false;
        DamageEvent e; e.amount = 8.0f; e.kind = DamageKind::Bite;
        const DamageResult r = Combat_Damage(LocalSession(), e);
        c.near("the player's poise (no armour)", survival.poise.max, 15.0, 1e-3);
        c.check("staggered", r.staggered && Poise_Staggered(survival.poise));
        survival.health = 100.0f;
        return true;
    });
    s.step("the stagger slows the walk, then wears off", 3.0, [](ScenarioCtx& c) {
        if (c.t < 0.1) return false;
        static bool slowChecked = false;
        if (!slowChecked) { slowChecked = true; c.lt("walk speed factor while staggered", player.statusSpeed, 0.5); }
        if (c.t < 0.8) return false;
        c.check("over after 0.6 s", !Poise_Staggered(survival.poise));
        c.gt("walk speed back", player.statusSpeed, 0.9);
        return true;
    });
    s.step("hyper-armour: the same bite in a heavy wind-up", 2.0, [](ScenarioCtx& c) {
        survival.poise.value = survival.poise.max;
        survival.hyperArmour = true; // as a heavy or two-handed wind-up sets it
        DamageEvent e; e.amount = 8.0f; e.kind = DamageKind::Bite;
        const DamageResult r = Combat_Damage(LocalSession(), e);
        survival.hyperArmour = false;
        c.check("the damage lands", r.dealt > 0.0f);
        c.check("no stagger", !r.staggered && !Poise_Staggered(survival.poise));
        survival.health = 100.0f;
        return true;
    });
}
