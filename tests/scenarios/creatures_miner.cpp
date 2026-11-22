// creatures/arena-miner (#1227): the lost miner (sandbox; seed 42 starts at night).
//   In the dark it sleeps through a player 5 m off (no light in hand). A torch in hand wakes it, and it comes and strikes.
//   Its plank shield: a slash at its face does almost nothing (85% taken, then its armour); a blunt hit at its face, or a
//   slash at its back, lands. Killed, it drops its pick and scrap.
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
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <cstring>
#include <string>

namespace {
struct St { uint64_t net = 0; Vector3 at{}; float d0 = 0.0f, health0 = 100.0f; int blocks0 = 0, scrap0 = 0; } S;
Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
int Idx() { return Animals_FindNet(S.net); }
float Dist(ScenarioCtx& c) { const int i = Idx(); return i < 0 ? -1.0f : Vector3Distance(Animals_Position(i), FeetLocal(c)); }
void Hand(const char* item) { // the selected hotbar cell (an empty one first) holds `item` ("" none)
    InvSlot& cell = LocalSession().inventory.grid()[Inventory::kGrid + 9];
    cell.id = item; cell.centi = *item ? 100 : 0;
    LocalSession().hotbar.selected = 9;
    terrainTool.heldTool = item;
}
float Hit(DamageType t, bool front) {
    DamageEvent e; e.amount = 10.0f; e.type = t; e.fromFront = front; e.fromBehind = !front;
    float dealt = 0.0f; Animals_HitEvent(Idx(), e, &dealt); return dealt;
}
}

SCENARIO("creatures/arena-miner") {
    s.tags = { "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;

    s.step("settle, empty-handed", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false); Animals_SetThreats(true); worldDifficulty = 1;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        Hand("");
        return true;
    });
    s.step("a lost miner 5 m off", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        Hand("");
        const Vector3 feet = FeetLocal(c), up = Vector3Normalize(feet);
        const Vector3 f = c.body->dirToLocal(player.forward);
        const Vector3 flat = Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
        const Vector3 dir = Vector3Normalize(Vector3Add(feet, Vector3Scale(flat, 5.0f)));
        S.at = Vector3Scale(dir, c.field->radius + c.field->elevation(dir));
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "lost_miner", S.at, 1, false, 0.01f);
        if (!c.check("spawned", Animals_Count() > before)) { c.fail("no miner"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        survival.health = 100.0f;
        return true;
    });
    s.step("in the dark it sleeps", 6.0, [](ScenarioCtx& c) {
        Hand("");
        if (c.t < 0.2) { S.d0 = Dist(c); return false; }
        if (c.t < 3.0) return false;
        c.check("asleep", std::strcmp(Animals_MoodName(Idx()), "sleep") == 0, Animals_MoodName(Idx()));
        c.gt("still about 5 m off (m)", Dist(c), S.d0 - 1.0f);
        c.near("no hit", 100.0f - survival.health, 0.0, 0.5);
        return true;
    });
    s.step("a torch in hand wakes it; it comes and strikes", 12.0, [](ScenarioCtx& c) {
        Hand("torch");
        survival.health = std::max(survival.health, 50.0f);
        static bool hit = false;
        if (survival.health < S.health0 - 0.5f) hit = true;
        S.health0 = survival.health;
        if (!hit) { c.waitingFor(std::string("its strike (mood ") + Animals_MoodName(Idx()) + ")"); return false; }
        c.check("awake", std::strcmp(Animals_MoodName(Idx()), "sleep") != 0, Animals_MoodName(Idx()));
        c.check("it struck", true);
        survival.health = 100.0f;
        return true;
    });
    s.step("the plank shield", 3.0, [](ScenarioCtx& c) {
        const int i = Idx();
        if (i < 0) { c.fail("the miner is gone"); return true; }
        S.blocks0 = Animals_ShieldBlocks();
        c.lt("a slash at its face (85% taken, then armour 3)", Hit(DamageType::Slash, true), 1.0);
        c.eq("the shield took it", Animals_ShieldBlocks() - S.blocks0, 1);
        c.near("a blunt hit at its face lands (10 - armour 3)", Hit(DamageType::Blunt, true), 7.0, 0.01);
        c.near("a slash at its back lands", Hit(DamageType::Slash, false), 7.0, 0.01);
        S.scrap0 = LocalSession().inventory.count("scrap");
        { DamageEvent e; e.amount = 999.0f; e.type = DamageType::Blunt; Animals_HitEvent(Idx(), e); }
        c.eq("killed, it drops scrap", LocalSession().inventory.count("scrap") - S.scrap0, 2);
        Hand("");
        return true;
    });
}
