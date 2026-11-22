// combat/weapons-composed (#1219, combat 5): a composed weapon from the real recipes, in hand, against the training dummy.
//   craft:  two bronze blades and planks hafted at the grid (recipes.txt) make one weapon_sword whose record names a bronze
//           head (the assembly hook, #1501); no plain sword.
//   hold:   it in the selected hotbar cell; the held copy's record is what ToolOf and the damage read (gToolDataOf).
//   hit:    a light tap at a hornhide held side-on 1.3 m ahead (as combat/melee): one hit, damage (3 + 3 x bronze's hardness 2)
//           = 9 on the hit's frame (bare hands do 2; the slash then bleeds on, as every slash does), the sword's light timing
//           (wind-up 0.30 s), and the copy's own wear rises by the swing.
#include "scenario.h"
#include "animals.h"
#include "composed.h"
#include "gameui.h"
#include "hotbar.h"
#include "inventory.h"
#include "itemdata.h"
#include "melee.h"
#include "damage.h"
#include "combat.h"
#include "player.h"
#include "playerinput.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
struct St {
    uint64_t net = 0; Vector3 at{}, face{};
    float hp0 = 0.0f, hpAtHit = -1.0f; int hits0 = 0, swings0 = 0, wear0 = 0, cell = -1;
    double tStart = -1, tHit = -1, tIdle = -1;
    long pressFrame = -1;
} S;

constexpr uint64_t kLmb = 1ull << (int)Act::Primary;

void Stand(ScenarioCtx& c) { player.flyMode = false; player.spawnOn(c.world, *c.body, c.spawnDirWorld()); }
Vector3 EyeLocal(ScenarioCtx& c) { return c.body->toLocal(c.world.camera.position); }
Vector3 GroundAhead(ScenarioCtx& c, float d, Vector3& right) {
    const Vector3 eye = EyeLocal(c), up = Vector3Normalize(eye);
    const Vector3 f = c.body->dirToLocal(player.forward);
    const Vector3 flat = Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
    right = Vector3Normalize(Vector3CrossProduct(flat, up));
    const Vector3 dir = Vector3Normalize(Vector3Add(eye, Vector3Scale(flat, d)));
    return Vector3Scale(dir, c.field->radius + c.field->elevation(dir));
}
int Idx() { return S.net ? Animals_FindNet(S.net) : -1; }
void Pin() { if (Idx() >= 0) Animals_TestPin(Idx(), S.at, S.face); }
void LookAtDummy(ScenarioCtx& c) {
    const Vector3 centre = Vector3Add(S.at, Vector3Scale(Vector3Normalize(S.at), 0.6f));
    const Vector3 to = Vector3Subtract(c.body->toWorld(centre), c.world.camera.position);
    if (Vector3Length(to) > 0.1f) player.forward = Vector3Normalize(to);
}
int SwordCell() { for (int i = 0; i < Inventory::kCells; ++i) if (inventory.grid()[i].id == "weapon_sword" && inventory.grid()[i].hasData()) return i; return -1; }
}

SCENARIO("combat/weapons-composed") {
    s.tags = { "combat" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.args = { "--day" };

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Stand(c);
        return true;
    });
    s.step("haft two bronze blades into a sword", 5.0, [](ScenarioCtx& c) {
        const Recipe* sword = nullptr;
        for (const Recipe& r : itemDB.recipes()) if (r.output == "weapon_sword" && !r.inputs.empty() && r.inputs[0].first == "bronze_blade") sword = &r;
        if (!c.check("a bronze sword recipe", sword != nullptr)) { c.fail("no recipe"); return true; }
        for (const auto& in : sword->inputs) inventory.add(ItemMatches(in.first, "oak_planks") ? "oak_planks" : in.first, in.second);
        c.check("crafted", inventory.craft(*sword));
        c.eq("composed swords", inventory.countInstances("weapon_sword"), 1);
        c.eq("plain swords", (int)inventory.count("weapon_sword"), 0);
        S.cell = SwordCell();
        if (!c.check("the sword in a cell", S.cell >= 0)) { c.fail("no sword"); return true; }
        const std::string d = inventory.grid()[S.cell].data;
        c.check("its head is bronze", ItemData::Get(d, "h") == "bronze", d);
        c.near("its damage from the head", Weapon_BaseDamage("weapon_sword", d), 9.0, 1e-4);
        c.check("it fights as a sword", Weapon_ClassOf("weapon_sword") == WeaponClass::Sword);
        return true;
    });
    s.step("a dummy ahead, the sword in hand", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        Vector3 right;
        S.at = GroundAhead(c, 1.3f, right);
        S.face = right;
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "hornhide", S.at, 1, false, 0.01f);
        if (!c.check("a hornhide for the dummy", Animals_Count() > before)) { c.fail("no dummy"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        Hotbar_Bind(3, "weapon_sword"); // key 4, as a player binds it
        gScenarioInput = true;
        gFakeSlot = 3;                  // and presses 4
        Pin(); LookAtDummy(c);
        return true;
    });
    s.step("the sword selected (key 4)", 3.0, [](ScenarioCtx& c) {
        Pin(); LookAtDummy(c);
        if (c.t < 0.3) return false;
        c.eq("hotbar slot selected", Hotbar_Selected(), 3);
        c.check("the sword in hand", terrainTool.heldTool == "weapon_sword", "held '" + terrainTool.heldTool + "'");
        S.cell = SwordCell();
        c.check("the selected cell holds the copy", S.cell == Inventory::kGrid + 3, "cell " + std::to_string(S.cell));
        return true;
    });
    s.step("a light swing with the sword", 6.0, [](ScenarioCtx& c) {
        Pin(); LookAtDummy(c);
        if (c.t < 0.3) return false;
        const MeleeState& m = Animals_Melee();
        if (S.pressFrame < 0) {
            S.hp0 = Animals_HealthOf(Idx()); S.hits0 = Animals_MeleeHits(); S.swings0 = m.swings;
            S.wear0 = ItemData::GetInt(inventory.grid()[S.cell].data, "w", 0);
            gFakeHeld |= kLmb; S.pressFrame = c.frame; c.cue();
            return false;
        }
        if (c.frame > S.pressFrame) gFakeHeld &= ~kLmb;
        if (S.tStart < 0 && m.swings > S.swings0) S.tStart = c.T;
        if (S.tStart >= 0 && S.tHit < 0 && m.phase == MeleePhase::Hit) S.tHit = c.T;
        if (S.tHit >= 0 && S.tIdle < 0 && m.phase == MeleePhase::Idle) S.tIdle = c.T;
        if (S.hpAtHit < 0 && Animals_MeleeHits() > S.hits0) S.hpAtHit = Animals_HealthOf(Idx()); // the hit's own damage (a slash bleeds on after it)
        if (S.tIdle < 0 || c.t < 1.5) { c.waitingFor("the sword's swing to finish"); return false; }
        c.check("the swing started", S.tStart >= 0);
        c.near("bronze sword light wind-up (s): the class's x its weight's time scale", S.tHit - S.tStart, Weapon_SpecOf("weapon_sword", inventory.grid()[S.cell].data).windUp, 1.5 / 60.0);
        c.eq("one hit", Animals_MeleeHits() - S.hits0, 1);
        c.near("light hit with a bronze sword (on the hit's frame)", S.hp0 - S.hpAtHit, 9.0, 1e-3);
        c.check("then the slash bleeds", S.hp0 - Animals_HealthOf(Idx()) > 9.0 + 0.1f,
                "total " + std::to_string(S.hp0 - Animals_HealthOf(Idx())) + " after " + std::to_string(c.t) + " s");
        const int cell = SwordCell();
        c.check("the sword is still held", cell >= 0);
        if (cell >= 0) c.check("its own wear rose", ItemData::GetInt(inventory.grid()[cell].data, "w", 0) > S.wear0,
                               "w " + std::to_string(S.wear0) + " -> " + inventory.grid()[cell].data);
        c.keepCapture("a composed bronze sword's light hit (#1219)");
        return true;
    });
    // (#1219) A ruby in its one socket: the next hit burns (instead of bleeding).
    s.step("a ruby set in the sword: the next hit burns", 8.0, [](ScenarioCtx& c) {
        Pin(); LookAtDummy(c);
        static long press = -1; static int hits0 = 0;
        if (press < 0) {
            if (c.t < 0.2) return false;
            const int cell = SwordCell();
            if (!c.check("the sword", cell >= 0)) { c.fail("no sword"); return true; }
            c.check("the ruby set", inventory.setInstanceData(cell, ItemData::Set(inventory.grid()[cell].data, "j", "ruby")));
            const StatusState* st0 = Animals_Status(Idx());
            if (st0) c.note("burn before: " + std::to_string(st0->burn));
            hits0 = Animals_MeleeHits();
            if (Animals_Melee().phase != MeleePhase::Idle) return false;
            gFakeHeld |= kLmb; press = c.frame;
            return false;
        }
        if (c.frame > press) gFakeHeld &= ~kLmb;
        if (Animals_MeleeHits() == hits0) { c.waitingFor("the second hit"); return false; }
        const StatusState* st = Animals_Status(Idx());
        c.check("it burns", st && st->burn > 0.0f, st ? "burn " + std::to_string(st->burn) : std::string("gone"));
        gScenarioInput = false;
        return true;
    });
    // (#1509) A client's hit on the server: the event the server builds from the session's held copy, through the typed net
    // path. A sapphire instead of the ruby: the dummy chills, and the damage is the sword's (as a local hit's).
    s.step("a client's hit on the server keeps its effect", 3.0, [](ScenarioCtx& c) {
        Pin();
        const int cell = SwordCell();
        if (!c.check("the sword", cell >= 0)) { c.fail("no sword"); return true; }
        const std::string d = ItemData::Set(inventory.grid()[cell].data, "j", "sapphire");
        const DamageEvent ev = Weapon_HitEvent("weapon_sword", d);
        c.check("the event chills", ev.status == StatusKind::Chill);
        c.near("the event's damage", ev.amount, 9.0, 1e-3);
        std::vector<std::pair<std::string, int>> drops;
        const float hp = Animals_HealthOf(Idx());
        float dealt = 0.0f;
        Animals_HitNetEvent(S.net, ev, drops, false, &dealt);
        const StatusState* st = Animals_Status(Idx());
        c.check("the dummy chills", st && st->chill > 0.0f, st ? "chill " + std::to_string(st->chill) : std::string("gone"));
        c.near("the health it took", hp - Animals_HealthOf(Idx()), dealt, 1e-3);
        c.check("typed damage", dealt > 0.0f, "dealt " + std::to_string(dealt));
        return true;
    });
    // (#1219, user 2026-11-22) Setting a gem in the inventory: a steel sword (two sockets) and three emeralds on the cursor,
    // clicked onto it through the slot rules: one is set each click; with both sockets taken, the third click swaps instead.
    s.step("emeralds dropped onto a steel sword: set, then swap when full", 3.0, [](ScenarioCtx& c) {
        const int before = inventory.countInstances("weapon_sword");
        const int cell = inventory.addInstance("weapon_sword", Composed_MakeData("steel"));
        if (!c.check("a steel sword", cell >= 0 && inventory.countInstances("weapon_sword") == before + 1)) { c.fail("no room"); return true; }
        const int durability = Composed_Stats("weapon_sword", inventory.grid()[cell].data).durability;
        InvSlot held; held.id = "emerald"; held.centi = 300;
        held = GameUI_TestClickHolding(held, cell);
        c.eq("first click: one emerald used", (int)(held.centi / 100), 2);
        held = GameUI_TestClickHolding(held, cell);
        c.eq("second click: another", (int)(held.centi / 100), 1);
        const ComposedStats s = Composed_Stats("weapon_sword", inventory.grid()[cell].data);
        c.eq("two gems set", (int)s.gems.size(), 2);
        c.eq("each emerald +25% durability", s.durability, (int)std::lround(std::lround(durability * 1.25f) * 1.25f));
        const InvSlot after = GameUI_TestClickHolding(held, cell);
        c.check("full: the third click swaps (the sword onto the cursor)", after.id == "weapon_sword", "cursor " + after.id);
        return true;
    });
}
