// input/rmb-* (#1716, playtest 2026-11-23: "right-click to interact on a machine or chest with an item with right-click
// capabilities in hand (telescope) opens the inventory and triggers the item"): an aimed interactable takes the RMB click, and the
// held item's own RMB (telescope, water bucket, food, shield block, ground placing, torch) fires only when nothing interactable
// is aimed. One scenario per held item; each runs the real input path (gFakeHeld RMB) against a chest, a furnace, a campfire and a
// bed (each in its own plot, aimed from 2.5 m on foot): the interactable's own effect happens (its screen opens, the bed becomes
// the spawn) and the item's does not; and a control, RMB at the bare ground with nothing aimed, where the item does fire.
// Not covered: a ladder (it is walked through and climbed; RMB on it does nothing, so no interactable takes the click there) and
// the fishing rod (it reads the raw mouse button, which a scenario cannot press; its handler takes the same Interact_RmbConsumed).
#include "scenario.h"
#include "buildtool.h"
#include "gameui.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "playersession.h"
#include "plots.h"
#include "survival.h"
#include "telescope.h"
#include "terraintool.h"
#include "world.h"
#include <functional>
#include <raymath.h>
#include <string>

namespace {
constexpr uint64_t kRmb = 1ull << (int)Act::Secondary;

struct Spec {
    const char* name;                       // the held item, for the check names
    const char* item;                       // put in the selected hotbar cell ("" bare hands)
    const char* offHand;                    // put in the off hand ("" none)
    int holdFrames;                         // RMB held this long at each case
    std::function<void(ScenarioCtx&)> reset; // before each case: put the item's state back
    std::function<bool()> fired;            // polled every frame of a case: did the item's RMB action happen?
};

struct Interactable { const char* name; PartType type; };
const Interactable kInteractables[] = { { "chest", PartType::Chest }, { "furnace", PartType::Furnace }, { "campfire", PartType::Campfire }, { "bed", PartType::Bed } };
constexpr int kN = (int)(sizeof(kInteractables) / sizeof(kInteractables[0]));

enum class Ph { Init, Spawn, Wait, Hold, Settle, Done };
struct Run {
    Ph ph = Ph::Init;
    int idx = 0;                // 0..kN-1 the interactables, kN the control (bare ground)
    long f = 0;                 // frames in the phase
    uint32_t plotId[kN + 1] = {};
    int part[kN + 1] = {};
    Vector3 at[kN + 1] = {};    // body-local: the part's centre, or the ground point of the control
    bool fired = false, screen = false, bed = false;
    int cell = 0;               // the hotbar cell holding the item
} R;

// Body-local position `k` x 25 m along the ground from the spawn: its direction and the ground radius there.
Vector3 Around(ScenarioCtx& c, int k, float& ground) {
    const Vector3 d0 = c.spawnDirLocal();
    const Vector3 ref = fabsf(d0.x) < 0.9f ? Vector3{ 1, 0, 0 } : Vector3{ 0, 0, 1 };
    const Vector3 t = Vector3Normalize(Vector3CrossProduct(d0, ref));
    const float r0 = c.field->radius + c.field->elevation(d0);
    const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(d0, r0), Vector3Scale(t, 25.0f * (float)k)));
    ground = c.field->radius + c.field->elevation(d);
    return d;
}
Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocalD(player.position).toFloat(); }

bool Drive(ScenarioCtx& c, const Spec& sp) {
    PlayerSession& ses = LocalSession();
    if (R.ph == Ph::Init) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        player.flyMode = false;
        terrainTool.creative = false;
        gScenarioInput = true;
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        for (int i = 0; i < kN; ++i) { // each interactable in a plot of its own, 25 m apart, a metre above the ground
            float ground = 0.0f;
            const Vector3 d = Around(c, i + 1, ground);
            Plot* plot = c.body->plots->plotFor(Vector3Scale(d, ground + 1.0f), 0.0f, true);
            R.plotId[i] = plot ? plot->id : 0;
            R.part[i] = plot ? plot->place({ kInteractables[i].type, { 0, 0, 0 }, 0 }) : -1;
            if (!c.check(std::string(kInteractables[i].name) + " placed", R.part[i] >= 0)) return true;
            R.at[i] = plot->toBody({ 0.5f, 0.5f, 0.5f });
        }
        { float ground = 0.0f; const Vector3 d = Around(c, 0, ground); R.at[kN] = Vector3Scale(d, ground); } // the control: bare ground
        // The held item in a free hotbar cell, selected (as the player's key press does); bare hands: an empty one.
        int cell = -1;
        for (int k = 0; k < 10 && cell < 0; ++k) if (ses.inventory.grid()[Inventory::kGrid + k].empty() && k != ses.hotbar.selected) cell = k;
        if (cell < 0) { c.fail("no free hotbar cell"); return true; }
        if (sp.item[0]) ses.inventory.grid()[Inventory::kGrid + cell] = InvSlot{ sp.item, 100 * 20, 0.0f };
        if (sp.offHand[0]) ses.inventory.grid()[Inventory::kOffHand] = InvSlot{ sp.offHand, 100, 0.0f };
        R.cell = cell;
        Hotbar_Press(cell, false);
        R.ph = Ph::Spawn; R.idx = 0; R.f = 0;
        return false;
    }
    if (R.idx > kN) return true;
    const bool control = R.idx == kN;
    const std::string what = std::string(sp.name) + " vs " + (control ? "bare ground" : kInteractables[R.idx].name);
    const Vector3 up = Vector3Normalize(R.at[R.idx]);
    const Vector3 side = Vector3Normalize(Vector3CrossProduct(up, { 0.0f, 0.0f, 1.0f }));
    if (R.ph == Ph::Spawn) {
        gFakeHeld &= ~kRmb;
        GameUI_CloseContainer();
        Telescope_SetActive(c.world, false);
        survival.hasBed = false;
        survival.food = 50.0f; survival.stamina = 100.0f;
        sp.reset(c);
        Hotbar_Press(R.cell, false); // (a used item, e.g. a poured bucket, is taken in hand again)
        // 2.5 m off the part, on the ground (the control: standing where its ground point is, looking at the ground ahead)
        const Vector3 from = control ? Vector3Add(Vector3Scale(up, Vector3Length(R.at[R.idx]) + 1.5f), Vector3Scale(side, -3.0f))
                                     : Vector3Add(Vector3Scale(up, Vector3Length(R.at[R.idx])), Vector3Scale(side, 2.5f));
        player.placeAt(c.world, *c.body, c.body->toWorld(from));
        R.ph = Ph::Wait; R.f = 0;
        return false;
    }
    // The view, every frame: at the part (the control: down at the ground 3 m ahead, away from the others)
    {
        const Vector3 eye = Vector3Add(FeetLocal(c), Vector3Scale(up, 1.6f));
        const Vector3 target = control ? Vector3Add(Vector3Scale(up, Vector3Length(R.at[R.idx]) - 0.2f), Vector3Scale(side, -6.0f)) : R.at[R.idx];
        player.forward = c.body->dirToWorld(Vector3Normalize(Vector3Subtract(target, eye)));
    }
    ++R.f;
    if (R.ph == Ph::Wait) {
        c.waitingFor("on the ground, " + std::string(control ? "the ground aimed, nothing interactable" : "the part aimed") + " (" + what + ")");
        Planetoid* hb = nullptr; Vector3 hl{}, hn{};
        const bool aimOk = control ? (!Interact_RmbConsumed() && TerrainTool_LastHit(hb, hl, hn)) : Interact_RmbConsumed();
        if (R.f > 45 && player.grounded && aimOk && !GameUI_ContainerOpen() && TerrainTool_Cooldown() <= 0.0f) {
            R.fired = R.screen = R.bed = false;
            R.ph = Ph::Hold; R.f = 0;
            gFakeHeld |= kRmb;
        } else if (R.f > 600) { c.check(what + ": aimed", false); return true; }
        return false;
    }
    R.fired = R.fired || sp.fired();
    R.screen = R.screen || GameUI_ContainerOpen();
    R.bed = R.bed || survival.hasBed;
    if (R.ph == Ph::Hold) {
        if (R.f >= sp.holdFrames) { gFakeHeld &= ~kRmb; R.ph = Ph::Settle; R.f = 0; }
        return false;
    }
    if (R.ph == Ph::Settle) {
        if (R.f < 6) return false;
        if (control) {
            c.check(what + ": the item's RMB fires", R.fired);
            c.check(what + ": no interactable screen", !R.screen);
        } else {
            const bool bedCase = kInteractables[R.idx].type == PartType::Bed;
            c.check(what + (bedCase ? ": the bed became the spawn" : ": its screen opened"), bedCase ? R.bed : R.screen);
            c.check(what + ": the item's RMB did not fire", !R.fired);
        }
        R.ph = Ph::Spawn; R.f = 0; ++R.idx;
        if (R.idx > kN) { gScenarioInput = false; gFakeHeld &= ~kRmb; GameUI_CloseContainer(); Telescope_SetActive(c.world, false); return true; }
    }
    return false;
}

}

SCENARIO("input/rmb-telescope") {
    s.tags = { "input", "rmb" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    static const Spec sp = { "telescope", "telescope", "", 4, [](ScenarioCtx&) {}, [] { return telescope.active; } };
    s.step("RMB with a telescope in hand: an aimed chest, furnace, campfire or bed takes it; bare ground raises the scope", 90.0, [](ScenarioCtx& c) { return Drive(c, sp); });
}

SCENARIO("input/rmb-bucket") {
    s.tags = { "input", "rmb" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    static int before = 0;
    static const Spec sp = { "water bucket", "water_bucket", "", 4,
        [](ScenarioCtx&) { before = (int)LocalSession().inventory.count("water_bucket"); if (before < 5) { LocalSession().inventory.add("water_bucket", 5 - before); before = 5; } },
        [] { return (int)LocalSession().inventory.count("water_bucket") < before; } };
    s.step("RMB with a full bucket: an aimed interactable takes it (nothing is poured); bare ground pours it", 90.0, [](ScenarioCtx& c) { return Drive(c, sp); });
}

SCENARIO("input/rmb-food") {
    s.tags = { "input", "rmb" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    static int eaten0 = 0;
    static const Spec sp = { "food", "cooked_meat", "", 100,
        [](ScenarioCtx&) { eaten0 = LocalSession().hotbar.eatenCount; },
        [] { return LocalSession().hotbar.eatProgress > 0.0f || LocalSession().hotbar.eatenCount > eaten0; } };
    s.step("RMB held with food in hand: an aimed interactable takes it (no bite starts); bare ground eats", 120.0, [](ScenarioCtx& c) { return Drive(c, sp); });
}

SCENARIO("input/rmb-shield") {
    s.tags = { "input", "rmb" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    static const Spec sp = { "shield", "", "wooden_shield", 20, [](ScenarioCtx&) {}, [] { return survival.blocking; } };
    s.step("RMB with a shield in the off hand: an aimed interactable takes it (no block); bare ground raises it", 90.0, [](ScenarioCtx& c) { return Drive(c, sp); });
}

SCENARIO("input/rmb-ground") {
    s.tags = { "input", "rmb" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    static int before = 0;
    static const Spec sp = { "dirt", "dirt", "", 4,
        [](ScenarioCtx&) { before = (int)LocalSession().inventory.count("dirt"); },
        [] { return (int)LocalSession().inventory.count("dirt") < before; } };
    s.step("RMB with dirt in hand: an aimed interactable takes it (nothing is placed); bare ground takes the dirt", 90.0, [](ScenarioCtx& c) { return Drive(c, sp); });
}

SCENARIO("input/rmb-torch") {
    s.tags = { "input", "rmb" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    static int before = 0;
    static const Spec sp = { "torch", "torch", "", 4,
        [](ScenarioCtx&) { before = (int)LocalSession().inventory.count("torch"); },
        [] { return (int)LocalSession().inventory.count("torch") < before; } };
    s.step("RMB with a torch in hand: an aimed interactable takes it (nothing is set down); bare ground places it", 90.0, [](ScenarioCtx& c) { return Drive(c, sp); });
}
