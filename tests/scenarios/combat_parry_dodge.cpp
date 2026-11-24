// combat/parry-dodge (#1218, combat step 4): parry timing, a dodge through an attack, and the winded state.
//   parry:  the wooden (round) shield raised (RMB through the real input path) a few steps before a hit from in front: stopped
//           whole, for no stamina. Held past the 330 ms window (#1218: was 150): an ordinary block (a quarter of the hit lands).
//   dodge:  a double-tap of D: a 0.45 s step of ~3 m to the right for 20 stamina; a hit inside its first 0.2 s passes through.
//   winded: run out of stamina: winded 1.5 s (walk x0.7), and a double-tap does not dodge.
#include "scenario.h"
#include "combat.h"
#include "defense.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "playersession.h"
#include "statusbar.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>

namespace {
struct St { Vector3 start{}; float stamina0 = 0.0f; int dodges0 = 0; long frame0 = 0; bool hitDuring = false; double tapAt = -1.0; } S;
constexpr uint64_t kRmb = 1ull << (int)Act::Secondary, kRight = 1ull << (int)Act::Right;

Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }

// A hit from 2 m in front of the player.
DamageResult HitFromFront(float amount, DamageKind kind = DamageKind::Melee) {
    DamageEvent e; e.amount = amount; e.kind = kind; e.hasSource = true;
    const Vector3 flat = Vector3Normalize(Vector3Subtract(player.forward, Vector3Scale(player.up, Vector3DotProduct(player.forward, player.up))));
    e.source = player.position + DVec3(Vector3Scale(flat, 2.0f));
    return Combat_Damage(LocalSession(), e);
}
}

SCENARIO("combat/parry-dodge") {
    s.tags = { "combat" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;

    s.step("settle, a shield in the off hand", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        InvSlot& off = LocalSession().inventory.grid()[Inventory::kOffHand];
        off.id = "wooden_shield"; off.centi = 100;
        // An empty hand: a tool or filler in the selected cell keeps RMB for itself (a shovel tills), and the shield won't block.
        for (int k = 0; k < 10; ++k) if (LocalSession().inventory.grid()[Inventory::kGrid + k].id.empty()) { LocalSession().hotbar.selected = k; break; }
        gScenarioInput = true;
        return true;
    });
    s.step("parry: raised just before the hit", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false; // standing
        if (!(gFakeHeld & kRmb)) { gFakeHeld |= kRmb; S.frame0 = c.frame; return false; }
        if (!survival.parrying) { c.waitingFor("the parry window to open (#1776: the press): rmbBlocks " + std::to_string(LocalSession().terrainTool.rmbBlocks) + " off " + LocalSession().inventory.offHand().id + " rmb held " + std::to_string((PlayerInput_Step().held >> (int)Act::Secondary) & 1) + " stamina " + std::to_string(survival.stamina)); return false; }
        if (c.frame < S.frame0 + 4) return false; // ~50 ms after the raise
        survival.health = 100.0f; survival.stamina = 100.0f;
        const int parries0 = survival.parries;
        c.check("the PARRY bar shows the window (#1777)", StatusBar_Active("parry"), TextFormat("%.2f left", StatusBar_Remaining("parry")));
        const DamageResult r = HitFromFront(20.0f);
        c.check("parried", r.parried);
        c.near("no damage", 100.0f - survival.health, 0.0, 1e-4);
        c.near("no stamina", 100.0f - survival.stamina, 0.0, 1e-4);
        c.eq("counted", survival.parries - parries0, 1);
        c.check("the bar goes when the parry lands", !StatusBar_Active("parry"));
        return true;
    });
    s.step("held past the window: an ordinary block", 3.0, [](ScenarioCtx& c) {
        if (c.t < 0.4) return false;
        survival.health = 100.0f; survival.stamina = 100.0f;
        const DamageResult r = HitFromFront(20.0f);
        c.check("not a parry", !r.parried);
        c.check("blocked", r.blocked);
        c.near("a quarter lands", r.dealt, 5.0, 1e-3);
        gFakeHeld &= ~kRmb;
        return true;
    });
    // (#1776, user) A tap parries (RMB already up again); a parry turns a blow whose live span overlaps its window, even
    // just after the window closed; a hit well after it lands; a projectile is parried too (no damage). Each tap waits out
    // the recovery first. `when` is the hit's time after the tap, as a share of the window plus seconds.
    struct TapCase { const char* name; float windowShare; float extra; DamageKind kind; bool parried; };
    static const TapCase kTaps[] = {
        { "a tap parries, RMB already up", 0.5f, 0.0f, DamageKind::Melee, true },
        { "the window just closed, the blow still live", 1.0f, 0.08f, DamageKind::Melee, true },
        { "a hit well after the window lands", 1.0f, 0.5f, DamageKind::Melee, false },
        { "a projectile is parried", 0.5f, 0.0f, DamageKind::Projectile, true },
    };
    for (const TapCase& tc : kTaps)
        s.step(tc.name, 5.0, [tc](ScenarioCtx& c) {
            if (c.t < 1.0) { S.tapAt = -1.0; return false; } // the recovery and stamina back
            if (S.tapAt < 0.0) { gFakeHeld |= kRmb; S.tapAt = SimTime(); S.frame0 = c.frame; return false; }
            if (c.frame == S.frame0 + 3) gFakeHeld &= ~kRmb; // a tap: up again after 3 frames
            if (SimTime() - S.tapAt < survival.shieldParryWindow * tc.windowShare + tc.extra) return false;
            gFakeHeld &= ~kRmb;
            survival.health = 100.0f; survival.stamina = 100.0f;
            const DamageResult r = HitFromFront(20.0f, tc.kind);
            c.check("parried as expected", r.parried == tc.parried, TextFormat("parried %d, %.2f s after the tap, window %.2f s", (int)r.parried, SimTime() - S.tapAt, survival.shieldParryWindow));
            c.check("no block (RMB is up)", !r.blocked);
            c.near("damage", 100.0f - survival.health, tc.parried ? 0.0 : 20.0, 1e-3);
            return true;
        });
    s.step("dodge: double-tap D, a hit passes through", 5.0, [](ScenarioCtx& c) {
        if (c.t < 1.2) { S.dodges0 = survival.dodges; S.hitDuring = false; return false; } // the block down, stamina back
        const long f = c.frame;
        static int phase = 0;
        if (phase == 0) { survival.stamina = 100.0f; S.start = FeetLocal(c); gFakeHeld |= kRight; S.frame0 = f; phase = 1; return false; }
        if (phase == 1 && f >= S.frame0 + 2) { gFakeHeld &= ~kRight; phase = 2; return false; }
        if (phase == 2 && f >= S.frame0 + 5) { gFakeHeld |= kRight; phase = 3; return false; }
        if (phase == 3 && survival.dodges > S.dodges0 && !S.hitDuring) {
            gFakeHeld &= ~kRight;
            S.hitDuring = true;
            survival.health = 100.0f;
            const DamageResult r = HitFromFront(20.0f);
            c.check("a hit inside the dodge passes", r.dodged && survival.health == 100.0f);
            S.frame0 = f; phase = 4;
            return false;
        }
        if (phase == 3 && f > S.frame0 + 30) { c.check("the double-tap dodged", false); gFakeHeld &= ~kRight; return true; }
        if (phase == 4 && f < S.frame0 + 40) return false;
        if (phase != 4) return false;
        c.eq("one dodge", survival.dodges - S.dodges0, 1);
        c.near("its stamina", 100.0f - survival.stamina, 20.0, 1.0);
        const Vector3 moved = Vector3Subtract(FeetLocal(c), S.start);
        const Vector3 up = Vector3Normalize(S.start);
        const float level = Vector3Length(Vector3Subtract(moved, Vector3Scale(up, Vector3DotProduct(moved, up))));
        c.between("dodged about 3 m (with the walk in the taps)", level, 2.4, 4.0);
        return true;
    });
    s.step("winded: out of stamina, no dodge", 5.0, [](ScenarioCtx& c) {
        static int phase = 0;
        static long f0 = 0;
        if (c.t < 0.6) return false;
        if (phase == 0) { survival.stamina = 0.0f; survival.staminaUsedAt = 1e9; phase = 1; f0 = c.frame; return false; }
        if (phase == 1 && c.frame < f0 + 3) return false;
        if (phase == 1) {
            c.gt("winded", survival.windedLeft, 1.0);
            c.near("walk x0.7", player.statusSpeed, 0.7, 0.05);
            S.dodges0 = survival.dodges;
            survival.stamina = 100.0f; // stamina back, still winded: the winded state alone stops the dodge
            gFakeHeld |= kRight; phase = 2; f0 = c.frame; return false;
        }
        if (phase == 2 && c.frame >= f0 + 2) { gFakeHeld &= ~kRight; phase = 3; return false; }
        if (phase == 3 && c.frame >= f0 + 5) { gFakeHeld |= kRight; phase = 4; return false; }
        if (phase == 4 && c.frame >= f0 + 8) { gFakeHeld &= ~kRight; phase = 5; return false; }
        if (phase != 5 || c.frame < f0 + 12) return false;
        c.eq("no dodge while winded", survival.dodges - S.dodges0, 0);
        survival.staminaUsedAt = -1e9; // regen again
        gFakeHeld = 0; gScenarioInput = false;
        return true;
    });
}
