// combat/melee (#1216, combat step 2): light and heavy swings in phases against a training dummy (a hornhide held in place,
// side-on, 1.3 m ahead), on the 60 Hz sim clock, with scripted LMB presses through the real input path.
//   light: a tap. Wind-up 250 ms, then the hit window (the blade sweeps the arc and meets the dummy once), back to idle
//          650 ms after the start (bare hands fight as "tools": 250/100/300 ms), 2 damage, 6 stamina.
//   heavy: LMB held through the wind-up: wind-up x1.6 (400 ms), damage x1.8 (3.6), stamina x1.7 (10.2).
//   chain: a second press in the first light's hit window: the second light starts a third of the way into the recovery
//          (100 ms after the hit window ends) and hits again; two swings, two hits.
//   lock-on: middle mouse with the dummy 5 m off to the right: the view eases onto it; pressed again, it lets go.
#include "scenario.h"
#include "animals.h"
#include "melee.h"
#include "player.h"
#include "playerinput.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
struct St {
    uint64_t net = 0; Vector3 at{}, face{};
    float hp0 = 0.0f, stamina0 = 0.0f; int hits0 = 0, swings0 = 0, heavies0 = 0;
    double tStart = -1, tHit = -1, tContact = -1, tIdle = -1, tSecond = -1, tHit2 = -1;
    long pressFrame = -1; int presses = 0;
} S;

constexpr uint64_t kLmb = 1ull << (int)Act::Primary, kMmb = 1ull << (int)Act::LockOn;

void Stand(ScenarioCtx& c) { player.flyMode = false; player.spawnOn(c.world, *c.body, c.spawnDirWorld()); }

Vector3 EyeLocal(ScenarioCtx& c) { return c.body->toLocal(c.world.camera.position); }

// The ground point d m ahead of the eye (deg to the right), and the view's flat ahead and right (body-local).
Vector3 GroundAhead(ScenarioCtx& c, float d, float deg, Vector3& flat, Vector3& right) {
    const Vector3 eye = EyeLocal(c), up = Vector3Normalize(eye);
    const Vector3 f = c.body->dirToLocal(player.forward);
    flat = Vector3Normalize(Vector3Subtract(f, Vector3Scale(up, Vector3DotProduct(f, up))));
    right = Vector3Normalize(Vector3CrossProduct(flat, up));
    const Vector3 p = Vector3Add(eye, Vector3Add(Vector3Scale(flat, d * cosf(deg * DEG2RAD)), Vector3Scale(right, d * sinf(deg * DEG2RAD))));
    const Vector3 dir = Vector3Normalize(p);
    return Vector3Scale(dir, c.field->radius + c.field->elevation(dir));
}

// The dummy's index now, by its stable id (an index shifts when another animal goes).
int Idx() { return S.net ? Animals_FindNet(S.net) : -1; }
void Pin() { if (Idx() >= 0) Animals_TestPin(Idx(), S.at, S.face); }

Vector3 DummyCentre() { return Vector3Add(S.at, Vector3Scale(Vector3Normalize(S.at), 0.6f)); } // half a hornhide's 1.2 m

// Look at the dummy's middle (each frame: the blade sweeps from the view).
void LookAtDummy(ScenarioCtx& c) {
    const Vector3 to = Vector3Subtract(c.body->toWorld(DummyCentre()), c.world.camera.position);
    if (Vector3Length(to) > 0.1f) player.forward = Vector3Normalize(to);
}

void ResetTrack() {
    const MeleeState& m = Animals_Melee();
    S.hp0 = Animals_HealthOf(Idx()); S.stamina0 = survival.stamina; S.hits0 = Animals_MeleeHits();
    S.swings0 = m.swings; S.heavies0 = m.heavies;
    S.tStart = S.tHit = S.tContact = S.tIdle = S.tSecond = S.tHit2 = -1; S.pressFrame = -1; S.presses = 0;
}

void Track(ScenarioCtx& c) {
    const MeleeState& m = Animals_Melee();
    if (S.tStart < 0 && m.swings > S.swings0) S.tStart = c.T;
    if (S.tStart >= 0 && S.tHit < 0 && m.phase == MeleePhase::Hit) S.tHit = c.T;
    if (S.tContact < 0 && Animals_MeleeHits() > S.hits0) S.tContact = c.T;
    if (S.tSecond < 0 && m.chain == 1 && m.phase == MeleePhase::WindUp) S.tSecond = c.T;
    if (S.tSecond >= 0 && S.tHit2 < 0 && m.phase == MeleePhase::Hit) S.tHit2 = c.T;
    if (S.tHit >= 0 && S.tIdle < 0 && m.phase == MeleePhase::Idle) S.tIdle = c.T;
}

// One frame of LMB down (a tap).
void Tap(ScenarioCtx& c) { gFakeHeld |= kLmb; S.pressFrame = c.frame; ++S.presses; }
void ReleaseAfterTap(ScenarioCtx& c) { if (S.pressFrame >= 0 && c.frame > S.pressFrame && (gFakeHeld & kLmb) && S.presses > 0) gFakeHeld &= ~kLmb; }
}

SCENARIO("combat/melee") {
    s.tags = { "combat" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.args = { "--day" }; // the dummy and the swing readable in the recording
    constexpr double kStep = 1.0 / 60.0;

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Stand(c);
        return true;
    });
    s.step("a dummy ahead", 5.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false; // standing
        Vector3 flat, right;
        S.at = GroundAhead(c, 1.3f, 0.0f, flat, right);
        S.face = right; // side-on
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "hornhide", S.at, 1, false, 0.01f);
        if (!c.check("a hornhide for the dummy", Animals_Count() > before)) { c.fail("no dummy"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        gScenarioInput = true; // scripted LMB and middle mouse through the real input path
        Pin(); LookAtDummy(c);
        return true;
    });
    s.step("light: a tap", 5.0, [](ScenarioCtx& c) {
        Pin(); LookAtDummy(c);
        if (c.t < 0.3) return false; // the view settles on it
        if (S.pressFrame < 0) { ResetTrack(); Tap(c); c.cue(); return false; }
        ReleaseAfterTap(c); Track(c);
        if (S.tIdle < 0) { const MeleeState& m = Animals_Melee(); c.waitingFor("the light swing to finish: swings " + std::to_string(m.swings) + " phase " + std::to_string((int)m.phase) + " aiming " + std::to_string(Animals_Aiming()) + " input " + std::to_string(PlayerInput_Frame().enabled) + " held " + std::to_string(PlayerInput_Frame().held) + " fly " + std::to_string(player.flyMode)); return false; }
        c.check("the swing started", S.tStart >= 0);
        c.near("light wind-up (s)", S.tHit - S.tStart, 0.25, 1.5 * kStep);
        c.check("contact in the hit window", S.tContact >= S.tHit - 1e-6 && S.tContact <= S.tHit + 0.10 + 1.5 * kStep,
                "contact " + std::to_string(S.tContact - S.tStart) + " s after the start");
        c.near("light swing, start to idle (s)", S.tIdle - S.tStart, 0.65, 1.5 * kStep);
        c.eq("one hit per target per swing", Animals_MeleeHits() - S.hits0, 1);
        c.near("light damage (bare hands)", S.hp0 - Animals_HealthOf(Idx()), 2.0, 1e-3);
        c.near("light stamina", S.stamina0 - survival.stamina, 6.0, 0.5);
        return true;
    });
    s.step("heavy: held through the wind-up", 5.0, [](ScenarioCtx& c) {
        Pin(); LookAtDummy(c);
        if (c.t < 1.2) { S.pressFrame = -1; return false; } // stamina regen waits; the numbers stay clean
        if (S.pressFrame < 0) { ResetTrack(); gFakeHeld |= kLmb; S.pressFrame = c.frame; return false; }
        Track(c);
        if (S.tHit >= 0) gFakeHeld &= ~kLmb;
        if (S.tIdle < 0) { c.waitingFor("the heavy swing to finish"); return false; }
        c.eq("it became a heavy", Animals_Melee().heavies - S.heavies0, 1);
        c.near("heavy wind-up (s)", S.tHit - S.tStart, 0.40, 1.5 * kStep);
        c.near("heavy damage (x1.8)", S.hp0 - Animals_HealthOf(Idx()), 3.6, 1e-3);
        c.near("heavy stamina (x1.7)", S.stamina0 - survival.stamina, 10.2, 0.5);
        return true;
    });
    s.step("two lights chain", 5.0, [](ScenarioCtx& c) {
        Pin(); LookAtDummy(c);
        if (c.t < 1.2) { S.pressFrame = -1; return false; }
        if (S.pressFrame < 0) { ResetTrack(); Tap(c); return false; }
        ReleaseAfterTap(c); Track(c);
        if (S.tHit >= 0 && S.presses == 1 && !(gFakeHeld & kLmb)) Tap(c); // the second press, in the hit window
        if (S.tIdle < 0 || S.tSecond < 0 || S.tIdle < S.tHit2) { c.waitingFor("both swings to finish"); return false; }
        c.eq("two swings", Animals_Melee().swings - S.swings0, 2);
        c.eq("two hits", Animals_MeleeHits() - S.hits0, 2);
        c.near("the second starts a third into the recovery (s after the first hit window)", S.tSecond - S.tHit, 0.10 + 0.10, 1.5 * kStep);
        return true;
    });
    s.step("lock-on: middle mouse", 6.0, [](ScenarioCtx& c) {
        static int stage = 0;
        static double t0 = 0.0;
        if (stage == 0) {
            // The dummy 5 m out, 25 degrees right; the view level and straight ahead.
            Vector3 flat, right;
            S.at = GroundAhead(c, 5.0f, 25.0f, flat, right);
            S.face = right;
            player.forward = c.body->dirToWorld(flat);
            Pin();
            gFakeHeld |= kMmb; stage = 1; t0 = c.T;
            return false;
        }
        Pin();
        if (stage == 1) { gFakeHeld &= ~kMmb; stage = 2; return false; }
        if (stage == 2 && c.T - t0 < 0.8) return false;
        if (stage == 2) {
            c.check("locked on", Animals_LockOn() != 0);
            const Vector3 to = Vector3Normalize(Vector3Subtract(c.body->toWorld(DummyCentre()), c.world.camera.position));
            const double deg = acos(std::fmin(1.0, (double)Vector3DotProduct(Vector3Normalize(player.forward), to))) * RAD2DEG;
            c.lt("the view eased onto it (degrees off)", deg, 6.0);
            gFakeHeld |= kMmb; stage = 3;
            return false;
        }
        if (stage == 3) { gFakeHeld &= ~kMmb; stage = 4; return false; }
        c.check("pressed again, it lets go", Animals_LockOn() == 0);
        gFakeHeld = 0; gScenarioInput = false;
        return true;
    });
}
