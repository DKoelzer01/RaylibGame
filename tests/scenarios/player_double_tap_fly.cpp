// player/double-tap-fly (#1562): in creative a double-tap of Space toggles flying (as F does); a single tap still jumps. Out of
// creative a double-tap is two jumps. Space is pressed through the real input path (gFakeHeld), one sim step a frame; the
// player's own view follows it (c.camera would put it in flight).
#include "scenario.h"
#include "player.h"
#include "playerinput.h"
#include "playersession.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
constexpr uint64_t kSpace = 1ull << (int)Act::Jump;
// Presses Space at frames f0, and (twice) f0 + 6: 2 frames down each (0.1 s apart at 60 Hz). Returns true when done.
bool Taps(ScenarioCtx& c, long f0, int taps) {
    const long f = c.frame - f0;
    const bool down = (f >= 0 && f < 2) || (taps == 2 && f >= 6 && f < 8);
    if (down) gFakeHeld |= kSpace; else gFakeHeld &= ~kSpace;
    return f >= 10;
}
long f0 = 0;
bool airborne = false;
}

SCENARIO("player/double-tap-fly") {
    s.tags = { "skills" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    s.step("settle on the ground, creative", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play, an idle world, on the ground");
        terrainTool.creative = true;
        player.flyMode = false;
        if (c.T < 1.0 || !WorldIdle() || !player.grounded) return false;
        f0 = c.frame + 5;
        return true;
    });
    s.step("creative: one tap jumps, no flying", 5.0, [](ScenarioCtx& c) {
        airborne = airborne || !player.grounded;
        if (!Taps(c, f0, 1) || c.frame < f0 + 45) return false;
        c.check("a single tap jumped", airborne);
        c.check("a single tap doesn't fly", !player.flyMode);
        f0 = c.frame + 30; // (landed by then)
        return true;
    });
    s.step("creative: a double-tap flies", 5.0, [](ScenarioCtx& c) {
        if (!Taps(c, f0, 2)) return false;
        c.check("a double-tap toggles flying on", player.flyMode);
        f0 = c.frame + 20;
        return true;
    });
    s.step("creative: another double-tap lands", 5.0, [](ScenarioCtx& c) {
        if (!Taps(c, f0, 2)) return false;
        c.check("a double-tap toggles flying off", !player.flyMode);
        f0 = c.frame + 90; // (down on the ground again)
        return true;
    });
    s.step("survival: a double-tap doesn't fly", 8.0, [](ScenarioCtx& c) {
        terrainTool.creative = false;
        if (!Taps(c, f0, 2) || c.frame < f0 + 20) return false;
        c.check("survival: a double-tap doesn't fly", !player.flyMode);
        gFakeHeld &= ~kSpace;
        return true;
    });
}
