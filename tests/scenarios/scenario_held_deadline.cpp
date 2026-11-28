// harness/held-deadline (#2119): a step that waits while the world streams doesn't spend its deadline on it. The headless
// sim runs many times real time while chunk builds and loading screens take wall time, so on a loaded runner a sim deadline
// ran out before the world was idle again (11 of the 20 flaky scenarios of #2107). Here: /tp across the sandbox, then wait
// for the world on a deadline of 0.5 s sim, far less than the streaming takes: before #2119 that step failed on any machine.
#include "scenario.h"
#include "commandline.h"
#include "playersession.h"
#include "terraintool.h"
#include "world.h"
#include <string>

namespace {
bool gSent = false;
}

SCENARIO("harness/held-deadline") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("settle: the start built", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        return c.T >= 1.0 && WorldIdle();
    });
    s.step("/tp across the world, wait for it on a 0.5 s deadline", 0.5, [](ScenarioCtx& c) {
        if (!gSent) {
            terrainTool.creative = true; // (/tp is creative's)
            c.note(("/tp -35 150 -> " + CommandLine_Run("/tp -35 150")).c_str());
            gSent = true;
        }
        c.waitingFor("the world streamed in at the far side (an idle world)");
        if (!WorldIdle() || !LocalSession().player.nearest) return false;
        c.note(TextFormat("idle after %.2f s sim, %.2f s wall; the deadline counted %.2f s of it", c.t, c.wall, c.clock));
        c.check("the world streamed for a while (the step had something to wait for)", c.t > c.clock + 0.3,
                TextFormat("%.2f s sim, %.2f s counted", c.t, c.clock));
        c.lt("the deadline counted only the frames outside the streaming (s)", c.clock, 0.5);
        return true;
    });
}
