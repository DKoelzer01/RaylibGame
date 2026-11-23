// ui/frame-cap (#1624, the user): the frame cap setting is a slider (Uncapped, then 1-240) with a number box that takes more;
// fpscap 0 is uncapped (it was "Monitor", which read 60 on a 144 Hz screen). Applied as closing Settings does: 0 leaves
// the frame pacer (#1609) uncapped, 500 (past the slider, typed in the box) caps it at 500, 144 at 144.
#include "scenario.h"
#include "framepacer.h"
#include "headless.h"
#include "menus.h"
#include "world.h"

SCENARIO("ui/frame-cap") {
    s.tags = { "world", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("0 is uncapped, 500 and 144 cap the pacer", 5.0, [](ScenarioCtx& c) {
        if (Headless()) { c.note("a headless run paces at the sim rate whatever the setting: nothing to check"); return true; }
        const int saved = settings.fpsCap;
        const bool vsync = settings.vsync;
        settings.vsync = false;
        for (const int cap : { 0, 500, 144 }) {
            settings.fpsCap = cap;
            Settings_Apply();
            c.eq(TextFormat("fpscap %d: the pacer's target", cap), FramePacer_Target(), cap);
        }
        c.gt("the monitor's refresh rate reads (Hz)", Display_MonitorHz(), 0);
        settings.fpsCap = saved; settings.vsync = vsync;
        Settings_Apply();
        return true;
    });
}
