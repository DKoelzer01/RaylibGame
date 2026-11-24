// ui/status-bar (#1777): the status-effect widget on the left edge. Two timed effects (1 s and 2 s) are added; their bars
// must drain with sim time and each entry must be gone once its time is up. Re-adding an effect restarts it full in its own
// slot; removing one ends it early. A grab of the HUD with both entries showing (RAYLIBGAME_SLIDE_SHOTS).
#include "scenario.h"
#include "statusbar.h"
#include "world.h"
#include <cmath>

namespace {
double t0 = 0.0; // sim time when the two effects were added
}

SCENARIO("ui/status-bar") {
    s.tags = { "world", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day", "--record-hud", "1" }; // the recording shows the HUD: the widget is what this scenario is about
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("two effects: 1 s and 2 s", 1.0, [](ScenarioCtx& c) {
        StatusBar_Clear();
        StatusBar_Add("test_short", 1.0f, "SHORT", GOLD);
        StatusBar_Add("test_long", 2.0f, "LONG", SKYBLUE, "telescope"); // with an item icon (one every build has)
        t0 = c.T;
        c.eq("two entries showing", StatusBar_Count(), 2);
        c.check("both start full", StatusBar_Remaining("test_short") > 0.97f && StatusBar_Remaining("test_long") > 0.97f);
        return true;
    });
    s.step("half a second in: both draining", 2.0, [](ScenarioCtx& c) {
        if (c.T - t0 < 0.5) return false;
        const double el = c.T - t0;
        c.near("the 1 s bar drains with sim time", StatusBar_Remaining("test_short"), 1.0 - el / 1.0, 0.05);
        c.near("the 2 s bar drains at half that rate", StatusBar_Remaining("test_long"), 1.0 - el / 2.0, 0.05);
        c.grab("two-effects"); // (a picture of the widget)
        return true;
    });
    s.step("past 1 s: the short one is gone", 2.0, [](ScenarioCtx& c) {
        if (c.T - t0 < 1.05) return false;
        c.check("the 1 s effect is removed on time", !StatusBar_Active("test_short"));
        c.check("the 2 s effect still shows", StatusBar_Active("test_long"));
        c.eq("one entry showing", StatusBar_Count(), 1);
        return true;
    });
    s.step("past 2 s: none left", 3.0, [](ScenarioCtx& c) {
        if (c.T - t0 < 2.05) return false;
        c.eq("no entries showing", StatusBar_Count(), 0);
        c.eq("an ended effect has nothing left", (double)StatusBar_Remaining("test_long"), 0.0);
        return true;
    });
    s.step("re-add restarts in place; remove ends early", 2.0, [](ScenarioCtx& c) {
        static double ta = -1.0;
        if (ta < 0.0) { StatusBar_Add("a", 2.0f, "A", GOLD); StatusBar_Add("b", 2.0f, "B", GOLD); ta = c.T; return false; }
        if (c.T - ta < 0.5) return false;
        StatusBar_Add("a", 2.0f, "A", GOLD); // restarted: full again
        c.check("a re-added effect starts full again", StatusBar_Remaining("a") > 0.97f);
        c.check("the other keeps draining", StatusBar_Remaining("b") < 0.8f);
        c.eq("still two entries (no duplicate)", StatusBar_Count(), 2);
        StatusBar_Remove("a");
        c.check("a removed effect is gone at once", !StatusBar_Active("a") && StatusBar_Count() == 1);
        StatusBar_Clear();
        return true;
    });
}
