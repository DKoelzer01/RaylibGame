// visual/held-items and visual/held-items-lamp (#1588): the hand and what it holds, on the sandbox: a composed bronze pickaxe
// and a steel sword, each in hand (hotbar keys, the real input path) with the player standing, a screenshot of each (the
// HUD and the viewmodel are in shots). The playtest saw every held item as a black silhouette in caves and at night: the
// hand took only the sky, the sun and the headlamp, not the lamps lighting the world around it.
//   held-items:       by day (--day): lit by the sun.
//   held-items-lamp:  at night with the light at the eye (RAYLIBGAME_CAMERA_LIGHT, as a lamp-lit cave): lit by the lamp.
// Checks the viewmodel's light value (finite, above 0.3) and keeps the shots for a person to judge.
#include "scenario.h"
#include "composed.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "terraintool.h"
#include "viewmodel.h"
#include "world.h"
#include <cmath>
#include <string>

namespace {
void Stand(ScenarioCtx& c) { player.flyMode = false; player.spawnOn(c.world, *c.body, c.spawnDirWorld()); }

void Build(ScenarioDef& s, bool lamp) {
    s.tags = { "visual" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 40.0;
    if (lamp) s.env = { { "RAYLIBGAME_CAMERA_LIGHT", "1" } }; // night (the sandbox's start), lit only by the light at the eye
    else s.args = { "--day" };                                  // a mid-morning sun over the start
    s.step("stand", 15.0, [](ScenarioCtx& c) {
        if (c.t < 0.05) { Stand(c); terrainTool.creative = false; }
        c.waitingFor("an idle world");
        return c.t > 1.0 && WorldIdle();
    });
    s.step("give and bind", 1.0, [](ScenarioCtx& c) {
        inventory.addInstance("tool_pickaxe", Composed_MakeData("bronze"));
        inventory.addInstance("weapon_sword", Composed_MakeData("steel"));
        Hotbar_Bind(2, "tool_pickaxe");
        Hotbar_Bind(3, "weapon_sword");
        gScenarioInput = true;
        gFakeSlot = 2;
        return true;
    });
    for (int k = 0; k < 2; ++k) {
        const std::string what = k == 0 ? "pickaxe" : "sword";
        s.step("hold the " + what, 5.0, [k, what, lamp](ScenarioCtx& c) {
            player.flyMode = false;
            if (c.t < 1.5) return false; // selected, raised, and the light eased in
            c.check("the " + what + " in hand", terrainTool.heldTool == (k == 0 ? "tool_pickaxe" : "weapon_sword"), "held '" + terrainTool.heldTool + "'");
            const float l = Viewmodel_Light();
            c.metric("viewmodel light", l, "");
            c.check("the hand's light is a number", std::isfinite(l), std::to_string(l));
            c.check(lamp ? "lit by the lamp at night (> 0.3)" : "lit by day (> 0.3)", l > 0.3f, std::to_string(l));
            c.shot(what, 100.0);
            if (k == 0) gFakeSlot = 3;
            return true;
        });
    }
    s.step("done", 1.0, [lamp](ScenarioCtx& c) { gScenarioInput = false; c.keepCapture(lamp ? "held items by a lamp at night (#1588)" : "held items by day (#1588)"); return true; });
}
}

SCENARIO("visual/held-items") { Build(s, false); }
SCENARIO("visual/held-items-lamp") { Build(s, true); }
