// ui/interface-size (#1153): changing Settings > Interface size away from Auto and closing Settings closed the game once in
// a playtest (not reproduced since). Closing Settings applies the settings (Settings_Apply; the file save is left out here).
// Each size in turn - Auto, 75 %, 100 %, 125 %, 150 %, back to Auto - is applied and then held for a second with the game
// window open and a second with it closed (the ImGui style and fonts rebuild at the new scale on the next frame). The game
// must still be running and report the scale it was given.
#include "scenario.h"
#include "gameui.h"
#include "menus.h"
#include "world.h"
#include <cmath>

namespace {
const float kSizes[] = { 0.0f, 0.75f, 1.0f, 1.25f, 1.5f, 0.0f };
constexpr int kCount = sizeof(kSizes) / sizeof(kSizes[0]);
}

SCENARIO("ui/interface-size") {
    s.tags = { "world", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    for (int k = 0; k < kCount; ++k) {
        s.step(TextFormat("size %d of %d (%.2f; 0 = Auto)", k + 1, kCount, kSizes[k]), 6.0, [k](ScenarioCtx& c) {
            if (c.t < 0.05) { settings.uiScale = kSizes[k]; Settings_Apply(); if (!GameUI_InventoryOpen()) GameUI_ToggleInventory(); }
            if (c.t > 0.6) c.grab(TextFormat("size%d", k)); // the window at this size (RAYLIBGAME_SLIDE_SHOTS: a picture to look at)
            if (c.t > 1.0 && c.t < 1.1 && GameUI_InventoryOpen()) GameUI_ToggleInventory(); // a second open, then closed
            if (c.t < 2.0) return false;
            const float ui = GameUI_Scale();
            static float autoScale = 0.0f; // the scale at Auto (the first size): the screen height's base (GameUI_Scale)
            if (k == 0) autoScale = ui;
            c.check(TextFormat("still running at size %.2f, the interface scale %.2f", kSizes[k], ui), std::isfinite(ui) && ui > 0.1f);
            c.near(TextFormat("the scale = Auto's x %.2f", kSizes[k] > 0.0f ? kSizes[k] : 1.0f), ui, autoScale * (kSizes[k] > 0.0f ? kSizes[k] : 1.0f), 0.01f);
            return true;
        });
    }
}
