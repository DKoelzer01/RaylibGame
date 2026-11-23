// ui/catalog-slots (#887, the user: "Catalog slots are larger than inventory slots"): the Catalog drew a fixed 48 px cell at
// the interface scale while the Inventory's shrinks to fit the screen (#1148). At two interface sizes the Inventory page and
// then the Catalog page are opened; the slot size each draws with must be the same.
#include "scenario.h"
#include "gameui.h"
#include "menus.h"
#include "terraintool.h"
#include "world.h"
#include <cmath>

SCENARIO("ui/catalog-slots") {
    s.tags = { "world", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.step("settle, creative", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        terrainTool.creative = true; // (the Catalog is creative's)
        return c.T >= 1.0 && WorldIdle();
    });
    for (const float size : { 1.0f, 1.5f }) {
        s.step(TextFormat("interface %.2f: the Inventory's and the Catalog's slots the same size", size), 6.0, [size](ScenarioCtx& c) {
            static float inv = 0.0f;
            if (c.t < 0.05) { settings.uiScale = size; Settings_Apply(); GameUI_ShowScreen(GameScreen::Inventory); return false; }
            if (c.t < 0.5) return false;
            float a = 0.0f, b = 0.0f;
            if (c.t < 0.55) { GameUI_TestDrawnCells(a, b); inv = a; GameUI_ShowScreen(GameScreen::Catalog); return false; }
            if (c.t < 1.0) return false;
            GameUI_TestDrawnCells(a, b);
            c.gt(TextFormat("interface %.2f: the Inventory drew its slots (px)", size), inv, 1.0f);
            c.near(TextFormat("interface %.2f: the Catalog's slot = the Inventory's (px)", size), b, inv, 0.01f);
            c.grab(TextFormat("catalog-%.0f", size * 100.0f));
            return true;
        });
    }
}
