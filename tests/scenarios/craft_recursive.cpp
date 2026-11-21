// craft/recursive-circuit (#928, the user's example): a basic circuit from raw materials only. It takes 2 vacuum tubes,
// 4 copper wire and a steel plate; a tube takes glass, 2 wire and a steel plate; wire comes in pairs from a copper ingot.
// Carrying 4 copper ingots, 2 glass and 3 steel plates (and knowing the recipes: wire and tubes held once), one request
// queues 4 wire, 2 tubes, 4 more wire, then the circuit; the queue runs them (0.3 s a craft) and the raw materials are used up.
// The crafting screen stays open on the circuit's page so the recording shows the queue.
#include "scenario.h"
#include "gameui.h"
#include "inventory.h"
#include "playersession.h"
#include "terraintool.h"

SCENARIO("craft/recursive-circuit") {
    s.tags = { "craft", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    static int circuit = -1, steps = 0;
    s.step("queue a circuit from raw materials", 10.0, [](ScenarioCtx& c) {
        if (c.T < 0.5) return false;
        terrainTool.creative = false;
        const auto& recipes = itemDB.recipes();
        for (size_t i = 0; i < recipes.size(); ++i) if (recipes[i].output == "basic_circuit") circuit = (int)i;
        if (!c.check("the basic circuit recipe exists", circuit >= 0)) return true;
        // Known recipes: their ingredients have been held once (RecipeBook_Discovered).
        for (const char* id : { "copper_wire", "vacuum_tube" }) { inventory.add(id, 1); inventory.remove(id, 1); }
        inventory.add("copper_ingot", 4); inventory.add("glass", 2); inventory.add("steel_plate", 3);
        GameUI_ShowScreen(GameScreen::Crafting);
        GameUI_SelectRecipe(circuit);
        const bool ok = GameUI_QueueRecursive(circuit, 1, &steps);
        c.check("the circuit can be planned from raw materials", ok);
        c.eq("plan steps (wire for the tubes, tubes, wire for the circuit, circuit)", steps, 4);
        c.eq("crafts queued (4 wire, 2 tubes, 1 circuit)", GameUI_CraftQueueLength(), 7);
        return true;
    });
    s.step("the queue runs: intermediates first, then the circuit", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("the queue to empty");
        static bool shot = false; // (the recording has no UI: one screenshot of the page and the queue mid-run)
        if (!shot && GameUI_CraftQueueLength() <= 5) { TakeScreenshot("bench_craft_recursive.png"); c.note("screenshot: bench_craft_recursive.png"); shot = true; }
        if (GameUI_CraftQueueLength() > 0) return false;
        c.eq("basic circuits", inventory.count("basic_circuit"), 1);
        c.eq("copper ingots left", inventory.count("copper_ingot"), 0);
        c.eq("glass left", inventory.count("glass"), 0);
        c.eq("steel plates left", inventory.count("steel_plate"), 0);
        c.eq("wire left (8 made, 8 used)", inventory.count("copper_wire"), 0);
        c.eq("tubes left", inventory.count("vacuum_tube"), 0);
        GameUI_ToggleInventory();
        return true;
    });
}
