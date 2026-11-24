// build/multiblock-ghost (#1429 PR 5): a large tank's controller alone in a plot; E on it shows the structure ghost and the
// parts list (15 bronze casings, 8 glass, a fluid in and a fluid out hatch). The parts go in where the ghost says, a layer at
// a time with the camera on them (the ghost shrinks as they land), the two hatch spots left empty; then I on the controller
// ("Craft missing") crafts the two hatches through the planner: the first from a casing carried, the second from a casing
// it crafts first (each plan sees what the one before leaves). The hatches placed, the tank validates; E hides the ghost.
#include "scenario.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "multiblock.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "survival.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
Int3 base{}, ctl{};
int ctlPart = -1;
Vector3 At(float x, float y, float z) { return plot->toBody({ x * kCellSize, y * kCellSize, z * kCellSize }); }
Vector3 Up() { return Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0))); }
// The part for an item the plan asks for (glass: the building block whose item is "glass").
bool PartFor(const std::string& item, Part& out) {
    PartType t; int m = 0;
    if (!Multiblock_ItemPart(item, t, m)) {
        if (item != "glass") return false;
        t = PartType::Cube; m = std::max(0, BuildMaterialForItem("glass"));
    }
    out = Part{ t, {}, 0, (uint8_t)m };
    return true;
}
MbPlan PlanNow(MbResult* res = nullptr) {
    const MbDef* d = nullptr;
    const MbResult r = Multiblock_CheckController(*plot, ctlPart, &d);
    if (res) *res = r;
    return d ? Multiblock_Plan(*d, Multiblock_All(), r) : MbPlan{};
}
// One key press through Machine_Update, as a player does: the camera in front of the controller looking at it (the machine
// aim casts from the camera), on foot for the aim.
void PressOnController(ScenarioCtx& c, Act a) {
    PlayerInput in; in.press(a);
    player.flyMode = false;
    Machine_Update(c.world, in);
    player.flyMode = true;
}
// (The controller's HUD can't be recorded here: the scenario camera flies, and the machine aim works only on foot.)
void LookAtController(ScenarioCtx& c) {
    c.camera(At(ctl.x + 0.5f, ctl.y + 1.5f, ctl.z + 4.5f), At(ctl.x + 0.5f, ctl.y + 0.5f, ctl.z + 0.5f), Up());
}
}

SCENARIO("build/multiblock-ghost") {
    s.tags = { "world", "build", "multiblock", "craft" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a large tank controller alone in a plot", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const MbDef* d = Multiblock_All().find("large_tank");
        if (!c.check("large_tank defined", d != nullptr)) return true;
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.x += 3; base.y += 1; // (beside the player, clear of the ground)
        ctl = { base.x + d->controller.x, base.y + d->controller.y, base.z + d->controller.z };
        Part p;
        if (!c.check("the controller is a part", PartFor("large_tank_controller", p))) return true;
        p.origin = ctl;
        ctlPart = plot->place(p);
        return c.check("controller placed", ctlPart >= 0);
    });
    s.step("E on the controller shows the ghost; the parts list", 10.0, [](ScenarioCtx& c) {
        static bool pressed = false;
        if (c.t < 0.05) pressed = false;
        LookAtController(c);
        if (c.t < 0.8) return false; // (the controller in view first)
        if (!pressed) {
            pressed = true;
            c.check("no ghost before E", !BuildTool_GhostShown(c.body->id, plot->id, ctl));
            PressOnController(c, Act::Use);
            c.check("E on the controller shows its ghost", BuildTool_GhostShown(c.body->id, plot->id, ctl));
            c.check("and does not open the game window", !GameUI_InventoryOpen());
            const MbPlan p = PlanNow();
            auto need = [&](const char* id) { const auto it = p.need.find(id); return it == p.need.end() ? 0 : it->second; };
            c.eq("bronze casings needed", need("bronze_casing"), 15);
            c.eq("glass needed", need("glass"), 8);
            c.eq("fluid input hatches needed", need("fluid_in_hatch"), 1);
            c.eq("fluid output hatches needed", need("fluid_out_hatch"), 1);
            c.eq("kinds of part needed", (int)p.need.size(), 4);
        }
        // (the camera back and above, the whole ghost in view, for the capture)
        c.camera(At(base.x + 4.5f, base.y + 4.0f, base.z + 7.0f), At(base.x + 1.5f, base.y + 1.2f, base.z + 1.5f), Up());
        return c.t > 2.5;
    });
    s.step("the parts placed where the ghost says, a layer at a time, the hatch spots left", 20.0, [](ScenarioCtx& c) {
        static int layer = 0;
        if (c.t < 0.05) layer = 0;
        c.camera(At(base.x + 4.5f, base.y + 4.0f, base.z + 7.0f), At(base.x + 1.5f, base.y + 1.2f, base.z + 1.5f), Up());
        if (c.t < 1.0 * (layer + 1)) return false; // (a layer a second, so the capture sees the ghost shrink)
        MbResult r;
        const MbPlan p = PlanNow(&r);
        int placed = 0;
        for (size_t i = 0; i < r.cells.size(); ++i) {
            const std::string& item = p.part[i];
            if (item.empty() || r.cells[i].layer != layer || item.find("_hatch") != std::string::npos) continue;
            Part q;
            if (!PartFor(item, q)) { c.check(("a part for " + item).c_str(), false); continue; }
            q.origin = { r.cells[i].world.x, r.cells[i].world.y, r.cells[i].world.z };
            placed += plot->place(q) >= 0;
        }
        c.note(TextFormat("layer %d: %d parts placed", layer, placed));
        if (++layer < 3) return false;
        const MbPlan left = PlanNow(&r);
        c.eq("left to place: the two hatches", (int)left.need.size(), 2);
        int open = 0;
        for (const MbCellResult& cell : r.cells) open += cell.status != MbCell::Ok;
        c.eq("cells still to build (the two hatch spots)", open, 2);
        return true;
    });
    s.step("I on the controller crafts the missing hatches", 25.0, [](ScenarioCtx& c) {
        static bool pressed = false;
        if (c.t < 0.05) pressed = false;
        LookAtController(c);
        if (c.t < 0.2) return false; // (past the reset frames above: one press)
        if (!pressed) {
            pressed = true;
            terrainTool.creative = false; // (creative crafting costs nothing: the plan would not be tested)
            QuestBook_SetTier(3);         // (the casings and hatches unlock at tier 3)
            for (const char* id : { "bronze_casing", "bronze_plate", "bronze_screw", "bronze_rod", "valve", "fluid_in_hatch", "fluid_out_hatch" })
                inventory.remove(id, inventory.count(id));
            // One casing carried, a casing's worth of parts, two valves and the hatches' plates (2 each + 4 for the casing).
            inventory.add("bronze_casing", 1); inventory.add("bronze_plate", 8); inventory.add("bronze_screw", 4);
            inventory.add("bronze_rod", 2); inventory.add("valve", 2);
            PressOnController(c, Act::Insert);
            c.eq("crafts queued (a hatch; a casing, then a hatch)", GameUI_CraftQueueLength(), 3);
        }
        c.waitingFor("the crafting queue to empty");
        if (GameUI_CraftQueueLength() > 0) return false;
        c.eq("fluid input hatches", inventory.count("fluid_in_hatch"), 1);
        c.eq("fluid output hatches", inventory.count("fluid_out_hatch"), 1);
        c.eq("casings left (one carried, one made, both used)", inventory.count("bronze_casing"), 0);
        c.eq("plates left", inventory.count("bronze_plate"), 0);
        return true;
    });
    s.step("the hatches placed, the tank validates; E hides the ghost", 10.0, [](ScenarioCtx& c) {
        static bool done = false;
        if (c.t < 0.05) done = false;
        c.camera(At(base.x + 4.5f, base.y + 4.0f, base.z + 7.0f), At(base.x + 1.5f, base.y + 1.2f, base.z + 1.5f), Up());
        if (c.t < 0.2) return false; // (past the reset frames above: once)
        if (!done) {
            done = true;
            MbResult r;
            const MbPlan p = PlanNow(&r);
            for (size_t i = 0; i < r.cells.size(); ++i) {
                Part q;
                if (p.part[i].empty() || !PartFor(p.part[i], q)) continue;
                if (inventory.count(p.part[i]) > 0) inventory.remove(p.part[i], 1);
                q.origin = { r.cells[i].world.x, r.cells[i].world.y, r.cells[i].world.z };
                plot->place(q);
            }
            PlanNow(&r);
            c.note(r.problem.empty() ? "validates" : r.problem.c_str());
            c.check("the tank validates", r.formed);
            c.check("nothing left on the parts list", PlanNow().need.empty());
        }
        if (c.t < 1.5) return false; // (the finished tank, no ghost left, in view)
        LookAtController(c);
        PressOnController(c, Act::Use);
        c.check("E again hides the ghost", !BuildTool_GhostShown(c.body->id, plot->id, ctl));
        return true;
    });
}
