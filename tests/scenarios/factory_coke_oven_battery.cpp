// factory/coke-oven-battery (#1429 PR 7): the first multiblock that runs recipes, end to end in the game. A coke oven battery
// (multiblocks.txt: a 3x3x3 block of bronze casings, an item input and an item output hatch, `machine coke_oven`,
// `parallel 4`) built whole beside the player; the wrench forms it; 8 coal into its input hatch and it bakes four loads at
// once: the first coke comes out as 4 together (a single coke oven makes 1 from 2 coal), all 8 coal taken; the controller
// reads working meanwhile. The camera holds the structure for the capture.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "hotbar.h"
#include "inventory.h"
#include "multiblock.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
Int3 base{}, ctl{};
int ctlPart = -1, inHatch = -1, outHatch = -1;
bool sawWorking = false;
Vector3 At(float x, float y, float z) { return plot->toBody({ x * kCellSize, y * kCellSize, z * kCellSize }); }
Vector3 Up() { return Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0))); }
// Front, to the east and above: the controller, both hatches and the top of the block in view.
void View(ScenarioCtx& c) { c.camera(At(ctl.x + 4.0f, ctl.y + 3.5f, ctl.z + 5.0f), At(ctl.x + 0.5f, ctl.y + 1.0f, ctl.z - 1.0f), Up()); }
int Count(const std::vector<std::pair<std::string, int>>& slots, const char* item) { int n = 0; for (const auto& [i, k] : slots) if (i == item) n += k; return n; }
MachineView ViewOf(ScenarioCtx& c, int part) { // (View follows one machine: the first call after a switch is stale)
    Factory::View(c.body->id, plot->id, part); Factory::WaitTicks(2);
    return Factory::View(c.body->id, plot->id, part);
}
}

SCENARIO("factory/coke-oven-battery") {
    s.tags = { "factory", "multiblock", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a whole coke oven battery beside the player", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const MbDef* d = Multiblock_All().find("coke_oven_battery");
        if (!c.check("coke_oven_battery defined", d != nullptr)) return true;
        c.check("its machine: the coke oven", d->machine == "coke_oven", d->machine.c_str());
        c.eq("parallel", d->parallel, 4);
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.x += 3; base.y += 1; // (beside the player, clear of the ground)
        ctl = { base.x + d->controller.x, base.y + d->controller.y, base.z + d->controller.z };
        PartType t; int m = 0;
        if (!c.check("the controller is a part", Multiblock_ItemPart("coke_oven_battery_controller", t, m))) return true;
        ctlPart = plot->place(Part{ t, ctl, 0, (uint8_t)m });
        if (!c.check("controller placed", ctlPart >= 0)) return true;
        const MbResult r = Multiblock_CheckController(*plot, ctlPart);
        const MbPlan plan = Multiblock_Plan(*d, Multiblock_All(), r);
        for (size_t i = 0; i < r.cells.size(); ++i) {
            if (plan.part[i].empty() || !Multiblock_ItemPart(plan.part[i], t, m)) continue;
            plot->place(Part{ t, { r.cells[i].world.x, r.cells[i].world.y, r.cells[i].world.z }, 0, (uint8_t)m });
        }
        const MbResult done = Multiblock_CheckController(*plot, ctlPart);
        c.check("it validates", done.formed, done.problem.c_str());
        for (int i = 0; i < (int)plot->parts.size(); ++i) {
            if (plot->parts[i].type != PartType::MbHatch) continue;
            const std::string kind = Multiblock_HatchKind(plot->parts[i].material);
            if (kind == "item_in") inHatch = i; else if (kind == "item_out") outHatch = i;
        }
        c.check("an input and an output hatch", inHatch >= 0 && outHatch >= 0);
        return true;
    });
    s.step("the wrench forms it", 10.0, [](ScenarioCtx& c) {
        static bool pressed = false;
        if (c.t < 0.05) pressed = false;
        c.camera(At(ctl.x + 0.5f, ctl.y + 1.5f, ctl.z + 4.5f), At(ctl.x + 0.5f, ctl.y + 0.5f, ctl.z + 0.5f), Up());
        if (c.t < 0.5) return false;
        if (!pressed) {
            pressed = true;
            if (inventory.count("wrench") <= 0) inventory.add("wrench", 1);
            Hotbar_Bind(0, "wrench");
            Hotbar_Press(0, false);
            PlayerInput in; in.press(Act::Secondary);
            player.flyMode = false; // (Machine_Update aims only on foot; c.camera flies the player to the eye)
            Machine_Update(c.world, in);
            player.flyMode = true;
            c.check("formed by the wrench", std::binary_search(plot->formed.begin(), plot->formed.end(), ctlPart));
        }
        View(c);
        if (c.t < 1.5) return false;
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        c.check("idle before coal", !Factory::Working(c.body->id, plot->id, ctlPart));
        Factory::Insert(c.body->id, plot->id, inHatch, "coal", 8); // (four loads' worth: all start at once)
        Factory::SetSpeed(10);
        sawWorking = false;
        return true;
    });
    s.step("four loads at once: the first coke comes out 4 together", 30.0, [](ScenarioCtx& c) {
        View(c);
        sawWorking |= Factory::Working(c.body->id, plot->id, ctlPart);
        const MachineView out = ViewOf(c, outHatch);
        const int coke = Count(out.outSlots, "coke");
        c.waitingFor(TextFormat("coke %d", coke));
        if (coke == 0) return false;
        Factory::SetSpeed(1);
        c.eq("the first coke: four loads together", coke, 4);
        c.eq("all 8 coal taken (4 x 2)", Count(ViewOf(c, inHatch).inSlots, "coal"), 0);
        c.check("the controller read working meanwhile", sawWorking);
        return c.t >= 1.0; // (a moment on the full output for the recording)
    });
}
