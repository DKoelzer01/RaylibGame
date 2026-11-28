// factory/large-macerator (#1429 PR 8): the first electric structure, end to end in the game. A large macerator
// (multiblocks.txt: a steel crusher, ore in at the hopper on top, crushed ore out at the back of the base, `machine
// macerator`, `parallel 4`) built whole beside the player, an LV energy hatch on its side cabled to a generator. The wrench
// forms it; 8 iron ore into its input hatch: with the generator unfuelled nothing is crushed (it runs only on its energy
// hatch); fuelled, it crushes four loads at once: the first crushed iron comes out as 12 together (a macerator makes 3 from
// 2 ore), all 8 ore taken; the controller reads working meanwhile. The camera holds the structure for the capture.
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
int ctlPart = -1, inHatch = -1, outHatch = -1, energyHatch = -1, gen = -1;
bool sawWorking = false;
Vector3 At(float x, float y, float z) { return plot->toBody({ x * kCellSize, y * kCellSize, z * kCellSize }); }
Vector3 Up() { return Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0))); }
// High in front, between the structure and the generator: the controller, the cable and the generator, the hopper on top.
void View(ScenarioCtx& c) { c.camera(At(ctl.x - 1.5f, ctl.y + 7.5f, ctl.z + 7.0f), At(ctl.x - 1.5f, ctl.y + 1.5f, ctl.z - 1.0f), Up()); }
int Count(const std::vector<std::pair<std::string, int>>& slots, const char* item) { int n = 0; for (const auto& [i, k] : slots) if (i == item) n += k; return n; }
MachineView ViewOf(ScenarioCtx& c, int part) { // (View follows one machine: the first call after a switch is stale)
    Factory::View(c.body->id, plot->id, part); Factory::WaitTicks(2);
    return Factory::View(c.body->id, plot->id, part);
}
}

SCENARIO("factory/large-macerator") {
    s.tags = { "factory", "multiblock", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a whole large macerator beside the player, cabled to a generator", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const MbDef* d = Multiblock_All().find("large_macerator");
        if (!c.check("large_macerator defined", d != nullptr)) return true;
        c.check("its machine: the macerator", d->machine == "macerator", d->machine.c_str());
        c.eq("parallel", d->parallel, 4);
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.x += 5; base.y += 1; // (beside the player, the generator between, clear of the ground)
        ctl = { base.x + d->controller.x, base.y + d->controller.y, base.z + d->controller.z };
        PartType t; int m = 0;
        if (!c.check("the controller is a part", Multiblock_ItemPart("large_macerator_controller", t, m))) return true;
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
            if (kind == "item_in") inHatch = i; else if (kind == "item_out") outHatch = i; else if (kind == "energy") energyHatch = i;
        }
        c.check("an input, an output and an energy hatch", inHatch >= 0 && outHatch >= 0 && energyHatch >= 0);
        if (energyHatch < 0) return true;
        // The generator (2x2x2) two cells west of the energy hatch, a cable between them.
        const Int3 eh = plot->parts[energyHatch].origin;
        gen = plot->place(Part{ PartType::Generator, { eh.x - 3, eh.y, eh.z }, 0 });
        const int cable = plot->place(Part{ PartType::Cable, { eh.x - 1, eh.y, eh.z }, 0 });
        c.check("generator and cable placed", gen >= 0 && cable >= 0);
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
        c.check("idle before ore", !Factory::Working(c.body->id, plot->id, ctlPart));
        Factory::Insert(c.body->id, plot->id, inHatch, "iron_ore", 8); // (four loads' worth: all start at once)
        return true;
    });
    s.step("no fuel in the generator: nothing crushed", 10.0, [](ScenarioCtx& c) {
        View(c);
        Factory::WaitTicks(400); // (over twice a load's 150 ticks)
        c.eq("no crushed iron without power", Count(ViewOf(c, outHatch).outSlots, "iron_crushed"), 0);
        Factory::InsertFuel(c.body->id, plot->id, gen, "coal", 16);
        Factory::SetSpeed(10);
        sawWorking = false;
        return true;
    });
    s.step("powered: four loads at once, the first crushed iron 12 together", 30.0, [](ScenarioCtx& c) {
        View(c);
        sawWorking |= Factory::Working(c.body->id, plot->id, ctlPart);
        const MachineView out = ViewOf(c, outHatch);
        const int crushed = Count(out.outSlots, "iron_crushed");
        c.waitingFor(TextFormat("crushed iron %d", crushed));
        if (crushed == 0) return false;
        Factory::SetSpeed(1);
        c.eq("the first crushed iron: four loads together", crushed, 12);
        c.eq("all 8 ore taken (4 x 2)", Count(ViewOf(c, inHatch).inSlots, "iron_ore"), 0);
        c.check("the controller read working meanwhile", sawWorking);
        return c.t >= 1.0; // (a moment on the full output for the recording)
    });
}
