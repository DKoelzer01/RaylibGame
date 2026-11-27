// build/multiblock-lamps (#1429 PR 6): a large tank built whole beside the player shows a dim red lamp on its controller
// (unformed); the wrench forms it: the lamp turns green, the two hatches light in their fluids' blues and a seam glow runs
// along the tank's edges; working (as a joined client hears it from its server: the tank runs no recipes yet) the lamps
// pulse; a casing broken unforms it and the lamp goes red again. The camera holds each state for the capture.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "headless.h"
#include "hotbar.h"
#include "inventory.h"
#include "machinemodels.h"
#include "multiblock.h"
#include "multiblockmodels.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
Int3 base{}, ctl{};
int ctlPart = -1;
Vector3 At(float x, float y, float z) { return plot->toBody({ x * kCellSize, y * kCellSize, z * kCellSize }); }
Vector3 Up() { return Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0))); }
bool PartFor(const std::string& item, Part& out) { // (glass: the building block whose item is "glass")
    PartType t; int m = 0;
    if (!Multiblock_ItemPart(item, t, m)) {
        if (item != "glass") return false;
        t = PartType::Cube; m = std::max(0, BuildMaterialForItem("glass"));
    }
    out = Part{ t, {}, 0, (uint8_t)m };
    return true;
}
MbLook LookNow() {
    std::vector<MbLook> looks;
    Multiblock_Looks(*plot, looks);
    return looks.size() == 1 ? looks[0] : MbLook{};
}
// Front and a little to the east, above: the controller's face, a hatch's lamp on the top and two edges' seams in view.
void View(ScenarioCtx& c) { c.camera(At(base.x + 4.0f, base.y + 4.2f, base.z + 6.5f), At(base.x + 1.5f, base.y + 1.2f, base.z + 2.0f), Up()); }
void CheckDrawn(ScenarioCtx& c) {
    if (!Headless()) c.check("the controller's lamp drawn", MultiblockModels_DrawnCount() >= 1);
}
}

SCENARIO("build/multiblock-lamps") {
    s.tags = { "world", "build", "multiblock", "render" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 60.0;
    s.step("a whole large tank beside the player, not formed: a red lamp", 20.0, [](ScenarioCtx& c) {
        static bool built = false;
        static double builtAt = 0.0;
        if (c.t < 0.05) built = false;
        if (!built) {
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
            if (!c.check("controller placed", ctlPart >= 0)) return true;
            MbResult r = Multiblock_CheckController(*plot, ctlPart);
            const MbPlan plan = Multiblock_Plan(*d, Multiblock_All(), r);
            for (size_t i = 0; i < r.cells.size(); ++i) {
                Part q;
                if (plan.part[i].empty() || !PartFor(plan.part[i], q)) continue;
                q.origin = { r.cells[i].world.x, r.cells[i].world.y, r.cells[i].world.z };
                plot->place(q);
            }
            c.check("it validates", Multiblock_CheckController(*plot, ctlPart).formed);
            built = true;
            builtAt = c.t;
        }
        View(c);
        if (c.t - builtAt < 2.0) return false;
        const MbLook l = LookNow();
        c.check("one controller, unformed", l.controller == ctlPart && !l.formed);
        c.eq("lamps: the controller's only", (int)l.lamps.size(), 1);
        CheckDrawn(c);
        return true;
    });
    s.step("the wrench forms it: a green lamp, the hatches' lamps and the seams", 10.0, [](ScenarioCtx& c) {
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
        if (c.t < 2.8) return false;
        const MbLook l = LookNow();
        c.check("formed", l.formed);
        c.eq("lamps: the controller and two hatches", (int)l.lamps.size(), 3);
        c.eq("the box: 3 cells wide", l.hi.x - l.lo.x + 1, 3);
        c.check("idle: no recipes (a tank)", !Factory::Working(c.body->id, plot->id, ctlPart));
        CheckDrawn(c);
        return true;
    });
    // The tank runs no recipes yet (its storage is the first-multiblock PR), so the working look comes the way a joined
    // client gets it: its server's working list, which the lamps read through Factory::Working.
    s.step("working: the lamps pulse", 10.0, [](ScenarioCtx& c) {
        if (c.t < 0.05) {
            Factory::SetRemote(true);
            Factory::SetRemoteWorking(c.body->id, plot->id, { ctlPart });
        }
        View(c);
        if (c.t < 3.0) return false;
        c.check("working (the server's view)", Factory::Working(c.body->id, plot->id, ctlPart));
        CheckDrawn(c);
        Factory::SetRemote(false);
        return true;
    });
    s.step("a casing broken: unformed, the lamp red again", 10.0, [](ScenarioCtx& c) {
        static bool broken = false;
        if (c.t < 0.05) broken = false;
        View(c);
        if (!broken) {
            broken = true;
            const int corner = plot->partAt(base); // (the tank's bottom north-west corner: a casing)
            if (!c.check("a casing at the corner", corner >= 0 && plot->parts[corner].type == PartType::MbCasing)) return true;
            plot->remove(corner);
            c.eq("the recheck drops the tank", (int)Multiblock_Recheck(*plot).size(), 1);
        }
        if (c.t < 2.0) return false;
        const MbLook l = LookNow();
        c.check("unformed", !l.formed);
        c.eq("lamps: the controller's only", (int)l.lamps.size(), 1);
        return true;
    });
}
