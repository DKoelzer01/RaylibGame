// build/multiblock-tank-parts (#1429 PR 2): the large tank's multiblock parts (bronze casings, a fluid in and a fluid out
// hatch, the controller) and its glass walls placed layer by layer in a plot, the camera on them; then the structure
// validates against multiblocks.txt through PartItem (the parts carry their catalog index in Part::material). (PR 3) The
// wrench on the controller forms it; breaking a casing unforms it.
#include "scenario.h"
#include "buildtool.h"
#include "hotbar.h"
#include "inventory.h"
#include "multiblock.h"
#include "player.h"
#include "playerinput.h"
#include "plots.h"
#include "world.h"
#include <algorithm>
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
Int3 base{};
Vector3 At(float x, float y, float z) { return plot->toBody({ x * kCellSize, y * kCellSize, z * kCellSize }); }
// A part for a legend key of the large tank (glass: the building block whose item is "glass").
bool PartFor(char key, int& hatchesPlaced, Part& out) {
    PartType t; int m = 0;
    switch (key) {
    case 'c': if (!Multiblock_ItemPart("bronze_casing", t, m)) return false; break;
    case 'C': if (!Multiblock_ItemPart("large_tank_controller", t, m)) return false; break;
    case 'h': if (!Multiblock_ItemPart(hatchesPlaced++ == 0 ? "fluid_in_hatch" : "fluid_out_hatch", t, m)) return false; break;
    case 'g': t = PartType::Cube; m = std::max(0, BuildMaterialForItem("glass")); break;
    default: return false;
    }
    out = Part{ t, {}, 0, (uint8_t)m };
    return true;
}
}

SCENARIO("build/multiblock-tank-parts") {
    s.tags = { "world", "build", "multiblock" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 60.0;
    s.step("a plot beside the player", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        base.x += 3; base.y += 1; // (beside the player, clear of the ground)
        return true;
    });
    s.step("the tank placed a layer at a time", 20.0, [](ScenarioCtx& c) {
        static int layer = 0, hatches = 0, placed = 0;
        if (c.t < 0.05) { layer = 0; hatches = 0; placed = 0; }
        const MbDef* d = Multiblock_All().find("large_tank");
        if (!c.check("large_tank defined", d != nullptr)) return true;
        const Vector3 up = Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0)));
        c.camera(At(base.x + 1.5f, base.y + 4.5f, base.z + 7.5f), At(base.x + 1.5f, base.y + 1.2f, base.z + 1.5f), up); // in front (+z: the controller's face), above
        if (c.t < 0.6 * (layer + 1)) return false; // a layer every 0.6 s, so the capture sees it build
        if (layer >= (int)d->layers.size()) {
            c.eq("parts placed (9 casings, hatch and controller below; 8 glass around an air cell; 9 on top)", placed, 9 + 8 + 9);
            return true;
        }
        for (int z = 0; z < (int)d->layers[layer].size(); ++z)
            for (int x = 0; x < (int)d->layers[layer][z].size(); ++x) {
                Part p;
                if (!PartFor(d->layers[layer][z][x], hatches, p)) continue;
                p.origin = { base.x + x, base.y + layer, base.z + z };
                placed += plot->place(p) >= 0;
            }
        ++layer;
        return false;
    });
    s.step("it validates as a formed large tank", 10.0, [](ScenarioCtx& c) {
        const MbCatalog& cat = Multiblock_All();
        const MbDef* d = cat.find("large_tank");
        const MbInt3 ctl{ base.x + d->controller.x, base.y + d->controller.y, base.z + d->controller.z };
        const MbResult r = Multiblock_Validate(*d, cat, [](MbInt3 w) {
            const int i = plot->partAt({ w.x, w.y, w.z });
            return i >= 0 ? std::string(PartItem(plot->parts[i])) : std::string();
        }, ctl, 0, 1);
        c.note(r.problem.empty() ? "formed" : r.problem.c_str());
        c.check("the placed parts read as their catalog ids and the tank validates", r.formed);
        const Vector3 up = Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0)));
        c.camera(At(base.x + 4.5f, base.y + 3.5f, base.z + 6.0f), At(base.x + 1.5f, base.y + 1.2f, base.z + 1.5f), up);
        return c.t > 1.5; // (hold the view for the capture)
    });
    // (#1429 PR 3) The wrench forms it, as a player does: the camera 2 m in front of the controller and looking at it (the
    // machine aim casts from the camera), the wrench in the hand, one RMB through Machine_Update.
    s.step("RMB with the wrench on the controller forms it", 10.0, [](ScenarioCtx& c) {
        const MbDef* d = Multiblock_All().find("large_tank");
        const Int3 ctl{ base.x + d->controller.x, base.y + d->controller.y, base.z + d->controller.z };
        const int idx = plot->partAt(ctl);
        if (!c.check("a part at the controller's cell", idx >= 0)) return true;
        const Vector3 up = Vector3Normalize(Vector3Subtract(At(0, 1, 0), At(0, 0, 0)));
        const Vector3 face = At(ctl.x + 0.5f, ctl.y + 0.5f, ctl.z + 0.5f);
        c.camera(At(ctl.x + 0.5f, ctl.y + 1.5f, ctl.z + 4.5f), face, up);
        static bool pressed = false;
        if (c.t < 0.05) pressed = false;
        if (c.t < 0.8) return false; // (the controller in view first)
        if (!pressed) {
            pressed = true;
            c.check("not formed before the wrench", !std::binary_search(plot->formed.begin(), plot->formed.end(), idx));
            if (inventory.count("wrench") <= 0) inventory.add("wrench", 1);
            Hotbar_Bind(0, "wrench");
            Hotbar_Press(0, false);
            PlayerInput in; in.press(Act::Secondary);
            player.flyMode = false; // (Machine_Update aims only on foot; c.camera flies the player to the eye)
            Machine_Update(c.world, in);
            player.flyMode = true;
            c.check("formed by the wrench", std::binary_search(plot->formed.begin(), plot->formed.end(), idx));
        }
        return c.t > 2.5; // (hold the view)
    });
    s.step("breaking a casing unforms it", 10.0, [](ScenarioCtx& c) {
        const MbDef* d = Multiblock_All().find("large_tank");
        const Int3 ctl{ base.x + d->controller.x, base.y + d->controller.y, base.z + d->controller.z };
        const int idx = plot->partAt(ctl);
        const int casing = plot->partAt({ base.x, base.y, base.z }); // (a corner casing of the bottom layer)
        if (!c.check("a casing at the corner", casing >= 0 && casing != idx)) return true;
        plot->remove(casing);
        const std::vector<int> dropped = Multiblock_Recheck(*plot);
        c.eq("the controller unformed (the one structure dropped)", (int)dropped.size(), 1);
        c.check("formed list empty", plot->formed.empty());
        return true;
    });
}
