// factory/disassembly (#1427 L7c): an LV assembler set to Disassemble, with the game's own data. In the sandbox, an LV
// generator cabled to an LV assembler; its screen open, switched to Disassemble. Fed two basic circuits (LV: 50% back) and an
// iron pickaxe (a composed tool, iron parts two tiers under LV: 100%), it gives back the circuits' direct inputs (the raw
// steel plate as a fraction: two halves make one) and the pickaxe's raw materials (3 iron ingots, 2 planks, 2 fibre).
// A raw item (a steel plate) is refused. The recording shows the screen with each input's yield.
#include "scenario.h"
#include "buildtool.h"
#include "composed.h"
#include "disassembly.h"
#include "factory.h"
#include "gameui.h"
#include "machineact.h"
#include "machinetiers.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int gen = -1, asmPart = -1;
bool Place(ScenarioCtx& c, PartType t, Int3 at) { return BuildTool_Place(*c.body, *plot, { t, at, 0 }, false); }
MachineView View() { Factory::View(bodyId, plot->id, asmPart); Factory::WaitTicks(2); return Factory::View(bodyId, plot->id, asmPart); }
int Count(const std::vector<std::pair<std::string, int>>& slots, const std::string& item) { int n = 0; for (const auto& s : slots) if (s.first == item) n += s.second; return n; }
}

SCENARIO("factory/disassembly") {
    s.tags = { "factory", "machines", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("an LV generator cabled to an LV assembler, its screen open", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 2.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        bool ok = Place(c, PartType::Generator, { 0, 0, 0 }) && Place(c, PartType::LvAssembler, { 3, 0, 0 });
        for (int x = 0; x <= 4 && ok; ++x) ok = Place(c, PartType::Cable, { x, 0, -1 });
        if (!c.check("generator, assembler and cable placed", ok)) return true;
        gen = plot->partAt({ 0, 0, 0 }); asmPart = plot->partAt({ 3, 0, 0 });
        c.camera(plot->toBody({ 2.0f, 2.5f, -4.5f }), plot->toBody({ 2.0f, 0.5f, 0.5f }), plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("the data: circuits are LV and taken apart, a steel plate is raw", 5.0, [](ScenarioCtx& c) {
        const DisassemblyEntry* e = Disassembly_Find("basic_circuit");
        if (!c.check("basic circuits are taken apart", e != nullptr)) return true;
        c.eq("a basic circuit is an LV item", e->tier, 0);
        c.eq("an LV assembler gives 50% of an LV item", MachineTiers_YieldPercent(0, e->tier), 50);
        c.check("a steel plate is raw: never taken", Disassembly_Find("steel_plate") == nullptr);
        c.check("research packs never are", Disassembly_Find("lv_pack") == nullptr);
        return true;
    });
    s.step("fuelled, its screen open, set to Disassemble", 5.0, [](ScenarioCtx& c) {
        static bool sent = false;
        Factory::Sync(c.world);
        if (c.T < 1.0) return false;
        if (!sent) {
            sent = true;
            Factory::InsertFuel(bodyId, plot->id, gen, "coal", 16);
            GameUI_OpenContainer(bodyId, plot->id, asmPart, (int)PartType::LvAssembler);
            MachineAct_SetMode(bodyId, plot->id, asmPart, 1);
        }
        const MachineView v = View();
        if (v.mode != 1) return false;
        c.check("a steel plate is refused", !Factory::Accepts(PartType::LvAssembler, "steel_plate", 1));
        return c.check("its screen is open", GameUI_ContainerOpen());
    });
    s.step("two circuits come back as their parts at 50%", 30.0, [](ScenarioCtx& c) {
        static bool fed = false;
        if (!fed) { fed = true; Factory::Insert(bodyId, plot->id, asmPart, "basic_circuit", 2); }
        Factory::Sync(c.world);
        const MachineView v = View();
        int left = 0; for (const auto& s : v.inSlots) left += s.second;
        const int steel = Count(v.outSlots, "steel_plate");
        c.waitingFor(TextFormat("inputs left %d, steel plates %d, %d/%d", left, steel, v.progress, v.duration));
        if (left > 0 || c.t < 1.5) return false; // (held a moment for the recording)
        c.eq("two circuits' steel plate at 50%: one (two halves)", steel, 1);
        const int wires = Count(v.outSlots, "copper_wire"), tubes = Count(v.outSlots, "vacuum_tube");
        c.note(TextFormat("rolled at 50%%: %d of 8 copper wires, %d of 4 vacuum tubes", wires, tubes));
        c.check("never more wires or tubes than went in", wires <= 8 && tubes <= 4);
        Factory::Extract(bodyId, plot->id, asmPart); // (its outputs to the inventory: room for the next)
        return true;
    });
    s.step("an iron pickaxe comes back as its parts' raw materials", 30.0, [](ScenarioCtx& c) {
        static bool fed = false;
        Factory::Sync(c.world);
        if (!fed) {
            if (!View().outSlots.empty() && Count(View().outSlots, "steel_plate") > 0) return false; // (the extract first)
            fed = true;
            Factory::Insert(bodyId, plot->id, asmPart, "tool_pickaxe", 1, 0, 0, Composed_MakeData("iron"));
        }
        const MachineView v = View();
        int left = 0; for (const auto& s : v.inSlots) left += s.second;
        const int ingots = Count(v.outSlots, "iron_ingot");
        c.waitingFor(TextFormat("inputs left %d, iron ingots %d, %d/%d", left, ingots, v.progress, v.duration));
        if (left > 0 || ingots == 0 || c.t < 1.5) return false;
        c.eq("its iron head: 3 ingots (iron is two tiers under LV: 100%)", ingots, 3);
        c.eq("its wood handle: 2 planks", Count(v.outSlots, "oak_planks"), 2);
        return c.eq("its fibre binding: 2 fibre", Count(v.outSlots, "fibre"), 2);
    });
}
