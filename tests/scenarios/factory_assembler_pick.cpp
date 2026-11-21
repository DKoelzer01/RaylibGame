// factory/assembler-pick (#1271 L6c): the LV assembler's screen. In the sandbox, an LV generator cabled to an LV assembler;
// its screen open, with "Makes" (the recipe picker) and the Assemble/Disassemble switch. Fed a circuit's parts (which an
// LV pack's recipe also fits), it makes what is picked: LV packs, then, picked back to basic circuits, a circuit. The pick
// is in the machine's view (what the screen and a joined player see). The recording shows the screen.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "gameui.h"
#include "machineact.h"
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
void Feed() {
    Factory::Insert(bodyId, plot->id, asmPart, "vacuum_tube", 2);
    Factory::Insert(bodyId, plot->id, asmPart, "copper_wire", 4);
    Factory::Insert(bodyId, plot->id, asmPart, "steel_plate", 1);
}
int Count(const std::vector<std::pair<std::string, int>>& slots, const std::string& item) { int n = 0; for (const auto& s : slots) if (s.first == item) n += s.second; return n; }
}

SCENARIO("factory/assembler-pick") {
    s.tags = { "factory", "machines", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("an LV generator cabled to an LV assembler, its screen open", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
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
    s.step("the generator fuelled; the screen open", 5.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        if (c.T < 1.0) return false;
        Factory::InsertFuel(bodyId, plot->id, gen, "coal", 16);
        GameUI_OpenContainer(bodyId, plot->id, asmPart, (int)PartType::LvAssembler);
        return c.check("its screen is open", GameUI_ContainerOpen());
    });
    s.step("picked LV packs: a circuit's parts make LV packs", 20.0, [](ScenarioCtx& c) {
        static bool fed = false;
        if (!fed) { fed = true; MachineAct_SetPick(bodyId, plot->id, asmPart, "lv_pack"); Feed(); }
        Factory::Sync(c.world);
        const MachineView v = View();
        c.waitingFor(TextFormat("pick '%s', lv packs %d", v.pick.c_str(), Count(v.outSlots, "lv_pack")));
        if (Count(v.outSlots, "lv_pack") < 2) return false;
        c.check("the view carries the pick", v.pick == "lv_pack");
        return c.eq("no circuit while LV packs are picked", Count(v.outSlots, "basic_circuit"), 0);
    });
    s.step("picked basic circuits: the rest makes a circuit", 25.0, [](ScenarioCtx& c) {
        static bool fed = false;
        if (!fed) { fed = true; MachineAct_SetPick(bodyId, plot->id, asmPart, "basic_circuit"); Feed(); }
        Factory::Sync(c.world);
        const MachineView v = View();
        c.waitingFor(TextFormat("pick '%s', circuits %d", v.pick.c_str(), Count(v.outSlots, "basic_circuit")));
        if (Count(v.outSlots, "basic_circuit") < 1) return false;
        return c.check("the view carries the pick", v.pick == "basic_circuit");
    });
    s.step("picked Any again", 5.0, [](ScenarioCtx& c) {
        static bool sent = false;
        if (!sent) { sent = true; MachineAct_SetPick(bodyId, plot->id, asmPart, ""); }
        Factory::Sync(c.world);
        const MachineView v = View();
        if (!v.pick.empty()) return false;
        return c.check("no pick: what its inputs allow", v.pick.empty());
    });
}
