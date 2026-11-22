// factory/tiers-ui (#1427 L7d): tiered machines on screen. In the sandbox, an LV generator cabled to an LV assembler and,
// beside it, an MV assembler (the same part, its variant step 1): the MV one is tinted toward MV's orange. Its screen shows
// its own name, "MV 128 V" and the circuit it makes overclocked (an LV recipe in an MV machine: twice as fast); then the
// Help page of the MV assembler lists its tiers and its recipes with their times, the MV pack's at its own tier.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "gameui.h"
#include "helpmanual.h"
#include "machineact.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int gen = -1, lvAsm = -1, mvAsm = -1;
bool Place(ScenarioCtx& c, PartType t, Int3 at, int step = 0) { return BuildTool_Place(*c.body, *plot, { t, at, 0, (uint8_t)step }, false); }
}

SCENARIO("factory/tiers-ui") {
    s.tags = { "factory", "machines", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 60.0;
    s.step("a generator, an LV and an MV assembler", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 2.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        bool ok = Place(c, PartType::Generator, { 0, 0, 0 }) && Place(c, PartType::LvAssembler, { 3, 0, 0 }) && Place(c, PartType::LvAssembler, { 6, 0, 0 }, 1);
        for (int x = 0; x <= 7 && ok; ++x) ok = Place(c, PartType::Cable, { x, 0, -1 });
        if (!c.check("placed", ok)) return true;
        gen = plot->partAt({ 0, 0, 0 }); lvAsm = plot->partAt({ 3, 0, 0 }); mvAsm = plot->partAt({ 6, 0, 0 });
        c.check("the MV assembler is the assembler part at step 1", plot->parts[mvAsm].type == PartType::LvAssembler && plot->parts[mvAsm].material == 1);
        c.check("its item is the MV assembler", std::string(PartItem(plot->parts[mvAsm])) == "mv_assembler");
        c.camera(plot->toBody({ 3.5f, 2.5f, -5.0f }), plot->toBody({ 3.5f, 0.5f, 0.5f }), plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("the machines from outside (the MV one tinted)", 3.0, [](ScenarioCtx& c) { Factory::Sync(c.world); return c.t > 1.5; });
    s.step("the MV assembler's screen, making a circuit overclocked", 20.0, [](ScenarioCtx& c) {
        static bool fed = false;
        if (!fed) {
            fed = true;
            Factory::InsertFuel(bodyId, plot->id, gen, "coal", 16);
            MachineAct_SetPick(bodyId, plot->id, mvAsm, "basic_circuit");
            Factory::Insert(bodyId, plot->id, mvAsm, "vacuum_tube", 4);
            Factory::Insert(bodyId, plot->id, mvAsm, "copper_wire", 8);
            Factory::Insert(bodyId, plot->id, mvAsm, "steel_plate", 2);
            GameUI_OpenContainer(bodyId, plot->id, mvAsm, (int)PartType::LvAssembler);
        }
        Factory::Sync(c.world);
        Factory::View(bodyId, plot->id, mvAsm); Factory::WaitTicks(2);
        const MachineView v = Factory::View(bodyId, plot->id, mvAsm);
        c.waitingFor(TextFormat("duration %d, status %s", v.duration, MachineStatusName(v.status)));
        if (v.duration <= 0 || c.t < 2.0) return false;
        return c.eq("an LV recipe in the MV assembler: half its 160 ticks", v.duration, 80);
    });
    s.step("the Help page of the MV assembler", 5.0, [](ScenarioCtx& c) {
        static bool opened = false;
        if (!opened) { opened = true; HelpManual_OpenItem("mv_assembler"); GameUI_ShowTab("Help"); }
        return c.t > 2.0;
    });
}
