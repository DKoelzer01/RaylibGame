// factory/anvil-repair (#1219, user 2026-11-22): an anvil works on a composed copy, through the machine paths a player uses
// (insert, strikes, take everything out), in the sandbox:
//   repair:  a worn bronze sword with a ruby and a bronze ingot; two strikes: back in the inventory with its wear down by
//            half its durability, the ruby kept, the ingot used;
//   a gem out: the same sword alone; two strikes: back without the ruby, and the ruby with it unless it broke (25%, GUESS).
#include "scenario.h"
#include "buildtool.h"
#include "composed.h"
#include "factory.h"
#include "inventory.h"
#include "itemdata.h"
#include "machineact.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int anvil = -1;
std::string worn, repaired;
int stage = 0;
double at = 0.0;
int Wear(const std::string& d) { return ItemData::GetInt(d, "w", 0); }
int SwordCell() { for (int i = 0; i < Inventory::kCells; ++i) if (inventory.grid()[i].id == "weapon_sword" && inventory.grid()[i].hasData()) return i; return -1; }
}

SCENARIO("factory/anvil-repair") {
    s.tags = { "factory" }; // (one runner step: a scenario matching two runs twice)
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 60.0;
    s.step("an anvil", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 p = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 2.0f)));
        const Vector3 ground = Vector3Scale(p, c.field->radius + c.field->elevation(p) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        if (!c.check("placed", BuildTool_Place(*c.body, *plot, { PartType::Anvil, { 0, 0, 0 }, 0, 0 }, false))) return true;
        anvil = plot->partAt({ 0, 0, 0 });
        c.camera(plot->toBody({ 0.5f, 1.8f, -2.2f }), plot->toBody({ 0.5f, 0.4f, 0.5f }), plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("repair: a worn bronze sword with a ruby, and a bronze ingot; two strikes", 20.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        if (stage == 0) {
            if (c.t < 0.5) return false;
            worn = ItemData::Set(ItemData::SetInt(Composed_MakeData("bronze"), "w", 300), "j", "ruby");
            MachineAct_Insert(bodyId, plot->id, anvil, "weapon_sword", 1, false, 0, worn);
            MachineAct_Insert(bodyId, plot->id, anvil, "bronze_ingot", 1);
            stage = 1; at = c.t; return false;
        }
        if (stage == 1 && c.t > at + 0.5) { Factory::Strike(bodyId, plot->id, anvil); stage = 2; at = c.t; return false; }
        if (stage == 2 && c.t > at + 1.0) { Factory::Strike(bodyId, plot->id, anvil); stage = 3; at = c.t; return false; }
        if (stage == 3 && c.t > at + 1.5) { MachineAct_Extract(bodyId, plot->id, anvil); stage = 4; at = c.t; return false; }
        if (stage == 4 && c.t > at + 0.5) {
            const ComposedStats w0 = Composed_Stats("weapon_sword", worn);
            const int cell = SwordCell();
            if (!c.check("the sword back in the inventory", cell >= 0)) { c.fail("no sword"); return true; }
            repaired = inventory.grid()[cell].data;
            c.eq("its wear down by half its durability", Wear(repaired), std::max(0, w0.wear - w0.durability / 2));
            c.check("the ruby kept", ItemData::Get(repaired, "j") == "ruby", repaired);
            c.eq("the ingot used (none came back)", (int)inventory.count("bronze_ingot"), 0);
            c.keepCapture("an anvil repairs a composed sword (#1219)");
            stage = 0;
            return true;
        }
        return false;
    });
    s.step("a gem out: the sword alone; two strikes", 20.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        if (stage == 0) {
            if (c.t < 0.3) return false;
            const int cell = SwordCell();
            if (!c.check("the sword", cell >= 0)) { c.fail("no sword"); return true; }
            inventory.takeInstance(cell);
            MachineAct_Insert(bodyId, plot->id, anvil, "weapon_sword", 1, false, 0, repaired);
            stage = 1; at = c.t; return false;
        }
        if (stage == 1 && c.t > at + 0.5) { Factory::Strike(bodyId, plot->id, anvil); stage = 2; at = c.t; return false; }
        if (stage == 2 && c.t > at + 1.0) { Factory::Strike(bodyId, plot->id, anvil); stage = 3; at = c.t; return false; }
        if (stage == 3 && c.t > at + 1.5) { MachineAct_Extract(bodyId, plot->id, anvil); stage = 4; at = c.t; return false; }
        if (stage == 4 && c.t > at + 0.5) {
            const int cell = SwordCell();
            if (!c.check("the sword back", cell >= 0)) { c.fail("no sword"); return true; }
            const std::string d = inventory.grid()[cell].data;
            c.check("no gem in it", ItemData::Get(d, "j").empty(), d);
            c.eq("its wear unchanged", Wear(d), Wear(repaired));
            c.note("the ruby came back: " + std::to_string((int)inventory.count("ruby")) + " (0: it broke, 25%)");
            c.check("at most the one ruby", inventory.count("ruby") <= 1);
            return true;
        }
        return false;
    });
}
