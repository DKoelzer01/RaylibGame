// factory/anvil-workpiece (#1613, the user's rules): RMB on an anvil without a hammer puts the held metal on it (a
// workpiece, drawn on its face) and, empty-handed, takes it off - the finished piece first. Through BuildTool_AnvilRmb,
// the call the RMB press makes, in the sandbox:
//   refused: an oak log (the anvil doesn't work it) stays in the inventory;
//   put on:  an iron bloom leaves the inventory and shows as the anvil's workpiece;
//   off:     empty-handed, it comes back;
//   forged:  on again, struck until it is an iron ingot: the finished piece is the workpiece, and RMB takes the ingot.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "inventory.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int anvil = -1;
int stage = 0;
double at = 0.0;
int strikes = 0;
}

SCENARIO("factory/anvil-workpiece") {
    s.tags = { "factory" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("an anvil, and a bloom and a log in the inventory", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = false;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 p = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 2.0f)));
        const Vector3 ground = Vector3Scale(p, c.field->radius + c.field->elevation(p) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        if (!c.check("placed", BuildTool_Place(*c.body, *plot, { PartType::Anvil, { 0, 0, 0 }, 0, 0 }, false))) return true;
        anvil = plot->partAt({ 0, 0, 0 });
        inventory.add("iron_bloom", 1);
        inventory.add("oak_log", 1);
        c.camera(plot->toBody({ 0.5f, 1.4f, -1.1f }), plot->toBody({ 0.5f, 0.6f, 0.5f }), plot->up); // close over its face
        Factory::Sync(c.world);
        Factory::View(bodyId, plot->id, anvil); // watched from now on, as an aimed anvil is (its panel)
        return true;
    });
    s.step("refused, then put on", 10.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        Factory::View(bodyId, plot->id, anvil);
        if (stage == 0) {
            if (c.t < 0.5) return false;
            c.check("an oak log is refused", BuildTool_AnvilRmb(bodyId, plot->id, anvil, "oak_log") == AnvilRmb::Refused);
            c.eq("and stays in the inventory", (int)inventory.count("oak_log"), 1);
            c.check("the bloom is put on", BuildTool_AnvilRmb(bodyId, plot->id, anvil, "iron_bloom") == AnvilRmb::Put);
            c.eq("and leaves the inventory", (int)inventory.count("iron_bloom"), 0);
            stage = 1; at = c.t; return false;
        }
        if (c.t < at + 0.5) return false;
        const std::string w = Factory::Workpiece(bodyId, plot->id, anvil);
        c.check("the anvil's workpiece is the bloom", w == "iron_bloom", w);
        c.keepCapture("a bloom on the anvil (#1613)");
        stage = 0;
        return true;
    });
    s.step("off again, empty-handed", 10.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        Factory::View(bodyId, plot->id, anvil);
        if (stage == 0) {
            if (c.t < 0.5) return false;
            c.check("taken off", BuildTool_AnvilRmb(bodyId, plot->id, anvil, "") == AnvilRmb::Took);
            stage = 1; at = c.t; return false;
        }
        if (c.t < at + 0.5) return false;
        c.eq("the bloom is back in the inventory", (int)inventory.count("iron_bloom"), 1);
        const std::string w = Factory::Workpiece(bodyId, plot->id, anvil);
        c.check("the anvil is bare", w.empty(), w);
        c.check("bare: RMB empty-handed takes nothing", BuildTool_AnvilRmb(bodyId, plot->id, anvil, "") == AnvilRmb::Nothing);
        stage = 0;
        return true;
    });
    s.step("forged: on, struck into an ingot, the ingot taken", 40.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        Factory::View(bodyId, plot->id, anvil);
        if (stage == 0) {
            if (c.t < 0.5) return false;
            c.check("on again", BuildTool_AnvilRmb(bodyId, plot->id, anvil, "iron_bloom") == AnvilRmb::Put);
            stage = 1; at = c.t; return false;
        }
        if (stage == 1) { // a strike every 0.5 s (the hold-RMB pace is 0.45 s) until the ingot is the workpiece
            if (Factory::Workpiece(bodyId, plot->id, anvil) == "iron_ingot") {
                c.note("struck " + std::to_string(strikes) + " times");
                c.keepCapture("the finished ingot on the anvil (#1613)");
                c.check("the finished piece first", BuildTool_AnvilRmb(bodyId, plot->id, anvil, "") == AnvilRmb::Took);
                stage = 2; at = c.t; return false;
            }
            if (c.t > at + 0.5) { Factory::Strike(bodyId, plot->id, anvil); ++strikes; at = c.t; }
            return false;
        }
        if (c.t < at + 0.5) return false;
        c.eq("the ingot is in the inventory", (int)inventory.count("iron_ingot"), 1);
        const std::string w = Factory::Workpiece(bodyId, plot->id, anvil);
        c.check("the anvil is bare", w.empty(), w);
        return true;
    });
}
