// machines/campfire-fire (#1103, the user): a campfire can be lit or put out whether or not it has something to cook.
// Lit with nothing to cook it burns its fuel (light, warmth); lit, a log put in chars at once; put out, it neither burns
// nor cooks; put out with fuel and given a log, it lights by itself and cooks (as before). Coal is the fuel (a fuel only),
// oak logs the food (a log is fuel too, but the fuel slot already holds coal, so it goes to the inputs).
#include "scenario.h"
#include "factory.h"
#include "plots.h"
#include "world.h"
#include <raymath.h>

namespace {
uint32_t plotId = 0;
int fire = -1;
MachineView ViewNow(ScenarioCtx& c) {
    Factory::View(c.body->id, plotId, fire); // ask for it, let the sim publish it
    Factory::WaitTicks(2);
    Factory::Sync(c.world);
    return Factory::View(c.body->id, plotId, fire);
}
int Out(const MachineView& v, const char* item) { int n = 0; for (const auto& s : v.outSlots) if (s.first == item) n += s.second; return n; }
}

SCENARIO("machines/campfire-fire") {
    s.tags = { "machines", "craft" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    s.step("a campfire with coal: lit with nothing to cook, it burns", 10.0, [](ScenarioCtx& c) {
        if (c.T < 0.5) return false;
        const Vector3 dir = c.spawnDirLocal();
        const float ground = c.field->radius + c.field->elevation(dir);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(Vector3Scale(dir, ground + 1.0f), 0.0f, true);
        plotId = plot->id;
        fire = plot->place({ PartType::Campfire, { 0, 0, 0 }, 0 });
        if (!c.check("campfire placed", fire >= 0)) return true;
        const Vector3 at = plot->toBody({ 0.5f, 0.5f, 0.5f });
        const Vector3 side = Vector3Normalize(Vector3CrossProduct(dir, { 0.0f, 1.0f, 0.0f }));
        c.camera(Vector3Add(Vector3Add(at, Vector3Scale(dir, 2.0f)), Vector3Scale(side, 3.0f)), at, dir); // the fire in view
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::InsertFuel(c.body->id, plotId, fire, "coal", 3);
        MachineView v = ViewNow(c);
        c.check("the campfire screen offers Light / Put out", v.canLight);
        c.check("with fuel and nothing to cook it stays out", v.status != MachineStatus::Working, MachineStatusName(v.status));
        c.eq("coal in the fuel slot", v.fuelSlot.second, 3);
        Factory::SetFire(c.body->id, plotId, fire, true);
        Factory::WaitTicks(20);
        v = ViewNow(c);
        c.check("lit with nothing to cook, it burns", v.status == MachineStatus::Working, MachineStatusName(v.status));
        c.check("it is burning a coal", v.fuelSlot.second == 2 && v.energy > 0, TextFormat("%d coal, %d ticks left", v.fuelSlot.second, v.energy));
        c.check("the world sees it working (light, smoke, stalkers)", Factory::Working(c.body->id, plotId, fire));
        return true;
    });
    s.step("lit, a log put in chars at once; the fire stays lit after", 10.0, [](ScenarioCtx& c) {
        Factory::Insert(c.body->id, plotId, fire, "oak_log", 1);
        Factory::WaitTicks(230); // the recipe takes 200 ticks
        const MachineView v = ViewNow(c);
        c.eq("charcoal made", Out(v, "charcoal"), 1);
        c.check("still lit with nothing left to cook", v.status == MachineStatus::Working, MachineStatusName(v.status));
        return true;
    });
    s.step("put out: no burn, no cooking", 10.0, [](ScenarioCtx& c) {
        Factory::SetFire(c.body->id, plotId, fire, false);
        Factory::WaitTicks(5);
        const MachineView a = ViewNow(c);
        Factory::WaitTicks(60);
        const MachineView b = ViewNow(c);
        c.check("put out, it is not working", b.status != MachineStatus::Working && b.fire == kFireOut, MachineStatusName(b.status));
        c.check("no fuel burnt while out", b.energy == a.energy && b.fuelSlot.second == a.fuelSlot.second,
                TextFormat("ticks left %d -> %d, coal %d -> %d", a.energy, b.energy, a.fuelSlot.second, b.fuelSlot.second));
        return true;
    });
    s.step("out with fuel, a log put in: it lights by itself and cooks", 10.0, [](ScenarioCtx& c) {
        Factory::Insert(c.body->id, plotId, fire, "oak_log", 1);
        Factory::WaitTicks(5);
        const MachineView a = ViewNow(c);
        c.check("given something to cook, it lights", a.status == MachineStatus::Working && a.fire == kFireAuto, MachineStatusName(a.status));
        Factory::WaitTicks(230);
        const MachineView b = ViewNow(c);
        c.eq("charcoal made in all", Out(b, "charcoal"), 2);
        c.check("lit by itself, it goes out when done", b.status != MachineStatus::Working, MachineStatusName(b.status));
        return true;
    });
}
