// quests/campfire-no-screen (#1728): tier 0's "get fired_crucible 1" did not register unless the player opened the campfire's
// screen and took the crucible out. A machine's finished recipe now counts for its owner's quest book whether or not anyone
// takes the output: place a campfire, put a crucible mould and fuel in it without opening any screen, wait for the firing,
// and the task is done while the fired crucible is still in the campfire (nothing was taken, nothing in the inventory).
#include "scenario.h"
#include "factory.h"
#include "inventory.h"
#include "plots.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>

namespace {
uint32_t plotId = 0;
int fire = -1;
constexpr int kTask = 6; // tiers.txt tier 0: "task | get fired_crucible 1" is its 7th task
MachineView ViewNow(ScenarioCtx& c) {
    Factory::View(c.body->id, plotId, fire); // ask for it, let the sim publish it
    Factory::WaitTicks(2);
    Factory::Sync(c.world);
    return Factory::View(c.body->id, plotId, fire);
}
int Out(const MachineView& v, const char* item) { int n = 0; for (const auto& s : v.outSlots) if (s.first == item) n += s.second; return n; }
}

SCENARIO("quests/campfire-no-screen") {
    s.tags = { "quests", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("a campfire, a crucible mould and fuel put in without a screen", 10.0, [](ScenarioCtx& c) {
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
        QuestBook_Update();
        c.check("the task is not done before", !QuestBook_TaskDone(0, kTask));
        Factory::InsertFuel(c.body->id, plotId, fire, "coal", 2);
        Factory::Insert(c.body->id, plotId, fire, "crucible_mold", 1);
        const MachineView v = ViewNow(c);
        c.check("the campfire is cooking the mould", v.status == MachineStatus::Working, MachineStatusName(v.status));
        return true;
    });
    s.step("fired: the task is done and nobody took the crucible", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("the firing (400 ticks)");
        Factory::WaitTicks(430);
        const MachineView v = ViewNow(c);
        QuestBook_Update();
        c.eq("the fired crucible is still in the campfire", Out(v, "fired_crucible"), 1);
        c.eq("nothing in the inventory", inventory.count("fired_crucible"), 0);
        c.eq("nothing taken: not obtained", inventory.obtainedTotal("fired_crucible"), 0.0);
        c.check("tier 0's \"get fired_crucible 1\" is done", QuestBook_TaskDone(0, kTask));
        return true;
    });
}
