// factory/rotary-sieve (#1734 PR 2): a rotary sieve on a hand crank, in the sandbox. It sifts 13 gravel into 13 sand; its
// outputs by chance are carried per machine in hundredths, so the 13th gravel (13 x 8% = 104) gives exactly one iron ore.
// The camera looks at the sieve from the front right, so the recording shows its drum turning.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "plots.h"
#include "rotationplot.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int crank = -1, sieve = -1;

bool Place(ScenarioCtx& c, PartType t, Int3 at, int orient = 0) { return BuildTool_Place(*c.body, *plot, RotationPlot_PartAt(t, at, orient), false); }
int Count(const char* item) { // what the sieve's output slots hold of an item
    Factory::View(bodyId, plot->id, sieve); Factory::Crank(bodyId, plot->id, crank); Factory::WaitTicks(2);
    int n = 0;
    for (const auto& s : Factory::View(bodyId, plot->id, sieve).outSlots) if (s.first == item) n += s.second;
    return n;
}
}

SCENARIO("factory/rotary-sieve") {
    s.tags = { "factory", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a rotary sieve on a crank, 13 gravel in it", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        const int xAxis = RotationPlot_OrientFor({ 1, 0, 0 });
        const bool ok = Place(c, PartType::HandCrank, { 0, 0, 0 }, xAxis) && Place(c, PartType::Shaft, { 1, 0, 0 }, xAxis) && Place(c, PartType::RotarySieve, { 2, 0, 0 });
        crank = plot->partAt({ 0, 0, 0 }); sieve = plot->partAt({ 2, 0, 0 });
        if (!c.check("a crank, a shaft and a rotary sieve placed", ok)) return true;
        Factory::Insert(bodyId, plot->id, sieve, "gravel", 13);
        c.camera(plot->toBody({ 5.0f, 2.6f, -2.5f }), plot->toBody({ 2.5f, 0.8f, 0.5f }), plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("cranked: it turns and sifts", 20.0, [](ScenarioCtx& c) {
        Factory::Crank(bodyId, plot->id, crank); Factory::Sync(c.world);
        if (c.t < 3.0) return false;
        Factory::View(bodyId, plot->id, sieve); Factory::Crank(bodyId, plot->id, crank); Factory::WaitTicks(2);
        const MachineView m = Factory::View(bodyId, plot->id, sieve);
        c.check("the sieve works", m.status == MachineStatus::Working, MachineStatusName(m.status));
        c.near("the sieve turns at the crank's 32 RPM", m.rpm, 32.0, 0.01);
        c.keepCapture("the rotary sieve's drum turning (#1734)");
        return true;
    });
    s.step("13 gravel: 13 sand and exactly one iron ore", 45.0, [](ScenarioCtx& c) {
        Factory::Crank(bodyId, plot->id, crank); Factory::Sync(c.world);
        const int sand = Count("sand"), iron = Count("iron_ore");
        c.waitingFor(TextFormat("%d sand, %d iron ore", sand, iron));
        if (sand < 13 && c.t < 40.0) return false; // 40 ticks a gravel at 32 RPM: 26 s
        c.eq("13 gravel sifted into 13 sand", sand, 13);
        c.eq("13 x 8% = 1.04: one iron ore", iron, 1);
        c.eq("13 x 6% = 0.78: no copper ore yet", Count("copper_ore"), 0);
        return true;
    });
}
