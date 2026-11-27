// factory/rotation-age (#1734): the rotation-age machines on a hand crank each, in the sandbox. A trip hammer works an iron
// bloom into an ingot; a mixer stirs copper and tin dust 3:1 into bronze dust. Each line: crank, shaft, machine, along x;
// one crank's 256 SU at 32 RPM drives one machine at its rated speed. The camera looks at both from the front, so the
// recording shows the camshaft and the mixer's paddles turning.
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
int crank[2] = { -1, -1 }, machine[2] = { -1, -1 };

int OrientFor(Int3 zTo) { return RotationPlot_OrientFor(zTo); }
bool Place(ScenarioCtx& c, PartType t, Int3 at, int orient = 0) { return BuildTool_Place(*c.body, *plot, RotationPlot_PartAt(t, at, orient), false); }
void CrankBoth() { for (int k : crank) Factory::Crank(bodyId, plot->id, k); }
}

SCENARIO("factory/rotation-age") {
    s.tags = { "factory", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 60.0;
    s.step("a trip hammer and a mixer, each on a crank", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        const PartType kinds[2] = { PartType::TripHammer, PartType::Mixer };
        bool ok = true;
        for (int i = 0; i < 2; ++i) { // z 0 and 3: crank at x 0, shaft at x 1, the machine from x 2
            const int z = 3 * i;
            ok = ok && Place(c, PartType::HandCrank, { 0, 0, z }, OrientFor({ 1, 0, 0 })) && Place(c, PartType::Shaft, { 1, 0, z }, OrientFor({ 1, 0, 0 }))
                 && Place(c, kinds[i], { 2, 0, z });
            crank[i] = plot->partAt({ 0, 0, z }); machine[i] = plot->partAt({ 2, 0, z });
        }
        if (!c.check("cranks, shafts, a trip hammer and a mixer placed", ok)) return true;
        Factory::Insert(bodyId, plot->id, machine[0], "iron_bloom", 4);
        Factory::Insert(bodyId, plot->id, machine[1], "copper_dust", 6);
        Factory::Insert(bodyId, plot->id, machine[1], "tin_dust", 2);
        // From the front right, a little above: both machines side by side, their models on top.
        c.camera(plot->toBody({ 6.5f, 3.0f, -1.5f }), plot->toBody({ 2.5f, 1.0f, 2.0f }), plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("cranked: both turn and work", 20.0, [](ScenarioCtx& c) {
        CrankBoth(); Factory::Sync(c.world);
        if (c.t < 3.0) return false;
        for (int i = 0; i < 2; ++i) {
            Factory::View(bodyId, plot->id, machine[i]); CrankBoth(); Factory::WaitTicks(2);
            const MachineView m = Factory::View(bodyId, plot->id, machine[i]);
            const char* name = i == 0 ? "the trip hammer" : "the mixer";
            c.check(std::string(name) + " works", m.status == MachineStatus::Working, MachineStatusName(m.status));
            c.near(std::string(name) + " turns at the crank's 32 RPM", m.rpm, 32.0, 0.01);
        }
        c.keepCapture("the trip hammer's camshaft and the mixer's paddles turning (#1734)");
        return true;
    });
    s.step("made: an ingot from every bloom, bronze dust from the copper and tin", 40.0, [](ScenarioCtx& c) {
        CrankBoth(); Factory::Sync(c.world);
        auto made = [&](int i, const char* item) { // what one machine's output slots hold of an item
            Factory::View(bodyId, plot->id, machine[i]); CrankBoth(); Factory::WaitTicks(2);
            int n = 0;
            for (const auto& s : Factory::View(bodyId, plot->id, machine[i]).outSlots) if (s.first == item) n += s.second;
            return n;
        };
        const int ingots = made(0, "iron_ingot"), dust = made(1, "bronze_dust");
        c.waitingFor(TextFormat("%d ingots, %d bronze dust", ingots, dust));
        if ((ingots < 4 || dust < 8) && c.t < 35.0) return false; // 4 blooms at 80 ticks; 2 mixes at 200 ticks
        c.eq("the trip hammer worked all 4 blooms", ingots, 4);
        c.eq("the mixer made 8 bronze dust (6 copper + 2 tin)", dust, 8);
        return true;
    });
}
