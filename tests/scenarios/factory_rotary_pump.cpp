// factory/rotary-pump (#1734 PR 3): a rotary pump on a hand crank, in the sandbox, over a pit with 3 m^3 of water poured in
// (a finite pond), with a boiler beside it. At the crank's 32 RPM (its rated speed) it lifts a water unit (0.5 m^3) every
// 25 s, 20 L/s: the boiler gets the first 500 L after about 25 s, and the pond loses what the boiler gained (no water made
// from nothing). Uncranked, it stops. The camera looks at the pump from the front right, so the recording shows its wheel.
// The poured water pools at the pit's rim rather than in it (#2005); it is within the intake's reach and the 3 m measure
// either way.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "fluid.h"
#include "plots.h"
#include "rotationplot.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int crank = -1, pump = -1, boiler = -1;
Vector3 pit{};       // body-local: the pit's middle, under the pump
float pit0 = 0.0f;   // m^3 in the pit when cranking starts

bool Place(ScenarioCtx& c, PartType t, Int3 at, int orient = 0) { return BuildTool_Place(*c.body, *plot, RotationPlot_PartAt(t, at, orient), false); }
// Each frame asks again; a view fills from the tick after it is asked. No WaitTicks here: it runs extra factory ticks, and
// called every frame it ran the pump ~6x faster than sim time, so its 25 s showed as 6.5 s.
MachineView Watch(int part) {
    static MachineView last[2]; // the pump's, the boiler's
    MachineView& l = last[part == pump ? 0 : 1];
    const MachineView v = Factory::View(bodyId, plot->id, part);
    if (v.valid) l = v;
    return l;
}
void LookAtPump(ScenarioCtx& c) { c.camera(plot->toBody({ 5.5f, 2.8f, -2.8f }), plot->toBody({ 2.0f, 0.6f, 0.8f }), plot->up); }
}

SCENARIO("factory/rotary-pump") {
    s.tags = { "factory", "machines", "fluid" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 150.0;
    s.step("a rotary pump on a crank over a pit, a boiler beside it", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("2 s of play and an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        // The pit first (as factory/rotation-sources digs its trench): a 2.2 m ball carved 1.6 m under the pump's floor (its
        // intake reaches 0.7 m down, then 2.5 m round).
        pit = plot->toBody({ 1.5f, -1.6f, 0.5f });
        ApplyTerrainEdit(*c.body, TerrainEdit{ pit, 2.2f, -1.0f, true });
        const int xAxis = RotationPlot_OrientFor({ 1, 0, 0 });
        // The pump's cells are x 2-3, z 0-1; the boiler's x 2-3, z 2-3: it touches the pump face to face.
        const bool ok = Place(c, PartType::HandCrank, { 0, 0, 0 }, xAxis) && Place(c, PartType::Shaft, { 1, 0, 0 }, xAxis) &&
                        Place(c, PartType::RotaryPump, { 2, 0, 0 }) && Place(c, PartType::Boiler, { 2, 0, 2 });
        crank = plot->partAt({ 0, 0, 0 }); pump = plot->partAt({ 2, 0, 0 }); boiler = plot->partAt({ 2, 0, 2 });
        if (!c.check("a crank, a shaft, a rotary pump and a boiler placed", ok && boiler >= 0)) return true;
        LookAtPump(c);
        Factory::Sync(c.world);
        return true;
    });
    s.step("pour 3 m^3 into the pit", 30.0, [](ScenarioCtx& c) {
        LookAtPump(c);
        c.waitingFor("the pit dug (an idle world)");
        if (c.t < 5.0 || !WorldIdle()) return false; // the fluid sim sees the carve a few seconds after the chunks
        Fluid::Pour(*c.body, plot->toBody({ 1.5f, -1.2f, -0.6f }), 3.0f); // on the open side of the pit, clear of every part
        return true;
    });
    s.step("the water settles in the pit", 20.0, [](ScenarioCtx& c) {
        LookAtPump(c);
        if (c.t < 6.0) return false;
        pit0 = Fluid::WaterAt(*c.body, pit, 3.0f);
        c.gt("the pit holds most of the 3 m^3 (m^3)", pit0, 2.0);
        return true;
    });
    s.step("cranked: it turns at 32 RPM", 10.0, [](ScenarioCtx& c) {
        LookAtPump(c);
        Factory::Crank(bodyId, plot->id, crank); Factory::Sync(c.world);
        const MachineView m = Watch(pump);
        if (c.t < 3.0) return false;
        c.check("the pump works", m.status == MachineStatus::Working, MachineStatusName(m.status));
        c.near("the pump turns at the crank's 32 RPM", m.rpm, 32.0, 0.01);
        c.keepCapture("the rotary pump's wheel turning (#1734)");
        return true;
    });
    s.step("the first 500 L reach the boiler after about 25 s", 45.0, [](ScenarioCtx& c) {
        LookAtPump(c);
        Factory::Crank(bodyId, plot->id, crank); Factory::Sync(c.world);
        const MachineView b = Watch(boiler);
        const double cranked = c.t + 3.0; // this step starts 3 s into cranking
        c.waitingFor(TextFormat("the boiler holds %d L after %.0f s of cranking", b.water, cranked));
        if (b.water <= 0 && c.t < 40.0) return false;
        c.metric("first water unit after", cranked, "s", false);
        c.eq("one water unit (L)", b.water, 500);
        c.gt("not before 20 L/s allows (s)", cranked, 23.0);
        c.lt("and not much later (s)", cranked, 30.0);
        return true;
    });
    s.step("the pit lost what the boiler gained", 10.0, [](ScenarioCtx& c) {
        LookAtPump(c);
        Factory::Crank(bodyId, plot->id, crank); Factory::Sync(c.world);
        if (c.t < 2.0) return false; // the take reaches the published water
        const float lost = pit0 - Fluid::WaterAt(*c.body, pit, 3.0f);
        c.note(TextFormat("the pit held %.2f m^3, lost %.2f m^3", pit0, lost));
        c.near("the pit lost one water unit (m^3)", lost, 0.5, 0.15);
        return true;
    });
    s.step("uncranked, it stops", 15.0, [](ScenarioCtx& c) {
        LookAtPump(c);
        const MachineView m = Watch(pump);
        if (c.t < 5.0) return false;
        c.near("the pump stands still (RPM)", m.rpm, 0.0, 0.01);
        c.check("and does not work", m.status != MachineStatus::Working, MachineStatusName(m.status));
        return true;
    });
}
