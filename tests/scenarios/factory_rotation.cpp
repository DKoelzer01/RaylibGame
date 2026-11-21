// factory/rotation (#1269, #159): mechanical power on the factory tick, in the sandbox. A hand crank (256 SU at 32 RPM) drives
// a line of gearboxes along x; each gearbox turns a vertical shaft down into a millstone grinding stone (128 SU each at
// 32 RPM). Two millstones: 100% stress, full speed. Four: 200%, the whole network at 50% (16 RPM: still working, slowed).
// Five: 250%, 40% (12.8 RPM): under the millstones' 16 RPM they stall ("Too slow") and make no progress. Let go of the crank:
// everything stops. The camera looks along the line, so the recording shows the shafts turning slower as the load grows.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int crank = -1, lineShaft = -1;
int mills[5] = { -1, -1, -1, -1, -1 };
int progressA = 0;

int OrientFor(Int3 zTo) { // an orientation that turns a part's local z (its axis) onto zTo
    for (int i = 0; i < 24; ++i) { const int* M = Orientation(i); if (M[2] == zTo.x && M[5] == zTo.y && M[8] == zTo.z) return i; }
    return 0;
}
bool Place(ScenarioCtx& c, PartType t, Int3 at, int orient = 0) { return BuildTool_Place(*c.body, *plot, { t, at, (uint8_t)orient }, false); }
// Millstone k at x = 3k (y 0..1), its vertical shaft at y 2, its gearbox at y 3 (axis z: ports on x and y), and the x shafts
// joining it to the previous gearbox.
bool AddMill(ScenarioCtx& c, int k) {
    const int x = 3 * k;
    bool ok = Place(c, PartType::Millstone, { x, 0, 0 }) && Place(c, PartType::Shaft, { x, 2, 0 }, OrientFor({ 0, 1, 0 })) && Place(c, PartType::Gearbox, { x, 3, 0 });
    if (k > 0) ok = ok && Place(c, k == 1 ? PartType::StressGauge : PartType::Shaft, { x - 2, 3, 0 }, OrientFor({ 1, 0, 0 })) // the first: the gauge
                       && Place(c, PartType::Shaft, { x - 1, 3, 0 }, OrientFor({ 1, 0, 0 }));
    mills[k] = plot->partAt({ x, 0, 0 });
    Factory::Insert(c.body->id, plot->id, mills[k], "stone", 40);
    return ok;
}
void Crank() { Factory::Crank(bodyId, plot->id, crank); }
RotationView Line() { return Factory::RotationAt(bodyId, plot->id, lineShaft); }
// View follows one machine at a time: switch, let two ticks pass (the crank keeps turning for 6), read.
MachineView Mill(int k) { Factory::View(bodyId, plot->id, mills[k]); Factory::WaitTicks(2); return Factory::View(bodyId, plot->id, mills[k]); }
int placed = 2;
}

SCENARIO("factory/rotation") {
    s.tags = { "factory", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 150.0;
    s.step("a crank, gearboxes and two millstones with stone", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        const bool ok = Place(c, PartType::HandCrank, { -1, 3, 0 }, OrientFor({ 1, 0, 0 })) && AddMill(c, 0) && AddMill(c, 1);
        if (!c.check("crank, gearboxes, shafts, a stress gauge and two millstones placed", ok)) return true;
        crank = plot->partAt({ -1, 3, 0 }); lineShaft = plot->partAt({ 1, 3, 0 }); // the stress gauge
        // Looking along the line from the front, a little above: the shafts, the gearboxes and the runner stones in view.
        const Vector3 mid = plot->toBody({ 3.0f * kCellSize * 2.0f, 1.0f, 0.5f });
        const Vector3 eye = plot->toBody({ 3.0f * kCellSize * 2.0f - 1.5f, 3.2f, -4.5f });
        c.camera(eye, mid, plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("cranked: two millstones at full speed", 20.0, [](ScenarioCtx& c) {
        Crank(); Factory::Sync(c.world);
        const RotationView v = Line();
        c.waitingFor(TextFormat("line %.1f RPM, stress %d%%", v.rpm, v.stressPercent));
        if (c.t < 3.0 || !v.valid) return false;
        c.near("the line turns at the crank's 32 RPM", fabsf(v.rpm), 32.0, 0.01);
        c.eq("stress 100% (2 x 128 SU on 256)", v.stressPercent, 100);
        c.eq("running at 100%", v.speedPercent, 100);
        const MachineView m = Mill(0);
        c.check("the millstone works", m.status == MachineStatus::Working, MachineStatusName(m.status));
        c.near("the millstone turns at 32 RPM", m.rpm, 32.0, 0.01);
        return true;
    });
    s.step("two more millstones: 200% stress, everything at half speed", 20.0, [](ScenarioCtx& c) {
        Crank();
        if (placed == 2) { placed = 4; if (!c.check("two more placed", AddMill(c, 2) && AddMill(c, 3))) return true; }
        Factory::Sync(c.world);
        const RotationView v = Line();
        c.waitingFor(TextFormat("line %.1f RPM, stress %d%%", v.rpm, v.stressPercent));
        if (c.t < 4.0 || !v.valid) return false;
        c.eq("stress 200% (4 x 128 SU on 256)", v.stressPercent, 200);
        c.eq("running at 50%", v.speedPercent, 50);
        c.near("the line turns at 16 RPM", fabsf(v.rpm), 16.0, 0.01);
        for (int k = 0; k < 4; ++k) {
            const MachineView m = Mill(k);
            c.check(TextFormat("millstone %d slowed but working", k), m.status == MachineStatus::Working && std::string(m.rotStatus) == "Slowed", m.rotStatus);
        }
        progressA = Mill(3).progress;
        return true;
    });
    s.step("a fifth: 250% stress, under 16 RPM the millstones stall", 20.0, [](ScenarioCtx& c) {
        Crank();
        if (placed == 4) { placed = 5; if (!c.check("a fifth placed", AddMill(c, 4))) return true; }
        Factory::Sync(c.world);
        const RotationView v = Line();
        c.waitingFor(TextFormat("line %.1f RPM, stress %d%%", v.rpm, v.stressPercent));
        if (c.t < 3.0 || !v.valid) return false;
        if (c.t < 3.1) { progressA = Mill(0).progress; return false; }
        if (c.t < 5.0) return false;
        c.eq("stress 250% (5 x 128 SU on 256)", v.stressPercent, 250);
        c.eq("running at 40%", v.speedPercent, 40);
        c.near("the line turns at 12.8 RPM", fabsf(v.rpm), 12.8, 0.01);
        const MachineView m = Mill(0);
        c.check("the millstone is too slow", std::string(m.rotStatus) == "Too slow", m.rotStatus);
        c.check("and shows no power", m.status == MachineStatus::NoPower, MachineStatusName(m.status));
        c.eq("its progress stands still", m.progress - progressA, 0);
        return true;
    });
    s.step("let go of the crank: everything stops", 10.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const RotationView v = Line();
        if (c.t < 1.0 || !v.valid) return false;
        c.eq("the line stands", fabsf(v.rpm), 0.0);
        c.check("the gauge says stopped", std::string(v.status) == "Stopped", v.status);
        const std::vector<std::pair<std::string, int>> out = Mill(0).outSlots;
        int gravel = 0;
        for (const auto& [item, n] : out) if (item == "gravel") gravel += n;
        c.gt("the first millstone ground gravel", gravel, 0);
        return true;
    });
}
