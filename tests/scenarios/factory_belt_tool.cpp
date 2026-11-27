// factory/belt-tool (#2054): the belt tool lays a run between two chests. A scripted drag (BeltTool_BenchDrag: the mouse
// can't be driven) from inside chest A to inside chest B, off their face centres: the plan snaps both ends to the faces'
// item ports and checks green; a drag through chest B checks red. Placing charges conveyor items from the inventory, the
// stone rides the belt from A to B, and taking the belt up refunds them. The camera looks down on the two chests with the
// build tool's ghost and outlines on (belts themselves are drawn by #2053).
#include "scenario.h"
#include "belttool.h"
#include "buildtool.h"
#include "factory.h"
#include "inventory.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int chestA = -1, chestB = -1, xa = 0, xb = 6;
int conveyorsBefore = 0, cost = 0;

int Stone(int part) {
    static int last[2] = { -1, -1 };
    const MachineView v = Factory::View(bodyId, plot->id, part);
    int& l = last[part == chestA ? 0 : 1];
    if (v.valid) {
        l = 0;
        for (const auto& [item, n] : v.inSlots) if (item == "stone") l += n;
    }
    return l;
}
void Look(ScenarioCtx& c) { c.camera(plot->toBody({ 2.0f, 3.5f, 4.5f }), plot->toBody({ 2.0f, 0.3f, 0.5f }), plot->up); }
}

SCENARIO("factory/belt-tool") {
    s.tags = { "factory", "belts", "ui" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("two chests, build mode with the conveyor", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("2 s of play and an idle world");
        if (c.T < 2.0 || !WorldIdle()) return false;
        terrainTool.creative = false;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        const bool ok = BuildTool_Place(*c.body, *plot, { PartType::Chest, { 0, 0, 0 }, 0 }, false)
                     && BuildTool_Place(*c.body, *plot, { PartType::Chest, { 6, 0, 0 }, 0 }, false);
        chestA = plot->partAt({ 0, 0, 0 }); chestB = plot->partAt({ 6, 0, 0 });
        if (!c.check("two chests placed", ok && chestA >= 0 && chestB >= 0)) return true;
        while (plot->partAt({ xa + 1, 0, 0 }) == chestA) ++xa;
        while (xb > 0 && plot->partAt({ xb - 1, 0, 0 }) == chestB) --xb;
        inventory.add("conveyor", 20);
        conveyorsBefore = inventory.count("conveyor");
        buildTool.active = true;
        buildTool.part = (int)PartType::Conveyor;
        c.check("the conveyor in build mode is the belt tool", BeltTool_On());
        Look(c);
        return true;
    });
    s.step("a drag through chest B is red", 3.0, [](ScenarioCtx& c) {
        Look(c);
        BeltTool_BenchDrag(c.body, plot, { 2, 0, 2 }, { (xb + 3) * 8 + 4, 0, 2 }); // into A's side, out past B's far side: through it
        if (c.t < 1.5) return false;
        const ConduitPlaceResult& r = BeltTool_PlannedCheck();
        c.note(TextFormat("why: %s", r.why.c_str()));
        c.check("refused", !r.ok && BeltTool_Planned().points.size() >= 2);
        return true;
    });
    s.step("a drag from inside A to inside B snaps to their ports and is green", 3.0, [](ScenarioCtx& c) {
        Look(c);
        BeltTool_BenchDrag(c.body, plot, { 2, 0, 2 }, { xb * 8 + 4, 0, 2 }); // off the face centres (z 2, not 4)
        if (c.t < 1.5) return false;
        const ConduitRun& run = BeltTool_Planned();
        const ConduitPlaceResult& r = BeltTool_PlannedCheck();
        if (!c.check("a straight run", run.points.size() == 2)) return true;
        c.check("it starts on A's east face", run.points[0] == CInt3{ xa * 8 + 6, 0, 4 });
        c.check("it ends on B's west face", run.points[1] == CInt3{ xb * 8, 0, 4 });
        c.note(TextFormat("why: %s", r.why.c_str()));
        c.check("the check passes", r.ok);
        c.eq("both ends on an item port", (int)r.portIdx.size(), 2);
        cost = BeltTool_Cost(run);
        c.eq("it costs a conveyor a cell, rounded up", cost, (int)((Conduit_Length(run) + 7) / 8));
        return true;
    });
    s.step("released: laid, charged", 5.0, [](ScenarioCtx& c) {
        Look(c);
        const ConduitRun run = BeltTool_Planned();
        BeltTool_BenchDrag(nullptr, nullptr, {}, {});
        std::string why;
        const bool ok = BeltTool_Place(*c.body, *plot, run, true, &why);
        c.note(why.c_str());
        c.check("placed", ok);
        c.eq("one belt on the plot", (int)plot->belts.size(), 1);
        c.eq("conveyors spent", conveyorsBefore - inventory.count("conveyor"), cost);
        BeltTool_BenchShow(c.body); // the laid belt's outline, as build mode shows it
        Factory::Sync(c.world);
        return true;
    });
    s.step("12 stone into the first chest", 10.0, [](ScenarioCtx& c) {
        Look(c);
        Factory::Sync(c.world);
        c.waitingFor("the factory has the plot (chest A's view)");
        if (!Factory::View(bodyId, plot->id, chestA).valid) return false;
        Factory::Insert(bodyId, plot->id, chestA, "stone", 12);
        return true;
    });
    s.step("the stone rides the belt into the second chest", 40.0, [](ScenarioCtx& c) {
        Look(c);
        Factory::Sync(c.world);
        const int b = Stone(chestB);
        c.waitingFor(TextFormat("chest B holds %d stone", b));
        if (b < 12 && c.t < 35.0) return false;
        c.metric("12 stone across after", c.t, "s", false);
        c.eq("chest B holds all 12 stone", b, 12);
        return true;
    });
    s.step("taken up: refunded", 3.0, [](ScenarioCtx& c) {
        Look(c);
        if (c.t < 1.0) return false;
        const int id = plot->belts.empty() ? 0 : plot->belts[0].id;
        c.eq("the belt under its middle is found", BeltTool_RunAt(*plot, { (xa * 8 + 6 + xb * 8) / 2 / 2 * 2, 0, 4 }), id);
        c.check("removed", BeltTool_Remove(*c.body, *plot, id, true));
        c.eq("no belts left", (int)plot->belts.size(), 0);
        BeltTool_BenchShow(nullptr);
        c.eq("conveyors back", inventory.count("conveyor"), conveyorsBefore);
        buildTool.active = false;
        return true;
    });
}
