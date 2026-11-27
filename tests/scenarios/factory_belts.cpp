// factory/belts (#1735 PR 3): a thin belt between two chests, through the game's own path: the run is added to the plot
// (Plot::addBelt), Factory::Sync carries it to the factory thread in the layout, the first chest loads the belt from its
// contents and the belt feeds the second chest. Belts aren't drawn or placeable yet (Lighting and UI/UX, on this data), so
// the camera shows the two chests; this checks the wiring, the unit tests (Factory.*Belt*) check the motion in detail.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int chestA = -1, chestB = -1;

int Stone(int part) {
    static int last[2] = { 0, 0 };
    const MachineView v = Factory::View(bodyId, plot->id, part);
    int& l = last[part == chestA ? 0 : 1];
    if (v.valid) {
        l = 0;
        for (const auto& [item, n] : v.inSlots) if (item == "stone") l += n;
    }
    return l;
}
void Look(ScenarioCtx& c) { c.camera(plot->toBody({ 4.0f, 4.0f, 6.0f }), plot->toBody({ 4.0f, 0.5f, 0.5f }), plot->up); }
}

SCENARIO("factory/belts") {
    s.tags = { "factory", "belts" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("two chests and a belt between them", 30.0, [](ScenarioCtx& c) {
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
        const bool ok = BuildTool_Place(*c.body, *plot, { PartType::Chest, { 0, 0, 0 }, 0 }, false)
                     && BuildTool_Place(*c.body, *plot, { PartType::Chest, { 6, 0, 0 }, 0 }, false);
        chestA = plot->partAt({ 0, 0, 0 }); chestB = plot->partAt({ 6, 0, 0 });
        if (!c.check("two chests placed", ok && chestA >= 0 && chestB >= 0)) return true;
        // From A's east face to B's west face along the cell row z = 0 (slot z 4, the face centre): A's last cell x and B's first.
        int xa = 0, xb = 6;
        while (plot->partAt({ xa + 1, 0, 0 }) == chestA) ++xa;
        while (xb > 0 && plot->partAt({ xb - 1, 0, 0 }) == chestB) --xb;
        ConduitRun run; run.material = "basic"; run.points = { { xa * 8 + 6, 0, 4 }, { xb * 8, 0, 4 } };
        const int id = plot->addBelt(run);
        c.check("the belt added to the plot", id > 0);
        Factory::Sync(c.world);
        Look(c);
        return true;
    });
    // An insert before the factory has the plot comes back to the player (#1849): wait for the chest's view first.
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
        const int a = Stone(chestA), b = Stone(chestB);
        c.waitingFor(TextFormat("chest A %d stone, chest B %d", a, b));
        if (b < 12 && c.t < 35.0) return false;
        c.metric("12 stone across after", c.t, "s", false);
        c.eq("chest B holds all 12 stone", b, 12);
        c.eq("chest A is empty", a, 0);
        // 12 items at 7.5 a second, plus the 2.4-cell belt's 1.3 s ride: about 3 s
        c.lt("at the basic belt's 7.5 items/s, a few seconds (s)", c.t, 6.0);
        c.gt("and no faster than the belt allows (s)", c.t, 2.0);
        return true;
    });
}
