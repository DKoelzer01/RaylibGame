// factory/alt-recipe (#2072 slice 2b): a ruin cache's alternate recipe, chosen for one machine through the real input
// path (E on the machine, on foot). Before the team has one, E on a coke oven is not the oven's; with packed_coke the
// panel offers it and E picks it; the oven then makes 3 coke from 5 coal (the base recipe: 1 from 2); E again goes back
// to the base recipes.
#include "scenario.h"
#include "buildtool.h"
#include "caches.h"
#include "factory.h"
#include "inventory.h"
#include "player.h"
#include "playerinput.h"
#include "playerrecord.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
uint32_t plotId = 0;
int oven = -1;
constexpr uint64_t kE = 1ull << (int)Act::Use;
const Int3 kOven{ 0, 0, 0 };

Plot* ThePlot(ScenarioCtx& c) { return c.body->plots ? c.body->plots->byId(plotId) : nullptr; }
Vector3 CellCentre(Plot& p, Int3 cell) { return p.toBody({ (cell.x + 0.5f) * kCellSize, (cell.y + 0.4f) * kCellSize, (cell.z + 0.5f) * kCellSize }); }
// Each frame: the player's view on the oven (the player stands; only the view turns).
void Aim(ScenarioCtx& c) {
    Plot* p = ThePlot(c);
    if (!p) return;
    const Vector3 at = CellCentre(*p, kOven);
    const Vector3 feet = c.body->toLocalD(player.position).toFloat();
    const Vector3 eye = Vector3Add(feet, Vector3Scale(Vector3Normalize(feet), 1.6f));
    player.forward = c.body->dirToWorld(Vector3Normalize(Vector3Subtract(at, eye)));
}
// One E press at the step's start (held 0.2 s), then `after` seconds.
bool PressE(ScenarioCtx& c, double after) {
    if (c.t < 0.2) { gFakeHeld |= kE; return false; }
    gFakeHeld &= ~kE;
    return c.t >= 0.2 + after;
}
MachineView OvenView(ScenarioCtx& c) { return Factory::View(c.body->id, plotId, oven); }
int Count(const std::vector<std::pair<std::string, int>>& slots, const char* item) { int n = 0; for (const auto& [i, k] : slots) if (i == item) n += k; return n; }
int Team() { int t = Teams_Of(PlayerId_Local()); return t ? t : Teams_CreateSolo(PlayerId_Local(), "Player"); }
} // namespace

SCENARIO("factory/alt-recipe") {
    s.tags = { "factory", "machines", "adventure" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a coke oven on the ground, the player on foot before it", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const float ground = c.field->radius + c.field->elevation(dir);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        Plot* plot = c.body->plots->plotFor(Vector3Scale(dir, ground + 0.05f), 0.0f, true);
        plotId = plot->id;
        if (!c.check("coke oven placed", BuildTool_Place(*c.body, *plot, { PartType::CokeOven, kOven, 0 }, false))) return true;
        oven = plot->partAt(kOven);
        const Vector3 mid = CellCentre(*plot, kOven);
        const Vector3 up = Vector3Normalize(mid), out = Vector3Normalize(Vector3CrossProduct(up, plot->east));
        player.flyMode = false;
        player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Add(Vector3Scale(up, Vector3Length(mid) + 0.3f), Vector3Scale(out, 3.0f))));
        Caches_ResetTeams();
        gScenarioInput = true; // (the scripted E is used: scenarios are benches, which take no input otherwise)
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        return true;
    });
    s.step("without the alternate, E on the oven is not the oven's", 10.0, [](ScenarioCtx& c) {
        Aim(c);
        c.waitingFor("grounded, the oven aimed");
        if (c.t < 1.0 || !player.grounded || BuildTool_AimedType() != PartType::CokeOven) return false;
        c.check("the oven has an alternate in the data", !Factory::AltsOf(PartType::CokeOven).empty());
        c.check("E is not the oven's", !BuildTool_AimingUsable());
        Caches_AddTeamAlt(Team(), "packed_coke"); // as a cache's alt reward does (Caches_ApplyTeam)
        c.check("E is the oven's once the team has packed coke", BuildTool_AimingUsable());
        return true;
    });
    s.step("E picks packed coke", 10.0, [](ScenarioCtx& c) {
        Aim(c);
        if (!PressE(c, 0.5)) return false;
        OvenView(c); Factory::WaitTicks(2); // (View follows one machine: the first call after a switch is stale)
        c.check("the oven runs packed coke", OvenView(c).alt == "packed_coke", OvenView(c).alt.c_str());
        Factory::Insert(c.body->id, plotId, oven, "coal", 10);
        Factory::SetSpeed(20);
        return true;
    });
    s.step("3 coke from 5 coal", 30.0, [](ScenarioCtx& c) {
        Aim(c);
        const MachineView v = OvenView(c);
        c.waitingFor(TextFormat("coke %d / 3", Count(v.outSlots, "coke")));
        if (Count(v.outSlots, "coke") < 3) return false;
        Factory::SetSpeed(1);
        Factory::WaitTicks(2);
        const MachineView w = OvenView(c);
        // (the second craft may have started: its coal is taken only when it finishes, #193)
        c.eq("3 coke from the first craft", Count(w.outSlots, "coke"), 3);
        c.eq("5 coal taken for it", Count(w.inSlots, "coal"), 5);
        return true;
    });
    s.step("E again: the base recipes", 10.0, [](ScenarioCtx& c) {
        Aim(c);
        if (!PressE(c, 0.5)) return false;
        Factory::WaitTicks(2);
        c.check("no alternate chosen", OvenView(c).alt.empty(), OvenView(c).alt.c_str());
        return true;
    });
}
