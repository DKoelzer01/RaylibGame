// factory/stats (#1199, --test-factory-stats): the production statistics count exactly what the factory did. A generator on a
// cable powers a macerator (stone 1 -> sand 2, 100 ticks at 4 EU/t) in the sandbox at 20x. 10 stone and 16 coal are put in by
// hand; when the macerator is idle again: 10 stone used and 20 sand made, by hand 10 stone + 16 coal in, 20 sand out once
// taken; the coal burnt plus what is left is what went in; EU used = 10 x 100 x 4 exactly and the generator made at least
// that; the ticks per status add up to the ticks that passed; the minute window's sparkline adds up to its rate; and the tier
// clock ran with the sim time.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "gameui.h"
#include "statstab.h"
#include "plots.h"
#include "survival.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <numeric>
#include <string>

namespace {
Plot* plot = nullptr;
int gen = -1, mac = -1;
uint64_t bodyId = 0;
long long ticksA = 0, statusA = 0;
double tierA = 0.0, simA = 0.0;
const StatLine* Line(const FactoryStatsView& v, const char* item) { for (const StatLine& l : v.items) if (l.item == item) return &l; return nullptr; }
const StatKind* Kind(const FactoryStatsView& v, PartType t) { for (const StatKind& k : v.kinds) if (k.type == (int)t) return &k; return nullptr; }
long long StatusTicks(const StatKind* k) { return k ? std::accumulate(std::begin(k->ticks), std::end(k->ticks), 0LL) : 0; }
FactoryStatsView PlotView() { return Factory::StatsView(StatWindow::Minute, true, bodyId, plot->id); }
}

SCENARIO("factory/stats") {
    s.tags = { "factory", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a generator, a cable and a macerator; stone and coal put in by hand", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        const bool placed = BuildTool_Place(*c.body, *plot, { PartType::Generator, { 0, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::Cable, { 2, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::Macerator, { 3, 0, 0 }, 0 }, false);
        if (!c.check("generator, cable and macerator placed", placed)) return true;
        gen = plot->partAt({ 0, 0, 0 }); mac = plot->partAt({ 3, 0, 0 }); bodyId = c.body->id;
        c.camera(Vector3Add(ground, Vector3Add(Vector3Scale(at, 2.5f), Vector3Scale(east, -3.0f))), Vector3Add(ground, Vector3Scale(at, 0.5f)), at);
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::SetSpeed(20);
        Factory::InsertFuel(c.body->id, plot->id, gen, "coal", 16);
        Factory::Insert(c.body->id, plot->id, mac, "stone", 10);
        Factory::WaitTicks(2);
        const FactoryStatsView v = PlotView();
        ticksA = Factory::Stats().ticks; statusA = StatusTicks(Kind(v, PartType::Macerator));
        tierA = QuestBook_TierSeconds(QuestBook_Tier()); simA = SimTime();
        return true;
    });
    s.step("the macerator works the stone through", 40.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const FactoryStatsView v = PlotView();
        const StatLine* sand = Line(v, "sand");
        c.waitingFor(TextFormat("sand made %lld / 20", sand ? sand->made : 0LL));
        return sand && sand->made >= 20 && Factory::View(c.body->id, plot->id, mac).status != MachineStatus::Working;
    });
    s.step("the counts match what happened", 10.0, [](ScenarioCtx& c) {
        Factory::Extract(c.body->id, plot->id, mac);
        Factory::WaitTicks(2);
        Factory::Sync(c.world);
        const FactoryStatsView v = PlotView();
        const long long ticksB = Factory::Stats().ticks;
        const StatLine* stone = Line(v, "stone");
        const StatLine* sand = Line(v, "sand");
        const StatLine* coal = Line(v, "coal");
        if (!c.check("stone, sand and coal are counted", stone && sand && coal)) return true;
        c.check("10 stone used", stone->used == 10, TextFormat("%lld", stone->used));
        c.check("20 sand made", sand->made == 20, TextFormat("%lld", sand->made));
        c.check("10 stone put in by hand", stone->handIn == 10, TextFormat("%lld", stone->handIn));
        c.check("16 coal put in by hand", coal->handIn == 16, TextFormat("%lld", coal->handIn));
        c.check("20 sand taken out by hand", sand->handOut == 20, TextFormat("%lld", sand->handOut));
        Factory::View(c.body->id, plot->id, gen); Factory::WaitTicks(2); Factory::Sync(c.world); // (View follows one machine: the first call after a switch is stale)
        const MachineView g = Factory::View(c.body->id, plot->id, gen);
        c.check("coal burnt + coal left = coal put in", coal->used + g.fuelSlot.second == 16, TextFormat("burnt %lld, left %d", coal->used, g.fuelSlot.second));
        const StatKind* macK = Kind(v, PartType::Macerator);
        const StatKind* genK = Kind(v, PartType::Generator);
        const double euUsed = macK ? macK->euUsedPerTick : 0.0, euMade = genK ? genK->euMadePerTick : 0.0;
        c.metric("EU used per tick (minute window)", euUsed, "EU/t");
        c.metric("EU made per tick (minute window)", euMade, "EU/t");
        // Totals over the minute: 10 crafts x 100 ticks x 4 EU/t, all inside the last minute of factory time (50 s of work).
        const double minuteTicks = v.coveredMinutes * 60.0 * kTicksPerSecond;
        c.between("EU used = 10 x 100 x 4", euUsed * minuteTicks, 3999.0, 4001.0);
        c.check("the generator made at least what was used", euMade * minuteTicks >= euUsed * minuteTicks - 0.5);
        const long long statusB = StatusTicks(macK);
        c.between("the macerator's ticks per status add up to the ticks passed", (double)(statusB - statusA) - (double)(ticksB - ticksA), -2.0, 2.0);
        float sparkSum = 0.0f;
        for (float x : sand->madeSpark) sparkSum += x;
        c.between("the minute's sparkline adds up to its rate", sparkSum - sand->madePerMin * v.coveredMinutes, -0.01, 0.01);
        const double tierRan = QuestBook_TierSeconds(QuestBook_Tier()) - tierA, simRan = SimTime() - simA;
        c.between("the tier clock ran with the sim time", tierRan / std::max(0.01, simRan), 0.95, 1.05);
        Factory::SetSpeed(1);
        GameUI_ShowTab("Statistics"); // (the tab with these figures in it: a screenshot with RAYLIBGAME_STATS_SHOT=<png>)
        return true;
    });
    s.step("the Statistics tab shows them", 5.0, [](ScenarioCtx& c) {
        if (c.t < 1.0) return false;
        if (const char* shot = getenv("RAYLIBGAME_STATS_SHOT")) TakeScreenshot(shot); // (raylib writes it in the working directory, by its file name)
        c.check("the tab lists the sand", [] { for (const StatLine& l : StatsTab_View().items) if (l.item == "sand") return true; return false; }());
        return true;
    });
}
