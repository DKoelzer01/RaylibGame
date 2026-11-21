// factory/summary-ingame (#1261): the summary sim in the running game. A chain (a chest of stone -> a macerator -> a chest)
// on a generator runs near for 90 s of factory time; the camera then goes 2.5 km away (the plot leaves "near": 60 s of
// measuring, then the summary advances it every 10 s), stays away for 180 s, and comes back. The sand made while away is
// at most what the near rate gives over that time (plus one craft) and at least 80% of it (the measuring minute takes the
// 20-tick step); while away the Statistics view marks the plot's sand as estimated. The factory runs at 20x.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
long long tick0 = 0, nearTicks = 0, awayTick = 0;
long long sandAtWarm = 0, sandAtLeave = 0;
bool estimatedSeen = false;
Vector3 ground0{}, up0{}, east0{};
long long Sand(bool* estimated = nullptr) {
    const FactoryStatsView v = Factory::StatsView(StatWindow::Minute, true, bodyId, plot->id);
    for (const StatLine& l : v.items) if (l.item == "sand") { if (estimated) *estimated = l.estimated; return l.made; }
    return 0;
}
long long Ticks() { return Factory::Stats().ticks; }
}

SCENARIO("factory/summary-ingame") {
    s.tags = { "factory", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("a chest of stone, a macerator and a chest on a generator", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 dir = c.spawnDirLocal();
        east0 = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east0, 6.0f)));
        ground0 = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        up0 = at;
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground0, 0.0f, true);
        bodyId = c.body->id;
        const bool placed = BuildTool_Place(*c.body, *plot, { PartType::Generator, { 0, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::Chest, { 3, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::Macerator, { 6, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::Chest, { 9, 0, 0 }, 0 }, false);
        bool lines = true;
        for (int x = 0; x <= 7; ++x) lines = lines && BuildTool_Place(*c.body, *plot, { PartType::Cable, { x, 0, -1 }, 0 }, false);
        for (int x = 3; x <= 10; ++x) lines = lines && BuildTool_Place(*c.body, *plot, { PartType::Conveyor, { x, 0, 2 }, 0 }, false);
        if (!c.check("the chain placed", placed && lines)) return true;
        c.camera(Vector3Add(ground0, Vector3Add(Vector3Scale(up0, 4.0f), Vector3Scale(east0, -4.0f))), ground0, up0);
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::SetSpeed(20);
        Factory::InsertFuel(bodyId, plot->id, plot->partAt({ 0, 0, 0 }), "coal", 64);
        Factory::Insert(bodyId, plot->id, plot->partAt({ 3, 0, 0 }), "stone", 640);
        tick0 = Ticks();
        return true;
    });
    s.step("near for 90 s of factory time (the last 30 s give the near rate)", 30.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const long long t = Ticks() - tick0;
        c.waitingFor(TextFormat("%lld / %d ticks near", t, 90 * kTicksPerSecond));
        if (t >= 60 * kTicksPerSecond && sandAtWarm == 0) { sandAtWarm = Sand(); nearTicks = Ticks(); }
        if (t < 90 * kTicksPerSecond) return false;
        const long long sand = Sand();
        nearTicks = Ticks() - nearTicks;
        sandAtLeave = sand;
        c.metric("near rate (sand a minute)", (double)(sand - sandAtWarm) * 60.0 * kTicksPerSecond / std::max(1LL, nearTicks), "/min", false);
        sandAtWarm = sand - sandAtWarm; // (now: sand made in the measured near ticks)
        // The camera leaves: 2.5 km east, 300 m up (the plot is "near" within 1000 m).
        c.camera(Vector3Add(ground0, Vector3Add(Vector3Scale(up0, 300.0f), Vector3Scale(east0, 2500.0f))), Vector3Add(ground0, Vector3Scale(east0, 2600.0f)), up0);
        awayTick = Ticks();
        return true;
    });
    s.step("away for 180 s of factory time (measuring, then the summary)", 40.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const long long t = Ticks() - awayTick;
        bool est = false;
        Sand(&est);
        estimatedSeen = estimatedSeen || est;
        c.waitingFor(TextFormat("%lld / %d ticks away", t, 180 * kTicksPerSecond));
        if (t < 180 * kTicksPerSecond) return false;
        c.check("while away the plot's sand is marked estimated (the summary sim ran it)", estimatedSeen);
        c.camera(Vector3Add(ground0, Vector3Add(Vector3Scale(up0, 4.0f), Vector3Scale(east0, -4.0f))), ground0, up0); // back
        return true;
    });
    s.step("back near: what it made while away", 10.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        if (c.t < 1.0) return false; // (the near message reaches the sim and advances the plot to now)
        const long long away = Ticks() - awayTick;
        const long long made = Sand() - sandAtLeave;
        const double expected = (double)sandAtWarm * away / std::max(1LL, nearTicks);
        c.metric("sand made away", (double)made, "", false);
        c.metric("at the near rate", expected, "", false);
        c.between("away / near rate (summary at or under the full sim, not much under)", made / std::max(1.0, expected), 0.80, 1.0 + 2.0 / std::max(1.0, expected));
        Factory::SetSpeed(1);
        return true;
    });
}
