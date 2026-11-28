// factory/drill-deposit (#2071): a miner works the finite deposit below it. In the sandbox the scenario looks for ground
// near the start with a deposit within the miner's 128 m reach (PlanetField::depositBelow), builds an MV generator, cables
// and two miners 3 m apart there: the first mines the deposit's main ore, the second is too close to it (deposits.txt:
// 8 m spacing) and says so. The camera looks at the row from the front.
#include "scenario.h"
#include "buildtool.h"
#include "deposits.h"
#include "factory.h"
#include "oreveins.h"
#include "plots.h"
#include "terraintool.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
uint64_t bodyId = 0;
int gen = -1, minerA = -1, minerB = -1;
std::string ore; // the deposit's main ore, as the miner's variant

bool Place(ScenarioCtx& c, PartType t, Int3 at) { return BuildTool_Place(*c.body, *plot, { t, at, 0 }, false); }
MachineView Look(int part) { Factory::View(bodyId, plot->id, part); Factory::WaitTicks(2); return Factory::View(bodyId, plot->id, part); }
int Count(const MachineView& m, const std::string& item) {
    int n = 0;
    for (const auto& s : m.outSlots) if (s.first == item) n += s.second;
    return n;
}
}

SCENARIO("factory/drill-deposit") {
    s.tags = { "factory", "machines", "underground" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("ground with a deposit in a miner's reach; an MV generator and two miners on it", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const DepositDrill* dd = Deposits_Drill(Deposits_All(), "miner");
        if (!c.check("deposits.txt has the miner's reach", dd != nullptr)) return true;
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir)), north = Vector3CrossProduct(dir, east);
        // Rings of 24 m steps out to ~480 m: the first ground point with a deposit below it.
        Vector3 ground = { 0, 0, 0 };
        bool found = false;
        for (int r = 1; r <= 20 && !found; ++r)
            for (int a = 0; a < 8 * r && !found; ++a) {
                const float ang = 6.2831853f * (float)a / (float)(8 * r), d = 24.0f * (float)r;
                const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Add(Vector3Scale(east, d * cosf(ang)), Vector3Scale(north, d * sinf(ang)))));
                const Vector3 g = Vector3Scale(at, c.field->radius + c.field->elevation(at));
                VeinDeposit dep; int cx, cy, cz;
                if (c.field->depositBelow(Vector3Subtract(g, Vector3Scale(at, 1.0f)), dd->reach - 8.0f, dep, cx, cy, cz)) { ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f); found = true; }
            }
        if (!c.check("a deposit within reach of ground near the start", found)) return true;
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        bool ok = Place(c, PartType::MvGenerator, { 0, 0, 0 }) && Place(c, PartType::Miner, { 3, 0, 0 }) && Place(c, PartType::Miner, { 6, 0, 0 });
        for (int x = 0; x <= 7 && ok; ++x) ok = Place(c, PartType::Cable, { x, 0, -1 });
        gen = plot->partAt({ 0, 0, 0 }); minerA = plot->partAt({ 3, 0, 0 }); minerB = plot->partAt({ 6, 0, 0 });
        if (!c.check("a generator, two miners and a cable row placed", ok)) return true;
        Factory::Insert(bodyId, plot->id, gen, "coal", 64);
        c.camera(plot->toBody({ 4.0f, 4.0f, 9.0f }), plot->toBody({ 4.0f, 0.5f, 0.0f }), plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("the first miner works the deposit, the second is too close", 20.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const MachineView a = Look(minerA), b = Look(minerB);
        c.waitingFor(TextFormat("miner 1 %s, miner 2 %s", MachineStatusName(a.status), MachineStatusName(b.status)));
        if ((a.status != MachineStatus::Working || b.status != MachineStatus::TooClose) && c.t < 15.0) return false;
        c.check("the first miner works", a.status == MachineStatus::Working, MachineStatusName(a.status));
        c.check("the second, 3 m away, is too close", b.status == MachineStatus::TooClose, MachineStatusName(b.status));
        c.keepCapture("a miner on a finite deposit and one too close (#2071)");
        return true;
    });
    s.step("5 s on: the first has mined its deposit's ore, the second nothing", 30.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const MachineView a = Look(minerA);
        for (const auto& s : a.outSlots) if (ore.empty() && s.second > 0) ore = s.first; // they hold only the deposit's main ore
        const int n = ore.empty() ? 0 : Count(a, ore);
        c.waitingFor(TextFormat("%d %s", n, ore.c_str()));
        if (n < 4 && c.t < 25.0) return false; // 1 unit a second at grade >= 1
        c.check("it mines an ore", ore.find("_ore") != std::string::npos || ore == "coal", ore);
        c.note("the deposit's ore: " + ore);
        c.gt("the first miner mined ore", n, 3);
        int other = 0;
        for (const auto& s : Look(minerB).outSlots) other += s.second;
        c.eq("the too-close miner mined nothing", other, 0);
        return true;
    });
}

namespace {
int rig = -1;
long long leftAtStart = -1;
// The first ground point near the start, in rings of 24 m steps, with a deposit within `reach` below it.
bool GroundOverDeposit(ScenarioCtx& c, float reach, Vector3& ground) {
    const Vector3 dir = c.spawnDirLocal();
    const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir)), north = Vector3CrossProduct(dir, east);
    for (int r = 1; r <= 30; ++r)
        for (int a = 0; a < 8 * r; ++a) {
            const float ang = 6.2831853f * (float)a / (float)(8 * r), d = 24.0f * (float)r;
            const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Add(Vector3Scale(east, d * cosf(ang)), Vector3Scale(north, d * sinf(ang)))));
            const Vector3 g = Vector3Scale(at, c.field->radius + c.field->elevation(at));
            VeinDeposit dep; int cx, cy, cz;
            if (c.field->depositBelow(Vector3Subtract(g, Vector3Scale(at, 1.0f)), reach, dep, cx, cy, cz)) { ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f); return true; }
        }
    return false;
}
}

// factory/borehole-rig (#2071 PR 2): an LV borehole rig over a deposit within its 64 m reach, cabled to an LV generator,
// drills its ore at 1 unit per 2 s, and the machine view counts the deposit's units down.
SCENARIO("factory/borehole-rig") {
    s.tags = { "factory", "machines", "underground" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 120.0;
    s.step("ground with a deposit within 64 m; an LV generator and an LV borehole rig on it", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const DepositDrill* dd = Deposits_Drill(Deposits_All(), "lv_borehole_rig");
        if (!c.check("deposits.txt has the LV rig's reach", dd != nullptr)) return true;
        Vector3 ground{};
        if (!c.check("a deposit within the LV rig's reach near the start", GroundOverDeposit(c, dd->reach - 8.0f, ground))) return true;
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        bodyId = c.body->id;
        bool ok = Place(c, PartType::Generator, { 0, 0, 0 }) && Place(c, PartType::LvBoreholeRig, { 3, 0, 0 });
        for (int x = 0; x <= 4 && ok; ++x) ok = Place(c, PartType::Cable, { x, 0, -1 });
        gen = plot->partAt({ 0, 0, 0 }); rig = plot->partAt({ 3, 0, 0 });
        if (!c.check("a generator, an LV rig and a cable row placed", ok)) return true;
        Factory::Insert(bodyId, plot->id, gen, "coal", 64);
        c.camera(plot->toBody({ 3.0f, 3.5f, 7.0f }), plot->toBody({ 3.5f, 1.0f, 0.0f }), plot->up);
        Factory::Sync(c.world);
        return true;
    });
    s.step("the rig works; its view shows the deposit's units", 20.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const MachineView v = Look(rig);
        c.waitingFor(TextFormat("rig %s, %lld left", MachineStatusName(v.status), v.depositLeft));
        if (v.status != MachineStatus::Working && c.t < 15.0) return false;
        c.check("the LV rig works", v.status == MachineStatus::Working, MachineStatusName(v.status));
        c.gt("the view counts the deposit's units", v.depositLeft, 0LL);
        leftAtStart = v.depositLeft;
        c.keepCapture("an LV borehole rig drilling a deposit (#2071)");
        return true;
    });
    s.step("10 s on: ore drilled, and the deposit down by as many units", 30.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        const MachineView v = Look(rig);
        int n = 0;
        for (const auto& sl : v.outSlots) n += sl.second;
        c.waitingFor(TextFormat("%d ore, %lld left", n, v.depositLeft));
        if (c.t < 10.0) return false;
        // (The rate itself is the unit test's: Look's WaitTicks runs the factory faster than scenario time.)
        c.gt("the rig drilled ore", n, 2);
        const long long drawn = leftAtStart - v.depositLeft;
        c.check("the deposit's units went down", drawn > 0, TextFormat("%lld -> %lld", leftAtStart, v.depositLeft));
        c.check("its ore is the units drawn times the grade (1 to 1.5)", n + 1 >= drawn && n <= drawn * 3 / 2 + 1, TextFormat("%d ore for %lld units", n, drawn));
        return true;
    });
}
