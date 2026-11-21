// research/labs (#1193, #1030): LV labs on a generator in the sandbox (the factory at 20x), with a small test node table so
// it runs in seconds. One lab finishes a 4-unit tier 3 node from steam packs; then two labs finish an 8-unit node in about
// the same time (linear, decided); and an LV lab never starts a tier 5 node (it idles and keeps its packs). The real table
// is put back at the end. Inert in play: labs only advance the team's queue; nothing is gated by research yet.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "plots.h"
#include "research.h"
#include "researchlabs.h"
#include "researchtab.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
int gen = -1, lab1 = -1, lab2 = -1;
double t0 = 0.0, oneLab = 0.0;
int unitsAt1s = 0;
ResearchState& Team() { return Research_Team(ResearchTab_Team()); }
int PartAt(Int3 c) { return plot ? plot->partAt(c) : -1; }
}

SCENARIO("research/labs") {
    s.tags = { "research", "machines" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };                     // (the labs in daylight)
    s.timeout = 120.0;
    s.step("an LV lab on a generator, a test node queued", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        std::string error;
        c.check("the test nodes load", Research_LoadText(
            "node t3_small | Test node | tier 3 | core | units 4 | time 10 | packs steam_pack | needs -\n"
            "node t3_big | Bigger test node | tier 3 | core | units 8 | time 10 | packs steam_pack | needs -\n"
            "node t5_far | Tier 5 test node | tier 5 | core | units 2 | time 20 | packs steam_pack lv_pack mv_pack | needs -\n", error), error);
        // The plot 6 m ahead on the ground: a generator (0,0,0), cables (2..3, 0) and (3, 1..2), a lab at (4,0,0).
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 6.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        const bool placed = BuildTool_Place(*c.body, *plot, { PartType::Generator, { 0, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::Cable, { 2, 0, 0 }, 0 }, false) && BuildTool_Place(*c.body, *plot, { PartType::Cable, { 3, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::Cable, { 3, 0, 1 }, 0 }, false) && BuildTool_Place(*c.body, *plot, { PartType::Cable, { 3, 0, 2 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::LvLab, { 4, 0, 0 }, 0 }, false);
        if (!c.check("generator, cables and a lab placed", placed)) return true;
        gen = PartAt({ 0, 0, 0 }); lab1 = PartAt({ 4, 0, 0 });
        c.camera(Vector3Add(ground, Vector3Add(Vector3Scale(at, 2.5f), Vector3Scale(east, -3.5f))), Vector3Add(ground, Vector3Scale(at, 0.5f)), at); // the labs in view
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::SetSpeed(20);
        Factory::InsertFuel(c.body->id, plot->id, gen, "coal", 64);
        Factory::Insert(c.body->id, plot->id, lab1, "steam_pack", 32);
        ResearchState& st = Team();
        st = ResearchState{};
        st.queue = { "t3_small" };
        ResearchLabs_Reset();
        t0 = c.T;
        return true;
    });
    s.step("one lab finishes the 4-unit node", 40.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        ResearchLabs_Update(c.world, c.T);
        c.waitingFor(TextFormat("t3_small: %d / 4 units", Team().progress.count("t3_small") ? Team().progress.at("t3_small") : 0));
        if (!Team().Done("t3_small")) return false;
        oneLab = c.T - t0;
        c.metric("one lab: 4 units", oneLab, "s");
        c.check("done and out of the queue", Team().queue.empty());
        // A second lab beside the cables; 8 units with two labs.
        if (!c.check("a second lab placed", BuildTool_Place(*c.body, *plot, { PartType::LvLab, { 4, 0, 2 }, 0 }, false))) return true;
        lab2 = PartAt({ 4, 0, 2 });
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::Insert(c.body->id, plot->id, lab2, "steam_pack", 32);
        Team().queue = { "t3_big" };
        t0 = c.T;
        return true;
    });
    s.step("two labs finish twice the units in about the same time", 40.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        ResearchLabs_Update(c.world, c.T);
        c.waitingFor(TextFormat("t3_big: %d / 8 units", Team().progress.count("t3_big") ? Team().progress.at("t3_big") : 0));
        if (!Team().Done("t3_big")) return false;
        const double twoLabs = c.T - t0;
        c.metric("two labs: 8 units", twoLabs, "s");
        c.between("two labs on twice the units take about as long as one lab on half (linear)", twoLabs / std::max(0.01, oneLab), 0.6, 1.5);
        // An LV lab never runs a tier 5 node: it idles and keeps its packs.
        Factory::Insert(c.body->id, plot->id, lab1, "lv_pack", 8);
        Factory::Insert(c.body->id, plot->id, lab1, "mv_pack", 8);
        Team().queue = { "t5_far" };
        t0 = c.T;
        return true;
    });
    s.step("an LV lab leaves a tier 5 node alone", 20.0, [](ScenarioCtx& c) {
        Factory::Sync(c.world);
        ResearchLabs_Update(c.world, c.T);
        // (a tier 3 unit made just before the switch may still be collected in the first second: it is spent, as decided)
        if (c.T - t0 < 1.0) { unitsAt1s = ResearchLabs_UnitsCounted(); return false; }
        if (c.T - t0 < 6.0) return false;
        c.check("no units made while only tier 5 is queued", ResearchLabs_UnitsCounted() == unitsAt1s);
        const MachineView v = Factory::View(c.body->id, plot->id, lab1);
        c.check("the LV lab is not working", v.status != MachineStatus::Working, MachineStatusName(v.status));
        c.check("t5_far untouched", !Team().progress.count("t5_far") && !Team().Done("t5_far"));
        Team() = ResearchState{};
        Factory::SetSpeed(1);
        Research_Load(); // (the real table back)
        return true;
    });
}
