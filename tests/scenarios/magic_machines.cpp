// magic/machines (#1252, #1327): with RAYLIBGAME_MAGIC on, in the sandbox, on the ley node nearest the start: a mana siphon
// on the node fills empty jars with the node's kind; a mana dynamo burns a jar; standing by a research table lets the Tome
// tab spend items (its reach). The camera watches the three machines on the node (the factory at 20x).
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "magic.h"
#include "magicmachines.h"
#include "manafield.h"
#include "plots.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>
#include <vector>

namespace {
Plot* plot = nullptr;
int siphon = -1, dynamo = -1;
std::string jar;
Vector3 nodeGround{}, nodeUp{}, east{};
int PartAt(Int3 c) { return plot ? plot->partAt(c) : -1; }
int Count(const std::vector<std::pair<std::string, int>>& slots, const std::string& item) { int n = 0; for (const auto& [id, k] : slots) if (id == item) n += k; return n; }
}

SCENARIO("magic/machines") {
    s.tags = { "magic", "machines", "factory" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.env = { { "RAYLIBGAME_MAGIC", "1" } };
    s.timeout = 90.0;
    s.step("a siphon, a dynamo and a research table on the nearest ley node", 15.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        const ManaField* f = Magic_Field(c.world);
        if (!c.check("the mana field (RAYLIBGAME_MAGIC)", f != nullptr)) return true;
        std::vector<ManaNode> nodes;
        const Vector3 start = c.spawnDirLocal();
        f->nodesNear(start, 3000.0f, nodes);
        if (!c.check("a node within 3 km of the start", !nodes.empty())) return true;
        const ManaNode* best = &nodes[0];
        for (const ManaNode& n : nodes) if (Vector3DotProduct(n.dir, start) > Vector3DotProduct(best->dir, start)) best = &n;
        nodeUp = best->dir;
        nodeGround = Vector3Scale(nodeUp, c.field->radius + c.field->elevation(nodeUp) + 0.05f);
        east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, nodeUp));
        jar = MagicMachines_SiphonJar(*f, nodeGround);
        c.metric("node strength", best->strength, "");
        c.check("the node has a jar", jar != "jar_none", jar);
        terrainTool.creative = true;
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(nodeGround, 0.0f, true);
        const bool placed = BuildTool_Place(*c.body, *plot, { PartType::ManaSiphon, { 0, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::ManaDynamo, { 4, 0, 0 }, 0 }, false)
            && BuildTool_Place(*c.body, *plot, { PartType::ResearchTable, { 0, 0, 4 }, 0 }, false);
        if (!c.check("siphon, dynamo and research table placed", placed)) return true;
        siphon = PartAt({ 0, 0, 0 }); dynamo = PartAt({ 4, 0, 0 });
        Factory::Sync(c.world);
        Factory::WaitTicks(3);
        Factory::SetSpeed(20);
        Factory::Insert(c.body->id, plot->id, siphon, "mana_jar", 4);
        return true;
    });
    s.step("the siphon fills jars with the node's kind", 30.0, [](ScenarioCtx& c) {
        // The camera 5 m off the node and 3 m up, the machines in view (and the player within the table's reach).
        c.camera(Vector3Add(nodeGround, Vector3Add(Vector3Scale(east, -4.0f), Vector3Scale(nodeUp, 3.0f))), Vector3Add(nodeGround, Vector3Scale(nodeUp, 0.5f)), nodeUp);
        Factory::Sync(c.world);
        MagicMachines_Update(c.world, c.T);
        const MachineView v = Factory::View(c.body->id, plot->id, siphon);
        const int filled = Count(v.outSlots, jar);
        c.waitingFor(TextFormat("%s: %d", jar.c_str(), filled));
        if (filled < 2) return false;
        c.check("only the node's kind", Count(v.outSlots, "jar_none") == 0);
        c.check("the research table is in reach", MagicMachines_NearResearchTable());
        Factory::InsertFuel(c.body->id, plot->id, dynamo, jar, 1);
        return true;
    });
    s.step("the dynamo takes a jar as fuel", 10.0, [](ScenarioCtx& c) {
        c.camera(Vector3Add(nodeGround, Vector3Add(Vector3Scale(east, -4.0f), Vector3Scale(nodeUp, 3.0f))), Vector3Add(nodeGround, Vector3Scale(nodeUp, 0.5f)), nodeUp);
        Factory::Sync(c.world);
        const MachineView v = Factory::View(c.body->id, plot->id, dynamo);
        c.waitingFor(TextFormat("dynamo fuel %s, energy %d", v.input.c_str(), v.energy));
        if (v.input.empty() && v.energy <= 0) return false;
        c.check("a jar of mana is the dynamo's fuel", true);
        Factory::SetSpeed(1);
        return true;
    });
}
