// saves/migrate (#1272, #1271 L5): an old save's state moves to the new tier ladder once. In the sandbox the world is made to
// look like one from before the ladder: the quest book at the old tier 2 with a boiler placed (so the new tier 2's milestone is
// granted: new tier 3), a forge hammer (a removed 8 V machine) standing on a plot, and research never
// migrated. The migration then: tier 3, the steam node granted (the boiler it already has), the forge hammer turned into a chest
// in its place holding the forge hammer itself, a toast (what a saved machine held loads into the chest: the part keeps its index, unit-tested). A second run changes nothing.
#include "scenario.h"
#include "buildtool.h"
#include "factory.h"
#include "inventory.h"
#include "laddermigrate.h"
#include "plots.h"
#include "research.h"
#include "researchtab.h"
#include "survival.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Plot* plot = nullptr;
int hammer = -1;
bool ChestHas(ScenarioCtx& c, const char* item) {
    for (const auto& [id, n] : Factory::ChestContents(c.body->id, plot->id, hammer)) if (id == item && n > 0) return true;
    return false;
}
}

SCENARIO("saves/migrate") {
    s.tags = { "factory", "saves" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    s.step("a world from before the ladder: old tier 2, a boiler placed, a forge hammer", 10.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        const Vector3 d = Vector3Normalize(c.spawnDirLocal());
        const Vector3 ground = Vector3Scale(d, c.field->radius + c.field->elevation(d) + 0.05f);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        if (!c.check("a forge hammer placed", BuildTool_Place(*c.body, *plot, { PartType::ForgeHammer, { 0, 0, 0 }, 0 }, false))) return true;
        hammer = plot->partAt({ 0, 0, 0 });
        c.camera(plot->toBody({ 3.0f, 2.5f, -3.0f }), plot->toBody({ 0.5f, 0.5f, 0.5f }), plot->up);
        Factory::Sync(c.world);
        Factory::WaitTicks(2);
        QuestBook_SetTier(2);
        QuestBook_NotePlaced("boiler");
        ResearchState& rs = Research_Team(ResearchTab_Team());
        rs = ResearchState{}; // never migrated
        Factory::WaitTicks(2);
        return true;
    });
    s.step("the migration: tier 3, the boiler's node, the forge hammer a chest", 10.0, [](ScenarioCtx& c) {
        static bool ran = false;
        if (!ran) {
            ran = true;
            c.eq("one team migrated", LadderMigrate_World(), 1);
            Factory::Sync(c.world);
            return false;
        }
        Factory::Sync(c.world);
        Factory::WaitTicks(2);
        Factory::ChestContents(c.body->id, plot->id, hammer); // (watched: filled from the next tick)
        if (c.t < 1.0) return false;
        c.eq("old tier 2 with a boiler placed: new tier 3", QuestBook_Tier(), 3);
        const ResearchState& rs = Research_Team(ResearchTab_Team());
        c.check("the boiler's node (steam power) granted", rs.Done("steam_power"));
        c.check("the team is marked migrated", rs.ladderMigrated);
        c.check("the forge hammer became a chest in its place", plot->parts[hammer].type == PartType::Chest);
        c.check("the chest holds the forge hammer", ChestHas(c, "forge_hammer"));
        return true;
    });
    s.step("a second load migrates nothing", 5.0, [](ScenarioCtx& c) {
        const int tier = QuestBook_Tier();
        const std::string before = Research_Team(ResearchTab_Team()).Text();
        c.eq("no team migrated again", LadderMigrate_World(), 0);
        c.eq("the tier stays", QuestBook_Tier(), tier);
        c.check("the research stays", Research_Team(ResearchTab_Team()).Text() == before);
        c.check("the chest stays a chest", plot->parts[hammer].type == PartType::Chest);
        return true;
    });
}
