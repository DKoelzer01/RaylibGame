#include "laddermigrate.h"
#include "factory.h"
#include "gameui.h"
#include "playerrecord.h"
#include "plots.h"
#include "research.h"
#include "researchgate.h"
#include "researchtab.h"
#include "survival.h"
#include "teamquests.h"
#include "worldtypes.h"
#include "headless.h"
#include <cstdio>

// Save migration to the new tier ladder (#1272). See laddermigrate.h.

int LadderMigrate_MapTier(int oldTier, bool steelMade, bool steamPlaced) {
    if (oldTier == 2) return steelMade || steamPlaced ? 3 : 2; // (decided: old 2 spans the new 2 and 3)
    return oldTier;                                            // 0, 1 and 3-7 keep their numbers; 8 and up are untouched
}

LadderMigratePlan LadderMigrate_Plan(const LadderMigrateInput& in) {
    LadderMigratePlan p;
    const bool steam = in.placed.count("boiler") > 0 || in.placed.count("forge_hammer") > 0;
    p.newTier = LadderMigrate_MapTier(in.oldTier, in.steelMade, steam);
    for (const ResearchNode& n : Research_Nodes()) {
        if (n.tier < p.newTier) { p.nodes.insert(n.id); continue; } // everything below: core and branch (generous on purpose)
        if (n.tier != p.newTier) continue;
        // Inside the tier: a node whose machine or recipe the team already placed or made.
        for (const std::string& u : n.unlocks) {
            const size_t colon = u.find(':');
            if (colon == std::string::npos) continue;
            const std::string kind = u.substr(0, colon), what = u.substr(colon + 1);
            if ((kind == "machine" && (in.placed.count(what) || (in.made && in.made(what)))) || (kind == "recipe" && in.made && in.made(what))) { p.nodes.insert(n.id); break; }
        }
    }
    for (int t = 3; t < p.newTier && t <= 7; ++t)
        if (const TierMilestone* m = ResearchGate_Milestone(t))
            for (int l = 0; l < (int)m->lines.size(); ++l) p.milestoneLines.insert({ t, l });
    return p;
}

bool LadderMigrate_Apply(QuestBookState& book, ResearchState& research, const LadderMigratePlan& plan, const std::function<int(int)>& taskCount) {
    if (research.ladderMigrated) return false;
    // (#1271 L5) the new tiers' task lists differ from the old ones: from the new tier up, tasks done are judged again
    for (auto it = book.tasksDone.begin(); it != book.tasksDone.end();) it = it->first >= plan.newTier ? book.tasksDone.erase(it) : std::next(it);
    if (plan.newTier > book.tier) {
        for (int t = 0; t < plan.newTier; ++t) {
            const int n = taskCount ? taskCount(t) : 0;
            for (int k = 0; k < n; ++k) book.tasksDone.insert({ t, k });
        }
        book.tier = plan.newTier;
    }
    research.done.insert(plan.nodes.begin(), plan.nodes.end());
    for (const std::string& d : plan.nodes) research.progress.erase(d);
    research.milestoneDone.insert(plan.milestoneLines.begin(), plan.milestoneLines.end());
    research.ladderMigrated = true;
    return true;
}

bool LadderMigrate_IsRemovedMachine(int t) {
    return t == (int)PartType::ForgeHammer || t == (int)PartType::AlloySmelter || t == (int)PartType::SteamMacerator;
}

void LadderMigrate_RefundIntoChest(Planetoid& body, Plot& plot, int part) {
    if (part < 0 || part >= (int)plot.parts.size()) return;
    const char* item = PartCostItem(plot.parts[part].type); // the machine itself goes back as an item
    const std::string machine = item ? item : "";
    if (!plot.retype(part, PartType::Chest)) return;      // a chest has the same 1 m footprint: what it held loads into it
    if (!machine.empty()) Factory::Insert(body.id, plot.id, part, machine, 1);
}

int LadderMigrate_World() {
    int migrated = 0;
    const std::function<int(int)> tasks = [](int t) { return QuestBook_TaskCount(t); };
    auto migrate = [&](int team, QuestBookState& book) {
        ResearchState& rs = Research_Team(team);
        if (rs.ladderMigrated) return;
        QuestBook_Select(&book);
        LadderMigrateInput in;
        in.oldTier = book.tier;
        in.steelMade = QuestBook_Met("get", "steel_ingot", 1.0);
        in.placed = book.placed;
        in.made = [](const std::string& item) { return QuestBook_Met("get", item, 1.0); };
        const LadderMigratePlan plan = LadderMigrate_Plan(in);
        QuestBook_Select(nullptr);
        if (LadderMigrate_Apply(book, rs, plan, tasks)) {
            ++migrated;
            if (plan.newTier != in.oldTier) printf("ladder: team %d moved from tier %d to tier %d (%zu research nodes granted)\n", team, in.oldTier, plan.newTier, plan.nodes.size());
        }
    };
    const int mine = Headless() ? 0 : ResearchTab_Team(); // the local player's team (a team of one to start with; a server has none)
    for (const TeamRecord& t : Teams_All()) {
        if (t.id == mine) migrate(t.id, QuestBook_LocalState());
        else migrate(t.id, TeamQuests_Book(t.id));
    }
    // The 8 V machines that are gone.
    if (!LadderMigrate_RemovedMachineHook) LadderMigrate_RemovedMachineHook = LadderMigrate_RefundIntoChest;
    int removed = 0;
    for (Planetoid* body : Planetoid::All()) {
        if (!body || !body->plots) continue;
        for (auto& pl : body->plots->plots)
            for (int i = 0; i < (int)pl->parts.size(); ++i)
                if (LadderMigrate_IsRemovedMachine((int)pl->parts[i].type)) { LadderMigrate_RemovedMachineHook(*body, *pl, i); ++removed; }
    }
    // (#1427 L7a) A machine's material is its variant step now: an older save's machines carry the build tool's shape material
    // (it leaked in when placed). Cleared once, on the load that migrates the book (a ladder save is never touched).
    if (migrated > 0)
        for (Planetoid* body : Planetoid::All()) {
            if (!body || !body->plots) continue;
            for (auto& pl : body->plots->plots) {
                bool changed = false;
                for (Part& p : pl->parts) if (IsMachine(p.type) && p.material != 0) { p.material = 0; changed = true; }
                if (changed) ++pl->version; // (the factory takes the layout again)
            }
        }
    if (removed > 0)
        GameUI_Toast(TextFormat("%d old steam machine%s became chests holding them: the press, millstone and steam engine do their work now", removed, removed == 1 ? "" : "s"));
    return migrated;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1272): the tier mapping, what is granted, and that it runs once ----
#include <gtest/gtest.h>

TEST(LadderMigrate, OldTiersMapToTheNewLadder) {
    for (int t : { 0, 1, 3, 4, 5, 6, 7, 8, 9, 10 }) EXPECT_EQ(LadderMigrate_MapTier(t, false, false), t) << "tier " << t;
    EXPECT_EQ(LadderMigrate_MapTier(2, false, false), 2) << "old Bronze and steam, no steel, no boiler: new 2";
    EXPECT_EQ(LadderMigrate_MapTier(2, true, false), 3) << "steel made: new 3";
    EXPECT_EQ(LadderMigrate_MapTier(2, false, true), 3) << "a boiler or forge hammer placed: new 2's milestone granted";
}

TEST(LadderMigrate, EverythingBelowIsGrantedAndOwnedMachinesInsideTheTier) {
    ASSERT_TRUE(Research_Load());
    ASSERT_TRUE(ResearchGate_LoadMilestones());
    LadderMigrateInput in;
    in.oldTier = 5;
    in.placed = { "macerator", "ore_washer" };
    in.made = [](const std::string& item) { return item == "aluminium_ingot"; };
    const LadderMigratePlan p = LadderMigrate_Plan(in);
    EXPECT_EQ(p.newTier, 5);
    for (const ResearchNode& n : Research_Nodes())
        if (n.tier < 5) EXPECT_TRUE(p.nodes.count(n.id)) << n.id << " (tier " << n.tier << ", " << (n.core ? "core" : "branch") << ")";
    bool aTier5Node = false, notOwned = false;
    for (const ResearchNode& n : Research_Nodes()) {
        if (n.tier != 5) continue;
        if (p.nodes.count(n.id)) aTier5Node = true; else notOwned = true;
    }
    EXPECT_TRUE(aTier5Node) << "a tier 5 node whose machine or recipe the team has";
    EXPECT_TRUE(notOwned) << "but not the whole tier";
    for (int t = 3; t < 5; ++t) {
        const TierMilestone* m = ResearchGate_Milestone(t);
        ASSERT_NE(m, nullptr);
        for (int l = 0; l < (int)m->lines.size(); ++l) EXPECT_TRUE(p.milestoneLines.count({ t, l })) << "tier " << t << " line " << l;
    }
    EXPECT_FALSE(p.milestoneLines.count({ 5, 0 })) << "the current tier's milestone is still to do";
}

TEST(LadderMigrate, ItRunsOnceAndMovesTheBookUp) {
    ASSERT_TRUE(Research_Load());
    QuestBookState book;
    book.tier = 2;
    book.tasksDone = { { 2, 0 }, { 3, 1 } }; // (3, 1): an old tier 3 task index, meaningless in the new tier 3
    ResearchState rs;
    LadderMigrateInput in;
    in.oldTier = 2;
    in.steelMade = true;
    const LadderMigratePlan p = LadderMigrate_Plan(in);
    ASSERT_EQ(p.newTier, 3);
    const auto tasks = [](int t) { return t + 2; }; // (any counts: tier t has t + 2 tasks here)
    EXPECT_TRUE(LadderMigrate_Apply(book, rs, p, tasks));
    EXPECT_EQ(book.tier, 3);
    for (int t = 0; t < 3; ++t) for (int k = 0; k < t + 2; ++k) EXPECT_TRUE(book.tasksDone.count({ t, k })) << t << "/" << k;
    EXPECT_FALSE(book.tasksDone.count({ 3, 1 })) << "the new tier's tasks are judged again";
    EXPECT_TRUE(rs.ladderMigrated);
    // A second time (a second load) changes nothing.
    const QuestBookState before = book;
    const std::string rsBefore = rs.Text();
    LadderMigratePlan more = p; more.newTier = 6; more.nodes.insert("arc_furnace");
    EXPECT_FALSE(LadderMigrate_Apply(book, rs, more, tasks));
    EXPECT_EQ(book.tier, before.tier);
    EXPECT_EQ(book.tasksDone, before.tasksDone);
    EXPECT_EQ(rs.Text(), rsBefore);
    // The mark is saved with the team's research.
    ResearchState back;
    ASSERT_TRUE(back.FromText(rs.Text()));
    EXPECT_TRUE(back.ladderMigrated);
}

TEST(LadderMigrate, OnlyTheEightVoltMachinesAreRemoved) {
    EXPECT_TRUE(LadderMigrate_IsRemovedMachine((int)PartType::ForgeHammer));
    EXPECT_TRUE(LadderMigrate_IsRemovedMachine((int)PartType::AlloySmelter));
    EXPECT_TRUE(LadderMigrate_IsRemovedMachine((int)PartType::SteamMacerator));
    EXPECT_FALSE(LadderMigrate_IsRemovedMachine((int)PartType::Boiler));
    EXPECT_FALSE(LadderMigrate_IsRemovedMachine((int)PartType::Macerator));
}

TEST(LadderMigrate, ARemovedMachineBecomesAChestInPlace) {
    Plot plot;
    plot.east = { 1, 0, 0 }; plot.up = { 0, 1, 0 }; plot.north = { 0, 0, 1 };
    const int hammer = plot.place({ PartType::ForgeHammer, { 0, 0, 0 }, 0 });
    ASSERT_GE(hammer, 0);
    const uint32_t v = plot.version;
    ASSERT_TRUE(plot.retype(hammer, PartType::Chest));
    EXPECT_EQ(plot.parts[hammer].type, PartType::Chest) << "the same index: the saved contents load into it";
    EXPECT_EQ(plot.partAt({ 1, 1, 1 }), hammer);
    EXPECT_GT(plot.version, v) << "the factory resyncs";
    EXPECT_FALSE(plot.retype(hammer, PartType::Shaft)) << "another footprint: refused";
}
#endif
