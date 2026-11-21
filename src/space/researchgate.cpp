#include "researchgate.h"
#include "researchtab.h"
#include "factory.h"
#include "net.h"
#include "simclock.h"
#include "terraintool.h"
#include "inventory.h"
#include "gameui.h"
#include "survival.h"
#include <algorithm>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>

// Research gating and tier milestones (#1271 L1). See researchgate.h.

namespace {
int gateOn = -1; // -1: not read yet (RAYLIBGAME_LADDER)
int selectedTeam = 0;

// recipe output / machine item -> the node that unlocks it, rebuilt when the node table changes
struct Index { const void* data = nullptr; size_t size = 0; std::unordered_map<std::string, const ResearchNode*> recipe, machine; } index;
const Index& Idx() {
    const std::vector<ResearchNode>& nodes = Research_Nodes();
    if (index.data != (const void*)nodes.data() || index.size != nodes.size()) {
        index = Index{};
        index.data = nodes.data(); index.size = nodes.size();
        for (const ResearchNode& n : nodes)
            for (const std::string& u : n.unlocks) {
                if (u.rfind("recipe:", 0) == 0) index.recipe.emplace(u.substr(7), &n);
                else if (u.rfind("machine:", 0) == 0) index.machine.emplace(u.substr(8), &n);
            }
    }
    return index;
}

std::map<int, TierMilestone> milestones;

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
const char* const kFinds[] = { "coal_seam", "crystal_sample", "vent_sample", "seafloor_survey", "moon_survey", "sounding_rocket" };
bool KnownFind(const std::string& f) { for (const char* k : kFinds) if (f == k) return true; return false; }
} // namespace

bool ResearchGate_On() {
    if (gateOn < 0) { const char* e = getenv("RAYLIBGAME_LADDER"); gateOn = e && e[0] == '0' ? 0 : 1; } // (#1271 L6b) on; RAYLIBGAME_LADDER=0: off
    return gateOn == 1;
}
void ResearchGate_SetOn(bool on) { gateOn = on ? 1 : 0; }

const ResearchNode* ResearchGate_NodeForRecipe(const std::string& output) {
    const Index& ix = Idx();
    auto it = ix.recipe.find(output);
    return it == ix.recipe.end() ? nullptr : it->second;
}

const ResearchNode* ResearchGate_NodeForPart(PartType t) {
    const Index& ix = Idx();
    for (const char* name : { PartCostItem(t), PartMachineName(t) }) {
        if (!name) continue;
        auto it = ix.machine.find(name);
        if (it != ix.machine.end()) return it->second;
    }
    return nullptr;
}

bool ResearchGate_RecipeAllowed(const ResearchState& s, const std::string& output) {
    if (!ResearchGate_On()) return true;
    const ResearchNode* n = ResearchGate_NodeForRecipe(output);
    if (!n) n = Idx().machine.count(output) ? Idx().machine.at(output) : nullptr; // crafting a machine's item: its machine node
    return !n || s.Done(n->id);
}

bool ResearchGate_PartAllowed(const ResearchState& s, PartType t) {
    if (!ResearchGate_On()) return true;
    const ResearchNode* n = ResearchGate_NodeForPart(t);
    return !n || s.Done(n->id);
}

void ResearchGate_SelectTeam(int team) { selectedTeam = team; }

const ResearchNode* ResearchGate_NodeForMachineItem(const std::string& item) {
    const Index& ix = Idx();
    auto it = ix.machine.find(item);
    return it == ix.machine.end() ? nullptr : it->second;
}

void CheckThroughput(); // (#1271 L6b-2, below)
void ResearchGate_Update() {
    if (gNetRole == NetRole::Client) return; // (the server's factory runs the machines)
    static double next = 0.0;
    static std::vector<std::string> sent;
    static bool sentOnce = false;
    const double now = SimTime();
    if (next - now > 1.0) next = now;
    if (now < next) return;
    next = now + 0.5;
    static double nextThroughput = 0.0;
    if (now >= nextThroughput || nextThroughput - now > 20.0) { nextThroughput = now + 10.0; if (ResearchGate_On()) CheckThroughput(); } // (#1271 L6b-2)
    std::vector<std::string> locked;
    if (ResearchGate_On() && !terrainTool.creative) { // (creative: every machine recipe)
        const ResearchState& team = ResearchGate_Team();
        for (const Factory::RecipeInfo& r : Factory::RecipeList())
            for (const auto& o : r.out)
                if (!ResearchGate_RecipeAllowed(team, o.first)) locked.push_back(o.first);
        std::sort(locked.begin(), locked.end());
        locked.erase(std::unique(locked.begin(), locked.end()), locked.end());
    }
    if (sentOnce && locked == sent) return;
    sentOnce = true;
    sent = locked;
    Factory::SetLockedOutputs(locked);
}

std::string ResearchGate_MachineLockNote(PartType t) {
    if (terrainTool.creative) return "";
    if (!ResearchGate_On()) return "";
    const ResearchState& team = ResearchGate_Team();
    int n = 0;
    const ResearchNode* first = nullptr;
    for (const Factory::RecipeInfo& r : Factory::RecipeList()) {
        if (r.machine != t || r.out.empty() || ResearchGate_RecipeAllowed(team, r.out[0].first)) continue;
        ++n;
        if (!first) first = ResearchGate_NodeForRecipe(r.out[0].first);
    }
    if (!n) return "";
    return std::to_string(n) + (n == 1 ? " recipe" : " recipes") + " locked: research " + (first ? first->name + " (tier " + std::to_string(first->tier) + ")" : std::string("it"));
}
const ResearchState& ResearchGate_Team() { return Research_Team(selectedTeam ? selectedTeam : ResearchTab_Team()); }

bool ResearchGate_LoadMilestonesText(const std::string& text, std::string& error) {
    std::map<int, TierMilestone> out;
    std::istringstream in(text);
    int lineNo = 0, tier = -1;
    for (std::string raw; std::getline(in, raw);) {
        ++lineNo;
        const std::string line = Trim(raw);
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> cols;
        { std::istringstream cs(line); for (std::string c; std::getline(cs, c, '|');) cols.push_back(Trim(c)); }
        if (cols[0] == "milestone") {
            if (tier < 0 || cols.size() < 2) { error = "tiers.txt line " + std::to_string(lineNo) + ": a milestone line outside a tier"; return false; }
            MilestoneLine ml;
            ml.text = cols.size() > 2 ? cols[2] : "";
            // Options separated by " or "; the first names the kind, later ones may leave it out ("throughput lv_pack 9 10").
            std::string spec = cols[1], kind;
            size_t at = 0;
            while (at <= spec.size()) {
                const size_t orAt = spec.find(" or ", at);
                std::istringstream os(spec.substr(at, orAt == std::string::npos ? std::string::npos : orAt - at));
                std::vector<std::string> w;
                for (std::string x; os >> x;) w.push_back(x);
                if (!w.empty() && (w[0] == "deliver" || w[0] == "build" || w[0] == "find")) { kind = w[0]; w.erase(w.begin()); }
                if (kind.empty() || w.empty()) { error = "tiers.txt line " + std::to_string(lineNo) + ": milestone '" + cols[1] + "': deliver, build or find what?"; return false; }
                MilestoneOption o;
                o.kind = kind; o.what = w[0];
                if (kind == "deliver") o.count = w.size() > 1 ? atoi(w[1].c_str()) : 0;
                if (kind == "build" && o.what == "throughput") {
                    o.item = w.size() > 1 ? w[1] : "";
                    o.count = w.size() > 2 ? atoi(w[2].c_str()) : 0;
                    o.minutes = w.size() > 3 ? atoi(w[3].c_str()) : 0;
                }
                ml.options.push_back(o);
                if (orAt == std::string::npos) break;
                at = orAt + 4;
            }
            TierMilestone& tm = out[tier];
            tm.tier = tier;
            tm.lines.push_back(ml);
        } else if (cols[0] != "task" && cols.size() >= 4) {
            tier = atoi(cols[0].c_str()); // a tier's header line: "3 | Steel | ..."
        }
    }
    milestones = std::move(out);
    return true;
}

void ResearchGate_InstallTierGate() { // (#1271 L6b-2) tiers 3-7 close on their research milestone as well as their tasks (the gate on)
    gQuestTierGate = [](int tier) { return !ResearchGate_On() || tier < 3 || tier > 7 || ResearchGate_MilestoneMet(ResearchGate_Team(), tier); };
}

bool ResearchGate_LoadMilestones(const char* path) {
    std::ifstream f(path);
    if (!f) return false;
    std::stringstream ss; ss << f.rdbuf();
    std::string error;
    if (!ResearchGate_LoadMilestonesText(ss.str(), error)) { printf("research: %s\n", error.c_str()); return false; }
    const std::string problems = ResearchGate_ValidateMilestones();
    if (!problems.empty()) printf("research: %s", problems.c_str());
    return true;
}

const TierMilestone* ResearchGate_Milestone(int tier) {
    auto it = milestones.find(tier);
    return it == milestones.end() ? nullptr : &it->second;
}

std::string ResearchGate_ValidateMilestones() {
    std::string out;
    for (const auto& [tier, m] : milestones)
        if (tier < 3 || tier > 7) out += "tier " + std::to_string(tier) + ": a milestone (only tiers 3-7 have them; tiers 0-2 end on tasks)\n";
    for (int tier = 3; tier <= 7; ++tier) {
        const TierMilestone* m = ResearchGate_Milestone(tier);
        if (!m) { out += "tier " + std::to_string(tier) + ": no milestone\n"; continue; }
        bool deliver = false, build = false, find = false;
        for (const MilestoneLine& l : m->lines)
            for (const MilestoneOption& o : l.options) {
                const std::string where = "tier " + std::to_string(tier) + " milestone '" + o.kind + " " + o.what + "': ";
                if (o.kind == "deliver") { deliver = true; if (o.count <= 0) out += where + "deliver how many?\n"; }
                else if (o.kind == "build") {
                    build = true;
                    if (o.what == "throughput") { if (o.item.empty() || o.count <= 0 || o.minutes <= 0) out += where + "throughput needs <item> <per minute> <minutes>\n"; }
                    else if (o.what != "lab_node" && o.what != "rocket_on_pad") out += where + "no such build (lab_node, throughput, rocket_on_pad)\n";
                } else if (o.kind == "find") { find = true; if (!KnownFind(o.what)) out += where + "no such find\n"; }
            }
        if (!deliver) out += "tier " + std::to_string(tier) + ": the milestone has no delivery\n";
        if (!build) out += "tier " + std::to_string(tier) + ": the milestone has no working build\n";
        if (!find) out += "tier " + std::to_string(tier) + ": the milestone has no find\n";
    }
    return out;
}

bool ResearchGate_CoreDone(const ResearchState& s, int tier) {
    for (const ResearchNode& n : Research_Nodes()) if (n.tier == tier && n.core && !s.Done(n.id)) return false;
    return true;
}

bool ResearchGate_MilestoneMet(const ResearchState& s, int tier) {
    const TierMilestone* m = ResearchGate_Milestone(tier);
    if (!m || !ResearchGate_CoreDone(s, tier)) return false;
    for (int l = 0; l < (int)m->lines.size(); ++l) if (!s.milestoneDone.count({ tier, l })) return false;
    return true;
}

// ---- (#1271 L6b-2) Working on the milestone: hand-ins at a lab, builds, throughput ----
namespace {
ResearchState& TeamMut() { return Research_Team(selectedTeam ? selectedTeam : ResearchTab_Team()); }
void Mark(ResearchState& s, int tier, int line, const std::string& what) {
    if (!s.milestoneDone.insert({ tier, line }).second) return;
    GameUI_Toast(TextFormat("Tier %d milestone: %s", tier, what.c_str()));
}
}

std::string ResearchGate_LineLabel(const MilestoneLine& l) {
    std::string s;
    for (const MilestoneOption& o : l.options) {
        if (!s.empty()) s += " or ";
        if (o.kind == "deliver") s += "deliver " + std::to_string(o.count) + " " + o.what;
        else if (o.kind == "build" && o.what == "throughput") s += "make " + std::to_string(o.count) + " " + o.item + " a minute for " + std::to_string(o.minutes) + " minutes";
        else s += o.kind + " " + o.what;
    }
    return s;
}

int ResearchGate_HandIn() {
    const int tier = QuestBook_Tier();
    const TierMilestone* m = tier >= 3 && tier <= 7 ? ResearchGate_Milestone(tier) : nullptr;
    if (!m) return 0;
    ResearchState& s = TeamMut();
    int done = 0;
    for (int l = 0; l < (int)m->lines.size(); ++l) {
        if (s.milestoneDone.count({ tier, l })) continue;
        const MilestoneLine& line = m->lines[l];
        for (const MilestoneOption& o : line.options) {
            if (o.kind == "deliver" && o.count > 0 && inventory.count(o.what) >= o.count) {
                inventory.remove(o.what, o.count);
                Mark(s, tier, l, line.text.empty() ? ResearchGate_LineLabel(line) : line.text); ++done; break;
            }
            if (o.kind == "find") { // its sample, or (the finds step is to come) a ruin's schematic standing in
                const std::string sample = o.what + "_sample";
                const char* use = inventory.count(sample) > 0 ? sample.c_str() : inventory.count("ruin_schematic") > 0 ? "ruin_schematic" : nullptr;
                if (use) { inventory.remove(use, 1); Mark(s, tier, l, line.text.empty() ? ResearchGate_LineLabel(line) : line.text); ++done; break; }
            }
        }
    }
    return done;
}

void ResearchGate_MarkBuild(const std::string& what) {
    ResearchState& s = TeamMut();
    for (int tier = 3; tier <= 7; ++tier)
        if (const TierMilestone* m = ResearchGate_Milestone(tier))
            for (int l = 0; l < (int)m->lines.size(); ++l)
                for (const MilestoneOption& o : m->lines[l].options)
                    if (o.kind == "build" && o.what == what) Mark(s, tier, l, m->lines[l].text.empty() ? ResearchGate_LineLabel(m->lines[l]) : m->lines[l].text);
}

bool ResearchGate_ThroughputMet(const std::vector<float>& b, int perMinute, int minutes) {
    if (perMinute <= 0 || minutes <= 0 || (int)b.size() < minutes) return false;
    double sum = 0.0;
    for (int i = (int)b.size() - minutes; i < (int)b.size(); ++i) {
        if (b[i] < 0.8f * perMinute) return false; // no minute under 80% (decided)
        sum += b[i];
    }
    return sum / minutes >= perMinute;
}

// Every 10 s: the current tier's throughput lines, from the statistics' per-minute buckets (machine crafts only).
void CheckThroughput() {
    const int tier = QuestBook_Tier();
    const TierMilestone* m = tier >= 3 && tier <= 7 ? ResearchGate_Milestone(tier) : nullptr;
    if (!m) return;
    ResearchState& s = TeamMut();
    bool any = false;
    for (int l = 0; l < (int)m->lines.size() && !any; ++l)
        if (!s.milestoneDone.count({ tier, l })) for (const MilestoneOption& o : m->lines[l].options) any = any || (o.kind == "build" && o.what == "throughput");
    if (!any) return;
    const FactoryStatsView v = Factory::StatsView(StatWindow::Hour);
    for (int l = 0; l < (int)m->lines.size(); ++l) {
        if (s.milestoneDone.count({ tier, l })) continue;
        for (const MilestoneOption& o : m->lines[l].options) {
            if (o.kind != "build" || o.what != "throughput") continue;
            for (const StatLine& sl : v.items)
                if (sl.item == o.item && v.coveredMinutes >= o.minutes && ResearchGate_ThroughputMet(sl.madeSpark, o.count, o.minutes)) {
                    Mark(s, tier, l, m->lines[l].text.empty() ? ResearchGate_LineLabel(m->lines[l]) : m->lines[l].text);
                    break;
                }
        }
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1271 L1): the gate (off by default), what it gates, and the tier milestones ----
#include <gtest/gtest.h>

namespace {
struct GateOn { bool was; GateOn() : was(ResearchGate_On()) { ResearchGate_SetOn(true); } ~GateOn() { ResearchGate_SetOn(was); } };
}

TEST(ResearchGate, OffItGatesNothing) {
    ASSERT_TRUE(Research_Load());
    const bool was = ResearchGate_On();
    ResearchGate_SetOn(false);
    ResearchState none;
    EXPECT_TRUE(ResearchGate_RecipeAllowed(none, "steel_ingot"));
    EXPECT_TRUE(ResearchGate_PartAllowed(none, PartType::BlastFurnace));
    ResearchGate_SetOn(was);
}

TEST(ResearchGate, OnItGatesWhatANodeUnlocksUntilTheNodeIsDone) {
    ASSERT_TRUE(Research_Load());
    GateOn on;
    ASSERT_NE(ResearchGate_NodeForRecipe("steel_ingot"), nullptr);
    EXPECT_EQ(ResearchGate_NodeForRecipe("steel_ingot")->id, "steelmaking");
    ASSERT_NE(ResearchGate_NodeForPart(PartType::BlastFurnace), nullptr) << "machine:blast_furnace";
    ASSERT_NE(ResearchGate_NodeForPart(PartType::Furnace), nullptr) << "machine:electric_furnace: the furnace part's item";
    ResearchState s;
    EXPECT_FALSE(ResearchGate_RecipeAllowed(s, "steel_ingot"));
    EXPECT_FALSE(ResearchGate_PartAllowed(s, PartType::BlastFurnace));
    EXPECT_FALSE(ResearchGate_RecipeAllowed(s, "blast_furnace")) << "crafting the machine's item: its node too";
    EXPECT_TRUE(ResearchGate_RecipeAllowed(s, "stone_pickaxe")) << "no node names it: not gated here";
    EXPECT_TRUE(ResearchGate_PartAllowed(s, PartType::Kiln));
    s.done.insert("steelmaking");
    EXPECT_TRUE(ResearchGate_RecipeAllowed(s, "steel_ingot"));
    EXPECT_TRUE(ResearchGate_PartAllowed(s, PartType::BlastFurnace));
}

TEST(ResearchGate, TheShippedMilestonesLoadAndValidate) {
    ASSERT_TRUE(ResearchGate_LoadMilestones()) << "resources/data/tiers.txt (tests run from the repo root)";
    EXPECT_EQ(ResearchGate_ValidateMilestones(), "");
    for (int t = 0; t <= 2; ++t) EXPECT_EQ(ResearchGate_Milestone(t), nullptr) << "tiers 0-2 end on tasks";
    const TierMilestone* t4 = ResearchGate_Milestone(4);
    ASSERT_NE(t4, nullptr);
    bool alternatives = false;
    for (const MilestoneLine& l : t4->lines) if (l.options.size() == 2 && l.options[1].what == "throughput") alternatives = true;
    EXPECT_TRUE(alternatives) << "T4's throughput target: circuits or LV packs";
}

TEST(ResearchGate, ParsingAndValidationCatchBadMilestones) {
    std::string error;
    ASSERT_TRUE(ResearchGate_LoadMilestonesText(
        "3 | Steam | - | text\n"
        "    task | get coal 1 | a task\n"
        "    milestone | deliver steel_plate 50 | fifty plates\n"
        "    milestone | build throughput lv_pack 9 10 or throughput basic_circuit 15 10 | a line\n"
        "    milestone | find coal_seam | a find\n"
        "4 | LV | - | text\n"
        "    milestone | deliver basic_circuit 0 | none\n"
        "    milestone | build windmill_farm | ?\n"
        "    milestone | find dragon_egg | ?\n"
        "1 | Copper | - | text\n"
        "    milestone | find coal_seam | early\n", error)) << error;
    const TierMilestone* t3 = ResearchGate_Milestone(3);
    ASSERT_NE(t3, nullptr);
    ASSERT_EQ(t3->lines.size(), 3u);
    ASSERT_EQ(t3->lines[1].options.size(), 2u);
    EXPECT_EQ(t3->lines[1].options[1].item, "basic_circuit");
    EXPECT_EQ(t3->lines[1].options[1].count, 15);
    EXPECT_EQ(t3->lines[1].options[1].minutes, 10);
    const std::string v = ResearchGate_ValidateMilestones();
    EXPECT_NE(v.find("deliver how many"), std::string::npos) << v;
    EXPECT_NE(v.find("no such build"), std::string::npos);
    EXPECT_NE(v.find("no such find"), std::string::npos);
    EXPECT_NE(v.find("tier 1: a milestone"), std::string::npos);
    EXPECT_NE(v.find("tier 5: no milestone"), std::string::npos);
    EXPECT_FALSE(ResearchGate_LoadMilestonesText("    milestone | find coal_seam | before any tier\n", error));
    ResearchGate_LoadMilestones(); // (the shipped table back)
}

TEST(ResearchGate, AMilestoneNeedsTheCoreResearchAndEveryLine) {
    ASSERT_TRUE(Research_Load());
    ASSERT_TRUE(ResearchGate_LoadMilestones());
    const TierMilestone* m = ResearchGate_Milestone(3);
    ASSERT_NE(m, nullptr);
    ResearchState s;
    for (int l = 0; l < (int)m->lines.size(); ++l) s.milestoneDone.insert({ 3, l });
    EXPECT_FALSE(ResearchGate_MilestoneMet(s, 3)) << "every line, but no core research";
    for (const ResearchNode& n : Research_Nodes()) if (n.tier == 3 && n.core) s.done.insert(n.id);
    EXPECT_TRUE(ResearchGate_CoreDone(s, 3));
    EXPECT_TRUE(ResearchGate_MilestoneMet(s, 3));
    s.milestoneDone.erase({ 3, 0 });
    EXPECT_FALSE(ResearchGate_MilestoneMet(s, 3)) << "a delivery still owed";
    // The lines done are kept with the team's research (save lines) and merged when a player joins.
    s.milestoneDone.insert({ 3, 0 });
    ResearchState back;
    ASSERT_TRUE(back.FromText(s.Text()));
    EXPECT_EQ(back.milestoneDone, s.milestoneDone);
    ResearchState team;
    Research_Merge(s, team);
    EXPECT_EQ(team.milestoneDone, s.milestoneDone);
}

// (#1271 L6b-2) A throughput target: a 10-minute mean of at least N a minute, no minute under 80% of N (decided).
TEST(ResearchGate, ThroughputNeedsTheMeanAndNoWeakMinute) {
    std::vector<float> steady(12, 15.0f);
    EXPECT_TRUE(ResearchGate_ThroughputMet(steady, 15, 10));
    EXPECT_FALSE(ResearchGate_ThroughputMet(steady, 16, 10)) << "the mean under the target";
    std::vector<float> dip = steady;
    dip[8] = 11.0f; // under 0.8 x 15 = 12
    EXPECT_FALSE(ResearchGate_ThroughputMet(dip, 15, 10)) << "one weak minute";
    std::vector<float> burst = steady;
    burst[0] = 0.0f; // outside the last 10 minutes: it doesn't count
    EXPECT_TRUE(ResearchGate_ThroughputMet(burst, 15, 10));
    EXPECT_FALSE(ResearchGate_ThroughputMet(std::vector<float>(5, 30.0f), 15, 10)) << "not 10 minutes of data yet";
    MilestoneLine l;
    l.options = { { "build", "throughput", "basic_circuit", 15, 10 }, { "build", "throughput", "lv_pack", 9, 10 } };
    EXPECT_EQ(ResearchGate_LineLabel(l), "make 15 basic_circuit a minute for 10 minutes or make 9 lv_pack a minute for 10 minutes");
}
#endif
