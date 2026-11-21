#include "research.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <functional>
#include <sstream>

namespace {
std::vector<ResearchNode> nodes;
std::map<int, ResearchState> teams;

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
std::vector<std::string> Words(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    for (std::string w; in >> w;) if (w != "-") out.push_back(w);
    return out;
}
}

const char* Research_PackOfTier(int tier) {
    static const char* kPacks[] = { "steam_pack", "lv_pack", "mv_pack", "hv_pack", "rocket_pack" };
    return tier >= 3 && tier <= 7 ? kPacks[tier - 3] : "";
}

bool Research_LoadText(const std::string& text, std::string& error) {
    std::vector<ResearchNode> out;
    std::istringstream in(text);
    int lineNo = 0;
    for (std::string line; std::getline(in, line);) {
        ++lineNo;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        std::vector<std::string> fields;
        { std::istringstream fs(line); for (std::string f; std::getline(fs, f, '|');) fields.push_back(Trim(f)); }
        if (fields.size() < 2 || fields[0].rfind("node ", 0) != 0) { error = "research.txt line " + std::to_string(lineNo) + ": not a node line"; return false; }
        ResearchNode n;
        n.id = Trim(fields[0].substr(5));
        n.name = fields[1];
        for (size_t k = 2; k < fields.size(); ++k) {
            const std::string& f = fields[k];
            std::istringstream fs(f);
            std::string key; fs >> key;
            std::string rest; std::getline(fs, rest); rest = Trim(rest);
            if (key == "tier") n.tier = atoi(rest.c_str());
            else if (key == "core") n.core = true;
            else if (key == "branch") { n.core = false; n.line = rest; }
            else if (key == "units") n.units = std::max(1, atoi(rest.c_str()));
            else if (key == "time") n.unitSeconds = (float)atof(rest.c_str());
            else if (key == "packs") n.packs = Words(rest);
            else if (key == "needs") n.needs = Words(rest);
            else if (key == "unlocks") n.unlocks = Words(rest);
            else if (key == "tome") n.tome = rest;
            else { error = "research.txt line " + std::to_string(lineNo) + ": unknown field '" + key + "'"; return false; }
        }
        if (n.id.empty()) { error = "research.txt line " + std::to_string(lineNo) + ": a node without an id"; return false; }
        out.push_back(std::move(n));
    }
    nodes = std::move(out);
    return true;
}

bool Research_Load(const char* path) {
    std::ifstream f(path);
    if (!f) return false;
    std::stringstream ss; ss << f.rdbuf();
    std::string error;
    if (!Research_LoadText(ss.str(), error)) { printf("research: %s\n", error.c_str()); return false; }
    const std::string problems = Research_Validate();
    if (!problems.empty()) printf("research: %s", problems.c_str());
    return true;
}

const std::vector<ResearchNode>& Research_Nodes() { return nodes; }
const ResearchNode* Research_Find(const std::string& id) {
    for (const ResearchNode& n : nodes) if (n.id == id) return &n;
    return nullptr;
}

std::string Research_Validate() {
    std::string out;
    std::set<std::string> ids;
    for (const ResearchNode& n : nodes) {
        if (!ids.insert(n.id).second) out += "node " + n.id + ": the id is used twice\n";
        if (n.tier < 3 || n.tier > 7) out += "node " + n.id + ": tier " + std::to_string(n.tier) + " (research is tiers 3-7)\n";
        if (n.core && !n.tome.empty()) out += "node " + n.id + ": a core node needs a tome page\n";
        bool ownPack = false;
        for (const std::string& p : n.packs) {
            int pt = 0;
            for (int t = 3; t <= 7; ++t) if (p == Research_PackOfTier(t)) pt = t;
            if (!pt) out += "node " + n.id + ": '" + p + "' is no research pack\n";
            else if (pt > n.tier) out += "node " + n.id + ": needs " + p + ", a later tier's pack\n";
            ownPack = ownPack || pt == n.tier;
        }
        if (!ownPack) out += "node " + n.id + ": its tier's pack (" + Research_PackOfTier(n.tier) + ") is not among its packs\n";
        for (const std::string& d : n.needs) {
            const ResearchNode* p = Research_Find(d);
            if (!p) out += "node " + n.id + ": needs '" + d + "', no such node\n";
            else if (p->tier > n.tier) out += "node " + n.id + ": needs " + d + " of a later tier\n";
        }
    }
    // Cycles: a depth-first walk over the prerequisites.
    std::map<std::string, int> mark; // 1 on the path, 2 finished
    std::function<bool(const ResearchNode&)> visit = [&](const ResearchNode& n) {
        int& m = mark[n.id];
        if (m == 2) return true;
        if (m == 1) return false;
        m = 1;
        for (const std::string& d : n.needs) if (const ResearchNode* p = Research_Find(d)) if (!visit(*p)) { out += "node " + n.id + ": a prerequisite cycle through " + d + "\n"; return false; }
        m = 2;
        return true;
    };
    for (const ResearchNode& n : nodes) visit(n);
    return out;
}

std::string ResearchState::Text() const {
    std::string s;
    if (!done.empty()) { s += "done"; for (const std::string& d : done) s += " " + d; s += "\n"; }
    if (!queue.empty()) { s += "queue"; for (const std::string& q : queue) s += " " + q; s += "\n"; }
    for (const auto& [id, units] : progress) if (units > 0) s += "progress " + id + " " + std::to_string(units) + "\n";
    if (insight > 0.0) { char b[64]; snprintf(b, sizeof b, "insight %.3f\n", insight); s += b; }
    return s;
}

bool ResearchState::FromText(const std::string& text) {
    *this = ResearchState{};
    std::istringstream in(text);
    bool any = false;
    for (std::string line; std::getline(in, line);) {
        std::istringstream ls(line);
        std::string key; if (!(ls >> key)) continue;
        if (key == "done") { for (std::string w; ls >> w;) done.insert(w); any = true; }
        else if (key == "queue") { for (std::string w; ls >> w;) queue.push_back(w); any = true; }
        else if (key == "progress") { std::string id; int u = 0; if (ls >> id >> u && u > 0) { progress[id] = u; any = true; } }
        else if (key == "insight") { ls >> insight; any = true; }
    }
    return any;
}

bool Research_CanStart(const ResearchState& s, const ResearchNode& n) {
    if (s.Done(n.id)) return false;
    for (const std::string& d : n.needs) if (!s.Done(d)) return false;
    return true;
}

ResearchSight Research_Sight(const ResearchState& s, const ResearchNode& n, int teamTier) {
    if (s.Done(n.id) || n.tier <= teamTier) return ResearchSight::Full;
    if (n.tier == teamTier + 1) return ResearchSight::NameOnly;
    return ResearchSight::Hidden;
}

void Research_Merge(const ResearchState& from, ResearchState& into) {
    into.done.insert(from.done.begin(), from.done.end());
    for (const auto& [id, units] : from.progress) {
        if (into.done.count(id)) continue;
        int& u = into.progress[id];
        u = std::max(u, units);
    }
    for (auto it = into.progress.begin(); it != into.progress.end();) it = into.done.count(it->first) ? into.progress.erase(it) : std::next(it);
    into.queue.erase(std::remove_if(into.queue.begin(), into.queue.end(), [&](const std::string& q) { return into.done.count(q) > 0; }), into.queue.end());
    into.insight = std::max(into.insight, from.insight);
}

ResearchState& Research_Team(int team) { return teams[team]; }
void Research_Reset() { teams.clear(); }
std::string Research_TeamText(int team) { auto it = teams.find(team); return it == teams.end() ? std::string() : it->second.Text(); }
void Research_SetTeamText(int team, const std::string& text) { teams[team].FromText(text); }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1030): the node table, the state's save lines, sight and the join rule ----
#include <gtest/gtest.h>

TEST(Research, TheShippedTableLoadsAndValidates) {
    ASSERT_TRUE(Research_Load()) << "resources/data/research.txt (tests run from the repo root)";
    EXPECT_GT(Research_Nodes().size(), 40u);
    EXPECT_EQ(Research_Validate(), "");
    for (int t = 3; t <= 7; ++t) {
        int core = 0;
        for (const ResearchNode& n : Research_Nodes()) core += n.tier == t && n.core;
        EXPECT_GT(core, 0) << "tier " << t << " has core nodes (its milestone needs them)";
    }
}

TEST(Research, ValidationCatchesBadTables) {
    std::string error;
    ASSERT_TRUE(Research_LoadText(
        "node a | A | tier 3 | core | packs steam_pack | needs b\n"
        "node b | B | tier 3 | core | packs steam_pack | needs a\n"            // a cycle
        "node c | C | tier 3 | core | packs lv_pack | needs d\n"              // a later pack, a missing node
        "node e | E | tier 4 | core | packs steam_pack lv_pack | tome page\n" // a core node with a tome page
        "node f | F | tier 3 | branch magic | packs steam_pack | needs g\n"
        "node g | G | tier 5 | core | packs mv_pack\n",                        // f needs a later tier
        error)) << error;
    const std::string v = Research_Validate();
    EXPECT_NE(v.find("cycle"), std::string::npos) << v;
    EXPECT_NE(v.find("a later tier's pack"), std::string::npos);
    EXPECT_NE(v.find("no such node"), std::string::npos);
    EXPECT_NE(v.find("tome page"), std::string::npos);
    EXPECT_NE(v.find("of a later tier"), std::string::npos);
    EXPECT_FALSE(Research_LoadText("node x | X | colour red\n", error)) << "an unknown field is refused";
    Research_Load(); // (the shipped table back for the other tests)
}

TEST(Research, StateRoundTripsThroughItsSaveLines) {
    ResearchState s;
    s.done = { "steam_power", "wire_drawing" };
    s.queue = { "vacuum_tube", "electricity" };
    s.progress["vacuum_tube"] = 12;
    s.insight = 3.5;
    ResearchState back;
    ASSERT_TRUE(back.FromText(s.Text()));
    EXPECT_EQ(back.done, s.done);
    EXPECT_EQ(back.queue, s.queue);
    EXPECT_EQ(back.progress.at("vacuum_tube"), 12);
    EXPECT_DOUBLE_EQ(back.insight, 3.5);
    Research_Reset();
    Research_SetTeamText(7, s.Text());
    EXPECT_EQ(Research_TeamText(7), s.Text()) << "a team's block survives the record";
}

TEST(Research, StartSightAndTheJoinRule) {
    ASSERT_TRUE(Research_Load());
    const ResearchNode* vac = Research_Find("vacuum_tube");
    const ResearchNode* circuits = Research_Find("basic_circuits");
    const ResearchNode* aluminium = Research_Find("aluminium");
    ASSERT_TRUE(vac && circuits && aluminium);
    ResearchState s;
    EXPECT_FALSE(Research_CanStart(s, *vac)) << "wire drawing first";
    s.done.insert("wire_drawing");
    EXPECT_TRUE(Research_CanStart(s, *vac));
    // Sight at tier 3 (decided): tier 3 in full, tier 4 names only, tier 5 hidden.
    EXPECT_EQ(Research_Sight(s, *vac, 3), ResearchSight::Full);
    EXPECT_EQ(Research_Sight(s, *circuits, 3), ResearchSight::NameOnly);
    EXPECT_EQ(Research_Sight(s, *aluminium, 3), ResearchSight::Hidden);
    // A joiner's research goes to the team: done the union, progress the larger (dropped once done), the team's queue kept.
    ResearchState team; team.done = { "steam_power" }; team.queue = { "wire_drawing", "steelmaking" }; team.progress["steelmaking"] = 4;
    ResearchState joiner; joiner.done = { "wire_drawing" }; joiner.progress["steelmaking"] = 9; joiner.progress["vacuum_tube"] = 2; joiner.insight = 2.0;
    Research_Merge(joiner, team);
    EXPECT_TRUE(team.Done("steam_power") && team.Done("wire_drawing"));
    EXPECT_EQ(team.progress.at("steelmaking"), 9);
    EXPECT_EQ(team.progress.at("vacuum_tube"), 2);
    EXPECT_EQ(team.queue, std::vector<std::string>{ "steelmaking" }) << "a node now done leaves the queue";
    EXPECT_DOUBLE_EQ(team.insight, 2.0);
}

#include "playerrecord.h"
// The state rides in the world file with its team (team_research ... team_research_end) and comes back on load.
TEST(Research, ATeamsResearchIsSavedWithTheTeam) {
    Teams_Reset();
    const std::string me = PlayerId_New();
    const int team = Teams_CreateSolo(me, "testers");
    ASSERT_GT(team, 0);
    ResearchState& s = Research_Team(team);
    s.done = { "steam_power" }; s.queue = { "steelmaking" }; s.progress["steelmaking"] = 5;
    const std::string want = s.Text();
    const std::string path = "unit_research_teams.tmp";
    { FILE* f = fopen(path.c_str(), "w"); ASSERT_TRUE(f); Teams_Save(f, "", nullptr); fclose(f); }
    Teams_Reset();
    EXPECT_EQ(Research_TeamText(team), "") << "a reset forgets it";
    FILE* f = fopen(path.c_str(), "r");
    ASSERT_TRUE(f);
    char tok[64];
    while (fscanf(f, " %63s", tok) == 1) Teams_LoadLine(tok, f, "");
    fclose(f);
    std::remove(path.c_str());
    EXPECT_EQ(Research_TeamText(team), want);
    Teams_Reset();
}
#endif
