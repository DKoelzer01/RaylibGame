#include "teamstore.h"

#include "research.h"

#include <cstring>
#include <sstream>

namespace {
// A section: its name (the save keyword is team_<name>) and the module's own text form of one team's state.
struct Section {
    const char* name;
    std::string (*get)(int team);
    void (*set)(int team, const std::string& text); // "" clears the team's state
};
const Section kSections[] = {
    {"research", Research_TeamText, Research_SetTeamText}, // (#1030) ResearchState::Text lines
};

const Section* FindSection(const std::string& name) {
    for (const Section& s : kSections) if (name == s.name) return &s;
    return nullptr;
}
} // namespace

std::string TeamStore_Text(int team) {
    std::string out;
    for (const Section& s : kSections) {
        std::string body = s.get(team);
        if (body.empty()) continue;
        if (body.back() != '\n') body += '\n';
        out += "section "; out += s.name; out += "\n"; out += body; out += "section_end\n";
    }
    return out;
}

bool TeamStore_SetText(int team, const std::string& text) {
    bool seen[sizeof(kSections) / sizeof(kSections[0])] = {};
    bool any = false;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.compare(0, 8, "section ") != 0) continue;
        any = true;
        const Section* s = FindSection(line.substr(8));
        std::string body;
        for (std::string l; std::getline(in, l);) {
            if (!l.empty() && l.back() == '\r') l.pop_back();
            if (l == "section_end") break;
            body += l; body += '\n';
        }
        if (!s) continue; // a section this build doesn't know
        s->set(team, body);
        seen[s - kSections] = true;
    }
    if (!any) return false;
    for (const Section& s : kSections) if (!seen[&s - kSections]) s.set(team, "");
    return true;
}

void TeamStore_SaveTeam(FILE* f, int team) {
    for (const Section& s : kSections) {
        std::string body = s.get(team);
        if (body.empty()) continue;
        if (body.back() != '\n') body += '\n';
        fprintf(f, "team_%s %d\n%steam_%s_end\n", s.name, team, body.c_str(), s.name);
    }
}

bool TeamStore_LoadLine(const char* tok, FILE* f) {
    if (strncmp(tok, "team_", 5) != 0) return false;
    const Section* s = FindSection(tok + 5);
    if (!s) return false;
    int id = 0;
    if (fscanf(f, " %d", &id) != 1) return true;
    const std::string end = std::string("team_") + s->name + "_end";
    // Read whole lines (a research "done" line outgrows any fixed buffer) up to the end keyword.
    std::string text, line;
    bool skippedRest = false; // the rest of the header line comes first
    for (int c; (c = fgetc(f)) != EOF;) {
        if (c != '\n') { line += (char)c; continue; }
        if (!skippedRest) { skippedRest = true; line.clear(); continue; }
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t a = line.find_first_not_of(" \t");
        if (a != std::string::npos && line.compare(a, end.size(), end) == 0) break;
        text += line; text += '\n';
        line.clear();
    }
    if (id > 0) s->set(id, text);
    return true;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1886): the store's text round trip, team separation and the save block ----
#include <gtest/gtest.h>

#include "playerrecord.h"

namespace {
struct TeamStoreTest : ::testing::Test {
    void SetUp() override { Teams_Reset(); }
    void TearDown() override { Teams_Reset(); }
};
} // namespace

TEST_F(TeamStoreTest, TextRoundTripsResearch) {
    ResearchState& a = Research_Team(3);
    a.done = {"basic_tools", "smelting"};
    a.queue = {"steam_power"};
    a.progress["steam_power"] = 12;
    const std::string text = TeamStore_Text(3);
    EXPECT_NE(text.find("section research\n"), std::string::npos) << text;
    const std::string want = Research_TeamText(3);

    Teams_Reset();
    ASSERT_TRUE(TeamStore_SetText(3, text));
    EXPECT_EQ(Research_TeamText(3), want) << "a client applying the host's text gets the same research";
    EXPECT_EQ(TeamStore_Text(3), text);
}

TEST_F(TeamStoreTest, TeamsStaySeparate) {
    Research_Team(1).done = {"smelting"};
    Research_Team(2).done = {"basic_tools"};
    const std::string one = TeamStore_Text(1);
    EXPECT_EQ(one.find("basic_tools"), std::string::npos) << "team 2's research isn't in team 1's text";
    ASSERT_TRUE(TeamStore_SetText(2, one));
    EXPECT_TRUE(Research_Team(2).done.count("smelting"));
    EXPECT_FALSE(Research_Team(2).done.count("basic_tools")) << "applying replaces, it doesn't merge";
    EXPECT_TRUE(Research_Team(1).done.count("smelting")) << "team 1 untouched";
}

TEST_F(TeamStoreTest, MissingSectionClearsUnknownIsSkipped) {
    Research_Team(4).done = {"smelting"};
    ASSERT_TRUE(TeamStore_SetText(4, "section future_thing\nwhatever 1\nsection_end\n"));
    EXPECT_EQ(Research_TeamText(4), "") << "the host sent no research: the team has none";
    EXPECT_FALSE(TeamStore_SetText(4, "")) << "no section at all is not a valid text";
    EXPECT_TRUE(TeamStore_Text(4).empty());
}

TEST_F(TeamStoreTest, SaveBlockLoadsBack) {
    ResearchState& a = Research_Team(5);
    a.done = {"smelting", "basic_tools"};
    a.insight = 2.5;
    const std::string want = Research_TeamText(5);
    FILE* f = tmpfile();
    ASSERT_NE(f, nullptr);
    TeamStore_SaveTeam(f, 5);
    fprintf(f, "after 1\n");
    rewind(f);
    Teams_Reset();
    char tok[64];
    ASSERT_EQ(fscanf(f, " %63s", tok), 1);
    EXPECT_STREQ(tok, "team_research") << "research keeps its #1030 save keyword";
    ASSERT_TRUE(TeamStore_LoadLine(tok, f));
    EXPECT_EQ(Research_TeamText(5), want);
    ASSERT_EQ(fscanf(f, " %63s", tok), 1);
    EXPECT_STREQ(tok, "after") << "the block stops at its end line";
    EXPECT_FALSE(TeamStore_LoadLine("team_member", f)) << "not a store keyword";
    fclose(f);
}
#endif
