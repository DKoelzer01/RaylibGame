#include "teamstore.h"

#include "blueprint.h"
#include "navigation.h"
#include "research.h"
#include "teamquests.h"

#include <cstring>
#include <sstream>

namespace {
// A section: its name (the save keyword is team_<name>) and the module's own text form of one team's state.
struct Section {
    const char* name;
    std::string (*get)(int team);
    void (*load)(int team, const std::string& text);  // the world save's block back (null: the module saves it itself)
    void (*apply)(int team, const std::string& text); // a client: the host's text ("" when the host had none)
};
void ApplyMap(int, const std::string& text) { Nav_MergeTilesText(text); } // added to the local map: the player keeps their own
const Section kSections[] = {
    {"research", Research_TeamText, Research_SetTeamText, Research_SetTeamText}, // (#1030) ResearchState::Text lines
    {"map", TeamQuests_TilesText, nullptr, ApplyMap}, // (#338) saved with the team's book (team_tiles); Nav_TilesText lines
    {"blueprints", Blueprint_TeamText, Blueprint_SetTeamText, Blueprint_ApplyTeamText}, // the team library: .bp texts
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
        s->apply(team, body);
        seen[s - kSections] = true;
    }
    if (!any) return false;
    for (const Section& s : kSections) if (!seen[&s - kSections]) s.apply(team, "");
    return true;
}

void TeamStore_SaveTeam(FILE* f, int team) {
    for (const Section& s : kSections) {
        if (!s.load) continue; // the module saves it itself
        std::string body = s.get(team);
        if (body.empty()) continue;
        if (body.back() != '\n') body += '\n';
        fprintf(f, "team_%s %d\n%steam_%s_end\n", s.name, team, body.c_str(), s.name);
    }
}

bool TeamStore_LoadLine(const char* tok, FILE* f) {
    if (strncmp(tok, "team_", 5) != 0) return false;
    const Section* s = FindSection(tok + 5);
    if (!s || !s->load) return false;
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
    if (id > 0) s->load(id, text);
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

namespace {
const Vector3 kSpotA = { 0.0f, 1.0f, 0.0f };
const Vector3 kSpotB = { 0.0f, 0.0f, 1.0f }; // a quarter of the way round: no tiles in common with A
constexpr float kR = 60000.0f;
struct TeamMapTest : TeamStoreTest {
    void SetUp() override { TeamStoreTest::SetUp(); TeamQuests_Reset(); Nav_Reset(); }
    void TearDown() override { TeamQuests_Reset(); Nav_Reset(); TeamStoreTest::TearDown(); }
};
} // namespace

TEST_F(TeamMapTest, ClientMapGetsTheTeamsTilesAndKeepsItsOwn) {
    const int team = 9101;
    const int tiles = TeamQuests_Explore(team, "terra", kSpotA, kR);
    ASSERT_GT(tiles, 0);
    const std::string text = TeamStore_Text(team);
    EXPECT_NE(text.find("section map\n"), std::string::npos) << text;

    NavTiles mine; // the client explored somewhere else itself
    Nav_RevealTiles(kSpotB, kR, mine["terra"]);
    Nav_MergeTilesText(Nav_TilesText(mine));
    const int own = Nav_ExploredTiles("terra");
    ASSERT_TRUE(TeamStore_SetText(team, text));
    EXPECT_EQ(Nav_ExploredTiles("terra"), own + tiles) << "the team's tiles are added to the client's own";
    ASSERT_TRUE(TeamStore_SetText(team, "section research\nsection_end\n"));
    EXPECT_EQ(Nav_ExploredTiles("terra"), own + tiles) << "a text without a map never clears the local one";
}

TEST_F(TeamMapTest, AMoverCarriesTheGroundTheyExplored) {
    const int oldTeam = 9102, newTeam = 9103;
    const std::string uuid = PlayerId_FromName("mover");
    const int mine = TeamQuests_Explore(oldTeam, "terra", kSpotA, kR, uuid);
    ASSERT_GT(mine, 0);

    FILE* f = tmpfile(); // the player's tiles survive a world save (a reset forgets the teams' unstored tiles too)
    ASSERT_NE(f, nullptr);
    TeamQuests_SavePlayers(f);
    rewind(f);
    TeamQuests_Reset();
    char tok[64];
    ASSERT_EQ(fscanf(f, " %63s", tok), 1);
    ASSERT_TRUE(TeamQuests_LoadPlayerLine(tok, f));
    EXPECT_EQ(fscanf(f, " %63s", tok), EOF) << "the block read to its end";
    fclose(f);

    TeamQuests_Explore(oldTeam, "terra", kSpotB, kR); // a teammate's ground: stays with the old team
    const int theirs = TeamQuests_Explore(newTeam, "luna", kSpotA, 20000.0f);
    EXPECT_EQ(TeamQuests_CarryPlayer(uuid, newTeam), mine) << "the mover's own tiles join the new team";
    NavTiles got;
    Nav_TilesFromText(TeamQuests_TilesText(newTeam), got);
    EXPECT_EQ((int)got["terra"].size(), mine) << "only the mover's own ground, not the old teammate's";
    EXPECT_EQ((int)got["luna"].size(), theirs) << "what the team had stays";
    EXPECT_GE(TeamQuests_Explored(newTeam), (double)(mine + theirs)) << "\"explored\" rises to the union";
    EXPECT_EQ(TeamQuests_CarryPlayer(uuid, newTeam), 0) << "carrying twice adds nothing";
    EXPECT_EQ(TeamQuests_CarryPlayer(PlayerId_FromName("nobody"), newTeam), 0);
}

TEST(NavTiles, TextRoundTrips) {
    NavTiles a;
    Nav_RevealTiles(kSpotA, kR, a["terra"]);
    Nav_RevealTiles(kSpotB, 20000.0f, a["luna"]);
    NavTiles b;
    Nav_TilesFromText(Nav_TilesText(a), b);
    EXPECT_EQ(a, b);
}

namespace {
const char* kTower = "# Blueprint\nname tower\nparts 2\n0 0 0 0 0 0\n0 0 1 0 0 0\n";
const char* kWall = "name wall\r\nparts 1\r\n0 1 0 0 0 0\r\n"; // CRLF, as a Windows file
} // namespace

TEST_F(TeamStoreTest, BlueprintLibraryRoundTripsAndSaves) {
    const int team = 9201;
    Blueprint_AddToTeam(team, kTower);
    Blueprint_AddToTeam(team, kWall);
    Blueprint_AddToTeam(team, kTower); // the same one twice is kept once
    Blueprint_AddToTeam(team, "name empty\nparts 0\n"); // no parts: not kept
    ASSERT_EQ(Blueprint_TeamCount(team), 2);
    const std::string lib = Blueprint_TeamText(team);
    const std::string text = TeamStore_Text(team);
    EXPECT_NE(text.find("section blueprints\n"), std::string::npos) << text;

    FILE* f = tmpfile(); // the world save's block
    ASSERT_NE(f, nullptr);
    TeamStore_SaveTeam(f, team);
    rewind(f);
    Teams_Reset();
    EXPECT_EQ(Blueprint_TeamCount(team), 0) << "a reset forgets the libraries";
    char tok[64];
    ASSERT_EQ(fscanf(f, " %63s", tok), 1);
    EXPECT_STREQ(tok, "team_blueprints");
    ASSERT_TRUE(TeamStore_LoadLine(tok, f));
    fclose(f);
    EXPECT_EQ(Blueprint_TeamText(team), lib) << "loaded back as saved";

    Teams_Reset();
    ASSERT_TRUE(TeamStore_SetText(team, text)); // a client gets its team's library
    EXPECT_EQ(Blueprint_LocalTeam(), team) << "a client's N shows the team it was sent";
    EXPECT_EQ(Blueprint_TeamText(team), lib);
    ASSERT_TRUE(TeamStore_SetText(team, "section research\nsection_end\n"));
    EXPECT_EQ(Blueprint_TeamCount(team), 0) << "the host sent no library: the team has none";
}

TEST_F(TeamStoreTest, NListsTheTeamsBlueprintsFirst) {
    const int team = 9202;
    ASSERT_TRUE(TeamStore_SetText(team, "section blueprints\nblueprint\n" + std::string(kTower) + "blueprint_end\nsection_end\n"));
    const std::vector<Blueprint>& list = Blueprint_List();
    ASSERT_FALSE(list.empty());
    EXPECT_EQ(list.front().name, "tower");
    EXPECT_EQ(list.front().team, team);
    EXPECT_EQ(list.front().parts.size(), 2u);
    for (size_t i = 1; i < list.size(); ++i) EXPECT_EQ(list[i].team, 0) << "then files: the player's and the presets";
}
#endif
