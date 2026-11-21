#include "tome.h"

#include "magicsight.h"
#include "research.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <sstream>

namespace {
std::vector<AspectDef> aspects;
std::map<std::string, std::map<std::string, int>> sources; // type -> aspect -> amount
std::vector<TomePage> pages;
std::map<int, TomeState> teams;

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r\n");
    if (a == std::string::npos) return "";
    return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}
std::vector<std::string> Fields(const std::string& line) {
    std::vector<std::string> out;
    std::istringstream fs(line);
    for (std::string f; std::getline(fs, f, '|');) out.push_back(Trim(f));
    return out;
}
std::vector<std::string> Words(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream in(s);
    for (std::string w; in >> w;) if (w != "-") out.push_back(w);
    return out;
}
std::string ReadFile(const char* path) {
    std::ifstream f(path);
    if (!f) return {};
    std::stringstream ss; ss << f.rdbuf();
    return ss.str();
}
}

bool Aspects_LoadText(const std::string& text, std::string& error) {
    std::vector<AspectDef> outA;
    std::map<std::string, std::map<std::string, int>> outS;
    std::istringstream in(text);
    int lineNo = 0;
    for (std::string line; std::getline(in, line);) {
        ++lineNo;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto f = Fields(line);
        if (f[0].rfind("aspect ", 0) == 0 && f.size() >= 2) {
            AspectDef a;
            a.id = Trim(f[0].substr(7));
            a.name = f[1];
            for (size_t k = 2; k < f.size(); ++k) {
                if (f[k].rfind("latin ", 0) == 0) a.latin = Trim(f[k].substr(6));
                else if (f[k].rfind("colour ", 0) == 0) a.colour = (unsigned)strtoul(f[k].c_str() + 7, nullptr, 16);
            }
            outA.push_back(a);
        } else if (f[0].rfind("of ", 0) == 0 && f.size() >= 2) {
            auto& m = outS[Trim(f[0].substr(3))];
            std::istringstream ps(f[1]);
            for (std::string part; std::getline(ps, part, ',');) {
                std::istringstream ws(part);
                std::string id; int n = 0;
                if (!(ws >> id >> n) || n <= 0) { error = "aspects.txt line " + std::to_string(lineNo) + ": '" + Trim(part) + "' is not '<aspect> <n>'"; return false; }
                m[id] += n;
            }
        } else { error = "aspects.txt line " + std::to_string(lineNo) + ": not an aspect or of line"; return false; }
    }
    aspects = std::move(outA);
    sources = std::move(outS);
    return true;
}

bool Tome_LoadText(const std::string& text, std::string& error) {
    std::vector<TomePage> out;
    std::istringstream in(text);
    int lineNo = 0;
    for (std::string line; std::getline(in, line);) {
        ++lineNo;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const auto f = Fields(line);
        if (f.size() < 2 || f[0].rfind("page ", 0) != 0) { error = "tome.txt line " + std::to_string(lineNo) + ": not a page line"; return false; }
        TomePage p;
        p.id = Trim(f[0].substr(5));
        p.name = f[1];
        for (size_t k = 2; k < f.size(); ++k) {
            std::istringstream fs(f[k]);
            std::string key; fs >> key;
            std::string rest; std::getline(fs, rest); rest = Trim(rest);
            if (key == "latin") p.latin = rest;
            else if (key == "branch") p.branch = rest;
            else if (key == "text") p.text = rest;
            else if (key == "needs") p.needs = Words(rest);
            else if (key == "research") p.research = Words(rest);
            else if (key == "unlocks") p.unlocks = Words(rest);
            else if (key == "cost") {
                const auto w = Words(rest);
                if (w.size() % 2) { error = "tome.txt line " + std::to_string(lineNo) + ": cost is '<aspect> <n> ...'"; return false; }
                for (size_t i = 0; i + 1 < w.size(); i += 2) p.cost[w[i]] += atoi(w[i + 1].c_str());
            } else if (key == "effect") {
                std::istringstream es(rest);
                std::string what; es >> what;
                if (what == "sight_tier") es >> p.sightTier;
                else if (what == "through_rock") p.throughRock = true;
                else { error = "tome.txt line " + std::to_string(lineNo) + ": unknown effect '" + what + "'"; return false; }
            } else { error = "tome.txt line " + std::to_string(lineNo) + ": unknown field '" + key + "'"; return false; }
        }
        if (p.id.empty()) { error = "tome.txt line " + std::to_string(lineNo) + ": a page without an id"; return false; }
        out.push_back(std::move(p));
    }
    pages = std::move(out);
    return true;
}

bool Tome_Load(const char* aspectsPath, const char* tomePath) {
    std::string error;
    const std::string a = ReadFile(aspectsPath), t = ReadFile(tomePath);
    if (a.empty() || t.empty()) { printf("tome: cannot read %s or %s\n", aspectsPath, tomePath); return false; }
    if (!Aspects_LoadText(a, error) || !Tome_LoadText(t, error)) { printf("tome: %s\n", error.c_str()); return false; }
    const std::string problems = Tome_Validate();
    if (!problems.empty()) printf("tome: %s", problems.c_str());
    return true;
}

const std::vector<AspectDef>& Aspects_All() { return aspects; }
const AspectDef* Aspects_Find(const std::string& id) { for (const AspectDef& a : aspects) if (a.id == id) return &a; return nullptr; }
const std::map<std::string, int>* Aspects_Of(const std::string& type) { auto it = sources.find(type); return it == sources.end() ? nullptr : &it->second; }
const std::vector<TomePage>& Tome_Pages() { return pages; }
const TomePage* Tome_Find(const std::string& id) { for (const TomePage& p : pages) if (p.id == id) return &p; return nullptr; }

std::string Tome_Validate() {
    std::string out;
    std::set<std::string> ids;
    for (const AspectDef& a : aspects) if (!ids.insert(a.id).second) out += "aspect " + a.id + ": the id is used twice\n";
    for (const auto& [type, m] : sources) {
        if (type.find(':') == std::string::npos) out += "of " + type + ": a type is <kind>:<id>\n";
        for (const auto& [id, n] : m) if (!Aspects_Find(id)) out += "of " + type + ": no aspect '" + id + "'\n";
    }
    ids.clear();
    for (const TomePage& p : pages) {
        if (!ids.insert(p.id).second) out += "page " + p.id + ": the id is used twice\n";
        for (const auto& [id, n] : p.cost) if (!Aspects_Find(id)) out += "page " + p.id + ": costs '" + id + "', no such aspect\n";
        for (const std::string& d : p.needs) if (!Tome_Find(d)) out += "page " + p.id + ": needs '" + d + "', no such page\n";
        if (!Research_Nodes().empty())
            for (const std::string& r : p.research) if (!Research_Find(r)) out += "page " + p.id + ": needs research '" + r + "', no such node\n";
    }
    for (const ResearchNode& n : Research_Nodes())
        if (!n.tome.empty() && !Tome_Find(n.tome)) out += "research " + n.id + ": needs tome page '" + n.tome + "', no such page\n";
    std::map<std::string, int> mark;
    std::function<bool(const TomePage&)> visit = [&](const TomePage& p) {
        int& m = mark[p.id];
        if (m == 2) return true;
        if (m == 1) return false;
        m = 1;
        for (const std::string& d : p.needs) if (const TomePage* q = Tome_Find(d)) if (!visit(*q)) { out += "page " + p.id + ": a cycle through " + d + "\n"; return false; }
        m = 2;
        return true;
    };
    for (const TomePage& p : pages) visit(p);
    return out;
}

std::string TomeState::Text() const {
    if (pages.empty()) return {};
    std::string s = "pages";
    for (const std::string& p : pages) s += " " + p;
    return s + "\n";
}
bool TomeState::FromText(const std::string& text) {
    pages.clear();
    std::istringstream in(text);
    bool any = false;
    for (std::string line; std::getline(in, line);) {
        std::istringstream ls(line);
        std::string key; if (!(ls >> key) || key != "pages") continue;
        for (std::string w; ls >> w;) { pages.insert(w); any = true; }
    }
    return any;
}
TomeState& Tome_Team(int team) { return teams[team]; }
void Tome_Reset() { teams.clear(); }
std::string Tome_TeamText(int team) { auto it = teams.find(team); return it == teams.end() ? std::string() : it->second.Text(); }
void Tome_SetTeamText(int team, const std::string& text) { teams[team].FromText(text); }

void Tome_GrantStudy(MagicPlayerState& me, const std::string& type, float value) {
    const auto* m = Aspects_Of(type);
    if (!m || value <= 0.0f) return;
    for (const auto& [id, n] : *m) me.aspects[id] += value * (float)n;
}

bool Tome_CanUnlock(const TomeState& team, const MagicPlayerState& me, const ResearchState& research, const TomePage& p, std::string* why) {
    auto no = [&](const std::string& w) { if (why) *why = w; return false; };
    if (team.Has(p.id)) return no("the team has it");
    for (const std::string& d : p.needs) if (!team.Has(d)) { const TomePage* q = Tome_Find(d); return no("needs the page " + (q ? q->name : d)); }
    for (const std::string& r : p.research) if (!research.Done(r)) { const ResearchNode* n = Research_Find(r); return no("needs the research " + (n ? n->name : r)); }
    for (const auto& [id, n] : p.cost) {
        const auto it = me.aspects.find(id);
        const float have = it == me.aspects.end() ? 0.0f : it->second;
        if (have + 1e-4f < (float)n) { const AspectDef* a = Aspects_Find(id); return no("needs " + std::to_string(n) + " " + (a ? a->name : id)); }
    }
    return true;
}

bool Tome_Unlock(TomeState& team, MagicPlayerState& me, const ResearchState& research, const std::string& page) {
    const TomePage* p = Tome_Find(page);
    if (!p || !Tome_CanUnlock(team, me, research, *p)) return false;
    for (const auto& [id, n] : p->cost) {
        float& have = me.aspects[id];
        have = std::max(0.0f, have - (float)n);
        if (have <= 1e-4f) me.aspects.erase(id);
    }
    team.pages.insert(p->id);
    return true;
}

void Tome_ApplyEffects(const TomeState& team, MagicPlayerState& me) {
    for (const std::string& id : team.pages)
        if (const TomePage* p = Tome_Find(id)) {
            if (p->sightTier >= 0) me.sightTier = std::max(me.sightTier, p->sightTier);
            if (p->throughRock) me.throughRock = true;
            if (p->id == "first_sight") me.sight = true; // a teammate's first find teaches the whole team (pages per team)
        }
}

bool Tome_ResearchAllowed(const ResearchNode& n, const TomeState& team) { return n.tome.empty() || team.Has(n.tome); }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1251): the shipped tables, studies give aspects, pages cost them, the research cross link ----
#include <gtest/gtest.h>
#include "unit_helpers.h"
#include "inventory.h"

TEST(Tome, TheShippedTablesLoadAndValidate) {
    ASSERT_TRUE(Research_Load());
    ASSERT_TRUE(Tome_Load()) << "resources/data/aspects.txt and tome.txt (tests run from the repo root)";
    EXPECT_EQ(Tome_Validate(), "");
    EXPECT_EQ(Aspects_All().size(), 16u) << "the doc's 16 aspects to start";
    ASSERT_NE(Tome_Find("first_sight"), nullptr) << "discovery's page, and research's resonance tuning needs it";
    EXPECT_NE(Aspects_Of("node:storm"), nullptr);
    for (const char* k : { "earth", "water", "fire", "air", "storm", "shard", "aether" })
        EXPECT_NE(Aspects_Of(std::string("node:") + k), nullptr) << "every first-slice node kind can be studied: " << k;
}

TEST(Tome, EveryItemWithAspectsExists) {
    UnitTest::EnsureItemDB();
    ASSERT_TRUE(Tome_Load());
    ItemDB withMagic;
    ASSERT_TRUE(withMagic.load("resources/data/items.txt", "resources/data/recipes.txt"));
    withMagic.loadExtra("resources/data/magic_items.txt");
    for (const char* t : { "item:iron_ingot", "item:quartz", "item:resonant_crystal", "item:oak_log" }) EXPECT_NE(Aspects_Of(t), nullptr) << t;
    // Every "of item:" line names a real item (magic items included).
    std::string text = ReadFile("resources/data/aspects.txt");
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        line = Trim(line);
        if (line.rfind("of item:", 0) != 0) continue;
        const std::string id = Trim(Fields(line)[0].substr(8));
        EXPECT_NE(withMagic.find(id), nullptr) << "aspects.txt names an unknown item: " << id;
    }
}

TEST(Tome, ValidationCatchesBadTables) {
    std::string error;
    ASSERT_TRUE(Aspects_LoadText("aspect a | A\naspect a | A again\nof item:x | a 1, b 2\nof nokind | a 1\n", error)) << error;
    ASSERT_TRUE(Tome_LoadText("page p | P | needs q\npage q | Q | needs p | cost zz 3\npage r | R | needs missing\n", error)) << error;
    const std::string v = Tome_Validate();
    EXPECT_NE(v.find("used twice"), std::string::npos) << v;
    EXPECT_NE(v.find("no aspect 'b'"), std::string::npos) << v;
    EXPECT_NE(v.find("<kind>:<id>"), std::string::npos) << v;
    EXPECT_NE(v.find("cycle"), std::string::npos) << v;
    EXPECT_NE(v.find("no such aspect"), std::string::npos) << v;
    EXPECT_NE(v.find("no such page"), std::string::npos) << v;
    EXPECT_FALSE(Tome_LoadText("page x | X | colour red\n", error)) << "an unknown field is refused";
    EXPECT_FALSE(Tome_LoadText("page x | X | cost a\n", error)) << "a cost needs amounts";
    EXPECT_FALSE(Aspects_LoadText("of item:x | a two\n", error));
    ASSERT_TRUE(Tome_Load()); // (the shipped tables back)
}

TEST(Tome, StudiesGiveAspectsAndPagesSpendThem) {
    ASSERT_TRUE(Research_Load());
    ASSERT_TRUE(Tome_Load());
    MagicPlayerState me;
    TomeState team;
    ResearchState research;
    team.pages.insert("first_sight");
    const TomePage* wider = Tome_Find("wider_sight");
    ASSERT_TRUE(wider);
    std::string why;
    EXPECT_FALSE(Tome_CanUnlock(team, me, research, *wider, &why));
    EXPECT_NE(why.find("needs"), std::string::npos);
    // Study air nodes and torches' light... until the cost is met: full value the first time, half after.
    for (int i = 0; i < 4; ++i) Tome_GrantStudy(me, "node:air", 1.0f);      // air 12, motion 4
    for (int i = 0; i < 2; ++i) Tome_GrantStudy(me, "item:torch", 1.0f);    // light 4, fire 2
    for (int i = 0; i < 2; ++i) Tome_GrantStudy(me, "item:wool", 1.0f);     // order 2
    EXPECT_FLOAT_EQ(me.aspects["air"], 12.0f);
    EXPECT_TRUE(Tome_CanUnlock(team, me, research, *wider, &why)) << why;
    EXPECT_TRUE(Tome_Unlock(team, me, research, "wider_sight"));
    EXPECT_TRUE(team.Has("wider_sight"));
    EXPECT_FLOAT_EQ(me.aspects["air"], 6.0f) << "the cost is spent";
    EXPECT_EQ(me.aspects.count("light"), 0u) << "spent to nothing: gone";
    EXPECT_FALSE(Tome_Unlock(team, me, research, "wider_sight")) << "once per team";
    // Effects reach every member through the team's pages.
    MagicPlayerState mate;
    Tome_ApplyEffects(team, mate);
    EXPECT_TRUE(mate.sight) << "the team's First Sight teaches its members";
    EXPECT_EQ(mate.sightTier, 1);
    // A page that needs research: Insight needs the LV generator (and research's mana dynamo needs Insight).
    const TomePage* insight = Tome_Find("insight");
    ASSERT_TRUE(insight);
    team.pages.insert("research_table");
    for (const auto& [id, n] : insight->cost) me.aspects[id] = 99.0f;
    EXPECT_FALSE(Tome_CanUnlock(team, me, research, *insight, &why));
    EXPECT_NE(why.find("research"), std::string::npos) << why;
    research.done.insert("lv_power");
    EXPECT_TRUE(Tome_CanUnlock(team, me, research, *insight, &why)) << why;
}

TEST(Tome, ResearchNeedsItsTomePage) {
    ASSERT_TRUE(Research_Load());
    ASSERT_TRUE(Tome_Load());
    const ResearchNode* tuning = Research_Find("resonance_tuning");
    ASSERT_TRUE(tuning);
    TomeState team;
    EXPECT_FALSE(Tome_ResearchAllowed(*tuning, team)) << "resonance tuning needs First Sight (the Factory doc)";
    team.pages.insert("first_sight");
    EXPECT_TRUE(Tome_ResearchAllowed(*tuning, team));
    const ResearchNode* steam = Research_Find("steam_power");
    ASSERT_TRUE(steam);
    EXPECT_TRUE(Tome_ResearchAllowed(*steam, TomeState{})) << "no page needed elsewhere";
}

TEST(Tome, ATeamsPagesRoundTrip) {
    TomeState a; a.pages = { "first_sight", "wider_sight" };
    TomeState b;
    ASSERT_TRUE(b.FromText(a.Text()));
    EXPECT_EQ(b.pages, a.pages);
    Tome_Reset();
    Tome_SetTeamText(4, a.Text());
    EXPECT_EQ(Tome_TeamText(4), a.Text());
    EXPECT_EQ(Tome_TeamText(5), "");
}
#endif
