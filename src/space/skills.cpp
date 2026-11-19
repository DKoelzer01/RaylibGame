#include "skills.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>

namespace {
    std::vector<SkillDef> skills;
    std::vector<PerkDef> perks;
    bool loaded = false;

    // A line's fields: words, and "quoted strings" as one field.
    std::vector<std::string> Fields(const std::string& line) {
        std::vector<std::string> out;
        size_t i = 0;
        while (i < line.size()) {
            while (i < line.size() && isspace((unsigned char)line[i])) ++i;
            if (i >= line.size() || line[i] == '#') break;
            if (line[i] == '"') { const size_t e = line.find('"', i + 1); out.push_back(line.substr(i + 1, (e == std::string::npos ? line.size() : e) - i - 1)); i = e == std::string::npos ? line.size() : e + 1; }
            else { size_t e = i; while (e < line.size() && !isspace((unsigned char)line[e])) ++e; out.push_back(line.substr(i, e - i)); i = e; }
        }
        return out;
    }
    SkillDef* FindSkillMut(const std::string& id) { for (SkillDef& d : skills) if (d.id == id) return &d; return nullptr; }
    PerkDef* FindPerkMut(const std::string& id) { for (PerkDef& p : perks) if (p.id == id) return &p; return nullptr; }
}

bool Skills_LoadText(const std::string& text) {
    skills.clear(); perks.clear();
    std::istringstream in(text);
    std::string line;
    bool ok = true;
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const std::vector<std::string> f = Fields(line);
        if (f.empty()) continue;
        if (f[0] == "skill" && f.size() >= 3) skills.push_back({ f[1], f[2], {} });
        else if (f[0] == "xp" && f.size() >= 4) { if (SkillDef* d = FindSkillMut(f[1])) d->xp[f[2]] = atof(f[3].c_str()); else ok = false; }
        else if (f[0] == "perk" && f.size() >= 8) {
            PerkDef p; p.id = f[1]; p.skill = f[2]; p.tier = std::clamp(atoi(f[3].c_str()), 1, 3); p.cost = std::max(1, atoi(f[4].c_str()));
            p.ranks = std::max(1, atoi(f[5].c_str())); p.name = f[6]; p.text = f[7];
            if (!FindSkillMut(p.skill)) ok = false; else perks.push_back(std::move(p));
        } else if (f[0] == "needs" && f.size() >= 3) { if (PerkDef* p = FindPerkMut(f[1])) p->needs.push_back(f[2]); else ok = false; }
        else if (f[0] == "effect" && f.size() >= 4) { if (PerkDef* p = FindPerkMut(f[1])) p->effects.push_back({ f[2], atof(f[3].c_str()) }); else ok = false; }
        else if (f[0] == "ability" && f.size() >= 3) { if (PerkDef* p = FindPerkMut(f[1])) p->ability = f[2]; else ok = false; }
        else ok = false;
    }
    loaded = true;
    return ok;
}

bool Skills_Load(const char* path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) { loaded = true; return false; }
    std::stringstream ss; ss << f.rdbuf();
    return Skills_LoadText(ss.str());
}

const std::vector<SkillDef>& Skills_All() { if (!loaded) Skills_Load(); return skills; }
const std::vector<PerkDef>& Skills_Perks() { if (!loaded) Skills_Load(); return perks; }
const SkillDef* Skills_Find(const std::string& id) { for (const SkillDef& d : Skills_All()) if (d.id == id) return &d; return nullptr; }
const PerkDef* Skills_FindPerk(const std::string& id) { for (const PerkDef& p : Skills_Perks()) if (p.id == id) return &p; return nullptr; }

double Skills_XpForLevel(int n) { return n < 1 ? 0.0 : 100.0 * std::pow((double)n, 1.4); }

int Skills_LevelOf(double xp) {
    int n = 0;
    double need = Skills_XpForLevel(1);
    while (xp >= need && n < 100000) { xp -= need; ++n; need = Skills_XpForLevel(n + 1); }
    return n;
}

double Skills_Progress(double xp) {
    int n = 0;
    double need = Skills_XpForLevel(1);
    while (xp >= need && n < 100000) { xp -= need; ++n; need = Skills_XpForLevel(n + 1); }
    return need > 0.0 ? xp / need : 0.0;
}

int Skills_AddXp(SkillState& s, const std::string& skill, double xp, double rate) {
    if (xp <= 0.0 || rate <= 0.0 || !Skills_Find(skill)) return 0;
    double& have = s.xp[skill];
    const int before = Skills_LevelOf(have);
    have += xp * rate;
    return Skills_LevelOf(have) - before;
}

int Skills_Event(SkillState& s, const std::string& event, double n, double rate, std::vector<std::string>* levelled) {
    int up = 0;
    for (const SkillDef& d : Skills_All()) {
        auto it = d.xp.find(event);
        if (it == d.xp.end()) continue;
        const int k = Skills_AddXp(s, d.id, it->second * n, rate);
        if (k > 0 && levelled) levelled->push_back(d.id);
        up += k;
    }
    return up;
}

namespace { constexpr double kFillWindow = 600.0; } // s: ground the player filled stays "theirs" this long
void Skills_NoteFill(SkillState& s, uint64_t cell, double now) {
    if (s.filledCells.size() > 4096) // (old ones out now and then: the map stays small)
        for (auto it = s.filledCells.begin(); it != s.filledCells.end();) { if (now - it->second > kFillWindow) it = s.filledCells.erase(it); else ++it; }
    s.filledCells[cell] = now;
}

double Skills_DigShare(const SkillState& s, uint64_t cell, double now) {
    auto it = s.filledCells.find(cell);
    return it != s.filledCells.end() && now - it->second < kFillWindow ? 0.1 : 1.0;
}

void Skills_OnDeath(SkillState& s, double share) {
    share = std::clamp(share, 0.0, 1.0);
    for (auto& [id, xp] : s.xp) {
        const int n = Skills_LevelOf(xp);
        double floor = 0.0; // the XP the current level started at
        for (int k = 1; k <= n; ++k) floor += Skills_XpForLevel(k);
        xp = floor + (xp - floor) * (1.0 - share);
    }
}

int Skills_Level(const SkillState& s, const std::string& skill) { auto it = s.xp.find(skill); return it == s.xp.end() ? 0 : Skills_LevelOf(it->second); }

int Skills_PointsSpent(const SkillState& s, const std::string& skill) {
    int n = 0;
    for (const auto& [id, r] : s.ranks) if (const PerkDef* p = Skills_FindPerk(id); p && p->skill == skill) n += r * p->cost;
    return n;
}

int Skills_PointsFree(const SkillState& s, const std::string& skill) { return Skills_Level(s, skill) - Skills_PointsSpent(s, skill); }

std::string Skills_CanBuy(const SkillState& s, const std::string& perk) {
    const PerkDef* p = Skills_FindPerk(perk);
    if (!p) return "No such perk";
    auto rit = s.ranks.find(perk);
    const int rank = rit == s.ranks.end() ? 0 : rit->second;
    if (rank >= p->ranks) return "Fully learnt";
    const int need = kSkillTierLevel[p->tier - 1];
    if (Skills_Level(s, p->skill) < need) return "Needs " + (Skills_Find(p->skill) ? Skills_Find(p->skill)->name : p->skill) + " level " + std::to_string(need);
    for (const std::string& q : p->needs) {
        auto it = s.ranks.find(q);
        if (it == s.ranks.end() || it->second < 1) { const PerkDef* qd = Skills_FindPerk(q); return "Needs " + (qd ? qd->name : q); }
    }
    if (Skills_PointsFree(s, p->skill) < p->cost) return p->cost == 1 ? "Needs a free point" : "Needs " + std::to_string(p->cost) + " free points";
    return "";
}

bool Skills_Buy(SkillState& s, const std::string& perk) {
    if (!Skills_CanBuy(s, perk).empty()) return false;
    ++s.ranks[perk];
    return true;
}

void Skills_Respec(SkillState& s, const std::string& skill) {
    for (auto it = s.ranks.begin(); it != s.ranks.end();) {
        const PerkDef* p = Skills_FindPerk(it->first);
        if (p && p->skill == skill) it = s.ranks.erase(it); else ++it;
    }
}

double Skills_Stat(const SkillState& s, const std::string& stat) {
    double m = 1.0;
    for (const auto& [id, r] : s.ranks) {
        const PerkDef* p = Skills_FindPerk(id);
        if (!p || r <= 0) continue;
        for (const auto& [st, per] : p->effects) if (st == stat) m += per * std::min(r, p->ranks); // capped by the perk's ranks
    }
    return std::max(0.1, m);
}

std::string Skills_SaveText(const SkillState& s) {
    std::string out;
    char buf[160];
    for (const auto& [id, xp] : s.xp) if (xp > 0.0) { snprintf(buf, sizeof(buf), "skill_xp %s %.3f\n", id.c_str(), xp); out += buf; }
    for (const auto& [id, r] : s.ranks) if (r > 0) { snprintf(buf, sizeof(buf), "perk %s %d\n", id.c_str(), r); out += buf; }
    return out;
}

bool Skills_LoadLine(SkillState& s, const char* line) {
    char id[96] = {};
    double xp = 0.0;
    int r = 0;
    if (sscanf(line, "skill_xp %95s %lf", id, &xp) == 2) { s.xp[id] = std::max(0.0, xp); return true; }
    if (sscanf(line, "perk %95s %d", id, &r) == 2) { if (r > 0) s.ranks[id] = r; return true; }
    return false;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#883): the curve, XP and events, anti-farming, death, perks and their rules, stats, saves ----
#include "unit_helpers.h"

namespace {
const char* kTestSkills = R"(
skill mining "Mining"
skill digging "Digging"
xp mining rock_stroke 1
xp mining ore_stroke 3
xp digging soil_stroke 1
xp digging rock_stroke 0.5
perk quick_pick mining 1 1 3 "Quick Pick" "Rock strokes 10% faster"
effect quick_pick rock_speed 0.10
perk light_hands mining 1 1 2 "Light Hands" "Picks wear 10% slower"
effect light_hands pick_wear -0.10
perk gem_eye mining 2 2 2 "Gem Eye" "Gems 25% more often"
needs gem_eye quick_pick
effect gem_eye gem_chance 0.25
)";
struct SkillsTest : ::testing::Test { void SetUp() override { ASSERT_TRUE(Skills_LoadText(kTestSkills)); } };
double XpToReach(int n) { double t = 0.0; for (int k = 1; k <= n; ++k) t += Skills_XpForLevel(k); return t; }
}

TEST_F(SkillsTest, CurveFollowsTheDesign) {
    EXPECT_NEAR(Skills_XpForLevel(1), 100.0, 1e-9);
    EXPECT_NEAR(Skills_XpForLevel(10), 2511.9, 0.5) << "level 10 ~2,500";
    EXPECT_NEAR(Skills_XpForLevel(30), 11690.0, 15.0) << "level 30 ~11,700";
    EXPECT_EQ(Skills_LevelOf(0.0), 0);
    EXPECT_EQ(Skills_LevelOf(99.9), 0);
    EXPECT_EQ(Skills_LevelOf(XpToReach(12)), 12);
    EXPECT_EQ(Skills_LevelOf(XpToReach(12) - 0.01), 11);
    EXPECT_NEAR(Skills_Progress(XpToReach(3) + Skills_XpForLevel(4) / 2.0), 0.5, 1e-9);
    EXPECT_EQ(Skills_LevelOf(XpToReach(200)), 200) << "no cap";
}

TEST_F(SkillsTest, EventsPayEverySkillThatCountsThemAtTheWorldRate) {
    SkillState s;
    std::vector<std::string> up;
    EXPECT_EQ(Skills_Event(s, "rock_stroke", 100.0, 1.0, &up), 1);
    EXPECT_NEAR(s.xp["mining"], 100.0, 1e-9);
    EXPECT_NEAR(s.xp["digging"], 50.0, 1e-9) << "a rock stroke counts for digging too";
    ASSERT_EQ(up.size(), 1u); EXPECT_EQ(up[0], "mining");
    Skills_Event(s, "ore_stroke", 10.0, 2.0);
    EXPECT_NEAR(s.xp["mining"], 160.0, 1e-9) << "3 XP x 10 x the 2x world rate";
    EXPECT_EQ(Skills_Event(s, "nothing_counts_this"), 0);
    EXPECT_EQ(Skills_AddXp(s, "no_such_skill", 1000.0), 0);
}

TEST_F(SkillsTest, GroundYouFilledPaysLittleWhenDugOutAgain) {
    SkillState s;
    EXPECT_DOUBLE_EQ(Skills_DigShare(s, 42, 0.0), 1.0) << "untouched ground pays in full";
    EXPECT_DOUBLE_EQ(Skills_DigShare(s, 42, 1.0), 1.0) << "digging the same spot again (a hole going down) still pays";
    Skills_NoteFill(s, 42, 10.0);
    EXPECT_DOUBLE_EQ(Skills_DigShare(s, 42, 70.0), 0.1) << "ground you filled, dug out again within 10 minutes";
    EXPECT_DOUBLE_EQ(Skills_DigShare(s, 43, 70.0), 1.0) << "another cell";
    EXPECT_DOUBLE_EQ(Skills_DigShare(s, 42, 700.0), 1.0) << "after 10 minutes it is ground like any other";
}

TEST_F(SkillsTest, DeathCostsProgressNeverALevel) {
    SkillState s;
    s.xp["mining"] = XpToReach(5) + Skills_XpForLevel(6) * 0.8;
    Skills_OnDeath(s, 0.25);
    EXPECT_EQ(Skills_Level(s, "mining"), 5);
    EXPECT_NEAR(Skills_Progress(s.xp["mining"]), 0.6, 1e-6);
    s.xp["digging"] = XpToReach(3); // exactly at a level: nothing to lose
    Skills_OnDeath(s, 0.5);
    EXPECT_EQ(Skills_Level(s, "digging"), 3);
}

TEST_F(SkillsTest, PerksNeedTheTierLevelPrerequisitesAndPoints) {
    SkillState s;
    EXPECT_EQ(Skills_CanBuy(s, "quick_pick"), "Needs Mining level 5");
    s.xp["mining"] = XpToReach(5);
    EXPECT_EQ(Skills_PointsFree(s, "mining"), 5);
    EXPECT_TRUE(Skills_Buy(s, "quick_pick"));
    EXPECT_TRUE(Skills_Buy(s, "quick_pick"));
    EXPECT_TRUE(Skills_Buy(s, "quick_pick"));
    EXPECT_EQ(Skills_CanBuy(s, "quick_pick"), "Fully learnt");
    EXPECT_EQ(Skills_PointsFree(s, "mining"), 2);
    EXPECT_EQ(Skills_CanBuy(s, "gem_eye"), "Needs Mining level 15") << "tier 2";
    s.xp["mining"] = XpToReach(15);
    SkillState fresh; fresh.xp["mining"] = XpToReach(15);
    EXPECT_EQ(Skills_CanBuy(fresh, "gem_eye"), "Needs Quick Pick");
    EXPECT_TRUE(Skills_Buy(s, "gem_eye"));
    EXPECT_EQ(Skills_PointsSpent(s, "mining"), 5) << "3 ranks x 1 + 1 rank x 2";
    EXPECT_EQ(Skills_PointsFree(s, "digging"), 0) << "points belong to their own skill";
    Skills_Respec(s, "mining");
    EXPECT_EQ(Skills_PointsFree(s, "mining"), 15);
    EXPECT_TRUE(s.ranks.empty());
}

TEST_F(SkillsTest, StatsMultiplyByRanksAndWearNeverDropsBelowATenth) {
    SkillState s;
    EXPECT_DOUBLE_EQ(Skills_Stat(s, "rock_speed"), 1.0);
    s.ranks["quick_pick"] = 2;
    EXPECT_NEAR(Skills_Stat(s, "rock_speed"), 1.2, 1e-9) << "+10% a rank";
    s.ranks["quick_pick"] = 9; // more than its ranks (a hand-edited save): capped
    EXPECT_NEAR(Skills_Stat(s, "rock_speed"), 1.3, 1e-9);
    s.ranks["light_hands"] = 2;
    EXPECT_NEAR(Skills_Stat(s, "pick_wear"), 0.8, 1e-9);
    const double suitModule = 1.5; // (Q8) a module's multiplier and a perk's stack: both multiply
    EXPECT_NEAR(suitModule * Skills_Stat(s, "rock_speed"), 1.95, 1e-9);
}

TEST_F(SkillsTest, SaveTextRoundTrips) {
    SkillState a;
    a.xp["mining"] = 1234.5; a.xp["digging"] = 7.25; a.ranks["quick_pick"] = 2;
    SkillState b;
    std::istringstream in(Skills_SaveText(a));
    std::string line;
    int lines = 0;
    while (std::getline(in, line)) { EXPECT_TRUE(Skills_LoadLine(b, line.c_str())); ++lines; }
    EXPECT_EQ(lines, 3);
    EXPECT_NEAR(b.xp["mining"], 1234.5, 1e-3);
    EXPECT_NEAR(b.xp["digging"], 7.25, 1e-3);
    EXPECT_EQ(b.ranks["quick_pick"], 2);
    EXPECT_FALSE(Skills_LoadLine(b, "team 3 something"));
}

TEST(SkillsData, TheShippedTableLoadsCleanly) {
    UnitTest::EnsureItemDB();
    EXPECT_TRUE(Skills_Load("resources/data/skills.txt"));
    EXPECT_GE(Skills_All().size(), 8u) << "slice 1's skills";
    for (const PerkDef& p : Skills_Perks()) {
        EXPECT_NE(Skills_Find(p.skill), nullptr) << p.id;
        for (const std::string& q : p.needs) EXPECT_NE(Skills_FindPerk(q), nullptr) << p.id << " needs " << q;
    }
}

#include <set>
TEST(SkillsData, EverySkillHasATreeAndEveryStatIsReadByTheGame) { // (#919)
    ASSERT_TRUE(Skills_Load("resources/data/skills.txt"));
    // The stats the game reads (terraintool, skillplay, buildtool, farming, suit, survival): a typo in skills.txt would
    // be a perk that does nothing.
    static const std::set<std::string> kRead = { "rock_speed", "ore_speed", "rock_radius", "ore_yield", "stone_yield", "gem_chance", "pick_wear",
        "soil_speed", "soil_radius", "shovel_wear", "chop_speed", "log_yield", "axe_wear", "fibre_yield", "stone_find", "crop_yield", "seed_yield",
        "break_speed", "walk_speed", "fall_safe", "hunger", "cold" };
    std::set<std::string> withTree;
    for (const PerkDef& p : Skills_Perks()) {
        withTree.insert(p.skill);
        for (const auto& [stat, per] : p.effects) EXPECT_TRUE(kRead.count(stat)) << p.id << ": " << stat << " is read nowhere";
        for (const std::string& q : p.needs) { const PerkDef* n = Skills_FindPerk(q); if (n) EXPECT_EQ(n->skill, p.skill) << p.id << " needs " << q << " of another skill"; }
    }
    for (const SkillDef& s : Skills_All()) EXPECT_TRUE(withTree.count(s.id)) << s.id << " has no perks";
}

TEST(SkillsData, TheNewTreesStackRankByRank) { // (#919) about +10% a rank, added across perks (then times the suit's)
    ASSERT_TRUE(Skills_Load("resources/data/skills.txt"));
    SkillState s;
    EXPECT_DOUBLE_EQ(Skills_Stat(s, "walk_speed"), 1.0) << "no perks: no change";
    s.ranks["long_stride"] = 3; s.ranks["marathon"] = 1;
    EXPECT_NEAR(Skills_Stat(s, "walk_speed"), 1.25, 1e-9) << "Long Stride 3 and Marathon";
    s.ranks["iron_stomach"] = 3; s.ranks["hardy"] = 1;
    EXPECT_LT(Skills_Stat(s, "hunger"), 0.8) << "hunger perks slow it";
    EXPECT_GT(Skills_Stat(s, "hunger"), 0.5);
}
#endif
