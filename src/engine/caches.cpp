#include "caches.h"
#include "inventory.h"

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>

// (#2072) The cache reward pool (resources/data/caches.txt) and the roll.
namespace {
struct KindDef { int rolls = 1; int weight[4] = { 1, 0, 0, 0 }; }; // weights by CacheReward
struct PartEntry { std::string item; int lo = 1, hi = 1; };
struct Pool {
    bool loaded = false;
    KindDef kinds[3];
    int rarity[3] = { 70, 25, 5 };
    std::map<std::pair<int, int>, std::vector<PartEntry>> parts; // (tier, rarity) -> entries
    int maxTier = 0;
    std::vector<CacheBlueprint> blueprints;
    std::vector<CacheAlt> alts;
    int boostShare = 25, boostMin = 5;
};
Pool pool;
std::map<int, std::set<std::string>> teamAlts;

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t\r"), b = s.find_last_not_of(" \t\r");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}
std::vector<std::string> Fields(const std::string& line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    for (std::string f; std::getline(ss, f, '|');) out.push_back(Trim(f));
    return out;
}
// "word value" -> value, for "rolls 2", "tier 7", "share 25"
bool Keyed(const std::string& f, const char* key, int& v) {
    std::istringstream in(f); std::string k;
    return (in >> k >> v) && k == key;
}

// A splitmix64 stream: the roll's only source of chance.
struct Rng {
    uint64_t s;
    uint64_t next() { uint64_t z = (s += 0x9E3779B97F4A7C15ull); z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull; z = (z ^ (z >> 27)) * 0x94D049BB133111EBull; return z ^ (z >> 31); }
    int below(int n) { return n <= 0 ? 0 : (int)(next() % (uint64_t)n); }
    int pick(const int* w, int n) {
        int total = 0;
        for (int i = 0; i < n; ++i) total += std::max(0, w[i]);
        int r = below(total);
        for (int i = 0; i < n; ++i) { if (r < std::max(0, w[i])) return i; r -= std::max(0, w[i]); }
        return 0;
    }
};

std::string Shown(const std::string& item) {
    if (const ItemDef* d = itemDB.find(item)) return d->name;
    std::string s = item; std::replace(s.begin(), s.end(), '_', ' ');
    return s;
}

void EnsureLoaded() { if (!pool.loaded) Caches_Load(); }

// A parts reward of `rarity`, from tier `tier` (or the nearest lower tier with a line).
CacheRoll RollParts(Rng& rng, int tier, CacheRarity rarity) {
    CacheRoll r; r.kind = CacheReward::Parts; r.rarity = rarity;
    tier = std::min(std::max(tier, 0), pool.maxTier);
    const std::vector<PartEntry>* line = nullptr;
    for (int t = tier; t >= 0 && !line; --t) {
        auto it = pool.parts.find({ t, (int)rarity });
        if (it != pool.parts.end() && !it->second.empty()) line = &it->second;
    }
    if (!line) return r;
    const PartEntry& e = (*line)[(size_t)rng.below((int)line->size())];
    const int n = e.lo + rng.below(e.hi - e.lo + 1);
    r.items = { { e.item, n } };
    r.text = n == 1 ? Shown(e.item) : std::to_string(n) + " " + Shown(e.item);
    return r;
}
} // namespace

bool Caches_LoadText(const std::string& text, std::string& error) {
    Pool p;
    error.clear();
    std::istringstream in(text);
    int lineNo = 0;
    auto bad = [&](const std::string& why) { error = "caches.txt:" + std::to_string(lineNo) + ": " + why; return false; };
    for (std::string line; std::getline(in, line);) {
        ++lineNo;
        line = Trim(line);
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> f = Fields(line);
        std::istringstream head(f[0]);
        std::string word; head >> word;
        if (word == "cache") {
            std::string kind; head >> kind;
            const int k = kind == "small" ? 0 : kind == "vault" ? 1 : kind == "deep_vault" ? 2 : -1;
            if (k < 0) return bad("unknown cache kind '" + kind + "'");
            KindDef d; d.weight[0] = 0;
            for (size_t i = 1; i < f.size(); ++i) {
                int v = 0;
                if (Keyed(f[i], "rolls", v)) d.rolls = v;
                else if (Keyed(f[i], "parts", v)) d.weight[(int)CacheReward::Parts] = v;
                else if (Keyed(f[i], "blueprint", v)) d.weight[(int)CacheReward::Blueprint] = v;
                else if (Keyed(f[i], "alt", v)) d.weight[(int)CacheReward::Alt] = v;
                else if (Keyed(f[i], "boost", v)) d.weight[(int)CacheReward::Boost] = v;
                else return bad("unknown field '" + f[i] + "'");
            }
            if (d.rolls < 1 || d.rolls > 8 || d.weight[0] <= 0) return bad("a cache needs 1-8 rolls and a parts weight");
            p.kinds[k] = d;
        } else if (word == "rarity") {
            std::string name; int w = -1; head >> name >> w;
            const int r = name == "common" ? 0 : name == "uncommon" ? 1 : name == "rare" ? 2 : -1;
            if (r < 0 || w < 0) return bad("rarity <common|uncommon|rare> <weight>");
            p.rarity[r] = w;
        } else if (word == "parts") {
            int tier = -1; std::string name; head >> tier >> name;
            const int r = name == "common" ? 0 : name == "uncommon" ? 1 : name == "rare" ? 2 : -1;
            if (tier < 0 || tier > 20 || r < 0) return bad("parts <tier> <rarity> | item min-max | ...");
            std::vector<PartEntry>& out = p.parts[{ tier, r }];
            for (size_t i = 1; i < f.size(); ++i) {
                std::istringstream e(f[i]); PartEntry pe; std::string range; char dash = 0;
                if (!(e >> pe.item >> range)) return bad("an entry is 'item min-max'");
                std::istringstream rs(range);
                if (!(rs >> pe.lo >> dash >> pe.hi) || dash != '-' || pe.lo < 1 || pe.hi < pe.lo) return bad("bad count '" + range + "'");
                out.push_back(pe);
            }
            p.maxTier = std::max(p.maxTier, tier);
        } else if (word == "blueprint" || word == "alt") {
            std::string id; head >> id;
            int tier = -1;
            if (id.empty() || f.size() != 3 || !Keyed(f[1], "tier", tier) || tier < 0 || f[2].empty()) return bad(word + " <id> | tier <n> | <shown name>");
            if (word == "blueprint") p.blueprints.push_back({ id, f[2], tier }); else p.alts.push_back({ id, f[2], tier });
        } else if (word == "boost") {
            for (size_t i = 1; i < f.size(); ++i) {
                int v = 0;
                if (Keyed(f[i], "share", v)) p.boostShare = v;
                else if (Keyed(f[i], "min", v)) p.boostMin = v;
                else return bad("unknown field '" + f[i] + "'");
            }
        } else {
            return bad("unknown line '" + word + "'");
        }
    }
    if (p.parts.empty()) { lineNo = 0; return bad("no parts lines"); }
    p.loaded = true;
    pool = std::move(p);
    return true;
}

bool Caches_Load(const char* path) {
    std::ifstream f(path, std::ios::binary);
    std::string error;
    if (!f) { printf("caches: can't read %s\n", path); pool.loaded = true; return false; } // (an empty pool: parts rolls give nothing)
    std::stringstream ss; ss << f.rdbuf();
    if (!Caches_LoadText(ss.str(), error)) { printf("%s\n", error.c_str()); pool.loaded = true; return false; }
    return true;
}

const std::vector<CacheBlueprint>& Caches_Blueprints() { EnsureLoaded(); return pool.blueprints; }
const std::vector<CacheAlt>& Caches_Alts() { EnsureLoaded(); return pool.alts; }
const CacheAlt* Caches_FindAlt(const std::string& id) {
    for (const CacheAlt& a : Caches_Alts()) if (a.id == id) return &a;
    return nullptr;
}
std::vector<std::string> Caches_PartItems() {
    EnsureLoaded();
    std::vector<std::string> out;
    for (const auto& [key, entries] : pool.parts) for (const PartEntry& e : entries) out.push_back(e.item);
    return out;
}

std::vector<CacheRoll> Caches_Roll(uint64_t worldSeed, uint64_t cacheId, CacheKind kind, int band, const CacheTeam& team) {
    EnsureLoaded();
    // The team's progress is its tier and how much it has researched: the same team at the same point rolls the same.
    Rng rng{ worldSeed * 0x9E3779B97F4A7C15ull ^ cacheId };
    rng.s ^= (uint64_t)team.tier * 0xD6E8FEB86659FD93ull ^ (uint64_t)team.researched.size() * 0xA0761D6478BD642Full;
    rng.next();
    const KindDef& def = pool.kinds[std::min((int)kind, 2)];
    const int deep = band >= 2 ? 1 : 0;
    const int T = team.tier + deep;
    std::set<std::string> gotBp(team.blueprints.begin(), team.blueprints.end()), gotAlt(team.alts.begin(), team.alts.end());
    std::set<std::string> boosted; // one boost per node per cache
    std::vector<CacheRoll> out;
    for (int slot = 0; slot < def.rolls; ++slot) {
        CacheRoll got;
        bool done = false;
        for (int attempt = 0; attempt < 4 && !done; ++attempt) { // a category with nothing left re-rolls, up to 3 times
            const CacheReward cat = (CacheReward)rng.pick(def.weight, 4);
            if (cat == CacheReward::Parts) {
                const CacheRarity rar = (CacheRarity)rng.pick(pool.rarity, 3);
                int tier = T;
                if (rar == CacheRarity::Common && rng.below(100) < 30) tier = T - 1;
                if (rar == CacheRarity::Rare && band >= 1 && !deep) tier = T + 1;
                got = RollParts(rng, tier, rar);
                done = !got.items.empty();
            } else if (cat == CacheReward::Blueprint || cat == CacheReward::Alt) {
                std::vector<std::pair<std::string, std::string>> can; // id, shown name
                if (cat == CacheReward::Blueprint) { for (const CacheBlueprint& b : pool.blueprints) if (b.tier <= T && !gotBp.count(b.id)) can.push_back({ b.id, b.name }); }
                else for (const CacheAlt& a : pool.alts) if (a.tier <= T && !gotAlt.count(a.id)) can.push_back({ a.id, a.name });
                if (can.empty()) continue;
                const auto& [id, name] = can[(size_t)rng.below((int)can.size())];
                got = {}; got.kind = cat; got.rarity = CacheRarity::Rare; got.id = id;
                got.text = (cat == CacheReward::Blueprint ? "Blueprint: " : "Alternate recipe: ") + name;
                (cat == CacheReward::Blueprint ? gotBp : gotAlt).insert(id);
                done = true;
            } else { // Boost
                std::vector<const CacheTeam::Node*> can;
                std::vector<int> w;
                for (const CacheTeam::Node& n : team.frontier)
                    if (n.tier <= T && !boosted.count(n.id)) { can.push_back(&n); w.push_back(n.core && n.tier == T ? 2 : 1); }
                if (can.empty()) continue;
                const CacheTeam::Node& n = *can[(size_t)rng.pick(w.data(), (int)w.size())];
                got = {}; got.kind = CacheReward::Boost; got.rarity = CacheRarity::Uncommon; got.id = n.id;
                got.units = std::min(n.units, std::max(pool.boostMin, n.units * pool.boostShare / 100));
                got.text = "Research: " + std::to_string(got.units) + " units of " + n.name;
                boosted.insert(n.id);
                done = true;
            }
        }
        if (!done) { got = RollParts(rng, T, CacheRarity::Common); if (got.items.empty()) continue; }
        out.push_back(got);
    }
    return out;
}

// ---- A team's alternate recipes ----
const std::set<std::string>& Caches_TeamAlts(int team) { return teamAlts[team]; }
bool Caches_TeamHasAlt(int team, const std::string& id) { auto it = teamAlts.find(team); return it != teamAlts.end() && it->second.count(id) > 0; }
void Caches_AddTeamAlt(int team, const std::string& id) { teamAlts[team].insert(id); }
std::string Caches_TeamAltsText(int team) {
    std::string out;
    auto it = teamAlts.find(team);
    if (it != teamAlts.end()) for (const std::string& id : it->second) out += id + "\n";
    return out;
}
void Caches_SetTeamAltsText(int team, const std::string& text) {
    std::set<std::string>& s = teamAlts[team];
    s.clear();
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) { line = Trim(line); if (!line.empty()) s.insert(line); }
}
void Caches_ResetTeams() { teamAlts.clear(); }

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>

namespace {
CacheTeam TeamAt(int tier) {
    CacheTeam t; t.tier = tier;
    t.frontier = { { "a", "Node A", tier, 40, true }, { "b", "Node B", tier, 10, false }, { "c", "Node C", tier + 1, 40, true } };
    return t;
}
}

// (#2072) The shipped pool loads, and every item it names is a real item.
TEST(Caches, ThePoolLoadsAndNamesRealItems) {
    ASSERT_TRUE(itemDB.load("resources/data/items.txt", "resources/data/recipes.txt"));
    std::ifstream f("resources/data/caches.txt", std::ios::binary);
    std::stringstream ss; ss << f.rdbuf();
    std::string error;
    ASSERT_TRUE(Caches_LoadText(ss.str(), error)) << error;
    for (const std::string& id : Caches_PartItems())
        EXPECT_TRUE(IsWildcardItem(id) || itemDB.find(id) != nullptr) << "caches.txt names an unknown item: " << id;
    for (const CacheBlueprint& b : Caches_Blueprints())
        EXPECT_TRUE(std::ifstream("resources/blueprints/" + b.id + ".bp").good()) << "no preset for blueprint " << b.id;
}

// (#2072) The same seed, cache and progress roll the same; more progress or another cache rolls differently; a small cache
// holds its 2 rewards and a deep vault its 4.
TEST(Caches, TheRollIsDeterministicAndSizedByKind) {
    std::string error;
    ASSERT_TRUE(Caches_LoadText(
        "cache small | rolls 2 | parts 75 | blueprint 10 | alt 5 | boost 10\n"
        "cache vault | rolls 3 | parts 60 | blueprint 15 | alt 10 | boost 15\n"
        "cache deep_vault | rolls 4 | parts 50 | blueprint 15 | alt 20 | boost 15\n"
        "parts 3 common | steel_plate 4-10 | coke 4-10\nparts 3 uncommon | vacuum_tube 1-3\nparts 3 rare | boiler 1-1\n"
        "parts 4 common | cable 8-16\nparts 4 uncommon | lv_hull 1-1\nparts 4 rare | macerator 1-1\n"
        "boost | share 25 | min 5\n", error)) << error;
    const CacheTeam t = TeamAt(3);
    auto same = [](const std::vector<CacheRoll>& a, const std::vector<CacheRoll>& b) {
        if (a.size() != b.size()) return false;
        for (size_t i = 0; i < a.size(); ++i) if (a[i].text != b[i].text) return false;
        return true;
    };
    const auto a = Caches_Roll(42, 7, CacheKind::Small, 0, t);
    EXPECT_EQ(a.size(), 2u);
    EXPECT_TRUE(same(a, Caches_Roll(42, 7, CacheKind::Small, 0, t)));
    EXPECT_EQ(Caches_Roll(42, 7, CacheKind::DeepVault, 0, t).size(), 4u);
    int differ = 0;
    for (uint64_t id = 100; id < 120; ++id) differ += same(Caches_Roll(42, id, CacheKind::Vault, 0, t), Caches_Roll(42, id + 1000, CacheKind::Vault, 0, t)) ? 0 : 1;
    EXPECT_GT(differ, 10) << "other caches roll other rewards";
    CacheTeam more = t; more.researched = { "x" };
    differ = 0;
    for (uint64_t id = 100; id < 120; ++id) differ += same(Caches_Roll(42, id, CacheKind::Vault, 0, t), Caches_Roll(42, id, CacheKind::Vault, 0, more)) ? 0 : 1;
    EXPECT_GT(differ, 10) << "the same cache rolls anew for a team further on";
}

// (#2072) Parts come from the team's tier (common one below now and then); the deep band rolls a tier up; a boost names a
// node the team could start, at most its tier, with 25% of the node's units (at least 5); blueprints and alternate recipes
// the team has are not rolled again, and with none left the reward is parts.
TEST(Caches, RewardsFollowTheTeamsProgress) {
    std::string error;
    ASSERT_TRUE(Caches_LoadText(
        "cache small | rolls 2 | parts 1 | blueprint 1 | alt 1 | boost 1\n"
        "cache deep_vault | rolls 4 | parts 1\n"
        "parts 2 common | iron_plate 1-1\nparts 3 common | steel_plate 1-1\nparts 3 uncommon | vacuum_tube 1-1\nparts 3 rare | boiler 1-1\n"
        "parts 4 common | cable 1-1\nparts 4 uncommon | lv_hull 1-1\nparts 4 rare | macerator 1-1\n"
        "blueprint bp3 | tier 3 | Steam block\nblueprint bp5 | tier 5 | MV line\nalt alt3 | tier 3 | Packed coke\n"
        "boost | share 25 | min 5\n", error)) << error;
    const CacheTeam t = TeamAt(3);
    std::map<std::string, int> items; int bp3 = 0, alt3 = 0, boosts = 0;
    for (uint64_t id = 0; id < 400; ++id)
        for (const CacheRoll& r : Caches_Roll(9, id, CacheKind::Small, 0, t)) {
            if (r.kind == CacheReward::Parts) { ASSERT_EQ(r.items.size(), 1u); items[r.items[0].first]++; }
            if (r.kind == CacheReward::Blueprint) { EXPECT_EQ(r.id, "bp3") << "no blueprint above the team's tier"; bp3++; }
            if (r.kind == CacheReward::Alt) { EXPECT_EQ(r.id, "alt3"); alt3++; }
            if (r.kind == CacheReward::Boost) {
                boosts++;
                EXPECT_NE(r.id, "c") << "no boost above the team's tier";
                EXPECT_EQ(r.units, r.id == "a" ? 10 : 5) << "25% of 40, or the minimum 5";
            }
        }
    EXPECT_GT(items["steel_plate"], 0); EXPECT_GT(items["iron_plate"], 0) << "common: sometimes a tier below";
    EXPECT_EQ(items["cable"] + items["lv_hull"], 0) << "nothing above the tier in the shallows but rare";
    EXPECT_GT(bp3, 0); EXPECT_GT(alt3, 0); EXPECT_GT(boosts, 0);
    int deepUp = 0;
    for (uint64_t id = 0; id < 50; ++id)
        for (const CacheRoll& r : Caches_Roll(9, id, CacheKind::DeepVault, 2, t)) deepUp += (r.items[0].first == "cable" || r.items[0].first == "lv_hull" || r.items[0].first == "macerator") ? 1 : 0;
    EXPECT_GT(deepUp, 100) << "the deep band rolls a tier up";
    CacheTeam had = t; had.blueprints = { "bp3" }; had.alts = { "alt3" }; had.frontier.clear();
    for (uint64_t id = 0; id < 200; ++id)
        for (const CacheRoll& r : Caches_Roll(9, id, CacheKind::Small, 0, had))
            EXPECT_EQ(r.kind, CacheReward::Parts) << "nothing new to give but parts";
    Caches_Load(); // the shipped pool back for later tests
}

// (#2072) A team's alternate recipes round-trip through their store text.
TEST(Caches, TeamAltsRoundTrip) {
    Caches_ResetTeams();
    Caches_AddTeamAlt(3, "packed_coke"); Caches_AddTeamAlt(3, "dust_steel"); Caches_AddTeamAlt(3, "packed_coke");
    const std::string text = Caches_TeamAltsText(3);
    EXPECT_EQ(text, "dust_steel\npacked_coke\n");
    Caches_ResetTeams();
    EXPECT_FALSE(Caches_TeamHasAlt(3, "dust_steel"));
    Caches_SetTeamAltsText(3, text);
    EXPECT_TRUE(Caches_TeamHasAlt(3, "dust_steel") && Caches_TeamHasAlt(3, "packed_coke"));
    EXPECT_EQ(Caches_TeamAltsText(4), "");
    Caches_ResetTeams();
}
#endif
