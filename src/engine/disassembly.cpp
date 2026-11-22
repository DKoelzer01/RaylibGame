#include "disassembly.h"
#include "composed.h"
#include "itemdata.h"

#include <algorithm>
#include <climits>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <unordered_map>

namespace {

struct PartLine { int ladderTier = 0; DisassemblyPart role[3]; bool has[3] = {}; }; // head, handle, binding

int lvLadder = 4;
int handTicks = 100, handEu = 8;
std::map<std::string, std::string> wildcards;
std::set<std::string> excluded;
std::vector<std::string> excludedSuffixes;
std::map<std::string, PartLine> partLines;
std::unordered_map<std::string, DisassemblyEntry> entries;
std::unordered_map<std::string, int> tiers; // voltage tier per item, memoised by Disassembly_Build
bool loaded = false;

std::string Trim(const std::string& s) {
    const size_t a = s.find_first_not_of(" \t"), b = s.find_last_not_of(" \t\r\n");
    return a == std::string::npos ? std::string() : s.substr(a, b - a + 1);
}
std::vector<std::string> Columns(const std::string& line) {
    std::vector<std::string> cols;
    std::stringstream ss(line);
    std::string c;
    while (std::getline(ss, c, '|')) cols.push_back(Trim(c));
    return cols;
}
std::string Resolve(const std::string& item) { auto it = wildcards.find(item); return it != wildcards.end() ? it->second : item; }

} // namespace

bool Disassembly_LoadText(const std::string& text) {
    wildcards.clear(); excluded.clear(); excludedSuffixes.clear(); partLines.clear();
    std::stringstream in(text);
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        const std::vector<std::string> c = Columns(line);
        if (c.size() < 2) continue;
        if (c[0] == "lv_ladder_tier") lvLadder = atoi(c[1].c_str());
        else if (c[0] == "hand" && c.size() >= 3) { handTicks = std::max(1, atoi(c[1].c_str())); handEu = std::max(0, atoi(c[2].c_str())); }
        else if (c[0] == "wildcard" && c.size() >= 3) wildcards[c[1]] = c[2];
        else if (c[0] == "exclude") excluded.insert(c[1]);
        else if (c[0] == "exclude_suffix") excludedSuffixes.push_back(c[1]);
        else if (c[0] == "part" && c.size() >= 6) {
            PartLine& p = partLines[c[1]];
            p.ladderTier = atoi(c[2].c_str());
            for (int r = 0; r < 3; ++r) {
                std::stringstream ps(c[3 + r]);
                std::string item; int n = 0;
                if (ps >> item >> n && item != "-" && n > 0) { p.role[r] = { item, n, true }; p.has[r] = true; }
            }
        }
    }
    loaded = true;
    return true;
}

bool Disassembly_Load(const char* path) {
    std::ifstream f(path);
    if (!f) return false;
    std::stringstream ss;
    ss << f.rdbuf();
    return Disassembly_LoadText(ss.str());
}

bool Disassembly_Excluded(const std::string& item) {
    if (excluded.count(item)) return true;
    for (const std::string& s : excludedSuffixes)
        if (item.size() > s.size() && item.compare(item.size() - s.size(), s.size(), s) == 0) return true;
    return gDisassemblyExclude && gDisassemblyExclude(item);
}

void Disassembly_Build(const std::vector<DisassemblyRecipe>& recipes, const std::function<int(const std::string&)>& ladderTier) {
    if (!loaded) Disassembly_Load();
    entries.clear(); tiers.clear();
    // Every way of making an item, wildcards resolved.
    std::map<std::string, std::vector<const DisassemblyRecipe*>> ways;
    for (const DisassemblyRecipe& r : recipes) if (!r.output.empty() && r.outCount > 0 && !r.in.empty()) ways[r.output].push_back(&r);
    for (const auto& [item, rs] : ways) {
        if (Disassembly_Excluded(item)) continue;
        // The first recipe gives the shape (its count, time and EU). An item made in several ways gives back only what every
        // way surely used, the least per item made (torches: the planks, not the coal or the charcoal one of them took), so
        // making it the cheap way and taking it apart never turns one material into another.
        const DisassemblyRecipe& first = *rs[0];
        DisassemblyEntry e;
        e.outCount = first.outCount;
        e.ticks = first.ticks > 0 ? first.ticks : handTicks;
        e.eu = first.ticks > 0 ? first.eu : handEu;
        std::map<std::string, int> need; // per craft of `first`
        for (const auto& [in, n] : first.in) need[Resolve(in)] += n;
        for (size_t k = 1; k < rs.size(); ++k) {
            std::map<std::string, int> other;
            for (const auto& [in, n] : rs[k]->in) other[Resolve(in)] += n;
            for (auto& [in, n] : need) {
                auto it = other.find(in);
                // that much per item made: n / first.outCount against m / outCount, in first's crafts (rounded down)
                const long long m = it == other.end() ? 0 : (long long)it->second * first.outCount / rs[k]->outCount;
                n = (int)std::min<long long>(n, m);
            }
        }
        for (const auto& [in, n] : need) if (n > 0 && in != item) e.parts.push_back({ in, n, true });
        if (e.parts.empty()) continue;
        entries[item] = std::move(e);
    }
    for (auto& [item, e] : entries) for (DisassemblyPart& p : e.parts) p.raw = !ways.count(p.item);
    // Tiers: the unlock's, else the ingredients' highest (a cycle counts as the lowest).
    const int lowest = -lvLadder;
    std::set<std::string> visiting;
    std::function<int(const std::string&)> tierOf = [&](const std::string& item) -> int {
        auto hit = tiers.find(item);
        if (hit != tiers.end()) return hit->second;
        int t = ladderTier ? ladderTier(item) : -1;
        if (t >= 0) t -= lvLadder;
        else {
            t = lowest;
            auto w = ways.find(item);
            if (w != ways.end() && visiting.insert(item).second) {
                for (const auto& [in, n] : w->second[0]->in) t = std::max(t, tierOf(Resolve(in)));
                visiting.erase(item);
            }
        }
        tiers[item] = t;
        return t;
    };
    for (auto& [item, e] : entries) e.tier = tierOf(item);
    for (const auto& [item, rs] : ways) tierOf(item);
}

const DisassemblyEntry* Disassembly_Find(const std::string& item) {
    auto it = entries.find(item);
    return it != entries.end() ? &it->second : nullptr;
}

std::vector<std::string> Disassembly_Items() {
    std::vector<std::string> v;
    for (const auto& [item, e] : entries) v.push_back(item);
    std::sort(v.begin(), v.end());
    return v;
}

int Disassembly_ItemTier(const std::string& item) {
    auto it = tiers.find(item);
    return it != tiers.end() ? it->second : -lvLadder;
}

bool Disassembly_Composed(const std::string& cls, const std::string& data, DisassemblyEntry& out) {
    if (!loaded) Disassembly_Load();
    const ComposedStats st = Composed_Stats(cls, data);
    if (!st.valid) return false;
    // Its full record: the default handle and binding filled in, as the item was made
    const std::string full = Composed_MakeData(st.head, st.handle, st.binding);
    const std::string mats[3] = { st.head, ItemData::Get(full, "g"), ItemData::Get(full, "b") };
    out = DisassemblyEntry{};
    out.ticks = handTicks; out.eu = handEu; out.outCount = 1;
    int ladder = 0;
    for (int r = 0; r < 3; ++r) {
        auto it = partLines.find(mats[r]);
        if (it == partLines.end()) continue;
        ladder = std::max(ladder, it->second.ladderTier);
        if (!it->second.has[r]) continue;
        const DisassemblyPart& p = it->second.role[r];
        auto same = std::find_if(out.parts.begin(), out.parts.end(), [&](const DisassemblyPart& q) { return q.item == p.item; });
        if (same != out.parts.end()) same->count += p.count; else out.parts.push_back(p);
    }
    out.tier = ladder - lvLadder;
    return !out.parts.empty();
}

uint32_t Disassembly_Next(uint32_t& rng) {
    if (rng == 0) rng = 0x9e3779b9u;
    rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5;
    return rng;
}

void Disassembly_Give(const DisassemblyEntry& e, int yieldPct, uint32_t& rng, std::vector<DisassemblyGive>& out) {
    out.clear();
    if (yieldPct <= 0) return;
    const int per = 100 * std::max(1, e.outCount);
    for (const DisassemblyPart& p : e.parts) {
        int centi = 0;
        if (p.raw) centi = (int)((long long)p.count * yieldPct / std::max(1, e.outCount)); // hundredths: count x yield% / outCount
        else for (int u = 0; u < p.count; ++u) if ((int)(Disassembly_Next(rng) % (uint32_t)per) < yieldPct) centi += 100;
        if (centi > 0) out.push_back({ p.item, centi });
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
#include "machinetiers.h"

namespace {
const char* kTestData =
    "lv_ladder_tier | 4\n"
    "hand | 100 | 8\n"
    "wildcard | any_planks | oak_planks\n"
    "exclude_suffix | _pack\n"
    "part | wood | 0 | oak_planks 3 | oak_planks 2 | -\n"
    "part | fibre | 0 | - | - | fibre 2\n"
    "part | iron | 2 | iron_ingot 3 | iron_ingot 1 | iron_ingot 1\n"
    "part | steel | 3 | steel_ingot 3 | steel_ingot 1 | steel_ingot 1\n";
std::vector<DisassemblyRecipe> TestRecipes() {
    return {
        { "copper_wire", 2, { { "copper_ingot", 1 } }, 0, 0 },
        { "basic_circuit", 1, { { "vacuum_tube", 2 }, { "copper_wire", 4 }, { "steel_plate", 1 } }, 160, 16 },
        { "vacuum_tube", 1, { { "glass", 1 }, { "copper_wire", 1 } }, 0, 0 },
        { "macerator_mv", 1, { { "mv_hull", 1 }, { "basic_circuit", 2 } }, 0, 0 },
        { "mv_hull", 1, { { "steel_plate", 8 } }, 0, 0 },
        { "torch", 4, { { "any_planks", 1 }, { "coal", 1 } }, 0, 0 },
        { "torch", 4, { { "any_planks", 1 }, { "charcoal", 1 } }, 0, 0 },
        { "lv_pack", 2, { { "steel_plate", 1 } }, 200, 12 },
    };
}
int TestTier(const std::string& item) { return item == "basic_circuit" ? 4 : item == "macerator_mv" ? 5 : -1; }
}

TEST(Disassembly, OneRecipeLevelDownWithRawAndComponentInputs) {
    ASSERT_TRUE(Disassembly_LoadText(kTestData));
    Disassembly_Build(TestRecipes(), TestTier);
    const DisassemblyEntry* c = Disassembly_Find("basic_circuit");
    ASSERT_NE(c, nullptr);
    EXPECT_EQ(c->tier, 0) << "unlocked at ladder tier 4: LV";
    EXPECT_EQ(c->ticks, 160); EXPECT_EQ(c->eu, 16) << "the forward recipe's time and EU";
    ASSERT_EQ(c->parts.size(), 3u);
    for (const DisassemblyPart& p : c->parts) {
        if (p.item == "steel_plate") { EXPECT_TRUE(p.raw); EXPECT_EQ(p.count, 1); }
        else { EXPECT_FALSE(p.raw) << p.item << " is made: a component"; }
    }
    const DisassemblyEntry* m = Disassembly_Find("macerator_mv");
    ASSERT_NE(m, nullptr);
    EXPECT_EQ(m->tier, 1);
    EXPECT_EQ(m->ticks, 100) << "a hand recipe: the file's hand time";
    ASSERT_EQ(m->parts.size(), 2u) << "its direct inputs, not the steel plates under the hull";
    EXPECT_EQ(Disassembly_ItemTier("mv_hull"), -4) << "no unlock, raw inputs: the lowest";
    EXPECT_EQ(Disassembly_Find("lv_pack"), nullptr) << "packs are never taken";
    EXPECT_EQ(Disassembly_Find("steel_plate"), nullptr) << "raw: nothing to take apart";
    const DisassemblyEntry* t = Disassembly_Find("torch");
    ASSERT_NE(t, nullptr);
    ASSERT_EQ(t->parts.size(), 1u) << "two ways to make it: only the planks both used";
    EXPECT_EQ(t->parts[0].item, "oak_planks") << "a wildcard comes back as the file's item";
}

TEST(Disassembly, EachTierGapsYieldAndRawFractions) {
    ASSERT_TRUE(Disassembly_LoadText(kTestData));
    Disassembly_Build(TestRecipes(), TestTier);
    const DisassemblyEntry& c = *Disassembly_Find("basic_circuit"); // LV (0)
    const int want[6][2] = { { -3, 0 }, { -2, 1 }, { -1, 10 }, { 0, 50 }, { 1, 75 }, { 2, 100 } }; // assembler tier: yield
    for (const auto& w : want) {
        const int y = MachineTiers_YieldPercent(w[0], c.tier);
        EXPECT_EQ(y, w[1]) << "assembler tier " << w[0];
        uint32_t rng = 1;
        std::vector<DisassemblyGive> g;
        Disassembly_Give(c, y, rng, g);
        int steel = 0;
        for (const DisassemblyGive& x : g) if (x.item == "steel_plate") steel = x.centi;
        EXPECT_EQ(steel, w[1]) << "1 steel plate at " << w[1] << "%: that many hundredths";
    }
    // A torch (4 a craft) at 50%: 1 plank / 4 x 50% = 12 hundredths, rounded down
    uint32_t rng = 1;
    std::vector<DisassemblyGive> g;
    Disassembly_Give(*Disassembly_Find("torch"), 50, rng, g);
    ASSERT_EQ(g.size(), 1u);
    EXPECT_EQ(g[0].centi, 12);
}

TEST(Disassembly, ChanceRollsAreSeededAndNearTheYield) {
    ASSERT_TRUE(Disassembly_LoadText(kTestData));
    Disassembly_Build(TestRecipes(), TestTier);
    const DisassemblyEntry& c = *Disassembly_Find("basic_circuit");
    auto run = [&](uint32_t seed, int n) {
        uint32_t rng = seed; long long wires = 0;
        std::vector<DisassemblyGive> g;
        for (int i = 0; i < n; ++i) { Disassembly_Give(c, 75, rng, g); for (const DisassemblyGive& x : g) if (x.item == "copper_wire") wires += x.centi / 100; }
        return wires;
    };
    EXPECT_EQ(run(12345, 50), run(12345, 50)) << "the same seed, the same rolls";
    const long long w = run(777, 2000); // 4 wires a circuit, 75% each: 6000 expected
    EXPECT_NEAR((double)w, 6000.0, 200.0);
    uint32_t rng = 5; std::vector<DisassemblyGive> g;
    Disassembly_Give(c, 100, rng, g);
    int all = 0; for (const DisassemblyGive& x : g) if (x.item == "copper_wire") all = x.centi;
    EXPECT_EQ(all, 400) << "100%: every unit";
}

TEST(Disassembly, CraftThenScrapNeverGains) {
    ASSERT_TRUE(Disassembly_LoadText(kTestData));
    const std::vector<DisassemblyRecipe> rs = TestRecipes();
    Disassembly_Build(rs, TestTier);
    // For every item and every way of making it, at the best yield (100%): what comes back per item is no more than that
    // way used of the same item per item made.
    for (const DisassemblyRecipe& r : rs) {
        const DisassemblyEntry* e = Disassembly_Find(r.output);
        if (!e) continue;
        for (const DisassemblyPart& p : e->parts) {
            int used = 0;
            for (const auto& [in, n] : r.in) if (in == p.item || (in == "any_planks" && p.item == "oak_planks")) used += n;
            // back per item: p.count / e->outCount; used per item: used / r.outCount
            EXPECT_LE((long long)p.count * r.outCount, (long long)used * e->outCount) << r.output << ": " << p.item;
        }
    }
}

TEST(Disassembly, AComposedToolGivesItsPartsRawMaterials) {
    ASSERT_TRUE(Composed_Load());
    ASSERT_TRUE(Disassembly_LoadText(kTestData));
    DisassemblyEntry e;
    ASSERT_TRUE(Disassembly_Composed("tool_pickaxe", Composed_MakeData("iron"), e));
    int ingots = 0, planks = 0, fibre = 0;
    for (const DisassemblyPart& p : e.parts) {
        EXPECT_TRUE(p.raw);
        if (p.item == "iron_ingot") ingots = p.count; else if (p.item == "oak_planks") planks = p.count; else if (p.item == "fibre") fibre = p.count;
    }
    EXPECT_EQ(ingots, 3) << "an iron head: 3 ingots";
    EXPECT_EQ(planks, 2) << "the default wood handle";
    EXPECT_EQ(fibre, 2) << "the default fibre binding (iron's)";
    EXPECT_EQ(e.tier, 2 - 4) << "iron parts: ladder tier 2, two under LV";
    ASSERT_TRUE(Disassembly_Composed("tool_pickaxe", Composed_MakeData("steel"), e));
    EXPECT_EQ(e.tier, -1) << "steel: ladder tier 3";
    EXPECT_FALSE(Disassembly_Composed("iron_ingot", "", e)) << "not a composed item";
}
#endif
