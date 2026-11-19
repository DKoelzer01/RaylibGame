#include "craftplan.h"
#include <algorithm>

namespace {
constexpr int kMaxDepth = 8; // (the search is bounded: the deepest hand chain today is 3)

struct Planner {
    const std::vector<Recipe>& recipes;
    const std::function<bool(int)>& known;
    const std::set<std::string>& machineMade;
    std::map<std::string, int> have;
    CraftPlan plan;
    std::vector<int> visiting; // recipes on the current path (cycles: a recipe never needs itself)

    // Take what is carried of `want` (an id or a group): the exact item first.
    int Take(const std::string& want, int qty) {
        auto take = [&](const std::string& id) {
            auto it = have.find(id);
            if (it == have.end() || it->second <= 0) return;
            const int n = std::min(qty, it->second);
            it->second -= n; qty -= n;
        };
        take(want);
        for (auto& [id, n] : have) { if (qty <= 0) break; if (id != want && ItemMatches(want, id)) take(id); }
        return qty;
    }
    bool MachineMakes(const std::string& want) const {
        for (const std::string& id : machineMade) if (ItemMatches(want, id)) return true;
        return false;
    }
    // `qty` of `want` into the plan: carried first, then hand-crafted (recursively), else missing.
    bool Need(const std::string& want, int qty, int depth) {
        qty = Take(want, qty);
        if (qty <= 0) return true;
        std::vector<int> cands;
        if (depth < kMaxDepth)
            for (int k = 0; k < (int)recipes.size(); ++k)
                if (ItemMatches(want, recipes[k].output) && known(k) && std::find(visiting.begin(), visiting.end(), k) == visiting.end())
                    cands.push_back(k);
        std::stable_sort(cands.begin(), cands.end(), [&](int a, int b) { return (recipes[a].output == want) > (recipes[b].output == want); });
        for (int k : cands) {
            const Recipe& r = recipes[k];
            const int out = std::max(1, r.outputCount), runs = (qty + out - 1) / out;
            const auto have0 = have; const size_t steps0 = plan.steps.size(); const auto left0 = plan.leftover;
            visiting.push_back(k);
            bool ok = true;
            for (const auto& [in, n] : r.inputs) ok = Need(in, n * runs, depth + 1) && ok;
            visiting.pop_back();
            if (ok) {
                plan.steps.push_back({ k, runs });
                if (runs * out > qty) { have[r.output] += runs * out - qty; plan.leftover[r.output] += runs * out - qty; }
                return true;
            }
            have = have0; plan.steps.resize(steps0); plan.leftover = left0; // (its shortfalls stay: they say what is missing)
        }
        if (cands.empty()) {
            plan.missing[want] += qty;
            if (MachineMakes(want)) plan.machineOnly.insert(want);
        }
        return false;
    }
};
}

CraftPlan CraftPlan_Make(const std::vector<Recipe>& recipes, const std::function<bool(int)>& known,
                         const std::map<std::string, int>& have, int target, int times, const std::set<std::string>& machineMade) {
    Planner p{ recipes, known, machineMade, have, {}, {} };
    if (target < 0 || target >= (int)recipes.size() || times < 1 || !known(target)) return p.plan;
    const Recipe& r = recipes[target];
    p.visiting.push_back(target);
    bool ok = true;
    for (const auto& [in, n] : r.inputs) ok = p.Need(in, n * times, 1) && ok;
    if (ok) p.plan.steps.push_back({ target, times });
    // Neighbouring steps of one recipe run as one.
    std::vector<CraftStep> merged;
    for (const CraftStep& s : p.plan.steps) {
        if (!merged.empty() && merged.back().recipe == s.recipe) merged.back().times += s.times;
        else merged.push_back(s);
    }
    p.plan.steps = merged;
    p.plan.ok = ok;
    if (!ok) p.plan.steps.clear();
    return p.plan;
}

int CraftPlan_Max(const std::vector<Recipe>& recipes, const std::function<bool(int)>& known,
                  const std::map<std::string, int>& have, int target, int cap, const std::set<std::string>& machineMade) {
    if (!CraftPlan_Make(recipes, known, have, target, 1, machineMade).ok) return 0;
    int lo = 1, hi = std::max(1, cap); // lo: makeable
    while (lo < hi) {
        const int mid = lo + (hi - lo + 1) / 2;
        if (CraftPlan_Make(recipes, known, have, target, mid, machineMade).ok) lo = mid; else hi = mid - 1;
    }
    return lo;
}

std::string CraftPlan_Describe(const std::vector<Recipe>& recipes, const CraftPlan& plan) {
    std::string s;
    for (size_t i = 0; i < plan.steps.size(); ++i) {
        const Recipe& r = recipes[plan.steps[i].recipe];
        const ItemDef* d = itemDB.find(r.output);
        if (i) s += i + 1 == plan.steps.size() ? ", then " : ", ";
        s += std::to_string(plan.steps[i].times * std::max(1, r.outputCount)) + " " + (d ? d->name : r.output);
    }
    return s;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#928): the resolver on a small recipe set ----
#include <gtest/gtest.h>

namespace {
std::vector<Recipe> TestRecipes() {
    auto R = [](const char* out, int n, std::vector<std::pair<std::string, int>> in) { Recipe r; r.output = out; r.outputCount = n; r.inputs = std::move(in); return r; };
    return {
        R("copper_wire", 2, { { "copper_ingot", 1 } }),                                  // 0: wire in pairs
        R("circuit", 1, { { "copper_wire", 3 }, { "steel_plate", 1 }, { "glass", 1 } }),  // 1
        R("steel_plate", 1, { { "steel_ingot", 1 } }),                                  // 2
        R("computer", 1, { { "circuit", 2 }, { "steel_plate", 2 } }),                   // 3: two levels down
        R("planks", 4, { { "any_log", 1 } }),                                           // 4
        R("a", 1, { { "b", 1 } }), R("b", 1, { { "a", 1 } }),                           // 5, 6: a cycle
        R("lamp", 1, { { "glass", 1 }, { "iron_ingot", 1 } }),                          // 7: iron from a furnace only
    };
}
const std::function<bool(int)> kAllKnown = [](int) { return true; };
}

TEST(CraftPlan, DirectWhenTheIngredientsAreCarried) {
    const auto rs = TestRecipes();
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, { { "copper_ingot", 2 } }, 0, 2);
    ASSERT_TRUE(p.ok);
    ASSERT_EQ(p.steps.size(), 1u);
    EXPECT_EQ(p.steps[0].recipe, 0); EXPECT_EQ(p.steps[0].times, 2);
}

TEST(CraftPlan, OneLevelFromRawMaterialsWithLeftoverWire) {
    const auto rs = TestRecipes();
    // The user's example: copper, steel and glass, all raw; a circuit needs 3 wire = 2 crafts of 2 (one wire left over).
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, { { "copper_ingot", 5 }, { "steel_ingot", 1 }, { "glass", 1 } }, 1, 1);
    ASSERT_TRUE(p.ok);
    ASSERT_EQ(p.steps.size(), 3u);
    EXPECT_EQ(p.steps.back().recipe, 1) << "the circuit last";
    int wireRuns = 0, plateRuns = 0;
    for (const CraftStep& s : p.steps) { if (s.recipe == 0) wireRuns += s.times; if (s.recipe == 2) plateRuns += s.times; }
    EXPECT_EQ(wireRuns, 2); EXPECT_EQ(plateRuns, 1);
    EXPECT_EQ(p.leftover.at("copper_wire"), 1);
}

TEST(CraftPlan, TwoLevelsAndTheMostThatCanBeMade) {
    const auto rs = TestRecipes();
    // A computer: 2 circuits (6 wire = 3 ingots, 2 plates, 2 glass) and 2 plates: 3 copper, 4 steel, 2 glass.
    const std::map<std::string, int> have{ { "copper_ingot", 7 }, { "steel_ingot", 9 }, { "glass", 4 } };
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, have, 3, 1);
    ASSERT_TRUE(p.ok);
    EXPECT_EQ(p.steps.back().recipe, 3);
    EXPECT_EQ(CraftPlan_Max(rs, kAllKnown, have, 3), 2) << "glass for 2 (copper: 6 wire a computer, 12 for 2 = 6 ingots)";
}

TEST(CraftPlan, CarriedIntermediatesAreUsedFirst) {
    const auto rs = TestRecipes();
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, { { "copper_wire", 3 }, { "steel_plate", 1 }, { "glass", 1 } }, 1, 1);
    ASSERT_TRUE(p.ok);
    ASSERT_EQ(p.steps.size(), 1u) << "nothing to make first";
}

TEST(CraftPlan, AShortageIsNamedAndNothingIsPlanned) {
    const auto rs = TestRecipes();
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, { { "copper_ingot", 1 }, { "steel_ingot", 1 }, { "glass", 1 } }, 1, 1);
    EXPECT_FALSE(p.ok);
    EXPECT_TRUE(p.steps.empty());
    EXPECT_EQ(p.missing.at("copper_ingot"), 1) << "3 wire need 2 ingots";
}

TEST(CraftPlan, UnknownRecipesAreNotUsed) {
    const auto rs = TestRecipes();
    const auto noWire = [](int k) { return k != 0; };
    const CraftPlan p = CraftPlan_Make(rs, noWire, { { "copper_ingot", 5 }, { "steel_ingot", 1 }, { "glass", 1 } }, 1, 1);
    EXPECT_FALSE(p.ok);
    EXPECT_EQ(p.missing.at("copper_wire"), 3);
}

TEST(CraftPlan, ACycleEnds) {
    const auto rs = TestRecipes();
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, {}, 5, 1);
    EXPECT_FALSE(p.ok);
}

TEST(CraftPlan, MachineOnlyIngredientsAreNamed) {
    const auto rs = TestRecipes();
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, { { "glass", 1 } }, 7, 1, { "iron_ingot" });
    EXPECT_FALSE(p.ok);
    EXPECT_EQ(p.machineOnly.count("iron_ingot"), 1u) << "the page can say: needs a machine";
}

TEST(CraftPlan, GroupsMatchAnyMember) {
    const auto rs = TestRecipes();
    const CraftPlan p = CraftPlan_Make(rs, kAllKnown, { { "birch_log", 2 } }, 4, 2);
    EXPECT_TRUE(p.ok);
}
#endif
