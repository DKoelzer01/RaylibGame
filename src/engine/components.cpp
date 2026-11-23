#include "components.h"

#include <algorithm>
#include <cstdlib>
#include <set>
#include <sstream>

namespace {
double CostOf(const CostModel& m, const std::string& item, CostRoute route, std::set<std::string>& onPath, std::vector<std::string>* unknown) {
    if (auto g = m.gatherSeconds.find(item); g != m.gatherSeconds.end()) return g->second; // raw: gathered
    bool anyHand = false;
    for (const CostRecipe& r : m.recipes) anyHand = anyHand || (r.hand && r.output == item);
    if (onPath.count(item)) { if (unknown) unknown->push_back(item + " (cycle)"); return m.unknownSeconds; }
    onPath.insert(item);
    double best = -1.0;
    for (const CostRecipe& r : m.recipes) {
        if (r.output != item || r.outCount < 1) continue;
        if (route == CostRoute::Hand && anyHand && !r.hand) continue; // the hand route takes the hand recipe where there is one
        double c = 0.0;
        for (const auto& [in, n] : r.inputs) c += n * CostOf(m, in, route, onPath, unknown);
        c /= r.outCount;
        if (best < 0.0 || c < best) best = c;
    }
    onPath.erase(item);
    if (best < 0.0) { if (unknown) unknown->push_back(item); return m.unknownSeconds; }
    return best;
}
} // namespace

double Components_Cost(const CostModel& m, const std::string& item, CostRoute route, std::vector<std::string>* unknown) {
    std::set<std::string> onPath;
    return CostOf(m, item, route, onPath, unknown);
}

std::map<std::string, double> Components_ParseGather(const std::string& text) {
    std::map<std::string, double> out;
    std::istringstream in(text);
    for (std::string line; std::getline(in, line);) {
        const size_t a = line.find_first_not_of(" \t");
        if (a == std::string::npos || line[a] == '#') continue;
        std::istringstream ls(line);
        std::string id; double s = 0.0;
        if (ls >> id >> s) out[id] = s;
    }
    return out;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

TEST(Components, CostExpandsToRawSecondsByHandAndBestRoute) {
    CostModel m;
    m.gatherSeconds = Components_ParseGather("# raw\niron_ore 6\ncoal 5\n");
    m.recipes = {
        { true, "iron_ingot", 1, { { "iron_ore", 1 }, { "coal", 1 } } },   // 11 s
        { true, "iron_plate", 1, { { "iron_ingot", 2 } } },                // by hand: 22 s
        { false, "iron_plate", 1, { { "iron_ingot", 1 } } },               // the press: 11 s
        { true, "iron_screw", 2, { { "iron_rod", 1 } } },
        { true, "iron_rod", 1, { { "iron_ingot", 1 } } },                  // a screw: 5.5 s
        { true, "loop_a", 1, { { "loop_b", 1 } } }, { true, "loop_b", 1, { { "loop_a", 1 } } },
    };
    EXPECT_DOUBLE_EQ(Components_Cost(m, "iron_plate", CostRoute::Hand), 22.0) << "the hand route takes the hand recipe";
    EXPECT_DOUBLE_EQ(Components_Cost(m, "iron_plate", CostRoute::Best), 11.0) << "the best route takes the press";
    EXPECT_DOUBLE_EQ(Components_Cost(m, "iron_screw", CostRoute::Hand), 5.5) << "per item made (2 a craft)";
    std::vector<std::string> unknown;
    EXPECT_DOUBLE_EQ(Components_Cost(m, "loop_a", CostRoute::Best, &unknown), m.unknownSeconds) << "a cycle counts as raw";
    EXPECT_FALSE(unknown.empty());
    unknown.clear();
    EXPECT_DOUBLE_EQ(Components_Cost(m, "mystery", CostRoute::Best, &unknown), m.unknownSeconds);
    ASSERT_EQ(unknown.size(), 1u);
    EXPECT_EQ(unknown[0], "mystery");
}
#endif
