#ifndef COMPONENTS_H
#define COMPONENTS_H

// (#1614, machine components) The cost report: an item's raw cost in gathering seconds, expanded through the recipes
// (hand and machine) down to what the world gives, each raw item weighted by its gathering time (resources/data/gather.txt).
// Report only (the user's answer to Q4): the bench prints it, nothing fails on it. Pure: recipes and times are passed in.

#include <map>
#include <string>
#include <utility>
#include <vector>

struct CostRecipe {
    bool hand = true;                                   // false: a machine's
    std::string output;
    int outCount = 1;
    std::vector<std::pair<std::string, int>> inputs;
};

enum class CostRoute { Hand, Best }; // Hand: hand recipes where an item has one (else its cheapest machine recipe); Best: the cheapest of all

struct CostModel {
    std::vector<CostRecipe> recipes;
    std::map<std::string, double> gatherSeconds;        // raw items
    double unknownSeconds = 10.0;                       // a raw item with no time (GUESS), counted in `unknown`
};

// Seconds to gather what one `item` takes by that route. A cycle or an item with no recipe and no time counts as raw at
// unknownSeconds; `unknown` (if given) collects those items.
double Components_Cost(const CostModel& m, const std::string& item, CostRoute route, std::vector<std::string>* unknown = nullptr);
// gather.txt: "<item> <seconds>" lines (comments: #)
std::map<std::string, double> Components_ParseGather(const std::string& text);

#endif
