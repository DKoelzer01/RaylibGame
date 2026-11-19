#pragma once
// Recursive hand crafting (#928, Factorio-style): a known recipe can be crafted when its ingredients are carried or can be
// hand-crafted from what is carried, through known hand recipes only. Machine recipes are never planned: an ingredient
// only a machine makes counts as missing (and is named, so the page can say "needs a furnace").
// Pure: the recipes, which are known, what is carried and which items machines make come in; the steps come out.
#include "inventory.h"
#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

struct CraftStep { int recipe = -1; int times = 0; };
struct CraftPlan {
    bool ok = false;
    std::vector<CraftStep> steps;            // intermediates first, the target last (each step's inputs exist when it runs)
    std::map<std::string, int> missing;      // raw shortfall (an id or a group such as any_planks)
    std::set<std::string> machineOnly;       // missing items a machine recipe makes (no hand recipe)
    std::map<std::string, int> leftover;     // made beyond the need (e.g. wire comes in pairs)
};

// Plan `times` crafts of recipes[target] from `have` (item id -> count carried).
CraftPlan CraftPlan_Make(const std::vector<Recipe>& recipes, const std::function<bool(int)>& known,
                         const std::map<std::string, int>& have, int target, int times,
                         const std::set<std::string>& machineMade = {});
// The most crafts of recipes[target] a plan can make (0..cap).
int CraftPlan_Max(const std::vector<Recipe>& recipes, const std::function<bool(int)>& known,
                  const std::map<std::string, int>& have, int target, int cap = 999,
                  const std::set<std::string>& machineMade = {});
// "6 copper wire, then 3 circuits": the steps in words (item names from itemDB when it has them).
std::string CraftPlan_Describe(const std::vector<Recipe>& recipes, const CraftPlan& plan);
