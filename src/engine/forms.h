#ifndef FORMS_H
#define FORMS_H

// (#1614, machine components, PR 1) The form table: how each component shape (plate, rod, wire, gear, screw, spring) is
// made of a metal, by hand and by machine, from resources/data/forms.txt, applied to every material in materials.txt that
// has the shape (its forms column). One line per shape and route, not one per metal. The item loader (ItemDB::load) adds
// the hand recipes and the factory (LoadRecipes) the machine recipes, each only where the data has no recipe for that
// output already (an explicit recipe wins: iron_rod 2 = iron_ingot 1 stays).

#include <string>
#include <utility>
#include <vector>

struct FormRecipe {
    std::string machine;                                  // "" by hand; else a machines.txt name
    std::string output;                                   // <material>_<shape>
    int outCount = 1;
    std::vector<std::pair<std::string, int>> inputs;      // <material>_<form> ids (or fixed items)
    int ticks = 0, eu = 0;                                // machine recipes
};

// The recipes the form table gives for the materials (forms.txt and materials.txt as text: tests pass their own).
std::vector<FormRecipe> Forms_Generate(const std::string& formsText, const std::string& materialsText);
// From resources/data (the game's files), loaded once.
const std::vector<FormRecipe>& Forms_All();
const std::vector<std::string>& Forms_Shapes(); // the shapes forms.txt names (plate, rod, ...)

#endif
