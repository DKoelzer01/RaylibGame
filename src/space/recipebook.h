#ifndef RECIPEBOOK_H
#define RECIPEBOOK_H

// Which recipes the player knows. A recipe is known once every one of its ingredients has been
// discovered (held or obtained at least once): finding iron ingots reveals iron tools and parts,
// and so on. resources/data/recipe_unlocks.txt can add a further condition per output (a task in
// the quest book's terms: "explorer_suit | do stalker 1"). Creative mode knows everything.

#include <string>

struct Recipe;

bool RecipeBook_Known(const Recipe& r);
bool RecipeBook_Discovered(const std::string& item);
// Search: a recipe matches when its output's name or id, or any ingredient's, contains the words.
bool RecipeBook_Matches(const Recipe& r, const std::string& lowerQuery);
// Once a second: toasts "New recipes: ..." for recipes that became known (not on the first check).
void RecipeBook_Update();
void RecipeBook_Reset();                 // a new world or a load: the next check is silent
int RecipeBook_KnownCount();

#endif
