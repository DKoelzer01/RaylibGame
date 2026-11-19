#ifndef RECIPEBOOK_H
#define RECIPEBOOK_H

// Which recipes the player knows. A recipe is known once every one of its ingredients has been
// discovered (held or obtained at least once): finding iron ingots reveals iron tools and parts,
// and so on. resources/data/recipe_unlocks.txt can add a further condition per output (a task in
// the quest book's terms: "explorer_suit | do stalker 1"). Creative mode knows everything.

#include <functional>
#include <set>
#include <string>

struct Recipe;

// One player's recipe book (#313, in their PlayerSession): outputs known at the last check.
struct RecipeBookState {
    std::set<std::string> known;
    bool primed = false;      // the first check after a reset is silent
    double lastCheck = -10.0;
};

bool RecipeBook_Known(const Recipe& r);
bool RecipeBook_Unlocked(const Recipe& r);  // (#338) its recipe_unlocks condition alone, on the selected quest book (a server: the team's)
// (#338) Craft recipe `index` (into itemDB.recipes()) `count` times: in single player here, from the inventory; in a joined game a request
// to the server (gCraftRequest; the result comes back as inventory deltas). The one entry point for the crafting screens.
bool Crafting_Request(int index, int count = 1);
// (#946) Joined: craft it with its missing intermediates first, the server planning them (#928); false when not joined (the caller
// queues the plan's steps here).
bool Crafting_RequestRecursive(int index, int count);
inline std::function<void(int index, const std::string& output, int count, bool useChests, bool recursive)> gCraftRequest; // useChests (#854): what the inventory lacks, from the chests near the player
bool RecipeBook_Discovered(const std::string& item);
// Search: a recipe matches when its output's name or id, or any ingredient's, contains the words.
bool RecipeBook_Matches(const Recipe& r, const std::string& lowerQuery);
// Once a second: toasts "New recipes: ..." for recipes that became known (not on the first check).
void RecipeBook_Update();
void RecipeBook_Reset();                 // a new world or a load: the next check is silent
int RecipeBook_KnownCount();

#endif
