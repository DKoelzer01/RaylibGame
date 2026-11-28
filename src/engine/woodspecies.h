#pragma once
// Tree species (#1995 phase 2): their identity and order. A tree variant is species x size class (trees.h), so the order
// is the variant layout; plants are not saved by variant (removals keep plant ids), so appending or reordering only moves
// meshes. Each species' form (its structure) is code in trees.cpp; its log, tints, bark texture and size come from
// resources/data/trees.txt (WoodSpeciesInfo). biomes.txt `species` names them.
#include <cstdint>

enum class WoodSpecies : uint8_t { Oak, Birch, Maple, Pine, Spruce, Acacia, Kapok, Mangrove, Krummholz, Snag, Baobab, TreeFern, Palm, Willow, Count };
constexpr int kWoodSpecies = (int)WoodSpecies::Count;
inline constexpr const char* kWoodSpeciesNames[kWoodSpecies] = {
    "oak", "birch", "maple", "pine", "spruce", "acacia", "kapok", "mangrove", "krummholz", "snag", "baobab", "tree_fern", "palm", "willow",
};
inline int WoodSpeciesByName(const char* name) {
    for (int i = 0; i < kWoodSpecies; ++i) {
        const char *a = kWoodSpeciesNames[i], *b = name;
        while (*a && *a == *b) { ++a; ++b; }
        if (!*a && !*b) return i;
    }
    return -1;
}
