#pragma once
// Biome variants (#1964, #1977): every biome type (Biome) splits into variants by a regional noise
// (PlanetField::variantAt); resources/data/biomes.txt lists them. A variant id indexes BiomeVariants().defs; the first
// variant of each type is its base (the biome as it was before variants) and the fallback where a gate fails.
#include "voxel.h"
#include <string>
#include <vector>

struct BiomeVariantDef {
    std::string key;          // "forest.birch_grove": the type's lower-case name, a dot, the variant's id
    std::string name;         // shown to players: "Birch grove"
    Biome type = Biome::Forest;
    float weight = 1.0f;      // its share of the type's area, relative to the type's other variants
    // Gates (optional): the variant only exists inside these ranges; elsewhere the point falls to the type's next band.
    float elevMin = -1e9f, elevMax = 1e9f;   // m above sea level
    float tempMin = -1e9f, tempMax = 1e9f;   // C (the climate's annual mean, PlanetField::Climate::temp)
    float moistMin = -1e9f, moistMax = 1e9f; // 0..1
    // Flora (#1978), applied over the type's own (trees.cpp FloraOf): density multipliers for trees, bushes and loose
    // stones, a tree size multiplier, an optional closed-canopy floor (least tree odds; < 0 = the type's), and optional
    // species weights in WoodSpecies order (oak, birch, maple, pine, spruce; all 0 = the type's broadleaf/conifer mix).
    static constexpr int kSpecies = 5;
    float trees = 1.0f, bushes = 1.0f, stones = 1.0f, size = 1.0f, cover = -1.0f;
    float species[kSpecies] = {};
    bool hasSpecies() const { for (float w : species) if (w > 0.0f) return true; return false; }
    bool gated() const { return elevMin > -1e8f || elevMax < 1e8f || tempMin > -1e8f || tempMax < 1e8f || moistMin > -1e8f || moistMax < 1e8f; }
    bool admits(float elev, float temp, float moist) const {
        return elev >= elevMin && elev <= elevMax && temp >= tempMin && temp <= tempMax && moist >= moistMin && moist <= moistMax;
    }
};

struct BiomeVariantTable {
    std::vector<BiomeVariantDef> defs;
    std::vector<int> byType[(int)Biome::Count]; // ids per type in file order, the base first
    std::vector<float> bandLo, bandHi;          // per id: its band of the variant axis (0..1) within its type
    int find(const std::string& key) const;     // -1 if none
    int base(Biome t) const { return byType[(int)t].empty() ? -1 : byType[(int)t][0]; }
};

// Parses biomes.txt's text into out (replacing it). false with a message naming the line on any error: an unknown key,
// type or value, a type without variants, or a gated base.
bool BiomeVariants_Parse(const std::string& text, BiomeVariantTable& out, std::string* err);
void BiomeVariants_Default(BiomeVariantTable& out); // one base variant per type, named after it
// resources/data/biomes.txt, loaded on first use (thread-safe); the defaults if it is missing or bad (logged).
const BiomeVariantTable& BiomeVariants();
const char* BiomeVariantName(int id); // its display name; "" for an unknown id
