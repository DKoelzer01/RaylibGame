#pragma once
// Biome variants (#1964, #1977): every biome type (Biome) splits into variants by a regional noise
// (PlanetField::variantAt); resources/data/biomes.txt lists them. A variant id indexes BiomeVariants().defs; the first
// variant of each type is its base (the biome as it was before variants) and the fallback where a gate fails.
#include "voxel.h"
#include "woodspecies.h"
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
    // species weights in WoodSpecies order (woodspecies.h; all 0 = the type's broadleaf/conifer mix).
    static constexpr int kSpecies = kWoodSpecies;
    float trees = 1.0f, bushes = 1.0f, stones = 1.0f, size = 1.0f, cover = -1.0f;
    float species[kSpecies] = {};
    // (#2012) Rubber: the percent of its felled trees that give one rubber (deterministic per tree:
    // BiomeVariants_RubberRoll). The jungle type gives it from every tree whatever this says.
    int rubber = 0;
    // (#2032) Colour multipliers (linear RGB) on the grass and on leaves (trees, bushes, the far canopy), over the climate
    // colour. A variant without its own takes its type's base's; < 0 = not given (resolved: BiomeVariantTable::grassTint).
    float grassTint[3] = { -1.0f, -1.0f, -1.0f }, leafTint[3] = { -1.0f, -1.0f, -1.0f };
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
    std::vector<Vector3> grassTint, leafTint;   // (#2032) per id, resolved: its own, else its base's, else 1
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
// (#2012) Whether the plant with this id is among a percent of trees: a hash of the id, so a host and its clients agree
// and a tree gives the same answer however often it is asked.
bool BiomeVariants_RubberRoll(int percent, uint64_t plantId);

// (#2032) Per-biome tints. The climate code the mesher puts in terrain vertices (and grass and trees carry in their up
// vector's length): 2 + the temperature step (0-62 over -30..40 C) + 64 x the variant id (-1 = none: 0) + moisture as
// the fraction; below 1.5 = no climate. The shaders (biometint.glsl) decode the id and look its tints up in uniform
// arrays BiomeTint_Bind sets.
constexpr int kMaxBiomeVariants = 64; // biometint.glsl's array size
float BiomeClimateCode(float tempC, float moist, int variant);
bool BiomeTintOn(); // RAYLIBGAME_BIOME_TINT=0: off (every tint 1, the climate gradient at full strength, as before #2032)
// The tint arrays and the climate gradient's weight, for a shader that includes biometint.glsl (locations per program).
void BiomeTint_Bind(const Shader& s);
// How much of the climate gradient stays inside a biome, around the temperate colour: 1 = all (tints off).
float BiomeClimateWeight();
