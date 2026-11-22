#ifndef OREFAMILIES_H
#define OREFAMILIES_H

// (#1467, the underground plan #1465) Ore vein families from resources/data/ore_families.txt: a deposit is a family (a
// main ore and its byproducts, GregTech style) found in its host rocks, depth band, biomes and underground regions, in one
// vein shape. PlanetField picks a family per deposit cell and a member per point (Lighting's side); this is the table.
// Members are ores.txt ids, resolved to ore kinds at load. Loaded once, on first use (thread-safe).
// Format agreed with Lighting (2026-11-22): hosts are hard eligibility (a non-host rock at the deposit centre's depth
// skips the family; no eligible family: no deposit); size is a range rolled per deposit; shape extents stay in code.

#include <string>
#include <vector>

enum class VeinShape { Blob = 0, Sheet, Pipe, Seam, Stockwork }; // blob: today's deposits; sheet: a dike or vein plane

struct OreFamilyMember { std::string ore; float weight = 1.0f; int kind = -1; }; // kind: its ores.txt row (-1 unknown)
struct OreFamilyBias { std::string name; float factor = 1.0f; };                 // a biome or region and its multiplier

struct OreFamily {
    std::string id;
    std::vector<OreFamilyMember> members;  // the main ore first
    std::vector<std::string> hosts;        // VoxelMaterial rock names (rock, granite, basalt, ...); empty: any
    float minDepth = 0.0f, maxDepth = 1e9f; // m under the original surface
    std::vector<OreFamilyBias> biomes;     // surface biome multipliers (unlisted x1); empty: any
    std::vector<OreFamilyBias> regions;    // underground region multipliers (#1547: none, karst, volcanic, ...); unlisted x1
    float rarity = 1.0f;                   // its weight among the families eligible in a cell
    VeinShape shape = VeinShape::Blob;
    float minSize = 4.0f, maxSize = 8.0f;  // m, rolled per deposit: blob or pipe radius, sheet or seam thickness, vein spacing
    float totalWeight = 0.0f;              // the members' weights summed
};

// Parse ore_families.txt lines (comments and blank lines skipped). Bad lines are skipped with a warning.
std::vector<OreFamily> OreFamilies_Parse(const std::string& text);
const std::vector<OreFamily>& OreFamilies_All();
const OreFamily* OreFamilies_Find(const std::string& id);
// A member by a uniform r in [0, 1): weighted, in row order. nullptr for a family without members.
const OreFamilyMember* OreFamilies_Member(const OreFamily& f, float r);
// The family's weight in a cell: rarity x its biome multiplier x its region multiplier. Names match case-insensitively
// (BiomeName() prints "Desert"; the file says "desert").
float OreFamilies_Rarity(const OreFamily& f, const std::string& biome, const std::string& region = "none");
bool OreFamilies_Hosts(const OreFamily& f, const std::string& rock); // its hosts take this rock (any: true)
const char* OreFamilies_ShapeName(VeinShape s);
// Problems, one a line ("" none): members not in ores.txt, weights not over 0, no members, depth min >= max, unknown hosts,
// biomes or regions, rarity not over 0, size not over 0 or min > max, duplicate ids.
std::string OreFamilies_Validate(const std::vector<OreFamily>& fams);

#endif
