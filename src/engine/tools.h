#ifndef TOOLS_H
#define TOOLS_H

// Hand tools: pickaxes (stone, ore; the tier sets which ores), shovels (soft ground: dirt, sand,
// gravel, snow) and axes (trees), in stone, iron, steel and titanium. Without a tool the player can
// only scrape soft ground a little at a time and fell trees slowly. Each tier digs a bigger sphere.
// Items are "<tier>_<type>" (stone_pickaxe ... titanium_axe); see items.txt and recipes.txt.

#include <string>
#include "voxel.h"

enum class ToolType { Hands, Pickaxe, Shovel, Axe, Hoe };

struct ToolSpec {
    ToolType type = ToolType::Hands;
    int tier = 0; // 0 hands, 1 stone, 2 iron, 3 steel, 4 titanium
};

ToolSpec ToolOf(const std::string& itemId);  // hands for anything that is not a tool
const char* ToolTierName(int tier);          // "stone", "iron"...
const char* ToolTypeName(ToolType t);        // "pickaxe"...

// What digging `material` (holding `item` when it is ore) with a tool does.
struct DigRule {
    bool can = false;
    float radius = 0.0f;           // metres
    float strokesPerSecond = 0.0f;
    std::string why;               // when it can't: what is needed
};
DigRule DigRuleFor(ToolSpec tool, VoxelMaterial material, const char* item);
bool IsSoft(VoxelMaterial m);
int PickTierFor(VoxelMaterial m, const char* item); // the pickaxe tier stone and each ore need

// Swings to fell a tree, relative to an axe-less base count (bushes are unaffected).
float ChopFactor(ToolSpec tool);

// Strokes a tool lasts before it breaks (hands: unlimited, 0): stone 150, iron 400, steel 900, titanium 2000.
int ToolDurability(ToolSpec tool);

#endif
