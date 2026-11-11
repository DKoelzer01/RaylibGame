#include "tools.h"

#include <cstring>

namespace {
    const char* const kTiers[] = { "hands", "stone", "iron", "steel", "titanium" };
    const char* const kTypes[] = { "hands", "pickaxe", "shovel", "axe" };
}

const char* ToolTierName(int tier) { return kTiers[tier < 0 ? 0 : tier > 4 ? 4 : tier]; }
const char* ToolTypeName(ToolType t) { return kTypes[(int)t]; }

ToolSpec ToolOf(const std::string& id) {
    ToolSpec s;
    const size_t u = id.find('_');
    if (u == std::string::npos) return s;
    const std::string tier = id.substr(0, u), type = id.substr(u + 1);
    for (int t = 1; t <= 4; ++t) if (tier == kTiers[t]) s.tier = t;
    if (s.tier == 0) return ToolSpec{};
    if (type == "pickaxe") s.type = ToolType::Pickaxe;
    else if (type == "shovel") s.type = ToolType::Shovel;
    else if (type == "axe") s.type = ToolType::Axe;
    else return ToolSpec{};
    return s;
}

bool IsSoft(VoxelMaterial m) {
    switch (m) {
        case VoxelMaterial::Regolith: case VoxelMaterial::Sand: case VoxelMaterial::Gravel:
        case VoxelMaterial::Grass: case VoxelMaterial::Snow: return true;
        default: return false;
    }
}

int PickTierFor(VoxelMaterial m, const char* item) {
    if (IsSoft(m)) return 0;
    if (m == VoxelMaterial::Ore && item) {
        if (strcmp(item, "bauxite_ore") == 0) return 2;  // an iron pickaxe
        if (strcmp(item, "ilmenite_ore") == 0) return 3; // steel
    }
    return 1; // stone, ice, the common ores: any pickaxe
}

DigRule DigRuleFor(ToolSpec tool, VoxelMaterial m, const char* item) {
    DigRule r;
    const int t = tool.tier;
    if (IsSoft(m)) {
        r.can = true;
        switch (tool.type) {
            case ToolType::Shovel:  r.radius = 0.55f + 0.3f * t; r.strokesPerSecond = 3.0f + 0.5f * t; break;
            case ToolType::Pickaxe: r.radius = 0.4f + 0.15f * t; r.strokesPerSecond = 2.5f; break;
            default:                r.radius = 0.35f; r.strokesPerSecond = 2.0f; break; // hands, or an axe
        }
        return r;
    }
    const int need = PickTierFor(m, item);
    if (tool.type != ToolType::Pickaxe) { r.why = "Needs a pickaxe"; return r; }
    if (t < need) { r.why = std::string("Needs a ") + kTiers[need] + " pickaxe or better"; return r; }
    r.can = true;
    r.radius = 0.45f + 0.25f * t;
    r.strokesPerSecond = 1.5f + 0.5f * t;
    return r;
}

float ChopFactor(ToolSpec tool) {
    if (tool.type != ToolType::Axe) return 2.0f; // bare hands (or the wrong tool): twice the swings
    static const float k[] = { 2.0f, 1.0f, 0.7f, 0.5f, 0.35f };
    return k[tool.tier < 0 ? 0 : tool.tier > 4 ? 4 : tool.tier];
}
