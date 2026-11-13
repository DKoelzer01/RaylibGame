#include "toolwear.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "tools.h"
#include "imgui.h"

#include <raylib.h>
#include <cstring>
#include <map>

namespace {
    std::map<std::string, int> used; // strokes on the tool in use, per item id
}

void ToolWear_OnStroke(const std::string& id, int strokes) {
    const ToolSpec t = ToolOf(id);
    const int max = ToolDurability(t);
    if (max <= 0 || strokes <= 0 || inventory.count(id) <= 0) return;
    int& u = used[id];
    u += strokes;
    if (u < max) {
        if (u == max - max / 10) GameUI_Toast(TextFormat("Your %s %s is nearly worn out (repair it: Tab)", ToolTierName(t.tier), ToolTypeName(t.type)));
        return;
    }
    u = 0;
    inventory.remove(id, 1);
    Audio_Play(Sfx::Error, 0.6f, 0.7f);
    GameUI_Toast(TextFormat("Your %s %s broke", ToolTierName(t.tier), ToolTypeName(t.type)));
}

int ToolWear_Used(const std::string& id) { auto it = used.find(id); return it == used.end() ? 0 : it->second; }

float ToolWear_Left(const std::string& id) {
    const int max = ToolDurability(ToolOf(id));
    if (max <= 0) return 1.0f;
    return 1.0f - (float)ToolWear_Used(id) / (float)max;
}

bool ToolWear_RepairCost(const std::string& id, std::string& item, int& count) {
    const ToolSpec t = ToolOf(id);
    if (ToolDurability(t) <= 0) return false;
    static const char* kHead[] = { "", "stone", "iron_ingot", "steel_ingot", "titanium_ingot" };
    item = kHead[t.tier];
    count = t.tier == 1 ? 2 : 1;
    return true;
}

bool ToolWear_Repair(const std::string& id) {
    std::string item; int n = 0;
    if (ToolWear_Used(id) <= 0 || !ToolWear_RepairCost(id, item, n)) return false;
    if (inventory.count(item) < n || !inventory.remove(item, n)) {
        GameUI_Toast(TextFormat("Repairing needs %d %s", n, item.c_str()));
        return false;
    }
    used[id] = 0;
    Audio_Play(Sfx::Craft, 0.5f);
    return true;
}

void ToolWear_DrawRepairUI() {
    bool header = false;
    for (const auto& kv : used) {
        if (kv.second <= 0 || inventory.count(kv.first) <= 0) continue;
        if (!header) { ImGui::TextUnformatted("Tool repair"); header = true; }
        std::string item; int n = 0;
        ToolWear_RepairCost(kv.first, item, n);
        const ItemDef* d = itemDB.find(kv.first);
        ImGui::PushID(kv.first.c_str());
        if (ImGui::SmallButton("Repair")) ToolWear_Repair(kv.first);
        ImGui::SameLine();
        ImGui::Text("%s  %.0f%% left  (%d %s)", d ? d->name.c_str() : kv.first.c_str(), ToolWear_Left(kv.first) * 100.0f, n, item.c_str());
        ImGui::PopID();
    }
    if (header) ImGui::Separator();
}

void ToolWear_Reset() { used.clear(); }

void ToolWear_Save(FILE* f) {
    for (const auto& kv : used) if (kv.second > 0) fprintf(f, "tool_wear %s %d\n", kv.first.c_str(), kv.second);
}

bool ToolWear_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "tool_wear") != 0) return false;
    char id[64] = {}; int n = 0;
    if (fscanf(f, " %63s %d", id, &n) == 2) used[id] = n;
    return true;
}
