#include "toolwear.h"
#include "playersession.h"
#include "audio.h"
#include "gameui.h"
#include "inventory.h"
#include "tools.h"
#include "composed.h"
#include "itemdata.h"
#include "hotbar.h"
#include "imgui.h"

#include <raylib.h>
#include <cstring>
#include <algorithm>
#include <map>

namespace {
    std::map<std::string, int>& Used() { return LocalSession().toolWear.used; } // the local player's (#313)
    // (#1449) A composed tool's copy in use: the selected hotbar cell if it is one, else the first held.
    int HeldCell(const std::string& id) {
        const auto& g = inventory.grid();
        const int p = inventory.preferred;
        if (p >= 0 && p < Inventory::kCells && g[p].id == id && g[p].hasData()) return p;
        for (int i = 0; i < Inventory::kCells; ++i) if (g[i].id == id && g[i].hasData()) return i;
        return -1;
    }
    // A composed item's repair: a little of its head's material.
    bool ComposedRepairCost(const std::string& head, std::string& item, int& count) {
        static const std::map<std::string, std::pair<const char*, int>> k = {
            { "stone", { "stone", 2 } }, { "flint", { "flint", 1 } }, { "bone", { "bone", 1 } }, { "bronze", { "bronze_ingot", 1 } },
            { "iron", { "iron_ingot", 1 } }, { "steel", { "steel_ingot", 1 } }, { "titanium", { "titanium_ingot", 1 } } };
        auto it = k.find(head);
        if (it == k.end()) return false;
        item = it->second.first; count = it->second.second;
        return true;
    }
}

std::string ToolWear_HeldData(const std::string& id) { const int c = HeldCell(id); return c < 0 ? std::string() : inventory.grid()[c].data; }

void ToolWear_OnStroke(const std::string& id, int strokes) {
    if (Composed_Is(id)) { // (#1449) wear in the copy's own record (w): each copy wears on its own
        const int c = HeldCell(id);
        if (c < 0 || strokes <= 0) return;
        const std::string data = inventory.grid()[c].data;
        const ComposedStats st = Composed_Stats(id, data);
        if (!st.valid || st.durability <= 0) return;
        const int u = st.wear + strokes, warn = st.durability - st.durability / 10;
        if (u < st.durability) {
            inventory.setInstanceData(c, ItemData::SetInt(data, "w", u));
            if (st.wear < warn && u >= warn) GameUI_Toast(TextFormat("Your %s is nearly worn out (repair it: Tab)", st.name.c_str()));
            return;
        }
        inventory.takeInstance(c);
        if (inventory.countInstances(id) <= 0 && inventory.count(id) <= 0) Hotbar_Unbind(id);
        Audio_Play(Sfx::Error, 0.6f, 0.7f);
        GameUI_Toast(TextFormat("Your %s broke", st.name.c_str()));
        return;
    }
    const ToolSpec t = ToolOf(id);
    const int max = ToolDurability(t);
    if (max <= 0 || strokes <= 0 || inventory.count(id) <= 0) return;
    int& u = Used()[id];
    u += strokes;
    if (u < max) {
        if (u == max - max / 10) GameUI_Toast(TextFormat("Your %s %s is nearly worn out (repair it: Tab)", ToolTierName(t.tier), ToolTypeName(t.type)));
        return;
    }
    u = 0;
    inventory.remove(id, 1);
    // The last of it broke: it leaves the hotbar (#209: the binding stayed, an empty slot showing the tool).
    if (inventory.count(id) <= 0)
        Hotbar_Unbind(id);
    Audio_Play(Sfx::Error, 0.6f, 0.7f);
    GameUI_Toast(TextFormat("Your %s %s broke", ToolTierName(t.tier), ToolTypeName(t.type)));
}

int ToolWear_Used(const std::string& id) {
    if (Composed_Is(id)) { const int c = HeldCell(id); return c < 0 ? 0 : std::max(0, ItemData::GetInt(inventory.grid()[c].data, "w", 0)); } // (#1449)
    auto it = Used().find(id); return it == Used().end() ? 0 : it->second;
}

float ToolWear_Left(const std::string& id) {
    const int max = Composed_Is(id) ? ToolDurabilityOf(id, ToolWear_HeldData(id)) : ToolDurability(ToolOf(id));
    if (max <= 0) return 1.0f;
    return 1.0f - (float)ToolWear_Used(id) / (float)max;
}

bool ToolWear_RepairCost(const std::string& id, std::string& item, int& count) {
    if (Composed_Is(id)) return ComposedRepairCost(Composed_Stats(id, ToolWear_HeldData(id)).head, item, count); // (#1449)
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
    if (Composed_Is(id)) { const int c = HeldCell(id); if (c >= 0) inventory.setInstanceData(c, ItemData::Set(inventory.grid()[c].data, "w", "")); } // (#1449)
    else Used()[id] = 0;
    Audio_Play(Sfx::Craft, 0.5f);
    return true;
}

void ToolWear_DrawRepairUI() {
    bool header = false;
    for (int i = 0; i < Inventory::kCells; ++i) { // (#1449) composed copies, each with its own wear
        const InvSlot& s = inventory.grid()[i];
        if (s.empty() || !s.hasData() || !Composed_Is(s.id)) continue;
        const ComposedStats st = Composed_Stats(s.id, s.data);
        if (!st.valid || st.wear <= 0 || st.durability <= 0) continue;
        if (!header) { ImGui::TextUnformatted("Tool repair"); header = true; }
        std::string item; int n = 0;
        ComposedRepairCost(st.head, item, n);
        ImGui::PushID(i);
        if (ImGui::SmallButton("Repair")) {
            if (inventory.count(item) >= n && inventory.remove(item, n)) { inventory.setInstanceData(i, ItemData::Set(s.data, "w", "")); Audio_Play(Sfx::Craft, 0.5f); }
            else GameUI_Toast(TextFormat("Repairing needs %d %s", n, item.c_str()));
        }
        ImGui::SameLine();
        ImGui::Text("%s  %.0f%% left  (%d %s)", st.name.c_str(), 100.0f * (1.0f - (float)st.wear / (float)st.durability), n, item.c_str());
        ImGui::PopID();
    }
    for (const auto& kv : Used()) {
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

void ToolWear_Reset() { Used().clear(); }

void ToolWear_Save(FILE* f) {
    for (const auto& kv : Used()) if (kv.second > 0) fprintf(f, "tool_wear %s %d\n", kv.first.c_str(), kv.second);
}

bool ToolWear_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "tool_wear") != 0) return false;
    char id[64] = {}; int n = 0;
    if (fscanf(f, " %63s %d", id, &n) == 2) Used()[id] = n;
    return true;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1449): a composed tool wears in its own record and breaks at its own durability ----
#include "unit_helpers.h"

TEST(ToolWear, ComposedToolsWearInTheirOwnRecord) {
    UnitTest::EnsureItemDB();
    ASSERT_TRUE(Composed_Load());
    gToolDataOf = ToolWear_HeldData;
    inventory.clear();
    const int flint = inventory.addInstance("tool_pickaxe", Composed_MakeData("flint"));
    const int iron = inventory.addInstance("tool_pickaxe", Composed_MakeData("iron"));
    ASSERT_GE(flint, 0); ASSERT_GE(iron, 0);
    inventory.preferred = iron;                     // the iron one in hand
    EXPECT_EQ(ToolOf("tool_pickaxe").tier, 2) << "the held copy's head sets the tier";
    ToolWear_OnStroke("tool_pickaxe", 10);
    EXPECT_EQ(ItemData::GetInt(inventory.grid()[iron].data, "w", 0), 10);
    EXPECT_EQ(ItemData::GetInt(inventory.grid()[flint].data, "w", 0), 0) << "the other copy is untouched";
    EXPECT_NEAR(ToolWear_Left("tool_pickaxe"), 1.0f - 10.0f / 400.0f, 1e-4f);
    inventory.preferred = flint;                    // switch to the flint one: tier 1, its own wear
    EXPECT_EQ(ToolOf("tool_pickaxe").tier, 1);
    const int life = ToolDurabilityOf("tool_pickaxe", inventory.grid()[flint].data);
    EXPECT_EQ(life, 125) << "flint: toughness 5 x 25";
    ToolWear_OnStroke("tool_pickaxe", life);         // worn out: it breaks, the iron one stays
    EXPECT_EQ(inventory.countInstances("tool_pickaxe"), 1);
    EXPECT_TRUE(inventory.grid()[flint].empty());
    EXPECT_EQ(ItemData::GetInt(inventory.grid()[iron].data, "w", 0), 10);
    inventory.clear(); inventory.preferred = -1;
    gToolDataOf = nullptr;
}
#endif
