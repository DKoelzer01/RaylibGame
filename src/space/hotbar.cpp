#include "hotbar.h"
#include "toolwear.h"
#include "buildtool.h"
#include "gameui.h"
#include "inventory.h"
#include "plots.h"
#include "survival.h"
#include "audio.h"
#include "terraintool.h"
#include "tools.h"
#include "animals.h"

#include <raylib.h>
#include <cstring>

const char* const kHotbarTerrainTool = "tool:terrain";

namespace {
    std::string slots[kHotbarSlots];
    int selected = 0;
    const int kKeys[kHotbarSlots] = { KEY_ONE, KEY_TWO, KEY_THREE, KEY_FOUR, KEY_FIVE, KEY_SIX, KEY_SEVEN, KEY_EIGHT, KEY_NINE, KEY_ZERO };

    // A part type placed from this item (machines, conveyors, cables), or Count.
    PartType PartForItem(const std::string& id) {
        for (int t = 0; t < (int)PartType::Count; ++t)
            if (IsShape((PartType)t)) continue; else
            if (id == PartCostItem((PartType)t)) return (PartType)t;
        return PartType::Count;
    }

    const char* ShortName(const std::string& id) {
        if (id.empty()) return "";
        if (id == kHotbarTerrainTool) return "Hands";
        if (id == "oxygen_canister") return "O2 can";
        if (id == "power_cell") return "Power cell";
        if (id == "cable") return "Cable";
        const int m = BuildMaterialForItem(id.c_str());
        if (m >= 0) return GetBuildMaterial(m).name;
        const ItemDef* d = itemDB.find(id);
        return d ? d->name.c_str() : id.c_str();
    }

    void Select(int slot, bool piloting) {
        const std::string& id = slots[slot];
        if (id.empty()) { // an empty slot is bare hands
            if (piloting) return;
            selected = slot; buildTool.active = false; terrainTool.heldTool.clear();
            return;
        }
        // Food and animal feed: the first press takes it in hand (animals that eat it follow you); the
        // next press eats it (food only).
        if ((Survival_IsFood(id) || Animals_IsFeed(id)) && !piloting) {
            if (selected == slot && buildTool.active == false && terrainTool.heldTool.empty()) { if (Survival_IsFood(id)) Survival_UseItem(id); return; }
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); return; }
            selected = slot; buildTool.active = false; terrainTool.heldTool.clear();
            GameUI_Toast(Survival_IsFood(id) ? TextFormat("Holding %s: press %d again to eat", ShortName(id), (slot + 1) % 10) : TextFormat("Holding %s: animals that eat it follow you", ShortName(id)));
            return;
        }
        if (Survival_IsUsable(id)) { Survival_UseItem(id); return; } // consumables don't change the selection
        if (piloting) return;
        selected = slot;
        if (id == kHotbarTerrainTool) { buildTool.active = false; terrainTool.heldTool.clear(); return; }
        if (KindOf(id) == ItemKind::Tool && ToolOf(id).type != ToolType::Hands) {
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
            buildTool.active = false;
            terrainTool.heldTool = id;
            return;
        }
        if (TerrainTool_IsFiller(id) && BuildMaterialForItem(id.c_str()) < 0) { // loose ground (gravel, sand): in hand, RMB places it
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
            buildTool.active = false;
            terrainTool.heldTool = id;
            return;
        }
        const int m = BuildMaterialForItem(id.c_str());
        if (m >= 0) {
            buildTool.active = true;
            buildTool.material = m;
            if (!IsShape((PartType)buildTool.part)) buildTool.part = (int)PartType::Cube;
            return;
        }
        const PartType t = PartForItem(id);
        if (t != PartType::Count) {
            buildTool.active = true;
            buildTool.part = (int)t;
            return;
        }
        GameUI_Toast(TextFormat("%s can't be used from the hotbar", ShortName(id)));
        Audio_Play(Sfx::Error, 0.3f);
    }
}

void Hotbar_Reset() {
    // Empty: an empty slot is bare hands; tools bind themselves as they are made, the rest by the player.
    for (int i = 0; i < kHotbarSlots; ++i) slots[i].clear();
    selected = 0;
}

void Hotbar_Bind(int slot, const std::string& id) {
    if (slot < 0 || slot >= kHotbarSlots) return;
    for (auto& s : slots) if (s == id) s.clear(); // one slot per item
    slots[slot] = id;
    GameUI_Toast(TextFormat("%s on key %d", ShortName(id), (slot + 1) % 10));
}

const std::string& Hotbar_Slot(int slot) { return slots[slot]; }
int Hotbar_Selected() { return selected; }
std::string Hotbar_HeldFood() {
    if (selected < 0 || selected >= kHotbarSlots || buildTool.active || !terrainTool.heldTool.empty()) return "";
    const std::string& id = slots[selected];
    return (Survival_IsFood(id) || Animals_IsFeed(id)) && inventory.count(id) > 0 ? id : "";
}

int Hotbar_KeySlot() {
    for (int i = 0; i < kHotbarSlots; ++i) if (IsKeyPressed(kKeys[i])) return i;
    return -1;
}

void Hotbar_Press(int slot, bool piloting) { if (slot >= 0 && slot < kHotbarSlots) Select(slot, piloting); }

void Hotbar_Update(Scene&, bool inputEnabled, bool piloting) {
    static bool initialised = false;
    if (!initialised) { Hotbar_Reset(); initialised = true; }
    // G (build mode on/off) is the build tool's; keep the selection in step with it.
    // Leaving build mode with a block selected: the hand goes empty (to an empty slot, if there is one).
    const std::string& cur = slots[selected];
    const bool inHand = cur.empty() || cur == terrainTool.heldTool || Survival_IsFood(cur) || Animals_IsFeed(cur);
    if (!buildTool.active && !inHand) {
        for (int i = 0; i < kHotbarSlots; ++i) if (slots[i].empty()) { selected = i; break; }
        terrainTool.heldTool.clear();
    }
    // A tool lost (or never carried) falls back to bare hands.
    if (!terrainTool.heldTool.empty() && inventory.count(terrainTool.heldTool) <= 0) terrainTool.heldTool.clear();
    // New tools go on the hotbar: the first empty slot, else the first whose item isn't carried.
    for (const InvSlot& s : inventory.grid()) {
        if (s.empty() || KindOf(s.id) != ItemKind::Tool || ToolOf(s.id).type == ToolType::Hands) continue;
        bool bound = false;
        for (const auto& b : slots) bound |= b == s.id;
        if (bound) continue;
        int free = -1;
        for (int i = 0; i < kHotbarSlots && free < 0; ++i) if (slots[i].empty()) free = i;
        for (int i = 0; i < kHotbarSlots && free < 0; ++i)
            if (i != selected && inventory.count(slots[i]) <= 0) free = i;
        if (free >= 0) Hotbar_Bind(free, s.id);
    }
    if (!inputEnabled) return;
    const int k = Hotbar_KeySlot();
    if (k >= 0) Select(k, piloting);
}

void Hotbar_DrawHUD(bool piloting) {
    const int n = kHotbarSlots, size = 64, gap = 6;
    const int w = n * size + (n - 1) * gap;
    const int x0 = (GetScreenWidth() - w) / 2, y0 = GetScreenHeight() - size - 44;
    for (int i = 0; i < n; ++i) {
        const int x = x0 + i * (size + gap);
        const std::string& id = slots[i];
        const bool sel = !piloting && i == selected;
        DrawRectangle(x, y0, size, size, Fade(BLACK, sel ? 0.65f : 0.45f));
        DrawRectangleLines(x, y0, size, size, sel ? SKYBLUE : Fade(RAYWHITE, 0.35f));
        DrawText(TextFormat("%d", (i + 1) % 10), x + 4, y0 + 3, 14, Fade(RAYWHITE, 0.7f));
        if (id.empty()) continue;
        {
            const ItemDef* d = itemDB.find(id);
            if (d) DrawRectangle(x + size / 2 - 9, y0 + 14, 18, 18, d->color);
            DrawText(TextFormat("%d", inventory.count(id)), x + size - 4 - MeasureText(TextFormat("%d", inventory.count(id)), 14), y0 + 3, 14, RAYWHITE);
            const float left = ToolWear_Left(id); // a worn tool: a bar along the bottom, green to red
            if (left < 1.0f && inventory.count(id) > 0) {
                DrawRectangle(x + 4, y0 + size - 4, size - 8, 3, Fade(BLACK, 0.6f));
                DrawRectangle(x + 4, y0 + size - 4, (int)((size - 8) * left), 3, ColorLerp(RED, GREEN, left));
            }
        }
        const char* name = ShortName(id);
        int fs = 12;
        while (fs > 8 && MeasureText(name, fs) > size - 4) --fs;
        DrawText(name, x + (size - MeasureText(name, fs)) / 2, y0 + size - fs - 6, fs, RAYWHITE);
    }
}

void Hotbar_Save(FILE* f) {
    for (int i = 0; i < kHotbarSlots; ++i) fprintf(f, "hotbar %d %s\n", i, slots[i].empty() ? "-" : slots[i].c_str());
}

bool Hotbar_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "hotbar") != 0) return false;
    int slot = -1;
    char id[64] = {};
    if (fscanf(f, " %d %63s", &slot, id) == 2 && slot >= 0 && slot < kHotbarSlots) slots[slot] = strcmp(id, "-") == 0 || strcmp(id, kHotbarTerrainTool) == 0 ? "" : id; // the old "Hands" slot is an empty one now
    return true;
}
