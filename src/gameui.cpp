#include "gameui.h"
#include "inventory.h"
#include "player.h"
#include "audio.h"
#include "world.h"

#include "imgui.h"
#include <string>

namespace {
    bool inventoryOpen = false;
    std::string toast;
    double toastUntil = 0.0;
}

bool GameUI_InventoryOpen() { return inventoryOpen; }
void GameUI_ToggleInventory() { inventoryOpen = !inventoryOpen; }
void GameUI_Toast(const char* text) { toast = text; toastUntil = GetTime() + 2.5; }

static ImVec4 ToImColor(Color c) { return ImVec4(c.r / 255.0f, c.g / 255.0f, c.b / 255.0f, 1.0f); }

void GameUI_Draw(Scene&) {
    if (!inventoryOpen) return;
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(620, 420), ImGuiCond_Appearing);
    if (ImGui::Begin("Inventory (Tab)", &inventoryOpen, ImGuiWindowFlags_NoCollapse)) {
        ImGui::Columns(2, "inv", true);
        ImGui::TextUnformatted("Items");
        ImGui::Separator();
        bool any = false;
        for (const ItemDef& d : itemDB.items()) {
            const int n = inventory.count(d.id);
            if (n <= 0) continue;
            any = true;
            ImGui::ColorButton(("##c" + d.id).c_str(), ToImColor(d.color), ImGuiColorEditFlags_NoTooltip, ImVec2(14, 14));
            ImGui::SameLine();
            ImGui::Text("%-18s x%d", d.name.c_str(), n);
            if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s", d.description.c_str());
        }
        if (!any) ImGui::TextDisabled("Empty. Dig with LMB to collect material.");

        ImGui::NextColumn();
        ImGui::TextUnformatted("Crafting");
        ImGui::Separator();
        for (size_t i = 0; i < itemDB.recipes().size(); ++i) {
            const Recipe& r = itemDB.recipes()[i];
            const ItemDef* out = itemDB.find(r.output);
            const bool can = inventory.canCraft(r);
            ImGui::PushID((int)i);
            ImGui::BeginDisabled(!can);
            if (ImGui::Button("Craft")) {
                if (inventory.craft(r)) Audio_Play(Sfx::Craft);
            }
            ImGui::EndDisabled();
            ImGui::SameLine();
            ImGui::Text("%d x %s", r.outputCount, out ? out->name.c_str() : r.output.c_str());
            std::string needs = "   needs ";
            for (size_t k = 0; k < r.inputs.size(); ++k) {
                const ItemDef* in = itemDB.find(r.inputs[k].first);
                needs += std::to_string(inventory.count(r.inputs[k].first)) + "/" + std::to_string(r.inputs[k].second) + " " + (in ? in->name : r.inputs[k].first);
                if (k + 1 < r.inputs.size()) needs += ", ";
            }
            ImGui::TextDisabled("%s", needs.c_str());
            ImGui::PopID();
        }
        ImGui::Columns(1);
    }
    ImGui::End();
}

void GameUI_DrawHUD(Scene&) {
    const int w = GetScreenWidth();
    // Jetpack fuel bar (top right).
    DrawRectangle(w - 220, 70, 200, 12, Fade(BLACK, 0.5f));
    DrawRectangle(w - 220, 70, (int)(200 * player.jetpackFuel), 12, SKYBLUE);
    DrawText("JET", w - 260, 68, 16, RAYWHITE);
    const char* mode = player.flyMode ? "FLY (F)" : player.grounded ? "WALK" : player.gravity > 0 ? "AIR" : "SPACE";
    DrawText(TextFormat("%s   g %.1f", mode, player.gravity), w - 220, 90, 16, RAYWHITE);
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
}
