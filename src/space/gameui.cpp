#include "gameui.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "audio.h"
#include "terraintool.h"
#include "survival.h"
#include "world.h"
#include "ship.h"

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

static const char* KindName(ItemKind k) { return k == ItemKind::Tool ? "tool" : k == ItemKind::Block ? "block" : "resource"; }

// Whole numbers plain; fractions to two places ("12.35"); thousands shortened ("8.1k").
static std::string AmountText(const InvSlot& s) {
    const double a = s.amount();
    if (a >= 1000.0) return TextFormat("%.1fk", a / 1000.0);
    if (s.centi % 100 == 0) return TextFormat("%lld", (long long)(s.centi / 100));
    return TextFormat("%.2f", a);
}

// The grid: click a cell to pick its stack up, click another to drop it there (swapping, or merging
// the same item); hover and press 1-0 to put the item on the hotbar.
static void DrawGrid() {
    static int held = -1;
    auto& grid = inventory.grid();
    const float cell = 54.0f, gap = 4.0f;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (int i = 0; i < Inventory::kSlots; ++i) {
        if (i % Inventory::kCols != 0) ImGui::SameLine(0.0f, gap);
        ImGui::PushID(i);
        const ImVec2 p = ImGui::GetCursorScreenPos();
        if (ImGui::InvisibleButton("cell", ImVec2(cell, cell))) {
            if (held < 0) { if (!grid[i].empty()) held = i; }
            else { inventory.swapSlots(held, i); held = -1; }
        }
        const bool hover = ImGui::IsItemHovered();
        const InvSlot& s = grid[i];
        const ImVec2 q(p.x + cell, p.y + cell);
        dl->AddRectFilled(p, q, IM_COL32(30, 32, 36, 220), 4.0f);
        if (!s.empty()) {
            const ItemDef* d = itemDB.find(s.id);
            const Color c = d ? d->color : GRAY;
            dl->AddRectFilled(ImVec2(p.x + 6, p.y + 6), ImVec2(q.x - 6, q.y - 20), IM_COL32(c.r, c.g, c.b, 255), 3.0f);
            const std::string name = d ? d->name.substr(0, 7) : s.id.substr(0, 7);
            dl->AddText(ImVec2(p.x + 4, q.y - 18), IM_COL32(220, 220, 220, 255), name.c_str());
            if (KindOf(s.id) != ItemKind::Tool) {
                const std::string amt = AmountText(s);
                const ImVec2 sz = ImGui::CalcTextSize(amt.c_str());
                dl->AddRectFilled(ImVec2(q.x - sz.x - 6, p.y + 2), ImVec2(q.x - 2, p.y + sz.y + 2), IM_COL32(0, 0, 0, 160), 2.0f);
                dl->AddText(ImVec2(q.x - sz.x - 4, p.y + 2), IM_COL32(255, 255, 255, 255), amt.c_str());
            }
            if (hover) {
                ImGui::SetTooltip("%s  (%s)\n%.2f of %d per slot\n%s\n(click to move; 1-0: hotbar slot)", d ? d->name.c_str() : s.id.c_str(),
                                  KindName(KindOf(s.id)), s.amount(), StackLimit(s.id), d ? d->description.c_str() : "");
                const int k = Hotbar_KeySlot();
                if (k >= 0) Hotbar_Bind(k, s.id);
            }
        }
        const ImU32 edge = i == held ? IM_COL32(255, 210, 80, 255) : hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 72, 78, 255);
        dl->AddRect(p, q, edge, 4.0f, 0, i == held ? 2.5f : 1.0f);
        ImGui::PopID();
    }
    if (held >= 0 && grid[held].empty()) held = -1;
    if (ImGui::Button("Sort")) { inventory.sort(); held = -1; }
    ImGui::SameLine();
    ImGui::TextDisabled(held >= 0 ? "Click a cell to drop the stack" : "Tools: 1 per slot, blocks 256, resources 8192");
}

// Hand crafting. In creative mode every recipe can be made and costs nothing.
static void DrawCrafting() {
    const bool creative = terrainTool.creative;
    if (creative) ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f), "Creative: crafting costs nothing");
    for (size_t i = 0; i < itemDB.recipes().size(); ++i) {
        const Recipe& r = itemDB.recipes()[i];
        const ItemDef* out = itemDB.find(r.output);
        const bool can = creative || inventory.canCraft(r);
        ImGui::PushID((int)i);
        ImGui::BeginDisabled(!can);
        if (ImGui::Button("Craft")) {
            if (creative) { if (inventory.addAmount(r.output, r.outputCount) < r.outputCount) Audio_Play(Sfx::Craft); }
            else if (inventory.craft(r)) Audio_Play(Sfx::Craft);
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
}

// Creative mode: free building, crafting and blueprints, and any item into the inventory.
static void DrawCreative() {
    ImGui::Checkbox("Creative mode (build, craft and place blueprints for free)", &terrainTool.creative);
    static char filter[64] = {};
    ImGui::SetNextItemWidth(220.0f);
    ImGui::InputTextWithHint("##filter", "search items", filter, sizeof(filter));
    ImGui::SameLine();
    if (ImGui::Button("Clear inventory")) inventory.clear();
    ImGui::Separator();
    std::string f = filter;
    for (char& c : f) c = (char)tolower(c);
    ImGui::BeginChild("creativeitems");
    int shown = 0;
    for (const ItemDef& d : itemDB.items()) {
        std::string name = d.name;
        for (char& c : name) c = (char)tolower(c);
        if (!f.empty() && name.find(f) == std::string::npos && d.id.find(f) == std::string::npos) continue;
        ImGui::PushID(d.id.c_str());
        ImGui::ColorButton("##c", ToImColor(d.color), ImGuiColorEditFlags_NoTooltip, ImVec2(14, 14));
        ImGui::SameLine();
        ImGui::Text("%-20s", d.name.c_str());
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("%s (%s, %d per slot)\n%s", d.id.c_str(), KindName(d.kind), StackLimit(d.id), d.description.c_str());
        ImGui::SameLine(260.0f);
        if (ImGui::SmallButton("+1")) inventory.add(d.id, 1);
        if (d.kind != ItemKind::Tool) {
            ImGui::SameLine();
            if (ImGui::SmallButton("+10")) inventory.add(d.id, 10);
            ImGui::SameLine();
            if (ImGui::SmallButton("stack")) inventory.add(d.id, StackLimit(d.id));
        }
        ImGui::SameLine();
        ImGui::TextDisabled("have %.2f", inventory.amount(d.id));
        ImGui::PopID();
        ++shown;
    }
    if (shown == 0) ImGui::TextDisabled("No items match");
    ImGui::EndChild();
}

void GameUI_Draw(Scene&) {
    if (!inventoryOpen) return;
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(1000, 480), ImGuiCond_Appearing);
    if (ImGui::Begin("Inventory (Tab)", &inventoryOpen, ImGuiWindowFlags_NoCollapse)) {
        ImGui::Columns(2, "inv", true);
        static bool widthSet = false;
        if (!widthSet) { ImGui::SetColumnWidth(0, Inventory::kCols * 58.0f + 24.0f); widthSet = true; }
        ImGui::TextUnformatted("Items");
        ImGui::Separator();
        DrawGrid();

        ImGui::NextColumn();
        if (ImGui::BeginTabBar("invtabs")) {
            if (ImGui::BeginTabItem("Crafting")) {
                DrawCrafting();
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Quests")) {
                // The main quest line: done tiers collapsed, the current one open, the next shown.
                const int cur = QuestBook_Tier();
                for (int i = 0; i < QuestBook_TierCount() && i <= cur + 1; ++i) {
                    ImGui::PushID(i);
                    const bool done = i < cur;
                    const char* name = QuestBook_TierName(i);
                    ImGui::SetNextItemOpen(i == cur, ImGuiCond_Appearing);
                    if (ImGui::CollapsingHeader(TextFormat("%s Tier %d: %s", done ? "[done]" : i == cur ? "[now]" : "[next]", i, name))) {
                        ImGui::TextWrapped("%s", QuestBook_TierText(i));
                        for (int k = 0; k < QuestBook_TaskCount(i); ++k) {
                            const bool kd = QuestBook_TaskDone(i, k);
                            ImGui::TextColored(kd ? ImVec4(0.45f, 0.85f, 0.45f, 1.0f) : ImVec4(0.95f, 0.95f, 0.95f, 1.0f), "%s %s", kd ? "[x]" : "[ ]", QuestBook_TaskText(i, k));
                        }
                    }
                    ImGui::PopID();
                }
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Creative")) {
                DrawCreative();
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
        }
        ImGui::Columns(1);
    }
    ImGui::End();
}

void GameUI_DrawHUD(Scene&) {
    const int w = GetScreenWidth();
    if (inventory.overflowed) { GameUI_Toast("Inventory full: some items were lost"); inventory.overflowed = false; }
    if (ship.piloted) {
        DrawText(TextFormat("SHIP  %.0f m/s", Vector3Length(ship.velocity)), w - 220, 70, 20, SKYBLUE);
        if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
        return;
    }
    // Jetpack fuel bar (top right).
    DrawRectangle(w - 220, 70, 200, 12, Fade(BLACK, 0.5f));
    DrawRectangle(w - 220, 70, (int)(200 * player.jetpackFuel), 12, SKYBLUE);
    DrawText("JET", w - 260, 68, 16, RAYWHITE);
    const char* mode = player.flyMode ? "FLY (F)" : player.grounded ? "WALK" : player.gravity >= 0.5f ? "AIR" : "SPACE";
    DrawText(TextFormat("%s   g %.1f   %.1f m/s", mode, player.gravity, Vector3Length(player.velocity)), w - 220, 90, 16, RAYWHITE);
    DrawText(TextFormat("dampeners %s (B)", player.dampeners ? "on" : "off"), w - 220, 110, 16, player.dampeners ? RAYWHITE : ORANGE);
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
}
