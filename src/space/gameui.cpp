#include "gameui.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "audio.h"
#include "terraintool.h"
#include "survival.h"
#include "suit.h"
#include "toolwear.h"
#include "questtree.h"
#include "recipebook.h"
#include "storage.h"
#include "factory.h"
#include <map>
#include "world.h"
#include "ship.h"

#include "imgui.h"
#include <string>

#include "buildtool.h"
#include "itemicons.h"
#include "plots.h"
#include "worldtypes.h"
#include <raymath.h>

namespace {
    bool inventoryOpen = false;
    std::string toast;
    double toastUntil = 0.0;
    struct Container { bool open = false; uint64_t body = 0; uint32_t plot = 0; int part = -1; PartType type = PartType::Count; } container;
}

bool GameUI_InventoryOpen() { return inventoryOpen; }
void GameUI_ToggleInventory() { inventoryOpen = !inventoryOpen; if (!inventoryOpen) container.open = false; }
void GameUI_OpenContainer(uint64_t body, uint32_t plot, int part, int partType) {
    container = { true, body, plot, part, (PartType)partType };
    inventoryOpen = true;
}
bool GameUI_ContainerOpen() { return container.open && inventoryOpen; }
void GameUI_CloseContainer() { container.open = false; }

// An item's icon: its colour in a shape by kind - tools a handle and head, blocks a shaded cube,
// machines a framed box, resources a rounded lump. `dim` greys it (a recipe you can't make yet).
static void DrawItemIcon(ImDrawList* dl, ImVec2 p, float size, const std::string& id, bool dim = false) {
    if (ItemIcons_Has(id)) {
        // The painted icon (itemicons.h); dimmed for recipes you can't make yet.
        float uv0[2], uv1[2];
        ItemIcons_UV(id, uv0, uv1);
        const float m = size * 0.04f;
        dl->AddImage((ImTextureID)(intptr_t)ItemIcons_Texture().id, ImVec2(p.x + m, p.y + m), ImVec2(p.x + size - m, p.y + size - m),
                     ImVec2(uv0[0], uv0[1]), ImVec2(uv1[0], uv1[1]), dim ? IM_COL32(110, 110, 115, 150) : IM_COL32(255, 255, 255, 255));
        return;
    }
    const ItemDef* d = itemDB.find(id);
    Color c = d ? d->color : GRAY;
    if (dim) c = ColorLerp(c, Color{ 60, 62, 66, 255 }, 0.6f);
    auto col = [&](float k, int a = 255) { return IM_COL32((int)std::min(255.0f, c.r * k), (int)std::min(255.0f, c.g * k), (int)std::min(255.0f, c.b * k), a); };
    const float m = size * 0.16f, x0 = p.x + m, y0 = p.y + m, x1 = p.x + size - m, y1 = p.y + size - m;
    const float cx = (x0 + x1) * 0.5f, cy = (y0 + y1) * 0.5f, s = x1 - x0;
    const ImU32 edge = IM_COL32(0, 0, 0, dim ? 90 : 160);
    if (IsMachineItem(id)) {
        dl->AddRectFilled(ImVec2(x0, y0), ImVec2(x1, y1), col(0.75f), 3.0f);
        dl->AddRectFilled(ImVec2(x0 + s * 0.22f, y0 + s * 0.22f), ImVec2(x1 - s * 0.22f, y1 - s * 0.22f), col(1.15f), 2.0f);
        dl->AddRect(ImVec2(x0, y0), ImVec2(x1, y1), edge, 3.0f, 0, 1.5f);
        return;
    }
    switch (KindOf(id)) {
    case ItemKind::Tool: {
        // A diagonal handle with the coloured head at its top end.
        dl->AddLine(ImVec2(x0 + s * 0.15f, y1 - s * 0.1f), ImVec2(x1 - s * 0.3f, y0 + s * 0.3f), IM_COL32(150, 105, 60, dim ? 120 : 255), s * 0.13f);
        dl->AddQuadFilled(ImVec2(x1 - s * 0.62f, y0 + s * 0.05f), ImVec2(x1, y0 + s * 0.12f), ImVec2(x1 - s * 0.1f, y0 + s * 0.6f), ImVec2(x1 - s * 0.45f, y0 + s * 0.35f), col(1.0f));
        dl->AddQuad(ImVec2(x1 - s * 0.62f, y0 + s * 0.05f), ImVec2(x1, y0 + s * 0.12f), ImVec2(x1 - s * 0.1f, y0 + s * 0.6f), ImVec2(x1 - s * 0.45f, y0 + s * 0.35f), edge, 1.2f);
        break;
    }
    case ItemKind::Block: {
        // An isometric cube: light top, mid left, dark right.
        const ImVec2 t(cx, y0), l(x0, y0 + s * 0.27f), r(x1, y0 + s * 0.27f), mid(cx, y0 + s * 0.54f), bl(x0, y1 - s * 0.2f), br(x1, y1 - s * 0.2f), b(cx, y1 + s * 0.06f);
        dl->AddQuadFilled(t, r, mid, l, col(1.2f));
        dl->AddQuadFilled(l, mid, b, bl, col(0.85f));
        dl->AddQuadFilled(mid, r, br, b, col(0.6f));
        dl->AddQuad(t, r, mid, l, edge, 1.0f);
        dl->AddQuad(l, mid, b, bl, edge, 1.0f);
        dl->AddQuad(mid, r, br, b, edge, 1.0f);
        break;
    }
    default: {
        // A lump: an uneven blob with a highlight.
        const ImVec2 pts[7] = { ImVec2(cx - s * 0.42f, cy + s * 0.05f), ImVec2(cx - s * 0.25f, cy - s * 0.35f), ImVec2(cx + s * 0.15f, cy - s * 0.42f),
                                ImVec2(cx + s * 0.45f, cy - s * 0.1f), ImVec2(cx + s * 0.38f, cy + s * 0.32f), ImVec2(cx - s * 0.02f, cy + s * 0.45f), ImVec2(cx - s * 0.35f, cy + s * 0.35f) };
        dl->AddConvexPolyFilled(pts, 7, col(1.0f));
        dl->AddPolyline(pts, 7, edge, ImDrawFlags_Closed, 1.2f);
        dl->AddCircleFilled(ImVec2(cx - s * 0.12f, cy - s * 0.15f), s * 0.09f, IM_COL32(255, 255, 255, dim ? 40 : 90));
        break;
    }
    }
}

// Whether a container takes an item from the inventory (a chest: raw materials, not tools or machines).
static bool ContainerTakes(const std::string& id) {
    if (!Factory::Accepts(container.type, id)) return false;
    if (container.type == PartType::Chest && (KindOf(id) == ItemKind::Tool || id == "telescope" || IsMachineItem(id))) return false;
    return true;
}
static std::string forceTab; // the tab to select when the window next draws
static Scene* uiWorld = nullptr; // the world the window was last drawn for (nearby chests)
void GameUI_ShowTab(const char* tab) { inventoryOpen = true; forceTab = tab; }
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
            if (GameUI_ContainerOpen() && !grid[i].empty()) {
                // A container is open: the stack goes into it (what doesn't fit comes back next tick).
                const std::string id = grid[i].id;
                const int n = (int)grid[i].amount();
                if (!ContainerTakes(id)) GameUI_Toast(TextFormat("The %s doesn't take %s", PartName(container.type), itemDB.find(id) ? itemDB.find(id)->name.c_str() : id.c_str()));
                else if (n > 0 && inventory.remove(id, n)) {
                    // Shift+click: into the fuel slot (#194; fuel that is also a recipe input otherwise fills only an empty fuel slot).
                    if (IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT)) Factory::InsertFuel(container.body, container.plot, container.part, id, n);
                    else Factory::Insert(container.body, container.plot, container.part, id, n);
                    Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
                }
                held = -1;
            } else if (held < 0) { if (!grid[i].empty()) held = i; }
            else { inventory.swapSlots(held, i); held = -1; }
        }
        const bool hover = ImGui::IsItemHovered();
        const InvSlot& s = grid[i];
        const ImVec2 q(p.x + cell, p.y + cell);
        dl->AddRectFilled(p, q, IM_COL32(30, 32, 36, 220), 4.0f);
        if (!s.empty()) {
            const ItemDef* d = itemDB.find(s.id);
            DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y), cell * 0.64f, s.id);
            const std::string name = d ? d->name.substr(0, 7) : s.id.substr(0, 7);
            dl->AddText(ImVec2(p.x + 4, q.y - 18), IM_COL32(220, 220, 220, 255), name.c_str());
            if (KindOf(s.id) != ItemKind::Tool) {
                const std::string amt = AmountText(s);
                const ImVec2 sz = ImGui::CalcTextSize(amt.c_str());
                dl->AddRectFilled(ImVec2(q.x - sz.x - 6, p.y + 2), ImVec2(q.x - 2, p.y + sz.y + 2), IM_COL32(0, 0, 0, 160), 2.0f);
                dl->AddText(ImVec2(q.x - sz.x - 4, p.y + 2), IM_COL32(255, 255, 255, 255), amt.c_str());
            }
            if (hover) {
                ImGui::SetTooltip("%s  (%s)\n%.2f of %d per slot\n%s\n(%s; 1-0: hotbar slot)", d ? d->name.c_str() : s.id.c_str(),
                                  KindName(KindOf(s.id)), s.amount(), StackLimit(s.id), d ? d->description.c_str() : "",
                                  GameUI_ContainerOpen() ? TextFormat("click to put it in the %s", PartName(container.type)) : "click to move");
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
    // Deposit (button or Q): items a nearby chest already holds go into it.
    const int nChests = uiWorld ? (int)Storage_NearbyChests(*uiWorld).size() : 0;
    ImGui::BeginDisabled(nChests == 0);
    if ((ImGui::Button("Deposit (Q)") || (nChests > 0 && IsKeyPressed(KEY_Q))) && uiWorld) {
        const int kinds = Storage_Deposit(*uiWorld);
        GameUI_Toast(kinds > 0 ? TextFormat("Deposited %d kind%s of item into nearby chests", kinds, kinds > 1 ? "s" : "") : "Nothing here matches what the nearby chests hold");
        if (kinds > 0) Audio_Play(Sfx::Pickup, 0.5f, 0.8f);
        held = -1;
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::TextDisabled(held >= 0 ? "Click a cell to drop the stack" : nChests > 0 ? TextFormat("%d chest%s within reach", nChests, nChests > 1 ? "s" : "") : "Tools 1, blocks 256, resources 8192");
}

// Hand crafting: the recipes you know (every ingredient discovered, and any unlock met), searchable,
// optionally only what you can make now. Nearby chests count as if their contents were in the
// inventory: a craft short of something pulls it from them first (it arrives next tick).
// In creative mode every recipe is known and costs nothing.
static void DrawCrafting() {
    const bool creative = terrainTool.creative;
    if (creative) ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f), "Creative: crafting costs nothing");
    ToolWear_DrawRepairUI();
    static char query[64] = {};
    static bool craftableOnly = false;
    static int pending = -1;           // a craft waiting for items pulled from chests
    static double pendingUntil = 0.0;
    ImGui::SetNextItemWidth(220.0f);
    ImGui::InputTextWithHint("##recipesearch", "search recipes", query, sizeof(query));
    ImGui::SameLine();
    ImGui::Checkbox("Only what I can make", &craftableOnly);
    std::string q = query;
    for (char& c : q) c = (char)tolower((unsigned char)c);
    // What nearby chests hold, once for this frame.
    std::map<std::string, int> stored;
    if (uiWorld)
        for (const NearChest& c : Storage_NearbyChests(*uiWorld))
            for (const auto& s : Factory::ChestContents(c.body, c.plot, c.part)) stored[s.first] += s.second;
    auto have = [&](const std::string& id) { int n = inventory.count(id); for (const auto& s : stored) if (ItemMatches(id, s.first)) n += s.second; return n; };
    auto canMake = [&](const Recipe& r) { if (creative) return true; for (const auto& in : r.inputs) if (have(in.first) < in.second) return false; return true; };
    // A pending craft completes once the pulled items have arrived.
    if (pending >= 0 && pending < (int)itemDB.recipes().size()) {
        const Recipe& r = itemDB.recipes()[pending];
        if (inventory.canCraft(r)) { if (inventory.craft(r)) Audio_Play(Sfx::Craft); pending = -1; }
        else if (GetTime() > pendingUntil) { GameUI_Toast("The chests no longer hold enough"); pending = -1; }
    }
    // A grid of the outputs' icons; hover for name, cost and description, click to craft. Recipes
    // you can't make yet are greyed.
    int hidden = 0, shown = 0;
    ImGui::BeginChild("recipes");
    const float cell = 60.0f, gap = 6.0f;
    const int cols = std::max(1, (int)((ImGui::GetContentRegionAvail().x + gap) / (cell + gap)));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (size_t i = 0; i < itemDB.recipes().size(); ++i) {
        const Recipe& r = itemDB.recipes()[i];
        if (!RecipeBook_Known(r)) { ++hidden; continue; }
        if (!RecipeBook_Matches(r, q)) continue;
        const bool can = canMake(r);
        if (craftableOnly && !can) continue;
        if (shown % cols != 0) ImGui::SameLine(0.0f, gap);
        ++shown;
        const ItemDef* out = itemDB.find(r.output);
        ImGui::PushID((int)i);
        const ImVec2 p = ImGui::GetCursorScreenPos(), qd(p.x + cell, p.y + cell);
        const bool clicked = ImGui::InvisibleButton("recipe", ImVec2(cell, cell));
        const bool hover = ImGui::IsItemHovered();
        dl->AddRectFilled(p, qd, can ? IM_COL32(34, 40, 36, 230) : IM_COL32(30, 31, 34, 200), 4.0f);
        DrawItemIcon(dl, ImVec2(p.x + 4, p.y + 2), cell - 8, r.output, !can);
        if (r.outputCount > 1) {
            const char* n = TextFormat("x%d", r.outputCount);
            dl->AddText(ImVec2(qd.x - ImGui::CalcTextSize(n).x - 4, qd.y - 17), IM_COL32(255, 255, 255, can ? 255 : 130), n);
        }
        dl->AddRect(p, qd, hover ? IM_COL32(220, 220, 220, 255) : can ? IM_COL32(90, 150, 95, 255) : IM_COL32(64, 66, 72, 255), 4.0f, 0, hover ? 2.0f : 1.0f);
        if (hover) {
            ImGui::BeginTooltip();
            ImGui::Text("%s%s", out ? out->name.c_str() : r.output.c_str(), r.outputCount > 1 ? TextFormat("  x%d", r.outputCount) : "");
            ImGui::TextDisabled("%s", KindName(KindOf(r.output)));
            if (out && !out->description.empty()) { ImGui::PushTextWrapPos(320.0f); ImGui::TextUnformatted(out->description.c_str()); ImGui::PopTextWrapPos(); }
            ImGui::Separator();
            ImGui::TextUnformatted("Cost:");
            for (const auto& in : r.inputs) {
                const ItemDef* id = itemDB.find(in.first);
                const int inv = inventory.count(in.first), all = have(in.first);
                const bool enough = creative || all >= in.second;
                ImGui::TextColored(enough ? ImVec4(0.55f, 0.9f, 0.55f, 1.0f) : ImVec4(0.95f, 0.45f, 0.4f, 1.0f), "  %d %s   (have %d%s)", in.second,
                                   id ? id->name.c_str() : in.first.c_str(), all, all > inv ? TextFormat(", %d in chests", all - inv) : "");
            }
            ImGui::TextDisabled(can ? (pending >= 0 ? "Waiting for items from the chests..." : "Click to craft") : "Missing items");
            ImGui::EndTooltip();
        }
        if (clicked && can && pending < 0) {
            if (creative) { if (inventory.addAmount(r.output, r.outputCount) < r.outputCount) Audio_Play(Sfx::Craft); }
            else if (inventory.canCraft(r)) { if (inventory.craft(r)) Audio_Play(Sfx::Craft); }
            else if (uiWorld) {
                bool ok = true;
                for (const auto& in : r.inputs) ok = Storage_Pull(*uiWorld, in.first, in.second - inventory.count(in.first)) && ok;
                if (ok) { pending = (int)i; pendingUntil = GetTime() + 2.0; }
            }
        }
        ImGui::PopID();
    }
    if (shown == 0) ImGui::TextDisabled(craftableOnly ? "Nothing you can make right now" : "No known recipe matches");
    if (hidden > 0) ImGui::TextDisabled("%d recipes still to discover: find new materials (and some need more)", hidden);
    ImGui::EndChild();
}

// A machine's or chest's screen: status, progress and charge, then its slots (input and output for
// machines, 16 for a chest). Click a slot to take its stack; click an inventory cell to put one in.
static void DrawSlots(const std::vector<std::pair<std::string, int>>& slots, int which, int cols) { // which: 0 inputs, 1 outputs, 2 fuel
    const bool output = which == 1;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float cell = 54.0f, gap = 4.0f;
    for (int i = 0; i < (int)slots.size(); ++i) {
        if (i % cols != 0) ImGui::SameLine(0.0f, gap);
        ImGui::PushID(i + which * 1000);
        const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + cell, p.y + cell);
        const bool clicked = ImGui::InvisibleButton("slot", ImVec2(cell, cell));
        const bool hover = ImGui::IsItemHovered();
        dl->AddRectFilled(p, q, which == 2 ? IM_COL32(52, 30, 20, 220) : output ? IM_COL32(40, 36, 28, 220) : IM_COL32(30, 32, 36, 220), 4.0f);
        const auto& s = slots[i];
        if (s.second > 0) {
            DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y), cell * 0.64f, s.first);
            const ItemDef* d = itemDB.find(s.first);
            const std::string name = d ? d->name.substr(0, 7) : s.first.substr(0, 7);
            dl->AddText(ImVec2(p.x + 4, q.y - 18), IM_COL32(220, 220, 220, 255), name.c_str());
            const char* amt = TextFormat("%d", s.second);
            const ImVec2 sz = ImGui::CalcTextSize(amt);
            dl->AddRectFilled(ImVec2(q.x - sz.x - 6, p.y + 2), ImVec2(q.x - 2, p.y + sz.y + 2), IM_COL32(0, 0, 0, 160), 2.0f);
            dl->AddText(ImVec2(q.x - sz.x - 4, p.y + 2), IM_COL32(255, 255, 255, 255), amt);
            if (hover) ImGui::SetTooltip("%s x%d\n%s\n(click to take)", d ? d->name.c_str() : s.first.c_str(), s.second, d ? d->description.c_str() : "");
            if (clicked) { Factory::TakeSlot(container.body, container.plot, container.part, i, which); Audio_Play(Sfx::Pickup, 0.4f, 1.2f); }
        }
        dl->AddRect(p, q, hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 72, 78, 255), 4.0f, 0, 1.0f);
        ImGui::PopID();
    }
}

static void DrawContainer() {
    const MachineView v = Factory::View(container.body, container.plot, container.part);
    ImGui::Text("%s", PartName(container.type));
    ImGui::SameLine();
    if (v.valid) ImGui::TextColored(v.status == MachineStatus::Working ? ImVec4(0.5f, 0.9f, 0.5f, 1.0f) : ImVec4(0.9f, 0.8f, 0.5f, 1.0f), "  %s", MachineStatusName(v.status));
    ImGui::SameLine(ImGui::GetContentRegionAvail().x - 60.0f);
    if (ImGui::Button("Close")) { container.open = false; return; }
    ImGui::Separator();
    if (!v.valid) { ImGui::TextDisabled("..."); return; }
    if (v.duration > 0) ImGui::ProgressBar((float)v.progress / (float)v.duration, ImVec2(300.0f, 0.0f), TextFormat("%d%%", 100 * v.progress / v.duration));
    if (v.energyMax > 0) ImGui::ProgressBar((float)v.energy / (float)v.energyMax, ImVec2(300.0f, 0.0f), TextFormat("%d / %d EU", v.energy, v.energyMax));
    else if (v.energy > 0 && container.type != PartType::Chest) ImGui::Text("Fuel / charge: %d", v.energy);
    const bool chest = container.type == PartType::Chest;
    ImGui::TextUnformatted(chest ? "Storage" : "Input");
    DrawSlots(v.inSlots, 0, chest ? 8 : 4);
    if (v.hasFuelSlot) { // #194: fuel has its own slot (a recipe input that is also fuel is never burnt)
        ImGui::TextUnformatted("Fuel");
        DrawSlots({ v.fuelSlot }, 2, 1);
    }
    if (v.waterMax > 0) { // #229: boilers need water and heat up before full steam
        ImGui::ProgressBar((float)v.water / (float)v.waterMax, ImVec2(300.0f, 0.0f), TextFormat("Water %d / %d L", v.water, v.waterMax));
        ImGui::ProgressBar((float)v.heat / (float)std::max(1, v.heatMax), ImVec2(300.0f, 0.0f), TextFormat("Heat %d%%", 100 * v.heat / std::max(1, v.heatMax)));
    }
    if (!chest) {
        ImGui::TextUnformatted("Output");
        DrawSlots(v.outSlots, 1, 4);
    }
    ImGui::Spacing();
    if (ImGui::Button("Put in all (I)") || IsKeyPressed(KEY_I)) {
        const auto all = inventory.all(); // copy: we remove while walking it
        int kinds = 0;
        for (const auto& kv : all) {
            if (kv.second <= 0 || !ContainerTakes(kv.first)) continue;
            bool bound = false;
            for (int s = 0; s < 10; ++s) bound = bound || Hotbar_Slot(s) == kv.first;
            if (chest && bound) continue; // keep what's on the hotbar
            if (inventory.remove(kv.first, kv.second)) { Factory::Insert(container.body, container.plot, container.part, kv.first, kv.second); ++kinds; }
        }
        if (kinds > 0) Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
    }
    ImGui::SameLine();
    if (ImGui::Button("Take all (O)") || IsKeyPressed(KEY_O)) {
        if (chest) for (int i = 0; i < (int)v.inSlots.size(); ++i) { if (v.inSlots[i].second > 0) Factory::TakeSlot(container.body, container.plot, container.part, i, false); }
        else Factory::Extract(container.body, container.plot, container.part);
        Audio_Play(Sfx::Pickup, 0.4f, 1.2f);
    }
    ImGui::TextDisabled(chest ? "Click an inventory stack to store it; click a slot to take it" : "Click an inventory stack to load it (Shift+click: into the fuel slot); click a slot to take it");
}

// The open container must still exist (not broken) and be within reach.
static void CheckContainer(Scene& world) {
    if (!container.open) return;
    Planetoid* b = Planetoid::find(container.body);
    Plot* pl = b && b->plots ? b->plots->byId(container.plot) : nullptr;
    if (!pl || container.part >= (int)pl->parts.size() || pl->parts[container.part].type != container.type) { container.open = false; return; }
    const Part& part = pl->parts[container.part];
    const Vector3 c = pl->toBody({ (part.origin.x + 0.5f) * kCellSize, (part.origin.y + 0.5f) * kCellSize, (part.origin.z + 0.5f) * kCellSize });
    const Vector3 me = b->toLocal((player.position - world.renderOrigin).toFloat());
    if (Vector3Distance(c, me) * b->scale > 9.0f) { container.open = false; GameUI_Toast("Too far from it"); }
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
        if (IsWildcardItem(d.id)) continue; // tags, not items (#269)
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

void GameUI_Draw(Scene& world) {
    uiWorld = &world;
    RecipeBook_Update(); // new recipes: a toast
    if (!inventoryOpen) { container.open = false; return; }
    CheckContainer(world);
    const ImGuiIO& io = ImGui::GetIO();
    ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.5f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    ImGui::SetNextWindowSize(ImVec2(1320, 640), ImGuiCond_Appearing);
    if (ImGui::Begin("Inventory (Tab)", &inventoryOpen, ImGuiWindowFlags_NoCollapse)) {
        ImGui::Columns(2, "inv", true);
        static bool widthSet = false;
        if (!widthSet) { ImGui::SetColumnWidth(0, Inventory::kCols * 58.0f + 24.0f); widthSet = true; }
        ImGui::TextUnformatted("Items");
        ImGui::Separator();
        DrawGrid();

        ImGui::NextColumn();
        if (container.open) DrawContainer();
        else if (ImGui::BeginTabBar("invtabs")) {
            if (ImGui::BeginTabItem("Crafting", nullptr, forceTab == "Crafting" ? ImGuiTabItemFlags_SetSelected : 0)) {
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
            if (ImGui::BeginTabItem("Quest tree", nullptr, forceTab == "Quest tree" ? ImGuiTabItemFlags_SetSelected : 0)) {
                QuestTree_DrawTab();
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Suit", nullptr, forceTab == "Suit" ? ImGuiTabItemFlags_SetSelected : 0)) {
                Suit_DrawTab();
                ImGui::EndTabItem();
            }
            if (ImGui::BeginTabItem("Creative")) {
                DrawCreative();
                ImGui::EndTabItem();
            }
            ImGui::EndTabBar();
            forceTab.clear();
        }
        ImGui::Columns(1);
    }
    ImGui::End();
}

void GameUI_DrawHUD(Scene& world) {
    // Speeds against the body you are near (its co-rotating surface frame, as the rocket HUD's "surface" speed): the
    // world velocity includes the planet's spin and orbit, which read as 60-90 m/s while hovering (#239).
    // On the ground the player's velocity is already relative to the surface (Player::stepPhysics): only
    // world-frame velocities (airborne, fly mode, the ship) take the surface's motion out.
    auto relSpeed = [&](Vector3 v, bool surfaceRelative = false) {
        if (surfaceRelative) return Vector3Length(v);
        const Vector3 at = (player.position - world.renderOrigin).toFloat();
        Planetoid* body = player.nearest ? player.nearest : NearestPlanetoid(world, at); // fly mode never sets nearest
        return body ? Vector3Length(Vector3Subtract(v, body->surfaceVelocity(at))) : Vector3Length(v);
    };
    const int w = GetScreenWidth();
    if (player.speedVignette) {
        // Speed vignette (#168, comfort): the screen edges darken with speed against the ground (from ~8 m/s) and with
        // fast turning (from ~120 deg/s), eased in and out; at most half dark.
        static Vector3 lastFwd{}; static float strength = 0.0f;
        const float dt = std::max(GetFrameTime(), 1e-4f);
        const Vector3 fwd = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        const float turn = Vector3LengthSqr(lastFwd) > 0.5f ? acosf(std::clamp(Vector3DotProduct(fwd, lastFwd), -1.0f, 1.0f)) * RAD2DEG / dt : 0.0f;
        lastFwd = fwd;
        const float speed = relSpeed(player.velocity, player.grounded);
        const float want = std::clamp(std::clamp((speed - 8.0f) / 20.0f, 0.0f, 1.0f) + std::clamp((turn - 120.0f) / 240.0f, 0.0f, 1.0f), 0.0f, 1.0f);
        strength += (want - strength) * (1.0f - expf(-dt * 4.0f));
        if (strength > 0.01f) {
            const int h = GetScreenHeight(), band = (int)(0.22f * std::min(w, h));
            const Color dark = Fade(BLACK, 0.5f * strength), clear = Fade(BLACK, 0.0f);
            DrawRectangleGradientH(0, 0, band, h, dark, clear);
            DrawRectangleGradientH(w - band, 0, band, h, clear, dark);
            DrawRectangleGradientV(0, 0, w, band, dark, clear);
            DrawRectangleGradientV(0, h - band, w, band, clear, dark);
        }
    }
    if (inventory.overflowed) { GameUI_Toast("Inventory full: some items were lost"); inventory.overflowed = false; }
    if (ship.piloted) {
        DrawText(TextFormat("SHIP  %.0f m/s", relSpeed(ship.velocity)), w - 220, 70, 20, SKYBLUE);
        if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
        return;
    }
    // Jetpack fuel bar (top right).
    DrawRectangle(w - 220, 70, 200, 12, Fade(BLACK, 0.5f));
    DrawRectangle(w - 220, 70, (int)(200 * player.jetpackFuel), 12, SKYBLUE);
    DrawText("JET", w - 260, 68, 16, RAYWHITE);
    const char* mode = player.flyMode ? "FLY (F)" : player.grounded ? "WALK" : player.gravity >= 0.5f ? "AIR" : "SPACE";
    DrawText(TextFormat("%s   g %.1f   %.1f m/s", mode, player.gravity, relSpeed(player.velocity, player.grounded)), w - 220, 90, 16, RAYWHITE);
    DrawText(TextFormat("dampeners %s (B)", player.dampeners ? "on" : "off"), w - 220, 110, 16, player.dampeners ? RAYWHITE : ORANGE);
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
}
