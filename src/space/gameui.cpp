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
#include <algorithm>
#include <cmath>

#include "buildtool.h"
#include "itemicons.h"
#include "plots.h"
#include "worldtypes.h"
#include "helpmanual.h"
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

// An item's name along the bottom of its cell (#219): whole if it fits; otherwise as much as fits, and while hovered the
// whole name scrolls through the cell (a marquee, clipped to it), pausing at the start.
static void CellName(ImDrawList* dl, ImVec2 p, ImVec2 q, const std::string& name, bool hover) {
    const float room = q.x - p.x - 8.0f, y = q.y - 18.0f;
    const float tw = ImGui::CalcTextSize(name.c_str()).x;
    if (tw <= room) { dl->AddText(ImVec2(p.x + 4, y), IM_COL32(220, 220, 220, 255), name.c_str()); return; }
    dl->PushClipRect(ImVec2(p.x + 3, y - 1), ImVec2(q.x - 3, q.y), true);
    if (hover) {
        const float gap = 24.0f, speed = 35.0f, pause = 0.6f;    // px, px/s, s at the start of each pass
        const float period = (tw + gap) / speed + pause;
        const float t = std::max(0.0f, fmodf((float)ImGui::GetTime(), period) - pause);
        const float x = p.x + 4 - t * speed;
        dl->AddText(ImVec2(x, y), IM_COL32(255, 255, 255, 255), name.c_str());
        dl->AddText(ImVec2(x + tw + gap, y), IM_COL32(255, 255, 255, 255), name.c_str()); // the next pass coming in
    } else {
        size_t n = name.size();
        while (n > 1 && ImGui::CalcTextSize(name.c_str(), name.c_str() + n).x > room) --n;
        dl->AddText(nullptr, 0.0f, ImVec2(p.x + 4, y), IM_COL32(220, 220, 220, 255), name.c_str(), name.c_str() + n);
    }
    dl->PopClipRect();
}

// The stack on the cursor (#192): picked out of its cell, it follows the mouse until it is put down; whatever is still
// held when the inventory closes goes back in.
static InvSlot cursorStack;
static void ReturnCursorStack() {
    if (cursorStack.empty()) return;
    ++inventory.muted; // not a pickup (#208)
    inventory.addAmount(cursorStack.id, cursorStack.amount());
    --inventory.muted;
    cursorStack = InvSlot{};
}
// Puts up to `centi` of the cursor stack into a cell (empty, or the same item up to its stack limit).
static void DropInto(InvSlot& cell, int64_t centi) {
    if (cursorStack.empty()) return;
    if (!cell.empty() && cell.id != cursorStack.id) return;
    const int64_t room = (int64_t)StackLimit(cursorStack.id) * 100 - (cell.empty() ? 0 : cell.centi);
    const int64_t put = std::min({ centi, cursorStack.centi, room });
    if (put <= 0) return;
    if (cell.empty()) { cell.id = cursorStack.id; cell.centi = 0; }
    cell.centi += put;
    cursorStack.centi -= put;
    if (cursorStack.centi <= 0) cursorStack = InvSlot{};
}
// Shift+click: the stack to the other part (grid <-> hotbar), onto the same item first, then an empty cell.
static void QuickMove(std::vector<InvSlot>& grid, int i) {
    if (grid[i].empty()) return;
    const bool toHotbar = i < Inventory::kGrid;
    const int from = toHotbar ? Inventory::kGrid : 0, to = toHotbar ? Inventory::kSlots : Inventory::kGrid;
    InvSlot moving = grid[i];
    grid[i] = InvSlot{};
    const int64_t cap = (int64_t)StackLimit(moving.id) * 100;
    for (int pass = 0; pass < 2 && moving.centi > 0; ++pass)
        for (int k = from; k < to && moving.centi > 0; ++k) {
            InvSlot& c = grid[k];
            if (pass == 0 ? (c.empty() || c.id != moving.id) : !c.empty()) continue;
            if (c.empty()) { c.id = moving.id; c.centi = 0; }
            const int64_t put = std::min(moving.centi, cap - c.centi);
            c.centi += put; moving.centi -= put;
        }
    if (moving.centi > 0) { if (grid[i].empty()) grid[i] = moving; else { cursorStack = moving; ReturnCursorStack(); } } // no room: it stays
}

// One click on cell i of the player's inventory (#192): left picks up / puts down / merges / swaps (Shift: quick-move to the other
// part), right takes half or, holding a stack, puts one down.
static void CellClick(int i, bool right, bool shift) {
    auto& grid = inventory.grid();
    if (right) {
        if (cursorStack.empty() && !grid[i].empty()) { // half of it (whole units; the larger half stays)
            const int64_t half = (grid[i].centi / 100 / 2) * 100;
            if (half > 0) { cursorStack.id = grid[i].id; cursorStack.centi = half; grid[i].centi -= half; }
        } else if (!cursorStack.empty()) DropInto(grid[i], std::min<int64_t>(100, cursorStack.centi)); // one unit (or what fraction is left)
        return;
    }
    if (cursorStack.empty() && shift) QuickMove(grid, i);
    else if (cursorStack.empty()) { if (!grid[i].empty()) { cursorStack = grid[i]; grid[i] = InvSlot{}; } } // picked up
    else if (grid[i].empty() || grid[i].id == cursorStack.id) DropInto(grid[i], cursorStack.centi); // put down / merged (what doesn't fit stays held)
    else std::swap(grid[i], cursorStack); // swapped: the cell's stack is held now
}

bool GameUI_TestCursorStack() {
    auto& g = inventory.grid();
    auto amt = [&](int i, const char* id) { return !g[i].empty() && g[i].id == id ? g[i].centi / 100 : (int64_t)0; };
    ReturnCursorStack();
    inventory.clear();
    ++inventory.muted; inventory.add("stone", 10); inventory.add("oak_log", 3); --inventory.muted; // cells 0 and 1
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: inventory cursor: %s: NO\n", what); ok = false; } };
    CellClick(0, false, false); check(cursorStack.id == "stone" && cursorStack.centi == 1000 && g[0].empty(), "click picks the stack up");
    CellClick(5, false, false); check(amt(5, "stone") == 10 && cursorStack.empty(), "click puts it down");
    CellClick(5, true, false);  check(amt(5, "stone") == 5 && cursorStack.centi == 500, "right-click takes half");
    CellClick(6, true, false);  check(amt(6, "stone") == 1 && cursorStack.centi == 400, "right-click puts one down");
    CellClick(6, false, false); check(amt(6, "stone") == 5 && cursorStack.empty(), "click merges the rest");
    CellClick(1, false, false); CellClick(6, false, false); check(amt(6, "oak_log") == 3 && cursorStack.id == "stone" && cursorStack.centi == 500, "click on another item swaps");
    CellClick(1, false, false); check(amt(1, "stone") == 5 && cursorStack.empty(), "...and puts the swapped stack down");
    CellClick(6, false, true);  check(g[6].empty() && amt(Inventory::kGrid, "oak_log") == 3, "Shift-click moves to the hotbar");
    CellClick(5, false, false); ReturnCursorStack(); check(cursorStack.empty() && inventory.count("stone") == 10, "a held stack goes back in when the window closes");
    printf("bench: inventory cursor (#192): %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// The grid (#192): click a cell to pick its stack up onto the cursor, click a cell to put it down (merging the same item,
// swapping another); right-click takes half, or with a stack held puts one down; Shift+click moves a stack between the grid
// and the hotbar. Hover and press 1-0 to put the item on the hotbar. With a container open a click sends the stack into it.
static void DrawGrid() {
    auto& grid = inventory.grid();
    const float cell = 54.0f, gap = 4.0f;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (int i = 0; i < Inventory::kSlots; ++i) {
        // The 9 x 5 grid, then the hotbar's 10 cells in a row under it (#209: real slots; drop stacks into them).
        if (i == Inventory::kGrid) { ImGui::Dummy(ImVec2(0.0f, 2.0f)); ImGui::TextDisabled("Hotbar (keys 1-0)"); }
        else if (i < Inventory::kGrid ? i % Inventory::kCols != 0 : true) ImGui::SameLine(0.0f, gap);
        ImGui::PushID(i);
        const ImVec2 p = ImGui::GetCursorScreenPos();
        const float cw = i < Inventory::kGrid ? cell : (Inventory::kCols * cell - gap) / (float)Inventory::kHotbar; // the hotbar row as wide as the grid
        if (ImGui::InvisibleButton("cell", ImVec2(cw, cell))) {
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
            } else if (!GameUI_ContainerOpen()) CellClick(i, false, IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT));
        }
        if (!GameUI_ContainerOpen() && ImGui::IsItemClicked(ImGuiMouseButton_Right)) CellClick(i, true, false);
        const bool hover = ImGui::IsItemHovered();
        const InvSlot& s = grid[i];
        const ImVec2 q(p.x + cw, p.y + cell);
        dl->AddRectFilled(p, q, IM_COL32(30, 32, 36, 220), 4.0f);
        if (!s.empty()) {
            const ItemDef* d = itemDB.find(s.id);
            DrawItemIcon(dl, ImVec2(p.x + (cw - cell * 0.64f) * 0.5f, p.y), cell * 0.64f, s.id);
            CellName(dl, p, q, d ? d->name : s.id, hover);
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
                if (k >= 0) inventory.swapSlots(i, Inventory::kGrid + k); // into that hotbar slot (its stack comes here)
            }
        }
        if (i >= Inventory::kGrid) dl->AddText(ImVec2(p.x + 3, p.y + 1), IM_COL32(150, 150, 150, 255), TextFormat("%d", (i - Inventory::kGrid + 1) % 10));
        const ImU32 edge = hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 72, 78, 255);
        dl->AddRect(p, q, edge, 4.0f, 0, 1.0f);
        ImGui::PopID();
    }
    if (ImGui::Button("Sort")) { ReturnCursorStack(); inventory.sort(); }
    ImGui::SameLine();
    // Deposit (button or Q): items a nearby chest already holds go into it.
    const int nChests = uiWorld ? (int)Storage_NearbyChests(*uiWorld).size() : 0;
    ImGui::BeginDisabled(nChests == 0);
    if ((ImGui::Button("Deposit (Q)") || (nChests > 0 && IsKeyPressed(KEY_Q))) && uiWorld) {
        const int kinds = Storage_Deposit(*uiWorld);
        GameUI_Toast(kinds > 0 ? TextFormat("Deposited %d kind%s of item into nearby chests", kinds, kinds > 1 ? "s" : "") : "Nothing here matches what the nearby chests hold");
        if (kinds > 0) Audio_Play(Sfx::Pickup, 0.5f, 0.8f);
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::TextDisabled(!cursorStack.empty() ? "Click: put it down  right-click: one"  : nChests > 0 ? TextFormat("%d chest%s within reach", nChests, nChests > 1 ? "s" : "") : "Tools 1, blocks 256, resources 8192");
    if (!cursorStack.empty()) { // the held stack follows the mouse, over everything (#192)
        ImDrawList* fg = ImGui::GetForegroundDrawList();
        const ImVec2 m = ImGui::GetMousePos();
        DrawItemIcon(fg, ImVec2(m.x - 18.0f, m.y - 18.0f), 36.0f, cursorStack.id);
        if (KindOf(cursorStack.id) != ItemKind::Tool) fg->AddText(ImVec2(m.x + 12.0f, m.y + 8.0f), IM_COL32(255, 255, 255, 255), AmountText(cursorStack).c_str());
    }
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
            CellName(dl, p, q, d ? d->name : s.first, hover);
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
    if (!inventoryOpen) { ReturnCursorStack(); container.open = false; return; } // a stack still held goes back in (#192)
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
            if (ImGui::BeginTabItem("Help", nullptr, forceTab == "Help" ? ImGuiTabItemFlags_SetSelected : 0)) { // the manual (#354)
                HelpManual_DrawTab();
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

// Pickup feed (#208): "+N item" lines in the bottom right, newest at the bottom, each fading out ~3 s after its last
// pickup; more of an item still on screen adds to its line (digging shows one growing "+12.35 Dirt").
namespace {
    struct Pickup { std::string id; double amount; double t; };
    std::vector<Pickup> pickups;
    constexpr double kPickupShow = 3.0, kPickupFade = 0.7;
    constexpr size_t kPickupLines = 6;
    void OnPickup(const std::string& id, double amount) {
        const double now = GetTime();
        for (size_t i = 0; i < pickups.size(); ++i)
            if (pickups[i].id == id && now - pickups[i].t < kPickupShow) {
                Pickup p = pickups[i];
                p.amount += amount; p.t = now;
                pickups.erase(pickups.begin() + (long)i);
                pickups.push_back(p); // the freshest line goes to the bottom
                return;
            }
        pickups.push_back({ id, amount, now });
        if (pickups.size() > 12) pickups.erase(pickups.begin());
    }
    void DrawPickups() {
        const double now = GetTime();
        pickups.erase(std::remove_if(pickups.begin(), pickups.end(), [&](const Pickup& p) { return now - p.t > kPickupShow; }), pickups.end());
        const int w = GetScreenWidth(), h = GetScreenHeight(), line = 26, icon = 22;
        const size_t first = pickups.size() > kPickupLines ? pickups.size() - kPickupLines : 0;
        int y = h - 140 - (int)(pickups.size() - first) * line;
        for (size_t i = first; i < pickups.size(); ++i, y += line) {
            const Pickup& p = pickups[i];
            const float a = (float)std::clamp((kPickupShow - (now - p.t)) / kPickupFade, 0.0, 1.0);
            const ItemDef* def = itemDB.find(p.id);
            const double r = std::round(p.amount);
            const char* txt = fabs(p.amount - r) < 0.005 ? TextFormat("+%.0f %s", r, def ? def->name.c_str() : p.id.c_str())
                                                          : TextFormat("+%.2f %s", p.amount, def ? def->name.c_str() : p.id.c_str());
            const int tw = MeasureText(txt, 20), x = w - 30 - tw;
            DrawRectangle(x - icon - 14, y - 2, tw + icon + 22, line - 2, Fade(BLACK, 0.35f * a));
            ItemIcons_Draw(p.id, { (float)(x - icon - 8), (float)y, (float)icon, (float)icon }, Fade(WHITE, a));
            DrawText(txt, x, y + 2, 20, Fade(RAYWHITE, a));
        }
    }
}
void GameUI_InitPickupFeed() { inventory.onAdded = OnPickup; }

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
    DrawPickups(); // #208
    if (ship.piloted) {
        DrawText(TextFormat("SHIP  %.0f m/s", relSpeed(ship.velocity)), w - 220, 70, 20, SKYBLUE);
        if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
        return;
    }
    // Jetpack fuel bar (top right).
    DrawRectangle(w - 220, 70, 200, 12, Fade(BLACK, 0.5f));
    DrawRectangle(w - 220, 70, (int)(200 * player.jetpackFuel), 12, SKYBLUE);
    DrawText("JET", w - 260, 68, 16, RAYWHITE);
    const char* mode = player.flyMode ? "FLY (F)" : player.grounded ? "WALK" : player.hovering ? "HOVER (H)" : player.gravity >= 0.5f ? "AIR" : "SPACE";
    DrawText(TextFormat("%s   g %.1f   %.1f m/s", mode, player.gravity, relSpeed(player.velocity, player.grounded)), w - 220, 90, 16, RAYWHITE);
    DrawText(TextFormat("dampeners %s (B)", player.dampeners ? "on" : "off"), w - 220, 110, 16, player.dampeners ? RAYWHITE : ORANGE);
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
}
