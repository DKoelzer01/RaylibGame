#include "gameui.h"
#include "machineact.h"
#include "playerinput.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "audio.h"
#include "terraintool.h"
#include "animals.h"
#include "weather.h"
#include "survival.h"
#include "suit.h"
#include "toolwear.h"
#include "questtree.h"
#include "recipebook.h"
#include "net.h"
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
#include "menus.h"
#include "playerrecord.h"
#include "combat.h"
#include "playersession.h"
#include "avatar.h"
#include "netclient.h"
#include "gamesave.h"
#include "slotrules.h"
#include <raymath.h>

namespace {
    bool inventoryOpen = false;
    std::string toast;
    double toastUntil = 0.0;
    struct Container { bool open = false; uint64_t body = 0; uint32_t plot = 0; int part = -1; PartType type = PartType::Count; } container;
    // The game window (#821): the screen shown; Quests' view; the database pop-up (#826: alt-click an item).
    GameScreen screen = GameScreen::Inventory;
    bool questTree = false;
    bool databaseOpen = false;
}

bool GameUI_InventoryOpen() { return inventoryOpen; }
void GameUI_ToggleInventory() { inventoryOpen = !inventoryOpen; if (!inventoryOpen) { container.open = false; databaseOpen = false; } }
void GameUI_ShowScreen(GameScreen s) { inventoryOpen = true; screen = s; }
// (#826) Alt + click on an item (a cell, a slot, a recipe): its database page opens in the pop-up, instead of the click. True
// when it opened (the caller skips its own click).
static bool OpenIfAlt(const std::string& id) {
    if (id.empty() || !ImGui::GetIO().KeyAlt) return false;
    HelpManual_OpenItem(id);
    databaseOpen = true;
    return true;
}
GameScreen GameUI_Screen() { return screen; }
bool GameUI_TextInputActive() { return ImGui::GetCurrentContext() && ImGui::GetIO().WantTextInput; }
float GameUI_Scale() {
    // Settings' interface size (0: auto, by the screen's height: 1080 p is 100%).
    if (settings.uiScale > 0.0f) return settings.uiScale;
    return std::clamp(GetScreenHeight() / 1080.0f, 0.75f, 2.0f);
}
void GameUI_OpenContainer(uint64_t body, uint32_t plot, int part, int partType) {
    container = { true, body, plot, part, (PartType)partType };
    inventoryOpen = true;
}
bool GameUI_ContainerOpen() { return container.open && inventoryOpen; }
void GameUI_CloseContainer() { container.open = false; }

// An item's icon: its colour in a shape by kind - tools a handle and head, blocks a shaded cube,
// machines a framed box, resources a rounded lump. `dim` greys it (a recipe you can't make yet).
// A hot or radiant item's glow (#580: Item_Glow), drawn under its icon: the cell washed in the glow's colour and a
// bright rim, both as strong as the glow (a fresh ingot yellow-white, a cooling one dull red).
static void DrawGlow(ImDrawList* dl, ImVec2 p, ImVec2 q, const ItemGlow& g) {
    if (g.strength <= 0.02f) return;
    dl->AddRectFilled(p, q, IM_COL32(g.color.r, g.color.g, g.color.b, (int)(140.0f * g.strength)), 4.0f);
    dl->AddRect(p, q, IM_COL32(g.color.r, g.color.g, g.color.b, (int)(235.0f * g.strength)), 4.0f, 0, 2.0f);
}

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

// Whether a container takes an item from the inventory: a chest takes every item (#568: tools and machines too; tool
// wear is kept per item kind, not per stack, so nothing is lost), a machine what it burns or processes.
static bool ContainerTakes(const std::string& id) { return Factory::Accepts(container.type, id); }
static Scene* uiWorld = nullptr; // the world the window was last drawn for (nearby chests)
void GameUI_ShowTab(const char* tab) { // (the older tab names, from the tests and tutorials)
    databaseOpen = false; // (a named screen: the pop-up only when asked for, below)
    const std::string t = tab;
    if (t == "Crafting") GameUI_ShowScreen(GameScreen::Crafting);
    else if (t == "Quest tree" || t == "Quests") { GameUI_ShowScreen(GameScreen::Quests); questTree = t == "Quest tree"; }
    else if (t == "Help" || t == "Database") { if (!inventoryOpen) GameUI_ShowScreen(GameScreen::Inventory); databaseOpen = true; }
    else if (t == "Creative" || t == "Catalog") GameUI_ShowScreen(terrainTool.creative ? GameScreen::Catalog : GameScreen::Inventory);
    else GameUI_ShowScreen(GameScreen::Inventory); // "Suit": beside the inventory
}
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
    inventory.addAmount(cursorStack.id, cursorStack.amount(), cursorStack.heat); // still as hot (#580)
    --inventory.muted;
    cursorStack = InvSlot{};
}
// Getting rid of things when the grid is full (#564). The bin deletes the held stack; the last one binned can be taken
// back out until another goes in or the window closes. A click outside the window drops the held stack (right-click:
// one) on the ground at the feet as a pile, picked up again by walking back to it after leaving (Combat_Drop).
static InvSlot binned;
static void BinClick() {
    if (!cursorStack.empty()) {
        binned = cursorStack; cursorStack = InvSlot{};
        const ItemDef* d = itemDB.find(binned.id);
        GameUI_Toast(TextFormat("Deleted %s (click the bin to take it back)", d ? d->name.c_str() : binned.id.c_str()));
    } else if (!binned.empty()) { cursorStack = binned; binned = InvSlot{}; }
}
static void DropCursor(int64_t centi) {
    if (cursorStack.empty()) return;
    InvSlot s = cursorStack;
    s.centi = std::min(centi, cursorStack.centi);
    if (!Combat_Drop(LocalSession(), s)) return;
    cursorStack.centi -= s.centi;
    if (cursorStack.centi <= 0) cursorStack = InvSlot{};
    Audio_Play(Sfx::Pickup, 0.4f, 0.6f);
}

// ---- The slot grids the item rules act on (#862: slotrules.h, Minecraft's rules + Mouse Tweaks') ----
// The player's cells (grid, hotbar, off hand). Moves among them change no totals, so they stay here even when joined; what
// goes to or comes from a container is a request (MachineAct_*), what is dropped a pile at the feet (Combat_Drop).
struct PlayerGrid : SlotGrid {
    int size() const override { return Inventory::kCells; }
    InvSlot at(int i) const override { return inventory.grid()[i]; }
    int64_t room(int i, const std::string& id) const override {
        const InvSlot& c = inventory.grid()[i];
        if (!c.empty() && c.id != id) return 0;
        return (int64_t)StackLimit(id) * 100 - (c.empty() ? 0 : c.centi);
    }
    int64_t put(int i, const InvSlot& held, int64_t centi) override {
        InvSlot& c = inventory.grid()[i];
        const int64_t k = std::min(centi, room(i, held.id));
        if (k <= 0) return 0;
        if (c.empty()) { c.id = held.id; c.centi = 0; c.heat = 0.0f; }
        c.heat = (float)((c.heat * c.centi + held.heat * k) / (double)(c.centi + k)); // by amount, as merges do (#580)
        c.centi += k;
        return k;
    }
    int64_t take(int i, int64_t centi, InvSlot& held) override {
        InvSlot& c = inventory.grid()[i];
        if (c.empty() || (!held.empty() && held.id != c.id)) return 0;
        const int64_t k = std::min(centi, c.centi);
        if (held.empty()) { held.id = c.id; held.centi = 0; held.heat = c.heat; }
        held.heat = (float)((held.heat * held.centi + c.heat * k) / (double)(held.centi + k));
        held.centi += k;
        c.centi -= k;
        if (c.centi <= 0) c = InvSlot{};
        return k;
    }
    // To the open container (machines hold whole units), else between the grid and the hotbar (the off hand: to the grid).
    void moveOut(int i, int64_t centi) override {
        InvSlot& c = inventory.grid()[i];
        if (c.empty()) return;
        if (GameUI_ContainerOpen()) {
            const std::string id = c.id;
            if (!ContainerTakes(id)) { GameUI_Toast(TextFormat("The %s doesn't take %s", PartName(container.type), itemDB.find(id) ? itemDB.find(id)->name.c_str() : id.c_str())); return; }
            const int n = (int)(std::min(centi, c.centi) / 100);
            if (n <= 0) return;
            c.centi -= (int64_t)n * 100;
            if (c.centi <= 0) c = InvSlot{};
            MachineAct_Insert(container.body, container.plot, container.part, id, n);
            Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
            return;
        }
        const bool toHotbar = i < Inventory::kGrid;
        const int from = toHotbar ? Inventory::kGrid : 0, to = toHotbar ? Inventory::kSlots : Inventory::kGrid;
        int64_t left = std::min(centi, c.centi);
        for (int pass = 0; pass < 2 && left > 0; ++pass) // onto the same item first, then an empty cell
            for (int k = from; k < to && left > 0; ++k) {
                const InvSlot& t = inventory.grid()[k];
                if (pass == 0 ? (t.empty() || t.id != c.id) : !t.empty()) continue;
                const int64_t n = put(k, c, std::min(left, room(k, c.id)));
                left -= n; c.centi -= n;
            }
        if (c.centi <= 0) c = InvSlot{};
    }
    // That item into cell i: from the open chest (it arrives on the next tick, merging with the item's cells), else from the
    // other part (grid <-> hotbar).
    void pullIn(int i, int64_t centi) override {
        const InvSlot c = inventory.grid()[i];
        if (c.empty()) return;
        if (GameUI_ContainerOpen()) {
            if (container.type == PartType::Chest) MachineAct_Take(container.body, container.plot, container.part, c.id, (int)std::max<int64_t>(1, centi / 100));
            return;
        }
        const bool fromHotbar = i < Inventory::kGrid;
        const int from = fromHotbar ? Inventory::kGrid : 0, to = fromHotbar ? Inventory::kSlots : Inventory::kGrid;
        for (int k = from; k < to; ++k) {
            InvSlot& s = inventory.grid()[k];
            if (s.empty() || s.id != c.id) continue;
            const int64_t n = std::min({ centi, s.centi, room(i, c.id) });
            if (n <= 0) return;
            InvSlot part = s; part.centi = n;
            put(i, part, n);
            s.centi -= n;
            if (s.centi <= 0) s = InvSlot{};
            return;
        }
    }
    void swapCell(int i, int cell) override { if (i != cell) inventory.swapSlots(i, cell); }
    void drop(int i, int64_t centi) override {
        InvSlot& c = inventory.grid()[i];
        if (c.empty()) return;
        InvSlot s = c; s.centi = std::min(centi, c.centi);
        if (!Combat_Drop(LocalSession(), s)) return;
        c.centi -= s.centi;
        if (c.centi <= 0) c = InvSlot{};
        Audio_Play(Sfx::Pickup, 0.4f, 0.6f);
    }
};
// One part of the open container's screen (which: 0 inputs / a chest's storage, 1 outputs, 2 fuel), as of the view taken
// this frame. Machines hold whole units and answer on the factory's next tick: a take comes onto the hand through the
// cursor sink (CursorSink), a put's leftovers back into the inventory.
struct ContainerGrid : SlotGrid {
    int which = 0;
    std::vector<std::pair<std::string, int>> slots;
    explicit ContainerGrid(int w) : which(w) {}
    int size() const override { return (int)slots.size(); }
    InvSlot at(int i) const override {
        InvSlot s;
        if (i >= 0 && i < (int)slots.size() && slots[i].second > 0) { s.id = slots[i].first; s.centi = (int64_t)slots[i].second * 100; }
        return s;
    }
    int64_t room(int i, const std::string& id) const override {
        if (which == 1 || id.empty()) return 0; // outputs only give
        const InvSlot s = at(i);
        if (which == 2) {
            const auto fuels = Factory::FuelsOf(container.type);
            if (std::find(fuels.begin(), fuels.end(), id) == fuels.end()) return 0;
        } else if (!ContainerTakes(id)) return 0;
        if (container.type == PartType::Chest || which == 2) { if (!s.empty() && s.id != id) return 0; return (int64_t)StackLimit(id) * 100 - s.centi; }
        return (int64_t)StackLimit(id) * 100; // a machine routes its inputs (Insert): any input slot click puts it in
    }
    int64_t put(int i, const InvSlot& held, int64_t centi) override {
        const int n = (int)(std::min(centi, held.centi) / 100);
        if (n <= 0 || room(i, held.id) <= 0) return 0;
        if (which == 2) MachineAct_Insert(container.body, container.plot, container.part, held.id, n, true);
        else if (container.type == PartType::Chest) MachineAct_InsertSlot(container.body, container.plot, container.part, held.id, n, i);
        else MachineAct_Insert(container.body, container.plot, container.part, held.id, n);
        Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
        return (int64_t)n * 100;
    }
    int64_t take(int i, int64_t centi, InvSlot& held) override {
        const InvSlot s = at(i);
        if (s.empty() || (!held.empty() && held.id != s.id)) return 0;
        const int n = (int)std::clamp<int64_t>(centi / 100, 1, s.centi / 100);
        MachineAct_TakeSlot(container.body, container.plot, container.part, i, which, n);
        Audio_Play(Sfx::Pickup, 0.4f, 1.2f);
        return 0; // (next tick, onto the hand)
    }
    void moveOut(int i, int64_t centi) override { // to the inventory
        const InvSlot s = at(i);
        if (s.empty()) return;
        if (centi < s.centi && container.type == PartType::Chest && which == 0) MachineAct_Take(container.body, container.plot, container.part, s.id, (int)std::max<int64_t>(1, centi / 100));
        else MachineAct_TakeSlot(container.body, container.plot, container.part, i, which); // (a machine's slot: all of it)
        Audio_Play(Sfx::Pickup, 0.4f, 1.2f);
    }
    void pullIn(int i, int64_t centi) override { // that item from the inventory into this slot
        const InvSlot s = at(i);
        if (s.empty() || which == 1) return;
        const int n = (int)std::min<double>(std::max<int64_t>(1, centi / 100), std::floor(inventory.count(s.id)));
        if (n <= 0 || !inventory.remove(s.id, n)) return;
        InvSlot h = s; h.centi = (int64_t)n * 100;
        if (put(i, h, h.centi) == 0) { ++inventory.muted; inventory.add(s.id, n); --inventory.muted; }
    }
    void swapCell(int i, int cell) override { // the slot's stack to the inventory, that cell's stack into the slot
        InvSlot& c = inventory.grid()[cell];
        const InvSlot s = at(i);
        if (!c.empty() && room(i, c.id) <= 0 && s.empty()) return;
        if (!s.empty()) MachineAct_TakeSlot(container.body, container.plot, container.part, i, which);
        if (!c.empty() && which != 1) {
            const InvSlot moving = c;
            c = InvSlot{};
            if (put(i, moving, moving.centi) < moving.centi) { ++inventory.muted; inventory.addAmount(moving.id, moving.amount(), moving.heat); --inventory.muted; }
        }
    }
};
// The creative Catalog: every item, free. A take gives a stack; a held stack put on it is deleted (as in Minecraft).
struct CatalogGrid : SlotGrid {
    std::vector<std::string> ids; // the items shown (the search applied)
    int size() const override { return (int)ids.size(); }
    InvSlot at(int i) const override { InvSlot s; s.id = ids[i]; s.centi = 100; return s; }
    int64_t room(int, const std::string&) const override { return 0; }
    int64_t put(int, const InvSlot& held, int64_t) override { return held.centi; }
    int64_t take(int i, int64_t centi, InvSlot& held) override {
        if (!held.empty() && held.id != ids[i]) return 0;
        const int64_t k = std::min(centi, capOf(ids[i]) - (held.empty() ? 0 : held.centi));
        if (k <= 0) return 0;
        if (held.empty()) { held = InvSlot{}; held.id = ids[i]; }
        held.centi += k;
        return k;
    }
    void moveOut(int i, int64_t) override { ++inventory.muted; inventory.add(ids[i], StackLimit(ids[i])); --inventory.muted; } // Shift-click: a stack into the inventory
    void pullIn(int, int64_t) override {}
    bool infinite() const override { return true; }
};
static PlayerGrid playerGrid;
static ContainerGrid containerIn(0), containerOut(1), containerFuel(2);
static CatalogGrid catalogGrid;
// The slot under the mouse this frame (set while the grids are drawn, read once after: GameUI_Draw).
static SlotGrid* hoverGrid = nullptr;
static int hoverSlot = -1;
static void NoteHover(SlotGrid* g, int i) {
    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem)) { hoverGrid = g; hoverSlot = i; }
    ImGui::SetItemKeyOwner(ImGuiKey_MouseWheelY); // the wheel moves items over a slot, not the page (#862)
}
// The drag preview (a split in progress): the share a slot would get, drawn on it.
static void DrawDragShare(ImDrawList* dl, const SlotGrid* g, int i, ImVec2 p, ImVec2 q) {
    for (const SlotDragPreview& d : SlotRules_DragPreview())
        if (d.grid == g && d.slot == i) {
            dl->AddRectFilled(p, q, IM_COL32(255, 220, 120, 50), 4.0f);
            const char* t = TextFormat("+%lld", (long long)(d.centi / 100));
            dl->AddText(ImVec2(p.x + 3, q.y - ImGui::GetTextLineHeight() - 2), IM_COL32(255, 220, 120, 255), t);
        }
}
// A whole click (press and release) through the rules: what the tests use, and what a click on a cell means.
static void RulesClick(SlotGrid* g, int i, bool right, bool shift) {
    SlotRules_SetHeld(&cursorStack);
    SlotInput p; p.shift = shift; (right ? p.rPress : p.lPress) = true; (right ? p.rDown : p.lDown) = true;
    SlotRules_Frame(p, g, i);
    SlotInput r; r.shift = shift; (right ? r.rRelease : r.lRelease) = true;
    SlotRules_Frame(r, g, i);
}
// One click on cell i of the player's inventory (#192, #862's rules).
static void CellClick(int i, bool right, bool shift) { SlotRules_SetGrids({ &playerGrid }); RulesClick(&playerGrid, i, right, shift); }

static void SlotClick(int i, int which, bool right, bool shift, int have);
bool GameUI_TestContainerCursor(Scene& world, uint64_t body, uint32_t plot, int part, int partType) {
    GameUI_OpenContainer(body, plot, part, partType);
    auto chest = [&]() { Factory::View(body, plot, part); Factory::WaitTicks(2); Factory::Sync(world); return Factory::View(body, plot, part); };
    auto settle = [&]() { Factory::WaitTicks(2); Factory::Sync(world); };
    auto inSlot = [&](const MachineView& v) { return v.inSlots.empty() ? 0 : v.inSlots[0].second; };
    ReturnCursorStack();
    inventory.clear();
    ++inventory.muted; inventory.add("stone", 20); --inventory.muted;
    bool ok = true;
    auto check = [&](bool c, const char* what) { if (!c) { printf("bench: container cursor: %s: NO\n", what); ok = false; } };
    chest(); // (the view follows this machine from now on)
    CellClick(0, false, false);                    check(cursorStack.centi == 2000, "picked the stack up from the grid");
    SlotClick(0, 0, false, false, 0); settle();    MachineView v = chest(); check(inSlot(v) == 20 && cursorStack.empty(), "click puts it into the chest");
    SlotClick(0, 0, true, false, inSlot(v)); settle(); v = chest(); check(cursorStack.id == "stone" && cursorStack.centi == 1000 && inSlot(v) == 10, "right-click takes half onto the cursor");
    SlotClick(0, 0, true, false, inSlot(v)); settle(); v = chest(); check(cursorStack.centi == 900 && inSlot(v) == 11, "right-click with a stack held puts one down");
    SlotClick(0, 0, false, true, inSlot(v)); settle(); v = chest(); check(inSlot(v) == 0 && inventory.count("stone") == 11, "Shift+click moves the slot to the inventory");
    CellClick(1, false, false);                    check(cursorStack.empty() && inventory.count("stone") == 20, "the held rest goes into the grid: nothing lost");
    // #566: a chest slot holds what an inventory cell does (500 stone in one slot, not stacks of 64).
    inventory.clear();
    ++inventory.muted; inventory.add("stone", 500); --inventory.muted; // one cell of 500
    { int c = -1; for (int k = 0; k < Inventory::kSlots; ++k) if (inventory.grid()[k].id == "stone") { c = k; break; }
      if (c >= 0) CellClick(c, false, false); }
    SlotClick(0, 0, false, false, 0); settle(); v = chest();
    check(inSlot(v) == 500 && cursorStack.empty(), "a chest slot takes 500 stone (#566)");
    // #568: tools and machines go in too.
    ++inventory.muted; inventory.add("stone_pickaxe", 1); inventory.add("chest", 2); --inventory.muted;
    for (const char* id : { "stone_pickaxe", "chest" }) {
        int c = -1; for (int k = 0; k < Inventory::kSlots; ++k) if (inventory.grid()[k].id == id) { c = k; break; }
        if (c >= 0) CellClick(c, false, false);
        SlotClick(id[0] == 's' ? 1 : 2, 0, false, false, 0); settle();
        check(cursorStack.empty() && inventory.count(id) == 0, TextFormat("the chest takes %s (#568)", id));
    }
    v = chest();
    check(v.inSlots.size() > 2 && v.inSlots[1].first == "stone_pickaxe" && v.inSlots[2].first == "chest" && v.inSlots[2].second == 2, "...and holds them");
    // #567: the look-at panel lists totals per item, wrapped and capped.
    Factory::Insert(body, plot, part, "stone", 300); settle(); v = chest(); // 800 stone: still one slot
    check(v.input.find("stone x800") != std::string::npos && v.input.find("stone x", v.input.find("stone x800") + 1) == std::string::npos, TextFormat("totals per item (#567): '%s'", v.input.c_str()));
    std::string many;
    for (int k = 0; k < 40; ++k) many += TextFormat("%sitem_number_%d x%d", k ? ", " : "", k, 64 * k);
    const std::vector<std::string> lines = WrapItemList("holds: ", many, 404, 16, 4);
    bool fits = lines.size() == 4 && lines.back().find("more") != std::string::npos;
    for (const std::string& l : lines) fits = fits && MeasureText(l.c_str(), 16) <= 404;
    check(fits, TextFormat("a 40-item list wraps into 4 lines inside the panel, the last saying what's left out (%d lines)", (int)lines.size()));
    Factory::Extract(body, plot, part); settle();
    GameUI_CloseContainer();
    printf("bench: container cursor (#192): %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

static void CursorSink(const std::string& id, int count, float heat);
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
    // #564: the bin deletes the held stack and gives the last one back; a drop leaves a pile at the feet that is picked
    // up only after walking away and coming back.
    auto stoneCell = [&]() { for (int k = 0; k < Inventory::kGrid; ++k) if (g[k].id == "stone") return k; return -1; };
    CellClick(stoneCell(), false, false); BinClick(); check(cursorStack.empty() && inventory.count("stone") == 0 && binned.id == "stone", "the bin deletes the held stack");
    BinClick(); check(cursorStack.id == "stone" && cursorStack.centi == 1000 && binned.empty(), "...and gives it back on a second click");
    PlayerSession& s = LocalSession();
    const int piles0 = Combat_DroppedPiles();
    DropCursor(100); DropCursor(cursorStack.centi);
    check(cursorStack.empty() && inventory.count("stone") == 0 && Combat_DroppedPiles() == piles0 + 1, "a drop (one, then the rest) makes one pile");
    if (s.player.nearest) {
        const DVec3 here = s.player.position;
        const PlayerInput none{};
        Combat_Step(s, 1.0f / 60.0f, none); check(inventory.count("stone") == 0, "standing on it, it stays on the ground");
        const DVec3 local = s.player.nearest->toLocalD(here);
        const Vector3 side = Vector3Normalize(Vector3Perpendicular(Vector3Normalize(local.toFloat())));
        s.player.position = s.player.nearest->toWorldD(local + Vector3Scale(side, 6.0f)); // 6 m along the ground
        Combat_Step(s, 1.0f / 60.0f, none); check(inventory.count("stone") == 0, "walking away leaves it");
        s.player.position = here;
        Combat_Step(s, 1.0f / 60.0f, none); check(inventory.count("stone") == 10 && Combat_DroppedPiles() == piles0, "coming back picks it all up");
    } else check(false, "the player is on a body");
    // #580: a stack taken hot out of a furnace stays hot on the cursor, and when put down; a cold one joining it averages.
    CursorSink("iron_ingot", 4, 1.0f);
    check(cursorStack.id == "iron_ingot" && cursorStack.heat > 0.99f, "a hot stack onto the cursor stays hot (#580)");
    CursorSink("iron_ingot", 4, 0.0f);
    check(fabsf(cursorStack.heat - 0.5f) < 0.01f, "a cold one joining it halves the heat");
    int empty = -1; for (int k = 0; k < Inventory::kGrid && empty < 0; ++k) if (g[k].empty()) empty = k;
    if (empty >= 0) CellClick(empty, false, false);
    check(empty >= 0 && g[empty].id == "iron_ingot" && fabsf(g[empty].heat - 0.5f) < 0.01f && cursorStack.empty(), "put down, the cell keeps the heat");
    printf("bench: inventory cursor (#192): %s\n", ok ? "PASS" : "FAIL");
    return ok;
}

// The grid (#192): click a cell to pick its stack up onto the cursor, click a cell to put it down (merging the same item,
// swapping another); right-click takes half, or with a stack held puts one down; Shift+click moves a stack between the grid
// and the hotbar. Hover and press 1-0 to put the item on the hotbar. With a container open a click sends the stack into it.
// One inventory cell (grid, hotbar or the off hand): click to pick up or put down (#192), right-click for one or half,
// shift for the quick move, middle-click to sort; its icon, name, amount and tooltip. cw: its width.
static void DrawCell(int i, float cw, float cell) {
    auto& grid = inventory.grid();
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const bool hotbar = i >= Inventory::kGrid && i < Inventory::kSlots;
    ImGui::PushID(i);
    const ImVec2 p = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("cell", ImVec2(cw, cell), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
    if (ImGui::IsItemClicked()) OpenIfAlt(grid[i].id); // (#826; the rules skip alt-clicks)
    NoteHover(&playerGrid, i); // clicks, drags, keys and the wheel: the item rules (#862), once a frame in GameUI_Draw
    if (ImGui::IsItemClicked(ImGuiMouseButton_Middle) && !terrainTool.creative) { ReturnCursorStack(); inventory.sort(); } // as the Sort button (#577; creative: a full stack, the rules)
    const bool hover = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
    const InvSlot& s = grid[i];
    const ImVec2 q(p.x + cw, p.y + cell);
    dl->AddRectFilled(p, q, IM_COL32(30, 32, 36, 220), 4.0f);
    DrawGlow(dl, p, q, Item_Glow(grid[i])); // hot or radiant (#580)
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
        if (hover && cursorStack.empty())
            ImGui::SetTooltip("%s  (%s)\n%.2f of %d per slot\n%s\n(Shift+click: %s; 1-0: swap with that hotbar slot; wheel: move one)", d ? d->name.c_str() : s.id.c_str(),
                              KindName(KindOf(s.id)), s.amount(), StackLimit(s.id), d ? d->description.c_str() : "",
                              GameUI_ContainerOpen() ? TextFormat("into the %s", PartName(container.type)) : i < Inventory::kGrid ? "to the hotbar" : "to the grid");
    }
    DrawDragShare(dl, &playerGrid, i, p, q);
    if (hotbar) dl->AddText(ImVec2(p.x + 3, p.y + 1), IM_COL32(150, 150, 150, 255), TextFormat("%d", (i - Inventory::kGrid + 1) % 10));
    const ImU32 edge = hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 72, 78, 255);
    dl->AddRect(p, q, edge, 4.0f, 0, 1.0f);
    ImGui::PopID();
}

static void DrawGrid() {
    const float cell = 54.0f, gap = 4.0f;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (int i = 0; i < Inventory::kSlots; ++i) {
        // The 10 x 5 grid, then the hotbar's 10 cells in a row under it (#209: real slots; drop stacks into them). The
        // off-hand cell is beside the player now (#822).
        const bool hotbar = i >= Inventory::kGrid;
        if (i == Inventory::kGrid) { ImGui::Dummy(ImVec2(0.0f, 2.0f)); ImGui::TextDisabled("Hotbar (keys 1-0)"); }
        else if (i < Inventory::kGrid ? i % Inventory::kCols != 0 : true) ImGui::SameLine(0.0f, gap);
        DrawCell(i, hotbar ? (Inventory::kCols * cell - gap) / (float)Inventory::kHotbar : cell, cell); // the hotbar row as wide as the grid
    }
    ImGui::Dummy(ImVec2(0.0f, 2.0f));
    { // the bin, under the hotbar (#564)
        const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + cell, p.y + cell);
        if (ImGui::InvisibleButton("bin", ImVec2(cell, cell))) BinClick();
        const bool hover = ImGui::IsItemHovered();
        dl->AddRectFilled(p, q, IM_COL32(48, 26, 26, 220), 4.0f);
        if (!binned.empty()) DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y + cell * 0.1f), cell * 0.64f, binned.id, true);
        const ImVec2 sz = ImGui::CalcTextSize("Bin");
        dl->AddText(ImVec2(p.x + (cell - sz.x) * 0.5f, q.y - sz.y - 2.0f), IM_COL32(220, 160, 150, 255), "Bin");
        dl->AddRect(p, q, hover ? IM_COL32(230, 150, 140, 255) : IM_COL32(110, 60, 60, 255), 4.0f, 0, 1.0f);
        if (hover) ImGui::SetTooltip(!cursorStack.empty() ? "Click: delete the held stack" : !binned.empty() ? "Click: take the deleted stack back" :
                                     "The bin: put a stack here to delete it\n(or click outside the window to drop it on the ground)");
    }
    if (ImGui::Button("Sort (middle-click)")) { ReturnCursorStack(); inventory.sort(); }
    ImGui::SameLine();
    // Deposit (button or Q): items a nearby chest already holds go into it.
    const int nChests = uiWorld ? (int)Storage_NearbyChests(*uiWorld).size() : 0;
    ImGui::BeginDisabled(nChests == 0);
    if ((ImGui::Button("Deposit (Q)") || (nChests > 0 && PlayerInput_BindingPressed(Act::Deposit))) && uiWorld) {
        const int kinds = Storage_Deposit(*uiWorld);
        GameUI_Toast(kinds > 0 ? TextFormat("Deposited %d kind%s of item into nearby chests", kinds, kinds > 1 ? "s" : "") : "Nothing here matches what the nearby chests hold");
        if (kinds > 0) Audio_Play(Sfx::Pickup, 0.5f, 0.8f);
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    ImGui::TextDisabled(!cursorStack.empty() ? "Click: put it down  right-click: one  outside the window: drop it" : nChests > 0 ? TextFormat("%d chest%s within reach", nChests, nChests > 1 ? "s" : "") : "Tools 1, blocks 256, resources 8192");
    if (!cursorStack.empty()) { // the held stack follows the mouse, over everything (#192)
        ImDrawList* fg = ImGui::GetForegroundDrawList();
        const ImVec2 m = ImGui::GetMousePos();
        DrawGlow(fg, ImVec2(m.x - 20.0f, m.y - 20.0f), ImVec2(m.x + 20.0f, m.y + 20.0f), Item_Glow(cursorStack)); // (#580)
        DrawItemIcon(fg, ImVec2(m.x - 18.0f, m.y - 18.0f), 36.0f, cursorStack.id);
        if (KindOf(cursorStack.id) != ItemKind::Tool) fg->AddText(ImVec2(m.x + 12.0f, m.y + 8.0f), IM_COL32(255, 255, 255, 255), AmountText(cursorStack).c_str());
    }
}

// The Crafting screen (#824): search and category chips over the recipes you know; a click selects one, and its page
// shows what it makes, what it takes (with what you have, nearby chests counted), Craft 1 / 5 / All, and what it's used
// by. Crafts go in a queue, one every kCraftEvery, each through Crafting_Request (here, or the server when joined). A
// single craft short of something pulls it from nearby chests first (it arrives next tick). In creative everything is
// known and free.
namespace {
    enum class CraftCat { All, CanMake, Tools, Building, Materials, Food, Suit, Count };
    const char* kCraftCatNames[] = { "All", "Can make", "Tools", "Building", "Materials", "Food", "Suit" };
    CraftCat CategoryOf(const std::string& id) {
        if (Suit_IsSuit(id) || Suit_IsModule(id)) return CraftCat::Suit;
        if (Survival_IsFood(id)) return CraftCat::Food;
        switch (KindOf(id)) {
            case ItemKind::Tool: return CraftCat::Tools;
            case ItemKind::Block: return CraftCat::Building;
            default: return CraftCat::Materials;
        }
    }
    struct QueuedCraft { int recipe; int left; };
    std::vector<QueuedCraft> craftQueue;
    double nextCraftAt = 0.0;
    constexpr double kCraftEvery = 0.3; // s between queued crafts (you see them come)
    int selectedRecipe = -1;
}

// The queue runs whether the window is open or not (GameUI_Draw calls it every frame).
static void RunCraftQueue() {
    if (craftQueue.empty() || GetTime() < nextCraftAt) return;
    nextCraftAt = GetTime() + kCraftEvery;
    QueuedCraft& q = craftQueue.front();
    const auto& recipes = itemDB.recipes();
    bool done = q.recipe < 0 || q.recipe >= (int)recipes.size();
    if (!done) {
        const Recipe& r = recipes[q.recipe];
        if (terrainTool.creative) { inventory.addAmount(r.output, r.outputCount); Audio_Play(Sfx::Craft); --q.left; }
        else if (gNetRole == NetRole::Client || inventory.canCraft(r)) { Crafting_Request(q.recipe, 1); --q.left; } // (#338: the server checks)
        else { const ItemDef* d = itemDB.find(r.output); GameUI_Toast(TextFormat("Ran out of what %s takes", d ? d->name.c_str() : r.output.c_str())); done = true; }
    }
    if (done || q.left <= 0) craftQueue.erase(craftQueue.begin());
}

static void DrawCrafting() {
    const bool creative = terrainTool.creative;
    if (creative) ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f), "Creative: crafting costs nothing");
    ToolWear_DrawRepairUI();
    static char query[64] = {};
    static CraftCat cat = CraftCat::All;
    static int pending = -1;           // a craft waiting for items pulled from chests
    static double pendingUntil = 0.0;
    const auto& recipes = itemDB.recipes();
    ImGui::SetNextItemWidth(260.0f);
    ImGui::InputTextWithHint("##recipesearch", "search recipes", query, sizeof(query));
    for (int c = 0; c < (int)CraftCat::Count; ++c) { // the category chips
        ImGui::SameLine();
        const bool on = (int)cat == c;
        if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0.88f, 0.66f, 0.25f, 1.0f)), ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.1f, 0.08f, 0.02f, 1.0f));
        if (ImGui::SmallButton(kCraftCatNames[c])) cat = (CraftCat)c;
        if (on) ImGui::PopStyleColor(2);
    }
    std::string q = query;
    for (char& c : q) c = (char)tolower((unsigned char)c);
    // What nearby chests hold, once for this frame.
    std::map<std::string, int> stored;
    if (uiWorld)
        for (const NearChest& c : Storage_NearbyChests(*uiWorld))
            for (const auto& s : Factory::ChestContents(c.body, c.plot, c.part)) stored[s.first] += s.second;
    auto have = [&](const std::string& id) { int n = inventory.count(id); for (const auto& s : stored) if (ItemMatches(id, s.first)) n += s.second; return n; };
    auto canMake = [&](const Recipe& r) { if (creative) return true; for (const auto& in : r.inputs) if (have(in.first) < in.second) return false; return true; };
    auto maxFromInventory = [&](const Recipe& r) {
        if (creative) return 99;
        int n = 1 << 20;
        for (const auto& in : r.inputs) n = std::min(n, inventory.count(in.first) / std::max(1, in.second));
        return n;
    };
    // A pending craft completes once the pulled items have arrived.
    if (pending >= 0 && pending < (int)recipes.size()) {
        const Recipe& r = recipes[pending];
        if (inventory.canCraft(r)) { Crafting_Request(pending, 1); pending = -1; } // (#338) here, or the server when joined
        else if (GetTime() > pendingUntil) { GameUI_Toast("The chests no longer hold enough"); pending = -1; }
    }
    if (selectedRecipe >= (int)recipes.size() || (selectedRecipe >= 0 && !RecipeBook_Known(recipes[selectedRecipe]))) selectedRecipe = -1;

    // Left: the recipes as icons (greyed: can't make now). Right: the selected recipe's page and the queue.
    const float pageW = std::min(420.0f, ImGui::GetContentRegionAvail().x * 0.48f);
    int hidden = 0, shown = 0;
    ImGui::BeginChild("recipes", ImVec2(ImGui::GetContentRegionAvail().x - pageW - 8.0f, 0.0f));
    const float cell = 60.0f, gap = 6.0f;
    const int cols = std::max(1, (int)((ImGui::GetContentRegionAvail().x + gap) / (cell + gap)));
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (size_t i = 0; i < recipes.size(); ++i) {
        const Recipe& r = recipes[i];
        if (!RecipeBook_Known(r)) { ++hidden; continue; }
        if (!RecipeBook_Matches(r, q)) continue;
        const bool can = canMake(r);
        if (cat == CraftCat::CanMake && !can) continue;
        if (cat != CraftCat::All && cat != CraftCat::CanMake && CategoryOf(r.output) != cat) continue;
        if (selectedRecipe < 0) selectedRecipe = (int)i;
        if (shown % cols != 0) ImGui::SameLine(0.0f, gap);
        ++shown;
        const ItemDef* out = itemDB.find(r.output);
        ImGui::PushID((int)i);
        const ImVec2 p = ImGui::GetCursorScreenPos(), qd(p.x + cell, p.y + cell);
        if (ImGui::InvisibleButton("recipe", ImVec2(cell, cell)) && !OpenIfAlt(r.output)) selectedRecipe = (int)i;
        const bool hover = ImGui::IsItemHovered(), sel = selectedRecipe == (int)i;
        if (hover && ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left) && can) craftQueue.push_back({ (int)i, 1 }); // double-click: craft one
        dl->AddRectFilled(p, qd, can ? IM_COL32(34, 40, 36, 230) : IM_COL32(30, 31, 34, 200), 4.0f);
        DrawItemIcon(dl, ImVec2(p.x + 4, p.y + 2), cell - 8, r.output, !can);
        if (r.outputCount > 1) {
            const char* n = TextFormat("x%d", r.outputCount);
            dl->AddText(ImVec2(qd.x - ImGui::CalcTextSize(n).x - 4, qd.y - 17), IM_COL32(255, 255, 255, can ? 255 : 130), n);
        }
        dl->AddRect(p, qd, sel ? IM_COL32(224, 169, 64, 255) : hover ? IM_COL32(220, 220, 220, 255) : can ? IM_COL32(90, 150, 95, 255) : IM_COL32(64, 66, 72, 255),
                    4.0f, 0, sel || hover ? 2.0f : 1.0f);
        if (hover) ImGui::SetTooltip("%s%s", out ? out->name.c_str() : r.output.c_str(), r.outputCount > 1 ? TextFormat("  x%d", r.outputCount) : "");
        ImGui::PopID();
    }
    if (shown == 0) ImGui::TextDisabled(cat == CraftCat::CanMake ? "Nothing you can make right now" : "No known recipe matches");
    if (hidden > 0) ImGui::TextDisabled("%d recipes still to discover: find new materials (and some need more)", hidden);
    ImGui::EndChild();
    ImGui::SameLine();

    ImGui::BeginChild("recipepage", ImVec2(pageW, 0.0f), true);
    if (selectedRecipe >= 0) {
        const Recipe& r = recipes[selectedRecipe];
        const ItemDef* out = itemDB.find(r.output);
        const ImVec2 p = ImGui::GetCursorScreenPos();
        ImGui::GetWindowDrawList()->AddRectFilled(p, ImVec2(p.x + 56, p.y + 56), IM_COL32(34, 40, 36, 230), 4.0f);
        DrawItemIcon(ImGui::GetWindowDrawList(), ImVec2(p.x + 4, p.y + 2), 48, r.output);
        ImGui::Dummy(ImVec2(56, 56));
        ImGui::SameLine();
        ImGui::BeginGroup();
        ImGui::TextColored(ImVec4(0.88f, 0.66f, 0.25f, 1.0f), "%s%s", out ? out->name.c_str() : r.output.c_str(), r.outputCount > 1 ? TextFormat("  x%d", r.outputCount) : "");
        ImGui::TextDisabled("%s (%s)", KindName(KindOf(r.output)), kCraftCatNames[(int)CategoryOf(r.output)]);
        ImGui::EndGroup();
        if (out && !out->description.empty()) ImGui::TextWrapped("%s", out->description.c_str());
        ImGui::Separator();
        ImGui::TextUnformatted("Takes:");
        for (const auto& in : r.inputs) {
            const ItemDef* id = itemDB.find(in.first);
            const int inv = inventory.count(in.first), all = have(in.first);
            const bool enough = creative || all >= in.second;
            ImGui::TextColored(enough ? ImVec4(0.55f, 0.9f, 0.55f, 1.0f) : ImVec4(0.95f, 0.45f, 0.4f, 1.0f), "  %d %s   (have %d%s)", in.second,
                               id ? id->name.c_str() : in.first.c_str(), all, all > inv ? TextFormat(", %d in chests", all - inv) : "");
            if (ImGui::IsItemClicked()) { HelpManual_OpenItem(in.first); databaseOpen = true; } // (#826) an ingredient: its page
        }
        const bool can = canMake(r);
        const int most = maxFromInventory(r);
        ImGui::Spacing();
        ImGui::BeginDisabled(!can || pending >= 0);
        if (ImGui::Button("Craft 1")) {
            if (creative || gNetRole == NetRole::Client || inventory.canCraft(r)) craftQueue.push_back({ selectedRecipe, 1 });
            else if (uiWorld) { // short in the inventory: pulled from nearby chests first (#338: not when joined)
                bool ok = true;
                for (const auto& in : r.inputs) ok = Storage_Pull(*uiWorld, in.first, in.second - inventory.count(in.first)) && ok;
                if (ok) { pending = selectedRecipe; pendingUntil = GetTime() + 2.0; }
            }
        }
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(most < 5);
        if (ImGui::Button("x5")) craftQueue.push_back({ selectedRecipe, 5 });
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(most < 1);
        if (ImGui::Button(creative ? "x20" : TextFormat("All (%d)", most))) craftQueue.push_back({ selectedRecipe, creative ? 20 : most });
        ImGui::EndDisabled();
        if (!can) ImGui::TextDisabled("Missing items");
        else if (pending >= 0) ImGui::TextDisabled("Waiting for items from the chests...");
        // What this is used by: the recipes that take it.
        ImGui::Separator();
        ImGui::TextUnformatted("Used by:");
        int uses = 0;
        for (size_t k = 0; k < recipes.size(); ++k) {
            const Recipe& u = recipes[k];
            if (!RecipeBook_Known(u)) continue;
            bool takes = false;
            for (const auto& in : u.inputs) takes = takes || ItemMatches(in.first, r.output);
            if (!takes) continue;
            const ItemDef* ud = itemDB.find(u.output);
            ImGui::PushID(10000 + (int)k);
            if (ImGui::SmallButton(ud ? ud->name.c_str() : u.output.c_str())) selectedRecipe = (int)k;
            ImGui::PopID();
            if (++uses % 3 != 0) ImGui::SameLine();
        }
        if (uses == 0) ImGui::TextDisabled("nothing you know yet");
        else if (uses % 3 != 0) ImGui::NewLine();
    } else ImGui::TextDisabled("Pick a recipe");
    // The queue: what's being crafted, next first.
    ImGui::Separator();
    ImGui::TextUnformatted(craftQueue.empty() ? "Queue: empty" : "Queue:");
    for (size_t k = 0; k < craftQueue.size(); ++k) {
        const QueuedCraft& c = craftQueue[k];
        const ItemDef* d = c.recipe >= 0 && c.recipe < (int)recipes.size() ? itemDB.find(recipes[c.recipe].output) : nullptr;
        ImGui::PushID(20000 + (int)k);
        if (ImGui::SmallButton("x")) { craftQueue.erase(craftQueue.begin() + (long)k); ImGui::PopID(); break; }
        ImGui::SameLine();
        ImGui::Text("%s  x%d%s", d ? d->name.c_str() : "?", c.left, k == 0 ? "  (crafting)" : "");
        ImGui::PopID();
    }
    ImGui::EndChild();
}

// A machine's or chest's screen: status, progress and charge, then its slots (input and output for
// machines, 16 for a chest). Click a slot to take its stack; click an inventory cell to put one in.
// One click on slot i of the open container (#192: the same slot behaviour as the player's inventory). which: 0 inputs,
// 1 outputs, 2 fuel; have: the count in it (as of the last tick). Machine slots change on the factory's next tick, so a
// pick-up arrives on the cursor a tick later (Factory::TakeSlotToCursor -> CursorSink).
static ContainerGrid& ContainerPart(int which) { return which == 1 ? containerOut : which == 2 ? containerFuel : containerIn; }
// (the tests: as a click on the screen, with the slots as the view has them now)
static void SlotClick(int i, int which, bool right, bool shift, int) {
    const MachineView v = Factory::View(container.body, container.plot, container.part);
    containerIn.slots = v.inSlots; containerOut.slots = v.outSlots; containerFuel.slots = { v.fuelSlot };
    SlotRules_SetGrids({ &playerGrid, &containerIn, &containerOut, &containerFuel });
    RulesClick(&ContainerPart(which), i, right, shift);
}

static void DrawSlots(const std::vector<std::pair<std::string, int>>& slots, int which, int cols, bool hot = false) { // which: 0 inputs, 1 outputs, 2 fuel; hot: its items glow (#580)
    const bool output = which == 1;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float cell = 54.0f, gap = 4.0f;
    ContainerGrid& grid = ContainerPart(which);
    grid.slots = slots; // (this frame's view: what the rules see)
    for (int i = 0; i < (int)slots.size(); ++i) {
        if (i % cols != 0) ImGui::SameLine(0.0f, gap);
        ImGui::PushID(i + which * 1000);
        const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + cell, p.y + cell);
        ImGui::InvisibleButton("slot", ImVec2(cell, cell), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
        if (ImGui::IsItemClicked()) OpenIfAlt(slots[i].second > 0 ? slots[i].first : std::string()); // (#826; the rules skip alt-clicks)
        NoteHover(&grid, i); // clicks, drags, keys and the wheel: the item rules (#862), once a frame in GameUI_Draw
        const bool hover = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
        dl->AddRectFilled(p, q, which == 2 ? IM_COL32(52, 30, 20, 220) : output ? IM_COL32(40, 36, 28, 220) : IM_COL32(30, 32, 36, 220), 4.0f);
        if (slots[i].second > 0) DrawGlow(dl, p, q, Item_GlowOf(slots[i].first, hot ? 1.0f : 0.0f)); // fresh out of a furnace, or radiant (#580)
        const auto& s = slots[i];
        if (s.second > 0) {
            DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y), cell * 0.64f, s.first);
            const ItemDef* d = itemDB.find(s.first);
            CellName(dl, p, q, d ? d->name : s.first, hover);
            const char* amt = TextFormat("%d", s.second);
            const ImVec2 sz = ImGui::CalcTextSize(amt);
            dl->AddRectFilled(ImVec2(q.x - sz.x - 6, p.y + 2), ImVec2(q.x - 2, p.y + sz.y + 2), IM_COL32(0, 0, 0, 160), 2.0f);
            dl->AddText(ImVec2(q.x - sz.x - 4, p.y + 2), IM_COL32(255, 255, 255, 255), amt);
            if (hover && cursorStack.empty()) ImGui::SetTooltip("%s x%d\n%s\n(click: pick up  right-click: half  Shift+click: to the inventory  wheel: move one)", d ? d->name.c_str() : s.first.c_str(), s.second, d ? d->description.c_str() : "");
        }
        DrawDragShare(dl, &grid, i, p, q);
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
    DrawSlots(v.inSlots, 0, chest ? Inventory::kCols : 4); // a chest lines up with the player's grid (#495)
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
        DrawSlots(v.outSlots, 1, 4, v.hotOutputs);
    }
    ImGui::Spacing();
    if (ImGui::Button("Put in all (I)") || PlayerInput_BindingPressed(Act::Insert)) {
        // (#332) joined: one request, the server puts in all it takes (one per item would run into the rate limit)
        if (MachineAct_Remote()) { gMachineRequest(container.body, container.plot, container.part, 2, std::string(), 0, 0, 0); Audio_Play(Sfx::Pickup, 0.4f, 0.8f); }
        const auto all = MachineAct_Remote() ? decltype(inventory.all()){} : inventory.all(); // copy: we remove while walking it
        int kinds = 0;
        for (const auto& kv : all) {
            if (kv.second <= 0 || !ContainerTakes(kv.first) || KindOf(kv.first) == ItemKind::Tool) continue; // tools stay with the player (one in the off hand too)
            bool bound = false;
            for (int s = 0; s < 10; ++s) bound = bound || Hotbar_Slot(s) == kv.first;
            if (chest && bound) continue; // keep what's on the hotbar
            if (inventory.remove(kv.first, kv.second)) { MachineAct_Insert(container.body, container.plot, container.part, kv.first, kv.second); ++kinds; }
        }
        if (kinds > 0) Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
    }
    ImGui::SameLine();
    if (ImGui::Button("Take all (O)") || PlayerInput_BindingPressed(Act::Extract)) {
        if (chest) for (int i = 0; i < (int)v.inSlots.size(); ++i) { if (v.inSlots[i].second > 0) MachineAct_TakeSlot(container.body, container.plot, container.part, i, 0); }
        else MachineAct_Extract(container.body, container.plot, container.part);
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
    ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f), "Creative: build, craft and place blueprints for free (/gamemode survival to leave)"); // (#823)
    static char filter[64] = {};
    ImGui::SetNextItemWidth(220.0f);
    ImGui::InputTextWithHint("##filter", "search items", filter, sizeof(filter));
    ImGui::SameLine();
    if (ImGui::Button("Clear inventory")) inventory.clear();
    {   // Creatures at the crosshair (#529): any species, on the ground the view's centre is on (within 200 m).
        static int kind = 0, count = 1;
        static const auto list = Animals_SpeciesList();
        if (!list.empty()) {
            kind = std::clamp(kind, 0, (int)list.size() - 1);
            ImGui::SetNextItemWidth(180.0f);
            if (ImGui::BeginCombo("##creature", list[kind].second.c_str())) {
                for (int i = 0; i < (int)list.size(); ++i)
                    if (ImGui::Selectable(list[i].second.c_str(), i == kind)) kind = i;
                ImGui::EndCombo();
            }
            ImGui::SameLine();
            ImGui::SetNextItemWidth(90.0f);
            ImGui::InputInt("##count", &count);
            count = std::clamp(count, 1, 20);
            ImGui::SameLine();
            if (ImGui::Button("Spawn at crosshair") && uiWorld) {
                Scene& w = *uiWorld;
                const Vector3 dir = Vector3Normalize(Vector3Subtract(w.camera.target, w.camera.position));
                TerrainHit th;
                if (RaycastTerrain(w, Ray{ w.camera.position, dir }, 200.0f, th) && th.planetoid) {
                    const int n = Animals_Spawn(w, *th.planetoid, list[kind].first, th.local, count, false, count > 1 ? 3.0f : 0.0f);
                    GameUI_Toast(TextFormat("Spawned %d %s", n, list[kind].second.c_str()));
                } else GameUI_Toast("Nothing under the crosshair within 200 m");
            }
        }
    }
    ImGui::Separator();
    // Every item as an icon grid (#862, as Minecraft's creative menu): click takes a full stack, right-click one, Shift+click a
    // stack into the inventory; a held stack put on the grid is deleted. Kinds as tabs over it; the search narrows them.
    static int kindTab = -1; // -1 all, else an ItemKind
    {
        const char* kinds[] = { "All", "Resources", "Blocks", "Tools" };
        for (int k = -1; k < 3; ++k) {
            if (k > -1) ImGui::SameLine();
            const bool on = kindTab == k;
            if (on) ImGui::PushStyleColor(ImGuiCol_Button, ImGui::GetStyleColorVec4(ImGuiCol_ButtonActive));
            if (ImGui::Button(kinds[k + 1])) kindTab = k;
            if (on) ImGui::PopStyleColor();
        }
    }
    std::string f = filter;
    for (char& c : f) c = (char)tolower(c);
    catalogGrid.ids.clear();
    for (const ItemDef& d : itemDB.items()) {
        if (IsWildcardItem(d.id)) continue; // tags, not items (#269)
        if (kindTab >= 0 && (int)d.kind != kindTab) continue;
        std::string name = d.name;
        for (char& c : name) c = (char)tolower(c);
        if (!f.empty() && name.find(f) == std::string::npos && d.id.find(f) == std::string::npos) continue;
        catalogGrid.ids.push_back(d.id);
    }
    ImGui::BeginChild("creativeitems");
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float cell = 54.0f, gap = 4.0f;
    const int cols = std::max(1, (int)((ImGui::GetContentRegionAvail().x + gap) / (cell + gap)));
    for (int i = 0; i < (int)catalogGrid.ids.size(); ++i) {
        const ItemDef* d = itemDB.find(catalogGrid.ids[i]);
        if (i % cols != 0) ImGui::SameLine(0.0f, gap);
        ImGui::PushID(i);
        const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + cell, p.y + cell);
        ImGui::InvisibleButton("item", ImVec2(cell, cell), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
        if (ImGui::IsItemClicked()) OpenIfAlt(catalogGrid.ids[i]); // (#826)
        NoteHover(&catalogGrid, i);
        const bool hover = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
        dl->AddRectFilled(p, q, IM_COL32(30, 32, 36, 220), 4.0f);
        DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y), cell * 0.64f, catalogGrid.ids[i]);
        CellName(dl, p, q, d ? d->name : catalogGrid.ids[i], hover);
        if (hover && cursorStack.empty() && d)
            ImGui::SetTooltip("%s (%s, %d per slot)\n%s\nhave %.2f\n(click: a stack  right-click: one  Shift+click: a stack into the inventory)",
                              d->name.c_str(), KindName(d->kind), StackLimit(d->id), d->description.c_str(), inventory.amount(d->id));
        else if (hover && !cursorStack.empty()) ImGui::SetTooltip("Click: delete the held stack");
        dl->AddRect(p, q, hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 72, 78, 255), 4.0f, 0, 1.0f);
        ImGui::PopID();
    }
    if (catalogGrid.ids.empty()) ImGui::TextDisabled("No items match");
    ImGui::EndChild();
}

// The team's progress (#825): quests and tiers are a team's (#346); a solo player is a team of one.
static void DrawTeamProgress() {
    const int cur = QuestBook_Tier(), tiersN = QuestBook_TierCount();
    std::string team = "Your team";
    const int id = Teams_Of(PlayerId_Local());
    for (const TeamRecord& t : Teams_All()) if (t.id == id && !t.name.empty()) team = t.name;
    int done = 0, all = cur < tiersN ? QuestBook_TaskCount(cur) : 0;
    for (int k = 0; k < all; ++k) done += QuestBook_TaskDone(cur, k);
    ImGui::TextColored(ImVec4(0.88f, 0.66f, 0.25f, 1.0f), "%s", team.c_str());
    ImGui::SameLine();
    if (cur < tiersN) ImGui::TextDisabled("tier %d of %d (%s): %d of %d tasks done   pinned to the HUD: %d of %d", cur, tiersN - 1, QuestBook_TierName(cur), done, all,
                                          (int)QuestBook_Pins().size(), kMaxPins);
    else ImGui::TextDisabled("every tier done");
}

// The main quest line: done tiers collapsed, the current one open, the next shown. Each open task has a pin: up to
// kMaxPins show on the HUD (#825).
static void DrawQuestList() {
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
                ImGui::PushID(k);
                if (!kd) {
                    bool pinned = QuestBook_Pinned(i, k);
                    const bool full = !pinned && (int)QuestBook_Pins().size() >= kMaxPins;
                    ImGui::BeginDisabled(full);
                    if (ImGui::Checkbox("##pin", &pinned)) QuestBook_Pin(i, k, pinned);
                    ImGui::EndDisabled();
                    if (ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) ImGui::SetTooltip(full ? "Three are pinned: unpin one first" : pinned ? "Pinned to the HUD" : "Pin to the HUD");
                    ImGui::SameLine();
                } else ImGui::Dummy(ImVec2(ImGui::GetFrameHeight(), 0.0f)), ImGui::SameLine();
                ImGui::TextColored(kd ? ImVec4(0.45f, 0.85f, 0.45f, 1.0f) : ImVec4(0.95f, 0.95f, 0.95f, 1.0f), "%s%s", QuestBook_TaskText(i, k), kd ? "   (done)" : ""); // (the pin box is the only box: no "[ ]" beside it)
                ImGui::PopID();
            }
        }
        ImGui::PopID();
    }
}

// The game window's large tab across the top: true when clicked. The selected one is lit, with an accent underline.
static bool ScreenTab(const char* label, const char* key, bool on, float w, float h) {
    const ImVec4 bg = on ? ImVec4(0.20f, 0.27f, 0.37f, 1.0f) : ImVec4(0.13f, 0.17f, 0.24f, 1.0f);
    ImGui::PushStyleColor(ImGuiCol_Button, bg);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.24f, 0.32f, 0.44f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.28f, 0.37f, 0.50f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_Text, on ? ImVec4(0.95f, 0.96f, 0.98f, 1.0f) : ImVec4(0.62f, 0.68f, 0.76f, 1.0f));
    const ImVec2 p = ImGui::GetCursorScreenPos();
    const bool clicked = ImGui::Button(key && *key ? TextFormat("%s  [%s]", label, key) : label, ImVec2(w, h));
    if (on) ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(p.x + 6, p.y + h - 4), ImVec2(p.x + w - 6, p.y + h), IM_COL32(224, 169, 64, 255));
    ImGui::PopStyleColor(4);
    return clicked;
}

// Placeholder screens for what comes later (their features are separate issues).
static void DrawComing(const char* title, const char* text) {
    ImGui::Dummy(ImVec2(0, 20 * GameUI_Scale()));
    ImGui::TextColored(ImVec4(0.88f, 0.66f, 0.25f, 1.0f), "%s", title);
    ImGui::PushTextWrapPos(ImGui::GetContentRegionAvail().x * 0.7f);
    ImGui::TextWrapped("%s", text);
    ImGui::PopTextWrapPos();
}

// The Team screen (#604). Joined to a server: the team as the server says (members, who's online, their health, the owner),
// with Invite / Leave / Kick / Rename / team chat as /team commands (the server decides), and an invite waiting to accept.
// Single player: a team of one.
static void DrawTeam() {
    const float s = GameUI_Scale();
    ImGui::Dummy(ImVec2(0, 12 * s));
    NetClient* c = gNetClient && gNetClient->state == NetClient::State::Joined ? gNetClient : nullptr;
    if (!c) {
        const int id = Teams_Of(PlayerId_Local());
        const TeamRecord* team = nullptr;
        for (const TeamRecord& t : Teams_All()) if (t.id == id) team = &t;
        ImGui::TextColored(ImVec4(0.88f, 0.66f, 0.25f, 1.0f), "%s", team ? team->name.c_str() : "Your team");
        ImGui::TextWrapped("You are a team of one. Quests and tiers belong to the team. In multiplayer, invite others to share them: Multiplayer on the main menu.");
        return;
    }
    const MsgTeamState& t = c->team;
    auto say = [&](const std::string& line) { NetClient_SendChat(*c, line); };
    bool youOwn = false;
    for (const auto& m : t.members) youOwn = youOwn || (m.you && m.owner);
    if (!t.invite.empty()) { // an invite waiting
        ImGui::TextColored(ImVec4(0.55f, 0.9f, 0.55f, 1.0f), "%s invites you. %s", t.invite.c_str(), t.members.size() <= 1
                           ? "Joining adds your quest progress to the team's."
                           : "Joining shares that team's quest progress; your current team keeps its own.");
        if (ImGui::Button("Accept")) say("/team accept");
        ImGui::SameLine();
        if (ImGui::Button("Decline")) say("/team decline");
        ImGui::Separator();
    }
    ImGui::TextColored(ImVec4(0.88f, 0.66f, 0.25f, 1.0f), "%s", t.name.empty() ? "Your team" : t.name.c_str());
    ImGui::SameLine();
    ImGui::TextDisabled("%d member%s; quests, tiers and recipes are shared", (int)t.members.size(), t.members.size() == 1 ? "" : "s");
    if (ImGui::BeginTable("members", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        for (const auto& m : t.members) {
            ImGui::TableNextRow();
            ImGui::TableNextColumn();
            ImGui::Text("%s%s%s", m.name.c_str(), m.you ? " (you)" : "", m.owner ? "  - owner" : "");
            ImGui::TableNextColumn();
            if (m.online) ImGui::ProgressBar(std::clamp(m.health / 100.0f, 0.0f, 1.0f), ImVec2(160 * s, 0), TextFormat("HP %.0f", m.health));
            else ImGui::TextDisabled("offline");
            ImGui::TableNextColumn();
            ImGui::TableNextColumn();
            if (youOwn && !m.you) { ImGui::PushID(m.name.c_str()); if (ImGui::SmallButton("Kick")) say("/team kick " + m.name); ImGui::PopID(); }
        }
        ImGui::EndTable();
    }
    ImGui::Spacing();
    static char invite[32] = {}, rename[32] = {}, teamSay[160] = {};
    ImGui::SetNextItemWidth(220 * s);
    ImGui::InputTextWithHint("##invite", "player name", invite, sizeof(invite));
    ImGui::SameLine();
    if (ImGui::Button("Invite") && invite[0]) { say(std::string("/team invite ") + invite); invite[0] = 0; }
    if (youOwn) {
        ImGui::SetNextItemWidth(220 * s);
        ImGui::InputTextWithHint("##rename", "a new team name", rename, sizeof(rename));
        ImGui::SameLine();
        if (ImGui::Button("Rename") && rename[0]) { say(std::string("/team name ") + rename); rename[0] = 0; }
    }
    if (t.members.size() > 1 && ImGui::Button("Leave the team")) say("/team leave");
    ImGui::Separator();
    ImGui::SetNextItemWidth(420 * s);
    if (ImGui::InputTextWithHint("##teamsay", "say to your team (Enter)", teamSay, sizeof(teamSay), ImGuiInputTextFlags_EnterReturnsTrue) && teamSay[0]) {
        say(std::string("/t ") + teamSay); teamSay[0] = 0;
    }
    ImGui::TextDisabled("Also from the command line: /team, /team invite NAME, /team accept, /team leave, /t TEXT");
}

// ---- The inventory screen's player (#822): a turning preview, the suit, its module slots and the off hand ----
namespace {
    RenderTexture2D preview{};
    float previewTurn = 0.0f;   // radians: the preview turns slowly; dragging over it turns it by hand
    bool previewDrawn = false;  // rendered this frame (GameUI_PrepareFrame)
    Color SuitColour() { const uint32_t c = settings.mpSuit; return c ? Color{ (unsigned char)(c >> 16), (unsigned char)(c >> 8), (unsigned char)c, 255 } : Color{ 0, 0, 0, 0 }; }
}

bool GameUI_PreviewDrawn() { return previewDrawn; }
void GameUI_PrepareFrame(Scene& world) {
    previewDrawn = false;
    if (!inventoryOpen || screen != GameScreen::Inventory || container.open) return;
    const Vector3 up = Vector3Normalize(player.up);
    if (Vector3Length(player.up) < 0.5f) return; // (not placed yet)
    const float s = GameUI_Scale();
    const int w = (int)(250 * s), h = (int)(360 * s);
    if (preview.id == 0 || preview.texture.width != w || preview.texture.height != h) {
        if (preview.id) UnloadRenderTexture(preview);
        preview = LoadRenderTexture(w, h);
    }
    previewTurn += GameFrameTime() * 0.5f;
    AvatarState a = Avatar_FromPlayer(player, 0xFFFFFF00u);
    a.velocity = { 0, 0, 0 }; a.name.clear();
    a.held = Hotbar_Selected() >= 0 ? Hotbar_Slot(Hotbar_Selected()) : std::string();
    a.headlamp = survival.headlamp;
    a.suitColor = SuitColour();
    // The camera stays in front of where the player faces; the avatar turns in front of it.
    Vector3 f = Vector3Subtract(player.forward, Vector3Scale(up, Vector3DotProduct(player.forward, up)));
    f = Vector3Length(f) > 0.01f ? Vector3Normalize(f) : Vector3Normalize(Vector3CrossProduct(up, { 1, 0, 0 }));
    const Vector3 r = Vector3CrossProduct(up, f);
    a.facing = Vector3Normalize(Vector3Add(Vector3Scale(f, cosf(previewTurn) * -1.0f), Vector3Scale(r, sinf(previewTurn))));
    const Vector3 feet = (player.position - world.renderOrigin).toFloat();
    const Vector3 at = Vector3Add(feet, Vector3Scale(up, 0.92f));
    const Vector3 eye = Vector3Add(Vector3Add(feet, Vector3Scale(up, 1.15f)), Vector3Scale(f, 3.6f));
    a.position = player.position;
    Avatar_DrawPreview(world, a, preview, eye, at, up);
    previewDrawn = true;
}

// The suit slot: the worn suit; holding a suit (picked up from the inventory), a click wears it.
static void SuitSlot(float cell) {
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + cell, p.y + cell);
    if (ImGui::InvisibleButton("suit", ImVec2(cell, cell)) && !OpenIfAlt(Suit_Worn()) && !cursorStack.empty() && Suit_IsSuit(cursorStack.id)) {
        const std::string id = cursorStack.id;
        ReturnCursorStack();
        if (Suit_Wear(id)) Audio_Play(Sfx::Build, 0.4f, 1.2f);
    }
    const bool hover = ImGui::IsItemHovered();
    dl->AddRectFilled(p, q, IM_COL32(30, 36, 48, 230), 4.0f);
    DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y + cell * 0.12f), cell * 0.64f, Suit_Worn());
    dl->AddRect(p, q, hover ? IM_COL32(220, 200, 150, 255) : IM_COL32(110, 100, 70, 255), 4.0f, 0, 1.0f);
    if (hover) {
        const ItemDef* d = itemDB.find(Suit_Worn());
        ImGui::SetTooltip("Worn: %s\n%s", d ? d->name.c_str() : Suit_Worn().c_str(), !cursorStack.empty() && Suit_IsSuit(cursorStack.id) ? "Click: wear the held suit" : "Pick up another suit and click here to wear it");
    }
}

// A module slot: its module; holding a module, a click fits it here (the one there comes out); empty-handed, a click takes it out.
static void ModuleSlot(int i, float cell) {
    ImGui::PushID(1000 + i);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const std::string fitted = i < (int)Suit_Modules().size() ? Suit_Modules()[i] : std::string();
    const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + cell, p.y + cell);
    if (ImGui::InvisibleButton("module", ImVec2(cell, cell)) && !OpenIfAlt(fitted)) {
        if (!cursorStack.empty() && Suit_IsModule(cursorStack.id)) { const std::string id = cursorStack.id; ReturnCursorStack(); Suit_FitInto(i, id); }
        else if (cursorStack.empty() && !fitted.empty()) Suit_Unfit(i);
    }
    const bool hover = ImGui::IsItemHovered();
    dl->AddRectFilled(p, q, IM_COL32(26, 32, 44, 230), 4.0f);
    if (!fitted.empty()) DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y + cell * 0.12f), cell * 0.64f, fitted);
    else { const ImVec2 sz = ImGui::CalcTextSize("+"); dl->AddText(ImVec2(p.x + (cell - sz.x) * 0.5f, p.y + (cell - sz.y) * 0.5f), IM_COL32(110, 120, 140, 255), "+"); }
    dl->AddRect(p, q, hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 76, 90, 255), 4.0f, 0, 1.0f);
    if (hover) {
        const ItemDef* d = fitted.empty() ? nullptr : itemDB.find(fitted);
        if (!fitted.empty()) ImGui::SetTooltip("%s\n%s\nClick: take it out", d ? d->name.c_str() : fitted.c_str(), Suit_ModuleText(fitted));
        else ImGui::SetTooltip("An empty module slot: pick up a suit module and click here");
    }
    ImGui::PopID();
}

static void DrawYou() {
    const float cell = 54.0f, big = 64.0f;
    ImGui::TextUnformatted("You");
    ImGui::Separator();
    ImGui::BeginGroup();
    ImGui::TextDisabled("Suit");
    SuitSlot(big);
    ImGui::Dummy(ImVec2(0.0f, 6.0f));
    ImGui::TextDisabled("Off hand");
    DrawCell(Inventory::kOffHand, big, big);
    ImGui::PushTextWrapPos(ImGui::GetCursorPosX() + 150.0f);
    ImGui::TextDisabled("A shield blocks (hold RMB); a torch lights; food or an oxygen canister: %s uses one", PlayerInput_KeyName(PlayerInput_Binding(Act::OffHandUse)));
    ImGui::PopTextWrapPos();
    ImGui::EndGroup();
    ImGui::SameLine(0.0f, 16.0f);
    // The player, turning; drag over it to turn it.
    if (previewDrawn) {
        const ImVec2 size((float)preview.texture.width, (float)preview.texture.height);
        ImGui::Image((ImTextureID)(intptr_t)preview.texture.id, size, ImVec2(0, 1), ImVec2(1, 0));
        if (ImGui::IsItemActive()) previewTurn += ImGui::GetIO().MouseDelta.x * 0.012f;
        if (ImGui::IsItemHovered()) ImGui::SetTooltip("Drag to turn");
    }
    const int slots = (int)Suit_Modules().size();
    int fitted = 0; for (const std::string& m : Suit_Modules()) fitted += m.empty() ? 0 : 1;
    ImGui::TextDisabled("Suit modules (%d of %d fitted)", fitted, slots);
    for (int i = 0; i < slots; ++i) { if (i) ImGui::SameLine(0.0f, 4.0f); ModuleSlot(i, cell); }
    ImGui::Separator();
    Suit_DrawCarried();
    ImGui::Separator();
    Suit_DrawStats();
}

void GameUI_Draw(Scene& world) {
    uiWorld = &world;
    RecipeBook_Update(); // new recipes: a toast
    RunCraftQueue(); // (#824) queued crafts, window open or not
    if (!inventoryOpen) { SlotRules_Cancel(); ReturnCursorStack(); binned = InvSlot{}; container.open = false; return; } // a stack still held goes back in (#192)
    if (cursorStack.heat > 0.0f) cursorStack.heat = std::max(0.0f, cursorStack.heat - GameFrameTime() / Inventory::kCoolSeconds); // held, it cools as in a cell (#580)
    CheckContainer(world);
    hoverGrid = nullptr; hoverSlot = -1; // (set again as the grids are drawn: #862)
    if (container.open && screen != GameScreen::Inventory) screen = GameScreen::Inventory; // a chest or machine opens beside your items
    const ImGuiIO& io = ImGui::GetIO();
    const float s = GameUI_Scale();
    // The whole screen but a margin, no title bar: the tabs are the window's header.
    const float margin = 24.0f * s;
    ImGui::SetNextWindowPos(ImVec2(margin, margin));
    ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x - 2 * margin, io.DisplaySize.y - 2 * margin));
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.07f, 0.10f, 0.14f, 0.97f));
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoCollapse |
                                   ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoScrollbar;
    if (ImGui::Begin("Game window", nullptr, flags)) {
        ImGui::SetWindowFontScale(s);
        // The tabs (#810's order; Catalog beside Inventory in creative only, #823), their keys from the bindings; at the right
        // end the temporary Help page (until #826) and Close.
        struct Tab { GameScreen s; const char* name; Act key; };
        static const Tab tabs[] = { { GameScreen::Inventory, "Inventory", Act::Inventory }, { GameScreen::Catalog, "Catalog", Act::Count }, { GameScreen::Crafting, "Crafting", Act::Crafting },
                                    { GameScreen::Quests, "Quests", Act::Count }, { GameScreen::Skills, "Skills", Act::Skills }, { GameScreen::Team, "Team", Act::Count } };
        const float tabH = 46.0f * s, tabW = 190.0f * s;
        if (screen == GameScreen::Catalog && !terrainTool.creative) screen = GameScreen::Inventory; // (left creative)
        for (const Tab& t : tabs) {
            if (t.s == GameScreen::Catalog && !terrainTool.creative) continue;
            const char* key = t.key == Act::Count ? "" : PlayerInput_KeyName(PlayerInput_Binding(t.key));
            if (ScreenTab(t.name, key, screen == t.s, tabW, tabH)) { screen = t.s; if (t.s != GameScreen::Inventory) container.open = false; }
            ImGui::SameLine(0.0f, 6.0f * s);
        }
        const float right = ImGui::GetWindowContentRegionMax().x;
        const float smallW = 120.0f * s;
        ImGui::SetCursorPosX(right - 2 * smallW - 6.0f * s);
        if (ScreenTab("Database", "", databaseOpen, smallW, tabH)) databaseOpen = !databaseOpen; // (#826: alt-click an item for its page)
        ImGui::SameLine(0.0f, 6.0f * s);
        if (ScreenTab("Close", "Esc", false, smallW, tabH)) inventoryOpen = false;
        ImGui::Separator();
        ImGui::BeginChild("screen", ImVec2(0, 0), false);
        // Inventory and Crafting show your items on the left; the others use the whole width.
        const bool items = screen == GameScreen::Inventory || screen == GameScreen::Crafting || screen == GameScreen::Catalog;
        if (items) {
            ImGui::Columns(2, "inv", true);
            static float colW = -1.0f;
            if (colW != s) { ImGui::SetColumnWidth(0, (Inventory::kCols * 58.0f + 24.0f) * s); colW = s; }
            ImGui::TextUnformatted("Items");
            ImGui::Separator();
            DrawGrid();
            ImGui::NextColumn();
        }
        switch (screen) {
            case GameScreen::Catalog: DrawCreative(); break; // (#823) every item, free: creative only
            case GameScreen::Inventory:
                if (container.open) DrawContainer();
                else DrawYou(); // (#822) the player, the suit and its modules, the off hand
                break;
            case GameScreen::Crafting: DrawCrafting(); break;
            case GameScreen::Quests:
                DrawTeamProgress(); // (#825)
                if (ImGui::RadioButton("List", !questTree)) questTree = false;
                ImGui::SameLine();
                if (ImGui::RadioButton("Tree", questTree)) questTree = true;
                ImGui::Separator();
                if (questTree) QuestTree_DrawTab(); else DrawQuestList();
                break;
            case GameScreen::Skills:
                DrawComing("Skills are coming (#167)", "Skills level as you use them: mining, building, farming, husbandry, combat, engineering, piloting and more. "
                                                     "Each level gives a point to spend in that skill's perk tree: faster digging, better yields, longer air, new abilities.");
                break;
            case GameScreen::Team: DrawTeam(); break;
            case GameScreen::Count: break;
        }
        if (items) ImGui::Columns(1);
        ImGui::EndChild();
    }
    ImGui::End();
    ImGui::PopStyleColor();
    // The database (#826): a pop-up over the game window, opened by alt-clicking an item (at its page) or the Database
    // button; browse on from there (links, search, Back). Like Factorio's Factoriopedia.
    if (databaseOpen) {
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x * 0.62f, io.DisplaySize.y * 0.7f), ImGuiCond_Appearing);
        ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.52f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.08f, 0.11f, 0.16f, 1.0f)); // opaque: the screen under it doesn't show through
        if (ImGui::Begin("Database (alt-click any item)", &databaseOpen, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings)) {
            ImGui::SetWindowFontScale(s);
            HelpManual_DrawTab();
        }
        ImGui::End();
        ImGui::PopStyleColor();
    }
    // The item rules (#862), once a frame: the mouse and keys with the slot under the mouse (Alt-clicks open the database).
    if (!io.KeyAlt && !databaseOpen) {
        SlotInput in;
        in.lPress = ImGui::IsMouseClicked(ImGuiMouseButton_Left); in.rPress = ImGui::IsMouseClicked(ImGuiMouseButton_Right);
        in.mPress = ImGui::IsMouseClicked(ImGuiMouseButton_Middle);
        in.lDown = ImGui::IsMouseDown(ImGuiMouseButton_Left); in.rDown = ImGui::IsMouseDown(ImGuiMouseButton_Right);
        in.lRelease = ImGui::IsMouseReleased(ImGuiMouseButton_Left); in.rRelease = ImGui::IsMouseReleased(ImGuiMouseButton_Right);
        in.lDouble = ImGui::IsMouseDoubleClicked(ImGuiMouseButton_Left);
        in.wheel = hoverGrid ? io.MouseWheel : 0.0f;
        in.shift = io.KeyShift; in.ctrl = io.KeyCtrl;
        in.hotbarKey = hoverGrid && !io.WantTextInput ? Hotbar_KeySlot() : -1;
        in.offHandKey = !io.WantTextInput && PlayerInput_BindingPressed(Act::SwapOffHand);
        in.dropKey = !io.WantTextInput && PlayerInput_BindingPressed(Act::DropItem);
        in.creative = terrainTool.creative;
        in.wheelInvert = settings.invertSlotWheel;
        std::vector<SlotGrid*> grids = { &playerGrid };
        if (container.open) { grids.push_back(&containerIn); grids.push_back(&containerOut); grids.push_back(&containerFuel); }
        if (screen == GameScreen::Catalog) grids.push_back(&catalogGrid);
        SlotRules_SetHeld(&cursorStack);
        SlotRules_SetGrids(grids);
        SlotRules_Frame(in, hoverGrid, hoverSlot);
    }
    // A click outside every window drops the held stack (right-click: one unit) at the feet (#564).
    if (!cursorStack.empty() && !ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow)) {
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left)) DropCursor(cursorStack.centi);
        else if (ImGui::IsMouseClicked(ImGuiMouseButton_Right)) DropCursor(100);
    }
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
// Items taken from a machine or chest slot onto the cursor (#192), a tick after the click: onto the held stack when it is
// empty or the same item, else into the inventory.
static void CursorSink(const std::string& id, int count, float heat) { // heat: fresh out of a furnace stays hot on the cursor (#580)
    if (count <= 0) return;
    if (cursorStack.empty()) { cursorStack.id = id; cursorStack.centi = (int64_t)count * 100; cursorStack.heat = heat; return; }
    if (cursorStack.id == id) {
        const int64_t add = (int64_t)count * 100;
        cursorStack.heat = (float)((cursorStack.heat * cursorStack.centi + heat * add) / (double)(cursorStack.centi + add));
        cursorStack.centi += add;
        return;
    }
    ++inventory.muted; inventory.addAmount(id, count, heat); --inventory.muted;
}
void GameUI_InitPickupFeed() { inventory.onAdded = OnPickup; Factory::SetCursorSink(CursorSink); }

// (#604) Your teammates on the HUD: each one's name and health (offline ones greyed), on the left under the quests; only
// when joined to a server and in a team of two or more (non-teammates show no health).
static void DrawTeamStrip() {
    const NetClient* c = gNetClient && gNetClient->state == NetClient::State::Joined ? gNetClient : nullptr;
    if (!c || c->team.members.size() < 2) return;
    const float s = GameUI_Scale();
    const int x = (int)(12 * s), w = (int)(150 * s), h = (int)(10 * s), fs = (int)(14 * s), row = fs + h + (int)(8 * s);
    int y = (int)(GetScreenHeight() * 0.36f);
    int n = 0; for (const auto& m : c->team.members) n += !m.you;
    DrawRectangleRounded({ (float)(x - 6 * s), (float)(y - 6 * s), (float)(w + 12 * s), (float)(n * row + 8 * s) }, 0.12f, 6, Fade(BLACK, 0.42f));
    for (const auto& m : c->team.members) {
        if (m.you) continue;
        DrawText(m.name.c_str(), x, y, fs, m.online ? RAYWHITE : Fade(RAYWHITE, 0.45f));
        if (m.online) {
            DrawRectangle(x, y + fs + (int)(2 * s), w, h, Fade(BLACK, 0.55f));
            DrawRectangle(x, y + fs + (int)(2 * s), (int)(w * std::clamp(m.health / 100.0f, 0.0f, 1.0f)), h, Color{ 230, 70, 70, 255 });
        } else DrawText("offline", x + MeasureText(m.name.c_str(), fs) + (int)(8 * s), y, fs, Fade(RAYWHITE, 0.4f));
        y += row;
    }
}

void GameUI_DrawHUD(Scene& world) {
    DrawTeamStrip();
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
    // Jetpack fuel bar (top right): only with a jetpack (a thruster pack and suit power, #217).
    if (player.jetpackAccel > 0.0f) {
        DrawRectangle(w - 220, 70, 200, 12, Fade(BLACK, 0.5f));
        DrawRectangle(w - 220, 70, (int)(200 * player.jetpackFuel), 12, SKYBLUE);
        DrawText("JET", w - 260, 68, 16, RAYWHITE);
    }
    // Under the minimap (#574): where you are (latitude, longitude and height over the body whose sphere of influence you
    // are in; world kilometres out in space), then speed and gravity, then the biome and the weather.
    {
        int ly = 384; // the minimap and its scale line end at ~372
        auto line = [&](const char* t, Color c) { DrawText(t, w - 20 - MeasureText(t, 16), ly, 16, c); ly += 20; }; // right-aligned with the minimap
        bool rotate = false;
        Planetoid* ref = player.nearest ? player.nearest : player.flyReference(world, rotate);
        const DVec3 local = ref ? ref->toLocalD(player.position) : DVec3{};
        const double r = sqrt(local.x * local.x + local.y * local.y + local.z * local.z);
        if (ref && r < 2.0 * ref->radius) {
            const double lat = asin(std::clamp(local.y / std::max(r, 1e-9), -1.0, 1.0)) * RAD2DEG, lon = atan2(local.z, local.x) * RAD2DEG;
            line(TextFormat("%s  %.3f %c  %.3f %c  %.0f m", ref->name.c_str(), fabs(lat), lat >= 0.0 ? 'N' : 'S', fabs(lon), lon >= 0.0 ? 'E' : 'W',
                                (r - ref->radius) * ref->scale), RAYWHITE);
        } else {
            const DVec3 p = player.position;
            line(TextFormat("%.1f  %.1f  %.1f km", p.x / 1000.0, p.y / 1000.0, p.z / 1000.0), RAYWHITE);
        }
        const char* mode = player.flyMode ? "FLY (F)" : player.grounded ? "WALK" : player.hovering ? "HOVER (H)" : player.gravity >= 0.5f ? "AIR" : "SPACE";
        line(TextFormat("%s   g %.1f   %.1f m/s", mode, player.gravity, relSpeed(player.velocity, player.grounded)), RAYWHITE);
        const char* biome = Survival_BiomeName();
        const bool weather = Weather_Here().kind != WeatherKind::Clear && biome[0];
        if (biome[0]) { line(weather ? TextFormat("%s   %s", biome, Weather_Name(Weather_Here().kind)) : biome, Fade(RAYWHITE, 0.8f)); }
        if (player.jetpackAccel > 0.0f) // only with a jetpack: nothing else uses them (#531)
            line(TextFormat("dampeners %s (B)", player.dampeners ? "on" : "off"), player.dampeners ? RAYWHITE : ORANGE);
    }
    if (GetTime() < toastUntil) DrawText(toast.c_str(), 10, GetScreenHeight() - 90, 20, YELLOW);
}

// (#824) the crafting queue from outside (tests; later, e.g. a recipe's "craft" from the database pop-up, #826)
void GameUI_QueueCraft(int recipe, int count) { if (count > 0) craftQueue.push_back({ recipe, count }); }
int GameUI_CraftQueueLength() { int n = 0; for (const QueuedCraft& c : craftQueue) n += c.left; return n; }
void GameUI_SelectRecipe(int recipe) { selectedRecipe = recipe; }

// (#826) the database pop-up from outside the window's clicks (tests; later a key or a link from quests)
void GameUI_OpenDatabase(const std::string& id) { inventoryOpen = true; if (!id.empty()) HelpManual_OpenItem(id); databaseOpen = true; }
bool GameUI_DatabaseOpen() { return inventoryOpen && databaseOpen; }

// (#827) The key hints line along the bottom: from the key bindings (so a rebound key shows as it is), at the interface
// size; Settings > Gameplay > Key hints: Auto (fades out after your first 2 hours played in this world), Always, Off.
static void DrawKeyHints(int x, int y) {
    float alpha = 0.7f;
    if (settings.keyHints == 2) return;
    if (settings.keyHints == 0) {
        const double t = GameSave_PlayTime();
        if (t > 2.0 * 3600.0) alpha *= (float)std::clamp(1.0 - (t - 2.0 * 3600.0) / 900.0, 0.0, 1.0); // over a quarter of an hour
        if (alpha <= 0.01f) return;
    }
    auto K = [](Act a) { return std::string(PlayerInput_KeyName(PlayerInput_Binding(a))); }; // (a copy: the name comes back in one shared buffer)
    std::string line = terrainTool.terraform ? "LMB dig  RMB fill  Ctrl+wheel brush" : "LMB mine/chop  RMB place";
    line += TextFormat("  %s build  1-0 hotbar  %s inventory  %s telescope  %s lamp", K(Act::BuildMode).c_str(), K(Act::Inventory).c_str(), K(Act::Telescope).c_str(), K(Act::Headlamp).c_str());
    if (terrainTool.creative) line += TextFormat("  %s fly", K(Act::Fly).c_str()); // fly is creative only (#531)
    line += TextFormat("  %s save  %s commands", K(Act::Save).c_str(), K(Act::Chat).c_str());
    const float s = GameUI_Scale();
    const int fs = (int)(20 * s);
    DrawText(line.c_str(), x, std::min(y, GetScreenHeight() - fs - (int)(8 * s)), fs, Fade(RAYWHITE, alpha)); // (larger text: kept on screen)
}
void GameUI_InstallKeyHints() { gKeyHintsHook = DrawKeyHints; }
