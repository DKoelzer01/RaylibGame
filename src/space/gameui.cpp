#include "gameui.h"
#include "lightlevel.h"
#include "claims.h"
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
#include "composed.h"
#include "questtree.h"
#include "researchtab.h"
#include "tometab.h" // (#1251)
#include "statstab.h"
#include "navigation.h"
#include "claims.h"
#include "recipebook.h"
#include "craftplan.h"
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
#include "machinetiers.h"
#include "disassembly.h"
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
    int gTestHoverCell = -1; // (#1449) GameUI_TestHover: a cell shown as hovered (its tooltip), for scenarios
    bool inventoryOpen = false;
    std::string toast;
    double toastUntil = 0.0;
    struct Container { bool open = false; uint64_t body = 0; uint32_t plot = 0; int part = -1; PartType type = PartType::Count; } container;
    // The game window (#821): the screen shown; Quests' view; the database pop-up (#826: alt-click an item).
    GameScreen screen = GameScreen::Inventory;
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
    // (#891, the user: "the font and icon scale needs to be at least doubled") Twice the old base at 1080 p, following the
    // screen's height (720 p: 1.33); the Interface size setting multiplies it (0: auto = 100%).
    const float base = 2.0f * std::clamp(GetScreenHeight() / 1080.0f, 0.6f, 2.0f);
    return base * (settings.uiScale > 0.0f ? settings.uiScale : 1.0f);
}
// (#925, the user: the 2x HUD "crowds the screen") The HUD outside the game window keeps the pre-#891 size: the screen's
// height over 1080 p times the Interface size setting, half the menus' scale (the pinned quests too).
float GameUI_HudScale() {
    return std::clamp(GetScreenHeight() / 1080.0f, 0.6f, 2.0f) * (settings.uiScale > 0.0f ? settings.uiScale : 1.0f);
}
// (#891) Every ImGui window at the interface scale: the fonts (FontScaleMain: rasterized at that size by the dynamic atlas,
// so they stay crisp; tooltips and pop-ups too) and the padding and spacing (the style's base sizes x the scale). Before
// each ImGui frame; it changes the style only when the scale does.
void GameUI_ApplyStyleScale() {
    static ImGuiStyle base;
    static bool haveBase = false;
    static float applied = -1.0f;
    if (!ImGui::GetCurrentContext()) return;
    if (!haveBase) { base = ImGui::GetStyle(); haveBase = true; }
    const float s = GameUI_Scale();
    gHudScale = GameUI_HudScale(); // (#925) the engine's HUD text (the crosshair) at the HUD's size
    if (std::fabs(s - applied) < 1e-4f) return;
    ImGuiStyle st = base;
    st.ScaleAllSizes(s);
    st.FontScaleMain = s;
    ImGui::GetStyle() = st;
    applied = s;
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
static bool ContainerTakes(const std::string& id) { // (#1427 L7c) an assembler set to disassemble: what it takes apart
    const int mode = container.type == PartType::LvAssembler ? Factory::View(container.body, container.plot, container.part).mode : 0;
    return Factory::Accepts(container.type, id, mode);
}
static Scene* uiWorld = nullptr; // the world the window was last drawn for (nearby chests)
void GameUI_ShowTab(const char* tab) { // (the older tab names, from the tests and tutorials)
    databaseOpen = false; // (a named screen: the pop-up only when asked for, below)
    const std::string t = tab;
    if (t == "Crafting") GameUI_ShowScreen(GameScreen::Crafting);
    else if (t == "Quest tree" || t == "Quests") GameUI_ShowScreen(GameScreen::Quests); // (#890: one view)
    else if (t == "Research") GameUI_ShowScreen(GameScreen::Research); // (#1192)
    else if (t == "Tome") GameUI_ShowScreen(TomeTab_Shown() ? GameScreen::Tome : GameScreen::Inventory); // (#1251)
    else if (t == "Statistics") GameUI_ShowScreen(GameScreen::Statistics); // (#1199)
    else if (t == "Team") GameUI_ShowScreen(GameScreen::Team);     // (#828)
    else if (t == "Skills") GameUI_ShowScreen(GameScreen::Skills); // (#828)
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

// (#891) Text inside a cell follows the cell, not the window: a quarter of its height, never above the window's font (the
// doubled font in a small cell ran out of it). Drawn with the font at that size (crisp: the dynamic atlas).
// (#889 follow-up) The items grid's cell this frame (DrawItemsPanel): the recipe icons and the You slots are sized from it.
static float gGridCell = 0.0f;
static float gDrawnCellInventory = 0.0f, gDrawnCellCatalog = 0.0f; // (#887) the cell each page last drew with (tests)
void GameUI_TestDrawnCells(float& inventoryCell, float& catalogCell) { inventoryCell = gDrawnCellInventory; catalogCell = gDrawnCellCatalog; }
// (#1148, the user prefers the Catalog page's cells) The one slot size of the game window's item pages (Inventory, Crafting,
// Catalog, a container's player grid; the suit, off-hand and module slots), from the screen area under the tabs: 54 px
// at the interface scale, at most 10 across half its width (the Catalog's column), and the grid, hotbar and bin rows fitting
// its height.
static float SlotCellFor(ImVec2 area) {
    const float gap = 4.0f, line = ImGui::GetTextLineHeightWithSpacing(), sy = ImGui::GetStyle().ItemSpacing.y;
    const float byWidth = (area.x * 0.5f - 24.0f - (Inventory::kCols - 1) * gap) / (float)Inventory::kCols;
    const float byHeight = (area.y - line * 2.0f - sy * 9.0f - 8.0f) / (float)(Inventory::kRows + 2);
    return std::max(24.0f, std::min({ 54.0f * GameUI_Scale(), byWidth, byHeight }));
}
static float CellFontSize(float cellH) { return std::min(ImGui::GetFontSize(), std::max(11.0f, cellH * 0.24f)); }
static ImVec2 CellTextSize(float fs, const char* t, const char* end = nullptr) { return ImGui::GetFont()->CalcTextSizeA(fs, FLT_MAX, 0.0f, t, end); }
static void CellText(ImDrawList* dl, float fs, ImVec2 at, ImU32 col, const char* t, const char* end = nullptr) { dl->AddText(ImGui::GetFont(), fs, at, col, t, end); }
// A count in the top right corner of a cell, on a dark backing.
static void CellCount(ImDrawList* dl, ImVec2 p, ImVec2 q, const char* t, ImU32 col = IM_COL32(255, 255, 255, 255)) {
    const float fs = CellFontSize(q.y - p.y);
    const ImVec2 sz = CellTextSize(fs, t);
    dl->AddRectFilled(ImVec2(q.x - sz.x - 6, p.y + 2), ImVec2(q.x - 2, p.y + sz.y + 2), IM_COL32(0, 0, 0, 160), 2.0f);
    CellText(dl, fs, ImVec2(q.x - sz.x - 4, p.y + 2), col, t);
}

// (#913 review) An item's name along the bottom of its cell, never cut mid-word: one line, a little smaller if that fits it
// (down to 70% of the cell's font); else two lines split at the space nearest the middle; a single long word shrinks to fit.
static void CellName(ImDrawList* dl, ImVec2 p, ImVec2 q, const std::string& name, bool) {
    const float fs = CellFontSize(q.y - p.y), room = q.x - p.x - 8.0f;
    const ImU32 col = IM_COL32(220, 220, 220, 255);
    const float tw = CellTextSize(fs, name.c_str()).x;
    if (tw <= room || tw * 0.7f <= room) { // one line
        const float f = tw <= room ? fs : fs * room / tw;
        CellText(dl, f, ImVec2(p.x + 4, q.y - f - 2.0f), col, name.c_str());
        return;
    }
    size_t split = std::string::npos; // the space nearest the middle
    for (size_t i = 0; i < name.size(); ++i)
        if (name[i] == ' ' && (split == std::string::npos || std::abs((int)i - (int)name.size() / 2) < std::abs((int)split - (int)name.size() / 2))) split = i;
    if (split == std::string::npos) { const float f = std::max(7.0f, fs * room / tw); CellText(dl, f, ImVec2(p.x + 4, q.y - f - 2.0f), col, name.c_str()); return; }
    const std::string a = name.substr(0, split), b = name.substr(split + 1);
    const float wide = std::max(CellTextSize(fs, a.c_str()).x, CellTextSize(fs, b.c_str()).x);
    const float f = std::max(7.0f, std::min(fs * 0.9f, fs * room / wide));
    const float y2 = q.y - f - 2.0f, y1 = y2 - f;
    dl->AddRectFilled(ImVec2(p.x + 2, y1 - 1), ImVec2(q.x - 2, q.y - 1), IM_COL32(0, 0, 0, 110), 2.0f); // (over the icon's foot: readable)
    CellText(dl, f, ImVec2(p.x + 4, y1), col, a.c_str());
    CellText(dl, f, ImVec2(p.x + 4, y2), col, b.c_str());
}

// The stack on the cursor (#192): picked out of its cell, it follows the mouse until it is put down; whatever is still
// held when the inventory closes goes back in.
static InvSlot cursorStack;
// (#1728) The stack came onto the cursor from a machine or chest (CursorSink), which counted it as obtained: coming back in, it
// must not count again. Cleared whenever the cursor is empty (each frame, and when it goes back in).
static bool cursorCredited = false;
static void ReturnCursorStack() {
    if (cursorStack.empty()) { cursorCredited = false; return; }
    ++inventory.muted; // not a pickup (#208)
    if (cursorStack.hasData()) inventory.addInstance(cursorStack.id, cursorStack.data, cursorStack.heat, !cursorCredited); // (#1446) one item, record and all
    else inventory.addAmount(cursorStack.id, cursorStack.amount(), cursorStack.heat, !cursorCredited); // still as hot (#580)
    --inventory.muted;
    cursorStack = InvSlot{};
    cursorCredited = false;
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
// Creative items onto the hand (the Catalog, a middle-click): here at once; joined, the server gives them (#866) and they
// move from the inventory onto the hand when its delta brings them. Returns the centi on the hand now.
static int64_t CreativeTake(const std::string& id, int64_t centi, InvSlot& held) {
    if (centi <= 0 || id.empty() || (!held.empty() && held.id != id)) return 0;
    if (gCreativeGive) {
        const int n = (int)std::max<int64_t>(1, centi / 100);
        gCreativeGive(id, n);
        MachineAct_ExpectOnCursor(id, n);
        return 0;
    }
    if (held.empty()) { held = InvSlot{}; held.id = id; }
    held.centi += centi;
    return centi;
}
// The player's cells (grid, hotbar, off hand). Moves among them change no totals, so they stay here even when joined; what
// goes to or comes from a container is a request (MachineAct_*), what is dropped a pile at the feet (Combat_Drop).
static int64_t ContainerPutRecord(const InvSlot& item); // (#1730) below, once the container's grid exists
struct PlayerGrid : SlotGrid {
    int size() const override { return Inventory::kCells; }
    InvSlot at(int i) const override { return inventory.grid()[i]; }
    int64_t room(int i, const std::string& id) const override {
        const InvSlot& c = inventory.grid()[i];
        if (!c.empty() && (c.id != id || c.hasData())) return 0; // (#1730) an item with its own record never takes more
        return (int64_t)StackLimit(id) * 100 - (c.empty() ? 0 : c.centi);
    }
    int64_t put(int i, const InvSlot& held, int64_t centi) override {
        InvSlot& c = inventory.grid()[i];
        int64_t k = std::min(centi, room(i, held.id));
        if (held.hasData()) { if (!c.empty()) return 0; k = std::min<int64_t>(k, 100); } // (#1730) one item, its record goes along, into an empty cell only
        if (k <= 0) return 0;
        if (c.empty()) { c.id = held.id; c.centi = 0; c.heat = 0.0f; c.data = held.data; }
        c.heat = (float)((c.heat * c.centi + held.heat * k) / (double)(c.centi + k)); // by amount, as merges do (#580)
        c.centi += k;
        return k;
    }
    int64_t take(int i, int64_t centi, InvSlot& held) override {
        InvSlot& c = inventory.grid()[i];
        if (c.empty() || (!held.empty() && (held.id != c.id || held.hasData() || c.hasData()))) return 0; // (#1730) items with a record never stack
        const int64_t k = std::min(centi, c.centi);
        if (held.empty()) { held.id = c.id; held.centi = 0; held.heat = c.heat; held.data = c.data; } // (#1730) the record comes with it
        held.heat = (float)((held.heat * held.centi + c.heat * k) / (double)(held.centi + k));
        held.centi += k;
        c.centi -= k;
        if (c.centi <= 0) c = InvSlot{};
        return k;
    }
    // (#1219, user 2026-11-22) A gem dropped onto a weapon with a free socket is set in it (one gem used). Joined: the server
    // sets it in the inventory it keeps and takes the gem from there, so the held gems go back into the inventory first (the
    // server knows no cursor) and its changes come back as inventory deltas.
    bool combine(int i, InvSlot& held) override {
        const InvSlot& c = inventory.grid()[i];
        if (c.empty() || !c.hasData() || held.centi < 100) return false;
        const std::string set = Composed_SetGem(c.id, c.data, held.id);
        if (set.empty()) return false;
        const std::string gem = held.id, name = Composed_DisplayName(c.id, c.data);
        if (gSetGemRequest) {
            inventory.addAmount(held.id, held.amount(), held.heat);
            held = InvSlot{};
            gSetGemRequest(c.id, c.data, gem);
        } else {
            inventory.setInstanceData(i, set);
            held.centi -= 100;
            if (held.centi <= 0) held = InvSlot{};
        }
        const ItemDef* gd = itemDB.find(gem);
        GameUI_Toast(TextFormat("%s set in the %s", gd ? gd->name.c_str() : gem.c_str(), name.c_str()));
        Audio_Play(Sfx::Craft, 0.5f, 1.3f);
        return true;
    }
    // To the open container (machines hold whole units), else between the grid and the hotbar (the off hand: to the grid).
    void moveOut(int i, int64_t centi) override {
        InvSlot& c = inventory.grid()[i];
        if (c.empty()) return;
        if (GameUI_ContainerOpen()) {
            if (c.hasData()) { if (ContainerPutRecord(c) > 0) c = InvSlot{}; return; } // (#1730) one item, record and all, into a free chest slot
            const std::string id = c.id;
            if (!ContainerTakes(id)) { GameUI_Toast(TextFormat("The %s doesn't take %s", PartName(container.type), itemDB.find(id) ? itemDB.find(id)->name.c_str() : id.c_str())); return; }
            const int64_t amount = std::min(centi, c.centi);
            // (#1102) a chest (here, not as a joined client) keeps the fraction too; a machine takes whole units
            const bool fractions = container.type == PartType::Chest && !MachineAct_Remote();
            const int n = (int)(amount / 100), frac = fractions ? (int)(amount % 100) : 0;
            if (n <= 0 && frac <= 0) return;
            c.centi -= (int64_t)n * 100 + frac;
            if (c.centi <= 0) c = InvSlot{};
            MachineAct_Insert(container.body, container.plot, container.part, id, n, false, frac);
            Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
            return;
        }
        const bool toHotbar = i < Inventory::kGrid;
        const int from = toHotbar ? Inventory::kGrid : 0, to = toHotbar ? Inventory::kSlots : Inventory::kGrid;
        int64_t left = std::min(centi, c.centi);
        for (int pass = 0; pass < 2 && left > 0; ++pass) // onto the same item first, then an empty cell
            for (int k = from; k < to && left > 0; ++k) {
                const InvSlot& t = inventory.grid()[k];
                if (pass == 0 ? (t.empty() || t.id != c.id || t.hasData() || c.hasData()) : !t.empty()) continue;
                const int64_t n = put(k, c, std::min(left, room(k, c.id)));
                left -= n; c.centi -= n;
            }
        if (c.centi <= 0) c = InvSlot{};
    }
    // That item into cell i: from the open chest (it arrives on the next tick, merging with the item's cells), else from the
    // other part (grid <-> hotbar).
    void pullIn(int i, int64_t centi) override {
        const InvSlot c = inventory.grid()[i];
        if (c.empty() || c.hasData()) return; // (#1730) an item with a record is not pulled in by id (one of a kind)
        if (GameUI_ContainerOpen()) {
            if (container.type == PartType::Chest) MachineAct_Take(container.body, container.plot, container.part, c.id, (int)std::max<int64_t>(1, centi / 100));
            return;
        }
        const bool fromHotbar = i < Inventory::kGrid;
        const int from = fromHotbar ? Inventory::kGrid : 0, to = fromHotbar ? Inventory::kSlots : Inventory::kGrid;
        for (int k = from; k < to; ++k) {
            InvSlot& s = inventory.grid()[k];
            if (s.empty() || s.id != c.id || s.hasData()) continue;
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
    void cloneStack(int i, InvSlot& held) override { const InvSlot s = at(i); if (!s.empty()) CreativeTake(s.id, capOf(s.id), held); } // (#866)
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
    std::vector<int> fracs; // (#1102) a chest slot's fraction more (hundredths), as the view has it
    std::vector<std::string> datas; // (#1446 2b) a slot's own record (a composed tool in a chest)
    explicit ContainerGrid(int w) : which(w) {}
    int size() const override { return (int)slots.size(); }
    InvSlot at(int i) const override {
        InvSlot s;
        const int frac = i >= 0 && i < (int)fracs.size() ? fracs[i] : 0;
        if (i >= 0 && i < (int)slots.size() && (slots[i].second > 0 || frac > 0)) { s.id = slots[i].first; s.centi = (int64_t)slots[i].second * 100 + frac; }
        if (!s.empty() && i < (int)datas.size()) s.data = datas[i];
        return s;
    }
    // (#1102) A chest here keeps fractions (0.4 dirt); machines hold whole units, and a joined client's moves are whole too
    // (the fraction stays in its inventory) until the container messages carry fractions.
    bool keepsFractions() const { return container.type == PartType::Chest && which == 0 && !MachineAct_Remote(); }
    int64_t room(int i, const std::string& id) const override {
        if (which == 1 || id.empty()) return 0; // outputs only give
        const InvSlot s = at(i);
        if (which == 2) {
            const auto fuels = Factory::FuelsOf(container.type);
            if (std::find(fuels.begin(), fuels.end(), id) == fuels.end()) return 0;
        } else if (!ContainerTakes(id)) return 0;
        if (container.type == PartType::Chest || which == 2) { if (!s.empty() && (s.id != id || s.hasData())) return 0; return (int64_t)StackLimit(id) * 100 - s.centi; } // (#1730) a slot with a record takes nothing more
        return (int64_t)StackLimit(id) * 100; // a machine routes its inputs (Insert): any input slot click puts it in
    }
    int64_t put(int i, const InvSlot& held, int64_t centi) override {
        if (held.hasData()) { // (#1446 2b) one item with its record: an empty chest slot, or a machine set to disassemble (L7c)
            if (which != 0) return 0;
            if (container.type == PartType::Chest) {
                if (!at(i).empty()) return 0;
                MachineAct_InsertSlot(container.body, container.plot, container.part, held.id, 1, i, 0, held.data);
            } else {
                const MachineView v = Factory::View(container.body, container.plot, container.part);
                const bool anvil = container.type == PartType::Anvil && Composed_Is(held.id); // (#1219) to repair it or take a gem out
                if (v.mode != 1 && !anvil) return 0; // only an assembler taking things apart, or an anvil, takes an item with its own record
                MachineAct_Insert(container.body, container.plot, container.part, held.id, 1, false, 0, held.data);
            }
            Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
            return 100;
        }
        const int64_t amount = std::min({ centi, held.centi, room(i, held.id) });
        const int n = (int)(amount / 100), frac = keepsFractions() ? (int)(amount % 100) : 0;
        if ((n <= 0 && frac <= 0) || room(i, held.id) <= 0) return 0;
        if (which == 2) MachineAct_Insert(container.body, container.plot, container.part, held.id, n, true);
        else if (container.type == PartType::Chest) MachineAct_InsertSlot(container.body, container.plot, container.part, held.id, n, i, frac);
        else MachineAct_Insert(container.body, container.plot, container.part, held.id, n);
        Audio_Play(Sfx::Pickup, 0.4f, 0.8f);
        return (int64_t)n * 100 + frac;
    }
    int64_t take(int i, int64_t centi, InvSlot& held) override {
        const InvSlot s = at(i);
        if (s.empty() || (!held.empty() && held.id != s.id)) return 0;
        if (s.hasData() || held.hasData()) { // (#1446 2b) an item with its record comes alone, onto an empty hand
            if (!held.empty()) return 0;
            MachineAct_TakeSlot(container.body, container.plot, container.part, i, which, 1);
            Audio_Play(Sfx::Pickup, 0.4f, 1.2f);
            return 0; // (next tick, onto the hand with its record)
        }
        if (s.centi < 100) { MachineAct_TakeSlot(container.body, container.plot, container.part, i, which); return 0; } // (#1102) only a fraction: to the inventory (the hand holds whole items here)
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
        if (s.empty() || s.hasData() || which == 1) return; // (#1730) a slot with a record is one of a kind
        const int n = (int)std::min<double>(std::max<int64_t>(1, centi / 100), std::floor(inventory.count(s.id)));
        if (n <= 0 || !inventory.remove(s.id, n)) return;
        InvSlot h = s; h.centi = (int64_t)n * 100;
        if (put(i, h, h.centi) == 0) { ++inventory.muted; inventory.add(s.id, n); --inventory.muted; }
    }
    void cloneStack(int i, InvSlot& held) override { const InvSlot s = at(i); if (!s.empty()) CreativeTake(s.id, capOf(s.id), held); } // (#866)
    void swapCell(int i, int cell) override { // the slot's stack to the inventory, that cell's stack into the slot
        InvSlot& c = inventory.grid()[cell];
        const InvSlot s = at(i);
        if (!c.empty() && room(i, c.id) <= 0 && s.empty()) return;
        if (c.hasData() && !s.empty()) return; // (#1730) a record can only go into an empty slot, and the slot's answer comes a tick later
        if (!s.empty()) MachineAct_TakeSlot(container.body, container.plot, container.part, i, which);
        if (!c.empty() && which != 1) {
            const InvSlot moving = c;
            c = InvSlot{};
            if (put(i, moving, moving.centi) < moving.centi) { ++inventory.muted; if (moving.hasData()) inventory.addInstance(moving.id, moving.data, moving.heat); else inventory.addAmount(moving.id, moving.amount(), moving.heat); --inventory.muted; }
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
        return CreativeTake(ids[i], std::min(centi, capOf(ids[i]) - (held.empty() ? 0 : held.centi)), held);
    }
    void moveOut(int i, int64_t) override { // Shift-click: a stack into the inventory (joined: the server's, #866)
        if (gCreativeGive) { gCreativeGive(ids[i], StackLimit(ids[i])); return; }
        ++inventory.muted; inventory.add(ids[i], StackLimit(ids[i])); --inventory.muted;
    }
    void pullIn(int, int64_t) override {}
    bool infinite() const override { return true; }
};
static PlayerGrid playerGrid;
static ContainerGrid containerIn(0), containerOut(1), containerFuel(2);
// (#1730) Shift-click of an item with its own record onto the open container: the first free chest slot (a machine routes it).
static int64_t ContainerPutRecord(const InvSlot& item) {
    int slot = 0;
    if (container.type == PartType::Chest) {
        slot = -1;
        for (int k = 0; k < containerIn.size() && slot < 0; ++k) if (containerIn.at(k).empty()) slot = k;
        if (slot < 0) return 0;
    }
    return containerIn.put(slot, item, 100);
}
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

static void CursorSink(const std::string& id, int count, float heat, const std::string& data);
void GameUI_TestHover(int cell) { gTestHoverCell = cell; }

// (#1219 tests) A left click on inventory cell `cell` holding `held` on the cursor, through the slot rules; the cursor after it.
InvSlot GameUI_TestClickHolding(const InvSlot& held, int cell) {
    cursorStack = held;
    CellClick(cell, false, false);
    const InvSlot after = cursorStack;
    if (!cursorStack.empty()) ReturnCursorStack();
    return after;
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
    CursorSink("iron_ingot", 4, 1.0f, std::string());
    check(cursorStack.id == "iron_ingot" && cursorStack.heat > 0.99f, "a hot stack onto the cursor stays hot (#580)");
    CursorSink("iron_ingot", 4, 0.0f, std::string());
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
    const bool hover = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem) || i == gTestHoverCell; // (#1449) a scenario shows a cell's tooltip
    const InvSlot& s = grid[i];
    const ImVec2 q(p.x + cw, p.y + cell);
    dl->AddRectFilled(p, q, IM_COL32(30, 32, 36, 220), 4.0f);
    DrawGlow(dl, p, q, Item_Glow(grid[i])); // hot or radiant (#580)
    if (!s.empty()) {
        const ItemDef* d = itemDB.find(s.id);
        DrawItemIcon(dl, ImVec2(p.x + (cw - cell * 0.64f) * 0.5f, p.y), cell * 0.64f, ItemIcons_Key(s.id, s.data));
        const std::string own = Composed_DisplayName(s.id, s.data); // (#1449) "Bronze pickaxe" from its head
        CellName(dl, p, q, !own.empty() ? own : d ? d->name : s.id, hover);
        if (KindOf(s.id) != ItemKind::Tool) {
            const std::string amt = AmountText(s);
            CellCount(dl, p, q, amt.c_str()); // (#891: at the cell's size)
        }
        if (hover && cursorStack.empty()) {
            const std::string parts = Composed_PartsText(s.id, s.data); // (#1449) its own parts, quality, sockets, wear
            if (!parts.empty())
                ImGui::SetTooltip("%s  (%s)\n%s\n%s\n(Shift+click: %s; 1-0: swap with that hotbar slot)", own.c_str(), KindName(KindOf(s.id)), parts.c_str(),
                                  d ? d->description.c_str() : "", GameUI_ContainerOpen() ? TextFormat("into the %s", PartName(container.type)) : i < Inventory::kGrid ? "to the hotbar" : "to the grid");
            else
                ImGui::SetTooltip("%s  (%s)\n%.2f of %d per slot\n%s\n(Shift+click: %s; 1-0: swap with that hotbar slot; wheel: move one)", d ? d->name.c_str() : s.id.c_str(),
                                  KindName(KindOf(s.id)), s.amount(), StackLimit(s.id), d ? d->description.c_str() : "",
                                  GameUI_ContainerOpen() ? TextFormat("into the %s", PartName(container.type)) : i < Inventory::kGrid ? "to the hotbar" : "to the grid");
        }
    }
    DrawDragShare(dl, &playerGrid, i, p, q);
    if (hotbar) CellText(dl, CellFontSize(cell), ImVec2(p.x + 3, p.y + 1), IM_COL32(150, 150, 150, 255), TextFormat("%d", (i - Inventory::kGrid + 1) % 10));
    const ImU32 edge = hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 72, 78, 255);
    dl->AddRect(p, q, edge, 4.0f, 0, 1.0f);
    ImGui::PopID();
}

// cell: a cell's size (px): 54 beside another screen; the Inventory screen sizes it to fill its area (#887).
static void DrawGrid(float cell = 54.0f) {
    const float gap = 4.0f, bin = cell; // (#1148, #1161: the bin and the off hand at the one slot size, side by side)
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (int i = 0; i < Inventory::kSlots; ++i) {
        // The 10 x 5 grid, then the hotbar's 10 cells in a row under it (#209: real slots; drop stacks into them); the off
        // hand and the bin under them (#1161).
        const bool hotbar = i >= Inventory::kGrid;
        if (i == Inventory::kGrid) { ImGui::Dummy(ImVec2(0.0f, 2.0f)); ImGui::TextDisabled("Hotbar (keys 1-0)"); }
        else if (i < Inventory::kGrid ? i % Inventory::kCols != 0 : true) ImGui::SameLine(0.0f, gap);
        DrawCell(i, hotbar ? (Inventory::kCols * cell - gap) / (float)Inventory::kHotbar : cell, cell); // the hotbar row as wide as the grid
    }
    ImGui::Dummy(ImVec2(0.0f, 2.0f));
    // (#1161, the user) The off hand beside the bin, under the hotbar: on every page with the items grid (it was on the
    // Inventory page's You panel only).
    DrawCell(Inventory::kOffHand, cell, cell);
    if (ImGui::IsItemHovered() && inventory.grid()[Inventory::kOffHand].empty())
        ImGui::SetTooltip("Off hand: a shield blocks (hold RMB); a torch lights;\nfood or an oxygen canister: %s uses one", PlayerInput_KeyName(PlayerInput_Binding(Act::OffHandUse)));
    ImGui::SameLine(0.0f, gap);
    { // the bin, under the hotbar (#564)
        const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + bin, p.y + bin);
        if (ImGui::InvisibleButton("bin", ImVec2(bin, bin))) BinClick();
        const bool hover = ImGui::IsItemHovered();
        dl->AddRectFilled(p, q, IM_COL32(48, 26, 26, 220), 4.0f);
        if (!binned.empty()) DrawItemIcon(dl, ImVec2(p.x + bin * 0.18f, p.y + bin * 0.1f), bin * 0.64f, ItemIcons_Key(binned.id, binned.data), true);
        const float fs = CellFontSize(bin); // (#1161: the cells' label size: the window's font ran out of a small cell)
        const ImVec2 sz = ImGui::GetFont()->CalcTextSizeA(fs, FLT_MAX, 0.0f, "Bin");
        dl->AddText(ImGui::GetFont(), fs, ImVec2(p.x + (bin - sz.x) * 0.5f, q.y - sz.y - 2.0f), IM_COL32(220, 160, 150, 255), "Bin");
        dl->AddRect(p, q, hover ? IM_COL32(230, 150, 140, 255) : IM_COL32(110, 60, 60, 255), 4.0f, 0, 1.0f);
        if (hover) ImGui::SetTooltip(!cursorStack.empty() ? "Click: delete the held stack" : !binned.empty() ? "Click: take the deleted stack back" :
                                     "The bin: put a stack here to delete it\n(or click outside the window to drop it on the ground)");
    }
    ImGui::SameLine(0.0f, 8.0f); // (#887: beside the bin, a row saved)
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
    // (#1161) Under the buttons, wrapped (beside them it was squeezed to a word a line once the grid took less width); the
    // stack sizes are the Help tab's, not repeated here.
    const char* note = !cursorStack.empty() ? "Click: put it down  right-click: one  outside the window: drop it" : nChests > 0 ? TextFormat("%d chest%s within reach", nChests, nChests > 1 ? "s" : "") : "";
    if (note[0]) {
        ImGui::PushStyleColor(ImGuiCol_Text, ImGui::GetStyleColorVec4(ImGuiCol_TextDisabled));
        ImGui::TextWrapped("%s", note);
        ImGui::PopStyleColor();
    }
    if (!cursorStack.empty()) { // the held stack follows the mouse, over everything (#192)
        ImDrawList* fg = ImGui::GetForegroundDrawList();
        const ImVec2 m = ImGui::GetMousePos();
        DrawGlow(fg, ImVec2(m.x - 17.0f * GameUI_Scale(), m.y - 17.0f * GameUI_Scale()), ImVec2(m.x + 17.0f * GameUI_Scale(), m.y + 17.0f * GameUI_Scale()), Item_Glow(cursorStack)); // (#580)
        { const float hs = 30.0f * GameUI_Scale(); DrawItemIcon(fg, ImVec2(m.x - hs * 0.5f, m.y - hs * 0.5f), hs, ItemIcons_Key(cursorStack.id, cursorStack.data)); } // (#891: at the interface scale)
        if (KindOf(cursorStack.id) != ItemKind::Tool) fg->AddText(ImVec2(m.x + 10.0f * GameUI_Scale(), m.y + 6.0f * GameUI_Scale()), IM_COL32(255, 255, 255, 255), AmountText(cursorStack).c_str());
    }
}

// The Crafting screen (#824): search and category chips over the recipes you know; a click selects one, and its page
// shows what it makes, what it takes (with what you have, nearby chests counted), Craft 1 / 5 / All, and what it's used
// by. Crafts go in a queue, one every kCraftEvery, each through Crafting_Request (here, or the server when joined). A
// single craft short of something pulls it from nearby chests first (it arrives next tick). In creative everything is
// known and free.
namespace {
    enum class CraftCat { All, CanMake, Tools, Building, Materials, Food, Suit, Count };
    const char* kCraftCatNames[] = { "All", "Craftable", "Tools", "Building", "Materials", "Food", "Suit" }; // (#1149: "Can make" was "Craftable")
    CraftCat CategoryOf(const std::string& id) {
        if (Suit_IsSuit(id) || Suit_IsModule(id)) return CraftCat::Suit;
        if (Survival_IsFood(id)) return CraftCat::Food;
        switch (KindOf(id)) {
            case ItemKind::Tool: return CraftCat::Tools;
            case ItemKind::Block: return CraftCat::Building;
            default: return CraftCat::Materials;
        }
    }
    // group (#928): the steps of one recursive craft (intermediates first) share it; cancelling one cancels them all, and a
    // step that runs short drops the rest. 0: a craft on its own.
    struct QueuedCraft { int recipe; int left; int group = 0; };
    std::vector<QueuedCraft> craftQueue;
    int nextCraftGroup = 1;
    double nextCraftAt = 0.0;
    constexpr double kCraftEvery = 0.3; // s between queued crafts (you see them come)
    int selectedRecipe = -1;
    // (#921) The selected recipe's page is a pop-up beside the recipe clicked: open or not (Esc and GameUI_SelectRecipe set
    // it), whether ImGui showed it last frame (so a click outside, which ImGui closes it on, is seen), a click on another
    // recipe that should keep it open on that one, and the clicked recipe's cell (the anchor) and the page's last height.
    bool recipePageOpen = false, recipePageShown = false, recipePageReopen = false;
    ImVec2 recipeAnchorMin{}, recipeAnchorMax{};
    float recipePageH = 0.0f;
    struct { ImVec2 firstCell{}; int firstRecipe = -1; ImVec2 pageMin{}, pageMax{}, craft1{}; } recipeTest; // (unit tests: where things were drawn)
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
    if (done && q.group) { // (#928) a step of a recursive craft ran short: the rest of it goes too (what was made stays)
        const int g = q.group;
        craftQueue.erase(std::remove_if(craftQueue.begin(), craftQueue.end(), [g](const QueuedCraft& c) { return c.group == g; }), craftQueue.end());
        return;
    }
    if (done || q.left <= 0) craftQueue.erase(craftQueue.begin());
}

// (#928) A plan's steps into the queue, intermediates first, as one group (a single step is a craft on its own).
static void QueuePlan(const CraftPlan& p) {
    const int g = p.steps.size() > 1 ? nextCraftGroup++ : 0;
    for (const CraftStep& s : p.steps) craftQueue.push_back({ s.recipe, s.times, g });
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
    ImGui::SetNextItemWidth(-1.0f); // (#889: the crafting column's width)
    ImGui::InputTextWithHint("##recipesearch", "search recipes", query, sizeof(query));
    const float chipsRight = ImGui::GetCursorScreenPos().x + ImGui::GetContentRegionAvail().x;
    // (#1149, the user) A category chip only when some known recipe is in it (All and Craftable always).
    bool catKnown[(int)CraftCat::Count] = { true, true };
    for (const Recipe& r : recipes) if (RecipeBook_Known(r)) catKnown[(int)CategoryOf(r.output)] = true;
    if (!catKnown[(int)cat]) cat = CraftCat::All;
    for (int c = 0; c < (int)CraftCat::Count; ++c) { // the category chips, wrapping in the column
        if (!catKnown[c]) continue;
        const float w = ImGui::CalcTextSize(kCraftCatNames[c]).x + ImGui::GetStyle().FramePadding.x * 2.0f;
        if (c > 0 && ImGui::GetItemRectMax().x + ImGui::GetStyle().ItemSpacing.x + w < chipsRight) ImGui::SameLine();
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
    // (#928) Factorio-style: a known recipe whose ingredients are not all carried can still be made when the missing ones
    // can be hand-crafted from what is carried (known hand recipes only; what only a machine makes stays missing). Every
    // recipe's plan, refreshed twice a second rather than each frame; which recipes are known, once per refresh.
    static std::vector<char> knownNow, planOk;
    static double planAt = -1.0;
    static std::set<std::string> machineMade;
    if (machineMade.empty()) for (const auto& ri : Factory::RecipeList()) for (const auto& o : ri.out) machineMade.insert(o.first);
    std::map<std::string, int> carried;
    for (const auto& kv : inventory.all()) carried[kv.first] = kv.second;
    const std::function<bool(int)> known = [&](int k) { return k >= 0 && k < (int)knownNow.size() && knownNow[k]; };
    if (GetTime() - planAt > 0.5 || planOk.size() != recipes.size()) {
        planAt = GetTime();
        knownNow.assign(recipes.size(), 0);
        for (size_t k = 0; k < recipes.size(); ++k) knownNow[k] = RecipeBook_Known(recipes[k]);
        planOk.assign(recipes.size(), 0);
        if (!creative)
            for (size_t k = 0; k < recipes.size(); ++k)
                if (knownNow[k] && !canMake(recipes[k])) planOk[k] = CraftPlan_Make(recipes, known, carried, (int)k, 1, machineMade).ok;
    }
    auto canMakeAt = [&](int k) { return canMake(recipes[k]) || (k >= 0 && k < (int)planOk.size() && planOk[k]); };
    // (#1149, the user) A click on a recipe opens its page; a click on its icon again (the page open) crafts one, its
    // intermediates first (#928; joined, the server plans it, #946). It was a double-click.
    auto craftOne = [&](int k) {
        if (k < 0 || k >= (int)recipes.size()) return;
        if (!canMakeAt(k)) { GameUI_Toast("Missing ingredients (the page lists them)"); Audio_Play(Sfx::Error, 0.4f); return; }
        if (creative || canMake(recipes[k])) craftQueue.push_back({ k, 1 });
        else if (!Crafting_RequestRecursive(k, 1)) QueuePlan(CraftPlan_Make(recipes, known, carried, k, 1, machineMade));
    };
    // A pending craft completes once the pulled items have arrived.
    if (pending >= 0 && pending < (int)recipes.size()) {
        const Recipe& r = recipes[pending];
        if (inventory.canCraft(r)) { Crafting_Request(pending, 1); pending = -1; } // (#338) here, or the server when joined
        else if (GetTime() > pendingUntil) { GameUI_Toast("The chests no longer hold enough"); pending = -1; }
    }
    if (selectedRecipe >= (int)recipes.size() || (selectedRecipe >= 0 && !RecipeBook_Known(recipes[selectedRecipe]))) selectedRecipe = -1;
    // (#921) ImGui closed the page last frame (a click outside it): closed, unless that click picked another recipe.
    if (recipePageShown && !ImGui::IsPopupOpen("##recipepage") && !recipePageReopen) recipePageOpen = false;
    recipePageReopen = false;
    if (selectedRecipe < 0) recipePageOpen = false;

    // (#889) A column beside the items: the recipes as icons, a few to a row (greyed: can't make now), filling the height
    // but for the queue under them. (#921) A click on one opens its page in a pop-up beside it.
    const int queueRows = std::min(6, 1 + (int)craftQueue.size());
    const float queueH = queueRows * ImGui::GetFrameHeightWithSpacing() + ImGui::GetStyle().ItemSpacing.y * 2.0f + 4.0f;
    struct CellHit { ImVec2 a, b; int recipe; };
    std::vector<CellHit> cellsDrawn; // (the pop-up's clicks on other recipes)
    int hidden = 0, shown = 0;
    ImGui::BeginChild("recipes", ImVec2(0.0f, std::max(ImGui::GetFrameHeight() * 3.0f, ImGui::GetContentRegionAvail().y - queueH)));
    const float gap = 6.0f;
    // (#889 follow-up, the user: "around the same size as the inventory slots, maybe a little smaller") 0.9 of the grid's cell.
    const float want = gGridCell > 0.0f ? gGridCell * 0.9f : 72.0f * GameUI_Scale();
    const float usable = ImGui::GetContentRegionAvail().x - ImGui::GetStyle().ScrollbarSize;
    const int cols = std::max(3, (int)((usable + gap) / (want + gap)));
    const float cell = std::min(want, (usable - (cols - 1) * gap) / (float)cols);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    for (size_t i = 0; i < recipes.size(); ++i) {
        const Recipe& r = recipes[i];
        const bool known = RecipeBook_Known(r), locked = !known && RecipeBook_LockedByResearch(r); // (#1271 L6b) found, but its research not done: shown locked
        if (!known && !locked) { ++hidden; continue; }
        if (!RecipeBook_Matches(r, q)) continue;
        const bool can = !locked && canMakeAt((int)i); // (#928: intermediates counted)
        if (cat == CraftCat::CanMake && !can) continue;
        if (cat != CraftCat::All && cat != CraftCat::CanMake && CategoryOf(r.output) != cat) continue;
        if (shown % cols != 0) ImGui::SameLine(0.0f, gap);
        ++shown;
        const ItemDef* out = itemDB.find(r.output);
        ImGui::PushID((int)i);
        const ImVec2 p = ImGui::GetCursorScreenPos(), qd(p.x + cell, p.y + cell);
        if (ImGui::InvisibleButton("recipe", ImVec2(cell, cell)) && !locked && !OpenIfAlt(r.output)) { selectedRecipe = (int)i; recipePageOpen = true; }
        const bool hover = ImGui::IsItemHovered(), sel = recipePageOpen && selectedRecipe == (int)i;
        if (ImGui::IsItemVisible()) cellsDrawn.push_back({ p, qd, (int)i });
        if (sel) { recipeAnchorMin = p; recipeAnchorMax = qd; } // (follows the cell as the grid scrolls)
        if (shown == 1) { recipeTest.firstCell = ImVec2((p.x + qd.x) * 0.5f, (p.y + qd.y) * 0.5f); recipeTest.firstRecipe = (int)i; }
        dl->AddRectFilled(p, qd, can ? IM_COL32(34, 40, 36, 230) : IM_COL32(30, 31, 34, 200), 4.0f);
        DrawItemIcon(dl, ImVec2(p.x + 4, p.y + 2), cell - 8, r.output, !can);
        if (r.outputCount > 1) {
            const char* n = TextFormat("x%d", r.outputCount);
            CellText(dl, CellFontSize(cell), ImVec2(qd.x - CellTextSize(CellFontSize(cell), n).x - 4, qd.y - CellFontSize(cell) - 3), IM_COL32(255, 255, 255, can ? 255 : 130), n); // (#891)
        }
        dl->AddRect(p, qd, sel ? IM_COL32(224, 169, 64, 255) : hover ? IM_COL32(220, 220, 220, 255) : can ? IM_COL32(90, 150, 95, 255) : IM_COL32(64, 66, 72, 255),
                    4.0f, 0, sel || hover ? 2.0f : 1.0f);
        if (locked) { CellText(dl, CellFontSize(cell), ImVec2(p.x + 4, p.y + 3), IM_COL32(255, 190, 70, 255), "locked"); } // (#1271 L6b)
        if (hover) ImGui::SetTooltip("%s%s%s%s", out ? out->name.c_str() : r.output.c_str(), r.outputCount > 1 ? TextFormat("  x%d", r.outputCount) : "", locked ? "\n" : "", locked ? RecipeBook_LockReason(r).c_str() : "");
        ImGui::PopID();
    }
    if (shown == 0) ImGui::TextDisabled(cat == CraftCat::CanMake ? "Nothing you can make right now" : "No known recipe matches");
    if (hidden > 0) { // (#891: wrapped in the column)
        ImGui::PushTextWrapPos(0.0f);
        ImGui::TextDisabled("%d recipes still to discover: find new materials (and some need more)", hidden);
        ImGui::PopTextWrapPos();
    }
    if (shown > 0) { ImGui::PushTextWrapPos(0.0f); ImGui::TextDisabled("Click a recipe for its page; click its icon again to craft one"); ImGui::PopTextWrapPos(); }
    ImGui::EndChild();

    // (#921) The selected recipe's page: a pop-up beside its cell (to the right while it stays in the crafting column, else
    // to the left, over the items; kept on screen), a fixed width in characters. Closes on Esc, a click outside or the same
    // recipe again; a click on another recipe moves it there.
    if (recipePageOpen && !ImGui::IsPopupOpen("##recipepage")) ImGui::OpenPopup("##recipepage");
    {
        const ImVec2 disp = ImGui::GetIO().DisplaySize;
        const float w = ImGui::GetFontSize() * 24.0f, m = 8.0f;
        const float colRight = ImGui::GetWindowPos().x + ImGui::GetWindowWidth(); // (the crafting column: this child window)
        float x = recipeAnchorMax.x + m;
        if (x + w > colRight) x = recipeAnchorMin.x - m - w;
        const float y = std::clamp(recipeAnchorMin.y, m, std::max(m, disp.y - recipePageH - m));
        ImGui::SetNextWindowPos(ImVec2(std::max(m, x), y), ImGuiCond_Always);
        ImGui::SetNextWindowSizeConstraints(ImVec2(w, 0.0f), ImVec2(w, disp.y - 2.0f * m));
    }
    recipePageShown = false;
    if (ImGui::BeginPopup("##recipepage")) {
      recipePageShown = true;
      if (!recipePageOpen || selectedRecipe < 0 || ImGui::IsKeyPressed(ImGuiKey_Escape)) { recipePageOpen = false; ImGui::CloseCurrentPopup(); }
      else {
        // A click outside the page: on another recipe, the page goes there (ImGui closes it this frame; it opens again on
        // that one next frame); on the same recipe it crafts one (#1149); anywhere else, it closes.
        if (ImGui::IsMouseClicked(ImGuiMouseButton_Left) && !ImGui::IsWindowHovered(ImGuiHoveredFlags_ChildWindows | ImGuiHoveredFlags_AllowWhenBlockedByActiveItem)) {
            const ImVec2 mp = ImGui::GetIO().MousePos;
            for (const CellHit& c : cellsDrawn)
                if (mp.x >= c.a.x && mp.x < c.b.x && mp.y >= c.a.y && mp.y < c.b.y) {
                    if (c.recipe == selectedRecipe) { craftOne(c.recipe); recipePageReopen = true; } // (#1149: its icon again crafts one; the page stays: ImGui closed it on the click, it opens again)
                    else { selectedRecipe = c.recipe; recipePageReopen = true; recipeAnchorMin = c.a; recipeAnchorMax = c.b; }
                    break;
                }
        }
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
        // (#928) The plan for one: what it makes first from what is carried (empty when the ingredients are all there), and
        // the most that can be made that way (the binary search: twice a second, for the recipe on the page).
        const bool direct = canMake(r);
        const CraftPlan plan1 = creative ? CraftPlan{} : CraftPlan_Make(recipes, known, carried, selectedRecipe, 1, machineMade);
        static int mostFor = -1, mostCached = 0;
        static double mostAt = -1.0;
        if (mostFor != selectedRecipe || GetTime() - mostAt > 0.5) { mostFor = selectedRecipe; mostAt = GetTime(); mostCached = creative ? 99 : CraftPlan_Max(recipes, known, carried, selectedRecipe, 999, machineMade); }
        const int most = mostCached;
        ImGui::TextUnformatted("Takes:");
        for (const auto& in : r.inputs) {
            const ItemDef* id = itemDB.find(in.first);
            const int inv = inventory.count(in.first), all = have(in.first);
            const bool enough = creative || all >= in.second, madeFirst = !enough && plan1.ok; // (#928: crafted from what you carry)
            ImGui::TextColored(enough ? ImVec4(0.55f, 0.9f, 0.55f, 1.0f) : madeFirst ? ImVec4(0.95f, 0.8f, 0.4f, 1.0f) : ImVec4(0.95f, 0.45f, 0.4f, 1.0f), "  %d %s   (have %d%s)%s", in.second,
                               id ? id->name.c_str() : in.first.c_str(), all, all > inv ? TextFormat(", %d in chests", all - inv) : "", madeFirst ? "  - made first" : "");
            if (ImGui::IsItemClicked()) { HelpManual_OpenItem(in.first); databaseOpen = true; } // (#826) an ingredient: its page
        }
        const bool can = direct || plan1.ok;
        // Queues n: in creative, or with everything carried, as it was; else the plan's steps, intermediates first.
        auto queueN = [&](int n) {
            if (creative || inventory.canCraft(r)) { craftQueue.push_back({ selectedRecipe, n }); return; }
            if (Crafting_RequestRecursive(selectedRecipe, n)) return; // (#946) joined: one request, the server plans and crafts it whole
            const CraftPlan p = n == 1 ? plan1 : CraftPlan_Make(recipes, known, carried, selectedRecipe, n, machineMade);
            if (p.ok) QueuePlan(p); else craftQueue.push_back({ selectedRecipe, n });
        };
        auto planTip = [&](int n) { // what a button will make first
            if (creative || !ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenDisabled)) return;
            const CraftPlan p = n == 1 ? plan1 : CraftPlan_Make(recipes, known, carried, selectedRecipe, n, machineMade);
            if (p.ok && p.steps.size() > 1) ImGui::SetTooltip("Crafts %s", CraftPlan_Describe(recipes, p).c_str());
        };
        ImGui::Spacing();
        ImGui::BeginDisabled(!can || pending >= 0);
        const bool craft1 = ImGui::Button("Craft 1");
        { const ImVec2 a = ImGui::GetItemRectMin(), b = ImGui::GetItemRectMax(); recipeTest.craft1 = ImVec2((a.x + b.x) * 0.5f, (a.y + b.y) * 0.5f); }
        planTip(1);
        if (craft1) {
            if (creative || gNetRole == NetRole::Client || inventory.canCraft(r) || plan1.ok) queueN(1);
            else if (uiWorld) { // short in the inventory: pulled from nearby chests first (#338: not when joined)
                bool ok = true;
                for (const auto& in : r.inputs) ok = Storage_Pull(*uiWorld, in.first, in.second - inventory.count(in.first)) && ok;
                if (ok) { pending = selectedRecipe; pendingUntil = GetTime() + 2.0; }
            }
        }
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(most < 5);
        if (ImGui::Button("x5")) queueN(5);
        planTip(5);
        ImGui::EndDisabled();
        ImGui::SameLine();
        ImGui::BeginDisabled(most < 1);
        if (ImGui::Button(creative ? "x20" : TextFormat("All (%d)", most))) queueN(creative ? 20 : most);
        planTip(most);
        ImGui::EndDisabled();
        if (!can) {
            if (!plan1.machineOnly.empty()) { // (#928) what a machine must make first: hand crafting can't
                std::string names;
                for (const std::string& m : plan1.machineOnly) { const ItemDef* d = itemDB.find(m); names += (names.empty() ? "" : ", ") + (d ? d->name : m); }
                ImGui::PushTextWrapPos(0.0f); ImGui::TextDisabled("Missing items; a machine makes %s", names.c_str()); ImGui::PopTextWrapPos();
            } else ImGui::TextDisabled("Missing items");
        } else if (pending >= 0) ImGui::TextDisabled("Waiting for items from the chests...");
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
      }
      recipePageH = ImGui::GetWindowHeight();
      recipeTest.pageMin = ImGui::GetWindowPos();
      recipeTest.pageMax = ImVec2(recipeTest.pageMin.x + ImGui::GetWindowWidth(), recipeTest.pageMin.y + recipePageH);
      ImGui::EndPopup();
    }
    // The queue: what's being crafted, next first.
    ImGui::BeginChild("craftqueue", ImVec2(0.0f, 0.0f), true);
    if (craftQueue.empty()) ImGui::TextDisabled("Queue: empty"); else ImGui::TextUnformatted("Queue:");
    ImGui::PushTextWrapPos(0.0f); // (#946) long step labels wrap in the column ("(craftin..." was cut)
    for (size_t k = 0; k < craftQueue.size(); ++k) {
        const QueuedCraft& c = craftQueue[k];
        const ItemDef* d = c.recipe >= 0 && c.recipe < (int)recipes.size() ? itemDB.find(recipes[c.recipe].output) : nullptr;
        // (#928) A recursive craft's intermediate steps sit indented under it (the last step of a group is the item asked
        // for); cancelling any step cancels the group (nothing is reserved: what was made, and every ingredient not yet
        // used, stays in the inventory).
        const bool step = c.group && k + 1 < craftQueue.size() && craftQueue[k + 1].group == c.group;
        ImGui::PushID(20000 + (int)k);
        if (ImGui::SmallButton("x")) {
            if (c.group) { const int g = c.group; craftQueue.erase(std::remove_if(craftQueue.begin(), craftQueue.end(), [g](const QueuedCraft& e) { return e.group == g; }), craftQueue.end()); }
            else craftQueue.erase(craftQueue.begin() + (long)k);
            ImGui::PopID(); break;
        }
        if (ImGui::IsItemHovered() && c.group) ImGui::SetTooltip("Cancels this craft and the steps it makes first");
        ImGui::SameLine();
        if (step) ImGui::TextDisabled("- %s  x%d (first)%s", d ? d->name.c_str() : "?", c.left * (c.recipe >= 0 && c.recipe < (int)recipes.size() ? std::max(1, recipes[c.recipe].outputCount) : 1), k == 0 ? "  (crafting)" : "");
        else ImGui::Text("%s  x%d%s", d ? d->name.c_str() : "?", c.left, k == 0 ? "  (crafting)" : "");
        ImGui::PopID();
    }
    ImGui::PopTextWrapPos();
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
    containerIn.slots = v.inSlots; containerIn.fracs = v.inCenti; containerIn.datas = v.inData; // (#1446 2b)
    containerOut.slots = v.outSlots; containerOut.datas = v.outData; containerFuel.slots = { v.fuelSlot }; // (#1499; #1446 2c)
    SlotRules_SetGrids({ &playerGrid, &containerIn, &containerOut, &containerFuel });
    RulesClick(&ContainerPart(which), i, right, shift);
}

static void DrawSlots(const std::vector<std::pair<std::string, int>>& slots, int which, int cols, bool hot = false, const std::vector<int>* fracs = nullptr, const std::vector<std::string>* datas = nullptr) { // fracs (#1102): a chest's fractions; datas (#1446): each slot's record // which: 0 inputs, 1 outputs, 2 fuel; hot: its items glow (#580)
    const bool output = which == 1;
    ImDrawList* dl = ImGui::GetWindowDrawList();
    const float gap = 4.0f, cell = std::min(54.0f * GameUI_Scale(), (ImGui::GetContentRegionAvail().x - (cols - 1) * gap) / (float)cols); // (#891: at the scale, fitting the column)
    ContainerGrid& grid = ContainerPart(which);
    grid.slots = slots; // (this frame's view: what the rules see)
    grid.fracs = fracs ? *fracs : std::vector<int>{};
    grid.datas = datas ? *datas : std::vector<std::string>{}; // (#1446 2c) the drawn screen too, not only SlotClick: an item with data takes as one
    for (int i = 0; i < (int)slots.size(); ++i) {
        if (i % cols != 0) ImGui::SameLine(0.0f, gap);
        ImGui::PushID(i + which * 1000);
        const ImVec2 p = ImGui::GetCursorScreenPos(), q(p.x + cell, p.y + cell);
        ImGui::InvisibleButton("slot", ImVec2(cell, cell), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonRight | ImGuiButtonFlags_MouseButtonMiddle);
        const InvSlot here = grid.at(i); // (#1102) whole items and a chest's fraction
        if (ImGui::IsItemClicked()) OpenIfAlt(!here.empty() ? here.id : std::string()); // (#826; the rules skip alt-clicks)
        NoteHover(&grid, i); // clicks, drags, keys and the wheel: the item rules (#862), once a frame in GameUI_Draw
        const bool hover = ImGui::IsItemHovered(ImGuiHoveredFlags_AllowWhenBlockedByActiveItem);
        dl->AddRectFilled(p, q, which == 2 ? IM_COL32(52, 30, 20, 220) : output ? IM_COL32(40, 36, 28, 220) : IM_COL32(30, 32, 36, 220), 4.0f);
        if (!here.empty()) DrawGlow(dl, p, q, Item_GlowOf(here.id, hot ? 1.0f : 0.0f)); // fresh out of a furnace, or radiant (#580)
        if (!here.empty()) {
            DrawItemIcon(dl, ImVec2(p.x + cell * 0.18f, p.y), cell * 0.64f, ItemIcons_Key(here.id, here.data));
            const ItemDef* d = itemDB.find(here.id);
            CellName(dl, p, q, d ? d->name : here.id, hover);
            // the amount as the player inventory shows it: "12", "2.5", "0.4"
            const char* amt = here.centi % 100 == 0 ? TextFormat("%d", (int)(here.centi / 100)) : here.centi % 10 == 0 ? TextFormat("%.1f", here.centi / 100.0) : TextFormat("%.2f", here.centi / 100.0);
            CellCount(dl, p, q, amt); // (#891)
            if (hover && cursorStack.empty()) ImGui::SetTooltip("%s x%s\n%s\n(click: pick up  right-click: half  Shift+click: to the inventory  wheel: move one)", d ? d->name.c_str() : here.id.c_str(), amt, d ? d->description.c_str() : "");
        }
        DrawDragShare(dl, &grid, i, p, q);
        dl->AddRect(p, q, hover ? IM_COL32(200, 200, 200, 255) : IM_COL32(70, 72, 78, 255), 4.0f, 0, 1.0f);
        ImGui::PopID();
    }
}

// (#1427 L7d) The open machine's part (its variant step is its material), its voltage tier (-1: not a tiered machine).
static const Part* ContainerPart() {
    Planetoid* b = Planetoid::find(container.body);
    Plot* pl = b && b->plots ? b->plots->byId(container.plot) : nullptr;
    return pl && container.part >= 0 && container.part < (int)pl->parts.size() && pl->parts[container.part].type == container.type ? &pl->parts[container.part] : nullptr;
}
static int ContainerTier() {
    const Part* p = ContainerPart();
    return p && MachineTiers_Tiered(p->type) ? MachineTiers_PartTier(p->type, p->material, 0) : -1;
}
static void DrawContainer() {
    const MachineView v = Factory::View(container.body, container.plot, container.part);
    { // (#1427 L7d) a variant by its own name, a tiered machine with its tier and voltage
        const Part* p = ContainerPart();
        const ItemDef* d = p && IsMachine(p->type) && p->material > 0 ? itemDB.find(PartItem(*p)) : nullptr;
        ImGui::Text("%s", d ? d->name.c_str() : PartName(container.type));
    }
    ImGui::SameLine();
    if (v.valid) ImGui::TextColored(v.status == MachineStatus::Working ? ImVec4(0.5f, 0.9f, 0.5f, 1.0f) : ImVec4(0.9f, 0.8f, 0.5f, 1.0f), "  %s", MachineStatusName(v.status));
    ImGui::SameLine(ImGui::GetContentRegionAvail().x - 60.0f);
    if (ImGui::Button("Close")) { container.open = false; return; }
    if (const int t = ContainerTier(); t >= 0) { // (#1427 L7d) its tier and voltage on a line of their own (the title row: name, status, Close)
        const Color c = MachineTiers_Color(t);
        ImGui::TextColored(ImVec4(c.r / 255.0f, c.g / 255.0f, c.b / 255.0f, 1.0f), "%s, %d V", MachineTiers_Name(t), MachineTiers_Volts(t));
    }
    ImGui::Separator();
    if (!v.valid) { ImGui::TextDisabled("..."); return; }
    if (v.duration > 0) ImGui::ProgressBar((float)v.progress / (float)v.duration, ImVec2(300.0f, 0.0f), TextFormat("%d%%", 100 * v.progress / v.duration));
    if (v.energyMax > 0) ImGui::ProgressBar((float)v.energy / (float)v.energyMax, ImVec2(300.0f, 0.0f), TextFormat("%d / %d EU", v.energy, v.energyMax));
    else if (v.energy > 0 && container.type != PartType::Chest) ImGui::Text("Fuel / charge: %d", v.energy);
    const bool chest = container.type == PartType::Chest;
    ImGui::TextUnformatted(chest ? "Storage" : "Input");
    DrawSlots(v.inSlots, 0, chest ? Inventory::kCols : 4, false, &v.inCenti, &v.inData); // a chest lines up with the player's grid (#495)
    if (v.hasFuelSlot) { // #194: fuel has its own slot (a recipe input that is also fuel is never burnt)
        ImGui::TextUnformatted("Fuel");
        DrawSlots({ v.fuelSlot }, 2, 1);
    }
    if (v.canLight) { // (#1103) a campfire lights or goes out whether or not it has something to cook
        const bool burning = v.status == MachineStatus::Working;
        const bool fuel = v.fuelSlot.second > 0 || v.energy > 0;
        if (!burning && !fuel) ImGui::BeginDisabled();
        if (ImGui::Button(burning ? "Put out" : "Light")) MachineAct_SetFire(container.body, container.plot, container.part, !burning);
        if (!burning && !fuel) { ImGui::EndDisabled(); ImGui::SameLine(); ImGui::TextDisabled("needs fuel"); }
        else if (!burning && v.fire == kFireOut) { ImGui::SameLine(); ImGui::TextDisabled("put out: lights again with something to cook"); }
    }
    if (const PartType next = MachineTiers_Successor(container.type); next != PartType::Count) { // (#1427 L7e) its job from LV on
        const ItemDef* d = itemDB.find(PartCostItem(next));
        ImGui::TextDisabled("From LV on, the %s does this job (and has tiers).", d ? d->name.c_str() : PartName(next));
    }
    if (container.type == PartType::LvAssembler) { // (#1271 L6c) the mode switch
        if (ImGui::RadioButton("Assemble", v.mode == 0) && v.mode != 0) MachineAct_SetMode(container.body, container.plot, container.part, 0);
        ImGui::SameLine();
        if (ImGui::RadioButton("Disassemble", v.mode == 1) && v.mode != 1) MachineAct_SetMode(container.body, container.plot, container.part, 1);
        if (v.mode == 1) { // (#1427 L7c) what it gives back for each input, at its tier's yield
            ImGui::PushTextWrapPos(0.0f); // (the panel is narrow at small windows: wrap, never clip)
            ImGui::TextDisabled("Takes an item apart into what its recipe used (a tool or weapon: its parts' materials), at a yield by the item's tier against this assembler's.");
            const int myTier = ContainerTier();
            for (size_t i = 0; i < v.inSlots.size(); ++i) {
                const std::string& id = v.inSlots[i].first;
                if (id.empty() || v.inSlots[i].second <= 0) continue;
                DisassemblyEntry e;
                const std::string data = i < v.inData.size() ? v.inData[i] : std::string();
                const bool ok = !data.empty() ? Disassembly_Composed(id, data, e) : Disassembly_Find(id) ? (e = *Disassembly_Find(id), true) : false;
                const ItemDef* d = itemDB.find(id);
                const char* name = d ? d->name.c_str() : id.c_str();
                if (!ok) { ImGui::Text("%s: can't be taken apart", name); continue; }
                const int y = MachineTiers_YieldPercent(std::max(0, myTier), e.tier);
                if (y <= 0) ImGui::Text("%s: too advanced for this assembler (0%%)", name);
                else if (e.keepPct < 100) ImGui::Text("%s: %d%% back (worn, %d%% left)", name, y * e.keepPct / 100, e.keepPct);
                else ImGui::Text("%s: %d%% back", name, y);
            }
            ImGui::PopTextWrapPos();
        }
    }
    { // (#1271 L6c) a machine of several recipes: what it makes ("Any": the first its inputs allow)
        // (#1427 L7d) each output with its recipe's tier: one above the machine's can't run here; one under it is overclocked
        std::vector<std::string> outs;
        std::vector<int> outTier;
        for (const Factory::RecipeInfo& r : Factory::RecipeList())
            if (r.machine == container.type && !r.out.empty() && std::find(outs.begin(), outs.end(), r.out[0].first) == outs.end()) { outs.push_back(r.out[0].first); outTier.push_back(r.tier); }
        const int myTier = ContainerTier();
        auto tierNote = [&](int rt) -> std::string {
            if (rt < 0 || myTier < 0) return "";
            const int steps = MachineTiers_Steps(myTier, rt);
            if (steps < 0) return TextFormat("  (%s: needs a %s machine)", MachineTiers_Name(rt), MachineTiers_Name(rt));
            if (steps == 0) return TextFormat("  (%s)", MachineTiers_Name(rt));
            int speed = 1; for (int s = 0; s < steps; ++s) speed *= MachineTiers_Overclock().speed;
            return TextFormat("  (%s, overclocked: x%d speed, x%d EU/t)", MachineTiers_Name(rt), speed, MachineTiers_Eu(1, steps));
        };
        const bool steered = container.type == PartType::Miner || (container.type >= PartType::LvLab && container.type <= PartType::HvLab)
                          || container.type == PartType::SteamLab || container.type == PartType::ManaSiphon; // (the game picks these)
        if (outs.size() > 1 && !steered && !chest) {
            auto label = [](const std::string& id) { const ItemDef* d = itemDB.find(id); return d ? d->name : id; };
            const std::string now = v.pick.empty() ? std::string("Any (what its inputs allow)") : label(v.pick);
            ImGui::SetNextItemWidth(std::min(300.0f, ImGui::GetContentRegionAvail().x - ImGui::CalcTextSize("Makes").x - ImGui::GetStyle().ItemInnerSpacing.x)); // ("Makes" stays in view on a narrow screen)
            if (ImGui::BeginCombo("Makes", now.c_str())) {
                if (ImGui::Selectable("Any (what its inputs allow)", v.pick.empty())) MachineAct_SetPick(container.body, container.plot, container.part, "");
                for (size_t k = 0; k < outs.size(); ++k) {
                    const std::string& o = outs[k];
                    const bool tooHigh = myTier >= 0 && outTier[k] >= 0 && MachineTiers_Steps(myTier, outTier[k]) < 0;
                    if (tooHigh) ImGui::BeginDisabled();
                    if (ImGui::Selectable((label(o) + tierNote(outTier[k]) + "##" + o).c_str(), v.pick == o)) MachineAct_SetPick(container.body, container.plot, container.part, o);
                    if (tooHigh) ImGui::EndDisabled();
                }
                ImGui::EndCombo();
            }
        }
    }
    if (v.waterMax > 0) { // #229: boilers need water and heat up before full steam
        ImGui::ProgressBar((float)v.water / (float)v.waterMax, ImVec2(300.0f, 0.0f), TextFormat("Water %d / %d L", v.water, v.waterMax));
        ImGui::ProgressBar((float)v.heat / (float)std::max(1, v.heatMax), ImVec2(300.0f, 0.0f), TextFormat("Heat %d%%", 100 * v.heat / std::max(1, v.heatMax)));
    }
    if (!chest) {
        ImGui::TextUnformatted("Output");
        DrawSlots(v.outSlots, 1, 4, v.hotOutputs, nullptr, &v.outData);
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
        if (chest) for (int i = 0; i < (int)v.inSlots.size(); ++i) { if (v.inSlots[i].second > 0 || (i < (int)v.inCenti.size() && v.inCenti[i] > 0)) MachineAct_TakeSlot(container.body, container.plot, container.part, i, 0); }
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
    const Vector3 me = b->toLocalOf(player.position, world.renderOrigin);
    if (Vector3Distance(c, me) * b->scale > 9.0f) { container.open = false; GameUI_Toast("Too far from it"); }
}

// Creative mode: free building, crafting and blueprints, and any item into the inventory.
static void DrawCreative() {
    ImGui::TextColored(ImVec4(1.0f, 0.8f, 0.3f, 1.0f), "Creative: build, craft and place blueprints for free (/gamemode survival to leave)"); // (#823)
    static char filter[64] = {};
    ImGui::SetNextItemWidth(220.0f);
    ImGui::InputTextWithHint("##filter", "search items", filter, sizeof(filter));
    ImGui::SameLine();
    if (ImGui::Button("Clear inventory")) { if (gCreativeGive) gCreativeGive("", 0); else inventory.clear(); } // (joined: the server's, #866)
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
    // (#887, the user: "Catalog slots are larger than inventory slots") the one slot size of the item pages (#1148,
    // SlotCellFor); a fixed 48 px at the interface scale was larger than the Inventory's wherever that fits the screen smaller
    const float cell = gGridCell > 0.0f ? gGridCell : 48.0f * GameUI_Scale(), gap = 4.0f;
    gDrawnCellCatalog = cell; // (#887 test)
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

// The game window's large tab across the top: true when clicked. The selected one is lit, with an accent underline.
// (#1206) the screen tabs' font: one size for the row, the largest at which every name fits (down to 70%; ScreenTabFont)
static float screenTabFont = 1.0f;
static void ScreenTabFont(const char* const* labels, int n, float w) {
    const float pad = ImGui::GetStyle().FramePadding.x * 2.0f;
    screenTabFont = 1.0f;
    for (int i = 0; i < n; ++i) screenTabFont = std::min(screenTabFont, (w - pad) / std::max(1.0f, ImGui::CalcTextSize(labels[i]).x));
    screenTabFont = std::max(0.7f, screenTabFont);
}
static bool ScreenTab(const char* label, const char* key, bool on, float w, float h, bool rowFont = false) {
    const ImVec4 bg = on ? ImVec4(0.20f, 0.27f, 0.37f, 1.0f) : ImVec4(0.13f, 0.17f, 0.24f, 1.0f);
    ImGui::PushStyleColor(ImGuiCol_Button, bg);
    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.24f, 0.32f, 0.44f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0.28f, 0.37f, 0.50f, 1.0f));
    ImGui::PushStyleColor(ImGuiCol_Text, on ? ImVec4(0.95f, 0.96f, 0.98f, 1.0f) : ImVec4(0.62f, 0.68f, 0.76f, 1.0f));
    const ImVec2 p = ImGui::GetCursorScreenPos();
    // (#891) with its key when both fit the tab, else the name alone (a large font in a narrow tab clipped "Skills [K")
    // (#1206) the row's font when its names don't all fit (8 tabs in creative): smaller, rather than clipped names
    const bool shrink = rowFont && screenTabFont < 1.0f;
    if (shrink) ImGui::PushFont(nullptr, ImGui::GetStyle().FontSizeBase * screenTabFont); // (FontSizeBase: GetFontSize is scaled already)
    const char* withKey = key && *key ? TextFormat("%s  [%s]", label, key) : label;
    const float pad = ImGui::GetStyle().FramePadding.x * 2.0f;
    const bool fits = ImGui::CalcTextSize(withKey).x + pad <= w;
    const bool clicked = ImGui::Button(TextFormat("%s##%s", fits ? withKey : label, label), ImVec2(w, h));
    if (shrink) ImGui::PopFont();
    if (on) ImGui::GetWindowDrawList()->AddRectFilled(ImVec2(p.x + 6, p.y + h - 4), ImVec2(p.x + w - 6, p.y + h), IM_COL32(224, 169, 64, 255));
    ImGui::PopStyleColor(4);
    return clicked;
}

// Placeholder screens for what comes later (their features are separate issues).
// The Skills screen (#883, the design on #167): every skill on the left (level, progress, free points), the selected one's
// perk tree on the right in its three tiers (opening at levels 5, 15 and 30), each perk with its ranks and what one rank
// does; Learn spends the skill's points. A retraining manual takes back every point of one tree (Q1).
static int skillSelected = 0;
static void DrawSkills() {
    SkillState& st = LocalSession().skills;
    const auto& all = Skills_All();
    if (all.empty()) { ImGui::TextDisabled("No skills loaded (resources/data/skills.txt)"); return; }
    skillSelected = std::clamp(skillSelected, 0, (int)all.size() - 1);
    const float fh = ImGui::GetFontSize();
    ImGui::BeginChild("skilllist", ImVec2(fh * 17.0f, 0.0f), true);
    for (int i = 0; i < (int)all.size(); ++i) {
        const SkillDef& d = all[i];
        const double xp = st.xp.count(d.id) ? st.xp.at(d.id) : 0.0;
        const int free = Skills_PointsFree(st, d.id);
        ImGui::PushID(i);
        const ImVec2 p = ImGui::GetCursorScreenPos();
        if (ImGui::Selectable("##skill", skillSelected == i, 0, ImVec2(0.0f, fh * 2.4f))) skillSelected = i;
        ImDrawList* dl = ImGui::GetWindowDrawList();
        const float w = ImGui::GetContentRegionAvail().x;
        dl->AddText(ImVec2(p.x + 4, p.y + 2), IM_COL32(230, 230, 230, 255), d.name.c_str());
        dl->AddText(ImVec2(p.x + w * 0.62f, p.y + 2), IM_COL32(255, 210, 90, 255), TextFormat("Lv %d", Skills_LevelOf(xp)));
        if (free > 0) dl->AddText(ImVec2(p.x + w * 0.8f, p.y + 2), IM_COL32(120, 230, 140, 255), TextFormat("+%d", free));
        const ImVec2 b0(p.x + 4, p.y + fh * 1.45f), b1(p.x + w - 4, p.y + fh * 1.45f + fh * 0.55f);
        dl->AddRectFilled(b0, b1, IM_COL32(40, 44, 52, 255), 3.0f);
        dl->AddRectFilled(b0, ImVec2(b0.x + (b1.x - b0.x) * (float)Skills_Progress(xp), b1.y), IM_COL32(230, 180, 60, 255), 3.0f);
        ImGui::PopID();
    }
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("skilltree", ImVec2(0.0f, 0.0f), true);
    const SkillDef& d = all[skillSelected];
    const double xp = st.xp.count(d.id) ? st.xp.at(d.id) : 0.0;
    const int level = Skills_LevelOf(xp);
    double into = xp; for (int k = 1; k <= level; ++k) into -= Skills_XpForLevel(k);
    ImGui::TextColored(ImVec4(1.0f, 0.82f, 0.35f, 1.0f), "%s  level %d", d.name.c_str(), level);
    ImGui::SameLine();
    ImGui::TextDisabled("  %.0f / %.0f XP to level %d", into, Skills_XpForLevel(level + 1), level + 1);
    { // (#828) the points on a line of their own (the header ran out of the pane at 1080p)
        const int pts = Skills_PointsFree(st, d.id);
        ImGui::TextColored(pts > 0 ? ImVec4(0.55f, 0.9f, 0.6f, 1.0f) : ImVec4(0.6f, 0.6f, 0.65f, 1.0f), "%d point%s to spend", pts, pts == 1 ? "" : "s");
        ImGui::SameLine(); ImGui::TextDisabled("   world XP rate %gx", worldXpRate);
    }
    std::string earns;
    for (const auto& [ev, x] : d.xp) earns += TextFormat("%s%s %g", earns.empty() ? "" : ", ", ev.c_str(), x);
    ImGui::TextDisabled("XP from: %s", earns.empty() ? "nothing yet" : earns.c_str());
    ImGui::Separator();
    // (#828) The perk tree as nodes, a row per tier (the plan's grid): learned (green), can learn now (amber: click it),
    // locked or not yet (grey). Each node: its name, ranks and effect; hover for its cost, what it needs and why not.
    bool any = false;
    const int kNodeCols = 3; // (4 left names and their ranks overlapping at 1080p)
    const float gap = fh * 0.6f, nodeW = (ImGui::GetContentRegionAvail().x - gap * (kNodeCols - 1)) / kNodeCols, nodeH = fh * 5.0f;
    for (int tier = 1; tier <= 3; ++tier) {
        const int need = kSkillTierLevel[tier - 1];
        ImGui::TextColored(level >= need ? ImVec4(0.6f, 0.85f, 1.0f, 1.0f) : ImVec4(0.5f, 0.5f, 0.55f, 1.0f), "Tier %d  (level %d)%s", tier, need, level >= need ? "" : "  - locked");
        int col = 0;
        for (const PerkDef& p : Skills_Perks()) {
            if (p.skill != d.id || p.tier != tier) continue;
            any = true;
            const int rank = st.ranks.count(p.id) ? st.ranks.at(p.id) : 0;
            const std::string why = Skills_CanBuy(st, p.id);
            const bool learned = rank > 0, full = rank >= p.ranks, canLearn = why.empty();
            if (col > 0) ImGui::SameLine(0.0f, gap);
            ImGui::PushID(p.id.c_str());
            const ImVec2 a = ImGui::GetCursorScreenPos(), b(a.x + nodeW, a.y + nodeH);
            const bool clicked = ImGui::InvisibleButton("##node", ImVec2(nodeW, nodeH));
            const bool hovered = ImGui::IsItemHovered();
            ImDrawList* dl = ImGui::GetWindowDrawList();
            const ImU32 edge = canLearn ? IM_COL32(230, 180, 60, 255) : learned ? IM_COL32(110, 210, 130, 255) : IM_COL32(80, 88, 100, 255);
            const ImU32 text = canLearn ? IM_COL32(240, 200, 110, 255) : learned ? IM_COL32(150, 235, 165, 255) : IM_COL32(150, 155, 165, 255);
            dl->AddRectFilled(a, b, hovered ? IM_COL32(34, 44, 58, 255) : IM_COL32(22, 30, 40, 255), 6.0f);
            dl->AddRect(a, b, edge, 6.0f, 0, canLearn || learned ? 2.0f : 1.0f);
            const float pad = fh * 0.4f, small = fh * 0.8f;
            dl->PushClipRect(ImVec2(a.x + 2, a.y + 2), ImVec2(b.x - 2, b.y - 2), true); // every line stays in its node
            dl->AddText(ImVec2(a.x + pad, a.y + pad), text, p.name.c_str());
            // the effect: at most two lines (the rest in the tooltip), above the status row
            dl->PushClipRect(ImVec2(a.x, a.y), ImVec2(b.x, b.y - pad - small * 1.15f), true);
            dl->AddText(ImGui::GetFont(), fh * 0.85f, ImVec2(a.x + pad, a.y + pad + fh * 1.25f), IM_COL32(190, 196, 206, 255), p.text.c_str(), nullptr, nodeW - 2.0f * pad);
            dl->PopClipRect();
            // the status row: what you can do (or why not: just "locked" under a locked tier, whose header says so) and the ranks
            const char* status = full ? "learned" : canLearn ? "can learn" : level < need ? "locked" : why.c_str(); // (click the node: the tooltip says so)
            const char* ranks = TextFormat("rank %d/%d", rank, p.ranks);
            const float rw = ImGui::GetFont()->CalcTextSizeA(small, FLT_MAX, 0.0f, ranks).x;
            dl->PushClipRect(ImVec2(a.x, a.y), ImVec2(b.x - pad - rw - pad, b.y), true);
            dl->AddText(ImGui::GetFont(), small, ImVec2(a.x + pad, b.y - pad - small), canLearn ? IM_COL32(230, 180, 60, 255) : IM_COL32(130, 136, 146, 255), status);
            dl->PopClipRect();
            dl->AddText(ImGui::GetFont(), small, ImVec2(b.x - pad - rw, b.y - pad - small), learned ? IM_COL32(150, 235, 165, 255) : IM_COL32(130, 136, 146, 255), ranks);
            dl->PopClipRect();
            if (hovered) {
                ImGui::BeginTooltip();
                ImGui::TextUnformatted(p.name.c_str());
                ImGui::TextDisabled("%s", p.text.c_str());
                ImGui::TextDisabled("%d point%s a rank, %d rank%s", p.cost, p.cost == 1 ? "" : "s", p.ranks, p.ranks == 1 ? "" : "s");
                if (canLearn) ImGui::TextColored(ImVec4(0.9f, 0.7f, 0.25f, 1.0f), "Click to learn a rank");
                if (!p.needs.empty()) {
                    std::string n;
                    for (const std::string& q : p.needs) { const PerkDef* qd = Skills_FindPerk(q); n += (n.empty() ? "" : ", ") + (qd ? qd->name : q); }
                    ImGui::TextDisabled("Needs %s", n.c_str());
                }
                if (!canLearn && !full) ImGui::TextColored(ImVec4(0.9f, 0.6f, 0.4f, 1.0f), "%s", why.c_str());
                ImGui::EndTooltip();
            }
            if (clicked && canLearn) { if (gSkillSpendRequest) gSkillSpendRequest(p.id); else Skills_Buy(st, p.id); } // (#920: joined, the server decides)
            ImGui::PopID();
            col = (col + 1) % kNodeCols;
        }
        ImGui::Spacing();
    }
    if (!any) ImGui::TextDisabled("%s's perks come in the next update; its XP and levels count already.", d.name.c_str());
    ImGui::Separator();
    const int spent = Skills_PointsSpent(st, d.id);
    const int manuals = (int)inventory.count("retraining_manual");
    if (spent > 0) {
        if (manuals > 0) {
            if (ImGui::Button(TextFormat("Retrain %s (uses a retraining manual: %d points back)", d.name.c_str(), spent))) {
                inventory.remove("retraining_manual", 1);
                Skills_Respec(st, d.id);
            }
        } else ImGui::TextDisabled("A retraining manual (crafted: cloth, charcoal, planks) gives back every point of this tree.");
    }
    ImGui::EndChild();
}

// (#828) The Team screen's helpers: how long ago, a compass word, where a teammate is from here.
namespace {
    std::string AgoText(uint32_t sec) {
        if (sec == 0) return "";
        if (sec < 120) return "just now";
        if (sec < 2 * 3600) return TextFormat("%u min ago", sec / 60);
        if (sec < 2 * 86400) return TextFormat("%u h ago", sec / 3600);
        return TextFormat("%u days ago", sec / 86400);
    }
    const char* CompassWord(float deg) {
        static const char* const w[8] = { "N", "NE", "E", "SE", "S", "SW", "W", "NW" };
        return deg < 0.0f ? "" : w[(int)std::floor(fmodf(deg + 22.5f, 360.0f) / 45.0f) % 8];
    }
    // "420 m NE", "on luna", "here"; "" when their position hasn't arrived
    std::string WhereText(const NetClient& c, uint32_t session) {
        uint64_t body = 0; DVec3 local{};
        if (!session || !NetClient_OtherView(c, session, body, local)) return "";
        const Planetoid* theirs = Planetoid::findNet(body);
        Planetoid* mine = player.nearest;
        if (!theirs) return "";
        if (theirs != mine || !mine) return "on " + theirs->name;
        const double d = (local - mine->toLocalD(player.position)).length() * mine->scale;
        if (d < 8.0) return "here";
        const float bearing = uiWorld ? Nav_Bearing(*uiWorld, local.toFloat()) : -1.0f;
        return d < 1000.0 ? TextFormat("%.0f m %s", d, CompassWord(bearing)) : TextFormat("%.1f km %s", d / 1000.0, CompassWord(bearing));
    }
    // The right pane: what the team shares - its tier, the recipes you know, its claims (of its budget) and who may enter them.
    void DrawTeamShared(int claimBudget, const std::vector<std::string>& allies) {
        const float fh = ImGui::GetFontSize();
        ImGui::TextColored(ImVec4(0.6f, 0.85f, 1.0f, 1.0f), "Shared");
        ImGui::Separator();
        const int cur = QuestBook_Tier(), tiersN = QuestBook_TierCount();
        if (cur < tiersN) ImGui::Text("Quests: tier %d of %d (%s)", cur, tiersN - 1, QuestBook_TierName(cur));
        else ImGui::Text("Quests: every tier done");
        ImGui::Text("Recipes you know: %d", (int)LocalSession().recipes.known.size());
        ImGui::Dummy(ImVec2(0, fh * 0.4f));
        // Claims: a joined game's from the server's list, single player's from the claims here
        struct Row { std::string body; Vector3 centre; float radius; };
        std::vector<Row> rows;
        if (gNetClient && gNetClient->state == NetClient::State::Joined) { for (const ClaimView& cv : Claims_Mirror()) if (cv.ours) rows.push_back({ cv.body, cv.centre, cv.radius }); }
        else for (const Claim& cl : Claims_OfTeam(Teams_Of(PlayerId_Local()))) rows.push_back({ cl.body, cl.centre, cl.radius });
        ImGui::Text("Claims: %d of %d", (int)rows.size(), claimBudget);
        if (rows.empty()) ImGui::TextDisabled("  none yet: place a claim banner to claim land");
        for (const Row& r : rows) {
            std::string where;
            if (player.nearest && player.nearest->name == r.body) {
                const double d = (DVec3(r.centre) - player.nearest->toLocalD(player.position)).length() * player.nearest->scale;
                where = d < r.radius ? "you are in it" : (d < 1000.0 ? TextFormat("%.0f m %s", d, CompassWord(uiWorld ? Nav_Bearing(*uiWorld, r.centre) : -1.0f))
                                                                     : TextFormat("%.1f km %s", d / 1000.0, CompassWord(uiWorld ? Nav_Bearing(*uiWorld, r.centre) : -1.0f)));
            } else where = "on " + r.body;
            ImGui::BulletText("%.0f m across, %s", r.radius * 2.0f, where.c_str());
        }
        ImGui::Dummy(ImVec2(0, fh * 0.4f));
        ImGui::TextColored(ImVec4(0.6f, 0.85f, 1.0f, 1.0f), "In our claims");
        ImGui::Separator();
        ImGui::TextWrapped("Our team%s may build, dig, open chests and machines, farm and use beds; anyone else is refused.", allies.empty() ? "" : " and our allies");
        if (!allies.empty()) {
            std::string a;
            for (const std::string& n : allies) a += (a.empty() ? "" : ", ") + n;
            ImGui::TextWrapped("Allies: %s", a.c_str());
        } else ImGui::TextDisabled("No allies.");
    }
}

// The Team screen (#604; the two panes of #828). Joined to a server: the team as the server says (members, who's online,
// where they are, their health, the owner; when the others were last seen), with Invite / Leave / Kick / Rename / team chat
// as /team commands (the server decides), and an invite waiting to accept. Single player: a team of one, with its claims.
static void DrawTeam() {
    const float s = GameUI_Scale();
    ImGui::Dummy(ImVec2(0, 12 * s));
    NetClient* c = gNetClient && gNetClient->state == NetClient::State::Joined ? gNetClient : nullptr;
    const float sharedW = std::max(ImGui::GetContentRegionAvail().x * 0.38f, 320.0f * s);
    if (!c) {
        const int id = Teams_Of(PlayerId_Local());
        const TeamRecord* team = nullptr;
        for (const TeamRecord& t : Teams_All()) if (t.id == id) team = &t;
        ImGui::BeginChild("teamleft", ImVec2(ImGui::GetContentRegionAvail().x - sharedW - 8.0f * s, 0.0f), true);
        ImGui::TextColored(ImVec4(0.88f, 0.66f, 0.25f, 1.0f), "%s", team ? team->name.c_str() : "Your team");
        ImGui::SameLine(); ImGui::TextDisabled("1 member");
        ImGui::Separator();
        ImGui::Text("%s (you)  - owner", settings.mpName.empty() ? "You" : settings.mpName.c_str());
        ImGui::Dummy(ImVec2(0, 8 * s));
        ImGui::TextWrapped("You are a team of one. Quests, tiers and claims belong to the team. In multiplayer, invite others to share them: Multiplayer on the main menu.");
        ImGui::EndChild();
        ImGui::SameLine();
        ImGui::BeginChild("teamshared", ImVec2(0.0f, 0.0f), true);
        std::vector<std::string> allies;
        for (const int a : Claims_AlliesOf(id)) if (const TeamRecord* at = Teams_Find(a)) allies.push_back(at->name);
        DrawTeamShared(Claims_Budget(std::max(1, (int)QuestBook_LocalState().visitedBodies.size())), allies);
        ImGui::EndChild();
        return;
    }
    const MsgTeamState& t = c->team;
    ImGui::BeginChild("teamleft", ImVec2(ImGui::GetContentRegionAvail().x - sharedW - 8.0f * s, 0.0f), true);
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
            ImGui::TableNextColumn(); // (#828) where they are from here, or when they were last seen
            if (m.online && !m.you) { const std::string w = WhereText(*c, m.session); if (!w.empty()) ImGui::TextUnformatted(w.c_str()); }
            else if (!m.online && m.lastSeen) ImGui::TextDisabled("last seen %s", AgoText(m.lastSeen).c_str());
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
    // (#828) Land claims (#343): the team's own, the other teams' it has heard of (and whether they let you in), and for the
    // owner, which teams may build on yours (/claim ally, /claim unally: the server keeps the list; /claims prints it).
    ImGui::Separator();
    ImGui::TextColored(ImVec4(0.88f, 0.66f, 0.25f, 1.0f), "Land claims");
    int ours = 0, theirs = 0;
    for (const ClaimView& cv : Claims_Mirror()) (cv.ours ? ours : theirs) += 1;
    if (ours == 0 && theirs == 0) ImGui::TextDisabled("None known."); // (how to claim: the Shared panel's Claims line)
    if ((ours > 0 || theirs > 0) && ImGui::BeginTable("claims", 4, ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        for (int pass = 0; pass < 2; ++pass) // ours first
            for (const ClaimView& cv : Claims_Mirror()) {
                if (cv.ours != (pass == 0)) continue;
                ImGui::TableNextRow();
                ImGui::TableNextColumn();
                const ImVec4 col((cv.colour >> 16 & 255) / 255.0f, (cv.colour >> 8 & 255) / 255.0f, (cv.colour & 255) / 255.0f, 1.0f);
                ImGui::TextColored(col, "%s", cv.ours ? "Your team" : cv.teamName.c_str());
                ImGui::TableNextColumn();
                ImGui::TextUnformatted(cv.body.c_str());
                ImGui::TableNextColumn();
                ImGui::Text("%.0f m across", 2.0f * cv.radius);
                ImGui::TableNextColumn();
                if (!cv.ours) { if (cv.allowed) ImGui::TextColored(ImVec4(0.55f, 0.9f, 0.55f, 1.0f), "lets you in"); else ImGui::TextDisabled("closed to you"); }
            }
        ImGui::EndTable();
    }
    if (youOwn) {
        static char ally[32] = {};
        ImGui::SetNextItemWidth(-1.0f); // the column's width: the two buttons go on the next line (side by side they were cut off)
        ImGui::InputTextWithHint("##ally", "another team's name", ally, sizeof(ally));
        if (ImGui::Button("Let them in") && ally[0]) { say(std::string("/claim ally ") + ally); ally[0] = 0; }
        ImGui::SameLine();
        if (ImGui::Button("Shut them out") && ally[0]) { say(std::string("/claim unally ") + ally); ally[0] = 0; }
    }
    if (ImGui::SmallButton("List claims and allies in chat")) say("/claims");
    ImGui::Separator();
    ImGui::SetNextItemWidth(420 * s);
    if (ImGui::InputTextWithHint("##teamsay", "say to your team (Enter)", teamSay, sizeof(teamSay), ImGuiInputTextFlags_EnterReturnsTrue) && teamSay[0]) {
        say(std::string("/t ") + teamSay); teamSay[0] = 0;
    }
    ImGui::TextDisabled("Also from the command line: /team, /team invite NAME, /team accept, /team leave, /t TEXT");
    ImGui::EndChild();
    ImGui::SameLine();
    ImGui::BeginChild("teamshared", ImVec2(0.0f, 0.0f), true);
    DrawTeamShared(t.claimBudget, t.allies);
    ImGui::EndChild();
}

// ---- The inventory screen's player (#822): a turning preview, the suit, its module slots and the off hand ----
namespace {
    RenderTexture2D preview{};
    float previewTurn = 0.0f;   // radians: the preview turns slowly; dragging over it turns it by hand
    bool previewDrawn = false;  // rendered this frame (GameUI_PrepareFrame)
    ImVec2 previewWant{ 250.0f, 360.0f }; // (#887) the size its panel gives it (px), from the last frame
    Color SuitColour() { const uint32_t c = settings.mpSuit; return c ? Color{ (unsigned char)(c >> 16), (unsigned char)(c >> 8), (unsigned char)c, 255 } : Color{ 0, 0, 0, 0 }; }
}

bool GameUI_PreviewDrawn() { return previewDrawn; }
void GameUI_PrepareFrame(Scene& world) {
    previewDrawn = false;
    if (!inventoryOpen || screen != GameScreen::Inventory || container.open) return;
    const Vector3 up = Vector3Normalize(player.up);
    if (Vector3Length(player.up) < 0.5f) return; // (not placed yet)
    const float s = GameUI_Scale();
    const int w = std::clamp((int)previewWant.x, (int)(120 * s), 2048), h = std::clamp((int)previewWant.y, (int)(160 * s), 2048); // (#887: its panel's size)
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

// The turning player preview (#822). An InvisibleButton takes the mouse and the image is drawn over it (#892: as a plain
// ImGui::Image it took no input, so IsItemActive was never true and dragging never turned the player). Dragging turns it
// 0.012 rad a pixel.
static void PreviewWidget(ImTextureID tex, ImVec2 size, float& turn) {
    const ImVec2 p = ImGui::GetCursorScreenPos();
    ImGui::InvisibleButton("##preview", size);
    ImGui::GetWindowDrawList()->AddImage(tex, p, ImVec2(p.x + size.x, p.y + size.y), ImVec2(0, 1), ImVec2(1, 0));
    if (ImGui::IsItemActive()) turn -= ImGui::GetIO().MouseDelta.x * 0.012f; // (#1073) the side you drag to is the side it turns to
    if (ImGui::IsItemHovered() && !ImGui::IsItemActive()) ImGui::SetTooltip("Drag to turn");
}

// The player panel (#822; #887: sized to its area, the top right of the Inventory screen): the suit and the off hand on the
// left, the turning preview filling the rest, the module slots under it, and what else is carried.
// (#913 review) The slots down the left (suit and off hand, modules under them; the off hand's help in its tooltip, spare suits
// and modules behind a button), and the preview, the panel's feature, in all the rest at its full height. Fits 1080 p.
static void DrawYou() {
    ImGui::TextUnformatted("You");
    ImGui::Separator();
    const float avW = ImGui::GetContentRegionAvail().x, avH = ImGui::GetContentRegionAvail().y;
    ImGui::PushFont(nullptr, ImGui::GetStyle().FontSizeBase * 0.75f); // (the column's labels small: at 150% they left the slots no room)
    const int slots = (int)Suit_Modules().size();
    const int modCols = slots > 1 ? 2 : 1;
    // (#1148) The slots exactly the items grid's size (suit and off hand side by side, the modules two to a row under them);
    // a short panel scrolls rather than shrinking them (a height fit had made them half the grid's).
    const float big = gGridCell > 0.0f ? gGridCell : 54.0f * GameUI_Scale();
    const float cell = big;
    ImGui::BeginGroup();
    ImGui::TextDisabled("Suit");
    SuitSlot(big); // (#1161: the off hand moved beside the bin, under the hotbar)
    int fitted = 0; for (const std::string& m : Suit_Modules()) fitted += m.empty() ? 0 : 1;
    ImGui::TextDisabled("Modules %d/%d", fitted, slots);
    for (int i = 0; i < slots; ++i) { if (i % modCols) ImGui::SameLine(0.0f, 4.0f); ModuleSlot(i, cell); }
    if (ImGui::SmallButton("Spare parts")) ImGui::OpenPopup("##spares");
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Spare suits and modules you carry: wear or fit them");
    if (ImGui::BeginPopup("##spares")) { Suit_DrawCarried(); ImGui::EndPopup(); }
    ImGui::EndGroup();
    ImGui::PopFont();
    const float leftW = ImGui::GetItemRectSize().x;
    ImGui::SameLine(0.0f, 12.0f);
    // The player, turning; drag over it to turn it: the rest of the panel, full height (the size is the render target's
    // next frame: GameUI_PrepareFrame).
    previewWant = ImVec2(std::max(120.0f, avW - leftW - 16.0f), std::max(160.0f, avH - 4.0f));
    if (previewDrawn) PreviewWidget((ImTextureID)(intptr_t)preview.texture.id, ImVec2((float)preview.texture.width, (float)preview.texture.height), previewTurn);
    else ImGui::Dummy(previewWant);
}

// The stats panel (#887, the bottom right of the Inventory screen): what changes the player's numbers. The suit's
// multipliers, the skill perks' (#883; they stack with the suit's), and room for team bonuses and timed effects.
static const char* PerkStatName(const std::string& stat) {
    static const struct { const char* id; const char* name; } kNames[] = {
        { "rock_speed", "Rock strokes" }, { "ore_speed", "Ore strokes" }, { "rock_radius", "Rock stroke reach" }, { "ore_yield", "Ore per stroke" },
        { "stone_yield", "Stone of the rock's kind" }, { "gem_chance", "Gem odds" }, { "pick_wear", "Pick wear" },
        // (#919) the other trees
        { "soil_speed", "Soil strokes" }, { "soil_radius", "Soil stroke reach" }, { "shovel_wear", "Shovel wear" },
        { "chop_speed", "Chops" }, { "log_yield", "Logs per tree" }, { "axe_wear", "Axe wear" },
        { "fibre_yield", "Fibre per bush" }, { "stone_find", "Loose stone" }, { "crop_yield", "Crop per harvest" },
        { "seed_yield", "Seeds per harvest" }, { "break_speed", "Taking parts down" }, { "walk_speed", "Walk" },
        { "fall_safe", "Safe landing speed" }, { "hunger", "Hunger" }, { "cold", "Cold felt" } };
    for (const auto& n : kNames) if (stat == n.id) return n.name;
    return stat.c_str();
}
static void DrawStats() {
    ImGui::TextUnformatted("Stats");
    ImGui::Separator();
    // (#913 review) One stat a cell, two to a row, a size smaller than the window's text: it fits without scrolling.
    ImGui::PushFont(nullptr, ImGui::GetStyle().FontSizeBase * 0.8f);
    ImGui::PushTextWrapPos(0.0f);
    ImGui::TextColored(ImVec4(0.6f, 0.85f, 1.0f, 1.0f), "Suit");
    const SuitStats& su = Suit_Stats();
    const std::pair<const char*, std::string> rows[] = {
        { "Oxygen", TextFormat("x%.2f", su.oxygen) }, { "Power", TextFormat("x%.2f", su.power) },
        { "Jet fuel", TextFormat("x%.2f", su.jetFuel) }, { "Thrust", TextFormat("x%.2f", su.jetThrust) },
        { "Walk", TextFormat("x%.2f", su.walk) }, { "Lamp", TextFormat("x%.2f", su.lamp) },
        { "Hunger", TextFormat("x%.2f", su.food) }, { "Underwater O2", TextFormat("x%.1f", su.underwater) },
        { "Cold felt", TextFormat("%.0f%%", su.cold * 100.0f) }, { "Heat", su.heatProof ? "regulated" : "hurts > 35 C" },
        { "Grip", su.grip ? "on" : "off" }, { "Ore scanner", su.scanner ? "on" : "off" } };
    // (#887, #1148) Labels in a column as wide as the longest one, never wrapped (stretched columns squeezed them to a letter
    // a line in a narrow panel); two stats to a row only where both pairs fit, else one.
    const int nRows = (int)(sizeof(rows) / sizeof(rows[0]));
    float labelW = 0.0f, valueW = 0.0f;
    for (const auto& r : rows) { labelW = std::max(labelW, ImGui::CalcTextSize(r.first).x); valueW = std::max(valueW, ImGui::CalcTextSize(r.second.c_str()).x); }
    const float cellPad = ImGui::GetStyle().CellPadding.x * 2.0f;
    const int pairs = ImGui::GetContentRegionAvail().x >= 2.0f * (labelW + valueW + 2.0f * cellPad) + 8.0f ? 2 : 1;
    ImGui::PushTextWrapPos(-1.0f); // (no wrapping inside the table)
    if (ImGui::BeginTable("suitstats", pairs * 2, ImGuiTableFlags_SizingFixedFit)) {
        for (int c = 0; c < pairs; ++c) {
            ImGui::TableSetupColumn("label", ImGuiTableColumnFlags_WidthFixed, labelW);
            ImGui::TableSetupColumn("value", ImGuiTableColumnFlags_WidthFixed, valueW + 8.0f);
        }
        for (int k = 0; k < nRows; ++k) {
            if (k % pairs == 0) ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextDisabled("%s", rows[k].first);
            ImGui::TableNextColumn(); ImGui::TextUnformatted(rows[k].second.c_str());
        }
        ImGui::EndTable();
    }
    ImGui::PopTextWrapPos();
    ImGui::Spacing();
    ImGui::TextColored(ImVec4(0.6f, 0.85f, 1.0f, 1.0f), "Skills");
    const SkillState& st = LocalSession().skills;
    std::vector<std::string> stats;
    for (const PerkDef& p : Skills_Perks())
        for (const auto& [stat, per] : p.effects) if (std::find(stats.begin(), stats.end(), stat) == stats.end()) stats.push_back(stat);
    int shown = 0;
    for (const std::string& stat : stats) {
        const double m = Skills_Stat(st, stat);
        if (std::fabs(m - 1.0) < 1e-6) continue;
        std::string from; // the perks that give it
        for (const auto& [id, r] : st.ranks)
            if (const PerkDef* p = Skills_FindPerk(id); p && r > 0)
                for (const auto& [s2, per] : p->effects) if (s2 == stat) from += (from.empty() ? "" : ", ") + p->name + (p->ranks > 1 ? TextFormat(" %d", r) : "");
        ImGui::Text("%s x%.2f", PerkStatName(stat), m);
        ImGui::SameLine(); ImGui::TextDisabled("(%s)", from.c_str());
        ++shown;
    }
    if (!shown) ImGui::TextDisabled("No perks yet (the Skills tab)");
    ImGui::Spacing();
    ImGui::TextColored(ImVec4(0.6f, 0.85f, 1.0f, 1.0f), "Team and effects");
    ImGui::TextDisabled("None yet");
    ImGui::PopTextWrapPos();
    ImGui::PopFont();
}

// The items panel of the Inventory and Crafting screens (#887, #889): the grid and hotbar at the shared slot size (#1148:
// SlotCellFor; it was sized to fill 72% of the window, larger than the Catalog page's cells).
static void DrawItemsPanel() {
    // The panel is as wide as its grid (#891); the other side gets what is left.
    const ImVec2 a = ImGui::GetContentRegionAvail();
    const float gap = 4.0f, pad = ImGui::GetStyle().WindowPadding.x * 2.0f;
    const float cell = gGridCell; // (#1148) the one slot size of every item page (SlotCellFor, set as the screen opens)
    gDrawnCellInventory = cell;   // (#887 test)
    const float panelW = std::max(Inventory::kCols * cell + (Inventory::kCols - 1) * gap + pad, ImGui::CalcTextSize("Sort (middle-click)  Deposit (Q)").x + 2.0f * cell + gap + pad); // (#1161: the off hand, the bin and the buttons)
    ImGui::BeginChild("invgrid", ImVec2(std::min(panelW, a.x * 0.72f), 0.0f), false, ImGuiWindowFlags_NoScrollbar | ImGuiWindowFlags_NoScrollWithMouse);
    ImGui::TextUnformatted("Items");
    DrawGrid(cell);
    ImGui::EndChild();
}

// The Crafting screen (#889, the user's layout): the items panel as on the Inventory screen, and the crafting column beside
// it (search and chips, the recipes as large icons, the selected recipe's page and the queue under them).
static void DrawCraftingScreen() {
    DrawItemsPanel();
    ImGui::SameLine();
    ImGui::BeginChild("craftcol", ImVec2(0.0f, 0.0f), false, ImGuiWindowFlags_NoScrollbar);
    DrawCrafting();
    ImGui::EndChild();
}

// The Inventory screen with no container open (#887, the user's layout): the items panel on the left, the player at the
// top right (half the height), the stats under it.
static void DrawInventoryScreen() {
    DrawItemsPanel();
    ImGui::SameLine();
    ImGui::BeginChild("invright", ImVec2(0.0f, 0.0f), false, ImGuiWindowFlags_NoScrollbar);
    const float h = ImGui::GetContentRegionAvail().y;
    ImGui::BeginChild("you", ImVec2(0.0f, h * 0.50f), true); // (#891: scrolls when the large text doesn't fit)
    DrawYou();
    ImGui::EndChild();
    ImGui::BeginChild("stats", ImVec2(0.0f, 0.0f), true);
    DrawStats();
    ImGui::EndChild();
    ImGui::EndChild();
}

// (#894) The pointer is over none of our windows. IsWindowHovered(AnyWindow) alone says no while an item is active, and a
// cell's InvisibleButton turns active on the very press that picks a stack up: that press counted as "outside the window"
// and dropped the stack it had just picked up. Blocked by an active item or a pop-up, a window under the pointer still counts.
static bool PointerOutsideWindows() {
    return !ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow | ImGuiHoveredFlags_AllowWhenBlockedByActiveItem | ImGuiHoveredFlags_AllowWhenBlockedByPopup);
}

static bool tabWindowJustOpened = true, tabWindowUnsaved = false; // (#1080) place it from settings on opening; a move or resize to save
void GameUI_Draw(Scene& world) {
    uiWorld = &world;
    RecipeBook_Update(); // new recipes: a toast
    RunCraftQueue(); // (#824) queued crafts, window open or not
    if (!inventoryOpen) { SlotRules_Cancel(); ReturnCursorStack(); binned = InvSlot{}; container.open = false; tabWindowJustOpened = true; return; } // a stack still held goes back in (#192); (#1080) placed again when it opens
    if (cursorStack.empty()) cursorCredited = false; // (#1728)
    if (cursorStack.heat > 0.0f) cursorStack.heat = std::max(0.0f, cursorStack.heat - GameFrameTime() / Inventory::kCoolSeconds); // held, it cools as in a cell (#580)
    CheckContainer(world);
    hoverGrid = nullptr; hoverSlot = -1; // (set again as the grids are drawn: #862)
    if (container.open && screen != GameScreen::Inventory) screen = GameScreen::Inventory; // a chest or machine opens beside your items
    const ImGuiIO& io = ImGui::GetIO();
    const float s = GameUI_Scale();
    // (#1080) No title bar (the tabs are the window's header). It opens centred at 70% of the screen by default, moves (drag an
    // empty part of it) and resizes (its corner) like the database window, and keeps its place and size in settings as
    // fractions of the screen, so another resolution or interface size keeps the layout.
    const ImVec2 D = io.DisplaySize;
    static ImVec2 lastScreen{ 0.0f, 0.0f };
    const bool place = tabWindowJustOpened || D.x != lastScreen.x || D.y != lastScreen.y;
    tabWindowJustOpened = false; lastScreen = D;
    if (place) {
        float* f = settings.tabWin;
        ImGui::SetNextWindowPos(ImVec2(f[0] * D.x, f[1] * D.y), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(f[2] * D.x, f[3] * D.y), ImGuiCond_Always);
    }
    ImGui::SetNextWindowSizeConstraints(ImVec2(std::min(D.x, 0.45f * D.x + 120.0f * s), std::min(D.y, 0.45f * D.y + 80.0f * s)), D); // (the panels need room)
    ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.07f, 0.10f, 0.14f, 0.7f)); // (#1079) the world shows through at 30%; what's in it stays opaque
    const ImGuiWindowFlags flags = ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoScrollbar;
    if (ImGui::Begin("Game window", nullptr, flags)) {
        { // (#1080) moved or resized: kept (fractions of the screen), saved once the drag ends
            const ImVec2 p = ImGui::GetWindowPos(), z = ImGui::GetWindowSize();
            const float now[4] = { p.x / D.x, p.y / D.y, z.x / D.x, z.y / D.y };
            bool moved = false;
            for (int k = 0; k < 4; ++k) moved = moved || fabsf(now[k] - settings.tabWin[k]) > 0.002f;
            if (moved) { for (int k = 0; k < 4; ++k) settings.tabWin[k] = std::clamp(now[k], 0.0f, 1.0f); tabWindowUnsaved = true; }
            if (tabWindowUnsaved && !ImGui::IsMouseDown(ImGuiMouseButton_Left)) { tabWindowUnsaved = false; Settings_Save(); }
        }
        // The tabs (#810's order; Catalog beside Inventory in creative only, #823), their keys from the bindings; at the right
        // end the temporary Help page (until #826) and Close.
        struct Tab { GameScreen s; const char* name; Act key; };
        static const Tab tabs[] = { { GameScreen::Inventory, "Inventory", Act::Inventory }, { GameScreen::Catalog, "Catalog", Act::Count }, { GameScreen::Crafting, "Crafting", Act::Crafting },
                                    { GameScreen::Quests, "Quests", Act::Count }, { GameScreen::Skills, "Skills", Act::Skills },
                                    { GameScreen::Research, "Research", Act::Count }, { GameScreen::Tome, "Tome", Act::Count }, { GameScreen::Team, "Team", Act::Count }, { GameScreen::Statistics, "Statistics", Act::Count } }; // (#1192: Research beside Skills; #1199 Statistics last, beside the Database; #1251 the Tome beside Research, once magic is found)
        auto tabShown = [](GameScreen s) { return (s != GameScreen::Catalog || terrainTool.creative) && (s != GameScreen::Tome || TomeTab_Shown()); };
        // (#891) tabs as large as the row allows at the doubled scale (6 at 190 x s overflowed 1080 p): the row less Database and Close
        const float smallW = 92.0f * s, tabH = 34.0f * s; // (#913 review: "Close [Esc]" keeps its key; #1206: 102 -> 92, room for 8 tabs)
        int nTabs = 0;
        for (const Tab& t : tabs) nTabs += tabShown(t.s);
        const float tabW = std::min(145.0f * s, (ImGui::GetContentRegionAvail().x - 2 * smallW - 6.0f * s * (nTabs + 2)) / (float)nTabs);
        if (!tabShown(screen)) screen = GameScreen::Inventory; // (left creative; the Tome without magic)
        {   // (#1206) one font for the row: the largest at which every shown name fits
            const char* names[16]; int n = 0;
            for (const Tab& t : tabs) if (n < 16 && tabShown(t.s)) names[n++] = t.name;
            ScreenTabFont(names, n, tabW);
        }
        for (const Tab& t : tabs) {
            if (!tabShown(t.s)) continue;
            const char* key = t.key == Act::Count ? "" : PlayerInput_KeyName(PlayerInput_Binding(t.key));
            if (ScreenTab(t.name, key, screen == t.s, tabW, tabH, true)) { screen = t.s; if (t.s != GameScreen::Inventory) container.open = false; }
            ImGui::SameLine(0.0f, 6.0f * s);
        }
        const float right = ImGui::GetWindowContentRegionMax().x;
        ImGui::SetCursorPosX(right - 2 * smallW - 6.0f * s);
        if (ScreenTab("Database", "", databaseOpen, smallW, tabH, true)) databaseOpen = !databaseOpen; // (#826: alt-click an item for its page)
        ImGui::SameLine(0.0f, 6.0f * s);
        if (ScreenTab("Close", "Esc", false, smallW, tabH, true)) inventoryOpen = false;
        ImGui::Separator();
        ImGui::BeginChild("screen", ImVec2(0, 0), false);
        gGridCell = SlotCellFor(ImGui::GetContentRegionAvail()); // (#1148) one slot size for the Inventory, Crafting and Catalog pages
        // Crafting, the Catalog and an open container show your items on the left; the Inventory screen alone lays itself out
        // (#887); the others use the whole width.
        const bool ownLayout = (screen == GameScreen::Inventory && !container.open) || screen == GameScreen::Crafting; // (#887, #889)
        const bool items = !ownLayout && (screen == GameScreen::Inventory || screen == GameScreen::Crafting || screen == GameScreen::Catalog);
        if (items) {
            ImGui::Columns(2, "inv", true);
            static float colW = -1.0f;
            { const float w = ImGui::GetWindowContentRegionMax().x * 0.5f; if (fabsf(colW - w) > 1.0f) { ImGui::SetColumnWidth(0, w); colW = w; } } // (#891: half, at any scale)
            ImGui::TextUnformatted("Items");
            ImGui::Separator();
            DrawGrid(gGridCell); // (#1148: the shared slot size: this column was its model)
            ImGui::NextColumn();
        }
        switch (screen) {
            case GameScreen::Catalog: DrawCreative(); break; // (#823) every item, free: creative only
            case GameScreen::Inventory:
                if (container.open) DrawContainer();
                else DrawInventoryScreen(); // (#887) the grid scaled up, the player (#822) and the stats
                break;
            case GameScreen::Crafting: DrawCraftingScreen(); break; // (#889) the items panel and the crafting column
            case GameScreen::Quests:
                DrawTeamProgress(); // (#825)
                ImGui::Separator();
                QuestTree_DrawTab(); // (#890: the tree is the Quests screen; a tier's objectives and pins open in a pop-out)
                break;
            case GameScreen::Skills: DrawSkills(); break; // (#883)
            case GameScreen::Research: ResearchTab_Draw(); break; // (#1192)
            case GameScreen::Tome: TomeTab_Draw(); break; // (#1251) magic, once found
            case GameScreen::Statistics: StatsTab_Draw(); break; // (#1199)
            case GameScreen::Team: DrawTeam(); break;
            case GameScreen::Count: break;
        }
        if (items) ImGui::Columns(1);
        ImGui::EndChild();
    }
    ImGui::End();
    ImGui::PopStyleColor();
    // (#889 follow-up) The HUD is hidden while the window is open, but the window's own messages (the bin, Deposit, crafting)
    // still show: over the window's bottom edge, on a dark backing.
    if (GetTime() < toastUntil) {
        ImDrawList* fg = ImGui::GetForegroundDrawList();
        const float fs = ImGui::GetFontSize();
        const ImVec2 sz = ImGui::GetFont()->CalcTextSizeA(fs, FLT_MAX, 0.0f, toast.c_str());
        const float margin = 24.0f * s; // (the screen's bottom margin, as the window had before #1080)
        const ImVec2 at((io.DisplaySize.x - sz.x) * 0.5f, io.DisplaySize.y - margin - sz.y * 0.5f);
        fg->AddRectFilled(ImVec2(at.x - 10, at.y - 4), ImVec2(at.x + sz.x + 10, at.y + sz.y + 4), IM_COL32(10, 14, 20, 235), 4.0f);
        fg->AddText(ImGui::GetFont(), fs, at, IM_COL32(255, 230, 90, 255), toast.c_str());
    }
    // The database (#826): a pop-up over the game window, opened by alt-clicking an item (at its page) or the Database
    // button; browse on from there (links, search, Back). Like Factorio's Factoriopedia.
    if (databaseOpen) {
        ImGui::SetNextWindowSize(ImVec2(io.DisplaySize.x * 0.62f, io.DisplaySize.y * 0.7f), ImGuiCond_Appearing);
        ImGui::SetNextWindowPos(ImVec2(io.DisplaySize.x * 0.5f, io.DisplaySize.y * 0.52f), ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
        ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0.08f, 0.11f, 0.16f, 1.0f)); // opaque: the screen under it doesn't show through
        if (ImGui::Begin("Database (alt-click any item)", &databaseOpen, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoSavedSettings)) {
            HelpManual_DrawTab();
        }
        ImGui::End();
        ImGui::PopStyleColor();
    }
    // The item rules (#862), once a frame: the mouse and keys with the slot under the mouse (Alt-clicks open the database).
    const bool heldBefore = !cursorStack.empty(); // (#894: only a stack already held when the click came may be dropped by it)
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
    // A click outside every window drops the held stack (right-click: one unit) at the feet (#564). Only a stack held before
    // this click, and only outside every window even while a cell's button is active (#894: a pickup click dropped it).
    if (heldBefore && !cursorStack.empty() && PointerOutsideWindows()) {
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
        const float m = GameUI_HudScale(); // (#925) at the HUD's size, above the hotbar
        const int w = GetScreenWidth(), h = GetScreenHeight(), line = (int)(26 * m), icon = (int)(22 * m), fs = (int)(20 * m);
        const size_t first = pickups.size() > kPickupLines ? pickups.size() - kPickupLines : 0;
        int y = h - (int)(140 * m) - (int)(pickups.size() - first) * line;
        for (size_t i = first; i < pickups.size(); ++i, y += line) {
            const Pickup& p = pickups[i];
            const float a = (float)std::clamp((kPickupShow - (now - p.t)) / kPickupFade, 0.0, 1.0);
            const ItemDef* def = itemDB.find(p.id);
            const double r = std::round(p.amount);
            const char* txt = fabs(p.amount - r) < 0.005 ? TextFormat("+%.0f %s", r, def ? def->name.c_str() : p.id.c_str())
                                                          : TextFormat("+%.2f %s", p.amount, def ? def->name.c_str() : p.id.c_str());
            const int tw = MeasureText(txt, fs), x = w - 30 - tw;
            DrawRectangle(x - icon - 14, y - 2, tw + icon + 22, line - 2, Fade(BLACK, 0.35f * a));
            ItemIcons_Draw(p.id, { (float)(x - icon - 8), (float)y, (float)icon, (float)icon }, Fade(WHITE, a));
            DrawText(txt, x, y + 2, fs, Fade(RAYWHITE, a));
        }
    }
}
// Items taken from a machine or chest slot onto the cursor (#192), a tick after the click: onto the held stack when it is
// empty or the same item, else into the inventory.
static void CursorSink(const std::string& id, int count, float heat, const std::string& data) { // heat: fresh out of a furnace stays hot on the cursor (#580)
    if (count <= 0) return;
    // (#1728) A take that goes straight onto the cursor never passes addAmount: count it here (a joined client's server already did)
    const bool count1728 = gNetRole != NetRole::Client;
    if (!data.empty()) { // (#1446) one item with its own record: onto an empty cursor, else into the inventory as it is
        if (cursorStack.empty()) {
            cursorStack = InvSlot{}; cursorStack.id = id; cursorStack.centi = 100; cursorStack.heat = heat; cursorStack.data = data;
            if (count1728) inventory.noteObtained(id, 1.0, data);
            cursorCredited = true;
            return;
        }
        ++inventory.muted; inventory.addInstance(id, data, heat); --inventory.muted;
        return;
    }
    if (cursorStack.empty()) {
        cursorStack.id = id; cursorStack.centi = (int64_t)count * 100; cursorStack.heat = heat;
        if (count1728) inventory.noteObtained(id, count);
        cursorCredited = true;
        return;
    }
    if (cursorStack.id == id && !cursorStack.hasData()) {
        const int64_t add = (int64_t)count * 100;
        cursorStack.heat = (float)((cursorStack.heat * cursorStack.centi + heat * add) / (double)(cursorStack.centi + add));
        cursorStack.centi += add;
        if (count1728) inventory.noteObtained(id, count);
        cursorCredited = true;
        return;
    }
    ++inventory.muted; inventory.addAmount(id, count, heat); --inventory.muted;
}
void GameUI_InitPickupFeed() {
    inventory.onAdded = OnPickup; Factory::SetCursorSink(CursorSink);
    // (#1728) a machine finishing a recipe counts for this game's own quest book, whether or not anyone takes the item (the local
    // player's world: single player, and a host's own player; a server's teams are credited by Host_Start's hook on top of this)
    gFactoryMade = [](uint64_t, uint32_t, const std::string& item, int count) { QuestBook_NoteMade(item, count); };
}

// (#604) Your teammates on the HUD: each one's name and health (offline ones greyed), on the left under the quests; only
// when joined to a server and in a team of two or more (non-teammates show no health).
static void DrawTeamStrip() {
    const NetClient* c = gNetClient && gNetClient->state == NetClient::State::Joined ? gNetClient : nullptr;
    if (!c || c->team.members.size() < 2) return;
    const float s = GameUI_HudScale(); // (#925) the HUD's size
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

// (#226) A tier reached: a banner across the top for 6 s (it slides in, then fades): the tier, its name and the first
// sentence of what it is about. The fanfare and the sparkles are the caller's (survival.cpp, where the tier advances).
namespace { int bannerTier = -1; double bannerAt = -1e9; }
void GameUI_TierBanner(int tier) { bannerTier = tier; bannerAt = GetTime(); }
static void DrawTierBanner() {
    const double age = GetTime() - bannerAt;
    if (bannerTier < 0 || age > 6.0) return;
    const float s = GameUI_HudScale();
    const float in = (float)std::min(1.0, age / 0.3), alpha = (float)std::clamp((6.0 - age) / 1.0, 0.0, 1.0);
    const int w = GetScreenWidth();
    const int fsBig = (int)(44 * s), fsName = (int)(30 * s), fsText = (int)(18 * s), pad = (int)(16 * s);
    std::string text = QuestBook_TierText(bannerTier);
    if (const size_t dot = text.find(". "); dot != std::string::npos) text = text.substr(0, dot + 1); // the first sentence
    const char* head = TextFormat("TIER %d REACHED", bannerTier);
    const char* name = QuestBook_TierName(bannerTier);
    int tw = std::max({ MeasureText(head, fsBig), MeasureText(name, fsName), MeasureText(text.c_str(), fsText) });
    tw = std::min(tw, w - 4 * pad);
    const int bw = tw + 2 * pad, bh = fsBig + fsName + fsText + 4 * pad;
    const int x = (w - bw) / 2, y = (int)(GetScreenHeight() * 0.16f - (1.0f - in) * 40.0f * s);
    DrawRectangle(x, y, bw, bh, Fade(Color{ 16, 20, 28, 255 }, 0.82f * alpha));
    DrawRectangleLinesEx(Rectangle{ (float)x, (float)y, (float)bw, (float)bh }, 2.0f * s, Fade(Color{ 236, 196, 92, 255 }, alpha));
    int ty = y + pad;
    DrawText(head, (w - MeasureText(head, fsBig)) / 2, ty, fsBig, Fade(Color{ 250, 214, 110, 255 }, alpha)); ty += fsBig + pad;
    DrawText(name, (w - MeasureText(name, fsName)) / 2, ty, fsName, Fade(RAYWHITE, alpha)); ty += fsName + pad;
    DrawText(text.c_str(), (w - std::min(MeasureText(text.c_str(), fsText), tw)) / 2, ty, fsText, Fade(Color{ 200, 206, 216, 255 }, alpha));
}

// (#1713) While the shield is up: a ring round the crosshair, gold with "PARRY" during the parry window (a melee hit or a bite
// landing now is parried: the attacker staggers), then grey while it only blocks. The window is the shield's (0.22-0.45 s).
static void DrawBlockRing() {
    if (!survival.blocking) return;
    const float s = GameUI_HudScale();
    const float cx = GetScreenWidth() * 0.5f, cy = GetScreenHeight() * 0.5f, r = 26.0f * s;
    const double since = SimTime() - survival.blockStartedAt;
    const bool parry = since >= 0.0 && since <= survival.shieldParryWindow;
    if (parry) {
        const float left = 1.0f - (float)(since / std::max(0.01f, survival.shieldParryWindow)); // the window running out: the arc shrinks
        DrawRing({ cx, cy }, r - 3.0f * s, r, -90.0f, -90.0f + 360.0f * left, 48, Color{ 255, 205, 70, 235 });
        const int fs = (int)(16 * s);
        const char* t = "PARRY";
        DrawText(t, (int)(cx - MeasureText(t, fs) * 0.5f), (int)(cy + r + 4.0f * s), fs, Color{ 255, 205, 70, 235 });
    } else {
        DrawRing({ cx, cy }, r - 2.0f * s, r, 0.0f, 360.0f, 48, Fade(LIGHTGRAY, 0.55f));
    }
}

void GameUI_DrawHUD(Scene& world) {
    DrawBlockRing(); // (#1713)
    DrawTeamStrip();
    DrawTierBanner(); // (#226)
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
        const float speed = Vector3Length(player.groundVelocity(world));
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
        if (GetTime() < toastUntil) { const int fs = (int)(20 * GameUI_HudScale()); DrawText(toast.c_str(), (GetScreenWidth() - MeasureText(toast.c_str(), fs)) / 2, (int)(GetScreenHeight() * 0.70f), fs, YELLOW); } // (#891: at the scale, above the hotbar)
        return;
    }
    // Jetpack fuel bar (top right): only with a jetpack (a thruster pack and suit power, #217).
    if (player.jetpackAccel > 0.0f) {
        DrawRectangle(w - 220, 70, 200, 12, Fade(BLACK, 0.5f));
        DrawRectangle(w - 220, 70, (int)(200 * player.jetpackFuel), 12, SKYBLUE);
        DrawText("JET", w - 260, 68, 16, RAYWHITE);
    }
    // Under the minimap (#574): where you are (latitude, longitude and height over the body whose sphere of influence you
    // are in; world kilometres out in space), then speed and gravity, then the biome and the weather. (#932, the user) Only
    // while the minimap is on screen (a Map carried): it reads as the map's caption.
    if (Nav_MinimapShown()) {
        int ly = 384; // the minimap and its scale line end at ~372
        const int lfs = (int)(16 * GameUI_HudScale()); auto line = [&](const char* t, Color c) { DrawText(t, w - 20 - MeasureText(t, lfs), ly, lfs, c); ly += lfs + 4; }; // (#891: at the scale) // right-aligned with the minimap
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
        const char* mode = player.flyMode && terrainTool.creative ? "FLY (F)" : player.grounded ? "WALK" : player.hovering ? "HOVER (H)" : player.gravity >= 0.5f ? "AIR" : "SPACE"; // (#932: fly is creative only)
        line(TextFormat("%s   g %.1f   %.1f m/s", mode, player.gravity, Vector3Length(player.groundVelocity(world))), RAYWHITE);
        const char* biome = Survival_BiomeName();
        const bool weather = Weather_Here().kind != WeatherKind::Clear && biome[0];
        if (biome[0]) { line(weather ? TextFormat("%s   %s", biome, Weather_Name(Weather_Here().kind)) : biome, Fade(RAYWHITE, 0.8f)); }
        if (terrainTool.creative && player.nearest && LightLevel_On()) { // (#1739) the light grid at head height (creative: a debugging readout)
            const Vector3 feet = player.nearest->toLocalOf(player.position, world.renderOrigin);
            const LightCell l = LightLevel_At(*player.nearest, Vector3Add(feet, Vector3Scale(Vector3Normalize(feet), 1.5f)));
            line(TextFormat("light %d (%.1f m)%s", l.level, l.level / (float)LightGrid::kStep, l.level <= 16 ? "  dark" : ""), Fade(RAYWHITE, 0.8f));
        }
        if (player.jetpackAccel > 0.0f) // only with a jetpack: nothing else uses them (#531)
            line(TextFormat("dampeners %s (B)", player.dampeners ? "on" : "off"), player.dampeners ? RAYWHITE : ORANGE);
    }
    if (GetTime() < toastUntil) { const int fs = (int)(20 * GameUI_HudScale()); DrawText(toast.c_str(), (GetScreenWidth() - MeasureText(toast.c_str(), fs)) / 2, (int)(GetScreenHeight() * 0.70f), fs, YELLOW); } // (#891: at the scale, above the hotbar)
}

// (#824) the crafting queue from outside (tests; later, e.g. a recipe's "craft" from the database pop-up, #826)
void GameUI_QueueCraft(int recipe, int count) { if (count > 0) craftQueue.push_back({ recipe, count }); }
bool GameUI_QueueRecursive(int recipe, int count, int* steps) { // (#928) as the page's buttons: its intermediates first
    const auto& recipes = itemDB.recipes();
    std::set<std::string> machineMade;
    for (const auto& ri : Factory::RecipeList()) for (const auto& o : ri.out) machineMade.insert(o.first);
    std::map<std::string, int> carried;
    for (const auto& kv : inventory.all()) carried[kv.first] = kv.second;
    const std::function<bool(int)> known = [&](int k) { return k >= 0 && k < (int)recipes.size() && RecipeBook_Known(recipes[k]); };
    const CraftPlan p = CraftPlan_Make(recipes, known, carried, recipe, count, machineMade);
    if (steps) *steps = (int)p.steps.size();
    if (p.ok && !Crafting_RequestRecursive(recipe, count)) QueuePlan(p); // (#946: joined, the server plans and crafts it)
    return p.ok;
}
int GameUI_CraftQueueLength() { int n = 0; for (const QueuedCraft& c : craftQueue) n += c.left; return n; }
void GameUI_SelectRecipe(int recipe) { selectedRecipe = recipe; recipePageOpen = recipe >= 0; } // (#921: its page open)
bool GameUI_ClosePopout() { // (#921) Esc: a pop-out closes before the window does
    if (inventoryOpen && recipePageOpen && screen == GameScreen::Crafting) { recipePageOpen = false; return true; }
    return inventoryOpen && screen == GameScreen::Quests && QuestTree_ClosePopout();
}

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
    // (#1713, playtest: "can't find how to parry") A weapon swings; a shield in the off hand blocks with RMB, and raised just
    // before a hit, parries (the ring round the crosshair shows the window: DrawBlockRing).
    const bool weapon = terrainTool.heldTool.rfind("weapon_", 0) == 0;
    std::string line = terrainTool.terraform ? "LMB dig  RMB fill  Ctrl+wheel brush"
                     : std::string(weapon ? "LMB swing" : "LMB mine/chop") + (terrainTool.rmbBlocks ? "  RMB block (just before a hit: parry)" : "  RMB place");
    line += TextFormat("  %s build  1-0 hotbar  %s inventory  %s telescope  %s lamp", K(Act::BuildMode).c_str(), K(Act::Inventory).c_str(), K(Act::Telescope).c_str(), K(Act::Headlamp).c_str());
    if (terrainTool.creative) line += TextFormat("  %s fly", K(Act::Fly).c_str()); // fly is creative only (#531)
    line += TextFormat("  %s save  %s commands", K(Act::Save).c_str(), K(Act::Chat).c_str());
    const float s = GameUI_HudScale(); // (#925) the HUD's size
    int fs = (int)(20 * s);
    while (fs > 12 && x + MeasureText(line.c_str(), fs) > GetScreenWidth() - (int)(8 * s)) fs -= 2; // (#891: as large as the screen's width allows)
    DrawText(line.c_str(), x, std::min(y, GetScreenHeight() - fs - (int)(8 * s)), fs, Fade(RAYWHITE, alpha)); // (larger text: kept on screen)
}
void GameUI_InstallKeyHints() { gKeyHintsHook = DrawKeyHints; }

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#892): dragging over the player preview turns it (a headless ImGui context, the mouse injected) ----
#include "unit_helpers.h"

TEST(GameUiPreview, DraggingAcrossItTurnsThePlayer) {
    ImGuiContext* prev = ImGui::GetCurrentContext();
    ImGuiContext* ctx = ImGui::CreateContext();
    ImGui::SetCurrentContext(ctx);
    ImGuiIO& io = ImGui::GetIO();
    io.DisplaySize = ImVec2(800.0f, 600.0f);
    io.DeltaTime = 1.0f / 60.0f;
    io.IniFilename = nullptr;
    io.BackendFlags |= ImGuiBackendFlags_RendererHasTextures; // (as rlImGui sets: the 1.92 font atlas builds as fonts are used)
    io.Fonts->AddFontDefault();
    float turn = 0.0f;
    auto frame = [&](float mx, float my, bool down) {
        io.AddMousePosEvent(mx, my);
        io.AddMouseButtonEvent(0, down);
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f));
        ImGui::SetNextWindowSize(ImVec2(800.0f, 600.0f));
        ImGui::Begin("preview test", nullptr, ImGuiWindowFlags_NoDecoration);
        ImGui::SetCursorScreenPos(ImVec2(100.0f, 100.0f));
        PreviewWidget((ImTextureID)1, ImVec2(200.0f, 300.0f), turn);
        ImGui::End();
        ImGui::Render();
    };
    frame(200.0f, 250.0f, false);                         // over it
    frame(200.0f, 250.0f, true);                          // pressed
    for (int k = 1; k <= 10; ++k) frame(200.0f + 10.0f * k, 250.0f, true); // dragged 100 px right
    frame(300.0f, 250.0f, false);                         // released
    EXPECT_NEAR(turn, -1.2f, 0.05f) << "100 px at 0.012 rad a pixel; dragging right lowers the turn (#1073: it had spun the other way)";
    const float after = turn;
    frame(500.0f, 250.0f, false); frame(600.0f, 250.0f, false); // moving without the button: no turn
    EXPECT_FLOAT_EQ(turn, after);
    ImGui::DestroyContext(ctx);
    ImGui::SetCurrentContext(prev);
}

// (#1078) A drag in the quest tree pans it; it must not move the Tab window. ImGui moves a movable window from a press on
// no item, so the tree's canvas is an InvisibleButton (questtree.cpp). Without it (the control) the window moves.
TEST(GameUiQuestTree, ADragOnTheCanvasDoesNotMoveTheWindow) {
    auto windowAfterDrag = [](bool canvas) {
        ImGuiContext* prev = ImGui::GetCurrentContext();
        ImGuiContext* ctx = ImGui::CreateContext();
        ImGui::SetCurrentContext(ctx);
        ImGuiIO& io = ImGui::GetIO();
        io.DisplaySize = ImVec2(1280.0f, 720.0f); io.DeltaTime = 1.0f / 60.0f; io.IniFilename = nullptr;
        io.BackendFlags |= ImGuiBackendFlags_RendererHasTextures;
        io.Fonts->AddFontDefault();
        ImVec2 pos;
        auto frame = [&](float mx, float my, bool down) {
            io.AddMousePosEvent(mx, my); io.AddMouseButtonEvent(0, down);
            ImGui::NewFrame();
            ImGui::SetNextWindowPos(ImVec2(100.0f, 100.0f), ImGuiCond_FirstUseEver);
            ImGui::SetNextWindowSize(ImVec2(600.0f, 400.0f), ImGuiCond_FirstUseEver);
            ImGui::Begin("tab window", nullptr, ImGuiWindowFlags_NoTitleBar); // movable, as the Tab window (#1080)
            ImGui::BeginChild("questtree", ImVec2(0, 0), true, ImGuiWindowFlags_HorizontalScrollbar);
            if (canvas) ImGui::InvisibleButton("canvas", ImVec2(1200.0f, 900.0f), ImGuiButtonFlags_MouseButtonLeft | ImGuiButtonFlags_MouseButtonMiddle);
            else ImGui::Dummy(ImVec2(1200.0f, 900.0f));
            ImGui::EndChild();
            pos = ImGui::GetWindowPos();
            ImGui::End();
            ImGui::Render();
        };
        frame(300.0f, 300.0f, false);
        frame(300.0f, 300.0f, true);
        for (int k = 1; k <= 8; ++k) frame(300.0f + 10.0f * k, 300.0f + 5.0f * k, true);
        frame(380.0f, 340.0f, false);
        ImGui::DestroyContext(ctx);
        ImGui::SetCurrentContext(prev);
        return pos;
    };
    const ImVec2 moved = windowAfterDrag(false), kept = windowAfterDrag(true);
    EXPECT_GT(moved.x, 150.0f) << "the control: a press on no item drags the window (the playtest bug)";
    EXPECT_FLOAT_EQ(kept.x, 100.0f) << "a drag starting on the tree's canvas leaves the Tab window where it is";
    EXPECT_FLOAT_EQ(kept.y, 100.0f);
}

// ---- Unit tests (#894): the press that picks a stack up is not a click outside the window ----
#include "unit_helpers.h"

namespace {
struct HeadlessImGui { // a context with a 1280 x 720 display and the mouse injected per frame
    ImGuiContext* prev = nullptr;
    ImGuiContext* ctx = nullptr;
    HeadlessImGui() {
        prev = ImGui::GetCurrentContext();
        ctx = ImGui::CreateContext();
        ImGui::SetCurrentContext(ctx);
        ImGuiIO& io = ImGui::GetIO();
        io.DisplaySize = ImVec2(1280.0f, 720.0f); io.DeltaTime = 1.0f / 60.0f; io.IniFilename = nullptr;
        io.BackendFlags |= ImGuiBackendFlags_RendererHasTextures; // (as rlImGui sets: the 1.92 font atlas builds as fonts are used)
        io.Fonts->AddFontDefault();
    }
    ~HeadlessImGui() { ImGui::DestroyContext(ctx); ImGui::SetCurrentContext(prev); }
    // One frame: a window at (0,0) 800 x 600 with a cell button at (100,100); after it, `after` runs (as GameUI_Draw's tail).
    template <class F> void Frame(float mx, float my, bool down, F after) {
        ImGuiIO& io = ImGui::GetIO();
        io.AddMousePosEvent(mx, my); io.AddMouseButtonEvent(0, down);
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f)); ImGui::SetNextWindowSize(ImVec2(800.0f, 600.0f));
        ImGui::Begin("game window", nullptr, ImGuiWindowFlags_NoDecoration);
        ImGui::SetCursorScreenPos(ImVec2(100.0f, 100.0f));
        ImGui::InvisibleButton("cell", ImVec2(54.0f, 54.0f));
        ImGui::End();
        after();
        ImGui::Render();
    }
};
}

TEST(GameUiCursor, ThePressThatPicksAStackUpIsNotAClickOutside) {
    HeadlessImGui ui;
    bool oldCheckSaysOutside = false, outside = true;
    ui.Frame(120.0f, 120.0f, false, [] {});
    ui.Frame(120.0f, 120.0f, true, [&] { // the press on the cell (it picked the stack up: its button is active now)
        oldCheckSaysOutside = !ImGui::IsWindowHovered(ImGuiHoveredFlags_AnyWindow);
        outside = PointerOutsideWindows();
    });
    EXPECT_TRUE(oldCheckSaysOutside) << "the cause of #894: an active button hides the window from IsWindowHovered(AnyWindow)";
    EXPECT_FALSE(outside) << "over the window: the stack stays on the cursor";
    ui.Frame(120.0f, 120.0f, false, [&] { EXPECT_FALSE(PointerOutsideWindows()) << "the release on the cell"; });
}

TEST(GameUiCursor, AClickOutsideEveryWindowIsOutside) {
    HeadlessImGui ui;
    bool outside = false;
    ui.Frame(1000.0f, 650.0f, false, [] {});
    ui.Frame(1000.0f, 650.0f, true, [&] { outside = PointerOutsideWindows(); });
    EXPECT_TRUE(outside) << "past the window: a held stack is dropped at the feet";
}

// ---- Unit tests (#921): the recipe page pop-up's clicks (creative: every recipe known, crafts free) ----
TEST(GameUiRecipePage, AClickOpensItBesideTheRecipeCraftOneQueuesAndEscClosesItFirst) {
    UnitTest::EnsureItemDB();
    HeadlessImGui ui;
    const bool creative0 = terrainTool.creative;
    terrainTool.creative = true;
    craftQueue.clear(); recipePageOpen = false; selectedRecipe = -1;
    GameUI_ShowScreen(GameScreen::Crafting);
    auto frame = [&](ImVec2 m, bool down) {
        ImGuiIO& io = ImGui::GetIO();
        io.AddMousePosEvent(m.x, m.y); io.AddMouseButtonEvent(0, down);
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f)); ImGui::SetNextWindowSize(ImVec2(1280.0f, 720.0f));
        ImGui::Begin("game window", nullptr, ImGuiWindowFlags_NoDecoration);
        ImGui::SetCursorPosX(640.0f); // the crafting column on the right half, as the Crafting screen lays it out
        ImGui::BeginChild("craftcol", ImVec2(0.0f, 0.0f));
        DrawCrafting();
        ImGui::EndChild();
        ImGui::End();
        ImGui::Render();
    };
    auto click = [&](ImVec2 m) { frame(m, false); frame(m, true); frame(m, false); frame(m, false); };
    const ImVec2 away(5.0f, 5.0f);
    frame(away, false); frame(away, false);
    EXPECT_FALSE(recipePageOpen) << "nothing open until a recipe is clicked";
    ASSERT_GE(recipeTest.firstRecipe, 0) << "creative: recipes are drawn";

    click(recipeTest.firstCell);
    EXPECT_TRUE(recipePageOpen);
    EXPECT_EQ(selectedRecipe, recipeTest.firstRecipe);
    EXPECT_TRUE(recipePageShown) << "the pop-up is drawn";
    const float cellX = recipeTest.firstCell.x;
    EXPECT_TRUE(recipeTest.pageMax.x <= cellX || recipeTest.pageMin.x >= cellX) << "beside the recipe, not over it";
    EXPECT_NEAR(recipeTest.pageMin.y, recipeAnchorMin.y, 2.0f) << "level with the recipe's cell";

    const int queued0 = GameUI_CraftQueueLength();
    click(recipeTest.craft1);
    EXPECT_EQ(GameUI_CraftQueueLength(), queued0 + 1) << "Craft 1 queues one";
    EXPECT_TRUE(recipePageOpen) << "crafting keeps the page open";

    click(recipeTest.firstCell); // (#1149) its icon again, the page open: one more craft
    EXPECT_EQ(GameUI_CraftQueueLength(), queued0 + 2) << "a click on the open recipe's icon crafts one (it closed the page)";
    EXPECT_TRUE(recipePageOpen) << "and the page stays";
    EXPECT_EQ(selectedRecipe, recipeTest.firstRecipe);

    EXPECT_TRUE(GameUI_ClosePopout()) << "Esc: the pop-up closes first";
    frame(away, false); frame(away, false);
    EXPECT_FALSE(recipePageOpen);
    EXPECT_FALSE(recipePageShown);
    EXPECT_TRUE(GameUI_InventoryOpen()) << "the window stays";
    EXPECT_FALSE(GameUI_ClosePopout()) << "Esc again: nothing open, the window closes (main.cpp toggles it)";

    craftQueue.clear(); selectedRecipe = -1;
    GameUI_ToggleInventory();
    terrainTool.creative = creative0;
}
#endif

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1730): an item's own record (InvSlot::data) survives every move ----
#include "unit_helpers.h"
#include <cstring>

namespace {
const char* const kRec = "h=bronze;g=wood;b=leather;q=112;s=1;w=37";
const char* const kRec2 = "h=iron;g=oak;b=leather;q=100;s=0;w=0";
const char* const kRec3 = "h=copper;g=pine;b=sinew;q=90;s=0;w=5";

struct InventoryMoves : ::testing::Test {
    void SetUp() override {
        UnitTest::EnsureItemDB();
        cursorStack = InvSlot{}; binned = InvSlot{};
        SlotRules_SetHeld(&cursorStack);
        inventory.clear();
    }
    void TearDown() override { cursorStack = InvSlot{}; binned = InvSlot{}; inventory.clear(); SlotRules_SetHeld(nullptr); SlotRules_SetGrids({}); }
    // Every place a copy can be: the cells, the cursor, the bin. A move must leave exactly one.
    static int Copies(const std::string& rec) {
        int n = 0;
        for (const InvSlot& s : inventory.grid()) if (!s.empty() && s.data == rec) ++n;
        if (cursorStack.data == rec && !cursorStack.empty()) ++n;
        if (binned.data == rec && !binned.empty()) ++n;
        return n;
    }
    static int Put(const char* id, const char* rec) { ++inventory.muted; const int c = inventory.addInstance(id, rec); --inventory.muted; return c; }
    static int Free() { for (int k = 0; k < Inventory::kGrid; ++k) if (inventory.grid()[k].empty()) return k; return -1; }
    static int FreeAfter(int not1, int not2 = -1) { for (int k = 0; k < Inventory::kGrid; ++k) if (inventory.grid()[k].empty() && k != not1 && k != not2) return k; return -1; }
    static int Holding(const std::string& rec) { for (int k = 0; k < Inventory::kCells; ++k) if (!inventory.grid()[k].empty() && inventory.grid()[k].data == rec) return k; return -1; }
};
}

TEST_F(InventoryMoves, CursorPickDropAndSwapKeepTheRecord) {
    const int a = Put("tool_pickaxe", kRec);
    ASSERT_GE(a, 0);
    CellClick(a, false, false);
    EXPECT_EQ(cursorStack.id, "tool_pickaxe");
    EXPECT_EQ(cursorStack.data, kRec) << "the pick-up took the record with it (#1730)";
    EXPECT_EQ(cursorStack.centi, 100);
    EXPECT_TRUE(inventory.grid()[a].empty());
    const int b = Free();
    CellClick(b, false, false);
    EXPECT_EQ(inventory.grid()[b].data, kRec);
    EXPECT_TRUE(cursorStack.empty());
    EXPECT_EQ(Copies(kRec), 1);

    // onto another item: they trade places
    ++inventory.muted; inventory.add("stone", 10); --inventory.muted;
    int c = -1; for (int k = 0; k < Inventory::kCells; ++k) if (inventory.grid()[k].id == "stone") c = k;
    ASSERT_GE(c, 0);
    CellClick(b, false, false); CellClick(c, false, false);
    EXPECT_EQ(cursorStack.id, "stone"); EXPECT_EQ(cursorStack.centi, 1000); EXPECT_TRUE(cursorStack.data.empty());
    EXPECT_EQ(inventory.grid()[c].data, kRec);
    CellClick(b, false, false);
    EXPECT_EQ(inventory.count("stone"), 10);
    EXPECT_EQ(Copies(kRec), 1);

    // onto another copy of the same item with another record: a swap, never a merge
    const int d = Put("tool_pickaxe", kRec2);
    ASSERT_GE(d, 0);
    CellClick(c, false, false); CellClick(d, false, false);
    EXPECT_EQ(cursorStack.data, kRec2); EXPECT_EQ(inventory.grid()[d].data, kRec);
    CellClick(c, false, false);
    EXPECT_EQ(Copies(kRec), 1); EXPECT_EQ(Copies(kRec2), 1);
    EXPECT_EQ(inventory.grid()[c].data, kRec2);
    EXPECT_EQ(inventory.countInstances("tool_pickaxe"), 2);
}

TEST_F(InventoryMoves, RightClickAndDragsNeverSplitOrDuplicateIt) {
    const int a = Put("tool_pickaxe", kRec);
    CellClick(a, true, false); // right-click: "half" of one unit is the whole item
    EXPECT_EQ(cursorStack.data, kRec); EXPECT_EQ(cursorStack.centi, 100);
    EXPECT_TRUE(inventory.grid()[a].empty());
    const int b = Free();
    CellClick(b, true, false); // right-click holding: one down
    EXPECT_EQ(inventory.grid()[b].data, kRec); EXPECT_TRUE(cursorStack.empty());

    // a left drag across three cells and a right drag: the one item lands once
    for (int pass = 0; pass < 2; ++pass) {
        const bool right = pass == 1;
        const int from = Holding(kRec);
        CellClick(from, false, false);
        ASSERT_EQ(cursorStack.data, kRec);
        const int f0 = Free(), f1 = FreeAfter(f0), f2 = FreeAfter(f0, f1);
        SlotRules_SetGrids({ &playerGrid });
        SlotInput p; (right ? p.rPress : p.lPress) = true; (right ? p.rDown : p.lDown) = true;
        SlotRules_Frame(p, &playerGrid, f0);
        SlotInput m; (right ? m.rDown : m.lDown) = true;
        SlotRules_Frame(m, &playerGrid, f1); SlotRules_Frame(m, &playerGrid, f2);
        SlotInput r; (right ? r.rRelease : r.lRelease) = true;
        SlotRules_Frame(r, &playerGrid, f2);
        EXPECT_EQ(Copies(kRec), 1) << (right ? "right drag" : "left drag");
        EXPECT_TRUE(cursorStack.empty());
        EXPECT_GE(Holding(kRec), 0);
    }
}

TEST_F(InventoryMoves, ItemsWithARecordNeverMergeOrGather) {
    const int s = Put("stone", kRec);
    ++inventory.muted; inventory.add("stone", 10); --inventory.muted;
    int plain = -1; for (int k = 0; k < Inventory::kCells; ++k) if (inventory.grid()[k].id == "stone" && !inventory.grid()[k].hasData()) plain = k;
    ASSERT_GE(s, 0); ASSERT_GE(plain, 0);
    // the record's stone onto the plain stack: a swap, not a stack of 11
    CellClick(s, false, false); CellClick(plain, false, false);
    EXPECT_EQ(cursorStack.data, ""); EXPECT_EQ(cursorStack.centi, 1000);
    EXPECT_EQ(inventory.grid()[plain].data, kRec); EXPECT_EQ(inventory.grid()[plain].centi, 100);
    CellClick(s, false, false); // the plain stack down again
    // the plain stack onto the record's stone, with right-click
    CellClick(s, false, false); CellClick(plain, true, false);
    EXPECT_EQ(cursorStack.data, kRec); EXPECT_EQ(cursorStack.centi, 100);
    EXPECT_EQ(inventory.grid()[plain].centi, 1000);
    EXPECT_TRUE(inventory.grid()[plain].data.empty());
    // a double click while holding it gathers nothing
    SlotRules_SetGrids({ &playerGrid });
    SlotInput dbl; dbl.lPress = dbl.lDown = dbl.lDouble = true;
    SlotRules_Frame(dbl, &playerGrid, plain);
    EXPECT_EQ(cursorStack.data, kRec); EXPECT_EQ(cursorStack.centi, 100);
    EXPECT_EQ(inventory.grid()[plain].centi, 1000);
    EXPECT_EQ(Copies(kRec), 1);
}

TEST_F(InventoryMoves, ShiftClickAndHotbarKeysKeepTheRecord) {
    const int tool = Put("tool_pickaxe", kRec);   // a tool lands on the hotbar
    ASSERT_GE(tool, Inventory::kGrid);
    CellClick(tool, false, true);                  // Shift-click: hotbar -> grid
    int g = Holding(kRec);
    ASSERT_GE(g, 0); EXPECT_LT(g, Inventory::kGrid);
    CellClick(g, false, true);                     // and back
    g = Holding(kRec);
    EXPECT_GE(g, Inventory::kGrid);
    EXPECT_EQ(Copies(kRec), 1);

    const int s = Put("stone", kRec2);
    ASSERT_GE(s, 0);
    CellClick(s, false, true);
    EXPECT_EQ(Copies(kRec2), 1);

    // the 1-0 keys and the off-hand key swap the cell with that cell
    SlotRules_SetGrids({ &playerGrid });
    const int from = Holding(kRec2);
    SlotInput k; k.hotbarKey = 3;
    SlotRules_Frame(k, &playerGrid, from);
    EXPECT_EQ(inventory.grid()[Inventory::kGrid + 3].data, kRec2);
    SlotInput o; o.offHandKey = true;
    SlotRules_Frame(o, &playerGrid, Inventory::kGrid + 3);
    EXPECT_EQ(inventory.grid()[Inventory::kOffHand].data, kRec2);
    EXPECT_EQ(Copies(kRec2), 1);
    EXPECT_EQ(Copies(kRec), 1);
}

TEST_F(InventoryMoves, TheBinTheWindowCloseAndAServerTakeKeepTheRecord) {
    const int a = Put("tool_pickaxe", kRec);
    CellClick(a, false, false);
    BinClick();
    EXPECT_TRUE(cursorStack.empty()); EXPECT_EQ(binned.data, kRec);
    BinClick();
    EXPECT_EQ(cursorStack.data, kRec); EXPECT_TRUE(binned.empty());
    ReturnCursorStack(); // the window closes
    EXPECT_TRUE(cursorStack.empty());
    EXPECT_EQ(Copies(kRec), 1);
    EXPECT_GE(Holding(kRec), 0);

    // a take from a chest or machine arrives through the cursor sink: with its record, onto an empty cursor
    CursorSink("tool_pickaxe", 1, 0.0f, kRec2);
    EXPECT_EQ(cursorStack.data, kRec2); EXPECT_EQ(cursorStack.centi, 100);
    CursorSink("tool_pickaxe", 1, 0.0f, kRec3); // the cursor is busy: into the inventory as it is
    EXPECT_EQ(cursorStack.data, kRec2);
    EXPECT_EQ(inventory.countInstances("tool_pickaxe"), 2);
    CellClick(Free(), false, false);
    EXPECT_EQ(Copies(kRec), 1); EXPECT_EQ(Copies(kRec2), 1); EXPECT_EQ(Copies(kRec3), 1);
}

TEST_F(InventoryMoves, DroppedItemsKeepTheirRecordOnTheGroundAndInTheSave) {
    Combat_Reset();
    InvSlot a; a.id = "tool_pickaxe"; a.centi = 100; a.data = kRec;
    InvSlot b = a; b.data = kRec2;
    InvSlot plain; plain.id = "stone"; plain.centi = 300;
    const DVec3 spot{ 0.0, 60000.0, 0.0 };
    ASSERT_TRUE(Combat_DropFor("terra", spot, a, ""));
    ASSERT_TRUE(Combat_DropFor("terra", spot, b, "")) << "same item, another record: its own entry in the pile";
    ASSERT_TRUE(Combat_DropFor("terra", spot, plain, ""));
    ASSERT_TRUE(Combat_DropFor("terra", spot, plain, ""));
    EXPECT_EQ(Combat_DroppedPiles(), 1);

    // the pile survives a save and a load
    const std::string path = UnitTest::ScratchDir() + "/pile.txt";
    FILE* f = fopen(path.c_str(), "w+b");
    ASSERT_NE(f, nullptr);
    Combat_Save(f);
    fflush(f); rewind(f);
    Combat_Reset();
    char tok[64];
    while (fscanf(f, " %63s", tok) == 1) Combat_LoadLine(tok, f);
    fclose(f);
    EXPECT_EQ(Combat_DroppedPiles(), 1);

    // walking away arms it, walking back takes everything: both records, one each, and the plain stone together
    Inventory inv;
    const DVec3 away{ 0.0, 60000.0, 50.0 };
    Combat_PickUpGravesFor("terra", away, inv, "", false, false);
    EXPECT_TRUE(Combat_PickUpGravesFor("terra", spot, inv, "", false, false));
    EXPECT_EQ(inv.countInstances("tool_pickaxe"), 2);
    EXPECT_GE(inv.findInstance("tool_pickaxe", kRec), 0);
    EXPECT_GE(inv.findInstance("tool_pickaxe", kRec2), 0);
    EXPECT_EQ(inv.count("stone"), 6);
    EXPECT_EQ(Combat_DroppedPiles(), 0);
    Combat_Reset();
}
#endif
