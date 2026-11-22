#include "hotbar.h"
#include "melee.h"
#include "machinetiers.h"
#include "uipalette.h"
#include "playersession.h"
#include "playerinput.h"
#include "world.h" // GameFrameTime
#include "toolwear.h"
#include "composed.h"
#include "buildtool.h"
#include "itemicons.h"
#include "gameui.h"
#include "inventory.h"
#include "plots.h"
#include "survival.h"
#include "audio.h"
#include "terraintool.h"
#include "tools.h"
#include "animals.h"
#include "telescope.h"
#include "worldtypes.h"
#include <raymath.h>
#include <memory>

#include <raylib.h>
#include <cstring>
#include <string>
#include <algorithm>

const char* const kHotbarTerrainTool = "tool:terrain";

namespace {
    HotbarState& H() { return LocalSession().hotbar; } // the local player's (#313)
    // Eating (#202): right-click held with food in hand.
    constexpr float kEatSeconds = 1.2f;
    void UpdateEating(float dt, bool held) {
        const std::string food = held ? Hotbar_HeldFood() : std::string();
        if (!held || food.empty() || !Survival_IsFood(food)) { H().eatProgress = 0.0f; H().eatBlocked = false; return; }
        if (H().eatBlocked) return;
        if (H().eatProgress == 0.0f && survival.food >= 98.0f) { GameUI_Toast("You are full"); H().eatBlocked = true; return; }
        gPrimaryLatch.consume(); // (#896) eating: LMB strikes nothing, and one held now waits for a release afterwards
        H().eatProgress += dt / kEatSeconds;
        if (H().eatProgress >= 1.0f) {
            if (Survival_UseItem(food)) ++H().eatenCount;
            H().eatProgress = 0.0f; // holding on eats the next one
        }
    }
    // The hotbar is the inventory's last 10 cells (#209): a slot holds a real stack; empty = bare hands.
    const std::string kNone;
    const std::string& CellId(int slot) {
        const InvSlot& s = inventory.grid()[Inventory::kGrid + slot];
        return s.empty() ? kNone : s.id;
    }
    const int kKeys[kHotbarSlots] = { KEY_ONE, KEY_TWO, KEY_THREE, KEY_FOUR, KEY_FIVE, KEY_SIX, KEY_SEVEN, KEY_EIGHT, KEY_NINE, KEY_ZERO };

    // A part type placed from this item (machines, conveyors, cables), or Count.
    PartType PartForItem(const std::string& id, int* step = nullptr) {
        if (step) *step = 0;
        { PartType vt; int vs = 0; if (MachineTiers_ItemPart(id, vt, vs)) { if (step) *step = vs; return vt; } } // (#1427 L7a) a machine variant
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
        const std::string id = CellId(slot); // a copy: using the item can empty the cell
        if (id.empty()) { // an empty slot is bare hands
            if (piloting) return;
            H().selected = slot; buildTool.active = false; terrainTool.heldTool.clear(); terrainTool.holdingItem = false;
            return;
        }
        // Food and animal feed: the first press takes it in hand (animals that eat it follow you); the
        // next press eats it (food only).
        if ((Survival_IsFood(id) || Animals_IsFeed(id)) && !piloting) {
            if (H().selected == slot && buildTool.active == false && terrainTool.heldTool.empty()) { // eating is a held right-click now (#202)
                if (Survival_IsFood(id)) GameUI_Toast(TextFormat("Hold right-click to eat the %s", ShortName(id)));
                return;
            }
            if (inventory.held(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); return; }
            H().selected = slot; buildTool.active = false; terrainTool.heldTool.clear(); terrainTool.holdingItem = true; // RMB: no terrain (#201)
            GameUI_Toast(Survival_IsFood(id) ? TextFormat("Holding %s: hold right-click to eat", ShortName(id)) : TextFormat("Holding %s: animals that eat it follow you", ShortName(id)));
            return;
        }
        if ((Hotbar_IsLight(id) || Hotbar_IsPlaceable(id)) && !piloting) { // a torch (#218): held for its light; a snare (#176): set down with RMB; LMB still digs with bare hands
            if (inventory.held(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); return; }
            H().selected = slot; buildTool.active = false; terrainTool.heldTool.clear(); terrainTool.holdingItem = true;
            return;
        }
        if (Survival_IsUsable(id)) { Survival_UseItem(id); return; } // consumables don't change the selection
        if (piloting) return;
        // The selection moves only once the choice is valid (#272): an error (not carried, can't be used) leaves the
        // hand and build mode as they were.
        auto take = [&]() { H().selected = slot; terrainTool.holdingItem = false; };
        if (id == kHotbarTerrainTool) { take(); buildTool.active = false; terrainTool.heldTool.clear(); return; }
        if (id == "bucket" || id == "water_bucket" || id == "fishing_rod") { // in hand: RMB scoops / pours water (terraintool), or casts (fishing.cpp, #176)
            if (inventory.held(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
            take();
            buildTool.active = false;
            terrainTool.heldTool = id;
            return;
        }
        if (KindOf(id) == ItemKind::Tool && (ToolOf(id).type != ToolType::Hands || Weapon_ClassOf(id) != WeaponClass::Tools)) { // (#1219) a weapon goes in the hand too (it digs as hands)
            if (inventory.held(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
            take();
            buildTool.active = false;
            terrainTool.heldTool = id;
            return;
        }
        if (TerrainTool_IsFiller(id) && BuildMaterialForItem(id.c_str()) < 0) { // ground (loose fillers and raw stone, #1098): in hand, RMB places it
            if (inventory.held(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
            take();
            buildTool.active = false;
            terrainTool.heldTool = id;
            return;
        }
        const int m = BuildMaterialForItem(id.c_str());
        if (m >= 0) {
            take();
            buildTool.active = true;
            buildTool.material = m;
            if (!IsShape((PartType)buildTool.part)) buildTool.part = (int)PartType::Cube;
            return;
        }
        int variantStep = 0;
        const PartType t = PartForItem(id, &variantStep);
        if (t != PartType::Count) {
            take();
            buildTool.active = true;
            buildTool.part = (int)t;
            buildTool.variant = variantStep; // (#1427 L7a) the variant: its step rides in the part's material
            return;
        }
        // (#534, the user) Anything else is just held: LMB punches (bare hands), RMB does nothing unless it has a use (the
        // telescope's RMB raises it: Telescope_Update).
        if (inventory.held(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
        take();
        buildTool.active = false;
        terrainTool.heldTool.clear();
        terrainTool.holdingItem = true;
    }
}

void Hotbar_Reset() {
    // The cells are the inventory's (cleared with it); an empty slot is bare hands, new tools land in empty slots.
    H().selected = 0;
    terrainTool.holdingItem = false;
}

void Hotbar_Bind(int slot, const std::string& id) {
    // Moves a stack of the item into the slot's cell (swapping what was there out), the main grid's first.
    if (slot < 0 || slot >= kHotbarSlots) return;
    auto& g = inventory.grid();
    const int cell = Inventory::kGrid + slot;
    if (!g[cell].empty() && g[cell].id == id) return;
    int from = -1;
    for (int i = 0; i < Inventory::kSlots && from < 0; ++i) if (i != cell && !g[i].empty() && g[i].id == id) from = i;
    if (from < 0) return; // not carried
    if (!g[cell].empty() && g[from].id != g[cell].id) std::swap(g[from], g[cell]); else inventory.swapSlots(from, cell);
    GameUI_Toast(TextFormat("%s on key %d", ShortName(id), (slot + 1) % 10));
}

const std::string& Hotbar_Slot(int slot) { return slot >= 0 && slot < kHotbarSlots ? CellId(slot) : kNone; }
int Hotbar_Selected() { return H().selected; }
std::string Hotbar_HeldFood() {
    if (H().selected < 0 || H().selected >= kHotbarSlots || buildTool.active || !terrainTool.heldTool.empty()) return "";
    const std::string& id = CellId(H().selected);
    return (Survival_IsFood(id) || Animals_IsFeed(id)) && inventory.count(id) > 0 ? id : "";
}

int Hotbar_KeySlot() {
    for (int i = 0; i < kHotbarSlots; ++i) if (IsKeyPressed(kKeys[i])) return i;
    return -1;
}

void Hotbar_Press(int slot, bool piloting) { if (slot >= 0 && slot < kHotbarSlots) Select(slot, piloting); }

void Hotbar_Update(Scene&, const PlayerInput& in, bool piloting) {
    static bool initialised = false;
    if (!initialised) { Hotbar_Reset(); initialised = true; }
    // G (build mode on/off) is the build tool's; keep the selection in step with it.
    // Leaving build mode with a block selected: the hand goes empty (to an empty slot, if there is one).
    inventory.preferred = Inventory::kGrid + H().selected; // use comes out of the slot in hand first (#209)
    const std::string cur = CellId(H().selected);
    // The last block or machine of the slot in hand placed (#554): build mode ends with it, as the hand is empty now. (Build
    // mode entered with G on an empty hand is left alone: the slot was empty before too.)
    static int lastSlot = -1;
    static std::string lastId;
    if (buildTool.active && cur.empty() && lastSlot == H().selected && !lastId.empty()) buildTool.active = false;
    // (#1158) The same for a tool, bucket or ground in hand: when its stack leaves the selected slot (the last one placed
    // while more waits in the grid, or moved to the grid or a chest) or another item takes the slot, the hand follows the
    // slot - what is there now if it is held that way, else bare hands. (A tool a test put in hand with no slot is left be.)
    if (lastSlot == H().selected && cur != lastId && !lastId.empty() && terrainTool.heldTool == lastId) {
        const bool heldLikeThat = !cur.empty() && inventory.held(cur) > 0 &&
            ((KindOf(cur) == ItemKind::Tool && ToolOf(cur).type != ToolType::Hands) || cur == "bucket" || cur == "water_bucket" ||
             cur == "fishing_rod" || (TerrainTool_IsFiller(cur) && BuildMaterialForItem(cur.c_str()) < 0));
        terrainTool.heldTool = heldLikeThat ? cur : std::string();
    }
    lastSlot = H().selected; lastId = cur;
    // (#534) An item with no use from the hand (not a tool, filler, block or part) stays in hand too: just held.
    auto plainHeldId = [](const std::string& id) {
        return !id.empty() && ToolOf(id).type == ToolType::Hands && !TerrainTool_IsFiller(id) && BuildMaterialForItem(id.c_str()) < 0 && PartForItem(id) == PartType::Count &&
               id != "bucket" && id != "water_bucket" && id != "fishing_rod" && !Survival_IsUsable(id);
    };
    const bool inHand = cur.empty() || cur == terrainTool.heldTool || Survival_IsFood(cur) || Animals_IsFeed(cur) || Hotbar_IsLight(cur) || Hotbar_IsPlaceable(cur) || plainHeldId(cur);
    if (!buildTool.active && !inHand) {
        for (int i = 0; i < kHotbarSlots; ++i) if (CellId(i).empty()) { H().selected = i; break; }
        terrainTool.heldTool.clear();
    }
    // A tool lost (or never carried) falls back to bare hands.
    if (!terrainTool.heldTool.empty() && inventory.held(terrainTool.heldTool) <= 0) terrainTool.heldTool.clear();
    // (New tools land in an empty slot by themselves: Inventory::addAmount, #209.)
    if (in.enabled) {
        const int k = in.slot;
        if (k >= 0) Select(k, piloting);
        // The mouse wheel scrolls the hotbar (#227): down = the next slot, as most games do. With Ctrl or Alt held it is the
        // mode's own action (build shape, terraform brush); the telescope keeps it for its zoom.
        const float wheel = in.wheel;
        const bool mod = in.down(Act::Down) || in.down(Act::FreeLook); // Ctrl or Alt
        if (wheel != 0.0f && !mod && !piloting && !telescope.active) Select((H().selected + (wheel < 0.0f ? 1 : kHotbarSlots - 1)) % kHotbarSlots, piloting);
    }
    // Food or feed really in hand (not just selected: the last one eaten or fed leaves bare hands): RMB
    // never places terrain with it (#201). Derived each frame so it cannot go stale.
    const bool torchHeld = !buildTool.active && (Hotbar_IsLight(CellId(H().selected)) || Hotbar_IsPlaceable(CellId(H().selected))) && inventory.count(CellId(H().selected)) > 0;
    // (#534, the user) Any other item in hand that is not a tool, a filler or a bucket is just held: RMB does nothing with it
    // (unless it has a use: the telescope's RMB, Telescope_Update), LMB punches (bare hands).
    const std::string& selId = CellId(H().selected);
    const bool plainItem = !buildTool.active && inventory.count(selId) > 0 && plainHeldId(selId);
    terrainTool.holdingItem = !Hotbar_HeldFood().empty() || torchHeld || plainItem || terrainTool.heldTool == "fishing_rod"; // a torch in hand: RMB places it, not terrain (#218); a rod casts (#176)
    if (torchHeld && !piloting && in.hit(Act::Secondary)) Hotbar_PlaceTorch();
    // Eating (#202): hold right-click with food in hand for kEatSeconds (a bar under the crosshair); letting go cancels.
    UpdateEating(GameFrameTime(), !piloting && in.down(Act::Secondary));
}

bool Hotbar_TestEatHold(float seconds) { // tests: hold right-click this long (in 1/60 s steps); true if a bite was eaten
    const int before = H().eatenCount;
    for (float t = 0.0f; t < seconds; t += 1.0f / 60.0f) UpdateEating(1.0f / 60.0f, true);
    UpdateEating(1.0f / 60.0f, false);
    return H().eatenCount > before;
}

// (#827) The off hand on the HUD: its item (a torch lights, a shield blocks, food and oxygen are used with its key), left of the
// hotbar, with that key in the corner.
static void DrawOffHand(int x, int y, int size, float sc) {
    const InvSlot& s = inventory.offHand();
    DrawRectangle(x, y, size, size, Fade(BLACK, 0.45f));
    DrawRectangleLines(x, y, size, size, Fade(Color{ 224, 169, 64, 255 }, 0.6f));
    DrawText(PlayerInput_KeyName(PlayerInput_Binding(Act::OffHandUse)), x + (int)(4 * sc), y + (int)(3 * sc), (int)(12 * sc), Fade(Color{ 224, 169, 64, 255 }, 0.85f));
    if (s.empty()) { const char* t = "off hand"; const int fs = (int)(11 * sc); DrawText(t, x + (size - MeasureText(t, fs)) / 2, y + size - fs - (int)(6 * sc), fs, Fade(RAYWHITE, 0.45f)); return; }
    const float is = (float)size - 28.0f * sc;
    if (ItemIcons_Has(s.id)) ItemIcons_Draw(s.id, { x + (size - is) * 0.5f, (float)y + 8.0f * sc, is, is });
    if (KindOf(s.id) != ItemKind::Tool) { const char* n = TextFormat("%d", (int)s.amount()); DrawText(n, x + size - (int)(4 * sc) - MeasureText(n, (int)(14 * sc)), y + (int)(3 * sc), (int)(14 * sc), RAYWHITE); }
    std::string name = ShortName(s.id);
    if (s.hasData()) { const std::string cn = Composed_DisplayName(s.id, s.data); if (!cn.empty()) name = cn; } // (#1219)
    int fs = (int)(12 * sc);
    while (name.size() > 1 && MeasureText(name.c_str(), fs) > size - 4) name.pop_back();
    DrawText(name.c_str(), x + (size - MeasureText(name.c_str(), fs)) / 2, y + size - fs - (int)(6 * sc), fs, RAYWHITE);
}

void Hotbar_DrawHUD(bool piloting) {
    const float sc = GameUI_HudScale(); // (#827, #925) the HUD's size
    const int n = kHotbarSlots, size = (int)(64 * sc), gap = (int)(6 * sc);
    const int w = n * size + (n - 1) * gap;
    const int x0 = (GetScreenWidth() - w) / 2, y0 = GetScreenHeight() - size - (int)(44 * sc);
    DrawOffHand(x0 - size - (int)(18 * sc), y0, size, sc); // (#827) the off hand, left of the hotbar
    for (int i = 0; i < n; ++i) {
        const int x = x0 + i * (size + gap);
        const std::string& id = CellId(i);
        const int held = (int)inventory.grid()[Inventory::kGrid + i].amount(); // this slot's own stack
        const bool sel = !piloting && i == H().selected;
        DrawRectangle(x, y0, size, size, Fade(BLACK, sel ? 0.65f : 0.45f));
        DrawRectangleLines(x, y0, size, size, sel ? SKYBLUE : Fade(RAYWHITE, 0.35f));
        DrawText(TextFormat("%d", (i + 1) % 10), x + (int)(4 * sc), y0 + (int)(3 * sc), (int)(14 * sc), Fade(RAYWHITE, 0.7f));
        if (id.empty()) continue;
        {
            // A hot or radiant item glows (#580): the slot washed in its colour with a bright rim, as in the inventory.
            const ItemGlow g = Item_Glow(inventory.grid()[Inventory::kGrid + i]);
            if (g.strength > 0.02f) {
                DrawRectangle(x, y0, size, size, Fade(g.color, 0.55f * g.strength));
                DrawRectangleLinesEx({ (float)x, (float)y0, (float)size, (float)size }, 2.0f, Fade(g.color, 0.92f * g.strength));
            }
            const ItemDef* d = itemDB.find(id);
            const float is = (float)size - 28.0f * sc; // the painted icon (itemicons.h), between the slot number and the name
            if (ItemIcons_Has(id)) ItemIcons_Draw(id, { x + (size - is) * 0.5f, (float)y0 + 8.0f * sc, is, is });
            else if (d) DrawRectangle(x + size / 2 - 9, y0 + 14, 18, 18, d->color);
            if (KindOf(id) != ItemKind::Tool) DrawText(TextFormat("%d", held), x + size - 4 - MeasureText(TextFormat("%d", held), (int)(14 * sc)), y0 + (int)(3 * sc), (int)(14 * sc), RAYWHITE);
            // a worn tool: a bar along the bottom, green to red; (#1449) a composed copy: its own wear, from its record
            const InvSlot& cellNow = inventory.grid()[Inventory::kGrid + i];
            const ComposedStats own = cellNow.hasData() ? Composed_Stats(cellNow.id, cellNow.data) : ComposedStats{};
            const float left = own.valid && own.durability > 0 ? 1.0f - (float)own.wear / (float)own.durability : ToolWear_Left(id);
            if (left < 1.0f && held > 0) {
                DrawRectangle(x + 4, y0 + size - 4, size - 8, 3, Fade(BLACK, 0.6f));
                DrawRectangle(x + 4, y0 + size - 4, (int)((size - 8) * left), 3, UiPalette_Lerp(left));
            }
        }
        std::string name = ShortName(id);
        { const InvSlot& own = inventory.grid()[Inventory::kGrid + i]; const std::string cn = own.hasData() ? Composed_DisplayName(own.id, own.data) : std::string(); if (!cn.empty()) name = cn; } // (#1219) "Bronze sword", as the inventory names it
        // (#925 follow-up) As the inventory's names: never cut mid-word. One line, shrunk to 70% at most; else two lines split at
        // the space nearest the middle; a single long word shrinks to fit (raylib's font has a 10 px floor).
        const int fs0 = (int)(12 * sc), room = size - 4;
        int fs = fs0;
        while (fs > (int)(fs0 * 0.7f) && MeasureText(name.c_str(), fs) > room) --fs;
        size_t split = std::string::npos;
        if (MeasureText(name.c_str(), fs) > room)
            for (size_t i = 0; i < name.size(); ++i)
                if (name[i] == ' ' && (split == std::string::npos || std::abs((int)i - (int)name.size() / 2) < std::abs((int)split - (int)name.size() / 2))) split = i;
        if (split == std::string::npos) {
            while (fs > 10 && MeasureText(name.c_str(), fs) > room) --fs;
            DrawText(name.c_str(), x + (size - MeasureText(name.c_str(), fs)) / 2, y0 + size - fs - 6, fs, RAYWHITE);
        } else {
            const std::string a = name.substr(0, split), b = name.substr(split + 1);
            fs = fs0;
            while (fs > 10 && std::max(MeasureText(a.c_str(), fs), MeasureText(b.c_str(), fs)) > room) --fs;
            DrawText(a.c_str(), x + (size - MeasureText(a.c_str(), fs)) / 2, y0 + size - 2 * fs - 7, fs, RAYWHITE);
            DrawText(b.c_str(), x + (size - MeasureText(b.c_str(), fs)) / 2, y0 + size - fs - 5, fs, RAYWHITE);
        }
    }
    // Eating (#202): its bar under the crosshair, where the swing bar sits.
    if (H().eatProgress > 0.0f) {
        const int cx = GetScreenWidth() / 2, cy = GetScreenHeight() / 2, bw = 80;
        DrawRectangle(cx - bw / 2, cy + 22, bw, 6, Fade(BLACK, 0.5f));
        DrawRectangle(cx - bw / 2, cy + 22, (int)(bw * std::min(H().eatProgress, 1.0f)), 6, { 120, 200, 90, 255 });
        DrawText("eating", cx - MeasureText("eating", 14) / 2, cy + 32, 14, Fade(RAYWHITE, 0.8f));
    }
}

void Hotbar_Save(FILE*) {} // the slots are inventory cells now (their "slot 45-54" lines, #209)

bool Hotbar_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "hotbar") != 0) return false;
    int slot = -1;
    char id[64] = {};
    // An older save's binding (#209: before the hotbar held stacks): move a stack of the item into the slot (the inventory's
    // "slot" lines come first in the file).
    if (fscanf(f, " %d %63s", &slot, id) == 2 && slot >= 0 && slot < kHotbarSlots && strcmp(id, "-") != 0 && strcmp(id, kHotbarTerrainTool) != 0) {
        auto& g = inventory.grid();
        for (int i = 0; i < Inventory::kGrid; ++i) if (!g[i].empty() && g[i].id == id) { std::swap(g[i], g[Inventory::kGrid + slot]); break; }
    }
    return true;
}

void Hotbar_Unbind(const std::string&) {
    // A broken tool's slot is already empty (its cell held the stack, #209): the hand goes bare if that was the slot in hand.
    if (CellId(H().selected).empty()) { terrainTool.heldTool.clear(); terrainTool.holdingItem = false; }
}

bool Hotbar_PlaceTorch() {
    // A torch in hand (#218): RMB stands one where the crosshair meets the ground (torches never burn out, #1093; fibre
    // torches are too flimsy to stand). It goes into the cell just off the surface, on the nearest plot.
    const std::string id = CellId(H().selected);
    if (id == "fibre_torch") { GameUI_Toast("A fibre torch is too flimsy to stand: craft a torch"); Audio_Play(Sfx::Error, 0.3f); return false; }
    if ((id != "torch" && id != "snare") || inventory.count(id) <= 0) return false;
    Planetoid* body = nullptr; Vector3 local{}, n{};
    if (!TerrainTool_LastHit(body, local, n) || !body) return false;
    if (!body->plots) body->plots = std::make_unique<PlotSet>();
    Plot* pl = body->plots->plotFor(local, 6.0f, false); // a build within 6 m takes it (as the build tool does); else a plot of its own
    if (!pl) pl = body->plots->plotFor(local, 0.0f, true);
    if (!pl) return false;
    const Vector3 nl = body->dirToLocal(n);
    const Part part{ id == "snare" ? PartType::Snare : PartType::Torch, pl->cellAt(pl->toPlot(Vector3Add(local, Vector3Scale(nl, 0.15f)))), 0, 0 };
    if (!pl->canPlace(part)) { GameUI_Toast(TextFormat("No room for a %s there", id.c_str())); Audio_Play(Sfx::Error, 0.3f); return false; }
    if (!BuildTool_Place(*body, *pl, part, true)) return false;
    BuildTool_NotePlacedByPlayer(*body, *pl, part); // (#919) Building XP once it has stood
    Audio_Play(Sfx::Build, 0.5f, 1.3f);
    return true;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1158): the hand follows the selected slot ----
#include "unit_helpers.h"

// Dirt in the selected slot, more in the grid: the slot's stack placed out (or moved away) leaves bare hands, not dirt still
// in hand over an empty slot; a shovel moved into the slot is taken in hand; a tool a test set with no slot item stays.
TEST(Hotbar, TheHandFollowsTheSelectedSlotWhenItsStackLeaves) {
    UnitTest::EnsureItemDB();
    inventory.clear();
    Hotbar_Reset();
    Scene& world = UnitTest::FieldOnlyScene();
    Hotbar_Update(world, PlayerInput::None(), false); // (its first call resets the bindings: get that out of the way)
    const int slot = 2, cell = Inventory::kGrid + slot;
    inventory.grid()[cell] = InvSlot{ "dirt", 500 };
    inventory.grid()[0] = InvSlot{ "dirt", 1000 };     // more dirt in the grid
    Hotbar_Press(slot, false);
    Hotbar_Update(world, PlayerInput::None(), false);
    ASSERT_EQ(terrainTool.heldTool, "dirt");
    inventory.grid()[cell] = InvSlot{};                 // the last of the slot placed (or moved to the grid)
    Hotbar_Update(world, PlayerInput::None(), false);
    EXPECT_EQ(terrainTool.heldTool, "") << "an empty selected slot is bare hands, though dirt is still carried";
    inventory.grid()[cell] = InvSlot{ "stone_shovel", 100 }; // a shovel dragged into the selected slot
    Hotbar_Update(world, PlayerInput::None(), false);
    EXPECT_EQ(terrainTool.heldTool, "") << "(nothing was in hand to follow: the next selection takes it)";
    Hotbar_Press(slot, false);
    Hotbar_Update(world, PlayerInput::None(), false);
    ASSERT_EQ(terrainTool.heldTool, "stone_shovel");
    inventory.grid()[cell] = InvSlot{ "dirt", 300 };    // swapped for dirt in the inventory screen
    Hotbar_Update(world, PlayerInput::None(), false);
    EXPECT_EQ(terrainTool.heldTool, "dirt") << "the hand takes what the slot holds now";
    inventory.grid()[cell] = InvSlot{};
    Hotbar_Update(world, PlayerInput::None(), false);
    terrainTool.heldTool = "stone_pickaxe";            // a bench puts a tool in hand with no slot item
    inventory.grid()[1] = InvSlot{ "stone_pickaxe", 100 };
    Hotbar_Update(world, PlayerInput::None(), false);
    EXPECT_EQ(terrainTool.heldTool, "stone_pickaxe") << "a tool set with no slot change is left be";
    inventory.clear(); terrainTool.heldTool.clear();
}
#endif
