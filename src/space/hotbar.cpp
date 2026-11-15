#include "hotbar.h"
#include "playersession.h"
#include "playerinput.h"
#include "world.h" // GameFrameTime
#include "toolwear.h"
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
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); return; }
            H().selected = slot; buildTool.active = false; terrainTool.heldTool.clear(); terrainTool.holdingItem = true; // RMB: no terrain (#201)
            GameUI_Toast(Survival_IsFood(id) ? TextFormat("Holding %s: hold right-click to eat", ShortName(id)) : TextFormat("Holding %s: animals that eat it follow you", ShortName(id)));
            return;
        }
        if ((Hotbar_IsLight(id) || Hotbar_IsPlaceable(id)) && !piloting) { // a torch (#218): held for its light; a snare (#176): set down with RMB; LMB still digs with bare hands
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); return; }
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
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
            take();
            buildTool.active = false;
            terrainTool.heldTool = id;
            return;
        }
        if (KindOf(id) == ItemKind::Tool && ToolOf(id).type != ToolType::Hands) {
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
            take();
            buildTool.active = false;
            terrainTool.heldTool = id;
            return;
        }
        if (TerrainTool_IsFiller(id) && BuildMaterialForItem(id.c_str()) < 0) { // loose ground (gravel, sand): in hand, RMB places it
            if (inventory.count(id) <= 0) { GameUI_Toast(TextFormat("No %s", ShortName(id))); Audio_Play(Sfx::Error, 0.3f); return; }
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
        const PartType t = PartForItem(id);
        if (t != PartType::Count) {
            take();
            buildTool.active = true;
            buildTool.part = (int)t;
            return;
        }
        GameUI_Toast(TextFormat("%s can't be used from the hotbar", ShortName(id)));
        Audio_Play(Sfx::Error, 0.3f);
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
    const bool inHand = cur.empty() || cur == terrainTool.heldTool || Survival_IsFood(cur) || Animals_IsFeed(cur) || Hotbar_IsLight(cur) || Hotbar_IsPlaceable(cur);
    if (!buildTool.active && !inHand) {
        for (int i = 0; i < kHotbarSlots; ++i) if (CellId(i).empty()) { H().selected = i; break; }
        terrainTool.heldTool.clear();
    }
    // A tool lost (or never carried) falls back to bare hands.
    if (!terrainTool.heldTool.empty() && inventory.count(terrainTool.heldTool) <= 0) terrainTool.heldTool.clear();
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
    terrainTool.holdingItem = !Hotbar_HeldFood().empty() || torchHeld || terrainTool.heldTool == "fishing_rod"; // a torch in hand: RMB places it, not terrain (#218); a rod casts (#176)
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

void Hotbar_DrawHUD(bool piloting) {
    const int n = kHotbarSlots, size = 64, gap = 6;
    const int w = n * size + (n - 1) * gap;
    const int x0 = (GetScreenWidth() - w) / 2, y0 = GetScreenHeight() - size - 44;
    for (int i = 0; i < n; ++i) {
        const int x = x0 + i * (size + gap);
        const std::string& id = CellId(i);
        const int held = (int)inventory.grid()[Inventory::kGrid + i].amount(); // this slot's own stack
        const bool sel = !piloting && i == H().selected;
        DrawRectangle(x, y0, size, size, Fade(BLACK, sel ? 0.65f : 0.45f));
        DrawRectangleLines(x, y0, size, size, sel ? SKYBLUE : Fade(RAYWHITE, 0.35f));
        DrawText(TextFormat("%d", (i + 1) % 10), x + 4, y0 + 3, 14, Fade(RAYWHITE, 0.7f));
        if (id.empty()) continue;
        {
            const ItemDef* d = itemDB.find(id);
            const float is = (float)size - 28.0f; // the painted icon (itemicons.h), between the slot number and the name
            if (ItemIcons_Has(id)) ItemIcons_Draw(id, { x + (size - is) * 0.5f, (float)y0 + 8.0f, is, is });
            else if (d) DrawRectangle(x + size / 2 - 9, y0 + 14, 18, 18, d->color);
            if (KindOf(id) != ItemKind::Tool) DrawText(TextFormat("%d", held), x + size - 4 - MeasureText(TextFormat("%d", held), 14), y0 + 3, 14, RAYWHITE);
            const float left = ToolWear_Left(id); // a worn tool: a bar along the bottom, green to red
            if (left < 1.0f && held > 0) {
                DrawRectangle(x + 4, y0 + size - 4, size - 8, 3, Fade(BLACK, 0.6f));
                DrawRectangle(x + 4, y0 + size - 4, (int)((size - 8) * left), 3, ColorLerp(RED, GREEN, left));
            }
        }
        std::string name = ShortName(id);
        int fs = 12;
        while (fs > 9 && MeasureText(name.c_str(), fs) > size - 4) --fs;
        while (name.size() > 1 && MeasureText(name.c_str(), fs) > size - 4) name.pop_back(); // still too long: cut (#219: names ran into the next slot)
        DrawText(name.c_str(), x + (size - MeasureText(name.c_str(), fs)) / 2, y0 + size - fs - 6, fs, RAYWHITE);
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
    // A torch in hand (#218): RMB stands one where the crosshair meets the ground (a placed torch never burns out; fibre
    // torches burn too fast to set down). It goes into the cell just off the surface, on the nearest plot.
    const std::string id = CellId(H().selected);
    if (id == "fibre_torch") { GameUI_Toast("Fibre torches burn too fast to set down: craft a torch"); Audio_Play(Sfx::Error, 0.3f); return false; }
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
    Audio_Play(Sfx::Build, 0.5f, 1.3f);
    return true;
}
