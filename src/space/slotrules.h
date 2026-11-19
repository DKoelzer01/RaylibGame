#pragma once

// Item movement rules for every slot grid (#862, the user: "Minecraft's logic + the Mouse Tweaks mod's logic"): the player's
// cells, a chest's or machine's slots and the creative Catalog all go through these, so a click means the same everywhere.
// The rules hold no UI: a screen reports the mouse and keys each frame (SlotInput) with the slot under the mouse, and each
// grid says what a put, take or move does there (a chest sends requests, the Catalog gives for free). Tests drive them with
// made-up input.
//
// Minecraft:
//   LMB on a stack: pick it all up; LMB holding a stack: put it all down / merge / swap with another item.
//   RMB on a stack: pick up half (rounded up); RMB holding a stack: put one down (another item: swap).
//   Shift+LMB: quick-move the stack to the other side (the open container, else grid <-> hotbar).
//   Double-click LMB: gather every stack of that item onto the cursor (up to a full stack).
//   LMB drag holding a stack: split it evenly over the slots dragged across (on release); RMB drag: one into each.
//   1-0 over a slot: swap it with that hotbar cell; the off-hand key: with the off hand. Q: drop one; Ctrl+Q: the stack.
//   Middle-click in creative: a full stack of that item onto the cursor.
// Mouse Tweaks:
//   RMB drag: one per slot (as above). LMB drag after picking a stack up: pick up every stack of that item dragged across.
//   Shift+LMB drag: quick-move every slot dragged across. Scroll wheel over a slot: down moves one (Shift: the stack) to the
//   other side, up pulls one of that item in from it (the direction can be swapped in the settings).

#include "inventory.h"
#include <cstdint>
#include <string>
#include <vector>

struct SlotGrid {
    virtual ~SlotGrid() = default;
    virtual int size() const = 0;
    virtual InvSlot at(int i) const = 0;                              // what is in slot i now (a container: its last view)
    virtual int64_t room(int i, const std::string& id) const = 0;     // centi of id that fits there now (0: none, or not here)
    // Up to `centi` of the held stack into slot i; returns what went (gone from the hand at once: a container sends it).
    virtual int64_t put(int i, const InvSlot& held, int64_t centi) = 0;
    // Up to `centi` of slot i onto the held stack (empty, or the same item); returns what arrived now (a container's comes on
    // the factory's next tick through the cursor sink, and returns 0 here).
    virtual int64_t take(int i, int64_t centi, InvSlot& held) = 0;
    virtual void moveOut(int i, int64_t centi) = 0;                   // to the other side (shift-click, the wheel down)
    virtual void pullIn(int i, int64_t centi) = 0;                    // that item from the other side into i (the wheel up)
    virtual void swapCell(int i, int cell) {}                          // with the player's inventory cell (1-0, the off hand)
    virtual void drop(int i, int64_t centi) {}                         // on the ground at the feet (Q)
    virtual bool infinite() const { return false; }                   // the Catalog: takes are free, a put deletes the stack
    virtual int64_t capOf(const std::string& id) const { return (int64_t)StackLimit(id) * 100; } // a full stack (centi)
};

struct SlotInput {
    bool lPress = false, rPress = false, mPress = false;   // went down this frame
    bool lDown = false, rDown = false;                     // held
    bool lRelease = false, rRelease = false;               // came up this frame
    bool lDouble = false;                                  // the press was a double-click
    float wheel = 0.0f;                                    // + up, - down (notches)
    bool shift = false, ctrl = false;
    int hotbarKey = -1;                                    // 0-9: the hotbar key pressed this frame
    bool offHandKey = false, dropKey = false;
    bool creative = false;
    bool wheelInvert = false;                              // the setting: down pulls in, up moves out
};

// The cursor stack the rules move (the inventory screen's).
void SlotRules_SetHeld(InvSlot* held);
// The grids on screen this frame, in gathering order (the player's first): a double-click gathers from them.
void SlotRules_SetGrids(const std::vector<SlotGrid*>& grids);
// One frame: the input and the slot under the mouse (g null: none). Returns true when it did something.
bool SlotRules_Frame(const SlotInput& in, SlotGrid* g, int i);
void SlotRules_Cancel();                                   // the window closed: a drag in progress ends without effect
// A drag in progress, for drawing its preview: the slots dragged across and what each would get (centi; a split's share).
struct SlotDragPreview { const SlotGrid* grid; int slot; int64_t centi; };
std::vector<SlotDragPreview> SlotRules_DragPreview();
