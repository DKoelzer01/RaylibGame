#include "machineact.h"
#include "factory.h"
#include "inventory.h"
#include "playersession.h"
#include <raylib.h>
#include <algorithm>
#include <vector>
#include <cstdio>

bool MachineAct_Remote() { return (bool)gMachineRequest; }

void MachineAct_Insert(uint64_t body, uint32_t plot, int part, const std::string& item, int count, bool fuel, int centi) {
    if (count <= 0 && centi <= 0) return;
    if (gMachineRequest) {
        inventory.addAmount(item, count); // back: the server takes it from the inventory it keeps, and its delta from ours
        gMachineRequest(body, plot, part, fuel ? 3 : 1, item, count, 0, 0);
        return;
    }
    if (fuel) Factory::InsertFuel(body, plot, part, item, count);
    else Factory::Insert(body, plot, part, item, count, 0, centi);
}

void MachineAct_InsertSlot(uint64_t body, uint32_t plot, int part, const std::string& item, int count, int slot, int centi) {
    if (count <= 0 && centi <= 0) return;
    if (gMachineRequest) {
        inventory.addAmount(item, count); // back, as Insert: the server takes it from the inventory it keeps
        gMachineRequest(body, plot, part, 8, item, count, slot, 0);
        return;
    }
    Factory::InsertSlot(body, plot, part, item, count, slot, 0, centi);
}

void MachineAct_Take(uint64_t body, uint32_t plot, int part, const std::string& item, int count) {
    if (count <= 0) return;
    if (gMachineRequest) { gMachineRequest(body, plot, part, 9, item, count, 0, 0); return; }
    Factory::Take(body, plot, part, item, count);
}

void MachineAct_Extract(uint64_t body, uint32_t plot, int part) {
    if (gMachineRequest) { gMachineRequest(body, plot, part, 4, std::string(), 0, 0, 0); return; }
    Factory::Extract(body, plot, part);
}
void MachineAct_SetFire(uint64_t body, uint32_t plot, int part, bool lit) {
    if (gMachineRequest) { gMachineRequest(body, plot, part, 10, std::string(), lit ? 1 : 0, 0, 0); return; }
    Factory::SetFire(body, plot, part, lit);
}
void MachineAct_SetMode(uint64_t body, uint32_t plot, int part, int mode) { // (#1271 L6c)
    if (gMachineRequest) { gMachineRequest(body, plot, part, 11, std::string(), mode, 0, 0); return; }
    Factory::SetMode(body, plot, part, mode);
}

// (#854) a joined client's take onto the cursor: the server puts it in the inventory it keeps; when its delta has brought
// it here, it moves from the inventory onto the cursor (the cursor is this screen's: the server never knew it).
namespace {
    struct CursorTake { std::string item; int count = 0, base = 0; double until = 0.0; };
    std::vector<CursorTake> cursorTakes;
    int cursorMoves = 0;
}

void MachineAct_Update() {
    for (size_t i = 0; i < cursorTakes.size();) {
        CursorTake& t = cursorTakes[i];
        auto sink = Factory::CursorSink();
        if (inventory.count(t.item) >= t.base + t.count) {
            if (sink && inventory.removeAmount(t.item, t.count)) { sink(t.item, t.count, 0.0f); ++cursorMoves; }
            cursorTakes.erase(cursorTakes.begin() + i);
            continue;
        }
        if (GetTime() > t.until) { cursorTakes.erase(cursorTakes.begin() + i); continue; } // (it stays in the inventory)
        ++i;
    }
}

void MachineAct_TakeSlot(uint64_t body, uint32_t plot, int part, int slot, int which, int cursorCount) {
    if (gMachineRequest) {
        if (cursorCount >= 0) { // onto the cursor: that much of the slot as we see it (the server's view)
            const MachineView v = Factory::View(body, plot, part);
            const auto& slots = which == 1 ? v.outSlots : v.inSlots;
            const std::pair<std::string, int> s = which == 2 ? v.fuelSlot : (slot >= 0 && slot < (int)slots.size() ? slots[slot] : std::pair<std::string, int>{});
            const int n = cursorCount > 0 ? std::min(cursorCount, s.second) : s.second;
            if (s.first.empty() || n <= 0) return;
            cursorTakes.push_back({ s.first, n, inventory.count(s.first), GetTime() + 3.0 });
            gMachineRequest(body, plot, part, 5, std::string(), n, slot, which);
            return;
        }
        gMachineRequest(body, plot, part, 5, std::string(), 0, slot, which);
        return;
    }
    if (cursorCount >= 0) Factory::TakeSlotToCursor(body, plot, part, slot, which, cursorCount);
    else Factory::TakeSlot(body, plot, part, slot, which);
}

void MachineAct_ExpectOnCursor(const std::string& item, int count) {
    if (count > 0 && !item.empty()) cursorTakes.push_back({ item, count, inventory.count(item), GetTime() + 3.0 });
}

int MachineAct_CursorMoves() { return cursorMoves; }
