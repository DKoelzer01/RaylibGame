#include "machineact.h"
#include "factory.h"
#include "inventory.h"
#include "playersession.h"

bool MachineAct_Remote() { return (bool)gMachineRequest; }

void MachineAct_Insert(uint64_t body, uint32_t plot, int part, const std::string& item, int count, bool fuel) {
    if (count <= 0) return;
    if (gMachineRequest) {
        inventory.addAmount(item, count); // back: the server takes it from the inventory it keeps, and its delta from ours
        gMachineRequest(body, plot, part, fuel ? 3 : 1, item, count, 0, 0);
        return;
    }
    if (fuel) Factory::InsertFuel(body, plot, part, item, count);
    else Factory::Insert(body, plot, part, item, count);
}

void MachineAct_Extract(uint64_t body, uint32_t plot, int part) {
    if (gMachineRequest) { gMachineRequest(body, plot, part, 4, std::string(), 0, 0, 0); return; }
    Factory::Extract(body, plot, part);
}

void MachineAct_TakeSlot(uint64_t body, uint32_t plot, int part, int slot, int which, int cursorCount) {
    if (gMachineRequest) { gMachineRequest(body, plot, part, 5, std::string(), 0, slot, which); return; }
    if (cursorCount >= 0) Factory::TakeSlotToCursor(body, plot, part, slot, which, cursorCount);
    else Factory::TakeSlot(body, plot, part, slot, which);
}
