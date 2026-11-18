#pragma once

// A player's actions on a machine (#332): in single player straight into the factory, as before; in a joined game a
// request to the server (gMachineRequest, NetMsg MachineAction), which owns the machines and the player's inventory and
// answers with inventory deltas. The UI calls these instead of Factory::Insert / Extract / TakeSlot / ... for what the
// player does (world generation and tests still call the factory directly).

#include <cstdint>
#include <functional>
#include <string>

// The caller has already taken `count` of `item` out of the inventory (or off the cursor), as before. Joined: it is
// given back here and the server takes it out of the inventory it keeps (the delta then removes it here too).
void MachineAct_Insert(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, bool fuel = false);
void MachineAct_Extract(uint64_t bodyId, uint32_t plotId, int part);
// One slot of a machine's screen: to the inventory, or (cursorCount >= 0) onto the cursor (0: all). Joined: always to the
// inventory for now.
void MachineAct_TakeSlot(uint64_t bodyId, uint32_t plotId, int part, int slot, int which, int cursorCount = -1);
bool MachineAct_Remote();                         // a joined game: actions go to the server
void MachineAct_Update();                        // once a frame (#854: a joined client's cursor takes, once their items arrive)
int MachineAct_CursorMoves();                    // (#854 tests) cursor takes completed (the items moved onto the cursor)

// Set while joined (main.cpp): body is the local Planetoid id; action as MsgMachineAction.
inline std::function<void(uint64_t bodyId, uint32_t plotId, int part, int action, const std::string& item, int count, int slot, int which)> gMachineRequest;
