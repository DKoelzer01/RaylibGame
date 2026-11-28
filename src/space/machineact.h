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
void MachineAct_Insert(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, bool fuel = false, int centi = 0); // centi (#1102): a fraction more (a local chest keeps it)
void MachineAct_Extract(uint64_t bodyId, uint32_t plotId, int part);
void MachineAct_SetFire(uint64_t bodyId, uint32_t plotId, int part, bool lit); // (#1103) light or put out a campfire. Joined: action 10.
// (#862) As Insert, into one storage slot of a chest (a machine routes it as Insert does). Joined: action 8.
void MachineAct_InsertSlot(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count, int slot, int centi = 0);
// (#862) `count` of item out of a chest into the inventory (the scroll wheel's pull). Joined: action 9.
void MachineAct_Take(uint64_t bodyId, uint32_t plotId, int part, const std::string& item, int count);
// One slot of a machine's screen: to the inventory, or (cursorCount >= 0) onto the cursor (0: all). Joined: always to the
// inventory for now.
void MachineAct_TakeSlot(uint64_t bodyId, uint32_t plotId, int part, int slot, int which, int cursorCount = -1);
bool MachineAct_Remote();                         // a joined game: actions go to the server
void MachineAct_Update();                        // once a frame (#854: a joined client's cursor takes, once their items arrive)
// (#866) `count` of item is coming into the inventory from the server (a creative take): move it onto the cursor once its
// delta has brought it (as the cursor takes above).
void MachineAct_ExpectOnCursor(const std::string& item, int count);
int MachineAct_CursorMoves();                    // (#854 tests) cursor takes completed (the items moved onto the cursor)

// (#1684) MsgMachineAction codes for a dungeon's parts (not machines): the host unlocks a vault door (taking a dungeon_key
// from the requester) or pulls a lever, and sends the part changes to everyone.
constexpr int kMachineActUnlockVault = 11, kMachineActPullLever = 12;
// (#2072) A ruin cache, its hold done: the host opens it (spending the requester's power cell, as its kind needs), rolls,
// gives the parts, marks it opened for everyone and sends the requester the result lines.
constexpr int kMachineActOpenCache = 13;

// Set while joined (main.cpp): body is the local Planetoid id; action as MsgMachineAction.
inline std::function<void(uint64_t bodyId, uint32_t plotId, int part, int action, const std::string& item, int count, int slot, int which)> gMachineRequest;
