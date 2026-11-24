#ifndef BELTSIM_H
#define BELTSIM_H

// (#1735, thin belts PR 2) Items riding belts, Factorio-style: one lane, items spaced along the belt, moving at the belt's
// speed and backing up when the end is blocked. Pure integer logic for the factory thread (deterministic: the same lines,
// inserts and sink answers give the same state on every machine, which the MP sync's per-belt checksums rely on). No game
// code calls it yet (conversion of the old conveyors is PR 3).
//
// Units: positions along a belt are in 1/64 cell (8 per conduit slot), measured from the run's first point.

#include "conduit.h"

#include <cstdint>
#include <functional>
#include <vector>

constexpr int kBeltSub = 8;    // 1/64-cell units per conduit slot
constexpr int kBeltGap = 16;   // the minimum distance between two items' positions: a quarter cell, 4 items per cell

struct BeltItem {
    int item = 0;              // the caller's item id
    int64_t pos = 0;           // 0 .. line length
};

enum class BeltEnd : uint8_t { None, Sink, Belt };

struct BeltLine {
    int id = 0;                // the run id (checksums and saves key on it)
    int64_t length = 0;        // 1/64 cell: Conduit_Length(run) * kBeltSub
    int speed = 0;             // 1/64 cell per tick (the material's speed)
    std::vector<BeltItem> items; // front first: pos strictly decreasing, neighbours at least kBeltGap apart
    BeltEnd end = BeltEnd::None;
    int sink = -1;             // BeltEnd::Sink: the caller's sink id (a machine's item_in port)
    int outLine = -1;          // BeltEnd::Belt: the index of the line this one feeds
    int64_t outPos = 0;        //   and the position on it where the front item goes in (a side merge or 0 end to end)
    bool asleep = false;       // nothing moved and the front could not leave: skipped until woken
};

// A sink takes the front item (true) or refuses it (false: the belt waits). Called in a fixed order within a step.
typedef std::function<bool(int sink, int item)> BeltSinkFn;

// A line for a belt run (kind Belt) with the given speed; empty, no end.
BeltLine Belt_MakeLine(const ConduitRun& run, int speed);

// Puts an item onto the line at pos if no item is closer than kBeltGap to it (an output port loading the belt uses pos 0).
// Wakes the line. False when there is no room or pos is outside 0..length.
bool Belt_TryInsert(BeltLine& line, int item, int64_t pos = 0);
bool Belt_HasRoom(const BeltLine& line, int64_t pos);

// Wakes a line (its sink may take something again).
void Belt_Wake(BeltLine& line);

// One factory tick over all lines in index order. Each awake line moves its items front to back: the front item by speed
// up to the end, every other item by speed but no closer than kBeltGap behind the one ahead. An item at the end leaves
// into the line's sink or merge target when it has room. A line where nothing moved and the front stayed sleeps; a line
// whose items moved or left wakes the lines feeding it (their merge may have room now). Sinks wake their lines themselves
// (Belt_Wake) when they free room.
void Belt_Step(std::vector<BeltLine>& lines, const BeltSinkFn& sink);

// FNV-1a over the id and every item (id, pos): the MP sync compares these per belt.
uint32_t Belt_Checksum(const BeltLine& line);

// The point at pos along the run's centreline, in 1/64-cell units (slot * kBeltSub); clamped to the ends.
CInt3 Belt_PointAt(const ConduitRun& run, int64_t pos);
// The position along the run (1/64 cell) of a centreline slot, or -1 when the slot is not on the centreline.
int64_t Belt_PosOf(const ConduitRun& run, CInt3 slot);

#endif
