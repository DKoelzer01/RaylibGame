#ifndef BELTSYNC_H
#define BELTSYNC_H

// (#1735 PR 4) Belts in multiplayer: a joined client runs no factory, but it replays the server's belt sim exactly. A plot's
// belts are a snapshot (every line with its items, ends and sleep flags) plus, per belt step, the events the client can't
// work out itself: the lines a machine woke, the items loaders put on, and the lines whose front a machine took. Movement
// and side merges follow from the state (beltsim is deterministic). A checksum after the steps catches a drift; the client
// then asks for a new snapshot. The layout's epoch changes whenever the lines are rebuilt (the client needs a snapshot).
//
// Pure logic: the factory records, the net layer carries (item ids are the factory's; the wire uses names), the client's
// replica replays.

#include "beltsim.h"

#include <cstdint>
#include <deque>
#include <map>
#include <vector>

struct BeltSyncStep {
    std::vector<int> wakes;                   // lines a machine woke since the step before (Belt_Wake), before the loads
    std::vector<std::pair<int, int>> loads;   // (line, item) put on at 0, in loader order, before the step
    std::vector<int> takes;                   // lines whose front went into its machine this step, in line order
    bool operator==(const BeltSyncStep& o) const { return wakes == o.wakes && loads == o.loads && takes == o.takes; }
};

struct BeltSyncSnapshot {
    uint32_t epoch = 0;
    int64_t step = 0;                         // the steps done so far: the next one is step + 1
    std::vector<ConduitRun> runs;             // parallel to lines
    std::vector<BeltLine> lines;
};

// FNV-1a over every line's checksum, in line order.
uint32_t Belt_PlotChecksum(const std::vector<BeltLine>& lines);

// A client's copy of one plot's belts.
class BeltReplica {
public:
    static constexpr int kStepsPerSecond = 20;  // the factory's tick rate: one belt step a tick
    static constexpr int kPlayout = 4;          // steps kept in hand behind the newest received (smooth over jitter)

    // Takes the snapshot's lines (sinks become the line's own index: the takes say which line gives up its front).
    void reset(const BeltSyncSnapshot& s);
    // Steps first .. first + n - 1 and the checksum after the last. Ignored when it isn't the next run of steps of this
    // epoch (a gap or a new epoch: the replica then needs a snapshot).
    void push(uint32_t epoch, int64_t first, const std::vector<BeltSyncStep>& steps, uint32_t checksum);
    // Plays the queued steps at kStepsPerSecond, keeping about kPlayout in hand (catching up when far behind). False once
    // a step didn't replay (a load found no room, a take had no front at the end, a checksum differed): needResync.
    bool advance(float dt);
    bool playAll();                              // every queued step now (tests)

    bool valid() const { return valid_; }
    bool needResync() const { return resync_; }
    uint32_t epoch() const { return epoch_; }
    int64_t step() const { return step_; }       // the steps played
    int64_t received() const { return step_ + (int64_t)queue_.size(); }
    float fraction() const { return acc_; }      // 0..1 of the way to the next step (drawing interpolates)
    const std::vector<ConduitRun>& runs() const { return runs_; }
    const std::vector<BeltLine>& lines() const { return lines_; }

private:
    bool play(int n);
    bool valid_ = false, resync_ = false;
    uint32_t epoch_ = 0;
    int64_t step_ = 0;
    float acc_ = 0.0f;
    std::vector<ConduitRun> runs_;
    std::vector<BeltLine> lines_;
    std::deque<BeltSyncStep> queue_;
    std::map<int64_t, uint32_t> sums_;           // step -> the server's checksum after it
};

// Server side: one plot's recent steps (a ring the net layer reads from).
struct BeltSyncRing {
    static constexpr size_t kMax = 200;           // 10 s of steps: a client further behind gets a snapshot
    uint32_t epoch = 0;
    int64_t first = 1;                            // the step number of steps[0]
    std::deque<BeltSyncStep> steps;
    std::deque<uint32_t> sums;                    // the checksum after each step
    void add(uint32_t ep, int64_t step, BeltSyncStep s, uint32_t sum);
    // The steps after `after` (at most max): false when they aren't all here any more (or the epoch moved on).
    bool since(uint32_t ep, int64_t after, size_t max, std::vector<BeltSyncStep>& out, int64_t& firstOut, uint32_t& sumOut) const;
};

#endif
