#include "beltsync.h"

#include <algorithm>

uint32_t Belt_PlotChecksum(const std::vector<BeltLine>& lines) {
    uint32_t h = 2166136261u;
    for (const BeltLine& l : lines) {
        uint32_t v = Belt_Checksum(l);
        for (int b = 0; b < 4; ++b) { h ^= v & 0xff; h *= 16777619u; v >>= 8; }
    }
    return h;
}

void BeltReplica::reset(const BeltSyncSnapshot& s) {
    valid_ = true; resync_ = false;
    epoch_ = s.epoch; step_ = s.step; acc_ = 0.0f;
    runs_ = s.runs; lines_ = s.lines;
    for (size_t i = 0; i < lines_.size(); ++i) {
        BeltLine& l = lines_[i];
        if (l.end == BeltEnd::Sink) l.sink = (int)i;
        l.packed = 0; // a cache: 0 is always right (beltsim recomputes it)
        int64_t pos = l.length;
        for (const BeltItem& it : l.items) pos -= it.gap;
        l.tailPos = pos;
    }
    queue_.clear(); sums_.clear();
}

void BeltReplica::push(uint32_t epoch, int64_t first, const std::vector<BeltSyncStep>& steps, uint32_t checksum) {
    if (!valid_ || steps.empty()) return;
    if (epoch != epoch_ || first != received() + 1) {
        if (epoch != epoch_ || first > received() + 1) resync_ = true; // a new layout, or steps missed
        return;                                                         // (an old repeat: already have it)
    }
    queue_.insert(queue_.end(), steps.begin(), steps.end());
    sums_[received()] = checksum;
}

bool BeltReplica::play(int n) {
    for (int k = 0; k < n && !queue_.empty() && !resync_; ++k) {
        const BeltSyncStep s = std::move(queue_.front());
        queue_.pop_front();
        const int nl = (int)lines_.size();
        for (int w : s.wakes) { if (w < 0 || w >= nl) { resync_ = true; break; } Belt_Wake(lines_[w]); }
        for (const auto& [li, item] : s.loads)
            if (resync_ || li < 0 || li >= nl || !Belt_TryInsert(lines_[li], item, 0)) { resync_ = true; break; }
        if (resync_) break;
        size_t next = 0; // the takes come in line order, as Belt_Step asks the sinks
        bool order = true;
        std::vector<int> taken;
        Belt_Step(lines_, [&](int line, int) {
            if (next < s.takes.size() && s.takes[next] == line) { ++next; return true; }
            for (size_t j = next; j < s.takes.size(); ++j) if (s.takes[j] == line) order = false; // asked out of order
            return false;
        }, &taken);
        ++step_;
        if (!order || taken != s.takes) { resync_ = true; break; }
        const auto sum = sums_.find(step_);
        if (sum != sums_.end()) {
            if (sum->second != Belt_PlotChecksum(lines_)) resync_ = true;
            sums_.erase(sums_.begin(), std::next(sum));
        }
    }
    return !resync_;
}

bool BeltReplica::advance(float dt) {
    if (!valid_ || resync_) return !resync_;
    acc_ += std::max(0.0f, dt) * kStepsPerSecond;
    int n = (int)acc_;
    acc_ -= (float)n;
    const int64_t ahead = (int64_t)queue_.size();
    if (ahead - n > 2 * kPlayout) n = (int)(ahead - kPlayout); // far behind (a burst, a stall): catch up to the playout
    if (n >= ahead) { n = (int)ahead; if (ahead == 0) acc_ = 0.0f; } // starved: wait at the last step (no drift ahead)
    return play(n);
}

bool BeltReplica::playAll() { return play((int)queue_.size()); }

void BeltSyncRing::add(uint32_t ep, int64_t step, BeltSyncStep s, uint32_t sum) {
    if (ep != epoch || step != first + (int64_t)steps.size()) { epoch = ep; first = step; steps.clear(); sums.clear(); }
    steps.push_back(std::move(s));
    sums.push_back(sum);
    while (steps.size() > kMax) { steps.pop_front(); sums.pop_front(); ++first; }
}

bool BeltSyncRing::since(uint32_t ep, int64_t after, size_t max, std::vector<BeltSyncStep>& out, int64_t& firstOut, uint32_t& sumOut) const {
    out.clear();
    if (ep != epoch || after + 1 < first) return false;
    const int64_t last = first + (int64_t)steps.size() - 1;
    if (after > last) return false;      // ahead of the ring: can't be (a reset)
    firstOut = after + 1;
    for (int64_t s = after + 1; s <= last && out.size() < max; ++s) out.push_back(steps[(size_t)(s - first)]);
    if (!out.empty()) sumOut = sums[(size_t)(firstOut + (int64_t)out.size() - 1 - first)];
    return true;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
#include <random>

namespace {
// A server plot: random lines (some ending in a sink, some merging into a lower or higher line), random loaders, sinks that
// take at random, and random machine wakes; it records each step's events as the factory does.
struct Server {
    std::vector<ConduitRun> runs;
    std::vector<BeltLine> lines;
    std::mt19937 rng;
    int64_t step = 0;
    explicit Server(unsigned seed) : rng(seed) {
        const int n = 3 + (int)(rng() % 8);
        for (int i = 0; i < n; ++i) {
            ConduitRun r; r.id = i + 1; r.kind = ConduitKind::Belt; r.material = "basic"; r.size = kBeltWidth;
            const int len = 2 * (4 + (int)(rng() % 40));
            r.points = { { 0, 0, 2 * i }, { len, 0, 2 * i } };
            BeltLine l = Belt_MakeLine(r, 2 + (int)(rng() % 10));
            if (i > 0 && rng() % 3 == 0) { l.end = BeltEnd::Belt; l.outLine = (int)(rng() % n); if (l.outLine == i) l.outLine = 0; l.outPos = 8 * (int64_t)(rng() % 8); }
            else if (rng() % 4 != 0) { l.end = BeltEnd::Sink; l.sink = 100 + i; }
            runs.push_back(r);
            lines.push_back(l);
        }
    }
    BeltSyncSnapshot snapshot() const { BeltSyncSnapshot s; s.epoch = 1; s.step = step; s.runs = runs; s.lines = lines; return s; }
    BeltSyncStep tick() {
        BeltSyncStep s;
        for (size_t i = 0; i < lines.size(); ++i) if (rng() % 10 == 0) { Belt_Wake(lines[i]); s.wakes.push_back((int)i); }
        for (size_t i = 0; i < lines.size(); ++i)
            if (rng() % 3 == 0 && Belt_HasRoom(lines[i], 0)) { const int item = 1 + (int)(rng() % 5); Belt_TryInsert(lines[i], item, 0); s.loads.push_back({ (int)i, item }); }
        Belt_Step(lines, [&](int, int) { return rng() % 3 != 0; }, &s.takes);
        ++step;
        return s;
    }
};
} // namespace

// The replica fed only the snapshot and the events stays identical to the server, every step, sleep flags too; whether the
// steps arrive one at a time or in bundles.
TEST(BeltSync, AReplicaReplaysTheServerExactly) {
    for (unsigned seed = 1; seed <= 30; ++seed) {
        Server srv(seed);
        for (int t = 0; t < 50; ++t) srv.tick(); // join a running plot
        BeltReplica rep;
        rep.reset(srv.snapshot());
        const int bundle = 1 + (int)(seed % 4);
        std::vector<BeltSyncStep> pending;
        int64_t first = srv.step + 1;
        for (int t = 0; t < 1500; ++t) {
            pending.push_back(srv.tick());
            if ((int)pending.size() < bundle) continue;
            rep.push(1, first, pending, Belt_PlotChecksum(srv.lines));
            pending.clear();
            first = srv.step + 1;
            ASSERT_TRUE(rep.playAll()) << "seed " << seed << " step " << srv.step;
            ASSERT_EQ(rep.step(), srv.step);
            ASSERT_EQ(Belt_PlotChecksum(rep.lines()), Belt_PlotChecksum(srv.lines)) << "seed " << seed << " step " << srv.step;
            for (size_t i = 0; i < srv.lines.size(); ++i) ASSERT_EQ(rep.lines()[i].asleep, srv.lines[i].asleep) << "seed " << seed << " line " << i;
        }
    }
}

// A drift (a missed take), a missed run of steps and a new epoch each ask for a snapshot; a snapshot clears it.
TEST(BeltSync, DriftGapsAndNewLayoutsAskForASnapshot) {
    Server srv(7);
    BeltReplica rep;
    rep.reset(srv.snapshot());
    BeltSyncStep s;
    int guard = 0;
    do { s = srv.tick(); } while (s.takes.empty() && ++guard < 1000);
    ASSERT_FALSE(s.takes.empty());
    // replay everything up to the step with a take, but drop that take
    BeltReplica bad;
    { Server again(7); bad.reset(again.snapshot()); }
    {
        Server again(7);
        for (int64_t k = 1; k < srv.step; ++k) { BeltSyncStep e = again.tick(); bad.push(1, k, { e }, Belt_PlotChecksum(again.lines)); }
        BeltSyncStep e = again.tick();
        e.takes.clear();
        bad.push(1, srv.step, { e }, Belt_PlotChecksum(again.lines));
    }
    EXPECT_FALSE(bad.playAll());
    EXPECT_TRUE(bad.needResync()) << "a missed take";
    // a gap
    const BeltSyncStep later = srv.tick();
    rep.push(1, 5, { later }, 0);
    EXPECT_TRUE(rep.needResync()) << "steps 1 .. 4 missing";
    rep.reset(srv.snapshot());
    EXPECT_FALSE(rep.needResync());
    const BeltSyncStep e = srv.tick();
    rep.push(2, srv.step, { e }, Belt_PlotChecksum(srv.lines));
    EXPECT_TRUE(rep.needResync()) << "a new epoch";
}

// Playout: at 20 steps a second it keeps a few steps in hand, never plays steps it hasn't got, and catches up after a burst.
TEST(BeltSync, PlaysAtTheTickRateBehindTheNewest) {
    Server srv(3);
    BeltReplica rep;
    rep.reset(srv.snapshot());
    auto next = [&] { const BeltSyncStep e = srv.tick(); rep.push(1, srv.step, { e }, Belt_PlotChecksum(srv.lines)); }; // (tick first: argument order is unspecified)
    for (int k = 0; k < 6; ++k) next();
    EXPECT_TRUE(rep.advance(0.1f));
    EXPECT_EQ(rep.step(), 2) << "0.1 s: two steps";
    EXPECT_TRUE(rep.advance(1.0f));
    EXPECT_EQ(rep.step(), 6) << "never past the newest";
    for (int k = 0; k < 40; ++k) next();
    EXPECT_TRUE(rep.advance(0.05f));
    EXPECT_EQ(rep.received() - rep.step(), BeltReplica::kPlayout) << "a burst: caught up to the playout";
}

TEST(BeltSync, TheRingGivesTheStepsSinceAClientsLast) {
    BeltSyncRing ring;
    for (int64_t s = 1; s <= 250; ++s) { BeltSyncStep e; e.takes = { (int)s }; ring.add(3, s, e, (uint32_t)s * 7); }
    std::vector<BeltSyncStep> out; int64_t first = 0; uint32_t sum = 0;
    EXPECT_FALSE(ring.since(3, 10, 100, out, first, sum)) << "step 11 fell out of the ring";
    ASSERT_TRUE(ring.since(3, 240, 100, out, first, sum));
    EXPECT_EQ(first, 241);
    ASSERT_EQ(out.size(), 10u);
    EXPECT_EQ(out[0].takes[0], 241);
    EXPECT_EQ(sum, 250u * 7);
    ASSERT_TRUE(ring.since(3, 250, 100, out, first, sum));
    EXPECT_TRUE(out.empty()) << "up to date";
    EXPECT_FALSE(ring.since(4, 240, 100, out, first, sum)) << "another epoch";
}
#endif
