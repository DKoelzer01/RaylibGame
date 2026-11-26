#include "beltsim.h"

#include <algorithm>
#include <cstdlib>

// Items are kept front first as gaps (see beltsim.h). In position terms a tick is: the front item moves by speed up to the
// end; each item behind moves by speed but stops kBeltGap behind the item ahead's new position. With d the move of the item
// ahead and s = gap - kBeltGap an item's slack, it moves min(speed, s + d), and its gap changes by d minus its own move. Once
// an item moves the full speed every item behind it does too (s >= 0), and their gaps don't change: the walk stops there.

BeltLine Belt_MakeLine(const ConduitRun& run, int speed) {
    BeltLine l;
    l.id = run.id;
    l.length = Conduit_Length(run) * kBeltSub;
    l.speed = speed;
    return l;
}

std::vector<int64_t> Belt_Positions(const BeltLine& line) {
    std::vector<int64_t> p;
    p.reserve(line.items.size());
    int64_t at = line.length;
    for (const BeltItem& it : line.items) { at -= it.gap; p.push_back(at); }
    return p;
}

namespace {
// Where an item at pos would go (the index of the first item behind it), or -1 when an item is closer than kBeltGap.
// *ahead: the position of the item ahead of it (the line's length when it would be the front). Behind the tail (the
// loader at 0, the usual case) it's O(1): the tail is the last and lowest position.
int InsertIndex(const BeltLine& line, int64_t pos, int64_t* ahead) {
    if (pos < 0 || pos > line.length) return -1;
    *ahead = line.length;
    if (line.items.empty()) return 0;
    if (pos <= line.tailPos) {
        if (pos > line.tailPos - kBeltGap) return -1;
        *ahead = line.tailPos;
        return (int)line.items.size();
    }
    int64_t at = line.length;
    for (size_t k = 0; k < line.items.size(); ++k) {
        at -= line.items[k].gap;
        if (std::llabs(at - pos) < kBeltGap) return -1;
        if (at < pos) return (int)k;
        *ahead = at;
    }
    return (int)line.items.size();
}
} // namespace

bool Belt_HasRoom(const BeltLine& line, int64_t pos) { int64_t a; return InsertIndex(line, pos, &a) >= 0; }

bool Belt_TryInsert(BeltLine& line, int item, int64_t pos) {
    int64_t ahead = 0;
    const int k = InsertIndex(line, pos, &ahead);
    if (k < 0) return false;
    // The gap of the new item (to the item ahead, or to the end) and of the item behind it (now to the new one).
    BeltItem b;
    b.item = item;
    b.gap = ahead - pos;
    if (k < (int)line.items.size()) line.items[k].gap -= b.gap;
    else line.tailPos = pos;
    line.items.insert(line.items.begin() + k, b);
    line.packed = std::min(line.packed, k); // a packed run ahead of k stays packed; past it the cache is recomputed
    line.asleep = false;
    return true;
}

void Belt_Wake(BeltLine& line) { line.asleep = false; }

namespace {
// One tick of movement; true when any item moved.
bool MoveItems(BeltLine& L) {
    const size_t n = L.items.size();
    const int64_t speed = std::max(0, L.speed);
    int64_t d = std::min(speed, L.items[0].gap);
    L.items[0].gap -= d;
    bool moved = d > 0;
    size_t i = 1;
    if (d == 0) i = std::max<size_t>(1, std::min<size_t>(L.packed, n)); // the packed run behind a stuck front can't move
    for (; i < n && d < speed; ++i) {
        const int64_t nd = std::min(speed, L.items[i].gap - kBeltGap + d);
        L.items[i].gap += d - nd;
        d = nd;
        if (nd > 0) moved = true;
    }
    const int64_t lastMove = i < n ? speed : d; // the walk stopped at full speed: everything behind moved the full speed
    if (lastMove > 0) moved = true;
    L.tailPos += lastMove;
    // The packed run: the front at the end and items exactly kBeltGap apart behind it (they stay put next tick).
    if (L.items[0].gap == 0) {
        int p = std::max(1, L.packed);
        while (p < (int)n && L.items[p].gap == kBeltGap) ++p;
        L.packed = p;
    } else {
        L.packed = 0;
    }
    return moved;
}
} // namespace

void Belt_Step(std::vector<BeltLine>& lines, const BeltSinkFn& sink) {
    std::vector<char> changed(lines.size(), 0);
    for (size_t i = 0; i < lines.size(); ++i) {
        BeltLine& L = lines[i];
        if (L.asleep) continue;
        if (L.items.empty()) { L.asleep = true; continue; }
        bool moved = false;

        // The front item leaves first when it is at the end, so the items behind it can close up this same tick.
        if (L.items.front().gap <= 0) {
            const int item = L.items.front().item;
            bool left = false;
            if (L.end == BeltEnd::Sink && sink) {
                left = sink(L.sink, item);
            } else if (L.end == BeltEnd::Belt && L.outLine >= 0 && (size_t)L.outLine < lines.size() && (size_t)L.outLine != i) {
                left = Belt_TryInsert(lines[L.outLine], item, L.outPos); // the target's own items have the right of way
                if (left) changed[L.outLine] = 1;
            }
            if (left) {
                L.items.pop_front(); // the next item's gap to it becomes its gap to the end (the front was at the end)
                L.packed = 0;
                moved = true;
            }
        }

        if (!L.items.empty() && MoveItems(L)) moved = true;
        if (moved) changed[i] = 1;
        else L.asleep = true;
    }
    // A line that changed may have made room for the lines merging into it.
    for (size_t i = 0; i < lines.size(); ++i) {
        BeltLine& L = lines[i];
        if (L.asleep && L.end == BeltEnd::Belt && L.outLine >= 0 && (size_t)L.outLine < lines.size() && changed[L.outLine])
            L.asleep = false;
    }
}

uint32_t Belt_Checksum(const BeltLine& line) {
    uint32_t h = 2166136261u;
    auto mix = [&](uint64_t v) {
        for (int b = 0; b < 8; ++b) { h ^= (uint32_t)(v & 0xff); h *= 16777619u; v >>= 8; }
    };
    mix((uint64_t)(int64_t)line.id);
    mix(line.items.size());
    for (const BeltItem& it : line.items) { mix((uint64_t)(int64_t)it.item); mix((uint64_t)it.gap); }
    return h;
}

CInt3 Belt_PointAt(const ConduitRun& run, int64_t pos) {
    if (run.points.empty()) return {};
    auto scaled = [](CInt3 p) { return CInt3{ p.x * kBeltSub, p.y * kBeltSub, p.z * kBeltSub }; };
    if (pos <= 0) return scaled(run.points.front());
    int64_t base = 0;
    for (size_t i = 1; i < run.points.size(); ++i) {
        const CInt3 a = scaled(run.points[i - 1]), b = scaled(run.points[i]);
        const int64_t len = (int64_t)std::abs(b.x - a.x) + std::abs(b.y - a.y) + std::abs(b.z - a.z);
        if (pos <= base + len) {
            const int64_t t = pos - base;
            auto step = [&](int from, int to) { return (int)(from + (to > from ? t : to < from ? -t : 0)); };
            return { step(a.x, b.x), step(a.y, b.y), step(a.z, b.z) };
        }
        base += len;
    }
    return scaled(run.points.back());
}

int64_t Belt_PosOf(const ConduitRun& run, CInt3 slot) {
    int64_t base = 0;
    for (size_t i = 1; i < run.points.size(); ++i) {
        const CInt3 a = run.points[i - 1], b = run.points[i];
        const bool inX = std::min(a.x, b.x) <= slot.x && slot.x <= std::max(a.x, b.x);
        const bool inY = std::min(a.y, b.y) <= slot.y && slot.y <= std::max(a.y, b.y);
        const bool inZ = std::min(a.z, b.z) <= slot.z && slot.z <= std::max(a.z, b.z);
        if (inX && inY && inZ)
            return (base + std::abs(slot.x - a.x) + std::abs(slot.y - a.y) + std::abs(slot.z - a.z)) * kBeltSub;
        base += (int64_t)std::abs(b.x - a.x) + std::abs(b.y - a.y) + std::abs(b.z - a.z);
    }
    return -1;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
#include <chrono>
#include <cstdio>

namespace {
BeltLine Line(int id, int64_t length, int speed) {
    BeltLine l;
    l.id = id;
    l.length = length;
    l.speed = speed;
    return l;
}
std::vector<int64_t> Positions(const BeltLine& l) { return Belt_Positions(l); }
bool Spaced(const BeltLine& l) {
    const std::vector<int64_t> p = Belt_Positions(l);
    for (size_t k = 1; k < p.size(); ++k) if (p[k - 1] - p[k] < kBeltGap) return false;
    return true;
}
const BeltSinkFn kRefuse = [](int, int) { return false; };

// The position-based belt step this module had before the gap encoding (#1735 PR 2), kept as the reference: the
// differential test runs both on the same inputs and expects the same items at the same positions every tick.
struct RefItem { int item; int64_t pos; };
struct RefLine { int64_t length = 0; int speed = 0; std::vector<RefItem> items; BeltEnd end = BeltEnd::None; int sink = -1, outLine = -1; int64_t outPos = 0; bool asleep = false; };
bool RefHasRoom(const RefLine& l, int64_t pos) {
    if (pos < 0 || pos > l.length) return false;
    for (const RefItem& it : l.items) if (std::llabs(it.pos - pos) < kBeltGap) return false;
    return true;
}
bool RefInsert(RefLine& l, int item, int64_t pos) {
    if (!RefHasRoom(l, pos)) return false;
    auto at = std::find_if(l.items.begin(), l.items.end(), [&](const RefItem& it) { return it.pos < pos; });
    l.items.insert(at, RefItem{ item, pos });
    l.asleep = false;
    return true;
}
void RefStep(std::vector<RefLine>& lines, const BeltSinkFn& sink) {
    std::vector<char> changed(lines.size(), 0);
    for (size_t i = 0; i < lines.size(); ++i) {
        RefLine& L = lines[i];
        if (L.asleep) continue;
        if (L.items.empty()) { L.asleep = true; continue; }
        bool moved = false;
        if (L.items.front().pos >= L.length) {
            const int item = L.items.front().item;
            bool left = false;
            if (L.end == BeltEnd::Sink && sink) left = sink(L.sink, item);
            else if (L.end == BeltEnd::Belt && L.outLine >= 0 && (size_t)L.outLine < lines.size() && (size_t)L.outLine != i) {
                left = RefInsert(lines[L.outLine], item, L.outPos);
                if (left) changed[L.outLine] = 1;
            }
            if (left) { L.items.erase(L.items.begin()); moved = true; }
        }
        int64_t limit = L.length;
        for (RefItem& it : L.items) {
            const int64_t to = std::max(it.pos, std::min(it.pos + L.speed, limit));
            if (to != it.pos) { it.pos = to; moved = true; }
            limit = it.pos - kBeltGap;
        }
        if (moved) changed[i] = 1;
        else L.asleep = true;
    }
    for (size_t i = 0; i < lines.size(); ++i) {
        RefLine& L = lines[i];
        if (L.asleep && L.end == BeltEnd::Belt && L.outLine >= 0 && (size_t)L.outLine < lines.size() && changed[L.outLine]) L.asleep = false;
    }
}
} // namespace

// (#1735) The gap-encoded step against the position-based reference: random networks (chains, side merges, end to end),
// loaders, sinks that take or refuse at random and wakes; every tick the same items at the same positions, the same sleep
// states and the same sink calls.
TEST(BeltSim, GapEncodingMatchesThePositionReference) {
    for (uint32_t seed = 1; seed <= 40; ++seed) {
        uint32_t rng = seed * 2654435761u;
        auto next = [&](uint32_t mod) { rng = rng * 1664525u + 1013904223u; return (rng >> 8) % mod; };
        const int n = 2 + (int)next(4);
        std::vector<BeltLine> lines;
        std::vector<RefLine> ref;
        for (int k = 0; k < n; ++k) {
            BeltLine l = Line(k + 1, 16 + next(300), 1 + (int)next(12));
            RefLine r; r.length = l.length; r.speed = l.speed;
            if (k + 1 < n && next(3) != 0) { // feeds a later line: end to end or into its side
                l.end = r.end = BeltEnd::Belt;
                l.outLine = r.outLine = k + 1 + (int)next(n - k - 1);
                l.outPos = r.outPos = -1; // set below once the target's length is known
            } else {
                l.end = r.end = BeltEnd::Sink;
                l.sink = r.sink = k;
            }
            lines.push_back(l);
            ref.push_back(r);
        }
        for (int k = 0; k < n; ++k) if (lines[k].outPos < 0) {
            const int64_t len = lines[lines[k].outLine].length;
            lines[k].outPos = ref[k].outPos = next(2) ? 0 : (int64_t)next((uint32_t)len + 1);
        }
        uint32_t sa = seed, sb = seed; // the two sinks answer from identical streams
        std::vector<std::pair<int, int>> callsA, callsB;
        const BeltSinkFn sinkA = [&](int s, int item) { callsA.push_back({ s, item }); sa = sa * 1103515245u + 12345u; return ((sa >> 16) & 3) != 0; };
        const BeltSinkFn sinkB = [&](int s, int item) { callsB.push_back({ s, item }); sb = sb * 1103515245u + 12345u; return ((sb >> 16) & 3) != 0; };
        for (int t = 0; t < 3000; ++t) {
            for (int k = 0; k < n; ++k) {
                if (next(3) == 0) { // a loader at the start, or now and then anywhere
                    const int item = (int)next(9);
                    const int64_t pos = next(8) == 0 ? (int64_t)next((uint32_t)lines[k].length + 1) : 0;
                    const bool a = Belt_TryInsert(lines[k], item, pos), b = RefInsert(ref[k], item, pos);
                    ASSERT_EQ(a, b) << "seed " << seed << " t " << t << " line " << k << " insert at " << pos;
                }
                if (next(40) == 0) { Belt_Wake(lines[k]); ref[k].asleep = false; }
            }
            Belt_Step(lines, sinkA);
            RefStep(ref, sinkB);
            ASSERT_EQ(callsA, callsB) << "seed " << seed << " t " << t;
            for (int k = 0; k < n; ++k) {
                std::vector<int64_t> rp; std::vector<int> ri, gi;
                for (const RefItem& it : ref[k].items) { rp.push_back(it.pos); ri.push_back(it.item); }
                for (const BeltItem& it : lines[k].items) gi.push_back(it.item);
                ASSERT_EQ(Belt_Positions(lines[k]), rp) << "seed " << seed << " t " << t << " line " << k;
                ASSERT_EQ(gi, ri) << "seed " << seed << " t " << t << " line " << k;
                ASSERT_EQ(lines[k].asleep, ref[k].asleep) << "seed " << seed << " t " << t << " line " << k;
                if (!rp.empty()) ASSERT_EQ(lines[k].tailPos, rp.back()) << "seed " << seed << " t " << t << " line " << k;
            }
        }
    }
}

// (#1735) What the gap encoding buys: 40 lines of 25 cells, 100 items each (4000 items), stepped 2000 ticks, free-flowing
// (sinks take) and backed up (sinks refuse), against the position reference. Prints only (timings on a shared machine).
TEST(BeltSim, GapEncodingTiming) {
    auto build = [](auto& lines, auto make) {
        for (int k = 0; k < 40; ++k) {
            auto l = make(64 * 25);
            l.end = BeltEnd::Sink; l.sink = k;
            lines.push_back(l);
        }
    };
    for (int blocked = 0; blocked < 2; ++blocked) {
        const BeltSinkFn sink = [&](int, int) { return blocked == 0; };
        std::vector<BeltLine> g;
        build(g, [](int64_t len) { return Line(1, len, 6); });
        std::vector<RefLine> r;
        build(r, [](int64_t len) { RefLine l; l.length = len; l.speed = 6; return l; });
        for (auto& l : g) for (int64_t p = 0; p <= l.length - 16; p += 16) Belt_TryInsert(l, 1, p);
        for (auto& l : r) for (int64_t p = 0; p <= l.length - 16; p += 16) RefInsert(l, 1, p);
        auto run = [&](auto step) {
            const auto t0 = std::chrono::steady_clock::now();
            for (int t = 0; t < 2000; ++t) step();
            return std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / 2000.0;
        };
        // keep both awake (a blocked line would sleep and cost nothing in either)
        const double gu = run([&] { for (auto& l : g) { Belt_TryInsert(l, 1, 0); Belt_Wake(l); } Belt_Step(g, sink); });
        const double ru = run([&] { for (auto& l : r) { RefInsert(l, 1, 0); l.asleep = false; } RefStep(r, sink); });
        printf("beltsim timing (%s, 4000 items on 40 lines): gaps %.1f us/tick, positions %.1f us/tick (%.1fx)\n",
               blocked ? "backed up" : "flowing", gu, ru, ru / std::max(0.001, gu));
    }
}

TEST(BeltSim, InsertNeedsAGap) {
    BeltLine l = Line(1, 64, 6);
    EXPECT_TRUE(Belt_TryInsert(l, 7, 0));
    EXPECT_FALSE(Belt_TryInsert(l, 7, 0));
    EXPECT_FALSE(Belt_TryInsert(l, 7, 15));
    EXPECT_TRUE(Belt_TryInsert(l, 8, 16));
    EXPECT_TRUE(Belt_TryInsert(l, 9, 40));
    EXPECT_FALSE(Belt_TryInsert(l, 9, 65));   // past the end
    EXPECT_FALSE(Belt_TryInsert(l, 9, -1));
    EXPECT_EQ(Positions(l), (std::vector<int64_t>{ 40, 16, 0 })); // front first
    EXPECT_EQ(l.items[0].item, 9);
}

TEST(BeltSim, MovesAndBacksUp) {
    // one cell of belt, speed 6, the end refuses: items stop at the end and queue at the gap behind it
    std::vector<BeltLine> lines = { Line(1, 64, 6) };
    lines[0].end = BeltEnd::Sink;
    lines[0].sink = 3;
    int ticks = 0;
    for (; ticks < 200; ++ticks) {
        Belt_TryInsert(lines[0], 1, 0);
        Belt_Step(lines, kRefuse);
    }
    EXPECT_EQ(Positions(lines[0]), (std::vector<int64_t>{ 64, 48, 32, 16, 0 })); // 5 items: 4 per cell plus one at the start
    EXPECT_TRUE(lines[0].asleep);
    // the sink opens: one item a tick at most, and the rest close up behind it
    int taken = 0;
    const BeltSinkFn take = [&](int s, int) { EXPECT_EQ(s, 3); ++taken; return true; };
    Belt_Step(lines, take);
    EXPECT_EQ(taken, 0); // asleep until woken
    Belt_Wake(lines[0]);
    Belt_Step(lines, take);
    EXPECT_EQ(taken, 1);
    EXPECT_EQ(Positions(lines[0]), (std::vector<int64_t>{ 54, 38, 22, 6 }));
}

TEST(BeltSim, Throughput) {
    // a long belt into a sink that always takes: items/s = speed * 20 / gap once it is full
    std::vector<BeltLine> lines = { Line(1, 64 * 8, 6) };
    lines[0].end = BeltEnd::Sink;
    int taken = 0;
    const BeltSinkFn take = [&](int, int) { ++taken; return true; };
    for (int t = 0; t < 400; ++t) { Belt_TryInsert(lines[0], 1, 0); Belt_Step(lines, take); }
    const int before = taken;
    for (int t = 0; t < 800; ++t) { Belt_TryInsert(lines[0], 1, 0); Belt_Step(lines, take); }
    // 800 ticks at 6/16 items a tick = 300 (the insert can only happen every ceil(16/6) = 3 ticks: 266)
    EXPECT_GE(taken - before, 260);
    EXPECT_LE(taken - before, 300);
}

TEST(BeltSim, EndToEndAndMergeYield) {
    // line 0 (main) runs into a refusing sink; line 1 merges into its side at 32; line 2 continues line 1's start... keep simple:
    std::vector<BeltLine> lines = { Line(10, 128, 4), Line(11, 64, 4) };
    lines[0].end = BeltEnd::Sink;
    lines[1].end = BeltEnd::Belt;
    lines[1].outLine = 0;
    lines[1].outPos = 32;
    // fill the side belt, nothing on the main: items merge at 32 and travel on
    for (int t = 0; t < 60; ++t) { Belt_TryInsert(lines[1], 2, 0); Belt_Step(lines, kRefuse); }
    ASSERT_FALSE(lines[0].items.empty());
    for (const BeltItem& it : lines[0].items) EXPECT_EQ(it.item, 2);
    // the main belt backs up past the merge point: the side belt can't put anything in and sleeps
    for (int t = 0; t < 400; ++t) { Belt_TryInsert(lines[1], 2, 0); Belt_Step(lines, kRefuse); }
    EXPECT_TRUE(lines[0].asleep);
    EXPECT_TRUE(lines[1].asleep);
    EXPECT_TRUE(Spaced(lines[0]));
    // the main's sink takes one: the main moves, which wakes the side belt
    bool once = true;
    const BeltSinkFn takeOne = [&](int, int) { if (!once) return false; once = false; return true; };
    Belt_Wake(lines[0]);
    Belt_Step(lines, takeOne);
    EXPECT_FALSE(lines[0].asleep);
    EXPECT_FALSE(lines[1].asleep);

    // yield: main items passing the merge point have the right of way. Main carries item 1 continuously; the side gets in
    // only where the main has a gap.
    std::vector<BeltLine> m = { Line(20, 256, 4), Line(21, 64, 4) };
    m[0].end = BeltEnd::Sink;
    m[1].end = BeltEnd::Belt;
    m[1].outLine = 0;
    m[1].outPos = 64;
    int fromSide = 0, fromMain = 0;
    const BeltSinkFn count = [&](int, int item) { (item == 2 ? fromSide : fromMain)++; return true; };
    for (int t = 0; t < 2000; ++t) {
        Belt_TryInsert(m[0], 1, 0); // saturated main
        Belt_TryInsert(m[1], 2, 0);
        Belt_Step(m, count);
        ASSERT_TRUE(Spaced(m[0]));
    }
    EXPECT_GT(fromMain, 400);
    EXPECT_EQ(fromSide, 0); // a full main belt never lets the side in
}

// A fixed network (two lines merging into a third, a sink that takes 3 in 4 at random) for 1500 ticks: every line's
// checksum each 100 ticks, then how many items the sink took.
static std::vector<uint32_t> FixedBeltRun() {
    std::vector<BeltLine> lines = { Line(1, 200, 6), Line(2, 120, 4), Line(3, 300, 6) };
    lines[0].end = BeltEnd::Belt; lines[0].outLine = 2; lines[0].outPos = 0;
    lines[1].end = BeltEnd::Belt; lines[1].outLine = 2; lines[1].outPos = 96;
    lines[2].end = BeltEnd::Sink; lines[2].sink = 0;
    uint32_t seed = 12345;
    int accepted = 0;
    const BeltSinkFn sink = [&](int, int) { seed = seed * 1664525u + 1013904223u; if (seed >> 30) { ++accepted; return true; } return false; };
    std::vector<uint32_t> sums;
    for (int t = 0; t < 1500; ++t) {
        Belt_TryInsert(lines[0], t % 5, 0);
        if (t % 3 == 0) Belt_TryInsert(lines[1], 100 + t % 7, 0);
        Belt_Step(lines, sink);
        if (t % 50 == 0) for (BeltLine& l : lines) Belt_Wake(l); // a sink freeing room now and then
        if (t % 100 == 99) for (const BeltLine& l : lines) sums.push_back(Belt_Checksum(l));
    }
    sums.push_back((uint32_t)accepted);
    return sums;
}

// (#1015) The determinism guard across platforms: the fixed run hashes to one constant, which the Windows (w64devkit) and
// Linux (GCC) CI legs both check. Only one leg failing means the platforms diverge (a float or libm call crept into the sim,
// or undefined behaviour): fix that, never the constant. Both legs failing alike after a deliberate change to the belt
// rules: put the new value here in the same PR and say why in its description.
TEST(BeltSim, FixedRunMatchesTheGoldenHash) {
    uint64_t h = 1469598103934665603ull; // FNV-1a over the run's checksums
    for (uint32_t s : FixedBeltRun()) h = (h ^ s) * 1099511628211ull;
    constexpr uint64_t kGolden = 15827766116492930527ull;
    EXPECT_EQ(h, kGolden) << "the fixed belt run changed (#1015): see the comment above before updating the constant";
}

TEST(BeltSim, Deterministic) {
    const std::vector<uint32_t> a = FixedBeltRun(), b = FixedBeltRun();
    EXPECT_EQ(a, b);
    EXPECT_GT(a.back(), 0u);
    // the checksum sees a moved item
    BeltLine l = Line(1, 64, 6);
    Belt_TryInsert(l, 1, 0);
    const uint32_t before = Belt_Checksum(l);
    l.items[0].gap -= 1;
    EXPECT_NE(Belt_Checksum(l), before);
}

TEST(BeltSim, Geometry) {
    ConduitRun r;
    r.id = 5; r.kind = ConduitKind::Belt; r.material = "basic"; r.size = kBeltWidth;
    r.points = { { 2, 0, 4 }, { 18, 0, 4 }, { 18, 0, 20 } };
    const BeltLine l = Belt_MakeLine(r, 6);
    EXPECT_EQ(l.id, 5);
    EXPECT_EQ(l.length, 32 * kBeltSub);
    EXPECT_EQ(Belt_PointAt(r, 0), (CInt3{ 16, 0, 32 }));
    EXPECT_EQ(Belt_PointAt(r, 8), (CInt3{ 24, 0, 32 }));
    EXPECT_EQ(Belt_PointAt(r, 16 * kBeltSub), (CInt3{ 144, 0, 32 }));      // the corner
    EXPECT_EQ(Belt_PointAt(r, 16 * kBeltSub + 3), (CInt3{ 144, 0, 35 }));
    EXPECT_EQ(Belt_PointAt(r, 9999), (CInt3{ 144, 0, 160 }));             // clamped
    EXPECT_EQ(Belt_PosOf(r, { 2, 0, 4 }), 0);
    EXPECT_EQ(Belt_PosOf(r, { 10, 0, 4 }), 8 * kBeltSub);
    EXPECT_EQ(Belt_PosOf(r, { 18, 0, 12 }), 24 * kBeltSub);
    EXPECT_EQ(Belt_PosOf(r, { 10, 0, 5 }), -1);
    for (int64_t p = 0; p <= l.length; p += 8) {                           // round trip on slot positions
        const CInt3 q = Belt_PointAt(r, p);
        EXPECT_EQ(Belt_PosOf(r, { q.x / kBeltSub, q.y / kBeltSub, q.z / kBeltSub }), p);
    }
}
#endif
