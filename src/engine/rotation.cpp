#include "rotation.h"
#include <algorithm>
#include <cstdlib>
#include <numeric>

// Rotation networks (#159, #1268). See rotation.h.

namespace Rotation {

namespace {
// A ratio's numerator or denominator past this: a runaway chain (jams the net). (#1858 M8: 2^24 jammed a 16-stage 3/2 chain,
// 3^16 = 43M, whose speed span (~657) is well inside kMaxRatioSpan; the span check is what limits gearing. 2^40 keeps every
// product in Build/Solve inside i128: rpm x den x num x 1000 < 2^14 x 2^40 x 2^40 x 2^10.)
constexpr int64_t kMaxTerm = int64_t(1) << 40;
using i128 = __int128;

int64_t Abs(int64_t v) { return v < 0 ? -v : v; }
// num / den reduced, den > 0. False if a term runs past kMaxTerm.
bool Reduce(i128 num, i128 den, int64_t& outNum, int64_t& outDen) {
    if (den < 0) { num = -num; den = -den; }
    i128 a = num < 0 ? -num : num, b = den;
    while (b != 0) { const i128 t = a % b; a = b; b = t; }
    if (a > 1) { num /= a; den /= a; }
    if ((num < 0 ? -num : num) > kMaxTerm || den > kMaxTerm) return false;
    outNum = (int64_t)num; outDen = (int64_t)den;
    return true;
}
// a / b rounded toward zero (b > 0), so a net turning the other way has the same magnitudes.
int64_t DivTrunc(i128 a, i128 b) { return (int64_t)(a / b); } // C++ integer division truncates toward zero
}

const char* JamName(Jam j) {
    switch (j) {
    case Jam::GearLoop: return "gears disagree";
    case Jam::RatioSpan: return "geared too far";
    case Jam::Direction: return "sources turn opposite ways";
    default: return "";
    }
}

const char* StatusName(Status s) {
    switch (s) {
    case Status::Running: return "Running";
    case Status::Slowed: return "Slowed";
    case Status::TooSlow: return "Too slow";
    case Status::Idle: return "Idle";
    case Status::Jammed: return "Jammed";
    default: return "Stopped";
    }
}

Layout Build(const Graph& g) {
    Layout L;
    const int n = (int)g.nodes.size();
    L.netOf.assign(n, -1);
    L.num.assign(n, 1);
    L.den.assign(n, 1);
    // Adjacency: speed(to) = speed(from) x num / den, both ways (the reverse of num/den is den/num).
    struct Edge { int to; int64_t num, den; };
    std::vector<std::vector<Edge>> adj(n);
    std::vector<char> selfJam(n, 0);
    for (const Link& k : g.links) {
        if (k.a < 0 || k.b < 0 || k.a >= n || k.b >= n || k.num == 0 || k.den <= 0) continue; // a malformed link joins nothing
        if (k.a == k.b) { if (k.num != k.den) selfJam[k.a] = 1; continue; }
        adj[k.a].push_back({ k.b, k.num, k.den });
        adj[k.b].push_back({ k.a, k.num < 0 ? -(int64_t)k.den : (int64_t)k.den, Abs(k.num) });
    }
    std::vector<int> queue;
    for (int start = 0; start < n; ++start) {
        if (L.netOf[start] >= 0) continue;
        const int id = (int)L.nets.size();
        L.nets.push_back({});
        Layout::Net& net = L.nets.back();
        L.netOf[start] = id;
        L.num[start] = 1; L.den[start] = 1;
        queue.assign(1, start);
        for (size_t qi = 0; qi < queue.size(); ++qi) {
            const int u = queue[qi];
            if (selfJam[u]) net.jam = Jam::GearLoop;
            for (const Edge& e : adj[u]) {
                int64_t rn = 1, rd = 1;
                if (!Reduce((i128)L.num[u] * e.num, (i128)L.den[u] * e.den, rn, rd)) {
                    if (net.jam == Jam::None) net.jam = Jam::RatioSpan;
                    rn = 1; rd = 1; // (any value: the net is jammed)
                }
                if (L.netOf[e.to] < 0) {
                    L.netOf[e.to] = id; L.num[e.to] = rn; L.den[e.to] = rd;
                    queue.push_back(e.to);
                } else if ((L.num[e.to] != rn || L.den[e.to] != rd) && net.jam == Jam::None) {
                    net.jam = Jam::GearLoop; // a loop whose gears disagree: in Create the parts would break; here nothing turns
                }
            }
        }
        net.nodes = queue;
        std::sort(net.nodes.begin(), net.nodes.end());
        // The span between the fastest and the slowest node, whichever node is the reference (so it does not depend on order).
        if (net.jam == Jam::None && net.nodes.size() > 1) {
            int fast = net.nodes[0], slow = net.nodes[0];
            for (int v : net.nodes) {
                if ((i128)Abs(L.num[v]) * L.den[fast] > (i128)Abs(L.num[fast]) * L.den[v]) fast = v;
                if ((i128)Abs(L.num[v]) * L.den[slow] < (i128)Abs(L.num[slow]) * L.den[v]) slow = v;
            }
            // |fast| / |slow| > kMaxRatioSpan
            if ((i128)Abs(L.num[fast]) * L.den[slow] > (i128)kMaxRatioSpan * Abs(L.num[slow]) * L.den[fast]) net.jam = Jam::RatioSpan;
        }
    }
    return L;
}

int StressRatio(int64_t capacity, int64_t demandMilli) {
    if (demandMilli <= 0) return 1000;
    if (capacity <= 0) return 0;
    const i128 r = (i128)capacity * 1000 * 1000 / demandMilli;
    return r >= 1000 ? 1000 : (int)r;
}

int StressPercent(const NetResult& n) {
    if (n.demandMilli <= 0) return 0;
    if (n.capacity <= 0) return 9999;
    const int64_t p = n.demandMilli / 10 / n.capacity; // demand x 100 / capacity (demand in milli)
    return (int)std::min<int64_t>(p, 9999);
}

int ProgressPermille(ConsumerSpec spec, int64_t speedMilli, Status status) {
    if ((status != Status::Running && status != Status::Slowed) || spec.ratedRpm <= 0) return 0;
    return (int)std::min<int64_t>(kMaxProgressPermille, Abs(speedMilli) / spec.ratedRpm);
}

void Step(const Graph& g, const Layout& L, Result& out) {
    const int n = (int)g.nodes.size();
    out.nominalMilli.assign(n, 0);
    out.speedMilli.assign(n, 0);
    out.status.assign(n, Status::Stopped);
    out.nets.assign(L.nets.size(), {});
    if ((int)L.netOf.size() != n) return; // a stale layout: nothing turns until it is rebuilt
    for (size_t ni = 0; ni < L.nets.size(); ++ni) {
        const Layout::Net& net = L.nets[ni];
        NetResult& r = out.nets[ni];
        r.jam = net.jam;
        // The reference speed: each turning source asks for rpm / its ratio; they must agree in direction, the fastest wins.
        bool any = false;
        i128 refNum = 0, refDen = 1; // the reference speed, refNum / refDen RPM (refDen > 0)
        int sign = 0;
        for (int v : net.nodes) {
            const Node& node = g.nodes[v];
            if (node.kind != NodeKind::Source || node.sourceRpm == 0) continue;
            const int rpm = std::clamp(node.sourceRpm, -kMaxRpm, kMaxRpm);
            r.capacity += (int64_t)std::max(0, node.suPerRpm) * std::abs(rpm);
            // rpm / (num / den) = rpm x den / num
            i128 sn = (i128)rpm * L.den[v], sd = L.num[v];
            if (sd < 0) { sn = -sn; sd = -sd; }
            const int s = sn > 0 ? 1 : -1;
            if (sign != 0 && s != sign && r.jam == Jam::None) r.jam = Jam::Direction;
            sign = s;
            const i128 absN = sn < 0 ? -sn : sn, absRef = refNum < 0 ? -refNum : refNum;
            if (!any || absN * refDen > absRef * sd) { refNum = sn; refDen = sd; }
            any = true;
        }
        if (r.jam != Jam::None) {
            r.capacity = 0;
            for (int v : net.nodes) out.status[v] = Status::Jammed;
            continue;
        }
        if (!any) continue; // nothing turns: every node Stopped
        for (int v : net.nodes) {
            int64_t m = DivTrunc(refNum * L.num[v] * 1000, refDen * L.den[v]);
            m = std::clamp<int64_t>(m, -(int64_t)kMaxRpm * 1000, (int64_t)kMaxRpm * 1000);
            out.nominalMilli[v] = m;
            const Node& node = g.nodes[v];
            if (node.kind == NodeKind::Consumer && node.running) r.demandMilli += (int64_t)std::max(0, node.stressPerRpm) * Abs(m);
        }
        r.ratio = StressRatio(r.capacity, r.demandMilli);
        for (int v : net.nodes) {
            const int64_t nominal = out.nominalMilli[v];
            const int64_t eff = DivTrunc((i128)nominal * r.ratio, 1000);
            out.speedMilli[v] = eff;
            r.fastestMilli = std::max(r.fastestMilli, Abs(eff));
            const Node& node = g.nodes[v];
            Status st;
            if (eff == 0) st = Status::Stopped;
            else if (node.kind == NodeKind::Consumer && !node.running) st = Status::Idle;
            else if (node.kind == NodeKind::Consumer && Abs(eff) < (int64_t)node.minRpm * 1000) st = Status::TooSlow;
            else st = r.ratio < 1000 ? Status::Slowed : Status::Running;
            out.status[v] = st;
        }
    }
}

} // namespace Rotation

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1268): speeds through gears, the stress ratio, stalls, jams, and order independence ----
#include <gtest/gtest.h>
#include <random>

using namespace Rotation;

namespace {
Result Solve(const Graph& g) { Result r; Step(g, Build(g), r); return r; }
}

TEST(Rotation, OneCrankTurnsAMillstoneAtFullSpeed) {
    Graph g;
    const int crank = g.addSource(kHandCrank, kHandCrank.rpm);
    const int shaft = g.addShaft();
    const int mill = g.addConsumer(kMillstone);
    g.link(crank, shaft); g.link(shaft, mill);
    const Result r = Solve(g);
    ASSERT_EQ(r.nets.size(), 1u);
    EXPECT_EQ(r.nets[0].capacity, 256);              // 8 SU/rpm x 32 rpm
    EXPECT_EQ(r.nets[0].demandMilli, 4 * 32000);     // 4 SU/rpm x 32 rpm
    EXPECT_EQ(r.nets[0].ratio, 1000);
    EXPECT_EQ(r.speedMilli[mill], 32000);
    EXPECT_EQ(r.speedMilli[shaft], 32000);
    EXPECT_EQ(r.status[mill], Status::Running);
    EXPECT_EQ(StressPercent(r.nets[0]), 50);
    EXPECT_EQ(ProgressPermille(kMillstone, r.speedMilli[mill], r.status[mill]), 1000);
}

TEST(Rotation, AnOverloadedNetSlowsEveryNodeByCapacityOverDemand) {
    Graph g;
    const int crank = g.addSource(kHandCrank, kHandCrank.rpm); // 256 SU
    int mills[4];
    for (int& m : mills) { m = g.addConsumer(kMillstone); g.link(crank, m); } // 4 x 128 SU = 512
    const Result r = Solve(g);
    EXPECT_EQ(r.nets[0].ratio, 500);
    EXPECT_EQ(StressPercent(r.nets[0]), 200);
    EXPECT_EQ(r.speedMilli[crank], 16000);           // the source turns slower too
    for (int m : mills) {
        EXPECT_EQ(r.speedMilli[m], 16000);
        EXPECT_EQ(r.status[m], Status::Slowed);       // exactly at the 16 rpm minimum: still works
        EXPECT_EQ(ProgressPermille(kMillstone, r.speedMilli[m], r.status[m]), 500);
    }
}

TEST(Rotation, AStalledConsumerStillCountsTowardTheDemand) {
    Graph g;
    const int crank = g.addSource(kHandCrank, kHandCrank.rpm); // 256 SU
    const int mixer = g.addConsumer(kMixer);                    // 128 SU at 32 rpm, needs 32 rpm
    const int a = g.addConsumer(kMillstone), b = g.addConsumer(kMillstone);
    g.link(crank, mixer); g.link(crank, a); g.link(crank, b);
    const Result r = Solve(g);
    EXPECT_EQ(r.nets[0].demandMilli, 3 * 128000);    // the mixer's stress stays in
    EXPECT_EQ(r.nets[0].ratio, 666);
    EXPECT_EQ(r.status[mixer], Status::TooSlow);
    EXPECT_EQ(ProgressPermille(kMixer, r.speedMilli[mixer], r.status[mixer]), 0);
    EXPECT_EQ(r.speedMilli[a], 21312);               // not sped up by the mixer's stall
    EXPECT_EQ(r.status[a], Status::Slowed);
}

TEST(Rotation, AGearboxDoublesTheSpeedAndTheStressBehindIt) {
    Graph g;
    const int crank = g.addSource(kHandCrank, kHandCrank.rpm);
    const int box = g.addShaft();
    const int fast = g.addShaft();
    const int mill = g.addConsumer(kMillstone);
    g.link(crank, box); g.link(box, fast, 2, 1); g.link(fast, mill);
    Result r = Solve(g);
    EXPECT_EQ(r.speedMilli[fast], 64000);
    EXPECT_EQ(r.nets[0].demandMilli, 4 * 64000);     // 256 SU: the crank's whole capacity
    EXPECT_EQ(r.nets[0].ratio, 1000);
    EXPECT_EQ(ProgressPermille(kMillstone, r.speedMilli[mill], r.status[mill]), 2000);
    // One more millstone on the slow side: 384 SU of demand on 256.
    const int slowMill = g.addConsumer(kMillstone);
    g.link(crank, slowMill);
    r = Solve(g);
    EXPECT_EQ(r.nets[0].ratio, 666);
    EXPECT_EQ(r.speedMilli[mill], 42624);
    EXPECT_EQ(r.speedMilli[slowMill], 21312);
}

TEST(Rotation, MeshedCogsTurnTheOtherWayAndSmallDrivesLargeAtHalfSpeed) {
    Graph g;
    const int crank = g.addSource(kHandCrank, kHandCrank.rpm);
    const int large = g.addShaft();
    g.link(crank, large, -1, 2);
    const int down = g.addShaft();
    g.link(large, down, 1, 3);                       // a 1:3 reduction: 16 / 3 rpm
    const Result r = Solve(g);
    EXPECT_EQ(r.speedMilli[large], -16000);
    EXPECT_EQ(r.speedMilli[down], -5333);            // rounded toward zero
}

TEST(Rotation, SeparateNetsDoNotAffectEachOther) {
    Graph g;
    const int c1 = g.addSource(kHandCrank, kHandCrank.rpm), m1 = g.addConsumer(kMillstone);
    const int c2 = g.addSource(kHandCrank, kHandCrank.rpm), m2 = g.addConsumer(kPress), m3 = g.addConsumer(kPress);
    g.link(c1, m1); g.link(c2, m2); g.link(c2, m3);  // the second: 2 x 256 on 256
    const Result r = Solve(g);
    ASSERT_EQ(r.nets.size(), 2u);
    EXPECT_EQ(r.nets[0].ratio, 1000);
    EXPECT_EQ(r.nets[1].ratio, 500);
    EXPECT_EQ(r.speedMilli[m1], 32000);
    EXPECT_EQ(r.speedMilli[m2], 16000);
}

TEST(Rotation, SourcesTurningOppositeWaysJamTheNet) {
    Graph g;
    const int a = g.addSource(kHandCrank, 32), b = g.addSource(kHandCrank, -32);
    const int mill = g.addConsumer(kMillstone);
    g.link(a, mill); g.link(mill, b);
    Result r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::Direction);
    EXPECT_EQ(r.status[mill], Status::Jammed);
    EXPECT_EQ(r.speedMilli[mill], 0);
    // Through meshed cogs the same two agree.
    g.links.clear();
    g.link(a, mill); g.link(mill, b, -1, 1);
    r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::None);
    EXPECT_EQ(r.speedMilli[mill], 32000);
    EXPECT_EQ(r.nets[0].capacity, 512);
}

TEST(Rotation, AGearLoopMustAgreeWithItself) {
    Graph g;
    const int crank = g.addSource(kHandCrank, 32);
    const int a = g.addShaft(), b = g.addShaft();
    g.link(crank, a); g.link(a, b, 2, 1); g.link(b, crank, 1, 2); // up then back down: consistent
    Result r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::None);
    EXPECT_EQ(r.speedMilli[b], 64000);
    g.links.back() = { b, crank, 1, 1 };                         // up, then 1:1 back: disagrees
    r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::GearLoop);
    EXPECT_EQ(r.status[a], Status::Jammed);
    EXPECT_EQ(r.nets[0].capacity, 0);
}

TEST(Rotation, ARunawayGearChainJamsAndTheTopSpeedIsCapped) {
    Graph g;
    int prev = g.addSource(kHandCrank, 32);
    for (int i = 0; i < 12; ++i) { const int s = g.addShaft(); g.link(prev, s, 2, 1); prev = s; } // x4096: allowed
    Result r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::None);
    EXPECT_EQ(r.speedMilli[prev], kMaxRpm * 1000);  // 32 x 4096 rpm, capped
    const int last = g.addShaft(); g.link(prev, last, 2, 1);                                        // x8192: too far
    r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::RatioSpan);
    EXPECT_EQ(r.speedMilli[last], 0);
}

// (#1858 M8) A long chain of 3:2 steps: the terms grow (3^16 = 43M) though the speed span stays small (1.5^16 ~ 657), so
// only the span may jam it.
TEST(Rotation, ALongOddRatioChainWithinTheSpanTurns) {
    Graph g;
    int prev = g.addSource(kHandCrank, 32);
    for (int i = 0; i < 16; ++i) { const int s = g.addShaft(); g.link(prev, s, 3, 2); prev = s; }
    Result r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::None) << "x657: inside kMaxRatioSpan";
    EXPECT_GT(r.speedMilli[prev], 0) << "it turns";
    for (int i = 0; i < 5; ++i) { const int s = g.addShaft(); g.link(prev, s, 3, 2); prev = s; } // 1.5^21 ~ 4988: past the span
    r = Solve(g);
    EXPECT_EQ(r.nets[0].jam, Jam::RatioSpan);
}

TEST(Rotation, NoTurningSourceMeansEverythingStopped) {
    Graph g;
    const int crank = g.addSource(kHandCrank, 0); // nobody holds it
    const int mill = g.addConsumer(kMillstone);
    g.link(crank, mill);
    const int lone = g.addConsumer(kSaw);          // on no shaft at all
    const Result r = Solve(g);
    EXPECT_EQ(r.status[mill], Status::Stopped);
    EXPECT_EQ(r.status[lone], Status::Stopped);
    EXPECT_EQ(r.nets[0].ratio, 1000);
    EXPECT_EQ(r.nets[0].capacity, 0);
}

TEST(Rotation, IdleConsumersAddNoStress) {
    Graph g;
    const int crank = g.addSource(kHandCrank, 32);
    const int a = g.addConsumer(kMillstone), b = g.addConsumer(kMillstone), idle = g.addConsumer(kMillstone, false);
    g.link(crank, a); g.link(crank, b); g.link(crank, idle);
    const Result r = Solve(g);
    EXPECT_EQ(r.nets[0].demandMilli, 256000);
    EXPECT_EQ(r.nets[0].ratio, 1000);
    EXPECT_EQ(r.status[idle], Status::Idle);
    EXPECT_EQ(r.speedMilli[idle], 32000);          // it still turns with its shaft
}

TEST(Rotation, TheFastestSourceSetsTheSpeedAndAllAddCapacity) {
    Graph g;
    const int crank = g.addSource(kHandCrank, 32), engine = g.addSource(kSteamEngine, 64);
    const int mill = g.addConsumer(kMillstone);
    g.link(crank, mill); g.link(mill, engine);
    const Result r = Solve(g);
    EXPECT_EQ(r.nets[0].capacity, 256 + 1024);
    EXPECT_EQ(r.speedMilli[mill], 64000);
    EXPECT_EQ(r.speedMilli[crank], 64000);
}

TEST(Rotation, TheStressRatioFormula) {
    EXPECT_EQ(StressRatio(256, 0), 1000);
    EXPECT_EQ(StressRatio(0, 1000), 0);
    EXPECT_EQ(StressRatio(256, 512000), 500);
    EXPECT_EQ(StressRatio(256, 384000), 666);       // rounded down
    EXPECT_EQ(StressRatio(1000, 1000), 1000);       // capped at full speed
}

// Listing the nodes and links in any order, and writing a link the other way round, gives the same speeds and statuses.
TEST(Rotation, TheResultDoesNotDependOnTheOrderOfNodesOrLinks) {
    Graph g;
    const int crank = g.addSource(kHandCrank, 32), engine = g.addSource(kWaterWheel, -16);
    int prev = crank;
    for (int i = 0; i < 8; ++i) { const int s = g.addShaft(); g.link(prev, s, i % 3 == 0 ? -1 : 1, i % 4 == 1 ? 2 : 1); prev = s; }
    g.link(prev, engine, 1, 1);
    for (int i = 0; i < 6; ++i) { const int c = g.addConsumer(i % 2 ? kPress : kMillstone, i != 4); g.link(2 + i, c, i == 3 ? 3 : 1, i == 3 ? 2 : 1); }
    const int other = g.addSource(kHandCrank, 32), otherMill = g.addConsumer(kMixer);
    g.link(other, otherMill, 1, 3);
    const Result want = Solve(g);
    ASSERT_EQ(want.nets[0].jam, Jam::None) << "a busy, consistent net";
    std::mt19937 rng(42);
    for (int trial = 0; trial < 20; ++trial) {
        const int n = (int)g.nodes.size();
        std::vector<int> to(n);
        std::iota(to.begin(), to.end(), 0);
        std::shuffle(to.begin(), to.end(), rng);
        Graph p;
        p.nodes.resize(n);
        for (int v = 0; v < n; ++v) p.nodes[to[v]] = g.nodes[v];
        for (const Link& k : g.links) {
            if (rng() & 1) p.links.push_back({ to[k.a], to[k.b], k.num, k.den });
            else p.links.push_back({ to[k.b], to[k.a], k.num < 0 ? -k.den : k.den, std::abs(k.num) }); // the same link, reversed
        }
        std::shuffle(p.links.begin(), p.links.end(), rng);
        const Result got = Solve(p);
        for (int v = 0; v < n; ++v) {
            EXPECT_EQ(got.speedMilli[to[v]], want.speedMilli[v]) << "trial " << trial << " node " << v;
            EXPECT_EQ(got.status[to[v]], want.status[v]) << "trial " << trial << " node " << v;
        }
    }
}
#endif
