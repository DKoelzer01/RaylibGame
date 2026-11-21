#ifndef ROTATION_H
#define ROTATION_H

// Rotation networks (#159, #1268; the Factory and Progression doc, section 5): Create-style mechanical power.
//
// This is the core only: a graph of nodes (shafts, gears, gearboxes, sources, consumers) joined by links that each carry a
// fixed, signed speed ratio, and the solver the factory tick and the summary sim both run on it. No parts, no geometry, no
// rendering: the factory builds a Graph from a plot's parts (step 2, #1269) and reads the speeds and statuses back.
//
// The rules (decided, doc 5.2):
// - A connected set of nodes is one net. Every node's speed is a fixed multiple of the net's reference speed (its ratio,
//   from the links). Links that disagree around a loop, or sources turning opposite ways, jam the net: nothing turns.
// - The fastest source sets the net's speed; each source adds its stress capacity at its own speed (SU per RPM x RPM).
// - Demand = the running consumers' stress per RPM x their RPM. ratio = min(1000, capacity x 1000 / demand), rounded down,
//   one division a tick. Every node turns at its nominal speed x ratio / 1000: an overloaded net turns slower, smoothly.
// - A consumer under its minimum speed is "too slow" (it makes no progress) and still counts toward the demand, so stalling
//   one never speeds the others up within a tick: no oscillation, and no dependence on the order nodes are visited.
// - Integers only (speeds in milli-RPM), nets and nodes in a fixed order: the same graph gives the same result everywhere,
//   and listing the nodes in another order gives the same speeds.

#include <cstdint>
#include <vector>

namespace Rotation {

// ---- Tuning (the doc's guesses, as named constants; tuned by playtest) ----
constexpr int kMaxRpm = 256;              // nothing turns faster (Create's cap): a node geared above it runs at it
constexpr int kMaxRatioSpan = 4096;       // fastest / slowest node of one net; beyond it the net jams (a runaway gear chain)
constexpr int kMaxProgressPermille = 4000; // a consumer above its rated speed works faster, up to 4x

struct SourceSpec { int rpm, suPerRpm; };              // nominal speed, stress capacity per RPM (capacity = suPerRpm x rpm)
struct ConsumerSpec { int stressPerRpm, minRpm, ratedRpm; }; // stress per RPM while working; stalls under minRpm; 1x at ratedRpm
constexpr SourceSpec kHandCrank   { 32, 8 };    // 256 SU while held
constexpr SourceSpec kWaterWheel  { 16, 16 };   // 256 SU at full current (step 3 scales the speed by the flow)
constexpr SourceSpec kWindmillSail{ 8, 32 };    // per sail: 256 SU at full wind (doc: 256 SU a sail, Create's scale)
constexpr SourceSpec kSteamEngine { 64, 16 };   // 1024 SU (tier 3)
constexpr ConsumerSpec kMillstone { 4, 16, 32 }; // doc: a millstone needs 16 RPM
constexpr ConsumerSpec kSaw       { 4, 16, 32 };
constexpr ConsumerSpec kPress     { 8, 16, 32 };
constexpr ConsumerSpec kBellows   { 2, 8, 16 };
constexpr ConsumerSpec kMixer     { 4, 32, 64 };
constexpr ConsumerSpec kMechPump  { 4, 8, 32 };

// ---- The graph (the factory rebuilds it when a plot's layout changes; per tick it updates the sources and `running`) ----
enum class NodeKind : uint8_t { Transmission, Source, Consumer };
struct Node {
    NodeKind kind = NodeKind::Transmission;
    // Source: its speed now, signed (the sign is its turning direction on the reference axis; 0: not turning: a crank no one
    // holds, a windmill in calm) and its stress capacity per RPM of that speed.
    int sourceRpm = 0, suPerRpm = 0;
    // Consumer: stress per RPM while running, the speed under which it stalls, and whether it has work this tick (an idle
    // consumer adds no stress, decided: only running consumers load the net).
    int stressPerRpm = 0, minRpm = 0;
    bool running = false;
};
// speed(b) = speed(a) x num / den. A shaft joint is 1/1, meshed cogs of one size -1/1, a small cog driving a large one -1/2,
// a 2:1 gearbox 2/1. num != 0, den > 0.
struct Link { int a = 0, b = 0, num = 1, den = 1; };
struct Graph {
    std::vector<Node> nodes;
    std::vector<Link> links;
    int add(const Node& n) { nodes.push_back(n); return (int)nodes.size() - 1; }
    int addShaft() { return add(Node{}); }
    int addSource(SourceSpec s, int rpm) { Node n; n.kind = NodeKind::Source; n.sourceRpm = rpm; n.suPerRpm = s.suPerRpm; return add(n); }
    int addConsumer(ConsumerSpec c, bool running = true) { Node n; n.kind = NodeKind::Consumer; n.stressPerRpm = c.stressPerRpm; n.minRpm = c.minRpm; n.running = running; return add(n); }
    void link(int a, int b, int num = 1, int den = 1) { links.push_back({ a, b, num, den }); }
};

enum class Jam : uint8_t { None, GearLoop, RatioSpan, Direction };
const char* JamName(Jam j); // "", "gears disagree", "geared too far", "sources turn opposite ways"

// What a layout gives (Build): the nets and each node's ratio to its net's reference. Kept until the layout changes.
struct Layout {
    struct Net {
        std::vector<int> nodes;        // ascending
        Jam jam = Jam::None;           // from the links (a loop that disagrees, a runaway span); Step adds Direction
    };
    std::vector<Net> nets;             // ordered by their lowest node
    std::vector<int> netOf;            // per node
    std::vector<int64_t> num, den;     // per node: speed = reference x num / den (reduced, den > 0)
};
Layout Build(const Graph& g);

enum class Status : uint8_t { Stopped, Running, Slowed, TooSlow, Idle, Jammed };
const char* StatusName(Status s);      // "Stopped", "Running", "Slowed", "Too slow", "Idle", "Jammed"

struct NetResult {
    Jam jam = Jam::None;
    int64_t capacity = 0;              // SU: the turning sources' suPerRpm x |rpm|
    int64_t demandMilli = 0;           // SU x 1000: the running consumers' stressPerRpm x |nominal milli-RPM|
    int ratio = 1000;                  // permille: every node turns at nominal x ratio / 1000
    int64_t fastestMilli = 0;          // the fastest node's |effective speed| (milli-RPM), for the tooltip
};
struct Result {
    std::vector<int64_t> nominalMilli; // per node: signed milli-RPM before the stress ratio (0 if jammed or no source turns)
    std::vector<int64_t> speedMilli;   // per node: signed effective milli-RPM (what turns, animates and works)
    std::vector<Status> status;        // per node (transmissions and sources: Stopped / Running / Slowed / Jammed)
    std::vector<NetResult> nets;       // per Layout net
};
// One tick: the sources' speeds and the consumers' `running` as they are now in g. g's links must be the ones Build saw.
void Step(const Graph& g, const Layout& layout, Result& out);

// The formula itself, for the summary sim (#1261) to apply to its rates: permille, rounded down, 1000 when nothing loads it.
int StressRatio(int64_t capacity, int64_t demandMilli);
// What the gauge shows: demand as a percentage of capacity ("stress 140% - running at 71%"); 0 with no demand, and a large
// number when there is demand and no capacity.
int StressPercent(const NetResult& n);
// A consumer's recipe progress this tick, permille of its rated speed (0 when too slow, idle or jammed; up to 4x).
int ProgressPermille(ConsumerSpec spec, int64_t speedMilli, Status status);

} // namespace Rotation

#endif
