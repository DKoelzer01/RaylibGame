#include "rotationplot.h"
#include <algorithm>
#include <set>
#include <tuple>

// Rotation on a plot (#1269). See rotationplot.h.

using Rotation::Cell; using Rotation::Node; using Rotation::NodeKind; using Rotation::Status; using Rotation::NetResult; using Rotation::Result;
using Rotation::MakeShaft; using Rotation::MakeCog; using Rotation::MakeGearbox; using Rotation::MakeSource; using Rotation::MakeConsumer;
using Rotation::FaceOf; using Rotation::Neighbour; using Rotation::BuildPartGraph; using Rotation::kHandCrank; using Rotation::kMillstone; using Rotation::kSaw; using Rotation::kPress;

Rotation::ConsumerSpec RotationPlot_Consumer(PartType t) {
    switch (t) {
    case PartType::Millstone: return kMillstone;
    case PartType::Saw: return kSaw;
    case PartType::Press: return kPress;
    case PartType::Dynamo: return Rotation::kDynamo; // (#1271 L3)
    case PartType::SteamLab: return Rotation::kSteamLab; // (#1271 L6)
    case PartType::Bellows: return Rotation::kBellows; // (#1271 L2)
    default: return { 0, 0, 0 };
    }
}

int RotationPlot_OrientFor(Int3 to) {
    for (int i = 0; i < 24; ++i) { const int* M = Orientation(i); if (M[2] == to.x && M[5] == to.y && M[8] == to.z) return i; }
    return 0;
}

Part RotationPlot_PartAt(PartType t, Int3 cell, int orient) {
    Part p{ t, { 0, 0, 0 }, (uint8_t)orient };
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(p, fp);
    Int3 lo = fp.empty() ? Int3{ 0, 0, 0 } : fp[0].first;
    for (const auto& c : fp) lo = { std::min(lo.x, c.first.x), std::min(lo.y, c.first.y), std::min(lo.z, c.first.z) };
    p.origin = { cell.x - lo.x, cell.y - lo.y, cell.z - lo.z };
    return p;
}

bool RotationPlot_Describe(const Part& part, int partIndex, RotationPartMsg& out) {
    if (!IsRotationPart(part.type)) return false;
    std::vector<std::pair<Int3, CellShape>> fp;
    Plot::Footprint(part, fp);
    if (fp.empty()) return false;
    std::vector<Cell> cells;
    for (const auto& c : fp) cells.push_back({ c.first.x, c.first.y, c.first.z });
    // The local axis, turned by the part's orientation (row-major 3x3).
    const Int3 local = (part.type == PartType::Millstone || part.type == PartType::SteamLab) ? Int3{ 0, 1, 0 } : (part.type == PartType::Saw || part.type == PartType::Press || part.type == PartType::SteamEngine || part.type == PartType::Dynamo || part.type == PartType::Motor) ? Int3{ 1, 0, 0 } : Int3{ 0, 0, 1 };
    const int* M = Orientation(part.orient);
    const int w[3] = { M[0] * local.x + M[1] * local.y + M[2] * local.z, M[3] * local.x + M[4] * local.y + M[5] * local.z, M[6] * local.x + M[7] * local.y + M[8] * local.z };
    const int axis = w[0] != 0 ? 0 : w[1] != 0 ? 1 : 2;
    out = RotationPartMsg{};
    out.part = partIndex;
    out.type = part.type;
    out.axisSign = w[axis] < 0 ? -1 : 1;
    const Cell c0 = cells[0];
    switch (part.type) {
    case PartType::Shaft: case PartType::StressGauge: out.geom = MakeShaft(c0, axis); break;
    case PartType::Cogwheel: out.geom = MakeCog(c0, axis, false); break;
    case PartType::LargeCogwheel: out.geom = MakeCog(c0, axis, true); break;
    case PartType::Gearbox: out.geom = MakeGearbox(c0, axis); break;
    case PartType::HandCrank: out.geom = MakeSource({ c0 }, axis, { { c0, FaceOf(axis, +1) }, { c0, FaceOf(axis, -1) } }, kHandCrank, 0); break;
    case PartType::Windmill: out.geom = MakeSource({ c0 }, axis, { { c0, FaceOf(axis, +1) }, { c0, FaceOf(axis, -1) } }, Rotation::kWindmill, 0); break; // (#1270)
    case PartType::WaterWheel: { // (#1270) shafts join at its hub (footprint cell 4), both sides
        const Cell hub = cells.size() > 4 ? cells[4] : c0;
        out.geom = MakeSource(cells, axis, { { hub, FaceOf(axis, +1) }, { hub, FaceOf(axis, -1) } }, Rotation::kWaterWheel, 0);
        break;
    }
    default: { // a machine: every outer face of its cells along its axis
        std::set<std::tuple<int, int, int>> in;
        for (const Cell& c : cells) in.insert({ c.x, c.y, c.z });
        std::vector<std::pair<Cell, int>> ports;
        for (const Cell& c : cells)
            for (int s : { -1, 1 }) {
                const Cell n = Neighbour(c, FaceOf(axis, s));
                if (!in.count({ n.x, n.y, n.z })) ports.push_back({ c, FaceOf(axis, s) });
            }
        if (part.type == PartType::SteamEngine) out.geom = MakeSource(cells, axis, ports, Rotation::kSteamEngine, 0); // (#1271 L3) its speed from its boiler's steam
        else if (part.type == PartType::Motor) out.geom = MakeSource(cells, axis, ports, Rotation::kMotor, 0); // its speed from the EU it gets
        else out.geom = MakeConsumer(cells, axis, ports, RotationPlot_Consumer(part.type), false);
        break;
    }
    }
    return true;
}

void RotationPlot::Set(std::vector<RotationPartMsg> msgs) {
    std::unordered_map<int, int> held; // cranks still turning keep turning
    for (size_t i = 0; i < parts.size(); ++i) if (i < crankTicks.size() && crankTicks[i] > 0) held[parts[i].part] = crankTicks[i];
    std::unordered_map<int, int> sampled; // (#1270) wheels and windmills keep their last sampled speed
    for (size_t i = 0; i < parts.size(); ++i) if (i < sourceRpm.size() && sourceRpm[i] != 0) sampled[parts[i].part] = sourceRpm[i];
    parts = std::move(msgs);
    byPart.clear();
    std::vector<Rotation::Part> geoms;
    geoms.reserve(parts.size());
    for (size_t i = 0; i < parts.size(); ++i) { byPart[parts[i].part] = (int)i; geoms.push_back(parts[i].geom); }
    pg = BuildPartGraph(geoms);
    layout = Rotation::Build(pg.graph);
    crankTicks.assign(parts.size(), 0);
    for (size_t i = 0; i < parts.size(); ++i) if (auto it = held.find(parts[i].part); it != held.end() && parts[i].type == PartType::HandCrank) crankTicks[i] = it->second;
    sourceRpm.assign(parts.size(), 0);
    for (size_t i = 0; i < parts.size(); ++i)
        if (auto it = sampled.find(parts[i].part); it != sampled.end() && RotationPlot_FedSource(parts[i].type)) sourceRpm[i] = it->second;
    res = Result{};
    Rotation::Step(pg.graph, layout, res);
}

namespace {
int NodeOf(const RotationPlot& rp, int i) { return rp.pg.firstNode[i] + rp.parts[i].geom.mainNode; }
}

void RotationPlot::Hold(int part) {
    auto it = byPart.find(part);
    if (it != byPart.end() && parts[it->second].type == PartType::HandCrank) crankTicks[it->second] = kCrankHoldTicks;
}

void RotationPlot::SetRunning(int part, bool running) {
    auto it = byPart.find(part);
    if (it == byPart.end()) return;
    Node& n = pg.graph.nodes[NodeOf(*this, it->second)];
    if (n.kind == NodeKind::Consumer) n.running = running;
}

bool RotationPlot::NetHasLoad(int part) const {
    auto it = byPart.find(part);
    if (it == byPart.end()) return false;
    const int n = NodeOf(*this, it->second);
    if (n >= (int)layout.netOf.size()) return false;
    for (int v : layout.nets[layout.netOf[n]].nodes) {
        const Node& node = pg.graph.nodes[v];
        if (node.kind == NodeKind::Consumer && node.running) return true;
    }
    return false;
}

void RotationPlot::SetSourceRpm(int part, int rpm) {
    auto it = byPart.find(part);
    if (it == byPart.end()) return;
    const PartType t = parts[it->second].type;
    if (!RotationPlot_FedSource(t)) return;
    const int cap = t == PartType::WaterWheel ? Rotation::kWaterWheel.rpm : t == PartType::Windmill ? Rotation::kWindmill.rpm : t == PartType::SteamEngine ? Rotation::kSteamEngine.rpm : Rotation::kMotor.rpm;
    sourceRpm[it->second] = std::clamp(rpm, -cap, cap);
}

void RotationPlot::Step(int dt) {
    for (size_t i = 0; i < parts.size(); ++i) {
        if (RotationPlot_FedSource(parts[i].type)) { // (#1270) as the main thread sampled it (#1271 L3: or the factory set it: steam engines, motors)
            pg.graph.nodes[NodeOf(*this, (int)i)].sourceRpm = sourceRpm[i];
            continue;
        }
        if (parts[i].type != PartType::HandCrank) continue;
        pg.graph.nodes[NodeOf(*this, (int)i)].sourceRpm = crankTicks[i] > 0 ? parts[i].axisSign * kHandCrank.rpm : 0;
        crankTicks[i] = std::max(0, crankTicks[i] - dt);
    }
    Rotation::Step(pg.graph, layout, res);
}

int64_t RotationPlot::SpeedMilli(int part) const {
    auto it = byPart.find(part);
    if (it == byPart.end()) return 0;
    const int n = NodeOf(*this, it->second);
    return n < (int)res.speedMilli.size() ? res.speedMilli[n] : 0;
}

Status RotationPlot::StatusOf(int part) const {
    auto it = byPart.find(part);
    if (it == byPart.end()) return Status::Stopped;
    const int n = NodeOf(*this, it->second);
    return n < (int)res.status.size() ? res.status[n] : Status::Stopped;
}

const NetResult* RotationPlot::NetOf(int part) const {
    auto it = byPart.find(part);
    if (it == byPart.end()) return nullptr;
    const int n = NodeOf(*this, it->second);
    if (n >= (int)layout.netOf.size()) return nullptr;
    const int net = layout.netOf[n];
    return net >= 0 && net < (int)res.nets.size() ? &res.nets[net] : nullptr;
}

int RotationPlot::ProgressPermille(int part) const {
    auto it = byPart.find(part);
    if (it == byPart.end()) return 0;
    return Rotation::ProgressPermille(RotationPlot_Consumer(parts[it->second].type), SpeedMilli(part), StatusOf(part));
}

uint64_t RotationPlot::Hash() const {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](long long v) { h = (h ^ (uint64_t)v) * 1099511628211ull; };
    for (int64_t s : res.speedMilli) mix(s);
    for (Status s : res.status) mix((int)s);
    for (int c : crankTicks) mix(c);
    for (int r : sourceRpm) mix(r);
    return h;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1269): plot parts described as rotation parts, oriented, and stepped as the factory does ----
#include <gtest/gtest.h>

namespace {
int OrientFor(Int3 zTo) { return RotationPlot_OrientFor(zTo); }
// The part whose footprint starts at `at` (a turned part's origin is elsewhere: Footprint turns about the origin corner).
RotationPartMsg Describe(PartType t, Int3 at, int orient = 0, int index = 0) {
    RotationPartMsg m;
    EXPECT_TRUE(RotationPlot_Describe(RotationPlot_PartAt(t, at, orient), index, m));
    return m;
}
}

TEST(RotationPlot, ATurnedPartLandsOnTheCellAskedFor) {
    for (int o = 0; o < 24; ++o) {
        const RotationPartMsg m = Describe(PartType::Shaft, { 4, 5, 6 }, o);
        ASSERT_EQ(m.geom.cells.size(), 1u);
        EXPECT_TRUE(m.geom.cells[0] == (Rotation::Cell{ 4, 5, 6 })) << "orientation " << o;
    }
}

TEST(RotationPlot, PartsLieOnTheirTurnedAxis) {
    const RotationPartMsg z = Describe(PartType::Shaft, { 3, 1, 2 });
    EXPECT_EQ(z.geom.axis, 2);
    EXPECT_EQ(z.axisSign, 1);
    const RotationPartMsg x = Describe(PartType::Shaft, { 0, 0, 0 }, OrientFor({ -1, 0, 0 }));
    EXPECT_EQ(x.geom.axis, 0);
    EXPECT_EQ(x.axisSign, -1);
    const RotationPartMsg mill = Describe(PartType::Millstone, { 0, 0, 0 });
    EXPECT_EQ(mill.geom.axis, 1);
    EXPECT_EQ(mill.geom.cells.size(), 8u);
    EXPECT_EQ(mill.geom.ports.size(), 8u) << "four cells on top (+y), four underneath (-y)";
    EXPECT_EQ(Describe(PartType::Press, { 0, 0, 0 }).geom.axis, 0);
    RotationPartMsg none;
    EXPECT_FALSE(RotationPlot_Describe(Part{ PartType::Cube, { 0, 0, 0 }, 0 }, 0, none));
}

TEST(RotationPlot, ACrankTurnsAMillstoneWhileHeld) {
    const int up = OrientFor({ 0, 1, 0 });
    std::vector<RotationPartMsg> msgs = {
        Describe(PartType::Millstone, { 0, 0, 0 }, 0, 10),
        Describe(PartType::Shaft, { 0, 2, 0 }, up, 11),
        Describe(PartType::HandCrank, { 0, 3, 0 }, up, 12),
    };
    RotationPlot rp;
    rp.Set(msgs);
    rp.SetRunning(10, true);
    rp.Step(1);
    EXPECT_EQ(rp.SpeedMilli(10), 0) << "nobody holds the crank";
    rp.Hold(12);
    rp.Step(1);
    EXPECT_EQ(std::abs(rp.SpeedMilli(10)), 32000);
    EXPECT_EQ(rp.StatusOf(10), Status::Running);
    EXPECT_EQ(rp.ProgressPermille(10), 1000);
    ASSERT_NE(rp.NetOf(11), nullptr);
    EXPECT_EQ(rp.NetOf(11)->capacity, 256);
    for (int t = 0; t < kCrankHoldTicks; ++t) rp.Step(1);
    EXPECT_EQ(rp.SpeedMilli(10), 0) << "let go: it stops";
    // A layout change while held keeps the crank turning.
    rp.Hold(12);
    rp.Set(msgs);
    rp.SetRunning(10, true);
    rp.Step(1);
    EXPECT_EQ(std::abs(rp.SpeedMilli(10)), 32000);
}

TEST(RotationPlot, FiveMillstonesOnOneCrankAreTooSlow) {
    // A horizontal crank line along x with millstones under it would not join (they turn about y); a vertical shaft down
    // through a stack of millstones does: each passes rotation through its top and bottom faces.
    const int up = OrientFor({ 0, 1, 0 });
    std::vector<RotationPartMsg> msgs;
    for (int k = 0; k < 5; ++k) msgs.push_back(Describe(PartType::Millstone, { 0, 2 * k, 0 }, 0, k));
    msgs.push_back(Describe(PartType::HandCrank, { 0, 10, 0 }, up, 9));
    RotationPlot rp;
    rp.Set(msgs);
    for (int k = 0; k < 5; ++k) rp.SetRunning(k, true);
    rp.Hold(9);
    rp.Step(1);
    ASSERT_NE(rp.NetOf(0), nullptr);
    EXPECT_EQ(rp.NetOf(0)->ratio, 400);                // 256 SU against 5 x 128
    EXPECT_EQ(std::abs(rp.SpeedMilli(0)), 12800);
    EXPECT_EQ(rp.StatusOf(0), Status::TooSlow);
    EXPECT_EQ(rp.ProgressPermille(0), 0);
    rp.SetRunning(3, false); rp.SetRunning(4, false);  // two go idle: 3 x 128 = 384 on 256
    rp.Hold(9);
    rp.Step(1);
    EXPECT_EQ(rp.StatusOf(0), Status::Slowed);
    EXPECT_EQ(rp.ProgressPermille(0), 666);
}

TEST(RotationPlot, AWaterWheelAndAWindmillTurnAtTheirSampledSpeed) {
    std::vector<RotationPartMsg> msgs = {
        Describe(PartType::WaterWheel, { 0, 0, 0 }, 0, 1),  // 3 x 3 cells in z = 0, its hub at (1, 1, 0)
        Describe(PartType::Shaft, { 1, 1, 1 }, 0, 2),       // on its hub, along z
        Describe(PartType::Windmill, { 5, 5, 5 }, 0, 3),
        Describe(PartType::Shaft, { 5, 5, 4 }, 0, 4),
    };
    EXPECT_EQ(msgs[0].geom.cells.size(), 9u);
    EXPECT_EQ(msgs[0].geom.ports.size(), 2u) << "shafts join only at the hub";
    RotationPlot rp;
    rp.Set(msgs);
    rp.Step(1);
    EXPECT_EQ(rp.SpeedMilli(2), 0) << "nothing sampled yet: still water, calm";
    rp.SetSourceRpm(1, 12); rp.SetSourceRpm(3, -16);
    rp.Step(1);
    EXPECT_EQ(rp.SpeedMilli(2), 12000);
    ASSERT_NE(rp.NetOf(2), nullptr);
    EXPECT_EQ(rp.NetOf(2)->capacity, Rotation::kWaterWheel.suPerRpm * 12);
    EXPECT_EQ(rp.SpeedMilli(4), -16000);
    EXPECT_EQ(rp.NetOf(4)->capacity, 1024) << "four sails of 256 SU at full wind";
    rp.SetSourceRpm(1, 40);                              // a flood: capped at the wheel's top speed
    rp.Step(1);
    EXPECT_EQ(rp.SpeedMilli(2), Rotation::kWaterWheel.rpm * 1000);
    rp.Set(msgs);                                        // a layout change keeps the sampled speeds
    rp.Step(1);
    EXPECT_EQ(rp.SpeedMilli(2), Rotation::kWaterWheel.rpm * 1000);
    rp.SetSourceRpm(2, 30);                              // a shaft is no source: ignored
    rp.Step(1);
    EXPECT_EQ(rp.SpeedMilli(2), Rotation::kWaterWheel.rpm * 1000);
}
#endif
