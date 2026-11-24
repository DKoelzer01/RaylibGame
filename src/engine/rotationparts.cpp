#include "rotationparts.h"
#include <cstddef>
#include <map>
#include <tuple>

// Rotation part geometry (#1269). See rotationparts.h.

namespace Rotation {

Cell Neighbour(Cell c, int face) {
    const int s = FaceSign(face);
    switch (FaceAxis(face)) {
    case 0: c.x += s; break;
    case 1: c.y += s; break;
    default: c.z += s; break;
    }
    return c;
}

namespace {
Cell Along(Cell c, int axis, int k) { if (axis == 0) c.x += k; else if (axis == 1) c.y += k; else c.z += k; return c; }
}

Part MakeShaft(Cell at, int axis, int length) {
    Part p;
    p.kind = PartKind::Shaft; p.axis = axis;
    if (length < 1) length = 1;
    for (int k = 0; k < length; ++k) p.cells.push_back(Along(at, axis, k));
    p.nodes.push_back(Node{});
    p.ports.push_back({ p.cells.front(), FaceOf(axis, -1), 0, 1 });
    p.ports.push_back({ p.cells.back(), FaceOf(axis, +1), 0, 1 });
    return p;
}

Part MakeCog(Cell at, int axis, bool large) {
    Part p = MakeShaft(at, axis, 1); // a cog sits on a shaft through it
    p.kind = large ? PartKind::LargeCog : PartKind::Cog;
    return p;
}

Part MakeGearbox(Cell at, int axis) {
    Part p;
    p.kind = PartKind::Gearbox; p.axis = axis;
    p.cells.push_back(at);
    p.nodes.push_back(Node{});
    const int pAxis = axis == 0 ? 1 : 0; // the lower of the two axes across its own
    const int canon = FaceOf(pAxis, +1);
    for (int face = 0; face < 6; ++face) {
        if (FaceAxis(face) == axis) continue;
        int factor;
        if (face == canon) factor = 1;
        else if (face == Opposite(canon)) factor = -1;             // opposite faces turn opposite ways
        else factor = FaceSign(face) > 0 ? -1 : 1;                 // bevel rule against the canonical (+) face
        p.ports.push_back({ at, face, 0, factor });
    }
    return p;
}

Part MakeSpeedChanger(Cell at, int axis, int num, int den) {
    Part p;
    p.kind = PartKind::SpeedChanger; p.axis = axis;
    p.cells.push_back(at);
    p.nodes.resize(2);
    p.innerLinks.push_back({ 0, 1, num, den });
    p.ports.push_back({ at, FaceOf(axis, -1), 0, 1 });
    p.ports.push_back({ at, FaceOf(axis, +1), 1, 1 });
    p.mainNode = 1;
    return p;
}

Part MakeSource(const std::vector<Cell>& cells, int axis, const std::vector<std::pair<Cell, int>>& ports, SourceSpec spec, int rpm) {
    Part p;
    p.kind = PartKind::Source; p.axis = axis; p.cells = cells;
    Node n; n.kind = NodeKind::Source; n.sourceRpm = rpm; n.suPerRpm = spec.suPerRpm;
    p.nodes.push_back(n);
    for (const auto& [c, face] : ports) p.ports.push_back({ c, face, 0, 1 });
    return p;
}

Part MakeConsumer(const std::vector<Cell>& cells, int axis, const std::vector<std::pair<Cell, int>>& ports, ConsumerSpec spec, bool running) {
    Part p;
    p.kind = PartKind::Consumer; p.axis = axis; p.cells = cells;
    Node n; n.kind = NodeKind::Consumer; n.stressPerRpm = spec.stressPerRpm; n.minRpm = spec.minRpm; n.running = running;
    p.nodes.push_back(n);
    for (const auto& [c, face] : ports) p.ports.push_back({ c, face, 0, 1 });
    return p;
}

PartGraph BuildPartGraph(const std::vector<Part>& parts) {
    PartGraph out;
    Graph& g = out.graph;
    out.firstNode.resize(parts.size());
    for (size_t i = 0; i < parts.size(); ++i) {
        out.firstNode[i] = (int)g.nodes.size();
        for (const Node& n : parts[i].nodes) g.nodes.push_back(n);
        for (const Link& k : parts[i].innerLinks) g.link(out.firstNode[i] + k.a, out.firstNode[i] + k.b, k.num, k.den);
    }
    // Ports: by (cell, face). A port joins the port facing it from the next cell.
    struct PortRef { int part, port; };
    std::map<std::tuple<int, int, int, int>, PortRef> byFace;
    for (size_t i = 0; i < parts.size(); ++i)
        for (size_t j = 0; j < parts[i].ports.size(); ++j) {
            const Port& pt = parts[i].ports[j];
            byFace.emplace(std::make_tuple(pt.cell.x, pt.cell.y, pt.cell.z, pt.face), PortRef{ (int)i, (int)j }); // first wins
        }
    for (const auto& [key, ref] : byFace) {
        const Port& a = parts[ref.part].ports[ref.port];
        const Cell nc = Neighbour(a.cell, a.face);
        auto it = byFace.find(std::make_tuple(nc.x, nc.y, nc.z, Opposite(a.face)));
        if (it == byFace.end() || it->second.part == ref.part) continue;
        if (std::make_tuple(nc.x, nc.y, nc.z, Opposite(a.face)) < key) continue; // each pair once
        const Port& b = parts[it->second.part].ports[it->second.port];
        g.link(out.firstNode[ref.part] + a.node, out.firstNode[it->second.part] + b.node, a.factor * b.factor, 1);
    }
    // Cog meshes, per axis: small-small side by side (reversed, 1:1), small-large diagonal (reversed, large at half speed).
    std::map<std::tuple<int, int, int, int>, int> cogAt; // (axis, cell) -> part
    for (size_t i = 0; i < parts.size(); ++i)
        if (parts[i].kind == PartKind::Cog || parts[i].kind == PartKind::LargeCog)
            cogAt.emplace(std::make_tuple(parts[i].axis, parts[i].cells[0].x, parts[i].cells[0].y, parts[i].cells[0].z), (int)i);
    for (const auto& [key, i] : cogAt) {
        if (parts[i].kind != PartKind::Cog) continue; // meshes are found from the small cog's side
        const int axis = parts[i].axis;
        const int p = axis == 0 ? 1 : 0, q = axis == 2 ? 1 : 2;
        const Cell c = parts[i].cells[0];
        const int self = out.firstNode[i] + parts[i].mainNode;
        auto other = [&](Cell at) -> int {
            auto it = cogAt.find(std::make_tuple(axis, at.x, at.y, at.z));
            return it == cogAt.end() ? -1 : it->second;
        };
        // Side by side: only the + neighbours, so each small pair links once.
        for (int ax : { p, q }) {
            const int j = other(Along(c, ax, 1));
            if (j >= 0 && parts[j].kind == PartKind::Cog) g.link(self, out.firstNode[j] + parts[j].mainNode, -1, 1);
        }
        for (int sp : { -1, 1 })
            for (int sq : { -1, 1 }) {
                const int j = other(Along(Along(c, p, sp), q, sq));
                if (j >= 0 && parts[j].kind == PartKind::LargeCog) g.link(self, out.firstNode[j] + parts[j].mainNode, -1, 2);
            }
    }
    return out;
}

} // namespace Rotation

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1269): shaft lines, cog meshes, gearboxes and speed changers turned into links, solved end to end ----
#include <gtest/gtest.h>

using namespace Rotation;

namespace {
struct Solved { PartGraph pg; Result r; };
Solved SolveParts(const std::vector<Part>& parts) { Solved s; s.pg = BuildPartGraph(parts); Step(s.pg.graph, Build(s.pg.graph), s.r); return s; }
int64_t SpeedOf(const Solved& s, const std::vector<Part>& parts, int part) { return s.r.speedMilli[s.pg.MainNode(parts, part)]; }
// A one-cell crank on the x axis at `at` (shafts on both its x faces), and a one-cell millstone the same way.
Part Crank(Cell at, int axis, int rpm) { return MakeSource({ at }, axis, { { at, FaceOf(axis, +1) }, { at, FaceOf(axis, -1) } }, kHandCrank, rpm); }
Part Mill(Cell at, int axis, bool running = true) { return MakeConsumer({ at }, axis, { { at, FaceOf(axis, +1) }, { at, FaceOf(axis, -1) } }, kMillstone, running); }
}

TEST(RotationParts, AShaftLineCarriesTheCrankToTheMillstone) {
    std::vector<Part> parts = { Crank({ 0, 0, 0 }, 0, 32), MakeShaft({ 1, 0, 0 }, 0, 3), Mill({ 4, 0, 0 }, 0) };
    const Solved s = SolveParts(parts);
    EXPECT_EQ(s.r.nets.size(), 1u);
    EXPECT_EQ(SpeedOf(s, parts, 2), 32000);
    EXPECT_EQ(s.r.status[s.pg.MainNode(parts, 2)], Status::Running);
}

TEST(RotationParts, AGapOrATurnedShaftBreaksTheLine) {
    std::vector<Part> gap = { Crank({ 0, 0, 0 }, 0, 32), MakeShaft({ 2, 0, 0 }, 0), Mill({ 3, 0, 0 }, 0) };
    EXPECT_EQ(SpeedOf(SolveParts(gap), gap, 2), 0);
    std::vector<Part> turned = { Crank({ 0, 0, 0 }, 0, 32), MakeShaft({ 1, 0, 0 }, 2), Mill({ 2, 0, 0 }, 0) }; // a z shaft between
    EXPECT_EQ(SpeedOf(SolveParts(turned), turned, 2), 0);
    std::vector<Part> beside = { Crank({ 0, 0, 0 }, 0, 32), Mill({ 0, 1, 0 }, 0) };                             // side by side, not end to end
    EXPECT_EQ(SpeedOf(SolveParts(beside), beside, 1), 0);
}

TEST(RotationParts, SmallCogsSideBySideTurnOppositeWays) {
    std::vector<Part> parts = { Crank({ 0, 0, 0 }, 2, 32), MakeCog({ 0, 0, 1 }, 2, false), MakeCog({ 1, 0, 1 }, 2, false), MakeShaft({ 1, 0, 2 }, 2) };
    const Solved s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 1), 32000);
    EXPECT_EQ(SpeedOf(s, parts, 2), -32000);
    EXPECT_EQ(SpeedOf(s, parts, 3), -32000);
    // A cog on another axis next to it does not mesh.
    std::vector<Part> skew = { Crank({ 0, 0, 0 }, 2, 32), MakeCog({ 0, 0, 1 }, 2, false), MakeCog({ 1, 0, 1 }, 0, false) };
    EXPECT_EQ(SpeedOf(SolveParts(skew), skew, 2), 0);
}

TEST(RotationParts, ASmallCogDrivesALargeOneDiagonallyAtHalfSpeed) {
    std::vector<Part> parts = { Crank({ 0, 0, 0 }, 1, 32), MakeCog({ 0, 1, 0 }, 1, false), MakeCog({ 1, 1, 1 }, 1, true) };
    Solved s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 2), -16000);
    parts[2] = MakeCog({ 1, 1, 0 }, 1, true); // beside it: no mesh
    s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 2), 0);
    parts[2] = MakeCog({ 1, 1, 1 }, 1, false); // a small cog diagonal: no mesh either
    s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 2), 0);
}

TEST(RotationParts, AGearboxTurnsShaftsOnItsFourSides) {
    // A y-axis gearbox at the origin, driven from +x; shafts on -x, +z and -z.
    std::vector<Part> parts = {
        Crank({ 2, 0, 0 }, 0, 32), MakeShaft({ 1, 0, 0 }, 0), MakeGearbox({ 0, 0, 0 }, 1),
        MakeShaft({ -1, 0, 0 }, 0), MakeShaft({ 0, 0, 1 }, 2), MakeShaft({ 0, 0, -1 }, 2),
        MakeShaft({ 0, 1, 0 }, 1), // on its own axis: not joined
    };
    const Solved s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 1), 32000);
    EXPECT_EQ(SpeedOf(s, parts, 3), -32000); // opposite face: reversed about x
    EXPECT_EQ(SpeedOf(s, parts, 4), -32000); // +x to +z: reversed
    EXPECT_EQ(SpeedOf(s, parts, 5), 32000);  // +x to -z: the same
    EXPECT_EQ(SpeedOf(s, parts, 6), 0);
}

TEST(RotationParts, FourGearboxesInARingAgree) {
    std::vector<Part> parts = {
        MakeGearbox({ 0, 0, 0 }, 1), MakeShaft({ 1, 0, 0 }, 0), MakeGearbox({ 2, 0, 0 }, 1), MakeShaft({ 2, 0, 1 }, 2),
        MakeGearbox({ 2, 0, 2 }, 1), MakeShaft({ 1, 0, 2 }, 0), MakeGearbox({ 0, 0, 2 }, 1), MakeShaft({ 0, 0, 1 }, 2),
        Crank({ 0, 1, 0 }, 1, 32), // drives nothing: a gearbox has no port on its own axis
        Crank({ 0, 0, -1 }, 2, 32), // on the first gearbox's -z face
    };
    const Solved s = SolveParts(parts);
    ASSERT_EQ(s.r.nets[Build(s.pg.graph).netOf[s.pg.MainNode(parts, 0)]].jam, Jam::None) << "the ring's bevels agree";
    EXPECT_EQ(SpeedOf(s, parts, 0), 32000);  // -z port: the same as the gearbox's +x
    EXPECT_EQ(SpeedOf(s, parts, 1), 32000);
    EXPECT_EQ(SpeedOf(s, parts, 3), 32000);
    EXPECT_EQ(SpeedOf(s, parts, 5), -32000);
    EXPECT_EQ(SpeedOf(s, parts, 7), -32000);
}

TEST(RotationParts, ASpeedChangerGearsUpOneWayAndDownTheOther) {
    std::vector<Part> parts = { Crank({ 0, 0, 0 }, 0, 32), MakeSpeedChanger({ 1, 0, 0 }, 0, 2, 1), Mill({ 2, 0, 0 }, 0) };
    Solved s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 2), 64000);
    // Turned round (the crank on its output side): the crank's side is the output, so the millstone gets half.
    parts = { Mill({ 0, 0, 0 }, 0), MakeSpeedChanger({ 1, 0, 0 }, 0, 2, 1), Crank({ 2, 0, 0 }, 0, 32) };
    s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 0), 16000);
}

TEST(RotationParts, MultiCellMachinesTakeRotationThroughTheirPortsOnly) {
    // A 2x2x2 press taking rotation from above at one cell's +y face.
    std::vector<Cell> cells;
    for (int x = 0; x < 2; ++x) for (int y = 0; y < 2; ++y) for (int z = 0; z < 2; ++z) cells.push_back({ x, y, z });
    std::vector<Part> parts = {
        MakeConsumer(cells, 1, { { { 0, 1, 0 }, FaceOf(1, +1) } }, kPress),
        MakeShaft({ 0, 2, 0 }, 1, 2), Crank({ 0, 4, 0 }, 1, 32),
        MakeShaft({ 1, 2, 1 }, 1), // over another cell of the press: no port there
    };
    const Solved s = SolveParts(parts);
    EXPECT_EQ(SpeedOf(s, parts, 0), 32000);
    EXPECT_EQ(s.r.nets[Build(s.pg.graph).netOf[s.pg.MainNode(parts, 0)]].demandMilli, 8 * 32000);
    EXPECT_EQ(SpeedOf(s, parts, 3), 0);
}

TEST(RotationParts, PartOrderDoesNotMatter) {
    std::vector<Part> parts = {
        Crank({ 2, 0, 0 }, 0, 32), MakeShaft({ 1, 0, 0 }, 0), MakeGearbox({ 0, 0, 0 }, 1), MakeShaft({ 0, 0, 1 }, 2),
        MakeCog({ 0, 0, 2 }, 2, false), MakeCog({ 1, 1, 2 }, 2, true), MakeSpeedChanger({ 1, 1, 3 }, 2, 1, 2), Mill({ 1, 1, 4 }, 2),
        Mill({ -1, 0, 0 }, 0),
    };
    const Solved a = SolveParts(parts);
    std::vector<Part> rev(parts.rbegin(), parts.rend());
    const Solved b = SolveParts(rev);
    const int n = (int)parts.size();
    for (int i = 0; i < n; ++i) EXPECT_EQ(SpeedOf(a, parts, i), SpeedOf(b, rev, n - 1 - i)) << "part " << i;
    EXPECT_EQ(SpeedOf(a, parts, 7), 8000); // 32 -> (gearbox +x to +z: -32 about z) -> cog -> large +16 -> /2 = 8
}
#endif
