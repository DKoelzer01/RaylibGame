#ifndef ROTATIONPARTS_H
#define ROTATIONPARTS_H

// Rotation part geometry (#159, #1269): turns a plot's rotation parts (shafts, cogwheels, gearboxes, speed changers, sources
// and consumers, each with its cells and axis) into the Graph the rotation solver runs on (rotation.h). No PartType here:
// the factory maps its part types and orientations to these descriptions (step 2), so this stays testable on its own.
//
// Conventions:
// - A node's speed is signed about the POSITIVE global direction of its axis (x, y or z).
// - A port is a face of one of a part's cells through which a shaft end passes. Two ports join when they face each other
//   from adjacent cells (A's cell + its face direction = B's cell, B's face the opposite); both express speed about the
//   same global axis, so their speeds are equal (ratio = A's factor x B's factor, factors are +-1).
// - Cogwheels also mesh, Create-style, with cogs of the same axis: two small cogs side by side in the plane across the
//   axis turn opposite ways at the same speed; a small cog diagonal to a large one (one cell off in both directions across
//   the axis) drives it at half its speed, the other way. A large cog does not mesh with a large one, nor a small cog
//   with a large one beside it (only diagonally).
// - A gearbox (Create's) joins shafts on the four faces around its axis: opposite faces turn opposite ways about their
//   shared axis; faces at right angles follow the bevel rule (both positive or both negative directions: reversed; mixed:
//   the same).
// - A speed changer (the doc's "gearbox" that trades speed for stress) takes its input on its -axis face and gives
//   num/den of that speed on its +axis face: two graph nodes and a link between them.

#include "rotation.h"
#include <vector>

namespace Rotation {

struct Cell { int x = 0, y = 0, z = 0; };
inline bool operator==(const Cell& a, const Cell& b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
// Faces: 0 +x, 1 -x, 2 +y, 3 -y, 4 +z, 5 -z. Axis of a face = face / 2; its sign = face % 2 ? -1 : +1.
constexpr int FaceOf(int axis, int sign) { return axis * 2 + (sign < 0 ? 1 : 0); }
constexpr int FaceAxis(int face) { return face / 2; }
constexpr int FaceSign(int face) { return face % 2 ? -1 : 1; }
constexpr int Opposite(int face) { return face ^ 1; }
Cell Neighbour(Cell c, int face);

enum class PartKind : uint8_t { Shaft, Cog, LargeCog, Gearbox, SpeedChanger, Source, Consumer };

struct Port { Cell cell; int face = 0; int node = 0; int factor = 1; }; // node: index into the part's own nodes
struct Part {
    PartKind kind = PartKind::Shaft;
    int axis = 0;                      // 0 x, 1 y, 2 z
    std::vector<Cell> cells;
    std::vector<Port> ports;
    std::vector<Node> nodes;           // 1, or 2 for a speed changer (input, output)
    std::vector<Link> innerLinks;      // between the part's own nodes (local indices)
    // The part's speed is reported from this node (a speed changer: its output).
    int mainNode = 0;
};

// Builders, one per kind. Cells are plot cells; a part along its axis spans `length` cells from `at` in the + direction.
Part MakeShaft(Cell at, int axis, int length = 1);
Part MakeCog(Cell at, int axis, bool large);
Part MakeGearbox(Cell at, int axis);
Part MakeSpeedChanger(Cell at, int axis, int num, int den);
// A source or consumer occupying `cells`, taking or giving rotation through the given shaft faces (`ports`: cell + face;
// their factor is +1, so the node turns about its axis as the shafts on those faces do). Its axis is the ports' axis.
Part MakeSource(const std::vector<Cell>& cells, int axis, const std::vector<std::pair<Cell, int>>& ports, SourceSpec spec, int rpm);
Part MakeConsumer(const std::vector<Cell>& cells, int axis, const std::vector<std::pair<Cell, int>>& ports, ConsumerSpec spec, bool running = true);

struct PartGraph {
    Graph graph;
    std::vector<int> firstNode;        // per part: its first node in graph.nodes (its nodes are consecutive)
    int MainNode(const std::vector<Part>& parts, int part) const { return firstNode[part] + parts[part].mainNode; }
};
// The graph of a plot's rotation parts: their nodes in part order, their inner links, then the port joins and the cog
// meshes. The same parts in any order give the same speeds (the solver is order-independent).
PartGraph BuildPartGraph(const std::vector<Part>& parts);

} // namespace Rotation

#endif
