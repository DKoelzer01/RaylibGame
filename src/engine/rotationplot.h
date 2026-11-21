#ifndef ROTATIONPLOT_H
#define ROTATIONPLOT_H

// Rotation on a plot (#159, #1269): the bridge between plot parts (PartType, orientation) and the rotation solver
// (rotation.h, rotationparts.h). The main thread describes each rotation part when a plot's layout changes; the factory
// keeps one RotationPlot per simulated plot and steps it each tick before power. Everything here is per plot (no shared
// state), so plots tick on any worker in any order with the same result.
//
// Local axes (before the part's orientation): shafts, cogs, the gearbox, the hand crank and the stress gauge turn about
// their local z (a beam's length); the millstone about its local y (a vertical stone); the saw and the press about their
// local x. A machine takes rotation through every outer face of its cells along that axis.

#include "plots.h"
#include "rotationparts.h"
#include <cstdint>
#include <unordered_map>
#include <vector>

// A part as the factory gets it from the layout (built on the main thread).
struct RotationPartMsg {
    int part = -1;                 // plot part index
    PartType type = PartType::Count;
    int axisSign = 1;              // its local axis points along + (1) or - (-1) of the plot axis it lies on
    Rotation::Part geom;
};
// false: not a rotation part.
bool RotationPlot_Describe(const Part& part, int partIndex, RotationPartMsg& out);
// The machine's load (stress per RPM, minimum and rated speed); {0,0,0} for parts that are not rotation consumers.
Rotation::ConsumerSpec RotationPlot_Consumer(PartType t);
// An orientation that turns a part's local z (a shaft's axis) onto `to` (a unit plot axis, either sign); 0 if none.
int RotationPlot_OrientFor(Int3 to);
// The part of that type and orientation whose footprint starts at `cell` (Plot::Footprint turns about the origin corner, so
// a turned one-cell part's origin is not its cell).
Part RotationPlot_PartAt(PartType t, Int3 cell, int orient);
constexpr int kCrankHoldTicks = 6; // a crank turns this long (0.3 s) after the last "held" from the player

struct RotationPlot {
    std::vector<RotationPartMsg> parts;          // in part order
    std::unordered_map<int, int> byPart;         // plot part -> index in parts
    Rotation::PartGraph pg;
    Rotation::Layout layout;
    Rotation::Result res;
    std::vector<int> crankTicks;                 // per part: ticks left turning (hand cranks)
    bool empty() const { return parts.empty(); }
    // A new layout: rebuilds the graph; cranks still held keep turning.
    void Set(std::vector<RotationPartMsg> msgs);
    void Hold(int part);                         // a player turns this crank (for kCrankHoldTicks)
    void SetRunning(int part, bool running);     // a consumer has a recipe this tick (only running ones load the net)
    void Step(int dt);                           // dt ticks: cranks count down, then one solve
    bool Has(int part) const { return byPart.count(part) > 0; }
    int64_t SpeedMilli(int part) const;          // signed milli-RPM now (0 if unknown)
    Rotation::Status StatusOf(int part) const;
    const Rotation::NetResult* NetOf(int part) const;
    int ProgressPermille(int part) const;        // a consumer's work rate this tick (permille of its rated speed)
    uint64_t Hash() const;                       // speeds, statuses and crank state (the determinism test)
};

#endif
