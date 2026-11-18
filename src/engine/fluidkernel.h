#ifndef FLUIDKERNEL_H
#define FLUIDKERNEL_H
#include "envflag.h"

// The fluid step on the CPU (#654): the same three passes as resources/compute/fluid_step.comp, over the same packed
// layout, so a machine without compute shaders (a headless dedicated server, FLUID_CPU=1) simulates water as the GPU does.
// Keep the two in step: a change to one is a change to the other.
//
// Cells: uint32 = fill (bits 0-7) | flags (8-15) | target (16-23) | ground share (24-31), slot * 4096 + x + 16 y + 256 z.
// Faces (#642): per cell P (vec4: momentum flows across +x / +y as applied, carried units) and N (two vec4: this substep's
// flows and carries, the momentum flows alone).

#include <cstdint>
#include <vector>
#include <cstdlib>
#include <atomic>

namespace FluidKernel {

struct Vec4 { float x = 0.0f, y = 0.0f, z = 0.0f, w = 0.0f; };

struct Params {
    int minSpread[8] = {}, visc[8] = {};
    float frictionLin[8] = {}, frictionQuad[8] = {};
    float dt = 1.0f / 120.0f;
    bool evaporate = false;
    // The cell edge (m, W4c #642): fill units per m^3 (units) and per metre of height (hu); the column depth looks this many
    // cells down (columnCells, ~8 m). Set with setCell.
    float cellSize = 1.0f, units = 255.0f, hu = 255.0f;
    int columnCells = 8;
    bool columns = EnvFlag("FLUID_COLUMNS"); // #680 prototype: lateral flow per pair of water columns
    bool advect = !EnvFlag("FLUID_NOADVECT"); // #680: column flows carry their momentum downstream
    void setCell(float cs, int column) { cellSize = cs; units = 255.0f / (cs * cs * cs); hu = 255.0f / cs; columnCells = column; }
};

struct Grid {
    std::vector<uint32_t> a, b;   // A (the state between substeps), B (after the vertical pass)
    std::vector<Vec4> faceP;      // a vec4 a cell
    std::vector<Vec4> faceN;      // two vec4 a cell
    std::vector<int> nbr;         // slot * 6 + face (-x +x -y +y -z +z), -1 none
    std::vector<uint8_t> processed; // per slot
    std::vector<uint32_t> flags;  // per slot: bit 0 changed, bits 1-6 face requests
    std::vector<int> delta;       // per cell: column-pair flow units this substep (#680; zeroed in the vertical pass)
    // The column pairs' limiter (#680): pass 1 plans each pair's units and adds them to its giver's and taker's demand;
    // the move pass scales every pair by min(1, avail / out, room / in), so no cell gives or takes more than it has.
    // Per cell: demand out / in, the smallest avail / room any pair saw (all reset in the vertical pass); per run bottom
    // cell and pair (x, y, z, w): the giver's top cell and the taker cell.
    std::vector<int> outDem, inDem, availOf, roomOf;
    std::vector<int> pairGive, pairTake;
    std::vector<uint32_t> sinkOut; // per slot: fill units its sinks took in (river flow, #656; the caller zeroes it)
    std::vector<uint8_t> dbgMark; // TEMP #680 (DBTRACE): per slot, trace the transfers into it
    void resize(int slots);       // keeps nothing
};

// One substep over the listed slots: vertical (A -> B), face flows (B -> N), apply (B, N -> A, P). Slots run in parallel
// on the job system (each pass writes only its own slot's cells, faces and flags).
void Substep(Grid& g, const Params& p, const std::vector<int>& list);

// TEMP #680 (DBTRACE): column transfers into the slots marked in Grid::dbgMark (cell = slot * 4096 + local).
struct DbgTransfer { int giver, taker, branch, units, takerShut, capped, runBottom; };
constexpr int kDbgTrace = 65536;
extern DbgTransfer dbgTrace[kDbgTrace];
struct DbgPairRec { int a, nb, dx, dy, why, depth; }; // a pair a marked run planned: its exit code (fluidkernel.cpp WHY), depth in units
extern DbgPairRec dbgPairs[kDbgTrace];
extern std::atomic<int> dbgPairsN;
extern std::atomic<int> dbgTraceN;

} // namespace FluidKernel

#endif
