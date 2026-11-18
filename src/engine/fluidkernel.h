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
    // The layout is the GPU's (fluid_step.comp), integer atomics only, so CPU and GPU agree whatever the order:
    //   dem: per cell, demand out << 16 | in (sums stay under 2040: no carry); cap: avail << 16 | room (one value a cell,
    //   OR-ed); both reset in the vertical pass.
    //   pairCode: per run bottom cell and pair (x, y, z, w), 16 bits (the GPU packs a cell's four into two words): aGives |
    //   (giver row + 64) << 1 (8 bits: -64..191) | (taker row - giver row + 110) << 9 (7 bits: -110..17: mostly a drop, but
    //   a pair's momentum can carry water some rows uphill, a run-up), rows counted from the run's bottom row; the move pass
    //   walks back to the cells (giver in A's column when A gives, the neighbour's otherwise; the taker in the other). A
    //   pair out of that range plans no move.
    std::vector<uint32_t> dem, cap;
    std::vector<uint16_t> pairCode;
    std::vector<uint32_t> sinkOut; // per slot: fill units its sinks took in (river flow, #656; the caller zeroes it)
    std::vector<uint8_t> dbgMark; // TEMP #680 (DBTRACE): per slot, trace the transfers into it
    // The dry-cell skip (#747, as the GPU's #675 fluid_groups.comp): per slot, 3 rotating substeps x 64 group bits (the
    // 64-cell groups that held water: group = z * 4 + y / 4), and the groups the next substep runs. The host sets
    // groupsFull at each tick's start: slots are handed out afresh each tick, so its first substep runs every cell.
    std::vector<uint32_t> gmask;  // slot * 6 + substep * 2 + word
    std::vector<uint64_t> gon;    // per slot
    int gcur = 0;
    bool groupsFull = true;
    void resize(int slots);       // keeps nothing
};

// One substep over the listed slots: vertical (A -> B), face flows (B -> N), apply (B, N -> A, P). Slots run in parallel
// on the job system (each pass writes only its own slot's cells, faces and flags).
void Substep(Grid& g, const Params& p, const std::vector<int>& list);
// A tick's start (#747), after the host packed `slots` slots into A and B: clears the per-cell scratch the slots carry from
// other bricks (face records of the processed slots, every slot's column deltas) and seeds the dry-cell skip's masks
// from the cells, so the first substep skips dry groups like the rest. Runs the slots on the job system.
void BeginTick(Grid& g, const std::vector<int>& list, int slots);

// TEMP #747 (DBSTEPPROF): wall time per phase (ms: group selection, passes 0 / 1 / 3 / 2) and groups run / offered,
// summed over substeps; the host prints and zeroes them.
extern double dbgPhaseMs[5];
extern long long dbgGroupsRun, dbgGroupsAll, dbgSubsteps, dbgFullSubsteps;

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
