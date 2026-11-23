#ifndef RIVERS_H
#define RIVERS_H

// River drainage network (plan 2, milestone 14). Built once per planet from a coarse height grid
// (~420 m cells on a cube-sphere): a priority flood from the sea gives every cell a downstream
// neighbour (depressions fill up to their spill level: lakes), flow accumulates downstream, and cells
// draining more than a threshold area carry a river. Water levels are a running minimum going
// downstream (never rising), widths and depths grow with the flow, every river ends in the sea (lakes
// on the way). Each river link is a curve (Catmull-Rom through the cell centres plus a meander)
// split into short pieces kept in a per-cell index, so terrain queries are cheap.

#include <functional>
#include <unordered_map>
#include <vector>
#include <atomic>
#include <memory>
#include <mutex>
#include <raylib.h>

class PlanetErosion;

class RiverNetwork {
public:
    // The network's grid and rivers. With the erosion bake (Plan 7 step 5, #359): the bake's 200 m cells, rivers from
    // ~2 km^2 of drainage (upland streams), little meander (the baked valleys already wind), and fineHeight: each
    // link's level is also held under the fine ground sampled along its curve (no bank fill holds it up any more).
    struct Params {
        float cell = 420.0f;            // m, target cell size
        float riverArea = 7.06e6f;      // m^2 drained before a river starts
        float meander = 0.18f;          // sideways swing, fraction of a link's length (at most 70 m)
        std::function<float(Vector3)> fineHeight; // optional: the full ground height (m), sampled along the links
        std::function<float(Vector3)> bankHeight; // optional (bake-native, #640): the ground beside a river (m), levels are held under it
        float maxReach = 300.0f;        // m: the valley carving's reach from a centreline (the index covers it)
        float lakeDepth = 6.0f;         // m: a filled depression deeper than this holds water (a lake, or with the bake a pond)
        std::function<float(Vector3)> runoff; // optional (#637): water per km^2 drained at a place, x the default (moisture x the water cycle)
        float valleyFill = 1.0f;        // (#637) 0..1: how much of its worn valley floor a river fills (0: the discharge width alone; the water cycle sets it)
        // optional (bake-native, #903): how far the runtime moves the bake's ground at a place of bake height h (m: biomes'
        // regional reshaping after the bake, e.g. swamps pulled down to a floor 2-5 m over the sea). Levels move with it.
        std::function<float(Vector3 dir, float h)> regionalShift;
    };
    // height(dir): smooth land height (m above sea level) for a unit direction.
    RiverNetwork(float radius, unsigned int seed, const std::function<float(Vector3)>& height, const Params& params);
    RiverNetwork(float radius, unsigned int seed, const std::function<float(Vector3)>& height);
    // Bake-native rivers (#547, #607-#609): the rivers ARE the bake's final drainage (PlanetErosion::hasFlow) - every
    // 50 m cell draining params.riverArea or more is a node, linked to the cell it drains to; the water level is the
    // bake's filled surface (a lake's level in its basins, else the worn channel's ground plus most of its incision), so
    // it never rises downstream and meets the river it joins at that river's level. No stepped levels, no valley carve.
    RiverNetwork(const PlanetErosion& bake, float radius, const Params& params);
    bool bakeNative() const { return bake != nullptr; }
    // Bake-native: the ground under the channels - h lowered to each piece's bed within its half width and to a 0.6 bank
    // slope just beyond (the lowest of all pieces near, so overlapping ends never leave a ridge, #608). Lakes carve nothing.
    float carve(Vector3 dir, float h) const;
    // Bake-native: the river head (a source) nearest to dir with water over minLevel, and its downstream direction (#625 looks).
    bool nearestHead(Vector3 dir, float minLevel, Vector3& head, Vector3& down) const;

    // The river closest to a point (body-local unit direction), within `reach` metres of its centreline.
    struct Hit {
        bool found = false;
        float dist = 1e30f;          // m from the centreline
        float level = 0.0f;          // water surface (m above sea level)
        float halfWidth = 0.0f, depth = 0.0f;
        Vector3 down{};              // unit tangent (body-local) pointing downstream (#361)
        float slope = 0.0f;          // the reach's fall per metre (its link, ~200 m: flow speed for the bed, #595)
        int piece = -1;              // the link it lies on (river flow state, #656)
    };
    Hit nearest(Vector3 dir, float reach) const;
    // Lakes: a whole filled basin of the flood (every cell of it, not only the deep ones) whose level
    // is held below the ground sampled all round it, so the shore is wherever the terrain rises past
    // the level and never a cell's straight edge. The water region is the basin's cells plus
    // kLakeMargin around them. Here: the lake's level and how far outside its cells (0 inside, m);
    // reported out to 2 * kLakeMargin (the bank fades out there), so water only where outside < kLakeMargin.
    bool lakeAt(Vector3 dir, float& level, float& outside) const;
    bool lakeAt(Vector3 dir, float& level, float& outside, int& id) const; // also the lake's id (lakeInfo)
    // Each lake's extent (#360: the small ones are finite, drained by digs): the centre, a radius covering its cells,
    // its level, depth (level - lowest cell) and size in cells.
    struct LakeInfo { Vector3 dir{ 0, 1, 0 }; float r = 0.0f, level = 0.0f, depth = 0.0f; int cells = 0; };
    const LakeInfo& lakeInfo(int id) const { return lakes[id]; }
    float cellSize() const { return cell; }
    bool lakeAt(Vector3 dir, float& level) const { float o; return lakeAt(dir, level, o) && o < kLakeMargin; } // lake water here
    bool lakeNear(Vector3 dir, float extent) const;  // cheap: a lake's water region within `extent` m may be here
    bool riverNear(Vector3 dir, float extent) const; // cheap: a river within `extent` m may be here
    bool waterNear(Vector3 dir, float extent) const { return lakeNear(dir, extent) || riverNear(dir, extent); }
    static constexpr float kLakeMargin = 60.0f;

    // Tests and stats.
    struct Stats { int cells = 0, landCells = 0, riverCells = 0, lakeCells = 0, mouths = 0, riverCellsAbove700 = 0, heads = 0, headFalls = 0; float cellSize = 0.0f, longestKm = 0.0f, highestSource = 0.0f; }; // heads: river sources; headFalls: those with a waterfall in their first 3 links (#625)
    Stats stats() const { return st; }
    // Every river link: the level never rises and the width never shrinks going downstream; returns
    // the number of links checked and the worst rise (m) / shrink (m) seen.
    int checkDownstream(float& worstRise, float& worstShrink) const;
    // Points along the longest river from its source to the sea (body-local unit directions), every
    // `step` metres, with the level the network gives there.
    void longestPath(float step, std::vector<Vector3>& dirs, std::vector<float>& levels) const;
    // Every river's centreline every `step` metres: the direction, the water level there (waterfalls step, #540) and the depth.
    void forEachCentrePoint(float step, const std::function<void(Vector3 dir, float level, float depth)>& fn) const;
    // Every link (tests, #638): its ends (body-local m), levels, full widths and depths at each end, and whether it lies in a lake.
    struct Link { Vector3 a, b; float la, lb, wa, wb, da, db; bool inLake; };
    void forEachLink(const std::function<void(const Link&)>& fn) const;
    // The waterfalls (#540): where each is (its lip), how far it drops (m), the downstream tangent, the river's half width.
    struct Fall { Vector3 lip{}, down{}; float drop = 0.0f, top = 0.0f, halfWidth = 0.0f; };
    std::vector<Fall> waterfalls() const;

    static constexpr float kMaxReach = 300.0f; // valley carving reaches at most this far from a centreline
    float reach() const { return reachM; }    // this network's reach (Params::maxReach)

    // River flow (#656, user: 'if a player fully dams up a river and diverts it elsewhere ... I want this to be possible'):
    // each link carries a share of its natural discharge (1: as generated). Less flow lowers and narrows its water (level
    // - depth x (1 - f^0.4), width x sqrt(f)); none leaves a dry bed (the channel stays: the terrain carve is unchanged).
    // A released reach's measured outflow sets the link at its downstream edge (setReachFlow); stepFlow carries changes
    // down the network at ~5 m/s, mixing at confluences by discharge; lakes keep their level (their whole basin feeds them).
    // Bake-native networks only. Lock-free for readers (nearest() reads a relaxed atomic per piece).
    float flowOf(int piece) const;
    void setReachFlow(int piece, float share) const; // share of the natural discharge leaving a released reach here
    void clearReachFlow(int piece) const;
    // The first link below `piece` (downstream) that starts more than `dist` m from `centre` (body-local m): a released
    // reach drives the flow from there, so the links under the reach keep their natural water (its sink classifies on them).
    int pieceBeyond(int piece, Vector3 centre, float dist) const;
    // While one lives on this thread, nearest() reports the natural river (ground materials: a dried bed keeps its stones).
    struct NaturalScope { NaturalScope(); ~NaturalScope(); bool prev; };
    // Advances the flow by dt seconds; appends the pieces whose water changed visibly since they were last reported.
    bool stepFlow(float dt, std::vector<int>& changed) const;
    struct PieceBounds { Vector3 mid; float radius; };
    PieceBounds pieceBounds(int piece) const;    // body-local centre and a radius covering its water (refreshes)
    // Save / load: the pieces off their natural flow and the reach overrides.
    void flowState(std::vector<std::pair<int, float>>& cur, std::vector<std::pair<int, float>>& set) const;
    void setFlowState(const std::vector<std::pair<int, float>>& cur, const std::vector<std::pair<int, float>>& set) const;

private:
    float R = 1.0f, cell = 420.0f, reachM = kMaxReach;
    int n = 1; // cells per face edge
    std::vector<float> h, F, L, hw, dep;
    std::vector<int> parent, mainChild;
    std::vector<float> acc;
    std::vector<unsigned char> river, lake;
    std::vector<float> lakeLevel;          // per lake cell: its basin's water level (m)
    std::vector<unsigned char> lakeSteps;  // cells to the nearest lake cell (capped at 255)
    std::vector<int> lakeIdOf;             // per cell: its lake's id, or -1
    std::vector<LakeInfo> lakes;
    struct Piece { Vector3 a, b; float la, lb, wa, wb, da, db; Vector3 mid{}; float halfLen = 0.0f; float fallT = -1.0f; float reachSlope = 0.0f; bool inLake = false; }; // positions in metres (body-local); mid/halfLen: a bounding sphere for nearest(); fallT >= 0: a waterfall at that fraction (la before it, lb after, #540)
    std::vector<Piece> pieces;
    std::vector<std::vector<int>> index; // cell -> pieces within reach
    Stats st;

    int cellOf(Vector3 dir) const;
    Vector3 centre(int c) const;
    void neighbours(int c, int out[8]) const;
    float distToCell(Vector3 dir, int c) const; // m from a direction to cell c's square (0 inside)
    // Bake-native (#547)
    const PlanetErosion* bake = nullptr;
    std::unordered_map<int, int> fineLake;  // bake cell -> lake id
    std::vector<int> nodeCell, nodeParent, nodePiece; // per node: its bake cell, the node it drains to (-1: the sea), its piece
    std::vector<float> nodeDisc; // per node: the discharge's half width (m), before the valley fill (#637); never shrinks downstream
    // River flow (#656): per node its discharge (as an area at the default runoff), the order sources first, and the state.
    std::vector<float> nodeQ, nodeChildQ; // discharge (as a drainage area), and its tributaries' together
    std::vector<int> nodeTopo;
    std::unique_ptr<std::atomic<float>[]> flowCur; // per piece (nodes and pieces share the index)
    mutable std::vector<float> flowSet, flowSent;       // override (< 0: none), the value last reported as changed
    mutable std::atomic<bool> flowAny{ false };         // some piece is off its natural flow (readers skip the lookup otherwise)
    mutable std::mutex flowMx;                     // writers (stepFlow, setReachFlow, setFlowState)
    bool fineLakeAt(Vector3 dir, float& level, float& outside, int& id) const;
};

#endif
