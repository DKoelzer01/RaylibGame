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
#include <vector>
#include <raylib.h>

class RiverNetwork {
public:
    // height(dir): smooth land height (m above sea level) for a unit direction.
    RiverNetwork(float radius, unsigned int seed, const std::function<float(Vector3)>& height);

    // The river closest to a point (body-local unit direction), within `reach` metres of its centreline.
    struct Hit {
        bool found = false;
        float dist = 1e30f;          // m from the centreline
        float level = 0.0f;          // water surface (m above sea level)
        float halfWidth = 0.0f, depth = 0.0f;
    };
    Hit nearest(Vector3 dir, float reach) const;
    // Lakes: a whole filled basin of the flood (every cell of it, not only the deep ones) whose level
    // is held below the ground sampled all round it, so the shore is wherever the terrain rises past
    // the level and never a cell's straight edge. The water region is the basin's cells plus
    // kLakeMargin around them. Here: the lake's level and how far outside its cells (0 inside, m);
    // reported out to 2 * kLakeMargin (the bank fades out there), so water only where outside < kLakeMargin.
    bool lakeAt(Vector3 dir, float& level, float& outside) const;
    bool lakeAt(Vector3 dir, float& level) const { float o; return lakeAt(dir, level, o) && o < kLakeMargin; } // lake water here
    bool lakeNear(Vector3 dir, float extent) const;  // cheap: a lake's water region within `extent` m may be here
    bool riverNear(Vector3 dir, float extent) const; // cheap: a river within `extent` m may be here
    bool waterNear(Vector3 dir, float extent) const { return lakeNear(dir, extent) || riverNear(dir, extent); }
    static constexpr float kLakeMargin = 60.0f;

    // Tests and stats.
    struct Stats { int cells = 0, landCells = 0, riverCells = 0, lakeCells = 0, mouths = 0; float cellSize = 0.0f, longestKm = 0.0f; };
    Stats stats() const { return st; }
    // Every river link: the level never rises and the width never shrinks going downstream; returns
    // the number of links checked and the worst rise (m) / shrink (m) seen.
    int checkDownstream(float& worstRise, float& worstShrink) const;
    // Points along the longest river from its source to the sea (body-local unit directions), every
    // `step` metres, with the level the network gives there.
    void longestPath(float step, std::vector<Vector3>& dirs, std::vector<float>& levels) const;

    static constexpr float kMaxReach = 300.0f; // valley carving reaches at most this far from a centreline

private:
    float R = 1.0f, cell = 420.0f;
    int n = 1; // cells per face edge
    std::vector<float> h, F, L, hw, dep;
    std::vector<int> parent, mainChild;
    std::vector<float> acc;
    std::vector<unsigned char> river, lake;
    std::vector<float> lakeLevel;          // per lake cell: its basin's water level (m)
    std::vector<unsigned char> lakeSteps;  // cells to the nearest lake cell (capped at 255)
    struct Piece { Vector3 a, b; float la, lb, wa, wb, da, db; }; // positions in metres (body-local)
    std::vector<Piece> pieces;
    std::vector<std::vector<int>> index; // cell -> pieces within reach
    Stats st;

    int cellOf(Vector3 dir) const;
    Vector3 centre(int c) const;
    void neighbours(int c, int out[8]) const;
    float distToCell(Vector3 dir, int c) const; // m from a direction to cell c's square (0 inside)
};

#endif
