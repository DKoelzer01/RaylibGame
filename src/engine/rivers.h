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
    // A lake here (in a lake cell or next to one): its level, else false.
    bool lakeAt(Vector3 dir, float& level) const;
    bool waterNear(Vector3 dir, float extent) const; // cheap: a river or lake within `extent` m may be here

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
    struct Piece { Vector3 a, b; float la, lb, wa, wb, da, db; }; // positions in metres (body-local)
    std::vector<Piece> pieces;
    std::vector<std::vector<int>> index; // cell -> pieces within reach
    Stats st;

    int cellOf(Vector3 dir) const;
    Vector3 centre(int c) const;
    void neighbours(int c, int out[8]) const;
};

#endif
