#ifndef CANYONS_H
#define CANYONS_H

// Graded canyon floors (#544, #547; Plan 9): a pass before the erosion bake that makes every canyon line drain end to end.
//
// The canyon lines (PlanetField::canyonFault) are noise zero crossings: they cross ridges and valleys alike, so a notch a
// fixed depth under the land is a chain of closed basins, and the bake's rivers spilled out of it sideways. Here, on a
// 100 m cube-sphere grid: the corridor cells (notch >= 0.5) are grouped into canyons; each drains to the end whose
// surroundings are lowest; along it (the distance from that end, in 100 m bins) the floor is the closest non-decreasing
// profile to the notch's own floor (weighted isotonic regression: it cuts the humps and fills the dips, by the least
// squares), at least 1 m per km, and at least `minDepth` under the land. The bake input then takes that floor along the
// line, with walls as wide as the depth needs (PlanetField::land).
//
// Deterministic: per-cell samples on threads (each cell's own value), everything else serial in index order.

#include <raylib.h>
#include <cstdint>
#include <functional>
#include <vector>

class CanyonGrade {
public:
    // notch(dir): 0..1 along the canyon lines; height(dir): the bake input there without this grading (m above sea).
    void build(float planetRadius, float cell, const std::function<float(Vector3)>& notch, const std::function<float(Vector3)>& height, int threads);
    // The graded floor (m) at dir, and how much of it there is (0 outside the canyons and past their ends, 1 inside).
    float floorAt(Vector3 dir, float& coverage) const;
    bool empty() const { return floor.empty(); }

    // Stats (logged, tests).
    int canyons = 0, corridorCells = 0, divided = 0; // divided: canyons draining both ways from a divide (kMaxCut)
    float seconds = 0.0f, meanCut = 0.0f, maxCut = 0.0f, meanFill = 0.0f, maxFill = 0.0f, lengthKm = 0.0f;

private:
    int n = 0;          // cells per face side
    float R = 0.0f;
    std::vector<float> floor; // per cell: the graded floor, NaN where none (the corridors and a margin around them)
    int cellOf(Vector3 dir) const;
    Vector3 dirOf(int f, int i, int j) const;
};

#endif
