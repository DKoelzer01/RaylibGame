#ifndef CAVEISLANDS_H
#define CAVEISLANDS_H

// (#1116, #1465) Floating rock in cave chambers: for each chamber within R m of the start, the density on a 1 m grid over a
// box +-50 m across and +-28 m up around its centre (local frame); solid cells connected (6-neighbour) to the box's boundary
// are the surrounding rock, any other solid component floats. The user's rule (#1465): no cave kind may add floating rock,
// so the caves/islands scenario holds the total to today's budget. `--cave-islands R` prints it (print = true).

class PlanetField;

struct CaveIslands {
    int chambers = 0;
    long long pieces = 0, volume = 0;  // floating components and their m^3, over all chambers
    long long byHeight[8] = {};        // m^3 by mean height from the chamber centre: -28..+28 m in 8 bins
    long long byDistance[8] = {};      // m^3 by mean distance across from it: 0..72 m in 8 bins
    float seconds = 0.0f;
};

CaveIslands CaveIslands_Measure(const PlanetField& pf, unsigned seed, float R, bool print);

#endif
