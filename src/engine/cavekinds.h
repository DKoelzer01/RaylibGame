#ifndef CAVEKINDS_H
#define CAVEKINDS_H

// (#1465 infrastructure, step 2) Cave segment kinds as a table: a segment's kind (s[15]) picks one row, and the CPU
// (voxel.cpp CaveGraphDistance) and the GPU (terrain.glsl caveGraphDistance) read the same parameters, so a new kind
// (sheet, cavern, ice, crystal) is a row here plus its segments, and goes through the same wall-noise code (the user's
// rule: every cave uses today's walls). The GLSL copy is generated from this table when a compute shader is loaded
// (CaveKinds_Glsl, the `#pragma cave_kinds` line in terrain.glsl), so the two can't drift.
//
// Today's four rows reproduce the old inline rules exactly (bit-identical density, --test-terrain-gpu, --cave-stats hash).

#include <string>

enum class CaveKind : int { Passage = 0, Mouth = 1, Solid = 2, Chamber = 3, Cavern = 4, Sheet = 5 };

struct CaveKindDef {
    const char* name;
    int solid;          // 1: a solid inside (aprons, lips, boulders, columns): solid = max(solid, r - distance), no bore
    float outline;      // the large-scale outline noise's weight on the bore (passages 0.2; mouths' noise-cut 0.4)
    int lopsided;       // 1: taller on one side (mouths)
    float rubble;       // rubble mounds on the floor, x max(0, rubble noise) (mouths 0.35)
    float fixedLevel;   // >= 0: the floor level as a share of the bore instead of the wandering noise level (mouths 0.45)
    float floorR;       // > 0: a level floor this share of the radius down (chambers 0.4: boulders sit on it), else level x bore
    int tilt;           // 1: the floor tilts with the tilt noise (passages); 0: level
    int scaleToPassage; // 1: swell and outline scaled by min(1, 8 m / r): at a passage's size on a big bore (chambers, #1116)
    float intrude;      // the wall noise's reach as a share of the bore (0.55; mouths 0.65)
    float outcropCap;   // > 0: wall outcrops at most this many metres (chambers 4, #1116: floating rock), else the bore's
    int wander;         // 1: the roof pockets and the sideways wander scaled by min(1, 8 m / r) too (caverns: at full size they
                        //    fold a big dome and leave floating rock, --cave-islands)
};

constexpr CaveKindDef kCaveKinds[] = {
    //  name       solid outline lop rubble level floorR tilt scale intrude cap   wander
    { "passage",   0,    0.2f,   0,  0.0f,  -1.0f, 0.0f,  1,   0,    0.55f,  0.0f, 0 },
    { "mouth",     0,    0.4f,   1,  0.35f, 0.45f, 0.0f,  0,   0,    0.65f,  0.0f, 0 },
    { "solid",     1,    0.2f,   0,  0.0f,  -1.0f, 0.0f,  0,   0,    0.55f,  0.0f, 0 },
    { "chamber",   0,    0.2f,   0,  0.0f,  -1.0f, 0.4f,  0,   1,    0.55f,  4.0f, 0 },
    // (#1465) a large cavern (r 50-150 m, rare and deep): the chamber's rules exactly (swell and outline at a passage's
    // size, a level floor, outcrops capped at 4 m: the floating-rock control), and its roof pockets and wander at that size too
    { "cavern",    0,    0.2f,   0,  0.0f,  -1.0f, 0.4f,  0,   1,    0.55f,  4.0f, 1 },
    // (#1465) a sheet cave: a wide, low bedding-plane passage (half-width 12-25 m, height x0.12-0.2) in layered rock: a level
    // floor at the noise level (no tilt: along the bedding), swell, outline, pockets and wander at a passage's size, outcrops
    // capped at 3 m (they stay in the low roof)
    { "sheet",     0,    0.2f,   0,  0.0f,  -1.0f, 0.0f,  0,   1,    0.55f,  3.0f, 1 },
};
constexpr int kCaveKindCount = (int)(sizeof(kCaveKinds) / sizeof(kCaveKinds[0]));

// A segment's stored kind (a float, s[15]) -> its row: the nearest whole kind, clamped (past the last: the last, as the
// old `kind > 2.5` read anything above as a chamber; below 0: a passage).
inline int CaveKind_Index(float k) {
    const int i = (int)(k + 0.5f);
    return k < 0.5f ? 0 : i >= kCaveKindCount ? kCaveKindCount - 1 : i;
}

// The table as GLSL const arrays (kCK_* indexed by kind), for terrain.glsl's `#pragma cave_kinds` line.
std::string CaveKinds_Glsl();

#endif
