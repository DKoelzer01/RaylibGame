#ifndef DEPTHTABLE_H
#define DEPTHTABLE_H

// #1952: the one underground depth table. Every depth is metres under the ground surface above (the planet's elevation
// there), unless named otherwise. Caves (voxel.cpp CaveCell), underground regions (underground.cpp) and ore (voxel.cpp
// oreKind) read their numbers here, so the depth bands (#2070), heat with depth (#1472) and the ore gradient (the user's
// ruling on #1952: pockets thin out with depth, deposits fewer but bigger and richer, ore kinds per band) change in one
// place. Creatures keep their own depth= ranges in creatures.txt (data, one per creature).
//
// Band edges: the user's ruling (#1465, 2026-11-27): shallow 25-250 m, middle 250-500 m, deep 500 m to the mantle. Heat
// with depth starts at the middle band (#1472). The other values are GUESSES.

enum class DepthBand { Surface = -1, Shallow = 0, Middle, Deep };

struct DepthTable {
    // Bands (#2070): shallow from the first caves down to `middle`, middle down to `deep`, deep down to the mantle (caves
    // stop `mantleClearance` above it; the mantle itself is unmineable).
    float shallow = 25.0f, middle = 250.0f, deep = 500.0f;
    float mantleClearance = 20.0f;   // m above the mantle radius (not under the ground)
    float mantleFade = 10.0f;        // m: caves fade out over this above the clearance
    // The noise caverns (#999) in the density field: from noiseCavernMin under the ground (faded in over noiseCavernFade)
    // down to noiseCavernMantle above the mantle (faded over noiseCavernFade). The GPU density (terrain.glsl) mirrors
    // these and the two above as literals: --test-terrain-gpu (bit-exact) fails if they part.
    float noiseCavernMin = 30.0f, noiseCavernFade = 20.0f, noiseCavernMantle = 40.0f;

    // Cave passages (#615, #999): a passage's axis stays at least passageMin + passagePerRadius x its radius deep (so it
    // never reaches the surface) and at most passageMax (today's clamp, 100 m short of the shallow band's floor; #2070
    // lifts it to the mantle); networks start passageMin +
    // passageStartSpan x U deep.
    float passageMin = 25.0f, passagePerRadius = 1.3f, passageMax = 150.0f, passageStartSpan = 95.0f;
    float chamberMin = 40.0f;        // #1058: chambers at passage vertices this deep
    float roomRoofClear = 10.0f;     // chambers, caverns, sheets: skipped where their dome comes this close to the surface
    float cavernMin = 100.0f;        // #1553
    float sheetMin = 30.0f, sheetMax = 120.0f; // #1565

    // Underground regions (#1547, GUESS thresholds).
    float floodedSeaMax = 200.0f;    // under the sea: flooded caves above this
    float frozenMax = 150.0f;        // the cold reaches this deep
    float crystalMin = 150.0f;       // rare crystal pockets below this
    float floodedWetMax = 100.0f;    // very wet ground floods its caves above this
    float fungalMin = 30.0f, fungalMax = 250.0f;

    // Ore (#504, #598). Each ore's and vein family's own depth range stays in its data (ores.txt minDepth, the families'
    // min/max depth, #1466, #1467).
    float oreMin = 3.0f;             // no ore this close to the surface
    float pocketFadeStart = 20.0f, pocketFadeSpan = 200.0f, pocketFadeThreshold = 0.12f; // pockets thin out with depth

    DepthBand bandAt(float depth) const {
        return depth < shallow ? DepthBand::Surface : depth < middle ? DepthBand::Shallow : depth < deep ? DepthBand::Middle : DepthBand::Deep;
    }
};

inline constexpr DepthTable kDepth{};

#endif
