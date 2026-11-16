#ifndef STRATA_H
#define STRATA_H

// Rock strata (Plan 9 step 3, #587, folds in #581). One data source for the rock under a place: the erosion bake reads
// its hardness and erodibility (cliff bands where a resistant layer crops out, caprock mesas, buttes, talus), the mesher
// its surface material where rock shows (cliffs, scree, cave walls), and digging its item (with its gem).
//
// A column (StrataColumn, from region-scale fields only, so the bake's cells and the mesher's vertices agree) holds the
// province and how its layers lie; Strata_RockAt gives the rock at a height. Heights are in the bake's simulation units
// (the bake rescales its result by PlanetErosion::heightScale: divide a final height by it first).
#include <raylib.h>
#include <cstdint>

enum class RockType : uint8_t { Soft = 0, Granite, Basalt, Sandstone, Limestone, Marble, Count };

struct StrataColumn {
    uint8_t province = 0;   // 0 layered sedimentary rock, 1 a granite core (mountains), 2 basalt (sea floor, volcanic)
    uint8_t hardShare = 115; // of 255: the share of resistant layers (sedimentary / basalt flows)
    float offset = 0.0f;    // m: the layers' warp and tilt here (added to the height before layering)
    float thick = 80.0f;    // m: layer thickness
};

// Province and layer geometry from the field's continent value at this place (c: < 0 ocean, > 0.3 well inland) and the
// surface point s (m). `seed` = the planet's terrain seed.
StrataColumn Strata_Column(Vector3 s, float continent, unsigned int seed);
RockType Strata_RockAt(const StrataColumn& c, float h, unsigned int seed);
float Strata_Hardness(RockType r);     // 0 soft .. 1 hardest: the bake's talus limit
float Strata_Erodibility(RockType r);  // the stream-power multiplier (soft ~1.5, granite 0.15)
const char* Strata_Item(RockType r);   // its items.txt id ("stone" for soft rock)
// Fault scarps (#590): metres a few long curved faults raise the land at surface point s (0 off them; the bake input).
float Strata_FaultUplift(Vector3 s, unsigned int seed);

#endif
