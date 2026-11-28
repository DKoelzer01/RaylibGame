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
#include <vector>

enum class RockType : uint8_t { Soft = 0, Granite, Basalt, Sandstone, Limestone, Marble, Count };

struct StrataColumn {
    uint8_t province = 0;   // 0 layered sedimentary rock, 1 a granite core (mountains), 2 basalt (sea floor, volcanic), 3 karst limestone
    uint8_t hardShare = 115; // of 255: the share of resistant layers (sedimentary / basalt flows)
    float offset = 0.0f;    // m: the layers' warp and tilt here (added to the height before layering)
    float thick = 80.0f;    // m: layer thickness
    // (#1990) The biome's bake multipliers here (BiomeVariantTable::BakeMul; 1 = as before): only the erosion bake reads
    // them, from its own strata callback (PlanetField::bakeStrataColumn).
    float kMul = 1.0f, talusMul = 1.0f, creepMul = 1.0f, glacialMul = 1.0f;
    float pinch = 0.5f;     // (#2064) 0..1 a lateral value (~3.5 km): lenticular beds pinch out where it is under their own threshold
    float caprock = 0.0f;   // (#1990 step 4) 0..1 aridity: resistant beds thicken over the soft ones (mesas and buttes)
};

// (#2064) The bed at a height: its index (which picks its rock), where in it the height lies (0 at its base, 1 at its
// top), its thickness (m) and whether it is a resistant one. Beds come in sets of four nominal layers whose thicknesses
// are skewed (many thin beds, a few thick ones), some lenticular; RAYLIBGAME_STRATA_BEDS=0: the even layers of #587.
struct StrataBed { int32_t index; float frac, thick; bool hard; };
StrataBed Strata_BedAt(const StrataColumn& c, float h, unsigned int seed);

// Province and layer geometry from the field's continent value at this place (c: < 0 ocean, > 0.3 well inland) and the
// surface point s (m). `seed` = the planet's terrain seed.
StrataColumn Strata_Column(Vector3 s, float continent, unsigned int seed);
RockType Strata_RockAt(const StrataColumn& c, float h, unsigned int seed);
float Strata_Hardness(RockType r);     // 0 soft .. 1 hardest: the bake's talus limit
float Strata_Erodibility(RockType r);  // the stream-power multiplier (soft ~1.5, granite 0.15)
const char* Strata_Item(RockType r);   // its items.txt id ("stone" for soft rock)
// Fault scarps (#590): metres a few long curved faults raise the land at surface point s (0 off them; the bake input).
float Strata_FaultUplift(Vector3 s, unsigned int seed);
// The fault field itself (#1248: Earth ley lines follow the faults): the signed value whose zero crossing is the fault line
// (> 0 on the raised side; ~4e-5 per metre across it), and the region weight (0 outside the rare fault regions).
float Strata_FaultValue(Vector3 s, unsigned int seed, float* region);
// Folded strata (#591): 0..1, the weight of the fold regions (cuestas, hogbacks) at surface point s.
float Strata_FoldWeight(Vector3 s, unsigned int seed);
// #591: 0..1 across a resistant sedimentary layer at height h (0 in soft beds): a ridge profile for the bake input.
float Strata_OutcropBump(const StrataColumn& c, float h, unsigned int seed);
// Karst (#592): 0..1 weight of the karst regions (~10 km), and the towers (+) and dolines (-) there, in metres.
float Strata_KarstWeight(Vector3 s, unsigned int seed);
float Strata_KarstRelief(Vector3 s, unsigned int seed);
// #1052: the dolines whose surface centre (on the sphere of radius R) lies in the box [lo, hi): centre direction, bowl
// radius and depth there (before the karst weight), and a stable id.
struct KarstDoline { Vector3 dir; float radius, depth; uint32_t id; };
void Strata_KarstDolines(Vector3 lo, Vector3 hi, float R, unsigned int seed, std::vector<KarstDoline>& out);
// Volcanoes (#593): metres a volcano adds at surface point s (cones 1-3 km, craters/calderas), 0 off them; weight =
// 0..1 proximity to the vent (the cone's basalt-and-ash province).
float Strata_Volcano(Vector3 s, float planetRadius, unsigned int seed, float* weight = nullptr);
inline bool gStrataFolds = true; // #591; RAYLIBGAME_FOLDS=0 for A/B (its own bake cache)
inline bool gStrataBeds = true;  // #2064: skewed, lenticular beds and fault throw; RAYLIBGAME_STRATA_BEDS=0 for A/B (its own bake cache)
inline bool gStrataClimate = true; // #1990 step 4: caprock where dry, karst where humid; RAYLIBGAME_STRATA_CLIMATE=0 for A/B

#endif
