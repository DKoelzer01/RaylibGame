#ifndef UNDERGROUND_H
#define UNDERGROUND_H

// (#1547, the underground plan #1465) The underground region map: each 256 m cave cell (the cave graph's cell) gets an
// underground biome from the surface above it, its depth and the seed. Read later by cave generation (ice, crystal and
// fungal caves: material, light, spawns) and by ore families (#1467). CPU-side, cached per cell and thread. All the
// thresholds are GUESS values (tuned when the first consumer lands).

#include <raylib.h>
#include <cstdint>

class PlanetField;

enum class UndergroundRegion : uint8_t { None = 0, Karst, Volcanic, Frozen, Crystal, Flooded, Fungal, Count };
const char* UndergroundRegion_Name(UndergroundRegion r);

// What the classification reads (the cell centre's): depth under the surface above it (m), the rock province (strata.h:
// 0 sedimentary, 1 granite, 2 basalt / volcanic, 3 karst), the surface climate and biome there, the cell's hash 0..1.
struct UndergroundInputs {
    float depth = 0.0f;
    int province = 0;
    float temp = 15.0f, moist = 0.5f;
    int biome = 0;       // Biome (voxel.h)
    bool ocean = false;  // under the sea
    float roll = 0.0f;   // the cell's own random 0..1 (crystal pockets)
};
UndergroundRegion UndergroundRegion_Classify(const UndergroundInputs& in);

// The region at a body-local point: its 256 m cell's (cached).
UndergroundRegion Underground_RegionAt(const PlanetField& pf, Vector3 p);
UndergroundInputs Underground_InputsAt(const PlanetField& pf, Vector3 p); // the cell's inputs (benches, tests)

#endif
