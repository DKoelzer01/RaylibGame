#ifndef WORLDGENPARAMS_H
#define WORLDGENPARAMS_H

// World generation parameters (#504): chosen on the New World screen, saved with the world in a sidecar
// (<save>.worldgen), applied to the start planet's PlanetField before its terrain is made. The defaults are the
// terrain as it was before (every scale 1, no shift): an old save without a sidecar generates as before.
#include <cstdint>
#include <string>

struct WorldGenParams {
    float ocean = 0.0f;           // -1..1: less .. more ocean (shifts the continents' threshold)
    float mountains = 1.0f;       // 0.3..2: mountain range height (the ridges, and the erosion bake's uplift)
    float roughness = 1.0f;       // 0.3..2: hills and fine detail
    float biomeSize = 1.0f;       // 0.5..2: the size of climate regions (biomes)
    float resources = 1.0f;       // 0.5..2: ore pockets and deposits
    bool erosion = true;          // the erosion bake (Plan 7; off: the noise terrain)
    float erosionStrength = 1.0f; // 0.5..2: how hard rivers cut (the bake's stream power)
    float cliffs = 1.0f;          // 0..2: cliffs and canyons (#589): the share of resistant rock layers (0: none; Plan 9)
    bool karst = true;            // karst regions: limestone towers and dolines (#592)
    bool volcanoes = true;        // volcanoes: cones, calderas, lava caps (#593)

    bool isDefault() const;
    uint64_t key() const;         // a hash of every field (the erosion bake's cache key includes it)
};

// The parameters for the world at `savePath`: its sidecar, or the defaults when there is none (old saves, benches).
WorldGenParams WorldGenParams_ForSave(const std::string& savePath);
// The sidecar's text ("key value" lines; unknown keys are skipped, missing ones keep their defaults). The multiplayer
// Welcome carries the same text (#757), so a new key here travels to clients with no protocol change.
std::string WorldGenParams_ToText(const WorldGenParams& p);
WorldGenParams WorldGenParams_FromText(const std::string& text);
// Writes the sidecar (only the New World screen does: benches reuse their save and never write one).
bool WorldGenParams_Save(const std::string& savePath, const WorldGenParams& p);

// The start planet of a world seed in the current system definition (worldGen.system): its terrain seed and radius,
// drawn exactly as generateWorld draws them. False when the system has no planet.
bool WorldGen_StartPlanet(unsigned int worldSeed, unsigned int& terrainSeed, float& radius);

#endif
