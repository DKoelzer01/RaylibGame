#ifndef OREVEINS_H
#define OREVEINS_H

// (#1467, the underground plan #1465) Vein generation: a 96 m deposit cell holds at most one deposit of an ore family
// (ore_families.txt), picked among the families whose hosts and depth band fit the deposit's centre, by rarity x biome x
// region; the deposit takes its family's shape, and each point inside it one of the family's ores by weight. Pure
// functions of their inputs (PlanetField::oreKind supplies the rock, depth, biome and noise); CPU only (ores don't carve).

#include "orefamilies.h"
#include <cstdint>
#include <functional>
#include <raylib.h>
#include <string>

constexpr float kVeinCell = 96.0f;  // m: a deposit cell (as the old deposit grid)
constexpr float kVeinReach = 22.0f; // m: every shape stays this close to its centre (the centre is >= 24 m from the cell's edges)

struct VeinDeposit {
    bool exists = false;
    Vector3 centre{};
    int family = -1;          // OreFamilies_All() index
    VeinShape shape = VeinShape::Blob;
    float size = 0.0f;        // m: blob or pipe radius, sheet or seam thickness, stockwork vein spacing
    Vector3 axis{ 0, 1, 0 };  // a sheet's normal; a pipe's and a seam's up (the radial direction)
    float envelope = 0.0f;    // m: the radius it fits in (miners' and the scanner's reach)
    uint32_t seed = 0;        // the member pick's
};

// What the deposit's centre sits in, for the family pick (PlanetField looks these up once per cell).
struct VeinSite { std::string rock = "rock", biome = "meadow", region = "none"; float depth = 0.0f; };

// The cell's deposit from its hash h (the old deposit's rolls in the same order: existence 35%, the centre; then the size,
// the family and the orientation). `site` is asked only for an existing deposit, with its centre.
VeinDeposit OreVeins_Make(int cx, int cy, int cz, uint32_t h, const std::function<VeinSite(Vector3 centre)>& site);
// The family among those that fit the site, by rarity x biome x region; r uniform in [0, 1). -1: none fits.
int OreVeins_PickFamily(const std::vector<OreFamily>& fams, const VeinSite& s, float r);
// Is p inside the deposit's shape? n1: lumpiness noise in about [-1, 1] at p (the old blob's); n2: the stockwork's vein
// field at p / size (only read for a stockwork; |n2| small is a vein).
bool OreVeins_Inside(const VeinDeposit& d, Vector3 p, float n1, float n2);
// The ore kind at p inside the deposit: its family's member by weight, per 2 m block (clumps, not speckle). -1: none.
int OreVeins_Member(const VeinDeposit& d, Vector3 p);
uint32_t OreVeins_Hash(int x, int y, int z, uint32_t seed);
float OreVeins_Unit(uint32_t& h); // the old deposit's LCG step: uniform [0, 1)

#endif
