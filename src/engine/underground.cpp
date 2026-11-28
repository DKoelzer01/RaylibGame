#include "underground.h"
#include "voxel.h"
#include "depthtable.h"

#include <raymath.h>
#include <cmath>
#include <unordered_map>

namespace {
constexpr float kCell = 256.0f; // voxel.cpp kCaveCell: the cave graph's cell

uint32_t CellHash(int x, int y, int z, uint32_t seed) {
    uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
    h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
    return h;
}
} // namespace

const char* UndergroundRegion_Name(UndergroundRegion r) {
    static const char* names[] = { "none", "karst", "volcanic", "frozen", "crystal", "flooded", "fungal" };
    const int i = (int)r;
    return i >= 0 && i < (int)UndergroundRegion::Count ? names[i] : "?";
}

// In priority order (GUESS thresholds): the rock decides first (karst limestone, volcanic basalt), then the cold (frozen
// near the surface), then rare deep crystal pockets, then water (under the sea, swamps, very wet ground near the surface),
// then humid forests' fungal caves.
UndergroundRegion UndergroundRegion_Classify(const UndergroundInputs& in) {
    if (in.depth < 0.0f) return UndergroundRegion::None;                                   // above the ground
    if (in.ocean) return in.depth < kDepth.floodedSeaMax ? UndergroundRegion::Flooded : UndergroundRegion::None;
    if (in.province == 3) return UndergroundRegion::Karst;
    if (in.province == 2) return UndergroundRegion::Volcanic;
    if (in.temp < -2.0f && in.depth < kDepth.frozenMax) return UndergroundRegion::Frozen;
    if (in.depth > kDepth.crystalMin && in.roll < 0.12f) return UndergroundRegion::Crystal;
    const Biome b = (Biome)in.biome;
    if (b == Biome::Swamp || (in.moist > 0.8f && in.depth < kDepth.floodedWetMax)) return UndergroundRegion::Flooded;
    if ((b == Biome::Jungle || b == Biome::Forest) && in.moist > 0.55f && in.depth > kDepth.fungalMin && in.depth < kDepth.fungalMax && in.roll < 0.45f) return UndergroundRegion::Fungal;
    return UndergroundRegion::None;
}

UndergroundInputs Underground_InputsAt(const PlanetField& pf, Vector3 p) {
    const int cx = (int)floorf(p.x / kCell), cy = (int)floorf(p.y / kCell), cz = (int)floorf(p.z / kCell);
    const Vector3 c = { (cx + 0.5f) * kCell, (cy + 0.5f) * kCell, (cz + 0.5f) * kCell };
    UndergroundInputs in;
    const float len = Vector3Length(c);
    if (len < 1.0f) { in.depth = pf.radius; return in; }
    const Vector3 dir = Vector3Scale(c, 1.0f / len);
    const float elev = pf.elevation(dir);
    in.depth = pf.radius + elev - len;
    in.ocean = elev < 0.0f;
    in.province = pf.strataColumn(dir).province;
    const PlanetField::Climate cl = pf.climateAt(dir, elev);
    in.temp = cl.temp; in.moist = cl.moist;
    in.biome = (int)pf.biomeAt(dir, elev);
    in.roll = (CellHash(cx, cy, cz, pf.oreSeed ^ 0x5ca1ab1eu) >> 8) * (1.0f / 16777216.0f);
    return in;
}

UndergroundRegion Underground_RegionAt(const PlanetField& pf, Vector3 p) {
    if (pf.moon) return UndergroundRegion::None;
    // Per thread, per cell (the classification reads the elevation, strata and climate once per 256 m cell).
    thread_local std::unordered_map<uint64_t, uint8_t>* cache = nullptr; // (never destroyed: MinGW thread_local, #1133)
    thread_local const void* owner = nullptr;
    thread_local unsigned ownerSeed = 0;
    if (!cache) cache = new std::unordered_map<uint64_t, uint8_t>();
    if (owner != &pf || ownerSeed != pf.oreSeed || cache->size() > 65536) { cache->clear(); owner = &pf; ownerSeed = pf.oreSeed; }
    const int cx = (int)floorf(p.x / kCell), cy = (int)floorf(p.y / kCell), cz = (int)floorf(p.z / kCell);
    const uint64_t key = ((uint64_t)(uint32_t)(cx & 0x1FFFFF) << 42) | ((uint64_t)(uint32_t)(cy & 0x1FFFFF) << 21) | (uint64_t)(uint32_t)(cz & 0x1FFFFF);
    if (auto it = cache->find(key); it != cache->end()) return (UndergroundRegion)it->second;
    const UndergroundRegion r = UndergroundRegion_Classify(Underground_InputsAt(pf, p));
    (*cache)[key] = (uint8_t)r;
    return r;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"
TEST(Underground, ClassifyFollowsTheRules) {
    UndergroundInputs in;
    in.depth = 80.0f; in.temp = 15.0f; in.moist = 0.5f; in.biome = (int)Biome::Meadow; in.roll = 0.9f;
    EXPECT_EQ(UndergroundRegion_Classify(in), UndergroundRegion::None) << "a temperate meadow: plain caves";
    { auto k = in; k.depth = -5.0f; k.province = 3; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::None) << "above the ground"; }
    { auto k = in; k.province = 3; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::Karst); }
    { auto k = in; k.province = 2; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::Volcanic); }
    { auto k = in; k.temp = -10.0f; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::Frozen); }
    { auto k = in; k.temp = -10.0f; k.depth = 300.0f; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::None) << "deep: the cold doesn't reach"; }
    { auto k = in; k.depth = 300.0f; k.roll = 0.05f; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::Crystal); }
    { auto k = in; k.ocean = true; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::Flooded); }
    { auto k = in; k.biome = (int)Biome::Swamp; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::Flooded); }
    { auto k = in; k.biome = (int)Biome::Jungle; k.moist = 0.7f; k.roll = 0.2f; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::Fungal); }
    { auto k = in; k.biome = (int)Biome::Jungle; k.moist = 0.7f; k.roll = 0.6f; EXPECT_EQ(UndergroundRegion_Classify(k), UndergroundRegion::None) << "fungal in patches"; }
    EXPECT_STREQ(UndergroundRegion_Name(UndergroundRegion::Fungal), "fungal");
}

// #1952: the depth table's bands, and the orders its readers rely on.
TEST(DepthTable, BandsAndTheirOrder) {
    const DepthTable& t = kDepth;
    EXPECT_EQ(t.bandAt(10.0f), DepthBand::Surface);
    EXPECT_EQ(t.bandAt(t.shallow), DepthBand::Shallow);
    EXPECT_EQ(t.bandAt(t.middle - 0.01f), DepthBand::Shallow);
    EXPECT_EQ(t.bandAt(t.middle), DepthBand::Middle);
    EXPECT_EQ(t.bandAt(t.deep), DepthBand::Deep);
    EXPECT_EQ(t.bandAt(5000.0f), DepthBand::Deep);
    EXPECT_LT(t.shallow, t.middle); EXPECT_LT(t.middle, t.deep);
    EXPECT_LE(t.passageMin, t.chamberMin); EXPECT_LE(t.chamberMin, t.passageMax) << "chambers hang off passage vertices";
    EXPECT_LE(t.cavernMin, t.passageMax) << "caverns hang off passage vertices";
    EXPECT_LT(t.sheetMin, t.sheetMax); EXPECT_LE(t.sheetMax, t.passageMax);
    EXPECT_LT(t.fungalMin, t.fungalMax);
    EXPECT_LT(t.mantleClearance, t.noiseCavernMantle);
}
#endif
