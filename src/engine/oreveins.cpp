#include "oreveins.h"
#include <raymath.h>

#include <algorithm>
#include <cmath>

uint32_t OreVeins_Hash(int x, int y, int z, uint32_t seed) {
    uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
    h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
    return h;
}
float OreVeins_Unit(uint32_t& h) { h = h * 1664525u + 1013904223u; return (h >> 8) * (1.0f / 16777216.0f); }

int OreVeins_PickFamily(const std::vector<OreFamily>& fams, const VeinSite& s, float r) {
    float total = 0.0f;
    std::vector<float> w(fams.size(), 0.0f);
    for (size_t i = 0; i < fams.size(); ++i) {
        const OreFamily& f = fams[i];
        if (s.depth < f.minDepth || s.depth > f.maxDepth || !OreFamilies_Hosts(f, s.rock)) continue;
        w[i] = std::max(0.0f, OreFamilies_Rarity(f, s.biome, s.region));
        total += w[i];
    }
    if (total <= 0.0f) return -1;
    float pick = std::clamp(r, 0.0f, 1.0f) * total, acc = 0.0f;
    int last = -1;
    for (size_t i = 0; i < fams.size(); ++i) {
        if (w[i] <= 0.0f) continue;
        acc += w[i]; last = (int)i;
        if (pick < acc) return (int)i;
    }
    return last;
}

VeinDeposit OreVeins_Make(int cx, int cy, int cz, uint32_t h, const std::function<VeinSite(Vector3 centre)>& site) {
    VeinDeposit d;
    if (OreVeins_Unit(h) > 0.35f) return d;
    // The centre as the old deposits had it: 24-72 m into the cell, so every shape (kVeinReach) stays inside it
    d.centre = { (cx + 0.25f + 0.5f * OreVeins_Unit(h)) * kVeinCell, (cy + 0.25f + 0.5f * OreVeins_Unit(h)) * kVeinCell, (cz + 0.25f + 0.5f * OreVeins_Unit(h)) * kVeinCell };
    const float sizeRoll = OreVeins_Unit(h), familyRoll = OreVeins_Unit(h);
    const auto& fams = OreFamilies_All();
    const VeinSite s = site ? site(d.centre) : VeinSite{};
    d.family = OreVeins_PickFamily(fams, s, familyRoll);
    if (d.family < 0) return d; // no family fits this rock and depth: no deposit
    const OreFamily& f = fams[d.family];
    d.exists = true;
    d.shape = f.shape;
    d.size = f.minSize + (f.maxSize - f.minSize) * sizeRoll;
    const Vector3 up = Vector3Normalize(d.centre);
    if (d.shape == VeinShape::Sheet) { // a dike or vein plane: steep (within ~35 degrees of vertical) at a random strike
        const float a = 6.2831853f * OreVeins_Unit(h), tilt = 0.6f * (OreVeins_Unit(h) - 0.5f);
        const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 e2 = Vector3CrossProduct(up, e1);
        const Vector3 horiz = Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a))); // the normal, before the tilt
        d.axis = Vector3Normalize(Vector3Add(horiz, Vector3Scale(up, tilt)));
    } else {
        d.axis = up; // a pipe runs along it, a seam lies across it
    }
    d.envelope = d.shape == VeinShape::Blob ? d.size : kVeinReach;
    d.seed = h;
    return d;
}

bool OreVeins_Inside(const VeinDeposit& d, Vector3 p, float n1, float n2) {
    if (!d.exists) return false;
    const Vector3 q = Vector3Subtract(p, d.centre);
    switch (d.shape) {
        case VeinShape::Blob: // today's: lumpy, 0.35 of the noise on the radius
            return Vector3Length(q) / d.size + 0.35f * n1 < 1.0f;
        case VeinShape::Sheet: { // a slab size thick, wavy faces, within kVeinReach of the centre
            const float off = Vector3DotProduct(q, d.axis);
            const float lateral = Vector3Length(Vector3Subtract(q, Vector3Scale(d.axis, off)));
            return fabsf(off) < 0.5f * d.size * (1.0f + 0.5f * n1) && lateral < kVeinReach * (1.0f + 0.15f * n1) - 2.0f;
        }
        case VeinShape::Pipe: { // vertical, size in radius, tapering a little toward its ends
            const float along = Vector3DotProduct(q, d.axis);
            const float radial = Vector3Length(Vector3Subtract(q, Vector3Scale(d.axis, along)));
            const float taper = 1.0f - 0.3f * fabsf(along) / kVeinReach;
            return fabsf(along) < kVeinReach && radial < d.size * taper * (1.0f + 0.35f * n1);
        }
        case VeinShape::Seam: { // a level bed size thick, undulating 2 m, out to kVeinReach
            const float along = Vector3DotProduct(q, d.axis) + 2.0f * n1;
            const float lateral = Vector3Length(Vector3Subtract(q, Vector3Scale(d.axis, Vector3DotProduct(q, d.axis))));
            return fabsf(along) < 0.5f * d.size && lateral < kVeinReach * (1.0f + 0.1f * n1) - 1.0f;
        }
        case VeinShape::Stockwork: // thin veins (where the vein field crosses zero) through a 20 m envelope
            return Vector3Length(q) < (kVeinReach - 2.0f) * (1.0f + 0.1f * n1) && fabsf(n2) < 0.09f;
    }
    return false;
}

int OreVeins_Member(const VeinDeposit& d, Vector3 p) {
    const auto& fams = OreFamilies_All();
    if (d.family < 0 || d.family >= (int)fams.size()) return -1;
    uint32_t h = OreVeins_Hash((int)floorf(p.x / 2.0f), (int)floorf(p.y / 2.0f), (int)floorf(p.z / 2.0f), d.seed);
    const OreFamilyMember* m = OreFamilies_Member(fams[d.family], OreVeins_Unit(h));
    return m ? m->kind : -1;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

namespace {
VeinDeposit Shaped(VeinShape s, float size) {
    VeinDeposit d;
    d.exists = true; d.centre = { 0, 1000, 0 }; d.shape = s; d.size = size; d.axis = { 0, 1, 0 };
    d.envelope = s == VeinShape::Blob ? size : kVeinReach;
    return d;
}
bool In(const VeinDeposit& d, float x, float y, float z, float n2 = 0.5f) { return OreVeins_Inside(d, { d.centre.x + x, d.centre.y + y, d.centre.z + z }, 0.0f, n2); }
}

TEST(OreVeins, EachShapeHoldsItsPointsAndNoneReachesPastItsCell) {
    const VeinDeposit blob = Shaped(VeinShape::Blob, 6.0f);
    EXPECT_TRUE(In(blob, 5, 0, 0)); EXPECT_FALSE(In(blob, 7, 0, 0));
    VeinDeposit sheet = Shaped(VeinShape::Sheet, 4.0f); sheet.axis = { 1, 0, 0 }; // a vertical plane facing x
    EXPECT_TRUE(In(sheet, 1.5f, 10, 10)) << "within its thickness, far along the plane";
    EXPECT_FALSE(In(sheet, 3, 0, 0)) << "past half its thickness";
    const VeinDeposit pipe = Shaped(VeinShape::Pipe, 5.0f);
    EXPECT_TRUE(In(pipe, 3, 18, 0)) << "high up its axis";
    EXPECT_FALSE(In(pipe, 8, 0, 0)) << "off the axis";
    const VeinDeposit seam = Shaped(VeinShape::Seam, 3.0f);
    EXPECT_TRUE(In(seam, 18, 1, 0)) << "far out along the bed";
    EXPECT_FALSE(In(seam, 0, 3, 0)) << "above the bed";
    const VeinDeposit stock = Shaped(VeinShape::Stockwork, 8.0f);
    EXPECT_TRUE(In(stock, 5, 5, 5, 0.02f)) << "on a vein";
    EXPECT_FALSE(In(stock, 5, 5, 5, 0.5f)) << "between veins";
    // Nothing past kVeinReach (so a deposit never crosses into the next cell), whatever the noise
    for (const VeinDeposit& d : { blob, sheet, pipe, seam, stock })
        for (float n1 : { -1.0f, 1.0f })
            for (const Vector3 dir : { Vector3{ 1, 0, 0 }, Vector3{ 0, 1, 0 }, Vector3{ 0, 0, 1 }, Vector3Normalize({ 1, 1, 1 }) })
                EXPECT_FALSE(OreVeins_Inside(d, Vector3Add(d.centre, Vector3Scale(dir, kVeinReach + 2.5f)), n1, 0.0f)) << OreFamilies_ShapeName(d.shape);
}

TEST(OreVeins, TheFamilyFitsItsRockAndDepthAndFollowsRarity) {
    const auto f = OreFamilies_Parse(
        "a | iron 1   | granite   | 10 100 | any       | any        | 1 | blob | 4\n"
        "b | copper 1 | granite   | 50 -   | desert 3  | volcanic 2 | 1 | blob | 4\n"
        "c | coal 1   | sandstone | 0 -    | any       | any        | 1 | seam | 3\n");
    VeinSite s; s.rock = "granite"; s.depth = 20;
    EXPECT_EQ(OreVeins_PickFamily(f, s, 0.99f), 0) << "only a fits granite at 20 m";
    s.rock = "basalt";
    EXPECT_EQ(OreVeins_PickFamily(f, s, 0.5f), -1) << "nothing fits basalt: no deposit";
    s.rock = "granite"; s.depth = 60; s.biome = "desert"; s.region = "volcanic";
    int n[2] = {};
    for (int i = 0; i < 700; ++i) n[OreVeins_PickFamily(f, s, (i + 0.5f) / 700.0f)]++;
    EXPECT_EQ(n[0], 100); EXPECT_EQ(n[1], 600) << "b: 1 x desert 3 x volcanic 2 = 6 to a's 1";
}

TEST(OreVeins, ADepositIsTheSameEveryTimeAndKeepsTheOldCentres) {
    auto site = [](Vector3) { VeinSite s; s.rock = "rock"; s.depth = 40; return s; };
    int made = 0;
    for (int i = 0; i < 400; ++i) {
        const uint32_t h = OreVeins_Hash(i, 3, -7, 42u);
        const VeinDeposit a = OreVeins_Make(i, 3, -7, h, site), b = OreVeins_Make(i, 3, -7, h, site);
        EXPECT_EQ(a.exists, b.exists); EXPECT_EQ(a.family, b.family); EXPECT_EQ(a.size, b.size);
        // The old deposit roll: exists under 0.35, then the centre (the same three draws)
        uint32_t o = h;
        const bool old = OreVeins_Unit(o) <= 0.35f;
        if (!old) { EXPECT_FALSE(a.exists); continue; }
        const float x = (i + 0.25f + 0.5f * OreVeins_Unit(o)) * kVeinCell;
        EXPECT_EQ(a.centre.x, x) << "the same place as before";
        if (a.exists) {
            ++made;
            EXPECT_GE(a.family, 0);
            const int k = OreVeins_Member(a, a.centre);
            EXPECT_EQ(k, OreVeins_Member(a, a.centre)) << "deterministic";
            EXPECT_GE(k, 0);
        }
    }
    EXPECT_GT(made, 50) << "about 35% of 400 cells, those a family fits";
}
#endif
