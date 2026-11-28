#include "strata.h"
#include "detmath.h"
#include <algorithm>
#include <cmath>

namespace {
    uint32_t Hash(int32_t x, int32_t y, int32_t z, uint32_t seed) {
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
        h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
        return h;
    }
    float Rand01(int32_t x, int32_t y, int32_t z, uint32_t seed) { return (Hash(x, y, z, seed) & 0xFFFFFFu) / 16777216.0f; }
    // Smooth 3D value noise in [0, 1] (cell size 1).
    float ValueNoise(Vector3 p, uint32_t seed) {
        const float fx = floorf(p.x), fy = floorf(p.y), fz = floorf(p.z);
        const int32_t ix = (int32_t)fx, iy = (int32_t)fy, iz = (int32_t)fz;
        float tx = p.x - fx, ty = p.y - fy, tz = p.z - fz;
        tx = tx * tx * (3.0f - 2.0f * tx); ty = ty * ty * (3.0f - 2.0f * ty); tz = tz * tz * (3.0f - 2.0f * tz);
        auto r = [&](int dx, int dy, int dz) { return Rand01(ix + dx, iy + dy, iz + dz, seed); };
        const float x00 = r(0, 0, 0) + (r(1, 0, 0) - r(0, 0, 0)) * tx, x10 = r(0, 1, 0) + (r(1, 1, 0) - r(0, 1, 0)) * tx;
        const float x01 = r(0, 0, 1) + (r(1, 0, 1) - r(0, 0, 1)) * tx, x11 = r(0, 1, 1) + (r(1, 1, 1) - r(0, 1, 1)) * tx;
        const float y0 = x00 + (x10 - x00) * ty, y1 = x01 + (x11 - x01) * ty;
        return y0 + (y1 - y0) * tz;
    }
    Vector3 Scaled(Vector3 s, float inv) { return { s.x * inv, s.y * inv, s.z * inv }; }
}

namespace {
    float FaultField(Vector3 s, unsigned int seed, float& region);
    float FaultStep(Vector3 s, float f, float region, unsigned int seed);
}

StrataColumn Strata_Column(Vector3 s, float continent, unsigned int seed) {
    StrataColumn c;
    // Provinces: basalt on the sea floor and in rare volcanic regions (~9 km), granite cores in the inland ranges (~15 km
    // regions where the continent value is high), layered sediments everywhere else.
    const float volcanic = ValueNoise(Scaled(s, 1.0f / 9000.0f), seed ^ 0xB5A17u);
    const float granite = ValueNoise(Scaled(s, 1.0f / 15000.0f), seed ^ 0x6A41Eu);
    if (continent < -0.02f || volcanic > 0.80f) c.province = 2;
    else if (continent > 0.30f && granite > 0.42f) c.province = 1;
    // Broad layers (40-150 m thick by region), a share of them resistant; warped over ~2.5 km and tilted over ~8 km by
    // up to ~150 m, so cliff bands wander and dip rather than following the contours.
    c.thick = 40.0f + 110.0f * ValueNoise(Scaled(s, 1.0f / 20000.0f), seed ^ 0x7A1C3u);
    c.hardShare = (uint8_t)(255.0f * (0.35f + 0.25f * ValueNoise(Scaled(s, 1.0f / 12000.0f), seed ^ 0x5AE2u)));
    c.offset = 300.0f * (ValueNoise(Scaled(s, 1.0f / 8000.0f), seed ^ 0x711Fu) - 0.5f)
             + 70.0f * (ValueNoise(Scaled(s, 1.0f / 2500.0f), seed ^ 0x3A9Bu) - 0.5f);
    // Cuestas and hogbacks (#591): in rare fold regions (~20 km) the layers dip 5-40% along a regional axis, so each
    // resistant layer crops out as a ridge along the strike (thick / dip apart: ~200-1600 m) that the bake carves
    // out of the soft beds: gentle dips give cuestas (a scarp and a long back slope), steep ones hogbacks. Basalt
    // flows and granite cores don't fold here.
    if (c.province == 0) {
        const float w = Strata_FoldWeight(s, seed);
        if (w > 0.0f) {
            const float dip = 0.05f + 0.35f * ValueNoise(Scaled(s, 1.0f / 30000.0f), seed ^ 0xF01D1u);
            const float a = 6.2831853f * ValueNoise(Scaled(s, 1.0f / 40000.0f), seed ^ 0xF01D2u);
            // Folds (anticlines and synclines ~16 km apart) across an axis in a fixed frame (on a ~20 km region the
            // sphere is near enough flat): limbs dip up to `dip`. A periodic offset keeps the fade-in at the region's
            // edge gentle (a planar tilt would sweep through layers metres apart there).
            const Vector3 axis = { det::Cosf(a), 0.5f * det::Sinf(a), det::Sinf(a) };
            const float lambda = 2500.0f;
            c.offset += w * dip * lambda * det::Sinf((s.x * axis.x + s.y * axis.y + s.z * axis.z) / lambda);
        }
    }
    // Fault scarps (#590): hard rock on the raised side near the line, so the escarpment stands while the bake wears it.
    float region = 0.0f;
    const float f = FaultField(s, seed, region);
    if (region > 0.0f && f > -0.006f && f < 0.05f && c.province == 0)
        c.hardShare = (uint8_t)std::max<float>(c.hardShare, 255.0f * (0.6f + 0.25f * region));
    if (gStrataBeds) {
        // #2064: the raised block carries its beds up with it, so they step across the fault by its throw (before, the
        // land stepped and the layers ran on level through the scarp).
        if (region > 0.0f) c.offset -= FaultStep(s, f, region, seed);
        c.pinch = ValueNoise(Scaled(s, 1.0f / 3500.0f), seed ^ 0xB3D5u);
    }
    return c;
}

StrataBed Strata_BedAt(const StrataColumn& c, float h, unsigned int seed) {
    // Basalt flows with softer ash beds between them (flows ~0.6 of a layer).
    const bool basalt = c.province == 2;
    const float thick = (basalt ? 0.6f : 1.0f) * c.thick;
    const float share = basalt ? std::max(c.hardShare / 255.0f, 0.65f) : c.hardShare / 255.0f;
    const float t = (h + c.offset) / thick;
    if (!gStrataBeds) { // #587's even layers
        const int32_t i = (int32_t)floorf(t);
        return { i, t - floorf(t), thick, Rand01(i, 0, basalt ? 11 : 0, seed) < share };
    }
    // #2064: sets of four nominal layers. In each the beds' shares are skewed (0.2 + 2.4 r^3: most beds thin, now and then
    // one three or four times the mean), every bed a quarter likely lenticular (it thins to nothing where the column's
    // lateral value is under its threshold: a pinch-out, which shifts the beds above it down), resistant beds thickened by
    // the caprock weight; the set keeps its total, so the beds can be walked in one pass. Equal layers read as built
    // terraces (#543). One hash per bed: this runs per bake cell and step.
    constexpr int kSet = 4;
    const int32_t j = (int32_t)floorf(t / kSet);
    const float u = t - (float)j * kSet;
    float w[kSet], sum = 0.0f;
    bool hard[kSet];
    for (int k = 0; k < kSet; ++k) {
        const uint32_t hb = Hash(j, k, basalt ? 0x2064B : 0x2064, seed);
        const float r = (hb & 0x3FFu) / 1024.0f, hr = ((hb >> 10) & 0x3FFu) / 1024.0f;
        const float lens = ((hb >> 20) & 0x3Fu) / 64.0f, th = ((hb >> 26) & 0x3Fu) / 64.0f;
        hard[k] = hr < share;
        float wk = 0.2f + 2.4f * r * r * r;
        if (lens < 0.25f) wk *= std::clamp((c.pinch - th) / 0.15f + 0.5f, 0.0f, 1.0f);
        if (hard[k]) wk *= 1.0f + 3.0f * c.caprock; // (GUESS)
        w[k] = wk; sum += wk;
    }
    const float scale = kSet / std::max(sum, 1e-3f);
    float acc = 0.0f;
    for (int k = 0; k < kSet; ++k) {
        const float wk = w[k] * scale;
        if (u < acc + wk || k == kSet - 1)
            return { j * kSet + k, std::clamp((u - acc) / std::max(wk, 1e-4f), 0.0f, 1.0f), wk * thick, hard[k] };
        acc += wk;
    }
    return { j * kSet, 0.0f, thick, false }; // (unreached)
}

RockType Strata_RockAt(const StrataColumn& c, float h, unsigned int seed) {
    if (c.province == 1) return RockType::Granite;
    if (c.province == 3) return RockType::Limestone; // karst: massive limestone (#592)
    const StrataBed b = Strata_BedAt(c, h, seed);
    if (c.province == 2) return b.hard ? RockType::Basalt : RockType::Soft;
    if (!b.hard) return RockType::Soft;
    const float kind = Rand01(b.index, 0, 1, seed);
    // Marble: a rare thin band at the foot of a resistant layer.
    if (kind < 0.08f) return b.frac < 0.35f ? RockType::Marble : RockType::Limestone;
    return kind < 0.55f ? RockType::Sandstone : RockType::Limestone;
}

namespace {
    // Fault field (#590): the region weight (0 outside the rare fault regions) and the signed value whose zero crossing
    // is the fault line (> 0 on the raised side). ~25 km cells give lines tens of km long; a ~6 km warp curves them.
    float FaultField(Vector3 s, unsigned int seed, float& region) {
        const float r = ValueNoise(Scaled(s, 1.0f / 40000.0f), seed ^ 0xFA170u);
        region = r < 0.6f ? 0.0f : std::min(1.0f, (r - 0.6f) / 0.08f);
        if (region <= 0.0f) return 0.0f;
        return ValueNoise(Scaled(s, 1.0f / 25000.0f), seed ^ 0xFA171u) - 0.5f
             + 0.04f * (ValueNoise(Scaled(s, 1.0f / 6000.0f), seed ^ 0xFA172u) - 0.5f);
    }
}

float Strata_FaultValue(Vector3 s, unsigned int seed, float* region) {
    float r = 0.0f;
    const float f = FaultField(s, seed, r);
    if (region) *region = r;
    return f;
}

float Strata_FaultUplift(Vector3 s, unsigned int seed) {
    float region = 0.0f;
    const float f = FaultField(s, seed, region);
    if (region <= 0.0f) return 0.0f;
    return FaultStep(s, f, region, seed);
}

namespace {
float FaultStep(Vector3 s, float f, float region, unsigned int seed) {
    // A step ~150 m wide (the field changes ~4e-5 per metre): three cells of the bake's 50 m grid, which wears it; steep
    // enough (~0.8-2) that the cliff profile (#588) gives the hard raised side a face.
    float t = std::clamp((f + 0.003f) / 0.006f, 0.0f, 1.0f);
    t = t * t * (3.0f - 2.0f * t);
    const float offset = 50.0f + 250.0f * ValueNoise(Scaled(s, 1.0f / 30000.0f), seed ^ 0xFA173u);
    return region * offset * t;
}
}

float Strata_Hardness(RockType r) {
    switch (r) {
        case RockType::Granite: return 1.0f;
        case RockType::Basalt: return 0.9f;
        case RockType::Marble: return 0.85f;
        case RockType::Sandstone: return 0.65f;
        case RockType::Limestone: return 0.55f;
        default: return 0.1f;
    }
}

float Strata_Erodibility(RockType r) { return 1.6f - 1.45f * Strata_Hardness(r); }

const char* Strata_Item(RockType r) {
    switch (r) {
        case RockType::Granite: return "granite";
        case RockType::Basalt: return "basalt";
        case RockType::Sandstone: return "sandstone";
        case RockType::Limestone: return "limestone";
        case RockType::Marble: return "marble";
        default: return "stone";
    }
}

float Strata_FoldWeight(Vector3 s, unsigned int seed) {
    if (!gStrataFolds) return 0.0f;
    const float fold = ValueNoise(Scaled(s, 1.0f / 20000.0f), seed ^ 0xF01D0u);
    return fold > 0.62f ? std::min(1.0f, (fold - 0.62f) / 0.06f) : 0.0f;
}

float Strata_OutcropBump(const StrataColumn& c, float h, unsigned int seed) {
    if (c.province != 0) return 0.0f;
    const StrataBed b = Strata_BedAt(c, h, seed);
    return b.hard ? det::Sinf(3.14159265f * b.frac) : 0.0f;
}

float Strata_KarstWeight(Vector3 s, unsigned int seed) {
    const float k = ValueNoise(Scaled(s, 1.0f / 10000.0f), seed ^ 0x4A257u);
    return k > 0.66f ? std::min(1.0f, (k - 0.66f) / 0.05f) : 0.0f;
}

float Strata_KarstRelief(Vector3 s, unsigned int seed) {
    // Jittered 3D Voronoi on a 300 m grid (the surface slices it into a 2D pattern): the nearest feature's cell decides.
    const float cell = 300.0f;
    const float gx = s.x / cell, gy = s.y / cell, gz = s.z / cell;
    const int32_t ix = (int32_t)floorf(gx), iy = (int32_t)floorf(gy), iz = (int32_t)floorf(gz);
    float best = 1e9f; int32_t bx = 0, by = 0, bz = 0;
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
                const int32_t cx = ix + dx, cy = iy + dy, cz = iz + dz;
                const float fx = cx + 0.15f + 0.7f * Rand01(cx, cy, cz, seed ^ 0x4A258u) - gx;
                const float fy = cy + 0.15f + 0.7f * Rand01(cx, cy, cz, seed ^ 0x4A259u) - gy;
                const float fz = cz + 0.15f + 0.7f * Rand01(cx, cy, cz, seed ^ 0x4A25Au) - gz;
                const float d2 = fx * fx + fy * fy + fz * fz;
                if (d2 < best) { best = d2; bx = cx; by = cy; bz = cz; }
            }
    const float dist = sqrtf(best) * cell;
    const float kind = Rand01(bx, by, bz, seed ^ 0x4A25Bu);
    const float size = Rand01(bx, by, bz, seed ^ 0x4A25Cu);
    auto smooth = [](float e0, float e1, float x) { float t = std::clamp((x - e0) / (e1 - e0), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); };
    if (kind < 0.55f) { // a tower: 60-200 m tall, 50-140 m in radius, steep sides
        const float r0 = 50.0f + 90.0f * size;
        return (60.0f + 140.0f * Rand01(bx, by, bz, seed ^ 0x4A25Du)) * (1.0f - smooth(0.6f * r0, r0, dist));
    }
    if (kind > 0.85f) { // a doline: a 15-40 m bowl, 40-90 m in radius
        const float r0 = 40.0f + 50.0f * size;
        return -(15.0f + 25.0f * size) * (1.0f - smooth(0.0f, r0, dist));
    }
    return 0.0f;
}

void Strata_KarstDolines(Vector3 lo, Vector3 hi, float R, unsigned int seed, std::vector<KarstDoline>& out) {
    // The doline features of Strata_KarstRelief (the same grid, jitter and kinds) whose surface centre (the feature
    // projected onto the sphere of radius R) lies in [lo, hi] and where the bowl really dips (a nearer tower can own the
    // spot instead): for the karst caves' sinks (#1052).
    const float cell = 300.0f;
    const int32_t x0 = (int32_t)floorf((lo.x - cell) / cell), x1 = (int32_t)floorf((hi.x + cell) / cell);
    const int32_t y0 = (int32_t)floorf((lo.y - cell) / cell), y1 = (int32_t)floorf((hi.y + cell) / cell);
    const int32_t z0 = (int32_t)floorf((lo.z - cell) / cell), z1 = (int32_t)floorf((hi.z + cell) / cell);
    for (int32_t cz = z0; cz <= z1; ++cz)
        for (int32_t cy = y0; cy <= y1; ++cy)
            for (int32_t cx = x0; cx <= x1; ++cx) {
                if (Rand01(cx, cy, cz, seed ^ 0x4A25Bu) <= 0.85f) continue;   // not a doline
                const Vector3 F = { (cx + 0.15f + 0.7f * Rand01(cx, cy, cz, seed ^ 0x4A258u)) * cell,
                                    (cy + 0.15f + 0.7f * Rand01(cx, cy, cz, seed ^ 0x4A259u)) * cell,
                                    (cz + 0.15f + 0.7f * Rand01(cx, cy, cz, seed ^ 0x4A25Au)) * cell };
                const float fl = sqrtf(F.x * F.x + F.y * F.y + F.z * F.z);
                if (fl < 1.0f) continue;
                const Vector3 dir = { F.x / fl, F.y / fl, F.z / fl };
                const Vector3 onR = { dir.x * R, dir.y * R, dir.z * R };
                if (onR.x < lo.x || onR.y < lo.y || onR.z < lo.z || onR.x >= hi.x || onR.y >= hi.y || onR.z >= hi.z) continue;
                const float depth = -Strata_KarstRelief(onR, seed);
                if (depth < 4.0f) continue;                                    // off the sphere's slice, or a tower's ground
                const float size = Rand01(cx, cy, cz, seed ^ 0x4A25Cu);
                const float r0 = 40.0f + 50.0f * size, off = fl - R;
                out.push_back({ dir, sqrtf(std::max(0.0f, r0 * r0 - off * off)), depth, (uint32_t)Hash(cx, cy, cz, seed ^ 0x4A25Eu) });
            }
}

float Strata_Volcano(Vector3 s, float planetRadius, unsigned int seed, float* weight) {
    // A jittered 20 km grid; 2% of cells hold a volcano at their feature point (projected onto the sphere): ~10 on a
    // 47 km planet.
    const float cell = 20000.0f;
    const float gx = s.x / cell, gy = s.y / cell, gz = s.z / cell;
    const int32_t ix = (int32_t)floorf(gx), iy = (int32_t)floorf(gy), iz = (int32_t)floorf(gz);
    const float len = sqrtf(s.x * s.x + s.y * s.y + s.z * s.z);
    float best = 0.0f, bestW = 0.0f;
    for (int dz = -1; dz <= 1; ++dz)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dx = -1; dx <= 1; ++dx) {
                const int32_t cx = ix + dx, cy = iy + dy, cz = iz + dz;
                if (Rand01(cx, cy, cz, seed ^ 0x7015Au) >= 0.02f) continue;
                Vector3 f = { (cx + Rand01(cx, cy, cz, seed ^ 0x7015Bu)) * cell, (cy + Rand01(cx, cy, cz, seed ^ 0x7015Cu)) * cell,
                              (cz + Rand01(cx, cy, cz, seed ^ 0x7015Du)) * cell };
                const float fl = sqrtf(f.x * f.x + f.y * f.y + f.z * f.z);
                if (fl < 1.0f) continue;
                f = Scaled(f, planetRadius / fl);
                const Vector3 sp = Scaled(s, planetRadius / std::max(len, 1.0f));
                const float dist = sqrtf((sp.x - f.x) * (sp.x - f.x) + (sp.y - f.y) * (sp.y - f.y) + (sp.z - f.z) * (sp.z - f.z));
                const float H = 1000.0f + 2000.0f * Rand01(cx, cy, cz, seed ^ 0x7015Eu);
                const float R = 4.0f * H;
                if (dist >= R) continue;
                const float q = dist / R;
                // Concave cone; the summit opens into a crater (small) or a caldera (a quarter of them, wide).
                const float qc = Rand01(cx, cy, cz, seed ^ 0x7015Fu) < 0.25f ? 0.16f : 0.07f;
                const float cone = H * (1.0f - q) * (1.0f - q);
                float h = cone;
                if (q < qc) {
                    const float rim = H * (1.0f - qc) * (1.0f - qc);
                    const float depth = (qc > 0.1f ? 0.25f : 0.12f) * H * (1.0f - (q / qc) * (q / qc));
                    h = rim - depth;
                }
                if (h > best) { best = h; bestW = 1.0f - q; }
            }
    if (weight) *weight = bestW;
    return best;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include "unit_helpers.h"

// #2064: the beds tile the column without gaps, their thicknesses are skewed (not #543's even treads), lenticular beds
// pinch out with the lateral value, the caprock weight thickens the resistant beds, and the even layers come back with
// RAYLIBGAME_STRATA_BEDS=0.
TEST(Strata, BedsSkewedLenticularAndTiled) {
    StrataColumn c; c.thick = 80.0f; c.hardShare = 115;
    struct Stats { int beds = 0; double sum = 0, sum2 = 0, maxT = 0, hard = 0, soft = 0; int nHard = 0, nSoft = 0; std::vector<float> t; };
    auto walk = [&](const StrataColumn& col) {
        Stats s;
        int32_t last = INT32_MIN; float lastFrac = 1.0f;
        for (float h = -4000.0f; h < 4000.0f; h += 0.25f) {
            const StrataBed b = Strata_BedAt(col, h, 42u);
            EXPECT_GE(b.frac, 0.0f); EXPECT_LE(b.frac, 1.0f);
            if (b.index == last) { EXPECT_GE(b.frac, lastFrac - 1e-4f) << "a bed's height runs backwards"; lastFrac = b.frac; continue; }
            if (last != INT32_MIN) { EXPECT_GT(b.index, last); EXPECT_LT(b.frac * b.thick, 0.3f) << "a gap under bed " << b.index; }
            last = b.index; lastFrac = b.frac;
            ++s.beds; s.sum += b.thick; s.sum2 += (double)b.thick * b.thick; s.maxT = std::max<double>(s.maxT, b.thick); s.t.push_back(b.thick);
            if (b.hard) { s.hard += b.thick; ++s.nHard; } else { s.soft += b.thick; ++s.nSoft; }
        }
        return s;
    };
    const Stats a = walk(c);
    ASSERT_GT(a.beds, 50);
    const double mean = a.sum / a.beds, cv = std::sqrt(std::max(0.0, a.sum2 / a.beds - mean * mean)) / mean;
    EXPECT_GT(cv, 0.5) << "the beds are too even (#543's terraces)";
    EXPECT_GT(a.maxT / mean, 2.5) << "no thick bed";

    StrataColumn p0 = c, p1 = c; p0.pinch = 0.0f; p1.pinch = 1.0f;
    int pinched = 0;
    for (float h = -4000.0f; h < 4000.0f; h += 40.0f) pinched += Strata_BedAt(p0, h, 42u).index != Strata_BedAt(p1, h, 42u).index;
    EXPECT_GT(pinched, 10) << "no bed pinches out across the lateral value";

    StrataColumn cap = c; cap.caprock = 1.0f;
    const Stats k = walk(cap);
    EXPECT_GT((k.hard / k.nHard) / (k.soft / k.nSoft), 1.8 * (a.hard / a.nHard) / (a.soft / a.nSoft)) << "the caprock weight didn't thicken the resistant beds";

    gStrataBeds = false;
    const Stats e = walk(c);
    gStrataBeds = true;
    EXPECT_NEAR(e.maxT, 80.0, 1e-3); EXPECT_NEAR(e.sum / e.beds, 80.0, 1e-3);
}

// #2064: across a fault the beds step by its throw: the raised side's offset is lowered by the uplift.
TEST(Strata, BedsStepAcrossAFault) {
    int found = 0;
    for (int i = 0; i < 20000 && found < 20; ++i) {
        const float y = 1.0f - 2.0f * (i + 0.5f) / 20000.0f, r = std::sqrt(std::max(0.0f, 1.0f - y * y)), a = 2.39996323f * i;
        const Vector3 s = { 60000.0f * r * std::cos(a), 60000.0f * y, 60000.0f * r * std::sin(a) };
        const float up = Strata_FaultUplift(s, 42u);
        if (up < 20.0f) continue;
        gStrataBeds = false;
        const StrataColumn even = Strata_Column(s, 0.1f, 42u);
        gStrataBeds = true;
        const StrataColumn beds = Strata_Column(s, 0.1f, 42u);
        EXPECT_NEAR(even.offset - beds.offset, up, 1e-3f);
        ++found;
    }
    EXPECT_GT(found, 0) << "no fault on the test sphere";
}
#endif
