#include "strata.h"
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

namespace { float FaultField(Vector3 s, unsigned int seed, float& region); }

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
    // Fault scarps (#590): hard rock on the raised side near the line, so the escarpment stands while the bake wears it.
    float region = 0.0f;
    const float f = FaultField(s, seed, region);
    if (region > 0.0f && f > -0.006f && f < 0.05f && c.province == 0)
        c.hardShare = (uint8_t)std::max<float>(c.hardShare, 255.0f * (0.6f + 0.25f * region));
    return c;
}

RockType Strata_RockAt(const StrataColumn& c, float h, unsigned int seed) {
    if (c.province == 1) return RockType::Granite;
    const float share = c.hardShare / 255.0f;
    if (c.province == 2) {
        // Basalt flows with softer ash beds between them (flows ~0.6 of a layer).
        const float t = (h + c.offset) / (0.6f * c.thick);
        const int32_t i = (int32_t)floorf(t);
        return Rand01(i, 0, 11, seed) < std::max(share, 0.65f) ? RockType::Basalt : RockType::Soft;
    }
    const float t = (h + c.offset) / c.thick;
    const int32_t i = (int32_t)floorf(t);
    if (Rand01(i, 0, 0, seed) >= share) return RockType::Soft;
    const float kind = Rand01(i, 0, 1, seed);
    // Marble: a rare thin band at the foot of a resistant layer.
    if (kind < 0.08f) return t - floorf(t) < 0.35f ? RockType::Marble : RockType::Limestone;
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

float Strata_FaultUplift(Vector3 s, unsigned int seed) {
    float region = 0.0f;
    const float f = FaultField(s, seed, region);
    if (region <= 0.0f) return 0.0f;
    // A step ~150 m wide (the field changes ~4e-5 per metre): three cells of the bake's 50 m grid, which wears it; steep
    // enough (~0.8-2) that the cliff profile (#588) gives the hard raised side a face.
    float t = std::clamp((f + 0.003f) / 0.006f, 0.0f, 1.0f);
    t = t * t * (3.0f - 2.0f * t);
    const float offset = 50.0f + 250.0f * ValueNoise(Scaled(s, 1.0f / 30000.0f), seed ^ 0xFA173u);
    return region * offset * t;
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
