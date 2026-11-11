#include "terraindeltas.h"
#include "voxel.h"
#include "tools.h"

#include <algorithm>
#include <cmath>

namespace {
    constexpr int kBias = 1 << 20; // chunk coordinates are packed as 21-bit unsigned values
    inline int FloorDiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
    inline int Idx(int lx, int ly, int lz) { return (lz * kDeltaChunk + ly) * kDeltaChunk + lx; }
}

TerrainDeltas::Key TerrainDeltas::KeyOf(int cx, int cy, int cz) {
    return ((Key)(uint32_t)(cx + kBias) << 42) | ((Key)(uint32_t)(cy + kBias) << 21) | (Key)(uint32_t)(cz + kBias);
}

void TerrainDeltas::Unpack(Key k, int& cx, int& cy, int& cz) {
    cx = (int)((k >> 42) & 0x1FFFFF) - kBias;
    cy = (int)((k >> 21) & 0x1FFFFF) - kBias;
    cz = (int)(k & 0x1FFFFF) - kBias;
}

TerrainDeltas::Key TerrainDeltas::RegionOf(Key chunk) {
    int x, y, z;
    Unpack(chunk, x, y, z);
    return KeyOf(FloorDiv(x, kRegionChunks), FloorDiv(y, kRegionChunks), FloorDiv(z, kRegionChunks));
}

float TerrainDeltas::at(int x, int y, int z) const {
    if (chunks.empty()) return 0.0f;
    const int cx = FloorDiv(x, kDeltaChunk), cy = FloorDiv(y, kDeltaChunk), cz = FloorDiv(z, kDeltaChunk);
    auto it = chunks.find(KeyOf(cx, cy, cz));
    if (it == chunks.end()) return 0.0f;
    return it->second->d[Idx(x - cx * kDeltaChunk, y - cy * kDeltaChunk, z - cz * kDeltaChunk)] * kDeltaStep;
}

float TerrainDeltas::sample(Vector3 p) const {
    if (chunks.empty()) return 0.0f;
    const int x0 = (int)floorf(p.x), y0 = (int)floorf(p.y), z0 = (int)floorf(p.z);
    const float fx = p.x - x0, fy = p.y - y0, fz = p.z - z0;
    const float c00 = at(x0, y0, z0) * (1 - fx) + at(x0 + 1, y0, z0) * fx;
    const float c10 = at(x0, y0 + 1, z0) * (1 - fx) + at(x0 + 1, y0 + 1, z0) * fx;
    const float c01 = at(x0, y0, z0 + 1) * (1 - fx) + at(x0 + 1, y0, z0 + 1) * fx;
    const float c11 = at(x0, y0 + 1, z0 + 1) * (1 - fx) + at(x0 + 1, y0 + 1, z0 + 1) * fx;
    return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz;
}

void TerrainDeltas::chunksIn(Vector3 lo, Vector3 hi, std::vector<std::pair<Key, const DeltaChunk*>>& out) const {
    out.clear();
    if (chunks.empty()) return;
    const int x0 = FloorDiv((int)floorf(lo.x), kDeltaChunk), x1 = FloorDiv((int)floorf(hi.x), kDeltaChunk);
    const int y0 = FloorDiv((int)floorf(lo.y), kDeltaChunk), y1 = FloorDiv((int)floorf(hi.y), kDeltaChunk);
    const int z0 = FloorDiv((int)floorf(lo.z), kDeltaChunk), z1 = FloorDiv((int)floorf(hi.z), kDeltaChunk);
    const long long boxChunks = (long long)(x1 - x0 + 1) * (y1 - y0 + 1) * (z1 - z0 + 1);
    if (boxChunks < (long long)chunks.size()) {
        for (int z = z0; z <= z1; ++z)
            for (int y = y0; y <= y1; ++y)
                for (int x = x0; x <= x1; ++x) {
                    auto it = chunks.find(KeyOf(x, y, z));
                    if (it != chunks.end()) out.push_back({ it->first, it->second.get() });
                }
    } else {
        for (const auto& kv : chunks) {
            int x, y, z;
            Unpack(kv.first, x, y, z);
            if (x >= x0 && x <= x1 && y >= y0 && y <= y1 && z >= z0 && z <= z1) out.push_back({ kv.first, kv.second.get() });
        }
    }
}

void TerrainDeltasBuilder::apply(const TerrainEdit& e, const DensityField& field, std::unordered_set<TerrainDeltas::Key>* touched) {
    next->strokes++;
    // A carve also lowers solid points just outside the sphere, so the meshed surface sits on it.
    const float r = e.carve ? e.radius + 2.0f : e.radius;
    const int x0 = (int)floorf(e.center.x - r), x1 = (int)ceilf(e.center.x + r);
    const int y0 = (int)floorf(e.center.y - r), y1 = (int)ceilf(e.center.y + r);
    const int z0 = (int)floorf(e.center.z - r), z1 = (int)ceilf(e.center.z + r);
    DeltaChunk* chunk = nullptr;
    TerrainDeltas::Key chunkKey = ~0ull;
    for (int z = z0; z <= z1; ++z)
        for (int y = y0; y <= y1; ++y)
            for (int x = x0; x <= x1; ++x) {
                const float dx = x - e.center.x, dy = y - e.center.y, dz = z - e.center.z;
                const float d2 = dx * dx + dy * dy + dz * dz;
                if (d2 >= r * r) continue;
                const Vector3 p = { (float)x, (float)y, (float)z };
                if (!field.editable(p)) continue; // e.g. a planet's mantle
                const float dist = sqrtf(d2);
                const float t = 1.0f - dist / r;
                float add = e.strength * t * t * (3.0f - 2.0f * t);
                float carveTo = 0.0f; // carve: the density this point may have at most
                if (e.carve) carveTo = ISO_LEVEL + DENSITY_GRADIENT * (dist - e.radius);
                const int cx = (int)floorf((float)x / kDeltaChunk), cy = (int)floorf((float)y / kDeltaChunk), cz = (int)floorf((float)z / kDeltaChunk);
                const TerrainDeltas::Key k = TerrainDeltas::KeyOf(cx, cy, cz);
                if (k != chunkKey) {
                    chunkKey = k;
                    auto own = owned.find(k);
                    if (own != owned.end()) chunk = own->second;
                    else {
                        // First touch in this build: copy (or create) the chunk so the published
                        // snapshot workers may hold is never written.
                        auto it = next->chunks.find(k);
                        auto copy = it != next->chunks.end() ? std::make_shared<DeltaChunk>(*it->second) : std::make_shared<DeltaChunk>();
                        chunk = copy.get();
                        next->chunks[k] = copy;
                        owned[k] = chunk;
                    }
                    if (touched) touched->insert(k);
                }
                int16_t& v = chunk->d[Idx(x - cx * kDeltaChunk, y - cy * kDeltaChunk, z - cz * kDeltaChunk)];
                if (e.carve) {
                    const float base = field.baseDensity(p);
                    if (e.softOnly && base >= ISO_LEVEL && !IsSoft(field.material(p, base))) continue;
                    if (base + v * kDeltaStep > carveTo) v = (int16_t)std::clamp((int)lroundf((carveTo - base) / kDeltaStep), -32767, 32767);
                } else {
                    v = (int16_t)std::clamp((int)lroundf(v + add / kDeltaStep), -32767, 32767);
                }
            }
}
