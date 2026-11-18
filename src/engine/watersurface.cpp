#include "watersurface.h"
#include "voxel.h"
#include <raymath.h>
#include <algorithm>
#include <atomic>
#include <cmath>
#include <list>
#include <mutex>
#include <unordered_map>
#include <vector>

namespace WaterSurface {

struct Tile {
    std::vector<float> draw; // (kTile + 1)^2 node levels (m over sea level, >= 0); empty: no inland water near the tile
};

namespace {

constexpr int kNodes = kTile + 1;
constexpr int kMaxLevel = 14;
constexpr long long kCapBytes = 64ll << 20;

// Cube faces: axis a = face / 2, sign + for even faces. A direction's face is its largest component; (u, v) are the next
// two components (cyclically) over that one's magnitude, in [-1, 1].
void FaceOf(Vector3 d, int& face, float& u, float& v) {
    const float c[3] = { d.x, d.y, d.z }, ac[3] = { fabsf(d.x), fabsf(d.y), fabsf(d.z) };
    const int a = ac[0] >= ac[1] && ac[0] >= ac[2] ? 0 : ac[1] >= ac[2] ? 1 : 2;
    face = a * 2 + (c[a] >= 0.0f ? 0 : 1);
    const float m = std::max(ac[a], 1e-12f);
    u = c[(a + 1) % 3] / m;
    v = c[(a + 2) % 3] / m;
}
Vector3 DirOf(int face, float u, float v) {
    const int a = face / 2;
    float c[3];
    c[a] = (face & 1) ? -1.0f : 1.0f;
    c[(a + 1) % 3] = u;
    c[(a + 2) % 3] = v;
    return Vector3Normalize({ c[0], c[1], c[2] });
}

struct Entry { std::shared_ptr<const Tile> tile; std::list<std::pair<const PlanetField*, uint64_t>>::iterator lru; };
std::mutex mx;
std::unordered_map<const PlanetField*, std::unordered_map<uint64_t, Entry>> cache;
std::list<std::pair<const PlanetField*, uint64_t>> lru; // most recent first
std::atomic<long long> bytes{ 0 }, built{ 0 };

uint64_t KeyOf(int level, int face, int ti, int tj) { return ((uint64_t)level << 59) | ((uint64_t)face << 56) | ((uint64_t)(uint32_t)ti << 28) | (uint64_t)(uint32_t)tj; }
long long BytesOf(const Tile& t) { return (long long)sizeof(Tile) + (long long)t.draw.size() * (long long)sizeof(float); }

// The tile's nodes. Raw levels over the tile and a one-node apron (a dry node's wet neighbours may lie in the next tile), then
// the drawn levels: wet nodes keep theirs; dry ones beside water take the lower of the highest neighbouring level and just
// under their own ground; others 0.
std::shared_ptr<const Tile> Build(const PlanetField& pf, int n, float cell, float fpEval, int face, int ti, int tj) {
    auto t = std::make_shared<Tile>();
    built.fetch_add(1);
    auto nodeDir = [&](int iu, int iv) {
        iu = std::clamp(iu, 0, n); iv = std::clamp(iv, 0, n);
        return DirOf(face, -1.0f + 2.0f * iu / n, -1.0f + 2.0f * iv / n);
    };
    const int u0 = ti * kTile, v0 = tj * kTile;
    // Cheap rejection: no inland water within the tile's reach (its half diagonal, a node of apron, and the dry-node margin).
    const float extent = 0.75f * kTile * cell + 2.0f * cell + 2.0f;
    if (!pf.inlandWaterNear(nodeDir(u0 + kTile / 2, v0 + kTile / 2), extent, fpEval)) return t;
    constexpr int A = kNodes + 2; // with the apron
    std::vector<float> raw((size_t)A * A);
    bool any = false;
    for (int j = 0; j < A; ++j)
        for (int i = 0; i < A; ++i) {
            const float w = pf.waterElevationMesh(nodeDir(u0 + i - 1, v0 + j - 1), fpEval);
            raw[(size_t)j * A + i] = w;
            any = any || w != PlanetField::kNoWater;
        }
    if (!any) return t;
    t->draw.resize((size_t)kNodes * kNodes);
    for (int j = 0; j < kNodes; ++j)
        for (int i = 0; i < kNodes; ++i) {
            float w = raw[(size_t)(j + 1) * A + (i + 1)];
            if (w == PlanetField::kNoWater) {
                float wet = PlanetField::kNoWater;
                for (int b = 0; b <= 2; ++b)
                    for (int a = 0; a <= 2; ++a) wet = std::max(wet, raw[(size_t)(j + b) * A + (i + a)]);
                w = wet == PlanetField::kNoWater ? 0.0f : std::min(wet, pf.elevation(nodeDir(u0 + i, v0 + j), fpEval) - std::max(0.5f, cell));
            }
            t->draw[(size_t)j * kNodes + i] = std::max(0.0f, w);
        }
    return t;
}

} // namespace

Reader::Reader(const PlanetField& f, float fp) : pf(f) {
    level = std::clamp((int)ceilf(log2f(std::max(fp, 0.5f) / 0.5f) - 1e-3f), 0, kMaxLevel);
    cell = 0.5f * (float)(1 << level);
    fpEval = cell;
    n = std::max(kTile, (int)ceilf(2.0f * pf.radius / cell));
}
Reader::~Reader() = default;

const Tile* Reader::tileFor(int face, int ti, int tj) {
    const uint64_t key = KeyOf(level, face, ti, tj);
    for (const Slot& s : slots) if (s.key == key) return s.tile.get();
    std::shared_ptr<const Tile> tile;
    {
        std::lock_guard<std::mutex> lock(mx);
        auto& m = cache[&pf];
        auto it = m.find(key);
        if (it != m.end()) { lru.splice(lru.begin(), lru, it->second.lru); tile = it->second.tile; }
    }
    if (!tile) {
        tile = Build(pf, n, cell, fpEval, face, ti, tj); // outside the lock: two threads may build one tile (same result)
        std::lock_guard<std::mutex> lock(mx);
        auto& m = cache[&pf];
        auto it = m.find(key);
        if (it != m.end()) tile = it->second.tile;
        else {
            lru.push_front({ &pf, key });
            m[key] = { tile, lru.begin() };
            bytes += BytesOf(*tile);
            while (bytes.load() > kCapBytes && lru.size() > 1) { // oldest first; readers keep theirs alive through shared_ptr
                const auto [f, k] = lru.back();
                auto& mf = cache[f];
                auto e = mf.find(k);
                if (e != mf.end()) { bytes -= BytesOf(*e->second.tile); mf.erase(e); }
                lru.pop_back();
            }
        }
    }
    Slot& s = slots[next];
    next = (next + 1) & 3;
    s.key = key; s.tile = tile;
    return s.tile.get();
}

float Reader::drawLevel(Vector3 dir) {
    int face; float u, v;
    FaceOf(dir, face, u, v);
    const float fu = (u + 1.0f) * 0.5f * n, fv = (v + 1.0f) * 0.5f * n;
    const int iu = std::clamp((int)floorf(fu), 0, n - 1), iv = std::clamp((int)floorf(fv), 0, n - 1);
    const float fx = std::clamp(fu - iu, 0.0f, 1.0f), fy = std::clamp(fv - iv, 0.0f, 1.0f);
    const int ti = iu / kTile, tj = iv / kTile;
    const Tile* t = tileFor(face, ti, tj);
    if (!t || t->draw.empty()) return 0.0f;
    const int x = iu - ti * kTile, y = iv - tj * kTile; // x + 1, y + 1 <= kTile: inside the tile's nodes
    auto at = [&](int a, int b) { return t->draw[(size_t)(y + b) * kNodes + (x + a)]; };
    return (at(0, 0) * (1 - fx) + at(1, 0) * fx) * (1 - fy) + (at(0, 1) * (1 - fx) + at(1, 1) * fx) * fy;
}

void Clear() {
    std::lock_guard<std::mutex> lock(mx);
    cache.clear();
    lru.clear();
    bytes = 0;
}
void Invalidate(const PlanetField& pf, Vector3 dir, float radius) {
    // The circle's centre and 8 points on its rim each give a face and a (u, v) box; every tile of every level touching a
    // box goes (a face's uv is stretched up to ~2x at its corners, hence the margin).
    const float ang = 1.5f * radius / pf.radius;
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(dir, t1);
    std::lock_guard<std::mutex> lock(mx);
    auto cit = cache.find(&pf);
    if (cit == cache.end()) return;
    auto& m = cit->second;
    for (int s = 0; s < 9; ++s) {
        const float a = s * (2.0f * PI / 8.0f);
        const Vector3 d = s == 0 ? dir : Vector3Normalize(Vector3Add(dir, Vector3Add(Vector3Scale(t1, ang * cosf(a)), Vector3Scale(t2, ang * sinf(a)))));
        int face; float u, v;
        FaceOf(d, face, u, v);
        const float du = 2.0f * ang; // in uv (a face is +-1 over +-45 degrees, stretched at its corners)
        for (int level = 0; level <= kMaxLevel; ++level) {
            const int n = std::max(kTile, (int)ceilf(2.0f * pf.radius / (0.5f * (float)(1 << level))));
            auto toTile = [&](float x) { return std::clamp((int)floorf((x + 1.0f) * 0.5f * n), 0, n - 1) / kTile; };
            for (int ti = toTile(u - du); ti <= toTile(u + du); ++ti)
                for (int tj = toTile(v - du); tj <= toTile(v + du); ++tj) {
                    auto e = m.find(KeyOf(level, face, ti, tj));
                    if (e == m.end()) continue;
                    bytes -= BytesOf(*e->second.tile);
                    lru.erase(e->second.lru);
                    m.erase(e);
                }
        }
    }
}
long long CacheBytes() { return bytes.load(); }
long long TilesBuilt() { return built.load(); }

} // namespace WaterSurface
