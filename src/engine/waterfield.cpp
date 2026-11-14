#include "waterfield.h"
#include "voxel.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <mutex>

namespace WaterField {

namespace {
    std::mutex mx;
    SnapshotPtr current = std::make_shared<const Snapshot>();
    std::vector<Dirty> dirtyQueue;
    constexpr float kNone = -1e9f;
    constexpr int kShallow = 64; // fills under 25 cm on the ground: drawn only as part of a sheet (openSurface)

}

void Site::brickBounds(int bi, int bj, int bk, Vector3& blo, Vector3& bhi) const {
    cellBounds(bi * kBrick, bj * kBrick, bk * kBrick, bi * kBrick + kBrick, bj * kBrick + kBrick, bk * kBrick + kBrick, blo, bhi);
}

void Site::cellBounds(int i0, int j0, int k0, int i1, int j1, int k1, Vector3& blo, Vector3& bhi) const {
    blo = { 1e30f, 1e30f, 1e30f }; bhi = { -1e30f, -1e30f, -1e30f };
    for (int c = 0; c < 8; ++c) {
        const float x = (float)((c & 1) ? i1 : i0), y = (float)((c & 2) ? j1 : j0), z = (float)((c & 4) ? k1 : k0);
        const Vector3 p = Vector3Add(origin, Vector3Add(Vector3Add(Vector3Scale(e, x), Vector3Scale(n, y)), Vector3Scale(u, z)));
        blo = Vector3Min(blo, p); bhi = Vector3Max(bhi, p);
    }
}

SnapshotPtr Current() {
    std::lock_guard<std::mutex> lock(mx);
    return current;
}

void Publish(SnapshotPtr snap, std::vector<Dirty>&& dirty) {
    std::lock_guard<std::mutex> lock(mx);
    current = std::move(snap);
    if (dirtyQueue.empty()) dirtyQueue = std::move(dirty);
    else dirtyQueue.insert(dirtyQueue.end(), dirty.begin(), dirty.end());
}

void Clear() {
    std::lock_guard<std::mutex> lock(mx);
    auto empty = std::make_shared<Snapshot>();
    empty->seq = current ? current->seq + 1 : 1;
    current = std::move(empty);
    dirtyQueue.clear();
}

void TakeDirty(std::vector<Dirty>& out) {
    std::lock_guard<std::mutex> lock(mx);
    out.swap(dirtyQueue);
    dirtyQueue.clear();
}

int CellFill(const Snapshot& s, uint64_t body, Vector3 local, int kind, float* frac, uint8_t* flagsOut) {
    for (const Site& site : s.sites) {
        if (site.body != body) continue;
        if (local.x < site.lo.x - 1 || local.y < site.lo.y - 1 || local.z < site.lo.z - 1 || local.x > site.hi.x + 1 || local.y > site.hi.y + 1 || local.z > site.hi.z + 1) continue;
        float fx, fy, fz;
        site.frame(local, fx, fy, fz);
        uint8_t fill, flags;
        if (!site.cell((int)floorf(fx), (int)floorf(fy), (int)floorf(fz), fill, flags)) continue;
        if (frac) *frac = fz - floorf(fz);
        if (flagsOut) *flagsOut = flags;
        if ((flags & kSolid) || CellKind(flags) != kind) return 0;
        return fill;
    }
    return -1;
}

// ---------------------------------------------------------------------------------------------
// Chunk sampler
// ---------------------------------------------------------------------------------------------

ChunkSampler::ChunkSampler(const Snapshot* snap, uint64_t body, Vector3 lo, Vector3 hi) {
    if (!snap) return;
    for (const Site& site : snap->sites) {
        if (site.body != body || site.bricks.empty()) continue;
        if (hi.x < site.lo.x - 2 || hi.y < site.lo.y - 2 || hi.z < site.lo.z - 2 || lo.x > site.hi.x + 2 || lo.y > site.hi.y + 2 || lo.z > site.hi.z + 2) continue;
        // The box in site cells (its corners), one cell of margin: the bricks that may reach it.
        float mn[3] = { 1e30f, 1e30f, 1e30f }, mxv[3] = { -1e30f, -1e30f, -1e30f };
        for (int c = 0; c < 8; ++c) {
            const Vector3 p = { (c & 1) ? hi.x : lo.x, (c & 2) ? hi.y : lo.y, (c & 4) ? hi.z : lo.z };
            float f[3];
            site.frame(p, f[0], f[1], f[2]);
            for (int a = 0; a < 3; ++a) { mn[a] = std::min(mn[a], f[a]); mxv[a] = std::max(mxv[a], f[a]); }
        }
        int b0[3], b1[3];
        for (int a = 0; a < 3; ++a) { b0[a] = FloorDiv((int)floorf(mn[a]) - 1, kBrick); b1[a] = FloorDiv((int)floorf(mxv[a]) + 1, kBrick); }
        const long long count = (long long)(b1[0] - b0[0] + 1) * (b1[1] - b0[1] + 1) * (b1[2] - b0[2] + 1);
        bool found = false;
        uint8_t kinds = 0;
        auto inRange = [&](int bi, int bj, int bk) { return bi >= b0[0] && bi <= b1[0] && bj >= b0[1] && bj <= b1[1] && bk >= b0[2] && bk <= b1[2]; };
        if (count <= (long long)site.bricks.size()) {
            for (int bk = b0[2]; bk <= b1[2]; ++bk)
                for (int bj = b0[1]; bj <= b1[1]; ++bj)
                    for (int bi = b0[0]; bi <= b1[0]; ++bi)
                        if (const Brick* b = site.brick(bi, bj, bk)) { found = true; kinds |= b->kinds; }
        } else {
            for (const auto& kv : site.bricks) {
                const int bi = (int)((kv.first >> 22) & 2047) - 1024, bj = (int)((kv.first >> 11) & 2047) - 1024, bk = (int)(kv.first & 2047) - 1024;
                if (inRange(bi, bj, bk)) { found = true; kinds |= kv.second->kinds; }
            }
        }
        if (!found) continue;
        sites.push_back(&site);
        kindMask |= kinds;
    }
    if (!sites.empty()) kindMask |= 1; // water: the published cells may hold none of it here and still decide (a drained pond)
}

float ChunkSampler::openSurface(int si, int ci, int cj, int ck, int kind, int depth) {
    const Site& s = *sites[si];
    uint8_t fill, flags;
    if (!s.cell(ci, cj, ck, fill, flags) || (flags & kSolid)) return kNone;
    const int f = CellKind(flags) == kind ? fill : 0;
    if (f > 0) {
        if (f < kShallow) {
            // Shallow water on no fluid (the ground is in this cell): drawn as part of a sheet - a side neighbour at the
            // same level, or deeper water - not as the films a draining body leaves cell by cell, each at its own height
            // over a bed the cell does not know (they read as a grid of spikes, #360).
            uint8_t fb, flb;
            const bool topsFluid = s.cell(ci, cj, ck - 1, fb, flb) && !(flb & kSolid) && fb >= 255 && CellKind(flb) == kind;
            if (!topsFluid) {
                static const int side[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
                bool sheet = false;
                for (const auto& d : side) {
                    uint8_t fo, flo;
                    if (!s.cell(ci + d[0], cj + d[1], ck, fo, flo) || (flo & kSolid) || fo == 0 || CellKind(flo) != kind) continue;
                    if (fo >= kShallow || abs((int)fo - f) <= 8) { sheet = true; break; }
                }
                if (!sheet) return kNone;
            }
        }
        if (f < 255) return ck + f / 255.0f;
        // Full: the surface is higher up if the cell above holds this fluid too.
        uint8_t fa, fla;
        if (depth < 8 && s.cell(ci, cj, ck + 1, fa, fla) && !(fla & kSolid) && fa > 0 && CellKind(fla) == kind) return openSurface(si, ci, cj, ck + 1, kind, depth + 1);
        return ck + 1.0f;
    }
    // Empty: over fluid, the surface below (so the interpolated surface runs level across the cell border).
    uint8_t fb, flb;
    if (s.cell(ci, cj, ck - 1, fb, flb) && !(flb & kSolid) && fb > 0 && CellKind(flb) == kind) return (ck - 1) + fb / 255.0f;
    return kNone;
}

float ChunkSampler::surface(int si, int ci, int cj, int ck, int kind) {
    const int64_t key = ((int64_t)si << 56) | ((int64_t)kind << 52) | ((int64_t)((ci + 8192) & 0x3FFF) << 28) | ((int64_t)((cj + 8192) & 0x3FFF) << 14) | (int64_t)((ck + 8192) & 0x3FFF);
    auto it = memo.find(key);
    if (it != memo.end()) return it->second;
    const Site& s = *sites[si];
    uint8_t fill, flags;
    float r = kNone;
    if (s.cell(ci, cj, ck, fill, flags)) {
        const bool solid = flags & kSolid;
        if (!solid) r = openSurface(si, ci, cj, ck, kind, 0);
        uint8_t fb, flb;
        // A solid cell (its centre is in the ground), or a dry open cell resting on the ground: the wettest side
        // neighbour's surface runs on into it, and the terrain (the chunk's own density) cuts it. So the water meets the
        // real shore instead of stopping at a cell's edge as a box (#191; shallow puddles read as a grid of spikes). Not
        // over a ledge (the cell below open): that water is falling, not lying there.
        if (r <= kNone * 0.5f && (solid || (s.cell(ci, cj, ck - 1, fb, flb) && (flb & kSolid)))) {
            static const int side[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
            for (const auto& d : side) r = std::max(r, openSurface(si, ci + d[0], cj + d[1], ck, kind, 0));
            // Only a surface inside this cell's height: deeper water beside it is under the ground here.
            if (r > ck + 1.0f) r = ck + 1.0f;
            if (r < (float)ck) r = kNone;
        }
        // A solid cell under water (the bed cuts through the cell above): the water goes on down into it, so the samples
        // below the bed's cell stay "in water" and no second, upside-down surface forms at the cell border (it showed
        // through shallow sheets as a regular grid of holes where the frame's cells cross the chunk's samples).
        if (solid) {
            const float up = openSurface(si, ci, cj, ck + 1, kind, 0);
            if (up >= ck + 1.0f) r = std::max(r, up);
        }
    }
    memo.emplace(key, r);
    return r;
}

bool ChunkSampler::density(Vector3 p, int kind, float& wd) {
    for (int si = 0; si < (int)sites.size(); ++si) {
        const Site& s = *sites[si];
        float fx, fy, fz;
        s.frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        if (!s.brick(FloorDiv(ci, kBrick), FloorDiv(cj, kBrick), FloorDiv(ck, kBrick))) continue;
        const float sf = surface(si, ci, cj, ck, kind);
        wd = sf <= kNone * 0.5f ? ISO_LEVEL - DENSITY_GRADIENT * 0.5f : ISO_LEVEL + DENSITY_GRADIENT * (sf - fz);
        return true;
    }
    return false;
}

bool ChunkSampler::covered(Vector3 p) const {
    for (const Site* s : sites) {
        float fx, fy, fz;
        s->frame(p, fx, fy, fz);
        if (s->brick(FloorDiv((int)floorf(fx), kBrick), FloorDiv((int)floorf(fy), kBrick), FloorDiv((int)floorf(fz), kBrick))) return true;
    }
    return false;
}

float ChunkSampler::depthBelow(Vector3 p, int kind) {
    for (const Site* s : sites) {
        float fx, fy, fz;
        s->frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy);
        int ck = (int)floorf(fz);
        uint8_t fill, flags;
        if (!s->cell(ci, cj, ck, fill, flags)) continue;
        float d = 0.0f;
        for (int k = 0; k < 8; ++k, --ck) {
            if (!s->cell(ci, cj, ck, fill, flags) || (flags & kSolid)) break;
            if (CellKind(flags) != kind || fill == 0) { if (k > 0) break; else continue; }
            d += fill / 255.0f;
        }
        return d;
    }
    return 0.0f;
}

} // namespace WaterField
