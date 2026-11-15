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
    blend = Vector3Distance(lo, hi) < 70.0f; // a box up to ~40 m across (levels 0-2): coarser samples are wider apart than a cell (#381)
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
        // Dense caches over the box (+2 cells: the blend reads the next column, surfaces read their side neighbours).
        Dense d;
        for (int a = 0; a < 3; ++a) { d.c0[a] = (int)floorf(mn[a]) - 2; d.n[a] = (int)floorf(mxv[a]) + 3 - d.c0[a]; }
        if ((long long)d.n[0] * d.n[1] * d.n[2] <= 262144) d.surf.assign((size_t)d.n[0] * d.n[1] * d.n[2], std::nanf(""));
        else d.n[0] = d.n[1] = d.n[2] = 0; // a coarse chunk: the memo map only
        for (int a = 0; a < 3; ++a) { d.b0[a] = FloorDiv(d.c0[a], kBrick); d.bn[a] = FloorDiv(d.c0[a] + std::max(1, d.n[a]) - 1, kBrick) - d.b0[a] + 1; }
        if ((long long)d.bn[0] * d.bn[1] * d.bn[2] <= 4096) {
            d.brick.assign((size_t)d.bn[0] * d.bn[1] * d.bn[2], nullptr);
            for (int bk = 0; bk < d.bn[2]; ++bk)
                for (int bj = 0; bj < d.bn[1]; ++bj)
                    for (int bi = 0; bi < d.bn[0]; ++bi)
                        d.brick[((size_t)bk * d.bn[1] + bj) * d.bn[0] + bi] = site.brick(d.b0[0] + bi, d.b0[1] + bj, d.b0[2] + bk);
        }
        dense.push_back(std::move(d));
    }
    if (!sites.empty()) kindMask |= 1; // water: the published cells may hold none of it here and still decide (a drained pond)
}

const Brick* ChunkSampler::brickPtr(int si, int bi, int bj, int bk) const {
    const Dense& d = dense[si];
    const int x = bi - d.b0[0], y = bj - d.b0[1], z = bk - d.b0[2];
    if (!d.brick.empty() && x >= 0 && y >= 0 && z >= 0 && x < d.bn[0] && y < d.bn[1] && z < d.bn[2])
        return d.brick[((size_t)z * d.bn[1] + y) * d.bn[0] + x];
    return sites[si]->brick(bi, bj, bk);
}

// Site::cell without the hash lookup inside the chunk's box (#381: openSurface reads ~10 cells per column).
bool ChunkSampler::cellAt(int si, int ci, int cj, int ck, uint8_t& fill, uint8_t& flags) const {
    const Brick* b = brickPtr(si, ci >> 4, cj >> 4, ck >> 4);
    if (!b) return false;
    const int i = CellIdx(ci & 15, cj & 15, ck & 15);
    fill = b->fill[i]; flags = b->flags[i];
    return true;
}

float ChunkSampler::openSurface(int si, int ci, int cj, int ck, int kind, int depth) {
    uint8_t fill, flags;
    if (!cellAt(si, ci, cj, ck, fill, flags) || (flags & kSolid)) return kNone;
    const int f = CellKind(flags) == kind ? fill : 0;
    if (f > 0) {
        if (f < kShallow) {
            // Shallow water on no fluid (the ground is in this cell): drawn as part of a sheet - a side neighbour at the
            // same level, or deeper water - not as the films a draining body leaves cell by cell, each at its own height
            // over a bed the cell does not know (they read as a grid of spikes, #360).
            uint8_t fb, flb;
            const bool topsFluid = cellAt(si, ci, cj, ck - 1, fb, flb) && !(flb & kSolid) && fb >= 255 && CellKind(flb) == kind;
            if (!topsFluid) {
                static const int side[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
                bool sheet = false;
                for (const auto& d : side) {
                    uint8_t fo, flo;
                    if (!cellAt(si, ci + d[0], cj + d[1], ck, fo, flo) || (flo & kSolid) || fo == 0 || CellKind(flo) != kind) continue;
                    if (fo >= kShallow || abs((int)fo - f) <= 8) { sheet = true; break; }
                }
                if (!sheet) return kNone;
            }
        }
        if (f < 255) return ck + f / 255.0f;
        // Full: the surface is higher up if the cell above holds this fluid too.
        uint8_t fa, fla;
        if (depth < 8 && cellAt(si, ci, cj, ck + 1, fa, fla) && !(fla & kSolid) && fa > 0 && CellKind(fla) == kind) return openSurface(si, ci, cj, ck + 1, kind, depth + 1);
        return ck + 1.0f;
    }
    // Empty: over fluid, the surface below (so the interpolated surface runs level across the cell border).
    uint8_t fb, flb;
    if (cellAt(si, ci, cj, ck - 1, fb, flb) && !(flb & kSolid) && fb > 0 && CellKind(flb) == kind) return (ck - 1) + fb / 255.0f;
    return kNone;
}

float ChunkSampler::surface(int si, int ci, int cj, int ck, int kind) {
    // Water in the chunk's box: the dense cache (#423); other kinds and cells outside it: the memo map.
    float* slot = nullptr;
    if (kind == 0 && !dense[si].surf.empty()) {
        const Dense& d = dense[si];
        const int x = ci - d.c0[0], y = cj - d.c0[1], z = ck - d.c0[2];
        if (x >= 0 && y >= 0 && z >= 0 && x < d.n[0] && y < d.n[1] && z < d.n[2]) {
            slot = &dense[si].surf[((size_t)z * d.n[1] + y) * d.n[0] + x];
            if (!std::isnan(*slot)) return *slot;
        }
    }
    const int64_t key = ((int64_t)si << 56) | ((int64_t)kind << 52) | ((int64_t)((ci + 8192) & 0x3FFF) << 28) | ((int64_t)((cj + 8192) & 0x3FFF) << 14) | (int64_t)((ck + 8192) & 0x3FFF);
    if (!slot) {
        auto it = memo.find(key);
        if (it != memo.end()) return it->second;
    }
    uint8_t fill, flags;
    float r = kNone;
    if (cellAt(si, ci, cj, ck, fill, flags)) {
        const bool solid = flags & kSolid;
        if (!solid) r = openSurface(si, ci, cj, ck, kind, 0);
        uint8_t fb, flb;
        // A solid cell (its centre is in the ground), or a dry open cell resting on the ground: the wettest side
        // neighbour's surface runs on into it, and the terrain (the chunk's own density) cuts it. So the water meets the
        // real shore instead of stopping at a cell's edge as a box (#191; shallow puddles read as a grid of spikes). Not
        // over a ledge (the cell below open): that water is falling, not lying there.
        if (r <= kNone * 0.5f && (solid || (cellAt(si, ci, cj, ck - 1, fb, flb) && (flb & kSolid)))) {
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
    if (slot) *slot = r; else memo.emplace(key, r);
    return r;
}

bool ChunkSampler::density(Vector3 p, int kind, float& wd) {
    for (int si = 0; si < (int)sites.size(); ++si) {
        const Site& s = *sites[si];
        float fx, fy, fz;
        s.frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        if (!brickAt(si, ci, cj, ck)) continue;
        // The density of one column's surface at this height (dry: as if the surface were half a cell under the sample).
        auto columnDensity = [&](int i, int j) {
            const float sf = surface(si, i, j, ck, kind);
            return sf <= kNone * 0.5f ? ISO_LEVEL - DENSITY_GRADIENT * 0.5f : ISO_LEVEL + DENSITY_GRADIENT * (sf - fz);
        };
        const float own = columnDensity(ci, cj);
        // Far from its column's surface (2.5 m+ over or under it), the blend can't put the sample on the other side
        // (only neighbours 2.5 m higher or lower could: a waterfall's step stays a little sharper): its own value does
        // (#381: the blend made a water rebuild 1.7x dearer; most samples are well above or below the surface).
        if (!blend || fabsf(own - ISO_LEVEL) > DENSITY_GRADIENT * 2.5f) { wd = own; return true; }
        // Blended across the four nearest columns (#423): each 1 m column was a flat-topped box, so the water showed 1 m
        // steps and square blobs. Linear weights from the cell centres slope the surface evenly between columns (smoothstep
        // weights put each column's whole step at its border: a 1 m grid of ripples on a draining pond), and a dry
        // neighbour pulls the edge down, rounding it off like surface tension. A column in a brick
        // that is not published (generated water beyond the simulated part) counts as this one, so no dip forms there.
        const float gx = fx - 0.5f, gy = fy - 0.5f;
        const int i0 = (int)floorf(gx), j0 = (int)floorf(gy);
        const float tx = gx - (float)i0, ty = gy - (float)j0;
        float sum = 0.0f;
        for (int b = 0; b < 2; ++b)
            for (int a = 0; a < 2; ++a) {
                const int i = i0 + a, j = j0 + b;
                const float w = (a ? tx : 1.0f - tx) * (b ? ty : 1.0f - ty);
                if (w <= 0.0f) continue;
                const bool same = i == ci && j == cj;
                const bool known = same || brickAt(si, i, j, ck);
                sum += w * (same ? own : known ? columnDensity(i, j) : own);
            }
        wd = sum;
        return true;
    }
    return false;
}

bool ChunkSampler::covered(Vector3 p) const {
    for (int si = 0; si < (int)sites.size(); ++si) {
        float fx, fy, fz;
        sites[si]->frame(p, fx, fy, fz);
        if (brickAt(si, (int)floorf(fx), (int)floorf(fy), (int)floorf(fz))) return true; // the dense brick cache (#381)
    }
    return false;
}

float ChunkSampler::depthBelow(Vector3 p, int kind) {
    for (int si = 0; si < (int)sites.size(); ++si) {
        float fx, fy, fz;
        sites[si]->frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy);
        int ck = (int)floorf(fz);
        uint8_t fill, flags;
        if (!cellAt(si, ci, cj, ck, fill, flags)) continue;
        float d = 0.0f;
        for (int k = 0; k < 8; ++k, --ck) {
            if (!cellAt(si, ci, cj, ck, fill, flags) || (flags & kSolid)) break;
            if (CellKind(flags) != kind || fill == 0) { if (k > 0) break; else continue; }
            d += fill / 255.0f;
        }
        return d;
    }
    return 0.0f;
}

} // namespace WaterField
