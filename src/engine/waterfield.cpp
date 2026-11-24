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
    const int kFilmLift = std::min(128, (int)lroundf(13.0f / CellSize())); // a moving film is drawn at least ~5 cm over its ground share (#649)

}

// The fluid cell edge (m): 0.25 by default (#714, user 2026-11-17: at 1 m a one-block hole in a block dam stayed shut, #709,
// and a dug trench's rounded edge spilled); RAYLIBGAME_FLUID_CELL=1 is the old 1 m grid, opt-in.
float CellSize() {
    static const float c = [] { const char* e = getenv("RAYLIBGAME_FLUID_CELL"); const float v = e ? (float)atof(e) : 0.25f; return v > 0.05f && v <= 4.0f ? v : 0.25f; }();
    return c;
}

void Site::brickBounds(int bi, int bj, int bk, Vector3& blo, Vector3& bhi) const {
    cellBounds(bi * kBrick, bj * kBrick, bk * kBrick, bi * kBrick + kBrick, bj * kBrick + kBrick, bk * kBrick + kBrick, blo, bhi);
}

void Site::cellBounds(int i0, int j0, int k0, int i1, int j1, int k1, Vector3& blo, Vector3& bhi) const {
    blo = { 1e30f, 1e30f, 1e30f }; bhi = { -1e30f, -1e30f, -1e30f };
    for (int c = 0; c < 8; ++c) {
        const float cs = CellSize();
        const float x = (float)((c & 1) ? i1 : i0) * cs, y = (float)((c & 2) ? j1 : j0) * cs, z = (float)((c & 4) ? k1 : k0) * cs;
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

ChunkSampler::ChunkSampler(const Snapshot* snap, uint64_t body, Vector3 lo, Vector3 hi, bool probeOnly) {
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
        if (probeOnly) return; // any() is all the caller asks (#381 review: the caches cost up to ~1 MB a call)
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
    fill = b->drawn[i]; flags = b->flags[i]; // the mesher draws the eased water (#366)
    return true;
}
bool ChunkSampler::cellAt(int si, int ci, int cj, int ck, uint8_t& fill, uint8_t& flags, uint8_t& shut) const {
    const Brick* b = brickPtr(si, ci >> 4, cj >> 4, ck >> 4);
    if (!b) return false;
    const int i = CellIdx(ci & 15, cj & 15, ck & 15);
    fill = b->drawn[i]; flags = b->flags[i]; shut = b->shut[i];
    return true;
}

float ChunkSampler::openSurface(int si, int ci, int cj, int ck, int kind, int depth) {
    uint8_t fill, flags, sh; // sh: the cell's ground share: its water stands at sh + fill (#454)
    if (!cellAt(si, ci, cj, ck, fill, flags, sh) || (flags & kSolid)) return kNone;
    const int f = CellKind(flags) == kind ? fill : 0;
    if (f > 0) {
        const Brick* own = brickPtr(si, ci >> 4, cj >> 4, ck >> 4);
        if (f < kShallow && own && own->moving) {
            // A moving film (#649): drawn wherever it is, at least kFilmLift over its ground share - a film running down a
            // trench or across dry ground showed nothing until it pooled (the water seemed to skip to where it settled).
            // Resting films keep the sheet rule below (a drained pond's films read as a grid of spikes, #360).
            return ck + (sh + std::max(f, kFilmLift)) / 255.0f;
        }
        if (f < kShallow) {
            // Shallow water on no fluid (the ground is in this cell): drawn as part of a sheet - a side neighbour at the
            // same level, or deeper water - not as the films a draining body leaves cell by cell, each at its own height
            // over a bed the cell does not know (they read as a grid of spikes, #360).
            uint8_t fb, flb, shb;
            const bool topsFluid = cellAt(si, ci, cj, ck - 1, fb, flb, shb) && !(flb & kSolid) && fb + shb >= 255 && CellKind(flb) == kind;
            if (!topsFluid) {
                static const int side[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
                bool sheet = false;
                for (const auto& d : side) {
                    uint8_t fo, flo, sho;
                    if (!cellAt(si, ci + d[0], cj + d[1], ck, fo, flo, sho) || (flo & kSolid) || fo == 0 || CellKind(flo) != kind) continue;
                    if (fo >= kShallow || abs(((int)sho + fo) - ((int)sh + f)) <= 8) { sheet = true; break; }
                }
                if (!sheet) return kNone;
            }
        }
        if (sh + f < 255) return ck + (sh + f) / 255.0f;
        // Full: the surface is higher up if the cell above holds this fluid too.
        uint8_t fa, fla;
        if (depth < 8 && cellAt(si, ci, cj, ck + 1, fa, fla) && !(fla & kSolid) && fa > 0 && CellKind(fla) == kind) return openSurface(si, ci, cj, ck + 1, kind, depth + 1);
        return ck + 1.0f;
    }
    // Empty: over fluid, the surface below (so the interpolated surface runs level across the cell border).
    uint8_t fb, flb, shb;
    if (cellAt(si, ci, cj, ck - 1, fb, flb, shb) && !(flb & kSolid) && fb > 0 && CellKind(flb) == kind) return (ck - 1) + (shb + fb) / 255.0f;
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

// RAYLIBGAME_SHORE_SIM=1: published generated water is drawn by the simulation's rules too (before #1091), for A/B.
static const bool kShoreSim = [] { const char* e = getenv("RAYLIBGAME_SHORE_SIM"); return e && e[0] == '1'; }();

// #1091: a dig by the sea pulls the shore's bricks into a site; their sea cells (kSource | kStatic, held at the sea's level)
// were then drawn with the simulation's rules for films and draining sheets - over a beach at sea level, a moving brick's
// films lifted over the sand and tilted toward dry neighbours (shards), a resting one's 15 cm skin and sheet rule cut the
// shallow band into jagged pieces, and the drawing switched rules at every brick edge. Water the simulation holds at the
// generated level is the generated water: the mesher draws it so (the same level, so it meets the simulated water - the
// flooded hole, dug cells - without a seam). A drained column (no wet cell) stays the simulation's (a pond emptied by bucket),
// except in a brick that holds generated water at its level (#1590): there a dry undug column is land the dig's bricks reached
// beside the sea, lake or river. Drawn as the simulation's (dry), it pulled the generated surface next to it down into walls and
// shards along the beach; the generated water's own rules (its level, clipped by the ground and the original bed) draw nothing
// there. Held water never drains, and pond bricks hold none (their water is finite), so an emptied pond stays empty.
bool ChunkSampler::generatedColumn(int si, int ci, int cj, int ck) const {
    if (kShoreSim) return false;
    const int bk = ck >> 4;
    const int64_t key = ((int64_t)si << 56) | ((int64_t)((ci + 8192) & 0x3FFF) << 28) | ((int64_t)((cj + 8192) & 0x3FFF) << 14) | (int64_t)((bk + 8192) & 0x3FFF);
    auto it = genMemo.find(key);
    if (it != genMemo.end()) return it->second;
    bool wet = false, gen = true;
    uint8_t fill, flags;
    for (int k = bk * 16; k < bk * 16 + 16 && gen; ++k) {
        if (!cellAt(si, ci, cj, k, fill, flags) || (flags & kSolid)) continue;
        if (flags & kDug) gen = false;                                   // an edit opened it: the simulation's water
        else if (fill > 0) { if ((flags & kStatic) && CellKind(flags) == 0) wet = true; else gen = false; }
    }
    if (gen && !wet) gen = holdsGenerated(si, ci >> 4, cj >> 4, bk);
    const bool r = gen;
    genMemo.emplace(key, r);
    return r;
}

bool ChunkSampler::holdsGenerated(int si, int bi, int bj, int bk) const {
    const int64_t key = ((int64_t)si << 56) | ((int64_t)((bi + 8192) & 0x3FFF) << 28) | ((int64_t)((bj + 8192) & 0x3FFF) << 14) | (int64_t)((bk + 8192) & 0x3FFF);
    auto it = heldMemo.find(key);
    if (it != heldMemo.end()) return it->second;
    bool held = false;
    if (const Brick* b = brickPtr(si, bi, bj, bk))
        for (int i = 0; i < kBrick * kBrick * kBrick && !held; ++i) held = b->fill[i] > 0 && (b->flags[i] & kStatic) && CellKind(b->flags[i]) == 0;
    heldMemo.emplace(key, held);
    return held;
}

bool ChunkSampler::density(Vector3 p, int kind, float& wd) {
    for (int si = 0; si < (int)sites.size(); ++si) {
        const Site& s = *sites[si];
        float fx, fy, fz;
        s.frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        if (!brickAt(si, ci, cj, ck)) continue;
        if (kind == 0 && generatedColumn(si, ci, cj, ck)) continue; // drawn as the generated water (#1091)
        // The density of one column's surface at this height (dry: as if the surface were half a cell under the sample).
        const float cs = CellSize();
        auto columnDensity = [&](int i, int j) {
            const float sf = surface(si, i, j, ck, kind);
            return sf <= kNone * 0.5f ? ISO_LEVEL - DENSITY_GRADIENT * 0.5f * cs : ISO_LEVEL + DENSITY_GRADIENT * (sf - fz) * cs; // cells -> m
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

bool ChunkSampler::moving(Vector3 p) const {
    for (int si = 0; si < (int)sites.size(); ++si) {
        float fx, fy, fz;
        sites[si]->frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        if (const Brick* b = brickPtr(si, ci >> 4, cj >> 4, ck >> 4)) return b->moving;
    }
    return false;
}

bool ChunkSampler::covered(Vector3 p) const {
    for (int si = 0; si < (int)sites.size(); ++si) {
        float fx, fy, fz;
        sites[si]->frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        // The dense brick cache (#381); generated water held at its level is the generated water's to draw (#1091).
        if (brickAt(si, ci, cj, ck) && !generatedColumn(si, ci, cj, ck)) return true;
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
        for (int k = 0; k < (int)(8.0f / CellSize()); ++k, --ck) {
            if (!cellAt(si, ci, cj, ck, fill, flags) || (flags & kSolid)) break;
            if (CellKind(flags) != kind || fill == 0) { if (k > 0) break; else continue; }
            d += fill / 255.0f * CellSize(); // m
        }
        return d;
    }
    return 0.0f;
}

bool ChunkSampler::slope(Vector3 p, int kind, Vector3& downhill, float& grade) {
    for (int si = 0; si < (int)sites.size(); ++si) {
        const Site& s = *sites[si];
        float fx, fy, fz;
        s.frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        if (!brickAt(si, ci, cj, ck)) continue;
        const float own = surface(si, ci, cj, ck, kind);
        if (own <= kNone * 0.5f) return false;
        // Surfaces in cells (heights and spacing alike, so the ratio is the grade), across ~1 m (k cells each side; a dry or
        // missing side: one-sided, then the next cell in). One fill unit of quantisation over one cell read as a 0.4% grade -
        // Manning on 2 m of pool made that 3 m/s, so pools swirled red in the debug view (#695); over 2k cells and less a
        // dead zone of 1.5 units, a level pool stands still and a real channel grade (~0.3% and up) stays.
        const int k = std::max(1, (int)lroundf(1.0f / CellSize()));
        auto diff = [&](int di, int dj) {
            for (int r = k; r >= 1; r = r > 1 ? std::max(1, r / 2) : 0) {
                const float a = surface(si, ci + di * r, cj + dj * r, ck, kind), b = surface(si, ci - di * r, cj - dj * r, ck, kind);
                const bool ha = a > kNone * 0.5f, hb = b > kNone * 0.5f;
                if (ha && hb) return (a - b) / (2.0f * r);
                if (ha) return (a - own) / r;
                if (hb) return (own - b) / r;
            }
            return 0.0f;
        };
        const float gx = diff(1, 0), gy = diff(0, 1);
        const float raw = sqrtf(gx * gx + gy * gy), noise = 1.5f / (255.0f * 2.0f * k);
        grade = std::max(0.0f, raw - noise);
        if (grade < 1e-5f) { downhill = { 0, 0, 0 }; grade = 0.0f; return true; }
        downhill = Vector3Scale(Vector3Add(Vector3Scale(s.e, -gx), Vector3Scale(s.n, -gy)), 1.0f / raw);
        return true;
    }
    return false;
}

} // namespace WaterField

namespace WaterField {
bool ChunkSampler::simFlow(Vector3 p, Vector3& dir, float& speed) const {
    for (int si = 0; si < (int)sites.size(); ++si) {
        const Site& s = *sites[si];
        float fx, fy, fz;
        s.frame(p, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        float best = -1.0f, bx = 0.0f, by = 0.0f;
        for (int dk = 0; dk >= -1; --dk) {
            const Brick* b = brickPtr(si, ci >> 4, cj >> 4, (ck + dk) >> 4);
            if (!b || b->flow.empty()) continue;
            const int i = (ci & 15) + 16 * (cj & 15); // per column (#677)
            const float vx = b->flow[2 * i] * 0.01f, vy = b->flow[2 * i + 1] * 0.01f, v = sqrtf(vx * vx + vy * vy);
            if (v > best) { best = v; bx = vx; by = vy; }
        }
        if (best < 0.0f) continue;
        speed = best;
        const Vector3 d = Vector3Add(Vector3Scale(s.e, bx), Vector3Scale(s.n, by));
        const float dl = Vector3Length(d);
        dir = dl > 1e-6f ? Vector3Scale(d, 1.0f / dl) : Vector3{ 0, 0, 0 };
        return true;
    }
    return false;
}
} // namespace WaterField

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
using namespace WaterField;

// #1091: a published brick's column of generated sea (kStatic sources) is left to the generated water; a column holding a
// dug (kDug) cell, simulated water (no kStatic) or none at all is the simulation's to draw.
TEST(WaterField, GeneratedSeaColumnsAreDrawnAsGeneratedWater) {
    auto b = std::make_shared<Brick>();
    const int sea = 8; // cells 0-7 under the sea's level
    for (int y = 0; y < kBrick; ++y)
        for (int x = 0; x < kBrick; ++x)
            for (int z = 0; z < sea; ++z) {
                const int i = CellIdx(x, y, z);
                uint8_t fl = kSource | kStatic;
                if (x == 5) fl = kDug;   // column 5: dug, flooded from the sea
                if (x == 9) fl = 0;      // column 9: water the simulation moved (no generated level)
                if (x == 12) { b->flags[i] = kSolid; continue; } // column 12: the beach, its ground over the sea's level (#1590)
                b->flags[i] = fl; b->fill[i] = b->drawn[i] = 255;
            }
    b->kinds = 1;
    auto pond = std::make_shared<Brick>(); // a pond's brick, emptied by bucket: finite water, none held, nothing left
    pond->kinds = 1;
    Snapshot snap;
    Site s;
    s.body = 7; s.origin = { 0, 0, 0 }; s.e = { 1, 0, 0 }; s.n = { 0, 0, 1 }; s.u = { 0, 1, 0 };
    const float cs = CellSize(), side = kBrick * cs;
    s.lo = { 0, 0, 0 }; s.hi = { 2 * side, side, side };
    s.bricks[BrickKey(0, 0, 0)] = b;
    s.bricks[BrickKey(1, 0, 0)] = pond;
    snap.sites.push_back(s);
    ChunkSampler sampler(&snap, 7, { 0, 0, 0 }, { 2 * side, side, side });
    auto at = [&](int x) { return Vector3{ (x + 0.5f) * cs, (sea - 0.5f) * cs, 3.5f * cs }; }; // a sample just under the sea's surface
    auto over = [&](int x) { return Vector3{ (x + 0.5f) * cs, (sea + 0.5f) * cs, 3.5f * cs }; }; // just over it
    float w = 0.0f;
    EXPECT_FALSE(sampler.density(at(2), 0, w)) << "generated sea: the generated water draws it";
    EXPECT_FALSE(sampler.covered(at(2)));
    EXPECT_TRUE(sampler.density(at(5), 0, w)) << "a dug, flooded column";
    EXPECT_TRUE(sampler.covered(at(5)));
    EXPECT_TRUE(sampler.density(at(9), 0, w)) << "simulated water";
    EXPECT_FALSE(sampler.covered(over(12))) << "dry land beside the sea: the generated water's rules (#1590)";
    EXPECT_TRUE(sampler.density(at(kBrick + 4), 0, w)) << "an emptied pond stays empty";
    EXPECT_TRUE(sampler.covered(at(kBrick + 4)));
}

// #1773: water in a corner pen, its two walls meeting only diagonally at the dry cell (6,6). The simulation blocks the
// flow there; the drawn surface must not cross the corner into that cell either (the blend pulled it in).
TEST(WaterField, SurfaceStopsAtADiagonalCorner) {
    auto b = std::make_shared<Brick>();
    for (int cj = 0; cj < kBrick; ++cj)
        for (int ci = 0; ci < kBrick; ++ci)
            for (int ck = 0; ck < 8; ++ck) {
                const int i = CellIdx(ci, cj, ck);
                const bool wall = (ci == 6 && cj < 6) || (cj == 6 && ci < 6);
                if (ck < 4 || wall) { b->flags[i] = kSolid; continue; } // the ground, and the walls (tiles) on it
                if (ci < 6 && cj < 6 && ck < 7) b->fill[i] = b->drawn[i] = ck < 6 ? 255 : 128; // the pen: its surface at 6.5 cells
            }
    b->kinds = 1;
    Snapshot snap;
    Site s;
    s.body = 7; s.origin = { 0, 0, 0 }; s.e = { 1, 0, 0 }; s.n = { 0, 0, 1 }; s.u = { 0, 1, 0 }; // cells: i = x, j = z, k = y
    const float cs = CellSize(), side = kBrick * cs;
    s.lo = { 0, 0, 0 }; s.hi = { side, side, side };
    s.bricks[BrickKey(0, 0, 0)] = b;
    snap.sites.push_back(s);
    ChunkSampler sampler(&snap, 7, { 0, 0, 0 }, { side, side, side });
    auto cell = [&](float i, float j, float k) { return Vector3{ i * cs, k * cs, j * cs }; };
    float w = 0.0f;
    ASSERT_TRUE(sampler.density(cell(6.05f, 6.05f, 6.1f), 0, w));
    EXPECT_LT(w, ISO_LEVEL) << "under the pen's level just past the corner: the dry cell stays dry";
    ASSERT_TRUE(sampler.density(cell(5.95f, 5.95f, 6.1f), 0, w));
    EXPECT_GT(w, ISO_LEVEL) << "the pen's own corner cell is still water";
    ASSERT_TRUE(sampler.density(cell(3.5f, 3.5f, 6.9f), 0, w));
    EXPECT_LT(w, ISO_LEVEL) << "over the pen's surface";
}
#endif
