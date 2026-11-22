#include "caveislands.h"
#include "voxel.h"

#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace {
struct BoxResult { long long pieces = 0, volume = 0, largest = 0; };

// The density on a grid of `step` m over a box +-hx across (both ways) and +-hy up around c (local frame); solid cells
// connected (6-neighbour) to the box's boundary are the surrounding rock, any other solid component floats. Volumes in
// m^3 (cells x step^3). hb / db (optional): m^3 by mean height (-28..+28 m, 8 bins) and distance across (0..72 m, 8 bins).
BoxResult MeasureBox(const PlanetField& pf, Vector3 c, int hx, int hy, int step, long long* hb, long long* db, double binH = 7.0, double binR = 9.0) {
    BoxResult out;
    const int NX = 2 * hx / step + 1, NY = 2 * hy / step + 1, NZ = NX;   // x, z across; y up
    const Vector3 up = Vector3Normalize(c);
    const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })), e2 = Vector3CrossProduct(up, e1);
    const float fx = (float)hx, fy = (float)hy, fs = (float)step;
    std::vector<uint8_t> solid((size_t)NX * NY * NZ), seen((size_t)NX * NY * NZ, 0);
    auto idx = [&](int x, int y, int z) { return ((size_t)y * NZ + z) * NX + x; };
    for (int y = 0; y < NY; ++y)
        for (int z = 0; z < NZ; ++z)
            for (int x = 0; x < NX; ++x) {
                // (step 1: x - hx exactly as the bench always computed it)
                const float ax = step == 1 ? x - fx : x * fs - fx, ay = step == 1 ? y - fy : y * fs - fy, az = step == 1 ? z - fx : z * fs - fx;
                const Vector3 q = Vector3Add(c, Vector3Add(Vector3Scale(e1, ax), Vector3Add(Vector3Scale(up, ay), Vector3Scale(e2, az))));
                solid[idx(x, y, z)] = pf.baseDensity(q) >= ISO_LEVEL;
            }
    std::vector<int> stack;
    double sumY = 0.0, sumR = 0.0; // the last flooded component's summed height (m from the centre) and distance across
    const double cx = (NX - 1) * 0.5, cy = (NY - 1) * 0.5;
    auto flood = [&](int sx, int sy, int sz) {                  // marks one solid component, returns its cells
        long long n = 0;
        stack.clear(); stack.push_back((int)idx(sx, sy, sz)); seen[idx(sx, sy, sz)] = 1;
        while (!stack.empty()) {
            const int k = stack.back(); stack.pop_back(); ++n;
            const int x = k % NX, z = (k / NX) % NZ, y = k / (NX * NZ);
            sumY += (y - cy) * step; sumR += sqrt((x - cx) * (x - cx) + (z - cx) * (z - cx)) * step;
            const int nb[6][3] = { { x - 1, y, z }, { x + 1, y, z }, { x, y - 1, z }, { x, y + 1, z }, { x, y, z - 1 }, { x, y, z + 1 } };
            for (const auto& b : nb) {
                if (b[0] < 0 || b[1] < 0 || b[2] < 0 || b[0] >= NX || b[1] >= NY || b[2] >= NZ) continue;
                const size_t j = idx(b[0], b[1], b[2]);
                if (solid[j] && !seen[j]) { seen[j] = 1; stack.push_back((int)j); }
            }
        }
        return n;
    };
    for (int y = 0; y < NY; ++y)                                // the surrounding rock: from every boundary cell
        for (int z = 0; z < NZ; ++z)
            for (int x = 0; x < NX; ++x)
                if ((x == 0 || y == 0 || z == 0 || x == NX - 1 || y == NY - 1 || z == NZ - 1) && solid[idx(x, y, z)] && !seen[idx(x, y, z)]) flood(x, y, z);
    const long long cell = (long long)step * step * step;
    for (int y = 0; y < NY; ++y)
        for (int z = 0; z < NZ; ++z)
            for (int x = 0; x < NX; ++x)
                if (solid[idx(x, y, z)] && !seen[idx(x, y, z)]) {
                    sumY = sumR = 0.0;
                    const long long n = flood(x, y, z), v = n * cell;
                    ++out.pieces; out.volume += v; out.largest = std::max(out.largest, v);
                    const double my = sumY / n, mr = sumR / n; // bins: height -4..+4 binH in 8, distance 0..8 binR in 8
                    if (hb) hb[std::clamp((int)((my + 4.0 * binH) / binH), 0, 7)] += v;
                    if (db) db[std::clamp((int)(mr / binR), 0, 7)] += v;
                }
    return out;
}
} // namespace

// (moved from main.cpp's --cave-islands, #1465: the same measurement, now also a scenario's; caverns added)
CaveIslands CaveIslands_Measure(const PlanetField& pf, unsigned seed, float R, bool print, float cavernRange) {
    CaveIslands out;
    const Vector3 home = pf.findSpawnDir(seed);
    std::vector<Vector3> chambers;
    pf.caveMouthsList(home, R, 3, chambers);
    out.chambers = (int)chambers.size();
    const auto t0 = std::chrono::steady_clock::now();
    for (size_t ci = 0; ci < chambers.size(); ++ci) {
        const BoxResult b = MeasureBox(pf, chambers[ci], 50, 28, 1, out.byHeight, out.byDistance);
        out.pieces += b.pieces; out.volume += b.volume;
        if (print) {
            const Vector3 up = Vector3Normalize(chambers[ci]);
            const float depth = pf.radius + pf.elevation(up) - Vector3Length(chambers[ci]);
            printf("bench: cave islands: chamber %zu (%.0f m down): %lld floating pieces, %lld m^3, largest %lld m^3\n", ci, depth, b.pieces, b.volume, b.largest);
        }
    }
    // (#1465) Large caverns (centre kind 4) within cavernRange: a box over the whole cavern (1.9 r across, the dome and floor
    // with 12 m to spare), on a 2 m grid (pieces under ~8 m^3 can slip through it).
    if (cavernRange > 0.0f) {
        std::vector<Vector3> caverns;
        std::vector<float> radii;
        pf.caveMouthsList(home, cavernRange, 4, caverns, &radii);
        out.caverns = (int)caverns.size();
        for (size_t ci = 0; ci < caverns.size(); ++ci) {
            const float r = radii[ci];
            const int hx = 2 * (int)ceilf(0.95f * r), hy = 2 * (int)ceilf(0.5f * (0.4f * r * 1.82f + 12.0f));
            static const int step = getenv("RAYLIBGAME_ISLANDS_CAVERN_STEP") ? std::max(1, atoi(getenv("RAYLIBGAME_ISLANDS_CAVERN_STEP"))) : 2; // (1: exact, slow)
            long long hb[8] = {}, db[8] = {};
            const BoxResult b = MeasureBox(pf, caverns[ci], hx, hy, step, hb, db, 0.1 * r, 0.125 * r);
            out.cavernPieces += b.pieces; out.cavernVolume += b.volume;
            if (print) {
                const Vector3 up = Vector3Normalize(caverns[ci]);
                const float depth = pf.radius + pf.elevation(up) - Vector3Length(caverns[ci]);
                const float away = Vector3Distance(Vector3Scale(up, pf.radius), Vector3Scale(home, pf.radius));
                printf("bench: cave islands: cavern %zu (r %.0f m, %.0f m down, %.0f m away): %lld floating pieces, %lld m^3, largest %lld m^3\n", ci, r, depth, away, b.pieces, b.volume, b.largest);
                printf("bench: cave islands: cavern %zu m^3 by height (x0.1 r from the centre, -0.4..+0.4 r):", ci); for (int k = 0; k < 8; ++k) printf(" %lld", hb[k]); printf("; by distance across (x0.125 r, 0..r):"); for (int k = 0; k < 8; ++k) printf(" %lld", db[k]); printf("\n");
            }
        }
        // (#1465) Sheet caves (kind 5, listed at their middles; radius = half their length + their half-width) within the same
        // range: a box over the whole sheet (1.1 x that across, +-30 m up: low, but it follows the ground's slope), 2 m grid.
        std::vector<Vector3> sheets;
        pf.caveMouthsList(home, cavernRange, 5, sheets, &radii);
        out.sheets = (int)sheets.size();
        for (size_t si = 0; si < sheets.size(); ++si) {
            const int hx = 2 * (int)ceilf(0.5f * (1.1f * radii[si] + 6.0f)), hy = 30;
            const BoxResult b = MeasureBox(pf, sheets[si], hx, hy, 2, nullptr, nullptr);
            out.sheetPieces += b.pieces; out.sheetVolume += b.volume;
            if (print) printf("bench: cave islands: sheet %zu (half-extent %.0f m): %lld floating pieces, %lld m^3, largest %lld m^3\n", si, radii[si], b.pieces, b.volume, b.largest);
        }
    }
    out.seconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
    if (print) {
        printf("bench: cave islands: %d chambers within %.0f m: %lld floating pieces, %lld m^3 in all (%.1f s)\n", out.chambers, R, out.pieces, out.volume, out.seconds);
        printf("bench: cave islands: floating m^3 by height from the centre (m):"); for (int b = 0; b < 8; ++b) printf(" %+d..%+d:%lld", -28 + 7 * b, -21 + 7 * b, out.byHeight[b]); printf("\n");
        printf("bench: cave islands: floating m^3 by distance across from the centre (m):"); for (int b = 0; b < 8; ++b) printf(" %d..%d:%lld", 9 * b, 9 * b + 9, out.byDistance[b]); printf("\n");
        if (cavernRange > 0.0f) printf("bench: cave islands: %d caverns within %.0f m: %lld floating pieces, %lld m^3 (2 m grid)\n", out.caverns, cavernRange, out.cavernPieces, out.cavernVolume);
        if (cavernRange > 0.0f) printf("bench: cave islands: %d sheet caves within %.0f m: %lld floating pieces, %lld m^3 (2 m grid)\n", out.sheets, cavernRange, out.sheetPieces, out.sheetVolume);
    }
    return out;
}
