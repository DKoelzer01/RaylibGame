#include "caveislands.h"
#include "voxel.h"

#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <vector>

// (moved from main.cpp's --cave-islands, #1465: the same measurement, now also a scenario's)
CaveIslands CaveIslands_Measure(const PlanetField& pf, unsigned seed, float R, bool print) {
    CaveIslands out;
    const Vector3 home = pf.findSpawnDir(seed);
    std::vector<Vector3> chambers;
    pf.caveMouthsList(home, R, 3, chambers);
    out.chambers = (int)chambers.size();
    const int NX = 101, NY = 57, NZ = 101;                          // x, z across; y up
    const auto t0 = std::chrono::steady_clock::now();
    for (size_t ci = 0; ci < chambers.size(); ++ci) {
        const Vector3 c = chambers[ci], up = Vector3Normalize(c);
        const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })), e2 = Vector3CrossProduct(up, e1);
        std::vector<uint8_t> solid((size_t)NX * NY * NZ), seen((size_t)NX * NY * NZ, 0);
        auto idx = [&](int x, int y, int z) { return ((size_t)y * NZ + z) * NX + x; };
        for (int y = 0; y < NY; ++y)
            for (int z = 0; z < NZ; ++z)
                for (int x = 0; x < NX; ++x) {
                    const Vector3 q = Vector3Add(c, Vector3Add(Vector3Scale(e1, x - 50.0f), Vector3Add(Vector3Scale(up, y - 28.0f), Vector3Scale(e2, z - 50.0f))));
                    solid[idx(x, y, z)] = pf.baseDensity(q) >= ISO_LEVEL;
                }
        std::vector<int> stack;
        double sumY = 0.0, sumR = 0.0; // the last flooded component's summed height (m from the centre) and distance across
        auto flood = [&](int sx, int sy, int sz) {                  // marks one solid component, returns its size
            long long n = 0;
            stack.clear(); stack.push_back((int)idx(sx, sy, sz)); seen[idx(sx, sy, sz)] = 1;
            while (!stack.empty()) {
                const int k = stack.back(); stack.pop_back(); ++n;
                const int x = k % NX, z = (k / NX) % NZ, y = k / (NX * NZ);
                sumY += y - 28.0; sumR += sqrt((x - 50.0) * (x - 50.0) + (z - 50.0) * (z - 50.0));
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
        long long islands = 0, vol = 0, largest = 0;
        for (int y = 0; y < NY; ++y)
            for (int z = 0; z < NZ; ++z)
                for (int x = 0; x < NX; ++x)
                    if (solid[idx(x, y, z)] && !seen[idx(x, y, z)]) {
                        sumY = sumR = 0.0;
                        const long long n = flood(x, y, z); ++islands; vol += n; largest = std::max(largest, n);
                        const double my = sumY / n, mr = sumR / n; // bins: height -28..+28 m in 8, distance 0..71 m in 8
                        out.byHeight[std::clamp((int)((my + 28.0) / 7.0), 0, 7)] += n; out.byDistance[std::clamp((int)(mr / 9.0), 0, 7)] += n;
                    }
        out.pieces += islands; out.volume += vol;
        if (print) {
            const float depth = pf.radius + pf.elevation(up) - Vector3Length(c);
            printf("bench: cave islands: chamber %zu (%.0f m down): %lld floating pieces, %lld m^3, largest %lld m^3\n", ci, depth, islands, vol, largest);
        }
    }
    out.seconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
    if (print) {
        printf("bench: cave islands: %d chambers within %.0f m: %lld floating pieces, %lld m^3 in all (%.1f s)\n", out.chambers, R, out.pieces, out.volume, out.seconds);
        printf("bench: cave islands: floating m^3 by height from the centre (m):"); for (int b = 0; b < 8; ++b) printf(" %+d..%+d:%lld", -28 + 7 * b, -21 + 7 * b, out.byHeight[b]); printf("\n");
        printf("bench: cave islands: floating m^3 by distance across from the centre (m):"); for (int b = 0; b < 8; ++b) printf(" %d..%d:%lld", 9 * b, 9 * b + 9, out.byDistance[b]); printf("\n");
    }
    return out;
}
