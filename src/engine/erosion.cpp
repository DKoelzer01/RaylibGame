#include "erosion.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <thread>

namespace {
    // The cube faces as PlanetField::prepareBlock lays them out (#251): face = axis * 2 + (positive ? 1 : 0).
    void FaceAxes(int f, Vector3& nrm, Vector3& t1, Vector3& t2) {
        const float s = (f & 1) ? 1.0f : -1.0f;
        switch (f >> 1) {
            case 0: nrm = { s, 0, 0 }; t1 = { 0, 1, 0 }; t2 = { 0, 0, 1 }; break;
            case 1: nrm = { 0, s, 0 }; t1 = { 0, 0, 1 }; t2 = { 1, 0, 0 }; break;
            default: nrm = { 0, 0, s }; t1 = { 1, 0, 0 }; t2 = { 0, 1, 0 }; break;
        }
    }
    int FaceOf(Vector3 d) {
        const float ax = fabsf(d.x), ay = fabsf(d.y), az = fabsf(d.z);
        if (ax >= ay && ax >= az) return d.x >= 0.0f ? 1 : 0;
        if (ay >= az) return d.y >= 0.0f ? 3 : 2;
        return d.z >= 0.0f ? 5 : 4;
    }
    float Hash3(int x, int y, int z, unsigned int seed) {
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
        h ^= h >> 13; h *= 0x5bd1e995u; h ^= h >> 15;
        return (h & 0xFFFFFF) / 16777215.0f;
    }
    // Low-frequency value noise on the sphere (rock regions): smooth, deterministic, seeded.
    float ValueNoise(Vector3 p, unsigned int seed) {
        const int x0 = (int)floorf(p.x), y0 = (int)floorf(p.y), z0 = (int)floorf(p.z);
        const float fx = p.x - x0, fy = p.y - y0, fz = p.z - z0;
        auto s = [](float t) { return t * t * (3.0f - 2.0f * t); };
        const float u = s(fx), v = s(fy), w = s(fz);
        float c[2][2][2];
        for (int k = 0; k < 2; ++k) for (int j = 0; j < 2; ++j) for (int i = 0; i < 2; ++i) c[k][j][i] = Hash3(x0 + i, y0 + j, z0 + k, seed);
        const float a = (c[0][0][0] * (1 - u) + c[0][0][1] * u) * (1 - v) + (c[0][1][0] * (1 - u) + c[0][1][1] * u) * v;
        const float b = (c[1][0][0] * (1 - u) + c[1][0][1] * u) * (1 - v) + (c[1][1][0] * (1 - u) + c[1][1][1] * u) * v;
        return a * (1 - w) + b * w;
    }
    int16_t Q(float h) { return (int16_t)std::clamp((int)lroundf(h * 10.0f), -32767, 32767); }

    // One simulation domain: cells with directions, areas, outlets and rock.
    struct Domain {
        size_t cells = 0;
        std::vector<Vector3> dir;
        std::vector<float> cellArea, h, U, rockK, hardness, strata;
        std::vector<uint8_t> outlet;
    };

    // Rock and uplift of a domain from its input heights (h0) and the planet's highest input (maxH0).
    void SetRock(Domain& d, const std::vector<float>& h0, float maxH0, float R, unsigned int seed, const ErosionParams& p) {
        d.U.resize(d.cells); d.rockK.resize(d.cells); d.hardness.resize(d.cells); d.strata.resize(d.cells);
        for (size_t c = 0; c < d.cells; ++c) {
            const Vector3 q = Vector3Scale(d.dir[c], R / 25000.0f);
            d.hardness[c] = ValueNoise(q, seed ^ 0xA11u);
            d.strata[c] = std::clamp((ValueNoise(Vector3Scale(q, 0.7f), seed ^ 0x57Au) - 0.55f) / 0.15f, 0.0f, 1.0f);
            d.rockK[c] = p.k * (1.6f - 1.2f * d.hardness[c]);
            // A floor of 5% for all land: uplift only in proportion to the input left wide low plains dead flat, with
            // straight 8-direction river traces across them (#355, seed 7).
            d.U[c] = h0[c] > 0.0f ? p.uplift * (0.05f + 0.95f * h0[c] / maxH0) : 0.0f;
        }
    }

    // Stream power + talus for `steps` steps (see erosion.h). nbOf(c, k) = the k-th of 8 neighbours or -1.
    template <class NB>
    void Simulate(Domain& d, NB&& nbOf, int steps, const ErosionParams& p, float R, unsigned int seed, float kScale) {
        const size_t cells = d.cells;
        const float stratumOffset = Hash3(1, 2, 3, seed) * 120.0f;
        auto capRock = [&](float hc) { return fmodf(hc + stratumOffset + 1200.0f, 120.0f) < 36.0f; }; // the hard 30% of each 120 m band
        auto erodibility = [&](size_t c, float hc) {
            const float k = d.rockK[c] * kScale;
            if (d.strata[c] <= 0.0f) return k;
            const float bandK = capRock(hc) ? 0.15f : 1.6f;
            return k * (1.0f + (bandK - 1.0f) * d.strata[c]);
        };
        std::vector<float> hf(cells), area(cells);
        std::vector<int> receiver(cells), order;
        std::vector<uint8_t> seen(cells);
        order.reserve(cells);
        using Item = std::pair<float, int>;
        int step = 0;
        auto flood = [&] {
            std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
            std::fill(seen.begin(), seen.end(), 0);
            order.clear();
            for (size_t c = 0; c < cells; ++c)
                if (d.outlet[c]) { hf[c] = d.h[c]; receiver[c] = (int)c; seen[c] = 1; pq.push({ d.h[c], (int)c }); }
            if (pq.empty()) { // no sea: the lowest cell drains everything
                size_t lo = 0; for (size_t c = 1; c < cells; ++c) if (d.h[c] < d.h[lo]) lo = c;
                hf[lo] = d.h[lo]; receiver[lo] = (int)lo; seen[lo] = 1; pq.push({ d.h[lo], (int)lo });
            }
            while (!pq.empty()) {
                const int c = pq.top().second; pq.pop();
                order.push_back(c);
                for (int k = 0; k < 8; ++k) {
                    const int o = nbOf((size_t)c, k);
                    if (o < 0 || seen[o]) continue;
                    seen[o] = 1;
                    hf[o] = std::max(d.h[o], hf[c] + 1e-3f * (1.0f + 4.0f * Hash3(o, step, 7, seed))); // jittered across flats
                    receiver[o] = c;
                    pq.push({ hf[o], o });
                }
            }
            for (int c : order) { // steepest descent on the filled surface (the flood parent where none is lower)
                if (receiver[c] == c) continue;
                float best = 0.0f;
                for (int k = 0; k < 8; ++k) {
                    const int o = nbOf((size_t)c, k);
                    if (o < 0) continue;
                    const float s = (hf[c] - hf[o]) / (Vector3Distance(d.dir[c], d.dir[o]) * R);
                    if (s > best) { best = s; receiver[c] = o; }
                }
            }
            for (size_t c = 0; c < cells; ++c) area[c] = d.cellArea[c];
            for (size_t k = order.size(); k-- > 0;) { const int c = order[k]; if (receiver[c] != c) area[receiver[c]] += area[c]; }
        };
        for (step = 0; step < steps; ++step) {
            if (step % std::max(1, p.floodEvery) == 0) flood();
            for (int c : order) { // implicit stream power, receivers first
                const int r = receiver[c];
                if (r == c || d.outlet[c]) continue;
                if (hf[c] > d.h[c] + 1e-3f) { d.h[c] += d.U[c]; continue; } // a filled depression: a lake bed, no incision
                const float L = Vector3Distance(d.dir[c], d.dir[r]) * R;
                const float f = erodibility(c, d.h[c]) * powf(area[c], p.m) / L;
                d.h[c] = std::max((d.h[c] + d.U[c] + f * d.h[r]) / (1.0f + f), d.h[r]);
            }
            if (step % 2 == 1) // talus: slopes steeper than the rock holds slump (half the excess moves down)
                for (size_t c = 0; c < cells; ++c) {
                    if (d.outlet[c]) continue;
                    const float hard = std::max(d.hardness[c], d.strata[c] * (capRock(d.h[c]) ? 1.0f : 0.0f));
                    const float maxSlope = p.talusSoft + (p.talusHard - p.talusSoft) * hard;
                    for (int k = 0; k < 8; ++k) {
                        const int o = nbOf(c, k);
                        if (o < 0) continue;
                        const float excess = d.h[c] - d.h[o] - maxSlope * Vector3Distance(d.dir[c], d.dir[o]) * R;
                        if (excess > 0.0f) { const float mv = 0.5f * excess; d.h[c] -= mv; if (!d.outlet[o]) d.h[o] += mv; }
                    }
                }
        }
    }

    // Bilinear on one face of a 6 x side^2 grid (values at cell centres); u, v in face coordinates, the side cells
    // spanning [-span, span] (span 1 for a plain face).
    float FaceBilinear(const std::vector<int16_t>& g, int f, int side, float u, float v, float span) {
        const float du = 2.0f * span / side;
        const float fi = std::clamp((u + span) / du - 0.5f, 0.0f, (float)(side - 1) - 1e-3f);
        const float fj = std::clamp((v + span) / du - 0.5f, 0.0f, (float)(side - 1) - 1e-3f);
        const int i = (int)fi, j = (int)fj;
        const float a = fi - i, b = fj - j;
        const int16_t* r0 = &g[((size_t)f * side + j) * side + i];
        const int16_t* r1 = r0 + side;
        return ((r0[0] * (1 - a) + r0[1] * a) * (1 - b) + (r1[0] * (1 - a) + r1[1] * a) * b) * 0.1f;
    }

    // Hillshaded map of one face of a side^2-per-face grid (#355 evidence).
    void DumpFace(const char* path, const float* h, int N, float cell) {
        Image img = GenImageColor(N, N, BLACK);
        auto H = [&](int i, int j) { i = std::clamp(i, 0, N - 1); j = std::clamp(j, 0, N - 1); return h[(size_t)j * N + i]; };
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
                const float gx = (H(i + 1, j) - H(i - 1, j)) / (2 * cell), gy = (H(i, j + 1) - H(i, j - 1)) / (2 * cell);
                const Vector3 nrm = Vector3Normalize({ -gx, 1.0f, -gy });
                const float shade = std::clamp(0.35f + 0.75f * Vector3DotProduct(nrm, Vector3Normalize({ -0.6f, 0.7f, -0.4f })), 0.0f, 1.0f);
                const float e = H(i, j);
                const Vector3 c = e <= 0.0f ? Vector3{ 0.15f, 0.3f, 0.6f } : e < 200.0f ? Vector3{ 0.35f, 0.6f, 0.3f } : e < 700.0f ? Vector3{ 0.5f, 0.5f, 0.35f } : Vector3{ 0.85f, 0.85f, 0.9f };
                ImageDrawPixel(&img, i, j, { (unsigned char)(255 * c.x * shade), (unsigned char)(255 * c.y * shade), (unsigned char)(255 * c.z * shade), 255 });
            }
        ExportImage(img, path);
        UnloadImage(img);
    }
    bool DumpOn() { const char* e = getenv("RAYLIBGAME_EROSION_DUMP"); return e && e[0] == '1'; }
}

void PlanetErosion::build(const std::function<float(Vector3, float)>& height, float planetRadius, unsigned int seedIn, const ErosionParams& p, int threads) {
    const auto t0 = std::chrono::steady_clock::now();
    R = planetRadius; seed = seedIn;
    const int T = std::max(1, threads);
    static const int di[8] = { 1, -1, 0, 0, 1, 1, -1, -1 }, dj[8] = { 0, 0, 1, -1, 1, -1, 1, -1 };
    // ---- The planet, coarse ----
    const int N = std::max(16, (int)lroundf(R * PI * 0.5f / p.cell));
    Domain d;
    d.cells = (size_t)6 * N * N;
    auto idx = [N](int f, int i, int j) { return ((size_t)f * N + j) * N + i; };
    d.dir.resize(d.cells); d.cellArea.resize(d.cells);
    const float du = 2.0f / N, side = R * du;
    for (int f = 0; f < 6; ++f) {
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
                const float u = -1.0f + (i + 0.5f) * du, v = -1.0f + (j + 0.5f) * du;
                d.dir[idx(f, i, j)] = Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v))));
                d.cellArea[idx(f, i, j)] = side * side / powf(1.0f + u * u + v * v, 1.5f);
            }
    }
    auto cellOf = [&](Vector3 dv) -> int {
        const int f = FaceOf(dv);
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        const float s = Vector3DotProduct(dv, nrm);
        const int i = std::clamp((int)floorf((Vector3DotProduct(dv, t1) / s + 1.0f) / du), 0, N - 1);
        const int j = std::clamp((int)floorf((Vector3DotProduct(dv, t2) / s + 1.0f) / du), 0, N - 1);
        return (int)idx(f, i, j);
    };
    std::vector<int> nb(d.cells * 8); // 8 neighbours per cell, across face edges by direction
    for (int f = 0; f < 6; ++f) {
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
                const size_t c = idx(f, i, j);
                for (int k = 0; k < 8; ++k) {
                    const int ii = i + di[k], jj = j + dj[k];
                    if (ii >= 0 && jj >= 0 && ii < N && jj < N) { nb[c * 8 + k] = (int)idx(f, ii, jj); continue; }
                    const float u = -1.0f + (ii + 0.5f) * du, v = -1.0f + (jj + 0.5f) * du;
                    const int o = cellOf(Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v)))));
                    nb[c * 8 + k] = o == (int)c ? -1 : o;
                }
            }
    }
    std::vector<float> h0(d.cells);
    {
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t) pool.emplace_back([&, t] { for (size_t c = t; c < d.cells; c += T) h0[c] = height(d.dir[c], p.cell); });
        for (auto& th : pool) th.join();
    }
    float maxH0 = 1.0f;
    for (float v : h0) maxH0 = std::max(maxH0, v);
    SetRock(d, h0, maxH0, R, seed, p);
    d.outlet.resize(d.cells); d.h.resize(d.cells);
    for (size_t c = 0; c < d.cells; ++c) { d.outlet[c] = h0[c] <= 0.0f; d.h[c] = h0[c] > 0.0f ? h0[c] * p.start : h0[c]; } // start low (uplift mode) or from the input (carve mode)
    Simulate(d, [&](size_t c, int k) { return nb[c * 8 + k]; }, p.steps, p, R, seed, 1.0f);
    // Back to the input's height range: the 99th percentile of land as high as the input's. (Matching only the highest
    // peak left the ranges ~30% lower than the input: p90 1041 -> 723 m on seed 42.)
    float scale = 1.0f;
    {
        std::vector<float> a, b;
        for (size_t c = 0; c < d.cells; ++c) if (!d.outlet[c]) { a.push_back(h0[c]); b.push_back(d.h[c]); }
        if (!a.empty()) {
            const size_t k = (size_t)(0.99 * (a.size() - 1));
            std::nth_element(a.begin(), a.begin() + k, a.end()); std::nth_element(b.begin(), b.begin() + k, b.end());
            scale = std::max(1.0f, a[k]) / std::max(1.0f, b[k]);
        }
    }
    n = N;
    coarse.resize(d.cells);
    maxCut = maxRaise = 0.0f;
    for (size_t c = 0; c < d.cells; ++c) {
        const float hb = d.outlet[c] ? h0[c] : std::max(0.5f, d.h[c] * scale);
        coarse[c] = Q(hb);
        maxCut = std::max(maxCut, h0[c] - hb); maxRaise = std::max(maxRaise, hb - h0[c]);
    }
    int dumpFace = 0;
    if (DumpOn()) {
        std::vector<float> a, b; // land height percentiles, input vs baked
        for (size_t c = 0; c < d.cells; ++c) if (h0[c] > 0.0f) { a.push_back(h0[c]); b.push_back(coarse[c] * 0.1f); }
        std::sort(a.begin(), a.end()); std::sort(b.begin(), b.end());
        if (!a.empty()) {
            printf("erosion: land height percentiles (input -> baked):");
            for (float q : { 0.5f, 0.75f, 0.9f, 0.97f, 0.99f, 0.999f }) printf(" p%.1f %.0f->%.0f", q * 100, a[(size_t)(q * (a.size() - 1))], b[(size_t)(q * (b.size() - 1))]);
            printf("\n");
        }
        int bestLand = -1;
        for (int fc = 0; fc < 6; ++fc) { int land = 0; for (int k = 0; k < N * N; ++k) land += h0[(size_t)fc * N * N + k] > 0.0f; if (land > bestLand) { bestLand = land; dumpFace = fc; } }
        std::vector<float> after((size_t)N * N);
        for (size_t k = 0; k < after.size(); ++k) after[k] = coarse[(size_t)dumpFace * N * N + k] * 0.1f;
        DumpFace("erosion_before.png", &h0[(size_t)dumpFace * N * N], N, p.cell);
        DumpFace("erosion_after.png", after.data(), N, p.cell);
    }
    std::vector<float> coarseUnscaled = std::move(d.h); // the refinement continues in the simulation's own units
    coarseSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
    d = Domain(); nb = std::vector<int>(); h0 = std::vector<float>(); // free the coarse working set before the faces

    // ---- Each face at `refine` x the resolution, with a margin past its edges ----
    const auto t1r = std::chrono::steady_clock::now();
    nf = fm = fe = 0; fine.clear(); fineMetres = 0.0f;
    if (p.refine >= 2) {
        nf = N * p.refine; fm = p.margin; fe = nf + 2 * fm;
        fineMetres = p.cell / p.refine;
        fine.assign((size_t)6 * fe * fe, 0);
        const int FE = fe;
        const float duf = 2.0f / nf, span = 1.0f + fm * duf, sideF = R * duf;
        auto coarseUnscaledAt = [&](Vector3 dv) {
            const int f = FaceOf(dv);
            Vector3 nrm, a1, a2; FaceAxes(f, nrm, a1, a2);
            const float s = Vector3DotProduct(dv, nrm), u = Vector3DotProduct(dv, a1) / s, v = Vector3DotProduct(dv, a2) / s;
            const float fi = std::clamp((u + 1.0f) / du - 0.5f, 0.0f, (float)(N - 1) - 1e-3f), fj = std::clamp((v + 1.0f) / du - 0.5f, 0.0f, (float)(N - 1) - 1e-3f);
            const int i = (int)fi, j = (int)fj; const float a = fi - i, b = fj - j;
            const float* r0 = &coarseUnscaled[((size_t)f * N + j) * N + i];
            const float* r1 = r0 + N;
            return (r0[0] * (1 - a) + r0[1] * a) * (1 - b) + (r1[0] * (1 - a) + r1[1] * a) * b;
        };
        // Finer cells drain smaller areas: K is scaled by refine^(2m) so the coarse steady state holds, and the finer
        // valleys appear on top of it.
        const float kScale = powf((float)(p.refine * p.refine), p.m);
        auto refineFace = [&, FE](int f, int sampleThreads) {
            Vector3 nrm, a1, a2; FaceAxes(f, nrm, a1, a2);
            Domain fd;
            fd.cells = (size_t)FE * FE;
            fd.dir.resize(fd.cells); fd.cellArea.resize(fd.cells);
            for (int j = 0; j < FE; ++j)
                for (int i = 0; i < FE; ++i) {
                    const float u = -span + (i + 0.5f) * duf, v = -span + (j + 0.5f) * duf;
                    fd.dir[(size_t)j * FE + i] = Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(a1, u), Vector3Scale(a2, v))));
                    fd.cellArea[(size_t)j * FE + i] = sideF * sideF / powf(1.0f + u * u + v * v, 1.5f);
                }
            std::vector<float> fh0(fd.cells);
            {
                std::vector<std::thread> pool;
                for (int t = 0; t < sampleThreads; ++t)
                    pool.emplace_back([&, t] { for (size_t c = t; c < fd.cells; c += sampleThreads) fh0[c] = height(fd.dir[c], fineMetres); });
                for (auto& th : pool) th.join();
            }
            SetRock(fd, fh0, maxH0, R, seed, p);
            fd.outlet.resize(fd.cells); fd.h.resize(fd.cells);
            for (int j = 0; j < FE; ++j)
                for (int i = 0; i < FE; ++i) {
                    const size_t c = (size_t)j * FE + i;
                    const bool rim = i == 0 || j == 0 || i == FE - 1 || j == FE - 1; // the margin's edge: held at the coarse result
                    fd.outlet[c] = fh0[c] <= 0.0f || rim;
                    fd.h[c] = fh0[c] <= 0.0f ? fh0[c] : coarseUnscaledAt(fd.dir[c]);
                }
            Simulate(fd, [FE](size_t c, int k) {
                const int i = (int)(c % FE) + di[k], j = (int)(c / FE) + dj[k];
                return (i < 0 || j < 0 || i >= FE || j >= FE) ? -1 : j * FE + i;
            }, p.refineSteps, p, R, seed, kScale);
            for (size_t c = 0; c < fd.cells; ++c)
                fine[(size_t)f * FE * FE + c] = Q(fh0[c] <= 0.0f ? fh0[c] : std::max(0.5f, fd.h[c] * scale));
            if (DumpOn() && f == dumpFace) {
                std::vector<float> out(fd.cells);
                for (size_t c = 0; c < fd.cells; ++c) out[c] = fine[(size_t)f * FE * FE + c] * 0.1f;
                DumpFace("erosion_fine.png", out.data(), FE, fineMetres);
            }
        };
        const int parallelFaces = std::clamp(T / 4, 1, 3); // a face's working set is ~200 MB
        for (int f0 = 0; f0 < 6; f0 += parallelFaces) {
            std::vector<std::thread> pool;
            for (int f = f0; f < std::min(6, f0 + parallelFaces); ++f) pool.emplace_back(refineFace, f, std::max(1, T / parallelFaces));
            for (auto& th : pool) th.join();
        }
    }
    refineSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t1r).count();
    computeHash();
    buildSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
}

float PlanetErosion::coarseAt(Vector3 dv) const {
    const int f = FaceOf(dv);
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float s = Vector3DotProduct(dv, nrm);
    return FaceBilinear(coarse, f, n, Vector3DotProduct(dv, t1) / s, Vector3DotProduct(dv, t2) / s, 1.0f);
}

float PlanetErosion::fineAt(Vector3 dv) const {
    // Every face whose extended grid holds dv, weighted down across its margin: the faces agree at their edges.
    const float duf = 2.0f / nf, span = 1.0f + fm * duf, inner = span - duf, band = std::max(duf, fm * duf - duf);
    {   // Fast path: the home face alone (weight 1 inside |u|, |v| <= 1) when no other face's grid reaches dv. A
        // neighbour face sees a point at home coordinate e (near its edge) at ~1/e, so it counts only from e > 1/inner.
        const int f = FaceOf(dv);
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        const float s = Vector3DotProduct(dv, nrm);
        const float u = Vector3DotProduct(dv, t1) / s, v = Vector3DotProduct(dv, t2) / s;
        if (std::max(fabsf(u), fabsf(v)) * inner < 0.999f) return FaceBilinear(fine, f, fe, u, v, span);
    }
    float sum = 0.0f, wsum = 0.0f;
    for (int f = 0; f < 6; ++f) {
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        const float s = Vector3DotProduct(dv, nrm);
        if (s <= 0.0f) continue;
        const float u = Vector3DotProduct(dv, t1) / s, v = Vector3DotProduct(dv, t2) / s;
        const float edge = std::max(fabsf(u), fabsf(v));
        if (edge >= inner) continue;
        float w = std::clamp((inner - edge) / band, 0.0f, 1.0f);
        w = w * w * (3.0f - 2.0f * w);
        if (w <= 0.0f) continue;
        sum += w * FaceBilinear(fine, f, fe, u, v, span);
        wsum += w;
    }
    return wsum > 0.0f ? sum / wsum : coarseAt(dv);
}

float PlanetErosion::heightAt(Vector3 dv, float footprint) const {
    if (n <= 0) return 0.0f;
    if (fe <= 0 || footprint >= 2.0f * fineMetres) return coarseAt(dv);
    return fineAt(dv);
}

void PlanetErosion::computeHash() {
    uint64_t hv = 1469598103934665603ull; // FNV-1a over the stored grids
    for (int16_t q : coarse) hv = (hv ^ (uint16_t)q) * 1099511628211ull;
    for (int16_t q : fine) hv = (hv ^ (uint16_t)q) * 1099511628211ull;
    hashValue = hv;
}

bool PlanetErosion::save(const std::string& path) const {
    if (n <= 0) return false;
    // Each value minus the one before it (smooth heights: small numbers), low bytes then high bytes: deflates ~5x better.
    std::vector<int16_t> all(coarse);
    all.insert(all.end(), fine.begin(), fine.end());
    std::vector<unsigned char> raw(all.size() * 2);
    int16_t prev = 0;
    for (size_t k = 0; k < all.size(); ++k) {
        const uint16_t dq = (uint16_t)(all[k] - prev); prev = all[k];
        raw[k] = (unsigned char)(dq & 0xFF); raw[all.size() + k] = (unsigned char)(dq >> 8);
    }
    int compSize = 0;
    unsigned char* comp = CompressData(raw.data(), (int)raw.size(), &compSize);
    if (!comp) return false;
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) { MemFree(comp); return false; }
    const char magic[4] = { 'E', 'R', 'O', '3' };
    fwrite(magic, 1, 4, f);
    fwrite(&n, sizeof(n), 1, f); fwrite(&nf, sizeof(nf), 1, f); fwrite(&fm, sizeof(fm), 1, f); fwrite(&fineMetres, sizeof(fineMetres), 1, f);
    fwrite(&R, sizeof(R), 1, f); fwrite(&seed, sizeof(seed), 1, f); fwrite(&hashValue, sizeof(hashValue), 1, f);
    fwrite(&compSize, sizeof(compSize), 1, f);
    fwrite(comp, 1, compSize, f);
    fclose(f);
    MemFree(comp);
    return true;
}

bool PlanetErosion::load(const std::string& path, float planetRadius, unsigned int wantSeed) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    char magic[4]; int N = 0, NF = 0, FM = 0, compSize = 0; float fmet = 0, r = 0; unsigned int s = 0; uint64_t hv = 0;
    bool ok = fread(magic, 1, 4, f) == 4 && memcmp(magic, "ERO3", 4) == 0 && fread(&N, sizeof(N), 1, f) == 1 && fread(&NF, sizeof(NF), 1, f) == 1 &&
              fread(&FM, sizeof(FM), 1, f) == 1 && fread(&fmet, sizeof(fmet), 1, f) == 1 && fread(&r, sizeof(r), 1, f) == 1 && fread(&s, sizeof(s), 1, f) == 1 &&
              fread(&hv, sizeof(hv), 1, f) == 1 && fread(&compSize, sizeof(compSize), 1, f) == 1 &&
              N > 0 && N < 4096 && NF >= 0 && NF < 16384 && FM >= 0 && FM < 1024 && compSize > 0 && s == wantSeed && fabsf(r - planetRadius) < 1.0f;
    std::vector<unsigned char> comp;
    if (ok) { comp.resize(compSize); ok = fread(comp.data(), 1, compSize, f) == (size_t)compSize; }
    fclose(f);
    if (!ok) return false;
    int rawSize = 0;
    unsigned char* raw = DecompressData(comp.data(), compSize, &rawSize);
    if (!raw) return false;
    const size_t nc = (size_t)6 * N * N, fe2 = NF > 0 ? (size_t)NF + 2 * FM : 0, nfine = (size_t)6 * fe2 * fe2;
    if ((size_t)rawSize != (nc + nfine) * 2) { MemFree(raw); return false; }
    coarse.resize(nc); fine.resize(nfine);
    const size_t total = nc + nfine;
    int16_t prev = 0;
    for (size_t k = 0; k < total; ++k) {
        prev = (int16_t)(prev + (int16_t)(raw[k] | (raw[total + k] << 8)));
        (k < nc ? coarse[k] : fine[k - nc]) = prev;
    }
    MemFree(raw);
    n = N; nf = NF; fm = FM; fe = (int)fe2; fineMetres = fmet; R = r; seed = s;
    computeHash();
    if (hashValue != hv) { n = 0; coarse.clear(); fine.clear(); fe = 0; return false; } // damaged
    return true;
}
