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
    // Low-frequency value noise on the sphere (rock regions): smooth, deterministic, seeded.
    float Hash3(int x, int y, int z, unsigned int seed) {
        uint32_t h = (uint32_t)x * 73856093u ^ (uint32_t)y * 19349663u ^ (uint32_t)z * 83492791u ^ seed * 2654435761u;
        h ^= h >> 13; h *= 0x5bd1e995u; h ^= h >> 15;
        return (h & 0xFFFFFF) / 16777215.0f;
    }
    // Hillshaded maps of one face (#355 evidence): the face with the most land, before and after, drainage in blue.
    void DumpFace(const char* path, const std::vector<float>& h, const std::vector<float>* area, int f, int N, float cell) {
        Image img = GenImageColor(N, N, BLACK);
        auto H = [&](int i, int j) { i = std::clamp(i, 0, N - 1); j = std::clamp(j, 0, N - 1); return h[((size_t)f * N + j) * N + i]; };
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
                const float gx = (H(i + 1, j) - H(i - 1, j)) / (2 * cell), gy = (H(i, j + 1) - H(i, j - 1)) / (2 * cell);
                const Vector3 nrm = Vector3Normalize({ -gx, 1.0f, -gy });
                const float shade = std::clamp(0.35f + 0.75f * Vector3DotProduct(nrm, Vector3Normalize({ -0.6f, 0.7f, -0.4f })), 0.0f, 1.0f);
                const float e = H(i, j);
                Vector3 c = e <= 0.0f ? Vector3{ 0.15f, 0.3f, 0.6f } : e < 200.0f ? Vector3{ 0.35f, 0.6f, 0.3f } : e < 700.0f ? Vector3{ 0.5f, 0.5f, 0.35f } : Vector3{ 0.85f, 0.85f, 0.9f };
                if (area && e > 0.0f && (*area)[((size_t)f * N + j) * N + i] > 2e7f) c = { 0.1f, 0.3f, 0.9f };
                ImageDrawPixel(&img, i, j, { (unsigned char)(255 * c.x * shade), (unsigned char)(255 * c.y * shade), (unsigned char)(255 * c.z * shade), 255 });
            }
        ExportImage(img, path);
        UnloadImage(img);
    }
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
}

void PlanetErosion::build(const std::function<float(Vector3)>& height, float planetRadius, unsigned int seedIn, const ErosionParams& p, int threads) {
    const auto t0 = std::chrono::steady_clock::now();
    R = planetRadius; seed = seedIn;
    const int N = std::max(16, (int)lroundf(R * PI * 0.5f / p.cell));
    const size_t cells = (size_t)6 * N * N;
    auto idx = [N](int f, int i, int j) { return ((size_t)f * N + j) * N + i; };
    // Cell centres (unit directions) and areas.
    std::vector<Vector3> dir(cells);
    std::vector<float> cellArea(cells);
    const float du = 2.0f / N, side = R * du;
    for (int f = 0; f < 6; ++f) {
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
                const float u = -1.0f + (i + 0.5f) * du, v = -1.0f + (j + 0.5f) * du;
                dir[idx(f, i, j)] = Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v))));
                cellArea[idx(f, i, j)] = side * side / powf(1.0f + u * u + v * v, 1.5f);
            }
    }
    auto cellOf = [&](Vector3 d) -> int {
        const int f = FaceOf(d);
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        const float s = Vector3DotProduct(d, nrm);
        const int i = std::clamp((int)floorf((Vector3DotProduct(d, t1) / s + 1.0f) / du), 0, N - 1);
        const int j = std::clamp((int)floorf((Vector3DotProduct(d, t2) / s + 1.0f) / du), 0, N - 1);
        return (int)idx(f, i, j);
    };
    // 8 neighbours per cell, across face edges by direction.
    std::vector<int> nb(cells * 8);
    static const int di[8] = { 1, -1, 0, 0, 1, 1, -1, -1 }, dj[8] = { 0, 0, 1, -1, 1, -1, 1, -1 };
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
    // Input heights, sampled on threads.
    const auto ts = std::chrono::steady_clock::now();
    std::vector<float> h0(cells);
    {
        const int T = std::max(1, threads);
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] { for (size_t c = t; c < cells; c += T) h0[c] = height(dir[c]); });
        for (auto& th : pool) th.join();
    }
    sampleSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - ts).count();
    // Rock: erodibility and the slope it holds. Hard regions (~25 km) keep steep faces; strata regions alternate hard
    // caprock and soft layers every 120 m of height (mesas, plateaus, cliff bands).
    std::vector<float> rockK(cells), hardness(cells), strata(cells);
    float maxH0 = 1.0f;
    for (size_t c = 0; c < cells; ++c) {
        const Vector3 q = Vector3Scale(dir[c], R / 25000.0f);
        hardness[c] = ValueNoise(q, seed ^ 0xA11u);
        strata[c] = std::clamp((ValueNoise(Vector3Scale(q, 0.7f), seed ^ 0x57Au) - 0.55f) / 0.15f, 0.0f, 1.0f);
        rockK[c] = p.k * (1.6f - 1.2f * hardness[c]);
        maxH0 = std::max(maxH0, h0[c]);
    }
    const float stratumOffset = Hash3(1, 2, 3, seed) * 120.0f;
    auto capRock = [&](float hc) { return fmodf(hc + stratumOffset + 1200.0f, 120.0f) < 36.0f; }; // the hard 30% of each band
    auto erodibility = [&](size_t c, float hc) {
        if (strata[c] <= 0.0f) return rockK[c];
        const float bandK = capRock(hc) ? 0.15f : 1.6f;
        return rockK[c] * (1.0f + (bandK - 1.0f) * strata[c]);
    };
    // Start low; the relief grows from uplift where the input is high.
    std::vector<float> h(cells), U(cells);
    std::vector<uint8_t> outlet(cells);
    for (size_t c = 0; c < cells; ++c) {
        outlet[c] = h0[c] <= 0.0f;
        U[c] = h0[c] > 0.0f ? p.uplift * h0[c] / maxH0 : 0.0f;
        h[c] = h0[c] > 0.0f ? h0[c] * 0.1f : h0[c];
    }
    std::vector<float> hf(cells), area(cells);
    std::vector<int> receiver(cells), order;
    std::vector<uint8_t> seen(cells);
    order.reserve(cells);
    using Item = std::pair<float, int>;
    auto flood = [&] {
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
        std::fill(seen.begin(), seen.end(), 0);
        order.clear();
        for (size_t c = 0; c < cells; ++c)
            if (outlet[c]) { hf[c] = h[c]; receiver[c] = (int)c; seen[c] = 1; pq.push({ h[c], (int)c }); }
        if (pq.empty()) { // no sea: the lowest cell drains the planet
            size_t lo = 0; for (size_t c = 1; c < cells; ++c) if (h[c] < h[lo]) lo = c;
            hf[lo] = h[lo]; receiver[lo] = (int)lo; seen[lo] = 1; pq.push({ h[lo], (int)lo });
        }
        while (!pq.empty()) {
            const int c = pq.top().second; pq.pop();
            order.push_back(c);
            for (int k = 0; k < 8; ++k) {
                const int o = nb[(size_t)c * 8 + k];
                if (o < 0 || seen[o]) continue;
                seen[o] = 1;
                hf[o] = std::max(h[o], hf[c] + 1e-3f);
                receiver[o] = c;
                pq.push({ hf[o], o });
            }
        }
        // Steepest descent on the filled surface (the flood parent where none is lower).
        for (int c : order) {
            if (receiver[c] == c) continue;
            float best = 0.0f;
            for (int k = 0; k < 8; ++k) {
                const int o = nb[(size_t)c * 8 + k];
                if (o < 0) continue;
                const float s = (hf[c] - hf[o]) / (Vector3Distance(dir[c], dir[o]) * R);
                if (s > best) { best = s; receiver[c] = o; }
            }
        }
        for (size_t c = 0; c < cells; ++c) area[c] = cellArea[c];
        for (size_t k = order.size(); k-- > 0;) { const int c = order[k]; if (receiver[c] != c) area[receiver[c]] += area[c]; }
    };
    for (int step = 0; step < p.steps; ++step) {
        if (step % std::max(1, p.floodEvery) == 0) flood();
        // Implicit stream power, receivers first (`order` rises on the filled surface).
        for (int c : order) {
            const int r = receiver[c];
            if (r == c || outlet[c]) continue;
            if (hf[c] > h[c] + 1e-3f) { h[c] += U[c]; continue; } // a filled depression: a lake bed, no incision
            const float L = Vector3Distance(dir[c], dir[r]) * R;
            const float f = erodibility(c, h[c]) * powf(area[c], p.m) / L;
            h[c] = std::max((h[c] + U[c] + f * h[r]) / (1.0f + f), h[r]);
        }
        // Talus: slopes steeper than the rock holds slump to it (half the excess moves down, mass kept).
        if (step % 2 == 1)
            for (size_t c = 0; c < cells; ++c) {
                if (outlet[c]) continue;
                const float hard = std::max(hardness[c], strata[c] * (capRock(h[c]) ? 1.0f : 0.0f));
                const float maxSlope = p.talusSoft + (p.talusHard - p.talusSoft) * hard;
                for (int k = 0; k < 8; ++k) {
                    const int o = nb[c * 8 + k];
                    if (o < 0) continue;
                    const float excess = h[c] - h[o] - maxSlope * Vector3Distance(dir[c], dir[o]) * R;
                    if (excess > 0.0f) { const float mv = 0.5f * excess; h[c] -= mv; if (!outlet[o]) h[o] += mv; }
                }
            }
    }
    // Back to the input's height range: the highest land as high as the input's.
    float maxH = 1.0f;
    for (size_t c = 0; c < cells; ++c) if (!outlet[c]) maxH = std::max(maxH, h[c]);
    delta.assign(cells, 0.0f);
    maxCut = maxRaise = 0.0f;
    for (size_t c = 0; c < cells; ++c) {
        if (outlet[c]) continue;
        const float hb = std::max(0.5f, h[c] * maxH0 / maxH);
        delta[c] = hb - h0[c];
        maxCut = std::max(maxCut, -delta[c]); maxRaise = std::max(maxRaise, delta[c]);
    }
    if (const char* d = getenv("RAYLIBGAME_EROSION_DUMP"); d && d[0] == '1') {
        debugBefore = h0; debugAfter.resize(cells); debugArea = area;
        for (size_t c = 0; c < cells; ++c) debugAfter[c] = h0[c] + delta[c];
        int bestFace = 0, bestLand = -1;
        for (int fc = 0; fc < 6; ++fc) { int land = 0; for (int k = 0; k < N * N; ++k) land += h0[(size_t)fc * N * N + k] > 0.0f; if (land > bestLand) { bestLand = land; bestFace = fc; } }
        DumpFace("erosion_before.png", debugBefore, nullptr, bestFace, N, p.cell);
        DumpFace("erosion_after.png", debugAfter, &debugArea, bestFace, N, p.cell);
    }
    n = N;
    computeHash();
    buildSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
}

float PlanetErosion::sample(Vector3 d) const {
    if (n <= 0) return 0.0f;
    const int f = FaceOf(d);
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float s = Vector3DotProduct(d, nrm);
    const float du = 2.0f / n;
    const float fi = std::clamp((Vector3DotProduct(d, t1) / s + 1.0f) / du - 0.5f, 0.0f, (float)(n - 1) - 1e-3f);
    const float fj = std::clamp((Vector3DotProduct(d, t2) / s + 1.0f) / du - 0.5f, 0.0f, (float)(n - 1) - 1e-3f);
    const int i = (int)fi, j = (int)fj;
    const float a = fi - i, b = fj - j;
    const float* r0 = &delta[((size_t)f * n + j) * n + i];
    const float* r1 = r0 + n;
    return (r0[0] * (1 - a) + r0[1] * a) * (1 - b) + (r1[0] * (1 - a) + r1[1] * a) * b;
}

void PlanetErosion::computeHash() {
    uint64_t hv = 1469598103934665603ull; // FNV-1a over the stored (quantised) deltas
    for (float v : delta) {
        const int16_t q = (int16_t)std::clamp((int)lroundf(v * 10.0f), -32767, 32767);
        hv = (hv ^ (uint16_t)q) * 1099511628211ull;
    }
    hashValue = hv;
}

bool PlanetErosion::save(const std::string& path) const {
    if (n <= 0) return false;
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return false;
    const char magic[4] = { 'E', 'R', 'O', '1' };
    fwrite(magic, 1, 4, f);
    fwrite(&n, sizeof(n), 1, f); fwrite(&R, sizeof(R), 1, f); fwrite(&seed, sizeof(seed), 1, f); fwrite(&hashValue, sizeof(hashValue), 1, f);
    std::vector<int16_t> q(delta.size());
    for (size_t c = 0; c < delta.size(); ++c) q[c] = (int16_t)std::clamp((int)lroundf(delta[c] * 10.0f), -32767, 32767);
    fwrite(q.data(), sizeof(int16_t), q.size(), f);
    fclose(f);
    return true;
}

bool PlanetErosion::load(const std::string& path, float planetRadius, unsigned int wantSeed) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    char magic[4]; int N = 0; float r = 0; unsigned int s = 0; uint64_t hv = 0;
    bool ok = fread(magic, 1, 4, f) == 4 && memcmp(magic, "ERO1", 4) == 0 && fread(&N, sizeof(N), 1, f) == 1 && fread(&r, sizeof(r), 1, f) == 1 &&
              fread(&s, sizeof(s), 1, f) == 1 && fread(&hv, sizeof(hv), 1, f) == 1 && N > 0 && N < 4096 && s == wantSeed && fabsf(r - planetRadius) < 1.0f;
    std::vector<int16_t> q;
    if (ok) { q.resize((size_t)6 * N * N); ok = fread(q.data(), sizeof(int16_t), q.size(), f) == q.size(); }
    fclose(f);
    if (!ok) return false;
    delta.resize(q.size());
    for (size_t c = 0; c < q.size(); ++c) delta[c] = q[c] / 10.0f;
    n = N; R = r; seed = s;
    computeHash();
    if (hashValue != hv) { n = 0; delta.clear(); return false; } // damaged
    return true;
}
