#include "erosion.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <queue>
#include <cstdlib>


// Prototype evidence (#352): hillshaded height maps before and after, drainage in blue (RAYLIBGAME_EROSION_DUMP=1).
static void DumpMap(const char* path, const std::vector<float>& h, const std::vector<float>* area, int N, float dx) {
    const int S = 3;
    Image img = GenImageColor(N * S, N * S, BLACK);
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) {
            auto H = [&](int a, int b) { a = std::clamp(a, 0, N - 1); b = std::clamp(b, 0, N - 1); return h[(size_t)b * N + a]; };
            const float gx = (H(i + 1, j) - H(i - 1, j)) / (2 * dx), gy = (H(i, j + 1) - H(i, j - 1)) / (2 * dx);
            const Vector3 nrm = Vector3Normalize({ -gx, 1.0f, -gy });
            const float shade = std::clamp(0.35f + 0.75f * Vector3DotProduct(nrm, Vector3Normalize({ -0.6f, 0.7f, -0.4f })), 0.0f, 1.0f);
            const float e = h[(size_t)j * N + i];
            Vector3 c = e <= 0.0f ? Vector3{ 0.15f, 0.3f, 0.6f } : e < 200.0f ? Vector3{ 0.35f, 0.6f, 0.3f } : e < 700.0f ? Vector3{ 0.5f, 0.5f, 0.35f } : Vector3{ 0.85f, 0.85f, 0.9f };
            if (area && e > 0.0f && (*area)[(size_t)j * N + i] > 4e6f) c = { 0.1f, 0.3f, 0.9f };
            const Color col = { (unsigned char)(255 * c.x * shade), (unsigned char)(255 * c.y * shade), (unsigned char)(255 * c.z * shade), 255 };
            ImageDrawRectangle(&img, i * S, j * S, S, S, col);
        }
    ExportImage(img, path);
    UnloadImage(img);
}
void ErosionBake::build(const std::function<float(Vector3)>& height, Vector3 centre, float planetRadius, const ErosionParams& p) {
    const auto t0 = std::chrono::steady_clock::now();
    nrm = Vector3Normalize(centre);
    t1 = Vector3Normalize(fabsf(nrm.y) < 0.9f ? Vector3CrossProduct(nrm, { 0, 1, 0 }) : Vector3CrossProduct(nrm, { 1, 0, 0 }));
    t2 = Vector3CrossProduct(nrm, t1);
    R = planetRadius; dx = p.dx; half = 0.5f * p.extent; fade = std::max(p.fade, dx);
    const int N = std::max(8, (int)(p.extent / p.dx));
    const size_t cells = (size_t)N * N;
    std::vector<float> h0(cells), h(cells), hf(cells);
    auto dirOf = [&](int i, int j) {
        const float u = -half + (i + 0.5f) * dx, v = -half + (j + 0.5f) * dx;
        return Vector3Normalize(Vector3Add(Vector3Scale(nrm, R), Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v))));
    };
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i) h0[(size_t)j * N + i] = height(dirOf(i, j));
    h = h0;
    std::vector<float> uplift;
    float maxH0 = 1.0f;
    if (p.upliftMode) { // start at 10%: the relief grows from uplift where the input is high and erosion carves it
        for (float v : h0) maxH0 = std::max(maxH0, v);
        uplift.resize(cells);
        for (size_t c = 0; c < cells; ++c) { uplift[c] = h0[c] > 0.0f ? h0[c] / maxH0 * 6.0f : 0.0f; h[c] = h0[c] > 0.0f ? h0[c] * 0.1f : h0[c]; }
    }
    // Outlets: the grid edge and the sea (below sea level). Never eroded.
    std::vector<uint8_t> outlet(cells, 0);
    for (int j = 0; j < N; ++j)
        for (int i = 0; i < N; ++i)
            if (i == 0 || j == 0 || i == N - 1 || j == N - 1 || h0[(size_t)j * N + i] <= 0.0f) outlet[(size_t)j * N + i] = 1;
    static const int di[8] = { 1, -1, 0, 0, 1, 1, -1, -1 }, dj[8] = { 0, 0, 1, -1, 1, -1, 1, -1 };
    static const float dl[8] = { 1, 1, 1, 1, 1.41421356f, 1.41421356f, 1.41421356f, 1.41421356f };
    std::vector<int> receiver(cells), order;
    std::vector<float> area(cells);
    std::vector<uint8_t> seen(cells);
    order.reserve(cells);
    using Item = std::pair<float, int>;
    for (int step = 0; step < p.steps; ++step) {
        // Priority flood from the outlets: filled heights (a hair of slope across lakes) and a processing order in
        // which every cell's receiver comes first.
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
        std::fill(seen.begin(), seen.end(), 0);
        order.clear();
        for (size_t c = 0; c < cells; ++c)
            if (outlet[c]) { hf[c] = h[c]; pq.push({ h[c], (int)c }); seen[c] = 1; receiver[c] = (int)c; }
        while (!pq.empty()) {
            const int c = pq.top().second; pq.pop();
            order.push_back(c);
            const int ci = c % N, cj = c / N;
            for (int k = 0; k < 8; ++k) {
                const int ni = ci + di[k], nj = cj + dj[k];
                if (ni < 0 || nj < 0 || ni >= N || nj >= N) continue;
                const int nb = nj * N + ni;
                if (seen[nb]) continue;
                seen[nb] = 1;
                hf[nb] = std::max(h[nb], hf[c] + 1e-3f);
                receiver[nb] = c; // the flood's parent: always lower (filled)
                pq.push({ hf[nb], nb });
            }
        }
        // Steepest descent on the filled surface where there is a lower neighbour (the flood parent otherwise).
        for (int c : order) {
            if (outlet[c]) continue;
            const int ci = c % N, cj = c / N;
            float best = 0.0f;
            for (int k = 0; k < 8; ++k) {
                const int ni = ci + di[k], nj = cj + dj[k];
                if (ni < 0 || nj < 0 || ni >= N || nj >= N) continue;
                const int nb = nj * N + ni;
                const float s = (hf[c] - hf[nb]) / (dl[k] * dx);
                if (s > best) { best = s; receiver[c] = nb; }
            }
        }
        // Receivers are strictly lower on hf, and `order` is by rising hf: receivers come first. Drainage area: donors
        // add into their receivers, last to first.
        std::fill(area.begin(), area.end(), dx * dx);
        for (size_t k = order.size(); k-- > 0;) { const int c = order[k]; if (receiver[c] != c) area[receiver[c]] += area[c]; }
        // Implicit stream power (n = 1), outlets first.
        for (int c : order) {
            if (outlet[c]) continue;
            const int r = receiver[c];
            const float relaxUp = p.upliftMode ? uplift[c] : p.relax * (h0[c] - h[c]);
            if (hf[c] > h[c] + 1e-3f) { h[c] += relaxUp; continue; } // in a filled depression: a lake bed, no incision
            const int ri = r % N, rj = r / N, ci = c % N, cj = c / N;
            const float L = dx * sqrtf((float)((ri - ci) * (ri - ci) + (rj - cj) * (rj - cj)));
            const float f = p.k * powf(area[c], p.m) / L;
            const float hn = (h[c] + relaxUp + f * h[r]) / (1.0f + f);
            h[c] = std::max(hn, h[r]); // never below its receiver
        }
        // Hillslope diffusion (explicit, 4 neighbours).
        hf = h;
        for (int j = 1; j < N - 1; ++j)
            for (int i = 1; i < N - 1; ++i) {
                const size_t c = (size_t)j * N + i;
                if (outlet[c]) continue;
                h[c] += p.diffusion * (hf[c - 1] + hf[c + 1] + hf[c - N] + hf[c + N] - 4.0f * hf[c]);
            }
    }
    if (p.upliftMode) { // back to the input's range: the highest land as high as the input's
        float maxH = 1.0f;
        for (size_t c = 0; c < cells; ++c) if (h0[c] > 0.0f) maxH = std::max(maxH, h[c]);
        for (size_t c = 0; c < cells; ++c) if (h0[c] > 0.0f) h[c] = std::max(0.5f, h[c] * maxH0 / maxH);
    }
    if (const char* d = getenv("RAYLIBGAME_EROSION_DUMP"); d && d[0] == '1') { DumpMap("erosion_before.png", h0, nullptr, N, dx); DumpMap("erosion_after.png", h, &area, N, dx); }
    delta.resize(cells);
    maxCut = maxFill = 0.0f;
    for (size_t c = 0; c < cells; ++c) {
        delta[c] = h[c] - h0[c];
        maxCut = std::max(maxCut, -delta[c]); maxFill = std::max(maxFill, delta[c]);
    }
    n = N;
    buildSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
}

bool ErosionBake::sample(Vector3 dir, float& d, float& weight) const {
    if (n <= 0) return false;
    const float s = Vector3DotProduct(dir, nrm);
    if (s <= 0.0f) return false;
    const float u = Vector3DotProduct(dir, t1) / s * R, v = Vector3DotProduct(dir, t2) / s * R;
    if (u <= -half || u >= half || v <= -half || v >= half) return false;
    const float fi = std::clamp((u + half) / dx - 0.5f, 0.0f, (float)(n - 1) - 1e-3f);
    const float fj = std::clamp((v + half) / dx - 0.5f, 0.0f, (float)(n - 1) - 1e-3f);
    const int i = (int)fi, j = (int)fj;
    const float a = fi - i, b = fj - j;
    const float* r0 = &delta[(size_t)j * n + i];
    const float* r1 = r0 + n;
    d = (r0[0] * (1 - a) + r0[1] * a) * (1 - b) + (r1[0] * (1 - a) + r1[1] * a) * b;
    const float edge = std::min({ u + half, half - u, v + half, half - v });
    const float t = std::clamp(edge / fade, 0.0f, 1.0f);
    weight = t * t * (3.0f - 2.0f * t);
    return true;
}
