#include "canyons.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <queue>
#include <thread>

namespace {
    // The cube faces as the erosion bake lays them out (erosion.cpp FaceAxes): face = axis * 2 + (positive ? 1 : 0).
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
    constexpr float kMinDepth = 20.0f;  // m: the graded floor stays at least this far under the land (still a canyon)
    constexpr float kMinGrade = 0.001f; // the floor falls at least 1 m per km toward its outlet
    constexpr float kNotchDepth = 60.0f; // the bake input's notch (PlanetField::land): land = input + this x notch
    constexpr float kMaxCut = 400.0f;   // m: deeper than this end to end, a canyon drains both ways from a divide instead
    constexpr int kMargin = 6;          // cells past the corridor that carry its floor (the walls of deep cuts are wide)
}

Vector3 CanyonGrade::dirOf(int f, int i, int j) const {
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float du = 2.0f / n, u = -1.0f + (i + 0.5f) * du, v = -1.0f + (j + 0.5f) * du;
    return Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v))));
}

int CanyonGrade::cellOf(Vector3 dv) const {
    const int f = FaceOf(dv);
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float s = Vector3DotProduct(dv, nrm), du = 2.0f / n;
    const int i = std::clamp((int)floorf((Vector3DotProduct(dv, t1) / s + 1.0f) / du), 0, n - 1);
    const int j = std::clamp((int)floorf((Vector3DotProduct(dv, t2) / s + 1.0f) / du), 0, n - 1);
    return (f * n + j) * n + i;
}

void CanyonGrade::build(float planetRadius, float cell, const std::function<float(Vector3)>& notch, const std::function<float(Vector3)>& height, int threads) {
    const auto t0 = std::chrono::steady_clock::now();
    R = planetRadius;
    n = std::max(16, (int)lroundf(R * PI * 0.5f / cell));
    const int cells = 6 * n * n, T = std::max(1, threads);
    auto dirC = [&](int c) { return dirOf(c / (n * n), c % n, (c / n) % n); };
    static const int di[8] = { 1, -1, 0, 0, 1, 1, -1, -1 }, dj[8] = { 0, 0, 1, -1, 1, -1, 1, -1 };
    auto nbOf = [&](int c, int k) {
        const int f = c / (n * n), i = c % n + di[k], j = (c / n) % n + dj[k];
        if (i >= 0 && j >= 0 && i < n && j < n) return (f * n + j) * n + i;
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        const float du = 2.0f / n, u = -1.0f + (i + 0.5f) * du, v = -1.0f + (j + 0.5f) * du;
        const int o = cellOf(Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v)))));
        return o == c ? -1 : o;
    };
    // The notch everywhere (cheap noise), on threads: each cell its own value.
    std::vector<uint8_t> q(cells, 0);
    {
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] { for (int c = t; c < cells; c += T) q[c] = (uint8_t)std::clamp((int)lroundf(notch(dirC(c)) * 255.0f), 0, 255); });
        for (auto& th : pool) th.join();
    }
    std::vector<int> corr;
    for (int c = 0; c < cells; ++c) if (q[c] >= 128) corr.push_back(c);
    // The input height in the corridors (land() is costly; only here).
    std::vector<float> h(cells, NAN);
    {
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] { for (size_t k = t; k < corr.size(); k += T) h[corr[k]] = height(dirC(corr[k])); });
        for (auto& th : pool) th.join();
    }
    auto hAt = [&](int c) { if (std::isnan(h[c])) h[c] = height(dirC(c)); return h[c]; }; // serial use only
    floor.assign(cells, NAN);
    std::vector<int> comp(cells, -1), members, hops(cells, -1), queueV;
    std::vector<float> dist(cells, 1e30f);
    double cutSum = 0.0, fillSum = 0.0; long cutN = 0;
    canyons = 0; corridorCells = 0; divided = 0; maxCut = maxFill = 0.0f; lengthKm = 0.0f;
    auto inCorr = [&](int c) { return c >= 0 && q[c] >= 128 && h[c] > 0.0f; };
    for (int s0 : corr) {
        if (comp[s0] >= 0 || !inCorr(s0)) continue;
        members.assign(1, s0); comp[s0] = s0;
        for (size_t k = 0; k < members.size(); ++k)
            for (int e = 0; e < 8; ++e) { const int o = nbOf(members[k], e); if (inCorr(o) && comp[o] < 0) { comp[o] = s0; members.push_back(o); } }
        if (members.size() < 8) continue; // a speck of a line: left as the plain notch
        // The ends: the farthest cell (in steps) from any cell, then the farthest from that one.
        auto farthest = [&](int from) {
            for (int c : members) hops[c] = -1;
            queueV.assign(1, from); hops[from] = 0;
            int last = from;
            for (size_t k = 0; k < queueV.size(); ++k) {
                const int c = queueV[k]; last = c;
                for (int e = 0; e < 8; ++e) { const int o = nbOf(c, e); if (o >= 0 && comp[o] == s0 && hops[o] < 0) { hops[o] = hops[c] + 1; queueV.push_back(o); } }
            }
            return last;
        };
        const int A = farthest(members[0]), B = farthest(A);
        // Each end's outlet: the lowest ground within 3 cells of it outside the canyon (the sea counts).
        auto outletLevel = [&](int e0) {
            std::vector<int> ring(1, e0), seen(1, e0);
            float lo = 1e30f;
            for (int r = 0; r < 3; ++r) {
                std::vector<int> next;
                for (int c : ring)
                    for (int e = 0; e < 8; ++e) {
                        const int o = nbOf(c, e);
                        if (o < 0 || std::find(seen.begin(), seen.end(), o) != seen.end()) continue;
                        seen.push_back(o); next.push_back(o);
                        if (comp[o] != s0) lo = std::min(lo, hAt(o));
                    }
                ring.swap(next);
            }
            return lo;
        };
        const int X = outletLevel(A) <= outletLevel(B) ? A : B;
        // Distance along the canyon from the outlet (Dijkstra in metres).
        using QE = std::pair<float, int>;
        std::priority_queue<QE, std::vector<QE>, std::greater<QE>> pq;
        for (int c : members) dist[c] = 1e30f;
        dist[X] = 0.0f; pq.push({ 0.0f, X });
        float smax = 0.0f;
        while (!pq.empty()) {
            const auto [dc, c] = pq.top(); pq.pop();
            if (dc > dist[c]) continue;
            smax = std::max(smax, dc);
            const Vector3 pc = dirC(c);
            for (int e = 0; e < 8; ++e) {
                const int o = nbOf(c, e);
                if (o < 0 || comp[o] != s0) continue;
                const float nd = dc + Vector3Distance(pc, dirC(o)) * R;
                if (nd < dist[o]) { dist[o] = nd; pq.push({ nd, o }); }
            }
        }
        // Per 100 m bin of that distance: the notch's floor (the lowest cell), the most it may rise (minDepth under the
        // land), and how many cells.
        const int nb = (int)(smax / cell) + 1;
        std::vector<float> m(nb, 1e30f), U(nb, 1e30f), w(nb, 0.0f);
        for (int c : members) {
            const int b = std::min(nb - 1, (int)(dist[c] / cell));
            m[b] = std::min(m[b], h[c]);
            U[b] = std::min(U[b], h[c] + kNotchDepth * q[c] / 255.0f - kMinDepth);
            w[b] += 1.0f;
        }
        // fit(lo, hi, outletAtLo): weighted isotonic regression (pool adjacent violators) over bins lo..hi of the floor less
        // kMinGrade x the distance from the outlet end: the floor never falls going away from the outlet. Returns the
        // squared error; with `write`, sets P there, then holds it minDepth under the land (lowered, and below that the
        // floor keeps falling toward the outlet) and above the sea (an inland floor under it would be an outlet).
        std::vector<float> P(nb, 0.0f);
        struct Block { double sw, swy, swyy; int first, last; };
        std::vector<Block> st;
        auto fit = [&](int lo, int hi, bool outletAtLo, bool write) {
            st.clear();
            const int n2 = hi - lo + 1;
            for (int k = 0; k < n2; ++k) {
                const int b = outletAtLo ? lo + k : hi - k;
                if (w[b] <= 0.0f) continue;
                const double y = m[b] - kMinGrade * (k + 0.5f) * cell;
                st.push_back({ w[b], w[b] * y, w[b] * y * y, k, k });
                while (st.size() > 1 && st[st.size() - 2].swy / st[st.size() - 2].sw > st.back().swy / st.back().sw) {
                    const Block top = st.back(); st.pop_back();
                    st.back().sw += top.sw; st.back().swy += top.swy; st.back().swyy += top.swyy; st.back().last = top.last;
                }
            }
            double sse = 0.0;
            for (const Block& bl : st) sse += bl.swyy - bl.swy * bl.swy / bl.sw;
            if (write) {
                for (const Block& bl : st) for (int k = bl.first; k <= bl.last; ++k) P[outletAtLo ? lo + k : hi - k] = (float)(bl.swy / bl.sw) + kMinGrade * (k + 0.5f) * cell;
                float above = 1e30f;
                for (int k = n2 - 1; k >= 0; --k) {
                    const int b = outletAtLo ? lo + k : hi - k;
                    if (w[b] <= 0.0f) continue;
                    P[b] = std::max(std::min({ P[b], U[b], above - kMinGrade * cell }), 1.0f);
                    above = P[b];
                }
            }
            return sse;
        };
        fit(0, nb - 1, true, true);
        // End to end, unless that cuts a gorge deeper than kMaxCut through a range: then the canyon drains both ways from
        // a divide, placed where the two fits fit best.
        float worstCut = 0.0f;
        for (int b = 0; b < nb; ++b) if (w[b] > 0.0f) worstCut = std::max(worstCut, m[b] - P[b]);
        if (worstCut > kMaxCut && nb >= 2) {
            int bestK = -1; double best = 1e300;
            for (int k = 0; k + 1 < nb; ++k) {
                const double sse = fit(0, k, true, false) + fit(k + 1, nb - 1, false, false);
                if (sse < best) { best = sse; bestK = k; }
            }
            fit(0, bestK, true, true);
            fit(bestK + 1, nb - 1, false, true);
            ++divided;
        }
        for (int c : members) {
            const float fl = P[std::min(nb - 1, (int)(dist[c] / cell))];
            floor[c] = fl;
            const float dlt = h[c] - fl;
            if (dlt > 0.0f) { cutSum += dlt; maxCut = std::max(maxCut, dlt); } else { fillSum -= dlt; maxFill = std::max(maxFill, -dlt); }
            ++cutN;
        }
        ++canyons; corridorCells += (int)members.size(); lengthKm += smax / 1000.0f;
    }
    meanCut = cutN ? (float)(cutSum / cutN) : 0.0f; meanFill = cutN ? (float)(fillSum / cutN) : 0.0f;
    // The margin: each cell near a corridor takes the floor of the first corridor cell that reaches it (index order).
    {
        std::vector<int> frontier, next;
        for (int c = 0; c < cells; ++c) if (!std::isnan(floor[c])) frontier.push_back(c);
        for (int r = 0; r < kMargin; ++r) {
            next.clear();
            for (int c : frontier)
                for (int e = 0; e < 8; ++e) { const int o = nbOf(c, e); if (o >= 0 && std::isnan(floor[o])) { floor[o] = floor[c]; next.push_back(o); } }
            frontier.swap(next);
        }
    }
    seconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
}

float CanyonGrade::floorAt(Vector3 dv, float& coverage) const {
    coverage = 0.0f;
    if (floor.empty()) return 0.0f;
    const int f = FaceOf(dv);
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float s = Vector3DotProduct(dv, nrm), du = 2.0f / n;
    const float fi = std::clamp((Vector3DotProduct(dv, t1) / s + 1.0f) / du - 0.5f, 0.0f, (float)(n - 1) - 1e-3f);
    const float fj = std::clamp((Vector3DotProduct(dv, t2) / s + 1.0f) / du - 0.5f, 0.0f, (float)(n - 1) - 1e-3f);
    const int i = (int)fi, j = (int)fj;
    const float a = fi - i, b = fj - j;
    const float wq[4] = { (1 - a) * (1 - b), a * (1 - b), (1 - a) * b, a * b };
    const int cq[4] = { (f * n + j) * n + i, (f * n + j) * n + i + 1, (f * n + j + 1) * n + i, (f * n + j + 1) * n + i + 1 };
    float sum = 0.0f;
    for (int k = 0; k < 4; ++k) {
        const float v = floor[cq[k]];
        if (std::isnan(v)) continue;
        sum += wq[k] * v; coverage += wq[k];
    }
    return coverage > 0.0f ? sum / coverage : 0.0f;
}
