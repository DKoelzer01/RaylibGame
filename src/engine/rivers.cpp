#include "rivers.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <queue>
#include <thread>

namespace {
    const Vector3 kN[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
    const Vector3 kU[6] = { { 0, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 0 } };
    const Vector3 kV[6] = { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 } };

    constexpr float kTargetCell = 420.0f; // m
    constexpr float kRiverArea = 40.0f;   // cells draining here before a river starts (~7 km^2)
    constexpr float kDrop = 0.8f;         // water below the smooth land at a river node
    constexpr float kMinLevel = 0.3f;     // at the sea
    constexpr float kLakeDepth = 6.0f;    // a filled depression this deep is a lake

    uint32_t Hash(uint32_t x, uint32_t seed) {
        uint32_t h = x * 2654435761u ^ seed * 0x9e3779b9u;
        h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
        return h;
    }

    float SegDist(Vector3 p, Vector3 a, Vector3 b, float& t) {
        const Vector3 ab = Vector3Subtract(b, a);
        const float l2 = Vector3DotProduct(ab, ab);
        t = l2 > 0.0f ? std::clamp(Vector3DotProduct(Vector3Subtract(p, a), ab) / l2, 0.0f, 1.0f) : 0.0f;
        return Vector3Distance(p, Vector3Add(a, Vector3Scale(ab, t)));
    }

    Vector3 CatmullRom(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3, float t) {
        const float t2 = t * t, t3 = t2 * t;
        Vector3 r = Vector3Scale(p1, 2.0f);
        r = Vector3Add(r, Vector3Scale(Vector3Subtract(p2, p0), t));
        r = Vector3Add(r, Vector3Scale(Vector3Add(Vector3Subtract(Vector3Scale(p0, 2.0f), Vector3Scale(p1, 5.0f)), Vector3Subtract(Vector3Scale(p2, 4.0f), p3)), t2));
        r = Vector3Add(r, Vector3Scale(Vector3Add(Vector3Subtract(Vector3Scale(p1, 3.0f), p0), Vector3Subtract(p3, Vector3Scale(p2, 3.0f))), t3));
        return Vector3Scale(r, 0.5f);
    }
}

int RiverNetwork::cellOf(Vector3 d) const {
    int face = 0; float best = -2.0f;
    for (int f = 0; f < 6; ++f) { const float v = Vector3DotProduct(d, kN[f]); if (v > best) { best = v; face = f; } }
    const float da = 0.5f * PI / n;
    const float a = atanf(Vector3DotProduct(d, kU[face]) / best), b = atanf(Vector3DotProduct(d, kV[face]) / best);
    const int i = std::clamp((int)floorf((a + 0.25f * PI) / da), 0, n - 1), j = std::clamp((int)floorf((b + 0.25f * PI) / da), 0, n - 1);
    return (face * n + i) * n + j;
}

Vector3 RiverNetwork::centre(int c) const {
    const int face = c / (n * n), i = (c / n) % n, j = c % n;
    const float da = 0.5f * PI / n;
    const float a = -0.25f * PI + (i + 0.5f) * da, b = -0.25f * PI + (j + 0.5f) * da;
    return Vector3Normalize(Vector3Add(kN[face], Vector3Add(Vector3Scale(kU[face], tanf(a)), Vector3Scale(kV[face], tanf(b)))));
}

void RiverNetwork::neighbours(int c, int out[8]) const {
    const int face = c / (n * n), i = (c / n) % n, j = c % n;
    const float da = 0.5f * PI / n;
    int k = 0;
    for (int di = -1; di <= 1; ++di)
        for (int dj = -1; dj <= 1; ++dj) {
            if (di == 0 && dj == 0) continue;
            const float a = -0.25f * PI + (i + di + 0.5f) * da, b = -0.25f * PI + (j + dj + 0.5f) * da;
            const Vector3 d = Vector3Normalize(Vector3Add(kN[face], Vector3Add(Vector3Scale(kU[face], tanf(a)), Vector3Scale(kV[face], tanf(b)))));
            const int m = cellOf(d);
            out[k++] = m == c ? -1 : m;
        }
}

RiverNetwork::RiverNetwork(float radius, unsigned int seed, const std::function<float(Vector3)>& height) {
    R = radius;
    n = std::max(4, (int)ceilf(0.5f * PI * R / kTargetCell));
    cell = 0.5f * PI * R / n;
    const int N = 6 * n * n;
    st.cells = N; st.cellSize = cell;
    h.assign(N, 0.0f);
    {   // Heights on worker threads (the land function is pure).
        const int T = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] { for (int c = t; c < N; c += T) h[c] = height(centre(c)); });
        for (auto& th : pool) th.join();
    }
    // Priority flood from the sea: each cell drains to the cell that reached it; depressions fill to
    // their spill level (+1 cm per cell so every filled cell still drains).
    F.assign(N, 0.0f);
    parent.assign(N, -1);
    std::vector<unsigned char> seen(N, 0);
    std::vector<int> order;
    order.reserve(N);
    using QE = std::pair<float, int>;
    std::priority_queue<QE, std::vector<QE>, std::greater<QE>> q;
    // The sea, and the polar ice shelves (frozen sea raised a couple of metres): rivers end there.
    for (int c = 0; c < N; ++c) if (h[c] < 0.0f || (fabsf(centre(c).y) > 0.86f && h[c] < 3.0f)) h[c] = std::min(h[c], -0.01f);
    for (int c = 0; c < N; ++c) if (h[c] < 0.0f) { F[c] = h[c]; seen[c] = 1; q.push({ h[c], c }); }
    if (q.empty()) { int lo = 0; for (int c = 1; c < N; ++c) if (h[c] < h[lo]) lo = c; F[lo] = h[lo]; seen[lo] = 1; q.push({ h[lo], lo }); }
    int nb[8];
    while (!q.empty()) {
        const int c = q.top().second; q.pop();
        order.push_back(c);
        neighbours(c, nb);
        for (int k = 0; k < 8; ++k) {
            const int m = nb[k];
            if (m < 0 || seen[m]) continue;
            seen[m] = 1;
            F[m] = std::max(h[m], F[c] + 0.01f);
            parent[m] = c;
            q.push({ F[m], m });
        }
    }
    // Flow: each land cell adds itself, downstream.
    acc.assign(N, 0.0f);
    for (int c = 0; c < N; ++c) if (h[c] >= 0.0f) acc[c] = 1.0f;
    for (int k = (int)order.size() - 1; k >= 0; --k) { const int c = order[k]; if (parent[c] >= 0) acc[parent[c]] += acc[c]; }
    river.assign(N, 0);
    lake.assign(N, 0);
    lakeLevel.assign(N, 0.0f);
    for (int c = 0; c < N; ++c) {
        if (h[c] < 0.0f) continue;
        ++st.landCells;
        if (acc[c] >= kRiverArea) { river[c] = 1; ++st.riverCells; }
    }
    // Lakes: each basin the flood filled (connected filled land cells) that is deep enough. Its
    // level is one value over the whole basin, and is held 3 m under the land just outside it (the
    // basin's edges, and kLakeMargin beyond them): the flood saw only cell centres, and between
    // them the real ground can fall lower, which stood the lake above its surroundings on walls.
    {
        std::vector<int> group(N, -1), cells;
        std::vector<unsigned char> filled(N, 0);
        for (int c = 0; c < N; ++c) filled[c] = h[c] >= 0.0f && F[c] - h[c] > 0.05f;
        for (int s = 0; s < N; ++s) {
            if (!filled[s] || group[s] >= 0) continue;
            cells.clear();
            cells.push_back(s);
            group[s] = s;
            float level = F[s], lowest = h[s];
            for (size_t i = 0; i < cells.size(); ++i) {
                const int c = cells[i];
                level = std::min(level, F[c]);
                lowest = std::min(lowest, h[c]);
                neighbours(c, nb);
                for (int k = 0; k < 8; ++k)
                    if (nb[k] >= 0 && filled[nb[k]] && group[nb[k]] < 0) { group[nb[k]] = s; cells.push_back(nb[k]); }
            }
            if (level - lowest <= kLakeDepth || level <= 1.0f) continue; // a shallow dip: no lake
            // The land just outside: on each edge to a cell outside the basin, and kLakeMargin beyond.
            float rim = 1e30f;
            for (int c : cells) {
                neighbours(c, nb);
                const Vector3 pc = Vector3Scale(centre(c), R);
                for (int k = 0; k < 8; ++k) {
                    const int m = nb[k];
                    if (m < 0 || group[m] == s) continue;
                    const Vector3 pm = Vector3Scale(centre(m), R);
                    const Vector3 cm = Vector3Subtract(pm, pc);
                    const float len = Vector3Length(cm);
                    if (len < 1.0f) continue;
                    const bool diagonal = len > 1.25f * cell;
                    const Vector3 side = Vector3Normalize(Vector3CrossProduct(pc, cm));
                    for (float out : { 0.5f * cell, 0.5f * cell + kLakeMargin }) {
                        const Vector3 base = Vector3Add(pc, Vector3Scale(cm, (diagonal ? out * 1.41f : out) / len));
                        for (float s2 : { -0.5f, -0.25f, 0.0f, 0.25f, 0.5f }) {
                            if (diagonal && s2 != 0.0f) continue;
                            rim = std::min(rim, height(Vector3Normalize(Vector3Add(base, Vector3Scale(side, s2 * cell)))));
                        }
                    }
                }
            }
            level = std::min(level, rim - 3.0f);
            if (level - lowest <= kLakeDepth || level <= 1.0f) continue;
            for (int c : cells) { lake[c] = 1; lakeLevel[c] = level; ++st.lakeCells; }
        }
    }
    // Steps from each cell to the nearest lake cell (breadth first), for cheap "lake near" tests.
    lakeSteps.assign(N, 255);
    {
        std::vector<int> frontier, next;
        for (int c = 0; c < N; ++c) if (lake[c]) { lakeSteps[c] = 0; frontier.push_back(c); }
        for (int d = 1; d < 255 && !frontier.empty(); ++d) {
            next.clear();
            for (int c : frontier) {
                neighbours(c, nb);
                for (int k = 0; k < 8; ++k) if (nb[k] >= 0 && lakeSteps[nb[k]] == 255) { lakeSteps[nb[k]] = (unsigned char)d; next.push_back(nb[k]); }
            }
            frontier.swap(next);
        }
    }
    // Levels: a running minimum downstream (children before parents in reverse flood order).
    L.assign(N, 0.0f); hw.assign(N, 0.0f); dep.assign(N, 0.0f);
    mainChild.assign(N, -1);
    // A floor: no river runs below the next lake downstream (it must arrive at the lake's surface).
    // The floor only falls going downstream (a lake's level is above anything it drains to), and so
    // does the running minimum, so their maximum does too.
    std::vector<float> floorLevel(N, kMinLevel);
    for (int c : order) {
        const float below = parent[c] >= 0 ? floorLevel[parent[c]] : kMinLevel;
        floorLevel[c] = lake[c] ? std::max(lakeLevel[c], below) : below;
    }
    std::vector<float> minChild(N, 1e30f);
    for (int k = (int)order.size() - 1; k >= 0; --k) {
        const int c = order[k];
        if (!river[c]) continue;
        L[c] = std::max(floorLevel[c], std::min(lake[c] ? lakeLevel[c] : F[c] - kDrop, minChild[c]));
        const float s = sqrtf(acc[c] / kRiverArea);
        hw[c] = std::min(3.0f + 2.5f * s, 45.0f);
        dep[c] = std::min(1.2f + 0.6f * s, 6.0f);
        const int p = parent[c];
        if (p >= 0) {
            minChild[p] = std::min(minChild[p], L[c]);
            if (mainChild[p] < 0 || acc[c] > acc[mainChild[p]]) mainChild[p] = c;
        }
    }
    // Curves: each link c -> parent as 8 pieces, Catmull-Rom through the neighbouring nodes plus a
    // meander that vanishes at the nodes (so links join up).
    for (int c = 0; c < N; ++c) {
        if (!river[c] || parent[c] < 0) continue;
        const int p = parent[c];
        const bool sea = !river[p];
        if (sea) ++st.mouths;
        const Vector3 P1 = Vector3Scale(centre(c), R), P2 = Vector3Scale(centre(p), R);
        const Vector3 P0 = mainChild[c] >= 0 ? Vector3Scale(centre(mainChild[c]), R) : Vector3Add(P1, Vector3Subtract(P1, P2));
        const Vector3 P3 = !sea && parent[p] >= 0 ? Vector3Scale(centre(parent[p]), R) : Vector3Add(P2, Vector3Subtract(P2, P1));
        const float la = L[c], lb = sea ? kMinLevel : L[p];
        const float wa = hw[c], wb = sea ? hw[c] * 1.2f : hw[p];
        const float da = dep[c], db = sea ? dep[c] : dep[p];
        const float segLen = Vector3Distance(P1, P2);
        const uint32_t hs = Hash((uint32_t)c, seed);
        const float amp = ((hs >> 8) * (1.0f / 16777216.0f) * 2.0f - 1.0f) * std::min(0.18f * segLen, 70.0f);
        constexpr int K = 8;
        Vector3 prev = P1;
        for (int k = 1; k <= K; ++k) {
            const float t = (float)k / K;
            Vector3 q2 = CatmullRom(P0, P1, P2, P3, t);
            const Vector3 up = Vector3Normalize(q2);
            const Vector3 tan = Vector3Subtract(CatmullRom(P0, P1, P2, P3, std::min(1.0f, t + 0.01f)), CatmullRom(P0, P1, P2, P3, std::max(0.0f, t - 0.01f)));
            const Vector3 side = Vector3Normalize(Vector3CrossProduct(up, tan));
            q2 = Vector3Add(q2, Vector3Scale(side, amp * sinf(PI * t)));
            q2 = Vector3Scale(Vector3Normalize(q2), R);
            const float t0 = (float)(k - 1) / K;
            pieces.push_back({ prev, q2, la + (lb - la) * t0, la + (lb - la) * t, wa + (wb - wa) * t0, wa + (wb - wa) * t, da + (db - da) * t0, da + (db - da) * t });
            prev = q2;
        }
    }
    index.assign(N, {});
    for (int i = 0; i < (int)pieces.size(); ++i) {
        const Vector3 mid = Vector3Scale(Vector3Add(pieces[i].a, pieces[i].b), 0.5f);
        const int m = cellOf(Vector3Normalize(mid));
        index[m].push_back(i);
        neighbours(m, nb);
        for (int k = 0; k < 8; ++k) if (nb[k] >= 0) index[nb[k]].push_back(i);
    }
    for (auto& v : index) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); }
    // The longest river (source to sea), for stats and tests.
    std::vector<float> len(N, 0.0f);
    for (int k = (int)order.size() - 1; k >= 0; --k) { const int c = order[k]; if (river[c] && parent[c] >= 0) len[parent[c]] = std::max(len[parent[c]], len[c] + cell); }
    for (int c = 0; c < N; ++c) if (river[c]) st.longestKm = std::max(st.longestKm, (len[c] + cell) / 1000.0f);
}

RiverNetwork::Hit RiverNetwork::nearest(Vector3 dir, float reach) const {
    Hit hit;
    const auto& list = index[cellOf(dir)];
    if (list.empty()) return hit;
    const Vector3 p = Vector3Scale(dir, R);
    for (int i : list) {
        const Piece& pc = pieces[i];
        float t;
        const float d = SegDist(p, pc.a, pc.b, t);
        if (d >= hit.dist) continue;
        hit.dist = d;
        hit.level = pc.la + (pc.lb - pc.la) * t;
        hit.halfWidth = pc.wa + (pc.wb - pc.wa) * t;
        hit.depth = pc.da + (pc.db - pc.da) * t;
    }
    hit.found = hit.dist < reach;
    return hit;
}

float RiverNetwork::distToCell(Vector3 d, int c) const {
    // Equiangular cube-sphere: angles in the cell's face frame are (very nearly) arc lengths.
    const int face = c / (n * n), i = (c / n) % n, j = c % n;
    const float dn = Vector3DotProduct(d, kN[face]);
    if (dn < 0.2f) return 1e30f;
    const float da = 0.5f * PI / n;
    const float a = atanf(Vector3DotProduct(d, kU[face]) / dn), b = atanf(Vector3DotProduct(d, kV[face]) / dn);
    const float a0 = -0.25f * PI + i * da, b0 = -0.25f * PI + j * da;
    const float ea = std::max({ a0 - a, 0.0f, a - (a0 + da) }), eb = std::max({ b0 - b, 0.0f, b - (b0 + da) });
    return sqrtf(ea * ea + eb * eb) * R;
}

bool RiverNetwork::lakeAt(Vector3 dir, float& level, float& outside) const {
    const int c = cellOf(dir);
    if (lake[c]) { level = lakeLevel[c]; outside = 0.0f; return true; }
    if (lakeSteps[c] != 1) return false;
    // Next to a lake: within the margin of a lake cell (the nearest; the highest on a tie).
    int nb[8];
    neighbours(c, nb);
    bool found = false;
    for (int k = 0; k < 8; ++k) {
        const int m = nb[k];
        if (m < 0 || !lake[m]) continue;
        const float d = distToCell(dir, m);
        if (d >= 2.0f * kLakeMargin) continue; // the water ends at the margin; the bank fades out by twice it
        if (!found || d < outside - 0.01f || (d < outside + 0.01f && lakeLevel[m] > level)) { level = lakeLevel[m]; outside = d; found = true; }
    }
    return found;
}

bool RiverNetwork::lakeNear(Vector3 dir, float extent) const {
    const int s = lakeSteps[cellOf(dir)];
    if (s == 255) return false;
    // A cell s steps away is at least (s - 1) cells from the lake's cells.
    return (float)std::max(0, s - 1) * cell <= extent + kLakeMargin;
}

bool RiverNetwork::riverNear(Vector3 dir, float extent) const {
    const int c = cellOf(dir);
    if (index[c].empty()) return false;
    const Hit h2 = nearest(dir, 60.0f + extent);
    return h2.found && h2.dist < h2.halfWidth * 1.3f + extent;
}

int RiverNetwork::checkDownstream(float& worstRise, float& worstShrink) const {
    worstRise = 0.0f; worstShrink = 0.0f;
    int links = 0;
    for (int c = 0; c < (int)river.size(); ++c) {
        if (!river[c] || parent[c] < 0 || !river[parent[c]]) continue;
        ++links;
        worstRise = std::max(worstRise, L[parent[c]] - L[c]);
        worstShrink = std::max(worstShrink, hw[c] - hw[parent[c]]);
    }
    for (const Piece& pc : pieces) { worstRise = std::max(worstRise, pc.lb - pc.la); worstShrink = std::max(worstShrink, pc.wa - pc.wb); }
    return links;
}

void RiverNetwork::longestPath(float step, std::vector<Vector3>& dirs, std::vector<float>& levels) const {
    dirs.clear(); levels.clear();
    // The source with the longest run to the sea.
    int best = -1, bestLen = 0;
    for (int c = 0; c < (int)river.size(); ++c) {
        if (!river[c] || mainChild[c] >= 0) continue;
        int len = 0;
        for (int x = c; x >= 0 && river[x]; x = parent[x]) ++len;
        if (len > bestLen) { bestLen = len; best = c; }
    }
    if (best < 0) return;
    // Walk its pieces (each link's 8 pieces were stored together, in order).
    std::vector<int> firstPiece(river.size(), -1);
    int k = 0;
    for (int c = 0; c < (int)river.size(); ++c) if (river[c] && parent[c] >= 0) { firstPiece[c] = k; k += 8; }
    for (int x = best; x >= 0 && river[x]; x = parent[x]) {
        if (firstPiece[x] < 0) break;
        for (int i = firstPiece[x]; i < firstPiece[x] + 8; ++i) {
            const Piece& pc = pieces[i];
            const float len = Vector3Distance(pc.a, pc.b);
            for (float s = 0.0f; s < len; s += step) {
                const float t = s / len;
                dirs.push_back(Vector3Normalize(Vector3Lerp(pc.a, pc.b, t)));
                levels.push_back(pc.la + (pc.lb - pc.la) * t);
            }
        }
    }
}
