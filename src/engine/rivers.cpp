#include "rivers.h"
#include "envflag.h"
#include "erosion.h"
#include <unordered_map>

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <queue>
#include <thread>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstdlib>

namespace {
    const Vector3 kN[6] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
    const Vector3 kU[6] = { { 0, 1, 0 }, { 0, 1, 0 }, { 0, 0, 1 }, { 0, 0, 1 }, { 1, 0, 0 }, { 1, 0, 0 } };
    const Vector3 kV[6] = { { 0, 0, 1 }, { 0, 0, -1 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 } };

    constexpr float kDrop = 0.8f;         // water below the smooth land at a river node
    constexpr float kFallMin = 2.0f;      // m: a piece dropping more than this is a waterfall (#540)
    constexpr float kMinLevel = 0.3f;     // at the sea

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

RiverNetwork::RiverNetwork(float radius, unsigned int seed, const std::function<float(Vector3)>& height, const Params& params) {
    R = radius;
    n = std::max(4, (int)ceilf(0.5f * PI * R / params.cell));
    cell = 0.5f * PI * R / n;
    const float riverCells = params.riverArea / (cell * cell); // cells draining here before a river starts
    auto tStage = std::chrono::steady_clock::now();
    const bool timeStages = EnvFlag("RAYLIBGAME_RIVER_TIMING");
    auto stage = [&](const char* what) { if (!timeStages) return; const auto t = std::chrono::steady_clock::now(); printf("rivers: before %s: %.2f s\n", what, std::chrono::duration<float>(t - tStage).count()); tStage = t; };
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
    // The 8 neighbours of every cell, once, on threads: neighbours() costs trigonometry per call, and the flood and
    // the lake distances below called it ~20 M times (2.9 s of a 3.6 s build on the 200 m grid, #359).
    std::vector<int> nbTable((size_t)N * 8);
    {
        const int T = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] { for (int c = t; c < N; c += T) neighbours(c, &nbTable[(size_t)c * 8]); });
        for (auto& th : pool) th.join();
    }
    auto nbOf = [&nbTable](int c, int out[8]) { std::copy_n(&nbTable[(size_t)c * 8], 8, out); };
    stage("heights + neighbours");
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
        nbOf(c, nb);
        for (int k = 0; k < 8; ++k) {
            const int m = nb[k];
            if (m < 0 || seen[m]) continue;
            seen[m] = 1;
            F[m] = std::max(h[m], F[c] + 0.01f);
            parent[m] = c;
            q.push({ F[m], m });
        }
    }
    stage("flood");
    // Flow: each land cell adds itself, downstream.
    acc.assign(N, 0.0f);
    for (int c = 0; c < N; ++c) if (h[c] >= 0.0f) acc[c] = 1.0f;
    for (int k = (int)order.size() - 1; k >= 0; --k) { const int c = order[k]; if (parent[c] >= 0) acc[parent[c]] += acc[c]; }
    river.assign(N, 0);
    lake.assign(N, 0);
    lakeIdOf.assign(N, -1); lakes.clear();
    lakeLevel.assign(N, 0.0f);
    for (int c = 0; c < N; ++c) {
        if (h[c] < 0.0f) continue;
        ++st.landCells;
        if (acc[c] >= riverCells) { river[c] = 1; ++st.riverCells; }
    }
    stage("flow");
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
                nbOf(c, nb);
                for (int k = 0; k < 8; ++k)
                    if (nb[k] >= 0 && filled[nb[k]] && group[nb[k]] < 0) { group[nb[k]] = s; cells.push_back(nb[k]); }
            }
            if (level - lowest <= params.lakeDepth || level <= 1.0f) continue; // a shallow dip: no lake
            // The land just outside: on each edge to a cell outside the basin, and kLakeMargin beyond.
            float rim = 1e30f;
            for (int c : cells) {
                nbOf(c, nb);
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
            if (level - lowest <= params.lakeDepth || level <= 1.0f) continue;
            {   // The lake's id and extent (#360: small ones are drainable).
                Vector3 sum{ 0, 0, 0 };
                for (int c : cells) sum = Vector3Add(sum, centre(c));
                LakeInfo li; li.dir = Vector3Normalize(sum); li.level = level; li.depth = level - lowest; li.cells = (int)cells.size();
                for (int c : cells) li.r = std::max(li.r, Vector3Distance(li.dir, centre(c)) * R + 0.75f * cell);
                for (int c : cells) { lake[c] = 1; lakeLevel[c] = level; lakeIdOf[c] = (int)lakes.size(); ++st.lakeCells; }
                lakes.push_back(li);
            }
        }
    }
    stage("lakes");
    // Steps from each cell to the nearest lake cell (breadth first), for cheap "lake near" tests.
    lakeSteps.assign(N, 255);
    {
        std::vector<int> frontier, next;
        for (int c = 0; c < N; ++c) if (lake[c]) { lakeSteps[c] = 0; frontier.push_back(c); }
        for (int d = 1; d < 255 && !frontier.empty(); ++d) {
            next.clear();
            for (int c : frontier) {
                nbOf(c, nb);
                for (int k = 0; k < 8; ++k) if (nb[k] >= 0 && lakeSteps[nb[k]] == 255) { lakeSteps[nb[k]] = (unsigned char)d; next.push_back(nb[k]); }
            }
            frontier.swap(next);
        }
    }
    stage("lake distances");
    // The largest tributary of each node (the curves pass through it).
    mainChild.assign(N, -1);
    for (int c = 0; c < N; ++c) {
        if (!river[c]) continue;
        const int p = parent[c];
        if (p >= 0 && (mainChild[p] < 0 || acc[c] > acc[mainChild[p]])) mainChild[p] = c;
    }
    // Curves: each link c -> parent as K pieces, Catmull-Rom through the neighbouring nodes plus a meander that
    // vanishes at the nodes (so links join up). The points first: the levels below sample the ground along them.
    constexpr int K = 8;
    std::vector<int> curveOf(N, -1);
    std::vector<Vector3> curvePts;
    for (int c = 0; c < N; ++c) {
        if (!river[c] || parent[c] < 0) continue;
        const int p = parent[c];
        const bool sea = !river[p];
        const Vector3 P1 = Vector3Scale(centre(c), R), P2 = Vector3Scale(centre(p), R);
        const Vector3 P0 = mainChild[c] >= 0 ? Vector3Scale(centre(mainChild[c]), R) : Vector3Add(P1, Vector3Subtract(P1, P2));
        const Vector3 P3 = !sea && parent[p] >= 0 ? Vector3Scale(centre(parent[p]), R) : Vector3Add(P2, Vector3Subtract(P2, P1));
        const float segLen = Vector3Distance(P1, P2);
        const uint32_t hs = Hash((uint32_t)c, seed);
        const float amp = ((hs >> 8) * (1.0f / 16777216.0f) * 2.0f - 1.0f) * std::min(params.meander * segLen, 70.0f);
        curveOf[c] = (int)curvePts.size();
        curvePts.push_back(P1);
        for (int k = 1; k <= K; ++k) {
            const float t = (float)k / K;
            Vector3 q2 = CatmullRom(P0, P1, P2, P3, t);
            const Vector3 up = Vector3Normalize(q2);
            const Vector3 tan = Vector3Subtract(CatmullRom(P0, P1, P2, P3, std::min(1.0f, t + 0.01f)), CatmullRom(P0, P1, P2, P3, std::max(0.0f, t - 0.01f)));
            const Vector3 side = Vector3Normalize(Vector3CrossProduct(up, tan));
            q2 = Vector3Add(q2, Vector3Scale(side, amp * sinf(PI * t)));
            curvePts.push_back(Vector3Scale(Vector3Normalize(q2), R));
        }
    }
    stage("curves");
    // The fine ground along each link (#359): with no bank fill, a level above it would stand on nothing.
    std::vector<float> fineMin(N, 1e30f), pointGround; // pointGround: the fine ground at each curve point (#190)
    if (params.fineHeight) pointGround.assign(curvePts.size(), 1e30f);
    if (params.fineHeight) {
        std::vector<int> links;
        for (int c = 0; c < N; ++c) if (curveOf[c] >= 0) links.push_back(c);
        const int T = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] {
                for (size_t i = t; i < links.size(); i += T) {
                    const int c = links[i];
                    float m = 1e30f;
                    for (int k = 0; k <= K; ++k) { const float g = params.fineHeight(Vector3Normalize(curvePts[curveOf[c] + k])); pointGround[curveOf[c] + k] = g; m = std::min(m, g); }
                    fineMin[c] = m;
                }
            });
        for (auto& th : pool) th.join();
    }
    stage("fine ground");
    // Levels: a running minimum downstream (children before parents in reverse flood order).
    L.assign(N, 0.0f); hw.assign(N, 0.0f); dep.assign(N, 0.0f);
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
        // The fine ground lowers a level by at most 8 m: a ravine or deep gully crossing the link (90 m) dragged the
        // level down, the running minimum carried it on and the valley carving dug a canyon to it (#359). Over such a
        // notch the water spans it.
        const float fineCap = std::min(F[c], std::max(fineMin[c], F[c] - 8.0f));
        L[c] = std::max(floorLevel[c], std::min(lake[c] ? lakeLevel[c] : fineCap - kDrop, minChild[c]));
        const float s = sqrtf(acc[c] * cell * cell / 7.06e6f); // widths by drained area (m^2): as before for the same river
        hw[c] = std::min(3.0f + 2.5f * s, 45.0f);
        dep[c] = std::min(1.2f + 0.6f * s, 6.0f);
        if (parent[c] >= 0) minChild[parent[c]] = std::min(minChild[parent[c]], L[c]);
        if (h[c] > 700.0f) ++st.riverCellsAbove700;
        if (mainChild[c] < 0) st.highestSource = std::max(st.highestSource, h[c]);
    }
    // The level at each curve point: interpolated along the link, and with the fine ground (#190) at the sea where the
    // ground there is already under it (the coarse cells put the mouth offshore of the real coast, so a river stood
    // 1.2 m over the sea at the shoreline: the flowing water overfilled and never settled). Upstream first, a running
    // minimum within each link and from its children, so the level never rises downstream.
    // With the fine ground (the bake, #540): the level follows the ground, not a straight line between the nodes - that
    // stood as a wall of water where the ground fell faster than the line (a cliff) and cut a trench where it rose. Each
    // point takes the highest ground over itself and the next kAhead points (~50 m) less kDrop: a notch narrower than that
    // is spanned (a gully crossing would otherwise drag the river down, #359), while past a cliff, where the ground stays
    // low, the level drops - a waterfall, the river going on at the ground below it.
    constexpr int kAhead = 2;
    static const bool stepped = !(getenv("RAYLIBGAME_RIVER_STEPS") && getenv("RAYLIBGAME_RIVER_STEPS")[0] == '0'); // 0: the old straight levels (A/B)
    std::vector<float> pointLevel(curvePts.size(), 0.0f), endCap(N, 1e30f);
    for (int k = (int)order.size() - 1; k >= 0; --k) {
        const int c = order[k];
        if (curveOf[c] < 0) continue;
        const int p = parent[c];
        const float la = L[c], lb = !river[p] ? kMinLevel : L[p];
        const bool follow = stepped && !pointGround.empty() && !lake[c];
        auto groundAt = [&](int i) { // point i of this link; past its end, on into the parent's curve
            if (i <= K) return pointGround[curveOf[c] + i];
            if (river[p] && curveOf[p] >= 0) return pointGround[curveOf[p] + std::min(i - K, K)];
            return pointGround[curveOf[c] + K];
        };
        float run = endCap[c];
        for (int i = 0; i <= K; ++i) {
            float lv = la + (lb - la) * (float)i / K;
            if (follow) {
                float ahead = -1e30f;
                for (int j = 0; j <= kAhead; ++j) ahead = std::max(ahead, groundAt(i + j));
                lv = std::max(floorLevel[c], ahead - kDrop);
            }
            // Only at a real coast: the link's cell is lowland by the sea (inland, ground under 0 m is a ravine, a gully
            // or a basin, and dropping the river there cut a canyon), and never under the next lake downstream.
            if (!pointGround.empty() && pointGround[curveOf[c] + i] < 0.0f && h[c] < 3.0f) lv = std::max(floorLevel[c], std::min(lv, kMinLevel));
            run = std::min(run, lv);
            pointLevel[curveOf[c] + i] = run;
        }
        endCap[p] = std::min(endCap[p], run);
    }
    // Every link ends at the level its parent starts with (#540): a tributary following its own ground arrived above the
    // main river at the confluence (the surface rose there walking downstream); now it falls into it, as a hanging valley.
    if (stepped)
        for (int c = 0; c < N; ++c) {
            const int p = parent[c];
            if (curveOf[c] >= 0 && p >= 0 && river[p] && curveOf[p] >= 0)
                pointLevel[curveOf[c] + K] = std::min(pointLevel[curveOf[c] + K], pointLevel[curveOf[p]]);
        }
    for (int c = 0; c < N; ++c) {
        if (curveOf[c] < 0) continue;
        const int p = parent[c];
        const bool sea = !river[p];
        if (sea) ++st.mouths;
        const float wa = hw[c], wb = sea ? hw[c] * 1.2f : hw[p];
        const float da = dep[c], db = sea ? dep[c] : dep[p];
        // The reach's slope over the whole link (~200 m; its falls included): how fast it runs (#595: the bed material).
        const float linkLen = std::max(1.0f, Vector3Distance(curvePts[curveOf[c]], curvePts[curveOf[c] + K]));
        const float reachSlope = std::max(0.0f, pointLevel[curveOf[c]] - pointLevel[curveOf[c] + K]) / linkLen;
        for (int k = 1; k <= K; ++k) {
            const float t = (float)k / K, t0 = (float)(k - 1) / K;
            const Vector3 a = curvePts[curveOf[c] + k - 1], b = curvePts[curveOf[c] + k];
            pieces.push_back({ a, b, pointLevel[curveOf[c] + k - 1], pointLevel[curveOf[c] + k], wa + (wb - wa) * t0, wa + (wb - wa) * t, da + (db - da) * t0, da + (db - da) * t });
            Piece& pc = pieces.back();
            pc.mid = Vector3Scale(Vector3Add(a, b), 0.5f); pc.halfLen = 0.5f * Vector3Distance(a, b);
            pc.reachSlope = reachSlope;
            // A waterfall (#540): the level holds to the cliff's edge - where the ground first falls under it - and drops
            // there, not in a ramp of water standing out over the drop.
            if (stepped && params.fineHeight && pc.la - pc.lb > kFallMin) {
                auto ground = [&](float u) { return params.fineHeight(Vector3Normalize(Vector3Lerp(a, b, u))); };
                float lo = 0.0f, hi = 1.0f;
                for (int s = 1; s <= 8; ++s) if (ground(s / 8.0f) < pc.la) { hi = s / 8.0f; lo = (s - 1) / 8.0f; break; }
                if (hi < 1.0f) for (int it = 0; it < 5; ++it) { const float m = 0.5f * (lo + hi); (ground(m) < pc.la ? hi : lo) = m; }
                pc.fallT = hi;
            }
        }
    }
    stage("levels + pieces");
    // Index: each piece in every cell within the carving reach of it (one ring of 420 m cells; two of 200 m ones).
    index.assign(N, {});
    reachM = params.maxReach;
    const int rings = std::max(1, (int)ceilf(reachM / cell));
    std::vector<int> near, frontier, next;
    for (int i = 0; i < (int)pieces.size(); ++i) {
        const Vector3 mid = Vector3Scale(Vector3Add(pieces[i].a, pieces[i].b), 0.5f);
        near.assign(1, cellOf(Vector3Normalize(mid)));
        frontier = near;
        for (int r = 0; r < rings; ++r) {
            next.clear();
            for (int c : frontier) { nbOf(c, nb); for (int k = 0; k < 8; ++k) if (nb[k] >= 0 && std::find(near.begin(), near.end(), nb[k]) == near.end()) { near.push_back(nb[k]); next.push_back(nb[k]); } }
            frontier.swap(next);
        }
        for (int c : near) index[c].push_back(i);
    }
    for (auto& v : index) { std::sort(v.begin(), v.end()); v.erase(std::unique(v.begin(), v.end()), v.end()); }
    stage("index");
    // The longest river (source to sea), for stats and tests.
    std::vector<float> len(N, 0.0f);
    for (int k = (int)order.size() - 1; k >= 0; --k) { const int c = order[k]; if (river[c] && parent[c] >= 0) len[parent[c]] = std::max(len[parent[c]], len[c] + cell); }
    for (int c = 0; c < N; ++c) if (river[c]) st.longestKm = std::max(st.longestKm, (len[c] + cell) / 1000.0f);
}

static thread_local bool tNatural = false; // a NaturalScope is open on this thread (#656)
RiverNetwork::NaturalScope::NaturalScope() : prev(tNatural) { tNatural = true; }
RiverNetwork::NaturalScope::~NaturalScope() { tNatural = prev; }

RiverNetwork::Hit RiverNetwork::nearest(Vector3 dir, float reach) const {
    Hit hit;
    const auto& list = index[cellOf(dir)];
    if (list.empty()) return hit;
    const Vector3 p = Vector3Scale(dir, R);
    // The nearest piece - unless the point lies inside more than one channel (a tributary running beside the main river
    // in one canyon, two reaches across a bend): then the lowest of those levels, as the higher water would spill into the
    // lower one (on the strata's cliff bands the nearest one rose 8 m walking down the longest river).
    Hit inside;
    inside.level = 1e30f;
    auto fill = [&](Hit& h, const Piece& pc, float d, float t, int i) {
        h.dist = d;
        h.level = pc.fallT < 0.0f ? pc.la + (pc.lb - pc.la) * t : t < pc.fallT ? pc.la : pc.lb; // a waterfall steps (#540)
        h.halfWidth = pc.wa + (pc.wb - pc.wa) * t;
        h.depth = pc.da + (pc.db - pc.da) * t;
        h.down = Vector3Normalize(pc.la >= pc.lb ? Vector3Subtract(pc.b, pc.a) : Vector3Subtract(pc.a, pc.b)); // levels fall downstream
        h.slope = pc.reachSlope;
        h.piece = i;
        // Less flow than its natural (#656): a lower, narrower river; none, the water at its bed (drawn as nothing: a dry bed).
        if (flowAny.load(std::memory_order_relaxed) && flowCur && !tNatural) {
            const float f = flowCur[i].load(std::memory_order_relaxed);
            if (f < 0.999f) { h.level -= h.depth * (1.0f - powf(std::max(f, 0.0f), 0.4f)) + 0.3f * (1.0f - std::min(1.0f, f * 10.0f)); h.halfWidth *= sqrtf(std::max(f, 0.0f)); h.depth *= powf(std::max(f, 0.0f), 0.4f); }
        }
    };
    for (int i : list) {
        const Piece& pc = pieces[i];
        // Bounding sphere first: no closer than its surface (most pieces in a cell are far away). A piece can hold the
        // point within its half width (45 m at most) even when another is nearer.
        const float dc = Vector3Distance(p, pc.mid) - pc.halfLen;
        if (dc >= std::max(hit.dist, 46.0f)) continue;
        float t;
        const float d = SegDist(p, pc.a, pc.b, t);
        if (d < hit.dist) fill(hit, pc, d, t, i);
        if (d < pc.wa + (pc.wb - pc.wa) * t) {
            Hit h;
            fill(h, pc, d, t, i);
            if (h.level < inside.level) inside = h;
        }
    }
    if (inside.level < 1e29f) hit = inside;
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

bool RiverNetwork::lakeAt(Vector3 dir, float& level, float& outside) const { int id; return lakeAt(dir, level, outside, id); }

bool RiverNetwork::lakeAt(Vector3 dir, float& level, float& outside, int& id) const {
    if (bake) return fineLakeAt(dir, level, outside, id);
    const int c = cellOf(dir);
    id = -1;
    if (lake[c]) { level = lakeLevel[c]; outside = 0.0f; id = lakeIdOf[c]; return true; }
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
        if (!found || d < outside - 0.01f || (d < outside + 0.01f && lakeLevel[m] > level)) { level = lakeLevel[m]; outside = d; id = lakeIdOf[m]; found = true; }
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
    if (bake) { // bake-native: every node draining into another, and junctions continuous (#607): each piece ends at its outlet's level
        for (int k = 0; k < (int)nodeParent.size(); ++k) {
            const int p = nodeParent[k];
            if (p < 0) continue;
            ++links;
            worstRise = std::max(worstRise, fabsf(pieces[nodePiece[k]].lb - pieces[nodePiece[p]].la)); // a step at the junction
        }
    }
    return links;
}

void RiverNetwork::longestPath(float step, std::vector<Vector3>& dirs, std::vector<float>& levels) const {
    dirs.clear(); levels.clear();
    if (bake) { // bake-native (#547): the node with the longest run to the sea, then its pieces node by node
        const int M = (int)nodeCell.size();
        std::vector<float> toSea(M, -1.0f);
        std::vector<int> chain;
        for (int k = 0; k < M; ++k) {
            chain.clear();
            int x = k;
            while (x >= 0 && toSea[x] < 0.0f) { chain.push_back(x); x = nodeParent[x]; }
            float d = x >= 0 ? toSea[x] : 0.0f;
            for (int i = (int)chain.size() - 1; i >= 0; --i) { d += 2.0f * pieces[nodePiece[chain[i]]].halfLen; toSea[chain[i]] = d; }
        }
        int best = -1;
        for (int k = 0; k < M; ++k) if (best < 0 || toSea[k] > toSea[best]) best = k;
        for (int x = best; x >= 0; x = nodeParent[x]) {
            const Piece& pc = pieces[nodePiece[x]];
            const float len = 2.0f * pc.halfLen;
            for (float s = 0.0f; s < len; s += step) {
                const float t = s / len;
                dirs.push_back(Vector3Normalize(Vector3Lerp(pc.a, pc.b, t)));
                levels.push_back(pc.fallT < 0.0f ? pc.la + (pc.lb - pc.la) * t : t < pc.fallT ? pc.la : pc.lb);
            }
        }
        return;
    }
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
                levels.push_back(pc.fallT < 0.0f ? pc.la + (pc.lb - pc.la) * t : t < pc.fallT ? pc.la : pc.lb);
            }
        }
    }
}

RiverNetwork::RiverNetwork(float radius, unsigned int seed, const std::function<float(Vector3)>& height) : RiverNetwork(radius, seed, height, Params()) {}

void RiverNetwork::forEachCentrePoint(float step, const std::function<void(Vector3, float, float)>& fn) const {
    for (const Piece& pc : pieces) {
        const int n = std::max(1, (int)ceilf(2.0f * pc.halfLen / step));
        for (int i = 0; i < n; ++i) {
            const float t = (i + 0.5f) / n;
            const float level = pc.fallT < 0.0f ? pc.la + (pc.lb - pc.la) * t : t < pc.fallT ? pc.la : pc.lb;
            fn(Vector3Normalize(Vector3Lerp(pc.a, pc.b, t)), level, pc.da + (pc.db - pc.da) * t);
        }
    }
}

std::vector<RiverNetwork::Fall> RiverNetwork::waterfalls() const {
    std::vector<Fall> out;
    for (const Piece& pc : pieces) {
        if (pc.fallT < 0.0f) continue;
        Fall f;
        f.lip = Vector3Normalize(Vector3Lerp(pc.a, pc.b, pc.fallT));
        f.down = Vector3Normalize(Vector3Subtract(pc.b, pc.a));
        f.drop = pc.la - pc.lb;
        f.top = pc.la;
        f.halfWidth = pc.wa + (pc.wb - pc.wa) * pc.fallT;
        out.push_back(f);
    }
    return out;
}

// ---- Bake-native rivers (#547, #607, #608, #609) ----

RiverNetwork::RiverNetwork(const PlanetErosion& ero, float radius, const Params& params) {
    bake = &ero;
    R = radius;
    n = std::max(4, (int)ceilf(0.5f * PI * R / params.cell)); // the index grid (as the old network's cells)
    cell = 0.5f * PI * R / n;
    reachM = params.maxReach;
    const int N = 6 * n * n, G = ero.flowCells();
    const float fineM = 0.5f * PI * R / sqrtf((float)G / 6.0f); // ~50 m
    st.cells = G; st.cellSize = fineM;
    // The nodes: every bake cell on land draining params.riverArea or more, in cell order.
    std::unordered_map<int, int> nodeOf;
    for (int g = 0; g < G; ++g) {
        if (ero.flowHeight(g) <= 0.0f) continue;
        ++st.landCells;
        if (ero.flowArea(g) >= params.riverArea) { nodeOf[g] = (int)nodeCell.size(); nodeCell.push_back(g); }
    }
    const int M = (int)nodeCell.size();
    st.riverCells = M;
    // The bake's heights as the runtime ground has them (#903): biomes reshape the ground after the bake (a swamp pulled 24 m
    // down left its rivers standing in the air), so every height a level comes from moves with it.
    auto shifted = [&](int g, float h) { return params.regionalShift ? h + params.regionalShift(ero.flowDir(g), h) : h; };
    auto groundOf = [&](int g) { return shifted(g, ero.flowHeight(g)); };
    std::unordered_map<int, float> fillLevel;
    for (const auto& fc : ero.filledCells()) fillLevel[fc.first] = shifted(fc.first, fc.second);
    // Lakes: the flood's filled basins deeper than params.lakeDepth (level: the basin's highest filled height - the flood's
    // jitter climbs across a basin's flat floor, and the lowest stood under the ground at its far side: dry lake beds).
    {
        std::unordered_map<int, int> comp;
        int nb[8];
        for (const auto& fc : ero.filledCells()) {
            if (comp.count(fc.first)) continue;
            std::vector<int> cells(1, fc.first);
            comp[fc.first] = fc.first;
            float level = fillLevel[fc.first], lowest = groundOf(fc.first);
            for (size_t i = 0; i < cells.size(); ++i) {
                ero.flowNeighbours(cells[i], nb);
                for (int k = 0; k < 8; ++k)
                    if (nb[k] >= 0 && !comp.count(nb[k]) && fillLevel.count(nb[k])) { comp[nb[k]] = fc.first; cells.push_back(nb[k]); }
                level = std::max(level, fillLevel[cells[i]]);
                lowest = std::min(lowest, groundOf(cells[i]));
            }
            if (level - lowest <= params.lakeDepth || level <= 1.0f) continue;
            LakeInfo li; Vector3 sum{ 0, 0, 0 };
            for (int c : cells) sum = Vector3Add(sum, ero.flowDir(c));
            li.dir = Vector3Normalize(sum); li.level = level; li.depth = level - lowest; li.cells = (int)cells.size();
            for (int c : cells) li.r = std::max(li.r, Vector3Distance(li.dir, ero.flowDir(c)) * R + fineM);
            const int id = (int)lakes.size();
            for (int c : cells) { fineLake[c] = id; fillLevel[c] = level; }
            st.lakeCells += (int)cells.size();
            lakes.push_back(li);
        }
    }
    // Links and levels: a node's own level is its lake's, or the filled surface, or its worn channel's ground plus 0.8 of
    // the incision the bake gave it (the water fills most of the channel). Then, upstream first, never above what drains
    // into it, and never under the next lake downstream (a river arrives at a lake's surface).
    // A river widens and deepens from nothing at its head (#625: one was born 8 m wide at a waterfall): 0 at riverArea, full at 4x.
    const float riverA = params.riverArea;
    auto headRamp = [riverA](float a) { return std::clamp((a / riverA - 1.0f) / 3.0f, 0.0f, 1.0f); };
    nodeParent.assign(M, -1);
    std::vector<float> own(M), area(M), L(M, 0.0f), floorL(M, kMinLevel);
    std::vector<int> indeg(M, 0), seaCell(M, -1), mainChildN(M, -1);
    for (int k = 0; k < M; ++k) {
        const int g = nodeCell[k];
        area[k] = ero.flowArea(g);
        const auto fl = fillLevel.find(g);
        // A head starts as a seep (#625): its water rises from the ground with the drainage (none at riverArea, the full 0.8 of
        // the incision at 4x it).
        own[k] = fl != fillLevel.end() ? fl->second : groundOf(g) + 0.8f * PlanetErosion::IncisionDepth(area[k]) * headRamp(area[k]);
        const int r = ero.flowReceiver(g);
        const auto it = r >= 0 ? nodeOf.find(r) : nodeOf.end();
        if (it != nodeOf.end()) { nodeParent[k] = it->second; ++indeg[it->second]; }
        else { seaCell[k] = r; ++st.mouths; }
    }
    std::vector<int> topo; topo.reserve(M); // sources first (Kahn, in node order)
    for (int k = 0; k < M; ++k) if (indeg[k] == 0) topo.push_back(k);
    for (size_t i = 0; i < topo.size(); ++i) { const int p = nodeParent[topo[i]]; if (p >= 0 && --indeg[p] == 0) topo.push_back(p); }
    for (int i = (int)topo.size() - 1; i >= 0; --i) { // downstream first: the next lake's level
        const int k = topo[i], p = nodeParent[k];
        const float below = p >= 0 ? floorL[p] : kMinLevel;
        const auto lk = fineLake.find(nodeCell[k]);
        floorL[k] = lk != fineLake.end() ? std::max(lakes[lk->second].level, below) : below;
    }
    for (int k : topo) { const int p = nodeParent[k]; if (p >= 0 && (mainChildN[p] < 0 || area[k] > area[mainChildN[p]])) mainChildN[p] = k; }
    // Positions: each node's cell centre eased toward its main tributary's and its outlet's (shared by every link that
    // meets there, so the curves join exactly and the 50 m grid's stairs soften).
    std::vector<Vector3> centreN(M), pos(M);
    for (int k = 0; k < M; ++k) centreN[k] = ero.flowDir(nodeCell[k]);
    for (int k = 0; k < M; ++k) {
        const Vector3 a = mainChildN[k] >= 0 ? centreN[mainChildN[k]] : centreN[k];
        const Vector3 b = nodeParent[k] >= 0 ? centreN[nodeParent[k]] : centreN[k];
        pos[k] = Vector3Normalize(Vector3Add(Vector3Scale(centreN[k], 0.5f), Vector3Add(Vector3Scale(a, 0.25f), Vector3Scale(b, 0.25f))));
    }
    // Widths and depths by drainage. A head's trickle is at least 2 m wide and 0.5 m deep (#631: a 1.6 m x 0.3 m one showed
    // as an empty trench - the water mesh draws nothing under 0.15 m, and coarser terrain never resolved the channel).
    // Discharge (#637, user: most rivers 10-20 m wide, the big ones where several join 30-50 m and several metres deep): each
    // node's own drainage times its runoff (params.runoff: wetter climates and a busier water cycle give more water per km^2),
    // plus all its tributaries bring - so a river never narrows flowing into drier country. Hydraulic geometry on it: full
    // width 6 + 8 s m (s = sqrt(discharge area / 7.06 km^2): ~14 m at s = 1, ~38 m at s = 4, at most 60 m), depth 1 + 0.7 s m
    // (at most 6 m); a head still ramps up from its trickle.
    std::vector<float> wetArea = area; // the discharge, as an area drained at the default runoff (m^2)
    if (params.runoff) {
        std::vector<float> runoff(M, 1.0f);
        const int T = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t) pool.emplace_back([&, t] { for (int k = t; k < M; k += T) runoff[k] = std::max(0.05f, params.runoff(pos[k])); });
        for (auto& th : pool) th.join();
        // Each node's own drainage (its area less its tributaries') at its runoff, accumulated downstream (sources first).
        std::vector<float> childArea(M, 0.0f);
        for (int k = 0; k < M; ++k) if (nodeParent[k] >= 0) childArea[nodeParent[k]] += area[k];
        std::vector<float> acc(M, 0.0f);
        for (int k = 0; k < M; ++k) acc[k] = std::max(0.0f, area[k] - childArea[k]) * runoff[k];
        for (int k : topo) { wetArea[k] = acc[k]; if (nodeParent[k] >= 0) acc[nodeParent[k]] += acc[k]; }
    }
    auto sOf = [&](int k) { return sqrtf(wetArea[k] / 7.06e6f); };
    auto widthOf = [&](int k) { return std::max(1.0f, std::min(3.0f + 4.0f * sOf(k), 30.0f) * headRamp(area[k])); }; // half width
    auto depthOf = [&](int k) { return std::max(0.5f, std::min(1.0f + 0.7f * sOf(k), 6.0f) * headRamp(area[k])); };
    // Banks (#640): the bake's level stood over the real ground beside a third of the rivers (the runtime adds detail the
    // bake's 50 m cells average away), and the water ended there as a wall. Each node's level is held 0.3 m under the bank
    // ground sampled 2 m and 6 m beyond both edges, at the node and halfway along its link, lowered by at most kMaxLower
    // (a deeper deficit is a slope the river runs along, not a bank: cutting that far made trenches). In parallel: ~6 field
    // samples a node.
    if (params.bankHeight) {
        constexpr float kMaxLower = 3.0f;
        const int T = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] {
                for (int k = t; k < M; k += T) {
                    if (fillLevel.count(nodeCell[k])) continue; // lakes and filled basins keep the flood's level
                    const Vector3 P = Vector3Scale(pos[k], R);
                    const Vector3 Q = nodeParent[k] >= 0 ? Vector3Scale(pos[nodeParent[k]], R) : P;
                    const Vector3 along = Vector3Length(Vector3Subtract(Q, P)) > 1.0f ? Vector3Normalize(Vector3Subtract(Q, P))
                                         : Vector3Normalize(Vector3CrossProduct(pos[k], fabsf(pos[k].y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                    const Vector3 side = Vector3Normalize(Vector3CrossProduct(along, pos[k]));
                    const float half = widthOf(k);
                    float cap = 1e30f;
                    const Vector3 at[2] = { P, Vector3Scale(Vector3Add(P, Q), 0.5f) };
                    for (int a = 0; a < 2; ++a)
                        for (float off : { half + 2.0f, half + 6.0f }) {
                            if (a == 1 && off > half + 2.0f) continue;
                            for (int s = -1; s <= 1; s += 2)
                                cap = std::min(cap, params.bankHeight(Vector3Normalize(Vector3Add(at[a], Vector3Scale(side, s * off)))) - 0.3f);
                        }
                    own[k] = std::max(std::min(own[k], cap), own[k] - kMaxLower);
                }
            });
        for (auto& th : pool) th.join();
    }
    std::vector<float> minChild(M, 1e30f);
    for (int k : topo) {
        L[k] = std::max(floorL[k], std::min(own[k], minChild[k]));
        const int p = nodeParent[k];
        if (p >= 0) minChild[p] = std::min(minChild[p], L[k]);
    }
    // Pieces: one per link, node to node (or into the sea cell).
    nodePiece.assign(M, -1);
    for (int k = 0; k < M; ++k) {
        const int p = nodeParent[k];
        const Vector3 A = Vector3Scale(pos[k], R);
        const Vector3 B = p >= 0 ? Vector3Scale(pos[p], R) : Vector3Scale(seaCell[k] >= 0 ? ero.flowDir(seaCell[k]) : pos[k], R);
        const float la = L[k], lb = p >= 0 ? L[p] : std::min(la, kMinLevel);
        const float wa = widthOf(k), wb = p >= 0 ? widthOf(p) : wa * 1.2f;
        const float da = depthOf(k), db = p >= 0 ? depthOf(p) : da;
        Piece pc{ A, B, la, lb, wa, wb, da, db };
        pc.mid = Vector3Scale(Vector3Add(A, B), 0.5f); pc.halfLen = 0.5f * Vector3Distance(A, B);
        pc.reachSlope = std::max(0.0f, la - lb) / std::max(1.0f, 2.0f * pc.halfLen);
        pc.inLake = fineLake.count(nodeCell[k]) && (p < 0 || fineLake.count(nodeCell[p]));
        // A waterfall only where the worn ground itself drops steeply (#540; more than 0.4 m per m over the link): there the level
        // holds to where the ground first falls under it. A gentler drop is a rapid, its level along the link: stepping
        // every steep 50 m link wore flat terraces into the slopes.
        // Never near a head (#625): a stream still forming runs down its slope as a rapid.
        if (params.fineHeight && headRamp(area[k]) >= 1.0f && la - lb > std::max(kFallMin, 0.8f * pc.halfLen)) {
            auto ground = [&](float u) { return params.fineHeight(Vector3Normalize(Vector3Lerp(A, B, u))); };
            float lo = 0.0f, hi = 1.0f;
            for (int s = 1; s <= 8; ++s) if (ground(s / 8.0f) < la) { hi = s / 8.0f; lo = (s - 1) / 8.0f; break; }
            if (hi < 1.0f) for (int it = 0; it < 5; ++it) { const float m = 0.5f * (lo + hi); (ground(m) < la ? hi : lo) = m; }
            pc.fallT = hi;
        }
        nodePiece[k] = (int)pieces.size();
        pieces.push_back(pc);
        if (mainChildN[k] < 0) st.highestSource = std::max(st.highestSource, ero.flowHeight(nodeCell[k]));
        if (ero.flowHeight(nodeCell[k]) > 700.0f) ++st.riverCellsAbove700;
    }
    // River flow state (#656): each node's discharge (its drainage) for mixing at confluences, the order sources first, and
    // every piece at its natural flow (node k's piece is piece k).
    nodeQ.assign(area.begin(), area.end());
    nodeChildQ.assign(nodeQ.size(), 0.0f);
    for (int k = 0; k < M; ++k) if (nodeParent[k] >= 0) nodeChildQ[nodeParent[k]] += std::max(nodeQ[k], 1.0f);
    nodeTopo = topo;
    flowCur.reset(new std::atomic<float>[pieces.size()]);
    for (size_t i = 0; i < pieces.size(); ++i) flowCur[i].store(1.0f, std::memory_order_relaxed);
    flowSet.assign(pieces.size(), -1.0f);
    flowSent.assign(pieces.size(), 1.0f);
    // Heads (#625): every source, and whether a waterfall comes within its first 3 links.
    for (int k = 0; k < M; ++k) {
        if (mainChildN[k] >= 0) continue;
        ++st.heads;
        bool fall = false;
        for (int x = k, prev = -1, n2 = 0; x >= 0 && n2 < 3 && !fall && (prev < 0 || mainChildN[x] == prev); prev = x, x = nodeParent[x], ++n2)
            fall = pieces[nodePiece[x]].fallT >= 0.0f; // its own stem only (a short tributary joining a falling river is not its head)
        st.headFalls += fall;
    }
    // The index (the old network's grid, ~200 m cells): each piece in its cell and the ring around it; lakes flagged.
    std::vector<int> nbTable((size_t)N * 8);
    {
        const int T = std::max(1u, std::min(16u, std::thread::hardware_concurrency()));
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t) pool.emplace_back([&, t] { for (int c = t; c < N; c += T) neighbours(c, &nbTable[(size_t)c * 8]); });
        for (auto& th : pool) th.join();
    }
    index.assign(N, {});
    const int rings = std::max(1, (int)ceilf(reachM / cell));
    std::vector<int> nearC, frontier, next;
    for (int i = 0; i < (int)pieces.size(); ++i) {
        nearC.assign(1, cellOf(Vector3Normalize(pieces[i].mid)));
        frontier = nearC;
        for (int r = 0; r < rings; ++r) {
            next.clear();
            for (int c : frontier) for (int k = 0; k < 8; ++k) { const int o = nbTable[(size_t)c * 8 + k]; if (o >= 0 && std::find(nearC.begin(), nearC.end(), o) == nearC.end()) { nearC.push_back(o); next.push_back(o); } }
            frontier.swap(next);
        }
        for (int c : nearC) index[c].push_back(i);
    }
    lake.assign(N, 0);
    for (const auto& fl : fineLake) lake[cellOf(ero.flowDir(fl.first))] = 1;
    lakeSteps.assign(N, 255);
    {
        std::vector<int> fr, nx;
        for (int c = 0; c < N; ++c) if (lake[c]) { lakeSteps[c] = 0; fr.push_back(c); }
        for (int d = 1; d < 255 && !fr.empty(); ++d) {
            nx.clear();
            for (int c : fr) for (int k = 0; k < 8; ++k) { const int o = nbTable[(size_t)c * 8 + k]; if (o >= 0 && lakeSteps[o] == 255) { lakeSteps[o] = (unsigned char)d; nx.push_back(o); } }
            fr.swap(nx);
        }
    }
    // The longest river (source to sea), for stats and tests.
    std::vector<float> len(M, 0.0f);
    for (int k : topo) {
        const int p = nodeParent[k];
        const float l = len[k] + 2.0f * pieces[nodePiece[k]].halfLen;
        if (p >= 0) len[p] = std::max(len[p], l); else st.longestKm = std::max(st.longestKm, l / 1000.0f);
    }
}

void RiverNetwork::forEachLink(const std::function<void(const Link&)>& fn) const {
    for (const Piece& pc : pieces) fn({ pc.a, pc.b, pc.la, pc.lb, 2.0f * pc.wa, 2.0f * pc.wb, pc.da, pc.db, pc.inLake }); // wa/wb are half widths
}

bool RiverNetwork::fineLakeAt(Vector3 dir, float& level, float& outside, int& id) const {
    const int g = bake->flowCellOf(dir);
    id = -1;
    if (const auto it = fineLake.find(g); it != fineLake.end()) { id = it->second; level = lakes[id].level; outside = 0.0f; return true; }
    if (lakeSteps[cellOf(dir)] > 1) return false; // no lake cell within ~200 m (cheap)
    // Beside a lake: the nearest lake cell within two cells (its edge, approximately: the centre less half a cell).
    const float half = 0.5f * st.cellSize;
    int nb[8], nb2[8];
    bool found = false;
    bake->flowNeighbours(g, nb);
    auto consider = [&](int m) {
        const auto it = fineLake.find(m);
        if (it == fineLake.end()) return;
        const float d = std::max(0.0f, Vector3Distance(dir, bake->flowDir(m)) * R - half);
        if (d >= 2.0f * kLakeMargin) return;
        if (!found || d < outside) { outside = d; level = lakes[it->second].level; id = it->second; found = true; }
    };
    for (int k = 0; k < 8; ++k) {
        if (nb[k] < 0) continue;
        consider(nb[k]);
        bake->flowNeighbours(nb[k], nb2);
        for (int k2 = 0; k2 < 8; ++k2) if (nb2[k2] >= 0) consider(nb2[k2]);
    }
    return found;
}

bool RiverNetwork::nearestHead(Vector3 dir, float minLevel, Vector3& head, Vector3& down) const {
    if (!bake) return false;
    std::vector<uint8_t> hasChild(nodeParent.size(), 0);
    for (int p : nodeParent) if (p >= 0) hasChild[p] = 1;
    float best = 1e30f; int bk = -1;
    for (int k = 0; k < (int)nodeParent.size(); ++k) {
        if (hasChild[k] || nodeParent[k] < 0) continue;
        const Piece& pc = pieces[nodePiece[k]];
        if (pc.la < minLevel) continue;
        const float d = Vector3Distance(Vector3Normalize(pc.a), dir);
        if (d < best) { best = d; bk = k; }
    }
    if (bk < 0) return false;
    const Piece& pc = pieces[nodePiece[bk]];
    head = Vector3Normalize(pc.a);
    down = Vector3Normalize(Vector3Subtract(pc.b, pc.a));
    return true;
}

float RiverNetwork::carve(Vector3 dir, float h) const {
    const auto& list = index[cellOf(dir)];
    if (list.empty()) return h;
    const Vector3 p = Vector3Scale(dir, R);
    constexpr float kBank = 0.6f, kBankReach = 15.0f; // the bank's slope just beyond the channel, and how far it reaches
    for (int i : list) {
        const Piece& pc = pieces[i];
        if (pc.inLake) continue;
        if (Vector3Distance(p, pc.mid) - pc.halfLen >= 45.0f + kBankReach) continue;
        float t;
        const float d = SegDist(p, pc.a, pc.b, t);
        const float half = pc.wa + (pc.wb - pc.wa) * t;
        if (d >= half + kBankReach) continue;
        const float level = pc.fallT < 0.0f ? pc.la + (pc.lb - pc.la) * t : t < pc.fallT ? pc.la : pc.lb;
        const float depth = pc.da + (pc.db - pc.da) * t;
        const float dep = std::min(depth, 1.0f); // light: deeper water only where the worn ground is lower (#609)
        const float bed = d < half ? level + 0.3f - (dep + 0.3f) * (1.0f - (d / half) * (d / half)) : level + 0.3f + (d - half) * kBank;
        h = std::min(h, bed);
    }
    return h;
}

// ---------------------------------------------------------------------------------------------
// River flow (#656)
// ---------------------------------------------------------------------------------------------

float RiverNetwork::flowOf(int piece) const {
    if (!flowCur || piece < 0 || piece >= (int)pieces.size()) return 1.0f;
    return flowCur[piece].load(std::memory_order_relaxed);
}

void RiverNetwork::setReachFlow(int piece, float share) const {
    std::lock_guard<std::mutex> lock(flowMx);
    if (!flowCur || piece < 0 || piece >= (int)pieces.size()) return;
    flowSet[piece] = std::clamp(share, 0.0f, 1.5f);
    flowAny.store(true, std::memory_order_relaxed);
}

void RiverNetwork::clearReachFlow(int piece) const {
    std::lock_guard<std::mutex> lock(flowMx);
    if (!flowCur || piece < 0 || piece >= (int)pieces.size()) return;
    flowSet[piece] = -1.0f;
}

bool RiverNetwork::stepFlow(float dt, std::vector<int>& changed) const {
    std::lock_guard<std::mutex> lock(flowMx);
    if (!flowCur || !flowAny.load(std::memory_order_relaxed)) return false;
    // Sources first: each node's target is the share of its discharge arriving - its own drainage at its natural rate plus
    // what its tributaries carry now - unless a released reach sets it; it eases there, ~95% of the way in the time a wave takes along its link at kSpeed (so a
    // change runs down the river at that speed). Lakes keep their level (their whole basin feeds them).
    constexpr float kSpeed = 5.0f; // m/s: a flow wave runs ahead of the water itself
    std::vector<float> inflow(nodeQ.size(), 0.0f);
    bool off = false;
    for (int k : nodeTopo) {
        const Piece& pc = pieces[k];
        const float q = std::max(nodeQ[k], 1.0f);
        // Against the node's own natural total (its tributaries' plus its own drainage), so an untouched node is exactly 1: a
        // node whose drainage is a little under its tributaries' sum read > 1 and drifted, and 250k links remeshed (#656).
        const float own = std::max(0.0f, q - nodeChildQ[k]);
        const float natural = nodeChildQ[k] + own;
        float target = flowSet[k] >= 0.0f ? flowSet[k] : std::clamp((inflow[k] + own) / std::max(natural, 1.0f), 0.0f, 1.5f);
        if (pc.inLake) target = 1.0f;
        float f = flowCur[k].load(std::memory_order_relaxed);
        f += (target - f) * std::min(1.0f, 3.0f * dt * kSpeed / std::max(1.0f, 2.0f * pc.halfLen));
        if (fabsf(f - 1.0f) < 0.002f && flowSet[k] < 0.0f) f = 1.0f;
        flowCur[k].store(f, std::memory_order_relaxed);
        off = off || f != 1.0f || flowSet[k] >= 0.0f;
        if (fabsf(f - flowSent[k]) > 0.05f || (f == 1.0f && flowSent[k] != 1.0f)) { changed.push_back(k); flowSent[k] = f; }
        const int p = nodeParent[k];
        if (p >= 0) inflow[p] += q * f;
    }
    flowAny.store(off, std::memory_order_relaxed);
    return !changed.empty();
}

RiverNetwork::PieceBounds RiverNetwork::pieceBounds(int piece) const {
    if (piece < 0 || piece >= (int)pieces.size()) return { Vector3{}, 0.0f }; // stale (a queued refresh from another world)
    const Piece& pc = pieces[piece];
    return { pc.mid, pc.halfLen + std::max(pc.wa, pc.wb) + 20.0f };
}

void RiverNetwork::flowState(std::vector<std::pair<int, float>>& cur, std::vector<std::pair<int, float>>& set) const {
    std::lock_guard<std::mutex> lock(flowMx);
    if (!flowCur) return;
    for (size_t i = 0; i < pieces.size(); ++i) {
        const float f = flowCur[i].load(std::memory_order_relaxed);
        if (f != 1.0f) cur.push_back({ (int)i, f });
        if (flowSet[i] >= 0.0f) set.push_back({ (int)i, flowSet[i] });
    }
}

void RiverNetwork::setFlowState(const std::vector<std::pair<int, float>>& cur, const std::vector<std::pair<int, float>>& set) const {
    std::lock_guard<std::mutex> lock(flowMx);
    if (!flowCur) return;
    bool any = false;
    for (const auto& [i, f] : cur) if (i >= 0 && i < (int)pieces.size()) { flowCur[i].store(f, std::memory_order_relaxed); flowSent[i] = f; any = true; }
    for (const auto& [i, f] : set) if (i >= 0 && i < (int)pieces.size()) { flowSet[i] = f; any = true; }
    flowAny.store(any || flowAny.load(), std::memory_order_relaxed);
}

int RiverNetwork::pieceBeyond(int piece, Vector3 centre, float dist) const {
    if (piece < 0 || piece >= (int)pieces.size()) return -1;
    for (int k = nodeParent[piece], j = 0; k >= 0 && k < (int)pieces.size() && j < 64; k = nodeParent[k], ++j) // from the next link: this one starts upstream
        if (Vector3Distance(pieces[k].a, centre) > dist) return k;
    return -1;
}
