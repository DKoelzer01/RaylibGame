#include "erosion.h"
#include "utils/logger.h"
#include <raymath.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <queue>
#include <thread>
#include <atomic>
#include <condition_variable>
#include <memory>
#include <mutex>

// Bake progress 0..1 for the loading screen (#362): the coarse planet pass is the first 30%, the six face refinements the rest.
static std::atomic<float> gBakeProgress{ 0.0f };
static constexpr int kPieces = 16; // the cache file deflates in this many pieces at once (#362)
float ErosionBakeProgress() { return gBakeProgress.load(std::memory_order_relaxed); }
static std::atomic<bool> gBaking{ false };
static std::atomic<int> gBakes{ 0 };
bool ErosionBaking() { return gBaking.load(std::memory_order_relaxed); }
int ErosionBakeCount() { return gBakes.load(); }
// #605 A/B: RAYLIBGAME_BAKE_OLD=1 runs the bake as before its speedups (same result, the old timing).
static const bool gBakeOld = getenv("RAYLIBGAME_BAKE_OLD") && getenv("RAYLIBGAME_BAKE_OLD")[0] == '1';

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
        std::vector<float> cellArea, h, U, rockK, hardness, strata, glacier;
        std::vector<StrataColumn> col;      // #587: the rock column per cell (empty: the old hardness noise and bands)
        unsigned int strataSeed = 0;
        std::vector<uint8_t> outlet;
    };

    // Rock and uplift of a domain from its input heights (h0) and the planet's highest input (maxH0).
    // Each cell's own values only, so `threads` (interleaved cells) never changes the result (#605).
    void SetRock(Domain& d, const std::vector<float>& h0, float maxH0, float R, unsigned int seed, const ErosionParams& p, const StrataFn* strata, int threads) {
        if (strata && *strata) { d.col.resize(d.cells); d.strataSeed = seed; } // #587
        d.U.resize(d.cells); d.rockK.resize(d.cells); d.hardness.resize(d.cells); d.strata.resize(d.cells); d.glacier.resize(d.cells);
        const int T = gBakeOld ? 1 : std::max(1, threads);
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t) pool.emplace_back([&, t] {
        for (size_t c = t; c < d.cells; c += T) {
            if (!d.col.empty()) d.col[c] = (*strata)(d.dir[c]);
            const Vector3 q = Vector3Scale(d.dir[c], R / 25000.0f);
            d.hardness[c] = ValueNoise(q, seed ^ 0xA11u);
            d.strata[c] = std::clamp((ValueNoise(Vector3Scale(q, 0.7f), seed ^ 0x57Au) - 0.55f) / 0.15f, 0.0f, 1.0f);
            d.rockK[c] = p.k * (1.6f - 1.2f * d.hardness[c]);
            // A floor of 5% for all land: uplift only in proportion to the input left wide low plains dead flat, with
            // straight 8-direction river traces across them (#355, seed 7).
            d.U[c] = h0[c] > 0.0f ? p.uplift * (0.05f + 0.95f * h0[c] / maxH0) : 0.0f;
            // Ice-age glaciers (#358): above a snowline of p.glacialLine m at the equator, falling to sea level by
            // |latitude| ~0.8; full strength 300 m above it.
            {
                const float lat = fabsf(d.dir[c].y);
                const float line = p.glacialLine * (1.0f - std::clamp((lat - 0.3f) / 0.5f, 0.0f, 1.0f));
                d.glacier[c] = p.glacial > 0.0f ? std::clamp((h0[c] - line) / 300.0f, 0.0f, 1.0f) : 0.0f;
            }
        }
        });
        for (auto& th : pool) th.join();
    }

    // Seconds per part of Simulate (#605): setup, flood, glacial, stream power, talus, wind.
    struct SimProfile { double s[6] = {}; void add(const SimProfile& o) { for (int k = 0; k < 6; ++k) s[k] += o.s[k]; } };
    const char* const kSimParts[6] = { "setup", "flood", "glacial", "stream", "talus", "wind" };
    // Faces refined at once (a face's working set is ~250 MB): 6, 3, 2 or 1 (whole rounds), at least 2 threads each. Each
    // face's simulation runs on one thread, so 3 at a time left 13 of 16 threads idle through it (#605).
    int kParallelFaces(int T) {
        if (gBakeOld) return std::clamp(T / 4, 1, 3);
        for (int k : { 6, 3, 2 }) if (T >= 2 * k) return k;
        return 1;
    }

    // A monotone radix heap for the priority flood (#605): it pops (height, cell) pairs in exactly the order the binary
    // heap of std::pair<float, int> did - the keys are the same total order (height first, then the cell index), so the
    // flood visits the cells in the same order and the bake is unchanged - but pushes are O(1) and pops cheap. Valid
    // because the flood never pushes a key below the last one popped (filled heights only rise outward).
    class RadixHeap {
    public:
        static uint64_t Key(float h, int c) {
            uint32_t b; memcpy(&b, &h, 4);
            b = (b & 0x80000000u) ? ~b : (b | 0x80000000u); // order-preserving float -> unsigned
            return ((uint64_t)b << 32) | (uint32_t)c;
        }
        void push(uint64_t k) { buckets[Bucket(k)].push_back(k); ++count; }
        bool empty() const { return count == 0; }
        uint64_t pop() {
            if (buckets[0].empty()) {
                int i = 1;
                while (buckets[i].empty()) ++i;
                last = *std::min_element(buckets[i].begin(), buckets[i].end());
                for (uint64_t k : buckets[i]) buckets[Bucket(k)].push_back(k);
                buckets[i].clear();
            }
            const uint64_t k = buckets[0].back(); buckets[0].pop_back(); --count;
            return k;
        }
        void reset() { for (auto& b : buckets) b.clear(); last = 0; count = 0; }
    private:
        std::vector<uint64_t> buckets[65];
        uint64_t last = 0;
        size_t count = 0;
        int Bucket(uint64_t k) const { return k == last ? 0 : 64 - __builtin_clzll(k ^ last); }
    };

    // A barrier for the stream pool: spins briefly, then yields.
    struct SpinBarrier {
        explicit SpinBarrier(int n) : n(n) {}
        void wait() {
            const int g = gen.load(std::memory_order_acquire);
            if (count.fetch_add(1, std::memory_order_acq_rel) == n - 1) { count.store(0, std::memory_order_relaxed); gen.fetch_add(1, std::memory_order_release); return; }
            for (int spins = 0; gen.load(std::memory_order_acquire) == g; ++spins) if (spins > 64) std::this_thread::yield(); // never sleeps: a sleep is a whole timer tick (up to 15.6 ms) per level, #621
        }
        const int n;
        std::atomic<int> count{ 0 }, gen{ 0 };
    };

    // The implicit stream power on several threads, with exactly the serial pass's result (#605). The serial pass walks
    // the flood order and reads each receiver's height as it stands: updated this step when the receiver came earlier in
    // the order, the step's starting height when later (a back edge). Here each cell gets a level one past its receiver's
    // (0 on a back edge: its receiver's starting height is saved first), the cells of one level run in parallel, and the
    // levels in turn; small levels run on the calling thread alone. Every cell computes what the serial pass computes.
    class StreamPool {
    public:
        explicit StreamPool(int threads) : T(threads), bar(threads) { for (int t = 1; t < T; ++t) workers.emplace_back([this, t] { loop(t); }); }
        ~StreamPool() { { std::lock_guard<std::mutex> l(m); quit = true; ++stepGen; } cv.notify_all(); for (auto& w : workers) w.join(); }
        void plan(const std::vector<int>& order, const std::vector<int>& recv, const std::vector<uint8_t>& outlet, const float* heights) {
            const size_t cells = recv.size();
            receiver = recv.data(); h = heights;
            pos.resize(cells); lvl.assign(cells, -1); back.assign(cells, 0); rv.resize(cells); backList.clear();
            for (size_t k = 0; k < order.size(); ++k) pos[order[k]] = (int)k;
            int maxLvl = -1;
            for (int c : order) {
                const int r = recv[c];
                if (r == c || outlet[c]) continue;
                if (pos[r] < pos[c]) lvl[c] = lvl[r] + 1; // an unchanging receiver (-1) puts it at level 0
                else { lvl[c] = 0; back[c] = 1; backList.push_back(c); }
                maxLvl = std::max(maxLvl, lvl[c]);
            }
            std::vector<int> start(maxLvl + 2, 0);
            for (int c : order) if (lvl[c] >= 0) ++start[lvl[c] + 1];
            for (int l = 0; l <= maxLvl; ++l) start[l + 1] += start[l];
            cellsByLvl.resize(start[maxLvl + 1]);
            { std::vector<int> at(start.begin(), start.end() - 1); for (int c : order) if (lvl[c] >= 0) cellsByLvl[at[lvl[c]]++] = c; }
            items.clear();
            for (int l = 0; l <= maxLvl; ++l) {
                const bool par = start[l + 1] - start[l] >= kParallelLevel;
                if (!par && !items.empty() && !items.back().par) items.back().end = start[l + 1];
                else items.push_back({ start[l], start[l + 1], par });
            }
        }
        template <class F> void stream(F& f) {
            fn = [](void* ctx, int c, float hr) { (*(F*)ctx)(c, hr); };
            ctx = &f;
            for (int c : backList) rv[c] = h[receiver[c]];
            { std::lock_guard<std::mutex> l(m); ++stepGen; } cv.notify_all(); // the workers start (they wait on the condition between steps, not spinning)
            work(0);
        }
    private:
        static constexpr int kParallelLevel = 4096;
        struct Item { int begin, end; bool par; };
        int T;
        SpinBarrier bar;
        std::vector<std::thread> workers;
        bool quit = false;
        std::mutex m; std::condition_variable cv; uint64_t stepGen = 0; // the step start (#621)
        const int* receiver = nullptr; const float* h = nullptr;
        std::vector<int> pos, lvl, cellsByLvl, backList;
        std::vector<uint8_t> back;
        std::vector<float> rv;
        std::vector<Item> items;
        void (*fn)(void*, int, float) = nullptr; void* ctx = nullptr;
        float recvH(int c) const { return back[c] ? rv[c] : h[receiver[c]]; }
        void work(int t) {
            bool prevSerial = true; // the start signal stands for the first item's barrier
            for (size_t k = 0; k < items.size(); ++k) {
                const Item& it = items[k];
                if (!it.par) { if (t == 0) for (int i = it.begin; i < it.end; ++i) fn(ctx, cellsByLvl[i], recvH(cellsByLvl[i])); prevSerial = true; continue; }
                if (prevSerial && k > 0) bar.wait(); // the serial levels before it are done
                for (int i = it.begin + t; i < it.end; i += T) fn(ctx, cellsByLvl[i], recvH(cellsByLvl[i]));
                bar.wait();                            // the whole level is done
                prevSerial = false;
            }
            bar.wait(); // the step is done on every thread before the caller moves on (and replans)
        }
        void loop(int t) {
            uint64_t seen = 0;
            for (;;) {
                { std::unique_lock<std::mutex> l(m); cv.wait(l, [&] { return stepGen != seen; }); seen = stepGen; if (quit) return; }
                work(t);
            }
        }
    };

    // Stream power + talus for `steps` steps (see erosion.h). nbOf(c, k) = the k-th of 8 neighbours or -1.
    template <class NB>
    void Simulate(Domain& d, NB&& nbOf, int steps, const ErosionParams& p, float R, unsigned int seed, float kScale, std::vector<float>* filledOut = nullptr, std::vector<float>* windOut = nullptr, SimProfile* prof = nullptr, int threads = 1, bool fine = false) {
        const size_t cells = d.cells;
        // #633: the talus limit varies (valley walls were planes at exactly one angle: perfect V gorges under every
        // fall). Per cell a factor 0.6-1.4 from value noise at ~260 and ~95 m, and per 30 m height band 0.8-1.25 (the
        // bands' heights wobble by cell): walls come out with steeper and gentler stretches, ledges and benches.
        std::vector<float> talusVar, talusBandOff;
        std::vector<int> wide; // #633: channel cells for the lateral erosion (fine grid), refreshed with each flood
        if (!gBakeOld) {
            talusVar.resize(cells); talusBandOff.resize(cells);
            for (size_t c = 0; c < cells; ++c) {
                const float n = 0.65f * ValueNoise(Vector3Scale(d.dir[c], R / 260.0f), seed ^ 0x7A1u) + 0.35f * ValueNoise(Vector3Scale(d.dir[c], R / 95.0f), seed ^ 0x7A2u);
                talusVar[c] = 0.6f + 0.8f * std::clamp((n - 0.25f) / 0.5f, 0.0f, 1.0f);
                talusBandOff[c] = 30.0f * ValueNoise(Vector3Scale(d.dir[c], R / 700.0f), seed ^ 0x7A3u);
            }
        }
        using Clock = std::chrono::steady_clock;
        auto tick = Clock::now();
        auto lap = [&](int k) { if (!prof) return; const auto t = Clock::now(); prof->s[k] += std::chrono::duration<double>(t - tick).count(); tick = t; }; // #605
        const float stratumOffset = Hash3(1, 2, 3, seed) * 120.0f;
        auto capRock = [&](float hc) { return fmodf(hc + stratumOffset + 1200.0f, 120.0f) < 36.0f; }; // the hard 30% of each 120 m band
        const bool layered = !d.col.empty(); // #587: the rock column's layer at the cell's current height
        auto erodibility = [&](size_t c, float hc) {
            if (layered) return p.k * kScale * Strata_Erodibility(Strata_RockAt(d.col[c], hc, d.strataSeed));
            const float k = d.rockK[c] * kScale;
            if (d.strata[c] <= 0.0f) return k;
            const float bandK = capRock(hc) ? 0.15f : 1.6f;
            return k * (1.0f + (bandK - 1.0f) * d.strata[c]);
        };
        std::vector<float> hf(cells), area(cells), areaM(gBakeOld ? 0 : cells), areaG(gBakeOld ? 0 : cells), recvL(gBakeOld ? 0 : cells);
        std::vector<int> receiver(cells), order;
        std::vector<float> hsDelta; // #1884: the hillslope pass's changes (applied after it, so the scan order doesn't matter)
        // Wind (#370): per cell its downwind and upwind neighbours (latitude bands: trade winds blowing west and toward the
        // equator below ~30 deg, westerlies east and poleward to ~60, polar easterlies beyond), the wind's speed there, and
        // how dry the band is (the subtropical desert belt around 25 deg, a weaker polar one; the wet equator and
        // mid-latitudes barely).
        std::vector<int> downwind, upwind;
        std::vector<float> windSpeed, dryLat, windDelta, windE; std::vector<int> windN, windDrop;
        if (p.wind > 0.0f) {
            downwind.assign(cells, -1); upwind.assign(cells, -1); windSpeed.assign(cells, 0.0f); dryLat.assign(cells, 0.0f); windDelta.assign(cells, 0.0f);
            for (size_t c = 0; c < cells; ++c) {
                const Vector3 up = d.dir[c];
                Vector3 east = Vector3CrossProduct({ 0, 1, 0 }, up);
                if (Vector3LengthSqr(east) < 1e-6f) continue; // the poles: calm
                east = Vector3Normalize(east);
                const Vector3 north = Vector3CrossProduct(up, east);
                const float latDeg = asinf(std::clamp(up.y, -1.0f, 1.0f)) * RAD2DEG, a = fabsf(latDeg), hemi = latDeg >= 0.0f ? 1.0f : -1.0f;
                Vector3 w; float speed;
                if (a < 30.0f) { w = Vector3Add(Vector3Scale(east, -1.0f), Vector3Scale(north, -0.35f * hemi)); speed = 1.0f; }      // trades: west, equatorward
                else if (a < 60.0f) { w = Vector3Add(east, Vector3Scale(north, 0.3f * hemi)); speed = 1.2f; }                         // westerlies: east, poleward
                else { w = Vector3Add(Vector3Scale(east, -1.0f), Vector3Scale(north, -0.2f * hemi)); speed = 0.8f; }                 // polar easterlies
                w = Vector3Normalize(w);
                float bestDown = -2.0f, bestUp = 2.0f;
                for (int k = 0; k < 8; ++k) {
                    const int o = nbOf(c, k);
                    if (o < 0) continue;
                    const float dot = Vector3DotProduct(Vector3Normalize(Vector3Subtract(d.dir[o], up)), w);
                    if (dot > bestDown) { bestDown = dot; downwind[c] = o; }
                    if (dot < bestUp) { bestUp = dot; upwind[c] = o; }
                }
                windSpeed[c] = speed;
                dryLat[c] = expf(-((a - 25.0f) * (a - 25.0f)) / (2.0f * 9.0f * 9.0f)) + 0.4f * expf(-((a - 75.0f) * (a - 75.0f)) / (2.0f * 10.0f * 10.0f));
            }
        }
        std::vector<uint8_t> seen(cells);
        std::unique_ptr<StreamPool> pool; // #605: the stream power on `threads` threads (the coarse pass; the faces run side by side)
        if (threads > 1 && !gBakeOld) pool = std::make_unique<StreamPool>(threads);
        order.reserve(cells);
        using Item = std::pair<float, int>;
        int step = 0;
        RadixHeap rh; // #605: the same pop order as the binary heap (gBakeOld), faster
        auto flood = [&] {
            std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
            rh.reset();
            auto push = [&](float v, int c) { if (gBakeOld) pq.push({ v, c }); else rh.push(RadixHeap::Key(v, c)); };
            auto empty = [&] { return gBakeOld ? pq.empty() : rh.empty(); };
            std::fill(seen.begin(), seen.end(), 0);
            order.clear();
            for (size_t c = 0; c < cells; ++c)
                if (d.outlet[c]) { hf[c] = d.h[c]; receiver[c] = (int)c; seen[c] = 1; push(d.h[c], (int)c); }
            if (empty()) { // no sea: the lowest cell drains everything
                size_t lo = 0; for (size_t c = 1; c < cells; ++c) if (d.h[c] < d.h[lo]) lo = c;
                hf[lo] = d.h[lo]; receiver[lo] = (int)lo; seen[lo] = 1; push(d.h[lo], (int)lo);
            }
            while (!empty()) {
                int c;
                if (gBakeOld) { c = pq.top().second; pq.pop(); } else c = (int)(uint32_t)rh.pop();
                order.push_back(c);
                for (int k = 0; k < 8; ++k) {
                    const int o = nbOf((size_t)c, k);
                    if (o < 0 || seen[o]) continue;
                    seen[o] = 1;
                    hf[o] = std::max(d.h[o], hf[c] + 1e-3f * (1.0f + 4.0f * Hash3(o, step, 7, seed))); // jittered across flats
                    receiver[o] = c;
                    push(hf[o], o);
                }
            }
            for (int c : order) { // the receiver on the filled surface (the flood parent where none is lower)
                if (receiver[c] == c) continue;
                if (p.randomReceivers) {
                    // #1884: one of the lower neighbours, odds by slope, rolled again each flood (seeded: the bake stays
                    // deterministic). Over the floods the flow spreads as over a smooth surface instead of along the
                    // 8 grid directions. A lower neighbour comes earlier in `order` (lower hf), as the steepest one did.
                    int lower[8]; float slope[8], sum = 0.0f; int n = 0;
                    for (int k = 0; k < 8; ++k) {
                        const int o = nbOf((size_t)c, k);
                        if (o < 0 || hf[o] >= hf[c]) continue;
                        const float s = (hf[c] - hf[o]) / (Vector3Distance(d.dir[c], d.dir[o]) * R);
                        slope[n] = s * s; // (squared: the steepest still wins most rolls, so channels stay single threads; GUESS)
                        sum += slope[n]; lower[n++] = o;
                    }
                    if (n == 0 || sum <= 0.0f) continue;
                    float roll = Hash3(c, step, 0x1884, seed) * sum;
                    int pick = n - 1;
                    for (int i = 0; i < n - 1; ++i) { roll -= slope[i]; if (roll < 0.0f) { pick = i; break; } }
                    receiver[c] = lower[pick];
                    continue;
                }
                float best = 0.0f;
                for (int k = 0; k < 8; ++k) {
                    const int o = nbOf((size_t)c, k);
                    if (o < 0 || hf[o] >= hf[c]) continue; // not lower: never steeper than 0 (#605: skips the distance)
                    const float s = (hf[c] - hf[o]) / (Vector3Distance(d.dir[c], d.dir[o]) * R);
                    if (s > best) { best = s; receiver[c] = o; }
                }
            }
            for (size_t c = 0; c < cells; ++c) area[c] = d.cellArea[c];
            for (size_t k = order.size(); k-- > 0;) { const int c = order[k]; if (receiver[c] != c) area[receiver[c]] += area[c]; }
            if (!gBakeOld) // #605: what the stream power needs per cell changes only here; once per flood, not every step
                for (size_t c = 0; c < cells; ++c) {
                    areaM[c] = powf(area[c], p.m);
                    areaG[c] = powf(area[c], 0.25f);
                    recvL[c] = Vector3Distance(d.dir[c], d.dir[receiver[c]]) * R;
                }
        };
        lap(0);
        for (step = 0; step < steps; ++step) {
            if (kScale == 1.0f) gBakeProgress.store(0.3f * (float)step / (float)std::max(1, steps), std::memory_order_relaxed); // the coarse pass
            if (step % std::max(1, p.floodEvery) == 0) {
                flood();
                if (pool) pool->plan(order, receiver, d.outlet, d.h.data());
                if (fine && !gBakeOld) { wide.clear(); for (int c : order) if (area[c] >= 5.0e6f && !d.outlet[c]) wide.push_back(c); } // #633: the channel cells, in flood order
                lap(1);
                // Glacial troughs (#358): ice fills each valley to a thickness growing with its drainage area (~45 m
                // at 1 km^2, at most 200 m); valley walls standing below that ice surface are ground back, so V
                // valleys widen into U troughs. Per flood interval, a share of the wall above the floor goes.
                if (p.glacial > 0.0f)
                    for (int c : order) {
                        const float g = d.glacier[c];
                        if (g <= 0.0f || d.outlet[c]) continue;
                        const float ice = g * std::min(200.0f, 1.4f * powf(area[c], 0.25f));
                        if (ice < 5.0f) continue;
                        const float top = d.h[c] + ice;
                        for (int k = 0; k < 8; ++k) {
                            const int o = nbOf((size_t)c, k);
                            if (o < 0 || d.outlet[o] || d.h[o] <= d.h[c] || d.h[o] >= top) continue;
                            d.h[o] -= 0.15f * g * (d.h[o] - d.h[c]) * (1.0f - (d.h[o] - d.h[c]) / ice);
                        }
                    }
                lap(2);
            }
            // Implicit stream power, receivers first. hr: the receiver's height as the serial pass in `order` reads it.
            auto streamCell = [&](int c, float hr) {
                if (hf[c] > d.h[c] + 1e-3f) { d.h[c] += d.U[c]; return; } // a filled depression: a lake bed, no incision
                const int r = receiver[c];
                const float L = gBakeOld ? Vector3Distance(d.dir[c], d.dir[r]) * R : recvL[c];
                const float ek = erodibility(c, d.h[c]);
                float f = ek * (gBakeOld ? powf(area[c], p.m) : areaM[c]) / L;
                if (const float g = d.glacier[c]; g > 0.0f) {
                    // Glacial: erosion grows more slowly with drainage area (A^0.25, matched to the fluvial rate at
                    // 1 km^2 and p.glacial x it), so it reaches up to the valley heads: cirques, steep headwalls.
                    const float fg = ek * p.glacial * 31.6f * (gBakeOld ? powf(area[c], 0.25f) : areaG[c]) / L;
                    f += (fg - f) * g;
                }
                d.h[c] = std::max((d.h[c] + d.U[c] + f * hr) / (1.0f + f), hr);
            };
            if (pool) pool->stream(streamCell);
            else
                for (int c : order) {
                    const int r = receiver[c];
                    if (r == c || d.outlet[c]) continue;
                    streamCell(c, d.h[r]);
                }
            // #633: lateral erosion on the fine grid. A river of some size wears its banks back as well as its bed, so
            // its valley gets a floor instead of a V: each step the side neighbours of a channel cell standing over it
            // lose a share of the difference, growing with the drainage (from 5 km^2; full at ~80 km^2). Only the
            // direct neighbours, so the floor is at most 3 cells (~150 m) wide; small streams keep their V-notches.
            if (fine && !gBakeOld && step % 2 == 0)
                for (int c : wide) {
                    if (hf[c] > d.h[c] + 1e-3f) continue; // lakes: no banks to wear
                    const float w = 0.03f * std::min(1.0f, 0.25f + 0.1875f * log2f(area[c] / 5.0e6f));
                    for (int k = 0; k < 4; ++k) {
                        const int o = nbOf((size_t)c, k);
                        if (o < 0 || d.outlet[o] || o == receiver[c] || area[o] >= area[c]) continue;
                        const float dh = d.h[o] - d.h[c];
                        if (dh > 0.5f) d.h[o] -= w * dh;
                    }
                }
            // #1884: hillslope diffusion (soil creep): every other step (the talus takes the others), each cell moves
            // toward its 4 side neighbours by D x 2 steps / L^2 of each difference (the 4-neighbour Laplacian), capped
            // at 0.1 a neighbour for stability. Crests round over and V shoulders soften; the sea floor stays put.
            if (p.hillslope > 0.0f && step % 2 == 0) {
                hsDelta.assign(cells, 0.0f);
                for (size_t c = 0; c < cells; ++c) {
                    if (d.outlet[c]) continue;
                    float dh = 0.0f;
                    for (int k = 0; k < 4; ++k) {
                        const int o = nbOf(c, k);
                        if (o < 0) continue;
                        const float L = Vector3Distance(d.dir[c], d.dir[o]) * R;
                        dh += std::min(2.0f * p.hillslope / (L * L), 0.1f) * (d.h[o] - d.h[c]);
                    }
                    // #1965: on the fine pass creep never fills a closed hollow's deep floor (the flood filled it 0.4 m+, half
                    // a bake pond's 0.8 m minimum): it filled ~2/3 of the bake's ponds (seed 42: 609 hollows with it off, 223
                    // on). Rims and shallow pits still round over: a 5 cm guard kept every pit on a slope, and their puddles
                    // never drained (#2018: fluid/slope-drain, fluid/flow).
                    if (fine && dh > 0.0f && hf[c] > d.h[c] + 0.4f) dh = 0.0f;
                    hsDelta[c] = dh;
                }
                for (size_t c = 0; c < cells; ++c) d.h[c] += hsDelta[c];
            }
            lap(3);
            if (step % 2 == 1) // talus: slopes steeper than the rock holds slump (half the excess moves down)
                for (size_t c = 0; c < cells; ++c) {
                    if (d.outlet[c]) continue;
                    const float hard = layered ? Strata_Hardness(Strata_RockAt(d.col[c], d.h[c], d.strataSeed))
                                               : std::max(d.hardness[c], d.strata[c] * (capRock(d.h[c]) ? 1.0f : 0.0f));
                    float maxSlope = p.talusSoft + (p.talusHard - p.talusSoft) * hard;
                    if (!talusVar.empty()) maxSlope *= talusVar[c] * (0.8f + 0.45f * Hash3((int)floorf((d.h[c] + talusBandOff[c]) / 30.0f), 0, 0, seed ^ 0x7A4u));
                    for (int k = 0; k < 8; ++k) {
                        const int o = nbOf(c, k);
                        if (o < 0) continue;
                        const float excess = d.h[c] - d.h[o] - maxSlope * Vector3Distance(d.dir[c], d.dir[o]) * R;
                        if (excess > 0.0f) { const float mv = 0.5f * excess; d.h[c] -= mv; if (!d.outlet[o]) d.h[o] += mv; }
                    }
                }
            lap(4);
            if (p.wind > 0.0f && step % 2 == 0) {
                // Aeolian erosion (#370): on dry land (the latitude belt, away from rivers: large drainage areas are wet) the
                // wind deflates ground, most on faces turned into it (yardang-like streaks along the wind), and carries
                // it to the downwind neighbour, where it settles in the lee of ridges and in basins (dune fields). Moved,
                // not lost (into the sea it is). Applied together after the pass, so the scan order doesn't matter.
                std::fill(windDelta.begin(), windDelta.end(), 0.0f);
                // One cell: how much it loses (e > 0, or false: none), and where it goes: dropAt, else spread over the n
                // cells downwind of it (n = 0: lost). Reads only; the heights change after the whole pass.
                auto windCell = [&](size_t c, float& e, int& n, int& dropAt) {
                    const int dn = downwind[c], un = upwind[c];
                    if (dn < 0 || un < 0 || d.outlet[c] || d.h[c] <= 0.0f) return false;
                    const float wet = std::clamp((logf(std::max(area[c], 1.0f)) - logf(2.0e6f)) / (logf(2.0e7f) - logf(2.0e6f)), 0.0f, 1.0f);
                    const float dry = dryLat[c] * (1.0f - wet);
                    if (dry < 0.02f) return false;
                    const float L = Vector3Distance(d.dir[c], d.dir[un]) * R;
                    const float windward = std::max(0.0f, (d.h[c] - d.h[un]) / L); // rising into the wind
                    const float lee = std::max(0.0f, (d.h[c] - d.h[dn]) / L);       // dropping away downwind: sheltered
                    e = 0.05f * p.wind * dry * windSpeed[c] * (0.3f + 4.0f * std::min(windward, 0.5f)) / (1.0f + 3.0f * lee);
                    // Saltation carries it several cells downwind (one hop left one-cell stripes: each neighbour gave it
                    // straight back the next step). It settles at the first sheltered cell: a drop behind a ridge (the lee)
                    // or a hollow lower than both its up- and downwind neighbours; with none within 8 cells it is spread
                    // along the path. Into the sea it is lost.
                    int at = (int)c; n = 0; dropAt = -1;
                    for (int hop = 0; hop < 8; ++hop) {
                        const int nx = downwind[at];
                        if (nx < 0 || d.outlet[nx]) break;
                        if (area[nx] > 2.0e7f) { n = 0; break; } // a river: it carries the sand away (it filled channels 74 m deep)
                        ++n;
                        const float Lh = Vector3Distance(d.dir[at], d.dir[nx]) * R;
                        const bool leeDrop = d.h[at] - d.h[nx] > 0.12f * Lh;
                        const int after = downwind[nx];
                        const bool hollow = after >= 0 && d.h[nx] < d.h[at] && d.h[nx] < d.h[after];
                        if (hop > 0 && (leeDrop || hollow)) { dropAt = nx; break; }
                        at = nx;
                    }
                    return true;
                };
                // Applied in cell order (the sums' order fixes the result): the path is the downwind chain from c.
                auto windApply = [&](size_t c, float e, int n, int dropAt) {
                    windDelta[c] -= e;
                    if (dropAt >= 0) windDelta[dropAt] += e;
                    else for (int k = 0, at = downwind[c]; k < n; ++k, at = downwind[at]) windDelta[at] += e / n;
                };
                if (threads > 1 && !gBakeOld) { // #605: the cells on threads (contiguous ranges), then applied in order
                    windE.resize(cells); windN.resize(cells); windDrop.resize(cells);
                    std::vector<std::thread> team;
                    const size_t per = (cells + threads - 1) / threads;
                    for (int t = 0; t < threads; ++t)
                        team.emplace_back([&, t] {
                            for (size_t c = per * t, e2 = std::min(cells, per * (t + 1)); c < e2; ++c)
                                if (!windCell(c, windE[c], windN[c], windDrop[c])) windE[c] = -1.0f;
                        });
                    for (auto& th : team) th.join();
                    for (size_t c = 0; c < cells; ++c) if (windE[c] >= 0.0f) windApply(c, windE[c], windN[c], windDrop[c]);
                } else
                    for (size_t c = 0; c < cells; ++c) { float e; int n, dropAt; if (windCell(c, e, n, dropAt)) windApply(c, e, n, dropAt); }
                for (size_t c = 0; c < cells; ++c) d.h[c] += windDelta[c];
                if (windOut) { windOut->resize(cells, 0.0f); for (size_t c = 0; c < cells; ++c) (*windOut)[c] += windDelta[c]; } // the wind's own net change (#370 evidence)
                lap(5);
            }
        }
        if (filledOut) { flood(); *filledOut = hf; } // the depressions left (ponds, #360): filled heights of the result
        lap(1);
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

    // Catmull-Rom weights of the 4 nodes around a cell at fraction t (they sum to 1).
    inline void CatmullRomWeights(float t, float w[4]) {
        const float t2 = t * t, t3 = t2 * t;
        w[0] = 0.5f * (-t3 + 2.0f * t2 - t);
        w[1] = 0.5f * (3.0f * t3 - 5.0f * t2 + 2.0f);
        w[2] = 0.5f * (-3.0f * t3 + 4.0f * t2 + t);
        w[3] = 0.5f * (t3 - t2);
    }

    // Catmull-Rom bicubic on one face (#678), the same layout as FaceBilinear: C1 across cell edges, so the slope no
    // longer jumps every cell (bilinear drew the grid as squares and diagonal bands in every hillshade). It interpolates
    // (each node keeps its baked value: river beds, lake shores and peaks stay where the bake put them; a cubic B-spline
    // would lower one-cell channels by a third and every peak); the price is a small overshoot beside a sharp step (at
    // most ~7% of the step). Nodes past the face's last row repeat it.
    float FaceCatmullRom(const std::vector<int16_t>& g, int f, int side, float u, float v, float span) {
        const float du = 2.0f * span / side;
        const float fi = std::clamp((u + span) / du - 0.5f, 0.0f, (float)(side - 1) - 1e-3f);
        const float fj = std::clamp((v + span) / du - 0.5f, 0.0f, (float)(side - 1) - 1e-3f);
        const int i = (int)fi, j = (int)fj;
        float wu[4], wv[4];
        CatmullRomWeights(fi - i, wu);
        CatmullRomWeights(fj - j, wv);
        int ci[4], cj[4];
        for (int k = 0; k < 4; ++k) { ci[k] = std::clamp(i - 1 + k, 0, side - 1); cj[k] = std::clamp(j - 1 + k, 0, side - 1); }
        const int16_t* base = &g[(size_t)f * side * side];
        float sum = 0.0f;
        for (int r = 0; r < 4; ++r) {
            const int16_t* row = base + (size_t)cj[r] * side;
            sum += wv[r] * (wu[0] * row[ci[0]] + wu[1] * row[ci[1]] + wu[2] * row[ci[2]] + wu[3] * row[ci[3]]);
        }
        return sum * 0.1f;
    }
    bool gBakeBilinear = [] { const char* e = getenv("RAYLIBGAME_BAKE_BILINEAR"); return e && e[0] == '1'; }(); // #678 A/B
    inline float FaceSample(const std::vector<int16_t>& g, int f, int side, float u, float v, float span) {
        return gBakeBilinear ? FaceBilinear(g, f, side, u, v, span) : FaceCatmullRom(g, f, side, u, v, span);
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

    // Before/after/difference of one face (#358 evidence): `beforePath` holds an earlier run's heights (int side, then
    // floats). Writes erosion_compare.png (the 512-cell window with the most change: before | after | difference, 2x)
    // and erosion_compare_face.png (the whole face, difference over the after hillshade). Difference: red = lowered,
    // blue = raised, full colour at 150 m.
    void DumpCompare(const char* beforePath, const float* after, int N, float cell, const Vector3* dirs = nullptr) {
        FILE* fi = fopen(beforePath, "rb");
        if (!fi) { printf("erosion: compare: cannot read %s\n", beforePath); return; }
        int n0 = 0;
        std::vector<float> before((size_t)N * N);
        const bool ok = fread(&n0, sizeof n0, 1, fi) == 1 && n0 == N && fread(before.data(), sizeof(float), before.size(), fi) == before.size();
        fclose(fi);
        if (!ok) { printf("erosion: compare: %s does not match this face (%d vs %d cells)\n", beforePath, n0, N); return; }
        auto shadeAt = [&](const float* h, int i, int j) {
            auto H = [&](int a, int b) { a = std::clamp(a, 0, N - 1); b = std::clamp(b, 0, N - 1); return h[(size_t)b * N + a]; };
            const float gx = (H(i + 1, j) - H(i - 1, j)) / (2 * cell), gy = (H(i, j + 1) - H(i, j - 1)) / (2 * cell);
            const Vector3 nrm = Vector3Normalize({ -gx, 1.0f, -gy });
            return std::clamp(0.35f + 0.75f * Vector3DotProduct(nrm, Vector3Normalize({ -0.6f, 0.7f, -0.4f })), 0.0f, 1.0f);
        };
        auto terrain = [&](const float* h, int i, int j) -> Color {
            const float e = h[(size_t)j * N + i], s = shadeAt(h, i, j);
            const Vector3 c = e <= 0.0f ? Vector3{ 0.15f, 0.3f, 0.6f } : e < 200.0f ? Vector3{ 0.35f, 0.6f, 0.3f } : e < 700.0f ? Vector3{ 0.5f, 0.5f, 0.35f } : Vector3{ 0.85f, 0.85f, 0.9f };
            return { (unsigned char)(255 * c.x * s), (unsigned char)(255 * c.y * s), (unsigned char)(255 * c.z * s), 255 };
        };
        auto diffColor = [&](int i, int j) -> Color {
            const float g = 0.25f + 0.6f * shadeAt(after, i, j); // grey hillshade under the colour
            const float dlt = after[(size_t)j * N + i] - before[(size_t)j * N + i];
            const float t = std::clamp(fabsf(dlt) / 150.0f, 0.0f, 1.0f), a = fabsf(dlt) < 2.0f ? 0.0f : 0.25f + 0.75f * t;
            const Vector3 tint = dlt < 0.0f ? Vector3{ 1.0f, 0.15f, 0.1f } : Vector3{ 0.15f, 0.45f, 1.0f };
            return { (unsigned char)(255 * (g * (1 - a) + tint.x * a)), (unsigned char)(255 * (g * (1 - a) + tint.y * a)), (unsigned char)(255 * (g * (1 - a) + tint.z * a)), 255 };
        };
        // Stats, and the window with the most change (summed |difference| on a 32-cell grid).
        double sum = 0; float maxCutD = 0, maxRaiseD = 0; size_t changed = 0, cutAt = 0;
        const int B = 32, nb = N / B;
        std::vector<double> block((size_t)nb * nb, 0.0);
        for (int j = 0; j < N; ++j)
            for (int i = 0; i < N; ++i) {
                const float dlt = after[(size_t)j * N + i] - before[(size_t)j * N + i];
                sum += fabsf(dlt); changed += fabsf(dlt) >= 2.0f;
                if (-dlt > maxCutD) { maxCutD = -dlt; cutAt = (size_t)j * N + i; }
                maxRaiseD = std::max(maxRaiseD, dlt);
                if (i / B < nb && j / B < nb) block[(size_t)(j / B) * nb + i / B] += fabsf(dlt);
            }
        const int W = std::min(512, N), wb = W / B;
        int bi = 0, bj = 0; double best = -1;
        for (int j = 0; j + wb <= nb; ++j)
            for (int i = 0; i + wb <= nb; ++i) {
                double s = 0;
                for (int b = 0; b < wb; ++b) for (int a = 0; a < wb; ++a) s += block[(size_t)(j + b) * nb + i + a];
                if (s > best) { best = s; bi = i * B; bj = j * B; }
            }
        printf("erosion: compare: %zu of %d cells changed by 2 m or more, mean |change| %.2f m, lowered up to %.0f m, raised up to %.0f m; window at cells (%d, %d) %.0f m wide\n",
               changed, N * N, sum / ((double)N * N), maxCutD, maxRaiseD, bi, bj, W * cell);
        if (dirs) { // where to look in the game: RAYLIBGAME_SPAWN_DIR takes these
            const Vector3 wc = dirs[(size_t)(bj + W / 2) * N + bi + W / 2], dc = dirs[cutAt];
            printf("erosion: compare: window centre dir %.6f,%.6f,%.6f; deepest cut (%.0f m, %.0f m high after) dir %.6f,%.6f,%.6f\n",
                   wc.x, wc.y, wc.z, maxCutD, after[cutAt], dc.x, dc.y, dc.z);
        }
        const int S = 2, gap = 8; // 2x zoom, panels side by side
        Image img = GenImageColor(3 * W * S + 2 * gap, W * S + 28, { 20, 20, 24, 255 });
        for (int j = 0; j < W; ++j)
            for (int i = 0; i < W; ++i) {
                const Color c[3] = { terrain(before.data(), bi + i, bj + j), terrain(after, bi + i, bj + j), diffColor(bi + i, bj + j) };
                for (int p = 0; p < 3; ++p) ImageDrawRectangle(&img, p * (W * S + gap) + i * S, 28 + j * S, S, S, c[p]);
            }
        ImageDrawText(&img, "BEFORE", 8, 6, 20, RAYWHITE);
        ImageDrawText(&img, "AFTER", W * S + gap + 8, 6, 20, RAYWHITE);
        ImageDrawText(&img, "DIFFERENCE  red = lowered  blue = raised  (full at 150 m)", 2 * (W * S + gap) + 8, 6, 20, RAYWHITE);
        ExportImage(img, "erosion_compare.png");
        UnloadImage(img);
        Image face = GenImageColor(N, N, BLACK);
        for (int j = 0; j < N; ++j) for (int i = 0; i < N; ++i) ImageDrawPixel(&face, i, j, diffColor(i, j));
        ImageDrawRectangleLines(&face, { (float)bi, (float)bj, (float)W, (float)W }, 3, YELLOW);
        ExportImage(face, "erosion_compare_face.png");
        UnloadImage(face);
    }
}

uint64_t ErosionParamsKey(const ErosionParams& p) {
    uint64_t h = 1469598103934665603ull;
    auto mix = [&h](const void* v, size_t n) { const unsigned char* b = (const unsigned char*)v; for (size_t i = 0; i < n; ++i) { h ^= b[i]; h *= 1099511628211ull; } };
    const uint32_t ver = kErosionBakeVersion; mix(&ver, sizeof ver);
    for (float f : { p.cell, p.k, p.m, p.uplift, p.talusSoft, p.talusHard, p.glacial, p.glacialLine, p.start, p.wind, p.hillslope }) mix(&f, sizeof f);
    for (int i : { p.steps, p.floodEvery, p.refine, p.refineSteps, p.margin, (int)p.randomReceivers }) mix(&i, sizeof i);
    if (p.landKey) mix(&p.landKey, sizeof p.landKey); // #504: a world's generation parameters (0: defaults, old caches still load)
    { const float span = NoiseSeedSpanUnits(); if (span != 64.0f) mix(&span, sizeof span); } // #532: the noise offsets (an A/B span bakes its own)
    if (!gStrataFolds) { const int nf = 0x591; mix(&nf, sizeof nf); }                                    // #591 A/B: RAYLIBGAME_FOLDS=0 bakes its own
    if (BakeRiversOld()) { const int nr = 0x547; mix(&nr, sizeof nr); } // #547 A/B: no final drainage pass
    if (const char* e = getenv("RAYLIBGAME_CANYON_GRADE"); e && e[0] == '0') { const int ng = 0x544; mix(&ng, sizeof ng); } // #544 A/B: ungraded canyons bake their own
    return h;
}

void PlanetErosion::build(const std::function<float(Vector3, float)>& height, float planetRadius, unsigned int seedIn, const ErosionParams& p, int threads, const StrataFn& strata) {
    gBakeProgress.store(0.0f, std::memory_order_relaxed);
    struct BakingFlag { BakingFlag() { gBaking.store(true); ++gBakes; } ~BakingFlag() { gBaking.store(false); } } baking; // the loading screen tells a bake from a cache read (#505)
    paramsKey = ErosionParamsKey(p);
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
    auto tp = std::chrono::steady_clock::now();
    auto pt = [&tp]() { const auto t = std::chrono::steady_clock::now(); const double s = std::chrono::duration<double>(t - tp).count(); tp = t; return s; }; // #605: per-pass seconds
    char pbuf[512];
    SimProfile coarseProf, faceProf[6];
    double faceT[6][4] = {}; // per face: sampling, rock, simulation, the rest (drift, ponds, dumps)
    const double tGrid = pt();
    std::vector<float> h0(d.cells);
    {
        std::vector<std::thread> pool;
        for (int t = 0; t < T; ++t) pool.emplace_back([&, t] { for (size_t c = t; c < d.cells; c += T) h0[c] = height(d.dir[c], p.cell); });
        for (auto& th : pool) th.join();
    }
    float maxH0 = 1.0f;
    for (float v : h0) maxH0 = std::max(maxH0, v);
    const double tSample = pt();
    SetRock(d, h0, maxH0, R, seed, p, &strata, T);
    const double tRock = pt();
    d.outlet.resize(d.cells); d.h.resize(d.cells);
    for (size_t c = 0; c < d.cells; ++c) { d.outlet[c] = h0[c] <= 0.0f; d.h[c] = h0[c] > 0.0f ? h0[c] * p.start : h0[c]; } // start low (uplift mode) or from the input (carve mode)
    Simulate(d, [&](size_t c, int k) { return nb[c * 8 + k]; }, p.steps, p, R, seed, 1.0f, nullptr, nullptr, &coarseProf, T);
    const double tSim = pt();
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
    heightScale = scale; // #587: the strata are laid in simulation heights (materials divide final heights by this)
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
    std::vector<Pond> facePonds[6]; // per face (faces run in parallel)
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
            auto ft = std::chrono::steady_clock::now();
            auto flap = [&](int k) { const auto t = std::chrono::steady_clock::now(); faceT[f][k] += std::chrono::duration<double>(t - ft).count(); ft = t; };
            std::vector<float> fh0(fd.cells), filled, windAcc; // windAcc: the wind's net change on the dump face (#370)
            {
                std::vector<std::thread> pool;
                for (int t = 0; t < sampleThreads; ++t)
                    pool.emplace_back([&, t] { for (size_t c = t; c < fd.cells; c += sampleThreads) fh0[c] = height(fd.dir[c], fineMetres); });
                for (auto& th : pool) th.join();
            }
            flap(0);
            SetRock(fd, fh0, maxH0, R, seed, p, &strata, sampleThreads);
            flap(1);
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
            }, p.refineSteps, p, R, seed, kScale, &filled, DumpOn() && f == dumpFace ? &windAcc : nullptr, &faceProf[f], sampleThreads, true); // #605: the face's share of the threads; fine (#633)
            flap(2);
            {   // Hold the fine grid to the coarse one at large scales (#360): the refinement's own uplift and erosion drift
                // it by tens of metres (70 m in seed 42's ranges), so the coarse LOD, and the rivers and lakes built on it,
                // disagreed with the near ground (a lake 130 m under it). Fine = coarse + the high-pass part of
                // (fine - coarse): the drift is box-blurred over ~6 fine cells and taken out.
                std::vector<float> drift(fd.cells, 0.0f), tmp(fd.cells, 0.0f);
                for (size_t c = 0; c < fd.cells; ++c) drift[c] = fh0[c] > 0.0f ? fd.h[c] - coarseUnscaledAt(fd.dir[c]) : 0.0f;
                const int rad = 3 * p.refine / 2;
                for (int pass = 0; pass < 2; ++pass) { // separable box blur: rows, then columns
                    const std::vector<float>& src = pass == 0 ? drift : tmp;
                    std::vector<float>& dst = pass == 0 ? tmp : drift;
                    std::vector<float> out(fd.cells);
                    for (int j = 0; j < FE; ++j)
                        for (int i = 0; i < FE; ++i) {
                            float s = 0.0f; int cnt = 0;
                            for (int k = -rad; k <= rad; ++k) {
                                const int ii = pass == 0 ? i + k : i, jj = pass == 0 ? j : j + k;
                                if (ii < 0 || jj < 0 || ii >= FE || jj >= FE) continue;
                                s += src[(size_t)jj * FE + ii]; ++cnt;
                            }
                            out[(size_t)j * FE + i] = s / cnt;
                        }
                    dst = std::move(out);
                }
                for (size_t c = 0; c < fd.cells; ++c)
                    if (fh0[c] > 0.0f && !fd.outlet[c]) { fd.h[c] -= drift[c]; filled[c] -= drift[c]; }
            }
            for (size_t c = 0; c < fd.cells; ++c)
                fine[(size_t)f * FE * FE + c] = Q(fh0[c] <= 0.0f ? fh0[c] : std::max(0.5f, fd.h[c] * scale));
            // Ponds (#360): the closed depressions the refinement left in this face's own cells (not its margin, which
            // is another face's): components of cells the flood filled by more than 0.4 m, at most 12 cells, 0.8 m+ deep.
            {
                const int lo = fm, hi = fm + nf;
                auto wet = [&](size_t q) { return !fd.outlet[q] && fh0[q] > 0.0f && (filled[q] - fd.h[q]) * scale > 0.4f; };
                std::vector<uint8_t> seenP(fd.cells, 0);
                std::vector<size_t> stack, members;
                for (int j = lo; j < hi; ++j)
                    for (int i = lo; i < hi; ++i) {
                        const size_t c0 = (size_t)j * FE + i;
                        if (seenP[c0] || !wet(c0)) continue;
                        members.clear(); stack.assign(1, c0); seenP[c0] = 1;
                        bool edge = false;
                        while (!stack.empty()) {
                            const size_t q = stack.back(); stack.pop_back();
                            members.push_back(q);
                            const int qi = (int)(q % FE), qj = (int)(q / FE);
                            for (int k = 0; k < 8; ++k) {
                                const int oi = qi + di[k], oj = qj + dj[k];
                                if (oi < 0 || oj < 0 || oi >= FE || oj >= FE) continue;
                                const size_t o = (size_t)oj * FE + oi;
                                if (seenP[o] || !wet(o)) continue;
                                if (oi < lo || oj < lo || oi >= hi || oj >= hi) { edge = true; continue; }
                                seenP[o] = 1; stack.push_back(o);
                            }
                        }
                        if (edge || members.size() > 12) continue;
                        float level = -1e9f, low = 1e9f; Vector3 sum{ 0, 0, 0 };
                        for (size_t q : members) { level = std::max(level, filled[q]); low = std::min(low, fd.h[q]); sum = Vector3Add(sum, fd.dir[q]); }
                        level *= scale; low *= scale;
                        if (level - low < 0.8f) continue;
                        const float r = std::clamp(sqrtf((float)members.size() * fineMetres * fineMetres / PI) + 6.0f, 12.0f, 60.0f);
                        facePonds[f].push_back({ Vector3Normalize(sum), r, level, level - low });
                    }
            }
            if (DumpOn() && f == dumpFace) {
                std::vector<float> out(fd.cells);
                for (size_t c = 0; c < fd.cells; ++c) out[c] = fine[(size_t)f * FE * FE + c] * 0.1f;
                DumpFace("erosion_fine.png", out.data(), FE, fineMetres);
                // The wind's own work on this face (#370), apart from the drainage changes it sets off: erosion_wind.png,
                // red = deflated, blue = deposited (full colour at 10 m), over the grey hillshade of the result.
                if (!windAcc.empty()) {
                    Image img = GenImageColor(FE, FE, BLACK);
                    double lost = 0.0, laid = 0.0; float maxLost = 0.0f, maxLaid = 0.0f;
                    auto H = [&](int i, int j) { i = std::clamp(i, 0, FE - 1); j = std::clamp(j, 0, FE - 1); return out[(size_t)j * FE + i]; };
                    for (int j = 0; j < FE; ++j)
                        for (int i = 0; i < FE; ++i) {
                            const float w = windAcc[(size_t)j * FE + i];
                            if (w < 0.0f) { lost -= w; maxLost = std::max(maxLost, -w); } else { laid += w; maxLaid = std::max(maxLaid, w); }
                            const float gx = (H(i + 1, j) - H(i - 1, j)) / (2 * fineMetres), gy = (H(i, j + 1) - H(i, j - 1)) / (2 * fineMetres);
                            const float g = 0.25f + 0.6f * std::clamp(0.35f + 0.75f * Vector3DotProduct(Vector3Normalize({ -gx, 1.0f, -gy }), Vector3Normalize({ -0.6f, 0.7f, -0.4f })), 0.0f, 1.0f);
                            const float t = std::clamp(fabsf(w) / 10.0f, 0.0f, 1.0f), a = fabsf(w) < 0.2f ? 0.0f : 0.25f + 0.75f * t;
                            const Vector3 tint = w < 0.0f ? Vector3{ 1.0f, 0.2f, 0.1f } : Vector3{ 0.2f, 0.5f, 1.0f };
                            const float sea = out[(size_t)j * FE + i] <= 0.0f ? 0.6f : 1.0f;
                            ImageDrawPixel(&img, i, j, { (unsigned char)(255 * sea * (g * (1 - a) + tint.x * a)), (unsigned char)(255 * sea * (g * (1 - a) + tint.y * a)), (unsigned char)(255 * (g * (1 - a) + tint.z * a)), 255 });
                        }
                    ExportImage(img, "erosion_wind.png");
                    UnloadImage(img);
                    LOG_INFO("Erosion", "wind on the dump face: %.0f m3/m2-cells deflated (deepest %.1f m), %.0f laid down (thickest %.1f m)\n", lost, maxLost, laid, maxLaid);
                }
                // A/B evidence (#358): RAYLIBGAME_EROSION_RAW=path keeps this face's heights; a later run with
                // RAYLIBGAME_EROSION_COMPARE=path draws that run against this one (before | after | difference).
                if (const char* rp = getenv("RAYLIBGAME_EROSION_RAW")) {
                    if (FILE* fo = fopen(rp, "wb")) { fwrite(&FE, sizeof FE, 1, fo); fwrite(out.data(), sizeof(float), out.size(), fo); fclose(fo); }
                }
                if (const char* cp = getenv("RAYLIBGAME_EROSION_COMPARE")) DumpCompare(cp, out.data(), FE, fineMetres, fd.dir.data());
            }
            flap(3);
        };
        const int parallelFaces = kParallelFaces(T);
        for (int f0 = 0; f0 < 6; f0 += parallelFaces) {
            std::vector<std::thread> pool;
            for (int f = f0; f < std::min(6, f0 + parallelFaces); ++f) pool.emplace_back(refineFace, f, std::max(1, T / parallelFaces));
            for (auto& th : pool) th.join();
            gBakeProgress.store(0.3f + 0.7f * (float)std::min(6, f0 + parallelFaces) / 6.0f, std::memory_order_relaxed);
        }
    }
    pondList.clear();
    for (const auto& fp2 : facePonds) pondList.insert(pondList.end(), fp2.begin(), fp2.end()); // face order: deterministic
    refineSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t1r).count();
    flow.clear(); filled.clear(); routeSeconds = 0.0f;
    if (fe > 0 && !BakeRiversOld()) routeFinal(T); // #547: drainage and channels on the final terrain
    computeHash();
    buildSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
    {   // #605: where the time went (wall seconds; the faces' parts are summed over the faces = thread seconds)
        SimProfile fs; double ft4[4] = {};
        for (int f = 0; f < 6; ++f) { fs.add(faceProf[f]); for (int k = 0; k < 4; ++k) ft4[k] += faceT[f][k]; }
        std::string s;
        snprintf(pbuf, sizeof pbuf, "coarse %.1f s (grid %.2f, input sampling %.2f on %d threads, rock %.2f, simulate %.2f [", coarseSeconds, tGrid, tSample, T, tRock, tSim); s += pbuf;
        for (int k = 0; k < 6; ++k) { snprintf(pbuf, sizeof pbuf, "%s%s %.2f", k ? ", " : "", kSimParts[k], coarseProf.s[k]); s += pbuf; }
        snprintf(pbuf, sizeof pbuf, "]); faces %.1f s wall, %d at a time; summed over faces: sampling %.2f, rock %.2f, simulate %.2f [", refineSeconds, kParallelFaces(T), ft4[0], ft4[1], ft4[2]); s += pbuf;
        for (int k = 0; k < 6; ++k) { snprintf(pbuf, sizeof pbuf, "%s%s %.2f", k ? ", " : "", kSimParts[k], fs.s[k]); s += pbuf; }
        snprintf(pbuf, sizeof pbuf, "], drift+ponds %.2f; final drainage %.1f s; total %.1f s", ft4[3], routeSeconds, buildSeconds); s += pbuf;
        profile = s;
    }
}

bool BakeRiversOld() { static const bool old = getenv("RAYLIBGAME_RIVERS_OLD") && getenv("RAYLIBGAME_RIVERS_OLD")[0] == '1'; return old; }

int PlanetErosion::flowCellOf(Vector3 dv) const {
    const int f = FaceOf(dv);
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float s = Vector3DotProduct(dv, nrm), duf = 2.0f / nf;
    const int i = std::clamp((int)floorf((Vector3DotProduct(dv, t1) / s + 1.0f) / duf), 0, nf - 1);
    const int j = std::clamp((int)floorf((Vector3DotProduct(dv, t2) / s + 1.0f) / duf), 0, nf - 1);
    return (f * nf + j) * nf + i;
}

Vector3 PlanetErosion::flowDir(int g) const {
    const int f = g / (nf * nf), r = g % (nf * nf), j = r / nf, i = r % nf;
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float duf = 2.0f / nf, u = -1.0f + (i + 0.5f) * duf, v = -1.0f + (j + 0.5f) * duf;
    return Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v))));
}

void PlanetErosion::flowNeighbours(int g, int out[8]) const {
    static const int di[8] = { 1, -1, 0, 0, 1, 1, -1, -1 }, dj[8] = { 0, 0, 1, -1, 1, -1, 1, -1 };
    const int f = g / (nf * nf), r = g % (nf * nf), j = r / nf, i = r % nf;
    for (int k = 0; k < 8; ++k) {
        const int ii = i + di[k], jj = j + dj[k];
        if (ii >= 0 && jj >= 0 && ii < nf && jj < nf) { out[k] = (f * nf + jj) * nf + ii; continue; }
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        const float duf = 2.0f / nf, u = -1.0f + (ii + 0.5f) * duf, v = -1.0f + (jj + 0.5f) * duf;
        const int o = flowCellOf(Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v)))));
        out[k] = o == g ? -1 : o;
    }
}

float PlanetErosion::flowHeight(int g) const {
    const int f = g / (nf * nf), r = g % (nf * nf), j = r / nf, i = r % nf;
    return fine[((size_t)f * fe + j + fm) * fe + i + fm] * 0.1f;
}

int PlanetErosion::flowReceiver(int g) const {
    const int k = flow[g] & 15;
    if (k >= 8) return -1;
    int nb[8];
    flowNeighbours(g, nb);
    return nb[k];
}

float PlanetErosion::flowArea(int g) const {
    const int cls = flow[g] >> 4;
    return cls == 0 ? 0.0f : kFlowArea0 * ldexpf(1.41421356f, cls - 1);
}

float PlanetErosion::IncisionDepth(float area) {
    // Light: 0.5 m where a stream starts (0.1 km^2), ~2 m at 10 km^2, at most 4 m.
    return area < kFlowArea0 ? 0.0f : std::min(4.0f, 0.5f * powf(area / kFlowArea0, 0.3f));
}

void PlanetErosion::routeFinal(int threads) {
    const auto t0 = std::chrono::steady_clock::now();
    const int G = 6 * nf * nf;
    std::vector<float> h(G), hf(G), area(G);
    std::vector<int> recv(G, -1), order;
    order.reserve(G);
    for (int g = 0; g < G; ++g) h[g] = flowHeight(g);
    // Neighbour tables would be 680 MB at this size: computed per use (across a face edge by direction).
    auto nbs = [&](int g, int out[8]) { flowNeighbours(g, out); };
    std::vector<uint8_t> seen(G, 0), outlet(G, 0);
    for (int g = 0; g < G; ++g) {
        // The sea, and the polar ice shelves (frozen sea a couple of metres up), as the rivers' network has it.
        const Vector3 d = flowDir(g);
        if (h[g] <= 0.0f || (fabsf(d.y) > 0.86f && h[g] < 3.0f)) outlet[g] = 1;
    }
    const float duf = 2.0f / nf, side = R * duf;
    std::vector<int8_t> dirK(G, 15);
    // One routing: priority flood over 8 neighbours (the bake's own connectivity: its valleys run diagonally too), steepest
    // descent on the filled surface - with preferD4, the steepest of the 4 side neighbours when one is lower, else of all
    // 8 - and the drainage area, downstream in reverse flood order.
    auto route = [&](bool preferD4) {
        RadixHeap heap;
        std::fill(seen.begin(), seen.end(), 0);
        std::fill(recv.begin(), recv.end(), -1);
        order.clear();
        for (int g = 0; g < G; ++g) if (outlet[g]) { hf[g] = h[g]; seen[g] = 1; heap.push(RadixHeap::Key(h[g], g)); }
        int nb[8];
        while (!heap.empty()) {
            const int c = (int)(uint32_t)heap.pop();
            order.push_back(c);
            nbs(c, nb);
            for (int k = 0; k < 8; ++k) {
                const int o = nb[k];
                if (o < 0 || seen[o]) continue;
                seen[o] = 1;
                hf[o] = std::max(h[o], hf[c] + 1e-3f * (1.0f + 4.0f * Hash3(o, 0, 13, seed))); // jittered across flats
                recv[o] = c;
                heap.push(RadixHeap::Key(hf[o], o));
            }
        }
        std::vector<std::thread> pool;
        const int T = std::max(1, threads);
        for (int t = 0; t < T; ++t)
            pool.emplace_back([&, t] {
                int nb2[8];
                for (int c = t; c < G; c += T) {
                    dirK[c] = 15;
                    if (outlet[c] || recv[c] < 0) continue;
                    nbs(c, nb2);
                    const Vector3 dc = flowDir(c);
                    float best = 0.0f; int bestK = -1;
                    for (int pass = preferD4 ? 0 : 1; pass < 2 && bestK < 0; ++pass)
                        for (int k = pass == 0 ? 0 : 0; k < (pass == 0 ? 4 : 8); ++k) {
                            const int o = nb2[k];
                            if (o < 0 || hf[o] >= hf[c]) continue;
                            const float s = (hf[c] - hf[o]) / Vector3Distance(dc, flowDir(o));
                            if (s > best) { best = s; bestK = k; }
                        }
                    if (bestK < 0) for (int k = 0; k < 8; ++k) if (nb2[k] == recv[c]) { bestK = k; break; }
                    dirK[c] = (int8_t)bestK;
                    recv[c] = nb2[bestK];
                }
            });
        for (auto& th : pool) th.join();
        for (int g = 0; g < G; ++g) {
            const int r = g % (nf * nf), j = r / nf, i = r % nf;
            const float u = -1.0f + (i + 0.5f) * duf, v = -1.0f + (j + 0.5f) * duf;
            area[g] = side * side / powf(1.0f + u * u + v * v, 1.5f);
        }
        for (size_t k = order.size(); k-- > 0;) { const int c = order[k]; if (!outlet[c] && recv[c] >= 0) area[recv[c]] += area[c]; }
    };
    route(false);
    // Diagonal channel links (#547): the line between two diagonal cells crosses the corner of the two cells beside them,
    // where the bilinear ground stands at the four cells' mean - up to half a V valley's side over the channel, which the
    // river then cut as a slot. The lower of those two side cells is worn down to the link's middle height, so the
    // channel steps through it on 4-neighbour links (the second routing prefers those) and the bilinear ground follows.
    const std::vector<float> h0 = h;
    {
        static const int kSideA[8] = { -1, -1, -1, -1, 0, 0, 1, 1 }, kSideB[8] = { -1, -1, -1, -1, 2, 3, 2, 3 }; // NE: E, N; SE: E, S; NW: W, N; SW: W, S
        int nb[8];
        for (int c : order) {
            const int k = dirK[c];
            if (k < 4 || k > 7 || area[c] < kFlowArea0 || hf[c] - h[c] > 0.05f) continue;
            nbs(c, nb);
            const int r = nb[k], a = nb[kSideA[k]], b = nb[kSideB[k]];
            if (r < 0 || a < 0 || b < 0) continue;
            const int s = h[a] <= h[b] ? a : b;
            if (outlet[s]) continue;
            h[s] = std::min(h[s], 0.5f * (h[c] + h[r]));
        }
    }
    route(true);
    // Classes, the filled cells, and the incision (channel cells off the filled surface, by their drainage).
    flow.assign(G, 0);
    filled.clear();
    std::vector<float> cut(G, 0.0f);
    for (int g = 0; g < G; ++g) {
        int cls = 0;
        if (area[g] >= kFlowArea0) cls = std::min(15, 1 + (int)floorf(log2f(area[g] / kFlowArea0)));
        flow[g] = (uint8_t)((outlet[g] ? 15 : (dirK[g] & 15)) | (cls << 4));
        if (outlet[g]) continue;
        cut[g] = h0[g] - h[g]; // the corners worn above
        if (hf[g] - h[g] > 0.05f) filled.push_back({ g, hf[g] });
        else if (cls > 0) cut[g] += IncisionDepth(flowArea(g));
    }
    // Into the stored heights: each face's own cells, and its margin cells from the cell of the face that owns them.
    for (int f = 0; f < 6; ++f) {
        Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
        const float span = 1.0f + fm * duf;
        for (int j = 0; j < fe; ++j)
            for (int i = 0; i < fe; ++i) {
                int g;
                if (i >= fm && j >= fm && i < fm + nf && j < fm + nf) g = (f * nf + (j - fm)) * nf + (i - fm);
                else {
                    const float u = -span + (i + 0.5f) * duf, v = -span + (j + 0.5f) * duf;
                    g = flowCellOf(Vector3Normalize(Vector3Add(nrm, Vector3Add(Vector3Scale(t1, u), Vector3Scale(t2, v)))));
                }
                if (cut[g] > 0.0f) {
                    int16_t& q = fine[((size_t)f * fe + j) * fe + i];
                    q = Q(std::max(0.5f, q * 0.1f - cut[g]));
                }
            }
    }
    routeSeconds = std::chrono::duration<float>(std::chrono::steady_clock::now() - t0).count();
}

float PlanetErosion::coarseAt(Vector3 dv) const {
    const int f = FaceOf(dv);
    Vector3 nrm, t1, t2; FaceAxes(f, nrm, t1, t2);
    const float s = Vector3DotProduct(dv, nrm);
    return FaceSample(coarse, f, n, Vector3DotProduct(dv, t1) / s, Vector3DotProduct(dv, t2) / s, 1.0f);
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
        if (std::max(fabsf(u), fabsf(v)) * inner < 0.999f) return FaceSample(fine, f, fe, u, v, span);
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
        sum += w * FaceSample(fine, f, fe, u, v, span);
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
    for (const Pond& pd : pondList) for (float v : { pd.dir.x, pd.dir.y, pd.dir.z, pd.r, pd.level, pd.depth }) { uint32_t b; memcpy(&b, &v, 4); hv = (hv ^ b) * 1099511628211ull; }
    for (uint8_t q : flow) hv = (hv ^ q) * 1099511628211ull; // #547
    for (const auto& fc : filled) { uint32_t b; memcpy(&b, &fc.second, 4); hv = ((hv ^ (uint32_t)fc.first) * 1099511628211ull ^ b) * 1099511628211ull; }
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
    // ERO6 (#362): deflated as kPieces pieces in parallel. One deflate of the ~50 MB took 40 s after a 60-80 s bake.
    std::vector<unsigned char*> comp(kPieces, nullptr);
    std::vector<int> compSize(kPieces, 0), rawLen(kPieces, 0);
    {
        const size_t per = (raw.size() + kPieces - 1) / kPieces;
        std::vector<std::thread> pool;
        for (int k = 0; k < kPieces; ++k) {
            const size_t a = std::min(raw.size(), per * k), b = std::min(raw.size(), per * (k + 1));
            rawLen[k] = (int)(b - a);
            pool.emplace_back([&, k, a] { if (rawLen[k] > 0) comp[k] = CompressData(raw.data() + a, rawLen[k], &compSize[k]); });
        }
        for (auto& t : pool) t.join();
    }
    auto freeAll = [&] { for (unsigned char* c : comp) if (c) MemFree(c); };
    for (int k = 0; k < kPieces; ++k) if (rawLen[k] > 0 && !comp[k]) { freeAll(); return false; }
    const std::string tmp = path + ".tmp"; // written whole, then renamed: a crash mid-save leaves the old cache (#355 review)
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) { freeAll(); return false; }
    const char magic[4] = { 'E', 'R', 'O', '8' }; // ERO7 (#587): + heightScale after the hash; ERO8 (#547): + the final drainage
    fwrite(magic, 1, 4, f);
    fwrite(&paramsKey, sizeof(paramsKey), 1, f);
    fwrite(&n, sizeof(n), 1, f); fwrite(&nf, sizeof(nf), 1, f); fwrite(&fm, sizeof(fm), 1, f); fwrite(&fineMetres, sizeof(fineMetres), 1, f);
    fwrite(&R, sizeof(R), 1, f); fwrite(&seed, sizeof(seed), 1, f); fwrite(&hashValue, sizeof(hashValue), 1, f); fwrite(&heightScale, sizeof(heightScale), 1, f);
    const int pieces = kPieces;
    fwrite(&pieces, sizeof(pieces), 1, f);
    for (int k = 0; k < kPieces; ++k) {
        fwrite(&rawLen[k], sizeof(int), 1, f); fwrite(&compSize[k], sizeof(int), 1, f);
        if (compSize[k] > 0) fwrite(comp[k], 1, compSize[k], f);
    }
    const int np = (int)pondList.size();
    fwrite(&np, sizeof(np), 1, f);
    if (np) fwrite(pondList.data(), sizeof(Pond), np, f);
    {   // ERO8 (#547): the final drainage - the flow bytes (deflated in pieces) and the filled cells
        const int flowLen = (int)flow.size();
        fwrite(&flowLen, sizeof(flowLen), 1, f);
        std::vector<unsigned char*> fc(kPieces, nullptr);
        std::vector<int> fcs(kPieces, 0), flen(kPieces, 0);
        const size_t per = (flow.size() + kPieces - 1) / kPieces;
        std::vector<std::thread> pool;
        for (int k = 0; k < kPieces; ++k) {
            const size_t a = std::min(flow.size(), per * k), b = std::min(flow.size(), per * (k + 1));
            flen[k] = (int)(b - a);
            pool.emplace_back([&, k, a] { if (flen[k] > 0) fc[k] = CompressData(flow.data() + a, flen[k], &fcs[k]); });
        }
        for (auto& t : pool) t.join();
        for (int k = 0; k < kPieces; ++k) {
            const int cs = fc[k] ? fcs[k] : 0, rl = fc[k] ? flen[k] : 0;
            fwrite(&rl, sizeof(int), 1, f); fwrite(&cs, sizeof(int), 1, f);
            if (cs > 0) fwrite(fc[k], 1, cs, f);
            if (fc[k]) MemFree(fc[k]);
        }
        const int nfill = (int)filled.size();
        fwrite(&nfill, sizeof(nfill), 1, f);
        for (const auto& fl : filled) { fwrite(&fl.first, sizeof(int), 1, f); fwrite(&fl.second, sizeof(float), 1, f); }
    }
    const bool written = !ferror(f);
    fclose(f);
    freeAll();
    if (!written) { remove(tmp.c_str()); return false; }
    remove(path.c_str()); // rename does not replace an existing file on Windows
    return rename(tmp.c_str(), path.c_str()) == 0;
}

bool PlanetErosion::load(const std::string& path, float planetRadius, unsigned int wantSeed, uint64_t wantKey) {
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    char magic[4]; int N = 0, NF = 0, FM = 0, pieces = 0; float fmet = 0, r = 0, hs = 1.0f; unsigned int s = 0; uint64_t hv = 0, key = 0;
    bool ok = fread(magic, 1, 4, f) == 4 && memcmp(magic, "ERO8", 4) == 0 && fread(&key, sizeof(key), 1, f) == 1 && key == wantKey && fread(&N, sizeof(N), 1, f) == 1 && fread(&NF, sizeof(NF), 1, f) == 1 &&
              fread(&FM, sizeof(FM), 1, f) == 1 && fread(&fmet, sizeof(fmet), 1, f) == 1 && fread(&r, sizeof(r), 1, f) == 1 && fread(&s, sizeof(s), 1, f) == 1 &&
              fread(&hv, sizeof(hv), 1, f) == 1 && fread(&hs, sizeof(hs), 1, f) == 1 && fread(&pieces, sizeof(pieces), 1, f) == 1 &&
              N > 0 && N < 4096 && NF >= 0 && NF < 16384 && FM >= 0 && FM < 1024 && pieces > 0 && pieces <= 256 && s == wantSeed && fabsf(r - planetRadius) < 1.0f;
    std::vector<std::vector<unsigned char>> comp(ok ? pieces : 0);
    std::vector<int> rawLen(comp.size(), 0);
    for (size_t k = 0; ok && k < comp.size(); ++k) {
        int cs = 0;
        ok = fread(&rawLen[k], sizeof(int), 1, f) == 1 && fread(&cs, sizeof(int), 1, f) == 1 && rawLen[k] >= 0 && cs >= 0;
        if (ok) { comp[k].resize(cs); ok = cs == 0 || fread(comp[k].data(), 1, cs, f) == (size_t)cs; }
    }
    std::vector<Pond> pl;
    if (ok) { int np = 0; ok = fread(&np, sizeof(np), 1, f) == 1 && np >= 0 && np < 10000000; if (ok && np) { pl.resize(np); ok = fread(pl.data(), sizeof(Pond), np, f) == (size_t)np; } }
    std::vector<uint8_t> fl;                       // ERO8 (#547): the final drainage
    std::vector<std::pair<int, float>> fill;
    if (ok) {
        int flowLen = 0;
        ok = fread(&flowLen, sizeof(flowLen), 1, f) == 1 && flowLen >= 0 && flowLen < 400000000;
        if (ok) fl.resize(flowLen);
        size_t at = 0;
        for (int k = 0; ok && k < kPieces; ++k) {
            int rl = 0, cs = 0;
            ok = fread(&rl, sizeof(int), 1, f) == 1 && fread(&cs, sizeof(int), 1, f) == 1 && rl >= 0 && cs >= 0 && at + (size_t)rl <= fl.size();
            if (!ok || cs == 0) continue;
            std::vector<unsigned char> buf(cs);
            ok = fread(buf.data(), 1, cs, f) == (size_t)cs;
            int got = 0;
            unsigned char* d = ok ? DecompressData(buf.data(), cs, &got) : nullptr;
            ok = ok && d && got == rl;
            if (ok) memcpy(fl.data() + at, d, (size_t)rl);
            if (d) MemFree(d);
            at += (size_t)rl;
        }
        int nfill = 0;
        ok = ok && fread(&nfill, sizeof(nfill), 1, f) == 1 && nfill >= 0 && nfill <= flowLen;
        if (ok) fill.resize(nfill);
        for (int k = 0; ok && k < nfill; ++k) ok = fread(&fill[k].first, sizeof(int), 1, f) == 1 && fread(&fill[k].second, sizeof(float), 1, f) == 1;
    }
    fclose(f);
    if (!ok) return false;
    flow = std::move(fl); filled = std::move(fill);
    const size_t nc = (size_t)6 * N * N, fe2 = NF > 0 ? (size_t)NF + 2 * FM : 0, nfine = (size_t)6 * fe2 * fe2;
    size_t rawSize = 0;
    for (int len : rawLen) rawSize += (size_t)len;
    if (rawSize != (nc + nfine) * 2) return false;
    std::vector<unsigned char> rawBuf(rawSize);
    {   // the pieces in parallel (#362)
        std::atomic<bool> good{ true };
        std::vector<std::thread> pool;
        size_t at = 0;
        for (size_t k = 0; k < comp.size(); ++k) {
            pool.emplace_back([&, k, at] {
                if (rawLen[k] == 0) return;
                int got = 0;
                unsigned char* d = DecompressData(comp[k].data(), (int)comp[k].size(), &got);
                if (!d || got != rawLen[k]) good = false; else memcpy(rawBuf.data() + at, d, (size_t)got);
                if (d) MemFree(d);
            });
            at += (size_t)rawLen[k];
        }
        for (auto& t : pool) t.join();
        if (!good) return false;
    }
    const unsigned char* raw = rawBuf.data();
    coarse.resize(nc); fine.resize(nfine);
    const size_t total = nc + nfine;
    int16_t prev = 0;
    for (size_t k = 0; k < total; ++k) {
        prev = (int16_t)(prev + (int16_t)(raw[k] | (raw[total + k] << 8)));
        (k < nc ? coarse[k] : fine[k - nc]) = prev;
    }
    n = N; nf = NF; fm = FM; fe = (int)fe2; fineMetres = fmet; R = r; seed = s; paramsKey = key; heightScale = hs;
    pondList = std::move(pl);
    computeHash();
    if (hashValue != hv) { n = 0; coarse.clear(); fine.clear(); fe = 0; return false; } // damaged
    return true;
}

PlanetErosion::Slope PlanetErosion::slopeAt(Vector3 dv, float footprint) const {
    Slope s{ heightAt(dv, footprint), { 0, 0, 0 } };
    if (n <= 0) return s;
    const bool fineGrid = fe > 0 && footprint < 2.0f * fineMetres;
    const float step = fineGrid ? fineMetres : 2.0f * R / n; // one cell (the face-centre size; near enough elsewhere)
    const Vector3 a = Vector3Normalize(Vector3CrossProduct(dv, fabsf(dv.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 b = Vector3CrossProduct(dv, a);
    const float k = step / R;
    auto at = [&](Vector3 t, float sgn) { return heightAt(Vector3Normalize(Vector3Add(dv, Vector3Scale(t, sgn * k))), footprint); };
    const float ga = (at(a, 1.0f) - at(a, -1.0f)) / (2.0f * step), gb = (at(b, 1.0f) - at(b, -1.0f)) / (2.0f * step);
    s.grad = Vector3Add(Vector3Scale(a, ga), Vector3Scale(b, gb));
    return s;
}

namespace {
    // sin and cos of 2 pi x by a refined parabola (error ~1e-3): cheap, and the same on every platform.
    inline void SinCos2Pi(float x, float& s, float& c) {
        auto sn = [](float u) { // u in [-0.5, 0.5]
            float y = 8.0f * u - 16.0f * u * fabsf(u);
            return 0.225f * (y * fabsf(y) - y) + y;
        };
        float u = x - floorf(x + 0.5f);
        s = sn(u);
        u += 0.25f; u -= floorf(u + 0.5f);
        c = sn(u);
    }
}

float ErosionGullies(Vector3 p, Vector3 dir, Vector3 grad, float footprint, unsigned int seed, float strength) {
    static const float kLambda[6] = { 900.0f, 420.0f, 190.0f, 85.0f, 38.0f, 17.0f };
    constexpr float kTwoPi = 6.28318531f;
    static const float kGullyAmp = getenv("RAYLIBGAME_GULLY_AMP") ? (float)atof(getenv("RAYLIBGAME_GULLY_AMP")) : 0.14f; // tuning (0.11 read as smooth ridges, user 2026-11-14)
    float h = 0.0f;
    if (strength <= 0.0f) return 0.0f;
    // Regions (~3 km) of smoother and more deeply cut ground: not every slope combed alike (variety, user 2026-11-14).
    strength *= 0.3f + 1.0f * ValueNoise(Vector3Scale(p, 1.0f / 3000.0f), seed ^ 0x6A11u);
    // Strength from the baked slope (the bent gradient below only steers): flat ground stays flat at every octave.
    // Grows with the slope, full only on steep (~42 deg) ground: gentle hills get mild ribs and stay walkable, mountain
    // sides get deep gullies (at 0.14 x and full strength from ~19 deg, 120 m hills reached 57 deg, #363).
    const float sfl = std::clamp(Vector3Length(grad) / 0.9f, 0.0f, 1.0f), sf = sfl * sqrtf(sfl);
    if (sf <= 0.0f) return 0.0f;
    float octaveAmp = kGullyAmp; // x 0.72 an octave: each finer octave adds less slope (0.88, 0.63, 0.46, ... at full strength)
    for (int o = 0; o < 6; ++o, octaveAmp *= 0.72f) {
        const float lam = kLambda[o];
        const float fade = std::clamp(lam / (2.0f * footprint) - 1.0f, 0.0f, 1.0f);
        if (fade <= 0.0f) break;
        const float slope = Vector3Length(grad);
        if (slope < 1e-4f) continue;
        const Vector3 across = Vector3Scale(Vector3CrossProduct(dir, grad), 1.0f / slope); // unit, along the contour
        const float inv = 1.0f / lam;
        const float fx = p.x * inv, fy = p.y * inv, fz = p.z * inv;
        const int x0 = (int)floorf(fx), y0 = (int)floorf(fy), z0 = (int)floorf(fz);
        float val = 0.0f, dval = 0.0f, wsum = 0.0f;
        for (int c = 0; c < 8; ++c) {
            const int cx = x0 + (c & 1), cy = y0 + ((c >> 1) & 1), cz = z0 + (c >> 2);
            const float ox = fx - cx, oy = fy - cy, oz = fz - cz;
            const float r2 = ox * ox + oy * oy + oz * oz;
            if (r2 >= 1.2f) continue;
            float w = 1.0f - r2 / 1.2f;
            w *= w;
            float sn, cs;
            SinCos2Pi(ox * across.x + oy * across.y + oz * across.z + Hash3(cx, cy, cz, seed + 31u * o), sn, cs);
            val += w * cs;
            dval -= w * sn;
            wsum += w;
        }
        if (wsum <= 0.0f) continue;
        val /= wsum; dval /= wsum;
        // Sharp crests, V gullies (user 2026-11-14: ridges were too smooth): x = 0 at a crest, 1 at a gully floor;
        // x^0.6 rises steeply from the crest. (The bending below keeps the plain stripes' slope.)
        { const float x = std::clamp(0.5f * (1.0f - val), 0.0f, 1.0f); val = 1.0f - 2.0f * powf(x, 0.6f); }
        const float amp = octaveAmp * lam * sf * fade * strength;
        h += amp * val;
        // The ribs' own slope across the contour bends the next octave: finer gullies run into these.
        grad = Vector3Add(grad, Vector3Scale(across, amp * dval * kTwoPi * inv));
    }
    return h;
}
