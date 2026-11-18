#include "seedsurvey.h"
#include "newworld.h"
#include "worldgenparams.h"

#include <raylib.h>
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstdarg>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <string>
#include <thread>
#include <vector>

// --seed-survey COUNT FIRSTSEED OUTDIR [THREADS] (#532): how alike are the start planets of different world seeds?
// Builds the New World preview (base terrain, default parameters, no erosion) of COUNT world seeds (half consecutive
// from FIRSTSEED, half random 32-bit ones), compares their land masks pairwise, and writes stats.txt, 16 sample maps
// and a 4x4 contact sheet to OUTDIR. No window. RAYLIBGAME_SEED_SPAN=U widens the noise seed offsets (voxel.cpp).

namespace {
    constexpr int W = 128, H = 64;
    // Latitude band |lat| < 60 deg (ice caps and polar biomes are latitude driven by design): rows 11..52.
    constexpr int kRow0 = 11, kRow1 = 52, kRows = kRow1 - kRow0 + 1;
    constexpr int kWords = kRows * 2; // 128 bits a row
    constexpr int kVariants = W * 2;  // 128 longitude rotations x {as is, N/S flip}

    struct Mask { uint64_t w[kWords]; int count; };

    int Pop(uint64_t v) { return __builtin_popcountll(v); }

    Mask MakeMask(const std::vector<unsigned char>& land) { // land: kRows x W
        Mask m{}; m.count = 0;
        for (int r = 0; r < kRows; ++r)
            for (int x = 0; x < W; ++x)
                if (land[r * W + x]) { m.w[r * 2 + (x >> 6)] |= 1ull << (x & 63); ++m.count; }
        return m;
    }
    // Variant v: rotate by v % W columns, flip rows when v >= W.
    Mask Variant(const std::vector<unsigned char>& land, int v) {
        std::vector<unsigned char> t(kRows * W);
        const int s = v % W; const bool flip = v >= W;
        for (int r = 0; r < kRows; ++r)
            for (int x = 0; x < W; ++x)
                t[r * W + x] = land[(flip ? kRows - 1 - r : r) * W + (x + s) % W];
        return MakeMask(t);
    }
    void Compare(const Mask& a, const Mask& b, int& inter, int& uni) {
        inter = 0; uni = 0;
        for (int k = 0; k < kWords; ++k) { inter += Pop(a.w[k] & b.w[k]); uni += Pop(a.w[k] | b.w[k]); }
    }
    constexpr int kCells = kRows * W;
    float Agree(int inter, int uni) { return 1.0f - (float)(uni - inter) / kCells; }
    float IoU(int inter, int uni) { return uni > 0 ? (float)inter / uni : 1.0f; }

    struct Dist { double mean, p50, p95, p99, max, min; };
    Dist Stats(std::vector<float> v) {
        Dist d{};
        if (v.empty()) return d;
        double s = 0; for (float f : v) s += f;
        std::sort(v.begin(), v.end());
        auto q = [&](double p) { return (double)v[std::min(v.size() - 1, (size_t)(p * (v.size() - 1) + 0.5))]; };
        d.mean = s / v.size(); d.p50 = q(0.5); d.p95 = q(0.95); d.p99 = q(0.99); d.max = v.back(); d.min = v.front();
        return d;
    }
    void Line(FILE* f, const char* name, const Dist& d) {
        char b[256];
        snprintf(b, sizeof b, "  %-44s mean %.3f  p50 %.3f  p95 %.3f  p99 %.3f  max %.3f  (min %.3f)\n", name, d.mean, d.p50, d.p95, d.p99, d.max, d.min);
        fputs(b, f); fputs(b, stdout);
    }
    void Say(FILE* f, const char* fmt, ...) __attribute__((format(printf, 2, 3)));
    void Say(FILE* f, const char* fmt, ...) {
        char b[1024]; va_list ap; va_start(ap, fmt); vsnprintf(b, sizeof b, fmt, ap); va_end(ap);
        fputs(b, f); fputs(b, stdout);
    }

    double Secs(std::chrono::steady_clock::time_point t) { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t).count(); }
    template <class F> void Parallel(int n, int threads, F&& f) {
        std::atomic<int> next{ 0 };
        std::vector<std::thread> ts;
        for (int t = 0; t < threads; ++t) ts.emplace_back([&] { for (int i; (i = next.fetch_add(1)) < n;) f(i); });
        for (auto& t : ts) t.join();
    }
}

int SeedSurvey_Run(int argc, char** argv) {
    if (argc < 5) { printf("usage: --seed-survey COUNT FIRSTSEED OUTDIR [THREADS] [SUBSET]\n"); return 1; }
    const int count = std::max(4, atoi(argv[2]));
    const unsigned int first = (unsigned int)strtoul(argv[3], nullptr, 10);
    std::string out = argv[4];
    if (!out.empty() && out.back() != '/' && out.back() != '\\') out += '/';
    const int threads = argc >= 6 ? std::max(1, atoi(argv[5])) : 4;
    const int subset = std::min(count, argc >= 7 ? atoi(argv[6]) : 400); // pairs searched over shifts/flips, and correlated
    SetTraceLogLevel(LOG_WARNING);
    MakeDirectory(out.c_str());
    MakeDirectory((out + "samples").c_str());

    // Seeds: half consecutive, half random 32-bit (as the New World screen rolls them). Order: s0,r0,s1,r1,... so any
    // prefix (the subset) has both kinds.
    const int half = count / 2;
    std::vector<unsigned int> seeds;
    std::vector<int> kind; // 0 consecutive, 1 random
    std::mt19937 roll(0x532u);
    for (int i = 0; i < count; ++i) {
        if ((i & 1) == 0 && i / 2 < half) { seeds.push_back(first + i / 2); kind.push_back(0); }
        else { unsigned int s = roll(); if (s == 0) s = 1; seeds.push_back(s); kind.push_back(1); }
    }
    unsigned int t0 = 0; float radius = 0.0f;
    if (!WorldGen_StartPlanet(seeds[0], t0, radius)) { printf("seed-survey: no planet in the system\n"); return 1; } // loads the system once
    const WorldGenParams params;
    const char* span = getenv("RAYLIBGAME_SEED_SPAN");
    printf("seed-survey: %d seeds (%d consecutive from %u, %d random), radius %.0f m, %d threads, seed span %s\n",
           count, half, first, count - half, radius, threads, span ? span : "default (10000 m)");

    // Sample every seed's map.
    std::vector<std::vector<float>> heights(count);
    std::vector<unsigned int> tseeds(count);
    const auto tStart = std::chrono::steady_clock::now();
    std::atomic<int> done{ 0 };
    Parallel(count, threads, [&](int i) {
        unsigned int ts = 0; float r = 0.0f;
        WorldGen_StartPlanet(seeds[i], ts, r);
        tseeds[i] = ts;
        std::vector<Color> col;
        NewWorld_RenderMap(ts, r, params, W, H, col, &heights[i]);
        const int d = ++done;
        if (d % 200 == 0) printf("seed-survey: %d / %d maps (%.1f s)\n", d, count, Secs(tStart));
    });

    // Masks over the band, the shuffled control, and the land frequency per cell.
    std::vector<std::vector<unsigned char>> land(count, std::vector<unsigned char>(kRows * W));
    std::vector<Mask> masks(count), shuffled(count);
    std::vector<float> landFrac(count);
    std::vector<int> landCount(kRows * W, 0);
    std::mt19937 shuf(7u);
    for (int i = 0; i < count; ++i) {
        for (int r = 0; r < kRows; ++r)
            for (int x = 0; x < W; ++x) land[i][r * W + x] = heights[i][(kRow0 + r) * W + x] > 0.0f;
        masks[i] = MakeMask(land[i]);
        landFrac[i] = (float)masks[i].count / kCells;
        for (int c = 0; c < kCells; ++c) landCount[c] += land[i][c];
        std::vector<unsigned char> s = land[i];
        std::shuffle(s.begin(), s.end(), shuf);
        shuffled[i] = MakeMask(s);
    }

    std::string statsPath = out + "stats.txt";
    FILE* f = fopen(statsPath.c_str(), "w");
    if (!f) { printf("seed-survey: cannot write %s\n", statsPath.c_str()); return 1; }
    Say(f, "Seed survey (#532): %d world seeds (%d consecutive from %u, %d random 32-bit), start planet radius %.0f m\n",
        count, half, first, count - half, radius);
    Say(f, "Preview sampling: PlanetField::elevationBase (no erosion), default WorldGenParams, %dx%d map, land = elevation > 0,\n", W, H);
    Say(f, "band |lat| < 60 deg (rows %d..%d, %d cells). Seed offset span: %s\n\n", kRow0, kRow1, kCells, span ? span : "default (10000 m)");
    Say(f, "Land fraction per seed:\n");
    Line(f, "land fraction", Stats(landFrac));

    // All pairs: plain IoU and agreement, baselines.
    const size_t nPairs = (size_t)count * (count - 1) / 2;
    std::vector<float> iou(nPairs), agree(nPairs), expAgree(nPairs), expIoU(nPairs), excess(nPairs), shIoU(nPairs), shAgree(nPairs);
    std::vector<size_t> rowStart(count);
    { size_t k = 0; for (int i = 0; i < count; ++i) { rowStart[i] = k; k += count - 1 - i; } }
    Parallel(count, threads, [&](int i) {
        size_t k = rowStart[i];
        for (int j = i + 1; j < count; ++j, ++k) {
            int in, un;
            Compare(masks[i], masks[j], in, un);
            iou[k] = IoU(in, un); agree[k] = Agree(in, un);
            const float p = landFrac[i], q = landFrac[j];
            expAgree[k] = p * q + (1 - p) * (1 - q);
            expIoU[k] = (p + q - p * q) > 0 ? p * q / (p + q - p * q) : 1.0f;
            excess[k] = agree[k] - expAgree[k];
            Compare(shuffled[i], shuffled[j], in, un);
            shIoU[k] = IoU(in, un); shAgree[k] = Agree(in, un);
        }
    });
    Say(f, "\nAll %zu pairs, same position (no shift):\n", nPairs);
    Line(f, "IoU of land", Stats(iou));
    Line(f, "baseline IoU (independent, same fractions)", Stats(expIoU));
    Line(f, "control IoU (cells shuffled per seed)", Stats(shIoU));
    Line(f, "agreement (1 - Hamming/cells)", Stats(agree));
    Line(f, "baseline agreement (independent)", Stats(expAgree));
    Line(f, "control agreement (cells shuffled)", Stats(shAgree));
    Line(f, "excess agreement (agreement - baseline)", Stats(excess));
    {   // by seed kind
        std::vector<float> ex[3];
        for (int i = 0; i < count; ++i) { size_t k = rowStart[i]; for (int j = i + 1; j < count; ++j, ++k) ex[kind[i] + kind[j]].push_back(excess[k]); }
        Line(f, "excess agreement, consecutive x consecutive", Stats(ex[0]));
        Line(f, "excess agreement, consecutive x random", Stats(ex[1]));
        Line(f, "excess agreement, random x random", Stats(ex[2]));
    }

    // Land frequency per cell: with independent seeds every cell is land in about mean(landFrac) of them.
    {
        double mf = 0; for (float v : landFrac) mf += v; mf /= count;
        int hi = 0, lo = 0; std::vector<float> freq(kCells);
        for (int c = 0; c < kCells; ++c) { freq[c] = (float)landCount[c] / count; hi += freq[c] > 0.9f; lo += freq[c] < 0.1f; }
        const Dist d = Stats(freq);
        Say(f, "\nLand frequency per cell over all seeds (independent seeds: every cell ~ %.3f +- %.3f):\n", mf, sqrt(mf * (1 - mf) / count));
        Line(f, "per-cell land frequency", d);
        Say(f, "  cells land in > 90%% of seeds: %d of %d; land in < 10%%: %d\n", hi, kCells, lo);
        Image img = GenImageColor(W, H, MAGENTA);
        for (int y = 0; y < H; ++y)
            for (int x = 0; x < W; ++x) {
                float fr = 0; for (int i = 0; i < count; ++i) fr += heights[i][y * W + x] > 0.0f;
                fr /= count;
                const unsigned char g = (unsigned char)(fr * 255.0f);
                ImageDrawPixel(&img, x, y, { g, g, g, 255 });
            }
        ImageResizeNN(&img, W * 4, H * 4);
        ExportImage(img, (out + "land_frequency.png").c_str());
        UnloadImage(img);
    }

    // Subset: best over 128 longitude rotations x N/S flip; height correlation.
    std::vector<std::vector<Mask>> var(subset, std::vector<Mask>(kVariants)), varSh(subset, std::vector<Mask>(kVariants));
    Parallel(subset, threads, [&](int j) {
        std::vector<unsigned char> s(kCells);
        for (int c = 0; c < kCells; ++c) s[c] = 0;
        // the shuffled mask's cells back to an array
        for (int r = 0; r < kRows; ++r) for (int x = 0; x < W; ++x) s[r * W + x] = (shuffled[j].w[r * 2 + (x >> 6)] >> (x & 63)) & 1;
        for (int v = 0; v < kVariants; ++v) { var[j][v] = Variant(land[j], v); varSh[j][v] = Variant(s, v); }
    });
    const size_t nSub = (size_t)subset * (subset - 1) / 2;
    std::vector<float> bestAgree(nSub), bestIoU(nSub), bestExcess(nSub), bestShAgree(nSub), bestShExcess(nSub), pear(nSub), pearBand(nSub);
    std::vector<int> bestShift(nSub);
    std::vector<size_t> subStart(subset);
    { size_t k = 0; for (int i = 0; i < subset; ++i) { subStart[i] = k; k += subset - 1 - i; } }
    // Per-seed standardised heights for the correlation (whole map, and the band only).
    auto standardise = [&](const std::vector<float>& h, int y0, int y1) {
        std::vector<float> z((size_t)(y1 - y0 + 1) * W);
        double s = 0, s2 = 0;
        for (int y = y0; y <= y1; ++y) for (int x = 0; x < W; ++x) { const float v = h[y * W + x]; s += v; s2 += (double)v * v; }
        const double n = z.size(), m = s / n, sd = sqrt(std::max(s2 / n - m * m, 1e-9));
        for (int y = y0; y <= y1; ++y) for (int x = 0; x < W; ++x) z[(y - y0) * W + x] = (float)((h[y * W + x] - m) / sd);
        return z;
    };
    std::vector<std::vector<float>> zAll(subset), zBand(subset);
    for (int i = 0; i < subset; ++i) { zAll[i] = standardise(heights[i], 0, H - 1); zBand[i] = standardise(heights[i], kRow0, kRow1); }
    Parallel(subset, threads, [&](int i) {
        size_t k = subStart[i];
        for (int j = i + 1; j < subset; ++j, ++k) {
            float ba = -1, bi = 0, bs = -1; int bv = 0;
            for (int v = 0; v < kVariants; ++v) {
                int in, un;
                Compare(masks[i], var[j][v], in, un);
                const float a = Agree(in, un);
                if (a > ba) { ba = a; bi = IoU(in, un); bv = v; }
                Compare(shuffled[i], varSh[j][v], in, un);
                bs = std::max(bs, Agree(in, un));
            }
            const float p = landFrac[i], q = landFrac[j], e = p * q + (1 - p) * (1 - q);
            bestAgree[k] = ba; bestIoU[k] = bi; bestShift[k] = bv; bestExcess[k] = ba - e; bestShAgree[k] = bs; bestShExcess[k] = bs - e;
            double s = 0; for (size_t c = 0; c < zAll[i].size(); ++c) s += zAll[i][c] * zAll[j][c];
            pear[k] = (float)(s / zAll[i].size());
            s = 0; for (size_t c = 0; c < zBand[i].size(); ++c) s += zBand[i][c] * zBand[j][c];
            pearBand[k] = (float)(s / zBand[i].size());
        }
    });
    Say(f, "\nSubset: the first %d seeds (both kinds interleaved), %zu pairs, best of 128 longitude shifts x N/S flip:\n", subset, nSub);
    Line(f, "best agreement over shifts/flips", Stats(bestAgree));
    Line(f, "IoU at that shift", Stats(bestIoU));
    Line(f, "excess best agreement (best - independent)", Stats(bestExcess));
    Line(f, "control: shuffled cells, best agreement", Stats(bestShAgree));
    Line(f, "control: shuffled cells, excess", Stats(bestShExcess));
    {
        int zeroShift = 0; for (size_t k = 0; k < nSub; ++k) zeroShift += bestShift[k] == 0;
        Say(f, "  pairs whose best variant is 'no shift, no flip': %d of %zu (%.1f%%; 1/256 = 0.4%% if unrelated)\n", zeroShift, nSub, 100.0 * zeroShift / nSub);
    }
    Say(f, "\nHeight-field Pearson correlation, same position (subset pairs):\n");
    Line(f, "Pearson r, whole map", Stats(pear));
    Line(f, "Pearson r, |lat| < 60 band", Stats(pearBand));

    // The most similar pairs (plain agreement excess, all pairs).
    {
        std::vector<size_t> idx(nPairs);
        for (size_t k = 0; k < nPairs; ++k) idx[k] = k;
        const size_t top = std::min<size_t>(12, nPairs);
        std::partial_sort(idx.begin(), idx.begin() + top, idx.end(), [&](size_t a, size_t b) { return excess[a] > excess[b]; });
        Say(f, "\nMost similar pairs (same position, by excess agreement):\n");
        for (size_t t = 0; t < top; ++t) {
            const size_t k = idx[t];
            int i = (int)(std::upper_bound(rowStart.begin(), rowStart.end(), k) - rowStart.begin()) - 1;
            const int j = i + 1 + (int)(k - rowStart[i]);
            Say(f, "  world seeds %10u / %10u (terrain %10u / %10u): agreement %.3f (baseline %.3f), IoU %.3f, land %.2f / %.2f\n",
                seeds[i], seeds[j], tseeds[i], tseeds[j], agree[k], expAgree[k], iou[k], landFrac[i], landFrac[j]);
        }
    }

    // 16 sample maps at 512x256 (8 consecutive, 8 random) and a 4x4 contact sheet.
    const int SW = 512, SH = 256;
    std::vector<int> pick;
    for (int i = 0; i < count && (int)pick.size() < 16; ++i) pick.push_back(i);
    std::vector<std::vector<Color>> maps(pick.size());
    Parallel((int)pick.size(), threads, [&](int n) {
        NewWorld_RenderMap(tseeds[pick[n]], radius, params, SW, SH, maps[n], nullptr);
    });
    Image sheet = GenImageColor(SW * 4 + 3 * 6, SH * 4 + 3 * 6, { 20, 20, 24, 255 });
    Say(f, "\nContact sheet (sheet.png, row-major 4x4):\n");
    for (size_t n = 0; n < pick.size(); ++n) {
        Image img = GenImageColor(SW, SH, BLACK);
        for (int y = 0; y < SH; ++y) for (int x = 0; x < SW; ++x) ImageDrawPixel(&img, x, y, maps[n][(size_t)y * SW + x]);
        const char* name = TextFormat("%ssamples/seed_%u.png", out.c_str(), seeds[pick[n]]);
        ExportImage(img, name);
        ImageDraw(&sheet, img, { 0, 0, (float)SW, (float)SH }, { (float)((n % 4) * (SW + 6)), (float)((n / 4) * (SH + 6)), (float)SW, (float)SH }, WHITE);
        UnloadImage(img);
        Say(f, "  [%zu,%zu] world seed %u (%s), terrain seed %u, land %.2f\n", n / 4, n % 4, seeds[pick[n]], kind[pick[n]] ? "random" : "consecutive", tseeds[pick[n]], landFrac[pick[n]]);
    }
    ExportImage(sheet, (out + "sheet.png").c_str());
    UnloadImage(sheet);
    Say(f, "\nTook %.1f s\n", Secs(tStart));
    fclose(f);
    printf("seed-survey: wrote %s\n", statsPath.c_str());
    return 0;
}
