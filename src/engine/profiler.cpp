// The zone ring and its run summary (#871); see profiler.h.
#include "profiler.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <map>
#include <mutex>
#include <string>
#include <vector>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#elif defined(__linux__)
#include <cstdio>
#endif

#ifdef TRACY_ENABLE
void Profiler_Enable(bool) {}
bool Profiler_On() { return false; }
void Profiler_Sample(const char*, double) {}
void Profiler_Counter(const char*, double) {}
std::string Profiler_SummaryJson() { return ""; }
void Profiler_PrintSummary() {}
void Profiler_MeasureZoneCost(double& onNs, double& offNs) { onNs = offNs = 0.0; }
void Profiler_Reset() {}
double Profiler_ZoneMeanMs(const char*) { return -1.0; }
void Profiler_SetHitchStamp(double (*)(), std::string (*)()) {}
void Profiler_SetHitchCap(size_t) {}
bool Profiler_WriteSession(const std::string&, std::string&) { return false; }
#else

extern "C" {
    typedef void (*PfnGetIntegerv)(unsigned int, int*);
    extern PfnGetIntegerv glad_glGetIntegerv;
}

namespace {

// ---- the rings ----
struct Event { const char* name; uint64_t t0, t1; };
constexpr uint32_t kRing = 1u << 14; // 16384 zones a thread between two drains (a frame); 384 KB
struct Ring {
    Event ev[kRing];
    std::atomic<uint32_t> head{ 0 }; // written by the owner only
    uint32_t tail = 0;               // read by the main thread only
    uint64_t lost = 0;
};
std::mutex gRingsMutex;            // only when a thread records its first zone, and when the main thread drains
std::vector<Ring*> gRings;
thread_local Ring* tRing = nullptr;

// ---- the statistics (main thread) ----
// Log histogram: bin = floor(log2(ns) * 8) -> 1/8 octave (~9%) resolution, 10 ns .. 70 s.
constexpr int kBins = 8 * 37;
int BinOf(double ns) { if (ns < 1.0) return 0; const int b = (int)(std::log2(ns) * 8.0); return std::clamp(b, 0, kBins - 1); }
double BinValue(int b) { return std::exp2((b + 0.5) / 8.0); }
struct Stat {
    uint64_t count = 0; double sumNs = 0.0, maxNs = 0.0;
    uint32_t hist[kBins] = {};
    void add(double ns) { ++count; sumNs += ns; maxNs = std::max(maxNs, ns); ++hist[BinOf(ns)]; }
    double pct(double q) const {
        if (!count) return 0.0;
        const uint64_t want = (uint64_t)std::ceil(q * (double)count);
        uint64_t acc = 0;
        for (int b = 0; b < kBins; ++b) { acc += hist[b]; if (acc >= want) return std::min(BinValue(b), maxNs); }
        return maxNs;
    }
};
std::map<std::string, Stat> gZones;              // by name (the same literal can have several addresses)
std::map<const char*, Stat*> gZoneByPtr;
struct Counter { uint64_t n = 0; double sum = 0.0, max = -1e300, last = 0.0; };
std::map<std::string, Counter> gCounters;
Stat gFrames;
struct Hitch { double atS, ms; std::vector<std::pair<std::string, double>> top; double gameS = -1.0; std::string context; };
// (#1372) a hitch is a frame over this (ms): RAYLIBGAME_PROFILE_HITCH_MS lowers it to name the zones of shorter spikes
const double gHitchMs = [] { const char* e = getenv("RAYLIBGAME_PROFILE_HITCH_MS"); const double v = e ? atof(e) : 0.0; return v > 0.0 ? v : 20.0; }();
double (*gGameClock)() = nullptr;          // (#886) stamped on each hitch: the game time, and what was going on
std::string (*gHitchContext)() = nullptr;
size_t gHitchCap = 200;
std::vector<Hitch> gHitches;
uint64_t gHitchesOver = 0;
std::vector<std::pair<const char*, double>> gFrameZones; // this frame's zone totals (ms)
uint64_t gTsc0 = 0; std::chrono::steady_clock::time_point gWall0;
double gNsPerTick = 1.0;
uint64_t gLastFrameTsc = 0;
double gPeakVramMB = 0.0; uint64_t gLastVram = 0;
long long gLostTotal = 0;

void Calibrate() {
    const uint64_t t = Prof::Now();
    const double wallNs = std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - gWall0).count();
    if (t > gTsc0 && wallNs > 1e6) gNsPerTick = wallNs / (double)(t - gTsc0);
}
Stat& ZoneStat(const char* name) {
    auto it = gZoneByPtr.find(name);
    if (it != gZoneByPtr.end()) return *it->second;
    Stat& s = gZones[name];
    gZoneByPtr[name] = &s;
    return s;
}
void AddFrameZone(const char* name, double ms) {
    for (auto& z : gFrameZones) if (z.first == name) { z.second += ms; return; }
    gFrameZones.push_back({ name, ms });
}
std::string Esc(const std::string& s) { std::string o; for (char c : s) { if (c == '"' || c == '\\') o += '\\'; o += c; } return o; }
std::string N(double v) { char b[48]; snprintf(b, sizeof b, "%.4g", v); return b; }

double PeakCommitMB() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) return pmc.PeakPagefileUsage / 1048576.0;
#elif defined(__linux__)
    if (FILE* f = fopen("/proc/self/status", "r")) {
        char line[256]; double kb = 0.0;
        while (fgets(line, sizeof line, f)) if (sscanf(line, "VmPeak: %lf kB", &kb) == 1) break;
        fclose(f);
        return kb / 1024.0;
    }
#endif
    return 0.0;
}
void SampleVram() {
    // NVX_gpu_memory_info (NVIDIA): total minus available, in KB. Other drivers: 0 (not reported).
    if (!glad_glGetIntegerv) return;
    int total = 0, avail = 0;
    glad_glGetIntegerv(0x9048, &total);
    glad_glGetIntegerv(0x9049, &avail);
    if (total > 0 && avail > 0 && total >= avail) gPeakVramMB = std::max(gPeakVramMB, (total - avail) / 1024.0);
}

struct AutoEnable { AutoEnable() { const char* e = getenv("RAYLIBGAME_PROFILE"); if (e && e[0] == '1') Profiler_Enable(true); } } gAutoEnable;

} // namespace

namespace Prof {
void Record(const char* name, uint64_t t0, uint64_t t1) {
    Ring* r = tRing;
    if (!r) {
        r = new Ring(); // once per thread; kept for the process (the main thread may still drain it)
        std::lock_guard<std::mutex> lock(gRingsMutex);
        gRings.push_back(r);
        tRing = r;
    }
    const uint32_t h = r->head.load(std::memory_order_relaxed);
    r->ev[h & (kRing - 1)] = { name, t0, t1 };
    r->head.store(h + 1, std::memory_order_release);
}

void Frame() {
    if (!gOn) return;
    const uint64_t now = Now();
    // (gFrameZones holds what Profiler_Sample added since the last frame: cleared at the end)
    {
        std::lock_guard<std::mutex> lock(gRingsMutex);
        for (Ring* r : gRings) {
            const uint32_t head = r->head.load(std::memory_order_acquire);
            uint32_t tail = r->tail;
            if (head - tail > kRing - 64) { const uint32_t skip = head - tail - (kRing - 64); r->lost += skip; gLostTotal += skip; tail += skip; } // lapped: the oldest are gone
            for (; tail != head; ++tail) {
                const Event e = r->ev[tail & (kRing - 1)];
                const double ns = (double)(e.t1 - e.t0) * gNsPerTick;
                ZoneStat(e.name).add(ns);
                AddFrameZone(e.name, ns / 1e6);
            }
            r->tail = tail;
        }
    }
    if (gLastFrameTsc) {
        const double ms = (double)(now - gLastFrameTsc) * gNsPerTick / 1e6;
        gFrames.add(ms * 1e6);
        if (ms > gHitchMs) {
            ++gHitchesOver;
            if (gHitches.size() < gHitchCap) {
                Hitch h{ (double)(now - gTsc0) * gNsPerTick / 1e9, ms, {} };
                if (gGameClock) h.gameS = gGameClock();
                if (gHitchContext) h.context = gHitchContext();
                std::sort(gFrameZones.begin(), gFrameZones.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
                for (size_t i = 0; i < gFrameZones.size() && i < 3; ++i) h.top.push_back({ gFrameZones[i].first, gFrameZones[i].second });
                gHitches.push_back(std::move(h));
            }
        }
    }
    gLastFrameTsc = now;
    gFrameZones.clear();
    if (gFrames.count % 60 == 0) { Calibrate(); SampleVram(); }
}
} // namespace Prof

void Profiler_Enable(bool on) {
    if (on && !Prof::gOn) {
        gTsc0 = Prof::Now(); gWall0 = std::chrono::steady_clock::now();
        // A first estimate of the tick rate over 2 ms; Frame() refines it over the run.
        while (std::chrono::steady_clock::now() - gWall0 < std::chrono::milliseconds(2)) {}
        Calibrate();
    }
    Prof::gOn = on;
}
bool Profiler_On() { return Prof::gOn; }
void Profiler_Reset() {
    { // drop what the rings hold
        std::lock_guard<std::mutex> lock(gRingsMutex);
        for (Ring* r : gRings) r->tail = r->head.load(std::memory_order_acquire);
    }
    for (auto& kv : gZones) kv.second = Stat{};
    gCounters.clear(); gFrames = Stat{}; gHitches.clear(); gHitchesOver = 0; gFrameZones.clear(); gLastFrameTsc = 0; gLostTotal = 0;
}
void Profiler_Sample(const char* name, double ms) {
    if (!Prof::gOn) return;
    static std::map<std::string, std::string> interned; // callers may pass built names: keep one stable copy each
    const char* key = interned.emplace(name, name).first->second.c_str();
    ZoneStat(key).add(ms * 1e6);
    AddFrameZone(key, ms);
}
void Profiler_Counter(const char* name, double v) {
    if (!Prof::gOn) return;
    Counter& c = gCounters[name];
    ++c.n; c.sum += v; c.max = std::max(c.max, v); c.last = v;
}

double Profiler_ZoneMeanMs(const char* name) {
    if (!Prof::gOn) return -1.0;
    const auto it = gZones.find(name);
    return it == gZones.end() || it->second.count == 0 ? -1.0 : it->second.sumNs / (double)it->second.count / 1e6;
}

std::string Profiler_SummaryJson() {
    if (!Prof::gOn || (gFrames.count == 0 && gZones.empty())) return "";
    Calibrate();
    std::string j = "{\n    \"frame\": {\"count\": " + std::to_string(gFrames.count) + ", \"mean_ms\": " + N(gFrames.count ? gFrames.sumNs / gFrames.count / 1e6 : 0.0) +
                    ", \"p50_ms\": " + N(gFrames.pct(0.5) / 1e6) + ", \"p95_ms\": " + N(gFrames.pct(0.95) / 1e6) + ", \"p99_ms\": " + N(gFrames.pct(0.99) / 1e6) +
                    ", \"max_ms\": " + N(gFrames.maxNs / 1e6) + ", \"over_20ms\": " + std::to_string(gHitchesOver) + "},\n";
    // The frame histogram: frames per bucket of [lo, hi) ms.
    const double edges[] = { 0, 4, 8, 12, 16.7, 20, 25, 33.3, 50, 100, 1e9 };
    std::vector<uint64_t> hb(10, 0);
    for (int b = 0; b < kBins; ++b) {
        const double ms = BinValue(b) / 1e6;
        for (int k = 0; k < 10; ++k) if (ms >= edges[k] && ms < edges[k + 1]) { hb[k] += gFrames.hist[b]; break; }
    }
    j += "    \"frame_histogram\": [";
    for (int k = 0; k < 10; ++k) j += std::string(k ? ", " : "") + "{\"lo_ms\": " + N(edges[k]) + ", \"hi_ms\": " + (k == 9 ? std::string("null") : N(edges[k + 1])) + ", \"frames\": " + std::to_string(hb[k]) + "}";
    j += "],\n    \"zones\": [";
    std::vector<std::pair<std::string, const Stat*>> zs;
    for (const auto& kv : gZones) if (kv.second.count) zs.push_back({ kv.first, &kv.second }); // (zones reset and not seen since: left out)
    std::sort(zs.begin(), zs.end(), [](const auto& a, const auto& b) { return a.second->sumNs > b.second->sumNs; });
    for (size_t i = 0; i < zs.size(); ++i) {
        const Stat& s = *zs[i].second;
        j += std::string(i ? ",\n      " : "\n      ") + "{\"name\": \"" + Esc(zs[i].first) + "\", \"count\": " + std::to_string(s.count) + ", \"mean_ms\": " + N(s.sumNs / std::max<uint64_t>(1, s.count) / 1e6) +
             ", \"p50_ms\": " + N(s.pct(0.5) / 1e6) + ", \"p95_ms\": " + N(s.pct(0.95) / 1e6) + ", \"p99_ms\": " + N(s.pct(0.99) / 1e6) + ", \"max_ms\": " + N(s.maxNs / 1e6) +
             ", \"per_frame_ms\": " + N(gFrames.count ? s.sumNs / gFrames.count / 1e6 : 0.0) + "}";
    }
    j += zs.empty() ? "],\n" : "\n    ],\n";
    j += "    \"hitches\": [";
    for (size_t i = 0; i < gHitches.size(); ++i) {
        const Hitch& h = gHitches[i];
        j += std::string(i ? ",\n      " : "\n      ") + "{\"at_s\": " + N(h.atS) + ", \"ms\": " + N(h.ms) +
             (h.gameS >= 0.0 ? ", \"game_s\": " + N(h.gameS) : std::string()) + (h.context.empty() ? std::string() : ", \"context\": \"" + Esc(h.context) + "\"") + ", \"top\": [";
        for (size_t k = 0; k < h.top.size(); ++k) j += std::string(k ? ", " : "") + "{\"zone\": \"" + Esc(h.top[k].first) + "\", \"ms\": " + N(h.top[k].second) + "}";
        j += "]}";
    }
    j += gHitches.empty() ? "],\n" : "\n    ],\n";
    j += "    \"counters\": [";
    bool first = true;
    for (const auto& kv : gCounters) {
        j += std::string(first ? "" : ", ") + "{\"name\": \"" + Esc(kv.first) + "\", \"mean\": " + N(kv.second.sum / std::max<uint64_t>(1, kv.second.n)) + ", \"max\": " + N(kv.second.max) + ", \"last\": " + N(kv.second.last) + "}";
        first = false;
    }
    j += "],\n    \"peak_commit_mb\": " + N(PeakCommitMB()) + ", \"peak_vram_mb\": " + N(gPeakVramMB) + ", \"lost_zones\": " + std::to_string(gLostTotal) + "\n  }";
    return j;
}

void Profiler_PrintSummary() {
    if (!Prof::gOn || gFrames.count == 0) return;
    Calibrate();
    printf("profile: %llu frames, mean %.2f ms, p50 %.2f, p95 %.2f, p99 %.2f, max %.2f ms; %llu over %g ms; peak commit %.0f MB, VRAM %.0f MB\n",
           (unsigned long long)gFrames.count, gFrames.sumNs / gFrames.count / 1e6, gFrames.pct(0.5) / 1e6, gFrames.pct(0.95) / 1e6, gFrames.pct(0.99) / 1e6, gFrames.maxNs / 1e6,
           (unsigned long long)gHitchesOver, gHitchMs, PeakCommitMB(), gPeakVramMB);
    std::vector<std::pair<std::string, const Stat*>> zs;
    for (const auto& kv : gZones) if (kv.second.count) zs.push_back({ kv.first, &kv.second }); // (zones reset and not seen since: left out)
    std::sort(zs.begin(), zs.end(), [](const auto& a, const auto& b) { return a.second->sumNs > b.second->sumNs; });
    for (size_t i = 0; i < zs.size() && i < 12; ++i)
        printf("profile: %-28s %7.3f ms/frame  n %-7llu p50 %.3f p95 %.3f p99 %.3f max %.3f ms\n", zs[i].first.c_str(), zs[i].second->sumNs / gFrames.count / 1e6,
               (unsigned long long)zs[i].second->count, zs[i].second->pct(0.5) / 1e6, zs[i].second->pct(0.95) / 1e6, zs[i].second->pct(0.99) / 1e6, zs[i].second->maxNs / 1e6);
}

void Profiler_SetHitchStamp(double (*gameClock)(), std::string (*context)()) { gGameClock = gameClock; gHitchContext = context; }
void Profiler_SetHitchCap(size_t cap) { gHitchCap = cap; }

// (#886) A play session's summary: <dir>/profile-<date-time>.json (Profiler_SummaryJson with the session around it) and
// .txt (the same, to read: frames, the zones by time per frame, every hitch with its game time and what was going on).
bool Profiler_WriteSession(const std::string& dir, std::string& path) {
    const std::string summary = Profiler_SummaryJson();
    if (summary.empty()) return false;
    Calibrate();
    char stamp[32]; const time_t now = time(nullptr); strftime(stamp, sizeof stamp, "%Y-%m-%d_%H%M%S", localtime(&now));
    const std::string base = dir + "/profile-" + stamp;
    const double wallS = std::chrono::duration<double>(std::chrono::steady_clock::now() - gWall0).count();
    if (FILE* f = fopen((base + ".json").c_str(), "w")) {
        fprintf(f, "{\n  \"kind\": \"play session\",\n  \"ended\": \"%s\",\n  \"wall_s\": %.1f,\n  \"profile\": %s\n}\n", stamp, wallS, summary.c_str());
        fclose(f);
    } else return false;
    FILE* f = fopen((base + ".txt").c_str(), "w");
    if (!f) return false;
    fprintf(f, "RaylibGame profile, a play session ending %s (%.0f s of wall time)\n\n", stamp, wallS);
    const double mean = gFrames.count ? gFrames.sumNs / gFrames.count / 1e6 : 0.0;
    fprintf(f, "Frames: %llu, mean %.2f ms (%.0f FPS), p50 %.2f, p95 %.2f, p99 %.2f, max %.2f ms; %llu frames over 20 ms\n",
            (unsigned long long)gFrames.count, mean, mean > 0.0 ? 1000.0 / mean : 0.0, gFrames.pct(0.5) / 1e6, gFrames.pct(0.95) / 1e6,
            gFrames.pct(0.99) / 1e6, gFrames.maxNs / 1e6, (unsigned long long)gHitchesOver);
    fprintf(f, "Peak commit %.0f MB, peak VRAM %.0f MB%s\n\n", PeakCommitMB(), gPeakVramMB, gPeakVramMB > 0.0 ? "" : " (not reported by this driver)");
    std::vector<std::pair<std::string, const Stat*>> zs;
    for (const auto& kv : gZones) if (kv.second.count) zs.push_back({ kv.first, &kv.second });
    std::sort(zs.begin(), zs.end(), [](const auto& a, const auto& b) { return a.second->sumNs > b.second->sumNs; });
    fprintf(f, "Zones, by time per frame (ms):\n  %-32s %9s %9s %9s %9s %9s %10s\n", "zone", "/frame", "p50", "p95", "p99", "max", "count");
    for (const auto& z : zs)
        fprintf(f, "  %-32s %9.3f %9.3f %9.3f %9.3f %9.2f %10llu\n", z.first.c_str(), gFrames.count ? z.second->sumNs / gFrames.count / 1e6 : 0.0,
                z.second->pct(0.5) / 1e6, z.second->pct(0.95) / 1e6, z.second->pct(0.99) / 1e6, z.second->maxNs / 1e6, (unsigned long long)z.second->count);
    const std::string kept = gHitchesOver > gHitches.size() ? " (the first " + std::to_string(gHitches.size()) + " of " + std::to_string(gHitchesOver) + ")" : "";
    fprintf(f, "\nHitches (frames over 20 ms)%s: when, how long, the zones that took most of that frame (ms summed over all threads: the chunk builds run on the workers, so they can add up past the frame), and what was going on\n", kept.c_str());
    for (const Hitch& h : gHitches) {
        char when[48];
        if (h.gameS >= 0.0) snprintf(when, sizeof when, "game %d:%04.1f", (int)(h.gameS / 60.0), fmod(h.gameS, 60.0));
        else snprintf(when, sizeof when, "wall %.1f s", h.atS);
        fprintf(f, "  %-16s %7.1f ms  ", when, h.ms);
        for (size_t k = 0; k < h.top.size(); ++k) fprintf(f, "%s%s %.1f", k ? ", " : "", h.top[k].first.c_str(), h.top[k].second);
        if (!h.context.empty()) fprintf(f, "   [%s]", h.context.c_str());
        fprintf(f, "\n");
    }
    if (gHitches.empty()) fprintf(f, "  none\n");
    fprintf(f, "\nCounters (mean / max / last):\n");
    for (const auto& kv : gCounters) fprintf(f, "  %-32s %12.4g %12.4g %12.4g\n", kv.first.c_str(), kv.second.sum / std::max<uint64_t>(1, kv.second.n), kv.second.max, kv.second.last);
    if (gCounters.empty()) fprintf(f, "  none\n");
    fclose(f);
    path = base + ".txt";
    return true;
}

void Profiler_MeasureZoneCost(double& onNs, double& offNs) {
    // A zone recorded 200k times on this thread (the ring drained as it goes), on and off.
    const bool was = Prof::gOn;
    auto run = [](bool on) {
        Prof::gOn = on;
        const auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < 200000; ++i) {
            { PROFILE_SCOPE("profiler/overhead probe"); }
            if ((i & 4095) == 4095 && tRing) tRing->tail = tRing->head.load(); // keep the ring from lapping (the probe isn't a real zone)
        }
        return std::chrono::duration<double, std::nano>(std::chrono::steady_clock::now() - t0).count() / 200000.0;
    };
    run(true); // warm up (the ring's first allocation)
    onNs = run(true);
    offNs = run(false);
    if (tRing) tRing->tail = tRing->head.load();
    Prof::gOn = was;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#871): zones from several threads reach the summary; percentiles; a reset ----
#include "unit_helpers.h"
#include <thread>

TEST(Profiler, ZonesFromThreadsReachTheSummary) {
    Profiler_Enable(true);
    Profiler_Reset();
    Prof::Frame();
    auto work = [] { for (int i = 0; i < 1000; ++i) { PROFILE_SCOPE("unit/worker zone"); volatile int x = 0; for (int k = 0; k < 50; ++k) x = x + k; } };
    std::thread a(work), b(work);
    a.join(); b.join();
    { PROFILE_SCOPE("unit/main zone"); std::this_thread::sleep_for(std::chrono::milliseconds(2)); }
    Profiler_Sample("gpu:unit pass", 1.5);
    Profiler_Counter("unit counter", 3.0);
    Profiler_Counter("unit counter", 5.0);
    Prof::Frame();
    const std::string j = Profiler_SummaryJson();
    Profiler_Enable(false);
    EXPECT_NE(j.find("\"name\": \"unit/worker zone\", \"count\": 2000"), std::string::npos) << j;
    EXPECT_NE(j.find("\"name\": \"unit/main zone\", \"count\": 1"), std::string::npos) << j;
    EXPECT_NE(j.find("\"name\": \"gpu:unit pass\", \"count\": 1"), std::string::npos) << j;
    EXPECT_NE(j.find("\"name\": \"unit counter\", \"mean\": 4, \"max\": 5"), std::string::npos) << j;
    EXPECT_NE(j.find("\"lost_zones\": 0"), std::string::npos) << j;
}

TEST(Profiler, PercentilesFromTheHistogramAreWithinABin) {
    Stat s;
    for (int i = 1; i <= 1000; ++i) s.add(i * 1000.0); // 1..1000 us
    EXPECT_NEAR(s.pct(0.5), 500e3, 500e3 * 0.1);
    EXPECT_NEAR(s.pct(0.99), 990e3, 990e3 * 0.1);
    EXPECT_DOUBLE_EQ(s.pct(1.0), 1000e3);
    EXPECT_EQ(s.count, 1000u);
}

TEST(Profiler, OffRecordsNothingAndCostsLittle) {
    Profiler_Enable(true);
    Profiler_Reset();
    Profiler_Enable(false);
    for (int i = 0; i < 100; ++i) { PROFILE_SCOPE("unit/off zone"); }
    Profiler_Enable(true);
    Prof::Frame();
    const std::string j = Profiler_SummaryJson();
    // (#1817) The best of 5 runs, and "off" judged against "on" from the same runs: a busy runner (12.5 ns off once against a
    // fixed 5 ns) slows both, while an off zone that still records or reads the clock costs a large share of an on one.
    double onNs = 1e9, offNs = 1e9;
    for (int k = 0; k < 5; ++k) {
        double on = 0.0, off = 0.0;
        Profiler_MeasureZoneCost(on, off);
        onNs = std::min(onNs, on);
        offNs = std::min(offNs, off);
    }
    Profiler_Enable(false);
    EXPECT_EQ(j.find("unit/off zone"), std::string::npos);
    EXPECT_LT(offNs, std::max(5.0, 0.25 * onNs)) << "off " << offNs << " ns, on " << onNs << " ns";
    EXPECT_LT(onNs, 200.0) << "a zone on: two rdtsc and a ring write (the plan: 20-40 ns; loose here, machines vary)";
}
#endif
#endif

// Memory now (tests: reloads, #902), in both builds: the process's commit charge and the GPU memory in use (NVIDIA's
// NVX_gpu_memory_info: total minus available; 0 elsewhere).
extern "C" {
    typedef void (*PfnGetIntegervMem)(unsigned int, int*);
    extern PfnGetIntegervMem glad_glGetIntegerv;
}
double Profiler_CommitMB() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS pmc{};
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof pmc)) return pmc.PagefileUsage / 1048576.0;
#endif
    return 0.0;
}
double Profiler_VramUsedMB() {
    if (!glad_glGetIntegerv) return 0.0;
    int total = 0, avail = 0;
    glad_glGetIntegerv(0x9048, &total);
    glad_glGetIntegerv(0x9049, &avail);
    return total > 0 && avail > 0 && total >= avail ? (total - avail) / 1024.0 : 0.0;
}
