#ifndef PROFILER_H
#define PROFILER_H

// Profiling zones. PROFILE_SCOPE("name") / PROFILE_FUNCTION() time the rest of the scope; PROFILE_FRAME() ends a frame.
//
// With Tracy (CMake -DENABLE_TRACY=ON or the "profile" preset) they are Tracy zones: run the Tracy GUI (tracy-profiler.exe,
// v0.11.1) and it connects to the game.
//
// Without Tracy (#871) they record into a fixed per-thread ring: a name pointer and two timestamps per zone (rdtsc), no
// allocation and no locks after a thread's first zone. Off by default: a zone then costs one branch. On with
// RAYLIBGAME_PROFILE=1 (scenarios turn it on, scenario.h) or Profiler_Enable(true). The main thread drains the rings at each
// PROFILE_FRAME() into per-zone statistics (count, mean, p50/p95/p99 from a log histogram, max), a frame-time histogram, a
// hitch list (frames over 20 ms with their top zones), named counters, gputimer's CPU and GPU passes (as "cpu:" and "gpu:"
// zones) and the peak commit and VRAM. Profiler_SummaryJson() is the run's summary (results.json's "profile").
#include <cstdint>
#include <string>

#ifdef TRACY_ENABLE
    #include <tracy/Tracy.hpp>
    #define PROFILE_CONCAT_(a, b) a##b
    #define PROFILE_CONCAT(a, b) PROFILE_CONCAT_(a, b)
    // Unique variable per line, so several zones can share one scope.
    #define PROFILE_SCOPE(name) ZoneNamedN(PROFILE_CONCAT(profileZone_, __LINE__), name, true)
    #define PROFILE_FUNCTION()  ZoneScoped
    #define PROFILE_FRAME()     FrameMark
#else
    #if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
        #include <x86intrin.h>
    #else
        #include <chrono>
    #endif
    namespace Prof {
        inline bool gOn = false;
        inline uint64_t Now() {
    #if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
            return __rdtsc();
    #else
            return (uint64_t)std::chrono::steady_clock::now().time_since_epoch().count();
    #endif
        }
        void Record(const char* name, uint64_t t0, uint64_t t1); // into this thread's ring
        void Frame();                                            // PROFILE_FRAME (main thread)
        struct Zone {
            const char* name; uint64_t t0;
            explicit Zone(const char* n) : name(gOn ? n : nullptr), t0(name ? Now() : 0) {}
            ~Zone() { if (name) Record(name, t0, Now()); }
            Zone(const Zone&) = delete; Zone& operator=(const Zone&) = delete;
        };
    }
    #define PROFILE_CONCAT_(a, b) a##b
    #define PROFILE_CONCAT(a, b) PROFILE_CONCAT_(a, b)
    #define PROFILE_SCOPE(name) Prof::Zone PROFILE_CONCAT(profileZone_, __LINE__)(name)
    #define PROFILE_FUNCTION()  Prof::Zone PROFILE_CONCAT(profileZone_, __LINE__)(__func__)
    #define PROFILE_FRAME()     Prof::Frame()
#endif

void Profiler_Enable(bool on);          // (no-op with Tracy)
bool Profiler_On();
void Profiler_Sample(const char* name, double ms);   // a time measured elsewhere (gputimer's passes), main thread
void Profiler_Counter(const char* name, double v);   // a named number once a frame (bricks awake, bytes sent), main thread
void Profiler_Reset();                   // forget everything so far (a timed window starts; main thread)
double Profiler_ZoneMeanMs(const char* name); // the mean of a zone since the last reset (gputimer passes are "cpu:NAME" and "gpu:NAME"); -1 when unseen
double Profiler_ZoneMaxMs(const char* name);  // (#1795) its longest run since the last reset; -1 when unseen
std::string Profiler_SummaryJson();      // "" when off or nothing was recorded
// (#886) Play sessions (RAYLIBGAME_PROFILE=1 outside a scenario): each hitch stamped with the game time and a short note of
// what was going on; how many hitches to keep (200; a long session: more); and the summary written on quit as
// <dir>/profile-<date-time>.json and .txt (path: the .txt). False when off or nothing was recorded.
void Profiler_SetHitchStamp(double (*gameClock)(), std::string (*context)());
void Profiler_SetHitchCap(size_t cap);
bool Profiler_WriteSession(const std::string& dir, std::string& path);
void Profiler_PrintSummary();            // "profile: ..." lines
// Microbenchmark for the overhead check: ns per zone (on), and per zone while off.
void Profiler_MeasureZoneCost(double& onNs, double& offNs);
double Profiler_CommitMB();               // the process's commit charge now (MB; tests, #902)
double Profiler_VramUsedMB();             // GPU memory in use now (MB; NVIDIA only, else 0)
double Profiler_ThreadCpuSeconds();       // CPU seconds the calling thread has run (#1914: below its wall time = starved)
void Profiler_MemoryHealth(double& freePhysMB, unsigned long long& pageFaults); // free RAM now, the process's page faults so far (#1914)

#endif
