#include "gputimer.h"

#include <rlgl.h>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <vector>

// glad's GL entry points, as raylib loads them (the same trick as viewmodel.cpp's glClear).
extern "C" {
    typedef void (*PfnGenQueries)(int, unsigned int*);
    typedef void (*PfnQueryCounter)(unsigned int, unsigned int);
    typedef void (*PfnGetQueryObjectiv)(unsigned int, unsigned int, int*);
    typedef void (*PfnGetQueryObjectui64v)(unsigned int, unsigned int, uint64_t*);
    extern PfnGenQueries glad_glGenQueries;
    extern PfnQueryCounter glad_glQueryCounter;
    extern PfnGetQueryObjectiv glad_glGetQueryObjectiv;
    extern PfnGetQueryObjectui64v glad_glGetQueryObjectui64v;
}

namespace {
    constexpr unsigned int kTimestamp = 0x8E28, kResult = 0x8866, kResultAvailable = 0x8867;
    constexpr int kFrames = 4, kMarks = 24;
    struct Slot { unsigned int q[kMarks + 1] = {}; const char* name[kMarks + 1] = {}; int used = 0; bool pending = false; };
    Slot slots[kFrames];
    int cur = -1;
    bool enabled = false, made = false;
    struct Sum { const char* name; double ms; };
    std::vector<Sum> sums;
    long frames = 0;
    double totalMs = 0.0;
    // CPU time between marks on the main thread (same names), for frames while enabled.
    std::vector<Sum> cpuSums;
    long cpuFrames = 0;
    std::chrono::steady_clock::time_point cpuLast;
    bool cpuStarted = false;
    void AddTo(std::vector<Sum>& v, const char* name, double ms) {
        for (Sum& x : v) if (strcmp(x.name, name) == 0) { x.ms += ms; return; }
        v.push_back({ name, ms });
    }

    bool Ready() { return glad_glGenQueries && glad_glQueryCounter && glad_glGetQueryObjectiv && glad_glGetQueryObjectui64v; }

    void Collect(Slot& s) {
        if (!s.pending) return;
        s.pending = false;
        if (s.used < 2) return;
        int avail = 0;
        glad_glGetQueryObjectiv(s.q[s.used - 1], kResultAvailable, &avail);
        if (!avail) return; // dropped rather than waited for
        std::vector<uint64_t> t(s.used);
        for (int i = 0; i < s.used; ++i) glad_glGetQueryObjectui64v(s.q[i], kResult, &t[i]);
        for (int i = 1; i < s.used; ++i) {
            const double ms = (double)(t[i] - t[i - 1]) / 1e6;
            Sum* hit = nullptr;
            for (Sum& x : sums) if (strcmp(x.name, s.name[i]) == 0) { hit = &x; break; }
            if (!hit) { sums.push_back({ s.name[i], 0.0 }); hit = &sums.back(); }
            hit->ms += ms;
        }
        static uint64_t prevEnd = 0;
        if (prevEnd && t[0] > prevEnd && t[0] - prevEnd < 100000000ull) AddTo(sums, "(between frames)", (double)(t[0] - prevEnd) / 1e6);
        prevEnd = t[s.used - 1];
        totalMs += (double)(t[s.used - 1] - t[0]) / 1e6;
        ++frames;
    }
}

void GpuTimer_Enable(bool on) { enabled = on && Ready(); }

void GpuTimer_BeginFrame() {
    if (!enabled) return;
    if (!made) { for (Slot& s : slots) glad_glGenQueries(kMarks + 1, s.q); made = true; }
    cur = (cur + 1) % kFrames;
    Collect(slots[cur]); // written kFrames frames ago
    Slot& s = slots[cur];
    s.used = 0;
    s.pending = true;
    GpuTimer_Mark("start");
}

void GpuTimer_CpuMark(const char* name);
void GpuTimer_Mark(const char* name) {
    if (!enabled || cur < 0) return;
    Slot& s = slots[cur];
    if (s.used > kMarks) return;
    GpuTimer_CpuMark(name);
    rlDrawRenderBatchActive(); // raylib's batched draws belong to the pass that made them
    s.name[s.used] = name;
    glad_glQueryCounter(s.q[s.used], kTimestamp);
    ++s.used;
}

void GpuTimer_Print() {
    if (frames == 0) { printf("bench: gpu: no timings (timer queries unavailable)\n"); return; }
    printf("bench: gpu per frame (avg of %ld):", frames);
    for (const Sum& x : sums) printf(" %s %.2f", x.name, x.ms / frames);
    printf(" | total %.2f ms\n", totalMs / frames);
}

void GpuTimer_CpuFrameStart() {
    if (!enabled) { cpuStarted = false; return; }
    cpuLast = std::chrono::steady_clock::now();
    cpuStarted = true;
    ++cpuFrames;
}

void GpuTimer_CpuMark(const char* name) {
    if (!enabled || !cpuStarted) return;
    const auto now = std::chrono::steady_clock::now();
    AddTo(cpuSums, name, std::chrono::duration<double, std::milli>(now - cpuLast).count());
    cpuLast = now;
}

void GpuTimer_PrintCpu() {
    if (cpuFrames == 0) return;
    double total = 0.0;
    printf("bench: cpu per frame (avg of %ld):", cpuFrames);
    for (const Sum& x : cpuSums) { printf(" %s %.2f", x.name, x.ms / cpuFrames); total += x.ms; }
    printf(" | total %.2f ms\n", total / cpuFrames);
}
