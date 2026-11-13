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
    typedef void (*PfnBeginQuery)(unsigned int, unsigned int);
    typedef void (*PfnEndQuery)(unsigned int);
    typedef void (*PfnGetQueryObjectiv)(unsigned int, unsigned int, int*);
    typedef void (*PfnGetQueryObjectui64v)(unsigned int, unsigned int, uint64_t*);
    extern PfnGenQueries glad_glGenQueries;
    extern PfnBeginQuery glad_glBeginQuery;
    extern PfnEndQuery glad_glEndQuery;
    extern PfnGetQueryObjectiv glad_glGetQueryObjectiv;
    extern PfnGetQueryObjectui64v glad_glGetQueryObjectui64v;
}

// Each segment between two marks is a GL_TIME_ELAPSED query: the GPU time the commands in it took to
// execute. (Timestamp queries, used before, were stamped as the commands were submitted on this driver:
// queued work showed up as "swap+poll" - 2.5 ms that glFinish proved was draw work, not presenting.)
namespace {
    constexpr unsigned int kTimeElapsed = 0x88BF, kResult = 0x8866, kResultAvailable = 0x8867;
    constexpr int kFrames = 4, kMarks = 24;
    struct Slot { unsigned int q[kMarks + 1] = {}; const char* name[kMarks + 1] = {}; int used = 0; bool pending = false; };
    Slot slots[kFrames];
    int cur = -1;
    bool enabled = false, made = false, open = false;
    struct Sum { const char* name; double ms; };
    std::vector<Sum> sums;
    long frames = 0;
    double totalMs = 0.0;
    std::vector<Sum> cpuSums;
    long cpuFrames = 0;
    std::chrono::steady_clock::time_point cpuLast;
    bool cpuStarted = false;
    void AddTo(std::vector<Sum>& v, const char* name, double ms) {
        for (Sum& x : v) if (strcmp(x.name, name) == 0) { x.ms += ms; return; }
        v.push_back({ name, ms });
    }

    bool Ready() { return glad_glGenQueries && glad_glBeginQuery && glad_glEndQuery && glad_glGetQueryObjectiv && glad_glGetQueryObjectui64v; }

    void Collect(Slot& s) {
        if (!s.pending) return;
        s.pending = false;
        if (s.used < 1) return;
        int avail = 0;
        glad_glGetQueryObjectiv(s.q[s.used - 1], kResultAvailable, &avail);
        if (!avail) return; // dropped rather than waited for
        double frame = 0.0;
        for (int i = 0; i < s.used; ++i) {
            uint64_t ns = 0;
            glad_glGetQueryObjectui64v(s.q[i], kResult, &ns);
            const double ms = (double)ns / 1e6;
            AddTo(sums, s.name[i], ms);
            frame += ms;
        }
        totalMs += frame;
        ++frames;
    }
}

void GpuTimer_Enable(bool on) { enabled = on && Ready(); }

void GpuTimer_BeginFrame() {
    if (!enabled) return;
    if (!made) { for (Slot& s : slots) glad_glGenQueries(kMarks + 1, s.q); made = true; }
    if (open && cur >= 0) { // the previous frame's last segment (after its last mark: the swap) ends here
        rlDrawRenderBatchActive();
        glad_glEndQuery(kTimeElapsed);
        open = false;
    }
    cur = (cur + 1) % kFrames;
    Collect(slots[cur]); // written kFrames frames ago
    Slot& s = slots[cur];
    s.used = 0;
    s.pending = true;
    s.name[0] = "start";
    glad_glBeginQuery(kTimeElapsed, s.q[0]);
    open = true;
}

void GpuTimer_CpuMark(const char* name);
void GpuTimer_Mark(const char* name) {
    if (!enabled || cur < 0 || !open) return;
    Slot& s = slots[cur];
    if (s.used >= kMarks) return;
    GpuTimer_CpuMark(name);
    rlDrawRenderBatchActive(); // raylib's batched draws belong to the pass that made them
    glad_glEndQuery(kTimeElapsed);
    s.name[s.used] = name;          // this segment ends at the mark: charged to its name
    ++s.used;
    s.name[s.used] = "(between frames)";
    glad_glBeginQuery(kTimeElapsed, s.q[s.used]); // the next one runs until the next mark (or frame)
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
