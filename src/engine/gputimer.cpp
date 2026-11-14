#include "gputimer.h"

#include <rlgl.h>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include "utils/logger.h"

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

namespace {
    // Slow frames (benches): a frame's CPU marks, printed when the frame took over 20 ms (at most 10 of them).
    std::vector<std::pair<const char*, double>> frameMarks;
    double frameMs = 0.0;
    int slowPrinted = 0;
    const auto gRunStart = std::chrono::steady_clock::now(); // program start (static init)
}

void GpuTimer_CpuFrameStart() {
    if (enabled && cpuStarted && frameMs > 20.0 && slowPrinted < 10) {
        ++slowPrinted;
        printf("bench: slow frame %.1f ms (cpu) at %.1f s, frame %ld:", frameMs, std::chrono::duration<double>(std::chrono::steady_clock::now() - gRunStart).count(), cpuFrames);
        for (const auto& m : frameMarks) if (m.second > 1.0) printf(" %s %.1f", m.first, m.second);
        printf("\n");
    }
    frameMarks.clear(); frameMs = 0.0;
    if (!enabled) { cpuStarted = false; return; }
    cpuLast = std::chrono::steady_clock::now();
    cpuStarted = true;
    ++cpuFrames;
}

void GpuTimer_CpuMark(const char* name) {
    if (!enabled || !cpuStarted) return;
    const auto now = std::chrono::steady_clock::now();
    const double ms = std::chrono::duration<double, std::milli>(now - cpuLast).count();
    AddTo(cpuSums, name, ms);
    frameMarks.push_back({ name, ms }); frameMs += ms;
    cpuLast = now;
}

void GpuTimer_PrintCpu() {
    if (cpuFrames == 0) return;
    double total = 0.0;
    printf("bench: cpu per frame (avg of %ld):", cpuFrames);
    for (const Sum& x : cpuSums) { printf(" %s %.2f", x.name, x.ms / cpuFrames); total += x.ms; }
    printf(" | total %.2f ms\n", total / cpuFrames);
}

void GpuTimer_LogAndReset(const char* tag) {
    // In-game perf log (#204): the same per-pass averages as the bench lines, to log.txt, then start over.
    std::string g = "[Perf] ", c = "[Perf] ";
    char buf[128];
    snprintf(buf, sizeof(buf), "%s gpu per frame (%ld):", tag, frames); g += buf;
    for (const Sum& x : sums) { snprintf(buf, sizeof(buf), " %s %.2f", x.name, frames ? x.ms / frames : 0.0); g += buf; }
    snprintf(buf, sizeof(buf), " | total %.2f ms\n", frames ? totalMs / frames : 0.0); g += buf;
    double total = 0.0;
    snprintf(buf, sizeof(buf), "%s cpu per frame (%ld):", tag, cpuFrames); c += buf;
    for (const Sum& x : cpuSums) { snprintf(buf, sizeof(buf), " %s %.2f", x.name, cpuFrames ? x.ms / cpuFrames : 0.0); c += buf; total += x.ms; }
    snprintf(buf, sizeof(buf), " | total %.2f ms\n", cpuFrames ? total / cpuFrames : 0.0); c += buf;
    logger.logf("%s", g.c_str());
    logger.logf("%s", c.c_str());
    sums.clear(); frames = 0; totalMs = 0.0;
    cpuSums.clear(); cpuFrames = 0;
}
