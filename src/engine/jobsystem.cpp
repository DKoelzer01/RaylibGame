#include "jobsystem.h"
#include "utils/logger.h"
#include <cstdio>

#include <algorithm>
#include <cstdlib>
#ifdef _WIN32
// Not <windows.h> (its names clash with raylib's): the two calls the workers need (#387).
extern "C" __declspec(dllimport) void* __stdcall GetCurrentThread();
extern "C" __declspec(dllimport) int __stdcall SetThreadPriority(void* thread, int priority);
#else
#include <sys/resource.h>
#include <unistd.h>
#endif

JobSystem::JobSystem(unsigned threadCount) {
    // #387: the workers run below normal priority, so a busy machine (compiles, browsers) takes their time, not the main
    // thread's (under a 12-way compile: 17 -> 11 ms frames while streaming, worst 39 -> 11 ms after; spawn preload not slower).
    // A/B: RAYLIBGAME_JOB_PRIORITY=normal, RAYLIBGAME_JOB_THREADS=N workers.
    if (const char* e = getenv("RAYLIBGAME_JOB_THREADS"); e && threadCount == 0) threadCount = (unsigned)std::max(1, atoi(e));
    if (const char* e = getenv("RAYLIBGAME_JOB_PRIORITY")) lowPriority = e[0] != 'n'; // default low (#387); =normal for A/B
    if (threadCount == 0) {
        const unsigned hw = std::thread::hardware_concurrency();
        threadCount = std::max(1u, hw > 1 ? hw - 1 : 1u); // leave a core for the main/render thread
    }
    workers.reserve(threadCount);
    for (unsigned i = 0; i < threadCount; ++i) workers.emplace_back([this, i] { char name[16]; snprintf(name, sizeof(name), "job%u", i); Log_SetThreadName(name); workerLoop(); }); // named in log lines (#514)
}

JobSystem::~JobSystem() {
    {
        std::lock_guard<std::mutex> lock(mutex);
        stopping = true;
        jobs = {};
    }
    cv.notify_all();
    for (std::thread& t : workers) t.join();
}

const char* JobTagName(JobTag t) {
    static const char* names[(int)JobTag::Count] = { "other", "chunk-prep", "chunk-mesh", "ring-collect", "water", "trees", "fluid", "farsphere", "navmap", "terrain-fold" };
    return (int)t < (int)JobTag::Count ? names[(int)t] : "?";
}

static thread_local JobTag tNextTag = JobTag::Other;
void JobSystem::tagNext(JobTag t) { tNextTag = t; }

void JobSystem::submit(float priority, std::function<void()> fn, bool background, JobTag tag) {
    if (tag == JobTag::Other) tag = tNextTag;
    tNextTag = JobTag::Other;
    const bool fg = !background && priority < backgroundFrom.load();
    {
        std::lock_guard<std::mutex> lock(mutex);
        if (fg) ++foreground; // under the lock with the push, so clearQueued's count stays in step
        jobs.push(Job{ priority, nextSeq++, std::move(fn), fg, tag });
        queuedBy[(int)tag]++; submittedBy[(int)tag]++;
    }
    cv.notify_one();
}

void JobSystem::clearQueued() {
    std::lock_guard<std::mutex> lock(mutex);
    while (!jobs.empty()) { if (jobs.top().foreground) --foreground; queuedBy[(int)jobs.top().tag]--; jobs.pop(); }
}

size_t JobSystem::queued() const {
    std::lock_guard<std::mutex> lock(mutex);
    return jobs.size();
}

void JobSystem::workerLoop() {
#ifdef _WIN32
    if (lowPriority) SetThreadPriority(GetCurrentThread(), -1); // THREAD_PRIORITY_BELOW_NORMAL
#elif defined(__linux__)
    if (lowPriority) (void)setpriority(PRIO_PROCESS, (id_t)gettid(), 5); // Linux: a thread's nice value is its own
#endif
    for (;;) {
        std::function<void()> fn;
        bool fg = false;
        JobTag tag = JobTag::Other;
        {
            std::unique_lock<std::mutex> lock(mutex);
            cv.wait(lock, [this] { return stopping || !jobs.empty(); });
            if (stopping) return;
            fn = std::move(const_cast<Job&>(jobs.top()).fn);
            fg = jobs.top().foreground;
            tag = jobs.top().tag;
            queuedBy[(int)tag]--;
            jobs.pop();
            active++; // under the lock, so idle() never sees queued==0 && running==0 mid-handoff
        }
        fn();
        ++finished; // before active drops: done + queued + running never dips mid-handoff
        active--;
        doneBy[(int)tag]++;
        if (fg) --foreground; // after the job, so a job that submits follow-ups keeps the count above zero
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
#include <chrono>

// (#1671) A pool job that fans out into parts must finish with a single worker (that worker is the job itself): the old
// submit-and-spin pattern waited forever for helpers no thread could run.
TEST(JobSystem, ParallelPartsFromInsideAJobWithOneWorker) {
    JobSystem js(1);
    std::atomic<int> sum{ 0 };
    std::atomic<bool> finished{ false };
    js.submit(0.0f, [&] {
        ParallelParts(js, 8, [&](int p) { sum.fetch_add(p + 1); }, 0.0f, JobTag::Other);
        finished = true;
    });
    const auto t0 = std::chrono::steady_clock::now();
    while (!finished.load() && std::chrono::steady_clock::now() - t0 < std::chrono::seconds(5)) std::this_thread::yield();
    EXPECT_TRUE(finished.load());
    EXPECT_EQ(sum.load(), 36); // every part ran exactly once
}

// Helpers dropped before they start (unloadWorld's clearQueued): the caller runs their parts.
TEST(JobSystem, ParallelPartsWhenHelpersAreDropped) {
    JobSystem js(1);
    std::atomic<bool> release{ false };
    js.submit(0.0f, [&] { while (!release.load()) std::this_thread::yield(); }); // keeps the only worker busy
    std::atomic<int> runs{ 0 };
    ParallelParts(js, 4, [&](int) { runs.fetch_add(1); }, 0.0f, JobTag::Other); // its helpers queue behind the busy job
    js.clearQueued();
    release = true;
    EXPECT_EQ(runs.load(), 4);
}
#endif
