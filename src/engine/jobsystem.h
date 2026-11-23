#ifndef JOBSYSTEM_H
#define JOBSYSTEM_H

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

// What a job is (#387: the per-second job mix in --perf-log).
enum class JobTag : uint8_t { Other = 0, ChunkPrep, ChunkMesh, RingCollect, WaterRebuild, Trees, Fluid, FarSphere, NavMap, Count };
const char* JobTagName(JobTag t);

// Fixed-size worker pool with a priority queue. Lower priority value runs first; equal
// priorities run in submission order. Jobs must not touch raylib/GL: workers produce CPU data and
// the main thread uploads it.
class JobSystem {
public:
    explicit JobSystem(unsigned threadCount = 0); // 0 = hardware threads - 1 (at least 1)
    ~JobSystem();                                 // drops queued jobs, waits for running ones

    JobSystem(const JobSystem&) = delete;
    JobSystem& operator=(const JobSystem&) = delete;

    // A job is background when `background` is set or its priority is at or above backgroundFrom (the preload
    // ring, #296); every other job is foreground and keeps foregroundIdle() false while queued or running.
    void submit(float priority, std::function<void()> fn, bool background = false, JobTag tag = JobTag::Other);
    // Tags the next submit from this thread (#387: lets a multi-line submit call stay as it is).
    static void tagNext(JobTag t);
    // Per tag: queued now, and submitted / finished since the start (#387).
    long long queuedOf(JobTag t) const { return queuedBy[(int)t].load(); }
    long long submittedOf(JobTag t) const { return submittedBy[(int)t].load(); }
    long long doneOf(JobTag t) const { return doneBy[(int)t].load(); }
    void setBackgroundFrom(float priority) { backgroundFrom = priority; }
    void clearQueued();                // discard jobs that have not started yet
    size_t queued() const;             // waiting in the queue
    int running() const { return active.load(); }
    bool idle() const { return queued() == 0 && running() == 0; }
    bool foregroundIdle() const { return foreground.load() == 0; } // no foreground job queued or running (#296 review)
    unsigned threadCount() const { return (unsigned)workers.size(); }
    uint64_t done() const { return finished.load(); } // jobs finished so far (a loading bar's work done, #494)

private:
    struct Job {
        float priority;
        uint64_t seq;
        std::function<void()> fn;
        bool foreground = true;
        JobTag tag = JobTag::Other;
        bool operator<(const Job& o) const { // std::priority_queue is a max-heap
            if (priority != o.priority) return priority > o.priority;
            return seq > o.seq;
        }
    };

    void workerLoop();

    std::vector<std::thread> workers;
    std::priority_queue<Job> jobs;
    mutable std::mutex mutex;
    std::condition_variable cv;
    std::atomic<int> active{0};
    std::atomic<int> foreground{0};    // foreground jobs queued + running
    std::atomic<uint64_t> finished{0};
    std::atomic<float> backgroundFrom{3.0e38f};
    std::atomic<long long> queuedBy[(int)JobTag::Count]{}, submittedBy[(int)JobTag::Count]{}, doneBy[(int)JobTag::Count]{};
    uint64_t nextSeq = 0;
    bool stopping = false;
    bool lowPriority = true;           // #387: workers below normal priority (RAYLIBGAME_JOB_PRIORITY=normal: off)
};

// Runs fn(0) .. fn(parts - 1), spread over the pool and the calling thread, and returns when all are done. Safe to call
// from inside a pool job (#1671): the caller and the helper jobs pull parts from one counter, so the caller does whatever
// no worker picks up - with a single worker busy running the caller, or after clearQueued dropped the helpers, it simply
// runs every part itself instead of waiting forever. A helper that starts late finds nothing left and touches only the
// shared counters (never the caller's stack). Parts must write disjoint data; which thread runs a part changes nothing.
inline void ParallelParts(JobSystem& js, int parts, const std::function<void(int)>& fn, float priority, JobTag tag) {
    if (parts <= 1) { if (parts == 1) fn(0); return; }
    struct Shared { std::atomic<int> next{ 0 }, done{ 0 }; const std::function<void(int)>* fn = nullptr; };
    auto sh = std::make_shared<Shared>();
    sh->fn = &fn;
    auto pull = [](Shared& s, int n) { for (int p; (p = s.next.fetch_add(1)) < n;) { (*s.fn)(p); s.done.fetch_add(1); } };
    for (int k = 1; k < parts; ++k) js.submit(priority, [sh, parts, pull] { pull(*sh, parts); }, false, tag);
    pull(*sh, parts);
    while (sh->done.load() < parts) std::this_thread::yield(); // only parts a worker is running right now remain
}

#endif
