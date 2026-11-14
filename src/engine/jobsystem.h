#ifndef JOBSYSTEM_H
#define JOBSYSTEM_H

#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

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
    void submit(float priority, std::function<void()> fn, bool background = false);
    void setBackgroundFrom(float priority) { backgroundFrom = priority; }
    void clearQueued();                // discard jobs that have not started yet
    size_t queued() const;             // waiting in the queue
    int running() const { return active.load(); }
    bool idle() const { return queued() == 0 && running() == 0; }
    bool foregroundIdle() const { return foreground.load() == 0; } // no foreground job queued or running (#296 review)
    unsigned threadCount() const { return (unsigned)workers.size(); }

private:
    struct Job {
        float priority;
        uint64_t seq;
        std::function<void()> fn;
        bool foreground = true;
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
    std::atomic<float> backgroundFrom{3.0e38f};
    uint64_t nextSeq = 0;
    bool stopping = false;
};

#endif
