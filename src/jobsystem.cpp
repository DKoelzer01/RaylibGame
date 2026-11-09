#include "jobsystem.h"

#include <algorithm>

JobSystem::JobSystem(unsigned threadCount) {
    if (threadCount == 0) {
        const unsigned hw = std::thread::hardware_concurrency();
        threadCount = std::max(1u, hw > 1 ? hw - 1 : 1u); // leave a core for the main/render thread
    }
    workers.reserve(threadCount);
    for (unsigned i = 0; i < threadCount; ++i) workers.emplace_back([this] { workerLoop(); });
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

void JobSystem::submit(float priority, std::function<void()> fn) {
    {
        std::lock_guard<std::mutex> lock(mutex);
        jobs.push(Job{ priority, nextSeq++, std::move(fn) });
    }
    cv.notify_one();
}

void JobSystem::clearQueued() {
    std::lock_guard<std::mutex> lock(mutex);
    jobs = {};
}

size_t JobSystem::queued() const {
    std::lock_guard<std::mutex> lock(mutex);
    return jobs.size();
}

void JobSystem::workerLoop() {
    for (;;) {
        std::function<void()> fn;
        {
            std::unique_lock<std::mutex> lock(mutex);
            cv.wait(lock, [this] { return stopping || !jobs.empty(); });
            if (stopping) return;
            fn = std::move(const_cast<Job&>(jobs.top()).fn);
            jobs.pop();
            active++; // under the lock, so idle() never sees queued==0 && running==0 mid-handoff
        }
        fn();
        active--;
    }
}
