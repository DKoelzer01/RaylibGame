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

const char* JobTagName(JobTag t) {
    static const char* names[(int)JobTag::Count] = { "other", "chunk-prep", "chunk-mesh", "ring-collect", "water", "trees", "fluid", "farsphere", "navmap" };
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
        active--;
        doneBy[(int)tag]++;
        if (fg) --foreground; // after the job, so a job that submits follow-ups keeps the count above zero
    }
}
