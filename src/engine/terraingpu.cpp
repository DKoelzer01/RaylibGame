#include "terraingpu.h"
#include "glcaps.h"
#include "gputimer.h"
#include "gpucompute.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <rlgl.h>
#include <raylib.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <atomic>
#include <memory>
#include <mutex>

namespace {
    struct BlockGpu { float o[4], n[4], t1[4], t2[4], s[4]; int a[4], b[4]; }; // terrain_density.comp Block
    static_assert(sizeof(BlockGpu) == 112, "std430 Block");

    constexpr int kMaxBlocks = 32;          // per batch
    constexpr int kSlots = 2;               // batches in flight
    constexpr int kMaxBackgroundBlocks = 8; // per batch of background blocks while throttled (#296)
    float gBackgroundFrom = 3.0e38f;
    bool gBackgroundThrottle = false;

    // Readback (#267): the samples are copied on the GPU into a persistently mapped buffer, and the workers read them
    // there; the main thread no longer copies a multi-MB batch out with glGetBufferSubData (~0.7 ms a frame while
    // streaming). A ring of them, each busy until the last worker drops its slice.
    struct Readback { GLuint buf = 0; float* ptr = nullptr; size_t cap = 0; std::atomic<bool> busy{ false }; };
    constexpr int kRing = 6;
    Readback ring[kRing];
    bool ringOk = false;

    // A free ring buffer of at least `bytes` (re-created bigger if needed), marked busy; -1 if none is free.
    int AcquireRing(size_t bytes) {
        for (int i = 0; i < kRing; ++i) {
            Readback& r = ring[i];
            if (r.busy.load()) continue;
            if (r.cap < bytes) {
                if (r.buf) glDeleteBuffers(1, &r.buf);
                size_t cap = 1u << 20;
                while (cap < bytes) cap <<= 1;
                glGenBuffers(1, &r.buf);
                glBindBuffer(GL_COPY_WRITE_BUFFER, r.buf);
                const GLbitfield flags = GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
                glBufferStorage(GL_COPY_WRITE_BUFFER, (GLsizeiptr)cap, nullptr, flags | GL_CLIENT_STORAGE_BIT);
                r.ptr = (float*)glMapBufferRange(GL_COPY_WRITE_BUFFER, 0, (GLsizeiptr)cap, flags);
                glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
                if (!r.ptr) { glDeleteBuffers(1, &r.buf); r.buf = 0; r.cap = 0; ringOk = false; logger.logf("[TerrainGpu] persistent map failed: readback by copy\n"); return -1; }
                r.cap = cap;
            }
            r.busy.store(true);
            return i;
        }
        return -1;
    }

    struct Slot {
        ComputeBuffer grid, islands, out, blocks;
        int rb = -1;                        // the ring buffer this batch is copied into (-1: read back by copy)
        GLsync fence = nullptr;
        std::vector<GpuBlockJob> jobs;
        std::vector<size_t> outOffset;      // per job (samples)
        size_t outCount = 0;
        double t0 = 0.0;
    };

    struct State {
        bool tried = false, ok = false;
        ComputeProgram prog;
        Slot slots[kSlots];
        std::mutex m;
        std::deque<GpuBlockJob> queue;
        const PlanetField* field = nullptr;   // whose uniforms are set
        const PlanetField* verified = nullptr; // checked against the CPU (its first batch)
        bool disabled = false;                 // a GPU / CPU mismatch turned the GPU path off
        long long blocks = 0, batches = 0;
        double latencyMs = 0.0;
    } S;

    void Ensure(ComputeBuffer& b, size_t bytes) {
        if (b.id && b.bytes >= bytes) return;
        size_t cap = std::max<size_t>(b.bytes, 4096);
        while (cap < bytes) cap *= 2;
        b.create((unsigned)cap, nullptr, true);
    }

    bool Init() {
        if (S.tried) return S.ok;
        S.tried = true;
        const char* env = getenv("RAYLIBGAME_GPU_TERRAIN");
        if (env && env[0] == '0') { logger.logf("[TerrainGpu] off (RAYLIBGAME_GPU_TERRAIN=0)\n"); return false; }
        if (!Compute_Available() || !glad_glFenceSync || !S.prog.load("resources/compute/terrain_density.comp")) {
            logger.logf("[TerrainGpu] unavailable: density blocks on the CPU\n");
            return false;
        }
        S.ok = true;
        // RAYLIBGAME_GPU_READBACK_COPY=1: the old copy readback (A/B).
        const char* rbEnv = getenv("RAYLIBGAME_GPU_READBACK_COPY");
        ringOk = GL_Caps().bufferStorage && glad_glCopyBufferSubData && !(rbEnv && rbEnv[0] == '1');
        logger.logf("[TerrainGpu] readback: %s\n", ringOk ? "persistent mapped ring" : "copy");
        logger.logf("[TerrainGpu] density blocks on the GPU\n");
        return true;
    }

    void SetField(const PlanetField& pf) {
        S.prog.use();
        // Uploaded on every dispatch (22 vec4s): a cache keyed on the field pointer went stale when a reloaded
        // world's new field landed at the old one's address (PR #179 review).
        std::vector<float> params;
        pf.gpuNoiseParams(params);
        rlSetUniform(S.prog.uniform("noiseParams"), params.data(), RL_SHADER_UNIFORM_VEC4, (int)params.size() / 4);
        S.prog.setFloat("planetRadius", pf.radius);
        S.prog.setInt("planetMoon", pf.isMoon() ? 1 : 0);
        S.prog.setFloat("mantleRadius", pf.mantleRadiusValue());
        S.prog.setFloat("outerR", pf.radius + pf.maxElevationValue());
        S.prog.setFloat("innerR", pf.radius + pf.minElevationValue());
        S.field = &pf;
    }

    // Uploads the slot's jobs (all of one field) and dispatches them; sets the fence.
    void Dispatch(Slot& sl, bool readback = false) { // readback: into the mapped ring when it can (the streaming pump)
        const PlanetField& pf = *sl.jobs.front().field;
        SetField(pf);
        std::vector<BlockGpu> blocks;
        std::vector<float> grid, isl;
        sl.outOffset.clear();
        size_t out = 0;
        int maxN = 0;
        for (const GpuBlockJob& j : sl.jobs) {
            const PlanetField::BlockGrid& g = j.grid;
            BlockGpu b{};
            b.o[0] = j.block.origin.x; b.o[1] = j.block.origin.y; b.o[2] = j.block.origin.z; b.o[3] = j.block.stride;
            b.n[0] = g.n.x; b.n[1] = g.n.y; b.n[2] = g.n.z; b.n[3] = j.block.footprint;
            b.t1[0] = g.t1.x; b.t1[1] = g.t1.y; b.t1[2] = g.t1.z; b.t1[3] = g.umin;
            b.t2[0] = g.t2.x; b.t2[1] = g.t2.y; b.t2[2] = g.t2.z; b.t2[3] = g.vmin;
            b.s[0] = g.su; b.s[1] = g.sv;
            b.a[0] = j.block.E; b.a[1] = g.nu; b.a[2] = g.nv; b.a[3] = (int)(grid.size() / 4);
            b.b[0] = (int)(isl.size() / 12); b.b[1] = (int)(j.islands.size() / 12); b.b[2] = (int)out;
            blocks.push_back(b);
            for (size_t k = 0; k < g.h.size(); ++k) { const float v[4] = { g.h[k], g.cols[k].mask, g.cols[k].f, 0.0f }; grid.insert(grid.end(), v, v + 4); }
            isl.insert(isl.end(), j.islands.begin(), j.islands.end());
            sl.outOffset.push_back(out);
            const int n = j.block.E * j.block.E * j.block.E;
            out += (size_t)n;
            maxN = std::max(maxN, n);
        }
        if (isl.empty()) isl.assign(12, 0.0f);
        sl.outCount = out;
        Ensure(sl.grid, grid.size() * 4);
        Ensure(sl.islands, isl.size() * 4);
        Ensure(sl.out, out * 4);
        Ensure(sl.blocks, blocks.size() * sizeof(BlockGpu));
        sl.grid.upload(grid.data(), (unsigned)(grid.size() * 4));
        sl.islands.upload(isl.data(), (unsigned)(isl.size() * 4));
        sl.blocks.upload(blocks.data(), (unsigned)(blocks.size() * sizeof(BlockGpu)));
        sl.grid.bind(0); sl.islands.bind(1); sl.out.bind(2); sl.blocks.bind(3);
        S.prog.use();
        rlComputeShaderDispatch((unsigned)((maxN + 63) / 64), (unsigned)blocks.size(), 1);
        glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
        sl.rb = readback && ringOk ? AcquireRing(out * 4) : -1;
        if (sl.rb >= 0) { // the samples into the mapped ring buffer, on the GPU
            glBindBuffer(GL_COPY_READ_BUFFER, sl.out.id);
            glBindBuffer(GL_COPY_WRITE_BUFFER, ring[sl.rb].buf);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, (GLsizeiptr)(out * 4));
            glBindBuffer(GL_COPY_READ_BUFFER, 0);
            glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
            // (coherent mapping: the fence wait makes the copy visible; no client-mapped barrier needed)
        }
        sl.fence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
        glFlush();
        sl.t0 = GetTime();
        ++S.batches;
        S.blocks += (long long)sl.jobs.size();
    }

    bool Done(Slot& sl, GLuint64 timeoutNs) {
        if (!sl.fence) return false;
        const GLenum r = glClientWaitSync(sl.fence, timeoutNs ? GL_SYNC_FLUSH_COMMANDS_BIT : 0, timeoutNs);
        return r == GL_ALREADY_SIGNALED || r == GL_CONDITION_SATISFIED;
    }
}

bool TerrainGpu_Available() { return Init() && !S.disabled; }

void TerrainGpu_Submit(GpuBlockJob&& job) {
    std::lock_guard<std::mutex> lock(S.m);
    // In priority order (lower first, FIFO among equals): a dig's remesh (negative) goes ahead of the
    // background LOD builds instead of waiting behind them.
    auto it = S.queue.end();
    while (it != S.queue.begin() && std::prev(it)->priority > job.priority) --it;
    S.queue.insert(it, std::move(job));
}

int TerrainGpu_Pending() {
    std::lock_guard<std::mutex> lock(S.m);
    int n = (int)S.queue.size();
    for (const Slot& sl : S.slots) n += (int)sl.jobs.size();
    return n;
}

void TerrainGpu_SetBackground(float from, bool throttle) { gBackgroundFrom = from; gBackgroundThrottle = throttle; }

void TerrainGpu_Pump() {
    if (!S.ok) return;
    // Finished batches: hand out their samples.
    for (Slot& sl : S.slots) {
        if (!sl.fence || !Done(sl, 0)) continue;
        glDeleteSync(sl.fence);
        sl.fence = nullptr;
        // A CPU buffer from a pool (a fresh multi-MB zero-filled vector per batch was allocated on the main thread);
        // it comes back when the last worker drops its slice.
        static std::mutex poolM;
        static std::vector<std::vector<float>*> pool;
        auto pooled = [&]() {
            std::vector<float>* buf = nullptr;
            { std::lock_guard<std::mutex> l(poolM); if (!pool.empty()) { buf = pool.back(); pool.pop_back(); } }
            if (!buf) buf = new std::vector<float>();
            if (buf->size() < sl.outCount) buf->resize(sl.outCount);
            std::shared_ptr<std::vector<float>> v(buf, [](std::vector<float>* p) { std::lock_guard<std::mutex> l(poolM); pool.push_back(p); });
            return v;
        };
        // #267: the batch was copied on the GPU into a mapped ring buffer: the workers read it there (released when
        // the last slice is done). Otherwise it is read back by copy.
        std::shared_ptr<const float> all;
        if (sl.rb >= 0) {
            Readback* r = &ring[sl.rb];
            all = std::shared_ptr<const float>(r->ptr, [r](const float*) { r->busy.store(false); });
            sl.rb = -1;
        } else {
            auto v = pooled();
            sl.out.download(v->data(), (unsigned)(sl.outCount * 4));
            all = std::shared_ptr<const float>(v, v->data());
        }
        S.latencyMs += (GetTime() - sl.t0) * 1000.0;
        // The first batch of each field is checked against the CPU once (one block, ~4 ms): the GPU path rests
        // on bit-exact results (collision reads the CPU field), which another driver might not give (FMA,
        // denormals, sqrt). On any difference the GPU path turns off and this batch is redone on the CPU.
        if (!sl.jobs.empty() && S.verified != sl.jobs.front().field) {
            const GpuBlockJob& j = sl.jobs.front();
            const size_t n = (size_t)j.block.E * j.block.E * j.block.E;
            std::vector<float> cpu(n);
            j.field->sampleBlock(j.block.origin, j.block.stride, j.block.E, j.block.footprint, j.grid, cpu.data());
            if (memcmp(cpu.data(), all.get() + sl.outOffset[0], n * sizeof(float)) != 0) {
                logger.logf("[TerrainGpu] GPU density differs from the CPU on this driver: density blocks back on the CPU\n");
                S.disabled = true;
            } else {
                S.verified = j.field;
            }
        }
        if (S.disabled) { // redone on the CPU, into a writable buffer (the ring is mapped read-only)
            auto v = pooled();
            for (size_t k = 0; k < sl.jobs.size(); ++k) {
                GpuBlockJob& j = sl.jobs[k];
                if (j.wanted && !j.wanted()) continue;
                j.field->sampleBlock(j.block.origin, j.block.stride, j.block.E, j.block.footprint, j.grid, v->data() + sl.outOffset[k]);
            }
            all = std::shared_ptr<const float>(v, v->data());
        }
        for (size_t k = 0; k < sl.jobs.size(); ++k) {
            GpuBlockJob& j = sl.jobs[k];
            if (j.wanted && !j.wanted()) continue;
            j.done(all, sl.outOffset[k]); // the worker reads its slice: no per-block copy on the main thread
        }
        sl.jobs.clear();
    }
    // Free slots: the next batch (blocks of one field, still wanted). The queue is in priority order, so background
    // blocks (#296) only come after the rest; while throttled, one small batch of them a frame.
    bool backgroundSent = false;
    for (Slot& sl : S.slots) {
        if (sl.fence || !sl.jobs.empty()) continue;
        {
            std::lock_guard<std::mutex> lock(S.m);
            const PlanetField* f = nullptr;
            int background = 0;
            for (auto it = S.queue.begin(); it != S.queue.end() && (int)sl.jobs.size() < kMaxBlocks;) {
                if (it->wanted && !it->wanted()) { it = S.queue.erase(it); continue; }
                if (gBackgroundThrottle && it->priority >= gBackgroundFrom) { // the rest of the queue is background too
                    if (backgroundSent || background >= kMaxBackgroundBlocks) break;
                    ++background;
                }
                if (f && it->field != f) { ++it; continue; }
                f = it->field;
                sl.jobs.push_back(std::move(*it));
                it = S.queue.erase(it);
            }
            if (background > 0) backgroundSent = true;
        }
        if (sl.jobs.empty()) break;
        if (S.disabled) {
            // Turned off after a mismatch: whatever was still queued is sampled on the CPU (once; new builds
            // no longer come here, TerrainGpu_Available is false).
            for (GpuBlockJob& j : sl.jobs) {
                if (j.wanted && !j.wanted()) continue;
                const size_t n = (size_t)j.block.E * j.block.E * j.block.E;
                auto cpu = std::make_shared<std::vector<float>>(n);
                j.field->sampleBlock(j.block.origin, j.block.stride, j.block.E, j.block.footprint, j.grid, cpu->data());
                j.done(std::shared_ptr<const float>(cpu, cpu->data()), 0);
            }
            sl.jobs.clear();
            continue;
        }
        Dispatch(sl, true);
    }
}

bool TerrainGpu_SampleBlock(const PlanetField& pf, Vector3 o, float stride, int SS, float fp,
                            const PlanetField::BlockGrid& g, const std::vector<float>& islands, float* out) {
    if (!Init()) return false;
    // Tests only, with nothing streaming: refuse (the caller samples on the CPU) rather than share slot 0 with
    // queued or in-flight chunk builds.
    {
        std::lock_guard<std::mutex> lock(S.m);
        if (!S.queue.empty()) return false;
    }
    for (const Slot& s : S.slots) if (s.fence || !s.jobs.empty()) return false;
    Slot& sl = S.slots[0];
    GpuBlockJob j;
    j.field = &pf;
    j.block.origin = o; j.block.stride = stride; j.block.E = SS; j.block.footprint = fp;
    j.grid = g;
    j.islands = islands;
    sl.jobs.push_back(std::move(j));
    Dispatch(sl);
    while (!Done(sl, 1000000000ull)) {}
    glDeleteSync(sl.fence);
    sl.fence = nullptr;
    sl.out.download(out, (unsigned)(SS * SS * SS * 4));
    sl.jobs.clear();
    return true;
}

TerrainGpuStats TerrainGpu_Stats() {
    TerrainGpuStats s;
    s.blocks = S.blocks; s.batches = S.batches;
    s.avgLatencyMs = S.batches ? S.latencyMs / (double)S.batches : 0.0;
    return s;
}

void TerrainGpu_Unload() {
    for (Slot& sl : S.slots) {
        if (sl.fence) glDeleteSync(sl.fence);
        sl.fence = nullptr;
        sl.jobs.clear();
        sl.outOffset.clear();
        sl.grid.destroy(); sl.islands.destroy(); sl.out.destroy(); sl.blocks.destroy();
    }
    S.prog.unload();
    std::lock_guard<std::mutex> lock(S.m);
    S.queue.clear();
    S.tried = S.ok = false;
    S.field = S.verified = nullptr;
    S.disabled = false;
}
