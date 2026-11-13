#include "terraingpu.h"
#include "gpucompute.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <rlgl.h>
#include <raylib.h>
#include <algorithm>
#include <cstdlib>
#include <deque>
#include <mutex>

namespace {
    struct BlockGpu { float o[4], n[4], t1[4], t2[4], s[4]; int a[4], b[4]; }; // terrain_density.comp Block
    static_assert(sizeof(BlockGpu) == 112, "std430 Block");

    constexpr int kMaxBlocks = 32;          // per batch
    constexpr int kSlots = 2;               // batches in flight

    struct Slot {
        ComputeBuffer grid, islands, out, blocks;
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
        logger.logf("[TerrainGpu] density blocks on the GPU\n");
        return true;
    }

    void SetField(const PlanetField& pf) {
        S.prog.use();
        if (S.field == &pf) return;
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
    void Dispatch(Slot& sl) {
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

bool TerrainGpu_Available() { return Init(); }

void TerrainGpu_Submit(GpuBlockJob&& job) {
    std::lock_guard<std::mutex> lock(S.m);
    S.queue.push_back(std::move(job));
}

int TerrainGpu_Pending() {
    std::lock_guard<std::mutex> lock(S.m);
    int n = (int)S.queue.size();
    for (const Slot& sl : S.slots) n += (int)sl.jobs.size();
    return n;
}

void TerrainGpu_Pump() {
    if (!S.ok) return;
    // Finished batches: hand out their samples.
    for (Slot& sl : S.slots) {
        if (!sl.fence || !Done(sl, 0)) continue;
        glDeleteSync(sl.fence);
        sl.fence = nullptr;
        auto all = std::make_shared<std::vector<float>>(sl.outCount);
        sl.out.download(all->data(), (unsigned)(sl.outCount * 4));
        S.latencyMs += (GetTime() - sl.t0) * 1000.0;
        for (size_t k = 0; k < sl.jobs.size(); ++k) {
            GpuBlockJob& j = sl.jobs[k];
            if (j.wanted && !j.wanted()) continue;
            j.done(all, sl.outOffset[k]); // the worker reads its slice: no per-block copy on the main thread
        }
        sl.jobs.clear();
    }
    // Free slots: the next batch (blocks of one field, still wanted).
    for (Slot& sl : S.slots) {
        if (sl.fence || !sl.jobs.empty()) continue;
        {
            std::lock_guard<std::mutex> lock(S.m);
            const PlanetField* f = nullptr;
            for (auto it = S.queue.begin(); it != S.queue.end() && (int)sl.jobs.size() < kMaxBlocks;) {
                if (it->wanted && !it->wanted()) { it = S.queue.erase(it); continue; }
                if (f && it->field != f) { ++it; continue; }
                f = it->field;
                sl.jobs.push_back(std::move(*it));
                it = S.queue.erase(it);
            }
        }
        if (sl.jobs.empty()) break;
        Dispatch(sl);
    }
}

bool TerrainGpu_SampleBlock(const PlanetField& pf, Vector3 o, float stride, int SS, float fp,
                            const PlanetField::BlockGrid& g, const std::vector<float>& islands, float* out) {
    if (!Init()) return false;
    // Wait for anything in flight, then use slot 0 synchronously.
    for (Slot& sl : S.slots) if (sl.fence) { while (!Done(sl, 1000000000ull)) {} }
    TerrainGpu_Pump(); // hands out the finished batches
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
        sl.grid.destroy(); sl.islands.destroy(); sl.out.destroy(); sl.blocks.destroy();
    }
    S.prog.unload();
    std::lock_guard<std::mutex> lock(S.m);
    S.queue.clear();
    S.tried = S.ok = false;
    S.field = nullptr;
}
