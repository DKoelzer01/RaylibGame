#include "gpucompute.h"
#include "cavekinds.h"
#include "shaderasync.h"
#include "utils/logger.h"

#include <raylib.h>
#include <rlgl.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <random>
#include <unordered_map>
#include <utility>
#include <vector>

extern "C" {
    typedef void (*PfnMemoryBarrier)(unsigned int);
    extern PfnMemoryBarrier glad_glMemoryBarrier;
    typedef void (*PfnBindBuffer)(unsigned int, unsigned int);
    extern PfnBindBuffer glad_glBindBuffer;
    typedef void (*PfnDispatchComputeIndirect)(long long);
    extern PfnDispatchComputeIndirect glad_glDispatchComputeIndirect;
    // (#1769) program binaries (GL 4.1 core)
    typedef void (*PfnGetProgramiv)(unsigned, unsigned, int*);
    extern PfnGetProgramiv glad_glGetProgramiv;
    typedef void (*PfnGetProgramBinary)(unsigned, int, int*, unsigned*, void*);
    extern PfnGetProgramBinary glad_glGetProgramBinary;
    typedef void (*PfnProgramBinary)(unsigned, unsigned, const void*, int);
    extern PfnProgramBinary glad_glProgramBinary;
    typedef unsigned (*PfnCreateProgram)(void);
    extern PfnCreateProgram glad_glCreateProgram;
    typedef void (*PfnDeleteProgram)(unsigned);
    extern PfnDeleteProgram glad_glDeleteProgram;
    typedef void (*PfnGetIntegerv)(unsigned, int*);
    extern PfnGetIntegerv glad_glGetIntegerv;
    typedef const unsigned char* (*PfnGetString)(unsigned);
    extern PfnGetString glad_glGetString;
}
#if RAYLIBGAME_HEADLESS_PLATFORM
extern "C" int RaylibHeadless_GLMode(void); // src/platform/rcore_headless.c: 2 = the null GL (no GPU)
#endif

bool Compute_Available() {
#if RAYLIBGAME_HEADLESS_PLATFORM
    // A Linux server without a GPU (#792) runs on stub GL functions: every GPU path (fluid, terrain density, arena,
    // tree culling) takes its CPU fallback.
    if (RaylibHeadless_GLMode() == 2) return false;
#endif
    return rlGetVersion() == RL_OPENGL_43 && glad_glMemoryBarrier != nullptr;
}

std::string LoadShaderSource(const std::string& path, std::vector<std::string>& includes) {
    struct Rec {
        static std::string Load(const std::string& p, std::vector<std::string>& inc, int depth) {
            char* text = LoadFileText(p.c_str());
            if (!text) return std::string();
            std::string src(text), out;
            UnloadFileText(text);
            const std::string dir = p.substr(0, p.find_last_of("/\\") + 1);
            size_t pos = 0;
            while (pos < src.size()) {
                size_t eol = src.find('\n', pos);
                if (eol == std::string::npos) eol = src.size();
                const std::string line = src.substr(pos, eol - pos);
                const size_t at = line.find("#include");
                const size_t q0 = line.find('"'), q1 = line.rfind('"');
                if (at != std::string::npos && line.find_first_not_of(" \t") == at && q0 != std::string::npos && q1 > q0 && depth < 4) {
                    const std::string file = dir + line.substr(q0 + 1, q1 - q0 - 1);
                    inc.push_back(file);
                    out += Load(file, inc, depth + 1);
                } else if (line.find("#pragma cave_kinds") != std::string::npos && line.find_first_not_of(" \t") == line.find('#')) {
                    out += CaveKinds_Glsl(); // (#1465) the cave segment kinds' table, from cavekinds.h
                } else {
                    out += line;
                }
                out += "\n";
                pos = eol + 1;
            }
            return out;
        }
    };
    return Rec::Load(path, includes, 0);
}

namespace {
// The source as compiled (includes expanded); empty if missing (logged).
std::string ComputeSource(const std::string& p) {
    std::vector<std::string> inc;
    std::string src = LoadShaderSource(p, inc);
    if (src.empty()) { LOG_ERROR("Compute", "compute shader %s: not found\n", p.c_str()); return src; }
    // (#1047) RAYLIBGAME_COLD_SHADERS=1: a comment unique to this run first, so the driver's shader cache misses, as on a fresh
    // machine or CI runner (comments may come before #version)
    static const bool cold = getenv("RAYLIBGAME_COLD_SHADERS") && getenv("RAYLIBGAME_COLD_SHADERS")[0] == '1';
    static int coldN = 0;
    if (cold) src = std::string(TextFormat("// cold %lld %d\n", (long long)time(nullptr), ++coldN)) + src;
    return src;
}
// (#1047) Compiles started early by Compute_Prewarm, by path; load() takes its program from here.
std::unordered_map<std::string, std::pair<AsyncProgram, double>> gPending; // the program, when it was started

// (#1769) The linked program's binary on disk, keyed by the expanded source and the driver: a cold driver cache made every
// scenario process (and the runner, and a fresh install) compile fluid_step.comp (~8 s) and terrain_density.comp (~5 s) again.
// A binary the driver refuses (a driver update) is deleted and the shader compiled. RAYLIBGAME_SHADER_CACHE = a directory, or
// "0" for none; default cache/shaders. RAYLIBGAME_COLD_SHADERS runs (#1047) neither read nor write it.
namespace BinCache {
constexpr unsigned kLinkStatus = 0x8B82, kBinaryLength = 0x8741, kNumBinaryFormats = 0x87FE;
constexpr unsigned kVendor = 0x1F00, kRenderer = 0x1F01, kVersion = 0x1F02;
constexpr uint32_t kMagic = 0x32424352; // "RCB2": magic, format, length, checksum, then the binary

std::string Dir() {
    static const std::string d = [] {
        const char* e = getenv("RAYLIBGAME_SHADER_CACHE");
        if (e && (!strcmp(e, "0") || !strcmp(e, "off") || !strcmp(e, "false"))) return std::string();
        if (getenv("RAYLIBGAME_COLD_SHADERS") && getenv("RAYLIBGAME_COLD_SHADERS")[0] == '1') return std::string();
        return std::string(e && e[0] ? e : "cache/shaders");
    }();
    return d;
}
bool Usable() {
    if (Dir().empty() || !glad_glGetProgramBinary || !glad_glProgramBinary || !glad_glGetProgramiv || !glad_glCreateProgram) return false;
    static const int formats = [] { int n = 0; if (glad_glGetIntegerv) glad_glGetIntegerv(kNumBinaryFormats, &n); return n; }();
    return formats > 0;
}
uint64_t Fnv(uint64_t h, const std::string& s) { for (unsigned char c : s) { h ^= c; h *= 1099511628211ull; } return h; }
uint32_t Checksum(const std::vector<char>& b) { uint64_t h = 14695981039346656037ull; for (char c : b) { h ^= (unsigned char)c; h *= 1099511628211ull; } return (uint32_t)(h ^ (h >> 32)); }
std::string GlString(unsigned name) { const unsigned char* s = glad_glGetString ? glad_glGetString(name) : nullptr; return s ? (const char*)s : ""; }
// The file for a shader: its name and a hash of the driver and the source as compiled.
std::string PathFor(const std::string& shaderPath, const std::string& src) {
    static const std::string driver = GlString(kVendor) + "|" + GlString(kRenderer) + "|" + GlString(kVersion);
    const uint64_t h = Fnv(Fnv(14695981039346656037ull, driver), src);
    const size_t slash = shaderPath.find_last_of("/\\");
    return Dir() + "/" + shaderPath.substr(slash == std::string::npos ? 0 : slash + 1) + "-" + TextFormat("%016llx", (unsigned long long)h) + ".bin";
}
unsigned Load(const std::string& file) {
    FILE* f = fopen(file.c_str(), "rb");
    if (!f) return 0;
    uint32_t head[4] = {}; // magic, format, length, checksum
    std::vector<char> bin;
    bool ok = fread(head, sizeof head, 1, f) == 1 && head[0] == kMagic && head[2] > 0 && head[2] < (64u << 20);
    if (ok) { bin.resize(head[2]); ok = fread(bin.data(), 1, bin.size(), f) == bin.size() && Checksum(bin) == head[3]; }
    fclose(f);
    if (!ok) { // short, damaged or an older layout: never given to the driver; deleted (a rename on Windows may not replace it)
        std::error_code ec;                                    // and compiled, then saved again
        std::filesystem::remove(file, ec);
        return 0;
    }
    const unsigned prog = glad_glCreateProgram();
    glad_glProgramBinary(prog, head[1], bin.data(), (int)bin.size());
    int linked = 0;
    glad_glGetProgramiv(prog, kLinkStatus, &linked);
    if (linked) return prog;
    if (glad_glDeleteProgram) glad_glDeleteProgram(prog);
    std::error_code ec;
    std::filesystem::remove(file, ec); // refused (another driver): compiled and written again
    return 0;
}
void Save(const std::string& file, unsigned prog) {
    int len = 0;
    glad_glGetProgramiv(prog, kBinaryLength, &len);
    if (len <= 0) return;
    std::vector<char> bin((size_t)len);
    unsigned format = 0;
    int got = 0;
    glad_glGetProgramBinary(prog, len, &got, &format, bin.data());
    if (got <= 0) return;
    std::error_code ec;
    std::filesystem::create_directories(Dir(), ec);
    bin.resize((size_t)got);
    // Written aside under a name of its own, then renamed over the file: scenario processes run side by side (#1768) and may
    // fill a cold cache with the same shader at once; a reader sees the old file, the new one or none, never half of one.
    // (Where a rename can't replace an existing file, the one already there stays: the same program.)
    const std::string tmp = file + TextFormat(".%08x%08x.tmp", (unsigned)std::random_device{}(), (unsigned)time(nullptr));
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f) return;
    const uint32_t head[4] = { kMagic, format, (uint32_t)got, Checksum(bin) };
    const bool ok = fwrite(head, sizeof head, 1, f) == 1 && fwrite(bin.data(), 1, bin.size(), f) == bin.size();
    if (fclose(f) != 0 || !ok) { std::filesystem::remove(tmp, ec); return; }
    std::filesystem::rename(tmp, file, ec);
    if (ec) std::filesystem::remove(tmp, ec);
}
} // namespace BinCache
}

bool ProgramCache_Usable() { return BinCache::Usable(); }
std::string ProgramCache_Path(const std::string& name, const std::string& src) { return BinCache::PathFor(name, src); }
unsigned ProgramCache_Load(const std::string& file) { return BinCache::Load(file); }
void ProgramCache_Save(const std::string& file, unsigned program) { BinCache::Save(file, program); }

int Compute_Prewarm(const std::vector<std::string>& paths) {
    if (!Compute_Available() || !ShaderAsync_Available()) return 0;
    int n = 0;
    for (const std::string& p : paths) {
        if (gPending.count(p)) continue;
        const std::string src = ComputeSource(p);
        if (src.empty()) continue;
        if (BinCache::Usable() && std::filesystem::exists(BinCache::PathFor(p, src))) continue; // (#1769) load() reads its binary
        AsyncProgram a = ShaderAsync_BeginCompute(src);
        if (a.started()) { gPending[p] = { a, GetTime() }; ++n; }
    }
    if (n) LOG_INFO("Compute", "%d compute shaders compiling in the background\n", n);
    return n;
}

int Compute_PrewarmPending() {
    int n = 0;
    for (const auto& kv : gPending) n += ShaderAsync_Done(kv.second.first) ? 0 : 1;
    return n;
}

bool ComputeProgram::load(const std::string& p) {
    if (!Compute_Available()) return false;
    const bool cache = BinCache::Usable();
    std::string src;
    if (cache && !gPending.count(p)) { // (#1769) a binary from an earlier run
        src = ComputeSource(p);
        if (src.empty()) return false;
        const double t0 = GetTime();
        if (const unsigned prog = BinCache::Load(BinCache::PathFor(p, src))) {
            if (const double ms = (GetTime() - t0) * 1000.0; ms > 50.0) TraceLog(LOG_WARNING, "Compute shader %s: cached binary loaded in %.0f ms", p.c_str(), ms);
            unload();
            id = prog;
            path = p;
            return true;
        }
    }
    if (auto it = gPending.find(p); it != gPending.end()) { // (#1047) started early: usually done by now
        const double t0 = GetTime();
        const bool wasDone = ShaderAsync_Done(it->second.first);
        const unsigned prog = ShaderAsync_Finish(it->second.first, p.c_str());
        const double waitedMs = (GetTime() - t0) * 1000.0, sinceMs = (t0 - it->second.second) * 1000.0;
        gPending.erase(it);
        if (waitedMs > 50.0 || !wasDone)
            TraceLog(LOG_WARNING, "Compute shader %s: background compile %s; waited %.0f ms for it (started %.0f ms before)", p.c_str(),
                     wasDone ? "done" : "not done yet", waitedMs, sinceMs);
        if (prog == 0) { TraceLog(LOG_WARNING, "Compute shader %s failed to compile", p.c_str()); return false; }
        if (cache) { if (src.empty()) src = ComputeSource(p); if (!src.empty()) BinCache::Save(BinCache::PathFor(p, src), prog); }
        unload();
        id = prog;
        path = p;
        return true;
    }
    if (src.empty()) src = ComputeSource(p);
    if (src.empty()) return false;
    const double t0 = GetTime();
    const unsigned int sh = rlCompileShader(src.c_str(), RL_COMPUTE_SHADER);
    if (sh == 0) { LOG_ERROR("Compute", "compute shader %s failed to compile\n", p.c_str()); TraceLog(LOG_WARNING, "Compute shader %s failed to compile", p.c_str()); return false; }
    const unsigned int prog = rlLoadComputeShaderProgram(sh);
    if (prog == 0) { LOG_ERROR("Compute", "compute shader %s failed to link\n", p.c_str()); return false; }
    // (#1047) a slow compile stalls its frame: a driver with a cold shader cache took seconds a shader
    if (const double ms = (GetTime() - t0) * 1000.0; ms > 50.0) TraceLog(LOG_WARNING, "Compute shader %s: compiled and linked in %.0f ms", p.c_str(), ms);
    if (cache) BinCache::Save(BinCache::PathFor(p, src), prog);
    unload();
    id = prog;
    path = p;
    return true;
}
void ComputeProgram::unload() { if (id) rlUnloadShaderProgram(id); id = 0; }
void ComputeProgram::use() const { rlEnableShader(id); }
int ComputeProgram::uniform(const char* name) const { return rlGetLocationUniform(id, name); }
void ComputeProgram::setInt(const char* name, int v) const { const int l = uniform(name); if (l >= 0) rlSetUniform(l, &v, RL_SHADER_UNIFORM_INT, 1); }
void ComputeProgram::setUint(const char* name, unsigned v) const { const int l = uniform(name); if (l >= 0) rlSetUniform(l, &v, RL_SHADER_UNIFORM_UINT, 1); }
void ComputeProgram::setFloat(const char* name, float v) const { const int l = uniform(name); if (l >= 0) rlSetUniform(l, &v, RL_SHADER_UNIFORM_FLOAT, 1); }
void ComputeProgram::setVec3(const char* name, float x, float y, float z) const { const float v[3] = { x, y, z }; const int l = uniform(name); if (l >= 0) rlSetUniform(l, v, RL_SHADER_UNIFORM_VEC3, 1); }
void ComputeProgram::setVec4(const char* name, float x, float y, float z, float w) const { const float v[4] = { x, y, z, w }; const int l = uniform(name); if (l >= 0) rlSetUniform(l, v, RL_SHADER_UNIFORM_VEC4, 1); }
void ComputeProgram::dispatch(unsigned gx, unsigned gy, unsigned gz) const {
    rlComputeShaderDispatch(gx, gy, gz);
    Compute_Barrier();
}

void ComputeProgram::dispatchIndirect(const ComputeBuffer& args, unsigned offset) const {
    constexpr unsigned kDispatchIndirectBuffer = 0x90EE; // GL_DISPATCH_INDIRECT_BUFFER
    glad_glBindBuffer(kDispatchIndirectBuffer, args.id);
    glad_glDispatchComputeIndirect((long long)offset);
    glad_glBindBuffer(kDispatchIndirectBuffer, 0);
    Compute_Barrier();
}

void ComputeBuffer::create(unsigned int size, const void* data, bool dynamic) {
    destroy();
    id = rlLoadShaderBuffer(size, data, dynamic ? RL_DYNAMIC_COPY : RL_STATIC_DRAW);
    bytes = size;
}
void ComputeBuffer::resize(unsigned int size) { create(size, nullptr, true); }
void ComputeBuffer::upload(const void* data, unsigned int size, unsigned int offset) { if (id) rlUpdateShaderBuffer(id, data, size, offset); }
void ComputeBuffer::download(void* dest, unsigned int size, unsigned int offset) const { if (id) rlReadShaderBuffer(id, dest, size, offset); }
void ComputeBuffer::bind(unsigned int index) const { rlBindShaderBuffer(id, index); }
void ComputeBuffer::destroy() { if (id) rlUnloadShaderBuffer(id); id = 0; bytes = 0; }

void Compute_Barrier() {
    // GL_SHADER_STORAGE_BARRIER_BIT 0x2000 | GL_BUFFER_UPDATE_BARRIER_BIT 0x200 | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT 0x1
    // | GL_COMMAND_BARRIER_BIT 0x40 | GL_ELEMENT_ARRAY_BARRIER_BIT 0x2
    if (glad_glMemoryBarrier) glad_glMemoryBarrier(0x2000 | 0x200 | 0x1 | 0x40 | 0x2);
}
