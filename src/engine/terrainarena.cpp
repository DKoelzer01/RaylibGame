#include "terrainarena.h"
#include "voxel.h"
#include "frustum.h"
#include "gpucompute.h"
#include "glcaps.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <rlgl.h>
#include <raymath.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <map>

extern "C" void* glfwGetProcAddress(const char* name); // raylib links GLFW

namespace {
    constexpr int kRecordLocation = 9;     // per-draw attributes 9-13, after raylib's default attributes (0-8)
    constexpr int kRecordTexels = 5;       // model rows 0-2, texturing origin, clip

    // Positions in their own stream (the shadow passes read only them); everything else interleaved.
    struct Attr { float nrm[3], uv[2], uv2[2], tan[4]; unsigned char col[4]; };
    static_assert(sizeof(Attr) == 48, "interleaved terrain attributes");
    constexpr size_t kPosBytes = 12;

    struct Command { GLuint count, instanceCount, firstIndex; GLint baseVertex; GLuint baseInstance; };
    static_assert(sizeof(Command) == 20, "indirect command");
    struct TableEntry { float origin[3], size; float clip[4]; GLuint iCount, iFirst, vFirst, pad; }; // chunk_cull.comp Chunk
    static_assert(sizeof(TableEntry) == 48, "chunk table entry");

    // First-fit range allocator over [0, cap), free ranges coalesced.
    struct Ranges {
        int cap = 0, used = 0;
        std::map<int, int> free; // offset -> size
        int alloc(int n) {
            for (auto it = free.begin(); it != free.end(); ++it) {
                if (it->second < n) continue;
                const int off = it->first, rest = it->second - n;
                free.erase(it);
                if (rest > 0) free[off + n] = rest;
                used += n;
                return off;
            }
            return -1;
        }
        void release(int off, int n) {
            used -= n;
            auto next = free.lower_bound(off);
            if (next != free.end() && off + n == next->first) { n += next->second; next = free.erase(next); }
            if (next != free.begin()) {
                auto prev = std::prev(next);
                if (prev->first + prev->second == off) { prev->second += n; return; }
            }
            free[off] = n;
        }
        void grow(int newCap) { const int u = used; release(cap, newCap - cap); used = u; cap = newCap; }
    };

    struct State {
        bool tried = false, ok = false;
        GLuint vao = 0, vaoDepth = 0, pos = 0, attr = 0, ibo = 0;
        GLuint table = 0, list = 0, records = 0, indirect = 0;
        // #286 indirect count (4.6 / ARB_indirect_parameters): the visible commands compacted in order, and their count.
        typedef void (__stdcall* PfnDrawCount)(GLenum, GLenum, const void*, GLintptr, GLsizei, GLsizei);
        PfnDrawCount drawCount = nullptr;
        GLuint compacted = 0, counts = 0;
        ComputeProgram compact;
        GLint cpCount = -1, cpBase = -1;
        Ranges verts, idx;
        int tableCap = 0;
        std::vector<int> freeSlots;
        int slotCount = 0;
        // The per-pass command / record / list ring: each flush takes the next n entries (wrapping), so a
        // compute pass never rewrites entries a draw issued just before it may still read.
        int ringCap = 0, ringPos = 0;
        std::vector<GLuint> queued;
        ComputeProgram cull;
        struct CullLocs { GLint count, base, rot, scale, offset, texPeriod, useFrustum, planes, useColumn, colPoint, colDir, colRadius, overlap; } cl{}; // looked up once
        struct DrawLocs { GLuint shader; GLint arena, viewProj, model; };
        std::vector<DrawLocs> drawLocs; // per draw shader (hot reload makes new ids)
        int draws = 0, flushes = 0;
    } S;

    void BindRecords() {
        // Per-draw records as instanced attributes: each command's baseInstance picks its record.
        glBindBuffer(GL_ARRAY_BUFFER, S.records);
        for (int k = 0; k < kRecordTexels; ++k) {
            glVertexAttribPointer(kRecordLocation + k, 4, GL_FLOAT, GL_FALSE, kRecordTexels * 16, (void*)(size_t)(k * 16));
            glEnableVertexAttribArray(kRecordLocation + k);
            glVertexAttribDivisor(kRecordLocation + k, 1);
        }
    }

    void BindAttributes() {
        // Full VAO (lit pass).
        glBindVertexArray(S.vao);
        glBindBuffer(GL_ARRAY_BUFFER, S.pos);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, (GLsizei)kPosBytes, nullptr);
        glBindBuffer(GL_ARRAY_BUFFER, S.attr);
        const GLsizei st = sizeof(Attr);
        glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, st, (void*)offsetof(Attr, uv));
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, st, (void*)offsetof(Attr, nrm));
        glVertexAttribPointer(3, 4, GL_UNSIGNED_BYTE, GL_TRUE, st, (void*)offsetof(Attr, col));
        glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, st, (void*)offsetof(Attr, tan));
        glVertexAttribPointer(5, 2, GL_FLOAT, GL_FALSE, st, (void*)offsetof(Attr, uv2));
        for (int a = 0; a <= 5; ++a) glEnableVertexAttribArray(a);
        BindRecords();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, S.ibo);
        // Depth VAO: positions and records only.
        glBindVertexArray(S.vaoDepth);
        glBindBuffer(GL_ARRAY_BUFFER, S.pos);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, (GLsizei)kPosBytes, nullptr);
        glEnableVertexAttribArray(0);
        BindRecords();
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, S.ibo);
        glBindVertexArray(0);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    // A bigger buffer with the old contents copied over.
    GLuint Regrow(GLuint old, size_t oldBytes, size_t newBytes, GLenum usage = GL_STATIC_DRAW) {
        GLuint b = 0;
        glGenBuffers(1, &b);
        glBindBuffer(GL_COPY_WRITE_BUFFER, b);
        glBufferData(GL_COPY_WRITE_BUFFER, (GLsizeiptr)newBytes, nullptr, usage);
        if (old && oldBytes) {
            glBindBuffer(GL_COPY_READ_BUFFER, old);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, (GLsizeiptr)oldBytes);
            glBindBuffer(GL_COPY_READ_BUFFER, 0);
        }
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
        if (old) glDeleteBuffers(1, &old);
        return b;
    }

    void EnsureRing(int n) {
        if (n <= S.ringCap) return;
        int cap = std::max(S.ringCap, 8192);
        while (cap < n * 4) cap *= 2; // room for a frame's passes before wrapping
        // Contents need not survive (each flush writes what it draws).
        GLuint old[] = { S.list, S.indirect, S.records, S.compacted, S.counts };
        glDeleteBuffers(5, old);
        S.compacted = S.counts = 0;
        S.list = Regrow(0, 0, (size_t)cap * sizeof(GLuint), GL_DYNAMIC_DRAW);
        S.indirect = Regrow(0, 0, (size_t)cap * sizeof(Command), GL_DYNAMIC_COPY);
        S.records = Regrow(0, 0, (size_t)cap * kRecordTexels * 16, GL_DYNAMIC_COPY);
        if (S.drawCount) {
            S.compacted = Regrow(0, 0, (size_t)cap * sizeof(Command), GL_DYNAMIC_COPY);
            S.counts = Regrow(0, 0, (size_t)cap * sizeof(GLuint), GL_DYNAMIC_COPY); // one per flush, at its ring base
        }
        S.ringCap = cap;
        S.ringPos = 0;
        BindAttributes();
    }

    bool Init() {
        if (S.tried) return S.ok;
        S.tried = true;
        const char* env = getenv("RAYLIBGAME_ARENA");
        if (env && env[0] == '0') { logger.logf("[Arena] off (RAYLIBGAME_ARENA=0)\n"); return false; }
        if (rlGetVersion() != RL_OPENGL_43 || !glad_glMultiDrawElementsIndirect || !glad_glCopyBufferSubData || !Compute_Available()) {
            logger.logf("[Arena] multi-draw indirect / compute unavailable: per-chunk draws\n");
            return false;
        }
        if (!S.cull.load("resources/compute/chunk_cull.comp")) { logger.logf("[Arena] chunk_cull.comp failed: per-chunk draws\n"); return false; }
        {
            const GLuint p = S.cull.id;
            auto L = [p](const char* n) { return glGetUniformLocation(p, n); };
            S.cl = { L("count"), L("base"), L("rot"), L("scale"), L("offset"), L("texPeriod"), L("useFrustum"), L("planes"),
                     L("useColumn"), L("colPoint"), L("colDir"), L("colRadius"), L("overlap") };
        }
        {   // #286: RAYLIBGAME_INDIRECT_COUNT=1 draws only the visible commands (compacted in order, indirect count draws).
            const char* ic = getenv("RAYLIBGAME_INDIRECT_COUNT");
            if (GL_Caps().indirectParameters && ic && ic[0] == '1') { // opt-in: measured no faster on an RTX 5070 (empty commands cost ~nothing there)
                S.drawCount = (State::PfnDrawCount)glfwGetProcAddress("glMultiDrawElementsIndirectCount");
                if (!S.drawCount) S.drawCount = (State::PfnDrawCount)glfwGetProcAddress("glMultiDrawElementsIndirectCountARB");
                if (S.drawCount && !S.compact.load("resources/compute/chunk_compact.comp")) S.drawCount = nullptr;
                if (S.drawCount) { S.cpCount = glGetUniformLocation(S.compact.id, "count"); S.cpBase = glGetUniformLocation(S.compact.id, "base"); }
            }
            logger.logf("[Arena] draws: %s\n", S.drawCount ? "indirect count (visible chunks compacted in order)" : "all listed commands (culled ones empty)");
        }
        const int v0 = 4 << 20, i0 = 24 << 20; // about a planet walk's working set (5.7 M vertices seen); grows as needed
        S.pos = Regrow(0, 0, (size_t)v0 * kPosBytes);
        S.attr = Regrow(0, 0, (size_t)v0 * sizeof(Attr));
        S.ibo = Regrow(0, 0, (size_t)i0 * sizeof(uint16_t));
        S.verts.grow(v0);
        S.idx.grow(i0);
        S.tableCap = 4096;
        S.table = Regrow(0, 0, (size_t)S.tableCap * sizeof(TableEntry));
        glGenVertexArrays(1, &S.vao);
        glGenVertexArrays(1, &S.vaoDepth);
        EnsureRing(1); // creates the ring buffers and binds the VAOs
        S.ok = true;
        logger.logf("[Arena] terrain mesh arena: %d vertices, %d indices, GPU chunk culling\n", v0, i0);
        return true;
    }

    // Grows by half (not double: a doubled multi-hundred-MB buffer would be mostly empty) until n fits.
    int AllocVerts(int n) {
        int off = S.verts.alloc(n);
        while (off < 0) {
            const int oldCap = S.verts.cap, newCap = std::max(oldCap + oldCap / 2, oldCap + n);
            S.pos = Regrow(S.pos, (size_t)oldCap * kPosBytes, (size_t)newCap * kPosBytes);
            S.attr = Regrow(S.attr, (size_t)oldCap * sizeof(Attr), (size_t)newCap * sizeof(Attr));
            S.verts.grow(newCap);
            BindAttributes();
            logger.logf("[Arena] vertices grown to %d\n", newCap);
            off = S.verts.alloc(n);
        }
        return off;
    }
    int AllocIndices(int n) {
        int off = S.idx.alloc(n);
        while (off < 0) {
            const int oldCap = S.idx.cap, newCap = std::max(oldCap + oldCap / 2, oldCap + n);
            S.ibo = Regrow(S.ibo, (size_t)oldCap * sizeof(uint16_t), (size_t)newCap * sizeof(uint16_t));
            S.idx.grow(newCap);
            BindAttributes();
            logger.logf("[Arena] indices grown to %d\n", newCap);
            off = S.idx.alloc(n);
        }
        return off;
    }
    int AllocSlot() {
        if (!S.freeSlots.empty()) { const int s = S.freeSlots.back(); S.freeSlots.pop_back(); return s; }
        if (S.slotCount >= S.tableCap) {
            const int newCap = S.tableCap * 2;
            S.table = Regrow(S.table, (size_t)S.tableCap * sizeof(TableEntry), (size_t)newCap * sizeof(TableEntry));
            S.tableCap = newCap;
        }
        return S.slotCount++;
    }
    void WriteTable(int slot, const void* data, size_t bytes, size_t offset = 0) {
        glBindBuffer(GL_COPY_WRITE_BUFFER, S.table);
        glBufferSubData(GL_COPY_WRITE_BUFFER, (GLintptr)(slot * sizeof(TableEntry) + offset), (GLsizeiptr)bytes, data);
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    }
}

bool TerrainArena_Available() { return Init(); }

ArenaAlloc TerrainArena_Upload(const ChunkMeshData& d) {
    ArenaAlloc a;
    if (!Init() || d.positions.empty() || d.indices.empty()) return a;
    const int nv = (int)(d.positions.size() / 3), ni = (int)d.indices.size();
    std::vector<Attr> v(nv);
    const bool uv = d.texcoords.size() >= (size_t)nv * 2, uv2 = d.texcoords2.size() >= (size_t)nv * 2, tan = d.tangents.size() >= (size_t)nv * 4;
    const bool col = d.colors.size() >= (size_t)nv * 4;
    for (int i = 0; i < nv; ++i) {
        Attr& o = v[i];
        memcpy(o.nrm, &d.normals[i * 3], 12);
        if (uv) memcpy(o.uv, &d.texcoords[i * 2], 8); else o.uv[0] = o.uv[1] = 0.0f;
        if (uv2) memcpy(o.uv2, &d.texcoords2[i * 2], 8); else o.uv2[0] = o.uv2[1] = 0.0f;
        if (tan) memcpy(o.tan, &d.tangents[i * 4], 16); else { o.tan[0] = o.tan[1] = o.tan[2] = 0.0f; o.tan[3] = 1.0f; }
        if (col) memcpy(o.col, &d.colors[i * 4], 4); else memset(o.col, 255, 4);
    }
    a.vFirst = AllocVerts(nv);
    a.vCount = nv;
    a.iFirst = AllocIndices(ni);
    a.iCount = ni;
    a.slot = AllocSlot();
    glBindBuffer(GL_ARRAY_BUFFER, S.pos);
    glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)a.vFirst * kPosBytes, (GLsizeiptr)nv * kPosBytes, d.positions.data());
    glBindBuffer(GL_ARRAY_BUFFER, S.attr);
    glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)a.vFirst * sizeof(Attr), (GLsizeiptr)nv * sizeof(Attr), v.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    // The element buffer binding is VAO state: upload through COPY_WRITE so no VAO is disturbed.
    glBindBuffer(GL_COPY_WRITE_BUFFER, S.ibo);
    glBufferSubData(GL_COPY_WRITE_BUFFER, (GLintptr)a.iFirst * sizeof(uint16_t), (GLsizeiptr)ni * sizeof(uint16_t), d.indices.data());
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    const Vector3 o = d.key.origin();
    TableEntry e{};
    e.origin[0] = o.x; e.origin[1] = o.y; e.origin[2] = o.z; e.size = d.key.size();
    e.clip[0] = e.clip[1] = e.clip[2] = 1e30f;
    e.iCount = (GLuint)ni; e.iFirst = (GLuint)a.iFirst; e.vFirst = (GLuint)a.vFirst;
    WriteTable(a.slot, &e, sizeof(e));
    a.clip = { 1e30f, 1e30f, 1e30f };
    return a;
}

void TerrainArena_Free(ArenaAlloc& a) {
    if (S.ok) {
        if (a.vFirst >= 0) S.verts.release(a.vFirst, a.vCount);
        if (a.iFirst >= 0) S.idx.release(a.iFirst, a.iCount);
        if (a.slot >= 0) S.freeSlots.push_back(a.slot);
    }
    a = ArenaAlloc{};
}

void TerrainArena_Queue(ArenaAlloc& a, Vector3 clip) {
    if (!a.valid()) return;
    if (clip.x != a.clip.x || clip.y != a.clip.y || clip.z != a.clip.z) {
        const float c[4] = { clip.x, clip.y, clip.z, 0.0f };
        WriteTable(a.slot, c, sizeof(c), offsetof(TableEntry, clip));
        a.clip = clip;
    }
    S.queued.push_back((GLuint)a.slot);
}

int TerrainArena_Flush(unsigned int shaderId, bool positionsOnly, const ArenaBody& body, const ArenaCull& cull) {
    const int n = (int)S.queued.size();
    if (n == 0) return 0;
    EnsureRing(n);
    if (S.ringPos + n > S.ringCap) S.ringPos = 0;
    const int base = S.ringPos;
    S.ringPos += n;

    // Cull + build the commands and records on the GPU.
    glBindBuffer(GL_COPY_WRITE_BUFFER, S.list);
    glBufferSubData(GL_COPY_WRITE_BUFFER, (GLintptr)(base * sizeof(GLuint)), (GLsizeiptr)(n * sizeof(GLuint)), S.queued.data());
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    S.cull.use();
    const GLuint prog = S.cull.id;
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, S.table);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, S.list); // whole buffer: a range bind needs GL_SHADER_STORAGE_BUFFER_OFFSET_ALIGNMENT-aligned offsets (it failed, and the pass culled a stale list)
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, S.indirect);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, S.records);
    glUniform1ui(S.cl.count, (GLuint)n);
    glUniform1ui(S.cl.base, (GLuint)base);
    glUniformMatrix3dv(S.cl.rot, 1, GL_TRUE, body.rot);
    glUniform1d(S.cl.scale, body.scale);
    glUniform3d(S.cl.offset, body.offset[0], body.offset[1], body.offset[2]);
    glUniform1d(S.cl.texPeriod, body.texPeriod);
    extern bool frustumCullingEnabled;
    const bool useFrustum = cull.frustum && frustumCullingEnabled;
    glUniform1i(S.cl.useFrustum, useFrustum ? 1 : 0);
    if (useFrustum) glUniform4fv(S.cl.planes, 6, (const float*)cull.frustum->planes);
    glUniform1i(S.cl.useColumn, cull.colRadius > 0.0f ? 1 : 0);
    glUniform3f(S.cl.colPoint, cull.colPoint.x, cull.colPoint.y, cull.colPoint.z);
    glUniform3f(S.cl.colDir, cull.colDir.x, cull.colDir.y, cull.colDir.z);
    glUniform1f(S.cl.colRadius, cull.colRadius);
    glUniform1f(S.cl.overlap, (CHUNK_SIZE + 1.0f) / CHUNK_SIZE);
    glDispatchCompute((GLuint)((n + 63) / 64), 1, 1);
    if (S.drawCount) {
        // #286: the visible commands compacted in order, and their count (chunk_compact.comp).
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        S.compact.use();
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, S.indirect);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, S.compacted);
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, S.counts);
        glUniform1ui(S.cpCount, (GLuint)n);
        glUniform1ui(S.cpBase, (GLuint)base);
        glDispatchCompute(1, 1, 1);
    }
    glMemoryBarrier(GL_COMMAND_BARRIER_BIT | GL_VERTEX_ATTRIB_ARRAY_BARRIER_BIT);

    // Draw.
    const Matrix viewProj = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
    const State::DrawLocs* dl = nullptr;
    for (const auto& d : S.drawLocs) if (d.shader == shaderId) dl = &d;
    if (!dl) {
        S.drawLocs.push_back({ shaderId, glGetUniformLocation(shaderId, "arenaDraw"), glGetUniformLocation(shaderId, "matViewProj"), glGetUniformLocation(shaderId, "matModel") });
        dl = &S.drawLocs.back();
    }
    const GLint locArena = dl->arena, locVP = dl->viewProj, locModel = dl->model;
    glUseProgram(shaderId);
    if (locArena >= 0) glUniform1i(locArena, 1);
    if (locVP >= 0) rlSetUniformMatrix(locVP, viewProj);
    if (locModel >= 0) rlSetUniformMatrix(locModel, body.rs); // fragment normals: body rotation
    glBindVertexArray(positionsOnly ? S.vaoDepth : S.vao);
    if (S.drawCount) {
        constexpr GLenum kParameterBuffer = 0x80EE; // GL_PARAMETER_BUFFER
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, S.compacted);
        glBindBuffer(kParameterBuffer, S.counts);
        S.drawCount(GL_TRIANGLES, GL_UNSIGNED_SHORT, (const void*)(size_t)(base * sizeof(Command)), (GLintptr)(base * sizeof(GLuint)), n, 0);
        glBindBuffer(kParameterBuffer, 0);
    } else {
        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, S.indirect);
        glMultiDrawElementsIndirect(GL_TRIANGLES, GL_UNSIGNED_SHORT, (const void*)(size_t)(base * sizeof(Command)), n, 0);
    }
    glBindVertexArray(0);
    if (locArena >= 0) glUniform1i(locArena, 0); // other meshes drawn with this shader use their own uniforms
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
    S.draws += n;
    S.flushes++;
    S.queued.clear();
    return n;
}

ArenaStats TerrainArena_Stats() {
    ArenaStats s;
    s.vertices = S.verts.used; s.vertexCapacity = S.verts.cap;
    s.indices = S.idx.used; s.indexCapacity = S.idx.cap;
    s.draws = S.draws; s.flushes = S.flushes;
    return s;
}

void TerrainArena_Unload() {
    if (!S.ok) { S = State{}; return; }
    GLuint bufs[] = { S.pos, S.attr, S.ibo, S.table, S.list, S.records, S.indirect, S.compacted, S.counts };
    glDeleteBuffers(9, bufs);
    GLuint vaos[] = { S.vao, S.vaoDepth };
    glDeleteVertexArrays(2, vaos);
    S.cull.unload();
    S.compact.unload();
    S = State{};
}
