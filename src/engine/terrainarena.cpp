#include "terrainarena.h"
#include "voxel.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <rlgl.h>
#include <raymath.h>
#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <map>

namespace {
    constexpr int kRecordLocation = 9;     // per-draw attributes 9-13, after raylib's default attributes (0-8)
    constexpr int kRecordTexels = 5;       // model rows 0-2, texturing origin, clip

    // Positions in their own stream (the shadow passes read only them); everything else interleaved.
    struct Attr { float nrm[3], uv[2], uv2[2], tan[4]; unsigned char col[4]; };
    static_assert(sizeof(Attr) == 48, "interleaved terrain attributes");
    constexpr size_t kPosBytes = 12;

    struct Command { GLuint count, instanceCount, firstIndex; GLint baseVertex; GLuint baseInstance; };

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
        GLuint vao = 0, vaoDepth = 0, pos = 0, attr = 0, ibo = 0, records = 0, indirect = 0;
        Ranges verts, idx;
        std::vector<float> rec;
        std::vector<Command> cmds;
        Matrix firstModel{};  // the first queued chunk: its rotation (shared by the body's chunks) is the shader's matModel
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
    GLuint Regrow(GLuint old, size_t oldBytes, size_t newBytes) {
        GLuint b = 0;
        glGenBuffers(1, &b);
        glBindBuffer(GL_COPY_WRITE_BUFFER, b);
        glBufferData(GL_COPY_WRITE_BUFFER, (GLsizeiptr)newBytes, nullptr, GL_STATIC_DRAW);
        if (old && oldBytes) {
            glBindBuffer(GL_COPY_READ_BUFFER, old);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, (GLsizeiptr)oldBytes);
            glBindBuffer(GL_COPY_READ_BUFFER, 0);
        }
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
        if (old) glDeleteBuffers(1, &old);
        return b;
    }

    bool Init() {
        if (S.tried) return S.ok;
        S.tried = true;
        const char* env = getenv("RAYLIBGAME_ARENA");
        if (env && env[0] == '0') { logger.logf("[Arena] off (RAYLIBGAME_ARENA=0)\n"); return false; }
        if (rlGetVersion() != RL_OPENGL_43 || !glad_glMultiDrawElementsIndirect || !glad_glCopyBufferSubData) {
            logger.logf("[Arena] multi-draw indirect unavailable: per-chunk draws\n");
            return false;
        }
        const int v0 = 4 << 20, i0 = 24 << 20; // about a planet walk's working set (5.7 M vertices seen); grows as needed
        S.pos = Regrow(0, 0, (size_t)v0 * kPosBytes);
        S.attr = Regrow(0, 0, (size_t)v0 * sizeof(Attr));
        S.ibo = Regrow(0, 0, (size_t)i0 * sizeof(uint16_t));
        S.verts.grow(v0);
        S.idx.grow(i0);
        glGenVertexArrays(1, &S.vao);
        glGenVertexArrays(1, &S.vaoDepth);
        glGenBuffers(1, &S.records);
        glGenBuffers(1, &S.indirect);
        BindAttributes();
        S.ok = true;
        logger.logf("[Arena] terrain mesh arena: %d vertices, %d indices\n", v0, i0);
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
    glBindBuffer(GL_ARRAY_BUFFER, S.pos);
    glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)a.vFirst * kPosBytes, (GLsizeiptr)nv * kPosBytes, d.positions.data());
    glBindBuffer(GL_ARRAY_BUFFER, S.attr);
    glBufferSubData(GL_ARRAY_BUFFER, (GLintptr)a.vFirst * sizeof(Attr), (GLsizeiptr)nv * sizeof(Attr), v.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    // The element buffer binding is VAO state: upload through COPY_WRITE so no VAO is disturbed.
    glBindBuffer(GL_COPY_WRITE_BUFFER, S.ibo);
    glBufferSubData(GL_COPY_WRITE_BUFFER, (GLintptr)a.iFirst * sizeof(uint16_t), (GLsizeiptr)ni * sizeof(uint16_t), d.indices.data());
    glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
    return a;
}

void TerrainArena_Free(ArenaAlloc& a) {
    if (a.vFirst >= 0 && S.ok) S.verts.release(a.vFirst, a.vCount);
    if (a.iFirst >= 0 && S.ok) S.idx.release(a.iFirst, a.iCount);
    a = ArenaAlloc{};
}

void TerrainArena_Add(const ArenaAlloc& a, const Matrix& m, Vector3 texOrigin, Vector3 clip) {
    if (!a.valid()) return;
    const float r[kRecordTexels * 4] = {
        m.m0, m.m4, m.m8, m.m12,  m.m1, m.m5, m.m9, m.m13,  m.m2, m.m6, m.m10, m.m14,
        texOrigin.x, texOrigin.y, texOrigin.z, 0.0f,  clip.x, clip.y, clip.z, 0.0f };
    const GLuint index = (GLuint)S.cmds.size();
    if (index == 0) S.firstModel = m;
    S.rec.insert(S.rec.end(), r, r + kRecordTexels * 4);
    S.cmds.push_back({ (GLuint)a.iCount, 1u, (GLuint)a.iFirst, (GLint)a.vFirst, index });
}

void TerrainArena_Flush(unsigned int shaderId, bool positionsOnly) {
    if (S.cmds.empty()) return;
    const int n = (int)S.cmds.size();
    glBindBuffer(GL_ARRAY_BUFFER, S.records); // new storage each flush (orphaned): the VAOs keep the buffer name
    glBufferData(GL_ARRAY_BUFFER, (GLsizeiptr)(S.rec.size() * sizeof(float)), S.rec.data(), GL_STREAM_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, S.indirect);
    glBufferData(GL_DRAW_INDIRECT_BUFFER, (GLsizeiptr)(n * sizeof(Command)), S.cmds.data(), GL_STREAM_DRAW);

    const Matrix viewProj = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
    const GLint locArena = glGetUniformLocation(shaderId, "arenaDraw");
    const GLint locVP = glGetUniformLocation(shaderId, "matViewProj");
    const GLint locModel = glGetUniformLocation(shaderId, "matModel");
    glUseProgram(shaderId);
    if (locArena >= 0) glUniform1i(locArena, 1);
    if (locVP >= 0) rlSetUniformMatrix(locVP, viewProj);
    if (locModel >= 0) rlSetUniformMatrix(locModel, S.firstModel); // fragment normals: body rotation
    glBindVertexArray(positionsOnly ? S.vaoDepth : S.vao);
    glMultiDrawElementsIndirect(GL_TRIANGLES, GL_UNSIGNED_SHORT, nullptr, n, 0);
    glBindVertexArray(0);
    if (locArena >= 0) glUniform1i(locArena, 0); // other meshes drawn with this shader use their own uniforms
    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, 0);
    S.draws += n;
    S.flushes++;
    S.rec.clear();
    S.cmds.clear();
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
    GLuint bufs[] = { S.pos, S.attr, S.ibo, S.records, S.indirect };
    glDeleteBuffers(5, bufs);
    GLuint vaos[] = { S.vao, S.vaoDepth };
    glDeleteVertexArrays(2, vaos);
    S = State{};
}
