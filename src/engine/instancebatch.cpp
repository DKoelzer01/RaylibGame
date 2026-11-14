#include "instancebatch.h"
#include "scene.h"
#include "world.h"
#include "worldtypes.h"
#include <algorithm>

#include <raymath.h>
#include <rlgl.h>
#include <cmath>
#include <vector>

namespace {
    struct Inst { float row[12]; float color[4]; float params[4]; }; // model rows (3x4), colour, params
    std::vector<Inst> cubes, cyls;
    int lastCount = 0;

    struct Shape { unsigned int vao = 0, vbo = 0, ibo = 0, inst = 0, instCap = 0; int indexCount = 0; };
    Shape cubeShape, cylShape;
    HotShader shader{ "resources/instanced.vs", "resources/instanced.fs" };
    bool ready = false;
    float buried = 0.0f; // IB_SetBuried (#284)

    void Push(std::vector<Inst>& to, const Matrix& m, Color c, float p0, bool glow) {
        Inst in;
        // Matrix (raymath, column-major fields m0..m15): rows of the 3x4 transform.
        const float r[12] = { m.m0, m.m4, m.m8, m.m12, m.m1, m.m5, m.m9, m.m13, m.m2, m.m6, m.m10, m.m14 };
        for (int k = 0; k < 12; ++k) in.row[k] = r[k];
        in.color[0] = c.r / 255.0f; in.color[1] = c.g / 255.0f; in.color[2] = c.b / 255.0f; in.color[3] = c.a / 255.0f;
        in.params[0] = p0; in.params[1] = glow ? 1.0f : 0.0f; in.params[2] = buried; in.params[3] = 0.0f;
        to.push_back(in);
    }

    void BuildShape(Shape& s, const std::vector<float>& posNrm, const std::vector<unsigned short>& idx) {
        s.vao = rlLoadVertexArray();
        rlEnableVertexArray(s.vao);
        s.vbo = rlLoadVertexBuffer(posNrm.data(), (int)(posNrm.size() * sizeof(float)), false);
        rlSetVertexAttribute(0, 3, RL_FLOAT, false, 6 * sizeof(float), 0);
        rlEnableVertexAttribute(0);
        rlSetVertexAttribute(2, 3, RL_FLOAT, false, 6 * sizeof(float), 3 * sizeof(float));
        rlEnableVertexAttribute(2);
        s.ibo = rlLoadVertexBufferElement(idx.data(), (int)(idx.size() * sizeof(unsigned short)), false);
        s.indexCount = (int)idx.size();
        s.instCap = 256;
        s.inst = rlLoadVertexBuffer(nullptr, (int)(s.instCap * sizeof(Inst)), true);
        for (int a = 0; a < 5; ++a) {
            rlSetVertexAttribute(10 + a, 4, RL_FLOAT, false, sizeof(Inst), a * 4 * (int)sizeof(float));
            rlEnableVertexAttribute(10 + a);
            rlSetVertexAttributeDivisor(10 + a, 1);
        }
        rlDisableVertexArray();
    }

    void Init() {
        if (ready) return;
        ready = true;
        shader.reloadIfChanged(true);
        // Cube [-0.5, 0.5]^3, flat normals.
        std::vector<float> v; std::vector<unsigned short> ix;
        const float n[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
        for (int f = 0; f < 6; ++f) {
            const Vector3 N = { n[f][0], n[f][1], n[f][2] };
            const Vector3 U = fabsf(N.y) > 0.5f ? Vector3{ 1, 0, 0 } : Vector3{ 0, 1, 0 };
            const Vector3 R = Vector3CrossProduct(U, N);
            const unsigned short base = (unsigned short)(v.size() / 6);
            for (int k = 0; k < 4; ++k) {
                const float a = (k == 1 || k == 2) ? 0.5f : -0.5f, b = k >= 2 ? 0.5f : -0.5f;
                const Vector3 p = Vector3Add(Vector3Scale(N, 0.5f), Vector3Add(Vector3Scale(R, a), Vector3Scale(U, b)));
                v.insert(v.end(), { p.x, p.y, p.z, N.x, N.y, N.z });
            }
            ix.insert(ix.end(), { base, (unsigned short)(base + 1), (unsigned short)(base + 2), base, (unsigned short)(base + 2), (unsigned short)(base + 3) });
        }
        BuildShape(cubeShape, v, ix);
        // Cylinder: radius 1, y 0..1, 12 sides with smooth normals, capped.
        v.clear(); ix.clear();
        constexpr int S = 12;
        for (int k = 0; k <= S; ++k) {
            const float a = 2.0f * PI * k / S, c = cosf(a), s = sinf(a);
            v.insert(v.end(), { c, 0.0f, s, c, 0.0f, s });
            v.insert(v.end(), { c, 1.0f, s, c, 0.0f, s });
        }
        for (int k = 0; k < S; ++k) {
            const unsigned short a = (unsigned short)(2 * k), b = (unsigned short)(2 * k + 2);
            ix.insert(ix.end(), { a, (unsigned short)(a + 1), (unsigned short)(b + 1), a, (unsigned short)(b + 1), b });
        }
        for (int cap = 0; cap < 2; ++cap) {
            const unsigned short centre = (unsigned short)(v.size() / 6);
            const float y = (float)cap, ny = cap ? 1.0f : -1.0f;
            v.insert(v.end(), { 0.0f, y, 0.0f, 0.0f, ny, 0.0f });
            for (int k = 0; k <= S; ++k) { const float a = 2.0f * PI * k / S; v.insert(v.end(), { cosf(a), y, sinf(a), 0.0f, ny, 0.0f }); }
            for (int k = 0; k < S; ++k) ix.insert(ix.end(), { centre, (unsigned short)(centre + 1 + k), (unsigned short)(centre + 2 + k) });
        }
        BuildShape(cylShape, v, ix);
    }

    void Draw(Shape& s, std::vector<Inst>& list, int shapeId, const Matrix& viewProj) {
        if (list.empty()) return;
        if ((unsigned)list.size() > s.instCap) {
            while (s.instCap < (unsigned)list.size()) s.instCap *= 2;
            rlEnableVertexArray(s.vao);
            rlUnloadVertexBuffer(s.inst);
            s.inst = rlLoadVertexBuffer(nullptr, (int)(s.instCap * sizeof(Inst)), true);
            for (int a = 0; a < 5; ++a) {
                rlSetVertexAttribute(10 + a, 4, RL_FLOAT, false, sizeof(Inst), a * 4 * (int)sizeof(float));
                rlEnableVertexAttribute(10 + a);
                rlSetVertexAttributeDivisor(10 + a, 1);
            }
            rlDisableVertexArray();
        }
        rlUpdateVertexBuffer(s.inst, list.data(), (int)(list.size() * sizeof(Inst)), 0);
        rlSetUniform(rlGetLocationUniform(shader.shader.id, "shape"), &shapeId, RL_SHADER_UNIFORM_INT, 1);
        rlSetUniformMatrix(rlGetLocationUniform(shader.shader.id, "matViewProj"), viewProj);
        rlEnableVertexArray(s.vao);
        rlDrawVertexArrayElementsInstanced(0, s.indexCount, 0, (int)list.size());
        rlDisableVertexArray();
        lastCount += (int)list.size();
        list.clear();
    }
}

void IB_Cube(Vector3 centre, Vector3 size, Color color, bool glow) {
    const Matrix m = MatrixMultiply(MatrixMultiply(MatrixScale(size.x, size.y, size.z), MatrixTranslate(centre.x, centre.y, centre.z)), rlGetMatrixTransform());
    Push(cubes, m, color, 1.0f, glow);
}

void IB_Cylinder(Vector3 base, float rTop, float rBottom, float height, Color color, bool glow) {
    const float rb = std::max(rBottom, 1e-4f);
    const Matrix m = MatrixMultiply(MatrixMultiply(MatrixScale(rb, height, rb), MatrixTranslate(base.x, base.y, base.z)), rlGetMatrixTransform());
    Push(cyls, m, color, rTop / rb, glow);
}

void IB_SetBuried(float b) { buried = std::clamp(b, 0.0f, 1.0f); }

void InstanceBatch_Flush(Scene& world) {
    lastCount = 0;
    if (cubes.empty() && cyls.empty()) return;
    Init();
    shader.reloadIfChanged();
    if (shader.shader.id == 0) { cubes.clear(); cyls.clear(); return; }
    rlDrawRenderBatchActive();
    const Matrix viewProj = MatrixMultiply(rlGetMatrixModelview(), rlGetMatrixProjection());
    rlEnableShader(shader.shader.id);
    const unsigned int id = shader.shader.id;
    auto v3 = [&](const char* n, Vector3 v) { rlSetUniform(rlGetLocationUniform(id, n), &v, RL_SHADER_UNIFORM_VEC3, 1); };
    v3("viewPos", world.camera.position);
    v3("sunDir", world.frameSun);
    v3("ambient", world.frameAmbient);
    v3("fogColor", world.frameFog);
    rlSetUniform(rlGetLocationUniform(id, "fogDensity"), &world.frameFogDensity, RL_SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(rlGetLocationUniform(id, "fogShell"), &world.frameFogShell, RL_SHADER_UNIFORM_VEC4, 1);
    {   // The night side (#284): the nearest planet or moon blocks the sun, as for parts (Planetoid::draw).
        float occ[4] = { 0, 0, 0, 0 };
        if (const Planetoid* p = NearestPlanetoid(world, world.camera.position); p && (p->kind == BodyKind::Planet || p->kind == BodyKind::Moon)) {
            const Vector3 c = (p->worldPos - world.renderOrigin).toFloat();
            occ[0] = c.x; occ[1] = c.y; occ[2] = c.z; occ[3] = std::max(0.0f, p->radius - 60.0f) * p->scale;
        }
        rlSetUniform(rlGetLocationUniform(id, "sunOccluder"), occ, RL_SHADER_UNIFORM_VEC4, 1);
    }
    if (world.frameLighting) { world.frameLighting(shader.shader); rlEnableShader(id); }
    rlDisableBackfaceCulling(); // instances may be mirrored (the plot frame is left-handed)
    Draw(cubeShape, cubes, 0, viewProj);
    Draw(cylShape, cyls, 1, viewProj);
    rlEnableBackfaceCulling();
    rlDisableShader();
}

int InstanceBatch_LastCount() { return lastCount; }

void InstanceBatch_Unload() {
    for (Shape* s : { &cubeShape, &cylShape }) {
        if (s->vao) rlUnloadVertexArray(s->vao);
        if (s->vbo) rlUnloadVertexBuffer(s->vbo);
        if (s->ibo) rlUnloadVertexBuffer(s->ibo);
        if (s->inst) rlUnloadVertexBuffer(s->inst);
        *s = Shape{};
    }
    shader.unload();
    ready = false;
}
