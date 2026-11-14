#include "scenetarget.h"
#include "glcaps.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <rlgl.h>
#include <cstdlib>
#include <cstring>

extern "C" void* glfwGetProcAddress(const char* name); // raylib links GLFW

namespace {
    typedef void (__stdcall* PfnClipControl)(GLenum origin, GLenum depth);
    constexpr GLenum kLowerLeft = 0x8CA1, kNegOneToOne = 0x935E, kZeroToOne = 0x935F;

    struct State {
        bool tried = false, ok = false;
        PfnClipControl clipControl = nullptr;
        GLuint fbo = 0, color = 0, depth = 0;     // the scene (multisampled when the window is)
        GLuint resolveFbo = 0, resolveTex = 0;    // single-sample copy, drawn to the window by a fullscreen pass
        GLuint prog = 0, vao = 0;
        int w = 0, h = 0, samples = 0;
    } S;

    // A fullscreen triangle sampling the resolved scene. Drawing (not blitting) to the window works whatever
    // its format and sample count: a multisampled blit needs identical formats and sample counts on both
    // sides, which only some drivers' default framebuffers happen to match (PR #156 review).
    const char* kVs = "#version 330\nout vec2 uv;\nvoid main() { vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2); uv = p; gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0); }\n";
    const char* kFs = "#version 330\nin vec2 uv;\nuniform sampler2D scene;\nout vec4 color;\nvoid main() { color = texture(scene, uv); }\n";

    void Release() {
        if (S.fbo) glDeleteFramebuffers(1, &S.fbo);
        if (S.resolveFbo) glDeleteFramebuffers(1, &S.resolveFbo);
        GLuint rb[] = { S.color, S.depth };
        if (S.color || S.depth) glDeleteRenderbuffers(2, rb);
        if (S.resolveTex) glDeleteTextures(1, &S.resolveTex);
        S.fbo = S.color = S.depth = S.resolveFbo = S.resolveTex = 0;
        S.w = S.h = 0;
    }

    bool Create(int w, int h) {
        Release();
        glGenFramebuffers(1, &S.fbo);
        glGenRenderbuffers(1, &S.color);
        glGenRenderbuffers(1, &S.depth);
        glBindRenderbuffer(GL_RENDERBUFFER, S.color);
        if (S.samples > 1) glRenderbufferStorageMultisample(GL_RENDERBUFFER, S.samples, GL_RGBA8, w, h);
        else glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA8, w, h);
        glBindRenderbuffer(GL_RENDERBUFFER, S.depth);
        if (S.samples > 1) glRenderbufferStorageMultisample(GL_RENDERBUFFER, S.samples, GL_DEPTH32F_STENCIL8, w, h);
        else glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH32F_STENCIL8, w, h);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, S.fbo);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, S.color);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, S.depth);
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        // The resolve target: the same RGBA8 format, single-sampled (a valid multisample resolve on any driver).
        glGenTextures(1, &S.resolveTex);
        glBindTexture(GL_TEXTURE_2D, S.resolveTex);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glBindTexture(GL_TEXTURE_2D, 0);
        glGenFramebuffers(1, &S.resolveFbo);
        glBindFramebuffer(GL_FRAMEBUFFER, S.resolveFbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, S.resolveTex, 0);
        if (status == GL_FRAMEBUFFER_COMPLETE) status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            logger.logf("[SceneTarget] framebuffer incomplete (0x%x): reverse-Z off\n", status);
            Release();
            return false;
        }
        S.w = w; S.h = h;
        logger.logf("[SceneTarget] %dx%d, %d samples, float depth, reverse-Z\n", w, h, S.samples);
        return true;
    }

    // glClipControl is core in 4.5 or ARB_clip_control: a non-null proc address alone does not prove the
    // context exposes it.
    bool ClipControlSupported() { return GL_Caps().clipControl; } // 4.5 or ARB_clip_control (glcaps, #286)

    bool Init() {
        if (S.tried) return S.ok;
        S.tried = true;
        const char* env = getenv("RAYLIBGAME_REVZ");
        if (env && env[0] == '0') { logger.logf("[SceneTarget] off (RAYLIBGAME_REVZ=0)\n"); return false; }
        if (rlGetVersion() != RL_OPENGL_43 || !ClipControlSupported()) { logger.logf("[SceneTarget] no clip control: reverse-Z off\n"); return false; }
        S.clipControl = (PfnClipControl)glfwGetProcAddress("glClipControl");
        if (!S.clipControl) S.clipControl = (PfnClipControl)glfwGetProcAddress("glClipControlARB");
        if (!S.clipControl) { logger.logf("[SceneTarget] no glClipControl: reverse-Z off\n"); return false; }
        const unsigned vs = rlCompileShader(kVs, RL_VERTEX_SHADER), fs = rlCompileShader(kFs, RL_FRAGMENT_SHADER);
        S.prog = (vs && fs) ? rlLoadShaderProgram(vs, fs) : 0;
        if (!S.prog) { logger.logf("[SceneTarget] resolve shader failed: reverse-Z off\n"); return false; }
        glGenVertexArrays(1, &S.vao);
        GLint samples = 0;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glGetIntegerv(GL_SAMPLES, &samples);
        S.samples = samples;
        S.ok = Create(GetRenderWidth(), GetRenderHeight());
        return S.ok;
    }
}

bool SceneTarget_Available() { return Init(); }

bool SceneTarget_Begin(Color clear) {
    if (!Init()) return false;
    rlDrawRenderBatchActive();
    const int w = GetRenderWidth(), h = GetRenderHeight();
    if ((w != S.w || h != S.h) && w > 0 && h > 0 && !Create(w, h)) { S.ok = false; return false; }
    glBindFramebuffer(GL_FRAMEBUFFER, S.fbo);
    S.clipControl(kLowerLeft, kZeroToOne);
    glDepthFunc(GL_GREATER);
    glDepthMask(GL_TRUE);
    glClearDepth(0.0);
    glClearStencil(0);
    glClearColor(clear.r / 255.0f, clear.g / 255.0f, clear.b / 255.0f, clear.a / 255.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
    return true;
}

void SceneTarget_ApplyProjection() {
    // From the standard projection BeginMode3D set: m0 = f / aspect, m5 = f, near = m14 / (m10 - 1).
    const Matrix p = rlGetMatrixProjection();
    const float n = p.m14 / (p.m10 - 1.0f);
    Matrix r = { 0 };
    r.m0 = p.m0;
    r.m5 = p.m5;
    r.m8 = p.m8; r.m9 = p.m9; // off-centre terms (none in practice)
    r.m14 = n;                // z' = n, w' = -z: depth n / -z, 1 at the near plane, 0 at infinity
    r.m11 = -1.0f;
    rlSetMatrixProjection(r);
}

void SceneTarget_End() {
    if (!S.fbo) return;
    rlDrawRenderBatchActive();
    // Resolve (or copy) into the single-sample texture: same RGBA8 format and size, valid everywhere.
    glBindFramebuffer(GL_READ_FRAMEBUFFER, S.fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, S.resolveFbo);
    glBlitFramebuffer(0, 0, S.w, S.h, 0, 0, S.w, S.h, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    S.clipControl(kLowerLeft, kNegOneToOne);
    glDepthFunc(GL_LEQUAL);
    glClearDepth(1.0);
    // Draw it to the window.
    glViewport(0, 0, S.w, S.h);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    glUseProgram(S.prog);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, S.resolveTex);
    glUniform1i(glGetUniformLocation(S.prog, "scene"), 0);
    glBindVertexArray(S.vao);
    glDrawArrays(GL_TRIANGLES, 0, 3);
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
    glUseProgram(0);
    // Back to rlgl's defaults (blend on, depth test as raylib leaves it after EndMode3D: off).
    glEnable(GL_BLEND);
    glEnable(GL_CULL_FACE);
}

void SceneTarget_Unload() {
    if (S.ok) Release();
    if (S.prog) rlUnloadShaderProgram(S.prog);
    if (S.vao) glDeleteVertexArrays(1, &S.vao);
    S = State{};
}
