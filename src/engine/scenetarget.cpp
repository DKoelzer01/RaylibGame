#include "scenetarget.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <rlgl.h>
#include <cstdlib>

extern "C" void* glfwGetProcAddress(const char* name); // raylib links GLFW

namespace {
    typedef void (__stdcall* PfnClipControl)(GLenum origin, GLenum depth);
    constexpr GLenum kLowerLeft = 0x8CA1, kNegOneToOne = 0x935E, kZeroToOne = 0x935F;

    struct State {
        bool tried = false, ok = false;
        PfnClipControl clipControl = nullptr;
        GLuint fbo = 0, color = 0, depth = 0;
        int w = 0, h = 0, samples = 0;
    } S;

    void Release() {
        if (S.fbo) glDeleteFramebuffers(1, &S.fbo);
        GLuint rb[] = { S.color, S.depth };
        if (S.color || S.depth) glDeleteRenderbuffers(2, rb);
        S.fbo = S.color = S.depth = 0;
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
        const GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
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

    bool Init() {
        if (S.tried) return S.ok;
        S.tried = true;
        const char* env = getenv("RAYLIBGAME_REVZ");
        if (env && env[0] == '0') { logger.logf("[SceneTarget] off (RAYLIBGAME_REVZ=0)\n"); return false; }
        if (rlGetVersion() != RL_OPENGL_43) return false;
        S.clipControl = (PfnClipControl)glfwGetProcAddress("glClipControl");
        if (!S.clipControl) S.clipControl = (PfnClipControl)glfwGetProcAddress("glClipControlARB");
        if (!S.clipControl) { logger.logf("[SceneTarget] no glClipControl: reverse-Z off\n"); return false; }
        GLint samples = 0;
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glGetIntegerv(GL_SAMPLES, &samples);
        S.samples = samples;
        S.ok = Create(GetRenderWidth(), GetRenderHeight());
        return S.ok;
    }
}

bool SceneTarget_Available() { return Init(); }

void SceneTarget_Begin(Color clear) {
    if (!Init()) return;
    rlDrawRenderBatchActive();
    const int w = GetRenderWidth(), h = GetRenderHeight();
    if ((w != S.w || h != S.h) && w > 0 && h > 0 && !Create(w, h)) { S.ok = false; return; }
    glBindFramebuffer(GL_FRAMEBUFFER, S.fbo);
    S.clipControl(kLowerLeft, kZeroToOne);
    glDepthFunc(GL_GREATER);
    glDepthMask(GL_TRUE);
    glClearDepth(0.0);
    glClearStencil(0);
    glClearColor(clear.r / 255.0f, clear.g / 255.0f, clear.b / 255.0f, clear.a / 255.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
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
    if (!S.ok) return;
    rlDrawRenderBatchActive();
    glBindFramebuffer(GL_READ_FRAMEBUFFER, S.fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
    glBlitFramebuffer(0, 0, S.w, S.h, 0, 0, S.w, S.h, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    S.clipControl(kLowerLeft, kNegOneToOne);
    glDepthFunc(GL_LEQUAL);
    glClearDepth(1.0);
}

void SceneTarget_Unload() {
    if (S.ok) Release();
    S = State{};
}
