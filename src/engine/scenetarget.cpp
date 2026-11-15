#include "scenetarget.h"
#include "glcaps.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <rlgl.h>
#include "world.h" // GameFrameTime (auto-exposure adapts on the simulation clock: recordings stay reproducible)
#include <algorithm>
#include <cstdio>
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
        // HDR (#131): an RGBA16F scene, so values over 1 survive to the composite, which tonemaps them (and auto-exposure
        // measures the scene, kCsExposure). RAYLIBGAME_HDR=0: the old RGBA8 target, copied as it is.
        bool hdr = true, autoExposure = true;
        int tonemap = 0;                       // tonemap: 0 knee (today's look up to the knee, highlights roll off), 1 ACES, 2 clip
        float knee = 0.8f, key = 0.35f, eMin = 0.35f, eMax = 1.0f, tau = 1.5f;
        GLuint expProg = 0, expBuf = 0;         // the exposure (and the measured luminance) in a 2-float SSBO
        GLint locTonemap = -1, locKnee = -1, locDt = -1, locKey = -1, locEMin = -1, locEMax = -1, locTau = -1;
    } S;
    constexpr GLuint kExposureBinding = 5;

    // A fullscreen triangle sampling the resolved scene. Drawing (not blitting) to the window works whatever
    // its format and sample count: a multisampled blit needs identical formats and sample counts on both
    // sides, which only some drivers' default framebuffers happen to match (PR #156 review).
    const char* kVs = "#version 330\nout vec2 uv;\nvoid main() { vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2); uv = p; gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0); }\n";
    const char* kFs = "#version 330\nin vec2 uv;\nuniform sampler2D scene;\nout vec4 color;\nvoid main() { color = texture(scene, uv); }\n";
    // HDR composite (#131). The world shaders still write display-encoded colour (pow 1/2.2, and the atmosphere, sky,
    // particles and raylib's default shader write display colour as they always did), now unclamped in RGBA16F. So: decode,
    // expose, tonemap, encode. Blending and MSAA stay where they were (display space), so nothing else changes.
    const char* kFsHdr = R"(#version 430
in vec2 uv;
uniform sampler2D scene;
uniform int tonemap;
uniform float knee;
uniform int nanDebug; // RAYLIBGAME_NAN_CHECK=1: non-finite pixels in magenta (#402)
layout(std430, binding = 5) readonly buffer Exposure { float exposure; float lum; float bad; };
out vec4 color;
vec3 kneeCurve(vec3 x) { // identity up to the knee, then an exponential shoulder towards 1 (slope 1 at the knee)
    vec3 over = max(x - knee, 0.0);
    return min(x, vec3(knee)) + (1.0 - knee) * (1.0 - exp(-over / (1.0 - knee)));
}
vec3 aces(vec3 x) { // Narkowicz's fit, with its usual 0.6 pre-scale
    x *= 0.6;
    return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}
void main() {
    vec3 c = texelFetch(scene, ivec2(gl_FragCoord.xy), 0).rgb;
    if (any(isnan(c)) || any(isinf(c))) { if (nanDebug == 1) { color = vec4(1.0, 0.0, 1.0, 1.0); return; } c = vec3(0.0); } // ~0.1% of pixels are NaN / Inf (#402); RGBA8 stored them as 0
    vec3 lin = pow(max(c, vec3(0.0)), vec3(2.2)) * (exposure > 0.0 ? exposure : 1.0);
    vec3 t = tonemap == 1 ? aces(lin) : tonemap == 2 ? min(lin, vec3(1.0)) : kneeCurve(lin);
    color = vec4(pow(t, vec3(1.0 / 2.2)), 1.0);
}
)";
    // Auto-exposure (#131): the log-average luminance of a 64 x 36 grid of scene samples (one workgroup) -> a target
    // exposure key / L, clamped to [eMin, eMax] (eMax 1: it never brightens, so nights stay as dark as tuned, #369),
    // approached with time constant tau. Not a mip chain: a few NaN / Inf pixels (written as 0 by the old RGBA8 target)
    // poisoned the whole average; here non-finite samples are skipped (and counted).
    const char* kCsExposure = R"(#version 430
layout(local_size_x = 256) in;
uniform sampler2D scene;
uniform float dt, key, eMin, eMax, tau;
layout(std430, binding = 5) buffer Exposure { float exposure; float lum; float bad; };
shared float sumLog[256];
shared float count[256];
shared float badCount[256];
void main() {
    const int GX = 64, GY = 36;
    ivec2 size = textureSize(scene, 0);
    float s = 0.0, n = 0.0, b = 0.0;
    for (int i = int(gl_LocalInvocationIndex); i < GX * GY; i += 256) {
        ivec2 q = ivec2((vec2(i % GX, i / GX) + 0.5) / vec2(GX, GY) * vec2(size));
        vec3 c = texelFetch(scene, q, 0).rgb;
        if (any(isnan(c)) || any(isinf(c))) { b += 1.0; continue; }
        vec3 lin = pow(max(c, vec3(0.0)), vec3(2.2));
        s += log(max(dot(lin, vec3(0.2126, 0.7152, 0.0722)), 1e-4));
        n += 1.0;
    }
    sumLog[gl_LocalInvocationIndex] = s;
    count[gl_LocalInvocationIndex] = n;
    badCount[gl_LocalInvocationIndex] = b;
    barrier();
    for (uint k = 128u; k > 0u; k >>= 1) {
        if (gl_LocalInvocationIndex < k) {
            sumLog[gl_LocalInvocationIndex] += sumLog[gl_LocalInvocationIndex + k];
            count[gl_LocalInvocationIndex] += count[gl_LocalInvocationIndex + k];
            badCount[gl_LocalInvocationIndex] += badCount[gl_LocalInvocationIndex + k];
        }
        barrier();
    }
    if (gl_LocalInvocationIndex == 0u) {
        float badN = badCount[0], goodN = count[0];
        float L = goodN > 0.0 ? exp(sumLog[0] / goodN) : key;
        float target = clamp(key / L, eMin, eMax);
        exposure = exposure <= 0.0 ? target : mix(target, exposure, exp(-dt / tau));
        lum = L;
        bad = badN;
    }
}
)";

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
        const GLenum colorFormat = S.hdr ? GL_RGBA16F : GL_RGBA8;
        if (S.samples > 1) glRenderbufferStorageMultisample(GL_RENDERBUFFER, S.samples, colorFormat, w, h);
        else glRenderbufferStorage(GL_RENDERBUFFER, colorFormat, w, h);
        glBindRenderbuffer(GL_RENDERBUFFER, S.depth);
        if (S.samples > 1) glRenderbufferStorageMultisample(GL_RENDERBUFFER, S.samples, GL_DEPTH32F_STENCIL8, w, h);
        else glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH32F_STENCIL8, w, h);
        glBindRenderbuffer(GL_RENDERBUFFER, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, S.fbo);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, S.color);
        glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, S.depth);
        GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
        // The resolve target: the same format, single-sampled (a valid multisample resolve on any driver). HDR: with a
        // mip chain, whose top level is the scene's average for auto-exposure.
        glGenTextures(1, &S.resolveTex);
        glBindTexture(GL_TEXTURE_2D, S.resolveTex);
        if (S.hdr) glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_HALF_FLOAT, nullptr);
        else glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
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
        logger.logf("[SceneTarget] %dx%d, %d samples, float depth, reverse-Z, %s\n", w, h, S.samples,
                    !S.hdr ? "RGBA8 (RAYLIBGAME_HDR=0)" : S.tonemap == 1 ? "HDR, ACES" : S.tonemap == 2 ? "HDR, clip" : "HDR, knee tonemap");
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
        // HDR switches (#131): RAYLIBGAME_HDR=0 (RGBA8, no tonemap), RAYLIBGAME_TONEMAP=knee|aces|clip,
        // RAYLIBGAME_AUTO_EXPOSURE=0, RAYLIBGAME_EXPOSURE_KEY / _MIN / _MAX, RAYLIBGAME_TONEMAP_KNEE.
        auto envF = [](const char* n, float d) { const char* v = getenv(n); return v && v[0] ? (float)atof(v) : d; };
        if (const char* h = getenv("RAYLIBGAME_HDR"); h && h[0] == '0') S.hdr = false;
        if (const char* t = getenv("RAYLIBGAME_TONEMAP")) S.tonemap = !strcmp(t, "aces") ? 1 : !strcmp(t, "clip") ? 2 : 0;
        if (const char* a = getenv("RAYLIBGAME_AUTO_EXPOSURE"); a && a[0] == '0') S.autoExposure = false;
        S.key = envF("RAYLIBGAME_EXPOSURE_KEY", S.key);
        S.eMin = envF("RAYLIBGAME_EXPOSURE_MIN", S.eMin);
        S.eMax = envF("RAYLIBGAME_EXPOSURE_MAX", S.eMax);
        S.knee = std::clamp(envF("RAYLIBGAME_TONEMAP_KNEE", S.knee), 0.1f, 0.99f);
        if (S.hdr) {
            const unsigned vs = rlCompileShader(kVs, RL_VERTEX_SHADER), fs = rlCompileShader(kFsHdr, RL_FRAGMENT_SHADER);
            S.prog = (vs && fs) ? rlLoadShaderProgram(vs, fs) : 0;
            const unsigned cs = rlCompileShader(kCsExposure, RL_COMPUTE_SHADER);
            S.expProg = cs ? rlLoadComputeShaderProgram(cs) : 0;
            if (!S.prog || !S.expProg) {
                logger.logf("[SceneTarget] HDR shaders failed: RGBA8 target\n");
                if (S.prog) rlUnloadShaderProgram(S.prog);
                if (S.expProg) rlUnloadShaderProgram(S.expProg);
                S.prog = S.expProg = 0;
                S.hdr = false;
            } else {
                S.locTonemap = glGetUniformLocation(S.prog, "tonemap");
                S.locKnee = glGetUniformLocation(S.prog, "knee");
                S.locDt = glGetUniformLocation(S.expProg, "dt");
                S.locKey = glGetUniformLocation(S.expProg, "key");
                S.locEMin = glGetUniformLocation(S.expProg, "eMin");
                S.locEMax = glGetUniformLocation(S.expProg, "eMax");
                S.locTau = glGetUniformLocation(S.expProg, "tau");
                const float init[3] = { S.autoExposure ? 0.0f : 1.0f, 0.0f, 0.0f }; // 0: the first frame takes its target at once
                glGenBuffers(1, &S.expBuf);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, S.expBuf);
                glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof init, init, GL_DYNAMIC_COPY);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
            }
        }
        if (!S.hdr) {
            const unsigned vs = rlCompileShader(kVs, RL_VERTEX_SHADER), fs = rlCompileShader(kFs, RL_FRAGMENT_SHADER);
            S.prog = (vs && fs) ? rlLoadShaderProgram(vs, fs) : 0;
        }
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
    if (S.hdr) {
        glBindBufferBase(GL_SHADER_STORAGE_BUFFER, kExposureBinding, S.expBuf);
        if (S.autoExposure) {
            // Auto-exposure: one workgroup measures the resolved scene and updates the exposure.
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, S.resolveTex);
            glUseProgram(S.expProg);
            glUniform1i(glGetUniformLocation(S.expProg, "scene"), 0);
            glUniform1f(S.locDt, std::clamp(GameFrameTime(), 0.0f, 0.25f));
            glUniform1f(S.locKey, S.key);
            glUniform1f(S.locEMin, S.eMin);
            glUniform1f(S.locEMax, S.eMax);
            glUniform1f(S.locTau, S.tau);
            glDispatchCompute(1, 1, 1);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            glBindTexture(GL_TEXTURE_2D, 0);
            // RAYLIBGAME_EXPOSURE_LOG=1: the exposure and the scene luminance every ~2 s (reads the buffer back: a stall).
            static const bool logExp = getenv("RAYLIBGAME_EXPOSURE_LOG") && getenv("RAYLIBGAME_EXPOSURE_LOG")[0] == '1';
            static int frame = 0;
            if (logExp && ++frame % 120 == 0) {
                float v[3] = { 0, 0, 0 };
                glMemoryBarrier(GL_BUFFER_UPDATE_BARRIER_BIT);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, S.expBuf);
                glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof v, v);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
                printf("exposure: %.3f (scene log-average luminance %.4f; %.0f of 2304 samples NaN/Inf)\n", v[0], v[1], v[2]);
            }
        }
    }
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
    if (S.hdr) {
        static const int nanCheck = getenv("RAYLIBGAME_NAN_CHECK") && getenv("RAYLIBGAME_NAN_CHECK")[0] == '1';
        glUniform1i(S.locTonemap, S.tonemap); glUniform1f(S.locKnee, S.knee);
        glUniform1i(glGetUniformLocation(S.prog, "nanDebug"), nanCheck);
    }
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
    if (S.expProg) rlUnloadShaderProgram(S.expProg);
    if (S.expBuf) glDeleteBuffers(1, &S.expBuf);
    if (S.vao) glDeleteVertexArrays(1, &S.vao);
    S = State{};
}
