#include "scenetarget.h"
#include "gpumem.h"
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
        // Contact shadows (#137): the scene's depth resolved to a texture; the composite marches it toward the sun.
        GLuint depthFbo = 0, depthTex = 0;
        bool contact = true;                      // RAYLIBGAME_CONTACT_SHADOWS=0: off
        float proj[3] = { 1.0f, 1.0f, 0.1f };     // the projection's m0, m5 and near plane (SceneTarget_ApplyProjection)
        float sunView[3] = { 0.0f, 1.0f, 0.0f };  // toward the sun, view space
        float contactStrength = 0.0f;             // 0: none this frame (night, overcast, no shadows)
        GLuint prog = 0, vao = 0;
        int w = 0, h = 0, samples = 0;
        // HDR (#131): an RGBA16F scene, so values over 1 survive to the composite, which tonemaps them (and auto-exposure
        // measures the scene, kCsExposure). RAYLIBGAME_HDR=0: the old RGBA8 target, copied as it is.
        bool hdr = true, autoExposure = true;
        int tonemap = 0;                       // tonemap: 0 knee (today's look up to the knee, highlights roll off), 1 ACES, 2 clip
        float knee = 0.8f, key = 0.35f, eMin = 0.35f, eMax = 1.0f, tau = 1.5f;
        GLuint expProg = 0, expBuf = 0;         // the exposure (and the measured luminance) in a 2-float SSBO
        GLint locTonemap = -1, locKnee = -1, locDt = -1, locKey = -1, locEMin = -1, locEMax = -1, locTau = -1;
        // Bloom (#140): a dual-filter (Kawase) chain on the linear scene, kBloomLevels from half size, mixed into the
        // composite before exposure and tonemapping (energy-conserving: lerp toward the blur by bloomStrength, no threshold).
        bool bloom = true;
        float bloomStrength = 0.04f;
        static constexpr int kBloomLevels = 6;
        GLuint bloomTex[kBloomLevels] = {}, bloomFbo[kBloomLevels] = {};
        int bloomW[kBloomLevels] = {}, bloomH[kBloomLevels] = {};
        GLuint downProg = 0, upProg = 0;
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
uniform sampler2D bloomTex; // #140: the blur chain's top level (linear), 0 strength = none
uniform float bloomStrength;
uniform float bloomLevels;
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
// Contact shadows (#137): a short march toward the sun through the resolved depth catches what the shadow maps' texels
// miss up close (stones on the ground, grass bases, part edges, feet). View-space positions from reverse-Z depth
// (d = near / -z); sun-facing pixels within 30 m only; darkens the lit colour by contactStrength (0: off).
uniform sampler2D depthTex;
uniform float contactStrength;
uniform vec3 sunView;   // toward the sun, view space
uniform vec3 projP;     // m0, m5, near
vec3 viewAt(vec2 uv, float d) { float z = -projP.z / max(d, 1e-9); vec2 ndc = uv * 2.0 - 1.0; return vec3(ndc.x * -z / projP.x, ndc.y * -z / projP.y, z); }
float contactShadow() {
    vec2 size = vec2(textureSize(depthTex, 0));
    vec2 uv = gl_FragCoord.xy / size;
    float d0 = texelFetch(depthTex, ivec2(gl_FragCoord.xy), 0).r;
    vec3 p = viewAt(uv, max(d0, 1e-6));
    vec3 dn = cross(dFdx(p), dFdy(p));                          // before any branch: derivatives need uniform flow
    if (d0 <= 0.0) return 1.0;                                  // sky
    float dist = -p.z;
    if (dist > 30.0 || dot(dn, dn) < 1e-20) return 1.0;
    vec3 n = normalize(dn);
    if (dot(n, -p) < 0.0) n = -n;                               // toward the camera
    float facing = dot(n, sunView);
    if (facing <= 0.03) return 1.0;                             // turned away: the sun term is dark already
    float len = mix(0.3, 0.9, clamp(dist / 20.0, 0.0, 1.0));
    float jitter = fract(52.9829189 * fract(dot(gl_FragCoord.xy, vec2(0.06711056, 0.00583715))));
    vec3 start = p + n * (0.01 + 0.002 * dist);
    float occ = 0.0;
    for (int i = 1; i <= 12; ++i) {
        vec3 q = start + sunView * ((float(i) - jitter) / 12.0 * len);
        if (q.z >= -projP.z) break;                             // behind the near plane
        vec2 s = vec2(q.x * projP.x, q.y * projP.y) / -q.z * 0.5 + 0.5;
        if (any(lessThan(s, vec2(0.0))) || any(greaterThan(s, vec2(1.0)))) break;
        float zs = -projP.z / max(texture(depthTex, s).r, 1e-9);
        float diff = zs - q.z;                                  // > 0: the scene there is in front of the ray point
        if (diff > 0.003 + 0.0015 * dist && diff < 0.12 + 0.01 * dist) { occ = 1.0; break; }
    }
    occ *= 1.0 - smoothstep(20.0, 30.0, dist);
    return 1.0 - occ * contactStrength * smoothstep(0.03, 0.25, facing);
}
void main() {
    vec3 c = texelFetch(scene, ivec2(gl_FragCoord.xy), 0).rgb;
    if (any(isnan(c)) || any(isinf(c))) { if (nanDebug == 1) { color = vec4(1.0, 0.0, 1.0, 1.0); return; } c = vec3(0.0); } // ~0.1% of pixels are NaN / Inf (#402); RGBA8 stored them as 0
    vec3 lin = pow(max(c, vec3(0.0)), vec3(2.2));
    if (contactStrength > 0.0) lin *= contactShadow();
    if (bloomStrength > 0.0) lin = mix(lin, texture(bloomTex, gl_FragCoord.xy / vec2(textureSize(scene, 0))).rgb / bloomLevels, bloomStrength); // #140
    lin = lin * (exposure > 0.0 ? exposure : 1.0);
    vec3 t = tonemap == 1 ? aces(lin) : tonemap == 2 ? min(lin, vec3(1.0)) : kneeCurve(lin);
    color = vec4(pow(t, vec3(1.0 / 2.2)), 1.0);
}
)";
    // Bloom blur passes (#140): the dual filter's downsample (centre x4 + four diagonal taps; the first level also decodes
    // the display-encoded scene to linear and drops NaN/Inf) and upsample (a tent of 8 taps), added level by level.
    const char* kFsDown = R"(#version 330
in vec2 uv;
uniform sampler2D src;
uniform vec2 texel;
uniform int decode;
out vec4 color;
vec3 fetch(vec2 p) {
    vec3 c = texture(src, p).rgb;
    if (any(isnan(c)) || any(isinf(c))) c = vec3(0.0);
    c = max(c, vec3(0.0));
    if (decode == 0) return c;
    c = pow(c, vec3(2.2));
    // Fireflies (#452): a few sub-pixel sun glints on water thousands of times brighter than the scene bloomed into
    // big white discs. Cap a pixel's luminance before it spreads.
    float l = dot(c, vec3(0.2126, 0.7152, 0.0722));
    return l > 16.0 ? c * (16.0 / l) : c;
}
// The first level also weights taps by 1 / (1 + luminance) (Karis average): a lone bright pixel can't dominate.
vec3 tap(vec2 p, float w, inout float wsum) {
    vec3 c = fetch(p);
    float k = decode == 1 ? w / (1.0 + dot(c, vec3(0.2126, 0.7152, 0.0722))) : w;
    wsum += k;
    return c * k;
}
void main() {
    vec2 h = texel * 0.5;
    float wsum = 0.0;
    vec3 s = tap(uv, 4.0, wsum) + tap(uv - h, 1.0, wsum) + tap(uv + h, 1.0, wsum) + tap(uv + vec2(h.x, -h.y), 1.0, wsum) + tap(uv - vec2(h.x, -h.y), 1.0, wsum);
    color = vec4(s / max(wsum, 1e-4), 1.0);
}
)";
    const char* kFsUp = R"(#version 330
in vec2 uv;
uniform sampler2D src;
uniform vec2 texel;
out vec4 color;
void main() {
    vec2 h = texel * 0.5;
    vec3 s = texture(src, uv + vec2(-2.0 * h.x, 0.0)).rgb + texture(src, uv + vec2(2.0 * h.x, 0.0)).rgb
           + texture(src, uv + vec2(0.0, -2.0 * h.y)).rgb + texture(src, uv + vec2(0.0, 2.0 * h.y)).rgb
           + 2.0 * (texture(src, uv + h).rgb + texture(src, uv - h).rgb + texture(src, uv + vec2(h.x, -h.y)).rgb + texture(src, uv + vec2(-h.x, h.y)).rgb);
    color = vec4(s / 12.0, 1.0);
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
        if (S.depthFbo) glDeleteFramebuffers(1, &S.depthFbo);
        if (S.depthTex) glDeleteTextures(1, &S.depthTex);
        S.fbo = S.color = S.depth = S.resolveFbo = S.resolveTex = S.depthFbo = S.depthTex = 0;
        S.w = S.h = 0;
        for (int i = 0; i < State::kBloomLevels; ++i) {
            if (S.bloomFbo[i]) glDeleteFramebuffers(1, &S.bloomFbo[i]);
            if (S.bloomTex[i]) glDeleteTextures(1, &S.bloomTex[i]);
            S.bloomFbo[i] = S.bloomTex[i] = 0;
        }
    }

    // The bloom levels (#140): half, quarter, ... size, RGBA16F, linear filtering.
    void CreateBloom(int w, int h) {
        for (int i = 0; i < State::kBloomLevels; ++i) {
            w = std::max(1, w / 2); h = std::max(1, h / 2);
            S.bloomW[i] = w; S.bloomH[i] = h;
            glGenTextures(1, &S.bloomTex[i]);
            glBindTexture(GL_TEXTURE_2D, S.bloomTex[i]);
            GpuMemTag gm196_("bloom"); // VRAM accounting (#302)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_HALF_FLOAT, nullptr);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glGenFramebuffers(1, &S.bloomFbo[i]);
            glBindFramebuffer(GL_FRAMEBUFFER, S.bloomFbo[i]);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, S.bloomTex[i], 0);
        }
        glBindTexture(GL_TEXTURE_2D, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
    }

    // The bloom chain (#140): down from the resolved scene (decoded to linear) through the levels, then back up adding
    // each smaller level into the next larger one; level 0 ends as the sum (the composite divides by the count).
    void RunBloom() {
        auto pass = [&](GLuint prog, GLuint src, int srcW, int srcH, int dst, bool decode) {
            glBindFramebuffer(GL_FRAMEBUFFER, S.bloomFbo[dst]);
            glViewport(0, 0, S.bloomW[dst], S.bloomH[dst]);
            glUseProgram(prog);
            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_2D, src);
            glUniform1i(glGetUniformLocation(prog, "src"), 0);
            glUniform2f(glGetUniformLocation(prog, "texel"), 1.0f / srcW, 1.0f / srcH);
            const GLint dl = glGetUniformLocation(prog, "decode");
            if (dl >= 0) glUniform1i(dl, decode ? 1 : 0);
            glDrawArrays(GL_TRIANGLES, 0, 3);
        };
        glBindVertexArray(S.vao);
        glDisable(GL_DEPTH_TEST);
        glDisable(GL_CULL_FACE);
        glDisable(GL_BLEND);
        glBindTexture(GL_TEXTURE_2D, S.resolveTex); // sampled with filtering here (the composite fetches texels)
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        pass(S.downProg, S.resolveTex, S.w, S.h, 0, true);
        for (int i = 1; i < State::kBloomLevels; ++i) pass(S.downProg, S.bloomTex[i - 1], S.bloomW[i - 1], S.bloomH[i - 1], i, false);
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE, GL_ONE);
        for (int i = State::kBloomLevels - 1; i > 0; --i) pass(S.upProg, S.bloomTex[i], S.bloomW[i], S.bloomH[i], i - 1, false);
        glDisable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); // rlgl's default
        glBindFramebuffer(GL_FRAMEBUFFER, 0);
        glBindVertexArray(0);
    }

    bool Create(int w, int h) {
        Release();
        glGenFramebuffers(1, &S.fbo);
        glGenRenderbuffers(1, &S.color);
        glGenRenderbuffers(1, &S.depth);
        glBindRenderbuffer(GL_RENDERBUFFER, S.color);
        const GLenum colorFormat = S.hdr ? GL_RGBA16F : GL_RGBA8;
        GpuMemTag gm249_("scene target"); // VRAM accounting (#302)
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
        GpuMemTag gm263_("scene target"); // VRAM accounting (#302)
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
        if (S.hdr && S.contact && status == GL_FRAMEBUFFER_COMPLETE) { // #137: the resolved depth (same format: a valid blit)
            glGenTextures(1, &S.depthTex);
            glBindTexture(GL_TEXTURE_2D, S.depthTex);
            GpuMemTag gmDepth_("scene target"); // VRAM accounting (#302)
            glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH32F_STENCIL8, w, h, 0, GL_DEPTH_STENCIL, GL_FLOAT_32_UNSIGNED_INT_24_8_REV, nullptr);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_COMPARE_MODE, GL_NONE);
            glBindTexture(GL_TEXTURE_2D, 0);
            glGenFramebuffers(1, &S.depthFbo);
            glBindFramebuffer(GL_FRAMEBUFFER, S.depthFbo);
            glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_TEXTURE_2D, S.depthTex, 0);
            glDrawBuffer(GL_NONE); glReadBuffer(GL_NONE);
            if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
                LOG_WARN("SceneTarget", "depth resolve target incomplete: no contact shadows\n");
                glDeleteFramebuffers(1, &S.depthFbo); glDeleteTextures(1, &S.depthTex); S.depthFbo = S.depthTex = 0;
            }
            glBindFramebuffer(GL_FRAMEBUFFER, 0);
        }
        if (status != GL_FRAMEBUFFER_COMPLETE) {
            LOG_WARN("SceneTarget", "framebuffer incomplete (0x%x): reverse-Z off\n", status);
            Release();
            return false;
        }
        S.w = w; S.h = h;
        if (S.hdr && S.bloom && S.downProg && S.upProg) CreateBloom(w, h); // #140
        LOG_INFO("SceneTarget", "%dx%d, %d samples, float depth, reverse-Z, %s\n", w, h, S.samples,
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
        if (env && env[0] == '0') { LOG_INFO("SceneTarget", "off (RAYLIBGAME_REVZ=0)\n"); return false; }
        if (rlGetVersion() != RL_OPENGL_43 || !ClipControlSupported()) { LOG_INFO("SceneTarget", "no clip control: reverse-Z off\n"); return false; }
        S.clipControl = (PfnClipControl)glfwGetProcAddress("glClipControl");
        if (!S.clipControl) S.clipControl = (PfnClipControl)glfwGetProcAddress("glClipControlARB");
        if (!S.clipControl) { LOG_INFO("SceneTarget", "no glClipControl: reverse-Z off\n"); return false; }
        // HDR switches (#131): RAYLIBGAME_HDR=0 (RGBA8, no tonemap), RAYLIBGAME_TONEMAP=knee|aces|clip,
        // RAYLIBGAME_AUTO_EXPOSURE=0, RAYLIBGAME_EXPOSURE_KEY / _MIN / _MAX, RAYLIBGAME_TONEMAP_KNEE.
        auto envF = [](const char* n, float d) { const char* v = getenv(n); return v && v[0] ? (float)atof(v) : d; };
        if (const char* h = getenv("RAYLIBGAME_HDR"); h && h[0] == '0') S.hdr = false;
        if (const char* c = getenv("RAYLIBGAME_CONTACT_SHADOWS"); c && c[0] == '0') S.contact = false; // #137
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
                LOG_WARN("SceneTarget", "HDR shaders failed: RGBA8 target\n");
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
                // Bloom (#140): RAYLIBGAME_BLOOM=0 off, RAYLIBGAME_BLOOM_STRENGTH (0.04).
                if (const char* b = getenv("RAYLIBGAME_BLOOM"); b && b[0] == '0') S.bloom = false;
                S.bloomStrength = std::clamp(envF("RAYLIBGAME_BLOOM_STRENGTH", S.bloomStrength), 0.0f, 1.0f);
                if (S.bloom) {
                    const unsigned bvs = rlCompileShader(kVs, RL_VERTEX_SHADER), dfs = rlCompileShader(kFsDown, RL_FRAGMENT_SHADER);
                    S.downProg = (bvs && dfs) ? rlLoadShaderProgram(bvs, dfs) : 0;
                    const unsigned uvs = rlCompileShader(kVs, RL_VERTEX_SHADER), ufs = rlCompileShader(kFsUp, RL_FRAGMENT_SHADER);
                    S.upProg = (uvs && ufs) ? rlLoadShaderProgram(uvs, ufs) : 0;
                    if (!S.downProg || !S.upProg) { LOG_WARN("SceneTarget", "bloom shaders failed: no bloom\n"); S.bloom = false; }
                }
                const float init[3] = { S.autoExposure ? 0.0f : 1.0f, 0.0f, 0.0f }; // 0: the first frame takes its target at once
                glGenBuffers(1, &S.expBuf);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, S.expBuf);
                GpuMemTag gm342_("exposure"); // VRAM accounting (#302)
                glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof init, init, GL_DYNAMIC_COPY);
                glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
            }
        }
        if (!S.hdr) {
            const unsigned vs = rlCompileShader(kVs, RL_VERTEX_SHADER), fs = rlCompileShader(kFs, RL_FRAGMENT_SHADER);
            S.prog = (vs && fs) ? rlLoadShaderProgram(vs, fs) : 0;
        }
        if (!S.prog) { LOG_WARN("SceneTarget", "resolve shader failed: reverse-Z off\n"); return false; }
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
    S.proj[0] = p.m0; S.proj[1] = p.m5; S.proj[2] = n; // for the contact shadows' depth reconstruction (#137)
    Matrix r = { 0 };
    r.m0 = p.m0;
    r.m5 = p.m5;
    r.m8 = p.m8; r.m9 = p.m9; // off-centre terms (none in practice)
    r.m14 = n;                // z' = n, w' = -z: depth n / -z, 1 at the near plane, 0 at infinity
    r.m11 = -1.0f;
    rlSetMatrixProjection(r);
}

void SceneTarget_SetContact(Vector3 sunView, float strength) {
    const float l = sqrtf(sunView.x * sunView.x + sunView.y * sunView.y + sunView.z * sunView.z);
    if (l > 1e-6f) { S.sunView[0] = sunView.x / l; S.sunView[1] = sunView.y / l; S.sunView[2] = sunView.z / l; }
    S.contactStrength = std::clamp(strength, 0.0f, 1.0f);
}

void SceneTarget_End() {
    if (!S.fbo) return;
    rlDrawRenderBatchActive();
    // Resolve (or copy) into the single-sample texture: same RGBA8 format and size, valid everywhere.
    glBindFramebuffer(GL_READ_FRAMEBUFFER, S.fbo);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, S.resolveFbo);
    glBlitFramebuffer(0, 0, S.w, S.h, 0, 0, S.w, S.h, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    const bool contactOn = S.depthFbo && S.contactStrength > 0.0f;
    if (contactOn) { // #137: the depth too (a multisampled depth resolves to one of its samples)
        glBindFramebuffer(GL_DRAW_FRAMEBUFFER, S.depthFbo);
        glBlitFramebuffer(0, 0, S.w, S.h, 0, 0, S.w, S.h, GL_DEPTH_BUFFER_BIT, GL_NEAREST);
    }
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
        if (S.bloom && S.bloomTex[0]) RunBloom(); // #140 (sets its own viewports; the composite below sets the window's)
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
        const bool bloomOn = S.bloom && S.bloomTex[0];
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, bloomOn ? S.bloomTex[0] : 0);
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(glGetUniformLocation(S.prog, "bloomTex"), 1);
        glUniform1f(glGetUniformLocation(S.prog, "bloomStrength"), bloomOn ? S.bloomStrength : 0.0f);
        glUniform1f(glGetUniformLocation(S.prog, "bloomLevels"), (float)State::kBloomLevels);
        // Contact shadows (#137).
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, contactOn ? S.depthTex : 0);
        glActiveTexture(GL_TEXTURE0);
        glUniform1i(glGetUniformLocation(S.prog, "depthTex"), 2);
        glUniform1f(glGetUniformLocation(S.prog, "contactStrength"), contactOn ? S.contactStrength : 0.0f);
        glUniform3fv(glGetUniformLocation(S.prog, "sunView"), 1, S.sunView);
        glUniform3fv(glGetUniformLocation(S.prog, "projP"), 1, S.proj);
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
    if (S.downProg) rlUnloadShaderProgram(S.downProg);
    if (S.upProg) rlUnloadShaderProgram(S.upProg);
    if (S.expProg) rlUnloadShaderProgram(S.expProg);
    if (S.expBuf) glDeleteBuffers(1, &S.expBuf);
    if (S.vao) glDeleteVertexArrays(1, &S.vao);
    S = State{};
}
