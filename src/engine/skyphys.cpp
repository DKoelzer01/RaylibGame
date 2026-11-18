#include "skyphys.h"
#include "gpumem.h"
#include "utils/logger.h"

#include "external/glad.h"
#include <raylib.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <string>

namespace {
    const char* kVs = "#version 430\nout vec2 uv;\nvoid main() { vec2 p = vec2((gl_VertexID << 1) & 2, gl_VertexID & 2); uv = p; gl_Position = vec4(p * 2.0 - 1.0, 0.0, 1.0); }\n";

    // Transmittance to the top of the air along (height, view zenith).
    const char* kFsTransmittance = R"(
in vec2 uv;
out vec4 color;
void main()
{
    float vh, mu;
    SkyTransmittanceParams(uv, vh, mu);
    vec3 ro = vec3(0.0, vh, 0.0), rd = vec3(sqrt(max(1.0 - mu * mu, 0.0)), mu, 0.0);
    float t = max(SkyRaySphere(ro, rd, skyRt), 0.0);
    const int N = 40;
    vec3 od = vec3(0.0);
    for (int i = 0; i < N; ++i) {
        vec3 p = ro + rd * (t * (float(i) + 0.5) / float(N));
        od += SkyMediumAt(length(p) - skyR).ext;
    }
    color = vec4(exp(-od * (t / float(N))), 1.0);
}
)";

    // Multiple scattering (Hillaire 2020, 5.5): second-order light at (sun zenith, height) over 64 directions, and the
    // transfer factor f_ms; the infinite series sums to L2 / (1 - f_ms). Ground albedo 0.3.
    const char* kFsMultiScatter = R"(
in vec2 uv;
out vec4 color;
uniform sampler2D transmittanceLut;
vec3 SunTransmittance(vec3 p, vec3 sun)
{
    float h = length(p);
    return texture(transmittanceLut, SkyTransmittanceUv(h, dot(p / h, sun))).rgb;
}
void main()
{
    float cs = uv.x * 2.0 - 1.0;
    float h = skyR + clamp(uv.y, 0.002, 0.998) * (skyRt - skyR);
    vec3 sun = vec3(0.0, cs, sqrt(max(1.0 - cs * cs, 0.0)));
    vec3 ro = vec3(0.0, h, 0.0);
    vec3 L = vec3(0.0), F = vec3(0.0);
    const int D = 8, N = 20;
    for (int j = 0; j < D; ++j)
    for (int i = 0; i < D; ++i) {
        float ph = 6.2831853 * (float(i) + 0.5) / float(D);
        float ct = 1.0 - 2.0 * (float(j) + 0.5) / float(D), st = sqrt(max(1.0 - ct * ct, 0.0));
        vec3 rd = vec3(st * cos(ph), ct, st * sin(ph));
        float tg = SkyRaySphere(ro, rd, skyR);
        float t = tg > 0.0 ? tg : max(SkyRaySphere(ro, rd, skyRt), 0.0);
        float dt = t / float(N);
        vec3 thr = vec3(1.0);
        for (int k = 0; k < N; ++k) {
            vec3 p = ro + rd * (dt * (float(k) + 0.5));
            SkyMedium m = SkyMediumAt(length(p) - skyR);
            vec3 Ts = SkyRaySphere(p, sun, skyR) > 0.0 ? vec3(0.0) : SunTransmittance(p, sun);
            vec3 stepT = exp(-m.ext * dt);
            vec3 inv = 1.0 / max(m.ext, vec3(1e-12));
            L += thr * (m.scat * Ts / (4.0 * kPi) - m.scat * Ts / (4.0 * kPi) * stepT) * inv;
            F += thr * (m.scat - m.scat * stepT) * inv;
            thr *= stepT;
        }
        if (tg > 0.0) { // the lit ground under it
            vec3 p = ro + rd * tg, n = normalize(p);
            L += thr * SunTransmittance(p, sun) * (0.3 / kPi) * max(dot(n, sun), 0.0);
        }
    }
    L /= float(D * D);
    F /= float(D * D);
    color = vec4(L / max(vec3(1.0) - F, vec3(1e-3)), 1.0);
}
)";

    // The sky view (Hillaire 2020, 5.3): in-scattered light (sun illuminance 1) for every view direction from one height.
    const char* kFsSkyView = R"(
in vec2 uv;
out vec4 color;
uniform sampler2D transmittanceLut;
uniform sampler2D multiScatterLut;
uniform float viewHeight;
uniform float sunCosZenith;
vec3 SunTransmittance(vec3 p, vec3 sun)
{
    float h = length(p);
    return texture(transmittanceLut, SkyTransmittanceUv(h, dot(p / h, sun))).rgb;
}
void main()
{
    float vh = clamp(viewHeight, skyR + 0.001, skyRt - 0.001);
    float vHorizon = sqrt(max(vh * vh - skyR * skyR, 0.0));
    float beta = acos(clamp(vHorizon / vh, -1.0, 1.0));
    float zenithHorizon = kPi - beta;
    float a;
    if (uv.y < 0.5) { float c = 1.0 - 2.0 * uv.y; a = zenithHorizon * (1.0 - c * c); }
    else { float c = uv.y * 2.0 - 1.0; a = zenithHorizon + beta * c * c; }
    float lightViewCos = 1.0 - 2.0 * uv.x * uv.x;
    float sa = sin(a), lvs = sqrt(max(1.0 - lightViewCos * lightViewCos, 0.0));
    vec3 rd = vec3(sa * lightViewCos, cos(a), sa * lvs);
    vec3 sun = vec3(sqrt(max(1.0 - sunCosZenith * sunCosZenith, 0.0)), sunCosZenith, 0.0);
    vec3 ro = vec3(0.0, vh, 0.0);
    float tg = SkyRaySphere(ro, rd, skyR);
    float t = tg > 0.0 ? tg : max(SkyRaySphere(ro, rd, skyRt), 0.0);
    float c = dot(rd, sun), pr = SkyPhaseRayleigh(c), pm = SkyPhaseMie(c);
    const int N = 32;
    vec3 L = vec3(0.0), thr = vec3(1.0);
    float prev = 0.0;
    for (int k = 0; k < N; ++k) {
        float f = (float(k) + 1.0) / float(N);
        float tk = t * f * f;            // denser near the camera
        float dt = tk - prev, tm = prev + 0.5 * dt;
        prev = tk;
        vec3 p = ro + rd * tm;
        float hp = length(p);
        SkyMedium m = SkyMediumAt(hp - skyR);
        float sunCos = dot(p / hp, sun);
        vec3 Ts = SkyRaySphere(p, sun, skyR) > 0.0 ? vec3(0.0) : SunTransmittance(p, sun);
        vec3 ms = texture(multiScatterLut, vec2(sunCos * 0.5 + 0.5, (hp - skyR) / (skyRt - skyR))).rgb;
        vec3 S = (m.rayleigh * pr + vec3(m.mie * pm)) * Ts + m.scat * ms;
        vec3 stepT = exp(-m.ext * dt);
        L += thr * (S - S * stepT) / max(m.ext, vec3(1e-12));
        thr *= stepT;
    }
    color = vec4(L, 1.0);
}
)";

    // The LUTs are built for an Earth-sized shell (km): see resources/skyphys.glsl.
    constexpr float kLutR = 6360.0f, kLutRt = 6460.0f;

    struct Lut { GLuint tex = 0, fbo = 0; int w = 0, h = 0; };
    struct State {
        bool tried = false, ok = false;
        GLuint vao = 0, progT = 0, progM = 0, progS = 0;
        Lut trans, multi, sky;
        float builtR = -1.0f, builtRt = -1.0f, skyH = -1.0f, skyMu = -2.0f;
        float tint = 0.35f, exposure = 10.0f;
    } S;

    bool EnvOn() { const char* e = getenv("RAYLIBGAME_SKY_PHYS"); return !(e && e[0] == '0'); } // on by default (user, 2026-11-18); =0 off

    GLuint Program(const std::string& common, const char* body, const char* name)
    {
        const std::string fs = "#version 430\n" + common + body;
        const unsigned v = rlCompileShader(kVs, RL_VERTEX_SHADER), f = rlCompileShader(fs.c_str(), RL_FRAGMENT_SHADER);
        const unsigned p = (v && f) ? rlLoadShaderProgram(v, f) : 0;
        if (!p) LOG_WARN("SkyPhys", "%s shader failed\n", name);
        return p;
    }

    void MakeLut(Lut& l, int w, int h)
    {
        l.w = w; l.h = h;
        glGenTextures(1, &l.tex);
        glBindTexture(GL_TEXTURE_2D, l.tex);
        GpuMemTag gm_("skyphys"); // VRAM accounting (#302)
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, w, h, 0, GL_RGBA, GL_HALF_FLOAT, nullptr);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        glGenFramebuffers(1, &l.fbo);
        glBindFramebuffer(GL_FRAMEBUFFER, l.fbo);
        glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, l.tex, 0);
    }

    bool Init()
    {
        S.tried = true;
        if (!EnvOn()) return false;
        if (const char* t = getenv("RAYLIBGAME_SKY_PHYS_TINT")) S.tint = std::clamp((float)atof(t), 0.0f, 1.0f);
        if (const char* e = getenv("RAYLIBGAME_SKY_PHYS_EXPOSURE")) S.exposure = std::max(0.0f, (float)atof(e));
        char* text = LoadFileText("resources/skyphys.glsl");
        if (!text) { LOG_WARN("SkyPhys", "resources/skyphys.glsl missing: physical sky off\n"); return false; }
        const std::string common = text;
        UnloadFileText(text);
        S.progT = Program(common, kFsTransmittance, "transmittance");
        S.progM = Program(common, kFsMultiScatter, "multi-scattering");
        S.progS = Program(common, kFsSkyView, "sky-view");
        if (!S.progT || !S.progM || !S.progS) return false;
        GLint prevFbo = 0;
        glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
        MakeLut(S.trans, 256, 64);
        MakeLut(S.multi, 32, 32);
        MakeLut(S.sky, 192, 108);
        glBindTexture(GL_TEXTURE_2D, 0);
        glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
        glGenVertexArrays(1, &S.vao);
        LOG_INFO("SkyPhys", "physical sky on (tint %.2f, exposure %.1f)\n", S.tint, S.exposure);
        return true;
    }

    void Pass(GLuint prog, const Lut& dst, float R, float Rt)
    {
        glBindFramebuffer(GL_FRAMEBUFFER, dst.fbo);
        glViewport(0, 0, dst.w, dst.h);
        glUseProgram(prog);
        glUniform1f(glGetUniformLocation(prog, "skyR"), R);
        glUniform1f(glGetUniformLocation(prog, "skyRt"), Rt);
        const GLint lt = glGetUniformLocation(prog, "transmittanceLut");
        if (lt >= 0) { glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, S.trans.tex); glUniform1i(lt, 0); }
        const GLint lm = glGetUniformLocation(prog, "multiScatterLut");
        if (lm >= 0) { glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, S.multi.tex); glUniform1i(lm, 1); }
        glDrawArrays(GL_TRIANGLES, 0, 3);
    }
}

bool SkyPhys_On()
{
    if (!S.tried) S.ok = Init();
    return S.ok;
}
float SkyPhys_Tint() { return S.tint; }
float SkyPhys_Exposure() { return S.exposure; }
unsigned SkyPhys_SkyViewTexture() { return S.ok && S.skyH > 0.0f ? S.sky.tex : 0; }
unsigned SkyPhys_TransmittanceTexture() { return S.ok && S.builtR > 0.0f ? S.trans.tex : 0; }

float SkyPhys_LutR() { return kLutR; }
float SkyPhys_LutRt() { return kLutRt; }
float SkyPhys_LutHeight(float planetRadius, float airTop, float camDist)
{
    // The camera's height in the air maps linearly onto the shell's (leaving terra's air at its top leaves the shell's).
    const float f = std::clamp((camDist - planetRadius) / std::max(airTop - planetRadius, 1.0f), 0.0f, 1.0f);
    return kLutR + f * (kLutRt - kLutR);
}

void SkyPhys_Update(float planetRadius, float airTop, float camDist, float sunCosZenith)
{
    if (!SkyPhys_On() || planetRadius <= 0.0f || airTop <= planetRadius) return;
    const float R = kLutR, Rt = kLutRt, viewHeight = SkyPhys_LutHeight(planetRadius, airTop, camDist);
    const bool planet = R != S.builtR || Rt != S.builtRt;
    // The sky view moves with the camera height (1 m, or 0.1% higher up) and the sun (~0.006 deg).
    const bool view = planet || fabsf(viewHeight - S.skyH) > std::max(0.001f, 0.001f * (viewHeight - R))
                   || fabsf(sunCosZenith - S.skyMu) > 1e-4f;
    if (!view) return;
    rlDrawRenderBatchActive();
    GLint prevFbo = 0, vp[4] = { 0, 0, 0, 0 }, prevProg = 0, prevVao = 0;
    glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFbo);
    glGetIntegerv(GL_VIEWPORT, vp);
    glGetIntegerv(GL_CURRENT_PROGRAM, &prevProg);
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &prevVao);
    const GLboolean depth = glIsEnabled(GL_DEPTH_TEST), blend = glIsEnabled(GL_BLEND), cull = glIsEnabled(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST); glDisable(GL_BLEND); glDisable(GL_CULL_FACE);
    glBindVertexArray(S.vao);
    if (planet) {
        Pass(S.progT, S.trans, R, Rt);
        Pass(S.progM, S.multi, R, Rt);
        S.builtR = R; S.builtRt = Rt;
    }
    glUseProgram(S.progS);
    glUniform1f(glGetUniformLocation(S.progS, "viewHeight"), viewHeight);
    glUniform1f(glGetUniformLocation(S.progS, "sunCosZenith"), sunCosZenith);
    Pass(S.progS, S.sky, R, Rt);
    S.skyH = viewHeight; S.skyMu = sunCosZenith;
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, 0);
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, 0);
    glBindVertexArray((GLuint)prevVao);
    glUseProgram((GLuint)prevProg);
    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)prevFbo);
    glViewport(vp[0], vp[1], vp[2], vp[3]);
    if (depth) glEnable(GL_DEPTH_TEST);
    if (blend) glEnable(GL_BLEND);
    if (cull) glEnable(GL_CULL_FACE);
}
