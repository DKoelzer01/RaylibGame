#include "wind.h"
#include "voxel.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>

// The wind field (#650): see wind.h.
namespace {
    constexpr int kTex = 256;           // texels per side
    constexpr float kPeriod = 512.0f;   // metres the texture spans (it tiles)

    bool enabled = true, envRead = false;
    const Planetoid* body = nullptr;
    Vector3 anchor{}, e1{ 1, 0, 0 }, e2{ 0, 0, 1 }; // body-local anchor on the surface near the player, and its tangent basis
    double shiftX = 0.0, shiftY = 0.0;                // how far the air has moved (m, along e1/e2), wrapped to the period
    Vector3 mean{};                                   // body-local mean wind (m/s) with the daily cycle and weather in
    float storm = 0.0f;
    std::vector<unsigned char> texels;                // RGBA8: waves, burst field, direction swing, 255
    Texture2D tex{};

    uint32_t Hash(uint32_t x) { x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16; return x; }
    // Tileable value noise: a lattice of `cells` per period, quintic blend, in [0, 1).
    float Tile(float u, float v, int cells, uint32_t seed) {
        const float x = u * cells, y = v * cells;
        const int x0 = (int)floorf(x), y0 = (int)floorf(y);
        const float fx = x - x0, fy = y - y0;
        auto q = [](float t) { return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f); };
        auto at = [&](int i, int j) { i = ((i % cells) + cells) % cells; j = ((j % cells) + cells) % cells; return (Hash((uint32_t)(i * 7919 + j * 104729) ^ seed) >> 8) * (1.0f / 16777216.0f); };
        const float a = at(x0, y0), b = at(x0 + 1, y0), c = at(x0, y0 + 1), d = at(x0 + 1, y0 + 1);
        const float sx = q(fx), sy = q(fy);
        return (a + (b - a) * sx) + ((c + (d - c) * sx) - (a + (b - a) * sx)) * sy;
    }
    void Build() {
        if (!texels.empty()) return;
        texels.resize((size_t)kTex * kTex * 4);
        for (int j = 0; j < kTex; ++j)
            for (int i = 0; i < kTex; ++i) {
                const float u = (i + 0.5f) / kTex, v = (j + 0.5f) / kTex;
                // Waves: ~64 m fronts (8 per period), ~16 m puffs (32), a little finer still (64).
                const float waves = 0.6f * Tile(u, v, 8, 0x1234u) + 0.3f * Tile(u, v, 32, 0x5678u) + 0.1f * Tile(u, v, 64, 0x9abcu);
                // Burst cells: ~100 m blobs (5 per period); the shader keeps only their peaks.
                const float burst = 0.7f * Tile(u, v, 5, 0x2468u) + 0.3f * Tile(u, v, 16, 0x1357u);
                const float swing = Tile(u, v, 12, 0xfeedu);
                unsigned char* t = &texels[((size_t)j * kTex + i) * 4];
                t[0] = (unsigned char)std::clamp((int)lroundf(waves * 255.0f), 0, 255);
                t[1] = (unsigned char)std::clamp((int)lroundf(burst * 255.0f), 0, 255);
                t[2] = (unsigned char)std::clamp((int)lroundf(swing * 255.0f), 0, 255);
                t[3] = 255;
            }
    }
    // Bilinear, wrapped: the same as the GPU's linear filter with REPEAT (texel centres at (i + 0.5) / kTex).
    void Sample(float u, float v, float out[3]) {
        const float x = u * kTex - 0.5f, y = v * kTex - 0.5f;
        const int x0 = (int)floorf(x), y0 = (int)floorf(y);
        const float fx = x - x0, fy = y - y0;
        auto at = [&](int i, int j, int c) { i = ((i % kTex) + kTex) % kTex; j = ((j % kTex) + kTex) % kTex; return texels[((size_t)j * kTex + i) * 4 + c] / 255.0f; };
        for (int c = 0; c < 3; ++c) {
            const float a = at(x0, y0, c), b = at(x0 + 1, y0, c), cc = at(x0, y0 + 1, c), d = at(x0 + 1, y0 + 1, c);
            out[c] = (a + (b - a) * fx) * (1.0f - fy) + (cc + (d - cc) * fx) * fy;
        }
    }
    float Smooth(float a, float b, float x) { const float t = std::clamp((x - a) / (b - a), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); }
    // The shared formula (grass.vs / trees.vs windAt mirror it): speed factor and direction swing from a texel.
    void Shape(const float s[3], float& factor, float& swing, float& burst) {
        const float waves = s[0] * 2.0f - 1.0f;
        const float thr = 0.80f - 0.14f * storm; // bursts: the burst field's peaks, more of them in a storm
        burst = Smooth(thr, thr + 0.1f, s[1]);
        factor = std::max(0.15f, 1.0f + 0.8f * waves + 1.0f * burst);
        swing = 0.35f * (s[2] * 2.0f - 1.0f);
    }
    bool ReadEnv() {
        if (!envRead) { envRead = true; if (const char* e = getenv("RAYLIBGAME_WIND")) enabled = e[0] != '0'; }
        return enabled;
    }
    void Basis(Vector3 up) {
        const Vector3 ref = fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 };
        e1 = Vector3Normalize(Vector3CrossProduct(ref, up));
        e2 = Vector3CrossProduct(up, e1);
    }
}

bool Wind_Enabled() { return ReadEnv(); }

float Wind_DailyFactor(Vector3 up, Vector3 sun, Vector3 axis) {
    // The sun's angle through the local day: 0 rising in the east, pi/2 noon, pi setting, then night.
    Vector3 east = Vector3CrossProduct(axis, up);
    if (Vector3LengthSqr(east) < 1e-6f) return 1.0f; // at a pole: no day cycle to speak of
    east = Vector3Normalize(east);
    float th = atan2f(Vector3DotProduct(sun, up), Vector3DotProduct(sun, east));
    if (th < 0.0f) th += 2.0f * PI;
    if (th <= PI) // day: calm at dawn, building into the afternoon, easing toward sunset
        return 0.45f + 0.85f * Smooth(0.0f, 0.7f * PI, th) - 0.4f * Smooth(0.75f * PI, PI, th);
    return 0.9f - 0.45f * Smooth(1.0f * PI, 2.0f * PI, th); // night: easing to the dawn calm
}

void Wind_SetFrame(const Planetoid* b, Vector3 playerLocal, Vector3 meanLocal, float s, float daily, float dt) {
    ReadEnv();
    Build();
    const float len = Vector3Length(playerLocal);
    if (!b || len < 1.0f) { body = nullptr; return; }
    const Vector3 up = Vector3Scale(playerLocal, 1.0f / len);
    // A new anchor when the player is 2 km from the old one (or on another body): the gust pattern restarts there.
    if (b != body || Vector3Distance(anchor, playerLocal) > 2000.0f) {
        body = b; anchor = playerLocal; Basis(up);
        shiftX = shiftY = 0.0;
    }
    storm = std::clamp(s, 0.0f, 1.0f);
    mean = Vector3Scale(meanLocal, daily * (1.0f + 0.5f * storm));
    if (dt > 0.0f) {
        shiftX = fmod(shiftX + Vector3DotProduct(mean, e1) * dt, (double)kPeriod);
        shiftY = fmod(shiftY + Vector3DotProduct(mean, e2) * dt, (double)kPeriod);
    }
}

Vector3 Wind_At(const Planetoid& b, Vector3 local, float heightAboveGround, bool withTerrain) {
    if (!ReadEnv() || &b != body || texels.empty()) return { 0, 0, 0 };
    const Vector3 rel = Vector3Subtract(local, anchor);
    const float u = (float)((Vector3DotProduct(rel, e1) - shiftX) / kPeriod), v = (float)((Vector3DotProduct(rel, e2) - shiftY) / kPeriod);
    float s[3]; Sample(u, v, s);
    float factor, swing, burst; Shape(s, factor, swing, burst);
    // Gameplay extras: a log profile over rough ground (z0 = 0.1 m; 10 m = 1) and the terrain's exposure.
    const float h = std::max(heightAboveGround, 0.5f);
    factor *= std::max(0.2f, logf(h / 0.1f) / logf(10.0f / 0.1f));
    if (withTerrain)
        if (const auto* pf = dynamic_cast<const PlanetField*>(b.field.get())) {
            const float r = Vector3Length(local);
            if (r > 1.0f) {
                const Vector3 up = Vector3Scale(local, 1.0f / r);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 })), t2 = Vector3CrossProduct(up, t1);
                const float here = pf->elevation(up);
                float ring = 0.0f;
                for (int k = 0; k < 4; ++k) {
                    const Vector3 off = Vector3Scale((k & 1) ? t2 : t1, (k & 2) ? -300.0f : 300.0f);
                    ring += pf->elevation(Vector3Normalize(Vector3Add(Vector3Scale(up, pf->radius), off)));
                }
                const float rise = here - 0.25f * ring; // > 0 on ridges and knolls, < 0 in hollows and valleys
                factor *= std::clamp(1.0f + rise / 150.0f, 0.5f, 1.3f);
            }
        }
    // The mean direction turned by the swing, scaled.
    const Vector3 m = mean;
    const float ms = Vector3Length(m);
    if (ms < 1e-4f) return { 0, 0, 0 };
    const Vector3 d = Vector3Scale(m, 1.0f / ms);
    const Vector3 up = Vector3Normalize(local);
    const Vector3 side = Vector3CrossProduct(up, d);
    const Vector3 dir = Vector3Add(Vector3Scale(d, cosf(swing)), Vector3Scale(side, sinf(swing)));
    return Vector3Scale(dir, ms * factor);
}

float Wind_Gust01(const Planetoid& b, Vector3 local) {
    if (!ReadEnv() || &b != body || texels.empty()) return 0.0f;
    const Vector3 rel = Vector3Subtract(local, anchor);
    float s[3]; Sample((float)((Vector3DotProduct(rel, e1) - shiftX) / kPeriod), (float)((Vector3DotProduct(rel, e2) - shiftY) / kPeriod), s);
    float factor, swing, burst; Shape(s, factor, swing, burst);
    return burst;
}

void Wind_Bind(Shader s, DVec3 renderOrigin, int slot) {
    const int onLoc = GetShaderLocation(s, "windOn");
    const float on = (ReadEnv() && body && !texels.empty()) ? 1.0f : 0.0f;
    if (onLoc >= 0) rlSetUniform(onLoc, &on, RL_SHADER_UNIFORM_FLOAT, 1);
    if (on == 0.0f) return;
    if (tex.id == 0) {
        Image img{ texels.data(), kTex, kTex, 1, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8 };
        tex = LoadTextureFromImage(img);
        SetTextureWrap(tex, TEXTURE_WRAP_REPEAT);
        SetTextureFilter(tex, TEXTURE_FILTER_BILINEAR);
    }
    rlActiveTextureSlot(slot);
    rlEnableTexture(tex.id);
    rlActiveTextureSlot(0);
    rlSetUniform(GetShaderLocation(s, "windTex"), &slot, RL_SHADER_UNIFORM_SAMPLER2D, 1);
    // The anchor and basis in render space; e1/e2 divided by the body's scale so dot products give body metres.
    const Vector3 o = (body->toWorldD(DVec3(anchor)) - renderOrigin).toFloat();
    const float inv = 1.0f / std::max(body->scale, 1e-6f);
    const Vector3 w1 = Vector3Scale(body->dirToWorld(e1), inv), w2 = Vector3Scale(body->dirToWorld(e2), inv);
    const float shift[2] = { (float)shiftX, (float)shiftY };
    const float m2[2] = { Vector3DotProduct(mean, e1), Vector3DotProduct(mean, e2) };
    const float period = kPeriod;
    rlSetUniform(GetShaderLocation(s, "windOrigin"), &o, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(s, "windE1"), &w1, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(s, "windE2"), &w2, RL_SHADER_UNIFORM_VEC3, 1);
    rlSetUniform(GetShaderLocation(s, "windShift"), shift, RL_SHADER_UNIFORM_VEC2, 1);
    rlSetUniform(GetShaderLocation(s, "windMean"), m2, RL_SHADER_UNIFORM_VEC2, 1);
    rlSetUniform(GetShaderLocation(s, "windStorm"), &storm, RL_SHADER_UNIFORM_FLOAT, 1);
    rlSetUniform(GetShaderLocation(s, "windPeriod"), &period, RL_SHADER_UNIFORM_FLOAT, 1);
}

void Wind_Unload() {
    if (tex.id) { UnloadTexture(tex); tex = Texture2D{}; }
}
