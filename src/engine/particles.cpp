#include "particles.h"
#include "scene.h"
#include "worldtypes.h"
#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace {

struct Particle {
    uint64_t body = 0;
    Vector3 pos{}, vel{};          // body-local
    float age = 0, life = 0, size = 0, groundR = 0, seed = 0;
    Color color = WHITE;
    ParticleKind kind = ParticleKind::Dust;
    bool bounced = false;
    bool alive() const { return age < life; }
};

constexpr size_t kMax = 6000;
std::vector<Particle> pool;
size_t nextSlot = 0;
uint32_t rng = 0x9E3779B9u;
Texture2D sprite{};

float Rand01() { rng ^= rng << 13; rng ^= rng >> 17; rng ^= rng << 5; return (rng & 0xFFFFFF) / 16777215.0f; }
float RandRange(float a, float b) { return a + (b - a) * Rand01(); }
Vector3 RandUnit() {
    const float z = RandRange(-1.0f, 1.0f), a = RandRange(0.0f, 6.2831853f), r = sqrtf(std::max(0.0f, 1.0f - z * z));
    return { r * cosf(a), r * sinf(a), z };
}

// Per kind: gravity (fraction of the body's surface gravity; negative rises), drag (1/s), additive.
struct KindInfo { float gravity, drag; bool additive; };
const KindInfo kInfo[(int)ParticleKind::Count] = {
    { -0.06f, 0.9f, false },  // Smoke
    { -0.12f, 1.4f, true },   // Fire
    { 0.03f, 1.6f, false },   // Dust
    { 1.0f, 0.15f, false },   // Debris
    { 0.6f, 0.4f, true },     // Spark
    { 1.0f, 0.3f, false },    // Splash
    { 0.0f, 1.5f, true },     // Firefly
};

Particle& Slot() {
    if (pool.size() < kMax) { pool.emplace_back(); return pool.back(); }
    // Reuse a dead slot near the cursor, else replace the oldest-placed one.
    for (size_t i = 0; i < 64; ++i) {
        Particle& p = pool[(nextSlot + i) % kMax];
        if (!p.alive()) { nextSlot = (nextSlot + i + 1) % kMax; return p; }
    }
    Particle& p = pool[nextSlot];
    nextSlot = (nextSlot + 1) % kMax;
    return p;
}

} // namespace

void Particles_Emit(const Planetoid& body, Vector3 local, const ParticleBurst& b) {
    const float r = Vector3Length(local);
    const Vector3 up = r > 1e-3f ? Vector3Scale(local, 1.0f / r) : Vector3{ 0, 1, 0 };
    const Vector3 dir = Vector3LengthSqr(b.dir) > 1e-8f ? Vector3Normalize(b.dir) : up;
    for (int i = 0; i < b.count; ++i) {
        Particle& p = Slot();
        p.body = body.id;
        p.kind = b.kind;
        const Vector3 d = Vector3Normalize(Vector3Add(dir, Vector3Scale(RandUnit(), b.spread * 1.5f)));
        p.pos = Vector3Add(local, Vector3Scale(RandUnit(), b.size * 0.5f));
        p.vel = Vector3Add(b.baseVel, Vector3Scale(d, b.speed * RandRange(0.6f, 1.2f)));
        p.age = 0.0f;
        p.life = b.life * RandRange(0.7f, 1.3f);
        p.size = b.size * RandRange(0.75f, 1.25f);
        p.groundR = r - 0.08f;
        p.seed = Rand01();
        p.bounced = false;
        // A little variation in colour, as real smoke and dust have.
        const float k = RandRange(0.88f, 1.08f);
        p.color = { (unsigned char)std::min(255.0f, b.color.r * k), (unsigned char)std::min(255.0f, b.color.g * k), (unsigned char)std::min(255.0f, b.color.b * k), b.color.a };
    }
}

void Particles_EmitWorld(const Planetoid& body, const DVec3& world, const ParticleBurst& b, Vector3 worldDir, Vector3 worldBaseVel) {
    ParticleBurst lb = b;
    lb.dir = Vector3LengthSqr(worldDir) > 1e-8f ? body.dirToLocal(worldDir) : Vector3{ 0, 0, 0 };
    lb.baseVel = body.dirToLocal(worldBaseVel);
    Particles_Emit(body, body.toLocalD(world).toFloat(), lb);
}

void Particles_Update(Scene&, float dt) {
    if (pool.empty() || dt <= 0.0f) return;
    dt = std::min(dt, 0.1f);
    std::unordered_map<uint64_t, float> gravity; // per body: surface gravity (m/s^2), Terra-style g = 9.81 R / 100 km
    for (Particle& p : pool) {
        if (!p.alive()) continue;
        p.age += dt;
        auto it = gravity.find(p.body);
        if (it == gravity.end()) {
            const Planetoid* b = Planetoid::find(p.body);
            it = gravity.emplace(p.body, b ? 9.81f * b->radius * b->scale / 100000.0f : -1.0f).first;
        }
        if (it->second < 0.0f) { p.age = p.life; continue; } // its body is gone
        const KindInfo& k = kInfo[(int)p.kind];
        const float r = Vector3Length(p.pos);
        const Vector3 up = r > 1e-3f ? Vector3Scale(p.pos, 1.0f / r) : Vector3{ 0, 1, 0 };
        Vector3 acc = Vector3Scale(up, -k.gravity * it->second);
        if (p.kind == ParticleKind::Firefly) acc = Vector3Add(acc, Vector3Scale(RandUnit(), 3.0f));
        p.vel = Vector3Scale(Vector3Add(p.vel, Vector3Scale(acc, dt)), expf(-k.drag * dt));
        p.pos = Vector3Add(p.pos, Vector3Scale(p.vel, dt));
        if ((p.kind == ParticleKind::Debris || p.kind == ParticleKind::Splash || p.kind == ParticleKind::Spark) && Vector3Length(p.pos) < p.groundR) {
            if (p.kind == ParticleKind::Debris && !p.bounced) {
                // One bounce: lose most of the speed, keep sliding a little, then lie for the rest of the life.
                const float vr = Vector3DotProduct(p.vel, up);
                p.vel = Vector3Scale(Vector3Subtract(p.vel, Vector3Scale(up, 1.3f * vr)), 0.35f);
                p.pos = Vector3Scale(up, p.groundR);
                p.bounced = true;
            } else if (p.kind == ParticleKind::Debris) {
                p.pos = Vector3Scale(up, p.groundR);
                p.vel = { 0, 0, 0 };
            } else {
                p.age = p.life;
            }
        }
    }
}

void Particles_Draw(Scene& world) {
    if (pool.empty()) return;
    if (sprite.id == 0) {
        Image img = GenImageGradientRadial(32, 32, 0.0f, WHITE, BLANK);
        sprite = LoadTextureFromImage(img);
        UnloadImage(img);
        SetTextureFilter(sprite, TEXTURE_FILTER_BILINEAR);
    }
    const Camera3D& cam = world.camera;
    const Vector3 fwd = Vector3Normalize(Vector3Subtract(cam.target, cam.position));
    const Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, cam.up));
    const Vector3 camUp = Vector3CrossProduct(right, fwd);
    const Vector3 sun = Vector3Normalize(world.sunDirection);

    // Render positions and the light at each body (day side bright, night dim) - per body once.
    struct BodyView { const Planetoid* b; };
    std::unordered_map<uint64_t, BodyView> bodies;
    struct Item { Vector3 at; float d2, size, light; const Particle* p; };
    static std::vector<Item> alpha, glow;
    alpha.clear(); glow.clear();
    for (const Particle& p : pool) {
        if (!p.alive()) continue;
        auto it = bodies.find(p.body);
        if (it == bodies.end()) it = bodies.emplace(p.body, BodyView{ Planetoid::find(p.body) }).first;
        const Planetoid* b = it->second.b;
        if (!b) continue;
        const Vector3 at = (b->toWorldD(DVec3(p.pos)) - world.renderOrigin).toFloat();
        const Vector3 rel = Vector3Subtract(at, cam.position);
        const float d2 = Vector3LengthSqr(rel);
        if (d2 > 1500.0f * 1500.0f || Vector3DotProduct(rel, fwd) < -2.0f) continue; // behind the camera
        const float t = p.age / p.life;
        float size = p.size;
        switch (p.kind) {
        case ParticleKind::Smoke: size *= 1.0f + 3.5f * t; break;
        case ParticleKind::Fire: size *= 1.0f - 0.6f * t; break;
        case ParticleKind::Dust: size *= 1.0f + 1.2f * t; break;
        default: break;
        }
        const float lr = Vector3Length(p.pos);
        const Vector3 upW = b->dirToWorld(lr > 1e-3f ? Vector3Scale(p.pos, 1.0f / lr) : Vector3{ 0, 1, 0 });
        const float light = std::clamp(Vector3DotProduct(upW, sun) * 2.5f + 0.35f, 0.12f, 1.0f);
        (kInfo[(int)p.kind].additive ? glow : alpha).push_back({ at, d2, size, light, &p });
    }
    if (alpha.empty() && glow.empty()) return;

    std::sort(alpha.begin(), alpha.end(), [](const Item& a, const Item& b) { return a.d2 > b.d2; }); // far first

    auto quads = [&](const std::vector<Item>& items) {
        rlSetTexture(sprite.id);
        rlBegin(RL_QUADS);
        for (const Item& it : items) {
            const Particle& p = *it.p;
            const float t = p.age / p.life;
            float a = 1.0f;
            Color c = p.color;
            switch (p.kind) {
            case ParticleKind::Smoke: a = 0.75f * powf(1.0f - t, 1.3f) * std::min(1.0f, t * 8.0f); break;
            case ParticleKind::Dust: a = 0.65f * (1.0f - t) * (1.0f - t); break;
            case ParticleKind::Debris: a = t > 0.8f ? (1.0f - t) * 5.0f : 1.0f; break;
            case ParticleKind::Splash: a = 0.75f * (1.0f - t); break;
            case ParticleKind::Fire: {
                const Color hot = { 255, 230, 140, 255 }, cool = { 140, 30, 10, 255 };
                c = t < 0.4f ? ColorLerp(hot, c, t / 0.4f) : ColorLerp(c, cool, (t - 0.4f) / 0.6f);
                a = (1.0f - t) * 0.6f;
                break;
            }
            case ParticleKind::Spark: a = 1.0f - t; break;
            case ParticleKind::Firefly: { const float s = 0.5f + 0.5f * sinf(p.age * 2.6f + p.seed * 40.0f); a = s * s * std::min(1.0f, t * 4.0f) * std::min(1.0f, (1.0f - t) * 4.0f); break; }
            default: break;
            }
            float k = kInfo[(int)p.kind].additive ? 1.0f : it.light;
            if (!kInfo[(int)p.kind].additive && p.kind != ParticleKind::Debris) k = std::max(k, 0.2f);
            rlColor4ub((unsigned char)(c.r * k), (unsigned char)(c.g * k), (unsigned char)(c.b * k), (unsigned char)std::clamp(a * c.a, 0.0f, 255.0f));
            const float h = it.size * 0.5f;
            const Vector3 r = Vector3Scale(right, h), u = Vector3Scale(camUp, h), q = it.at;
            rlTexCoord2f(0, 0); rlVertex3f(q.x - r.x + u.x, q.y - r.y + u.y, q.z - r.z + u.z);
            rlTexCoord2f(0, 1); rlVertex3f(q.x - r.x - u.x, q.y - r.y - u.y, q.z - r.z - u.z);
            rlTexCoord2f(1, 1); rlVertex3f(q.x + r.x - u.x, q.y + r.y - u.y, q.z + r.z - u.z);
            rlTexCoord2f(1, 0); rlVertex3f(q.x + r.x + u.x, q.y + r.y + u.y, q.z + r.z + u.z);
        }
        rlEnd();
        rlSetTexture(0);
    };
    rlDrawRenderBatchActive();
    rlDisableDepthMask();
    rlDisableBackfaceCulling();
    if (!alpha.empty()) quads(alpha);
    if (!glow.empty()) { rlDrawRenderBatchActive(); BeginBlendMode(BLEND_ADDITIVE); quads(glow); rlDrawRenderBatchActive(); EndBlendMode(); }
    rlDrawRenderBatchActive();
    rlEnableBackfaceCulling();
    rlEnableDepthMask();
}

void Particles_Clear() { pool.clear(); nextSlot = 0; }
int Particles_Count() { int n = 0; for (const Particle& p : pool) n += p.alive() ? 1 : 0; return n; }
void Particles_Unload() { Particles_Clear(); if (sprite.id > 0) UnloadTexture(sprite); sprite = Texture2D{}; }
