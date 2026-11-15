#include "events.h"
#include "audio.h"
#include "buildtool.h"
#include "gameui.h"
#include "particles.h"
#include "player.h"
#include "plots.h"
#include "scene.h"
#include "trees.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

namespace {
    constexpr double kBucket = 600.0, kShowerLength = 120.0;
    constexpr float kShowerOdds = 0.15f;
    long long lastBucket = -1;
    double showerEnd = -1.0, showerStart = -1.0;
    int impactsDone = 0;
    bool force = false;
    int streaksDrawn = 0;

    struct Streak { Vector3 from, dir; float life, len; };
    std::vector<Streak> streaks;
    float nextStreak = 0.0f;
    struct Glow { std::string body; Vector3 local; double at; };    // impact flashes and hot meteorites (real seconds)
    std::vector<Glow> glows;

    uint32_t Hash(uint64_t x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdull; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ull; x ^= x >> 33; return (uint32_t)x; }
    float Rand() { return (float)GetRandomValue(0, 1000000) / 1000000.0f; }

    bool DarkAt(Scene& world, const Planetoid& body, Vector3 local) {
        return Vector3DotProduct(Vector3Normalize(world.sunDirection), body.dirToWorld(Vector3Normalize(local))) < -0.1f;
    }

    bool Impact(Scene& world, Planetoid& body, float dist, float bearing) {
        const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
        if (!pf) return false;
        const Vector3 me = body.toLocal((player.position - world.renderOrigin).toFloat());
        const Vector3 up = Vector3Normalize(me);
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 t2 = Vector3CrossProduct(t1, up);
        const Vector3 d = Vector3Normalize(Vector3Add(up, Vector3Scale(Vector3Add(Vector3Scale(t1, cosf(bearing)), Vector3Scale(t2, sinf(bearing))), dist / pf->radius)));
        const float e = pf->elevation(d);
        if (e < 1.0f) return false;                                          // into the sea
        if (pf->inlandWaterNear(d, 0.0f)) { const float w = pf->waterElevation(d); if (w != PlanetField::kNoWater && w > e - 0.5f) return false; }
        const Vector3 ground = Vector3Scale(d, pf->radius + e);
        // The crater: a carved bowl, its centre 1.2 m under the ground (4.7 m deep at the middle, ~7 m across).
        TerrainEdit c{};
        c.center = Vector3Subtract(ground, Vector3Scale(d, 1.2f));
        c.radius = 3.5f;
        c.carve = true;
        ApplyTerrainEdit(body, c);
        // The meteorite at the bottom, half sunk.
        const Vector3 bottom = Vector3Subtract(ground, Vector3Scale(d, 1.2f + 3.5f - 0.6f));
        if (!body.plots) body.plots = std::make_unique<PlotSet>();
        Plot* plot = body.plots->plotFor(bottom, 0.0f, true);
        if (plot) {
            const Part m{ PartType::Meteorite, plot->cellAt(plot->toPlot(bottom)), (uint8_t)GetRandomValue(0, 23), 0 };
            if (plot->place(m) >= 0) BuildTool_LayFoundation(body, *plot, m);
        }
        if (body.trees) body.trees->addClearing(ground, 6.0f);
        ParticleBurst dust;
        dust.kind = ParticleKind::Dust; dust.color = { 110, 95, 80, 220 }; dust.count = 60;
        dust.size = 1.2f; dust.life = 4.0f; dust.speed = 7.0f; dust.spread = 0.8f;
        Particles_Emit(body, ground, dust);
        ParticleBurst sparks;
        sparks.kind = ParticleKind::Debris; sparks.color = { 255, 170, 80, 255 }; sparks.count = 40;
        sparks.size = 0.12f; sparks.life = 2.0f; sparks.speed = 12.0f; sparks.spread = 0.7f;
        Particles_Emit(body, ground, sparks);
        glows.push_back({ body.name, Vector3Add(bottom, Vector3Scale(d, 0.8f)), GetTime() });
        const float near = Vector3Distance(ground, me);
        Audio_Play(Sfx::Dig, std::clamp(1.4f - near / 250.0f, 0.3f, 1.0f), 0.35f);
        GameUI_Toast(TextFormat("A meteor came down %.0f m away", near));
        // A last streak, down onto it.
        const Vector3 g = body.toWorld(ground);
        const Vector3 dirDown = Vector3Normalize(Vector3Subtract(g, Vector3Add(g, Vector3Add(Vector3Scale(body.dirToWorld(d), 500.0f), Vector3Scale(body.dirToWorld(t1), 300.0f)))));
        streaks.push_back({ Vector3Subtract(g, Vector3Scale(dirDown, 520.0f)), dirDown, 0.6f, 160.0f });
        return true;
    }
}

void Events_Reset() { lastBucket = -1; showerEnd = showerStart = -1.0; impactsDone = 0; streaks.clear(); glows.clear(); force = false; }
bool Events_ShowerOn() { return WorldTime() < showerEnd; }
void Events_ForceShower() { force = true; }
int Events_Streaks() { return streaksDrawn; }

void Events_Update(Scene& world, float dt) {
    Planetoid* body = player.nearest;
    const auto* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
    if (!pf || player.flyMode) return;
    const Vector3 me = body->toLocal((player.position - world.renderOrigin).toFloat());
    if (Vector3Length(me) - pf->radius > 400.0f) return;                    // on the ground (or near it) only
    const double now = WorldTime();
    const long long bucket = (long long)floor(now / kBucket);
    if (bucket != lastBucket || force) {
        const bool fresh = lastBucket >= 0 || force;                        // not on the first look after loading
        lastBucket = bucket;
        const uint32_t h = Hash((uint64_t)bucket * 0x9e3779b97f4a7c15ull ^ (uint64_t)pf->foliageSeed ^ 0x3e7e02u);
        if (force || (fresh && DarkAt(world, *body, me) && (h & 0xFFFF) / 65535.0f < kShowerOdds)) {
            showerStart = now; showerEnd = now + kShowerLength; impactsDone = 0;
            GameUI_Toast("A meteor shower lights up the sky");
        }
        force = false;
    }
    if (now < showerEnd) {
        const double f = (now - showerStart) / kShowerLength;
        if ((impactsDone == 0 && f > 0.33) || (impactsDone == 1 && f > 0.72)) {
            ++impactsDone;
            for (int tries = 0; tries < 6; ++tries) if (Impact(world, *body, 80.0f + 100.0f * Rand(), Rand() * 2.0f * PI)) break;
        }
        if ((nextStreak -= dt) <= 0.0f) { // a streak high over the player, crossing the sky
            nextStreak = 0.12f + 0.35f * Rand();
            const Vector3 up = body->dirToWorld(Vector3Normalize(me));
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const Vector3 t2 = Vector3CrossProduct(t1, up);
            const float a = Rand() * 2.0f * PI, b = Rand() * 2.0f * PI;
            const Vector3 from = Vector3Add(world.camera.position, Vector3Add(Vector3Scale(up, 500.0f + 300.0f * Rand()), Vector3Add(Vector3Scale(t1, 900.0f * cosf(a)), Vector3Scale(t2, 900.0f * sinf(a)))));
            const Vector3 dir = Vector3Normalize(Vector3Add(Vector3Scale(up, -0.35f), Vector3Add(Vector3Scale(t1, cosf(b)), Vector3Scale(t2, sinf(b)))));
            streaks.push_back({ from, dir, 0.5f + 0.4f * Rand(), 90.0f + 120.0f * Rand() });
        }
    }
    for (Streak& s : streaks) { s.from = Vector3Add(s.from, Vector3Scale(s.dir, 900.0f * dt)); s.life -= dt; }
    streaks.erase(std::remove_if(streaks.begin(), streaks.end(), [](const Streak& s) { return s.life <= 0.0f; }), streaks.end());
    glows.erase(std::remove_if(glows.begin(), glows.end(), [](const Glow& g) { return GetTime() - g.at > 90.0; }), glows.end());
}

void Events_Draw3D(Scene& world) {
    (void)world;
    if (streaks.empty()) return;
    rlDrawRenderBatchActive();
    for (const Streak& s : streaks) {
        const Vector3 tail = Vector3Subtract(s.from, Vector3Scale(s.dir, s.len));
        const unsigned char a = (unsigned char)std::clamp(s.life * 400.0f, 0.0f, 255.0f);
        DrawLine3D(tail, s.from, Color{ 255, 240, 200, a });
        DrawLine3D(Vector3Add(tail, { 0.0f, 0.6f, 0.0f }), Vector3Add(s.from, { 0.0f, 0.6f, 0.0f }), Color{ 255, 200, 140, (unsigned char)(a / 2) });
        ++streaksDrawn;
    }
    rlDrawRenderBatchActive();
}

void Events_AddLights(Scene& world) {
    Planetoid* body = player.nearest;
    if (!body) return;
    for (const Glow& g : glows) {
        if (g.body != body->name) continue;
        const float t = (float)(GetTime() - g.at);
        LocalLight L;
        L.position = body->toWorld(g.local);
        if (t < 1.5f) { L.range = 40.0f; L.color = { 1.0f, 0.85f, 0.6f }; L.intensity = 4.0f * (1.0f - t / 1.5f) + 1.0f; } // the flash
        else { L.range = 9.0f; L.color = { 1.0f, 0.45f, 0.15f }; L.intensity = 1.6f * (1.0f - t / 90.0f); }               // still hot
        world.localLights.push_back(L);
    }
}

void Events_Save(FILE* f) { fprintf(f, "events_meteor %lld\n", lastBucket); }
bool Events_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "events_meteor") != 0) return false;
    long long b = -1;
    if (fscanf(f, " %lld", &b) == 1) lastBucket = b;
    return true;
}

bool Events_TestImpact(Scene& world, Planetoid& body, float dist, float bearing) { return Impact(world, body, dist, bearing); }
