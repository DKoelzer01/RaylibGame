#include "weather.h"
#include "audio.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include "SimplexNoise.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <cmath>
#include <random>

namespace {
    // The cloud field: a 3D noise over the sphere (features ~6 km) that drifts with the wind.
    struct Fields { SimplexNoise cloud, humid, wind; bool ready = false; };
    Fields fields;
    void Init(unsigned int seed) {
        if (fields.ready) return;
        std::mt19937 rng(seed ^ 0x5ea7u);
        std::uniform_real_distribution<float> unit(0.0f, 1.0f);
        auto make = [&](float f) { return SimplexNoise(f, 1.0f, 2.0f, 0.5f, unit(rng) * 10000.0f, unit(rng) * 10000.0f, unit(rng) * 10000.0f); };
        fields.cloud = make(1.0f / 6000.0f);
        fields.humid = make(1.0f / 9000.0f);
        fields.wind = make(1.0f / 30000.0f);
        fields.ready = true;
    }

    WeatherState here;            // smoothed at the player
    int forced = -1;
    double lastFlash = 0.0, nextFlash = 0.0, thunderAt = -1.0;
    float flash = 0.0f;           // lightning overlay 0..1

    Planetoid* StartBody(Scene& world) { return world.objects.empty() ? nullptr : dynamic_cast<Planetoid*>(world.objects[0].get()); }

    // The air temperature at a direction (survival's formula, for any place on the planet).
    float TempAt(Scene& world, const Planetoid& p, const PlanetField& pf, Vector3 dir) {
        const float lat = fabsf(dir.y);
        const float h = std::max(0.0f, pf.elevation(dir, 64.0f));
        const float sun = Vector3DotProduct(Vector3Normalize(world.sunDirection), p.dirToWorld(dir));
        return 32.0f - 57.0f * powf(lat, 1.5f) - 6.5f * h / 1000.0f + 7.0f * sun;
    }
}

const char* Weather_Name(WeatherKind k) {
    switch (k) {
        case WeatherKind::Clear: return "clear";
        case WeatherKind::Cloudy: return "overcast";
        case WeatherKind::Rain: return "rain";
        case WeatherKind::Storm: return "storm";
        case WeatherKind::Snow: return "snow";
        case WeatherKind::Fog: return "fog";
    }
    return "";
}

WeatherState Weather_At(Scene& world, Vector3 dir, double time) {
    WeatherState w;
    Planetoid* p = StartBody(world);
    const auto* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
    if (!p || !pf || pf->moon || !world.atmosphere.enabled) return w;
    Init(p->seed);
    const float R = pf->radius;
    // The wind: a slowly turning direction per region, 4-12 m/s; the cloud field moves with it.
    const float wa = fields.wind.fractal(2, dir.x * R, dir.y * R, dir.z * R + (float)(time * 0.02)) * PI * 2.0f;
    const float ws = 8.0f + 4.0f * fields.wind.fractal(1, dir.z * R, dir.x * R, dir.y * R);
    Vector3 east = Vector3CrossProduct({ 0, 1, 0 }, dir);
    east = Vector3Length(east) < 1e-3f ? Vector3{ 1, 0, 0 } : Vector3Normalize(east);
    const Vector3 north = Vector3CrossProduct(dir, east);
    w.wind = Vector3Add(Vector3Scale(east, cosf(wa) * ws), Vector3Scale(north, sinf(wa) * ws));
    const Vector3 drift = Vector3Scale(w.wind, (float)fmod(time, 1.0e6)); // metres the air has moved
    const Vector3 s = Vector3Subtract(Vector3Scale(dir, R), drift);
    // Cloud: wetter near the equator and the coasts' latitudes, with slow daily swings.
    const float base = 0.5f + 0.5f * fields.cloud.fractal(4, s.x, s.y, s.z);
    const float humid = 0.5f + 0.5f * fields.humid.fractal(2, dir.x * R, dir.y * R, dir.z * R + (float)(time * 0.05));
    w.cloud = std::clamp(base * 0.85f + humid * 0.3f - 0.1f, 0.0f, 1.0f);
    const float T = TempAt(world, *p, *pf, dir);
    if (w.cloud > 0.62f) {
        w.precipitation = std::clamp((w.cloud - 0.62f) / 0.25f, 0.0f, 1.0f);
        w.kind = T < 1.0f ? WeatherKind::Snow : (w.cloud > 0.78f && T > 12.0f ? WeatherKind::Storm : WeatherKind::Rain);
    } else if (w.cloud > 0.45f) {
        w.kind = WeatherKind::Cloudy;
    }
    // Fog: around dawn (the sun just below to a little above the horizon), humid and calm, low ground.
    const float sunUp = Vector3DotProduct(Vector3Normalize(world.sunDirection), p->dirToWorld(dir));
    if (w.kind != WeatherKind::Storm && sunUp > -0.15f && sunUp < 0.2f && humid > 0.62f && ws < 9.5f && pf->elevation(dir, 64.0f) < 300.0f) {
        w.fog = std::clamp((humid - 0.62f) / 0.15f, 0.0f, 1.0f) * (1.0f - fabsf(sunUp - 0.02f) / 0.2f);
        if (w.fog > 0.3f && w.kind != WeatherKind::Rain && w.kind != WeatherKind::Snow) w.kind = WeatherKind::Fog;
    }
    return w;
}

void Weather_Force(int kind) { forced = kind; }
const WeatherState& Weather_Here() { return here; }
bool Weather_StormHere() { return here.kind == WeatherKind::Storm; }

bool Weather_RainingAt(Scene& world, Vector3 dir) {
    if (forced >= 0) return forced == (int)WeatherKind::Rain || forced == (int)WeatherKind::Storm;
    const WeatherState w = Weather_At(world, dir, WorldTime());
    return w.kind == WeatherKind::Rain || w.kind == WeatherKind::Storm;
}

void Weather_Update(Scene& world, float dt) {
    Planetoid* p = StartBody(world);
    WeatherState target;
    if (p && player.nearest == p) {
        const Vector3 local = p->toLocal((player.position - world.renderOrigin).toFloat());
        if (Vector3Length(local) > 1.0f) target = Weather_At(world, Vector3Normalize(local), WorldTime());
    }
    if (forced >= 0) {
        target.kind = (WeatherKind)forced;
        target.cloud = forced == 0 ? 0.1f : forced == 1 ? 0.55f : forced == 5 ? 0.4f : forced == 3 ? 0.95f : 0.75f;
        target.precipitation = forced >= 2 && forced <= 4 ? 0.8f : 0.0f;
        target.fog = forced == 5 ? 0.9f : 0.0f;
    }
    // Ease toward it (a front takes a minute to come over).
    const float k = std::min(1.0f, dt / 20.0f);
    here.cloud += (target.cloud - here.cloud) * k;
    here.precipitation += (target.precipitation - here.precipitation) * k;
    here.fog += (target.fog - here.fog) * k;
    here.wind = Vector3Lerp(here.wind, target.wind, k);
    here.kind = target.kind;
    if (forced >= 0) { here.cloud = target.cloud; here.precipitation = target.precipitation; here.fog = target.fog; }
    // Above the air or underground the sky is not overhead.
    world.weatherCloud = here.cloud;
    world.weatherFog = here.fog;
    // Lightning in a storm: a flash every 6-20 s, thunder a moment later.
    if (here.kind == WeatherKind::Storm) {
        if (GetTime() > nextFlash) {
            flash = 1.0f;
            lastFlash = GetTime();
            nextFlash = GetTime() + 6.0 + 14.0 * (double)GetRandomValue(0, 1000) / 1000.0;
            thunderAt = GetTime() + 0.6 + 1.8 * (double)GetRandomValue(0, 1000) / 1000.0;
        }
    }
    if (thunderAt > 0.0 && GetTime() > thunderAt) { Audio_Play(Sfx::Dig, 1.0f, 0.25f); Audio_Play(Sfx::Build, 0.7f, 0.2f); thunderAt = -1.0; }
    flash = std::max(0.0f, flash - dt * 3.0f);
}

void Weather_Draw3D(Scene& world) {
    if (here.precipitation < 0.05f) return;
    Planetoid* p = StartBody(world);
    if (!p || player.nearest != p) return;
    const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
    const Vector3 local = p->toLocal((player.position - world.renderOrigin).toFloat());
    const float d = Vector3Length(local);
    if (!pf || d < 1.0f) return;
    const Vector3 up = Vector3Scale(local, 1.0f / d);
    if (d < pf->radius + pf->elevation(up) - 3.0f) return; // underground: no rain in the tunnels
    const Vector3 cam = world.camera.position;
    const Vector3 upW = p->dirToWorld(up);
    Vector3 e = Vector3CrossProduct(upW, { 0, 0, 1 });
    if (Vector3Length(e) < 1e-3f) e = Vector3CrossProduct(upW, { 1, 0, 0 });
    e = Vector3Normalize(e);
    const Vector3 n = Vector3CrossProduct(upW, e);
    const bool snow = here.kind == WeatherKind::Snow;
    const int count = (int)(snow ? 1600 * here.precipitation : 1600 * here.precipitation);
    const float fall = snow ? 1.2f : 9.0f, span = 24.0f, height = 16.0f;
    const Vector3 windW = p->dirToWorld(here.wind);
    const double t = GetTime();
    rlBegin(RL_LINES);
    for (int i = 0; i < count; ++i) {
        // Each drop: a fixed column around the camera, falling and wrapping (no state kept).
        const uint32_t h = (uint32_t)i * 2654435761u;
        const float ax = ((h & 1023) / 1023.0f - 0.5f) * 2.0f * span, az = (((h >> 10) & 1023) / 1023.0f - 0.5f) * 2.0f * span;
        const float phase = ((h >> 20) & 1023) / 1023.0f;
        float y = height - (float)fmod(phase * height + t * fall, (double)height) - 4.0f;
        Vector3 pos = Vector3Add(cam, Vector3Add(Vector3Scale(e, ax), Vector3Add(Vector3Scale(n, az), Vector3Scale(upW, y))));
        if (snow) {
            pos = Vector3Add(pos, Vector3Scale(e, 0.4f * sinf((float)t * 1.3f + phase * 20.0f)));
            const float s = 0.07f;
            rlColor4ub(240, 244, 250, 220);
            rlVertex3f(pos.x - e.x * s, pos.y - e.y * s, pos.z - e.z * s); rlVertex3f(pos.x + e.x * s, pos.y + e.y * s, pos.z + e.z * s);
            rlVertex3f(pos.x - n.x * s, pos.y - n.y * s, pos.z - n.z * s); rlVertex3f(pos.x + n.x * s, pos.y + n.y * s, pos.z + n.z * s);
        } else {
            const Vector3 streak = Vector3Add(Vector3Scale(upW, 0.7f), Vector3Scale(windW, 0.04f));
            rlColor4ub(185, 200, 225, 190);
            rlVertex3f(pos.x, pos.y, pos.z);
            rlVertex3f(pos.x + streak.x, pos.y + streak.y, pos.z + streak.z);
        }
    }
    rlEnd();
}

void Weather_DrawHUD() {
    if (flash > 0.0f) DrawRectangle(0, 0, GetScreenWidth(), GetScreenHeight(), Fade(WHITE, 0.55f * flash));
}
