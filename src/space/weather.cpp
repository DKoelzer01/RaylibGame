#include "weather.h"
#include "wind.h"
#include "audio.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"
#include "SimplexNoise.h"
#include "menus.h"    // settings.reduceFlashing (#1162)
#include "soundscape.h" // (#1241) thunder

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
    float flash = 0.0f;           // lightning now 0..1 (#1162: lights the sky, Scene::lightning)
    float flashAz = 0.0f, flashEl = 0.3f; // the strike: around the player's horizon and how high (radians)

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
    // Plan 5 B: the biome climate's moisture - deserts seldom cloud over, jungles and swamps often do.
    const float moist = pf->climateAt(dir, std::max(pf->elevation(dir, 64.0f), 0.0f)).moist;
    w.cloud = std::clamp(base * 0.85f + humid * 0.3f - 0.1f + (moist - 0.5f) * 0.45f, 0.0f, 1.0f);
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

void Weather_Update(Scene& world, float dt, float realDt) {
    Planetoid* p = StartBody(world);
    WeatherState target;
    if (p && player.nearest == p) {
        const Vector3 local = p->toLocalOf(player.position, world.renderOrigin);
        if (Vector3Length(local) > 1.0f) target = Weather_At(world, Vector3Normalize(local), WorldTime());
    }
    { static const int envForce = getenv("RAYLIBGAME_WEATHER_FORCE") ? atoi(getenv("RAYLIBGAME_WEATHER_FORCE")) : -1; if (envForce >= 0) forced = envForce; } // (#1083 captures) 0 clear .. 3 storm (WeatherKind)
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
    static const float cloudOverride = getenv("RAYLIBGAME_WEATHER_CLOUD") ? (float)atof(getenv("RAYLIBGAME_WEATHER_CLOUD")) : -1.0f; // (#933 captures) cloud cover, no rain
    if (cloudOverride >= 0.0f) { here.cloud = cloudOverride; here.precipitation = 0.0f; }
    // Above the air or underground the sky is not overhead.
    world.weatherCloud = here.cloud;
    world.weatherFog = here.fog;
    world.cloudDrift = Vector3Add(world.cloudDrift, Vector3Scale(here.wind, -dt * 3.0f)); // clouds move with the wind (a little faster than at the ground)
    // The wind field (#650): this mean wind, storminess (rain a little, storms fully) and the daily cycle at the player.
    if (p && player.nearest == p) {
        const Vector3 local = p->toLocalOf(player.position, world.renderOrigin);
        const float len = Vector3Length(local);
        static float stormEased = 0.0f;
        const float stormTarget = here.kind == WeatherKind::Storm ? 1.0f : here.kind == WeatherKind::Rain ? 0.3f : 0.0f;
        stormEased += (stormTarget - stormEased) * k;
        if (forced >= 0) stormEased = stormTarget;
        const float daily = len > 1.0f ? Wind_DailyFactor(Vector3Scale(local, 1.0f / len), Vector3Normalize(p->dirToLocal(world.sunDirection)), { 0, 1, 0 }) : 1.0f;
        // (#888) at the ground the wind is a fraction of the weather's (that moves the clouds and fronts): 8-12 m/s made every
        // day a fresh breeze with gale bursts, the grass pumping; the surface sees ~2-5 m/s, storms still blow hard (wind.cpp).
        const float kSurfaceWind = Wind_Old888() ? 1.0f : 0.45f; // (RAYLIBGAME_WIND_888=0: the old full wind)
        Wind_SetFrame(p, local, Vector3Scale(here.wind, kSurfaceWind), stormEased, daily, realDt >= 0.0f ? realDt : dt);
    } else Wind_SetFrame(nullptr, {}, {}, 0.0f, 1.0f, realDt >= 0.0f ? realDt : dt);
    // Lightning in a storm: a flash every 6-20 s, thunder a moment later. (#1162) It lights the sky and the clouds toward the
    // strike (atmosphere.fs) and lifts the open-sky ambient briefly (scene.cpp; fragSky keeps caves and interiors dark), with
    // a return stroke's second flicker; it no longer whites out the screen (a 55% overlay, user: "They should light up the
    // sky, not the whole screen"). Reduce flashing (settings): one soft glow, about a third as bright.
    if (here.kind == WeatherKind::Storm && p) { // (also with no nearest body set: the fly camera)
        if (GameClock() > nextFlash) {
            lastFlash = GameClock();
            static const double fixedPeriod = getenv("RAYLIBGAME_LIGHTNING_PERIOD") ? atof(getenv("RAYLIBGAME_LIGHTNING_PERIOD")) : 0.0; // captures: a flash every N s
            nextFlash = GameClock() + (fixedPeriod > 0.0 ? fixedPeriod : 6.0 + 14.0 * (double)GetRandomValue(0, 1000) / 1000.0);
            thunderAt = GameClock() + 0.6 + 1.8 * (double)GetRandomValue(0, 1000) / 1000.0;
            const float az = 6.2831853f * GetRandomValue(0, 1000) / 1000.0f;            // somewhere around the horizon,
            const float el = 0.17f + 0.45f * GetRandomValue(0, 1000) / 1000.0f;          // 10-36 degrees up
            flashAz = az; flashEl = el;
        }
    }
    if (thunderAt > 0.0 && GameClock() > thunderAt) {
        if (Audio_V2()) Soundscape_Thunder(world, (float)(thunderAt - lastFlash)); // (#1241) a rumble from the strike's side
        else { Audio_Play(Sfx::Dig, 1.0f, 0.25f); Audio_Play(Sfx::Build, 0.7f, 0.2f); }
        thunderAt = -1.0;
    }
    {
        const float ts = (float)(GameClock() - lastFlash);
        float e = 0.0f;
        if (lastFlash > 0.0 && ts >= 0.0f && ts < 2.0f) {
            static const bool reduceEnv = getenv("RAYLIBGAME_REDUCE_FLASHING") && getenv("RAYLIBGAME_REDUCE_FLASHING")[0] == '1'; // captures (benches skip settings.txt)
            if (settings.reduceFlashing || reduceEnv) e = 0.35f * std::min(1.0f, ts / 0.15f) * expf(-std::max(0.0f, ts - 0.15f) * 2.5f); // eased in and out
            else {
                e = expf(-ts * 12.0f);                                                    // the stroke
                if (ts > 0.14f) e = std::max(e, 0.75f * expf(-(ts - 0.14f) * 9.0f));     // the return stroke
            }
        }
        // Captures: RAYLIBGAME_LIGHTNING_HOLD=1 holds a flash at full strength (the recorder's frames rarely land on one).
        static const bool hold = getenv("RAYLIBGAME_LIGHTNING_HOLD") && getenv("RAYLIBGAME_LIGHTNING_HOLD")[0] == '1';
        if (hold && here.kind == WeatherKind::Storm) e = 1.0f;
        flash = e;
        world.lightning = e;
        if (p && e > 0.0f) {
            const Vector3 local = p->toLocal(world.camera.position);
            const float d = Vector3Length(local);
            if (d > 1.0f) {
                const Vector3 up = Vector3Scale(local, 1.0f / d);
                const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
                const Vector3 t2 = Vector3CrossProduct(up, t1);
                const Vector3 dl = Vector3Add(Vector3Scale(up, sinf(flashEl)), Vector3Scale(Vector3Add(Vector3Scale(t1, cosf(flashAz)), Vector3Scale(t2, sinf(flashAz))), cosf(flashEl)));
                world.lightningDir = Vector3Normalize(p->dirToWorld(dl));
            }
        }
    }
}

void Weather_Draw3D(Scene& world) {
    if (here.precipitation < 0.05f) return;
    Planetoid* p = StartBody(world);
    if (!p || player.nearest != p) return;
    const auto* pf = dynamic_cast<const PlanetField*>(p->field.get());
    const Vector3 local = p->toLocalOf(player.position, world.renderOrigin);
    const float d = Vector3Length(local);
    if (!pf || d < 1.0f) return;
    const Vector3 up = Vector3Scale(local, 1.0f / d);
    if (d < pf->radius + pf->elevation(up) - 3.0f) return; // underground: no rain in the tunnels
    // Above the cloud layer (and in orbit): none (#154). High peaks keep their snow: the cap follows the ground up.
    if (d > pf->radius + std::max(pf->elevation(up), 0.0f) + world.atmosphere.cloudHeight) return;
    const Vector3 cam = world.camera.position;
    const Vector3 upW = p->dirToWorld(up);
    Vector3 e = Vector3CrossProduct(upW, { 0, 0, 1 });
    if (Vector3Length(e) < 1e-3f) e = Vector3CrossProduct(upW, { 1, 0, 0 });
    e = Vector3Normalize(e);
    const Vector3 n = Vector3CrossProduct(upW, e);
    const bool snow = here.kind == WeatherKind::Snow;
    const int count = (int)(snow ? 1600 * here.precipitation : 1600 * here.precipitation);
    const float fall = snow ? 1.2f : 24.0f, span = 24.0f, height = 16.0f; // rain read slow at 9 m/s on screen (#33)
    const Vector3 windW = p->dirToWorld(here.wind);
    const double t = GameClock();
    // The camera along e and n in the body's frame: drop columns are fixed to the ground (wrapped into the box around the
    // camera), so walking doesn't carry the rain along.
    const DVec3 camL = p->toLocalD(world.renderOrigin + cam);
    const Vector3 eL = p->dirToLocal(e), nL = p->dirToLocal(n);
    const double ce = camL.x * eL.x + camL.y * eL.y + camL.z * eL.z, cn = camL.x * nL.x + camL.y * nL.y + camL.z * nL.z;
    auto mix = [](uint32_t x) { x ^= x >> 16; x *= 0x7feb352du; x ^= x >> 15; x *= 0x846ca68bu; x ^= x >> 16; return x; };
    auto wrap = [&](double v) { const double w = 2.0 * span; return (float)(v - w * floor(v / w)) - span; }; // into [-span, span)
    rlBegin(RL_LINES);
    for (int i = 0; i < count; ++i) {
        // Each drop (no state kept): its own fall speed (+-15%) and start phase, and a new column every time it wraps
        // (#616: every drop kept one column and all wrapped every 16 m / 24 m/s, so the whole pattern repeated each 0.67 s).
        const uint32_t h = mix((uint32_t)i * 2654435761u + 0x9e3779b9u);
        const float phase = (h & 1023) / 1023.0f;
        const double fallI = fall * (0.85 + 0.3 * (((h >> 10) & 1023) / 1023.0));
        const double travel = phase * height + t * fallI;
        const uint32_t cycle = (uint32_t)(int64_t)floor(travel / height);
        const uint32_t c = mix(h ^ (cycle * 0x85ebca6bu));
        const float ax = wrap((c & 0xFFFF) / 65535.0 * 2.0 * span - ce), az = wrap((c >> 16) / 65535.0 * 2.0 * span - cn);
        float y = height - (float)(travel - floor(travel / height) * height) - 4.0f;
        Vector3 pos = Vector3Add(cam, Vector3Add(Vector3Scale(e, ax), Vector3Add(Vector3Scale(n, az), Vector3Scale(upW, y))));
        if (snow) {
            pos = Vector3Add(pos, Vector3Scale(e, 0.4f * sinf((float)t * 1.3f + phase * 20.0f)));
            const float s = 0.07f;
            rlColor4ub(240, 244, 250, 220);
            rlVertex3f(pos.x - e.x * s, pos.y - e.y * s, pos.z - e.z * s); rlVertex3f(pos.x + e.x * s, pos.y + e.y * s, pos.z + e.z * s);
            rlVertex3f(pos.x - n.x * s, pos.y - n.y * s, pos.z - n.z * s); rlVertex3f(pos.x + n.x * s, pos.y + n.y * s, pos.z + n.z * s);
        } else {
            // A motion-blur streak along the drop's velocity (down at `fall`, carried by the wind), ~50 ms long.
            const Vector3 streak = Vector3Add(Vector3Scale(upW, fall * 0.05f), Vector3Scale(windW, -0.05f));
            rlColor4ub(185, 200, 225, 190);
            rlVertex3f(pos.x, pos.y, pos.z);
            rlVertex3f(pos.x + streak.x, pos.y + streak.y, pos.z + streak.z);
        }
    }
    rlEnd();
}

void Weather_DrawHUD() {
    // (#1162) Nothing: lightning lights the sky and the open ground in the 3D scene instead of a white overlay.
}
