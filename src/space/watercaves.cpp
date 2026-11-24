#include "watercaves.h"

#include "audio.h"
#include "fluid.h"
#include "player.h"
#include "profiler.h"
#include "rivers.h"
#include "scene.h"
#include "soundscape.h"
#include "voxel.h"
#include "weather.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <random>

float WaterCaves_RiverLoudness(float halfWidth, float depth, float slope, float flowShare) {
    // The discharge's size (width x depth) and how fast it runs (its slope): a 2 m brook murmurs, a 20 m river on a steep
    // reach roars. Square roots keep big rivers from drowning the rest.
    const float size = std::clamp(sqrtf(std::max(0.0f, halfWidth * 2.0f * depth)) / 8.0f, 0.0f, 1.0f);
    const float rush = std::clamp(0.4f + sqrtf(std::max(0.0f, slope)) * 6.0f, 0.4f, 1.0f);
    return std::clamp(flowShare, 0.0f, 1.0f) * std::clamp(0.25f + 0.75f * size, 0.0f, 1.0f) * rush;
}

float WaterCaves_FlowLoudness(float discharge) {
    if (discharge < 0.02f) return 0.0f;                         // still water is silent
    return std::clamp(sqrtf(discharge) * 0.6f, 0.0f, 1.0f);    // ~2.8 m^3/s: full
}

float WaterCaves_DripInterval(float depth, float random01) {
    const float wetness = std::clamp((depth - 3.0f) / 30.0f, 0.0f, 1.0f);   // deeper caves are wetter, up to 33 m
    const float mean = 4.0f - 2.8f * wetness;                              // 4 s near a mouth .. 1.2 s deep down
    return mean * (0.4f + 1.2f * std::clamp(random01, 0.0f, 1.0f));
}

namespace {
    constexpr int kFlowSlots = 7;      // the listener's 8 m cell and its six neighbours
    int riverEm[3] = { -1, -1, -1 };   // along the river's course: upstream, nearest, downstream
    int flowEm[kFlowSlots] = { -1, -1, -1, -1, -1, -1, -1 };
    int mouthEm = -1, rumbleEm = -1;
    double nextScan = 0.0, nextCave = 0.0, nextDrip = 0.0;
    std::mt19937 rng(99);
    std::uniform_real_distribution<float> u01(0.0f, 1.0f);

    int Sound(const char* name) { const int s = Audio_Sound(name); return s; }

    int Emitter(int& slot, const char* sound, AudioBus bus, const AudioAttenuation& att, bool positional = true) {
        if (slot >= 0) return slot;
        const int snd = Sound(sound);
        if (snd >= 0) slot = Audio_EmitterCreate(snd, bus, AudioPriority::Low, att, positional);
        return slot;
    }

    void Rivers(Planetoid& p, const PlanetField& pf, const DVec3& local) {
        const RiverNetwork* net = pf.rivers();
        float loud = 0.0f;
        Vector3 centre{}, down{};
        float halfWidth = 0.0f;
        const float r = (float)local.length();
        const Vector3 lf = local.toFloat();
        const Vector3 dir = Vector3Scale(lf, 1.0f / r);
        if (net) {
            const RiverNetwork::Hit h = net->nearest(dir, 150.0f);
            if (h.found && h.piece >= 0 && h.depth > 0.05f && !Fluid::Simulated(p, Vector3Scale(dir, pf.radius + h.level - 0.3f))) {
                loud = WaterCaves_RiverLoudness(h.halfWidth, h.depth, h.slope, net->flowOf(h.piece));
                // The centreline point: across the river from here, on whichever side its distance shrinks.
                down = Vector3Normalize(Vector3Subtract(h.down, Vector3Scale(dir, Vector3DotProduct(h.down, dir))));
                const Vector3 across = Vector3CrossProduct(dir, down);
                const RiverNetwork::Hit a = net->nearest(Vector3Normalize(Vector3Add(lf, Vector3Scale(across, 3.0f))), 150.0f);
                const float side = a.found && a.dist < h.dist ? 1.0f : -1.0f;
                centre = Vector3Add(Vector3Scale(dir, pf.radius + h.level), Vector3Scale(across, side * h.dist));
                halfWidth = h.halfWidth;
            }
        }
        AudioAttenuation att;
        att.refDistance = std::max(3.0f, halfWidth);
        att.maxDistance = 150.0f;
        att.lowpassDistance = 90.0f;
        for (int k = 0; k < 3; ++k) {
            if (Emitter(riverEm[k], "water", AudioBus::Ambience, att) < 0) return;
            if (loud > 0.0f)
                Audio_EmitterSetPosition(riverEm[k], p.toWorldD(DVec3(Vector3Add(centre, Vector3Scale(down, 30.0f * (k - 1))))));
            Audio_EmitterSetVolume(riverEm[k], loud * (k == 1 ? 1.0f : 0.7f));
            Audio_EmitterSetPitch(riverEm[k], 1.0f - 0.3f * std::clamp(halfWidth / 20.0f, 0.0f, 1.0f) + 0.04f * (k - 1)); // big rivers deeper
            Audio_EmitterSetPlaying(riverEm[k], loud > 0.01f);
        }
    }

    void SimulatedFlow(Planetoid& p, const DVec3& local) {
        // Only where the simulation holds water nearby (a cheap look first).
        const Vector3 lf = local.toFloat();
        const bool any = Fluid::AnyWaterNear(p, lf, 16.0f, 0.2f); // (WaterAt at 16 m read 2 M cells a site, #1736)
        AudioAttenuation att;
        att.refDistance = 2.0f;
        att.maxDistance = 40.0f;
        att.lowpassDistance = 30.0f;
        const float cell = 8.0f;
        const Vector3 c0 = { (floorf(lf.x / cell) + 0.5f) * cell, (floorf(lf.y / cell) + 0.5f) * cell, (floorf(lf.z / cell) + 0.5f) * cell };
        static const Vector3 kOff[kFlowSlots] = { { 0, 0, 0 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
        for (int k = 0; k < kFlowSlots; ++k) {
            float loud = 0.0f;
            const Vector3 at = Vector3Add(c0, Vector3Scale(kOff[k], cell));
            if (any) loud = WaterCaves_FlowLoudness(Vector3Length(Fluid::FlowAt(p, at, 5.0f))); // (keeps the spot watched)
            if (flowEm[k] < 0 && loud <= 0.0f) continue;
            if (Emitter(flowEm[k], "water", AudioBus::Ambience, att) < 0) return;
            Audio_EmitterSetPosition(flowEm[k], p.toWorldD(DVec3(at)));
            Audio_EmitterSetVolume(flowEm[k], loud);
            Audio_EmitterSetPitch(flowEm[k], 1.1f + 0.03f * k);
            Audio_EmitterSetPlaying(flowEm[k], loud > 0.01f);
        }
    }

    void Caves(Scene& world, Planetoid& p, const PlanetField& pf, const DVec3& local, float depth, float air) {
        const Vector3 lf = local.toFloat();
        const float r = Vector3Length(lf);
        const Vector3 dir = Vector3Scale(lf, 1.0f / r);
        const double now = GetTime();
        if (now >= nextCave) {
            nextCave = now + 1.0;
            // Wind at the nearest cave mouth within 80 m: a hollow, higher wind, stronger with the weather's.
            Vector3 at{}, out{};
            const bool mouth = air > 0.1f && pf.nearestCaveMouth(dir, 80.0f, at, out);
            AudioAttenuation att;
            att.refDistance = 4.0f;
            att.maxDistance = 70.0f;
            att.lowpassDistance = 50.0f;
            if (Emitter(mouthEm, "wind", AudioBus::Ambience, att) >= 0) {
                const float wind = Vector3Length(Weather_Here().wind);
                const float v = mouth ? air * (0.25f + 0.5f * std::clamp(wind / 10.0f, 0.0f, 1.0f)) : 0.0f;
                if (mouth) Audio_EmitterSetPosition(mouthEm, p.toWorldD(DVec3(at)));
                Audio_EmitterSetVolume(mouthEm, v);
                Audio_EmitterSetPitch(mouthEm, 1.35f);
                Audio_EmitterSetPlaying(mouthEm, v > 0.01f);
            }
            // Under a volcano, below ground: the deep rumble (lava tubes, #1063).
            const float under = std::clamp((depth - 3.0f) / 12.0f, 0.0f, 1.0f);
            const float volcano = under > 0.0f ? pf.volcanoWeight(Vector3Scale(dir, pf.radius)) : 0.0f;
            if (Emitter(rumbleEm, "bed_rumble", AudioBus::Ambience, {}, false) >= 0) {
                Audio_EmitterSetVolume(rumbleEm, 0.8f * under * volcano);
                Audio_EmitterSetPlaying(rumbleEm, under * volcano > 0.01f);
            }
        }
        // Drips from the ceiling, underground and in air.
        if (depth > 3.0f && air > 0.1f && now >= nextDrip) {
            nextDrip = now + WaterCaves_DripInterval(depth, u01(rng));
            static int drip = -2;
            if (drip == -2) drip = Sound("drip");
            if (drip < 0) return;
            // A spot 1.5-7 m around the listener, the ceiling straight above it (render space).
            const Vector3 upW = p.dirToWorld(dir);
            const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(upW, fabsf(upW.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
            const Vector3 t2 = Vector3CrossProduct(upW, t1);
            const float a = 6.2831853f * u01(rng), rr = 1.5f + 5.5f * u01(rng);
            const Vector3 from = Vector3Add(world.camera.position, Vector3Add(Vector3Scale(t1, cosf(a) * rr), Vector3Scale(t2, sinf(a) * rr)));
            TerrainHit hit;
            if (!RaycastTerrain(world, Ray{ from, upW }, 12.0f, hit)) return;
            AudioAttenuation att;
            att.refDistance = 1.5f;
            att.maxDistance = 30.0f;
            att.lowpassDistance = 25.0f;
            Audio_PlaySoundAt(drip, world.renderOrigin + hit.point, AudioBus::Ambience, AudioPriority::Low, 0.35f + 0.3f * u01(rng),
                              0.8f + 0.5f * u01(rng), att);
        }
    }
}

void WaterCaves_Reset() {
    for (int& e : riverEm) { Audio_EmitterDestroy(e); e = -1; }
    for (int& e : flowEm) { Audio_EmitterDestroy(e); e = -1; }
    Audio_EmitterDestroy(mouthEm); mouthEm = -1;
    Audio_EmitterDestroy(rumbleEm); rumbleEm = -1;
}

void WaterCaves_Update(Scene& world, float dt) {
    (void)dt;
    if (!Audio_V2()) return;
    Planetoid* p = player.nearest;
    const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
    if (!pf) return;
    PROFILE_SCOPE("audio:water and caves");
    const DVec3 local = p->toLocalD(Audio_GetListener().pos);
    if (local.length() < 1.0) return;
    const double now = GetTime();
    if (now >= nextScan) { // rivers and simulated flow twice a second
        nextScan = now + 0.5;
        { PROFILE_SCOPE("audio:rivers"); Rivers(*p, *pf, local); }
        { PROFILE_SCOPE("audio:flow"); SimulatedFlow(*p, local); }
    }
    { PROFILE_SCOPE("audio:caves"); Caves(world, *p, *pf, local, -Soundscape_HeightAboveGround(), Soundscape_Air()); }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1242): river and flow loudness, drip timing ----
#include "unit_helpers.h"

TEST(WaterCaves, RiversByTheirFlowAndADivertedRiverIsQuiet) {
    const float brook = WaterCaves_RiverLoudness(1.0f, 0.4f, 0.002f, 1.0f);
    const float river = WaterCaves_RiverLoudness(10.0f, 3.0f, 0.002f, 1.0f);
    const float rapids = WaterCaves_RiverLoudness(10.0f, 3.0f, 0.03f, 1.0f);
    EXPECT_GT(brook, 0.0f);
    EXPECT_GT(river, brook);
    EXPECT_GT(rapids, river);
    EXPECT_LE(rapids, 1.0f);
    EXPECT_FLOAT_EQ(WaterCaves_RiverLoudness(10.0f, 3.0f, 0.03f, 0.0f), 0.0f);           // dammed dry
    EXPECT_NEAR(WaterCaves_RiverLoudness(10.0f, 3.0f, 0.03f, 0.25f), 0.25f * rapids, 1e-6f); // a quarter of its water
}

TEST(WaterCaves, FlowLoudnessAndDrips) {
    EXPECT_FLOAT_EQ(WaterCaves_FlowLoudness(0.0f), 0.0f);
    EXPECT_FLOAT_EQ(WaterCaves_FlowLoudness(0.01f), 0.0f);
    EXPECT_GT(WaterCaves_FlowLoudness(0.5f), WaterCaves_FlowLoudness(0.1f));
    EXPECT_FLOAT_EQ(WaterCaves_FlowLoudness(100.0f), 1.0f);
    EXPECT_LT(WaterCaves_DripInterval(40.0f, 0.5f), WaterCaves_DripInterval(4.0f, 0.5f)); // deeper drips more
    for (float r : { 0.0f, 0.5f, 1.0f }) {
        EXPECT_GT(WaterCaves_DripInterval(10.0f, r), 0.3f);
        EXPECT_LT(WaterCaves_DripInterval(10.0f, r), 7.0f);
    }
}
#endif
