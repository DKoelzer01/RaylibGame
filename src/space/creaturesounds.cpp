#include "creaturesounds.h"

#include "audio.h"
#include "music.h"
#include "worldtypes.h"

#include <algorithm>
#include <cmath>
#include <functional>

const char* CreatureSound_Voice(const std::string& f) {
    for (const char* v : { "beast", "lurker", "scavenger", "ooze", "undead", "construct" })
        if (f == v) return v;
    return "beast";
}

std::string CreatureSound_Fallback(CreatureEvent e, const std::string& faction) {
    return std::string("creature_") + CreatureEvent_Name(e) + "_" + CreatureSound_Voice(faction);
}

AudioPriority CreatureSound_Priority(CreatureEvent e) {
    switch (e) {
        case CreatureEvent::Telegraph: return AudioPriority::Critical; // an attack coming: never dropped (decided)
        case CreatureEvent::Alert: return AudioPriority::High;
        case CreatureEvent::Idle: return AudioPriority::Low;
        default: return AudioPriority::Normal;
    }
}

float CreatureSound_Pitch(const std::string& species) {
    const size_t h = std::hash<std::string>{}(species);
    return 0.88f + 0.24f * (float)(h % 1000) / 999.0f;
}

float CreatureSound_Danger(CreatureEvent e, float dist) {
    float level = 0.0f, range = 0.0f;
    switch (e) {
        case CreatureEvent::Telegraph: level = 1.0f; range = 30.0f; break; // an attack winding up near you
        case CreatureEvent::Hurt: level = 0.8f; range = 30.0f; break;      // a fight
        case CreatureEvent::Alert: level = 0.6f; range = 40.0f; break;     // it has noticed something
        default: return 0.0f;
    }
    if (dist >= range) return 0.0f;
    return level * std::clamp((range - dist) / (0.5f * range), 0.0f, 1.0f); // full within half the range
}

namespace {
    float threat = 0.0f;

    void OnEvent(const std::string& species, CreatureEvent e, uint64_t bodyId, Vector3 local, const std::string& rowSound) {
        if (!Audio_V2()) return;
        Planetoid* body = Planetoid::find(bodyId);
        if (!body) return;
        const DVec3 at = body->toWorldD(DVec3(local));
        const float dist = (float)(at - Audio_GetListener().pos).length();
        threat = std::max(threat, CreatureSound_Danger(e, dist));
        const CreatureSpec* spec = Creature_Spec(species);
        int snd = rowSound.empty() ? -1 : Audio_Sound(rowSound.c_str()); // the row's own sound when the asset slot has it
        if (snd < 0) snd = Audio_Sound(CreatureSound_Fallback(e, spec ? spec->faction : std::string()).c_str());
        if (snd < 0) return;
        AudioAttenuation att;
        att.refDistance = 3.0f;
        att.maxDistance = e == CreatureEvent::Idle ? 30.0f : e == CreatureEvent::Telegraph ? 50.0f : 45.0f;
        att.lowpassDistance = 35.0f;
        static const float kVolume[5] = { 0.5f, 0.8f, 1.0f, 0.8f, 0.9f }; // idle, alert, telegraph, hurt, death
        Audio_PlaySoundAt(snd, at, AudioBus::Creatures, CreatureSound_Priority(e), kVolume[(int)e], CreatureSound_Pitch(species), att);
    }
}

void CreatureSounds_Init() {
    static bool done = false;
    if (done) return;
    done = true;
    auto previous = gCreatureEvent; // (another system may listen too: both are called)
    gCreatureEvent = [previous](const std::string& species, CreatureEvent e, uint64_t body, Vector3 local, const std::string& sound) {
        if (previous) previous(species, e, body, local, sound);
        OnEvent(species, e, body, local, sound);
    };
}

void CreatureSounds_Update(float dt) {
    if (!Audio_V2()) return;
    if (threat > 0.0f) Music_SetDanger(threat);
    threat = std::max(0.0f, threat - dt / 10.0f); // a fight's tension lingers ~10 s after the last call
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1246): voices, names, priorities, pitch, danger ----
#include "unit_helpers.h"

TEST(CreatureSounds, NamesVoicesAndPriorities) {
    EXPECT_STREQ(CreatureSound_Voice("lurker"), "lurker");
    EXPECT_STREQ(CreatureSound_Voice("construct"), "construct");
    EXPECT_STREQ(CreatureSound_Voice("wild"), "beast");
    EXPECT_STREQ(CreatureSound_Voice(""), "beast");
    EXPECT_EQ(CreatureSound_Fallback(CreatureEvent::Telegraph, "lurker"), "creature_telegraph_lurker"); // Lighting's naming
    EXPECT_EQ(CreatureSound_Fallback(CreatureEvent::Death, "ooze"), "creature_death_ooze");
    EXPECT_EQ(CreatureSound_Priority(CreatureEvent::Telegraph), AudioPriority::Critical);
    EXPECT_EQ(CreatureSound_Priority(CreatureEvent::Idle), AudioPriority::Low);
    const float a = CreatureSound_Pitch("thornback"), b = CreatureSound_Pitch("dusk_stalker");
    EXPECT_GE(a, 0.88f); EXPECT_LE(a, 1.12f);
    EXPECT_EQ(a, CreatureSound_Pitch("thornback")); // stable
    EXPECT_NE(a, b);
}

TEST(CreatureSounds, DangerNearOnly) {
    EXPECT_FLOAT_EQ(CreatureSound_Danger(CreatureEvent::Telegraph, 5.0f), 1.0f);
    EXPECT_FLOAT_EQ(CreatureSound_Danger(CreatureEvent::Telegraph, 40.0f), 0.0f);
    EXPECT_GT(CreatureSound_Danger(CreatureEvent::Telegraph, 20.0f), 0.0f);
    EXPECT_LT(CreatureSound_Danger(CreatureEvent::Telegraph, 20.0f), 1.0f);
    EXPECT_FLOAT_EQ(CreatureSound_Danger(CreatureEvent::Alert, 10.0f), 0.6f);
    EXPECT_FLOAT_EQ(CreatureSound_Danger(CreatureEvent::Idle, 1.0f), 0.0f);
    EXPECT_FLOAT_EQ(CreatureSound_Danger(CreatureEvent::Death, 1.0f), 0.0f);
}
#endif
