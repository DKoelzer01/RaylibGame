#ifndef AUDIOSPATIAL_H
#define AUDIOSPATIAL_H

// Positional sound and looping emitters (#1239, audio plan #1003 step 2). Pure: no raylib audio calls, unit tested at
// the end of audiospatial.cpp; audio.cpp plays what it decides. Positions are global doubles (planet scale); only the
// difference to the listener goes to float, so a sound 10 000 km from the origin pans as well as one at the origin.
// Simple model (the fallback once Steam Audio, #1020, is in): distance attenuation, constant-power-ish panning from the
// listener's right axis, a quieter and duller sound behind, and a one-pole low-pass that closes with distance.

#include "audiomixer.h"
#include "dvec3.h"

#include <cstdint>
#include <vector>

constexpr float kAudioNoFilterHz = 20000.0f;

struct AudioAttenuation {
    float refDistance = 2.0f;      // full volume inside this (m)
    float maxDistance = 60.0f;     // silent beyond it; fades out over its last 20%
    float rolloff = 1.0f;          // inverse distance: ref / (ref + rolloff (d - ref))
    float lowpassDistance = 40.0f; // the cutoff falls from 20 kHz at ref to 2 kHz here (air absorbs the highs)
};

struct AudioListener {
    DVec3 pos{};                   // global
    Vector3 forward{ 0, 0, -1 }, up{ 0, 1, 0 };
    float muffleHz = 0.0f;         // > 0: every sound low-passed to this (vacuum through the suit, underwater; step 5)
};

struct AudioSpatial {
    float gain = 1.0f;
    float pan = 0.5f;              // raylib's convention: 1 = left, 0 = right, 0.5 = centre (raudio's MixAudioFrames)
    float cutoffHz = kAudioNoFilterHz;
    float distance = 0.0f;
};

float AudioDistanceGain(float d, const AudioAttenuation& a);
AudioSpatial Audio_Spatialize(const AudioListener& l, const DVec3& src, const AudioAttenuation& a);
AudioSpatial Audio_Spatialize2D(const AudioListener& l); // non-positional: centred, only the listener's muffle
// One-pole low-pass coefficient (y += a (x - y)) for a cutoff at a sample rate; 1 = no filtering.
float AudioOnePoleAlpha(float cutoffHz, float sampleRate);

// Looping (or long) sounds in the world: machine hums, wind, rain, water. Emitters are virtual: one holds a voice only
// while it is audible and a voice is free (or a lower priority one can be taken); out of range it fades out and gives
// the voice back. With more audible emitters than voices, a louder emitter takes the voice of one at the same priority
// that is less than half as loud. Fades of kFadeSeconds avoid clicks.
class AudioEmitters {
public:
    static constexpr int kMax = 512;
    static constexpr float kFadeSeconds = 0.25f;
    static constexpr float kAudible = 0.001f;   // -60 dB

    struct Emitter {
        bool used = false;
        int sound = -1;
        AudioBus bus = AudioBus::Ambience;
        AudioPriority priority = AudioPriority::Low;
        AudioAttenuation att{};
        bool positional = true;
        DVec3 pos{};
        bool playing = false;               // the owner wants it heard (state: a machine running, rain falling)
        float volume = 1.0f, pitch = 1.0f;
        // Derived each update:
        float fade = 0.0f;                  // 0..1
        float audible = 0.0f;               // bus gain x volume x distance gain (0 when not playing)
        AudioSpatial spatial{};
        int voice = -1;
        uint32_t voiceSeq = 0;              // the mixer's seq when the voice was taken (a different one: it was stolen)
    };
    struct Action {
        enum Kind : uint8_t { Start, Stop } kind;
        int emitter, voice;
    };

    int create(int sound, AudioBus bus, AudioPriority priority, const AudioAttenuation& att, bool positional);
    // Frees the emitter; returns the voice it held (released in the mixer; the backend stops it) or -1.
    int destroy(AudioMixer& m, int e);
    Emitter* get(int e) { return e >= 0 && e < kMax && list[e].used ? &list[e] : nullptr; }
    const Emitter* get(int e) const { return e >= 0 && e < kMax && list[e].used ? &list[e] : nullptr; }
    // Spatialises every emitter, moves the fades, takes and gives back voices (appended to `out`), sets the voices'
    // volumes in the mixer. Call before AudioMixer::update.
    void update(AudioMixer& m, const AudioListener& l, double now, float dt, std::vector<Action>& out);
    int count() const;
    int voiced() const;

private:
    Emitter list[kMax];
    int highWater = 0;                      // one past the last used slot
};

#endif
