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
    bool hrtf = false;             // Steam Audio places the sounds (#1292): no pan, no behind model, only dir and blend
};

struct AudioSpatial {
    float gain = 1.0f;
    float pan = 0.5f;              // raylib's convention: 1 = left, 0 = right, 0.5 = centre (raudio's MixAudioFrames)
    float cutoffHz = kAudioNoFilterHz;
    float distance = 0.0f;
    Vector3 dir{ 0, 0, -1 };       // hrtf: Steam Audio's listener space (AudioSteamDirection)
    float blend = 0.0f;            // hrtf: spatial blend (0 for a centred, non-positional sound)
};

float AudioDistanceGain(float d, const AudioAttenuation& a);
AudioSpatial Audio_Spatialize(const AudioListener& l, const DVec3& src, const AudioAttenuation& a);
AudioSpatial Audio_Spatialize2D(const AudioListener& l); // non-positional: centred, only the listener's muffle
// One-pole low-pass coefficient (y += a (x - y)) for a cutoff at a sample rate; 1 = no filtering.
float AudioOnePoleAlpha(float cutoffHz, float sampleRate);

// ---- Steam Audio (#1292) helpers that need no DLL ----
// The direction to a source in Steam Audio's listener space (x right, y up, -z ahead), unit length; straight ahead when
// the source is at the listener. Computed from the double difference (iplCalculateRelativeDirection takes floats).
Vector3 AudioSteamDirection(const AudioListener& l, const DVec3& src);
// Air absorption per band (low, mid, high) over d metres: Steam Audio's default exponential model, its coefficients
// times `scale` (> 1: game-scale dulling, so 40 m sounds farther than real air would make it).
void AudioAirAbsorption(float d, float scale, float out[3]);

// Stereo frames rendered in fixed blocks (Steam Audio's frame size) and handed out in whatever counts raylib's stream
// callback asks for.
class AudioBlockFifo {
public:
    static constexpr int kBlock = 256;
    int available() const { return count - pos; }
    // Takes a rendered block (left and right, kBlock frames each); only when empty.
    void fill(const float* left, const float* right) {
        for (int i = 0; i < kBlock; ++i) { buf[2 * i] = left[i]; buf[2 * i + 1] = right[i]; }
        pos = 0;
        count = kBlock;
    }
    // Copies up to `frames` interleaved stereo frames to out; returns how many.
    int pop(float* out, int frames) {
        const int n = frames < available() ? frames : available();
        for (int i = 0; i < 2 * n; ++i) out[i] = buf[2 * pos + i];
        pos += n;
        return n;
    }
    void clear() { pos = count = 0; }

private:
    float buf[2 * kBlock] = {};
    int pos = 0, count = 0;
};

// (#1294) The reverb's send bus: voices add their blocks at their own write cursor, the reverb reads (and clears) at its
// read cursor. Everything runs on the audio thread (raylib calls every stream's callback from miniaudio's one mixing
// callback), so there are no locks. A writer that falls behind the reader or runs too far ahead (a voice at another
// pitch consumes frames at another rate) is put back kLead frames ahead of it: a dropped or doubled block in a reverb
// tail is inaudible.
class AudioSendRing {
public:
    static constexpr int kSize = 8192, kLead = 512;
    void add(uint64_t& cursor, const float* x, int n, float gain) {
        if (cursor < readPos || cursor + (uint64_t)n > readPos + kSize - 256) cursor = readPos + kLead;
        for (int i = 0; i < n; ++i) buf[(cursor + i) % kSize] += x[i] * gain;
        cursor += n;
    }
    void read(float* out, int n) {
        for (int i = 0; i < n; ++i) { float& s = buf[(readPos + i) % kSize]; out[i] = s; s = 0.0f; }
        readPos += n;
    }
    uint64_t readPosition() const { return readPos; }

private:
    float buf[kSize] = {};
    uint64_t readPos = 0;
};

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
