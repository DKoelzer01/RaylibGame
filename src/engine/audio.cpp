#include "audio.h"

#include <raylib.h>
#include <cmath>
#include <random>

namespace {
    bool ready = false;
    Sound sounds[(int)Sfx::Count] = {};

    // f(t, u): t = seconds, u = 0..1 progress; returns -1..1.
    template <typename F>
    Sound Synth(float seconds, F f) {
        const int rate = 44100;
        const int frames = (int)(seconds * rate);
        short* data = (short*)MemAlloc(frames * sizeof(short));
        for (int i = 0; i < frames; ++i) {
            const float t = (float)i / rate, u = (float)i / frames;
            const float v = f(t, u);
            data[i] = (short)(std::fmax(-1.0f, std::fmin(1.0f, v)) * 30000.0f);
        }
        Wave w = { (unsigned int)frames, (unsigned int)rate, 16, 1, data };
        Sound s = LoadSoundFromWave(w);
        UnloadWave(w);
        return s;
    }
}

void Audio_Init() {
    InitAudioDevice();
    if (!IsAudioDeviceReady()) return;
    std::mt19937 rng(1234);
    std::uniform_real_distribution<float> n(-1.0f, 1.0f);
    float lp = 0.0f;
    // Crunchy low-passed noise burst.
    sounds[(int)Sfx::Dig] = Synth(0.14f, [&](float, float u) { lp += (n(rng) - lp) * 0.25f; return lp * 1.6f * (1.0f - u) * (1.0f - u); });
    // Soft thump: falling sine plus a little noise.
    sounds[(int)Sfx::Build] = Synth(0.18f, [&](float t, float u) { return (sinf(6.2832f * (140.0f - 90.0f * u) * t) * 0.8f + n(rng) * 0.1f) * (1.0f - u); });
    // Two-tone chime.
    sounds[(int)Sfx::Craft] = Synth(0.3f, [&](float t, float u) { const float f = u < 0.4f ? 660.0f : 990.0f; return sinf(6.2832f * f * t) * 0.5f * (1.0f - u); });
    // Rising whoosh.
    sounds[(int)Sfx::Jump] = Synth(0.2f, [&](float t, float u) { lp += (n(rng) - lp) * (0.05f + 0.3f * u); return lp * 0.9f * sinf(3.1416f * u) + sinf(6.2832f * (200 + 300 * u) * t) * 0.1f; });
    sounds[(int)Sfx::Pickup] = Synth(0.08f, [&](float t, float u) { return sinf(6.2832f * 1320.0f * t) * 0.35f * (1.0f - u); });
    sounds[(int)Sfx::Error] = Synth(0.2f, [&](float t, float u) { return (fmodf(t * 180.0f, 1.0f) < 0.5f ? 0.3f : -0.3f) * (1.0f - u); });
    ready = true;
}

void Audio_Shutdown() {
    if (!ready) return;
    for (Sound& s : sounds) UnloadSound(s);
    CloseAudioDevice();
    ready = false;
}

void Audio_Play(Sfx s, float volume, float pitch) {
    if (!ready) return;
    Sound& snd = sounds[(int)s];
    SetSoundVolume(snd, volume);
    SetSoundPitch(snd, pitch);
    PlaySound(snd);
}
