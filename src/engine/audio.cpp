#include "audio.h"

#include "headless.h"
#include "utils/logger.h"
#include "profiler.h"

#include <raylib.h>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <random>

namespace {
    bool ready = false;
    float masterVolume = 0.35f;
    Sound sounds[(int)Sfx::Count] = {};
    float soundSeconds[(int)Sfx::Count] = {};    // length at pitch 1 (the mixer's expiry)
    bool fromFile[(int)Sfx::Count] = {};

    // v2 (#1238): the bookkeeping and one alias per voice (an alias shares the sample data, so the same effect can play
    // on several voices at once; the old engine restarted the one Sound).
    AudioMixer mixer;
    Sound voiceAlias[AudioMixer::kVoices] = {};
    int voiceAliasOf[AudioMixer::kVoices];       // the Sfx its alias was made from (-1 none)
    float voiceGain[AudioMixer::kVoices] = {};   // last volume handed to raylib (re-applied when the bus gains move)

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

    void MakeProceduralSounds() {
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
        // (#226) A tier reached: a rising arpeggio (C E G, then a long high C), each note a struck bell (fast attack, decay).
        sounds[(int)Sfx::Fanfare] = Synth(1.6f, [&](float t, float) {
            static const float kNotes[4] = { 523.25f, 659.25f, 783.99f, 1046.5f };
            const int i = std::min(3, (int)(t / 0.18f));
            const float s = t - i * 0.18f;                        // time into this note
            const float env = std::min(1.0f, s * 200.0f) * expf(-s * (i == 3 ? 2.2f : 6.0f));
            const float f = kNotes[i];
            return (sinf(6.2832f * f * t) * 0.45f + sinf(6.2832f * 2.0f * f * t) * 0.12f) * env;
        });
    }

    // The asset slot (v2): resources/audio/<name>.<ext> replaces the procedural sound of that name.
    void LoadAssetOverrides() {
        for (int i = 0; i < (int)Sfx::Count; ++i) {
            const std::string path = AudioFindAsset("resources/audio", Audio_SfxName((Sfx)i));
            if (path.empty()) continue;
            Sound s = LoadSound(path.c_str());
            if (!IsSoundValid(s)) { LOG_WARN("Audio", "could not load %s: keeping the procedural sound\n", path.c_str()); continue; }
            UnloadSound(sounds[i]);
            sounds[i] = s;
            fromFile[i] = true;
            LOG_INFO("Audio", "asset slot: %s replaces the procedural '%s'\n", path.c_str(), Audio_SfxName((Sfx)i));
        }
    }

    void StopVoice(int v) {
        if (voiceAliasOf[v] >= 0 && IsSoundPlaying(voiceAlias[v])) StopSound(voiceAlias[v]);
    }
}

bool Audio_V2() {
    static const bool on = [] { const char* e = getenv("RAYLIBGAME_AUDIO_V2"); return e && e[0] == '1'; }();
    return on;
}

const char* Audio_SfxName(Sfx s) {
    static const char* kNames[] = { "dig", "build", "craft", "jump", "pickup", "error", "fanfare" };
    return (int)s < (int)Sfx::Count ? kNames[(int)s] : "?";
}

AudioBus Audio_SfxBus(Sfx s) {
    switch (s) {
        case Sfx::Craft: case Sfx::Error: case Sfx::Fanfare: return AudioBus::UI;
        default: return AudioBus::Effects; // dig, build, jump, pickup: tool and world feedback
    }
}

AudioPriority Audio_SfxPriority(Sfx s) {
    switch (s) {
        case Sfx::Dig: return AudioPriority::Low;      // many a second while digging
        case Sfx::Fanfare: return AudioPriority::High;
        default: return AudioPriority::Normal;
    }
}

void Audio_Init() {
#if RAYLIBGAME_HEADLESS_PLATFORM
    return; // the Linux dedicated server build: no audio device at all (#1238)
#endif
    if (Headless()) return; // --host / --headless: a server never plays a sound
    InitAudioDevice();
    if (!IsAudioDeviceReady()) return;
    MakeProceduralSounds();
    if (Audio_V2()) {
        LoadAssetOverrides();
        SetMasterVolume(1.0f); // the master volume is a bus gain in v2 (the mixer applies it per voice)
        for (int v = 0; v < AudioMixer::kVoices; ++v) voiceAliasOf[v] = -1;
        LOG_INFO("Audio", "audio v2: %d voices, buses and the asset slot (RAYLIBGAME_AUDIO_V2)\n", AudioMixer::kVoices);
    } else {
        SetMasterVolume(masterVolume);
    }
    for (int i = 0; i < (int)Sfx::Count; ++i)
        soundSeconds[i] = sounds[i].stream.sampleRate ? (float)sounds[i].frameCount / (float)sounds[i].stream.sampleRate : 0.0f;
    mixer.setBusVolume(AudioBus::Master, masterVolume);
    ready = true;
}

void Audio_Shutdown() {
    if (!ready) return;
    if (Audio_V2())
        for (int v = 0; v < AudioMixer::kVoices; ++v)
            if (voiceAliasOf[v] >= 0) { StopSound(voiceAlias[v]); UnloadSoundAlias(voiceAlias[v]); voiceAliasOf[v] = -1; }
    for (Sound& s : sounds) UnloadSound(s);
    CloseAudioDevice();
    ready = false;
}

void Audio_Update() {
    if (!ready || !Audio_V2()) return;
    PROFILE_SCOPE("audio:update");
    const double now = GetTime();
    static double last = now;
    const float dt = (float)std::min(0.25, std::max(0.0, now - last));
    last = now;
    // A voice raylib finished early (a file shorter than measured, a stop elsewhere) is free again.
    for (int v = 0; v < AudioMixer::kVoices; ++v) {
        const AudioMixer::Voice& vo = mixer.voice(v);
        if (vo.active && now - vo.start > 0.05 && !IsSoundPlaying(voiceAlias[v])) mixer.release(v);
    }
    mixer.update(now, dt);
    for (int v = 0; v < AudioMixer::kVoices; ++v) {
        if (!mixer.voice(v).active) continue;
        const float g = mixer.gain(v);
        if (fabsf(g - voiceGain[v]) > 1e-4f) { SetSoundVolume(voiceAlias[v], g); voiceGain[v] = g; }
    }
}

void Audio_PlayOn(Sfx s, AudioBus bus, AudioPriority priority, float volume, float pitch) {
    if (!ready || (int)s < 0 || (int)s >= (int)Sfx::Count) return;
    if (!Audio_V2()) {
        Sound& snd = sounds[(int)s];
        SetSoundVolume(snd, volume);
        SetSoundPitch(snd, pitch);
        PlaySound(snd);
        return;
    }
    PROFILE_SCOPE("audio:play");
    int stolen = -1;
    const int v = mixer.play(GetTime(), (int)s, bus, priority, volume, pitch, soundSeconds[(int)s], &stolen);
    if (v < 0) return; // dropped: every voice holds a higher priority
    if (stolen >= 0) StopVoice(v);
    if (voiceAliasOf[v] != (int)s) {
        if (voiceAliasOf[v] >= 0) { StopSound(voiceAlias[v]); UnloadSoundAlias(voiceAlias[v]); }
        voiceAlias[v] = LoadSoundAlias(sounds[(int)s]);
        voiceAliasOf[v] = (int)s;
    }
    voiceGain[v] = mixer.gain(v);
    SetSoundVolume(voiceAlias[v], voiceGain[v]);
    SetSoundPitch(voiceAlias[v], pitch);
    SetSoundPan(voiceAlias[v], 0.5f);
    PlaySound(voiceAlias[v]);
}

void Audio_Play(Sfx s, float volume, float pitch) {
    Audio_PlayOn(s, Audio_SfxBus(s), Audio_SfxPriority(s), volume, pitch);
}

void Audio_SetMasterVolume(float v) {
    masterVolume = v < 0.0f ? 0.0f : (v > 1.0f ? 1.0f : v);
    mixer.setBusVolume(AudioBus::Master, masterVolume);
    if (ready && !Audio_V2()) SetMasterVolume(masterVolume);
}

float Audio_GetMasterVolume() { return masterVolume; }

void Audio_SetBusVolume(AudioBus b, float v) {
    if (b == AudioBus::Master) Audio_SetMasterVolume(v);
    else mixer.setBusVolume(b, v);
}

float Audio_GetBusVolume(AudioBus b) { return b == AudioBus::Master ? masterVolume : mixer.busVolume(b); }

void Audio_DuckFor(double seconds) {
    if (ready && Audio_V2()) mixer.duckFor(GetTime(), seconds);
}

AudioMixer::Stats Audio_Stats() { return mixer.stats(); }
const AudioMixer& Audio_Mixer() { return mixer; }
