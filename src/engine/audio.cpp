#include "audio.h"

#include "headless.h"
#include "utils/logger.h"
#include "profiler.h"
#include "steamaudio.h"
#include "audioscene.h"
#include "audiosynth.h"

#include <raylib.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <functional>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {
    constexpr int kRate = 44100;
    bool ready = false;
    float masterVolume = 0.35f;

    // ---- The old engine: one Sound per effect (a replay restarts it) ----
    Sound sounds[(int)Sfx::Count] = {};

    // ---- v2 (#1238, #1239) ----
    // Every sound is mono float PCM at 44.1 kHz (or a Music stream for long files). Every voice is its own AudioStream
    // whose callback (on miniaudio's thread) reads the sound it was given, loops it or ends it, and runs the one-pole
    // low-pass; raylib applies its volume, pan and pitch. (A static raylib Sound can neither loop nor be filtered.)
    struct SoundDef {
        std::string name;
        std::vector<float> pcm;
        bool stream = false;        // a Music (long beds, music pieces)
        Music music{};
        float seconds = 0.0f;
        bool fromFile = false, missing = false;
        int musicVoice = -1;        // the voice playing the stream (one at a time)
    };
    std::vector<SoundDef> table;    // index = sound id; the Sfx first

    struct VoiceDsp {               // shared with the audio thread
        const float* data = nullptr;    // written only while the stream is stopped (raylib's lock orders it)
        uint32_t len = 0;
        bool loop = false;
        uint32_t cursor = 0;            // audio thread
        float lp = 0.0f, lpR = 0.0f;    // audio thread
        std::atomic<float> alpha{ 1.0f };
        std::atomic<bool> done{ false };
        // Steam Audio (#1292): stereo streams; effects on 256-frame blocks, handed out through the FIFO.
        SteamAudio::VoiceFx* fx = nullptr;
        AudioBlockFifo fifo;
        float block[AudioBlockFifo::kBlock], left[AudioBlockFifo::kBlock], right[AudioBlockFifo::kBlock];
        std::atomic<float> dir[3] = { 0.0f, 0.0f, -1.0f };
        std::atomic<float> blend{ 0.0f };
        std::atomic<float> air[3] = { 1.0f, 1.0f, 1.0f };
        std::atomic<float> send{ 0.0f };    // (#1294) how much of the dry sound goes to the reverb
        uint64_t sendCursor = 0;            // audio thread: where its next block lands in the send ring
    };
    VoiceDsp dsp[AudioMixer::kVoices];
    AudioStream voiceStream[AudioMixer::kVoices] = {};
    std::shared_ptr<const AudioScene> worldScene; // (#1293) the ground and parts around the listener (Audio_SetWorld)
    static_assert(SteamAudio::kMaxVoices == AudioMixer::kVoices, "one Steam Audio source per voice");
    bool steam = false;                 // Steam Audio loaded: every voice is stereo and goes through its effects
    // (#1294) The reverb: the voices' sends summed in a ring, a stream of its own running Steam Audio's reverb on them.
    bool reverbLive = false;
    AudioSendRing sendRing;
    AudioStream reverbStream = {};
    AudioBlockFifo reverbFifo;
    float reverbBlock[AudioBlockFifo::kBlock], reverbL[AudioBlockFifo::kBlock], reverbR[AudioBlockFifo::kBlock];
    void ReverbCallback(void* out, unsigned int frames) {
        float* o = (float*)out;
        unsigned int done = 0;
        while (done < frames) {
            if (reverbFifo.available() == 0) {
                sendRing.read(reverbBlock, AudioBlockFifo::kBlock);
                SteamAudio::ReverbProcess(reverbBlock, reverbL, reverbR);
                reverbFifo.fill(reverbL, reverbR);
            }
            done += (unsigned int)reverbFifo.pop(o + 2 * done, (int)(frames - done));
        }
    }
    constexpr float kAirScale = 8.0f;   // air absorption x8: game-scale dulling with distance (guess, tune by ear)

    enum class VKind : uint8_t { None, Pcm, Stream };
    struct VoiceState {
        VKind kind = VKind::None;
        int sound = -1;
        int emitter = -1;               // the emitter it plays for, -1 a one-shot
        bool positional = false;
        DVec3 pos{};
        AudioAttenuation att{};
        float gain = -1.0f, pan = -1.0f, alpha = -1.0f, pitch = -1.0f; // last values handed to raylib
    };
    VoiceState vs[AudioMixer::kVoices];

    AudioMixer mixer;
    AudioEmitters emitters;
    AudioListener listener;
    bool emFramed[AudioEmitters::kMax] = {};
    FramePos emFrame[AudioEmitters::kMax] = {};
    std::vector<AudioEmitters::Action> actions;

    void FillVoice(VoiceDsp& d, float* out, unsigned int frames) {
        const float a = d.alpha.load(std::memory_order_relaxed);
        for (unsigned int i = 0; i < frames; ++i) {
            float x = 0.0f;
            if (d.cursor >= d.len && d.loop && d.len) d.cursor = 0;
            if (d.cursor < d.len) x = d.data[d.cursor++];
            else d.done.store(true, std::memory_order_relaxed);
            d.lp += a * (x - d.lp);
            out[i] = d.lp;
        }
    }
    // Steam Audio: one block of the sound, through the voice's effects (HRTF and air absorption; bypassed for a
    // non-positional sound), then the listener's muffle, into the FIFO.
    void RenderBlock(VoiceDsp& d) {
        constexpr int n = AudioBlockFifo::kBlock;
        for (int i = 0; i < n; ++i) {
            float x = 0.0f;
            if (d.cursor >= d.len && d.loop && d.len) d.cursor = 0;
            if (d.cursor < d.len) x = d.data[d.cursor++];
            else d.done.store(true, std::memory_order_relaxed);
            d.block[i] = x;
        }
        if (reverbLive) { // (#1294) the dry sound into the reverb's send ring
            const float s = d.send.load(std::memory_order_relaxed);
            if (s > 0.0f) sendRing.add(d.sendCursor, d.block, n, s);
        }
        const float blend = d.blend.load(std::memory_order_relaxed);
        if (d.fx && blend > 0.0f) {
            SteamAudio::Params p;
            for (int k = 0; k < 3; ++k) { p.dir[k] = d.dir[k].load(std::memory_order_relaxed); p.air[k] = d.air[k].load(std::memory_order_relaxed); }
            p.spatialBlend = blend;
            SteamAudio::Process(d.fx, d.block, d.left, d.right, p);
        } else {
            memcpy(d.left, d.block, sizeof d.left);
            memcpy(d.right, d.block, sizeof d.right);
        }
        const float a = d.alpha.load(std::memory_order_relaxed);
        if (a < 1.0f)
            for (int i = 0; i < n; ++i) { d.lp += a * (d.left[i] - d.lp); d.left[i] = d.lp; d.lpR += a * (d.right[i] - d.lpR); d.right[i] = d.lpR; }
        d.fifo.fill(d.left, d.right);
    }
    void FillVoiceStereo(VoiceDsp& d, float* out, unsigned int frames) {
        unsigned int done = 0;
        while (done < frames) {
            if (d.fifo.available() == 0) RenderBlock(d);
            done += (unsigned int)d.fifo.pop(out + 2 * done, (int)(frames - done));
        }
    }
    // raylib's stream callbacks carry no user pointer: one function per voice.
    template <int V> void VoiceCallback(void* out, unsigned int frames) {
        if (steam) FillVoiceStereo(dsp[V], (float*)out, frames);
        else FillVoice(dsp[V], (float*)out, frames);
    }
    template <int... I> constexpr std::array<AudioCallback, sizeof...(I)> MakeCallbacks(std::integer_sequence<int, I...>) { return { &VoiceCallback<I>... }; }
    constexpr auto kCallbacks = MakeCallbacks(std::make_integer_sequence<int, AudioMixer::kVoices>{});

    // f(t, u): t = seconds, u = 0..1 progress; returns -1..1. Scaled as the old 16-bit sounds were (x 30000 / 32768).
    template <typename F>
    std::vector<float> Synth(float seconds, F f) {
        const int frames = (int)(seconds * kRate);
        std::vector<float> pcm(frames);
        for (int i = 0; i < frames; ++i) {
            const float t = (float)i / kRate, u = (float)i / frames;
            pcm[i] = std::fmax(-1.0f, std::fmin(1.0f, f(t, u))) * (30000.0f / 32768.0f);
        }
        return pcm;
    }

    Sound SoundFromPcm(const std::vector<float>& pcm) { // the old engine's 16-bit Sound (the same samples as before, to an LSB)
        short* data = (short*)MemAlloc((unsigned int)(pcm.size() * sizeof(short)));
        for (size_t i = 0; i < pcm.size(); ++i) data[i] = (short)(pcm[i] * 32768.0f);
        Wave w = { (unsigned int)pcm.size(), (unsigned int)kRate, 16, 1, data };
        Sound s = LoadSoundFromWave(w);
        UnloadWave(w);
        return s;
    }

    // The seven effects, in the order the old engine made them (the same noise).
    void MakeProceduralSfx(std::vector<float> out[(int)Sfx::Count]) {
        std::mt19937 rng(1234);
        std::uniform_real_distribution<float> n(-1.0f, 1.0f);
        float lp = 0.0f;
        // Crunchy low-passed noise burst.
        out[(int)Sfx::Dig] = Synth(0.14f, [&](float, float u) { lp += (n(rng) - lp) * 0.25f; return lp * 1.6f * (1.0f - u) * (1.0f - u); });
        // Soft thump: falling sine plus a little noise.
        out[(int)Sfx::Build] = Synth(0.18f, [&](float t, float u) { return (sinf(6.2832f * (140.0f - 90.0f * u) * t) * 0.8f + n(rng) * 0.1f) * (1.0f - u); });
        // Two-tone chime.
        out[(int)Sfx::Craft] = Synth(0.3f, [&](float t, float u) { const float f = u < 0.4f ? 660.0f : 990.0f; return sinf(6.2832f * f * t) * 0.5f * (1.0f - u); });
        // Rising whoosh.
        out[(int)Sfx::Jump] = Synth(0.2f, [&](float t, float u) { lp += (n(rng) - lp) * (0.05f + 0.3f * u); return lp * 0.9f * sinf(3.1416f * u) + sinf(6.2832f * (200 + 300 * u) * t) * 0.1f; });
        out[(int)Sfx::Pickup] = Synth(0.08f, [&](float t, float u) { return sinf(6.2832f * 1320.0f * t) * 0.35f * (1.0f - u); });
        out[(int)Sfx::Error] = Synth(0.2f, [&](float t, float u) { return (fmodf(t * 180.0f, 1.0f) < 0.5f ? 0.3f : -0.3f) * (1.0f - u); });
        // (#226) A tier reached: a rising arpeggio (C E G, then a long high C), each note a struck bell (fast attack, decay).
        out[(int)Sfx::Fanfare] = Synth(1.6f, [&](float t, float) {
            static const float kNotes[4] = { 523.25f, 659.25f, 783.99f, 1046.5f };
            const int i = std::min(3, (int)(t / 0.18f));
            const float s = t - i * 0.18f;                        // time into this note
            const float env = std::min(1.0f, s * 200.0f) * expf(-s * (i == 3 ? 2.2f : 6.0f));
            const float f = kNotes[i];
            return (sinf(6.2832f * f * t) * 0.45f + sinf(6.2832f * 2.0f * f * t) * 0.12f) * env;
        });
    }

    // A seamless loop: `seconds` of sound whose last `xfade` seconds were blended into its start.
    template <typename F>
    std::vector<float> SynthLoop(float seconds, float xfade, F f) {
        std::vector<float> a = Synth(seconds + xfade, f);
        const int len = (int)(seconds * kRate), x = (int)a.size() - len;
        for (int i = 0; i < x; ++i) { const float w = (float)i / x; a[i] = a[i] * w + a[len + i] * (1.0f - w); }
        a.resize(len);
        return a;
    }

    // (#1240) Footsteps, one short sound per surface (footsteps.cpp plays them with a little pitch jitter): noise shaped by
    // a fast attack and a decay, filtered to the surface's colour, with its own extra (a click, a thump, grains, a ring).
    bool MakeProceduralStep(const std::string& surface, std::vector<float>& out) {
        std::mt19937 rng(std::hash<std::string>{}(surface) & 0xffffffffu);
        std::uniform_real_distribution<float> n(-1.0f, 1.0f), u01(0.0f, 1.0f);
        float lp = 0.0f, lp2 = 0.0f, grain = 0.0f;
        auto env = [](float t, float attack, float decay) { return std::min(1.0f, t / attack) * expf(-t * decay); };
        if (surface == "stone") { // a hard tap: bright noise and a click
            out = Synth(0.12f, [&](float t, float) { const float x = n(rng); lp += (x - lp) * 0.3f;
                return ((x - lp) * 0.9f * env(t, 0.001f, 55.0f)) + sinf(6.2832f * 2300.0f * t) * 0.25f * expf(-t * 220.0f); });
        } else if (surface == "dirt") { // a dull thud
            out = Synth(0.15f, [&](float t, float) { lp += (n(rng) - lp) * 0.08f;
                return lp * 2.2f * env(t, 0.003f, 35.0f) + sinf(6.2832f * 90.0f * t) * 0.35f * expf(-t * 40.0f); });
        } else if (surface == "grass") { // two quick rustles
            out = Synth(0.2f, [&](float t, float) { const float x = n(rng); lp += (x - lp) * 0.35f; lp2 += (x - lp2) * 0.06f;
                const float band = lp - lp2;
                return band * 1.4f * (env(t, 0.004f, 30.0f) + 0.6f * (t > 0.045f ? env(t - 0.045f, 0.004f, 35.0f) : 0.0f)); });
        } else if (surface == "sand") { // a soft, longer hiss
            out = Synth(0.22f, [&](float t, float) { lp += (n(rng) - lp) * 0.4f;
                return lp * 0.7f * env(t, 0.012f, 22.0f); });
        } else if (surface == "gravel") { // a crunch: dense little clicks over low noise
            out = Synth(0.2f, [&](float t, float) {
                if (u01(rng) < 500.0f / kRate) grain = 0.6f + 0.4f * u01(rng);
                grain *= 0.985f;
                lp += (n(rng) - lp) * 0.15f;
                return (grain * n(rng) * 0.8f + lp * 0.8f) * env(t, 0.002f, 20.0f); });
        } else if (surface == "snow") { // a squeaky crunch, dull
            out = Synth(0.24f, [&](float t, float) {
                if (u01(rng) < 250.0f / kRate) grain = 1.0f;
                grain *= 0.993f;
                lp += (n(rng) * grain - lp) * 0.12f;
                return lp * 2.5f * env(t, 0.01f, 16.0f); });
        } else if (surface == "wood") { // a hollow knock: a resonance near 170 Hz and its octave
            out = Synth(0.2f, [&](float t, float) { lp += (n(rng) - lp) * 0.25f;
                return (sinf(6.2832f * 170.0f * t) * 0.5f + sinf(6.2832f * 340.0f * t) * 0.2f) * expf(-t * 28.0f) + lp * 0.5f * expf(-t * 90.0f); });
        } else if (surface == "metal") { // a ringing clank: inharmonic partials
            out = Synth(0.35f, [&](float t, float) { lp += (n(rng) - lp) * 0.5f;
                return (sinf(6.2832f * 520.0f * t) * 0.3f * expf(-t * 12.0f) + sinf(6.2832f * 1310.0f * t) * 0.2f * expf(-t * 16.0f) +
                        sinf(6.2832f * 2150.0f * t) * 0.12f * expf(-t * 22.0f)) + lp * 0.4f * expf(-t * 120.0f); });
        } else if (surface == "water") { // a splash: noise and a few bubbles
            float bubF = 0.0f, bubT = 1.0f, bubPh = 0.0f;
            out = Synth(0.3f, [&](float t, float) { lp += (n(rng) - lp) * 0.3f;
                if (u01(rng) < 60.0f / kRate) { bubF = 400.0f + 700.0f * u01(rng); bubT = 0.0f; }
                float b = 0.0f;
                if (bubT < 0.03f) { bubPh += 6.2832f * bubF * (1.0f + 10.0f * bubT) / kRate; b = sinf(bubPh) * (1.0f - bubT / 0.03f) * 0.3f; bubT += 1.0f / kRate; }
                return lp * 1.2f * env(t, 0.004f, 12.0f) + b * expf(-t * 6.0f); });
        } else return false;
        return true;
    }

    // Procedural loops for emitters (later steps tune them and add more; files in resources/audio replace them).
    bool MakeProceduralLoop(const std::string& name, std::vector<float>& out) {
        std::mt19937 rng(std::hash<std::string>{}(name) & 0xffffffffu);
        std::uniform_real_distribution<float> n(-1.0f, 1.0f);
        if (name == "hum") { // a machine hum: 55 Hz and its harmonics, a slow wobble; whole cycles in 2 s, so it loops as is
            out = Synth(2.0f, [&](float t, float) {
                float s = 0.0f;
                for (int k = 1; k <= 6; ++k) s += sinf(6.2832f * 55.0f * k * t) / k * (k % 2 ? 1.0f : 0.5f);
                return 0.3f * s * (1.0f + 0.1f * sinf(6.2832f * 0.5f * t));
            });
            return true;
        }
        if (name == "wind") { // low-passed noise with a slowly moving cutoff and gusts
            float lp = 0.0f, lp2 = 0.0f;
            out = SynthLoop(6.0f, 1.0f, [&](float t, float) {
                const float g = 0.5f + 0.25f * sinf(6.2832f * 0.21f * t) + 0.25f * sinf(6.2832f * 0.13f * t + 1.0f);
                lp += (n(rng) - lp) * (0.01f + 0.04f * g);
                lp2 += (lp - lp2) * 0.2f;
                return lp2 * 3.0f * (0.4f + 0.6f * g);
            });
            return true;
        }
        if (name == "water") { // a brook: brown noise and bubbles (short rising chirps)
            float brown = 0.0f, hp = 0.0f, bubF = 0.0f, bubT = 1.0f, bubPh = 0.0f;
            std::uniform_real_distribution<float> u01(0.0f, 1.0f);
            out = SynthLoop(4.0f, 0.5f, [&](float, float) {
                brown = brown * 0.995f + n(rng) * 0.05f;
                hp += (brown - hp) * 0.02f;
                if (u01(rng) < 30.0f / kRate) { bubF = 300.0f + 600.0f * u01(rng); bubT = 0.0f; }
                float b = 0.0f;
                if (bubT < 0.04f) { bubPh += 6.2832f * bubF * (1.0f + 8.0f * bubT) / kRate; b = sinf(bubPh) * (1.0f - bubT / 0.04f) * 0.25f; bubT += 1.0f / kRate; }
                return (brown - hp) * 1.5f + b;
            });
            return true;
        }
        if (name.rfind("step_", 0) == 0) return MakeProceduralStep(name.substr(5), out); // (#1240)
        // (#1241) ambience beds (12 s loops of layered wind, birds, crickets, frogs, insects, surf, rain) and thunder
        const unsigned seed = (unsigned)(std::hash<std::string>{}(name) & 0xffffffffu);
        if (AudioBedRecipe r; AudioSynth_BedRecipe(name, r)) { out = AudioSynth_Bed(r, 12.0f, seed); return true; }
        if (name == "thunder_near") { out = AudioSynth_Thunder(1.0f, seed); return true; }
        if (name == "thunder_far") { out = AudioSynth_Thunder(0.0f, seed); return true; }
        return false;
    }

    bool LoadPcmFile(const std::string& path, std::vector<float>& out) {
        Wave w = LoadWave(path.c_str());
        if (!IsWaveValid(w)) return false;
        WaveFormat(&w, kRate, 32, 1);
        out.assign((const float*)w.data, (const float*)w.data + w.frameCount);
        UnloadWave(w);
        return !out.empty();
    }

    bool WantsStream(const std::string& name, const std::string& path) {
        std::string ext = std::filesystem::path(path).extension().string();
        for (char& c : ext) c = (char)tolower((unsigned char)c);
        if (ext == ".wav") return false;
        if (name.rfind("music_", 0) == 0 || name.rfind("bed_", 0) == 0) return true;
        std::error_code ec;
        return std::filesystem::file_size(path, ec) > (1u << 20) && !ec;
    }

    // The asset slot for one entry: a file replaces the procedural sound (or makes a sound that has none).
    void LoadFromFile(SoundDef& d, bool allowStream) {
        const std::string path = AudioFindAsset("resources/audio", d.name.c_str());
        if (path.empty()) return;
        if (allowStream && WantsStream(d.name, path)) {
            Music m = LoadMusicStream(path.c_str());
            if (!IsMusicValid(m)) { LOG_WARN("Audio", "could not stream %s\n", path.c_str()); return; }
            d.music = m;
            d.stream = true;
            d.pcm.clear();
            d.seconds = GetMusicTimeLength(m);
        } else {
            std::vector<float> pcm;
            if (!LoadPcmFile(path, pcm)) { LOG_WARN("Audio", "could not load %s: keeping the procedural sound\n", path.c_str()); return; }
            d.pcm = std::move(pcm);
            d.seconds = (float)d.pcm.size() / kRate;
        }
        d.fromFile = true;
        LOG_INFO("Audio", "asset slot: %s is '%s'%s\n", path.c_str(), d.name.c_str(), d.stream ? " (streamed)" : "");
    }

    void StopVoiceRaw(int v) {
        VoiceState& s = vs[v];
        if (s.kind == VKind::Pcm) StopAudioStream(voiceStream[v]);
        else if (s.kind == VKind::Stream && s.sound >= 0) {
            SoundDef& d = table[s.sound];
            if (d.musicVoice == v) { StopMusicStream(d.music); d.musicVoice = -1; }
        }
        s = VoiceState{};
    }

    void StartVoice(int v, int sound, bool loop, float pitch) {
        StopVoiceRaw(v);
        SoundDef& d = table[sound];
        VoiceState& s = vs[v];
        s.sound = sound;
        if (d.stream) {
            if (d.musicVoice >= 0 && d.musicVoice != v) { // one stream at a time: the newer play wins
                const int old = d.musicVoice;
                mixer.release(old);
                StopVoiceRaw(old);
            }
            d.music.looping = loop;
            StopMusicStream(d.music);  // rewinds
            PlayMusicStream(d.music);
            d.musicVoice = v;
            s.kind = VKind::Stream;
            return;
        }
        VoiceDsp& q = dsp[v];
        q.data = d.pcm.data();
        q.len = (uint32_t)d.pcm.size();
        q.loop = loop;
        q.cursor = 0;
        q.lp = q.lpR = 0.0f;
        q.done.store(false);
        q.alpha.store(1.0f);
        q.blend.store(0.0f);
        q.send.store(0.0f);
        q.sendCursor = 0; // (behind the reader: its first block resyncs it)
        q.fifo.clear();
        SteamAudio::ResetVoice(q.fx); // no HRTF tail from the last sound
        SetAudioStreamVolume(voiceStream[v], 0.0f); // set properly right after (ApplyVoice)
        SetAudioStreamPitch(voiceStream[v], pitch);
        PlayAudioStream(voiceStream[v]);
        s.kind = VKind::Pcm;
    }

    void ApplyVoice(int v, float gainIn, const AudioSpatial& sp, float pitch) {
        VoiceState& s = vs[v];
        // (#1310) A one-shot at a place gets its distance gain here (an emitter's voice volume has it already).
        const float gain = s.emitter >= 0 ? gainIn : gainIn * sp.gain;
        const float alpha = AudioOnePoleAlpha(sp.cutoffHz, (float)kRate);
        if (reverbLive && s.kind == VKind::Pcm) {
            // (#1294) The reverb send: the room's wet level, the sound's level with only half its distance loss (a
            // room's reverb is about as loud anywhere in it). Not for UI and music, nor for centred ambience beds.
            const AudioBus bus = mixer.voice(v).bus;
            const bool placed = s.positional || (s.emitter >= 0 && emitters.get(s.emitter) && emitters.get(s.emitter)->positional);
            const bool sends = bus != AudioBus::UI && bus != AudioBus::Music && (placed || bus != AudioBus::Ambience);
            const float dist = placed ? std::max(sp.gain, 1e-3f) : 1.0f;
            dsp[v].send.store(sends ? std::min(1.0f, SteamAudio::ReverbWet() * gain / sqrtf(dist)) : 0.0f, std::memory_order_relaxed);
        }
        if (s.kind == VKind::Pcm) {
            if (fabsf(gain - s.gain) > 1e-4f) SetAudioStreamVolume(voiceStream[v], gain);
            if (fabsf(sp.pan - s.pan) > 1e-3f) SetAudioStreamPan(voiceStream[v], sp.pan);
            if (fabsf(pitch - s.pitch) > 1e-4f) SetAudioStreamPitch(voiceStream[v], pitch);
            if (fabsf(alpha - s.alpha) > 1e-4f) dsp[v].alpha.store(alpha, std::memory_order_relaxed);
            if (steam) { // (#1292) the direction, how much is spatialised, air absorption by distance
                VoiceDsp& q = dsp[v];
                const bool placed = s.positional || (s.emitter >= 0 && emitters.get(s.emitter) && emitters.get(s.emitter)->positional);
                q.dir[0].store(sp.dir.x, std::memory_order_relaxed);
                q.dir[1].store(sp.dir.y, std::memory_order_relaxed);
                q.dir[2].store(sp.dir.z, std::memory_order_relaxed);
                q.blend.store(placed ? sp.blend : 0.0f, std::memory_order_relaxed);
                float air[3];
                AudioAirAbsorption(placed ? sp.distance : 0.0f, kAirScale, air);
                for (int k = 0; k < 3; ++k) q.air[k].store(air[k], std::memory_order_relaxed);
            }
        } else if (s.kind == VKind::Stream) {
            Music& m = table[s.sound].music;
            if (fabsf(gain - s.gain) > 1e-4f) SetMusicVolume(m, gain);
            if (fabsf(pitch - s.pitch) > 1e-4f) SetMusicPitch(m, pitch);
        }
        s.gain = gain; s.pan = sp.pan; s.pitch = pitch; s.alpha = alpha;
    }

    AudioSpatial VoiceSpatial(int v) {
        const VoiceState& s = vs[v];
        if (s.emitter >= 0) if (const AudioEmitters::Emitter* e = emitters.get(s.emitter)) return e->spatial;
        return s.positional ? Audio_Spatialize(listener, s.pos, s.att) : Audio_Spatialize2D(listener);
    }

    void PlayOneShot(int sound, AudioBus bus, AudioPriority priority, float volume, float pitch, const DVec3* pos, const AudioAttenuation& att) {
        if (sound < 0 || sound >= (int)table.size() || table[sound].missing) return;
        if (pos && Audio_Spatialize(listener, *pos, att).gain <= 0.0f) return; // out of earshot: no voice spent
        const SoundDef& d = table[sound];
        const int v = mixer.play(GetTime(), sound, bus, priority, volume, pitch, d.seconds);
        if (v < 0) return; // dropped: every voice holds a higher priority
        StartVoice(v, sound, false, pitch);
        vs[v].emitter = -1;
        vs[v].positional = pos != nullptr;
        if (pos) { vs[v].pos = *pos; vs[v].att = att; }
        ApplyVoice(v, mixer.gain(v), VoiceSpatial(v), mixer.voice(v).pitch);
    }

    bool V2Ready() { return ready && Audio_V2(); }
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
    if ((true)) return; // the Linux dedicated server build: no audio device at all (#1238)
#endif
    if (Headless()) return; // --host / --headless: a server never plays a sound
    InitAudioDevice();
    if (!IsAudioDeviceReady()) return;
    std::vector<float> pcm[(int)Sfx::Count];
    MakeProceduralSfx(pcm);
    if (Audio_V2()) {
        table.clear();
        for (int i = 0; i < (int)Sfx::Count; ++i) {
            SoundDef d;
            d.name = Audio_SfxName((Sfx)i);
            d.pcm = std::move(pcm[i]);
            LoadFromFile(d, false);
            d.seconds = (float)d.pcm.size() / kRate;
            table.push_back(std::move(d));
        }
        // (#1292) Steam Audio when its library is there: every voice stereo, through HRTF; else our own panning.
        steam = SteamAudio::Load();
        for (int v = 0; v < AudioMixer::kVoices && steam; ++v)
            if (!(dsp[v].fx = SteamAudio::CreateVoice(v))) {
                LOG_WARN("Audio", "steam audio: a voice's effects could not be made: using the fallback panning\n");
                for (int k = 0; k < v; ++k) { SteamAudio::DestroyVoice(dsp[k].fx); dsp[k].fx = nullptr; }
                SteamAudio::Shutdown();
                steam = false;
            }
        listener.hrtf = steam;
        if (steam) SteamAudio::StartSimulation(); // (#1293) occlusion and transmission (RAYLIBGAME_AUDIO_OCCLUSION=0: off)
        for (int v = 0; v < AudioMixer::kVoices; ++v) {
            voiceStream[v] = LoadAudioStream(kRate, 32, steam ? 2 : 1);
            SetAudioStreamCallback(voiceStream[v], kCallbacks[v]);
        }
        if (steam && SteamAudio::ReverbOn()) { // (#1294) the reverb's own stream, always playing (silent when dry)
            reverbStream = LoadAudioStream(kRate, 32, 2);
            SetAudioStreamCallback(reverbStream, ReverbCallback);
            PlayAudioStream(reverbStream);
            reverbLive = true;
        }
        SetMasterVolume(1.0f); // the master volume is a bus gain in v2 (the mixer applies it per voice)
        LOG_INFO("Audio", "audio v2: %d voices, buses, positional sound, emitters and the asset slot (RAYLIBGAME_AUDIO_V2); steam audio: %s\n",
                 AudioMixer::kVoices, SteamAudio::Status());
    } else {
        for (int i = 0; i < (int)Sfx::Count; ++i) sounds[i] = SoundFromPcm(pcm[i]);
        SetMasterVolume(masterVolume);
    }
    mixer.setBusVolume(AudioBus::Master, masterVolume);
    ready = true;
}

void Audio_Shutdown() {
    if (!ready) return;
    if (Audio_V2()) {
        for (int v = 0; v < AudioMixer::kVoices; ++v) { StopVoiceRaw(v); UnloadAudioStream(voiceStream[v]); }
        if (reverbLive) { StopAudioStream(reverbStream); UnloadAudioStream(reverbStream); reverbLive = false; }
        for (SoundDef& d : table) if (d.stream) UnloadMusicStream(d.music);
        table.clear();
        for (VoiceDsp& q : dsp) { SteamAudio::DestroyVoice(q.fx); q.fx = nullptr; } // the streams are gone: no callback runs
        SteamAudio::Shutdown();
        worldScene.reset();
        steam = false;
    } else {
        for (Sound& s : sounds) UnloadSound(s);
    }
    CloseAudioDevice();
    ready = false;
}

void Audio_Update() {
    if (!V2Ready()) return;
    PROFILE_SCOPE("audio:update");
    const double now = GetTime();
    static double last = now;
    const float dt = (float)std::min(0.25, std::max(0.0, now - last));
    last = now;
    // Finished one-shots (the callback reached the end; a stream stopped) are free again; streams are fed.
    for (int v = 0; v < AudioMixer::kVoices; ++v) {
        const AudioMixer::Voice& vo = mixer.voice(v);
        if (!vo.active) continue;
        if (vs[v].kind == VKind::Pcm && dsp[v].done.load(std::memory_order_relaxed)) mixer.release(v);
        else if (vs[v].kind == VKind::Stream) {
            Music& m = table[vs[v].sound].music;
            UpdateMusicStream(m);
            if (!IsMusicStreamPlaying(m) && now - vo.start > 0.1) mixer.release(v);
        }
    }
    // Emitters: places in a body's frame follow the body; then voices are taken and given back.
    for (int e = 0; e < AudioEmitters::kMax; ++e) {
        if (!emFramed[e]) continue;
        AudioEmitters::Emitter* em = emitters.get(e);
        DVec3 pg, vg;
        if (em && Frames_ToGlobal(emFrame[e], pg, vg)) em->pos = pg;
    }
    actions.clear();
    emitters.update(mixer, listener, now, dt, actions);
    for (const AudioEmitters::Action& a : actions) {
        if (a.kind == AudioEmitters::Action::Stop) { if (vs[a.voice].emitter == a.emitter) StopVoiceRaw(a.voice); }
        else if (const AudioEmitters::Emitter* em = emitters.get(a.emitter)) {
            StartVoice(a.voice, em->sound, true, em->pitch);
            vs[a.voice].emitter = a.emitter;
        }
    }
    mixer.update(now, dt);
    for (int v = 0; v < AudioMixer::kVoices; ++v) {
        if (!mixer.voice(v).active) { if (vs[v].kind != VKind::None) StopVoiceRaw(v); continue; }
        float pitch = mixer.voice(v).pitch;
        if (vs[v].emitter >= 0) if (const AudioEmitters::Emitter* em = emitters.get(vs[v].emitter)) pitch = em->pitch;
        ApplyVoice(v, mixer.gain(v), VoiceSpatial(v), pitch);
    }
    // (#1293) The occlusion simulation: where each placed voice is, the listener, the world snapshot.
    if (SteamAudio::Simulating() && worldScene) {
        PROFILE_SCOPE("audio:occlusion inputs");
        SteamAudio::SimVoice sv[SteamAudio::kMaxVoices];
        for (int v = 0; v < AudioMixer::kVoices; ++v) {
            if (!mixer.voice(v).active || vs[v].kind != VKind::Pcm) continue;
            if (vs[v].emitter >= 0) {
                const AudioEmitters::Emitter* em = emitters.get(vs[v].emitter);
                if (!em || !em->positional) continue;
                sv[v] = { true, em->pos };
            } else if (vs[v].positional) sv[v] = { true, vs[v].pos };
        }
        const float ahead[3] = { listener.forward.x, listener.forward.y, listener.forward.z }, up[3] = { listener.up.x, listener.up.y, listener.up.z };
        SteamAudio::SetSimInputs(worldScene, listener.pos, ahead, up, sv);
    }
}

void Audio_SetWorld(Scene& world) {
    // (#1293) A fresh snapshot of the ground and parts around the listener four times a second, or after 4 m.
    if (!V2Ready() || !SteamAudio::Simulating()) return;
    static double lastTime = -1e9;
    static DVec3 lastPos{ 1e30, 1e30, 1e30 };
    const double now = GetTime();
    if (now - lastTime < 0.25 && (listener.pos - lastPos).length() < 4.0) return;
    PROFILE_SCOPE("audio:scene snapshot");
    lastTime = now;
    lastPos = listener.pos;
    worldScene = AudioScene_Snapshot(world, listener.pos);
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
    PlayOneShot((int)s, bus, priority, volume, pitch, nullptr, {});
}

void Audio_Play(Sfx s, float volume, float pitch) {
    Audio_PlayOn(s, Audio_SfxBus(s), Audio_SfxPriority(s), volume, pitch);
}

// ---- Positional sound and emitters (#1239) ----

void Audio_SetListener(const DVec3& posGlobal, Vector3 forward, Vector3 up) {
    listener.pos = posGlobal;
    listener.forward = forward;
    listener.up = up;
}

const AudioListener& Audio_GetListener() { return listener; }
void Audio_SetMuffle(float cutoffHz) { listener.muffleHz = std::max(0.0f, cutoffHz); }

int Audio_Sound(const char* name) {
    if (!V2Ready() || !name || !*name) return -1;
    for (int i = 0; i < (int)table.size(); ++i)
        if (table[i].name == name) return table[i].missing ? -1 : i;
    SoundDef d;
    d.name = name;
    LoadFromFile(d, true);
    if (!d.fromFile) {
        if (MakeProceduralLoop(d.name, d.pcm)) d.seconds = (float)d.pcm.size() / kRate;
        else { d.missing = true; LOG_WARN("Audio", "no sound '%s' (no file in resources/audio, no procedural one)\n", name); }
    }
    table.push_back(std::move(d));
    return table.back().missing ? -1 : (int)table.size() - 1;
}

void Audio_PlayAt(Sfx s, const DVec3& posGlobal, float volume, float pitch, const AudioAttenuation& att) {
    if (!ready || (int)s < 0 || (int)s >= (int)Sfx::Count) return;
    if (!Audio_V2()) { // the old engine: 2D, quieter with distance
        const float g = AudioDistanceGain((float)(posGlobal - listener.pos).length(), att);
        if (g > 0.0f) Audio_PlayOn(s, Audio_SfxBus(s), Audio_SfxPriority(s), volume * g, pitch);
        return;
    }
    PROFILE_SCOPE("audio:play");
    PlayOneShot((int)s, Audio_SfxBus(s), Audio_SfxPriority(s), volume, pitch, &posGlobal, att);
}

void Audio_PlaySound(int sound, AudioBus bus, AudioPriority priority, float volume, float pitch) {
    if (V2Ready()) PlayOneShot(sound, bus, priority, volume, pitch, nullptr, {});
}

void Audio_PlaySoundAt(int sound, const DVec3& posGlobal, AudioBus bus, AudioPriority priority, float volume, float pitch, const AudioAttenuation& att) {
    if (V2Ready()) PlayOneShot(sound, bus, priority, volume, pitch, &posGlobal, att);
}

int Audio_EmitterCreate(int sound, AudioBus bus, AudioPriority priority, const AudioAttenuation& att, bool positional) {
    if (!V2Ready() || sound < 0 || sound >= (int)table.size() || table[sound].missing) return -1;
    const int e = emitters.create(sound, bus, priority, att, positional && !table[sound].stream); // streams play centred
    if (e >= 0) emFramed[e] = false;
    return e;
}

void Audio_EmitterDestroy(int e) {
    if (!V2Ready() || e < 0) return;
    const int v = emitters.destroy(mixer, e);
    if (v >= 0 && vs[v].emitter == e) StopVoiceRaw(v);
    if (e < AudioEmitters::kMax) emFramed[e] = false;
}

void Audio_EmitterSetPosition(int e, const DVec3& posGlobal) {
    if (AudioEmitters::Emitter* em = V2Ready() ? emitters.get(e) : nullptr) { em->pos = posGlobal; emFramed[e] = false; }
}

void Audio_EmitterSetFrame(int e, const FramePos& p) {
    if (AudioEmitters::Emitter* em = V2Ready() ? emitters.get(e) : nullptr) {
        emFrame[e] = p;
        emFramed[e] = true;
        DVec3 pg, vg;
        if (Frames_ToGlobal(p, pg, vg)) em->pos = pg;
    }
}

void Audio_EmitterSetPlaying(int e, bool playing) {
    if (AudioEmitters::Emitter* em = V2Ready() ? emitters.get(e) : nullptr) em->playing = playing;
}

void Audio_EmitterSetVolume(int e, float volume) {
    if (AudioEmitters::Emitter* em = V2Ready() ? emitters.get(e) : nullptr) em->volume = std::clamp(volume, 0.0f, 4.0f);
}

void Audio_EmitterSetPitch(int e, float pitch) {
    if (AudioEmitters::Emitter* em = V2Ready() ? emitters.get(e) : nullptr) em->pitch = std::clamp(pitch, 0.05f, 8.0f);
}

void Audio_EmitterStats(int& count, int& voiced) {
    count = emitters.count();
    voiced = emitters.voiced();
}

// ---- Volumes ----

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
    if (V2Ready()) mixer.duckFor(GetTime(), seconds);
}

AudioMixer::Stats Audio_Stats() { return mixer.stats(); }
const AudioMixer& Audio_Mixer() { return mixer; }

const char* Audio_SpatialBackend() {
    static std::string s;
    if (!V2Ready()) return "off (the old engine)";
    s = (steam ? std::string("steam audio ") : std::string("own panning; steam audio ")) + SteamAudio::Status();
    if (SteamAudio::Simulating()) { // (#1293)
        const SteamAudio::SimStats st = SteamAudio::GetSimStats();
        char buf[160];
        int occluded = 0;
        for (float o : st.occlusion) occluded += o < 0.99f ? 1 : 0;
        snprintf(buf, sizeof buf, "\nocclusion: %lld runs, %.2f ms (max %.2f), %d sources, %d occluded", st.ticks, st.lastMs, st.maxMs, st.sources, occluded);
        s += buf;
        if (reverbLive) { // (#1294)
            snprintf(buf, sizeof buf, "\nreverb: enclosed %.0f%%, mean %.1f m, RT60 %.2f/%.2f/%.2f s, wet %.2f",
                     st.enclosure * 100.0f, st.meanDistance, st.rt60[0], st.rt60[1], st.rt60[2], st.wet);
            s += buf;
        }
    }
    return s.c_str();
}
