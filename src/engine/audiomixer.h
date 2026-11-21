#ifndef AUDIOMIXER_H
#define AUDIOMIXER_H

// The audio engine's bookkeeping (#1238, audio plan #1003 step 1): buses with their own volume, ducking of music and
// ambience under alarms and combat, and a fixed pool of voices with priorities. Pure data: no raylib or miniaudio calls,
// so it is unit tested (the end of audiomixer.cpp) and carries over when step 3 (#1020, Steam Audio) replaces how voices
// are actually mixed. audio.cpp is the backend that plays what this decides.

#include <cstdint>
#include <string>

enum class AudioBus : uint8_t { Master, Effects, Creatures, Machines, Ambience, Music, UI, Count };
// Critical: alarms and creature telegraphs. They are never dropped (they steal any other voice, the oldest critical one
// last) and they duck the music and ambience while they play.
enum class AudioPriority : uint8_t { Low, Normal, High, Critical };

const char* AudioBus_Name(AudioBus b);   // "Master", "Effects", ... (Settings labels)
const char* AudioBus_Key(AudioBus b);    // "master", "effects", ... (settings.txt keys)

class AudioMixer {
public:
    static constexpr int kVoices = 32;
    static constexpr float kDuckGain = 0.5f;      // -6 dB on the ducked buses (guess, tune by ear)
    static constexpr float kDuckAttack = 0.15f;   // seconds to reach the ducked level
    static constexpr float kDuckRelease = 0.8f;   // seconds back to full

    struct Voice {
        bool active = false;
        int sound = -1;                 // the backend's sound id
        AudioBus bus = AudioBus::Effects;
        AudioPriority priority = AudioPriority::Normal;
        float volume = 1.0f, pitch = 1.0f;
        double start = 0.0, end = 0.0;  // end: when it stops by itself (start + length / pitch); infinite for loops
        uint32_t seq = 0;               // play order (the oldest is stolen first)
    };
    struct Stats { int active = 0; long long played = 0, stolen = 0, dropped = 0; };

    AudioMixer();

    // A voice for a new sound, or -1 when it was dropped (every voice holds a higher priority). `stolenSound` gets the
    // sound id the chosen voice was playing when it was stolen from a still-playing sound (-1 otherwise): the backend
    // stops it. `lengthSeconds` < 0 = loops until released. Critical sounds duck for their length.
    int play(double now, int sound, AudioBus bus, AudioPriority priority, float volume, float pitch, float lengthSeconds, int* stolenSound = nullptr);
    void release(int voice);                      // stopped (by the backend, or a loop ended)
    void update(double now, float dt);            // expires finished voices, moves the duck envelope
    void duckFor(double now, double seconds);     // duck music and ambience until now + seconds (alarms, combat)

    void setBusVolume(AudioBus b, float v);       // clamped to 0..1
    float busVolume(AudioBus b) const { return busVol[(int)b]; }
    static bool ducks(AudioBus b) { return b == AudioBus::Music || b == AudioBus::Ambience; }
    float duckLevel() const { return duck; }      // 1 = not ducked, kDuckGain = fully ducked
    float busGain(AudioBus b) const;              // master x bus x duck
    float gain(int voice) const;                  // busGain(bus) x the voice's own volume
    const Voice& voice(int v) const { return voices[v]; }
    Stats stats() const;

private:
    Voice voices[kVoices];
    float busVol[(int)AudioBus::Count];
    float duck = 1.0f;
    double duckUntil = -1.0;
    uint32_t nextSeq = 1;
    long long played = 0, stolen = 0, dropped = 0;
    void expire(double now);
};

// The asset slot: resources/audio/<name>.wav|.ogg|.flac|.mp3 replaces the procedural sound of that name. Returns the first
// file found in `dir` (in that order), or "" when there is none.
std::string AudioFindAsset(const std::string& dir, const char* name);

#endif
