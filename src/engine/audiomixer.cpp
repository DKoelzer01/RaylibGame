#include "audiomixer.h"

#include <algorithm>
#include <filesystem>
#include <limits>

const char* AudioBus_Name(AudioBus b) {
    static const char* kNames[] = { "Master", "Effects", "Creatures", "Machines", "Ambience", "Music", "UI" };
    return (int)b < (int)AudioBus::Count ? kNames[(int)b] : "?";
}

const char* AudioBus_Key(AudioBus b) {
    static const char* kKeys[] = { "master", "effects", "creatures", "machines", "ambience", "music", "ui" };
    return (int)b < (int)AudioBus::Count ? kKeys[(int)b] : "?";
}

AudioMixer::AudioMixer() {
    for (float& v : busVol) v = 1.0f;
}

void AudioMixer::expire(double now) {
    for (Voice& v : voices)
        if (v.active && now >= v.end) v.active = false;
}

int AudioMixer::play(double now, int sound, AudioBus bus, AudioPriority priority, float volume, float pitch, float lengthSeconds, int* stolenSound, bool stealEqual) {
    if (stolenSound) *stolenSound = -1;
    expire(now);
    int pick = -1;
    for (int i = 0; i < kVoices && pick < 0; ++i)
        if (!voices[i].active) pick = i;
    if (pick < 0) {
        // Steal: the lowest priority that is not above the new sound's, the oldest of those.
        for (int i = 0; i < kVoices; ++i) {
            const Voice& v = voices[i];
            if (v.priority > priority || (!stealEqual && v.priority == priority)) continue;
            if (pick < 0 || v.priority < voices[pick].priority || (v.priority == voices[pick].priority && v.seq < voices[pick].seq)) pick = i;
        }
        if (pick < 0) { if (stealEqual) ++dropped; return -1; } // (an emitter waiting for a voice retries: not a drop)
        ++stolen;
        if (stolenSound) *stolenSound = voices[pick].sound;
    }
    Voice& v = voices[pick];
    v.active = true;
    v.sound = sound;
    v.bus = bus == AudioBus::Master || bus >= AudioBus::Count ? AudioBus::Effects : bus; // Master is a volume, not a bus to play on
    v.priority = priority;
    v.volume = std::clamp(volume, 0.0f, 4.0f);
    v.pitch = pitch > 0.01f ? pitch : 0.01f;
    v.start = now;
    v.end = lengthSeconds < 0.0f ? std::numeric_limits<double>::infinity() : now + lengthSeconds / v.pitch;
    v.seq = nextSeq++;
    ++played;
    if (priority == AudioPriority::Critical && lengthSeconds > 0.0f) duckFor(now, v.end - now);
    return pick;
}

void AudioMixer::release(int voice) {
    if (voice >= 0 && voice < kVoices) voices[voice].active = false;
}

void AudioMixer::setVolume(int voice, float volume) {
    if (voice >= 0 && voice < kVoices) voices[voice].volume = std::clamp(volume, 0.0f, 4.0f);
}

void AudioMixer::duckFor(double now, double seconds) {
    duckUntil = std::max(duckUntil, now + seconds);
}

void AudioMixer::update(double now, float dt) {
    expire(now);
    bool ducked = now < duckUntil;
    for (const Voice& v : voices)
        if (v.active && v.priority == AudioPriority::Critical) ducked = true; // a looping alarm ducks while it plays
    const float target = ducked ? kDuckGain : 1.0f;
    const float rate = (1.0f - kDuckGain) / (ducked ? kDuckAttack : kDuckRelease); // level per second
    if (duck < target) duck = std::min(target, duck + rate * dt);
    else if (duck > target) duck = std::max(target, duck - rate * dt);
}

void AudioMixer::setBusVolume(AudioBus b, float v) {
    if ((int)b < (int)AudioBus::Count) busVol[(int)b] = std::clamp(v, 0.0f, 1.0f);
}

float AudioMixer::busGain(AudioBus b) const {
    if (b == AudioBus::Master) return busVol[0];
    return busVol[0] * busVol[(int)b] * (ducks(b) ? duck : 1.0f);
}

float AudioMixer::gain(int voice) const {
    if (voice < 0 || voice >= kVoices || !voices[voice].active) return 0.0f;
    return busGain(voices[voice].bus) * voices[voice].volume;
}

AudioMixer::Stats AudioMixer::stats() const {
    Stats s;
    for (const Voice& v : voices) s.active += v.active ? 1 : 0;
    s.played = played; s.stolen = stolen; s.dropped = dropped;
    return s;
}

std::string AudioFindAsset(const std::string& dir, const char* name) {
    static const char* kExt[] = { ".wav", ".ogg", ".flac", ".mp3" };
    std::error_code ec;
    for (const char* e : kExt) {
        const std::filesystem::path p = std::filesystem::path(dir) / (std::string(name) + e);
        if (std::filesystem::is_regular_file(p, ec)) return p.generic_string();
    }
    return "";
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1238): the voice pool, priorities, bus gains, ducking, the asset slot ----
#include "unit_helpers.h"
#include <cstdio>

namespace {
int Fill(AudioMixer& m, double now, AudioPriority p, float len = 1.0f) {
    int n = 0;
    for (int i = 0; i < AudioMixer::kVoices; ++i) n += m.play(now, i, AudioBus::Effects, p, 1.0f, 1.0f, len) >= 0;
    return n;
}
} // namespace

TEST(AudioMixer, FreeVoicesFirstThenStealOldestOfLowestPriority) {
    AudioMixer m;
    EXPECT_EQ(Fill(m, 0.0, AudioPriority::Normal), AudioMixer::kVoices);
    EXPECT_EQ(m.stats().active, AudioMixer::kVoices);
    EXPECT_EQ(m.stats().stolen, 0);
    int stolenSound = -2;
    const int v = m.play(0.1, 100, AudioBus::UI, AudioPriority::Normal, 1.0f, 1.0f, 1.0f, &stolenSound);
    ASSERT_GE(v, 0);
    EXPECT_EQ(stolenSound, 0);                    // the oldest normal one (sound 0) went
    EXPECT_EQ(m.voice(v).sound, 100);
    EXPECT_EQ(m.stats().stolen, 1);
    // A low-priority voice is stolen before any older normal one.
    m.release(5);
    ASSERT_EQ(m.play(0.2, 200, AudioBus::Effects, AudioPriority::Low, 1.0f, 1.0f, 1.0f), 5);
    m.play(0.3, 300, AudioBus::Effects, AudioPriority::High, 1.0f, 1.0f, 1.0f, &stolenSound);
    EXPECT_EQ(stolenSound, 200);
}

TEST(AudioMixer, LowerPriorityIsDroppedWhenEveryVoiceIsHigher) {
    AudioMixer m;
    Fill(m, 0.0, AudioPriority::High);
    EXPECT_EQ(m.play(0.0, 99, AudioBus::Effects, AudioPriority::Normal, 1.0f, 1.0f, 1.0f), -1);
    EXPECT_EQ(m.stats().dropped, 1);
    EXPECT_GE(m.play(0.0, 98, AudioBus::Effects, AudioPriority::High, 1.0f, 1.0f, 1.0f), 0); // same priority: steals
}

TEST(AudioMixer, CriticalIsNeverDropped) {
    AudioMixer m;
    Fill(m, 0.0, AudioPriority::Critical, -1.0f); // 32 looping alarms
    for (int k = 0; k < 10; ++k) EXPECT_GE(m.play(1.0, 500 + k, AudioBus::UI, AudioPriority::Critical, 1.0f, 1.0f, 0.5f), 0);
    EXPECT_EQ(m.stats().dropped, 0);
    EXPECT_EQ(m.play(1.0, 600, AudioBus::UI, AudioPriority::High, 1.0f, 1.0f, 0.5f), -1);
}

TEST(AudioMixer, VoicesExpireByLengthOverPitch) {
    AudioMixer m;
    const int v = m.play(10.0, 1, AudioBus::Effects, AudioPriority::Normal, 1.0f, 2.0f, 1.0f); // 1 s at double speed
    m.update(10.49, 0.016f);
    EXPECT_TRUE(m.voice(v).active);
    m.update(10.51, 0.016f);
    EXPECT_FALSE(m.voice(v).active);
    const int loop = m.play(0.0, 2, AudioBus::Ambience, AudioPriority::Low, 1.0f, 1.0f, -1.0f);
    m.update(1e6, 0.016f);
    EXPECT_TRUE(m.voice(loop).active);           // loops play until released
    m.release(loop);
    EXPECT_EQ(m.stats().active, 0);
}

TEST(AudioMixer, BusGainsMultiplyAndClamp) {
    AudioMixer m;
    m.setBusVolume(AudioBus::Master, 0.5f);
    m.setBusVolume(AudioBus::Effects, 0.4f);
    m.setBusVolume(AudioBus::UI, 3.0f);           // clamped to 1
    EXPECT_FLOAT_EQ(m.busVolume(AudioBus::UI), 1.0f);
    const int v = m.play(0.0, 1, AudioBus::Effects, AudioPriority::Normal, 0.5f, 1.0f, 1.0f);
    EXPECT_FLOAT_EQ(m.gain(v), 0.5f * 0.4f * 0.5f);
    const int u = m.play(0.0, 2, AudioBus::Master, AudioPriority::Normal, 1.0f, 1.0f, 1.0f);
    EXPECT_EQ(m.voice(u).bus, AudioBus::Effects); // nothing plays on Master itself
    m.release(v);
    EXPECT_FLOAT_EQ(m.gain(v), 0.0f);
}

TEST(AudioMixer, CriticalDucksMusicAndAmbienceOnly) {
    AudioMixer m;
    double t = 0.0;
    m.play(t, 1, AudioBus::UI, AudioPriority::Critical, 1.0f, 1.0f, 2.0f); // a 2 s alarm
    for (int i = 0; i < 30; ++i) m.update(t += 0.01, 0.01f);                // 0.3 s: past the attack
    EXPECT_FLOAT_EQ(m.duckLevel(), AudioMixer::kDuckGain);
    EXPECT_FLOAT_EQ(m.busGain(AudioBus::Music), AudioMixer::kDuckGain);
    EXPECT_FLOAT_EQ(m.busGain(AudioBus::Ambience), AudioMixer::kDuckGain);
    EXPECT_FLOAT_EQ(m.busGain(AudioBus::Effects), 1.0f);
    for (int i = 0; i < 170; ++i) m.update(t += 0.01, 0.01f);               // the alarm ends at 2 s
    EXPECT_LT(m.duckLevel(), 1.0f);                                          // releasing
    for (int i = 0; i < 100; ++i) m.update(t += 0.01, 0.01f);               // + 1 s > release time
    EXPECT_FLOAT_EQ(m.duckLevel(), 1.0f);
    // duckFor (combat) without a sound
    m.duckFor(t, 1.0);
    for (int i = 0; i < 30; ++i) m.update(t += 0.01, 0.01f);
    EXPECT_FLOAT_EQ(m.duckLevel(), AudioMixer::kDuckGain);
}

TEST(AudioMixer, AssetSlotFindsFilesInOrder) {
    const std::string dir = UnitTest::ScratchDir();
    EXPECT_EQ(AudioFindAsset(dir, "dig"), "");
    auto touch = [&](const char* f) { FILE* fp = fopen((dir + "/" + f).c_str(), "wb"); ASSERT_NE(fp, nullptr); fputs("x", fp); fclose(fp); };
    touch("dig.ogg");
    EXPECT_NE(AudioFindAsset(dir, "dig").find("dig.ogg"), std::string::npos);
    touch("dig.wav");
    EXPECT_NE(AudioFindAsset(dir, "dig").find("dig.wav"), std::string::npos); // .wav first
    EXPECT_EQ(AudioFindAsset(dir, "build"), "");
}
#endif
