#include "audiospatial.h"

#include <algorithm>
#include <cmath>

float AudioDistanceGain(float d, const AudioAttenuation& a) {
    const float ref = std::max(0.01f, a.refDistance);
    if (d >= a.maxDistance) return 0.0f;
    float g = d <= ref ? 1.0f : ref / (ref + a.rolloff * (d - ref));
    const float fadeStart = 0.8f * a.maxDistance;
    if (d > fadeStart) g *= (a.maxDistance - d) / (a.maxDistance - fadeStart);
    return g;
}

AudioSpatial Audio_Spatialize2D(const AudioListener& l) {
    AudioSpatial s;
    if (l.muffleHz > 0.0f) s.cutoffHz = std::min(s.cutoffHz, l.muffleHz);
    return s;
}

AudioSpatial Audio_Spatialize(const AudioListener& l, const DVec3& src, const AudioAttenuation& a) {
    AudioSpatial s = Audio_Spatialize2D(l);
    const DVec3 relD = src - l.pos;                    // in double: planet-scale positions cancel exactly enough
    const float rx = (float)relD.x, ry = (float)relD.y, rz = (float)relD.z;
    const float d = sqrtf(rx * rx + ry * ry + rz * rz);
    s.distance = d;
    s.gain = AudioDistanceGain(d, a);
    if (l.hrtf) { // Steam Audio places it (#1292): only the direction and how much of it is spatialised
        s.dir = AudioSteamDirection(l, src);
        s.blend = std::clamp(d / std::max(0.01f, a.refDistance), 0.0f, 1.0f); // close sounds surround
        return s;
    }
    if (d < 1e-3f) return s;
    const Vector3 dir = { rx / d, ry / d, rz / d };
    auto norm = [](Vector3 v) { const float n = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); return n > 1e-6f ? Vector3{ v.x / n, v.y / n, v.z / n } : Vector3{ 0, 0, 0 }; };
    const Vector3 f = norm(l.forward);
    const Vector3 right = norm({ f.y * l.up.z - f.z * l.up.y, f.z * l.up.x - f.x * l.up.z, f.x * l.up.y - f.y * l.up.x }); // forward x up
    const float x = dir.x * right.x + dir.y * right.y + dir.z * right.z;    // +1 right, -1 left
    const float z = dir.x * f.x + dir.y * f.y + dir.z * f.z;                // +1 ahead, -1 behind
    const float spread = std::clamp(d / std::max(0.01f, a.refDistance), 0.0f, 1.0f); // close sounds surround: less pan
    s.pan = 0.5f - 0.45f * x * spread;                 // never fully in one ear
    const float behind = std::max(0.0f, -z);
    s.gain *= 1.0f - 0.3f * behind;                    // the head's shadow
    const float ref = std::max(0.01f, a.refDistance);
    const float t = a.lowpassDistance > ref ? std::clamp((d - ref) / (a.lowpassDistance - ref), 0.0f, 1.0f) : 0.0f;
    float cutoff = kAudioNoFilterHz * powf(0.1f, t);
    cutoff *= 1.0f - 0.5f * behind;
    s.cutoffHz = std::min(s.cutoffHz, cutoff);
    return s;
}

float AudioOnePoleAlpha(float cutoffHz, float sampleRate) {
    if (cutoffHz >= kAudioNoFilterHz || sampleRate <= 0.0f) return 1.0f;
    return 1.0f - expf(-6.2831853f * std::max(10.0f, cutoffHz) / sampleRate);
}

Vector3 AudioSteamDirection(const AudioListener& l, const DVec3& src) {
    const DVec3 relD = src - l.pos;
    const float rx = (float)relD.x, ry = (float)relD.y, rz = (float)relD.z;
    const float d = sqrtf(rx * rx + ry * ry + rz * rz);
    if (d < 1e-4f) return { 0, 0, -1 };
    auto norm = [](Vector3 v) { const float n = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); return n > 1e-6f ? Vector3{ v.x / n, v.y / n, v.z / n } : Vector3{ 0, 0, 0 }; };
    const Vector3 f = norm(l.forward);
    const Vector3 r = norm({ f.y * l.up.z - f.z * l.up.y, f.z * l.up.x - f.x * l.up.z, f.x * l.up.y - f.y * l.up.x });
    const Vector3 u = { r.y * f.z - r.z * f.y, r.z * f.x - r.x * f.z, r.x * f.y - r.y * f.x }; // right x forward: the true up
    const Vector3 dir = { rx / d, ry / d, rz / d };
    return { dir.x * r.x + dir.y * r.y + dir.z * r.z, dir.x * u.x + dir.y * u.y + dir.z * u.z, -(dir.x * f.x + dir.y * f.y + dir.z * f.z) };
}

void AudioAirAbsorption(float d, float scale, float out[3]) {
    static const float kCoef[3] = { 0.0002f, 0.0017f, 0.0182f }; // Steam Audio's default model (1/m)
    for (int b = 0; b < 3; ++b) out[b] = expf(-kCoef[b] * scale * std::max(0.0f, d));
}

// ---- Emitters ----

int AudioEmitters::create(int sound, AudioBus bus, AudioPriority priority, const AudioAttenuation& att, bool positional) {
    for (int i = 0; i < kMax; ++i) {
        if (list[i].used) continue;
        Emitter& e = list[i];
        e = Emitter{};
        e.used = true;
        e.sound = sound;
        e.bus = bus;
        e.priority = priority;
        e.att = att;
        e.positional = positional;
        highWater = std::max(highWater, i + 1);
        return i;
    }
    return -1;
}

int AudioEmitters::destroy(AudioMixer& m, int e) {
    Emitter* em = get(e);
    if (!em) return -1;
    int v = -1;
    if (em->voice >= 0 && m.voice(em->voice).active && m.voice(em->voice).seq == em->voiceSeq) { v = em->voice; m.release(v); }
    *em = Emitter{};
    while (highWater > 0 && !list[highWater - 1].used) --highWater;
    return v;
}

void AudioEmitters::update(AudioMixer& m, const AudioListener& l, double now, float dt, std::vector<Action>& out) {
    // Spatialise, fade, give back voices no longer wanted.
    int want[kMax];
    int nWant = 0;
    for (int i = 0; i < highWater; ++i) {
        Emitter& e = list[i];
        if (!e.used) continue;
        if (e.voice >= 0 && !(m.voice(e.voice).active && m.voice(e.voice).seq == e.voiceSeq)) e.voice = -1; // stolen
        e.spatial = e.positional ? Audio_Spatialize(l, e.pos, e.att) : Audio_Spatialize2D(l);
        e.audible = e.playing ? m.busGain(e.bus) * e.volume * e.spatial.gain : 0.0f;
        const bool wanted = e.audible > kAudible;
        const float step = dt / kFadeSeconds;
        e.fade = wanted ? std::min(1.0f, e.fade + step) : std::max(0.0f, e.fade - step);
        if (!wanted && e.fade <= 0.0f && e.voice >= 0) {
            m.release(e.voice);
            out.push_back({ Action::Stop, i, e.voice });
            e.voice = -1;
        }
        if (wanted && e.voice < 0) want[nWant++] = i;
    }
    // Take voices, the loudest first.
    std::sort(want, want + nWant, [&](int a, int b) { return list[a].audible > list[b].audible; });
    for (int k = 0; k < nWant; ++k) {
        Emitter& e = list[want[k]];
        int v = m.play(now, e.sound, e.bus, e.priority, 0.0f, e.pitch, -1.0f, nullptr, false);
        if (v < 0) {
            // Every voice is busy at this priority or above: take the quietest same-priority emitter's voice if it is
            // under half as loud as this one.
            int quiet = -1;
            for (int i = 0; i < highWater; ++i) {
                const Emitter& o = list[i];
                if (!o.used || o.voice < 0 || o.priority != e.priority) continue;
                if (quiet < 0 || o.audible < list[quiet].audible) quiet = i;
            }
            if (quiet < 0 || list[quiet].audible >= 0.5f * e.audible) continue;
            Emitter& q = list[quiet];
            m.release(q.voice);
            out.push_back({ Action::Stop, quiet, q.voice });
            q.voice = -1;
            q.fade = 0.0f;
            v = m.play(now, e.sound, e.bus, e.priority, 0.0f, e.pitch, -1.0f, nullptr, false);
            if (v < 0) continue;
        }
        e.voice = v;
        e.voiceSeq = m.voice(v).seq;
        e.fade = std::min(e.fade, dt / kFadeSeconds); // fade in from (near) silence on a fresh voice
        out.push_back({ Action::Start, want[k], v });
    }
    // The voices' own volume: distance gain x the emitter's volume x the fade (the mixer adds the bus gains).
    for (int i = 0; i < highWater; ++i) {
        const Emitter& e = list[i];
        if (e.used && e.voice >= 0) m.setVolume(e.voice, e.volume * e.spatial.gain * e.fade);
    }
}

int AudioEmitters::count() const {
    int n = 0;
    for (int i = 0; i < highWater; ++i) n += list[i].used ? 1 : 0;
    return n;
}

int AudioEmitters::voiced() const {
    int n = 0;
    for (int i = 0; i < highWater; ++i) n += list[i].used && list[i].voice >= 0 ? 1 : 0;
    return n;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1239): attenuation, panning, the low-pass, emitter virtualisation ----
#include "unit_helpers.h"

namespace {
AudioListener ListenerAt(DVec3 p) {
    AudioListener l;
    l.pos = p;
    l.forward = { 0, 0, -1 };   // looking down -z, y up: +x is to the right
    l.up = { 0, 1, 0 };
    return l;
}
} // namespace

TEST(AudioSpatial, DistanceGainCurve) {
    AudioAttenuation a; // ref 2, max 60
    EXPECT_FLOAT_EQ(AudioDistanceGain(0.0f, a), 1.0f);
    EXPECT_FLOAT_EQ(AudioDistanceGain(2.0f, a), 1.0f);
    EXPECT_FLOAT_EQ(AudioDistanceGain(4.0f, a), 0.5f);
    EXPECT_NEAR(AudioDistanceGain(20.0f, a), 0.1f, 1e-6f);
    EXPECT_LT(AudioDistanceGain(55.0f, a), 2.0f / 55.0f);  // fading out over the last 20%
    EXPECT_FLOAT_EQ(AudioDistanceGain(60.0f, a), 0.0f);
    EXPECT_FLOAT_EQ(AudioDistanceGain(1e6f, a), 0.0f);
    float prev = 2.0f;
    for (float d = 0.0f; d < 70.0f; d += 0.5f) { const float g = AudioDistanceGain(d, a); EXPECT_LE(g, prev); prev = g; }
}

TEST(AudioSpatial, PanFollowsTheListenersRightAndWorksAtPlanetScale) {
    for (const DVec3 origin : { DVec3{ 0, 0, 0 }, DVec3{ 6.4e6, -3.1e6, 2.2e7 } }) {
        const AudioListener l = ListenerAt(origin);
        const AudioSpatial right = Audio_Spatialize(l, origin + DVec3{ 10, 0, 0 }, {});
        const AudioSpatial left = Audio_Spatialize(l, origin + DVec3{ -10, 0, 0 }, {});
        const AudioSpatial ahead = Audio_Spatialize(l, origin + DVec3{ 0, 0, -10 }, {});
        EXPECT_LT(right.pan, 0.1f);                     // raylib: 0 = right
        EXPECT_GT(left.pan, 0.9f);
        EXPECT_NEAR(ahead.pan, 0.5f, 1e-4f);
        EXPECT_NEAR(ahead.distance, 10.0f, 1e-3f);
        EXPECT_NEAR(right.gain, left.gain, 1e-6f);
    }
}

TEST(AudioSpatial, BehindIsQuieterAndDuller) {
    const AudioListener l = ListenerAt({ 0, 0, 0 });
    const AudioSpatial ahead = Audio_Spatialize(l, { 0, 0, -10 }, {});
    const AudioSpatial behind = Audio_Spatialize(l, { 0, 0, 10 }, {});
    EXPECT_LT(behind.gain, ahead.gain);
    EXPECT_LT(behind.cutoffHz, ahead.cutoffHz);
}

TEST(AudioSpatial, LowpassClosesWithDistanceAndMuffle) {
    AudioListener l = ListenerAt({ 0, 0, 0 });
    const AudioSpatial nearS = Audio_Spatialize(l, { 0, 0, -1 }, {});
    const AudioSpatial farS = Audio_Spatialize(l, { 0, 0, -40 }, {});
    EXPECT_FLOAT_EQ(nearS.cutoffHz, kAudioNoFilterHz);
    EXPECT_NEAR(farS.cutoffHz, 2000.0f, 1.0f);
    EXPECT_FLOAT_EQ(AudioOnePoleAlpha(kAudioNoFilterHz, 44100.0f), 1.0f);
    EXPECT_LT(AudioOnePoleAlpha(2000.0f, 44100.0f), AudioOnePoleAlpha(8000.0f, 44100.0f));
    l.muffleHz = 600.0f;                                // the suit in vacuum
    EXPECT_FLOAT_EQ(Audio_Spatialize(l, { 0, 0, -1 }, {}).cutoffHz, 600.0f);
    EXPECT_FLOAT_EQ(Audio_Spatialize2D(l).cutoffHz, 600.0f);
}

TEST(AudioEmitters, VoiceOnlyWhileAudibleWithFades) {
    AudioMixer m;
    AudioEmitters em;
    AudioListener l = ListenerAt({ 0, 0, 0 });
    std::vector<AudioEmitters::Action> acts;
    const int e = em.create(7, AudioBus::Machines, AudioPriority::Normal, {}, true);
    em.get(e)->pos = { 0, 0, -5 };
    em.update(m, l, 0.0, 0.016f, acts);
    EXPECT_TRUE(acts.empty());                          // not playing yet
    em.get(e)->playing = true;
    em.update(m, l, 0.016, 0.016f, acts);
    ASSERT_EQ(acts.size(), 1u);
    EXPECT_EQ(acts[0].kind, AudioEmitters::Action::Start);
    const int v = acts[0].voice;
    EXPECT_EQ(m.voice(v).sound, 7);
    EXPECT_LT(m.voice(v).volume, 0.1f);                 // fading in
    double t = 0.016;
    for (int i = 0; i < 30; ++i) { acts.clear(); em.update(m, l, t += 0.016, 0.016f, acts); m.update(t, 0.016f); }
    EXPECT_TRUE(acts.empty());
    EXPECT_NEAR(m.voice(v).volume, AudioDistanceGain(5.0f, {}), 1e-5f);
    EXPECT_TRUE(m.voice(v).active);                     // loops never expire
    // Walk out of range: it fades and gives the voice back.
    l.pos = { 0, 0, 500 };
    int stops = 0;
    for (int i = 0; i < 30; ++i) { acts.clear(); em.update(m, l, t += 0.016, 0.016f, acts); for (auto& a : acts) stops += a.kind == AudioEmitters::Action::Stop; }
    EXPECT_EQ(stops, 1);
    EXPECT_FALSE(m.voice(v).active);
    EXPECT_EQ(em.voiced(), 0);
    EXPECT_EQ(em.count(), 1);
    EXPECT_EQ(em.destroy(m, e), -1);
    EXPECT_EQ(em.count(), 0);
}

TEST(AudioEmitters, LouderTakesQuieterSamePriorityVoiceAndStolenVoicesAreNoticed) {
    AudioMixer m;
    AudioEmitters em;
    const AudioListener l = ListenerAt({ 0, 0, 0 });
    std::vector<AudioEmitters::Action> acts;
    // 32 far emitters fill the pool.
    for (int i = 0; i < AudioMixer::kVoices; ++i) {
        const int e = em.create(1, AudioBus::Ambience, AudioPriority::Low, {}, true);
        em.get(e)->pos = { 0, 0, -40.0 - i * 0.1 };
        em.get(e)->playing = true;
    }
    em.update(m, l, 0.0, 0.016f, acts);
    EXPECT_EQ(em.voiced(), AudioMixer::kVoices);
    // A near one of the same priority takes the quietest's voice (the farthest one).
    const int nearE = em.create(2, AudioBus::Ambience, AudioPriority::Low, {}, true);
    em.get(nearE)->pos = { 0, 0, -3 };
    em.get(nearE)->playing = true;
    acts.clear();
    em.update(m, l, 0.016, 0.016f, acts);
    EXPECT_GE(em.get(nearE)->voice, 0);
    EXPECT_EQ(em.get(AudioMixer::kVoices - 1)->voice, -1);
    EXPECT_EQ(em.voiced(), AudioMixer::kVoices);
    // A similar-loudness newcomer does not (no thrashing at equal loudness).
    const int twin = em.create(3, AudioBus::Ambience, AudioPriority::Low, {}, true);
    em.get(twin)->pos = { 0, 0, -40.0 };
    em.get(twin)->playing = true;
    acts.clear();
    em.update(m, l, 0.032, 0.016f, acts);
    EXPECT_EQ(em.get(twin)->voice, -1);
    // A one-shot alarm steals an emitter's voice: the emitter notices and goes virtual.
    int stolen = -1;
    const int alarmV = m.play(0.05, 99, AudioBus::UI, AudioPriority::Critical, 1.0f, 1.0f, 1.0f, &stolen);
    ASSERT_GE(alarmV, 0);
    acts.clear();
    em.update(m, l, 0.05, 0.016f, acts);
    int holders = 0;
    for (int i = 0; i < AudioEmitters::kMax; ++i) if (const auto* e = em.get(i); e && e->voice == alarmV) ++holders;
    EXPECT_EQ(holders, 0);
}

TEST(AudioSteam, DirectionInListenerSpace) {
    // Listener looking down -z, y up (Steam Audio's own convention, so directions come out unchanged) and at planet scale.
    for (const DVec3 o : { DVec3{ 0, 0, 0 }, DVec3{ 6.4e6, 1.0e5, -2.0e7 } }) {
        const AudioListener l = ListenerAt(o);
        const Vector3 r = AudioSteamDirection(l, o + DVec3{ 5, 0, 0 });
        EXPECT_NEAR(r.x, 1.0f, 1e-4f); EXPECT_NEAR(r.y, 0.0f, 1e-4f); EXPECT_NEAR(r.z, 0.0f, 1e-4f);
        const Vector3 a = AudioSteamDirection(l, o + DVec3{ 0, 0, -5 });
        EXPECT_NEAR(a.z, -1.0f, 1e-4f);
        const Vector3 up = AudioSteamDirection(l, o + DVec3{ 0, 5, 0 });
        EXPECT_NEAR(up.y, 1.0f, 1e-4f);
    }
    // A listener turned to look along +x with +z up: a source along +x is ahead, along -y is to its right.
    AudioListener t = ListenerAt({ 0, 0, 0 });
    t.forward = { 1, 0, 0 };
    t.up = { 0, 0, 1 };
    EXPECT_NEAR(AudioSteamDirection(t, { 3, 0, 0 }).z, -1.0f, 1e-4f);
    EXPECT_NEAR(AudioSteamDirection(t, { 0, -3, 0 }).x, 1.0f, 1e-4f);
    EXPECT_NEAR(AudioSteamDirection(t, { 0, 0, 3 }).y, 1.0f, 1e-4f);
    const Vector3 same = AudioSteamDirection(t, { 0, 0, 0 });
    EXPECT_FLOAT_EQ(same.z, -1.0f);
}

TEST(AudioSteam, HrtfListenerLeavesPanningToSteamAudio) {
    AudioListener l = ListenerAt({ 0, 0, 0 });
    l.hrtf = true;
    const AudioSpatial behind = Audio_Spatialize(l, { 0, 0, 10 }, {});
    EXPECT_FLOAT_EQ(behind.pan, 0.5f);
    EXPECT_FLOAT_EQ(behind.gain, AudioDistanceGain(10.0f, {}));   // no behind model
    EXPECT_FLOAT_EQ(behind.cutoffHz, kAudioNoFilterHz);
    EXPECT_NEAR(behind.dir.z, 1.0f, 1e-4f);
    EXPECT_FLOAT_EQ(behind.blend, 1.0f);
    EXPECT_FLOAT_EQ(Audio_Spatialize(l, { 0, 0, -1 }, {}).blend, 0.5f); // inside the reference distance: partly dry
    EXPECT_FLOAT_EQ(Audio_Spatialize2D(l).blend, 0.0f);
}

TEST(AudioSteam, AirAbsorptionFallsFasterForHighBands) {
    float a[3];
    AudioAirAbsorption(0.0f, 8.0f, a);
    for (float v : a) EXPECT_FLOAT_EQ(v, 1.0f);
    AudioAirAbsorption(40.0f, 8.0f, a);
    EXPECT_GT(a[0], a[1]);
    EXPECT_GT(a[1], a[2]);
    EXPECT_NEAR(a[2], expf(-0.0182f * 8.0f * 40.0f), 1e-6f);
}

TEST(AudioSteam, BlockFifoHandsOutAnyCounts) {
    AudioBlockFifo f;
    float l[AudioBlockFifo::kBlock], r[AudioBlockFifo::kBlock];
    for (int i = 0; i < AudioBlockFifo::kBlock; ++i) { l[i] = (float)i; r[i] = -(float)i; }
    EXPECT_EQ(f.available(), 0);
    f.fill(l, r);
    float out[2 * 300];
    EXPECT_EQ(f.pop(out, 100), 100);
    EXPECT_FLOAT_EQ(out[0], 0.0f); EXPECT_FLOAT_EQ(out[1], -0.0f);
    EXPECT_FLOAT_EQ(out[2 * 99], 99.0f); EXPECT_FLOAT_EQ(out[2 * 99 + 1], -99.0f);
    EXPECT_EQ(f.pop(out, 300), AudioBlockFifo::kBlock - 100);
    EXPECT_FLOAT_EQ(out[0], 100.0f);
    EXPECT_EQ(f.available(), 0);
    EXPECT_EQ(f.pop(out, 10), 0);
}

TEST(AudioSendRing, WritersSumAtTheirCursorsAndResync) {
    static AudioSendRing r; // (32 KB: not on the test's stack)
    float one[256], two[256], out[1024];
    for (int i = 0; i < 256; ++i) { one[i] = 1.0f; two[i] = 2.0f; }
    uint64_t a = 0, b = 0;                       // fresh writers start kLead ahead of the reader
    r.add(a, one, 256, 0.5f);
    r.add(b, two, 256, 1.0f);
    EXPECT_EQ(a, (uint64_t)AudioSendRing::kLead + 256);
    EXPECT_EQ(b, a);
    r.read(out, AudioSendRing::kLead);           // the lead: silence
    for (int i = 0; i < AudioSendRing::kLead; ++i) ASSERT_FLOAT_EQ(out[i], 0.0f);
    r.read(out, 256);                            // both writers' blocks, summed with their gains
    for (int i = 0; i < 256; ++i) ASSERT_FLOAT_EQ(out[i], 2.5f);
    r.read(out, 256);                            // read samples were cleared: nothing comes round again
    for (int i = 0; i < 256; ++i) ASSERT_FLOAT_EQ(out[i], 0.0f);
    // A writer the reader has passed is put back ahead of it.
    r.add(a, one, 256, 1.0f);
    EXPECT_EQ(a, r.readPosition() + AudioSendRing::kLead + 256);
    // One that ran too far ahead too.
    uint64_t far = r.readPosition() + AudioSendRing::kSize;
    r.add(far, one, 256, 1.0f);
    EXPECT_EQ(far, r.readPosition() + AudioSendRing::kLead + 256);
}
#endif
