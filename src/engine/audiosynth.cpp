#include "audiosynth.h"

#include <algorithm>
#include <cmath>
#include <random>

namespace {
    constexpr int kRate = 44100;
    constexpr float kTau = 6.2831853f;

    struct Chirp { float t = 0.0f, dur = 0.0f, f0 = 0.0f, sweep = 0.0f, vib = 0.0f, phase = 0.0f, amp = 0.0f; };
    struct Croak { float t = 0.0f, dur = 0.0f, f = 0.0f, phase = 0.0f, amp = 0.0f; };

    // Peak-normalise to `peak` (the layers' sum varies by recipe).
    void Normalise(std::vector<float>& a, float peak) {
        float m = 0.0f;
        for (float v : a) m = std::max(m, fabsf(v));
        if (m > 1e-6f) for (float& v : a) v *= peak / m;
    }
}

bool AudioSynth_BedRecipe(const std::string& name, AudioBedRecipe& r) {
    r = AudioBedRecipe{};
    if (name == "bed_forest_day") { r.wind = 0.4f; r.birds = 1.0f; r.insects = 0.1f; }
    else if (name == "bed_forest_night") { r.wind = 0.3f; r.crickets = 1.0f; }
    else if (name == "bed_meadow_day") { r.wind = 0.5f; r.birds = 0.6f; r.insects = 0.3f; }
    else if (name == "bed_jungle_day") { r.wind = 0.1f; r.birds = 1.6f; r.insects = 0.9f; }
    else if (name == "bed_jungle_night") { r.wind = 0.1f; r.crickets = 1.2f; r.frogs = 0.8f; r.insects = 0.5f; }
    else if (name == "bed_swamp_day") { r.wind = 0.2f; r.birds = 0.5f; r.insects = 0.6f; r.frogs = 0.4f; }
    else if (name == "bed_swamp_night") { r.wind = 0.1f; r.frogs = 1.3f; r.crickets = 0.6f; }
    else if (name == "bed_savanna_day") { r.wind = 0.7f; r.birds = 0.3f; r.insects = 0.5f; }
    else if (name == "bed_savanna_night") { r.wind = 0.5f; r.crickets = 0.8f; }
    else if (name == "bed_desert") { r.wind = 1.0f; r.insects = 0.1f; }
    else if (name == "bed_desert_night") { r.wind = 0.8f; r.crickets = 0.25f; }
    else if (name == "bed_taiga") { r.wind = 0.8f; r.birds = 0.25f; }
    else if (name == "bed_taiga_night") { r.wind = 0.8f; r.crickets = 0.1f; }
    else if (name == "bed_cold") { r.wind = 1.0f; }
    else if (name == "bed_shore") { r.wind = 0.4f; r.surf = 1.0f; r.birds = 0.15f; }
    else if (name == "bed_rain") { r.rain = 1.0f; }
    else if (name == "bed_rain_leaves") { r.leaves = 1.0f; r.rain = 0.3f; }
    else return false;
    return true;
}

std::vector<float> AudioSynth_Bed(const AudioBedRecipe& r, float seconds, unsigned seed) {
    const float xfade = std::min(1.5f, seconds * 0.25f);
    const int len = (int)(seconds * kRate), x = (int)(xfade * kRate), total = len + x;
    std::vector<float> a(total, 0.0f);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> n(-1.0f, 1.0f), u(0.0f, 1.0f);
    const float dt = 1.0f / kRate;
    // Layer state.
    float wLp = 0.0f, wLp2 = 0.0f;
    std::vector<Chirp> chirps;
    int phraseLeft = 0;
    float phraseGap = 0.0f, nextPhrase = 0.2f;
    struct Cricket { float f, period, next, phase; int pulses; float pulseT; };
    std::vector<Cricket> crickets;
    for (int k = 0; k < 3; ++k) crickets.push_back({ 4200.0f + 800.0f * u(rng), 0.45f + 0.35f * u(rng), u(rng), 0.0f, 0, 0.0f });
    std::vector<Croak> croaks;
    float nextCroak = 0.3f;
    float iLp = 0.0f, iLp2 = 0.0f;
    float sLp = 0.0f;
    float rLp = 0.0f, drop = 0.0f, leaf = 0.0f, leafLp = 0.0f;
    for (int i = 0; i < total; ++i) {
        const float t = i * dt;
        float s = 0.0f;
        if (r.wind > 0.0f) { // low-passed noise with slow gusts
            const float g = 0.5f + 0.25f * sinf(kTau * 0.21f * t) + 0.25f * sinf(kTau * 0.13f * t + 1.0f);
            wLp += (n(rng) - wLp) * (0.01f + 0.04f * g);
            wLp2 += (wLp - wLp2) * 0.2f;
            s += r.wind * wLp2 * 3.0f * (0.4f + 0.6f * g);
        }
        if (r.birds > 0.0f) { // phrases of 2-5 chirps, a few a second at full
            if (phraseLeft == 0 && t >= nextPhrase) { phraseLeft = 2 + (int)(u(rng) * 4.0f); phraseGap = 0.0f; nextPhrase = t + (0.6f + 1.6f * u(rng)) / r.birds; }
            if (phraseLeft > 0 && (phraseGap -= dt) <= 0.0f) {
                Chirp c;
                c.dur = 0.05f + 0.15f * u(rng);
                c.f0 = 2200.0f + 2400.0f * u(rng);
                c.sweep = (u(rng) - 0.4f) * 2000.0f;
                c.vib = 15.0f + 25.0f * u(rng);
                c.amp = 0.12f + 0.12f * u(rng);
                chirps.push_back(c);
                --phraseLeft;
                phraseGap = c.dur + 0.04f + 0.1f * u(rng);
            }
            for (Chirp& c : chirps) {
                const float uu = c.t / c.dur;
                const float f = c.f0 + c.sweep * uu + 150.0f * sinf(kTau * c.vib * c.t);
                c.phase += kTau * f * dt;
                s += r.birds * c.amp * sinf(c.phase) * sinf(3.14159f * uu);
                c.t += dt;
            }
            chirps.erase(std::remove_if(chirps.begin(), chirps.end(), [](const Chirp& c) { return c.t >= c.dur; }), chirps.end());
        }
        if (r.crickets > 0.0f) { // each cricket: three 25 Hz pulses, then a pause
            for (Cricket& c : crickets) {
                c.phase += kTau * c.f * dt;
                if (c.pulses == 0 && t >= c.next) { c.pulses = 3; c.pulseT = 0.0f; c.next = t + c.period; }
                if (c.pulses > 0) {
                    const float within = fmodf(c.pulseT, 0.04f);
                    const float env = within < 0.02f ? sinf(3.14159f * within / 0.02f) : 0.0f;
                    s += r.crickets * 0.06f * env * sinf(c.phase);
                    c.pulseT += dt;
                    if (c.pulseT >= 0.04f * 3) c.pulses = 0;
                }
            }
        }
        if (r.frogs > 0.0f) { // croaks: a pulsed low tone with harmonics
            if (t >= nextCroak) { croaks.push_back({ 0.0f, 0.18f + 0.15f * u(rng), 230.0f + 220.0f * u(rng), 0.0f, 0.2f + 0.15f * u(rng) }); nextCroak = t + (0.4f + 1.2f * u(rng)) / r.frogs; }
            for (Croak& c : croaks) {
                c.phase += kTau * c.f * dt;
                const float am = 0.5f + 0.5f * sinf(kTau * 32.0f * c.t);
                const float env = sinf(3.14159f * c.t / c.dur);
                s += r.frogs * c.amp * env * am * (sinf(c.phase) + 0.4f * sinf(2.0f * c.phase) + 0.2f * sinf(3.0f * c.phase));
                c.t += dt;
            }
            croaks.erase(std::remove_if(croaks.begin(), croaks.end(), [](const Croak& c) { return c.t >= c.dur; }), croaks.end());
        }
        if (r.insects > 0.0f) { // a high band of noise, buzzing and swelling
            const float x0 = n(rng);
            iLp += (x0 - iLp) * 0.55f;
            iLp2 += (x0 - iLp2) * 0.35f;
            const float band = iLp - iLp2;
            s += r.insects * band * 0.5f * (0.6f + 0.4f * sinf(kTau * 0.3f * t)) * (0.7f + 0.3f * sinf(kTau * 110.0f * t));
        }
        if (r.surf > 0.0f) { // waves: swells of low-passed noise every ~8 s
            const float sw = std::max(0.0f, sinf(kTau * t / 7.7f));
            const float env = 0.25f + 0.75f * sw * sw;
            sLp += (n(rng) - sLp) * (0.02f + 0.08f * sw);
            s += r.surf * sLp * 2.5f * env;
        }
        if (r.rain > 0.0f) { // hiss and drops on the ground
            const float x0 = n(rng);
            rLp += (x0 - rLp) * 0.2f;
            if (u(rng) < 220.0f / kRate) drop = 0.3f + 0.7f * u(rng);
            drop *= 0.992f;
            s += r.rain * ((x0 - rLp) * 0.25f + drop * n(rng) * 0.35f);
        }
        if (r.leaves > 0.0f) { // many small ticks, duller (drops on leaves)
            if (u(rng) < 700.0f / kRate) leaf = 0.4f + 0.6f * u(rng);
            leaf *= 0.985f;
            leafLp += (leaf * n(rng) - leafLp) * 0.35f;
            s += r.leaves * leafLp * 0.9f;
        }
        a[i] = s;
    }
    for (int i = 0; i < x; ++i) { const float w = (float)i / x; a[i] = a[i] * w + a[len + i] * (1.0f - w); }
    a.resize(len);
    Normalise(a, 0.7f);
    return a;
}

std::vector<float> AudioSynth_Thunder(float near, unsigned seed) {
    near = std::clamp(near, 0.0f, 1.0f);
    const float seconds = 5.5f;
    const int len = (int)(seconds * kRate);
    std::vector<float> a(len);
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> n(-1.0f, 1.0f), u(0.0f, 1.0f);
    float lp = 0.0f, lp2 = 0.0f, crackLp = 0.0f, swell = 0.0f, nextSwell = 0.2f;
    for (int i = 0; i < len; ++i) {
        const float t = (float)i / kRate;
        if (t >= nextSwell) { swell = 0.5f + 0.5f * u(rng); nextSwell = t + 0.3f + 0.7f * u(rng); }
        swell *= 0.99995f;
        const float x0 = n(rng);
        lp += (x0 - lp) * (0.01f + 0.03f * near);
        lp2 += (lp - lp2) * 0.05f;
        const float env = std::min(1.0f, t / (0.35f - 0.25f * near)) * expf(-t * 0.6f);
        float s = lp2 * 14.0f * env * (0.5f + swell);
        if (t < 0.25f && near > 0.0f) { // the crack of a near strike
            crackLp += (x0 - crackLp) * 0.6f;
            s += near * crackLp * 1.5f * expf(-t * 18.0f);
        }
        a[i] = s;
    }
    Normalise(a, 0.9f);
    return a;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1241): the beds are seamless, the same each start, and every named one exists ----
#include "unit_helpers.h"

TEST(AudioSynth, BedsLoopSeamlesslyAndRepeat) {
    AudioBedRecipe r;
    ASSERT_TRUE(AudioSynth_BedRecipe("bed_forest_day", r));
    const std::vector<float> a = AudioSynth_Bed(r, 4.0f, 7u), b = AudioSynth_Bed(r, 4.0f, 7u);
    ASSERT_EQ(a.size(), (size_t)(4.0f * 44100));
    EXPECT_EQ(a, b);                                              // deterministic
    float peak = 0.0f, meanStep = 0.0f;
    for (size_t i = 0; i < a.size(); ++i) { peak = std::max(peak, fabsf(a[i])); if (i) meanStep += fabsf(a[i] - a[i - 1]); }
    meanStep /= (float)(a.size() - 1);
    EXPECT_NEAR(peak, 0.7f, 1e-4f);
    EXPECT_LT(fabsf(a.front() - a.back()), 0.2f + 8.0f * meanStep); // the wrap is no bigger a jump than the sound makes
    for (const char* name : { "bed_forest_night", "bed_jungle_day", "bed_jungle_night", "bed_swamp_day", "bed_swamp_night",
                              "bed_savanna_day", "bed_savanna_night", "bed_desert", "bed_desert_night", "bed_taiga",
                              "bed_taiga_night", "bed_cold", "bed_shore", "bed_rain", "bed_rain_leaves", "bed_meadow_day" })
        EXPECT_TRUE(AudioSynth_BedRecipe(name, r)) << name;
    EXPECT_FALSE(AudioSynth_BedRecipe("bed_nowhere", r));
}

TEST(AudioSynth, ThunderIsBounded) {
    const std::vector<float> t = AudioSynth_Thunder(1.0f, 3u);
    float peak = 0.0f;
    for (float v : t) peak = std::max(peak, fabsf(v));
    EXPECT_NEAR(peak, 0.9f, 1e-4f);
    EXPECT_GT(t.size(), (size_t)44100 * 4);
}
#endif
