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
    else if (name == "bed_waterfall") { r.rain = 1.0f; r.surf = 0.8f; } // (#1242) a roar: dense hiss over heavy swells
    else if (name == "bed_rumble") { r.rumble = 1.0f; }                 // (#1242) deep under a volcano
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
    float rbBrown = 0.0f, rbLp = 0.0f;
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
        if (r.rumble > 0.0f) { // brown noise far down, heaving slowly, with a 31 Hz drone
            rbBrown = rbBrown * 0.999f + n(rng) * 0.02f;
            rbLp += (rbBrown - rbLp) * 0.004f;
            s += r.rumble * (rbLp * 30.0f * (0.6f + 0.4f * sinf(kTau * 0.11f * t)) + 0.15f * sinf(kTau * 31.0f * t));
        }
        a[i] = s;
    }
    for (int i = 0; i < x; ++i) { const float w = (float)i / x; a[i] = a[i] * w + a[len + i] * (1.0f - w); }
    a.resize(len);
    Normalise(a, 0.7f);
    return a;
}

std::vector<float> AudioSynth_Drip(unsigned seed) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> n(-1.0f, 1.0f), u(0.0f, 1.0f);
    const int len = (int)(0.18f * kRate);
    std::vector<float> a(len);
    const float f0 = 700.0f + 500.0f * u(rng); // the bubble's ring rises as it closes (a "plink")
    float phase = 0.0f, lp = 0.0f;
    for (int i = 0; i < len; ++i) {
        const float t = (float)i / kRate;
        phase += kTau * f0 * (1.0f + 2.5f * t / 0.18f) / kRate;
        lp += (n(rng) - lp) * 0.4f;
        a[i] = sinf(phase) * expf(-t * 35.0f) * std::min(1.0f, t * 800.0f) + lp * 0.25f * expf(-t * 120.0f);
    }
    Normalise(a, 0.8f);
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

namespace {
    // A seamless loop of `seconds`: f(t) for t over seconds + xfade, the tail blended into the head.
    template <typename F>
    std::vector<float> Loop(float seconds, float xfade, F f) {
        const int len = (int)(seconds * kRate), x = (int)(xfade * kRate);
        std::vector<float> a(len + x);
        for (int i = 0; i < len + x; ++i) a[i] = f((float)i / kRate);
        for (int i = 0; i < x; ++i) { const float w = (float)i / x; a[i] = a[i] * w + a[len + i] * (1.0f - w); }
        a.resize(len);
        return a;
    }
    // A mains-like hum: a fundamental and its harmonics (odd ones stronger), whole cycles in the loop's length.
    float Hum(float t, float f0, int harmonics, float evenShare) {
        float s = 0.0f;
        for (int k = 1; k <= harmonics; ++k) s += sinf(kTau * f0 * k * t) / k * (k % 2 ? 1.0f : evenShare);
        return s;
    }
}

bool AudioSynth_Machine(const std::string& name, unsigned seed, std::vector<float>& out) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> n(-1.0f, 1.0f), u(0.0f, 1.0f);
    float lp = 0.0f, lp2 = 0.0f, pop = 0.0f, grain = 0.0f, clank = 0.0f, clankT = 1.0f, clankF = 0.0f;
    if (name == "hum_fire") { // a roaring fire: low noise, flickering, with pops (campfires, kilns, coke ovens, blast furnaces)
        out = Loop(6.0f, 1.0f, [&](float t) {
            lp += (n(rng) - lp) * 0.03f;
            if (u(rng) < 6.0f / kRate) pop = 0.6f + 0.4f * u(rng);
            pop *= 0.995f;
            return lp * 4.0f * (0.7f + 0.3f * sinf(kTau * 0.7f * t)) + pop * n(rng) * 0.5f;
        });
    } else if (name == "hum_steam") { // a steam engine: a chug 2.5 times a second over a hiss
        out = Loop(4.0f, 0.5f, [&](float t) {
            const float x = n(rng);
            lp += (x - lp) * 0.2f;
            const float ph = fmodf(t * 2.5f, 1.0f);
            const float chug = expf(-ph * 9.0f);
            lp2 += (x - lp2) * 0.05f;
            return (x - lp) * 0.25f + lp2 * 3.0f * chug + 0.2f * sinf(kTau * 45.0f * t) * chug;
        });
    } else if (name == "hum_lv") { // low voltage: a 50 Hz hum and a soft mechanical whir
        out = Loop(2.0f, 0.3f, [&](float t) {
            lp += (n(rng) - lp) * 0.08f;
            return 0.35f * Hum(t, 50.0f, 6, 0.5f) + lp * 0.8f * (0.8f + 0.2f * sinf(kTau * 3.0f * t));
        });
    } else if (name == "hum_mv") { // medium voltage: higher and cleaner, a faint whine
        out = Loop(2.0f, 0.3f, [&](float t) {
            lp += (n(rng) - lp) * 0.1f;
            return 0.3f * Hum(t, 100.0f, 5, 0.6f) + 0.06f * sinf(kTau * 1200.0f * t) + lp * 0.4f;
        });
    } else if (name == "hum_hv") { // high voltage: a buzzing 150 Hz and a whine
        out = Loop(2.0f, 0.3f, [&](float t) {
            const float buzz = Hum(t, 150.0f, 9, 0.2f);
            return 0.25f * buzz + 0.07f * sinf(kTau * 2400.0f * t) * (0.7f + 0.3f * sinf(kTau * 0.5f * t));
        });
    } else if (name == "hum_ev") { // fusion, mass drivers: a deep throb and a shimmer
        out = Loop(4.0f, 0.5f, [&](float t) {
            const float throb = 0.5f + 0.5f * sinf(kTau * 1.5f * t);
            return 0.4f * sinf(kTau * 40.0f * t) * throb + 0.05f * sinf(kTau * 3000.0f * t) * (0.5f + 0.5f * sinf(kTau * 7.0f * t));
        });
    } else if (name == "factory_bed") { // many machines: hums, conveyor rattle and the odd clank, blended
        out = Loop(10.0f, 1.5f, [&](float t) {
            lp += (n(rng) - lp) * 0.06f;
            if (u(rng) < 40.0f / kRate) grain = 0.5f + 0.5f * u(rng);
            grain *= 0.99f;
            if (clankT >= 0.25f && u(rng) < 0.7f / kRate) { clankT = 0.0f; clankF = 400.0f + 600.0f * u(rng); clank = 0.4f + 0.4f * u(rng); }
            float c = 0.0f;
            if (clankT < 0.25f) { c = clank * expf(-clankT * 18.0f) * (sinf(kTau * clankF * clankT) + 0.5f * sinf(kTau * clankF * 2.71f * clankT)); clankT += 1.0f / kRate; }
            return 0.25f * Hum(t, 50.0f, 5, 0.5f) + 0.1f * Hum(t, 100.0f, 3, 0.6f) + lp * 0.9f + grain * n(rng) * 0.2f + c * 0.5f;
        });
    // (#1269, #1270) Mechanical power, made at 64 RPM (machinesounds.cpp plays them faster or slower by the shaft's speed).
    } else if (name == "mech_gears") { // wooden cogs: teeth clacking 8 times a second over a low rumble
        out = Loop(2.0f, 0.25f, [&](float t) {
            const float ph = fmodf(t * 8.0f, 1.0f);
            lp += (n(rng) - lp) * 0.05f;
            lp2 += (n(rng) - lp2) * 0.4f;
            return lp2 * 0.9f * expf(-ph * 40.0f) + 0.15f * sinf(kTau * 180.0f * t) * expf(-ph * 30.0f) + lp * 0.6f;
        });
    } else if (name == "mech_millstone") { // stone grinding on stone, a rumble swelling once a turn
        out = Loop(3.75f, 0.5f, [&](float t) {
            lp += (n(rng) - lp) * 0.04f;
            lp2 += (n(rng) - lp2) * 0.25f;
            const float turn = 0.7f + 0.3f * sinf(kTau * t * (64.0f / 60.0f));
            return (lp * 3.0f + lp2 * 0.3f) * turn;
        });
    } else if (name == "mech_saw") { // a blade cutting: a buzzing saw tooth with chips of noise
        out = Loop(2.0f, 0.25f, [&](float t) {
            const float saw = 2.0f * fmodf(t * 110.0f, 1.0f) - 1.0f;
            lp += (n(rng) - lp) * 0.5f;
            return saw * 0.35f * (0.8f + 0.2f * sinf(kTau * 4.0f * t)) + lp * 0.25f;
        });
    } else if (name == "mech_press") { // a press stamping about once a second: a heavy thump and a clank
        out = Loop(1.875f, 0.0f, [&](float t) {
            const float ph = fmodf(t * (64.0f / 60.0f), 1.0f) / (64.0f / 60.0f); // seconds into this stroke
            lp += (n(rng) - lp) * 0.1f;
            return sinf(kTau * 55.0f * ph) * expf(-ph * 14.0f) * 0.9f + sinf(kTau * 700.0f * ph) * expf(-ph * 30.0f) * 0.2f + lp * 0.2f * expf(-ph * 10.0f);
        });
    } else if (name == "mech_waterwheel") { // paddles slapping the water and a slow creak
        out = Loop(3.75f, 0.5f, [&](float t) {
            const float paddle = fmodf(t * 8.0f * (64.0f / 60.0f) / 4.0f, 1.0f); // a paddle every ~0.47 s
            lp += (n(rng) - lp) * 0.3f;
            const float creak = sinf(kTau * (220.0f + 40.0f * sinf(kTau * 3.0f * t)) * t) * 0.08f * (0.5f + 0.5f * sinf(kTau * t * (64.0f / 60.0f)));
            return lp * 1.2f * expf(-paddle * 6.0f) + lp * 0.25f + creak;
        });
    } else if (name == "mech_windmill") { // sails sweeping past (four a turn) and the frame's creak
        out = Loop(3.75f, 0.5f, [&](float t) {
            const float sail = 0.5f + 0.5f * cosf(kTau * t * 4.0f * (64.0f / 60.0f) / 4.0f);
            lp += (n(rng) - lp) * 0.06f;
            const float creak = sinf(kTau * 160.0f * t) * 0.06f * (0.5f + 0.5f * sinf(kTau * t * (64.0f / 60.0f)));
            return lp * 2.5f * sail * sail + creak;
        });
    } else if (name == "machine_stall") { // winding down: a falling tone and a click
        const int len = (int)(0.6f * kRate);
        out.resize(len);
        float ph = 0.0f;
        for (int i = 0; i < len; ++i) {
            const float t = (float)i / kRate;
            ph += kTau * (600.0f - 450.0f * (t / 0.6f)) / kRate;
            out[i] = sinf(ph) * 0.5f * (1.0f - t / 0.6f) + (t < 0.02f ? n(rng) * 0.6f * (1.0f - t / 0.02f) : 0.0f);
        }
    } else if (name == "machine_click") { // a relay's click
        const int len = (int)(0.04f * kRate);
        out.resize(len);
        for (int i = 0; i < len; ++i) { const float t = (float)i / kRate; lp += (n(rng) - lp) * 0.5f; out[i] = (n(rng) - lp) * expf(-t * 150.0f) + 0.3f * sinf(kTau * 2000.0f * t) * expf(-t * 200.0f); }
    } else if (name == "machine_done") { // a batch done: two soft bell notes
        const int len = (int)(0.7f * kRate);
        out.resize(len);
        for (int i = 0; i < len; ++i) {
            const float t = (float)i / kRate, t2 = t - 0.12f;
            float s = sinf(kTau * 880.0f * t) * expf(-t * 7.0f) * std::min(1.0f, t * 400.0f);
            if (t2 > 0.0f) s += sinf(kTau * 1320.0f * t2) * expf(-t2 * 6.0f) * std::min(1.0f, t2 * 400.0f);
            out[i] = s * 0.4f;
        }
    } else return false;
    Normalise(out, 0.7f);
    return true;
}

namespace {
    float Midi(float m) { return 440.0f * powf(2.0f, (m - 69.0f) / 12.0f); }
    // sin(2 pi cycles) from a table, the phase kept in double (a minute of a 2 kHz tone is 120 000 cycles): music is
    // tens of sines a sample over tens of seconds, too slow with sinf when it is made on first use.
    float SinCycles(double cycles) {
        static float table[4097];
        static bool made = false;
        if (!made) { for (int i = 0; i <= 4096; ++i) table[i] = sinf(kTau * i / 4096.0f); made = true; }
        const double f = cycles - floor(cycles);
        const float x = (float)(f * 4096.0);
        const int i = (int)x;
        return table[i] + (table[i + 1] - table[i]) * (x - i);
    }
    // One sustained tone of a pad: a soft organ-like stack and a slightly detuned twin (a slow chorus).
    float PadTone(float f, float t, float bright) {
        const double c = (double)f * t;
        return SinCycles(c) + bright * (0.3f * SinCycles(2.0 * c) + 0.12f * SinCycles(3.0 * c)) + 0.6f * SinCycles(c * 1.004);
    }
    // An echo: y[i] += fb * y[i - delay], in place (generated, then looped: the tail carries into the head by the crossfade).
    void Echo(std::vector<float>& a, float seconds, float fb) {
        const int d = (int)(seconds * kRate);
        for (size_t i = (size_t)d; i < a.size(); ++i) a[i] += fb * a[i - d];
    }
    // Crossfade the last `x` samples into the head and keep `len` (a seamless loop).
    void CloseLoop(std::vector<float>& a, int len) {
        const int x = (int)a.size() - len;
        for (int i = 0; i < x; ++i) { const float w = (float)i / x; a[i] = a[i] * w + a[len + i] * (1.0f - w); }
        a.resize(len);
    }
    struct Pluck { float t0, f, amp; };
}

bool AudioSynth_Music(const std::string& name, unsigned seed, std::vector<float>& out) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> u(0.0f, 1.0f);
    // Calm pads: four chords of `seg` seconds each, a plucked melody over them drawn from a scale, an echo.
    auto pad = [&](const std::vector<std::vector<float>>& chords, float seg, const std::vector<float>& scale, float pluckRate,
                   float bright, float echoS, float echoFb, float drone) {
        const float seconds = seg * chords.size(), xfade = 2.0f;
        const int len = (int)(seconds * kRate), total = (int)((seconds + xfade) * kRate);
        out.assign(total, 0.0f);
        std::vector<Pluck> plucks;
        for (float t = 0.5f; t < seconds + xfade; t += (0.5f + u(rng)) / std::max(0.05f, pluckRate))
            plucks.push_back({ t, Midi(scale[(size_t)(u(rng) * scale.size()) % scale.size()]), 0.25f + 0.2f * u(rng) });
        size_t firstPluck = 0; // (plucks are in time order: only those started within 3 s sound)
        for (int i = 0; i < total; ++i) {
            const float t = (float)i / kRate;
            while (firstPluck < plucks.size() && t - plucks[firstPluck].t0 > 3.0f) ++firstPluck;
            const int c = (int)(t / seg) % (int)chords.size(), prev = (c + (int)chords.size() - 1) % (int)chords.size();
            const float into = fmodf(t, seg);
            const float w = std::min(1.0f, into / 1.5f); // chords glide into each other over 1.5 s
            float a = 0.0f, b = 0.0f;
            for (float m : chords[c]) a += PadTone(Midi(m), t, bright);
            for (float m : chords[prev]) b += PadTone(Midi(m), t, bright);
            float s = (b * (1.0f - w) + a * w) * 0.12f * (0.85f + 0.15f * sinf(kTau * 0.07f * t));
            if (drone > 0.0f) s += drone * (SinCycles((double)Midi(38) * t) + 0.5f * SinCycles((double)Midi(45) * t)) * 0.2f;
            for (size_t k = firstPluck; k < plucks.size() && plucks[k].t0 <= t; ++k) {
                const Pluck& p = plucks[k];
                const float tau = t - p.t0;
                s += p.amp * SinCycles((double)p.f * tau) * expf(-tau * 2.5f) * std::min(1.0f, tau * 300.0f);
            }
            out[i] = s;
        }
        if (echoFb > 0.0f) Echo(out, echoS, echoFb);
        CloseLoop(out, len);
    };
    if (name == "music_pad_day") // D major, warm: D, Bm, G, A; a pentatonic melody
        pad({ { 50, 57, 62, 66 }, { 47, 54, 59, 62 }, { 43, 55, 59, 62 }, { 45, 52, 57, 61 } }, 5.0f, { 74, 76, 78, 81, 83, 86 }, 0.9f, 0.8f, 0.3f, 0.25f, 0.0f);
    else if (name == "music_pad_night") // A minor, soft: Am, F, C, G; a sparse melody with a longer echo
        pad({ { 45, 52, 57, 60 }, { 41, 53, 57, 60 }, { 48, 55, 60, 64 }, { 43, 55, 59, 62 } }, 5.0f, { 69, 72, 74, 76, 79, 81 }, 0.45f, 0.4f, 0.45f, 0.4f, 0.0f);
    else if (name == "music_pad_cave") // a low drone, faint chords, high bells far apart with long echoes
        pad({ { 62, 69 }, { 60, 67 }, { 62, 69 }, { 57, 64 } }, 5.0f, { 74, 77, 79, 81, 84, 86 }, 0.3f, 0.2f, 0.6f, 0.55f, 1.0f);
    else if (name == "music_pad_space") // two slow, open clusters, shimmering
        pad({ { 62, 69, 76, 81 }, { 60, 67, 74, 79 } }, 10.0f, { 86, 88, 93, 95 }, 0.2f, 0.25f, 0.7f, 0.5f, 0.0f);
    else if (name == "music_danger") { // 100 BPM: a low pulse and a tense minor second, tremolo
        const int len = 12 * kRate;
        out.assign(len, 0.0f);
        for (int i = 0; i < len; ++i) {
            const float t = (float)i / kRate, beat = fmodf(t, 0.6f);
            const float kick = sinf(kTau * (55.0f + 40.0f * expf(-beat * 30.0f)) * beat) * expf(-beat * 8.0f);
            const float drone = (sinf(kTau * Midi(40) * t) + sinf(kTau * Midi(41) * t)) * 0.25f * (0.7f + 0.3f * sinf(kTau * 5.0f * t));
            out[i] = kick * 0.8f + drone;
        }
    } else if (name == "music_discovery" || name == "music_achievement") { // stings: bells rising, over a chord
        const bool disc = name == "music_discovery";
        const float seconds = disc ? 6.0f : 4.0f;
        const int len = (int)(seconds * kRate);
        out.assign(len, 0.0f);
        const std::vector<float> notes = disc ? std::vector<float>{ 74, 78, 81, 86 } : std::vector<float>{ 72, 76, 79, 84, 88 };
        const std::vector<float> chord = disc ? std::vector<float>{ 50, 57, 62, 66 } : std::vector<float>{ 48, 55, 60, 64 };
        const float step = disc ? 0.4f : 0.12f;
        for (int i = 0; i < len; ++i) {
            const float t = (float)i / kRate;
            float s = 0.0f;
            for (size_t k = 0; k < notes.size(); ++k) {
                const float tau = t - step * k;
                if (tau > 0.0f) s += 0.3f * sinf(kTau * Midi(notes[k]) * tau) * expf(-tau * (disc ? 1.5f : 2.5f)) * std::min(1.0f, tau * 300.0f);
            }
            const float env = std::min(1.0f, t / 0.5f) * std::max(0.0f, 1.0f - t / seconds);
            for (float m : chord) s += PadTone(Midi(m), t, 0.5f) * 0.05f * env;
            out[i] = s;
        }
        Echo(out, 0.35f, 0.3f);
    } else return false;
    Normalise(out, 0.6f);
    return true;
}

bool AudioSynth_Creature(const std::string& event, const std::string& voice, unsigned seed, std::vector<float>& out) {
    // The event shapes the call: how long, the pitch from start to end, the loudness over time.
    struct Shape { float seconds, p0, p1, attack, decay; bool crescendo; };
    Shape sh;
    if (event == "idle") sh = { 0.45f, 1.0f, 0.95f, 0.05f, 6.0f, false };
    else if (event == "alert") sh = { 0.6f, 1.1f, 1.0f, 0.03f, 4.0f, false };
    else if (event == "telegraph") sh = { 0.5f, 0.9f, 1.35f, 0.0f, 0.0f, true }; // a wind-up: swelling and rising
    else if (event == "hurt") sh = { 0.3f, 1.6f, 1.1f, 0.01f, 9.0f, false };
    else if (event == "death") sh = { 1.1f, 1.0f, 0.55f, 0.05f, 2.2f, false };
    else return false;
    enum V { Beast, Lurker, Scavenger, Ooze, Undead, Construct } v;
    if (voice == "beast") v = Beast;
    else if (voice == "lurker") v = Lurker;
    else if (voice == "scavenger") v = Scavenger;
    else if (voice == "ooze") v = Ooze;
    else if (voice == "undead") v = Undead;
    else if (voice == "construct") v = Construct;
    else return false;
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> n(-1.0f, 1.0f), u(0.0f, 1.0f);
    const int len = (int)(sh.seconds * kRate);
    out.assign(len, 0.0f);
    double ph = 0.0, ph2 = 0.0;
    float lp = 0.0f, lp2 = 0.0f, click = 0.0f, bub = 0.0f, bubF = 0.0f, bubT = 1.0f;
    for (int i = 0; i < len; ++i) {
        const float t = (float)i / kRate, x = t / sh.seconds;
        const float p = sh.p0 + (sh.p1 - sh.p0) * x;
        const float env = sh.crescendo ? std::min(1.0f, 0.2f + x * 1.2f) * std::min(1.0f, (1.0f - x) * 25.0f)
                                       : std::min(1.0f, sh.attack > 0.0f ? t / sh.attack : 1.0f) * expf(-t * sh.decay);
        float s = 0.0f;
        switch (v) {
            case Beast: { // a growl: a rough saw tooth near 110 Hz, low-passed, rumbling
                const float f = 110.0f * p * (1.0f + 0.03f * n(rng));
                ph += f / kRate;
                const float saw = 2.0f * (float)(ph - floor(ph)) - 1.0f;
                lp += (saw + 0.3f * n(rng) - lp) * 0.18f;
                s = lp * (1.0f + 0.4f * sinf(kTau * 32.0f * t));
                break;
            }
            case Lurker: { // a rattle: a train of clicks, faster as it winds up, over a hiss
                const float rate = 30.0f * p * (sh.crescendo ? 1.0f + 1.5f * x : 1.0f);
                ph += rate / kRate;
                if (ph >= 1.0) { ph -= 1.0; click = 1.0f; }
                click *= 0.93f;
                lp += (n(rng) - lp) * 0.6f;
                s = click * n(rng) * 1.2f + (n(rng) - lp) * 0.15f;
                break;
            }
            case Scavenger: { // a chittering screech: a frequency-modulated whistle near 1.1 kHz
                const float f = 1100.0f * p;
                ph2 += 37.0f / kRate;
                ph += (f * (1.0f + 0.15f * SinCycles(ph2))) / kRate;
                s = SinCycles(ph) * 0.8f + 0.2f * n(rng) * (0.5f + 0.5f * SinCycles(ph2 * 0.5));
                break;
            }
            case Ooze: { // a gurgle: low bubbles over a wet, low-passed noise
                if (u(rng) < (sh.crescendo ? 40.0f + 120.0f * x : 50.0f) / kRate) { bubF = (120.0f + 220.0f * u(rng)) * p; bubT = 0.0f; }
                if (bubT < 0.06f) { ph += bubF * (1.0f + 4.0f * bubT) / kRate; bub = (float)SinCycles(ph) * (1.0f - bubT / 0.06f); bubT += 1.0f / kRate; } else bub = 0.0f;
                lp += (n(rng) - lp) * 0.05f;
                s = bub * 0.8f + lp * 1.5f;
                break;
            }
            case Undead: { // a moan: a voiced tone with vibrato and breath
                const float f = 130.0f * p * (1.0f + 0.02f * sinf(kTau * 5.0f * t));
                ph += f / kRate;
                lp2 += (n(rng) - lp2) * 0.3f;
                s = (float)(SinCycles(ph) + 0.5 * SinCycles(2.0 * ph) + 0.25 * SinCycles(3.0 * ph)) * 0.6f + (n(rng) - lp2) * 0.3f;
                break;
            }
            case Construct: { // tones: a square-ish buzz and an octave over it, charging up as it winds up
                const float f = 330.0f * p * (sh.crescendo ? 1.0f + 2.0f * x * x : 1.0f);
                ph += f / kRate;
                const float sq = (ph - floor(ph)) < 0.5 ? 0.5f : -0.5f;
                lp += (sq - lp) * 0.3f;
                s = lp + 0.3f * (float)SinCycles(2.0 * ph);
                break;
            }
        }
        out[i] = s * env;
    }
    Normalise(out, 0.8f);
    return true;
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

TEST(AudioSynth, MusicPadsLoopAndStingsAreBounded) {
    std::vector<float> a, b;
    ASSERT_TRUE(AudioSynth_Music("music_pad_day", 5u, a));
    ASSERT_TRUE(AudioSynth_Music("music_pad_day", 5u, b));
    EXPECT_EQ(a, b);
    EXPECT_EQ(a.size(), (size_t)(20 * 44100));
    float peak = 0.0f, meanStep = 0.0f;
    for (size_t i = 0; i < a.size(); ++i) { peak = std::max(peak, fabsf(a[i])); if (i) meanStep += fabsf(a[i] - a[i - 1]); }
    meanStep /= (float)(a.size() - 1);
    EXPECT_NEAR(peak, 0.6f, 1e-4f);
    EXPECT_LT(fabsf(a.front() - a.back()), 0.1f + 8.0f * meanStep); // seamless
    for (const char* n : { "music_pad_night", "music_pad_cave", "music_pad_space", "music_danger", "music_discovery", "music_achievement" }) {
        std::vector<float> m;
        ASSERT_TRUE(AudioSynth_Music(n, 1u, m)) << n;
        float p = 0.0f;
        for (float v : m) p = std::max(p, fabsf(v));
        EXPECT_NEAR(p, 0.6f, 1e-4f) << n;
    }
    std::vector<float> none;
    EXPECT_FALSE(AudioSynth_Music("music_unknown", 1u, none));
}

TEST(AudioSynth, CreatureVoicesForEveryEvent) {
    for (const char* voice : { "beast", "lurker", "scavenger", "ooze", "undead", "construct" })
        for (const char* ev : { "idle", "alert", "telegraph", "hurt", "death" }) {
            std::vector<float> s;
            ASSERT_TRUE(AudioSynth_Creature(ev, voice, 3u, s)) << ev << " " << voice;
            float peak = 0.0f;
            for (float v : s) peak = std::max(peak, fabsf(v));
            EXPECT_NEAR(peak, 0.8f, 1e-4f) << ev << " " << voice;
            if (std::string(ev) == "telegraph") {
                EXPECT_EQ(s.size(), (size_t)(0.5f * 44100)) << voice; // the wind-up fits a 0.5 s telegraph
                // It swells: the second half is louder than the first.
                float a = 0.0f, b = 0.0f;
                for (size_t i = 0; i < s.size() / 2; ++i) a += s[i] * s[i];
                for (size_t i = s.size() / 2; i < s.size(); ++i) b += s[i] * s[i];
                EXPECT_GT(b, a) << voice;
            }
        }
    std::vector<float> s;
    EXPECT_FALSE(AudioSynth_Creature("sing", "beast", 1u, s));
    EXPECT_FALSE(AudioSynth_Creature("idle", "dragon", 1u, s));
}
#endif
