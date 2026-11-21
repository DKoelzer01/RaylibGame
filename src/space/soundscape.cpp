#include "soundscape.h"

#include "audio.h"
#include "player.h"
#include "profiler.h"
#include "scene.h"
#include "weather.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace {
    const char* kBiomeKeys[(int)Biome::Count] = { "ocean", "jungle", "savanna", "desert", "badlands", "forest", "meadow",
                                                  "swamp", "taiga", "tundra", "alpine", "polar" };
    float Smooth(float a, float b, float x) { const float t = std::clamp((x - a) / (b - a), 0.0f, 1.0f); return t * t * (3.0f - 2.0f * t); }
}

bool Ambience_Parse(const std::string& text, AmbienceTable& t, std::string* error) {
    std::istringstream in(text);
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.resize(hash);
        std::istringstream ls(line);
        std::string key, a, b;
        float vol = 0.0f;
        if (!(ls >> key)) continue; // blank or comment
        auto fail = [&](const char* why) { if (error) *error = "line " + std::to_string(lineNo) + ": " + why; return false; };
        if (key[0] == '@') {
            if (!(ls >> a)) return fail("a layer needs a sound name");
            if (key == "@rain") t.rain = a;
            else if (key == "@leaves") t.leaves = a;
            else if (key == "@wind") t.wind = a;
            else return fail("unknown layer (@rain, @leaves, @wind)");
            continue;
        }
        int b0 = -1;
        for (int i = 0; i < (int)Biome::Count; ++i) if (key == kBiomeKeys[i]) b0 = i;
        if (b0 < 0) return fail("unknown biome");
        if (!(ls >> a >> b >> vol)) return fail("expected: biome day night volume");
        t.biome[b0] = { a, b, std::clamp(vol, 0.0f, 2.0f) };
    }
    return true;
}

void Ambience_Targets(const AmbienceTable& t, const AmbienceInputs& in, std::vector<std::pair<std::string, float>>& out) {
    out.clear();
    auto add = [&](const std::string& name, float v) {
        if (name.empty()) return;
        for (auto& kv : out) if (kv.first == name) { kv.second += v; return; }
        out.push_back({ name, v });
    };
    const float day = Smooth(-8.0f, 4.0f, in.sunElevation);
    const float under = Smooth(3.0f, 15.0f, -in.heightAboveGround);   // below the surface (caves: their own sounds, #1242)
    const float outside = 1.0f - under;
    const float lowEnough = 1.0f - Smooth(60.0f, 400.0f, in.heightAboveGround); // the ground's sounds thin out high up
    const float air = std::clamp(in.air, 0.0f, 1.0f);
    const float wet = in.underwater ? 0.15f : 1.0f;
    const float bedGain = outside * lowEnough * air * wet * (1.0f - 0.5f * std::clamp(in.snow, 0.0f, 1.0f)); // snow hushes
    for (int b = 0; b < (int)Biome::Count; ++b) {
        const float w = in.biomeWeight[b] * t.biome[b].volume * bedGain;
        add(t.biome[b].day, w * day);
        add(t.biome[b].night, w * (1.0f - day));
    }
    const float leafy = in.biomeWeight[(int)Biome::Forest] + in.biomeWeight[(int)Biome::Jungle] +
                        in.biomeWeight[(int)Biome::Taiga] + in.biomeWeight[(int)Biome::Swamp];
    const float rain = std::clamp(in.rain, 0.0f, 1.0f) * outside * air * wet * 0.8f;
    add(t.rain, rain * (1.0f - 0.6f * leafy));
    add(t.leaves, rain * leafy);
    const float wind = air * outside * wet * std::clamp(0.15f + in.wind / 15.0f, 0.0f, 1.0f) * (0.35f + 0.65f * Smooth(30.0f, 300.0f, in.heightAboveGround));
    add(t.wind, wind * 0.6f);
}

float Ambience_Muffle(const AmbienceInputs& in) {
    if (in.underwater) return 450.0f;
    if (in.air >= 0.5f) return 0.0f;
    const float a = std::max(0.0f, in.air) / 0.5f; // 0 vacuum .. 1
    return 500.0f + a * a * 8000.0f;              // through the suit: dull in vacuum, opening up as the air thickens
}

namespace {
    AmbienceTable table;
    bool loaded = false;
    struct Bed { int emitter = -1; float volume = 0.0f; bool tried = false; };
    std::unordered_map<std::string, Bed> beds;
    AmbienceInputs inputs;
    double nextBiomeSample = 0.0;
    std::vector<std::pair<std::string, float>> targets;

    void Load() {
        loaded = true;
        std::ifstream f("resources/data/ambience.txt");
        if (!f) { fprintf(stderr, "soundscape: resources/data/ambience.txt is missing: no ambience beds\n"); return; }
        std::stringstream ss;
        ss << f.rdbuf();
        std::string err;
        if (!Ambience_Parse(ss.str(), table, &err)) fprintf(stderr, "soundscape: ambience.txt %s\n", err.c_str());
    }

    // Where the listener is: the biomes around it, its height over the ground, the air.
    void Gather(Scene& world) {
        const DVec3 cam = world.cameraWorld();
        const Vector3 camR = world.camera.position;
        // Air: the start planet's atmosphere (Scene::atmosphere); anywhere else is vacuum.
        if (world.atmosphere.enabled) {
            const float alt = Vector3Distance(camR, world.atmosphere.centre) - world.atmosphere.radius;
            inputs.air = std::clamp(1.0f - alt / std::max(1.0f, world.atmosphere.height), 0.0f, 1.0f);
        } else inputs.air = 0.0f;
        inputs.underwater = IsUnderwater(world, camR);
        inputs.sunElevation = SunElevationAtCamera(world);
        const WeatherState& w = Weather_Here();
        inputs.snow = w.kind == WeatherKind::Snow ? w.precipitation : 0.0f;
        inputs.rain = w.kind == WeatherKind::Rain || w.kind == WeatherKind::Storm ? w.precipitation : 0.0f;
        inputs.wind = Vector3Length(w.wind);
        Planetoid* p = player.nearest;
        const PlanetField* pf = p ? dynamic_cast<const PlanetField*>(p->field.get()) : nullptr;
        if (!pf) {
            for (float& b : inputs.biomeWeight) b = 0.0f;
            inputs.heightAboveGround = 1e6f;
            return;
        }
        const Vector3 local = p->toLocalD(cam).toFloat();
        const float r = Vector3Length(local);
        if (r < 1.0f) return;
        const Vector3 dir = Vector3Scale(local, 1.0f / r);
        inputs.heightAboveGround = r - (pf->radius + pf->elevation(dir));
        if (GetTime() < nextBiomeSample) return;
        nextBiomeSample = GetTime() + 0.5;
        // The biome here and at six points 80 m around: neighbouring biomes blend in near a border (a shore's surf too).
        const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(dir, fabsf(dir.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 t2 = Vector3CrossProduct(dir, t1);
        float wsum[(int)Biome::Count] = {};
        for (int k = 0; k < 7; ++k) {
            Vector3 d = dir;
            if (k > 0) {
                const float a = 6.2831853f * (k - 1) / 6.0f;
                d = Vector3Normalize(Vector3Add(local, Vector3Scale(Vector3Add(Vector3Scale(t1, cosf(a)), Vector3Scale(t2, sinf(a))), 80.0f)));
            }
            const int b = (int)pf->biomeAt(d, pf->elevation(d));
            if (b >= 0 && b < (int)Biome::Count) wsum[b] += k == 0 ? 2.0f : 1.0f; // (where you stand counts double)
        }
        for (int b = 0; b < (int)Biome::Count; ++b) inputs.biomeWeight[b] = wsum[b] / 8.0f;
    }
}

void Soundscape_Update(Scene& world, float dt) {
    if (!Audio_V2()) return;
    if (!loaded) Load();
    PROFILE_SCOPE("audio:soundscape");
    Gather(world);
    Ambience_Targets(table, inputs, targets);
    Audio_SetMuffle(Ambience_Muffle(inputs));
    bool madeOne = false; // at most one new bed a frame (a procedural one is synthesised on first use)
    for (const auto& [name, target] : targets) {
        Bed& b = beds[name];
        if (b.emitter < 0) {
            if (b.tried || target < 0.003f || madeOne) continue;
            b.tried = true;
            madeOne = true;
            const int snd = Audio_Sound(name.c_str());
            if (snd >= 0) b.emitter = Audio_EmitterCreate(snd, AudioBus::Ambience, AudioPriority::Normal, {}, false);
            if (b.emitter < 0) continue;
        }
        b.volume += (target - b.volume) * std::min(1.0f, dt / 1.5f); // beds glide (a biome border, dusk)
        Audio_EmitterSetVolume(b.emitter, b.volume);
        Audio_EmitterSetPlaying(b.emitter, b.volume > 0.003f);
    }
}

float Soundscape_Air() { return inputs.air; }
const AmbienceInputs& Soundscape_Inputs() { return inputs; }
float Soundscape_HeightAboveGround() { return inputs.heightAboveGround; }

void Soundscape_PreloadStep() {
    if (!Audio_V2()) return;
    if (!loaded) Load();
    // Every sound the table names, then thunder: the first not yet looked up (Audio_Sound caches it, a file or synthesised).
    static std::vector<std::string> names;
    static size_t next = 0;
    if (names.empty()) {
        for (const AmbienceTable::Bed& b : table.biome) { names.push_back(b.day); names.push_back(b.night); }
        names.insert(names.end(), { table.rain, table.leaves, table.wind, "thunder_near", "thunder_far",
                                    "water", "bed_waterfall", "bed_rumble", "drip", // (#1242 water and caves too)
                                    "music_pad_day", "music_pad_night", "music_pad_cave", "music_pad_space", "music_danger",
                                    "music_discovery", "music_achievement" }); // (#1245 the procedural music)
    }
    while (next < names.size()) {
        const std::string& n = names[next++];
        if (n.empty()) continue;
        bool seen = false;
        for (size_t k = 0; k + 1 < next && !seen; ++k) seen = names[k] == n;
        if (seen) continue;
        PROFILE_SCOPE("audio:preload bed");
        Audio_Sound(n.c_str());
        return; // one a frame
    }
}

void Soundscape_Thunder(Scene& world, float delay) {
    const float dist = std::max(0.0f, delay) * 343.0f;
    const bool near = dist < 450.0f;
    static int nearId = -2, farId = -2;
    if (nearId == -2) { nearId = Audio_Sound("thunder_near"); farId = Audio_Sound("thunder_far"); }
    const int snd = near ? nearId : farId;
    if (snd < 0) return;
    // From the strike's direction, as a sound 60 m away that does not fade with that distance (it is kilometres off).
    const Vector3 dir = Vector3LengthSqr(world.lightningDir) > 0.5f ? world.lightningDir : Vector3{ 0, 1, 0 };
    AudioAttenuation att;
    att.refDistance = 100.0f;
    att.maxDistance = 1000.0f;
    att.lowpassDistance = 1000.0f;
    const float volume = std::clamp(1.2f - dist / 1500.0f, 0.3f, 1.0f) * (inputs.underwater ? 0.3f : 1.0f);
    Audio_PlaySoundAt(snd, Audio_GetListener().pos + DVec3(Vector3Scale(dir, 60.0f)), AudioBus::Ambience, AudioPriority::High, volume,
                      near ? 1.0f : 0.85f, att);
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1241): the table, day and night, underground, vacuum, rain on leaves, the muffle ----
#include "unit_helpers.h"

namespace {
float TargetOf(const std::vector<std::pair<std::string, float>>& v, const char* name) {
    for (const auto& kv : v) if (kv.first == name) return kv.second;
    return 0.0f;
}
}

TEST(Soundscape, ShippedTableParsesAndNamesEveryBiome) {
    std::ifstream f("resources/data/ambience.txt");
    ASSERT_TRUE(f.good());
    std::stringstream ss;
    ss << f.rdbuf();
    AmbienceTable t;
    std::string err;
    ASSERT_TRUE(Ambience_Parse(ss.str(), t, &err)) << err;
    for (int b = 0; b < (int)Biome::Count; ++b) {
        EXPECT_FALSE(t.biome[b].day.empty()) << b;
        EXPECT_GT(t.biome[b].volume, 0.0f) << b;
    }
    AmbienceTable bad;
    EXPECT_FALSE(Ambience_Parse("forest bed_a\n", bad, &err));
    EXPECT_FALSE(Ambience_Parse("narnia a b 1\n", bad, &err));
    EXPECT_FALSE(Ambience_Parse("@thunder x\n", bad, &err));
}

TEST(Soundscape, DayNightUndergroundAndVacuum) {
    AmbienceTable t;
    ASSERT_TRUE(Ambience_Parse("forest fday fnight 1\nmeadow mday fnight 0.5\n@rain r\n@leaves l\n@wind w\n", t));
    AmbienceInputs in;
    in.biomeWeight[(int)Biome::Forest] = 1.0f;
    std::vector<std::pair<std::string, float>> out;
    in.sunElevation = 40.0f;
    Ambience_Targets(t, in, out);
    EXPECT_FLOAT_EQ(TargetOf(out, "fday"), 1.0f);
    EXPECT_FLOAT_EQ(TargetOf(out, "fnight"), 0.0f);
    in.sunElevation = -30.0f;
    Ambience_Targets(t, in, out);
    EXPECT_FLOAT_EQ(TargetOf(out, "fday"), 0.0f);
    EXPECT_FLOAT_EQ(TargetOf(out, "fnight"), 1.0f);
    // A shared night bed sums its biomes' shares.
    in.biomeWeight[(int)Biome::Forest] = 0.5f;
    in.biomeWeight[(int)Biome::Meadow] = 0.5f;
    Ambience_Targets(t, in, out);
    EXPECT_FLOAT_EQ(TargetOf(out, "fnight"), 0.75f);
    // 20 m underground: silent beds, no wind, no rain.
    in.heightAboveGround = -20.0f;
    in.rain = 1.0f;
    Ambience_Targets(t, in, out);
    for (const auto& kv : out) EXPECT_FLOAT_EQ(kv.second, 0.0f) << kv.first;
    // Vacuum: nothing, and muffled through the suit.
    in.heightAboveGround = 0.0f;
    in.air = 0.0f;
    Ambience_Targets(t, in, out);
    for (const auto& kv : out) EXPECT_FLOAT_EQ(kv.second, 0.0f) << kv.first;
    EXPECT_FLOAT_EQ(Ambience_Muffle(in), 500.0f);
    in.air = 1.0f;
    EXPECT_FLOAT_EQ(Ambience_Muffle(in), 0.0f);
    in.underwater = true;
    EXPECT_GT(Ambience_Muffle(in), 0.0f);
}

TEST(Soundscape, RainOnLeavesInForestsAndWindWithHeight) {
    AmbienceTable t;
    ASSERT_TRUE(Ambience_Parse("forest a b 1\ndesert c d 1\n@rain r\n@leaves l\n@wind w\n", t));
    AmbienceInputs in;
    in.sunElevation = 30.0f;
    in.rain = 1.0f;
    std::vector<std::pair<std::string, float>> out;
    in.biomeWeight[(int)Biome::Forest] = 1.0f;
    Ambience_Targets(t, in, out);
    EXPECT_GT(TargetOf(out, "l"), TargetOf(out, "r"));
    in.biomeWeight[(int)Biome::Forest] = 0.0f;
    in.biomeWeight[(int)Biome::Desert] = 1.0f;
    Ambience_Targets(t, in, out);
    EXPECT_FLOAT_EQ(TargetOf(out, "l"), 0.0f);
    EXPECT_GT(TargetOf(out, "r"), 0.0f);
    const float windLow = TargetOf(out, "w");
    in.heightAboveGround = 500.0f; // on a peak or flying: the beds thin out, the wind rises
    Ambience_Targets(t, in, out);
    EXPECT_GT(TargetOf(out, "w"), windLow);
    EXPECT_FLOAT_EQ(TargetOf(out, "c"), 0.0f);
}
#endif
