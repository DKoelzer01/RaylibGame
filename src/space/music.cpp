#include "music.h"

#include "audio.h"
#include "profiler.h"
#include "soundscape.h"

#include <raylib.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <random>
#include <sstream>
#include <unordered_set>

bool Music_Parse(const std::string& text, MusicTable& t, std::string* error) {
    t.pieces.clear();
    std::istringstream in(text);
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.resize(hash);
        std::istringstream ls(line);
        std::string kind, a, b;
        if (!(ls >> kind)) continue;
        auto fail = [&](const char* why) { if (error) *error = "line " + std::to_string(lineNo) + ": " + why; return false; };
        if (kind == "piece") {
            MusicTable::Piece p;
            if (!(ls >> p.name)) return fail("a piece needs a name");
            while (ls >> a) p.places.push_back(a);
            if (p.places.empty()) return fail("a piece needs at least one place (or any)");
            t.pieces.push_back(p);
        } else if (kind == "layer") {
            if (!(ls >> a >> b)) return fail("expected: layer danger|discovery|achievement <sound>");
            if (a == "danger") t.danger = b;
            else if (a == "discovery") t.discovery = b;
            else if (a == "achievement") t.achievement = b;
            else return fail("unknown layer");
        } else return fail("expected piece or layer");
    }
    return true;
}

int Music_Choose(const MusicTable& t, const std::vector<std::string>& tags, int last, float random01) {
    std::vector<int> best;
    int bestScore = 0;
    // A place a piece names (a cave, space, a biome) counts most, the time of day less, "any" only makes it eligible.
    auto scoreOf = [&](const MusicTable::Piece& p) {
        int s = 0;
        for (const std::string& pl : p.places) {
            if (pl == "any") { s += 1; continue; }
            if (std::find(tags.begin(), tags.end(), pl) != tags.end()) s += pl == "day" || pl == "night" ? 2 : 6;
        }
        return s;
    };
    for (int pass = 0; pass < 2 && best.empty(); ++pass) // the last piece again only if nothing else fits
        for (int i = 0; i < (int)t.pieces.size(); ++i) {
            if (pass == 0 && i == last) continue;
            const int s = scoreOf(t.pieces[i]);
            if (s <= 0 || s < bestScore) continue;
            if (s > bestScore) { best.clear(); bestScore = s; }
            best.push_back(i);
        }
    if (best.empty()) return -1;
    return best[std::min(best.size() - 1, (size_t)(std::clamp(random01, 0.0f, 0.9999f) * best.size()))];
}

double Music_FirstAt(double now, float r) { return now + 180.0 + 180.0 * std::clamp(r, 0.0f, 1.0f); }
double Music_NextAfter(double ended, float r) { return ended + 600.0 + 600.0 * std::clamp(r, 0.0f, 1.0f); }

namespace {
    const char* kBiomeTags[] = { "ocean", "jungle", "savanna", "desert", "badlands", "forest", "meadow", "swamp", "taiga", "tundra", "alpine", "polar" };
    constexpr double kProceduralSeconds = 90.0, kFade = 8.0;
    MusicTable table;
    bool loaded = false;
    std::mt19937 rng((unsigned)time(nullptr));
    std::uniform_real_distribution<float> u01(0.0f, 1.0f);
    double nextAt = -1.0;
    int last = -1, current = -1, pieceEmitter = -1, dangerEmitter = -1;
    double startedAt = 0.0, endsAt = 0.0;
    float pieceVol = 0.0f, danger = 0.0f, dangerTarget = 0.0f;
    double lastSting = -1e9;
    std::unordered_set<int> biomesSeen;
    bool wentUnder = false;
    char status[160] = "off";

    void Load() {
        loaded = true;
        std::ifstream f("resources/data/music.txt");
        if (!f) { fprintf(stderr, "music: resources/data/music.txt is missing: no music\n"); return; }
        std::stringstream ss;
        ss << f.rdbuf();
        std::string err;
        if (!Music_Parse(ss.str(), table, &err)) fprintf(stderr, "music.txt %s\n", err.c_str());
    }

    std::vector<std::string> PlaceTags() {
        const AmbienceInputs& in = Soundscape_Inputs();
        std::vector<std::string> t;
        t.push_back(in.sunElevation > 0.0f ? "day" : "night");
        if (in.air < 0.05f) t.push_back("space");
        if (in.heightAboveGround < -8.0f) t.push_back("cave");
        if (in.heightAboveGround > 300.0f) t.push_back("high");
        int dom = -1;
        float w = 0.0f;
        for (int b = 0; b < (int)Biome::Count; ++b) if (in.biomeWeight[b] > w) { w = in.biomeWeight[b]; dom = b; }
        if (dom >= 0) t.push_back(kBiomeTags[dom]);
        if (dom == (int)Biome::Ocean) t.push_back("sea");
        return t;
    }

    void Stop() {
        if (pieceEmitter >= 0) { Audio_EmitterDestroy(pieceEmitter); pieceEmitter = -1; }
        current = -1;
    }

    bool Start(double now) {
        const int i = Music_Choose(table, PlaceTags(), last, u01(rng));
        if (i < 0) return false;
        const int snd = Audio_Sound(table.pieces[i].name.c_str());
        if (snd < 0) return false;
        Stop();
        pieceEmitter = Audio_EmitterCreate(snd, AudioBus::Music, AudioPriority::High, {}, false);
        if (pieceEmitter < 0) return false;
        current = last = i;
        startedAt = now;
        // A real piece plays once through; a procedural pad loops for kProceduralSeconds.
        endsAt = now + (Audio_SoundFromFile(snd) ? std::max(10.0, (double)Audio_SoundSeconds(snd)) : kProceduralSeconds);
        pieceVol = 0.0f;
        Audio_EmitterSetVolume(pieceEmitter, 0.0f);
        Audio_EmitterSetPlaying(pieceEmitter, true);
        return true;
    }

    void Sting(const std::string& name) {
        const double now = GetTime();
        if (now - lastSting < 60.0) return;
        const int snd = Audio_Sound(name.c_str());
        if (snd < 0) return;
        lastSting = now;
        Audio_PlaySound(snd, AudioBus::Music, AudioPriority::Normal, 0.7f, 1.0f);
    }
}

void Music_Reset() {
    Stop();
    if (dangerEmitter >= 0) { Audio_EmitterDestroy(dangerEmitter); dangerEmitter = -1; }
    nextAt = -1.0;
    biomesSeen.clear();
    wentUnder = false;
    danger = dangerTarget = 0.0f;
}

void Music_Event(MusicEvent e) {
    if (!Audio_V2()) return;
    if (!loaded) Load();
    Sting(e == MusicEvent::Discovery ? table.discovery : table.achievement);
}

void Music_SetDanger(float level) { dangerTarget = std::max(dangerTarget, std::clamp(level, 0.0f, 1.0f)); }

void Music_PlayNow() { if (Audio_V2()) { if (!loaded) Load(); Start(GetTime()); } }

const char* Music_Status() { return status; }

void Music_Update(Scene& world, float dt) {
    (void)world;
    if (!Audio_V2()) return;
    if (!loaded) Load();
    PROFILE_SCOPE("audio:music");
    const double now = GetTime();
    if (nextAt < 0.0) {
        static const bool soon = getenv("RAYLIBGAME_MUSIC_SOON") && getenv("RAYLIBGAME_MUSIC_SOON")[0] == '1'; // playtests: a piece 10 s in
        nextAt = soon ? now + 10.0 : Music_FirstAt(now, u01(rng));
    }
    // Discovery: a biome first reached this session (not the one you start in), the first time well underground.
    const AmbienceInputs& in = Soundscape_Inputs();
    int dom = -1;
    float w = 0.4f;
    for (int b = 0; b < (int)Biome::Count; ++b) if (in.biomeWeight[b] > w) { w = in.biomeWeight[b]; dom = b; }
    if (dom >= 0 && !biomesSeen.count(dom)) { if (!biomesSeen.empty()) Music_Event(MusicEvent::Discovery); biomesSeen.insert(dom); }
    if (!wentUnder && in.heightAboveGround < -10.0f) { wentUnder = true; Music_Event(MusicEvent::Discovery); }
    // Danger: what was set this frame, gliding (in over a second, out over four); it pushes a calm piece out.
    const float rate = dangerTarget > danger ? dt / 1.0f : dt / 4.0f;
    danger += std::clamp(dangerTarget - danger, -rate, rate);
    dangerTarget = 0.0f;
    if (dangerEmitter < 0 && danger > 0.01f) {
        const int snd = Audio_Sound(table.danger.c_str());
        if (snd >= 0) dangerEmitter = Audio_EmitterCreate(snd, AudioBus::Music, AudioPriority::High, {}, false);
    }
    if (dangerEmitter >= 0) { Audio_EmitterSetVolume(dangerEmitter, 0.6f * danger); Audio_EmitterSetPlaying(dangerEmitter, danger > 0.01f); }
    // The calm piece: fades in, plays, fades out by its end (or under danger); then silence until the next.
    if (current < 0 && now >= nextAt && danger < 0.1f) {
        if (!Start(now)) nextAt = now + 60.0; // nothing suits this place: look again in a minute
    }
    if (current >= 0) {
        const double fadeIn = std::clamp((now - startedAt) / kFade, 0.0, 1.0), fadeOut = std::clamp((endsAt - now) / kFade, 0.0, 1.0);
        pieceVol = (float)std::min(fadeIn, fadeOut) * (1.0f - danger);
        Audio_EmitterSetVolume(pieceEmitter, 0.7f * pieceVol);
        if (now >= endsAt) { Stop(); nextAt = Music_NextAfter(now, u01(rng)); }
    }
    snprintf(status, sizeof status, "%s; next in %.0f s; danger %.2f", current >= 0 ? table.pieces[current].name.c_str() : "silence",
             std::max(0.0, nextAt - now), danger);
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1245): the table, choosing by place, the schedule ----
#include "unit_helpers.h"

TEST(Music, ShippedTableAndChoosingByPlace) {
    std::ifstream f("resources/data/music.txt");
    ASSERT_TRUE(f.good());
    std::stringstream ss;
    ss << f.rdbuf();
    MusicTable t;
    std::string err;
    ASSERT_TRUE(Music_Parse(ss.str(), t, &err)) << err;
    ASSERT_GE(t.pieces.size(), 4u);
    auto nameOf = [&](int i) { return i >= 0 ? t.pieces[i].name : std::string("-"); };
    EXPECT_EQ(nameOf(Music_Choose(t, { "night", "cave" }, -1, 0.0f)), "music_pad_cave");
    EXPECT_EQ(nameOf(Music_Choose(t, { "night", "forest" }, -1, 0.0f)), "music_pad_night");
    EXPECT_EQ(nameOf(Music_Choose(t, { "day", "space" }, -1, 0.0f)), "music_pad_space");
    EXPECT_EQ(nameOf(Music_Choose(t, { "day", "desert" }, -1, 0.0f)), "music_pad_day"); // (only "any" and "day" fit)
}

TEST(Music, AvoidsTheLastPieceUnlessItIsTheOnlyOne) {
    MusicTable t;
    ASSERT_TRUE(Music_Parse("piece a day\npiece b day\npiece c cave\n", t));
    for (float r : { 0.0f, 0.5f, 0.99f }) EXPECT_EQ(Music_Choose(t, { "day" }, 0, r), 1);
    EXPECT_EQ(Music_Choose(t, { "cave" }, 2, 0.3f), 2);      // the only fit, again
    EXPECT_EQ(Music_Choose(t, { "space" }, -1, 0.3f), -1);   // nothing fits
    EXPECT_FALSE(Music_Parse("piece lonely\n", t));
    EXPECT_FALSE(Music_Parse("layer loud x\n", t));
    EXPECT_FALSE(Music_Parse("song x day\n", t));
}

TEST(Music, RareSchedule) {
    EXPECT_GE(Music_FirstAt(100.0, 0.0f), 280.0);
    EXPECT_LE(Music_FirstAt(100.0, 1.0f), 460.0);
    EXPECT_GE(Music_NextAfter(1000.0, 0.0f), 1600.0); // at least 10 minutes of silence after a piece
    EXPECT_LE(Music_NextAfter(1000.0, 1.0f), 2200.0); // at most 20
}
#endif
