#include "subtitles.h"

#include <raylib.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <unordered_map>

bool SubtitleTable::parse(const std::string& text, std::string* error) {
    rows.clear();
    std::istringstream in(text);
    std::string line;
    int lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t hash = line.find('#');
        if (hash != std::string::npos) line.resize(hash);
        std::istringstream ls(line);
        std::string key, level;
        if (!(ls >> key)) continue;
        std::string rest;
        std::getline(ls >> level, rest);
        const size_t a = rest.find_first_not_of(" \t");
        rest = a == std::string::npos ? "" : rest.substr(a);
        while (!rest.empty() && (rest.back() == ' ' || rest.back() == '\t')) rest.pop_back();
        if ((level != "important" && level != "all") || rest.empty()) {
            if (error) *error = "line " + std::to_string(lineNo) + ": expected: sound important|all text";
            return false;
        }
        Row r;
        r.prefix = key.back() == '*';
        r.key = r.prefix ? key.substr(0, key.size() - 1) : key;
        r.cap = { rest, level == "important" };
        rows.push_back(r);
    }
    return true;
}

const SubtitleCaption* SubtitleTable::find(const std::string& sound) const {
    const SubtitleCaption* best = nullptr;
    size_t bestLen = 0;
    for (const Row& r : rows) {
        if (!r.prefix) { if (r.key == sound) return &r.cap; continue; }
        if (sound.rfind(r.key, 0) == 0 && (!best || r.key.size() > bestLen)) { best = &r.cap; bestLen = r.key.size(); }
    }
    return best;
}

float Subtitle_Angle(const AudioListener& l, const DVec3& pos) {
    const DVec3 d = pos - l.pos;
    const Vector3 rel = { (float)d.x, (float)d.y, (float)d.z };
    auto norm = [](Vector3 v) { const float n = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z); return n > 1e-6f ? Vector3{ v.x / n, v.y / n, v.z / n } : Vector3{ 0, 0, 0 }; };
    const Vector3 f = norm(l.forward);
    const Vector3 r = norm({ f.y * l.up.z - f.z * l.up.y, f.z * l.up.x - f.x * l.up.z, f.x * l.up.y - f.y * l.up.x });
    const float x = rel.x * r.x + rel.y * r.y + rel.z * r.z, z = rel.x * f.x + rel.y * f.y + rel.z * f.z;
    if (fabsf(x) < 1e-6f && fabsf(z) < 1e-6f) return 0.0f;
    return atan2f(x, z) * 57.29578f;
}

void Subtitle_Arrows(float a, std::string& left, std::string& right) {
    left.clear();
    right.clear();
    if (fabsf(a) <= 35.0f) return;                  // ahead: in view
    if (fabsf(a) >= 145.0f) { left = "v "; return; } // behind
    if (a > 0.0f) right = " >"; else left = "< ";
}

void SubtitleLog::add(const SubtitleCaption& c, bool directional, float angleDeg, double now) {
    for (Line& l : lines)
        if (l.text == c.text) { l.last = now; l.directional = directional; l.angle = angleDeg; l.important = l.important || c.important; return; }
    Line l;
    l.text = c.text;
    l.important = c.important;
    l.directional = directional;
    l.angle = angleDeg;
    l.first = l.last = now;
    lines.push_back(l);
}

void SubtitleLog::expire(double now) {
    lines.erase(std::remove_if(lines.begin(), lines.end(), [&](const Line& l) { return now - l.last > kShowSeconds; }), lines.end());
}

void SubtitleLog::visible(SubtitleLevel level, std::vector<const Line*>& out) const {
    out.clear();
    if (level == SubtitleLevel::Off) return;
    for (const Line& l : lines)
        if (level == SubtitleLevel::All || l.important) out.push_back(&l);
    std::sort(out.begin(), out.end(), [](const Line* a, const Line* b) { return a->last > b->last; });
    if ((int)out.size() > kMaxLines) out.resize(kMaxLines);
}

namespace {
    SubtitleTable table;
    SubtitleLog logLines;
    SubtitleLevel level = SubtitleLevel::Important;
    bool loaded = false;
    std::unordered_map<std::string, SubtitleCaption> overrides; // (#1635) Subtitles_SetCaption: exact sound names, before the table
}

void Subtitles_Load() {
    if (loaded) return;
    loaded = true;
    std::ifstream f("resources/data/subtitles.txt");
    if (!f) return;
    std::stringstream ss;
    ss << f.rdbuf();
    std::string err;
    if (!table.parse(ss.str(), &err)) fprintf(stderr, "subtitles.txt %s\n", err.c_str());
}

void Subtitles_SetLevel(SubtitleLevel l) { level = l; }
SubtitleLevel Subtitles_Level() { return level; }

void Subtitles_SetCaption(const std::string& sound, const std::string& text, bool important) {
    if (sound.empty() || text.empty()) return;
    SubtitleCaption c; c.text = text; c.important = important;
    overrides[sound] = c;
}

void Subtitles_Event(const std::string& text, bool important, const DVec3* pos, const AudioListener& listener) {
    if (level == SubtitleLevel::Off || text.empty() || (level == SubtitleLevel::Important && !important)) return;
    SubtitleCaption c; c.text = text; c.important = important;
    logLines.add(c, pos != nullptr, pos ? Subtitle_Angle(listener, *pos) : 0.0f, GetTime());
}

void Subtitles_OnSound(const std::string& sound, const DVec3* pos, const AudioListener& listener) {
    if (level == SubtitleLevel::Off) return;
    const auto ov = overrides.find(sound); // (#1635) a creature's own words first
    const SubtitleCaption* c = ov != overrides.end() ? &ov->second : table.find(sound);
    if (!c || (level == SubtitleLevel::Important && !c->important)) return;
    logLines.add(*c, pos != nullptr, pos ? Subtitle_Angle(listener, *pos) : 0.0f, GetTime());
}

void Subtitles_Draw(float scale) {
    const double now = GetTime();
    logLines.expire(now);
    std::vector<const SubtitleLog::Line*> shown;
    logLines.visible(level, shown);
    if (shown.empty()) return;
    const int fs = (int)(18 * scale), pad = (int)(6 * scale), lineH = fs + 2 * pad;
    const int sw = GetScreenWidth(), sh = GetScreenHeight();
    int y = sh - (int)(150 * scale);
    std::string l, r;
    for (const SubtitleLog::Line* line : shown) { // newest at the bottom
        Subtitle_Arrows(line->directional ? line->angle : 0.0f, l, r);
        const std::string s = l + line->text + r;
        const int w = MeasureText(s.c_str(), fs);
        const float fade = (float)std::clamp((SubtitleLog::kShowSeconds - (now - line->last)) / 0.5, 0.0, 1.0);
        const int x = sw - w - 2 * pad - (int)(20 * scale);
        DrawRectangle(x - pad, y - lineH + pad, w + 2 * pad, lineH, Fade(BLACK, 0.55f * fade));
        DrawText(s.c_str(), x, y - lineH + 2 * pad, fs, Fade(line->important ? Color{ 255, 220, 120, 255 } : RAYWHITE, fade));
        y -= lineH + (int)(2 * scale);
        // Important sounds out of view: an arc at the screen's edge toward them (0 = ahead = the top of the screen).
        if (line->important && line->directional && fabsf(line->angle) > 35.0f && now - line->last < 1.5) {
            const float a = line->angle - 90.0f;
            const float rad = 0.42f * (float)std::min(sw, sh);
            DrawRing({ sw * 0.5f, sh * 0.5f }, rad - 6.0f * scale, rad, a - 18.0f, a + 18.0f, 24, Fade(Color{ 255, 90, 60, 255 }, 0.7f * fade));
        }
    }
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#1244): the table, arrows, the log ----
#include "unit_helpers.h"

TEST(Subtitles, TableExactThenLongestPrefix) {
    SubtitleTable t;
    std::string err;
    ASSERT_TRUE(t.parse("# c\nmachine_stall important Machine stalled\nhum_* all Machine humming\nhum_ev* all Fusion humming\n", &err)) << err;
    ASSERT_NE(t.find("machine_stall"), nullptr);
    EXPECT_TRUE(t.find("machine_stall")->important);
    EXPECT_EQ(t.find("hum_lv")->text, "Machine humming");
    EXPECT_EQ(t.find("hum_ev")->text, "Fusion humming");
    EXPECT_EQ(t.find("dig"), nullptr);
    EXPECT_FALSE(t.parse("x sometimes text\n", &err));
    EXPECT_FALSE(t.parse("x all\n", &err));
}

TEST(Subtitles, ShippedTableParses) {
    std::ifstream f("resources/data/subtitles.txt");
    ASSERT_TRUE(f.good());
    std::stringstream ss;
    ss << f.rdbuf();
    SubtitleTable t;
    std::string err;
    ASSERT_TRUE(t.parse(ss.str(), &err)) << err;
    ASSERT_NE(t.find("machine_stall"), nullptr);
    EXPECT_TRUE(t.find("machine_stall")->important);
    EXPECT_NE(t.find("step_grass"), nullptr);
}

TEST(Subtitles, DirectionArrows) {
    AudioListener l;
    l.pos = { 0, 0, 0 };
    l.forward = { 0, 0, -1 };
    l.up = { 0, 1, 0 };
    EXPECT_NEAR(Subtitle_Angle(l, { 0, 0, -10 }), 0.0f, 1e-3f);
    EXPECT_NEAR(Subtitle_Angle(l, { 10, 0, 0 }), 90.0f, 1e-3f);
    EXPECT_NEAR(Subtitle_Angle(l, { -10, 0, 0 }), -90.0f, 1e-3f);
    EXPECT_NEAR(fabsf(Subtitle_Angle(l, { 0, 0, 10 })), 180.0f, 1e-3f);
    std::string a, b;
    Subtitle_Arrows(0.0f, a, b);   EXPECT_EQ(a + b, "");
    Subtitle_Arrows(90.0f, a, b);  EXPECT_EQ(b, " >");
    Subtitle_Arrows(-90.0f, a, b); EXPECT_EQ(a, "< ");
    Subtitle_Arrows(170.0f, a, b); EXPECT_EQ(a, "v ");
}

TEST(Subtitles, LogMergesFiltersAndExpires) {
    SubtitleLog log;
    log.add({ "Footsteps", false }, true, 90.0f, 1.0);
    log.add({ "Machine stalled", true }, true, -90.0f, 1.2);
    log.add({ "Footsteps", false }, true, 0.0f, 2.0);  // the same caption again: one line, refreshed
    std::vector<const SubtitleLog::Line*> v;
    log.visible(SubtitleLevel::All, v);
    ASSERT_EQ(v.size(), 2u);
    EXPECT_EQ(v[0]->text, "Footsteps");                 // newest first
    EXPECT_FLOAT_EQ(v[0]->angle, 0.0f);
    log.visible(SubtitleLevel::Important, v);
    ASSERT_EQ(v.size(), 1u);
    EXPECT_EQ(v[0]->text, "Machine stalled");
    log.visible(SubtitleLevel::Off, v);
    EXPECT_TRUE(v.empty());
    log.expire(1.2 + SubtitleLog::kShowSeconds + 0.1);  // the stall is gone, the steps (2.0) stay
    log.visible(SubtitleLevel::All, v);
    ASSERT_EQ(v.size(), 1u);
    EXPECT_EQ(v[0]->text, "Footsteps");
    for (int i = 0; i < 10; ++i) log.add({ "S" + std::to_string(i), false }, false, 0.0f, 3.0 + i * 0.01);
    log.visible(SubtitleLevel::All, v);
    EXPECT_EQ((int)v.size(), SubtitleLog::kMaxLines);
}
#endif
