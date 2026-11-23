#ifndef SUBTITLES_H
#define SUBTITLES_H

// Sound subtitles and direction cues (#1244, audio plan #1003 step 8). Sounds are captioned by name from
// resources/data/subtitles.txt (audio.cpp reports every one-shot it plays and every emitter while it is heard). A list at
// the bottom right shows the captions of the last few seconds with an arrow toward each source ("< Machine stalled");
// important sounds out of view also get an arc at the screen's edge pointing at them. Settings: Off, Important (the
// default, user decision), All. Audio v2 only (the old engine's sounds have no names).

#include "audiospatial.h" // AudioListener, DVec3

#include <string>
#include <vector>

enum class SubtitleLevel : uint8_t { Off = 0, Important = 1, All = 2 };

struct SubtitleCaption { std::string text; bool important = false; };

// The caption table (pure).
class SubtitleTable {
public:
    // "name level text..." lines ('#' comments); a name ending in '*' is a prefix. False on a malformed line.
    bool parse(const std::string& text, std::string* error = nullptr);
    // The caption of a sound name (exact names first, then the longest matching prefix), or nullptr.
    const SubtitleCaption* find(const std::string& sound) const;

private:
    struct Row { std::string key; bool prefix = false; SubtitleCaption cap; };
    std::vector<Row> rows;
};

// The direction of a source in the listener's view: degrees, 0 ahead, +90 right, -90 left, +-180 behind.
float Subtitle_Angle(const AudioListener& l, const DVec3& pos);
// The arrows around a caption for that direction: {left, right} ("<", "") / ("", ">") / ("^ ", "") ahead / ("v ", "") behind.
void Subtitle_Arrows(float angleDeg, std::string& left, std::string& right);

// The lines on screen (pure): one per caption text, refreshed while its sound repeats or plays on.
class SubtitleLog {
public:
    static constexpr double kShowSeconds = 2.5;
    static constexpr int kMaxLines = 6;
    struct Line { std::string text; bool important = false, directional = false; float angle = 0.0f; double last = 0.0, first = 0.0; };
    void add(const SubtitleCaption& c, bool directional, float angleDeg, double now);
    void expire(double now);
    // The lines a level shows, newest first, at most kMaxLines.
    void visible(SubtitleLevel level, std::vector<const Line*>& out) const;

private:
    std::vector<Line> lines;
};

void Subtitles_Load();                          // the table (once; audio.cpp calls it on start)
void Subtitles_SetLevel(SubtitleLevel level);
SubtitleLevel Subtitles_Level();
// A sound heard (audio.cpp): its name, where it is (nullptr: everywhere, no direction), the listener.
// (#1635) A caption for one sound name, ahead of subtitles.txt (a creature's telegraph sound: its own wind-up words).
void Subtitles_SetCaption(const std::string& sound, const std::string& text, bool important);
// (#1635) A caption with no sound behind it (any audio engine: a creature's wind-up or hit). Same text as a sound's caption: one line.
void Subtitles_Event(const std::string& text, bool important, const DVec3* pos, const AudioListener& listener);
void Subtitles_OnSound(const std::string& sound, const DVec3* pos, const AudioListener& listener);
// The HUD: the list and the edge arcs. `scale`: the interface size.
void Subtitles_Draw(float scale);

#endif
