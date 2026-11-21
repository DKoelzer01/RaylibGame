#ifndef MUSIC_H
#define MUSIC_H

// Music (#1245, audio plan #1003 step 10): rare calm pieces chosen by place (resources/data/music.txt; the user's decision:
// rarely, with the Music slider), a danger layer that fades in under a threat, and stings for discovery and achievement.
// No constant soundtrack: silence is part of the design. Real pieces go in resources/audio (streamed); until then the
// procedural pads stand in. Audio v2 only, on the Music bus (ducked under alarms by the mixer).

#include <string>
#include <vector>

class Scene;

struct MusicTable {
    struct Piece { std::string name; std::vector<std::string> places; };
    std::vector<Piece> pieces;
    std::string danger = "music_danger", discovery = "music_discovery", achievement = "music_achievement";
};
bool Music_Parse(const std::string& text, MusicTable& t, std::string* error = nullptr);
// The piece that best suits a place (its tags: "day", "forest", "cave", ...): most tags matched, a piece naming the place
// itself over "any"; ties by `random01`; the last piece played only when nothing else fits. -1 when nothing does.
int Music_Choose(const MusicTable& t, const std::vector<std::string>& placeTags, int last, float random01);
// When the next piece may start (s of play time): the first 3-6 min in, then 10-20 min after the last one ended.
double Music_FirstAt(double now, float random01);
double Music_NextAfter(double ended, float random01);

enum class MusicEvent { Discovery, Achievement };
void Music_Event(MusicEvent e);       // a sting (at most one a minute)
void Music_SetDanger(float level);    // 0..1 now (creatures, combat): the danger layer follows it; a calm piece fades out
void Music_Update(Scene& world, float dt);
void Music_PlayNow();                 // debug: start a piece now
const char* Music_Status();           // debug: what plays, when the next one may
void Music_Reset();                   // a new or unloaded world

#endif
