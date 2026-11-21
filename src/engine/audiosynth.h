#ifndef AUDIOSYNTH_H
#define AUDIOSYNTH_H

// Procedural ambience (#1241): seamless loops built from layers (wind, birdsong, crickets, frogs, insects, surf, rain on the
// ground, rain on leaves), and thunder. The fallback for the beds the asset slot does not have yet: CC0 recordings named
// the same in resources/audio replace them (user decision, plan #1003 section 14). Mono float PCM at 44.1 kHz, roughly -1..1.

#include <string>
#include <vector>

struct AudioBedRecipe {
    float wind = 0.0f, birds = 0.0f, crickets = 0.0f, frogs = 0.0f, insects = 0.0f, surf = 0.0f, rain = 0.0f, leaves = 0.0f;
    float rumble = 0.0f; // (#1242) a deep, slowly heaving rumble (lava under a volcano)
};
// (#1242) A water drop falling into a pool from a cave ceiling: a short rising "plink" and a little splash.
std::vector<float> AudioSynth_Drip(unsigned seed);
// The recipe of a procedural bed by name ("bed_forest_day", ...); false: no such bed.
bool AudioSynth_BedRecipe(const std::string& name, AudioBedRecipe& r);
// A seamless loop of `seconds` mixing the recipe's layers (seeded by name: the same bed every start).
std::vector<float> AudioSynth_Bed(const AudioBedRecipe& r, float seconds, unsigned seed);
// A thunderclap: a crack (near strikes) and a long rumble. near: 0 (far, rumble only) .. 1 (overhead).
std::vector<float> AudioSynth_Thunder(float near, unsigned seed);
// (#1243) Machines: hums by class ("hum_fire", "hum_steam", "hum_lv", "hum_mv", "hum_hv", "hum_ev": seamless loops), the
// factory bed ("factory_bed": many machines at once) and the state cues ("machine_stall" a falling tone, "machine_click",
// "machine_done" a soft chime). False: not a machine sound.
bool AudioSynth_Machine(const std::string& name, unsigned seed, std::vector<float>& out);
// (#1245) Generative music, the fallback until real pieces are in resources/audio: calm pads as 20 s seamless loops
// ("music_pad_day", "music_pad_night", "music_pad_cave", "music_pad_space"), the danger layer ("music_danger", a 12 s
// loop) and two stings ("music_discovery", "music_achievement"). False: not one of them.
bool AudioSynth_Music(const std::string& name, unsigned seed, std::vector<float>& out);

#endif
