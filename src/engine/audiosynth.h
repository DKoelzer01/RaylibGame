#ifndef AUDIOSYNTH_H
#define AUDIOSYNTH_H

// Procedural ambience (#1241): seamless loops built from layers (wind, birdsong, crickets, frogs, insects, surf, rain on the
// ground, rain on leaves), and thunder. The fallback for the beds the asset slot does not have yet: CC0 recordings named
// the same in resources/audio replace them (user decision, plan #1003 section 14). Mono float PCM at 44.1 kHz, roughly -1..1.

#include <string>
#include <vector>

struct AudioBedRecipe {
    float wind = 0.0f, birds = 0.0f, crickets = 0.0f, frogs = 0.0f, insects = 0.0f, surf = 0.0f, rain = 0.0f, leaves = 0.0f;
};
// The recipe of a procedural bed by name ("bed_forest_day", ...); false: no such bed.
bool AudioSynth_BedRecipe(const std::string& name, AudioBedRecipe& r);
// A seamless loop of `seconds` mixing the recipe's layers (seeded by name: the same bed every start).
std::vector<float> AudioSynth_Bed(const AudioBedRecipe& r, float seconds, unsigned seed);
// A thunderclap: a crack (near strikes) and a long rumble. near: 0 (far, rumble only) .. 1 (overhead).
std::vector<float> AudioSynth_Thunder(float near, unsigned seed);

#endif
