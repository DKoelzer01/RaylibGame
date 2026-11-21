#ifndef SOUNDSCAPE_H
#define SOUNDSCAPE_H

// Ambience (#1241, audio plan #1003 step 5): beds by biome, crossfaded by the sun's height and the biomes around the
// listener (resources/data/ambience.txt); rain on the ground and on leaves; wind that strengthens with the weather and
// with height; the beds thin out high up, fade underground and are gone in vacuum, where everything is muffled through
// the suit (user decision); underwater muffles too. Thunder follows lightning by the sound's travel time (weather.cpp).
// Audio v2 only (the old engine has no beds). The beds are CC0 recordings when resources/audio has them, else procedural.

#include "voxel.h" // Biome

#include <string>
#include <utility>
#include <vector>

class Scene;

struct AmbienceTable {
    struct Bed { std::string day, night; float volume = 0.0f; };
    Bed biome[(int)Biome::Count];
    std::string rain = "bed_rain", leaves = "bed_rain_leaves", wind = "wind";
};
// Parses ambience.txt; false (and the line in `error`) on a malformed line. Unknown biomes are errors too.
bool Ambience_Parse(const std::string& text, AmbienceTable& t, std::string* error = nullptr);

struct AmbienceInputs {
    float biomeWeight[(int)Biome::Count] = {}; // the biomes around the listener (sums to 1; all 0 off a planet)
    float sunElevation = 0.0f;   // degrees over the horizon
    float heightAboveGround = 0.0f; // m (negative: underground)
    float air = 1.0f;            // 0 vacuum .. 1 sea level
    float rain = 0.0f, snow = 0.0f; // 0..1
    float wind = 0.0f;           // m/s
    bool underwater = false;
};
// The target volume of every sound named by the table (pure).
void Ambience_Targets(const AmbienceTable& t, const AmbienceInputs& in, std::vector<std::pair<std::string, float>>& out);
// The listener's muffle (Hz, 0 = none): through the suit in thin air or vacuum, underwater.
float Ambience_Muffle(const AmbienceInputs& in);

void Soundscape_Update(Scene& world, float dt);   // once a frame (v2)
// As of the last update (other sound systems, #1242): the air (0 vacuum .. 1) and the listener's height over the ground
// (m, negative underground; large off a planet).
float Soundscape_Air();
float Soundscape_HeightAboveGround();
// One bed a call, made ahead (a procedural bed takes a moment to synthesise): call it each loading-screen frame.
void Soundscape_PreloadStep();
// Thunder `delay` seconds after its flash (the strike is delay x 343 m away, toward world.lightningDir).
void Soundscape_Thunder(Scene& world, float delay);

#endif
