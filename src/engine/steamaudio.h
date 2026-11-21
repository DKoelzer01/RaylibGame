#ifndef STEAMAUDIO_H
#define STEAMAUDIO_H

// Steam Audio (#1020): HRTF and air absorption per voice (#1292), occlusion and transmission against the world (#1293).
// Valve's C library (Apache-2.0, headers in third_party/steamaudio) loaded at run time from phonon.dll / libphonon.so:
// no link-time dependency, so a build or a machine without it, CI and the dedicated server are unaffected, and audio.cpp
// falls back to its own panning and low-pass (audiospatial.h). tools/fetch_steamaudio.sh puts the DLL in
// deps/steamaudio; CMake copies it next to game.exe. RAYLIBGAME_STEAM_AUDIO=0 keeps it unloaded (A/B);
// RAYLIBGAME_AUDIO_OCCLUSION=0 leaves out the occlusion simulation.

#include "dvec3.h"

#include <memory>

struct AudioScene;

namespace SteamAudio {

constexpr int kBlock = 256;         // frames per effect call (AudioBlockFifo::kBlock)
constexpr int kSampleRate = 44100;
constexpr int kMaxVoices = 32;      // AudioMixer::kVoices

// Loads the library, creates the context and the default HRTF. False (and the fallback) when it isn't there.
bool Load();
void Shutdown();
bool Active();
const char* Status();               // "4.8.1, HRTF" / "not found (phonon.dll)" / "off (RAYLIBGAME_STEAM_AUDIO=0)"

// Per-voice effects. Process runs on the audio thread, one voice's effects never on two threads at once.
struct VoiceFx;
VoiceFx* CreateVoice(int index);    // index: the voice (0..kMaxVoices-1), whose occlusion results it applies
void DestroyVoice(VoiceFx* v);
void ResetVoice(VoiceFx* v);        // a new sound on the voice: no tail from the last one, no occlusion until simulated
struct Params {
    float dir[3] = { 0, 0, -1 };    // listener space (x right, y up, -z ahead), unit
    float spatialBlend = 1.0f;      // 0 = dry stereo, 1 = fully binaural
    float air[3] = { 1, 1, 1 };     // air absorption per band
};
// mono (kBlock) -> left, right (kBlock each)
void Process(VoiceFx* v, const float* mono, float* left, float* right, const Params& p);

// ---- Occlusion and transmission (#1293) ----
// A custom Steam Audio scene whose rays are AudioScene's (the density field and part cells: no meshes, no Embree),
// simulated ~10 times a second on its own thread for the voices that play at a place. The results are smoothed per
// block on the audio thread (a ~25 ms glide), so they never click.
bool StartSimulation();             // after Load; false: not available (then nothing is occluded)
void StopSimulation();
bool Simulating();
struct SimVoice { bool active = false; DVec3 pos{}; };  // a positional voice and where it is (global)
// Main thread, once a frame: the world snapshot, the listener (global, its axes) and every voice.
void SetSimInputs(std::shared_ptr<const AudioScene> scene, const DVec3& listener, const float ahead[3], const float up[3],
                  const SimVoice (&voices)[kMaxVoices]);
struct SimStats { long long ticks = 0; float lastMs = 0.0f, maxMs = 0.0f; int sources = 0; float occlusion[kMaxVoices] = {}; };
SimStats GetSimStats();

} // namespace SteamAudio

#endif
