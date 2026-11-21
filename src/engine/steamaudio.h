#ifndef STEAMAUDIO_H
#define STEAMAUDIO_H

// Steam Audio (#1020, slice #1292): HRTF and air absorption per voice. Valve's C library (Apache-2.0, headers in
// third_party/steamaudio) loaded at run time from phonon.dll / libphonon.so: no link-time dependency, so a build or a
// machine without it, CI and the dedicated server are unaffected, and audio.cpp falls back to its own panning and
// low-pass (audiospatial.h). tools/fetch_steamaudio.sh puts the DLL in deps/steamaudio; CMake copies it next to game.exe.
// RAYLIBGAME_STEAM_AUDIO=0 keeps it unloaded (A/B).

namespace SteamAudio {

constexpr int kBlock = 256;         // frames per effect call (AudioBlockFifo::kBlock)
constexpr int kSampleRate = 44100;

// Loads the library, creates the context and the default HRTF. False (and the fallback) when it isn't there.
bool Load();
void Shutdown();
bool Active();
const char* Status();               // "4.8.1, HRTF" / "not found (phonon.dll)" / "off (RAYLIBGAME_STEAM_AUDIO=0)"

// Per-voice effects. Process runs on the audio thread, one voice's effects never on two threads at once.
struct VoiceFx;
VoiceFx* CreateVoice();
void DestroyVoice(VoiceFx* v);
void ResetVoice(VoiceFx* v);        // a new sound on the voice: no tail from the last one
struct Params {
    float dir[3] = { 0, 0, -1 };    // listener space (x right, y up, -z ahead), unit
    float spatialBlend = 1.0f;      // 0 = dry stereo, 1 = fully binaural
    float air[3] = { 1, 1, 1 };     // air absorption per band
};
// mono (kBlock) -> left, right (kBlock each)
void Process(VoiceFx* v, const float* mono, float* left, float* right, const Params& p);

} // namespace SteamAudio

#endif
