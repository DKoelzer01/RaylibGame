#ifndef AUDIO_H
#define AUDIO_H

#include "audiomixer.h"

// Tiny procedural sound-effect set: generated at startup, no audio assets needed yet.
enum class Sfx { Dig, Build, Craft, Jump, Pickup, Error, Fanfare, Count }; // Fanfare: a tier reached (#226)

// The audio engine v2 (#1238, plan #1003): buses, a 32-voice pool with priorities, ducking and the asset slot
// (resources/audio/<name>.wav|ogg|flac|mp3 replaces a procedural sound). RAYLIBGAME_AUDIO_V2=1 turns it on; off (the
// default until step 3, #1020) keeps the old engine: one Sound per effect, the master volume only.
bool Audio_V2();

void Audio_Init();       // never opens a device on a dedicated server (--headless) or the headless platform build
void Audio_Shutdown();
void Audio_Update();     // once a frame: ducking, finished voices, gains (v2; nothing in the old engine)
// Plays on the effect's default bus and priority (Audio_SfxBus / Audio_SfxPriority).
void Audio_Play(Sfx s, float volume = 1.0f, float pitch = 1.0f);
void Audio_PlayOn(Sfx s, AudioBus bus, AudioPriority priority, float volume = 1.0f, float pitch = 1.0f);
AudioBus Audio_SfxBus(Sfx s);
AudioPriority Audio_SfxPriority(Sfx s);
const char* Audio_SfxName(Sfx s);   // "dig", "build", ...: the asset slot's file name

// Master volume 0..1 (applies to every sound). Default 0.35: the procedural effects are loud at full.
void Audio_SetMasterVolume(float v);
float Audio_GetMasterVolume();
// Per-bus volume 0..1 (v2; Master is Audio_SetMasterVolume). Kept in the old engine too, so Settings round-trips.
void Audio_SetBusVolume(AudioBus b, float v);
float Audio_GetBusVolume(AudioBus b);
void Audio_DuckFor(double seconds); // music and ambience -6 dB for a while (combat; alarms duck by themselves)
AudioMixer::Stats Audio_Stats();
const AudioMixer& Audio_Mixer();

#endif
