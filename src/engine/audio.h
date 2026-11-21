#ifndef AUDIO_H
#define AUDIO_H

#include "audiomixer.h"
#include "audiospatial.h"
#include "frames.h"

// Tiny procedural sound-effect set: generated at startup, no audio assets needed yet.
enum class Sfx { Dig, Build, Craft, Jump, Pickup, Error, Fanfare, Count }; // Fanfare: a tier reached (#226)

// The audio engine v2 (#1238, plan #1003): buses, a 32-voice pool with priorities, ducking and the asset slot
// (resources/audio/<name>.wav|ogg|flac|mp3 replaces a procedural sound). RAYLIBGAME_AUDIO_V2=1 turns it on; off (the
// default until step 3, #1020) keeps the old engine: one Sound per effect, the master volume only.
bool Audio_V2();

void Audio_Init();       // never opens a device on a dedicated server (--headless) or the headless platform build
void Audio_Shutdown();
void Audio_Update();     // once a frame: emitters, ducking, finished voices, gains, panning (v2; nothing in the old engine)
// Plays on the effect's default bus and priority (Audio_SfxBus / Audio_SfxPriority).
void Audio_Play(Sfx s, float volume = 1.0f, float pitch = 1.0f);
void Audio_PlayOn(Sfx s, AudioBus bus, AudioPriority priority, float volume = 1.0f, float pitch = 1.0f);
AudioBus Audio_SfxBus(Sfx s);
AudioPriority Audio_SfxPriority(Sfx s);
const char* Audio_SfxName(Sfx s);   // "dig", "build", ...: the asset slot's file name

// ---- Positional sound and emitters (#1239, v2) ----
// The listener: the camera, in global coordinates (Scene::cameraWorld) and its axes. Set once a frame before Audio_Update.
void Audio_SetListener(const DVec3& posGlobal, Vector3 forward, Vector3 up);
const AudioListener& Audio_GetListener();
void Audio_SetMuffle(float cutoffHz);   // > 0: everything low-passed (vacuum through the suit, underwater); 0 = off

// Sounds by name (v2): the Sfx names, the procedural loops ("hum", "wind", "water"), or any file in resources/audio
// (a long .ogg/.mp3/.flac over 1 MB, or one named music_* / bed_*, is streamed: centred, no low-pass, one at a time).
// -1 when there is no such sound or audio is off. Sfx s is sound id (int)s.
int Audio_Sound(const char* name);
// One-shots at a place (the old engine: a 2D sound with the distance's gain).
void Audio_PlayAt(Sfx s, const DVec3& posGlobal, float volume = 1.0f, float pitch = 1.0f, const AudioAttenuation& att = {});
void Audio_PlaySound(int sound, AudioBus bus, AudioPriority priority, float volume = 1.0f, float pitch = 1.0f);
void Audio_PlaySoundAt(int sound, const DVec3& posGlobal, AudioBus bus, AudioPriority priority, float volume = 1.0f,
                       float pitch = 1.0f, const AudioAttenuation& att = {});

// Looping emitters (v2; -1 and no-ops otherwise). Created stopped: set a position and Audio_EmitterSetPlaying(e, true).
// A positional emitter's place is either global (SetPosition, kept as is) or in a body's frame (SetFrame, converted to
// global every frame, so a machine on a spinning planet stays put). positional = false: a centred bed (wind, rain).
int Audio_EmitterCreate(int sound, AudioBus bus, AudioPriority priority, const AudioAttenuation& att = {}, bool positional = true);
void Audio_EmitterDestroy(int e);
void Audio_EmitterSetPosition(int e, const DVec3& posGlobal);
void Audio_EmitterSetFrame(int e, const FramePos& p);
void Audio_EmitterSetPlaying(int e, bool playing);
void Audio_EmitterSetVolume(int e, float volume);
void Audio_EmitterSetPitch(int e, float pitch);
void Audio_EmitterStats(int& count, int& voiced);
// (#1292) How sounds are placed: "steam audio 4.8.1, HRTF", or our own panning and why (not found, off, ...).
const char* Audio_SpatialBackend();

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
