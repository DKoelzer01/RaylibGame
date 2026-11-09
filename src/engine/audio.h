#ifndef AUDIO_H
#define AUDIO_H

// Tiny procedural sound-effect set: generated at startup, no audio assets needed yet.
enum class Sfx { Dig, Build, Craft, Jump, Pickup, Error, Count };

void Audio_Init();
void Audio_Shutdown();
void Audio_Play(Sfx s, float volume = 1.0f, float pitch = 1.0f);

#endif
