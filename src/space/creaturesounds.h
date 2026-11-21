#ifndef CREATURESOUNDS_H
#define CREATURESOUNDS_H

// Creature sounds (#1246, audio plan #1003 step 11): the creature framework's events (#1223: gCreatureEvent - idle, alert,
// the attack's telegraph at its wind-up, hurt, death) played where the creature is, on the Creatures bus. A row's own
// sound name (creatures.txt sound_<event>=) plays when the asset slot has it; otherwise the procedural voice of its
// faction ("creature_<event>_<faction>": beast growls, lurker rattles, scavenger screeches, ooze gurgles, undead moans,
// construct tones), pitched a little per species. Telegraphs are Critical (never dropped, captioned as important with an
// edge arc, #1244); creatures alerted to or fighting the player near by raise the music's danger layer (#1245).
// Audio v2 only; no change to the creature code (it only calls the hook).

#include "audiomixer.h" // AudioPriority
#include "creatures.h"

#include <string>

// Pure helpers (tested).
// The voice for a faction: one of the six procedural voices ("beast" for factions without their own).
const char* CreatureSound_Voice(const std::string& faction);
// The procedural fallback's name for an event and a faction ("creature_telegraph_lurker").
std::string CreatureSound_Fallback(CreatureEvent e, const std::string& faction);
AudioPriority CreatureSound_Priority(CreatureEvent e);
// A species' pitch (0.88..1.12, from its id): two species of one faction sound apart.
float CreatureSound_Pitch(const std::string& species);
// How much an event near the player (at `dist` m) raises the danger layer (0..1).
float CreatureSound_Danger(CreatureEvent e, float dist);

void CreatureSounds_Init();               // installs the hook (keeps any hook already set and calls it too)
void CreatureSounds_Update(float dt);     // once a frame: the danger layer, decaying

#endif
