#ifndef TRAPPING_H
#define TRAPPING_H

// Trapping (#176): snares set on the ground (hold one, RMB) catch small game now and then: about every 20 minutes of
// world time a snare gets a chance, better where small game lives (grass and woods) than in desert, ice or swamp.
// A snare holding a catch shows it; walking up to it (within 2.5 m) takes it: raw meat, and often a hide. Snares are
// checked near the player only, but on world time, so ones left alone have caught something on your return.

#include <cstdio>
#include <cstdint>

class Scene;
class Planetoid;

void Trapping_Update(Scene& world);   // once a second inside
void Trapping_Reset();
void Trapping_Save(FILE* f);
bool Trapping_LoadLine(const char* tok, FILE* f);
bool Trapping_HasCatch(uint64_t bodyId, uint32_t plotId, int part); // the model shows the catch
// Tests (--test-trapping): every snare near the player set back `seconds` of world time (a catch is due) and forced
// to catch; returns how many snares were found.
int Trapping_TestAge(Scene& world, double seconds);

#endif
