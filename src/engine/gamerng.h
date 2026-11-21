#ifndef GAMERNG_H
#define GAMERNG_H

// Gameplay randomness (#316, multiplayer phase 1): one seeded stream per world. Outcomes that matter to the game (gem
// and seed drops, harvest yields, food sickness, where night hunters appear, meteor impact sites and meteorite ore) roll
// here; in multiplayer only the server rolls and clients receive the outcomes (user, 2026-11-15: every gameplay roll is
// the authority's). The stream is seeded from the world seed by generateWorld and saved with the world, so a seed plus
// the same actions gives the same results. Other gameplay choices are already derived from seeded state that the
// authority simulates: herds and their decisions (per-animal LCG seeded from the spawn hash), snare catches (hash of
// plot, part and set time), structure loot (site seed), the weather schedule and meteor-shower nights (world seed hashes).
// Cosmetic randomness (sound pitch, particles, ambience, meteor streaks, lightning) stays on raylib's GetRandomValue and
// may differ per client.

#include <cstdint>
#include <cstdio>

void GameRng_Seed(uint64_t seed);
int GameRng_Int(int min, int max);   // inclusive, like GetRandomValue
float GameRng_Unit();                // [0, 1)
uint64_t GameRng_State();
void GameRng_Save(FILE* f);          // "game_rng N"
bool GameRng_LoadLine(const char* tok, FILE* f);

#endif
