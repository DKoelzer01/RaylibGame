#include "gamerng.h"

#include <cstring>

namespace {
uint64_t state = 0x9E3779B97F4A7C15ull;

uint64_t Next() { // splitmix64: tiny, fast, good enough for drops
    uint64_t z = (state += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
}
}

void GameRng_Seed(uint64_t seed) { state = seed * 0x2545F4914F6CDD1Dull + 0x6A09E667F3BCC909ull; }

int GameRng_Int(int min, int max) {
    if (max <= min) return min;
    const uint64_t span = (uint64_t)((int64_t)max - (int64_t)min) + 1u;
    return (int)((int64_t)min + (int64_t)(Next() % span));
}

float GameRng_Unit() { return (float)((Next() >> 40) * (1.0 / 16777216.0)); }
uint64_t GameRng_State() { return state; }
void GameRng_Save(FILE* f) { fprintf(f, "game_rng %llu\n", (unsigned long long)state); }

bool GameRng_LoadLine(const char* tok, FILE* f) {
    if (strcmp(tok, "game_rng") != 0) return false;
    unsigned long long s = 0;
    if (fscanf(f, " %llu", &s) == 1) state = s;
    return true;
}
