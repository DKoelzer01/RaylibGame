#ifndef MPBOTS_H
#define MPBOTS_H

// Multiplayer test harness (#347): N bot clients in this process, joined to the --host server running here over loopback,
// under simulated network conditions (net.h NetConditions: lag, jitter, loss, reorder; seeded, so a run replays).
// `--test-mp-bots N` (with --host <save>): every bot joins; they ping for a while (round trips under the conditions);
// then each walks a scripted path for 8 s predicting itself (#327: circles, sprint + jumps, crouch zigzag; the prediction
// at each acknowledged input within 0.5 m, at most 1 correction in 10 states, bot1 sees the others move smoothly), then
// join/leave churn - half the bots leave and come back, three rounds - each must get its own session back, and the
// server's online count must follow. Prints per-bot join and round-trip timings (a base for #348's numbers).
// RAYLIBGAME_MP_CSV=<file>: a row per own state of bot1 (ticks, server and client body-local feet, error, corrected).
// Test-only: single player never creates one.

#include <cstdint>

class Scene;

void MpBots_Start(int count, uint16_t port);
// Per frame after Host_Update. Returns true once finished (the result line has printed).
bool MpBots_Update(Scene& world);
bool MpBots_Active();

#endif
