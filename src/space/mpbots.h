#ifndef MPBOTS_H
#define MPBOTS_H

// Multiplayer test harness (#347): N bot clients in this process, joined to the --host server running here over loopback,
// under simulated network conditions (net.h NetConditions: lag, jitter, loss, reorder; seeded, so a run replays).
// `--test-mp-bots N` (with --host <save>): every bot joins; they ping for a while (round trips under the conditions);
// then join/leave churn - half the bots leave and come back, three rounds - each must get its own session back, and the
// server's online count must follow. Prints per-bot join and round-trip timings (a base for #348's numbers).
// Test-only: single player never creates one.

#include <cstdint>

class Scene;

void MpBots_Start(int count, uint16_t port);
// Per frame after Host_Update. Returns true once finished (the result line has printed).
bool MpBots_Update(Scene& world);
bool MpBots_Active();

#endif
