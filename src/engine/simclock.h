#ifndef SIMCLOCK_H
#define SIMCLOCK_H

#include <algorithm>

// The simulation clock (#315, multiplayer phase 1): gameplay runs in fixed steps of SimDt() = 1/60 s (kSimHz, the
// server's tick rate; the player's 120 Hz collision substeps run inside each step). Each frame main.cpp hands the frame's
// time (GameFrameTime(): fixed while recording, so captures stay reproducible) to SimClock_BeginFrame, runs the gameplay
// block as many times as it returns (at most kSimMaxSteps: a hitch drops the rest, as the old 0.1 s frame clamp did) and
// calls SimClock_EndStep after each step. SimTime() is the sum of the steps run: gameplay timers (cooldowns, AI and spawn
// intervals, respawn delays) read it, never the wall clock. SimAlpha() is how far the frame is into the next step (0..1):
// the camera interpolates the player between the last two steps by it. Benches, streaming and UI timing keep GetTime().
constexpr int kSimHz = 60;
constexpr int kSimMaxSteps = 6;
inline double gSimAccum = 0.0;   // frame time not yet stepped (s)
inline long long gSimTick = 0;   // steps run since start
inline float SimDt() { return 1.0f / (float)kSimHz; }
inline double SimTime() { return (double)gSimTick / (double)kSimHz; }
inline long long SimTick() { return gSimTick; }
inline float SimAlpha() { return (float)std::clamp(gSimAccum * kSimHz, 0.0, 1.0); }
inline int SimClock_BeginFrame(float frameDt) {
    gSimAccum += frameDt > 0.0f ? frameDt : 0.0f;
    int steps = (int)(gSimAccum * kSimHz + 1e-6);
    if (steps > kSimMaxSteps) { steps = kSimMaxSteps; gSimAccum = (double)steps / kSimHz; }
    return steps;
}
inline void SimClock_EndStep() { gSimAccum = std::max(0.0, gSimAccum - 1.0 / kSimHz); ++gSimTick; }

#endif
