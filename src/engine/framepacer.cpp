// (#1583) Frame pacing: see framepacer.h. No raylib header here (windows.h clashes with it); raylib's SetTargetFPS is
// declared by hand.
#include "framepacer.h"

#include <chrono>
#include <cstdlib>
#include <thread>

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#ifndef CREATE_WAITABLE_TIMER_HIGH_RESOLUTION
#define CREATE_WAITABLE_TIMER_HIGH_RESOLUTION 0x00000002
#endif
#endif

extern "C" void SetTargetFPS(int fps); // raylib

namespace {
using Clock = std::chrono::steady_clock;
int gTarget = 0;
Clock::time_point gNext{};
bool gHaveNext = false;

void SleepUntil(Clock::time_point deadline) {
    // Sleep to ~0.5 ms before the deadline, then spin (a sleep can wake late; the spin is exact).
    const auto spinFrom = deadline - std::chrono::microseconds(500);
#if defined(_WIN32)
    static HANDLE timer = CreateWaitableTimerExW(nullptr, nullptr, CREATE_WAITABLE_TIMER_HIGH_RESOLUTION, TIMER_ALL_ACCESS);
    if (!timer) timer = CreateWaitableTimerExW(nullptr, nullptr, 0, TIMER_ALL_ACCESS); // older Windows: the plain timer
    const auto now = Clock::now();
    if (timer && spinFrom > now) {
        LARGE_INTEGER due;
        due.QuadPart = -(LONGLONG)(std::chrono::duration_cast<std::chrono::nanoseconds>(spinFrom - now).count() / 100); // relative, 100 ns
        if (due.QuadPart < 0 && SetWaitableTimer(timer, &due, 0, nullptr, nullptr, FALSE)) WaitForSingleObject(timer, INFINITE);
    }
#else
    if (spinFrom > Clock::now()) std::this_thread::sleep_until(spinFrom);
#endif
    while (Clock::now() < deadline) std::this_thread::yield();
}
} // namespace

void FramePacer_SetTarget(int fps) {
    gHaveNext = false;
    // RAYLIBGAME_RAYLIB_PACING=1: raylib's own cap inside EndDrawing, as before (an A/B, and the way back)
    static const bool raylibPacing = [] { const char* v = std::getenv("RAYLIBGAME_RAYLIB_PACING"); return v && *v && *v != '0'; }();
    if (raylibPacing) { gTarget = 0; SetTargetFPS(fps > 0 ? fps : 0); return; }
    gTarget = fps > 0 ? fps : 0;
    SetTargetFPS(0); // the game paces; raylib doesn't wait
}

int FramePacer_Target() { return gTarget; }

void FramePacer_Wait() {
    if (gTarget <= 0) { gHaveNext = false; return; }
    const auto period = std::chrono::nanoseconds(1000000000LL / gTarget);
    const auto now = Clock::now();
    // A frame that ran over (or the first): start a new rhythm from now, no catching up with short frames.
    if (!gHaveNext || now > gNext + period) { gNext = now + period; gHaveNext = true; return; }
    if (now < gNext) SleepUntil(gNext);
    gNext += period;
}
