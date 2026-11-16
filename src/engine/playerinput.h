#ifndef PLAYERINPUT_H
#define PLAYERINPUT_H

#include <raylib.h>
#include <cstdint>

// One player's input (#312, multiplayer phase 1): the gameplay systems (player, terrain tool, build tool, machines,
// hotbar, survival, animals, ship, rocket) read only this struct, never the keyboard or mouse, so a server can feed
// them a remote player's input from the network. The local client fills it from raylib (PlayerInput_BeginFrame).
// UI-only input (menus, the inventory screen, the map, the telescope zoom, the debug UI, the free-fly camera) stays
// client-side and keeps reading raylib.
enum class Act : int {
    Forward, Back, Left, Right, Sprint, Jump, Down,  // W S A D, Shift, Space (jump, jetpack up, swim up), Ctrl
    Primary, Secondary,                              // LMB (dig, break, hit), RMB (place, use, eat)
    Fly, Dampeners, Hover, Headlamp,                 // F, B, H, L
    Use, Insert, Extract, Rotate, Flip, BuildMode,   // E (board, leave), I, O, R, V, G
    Brake, Map, Autopilot, WarpUp, WarpDown,         // X, M, P, '.', ','
    Dock, Transfer, Undock, RollLeft, RollRight,     // K, T, U, Z, C
    FreeLook, Save, Load,                            // Alt (rocket camera), F5, F9
    Count
};

struct PlayerInput {
    bool enabled = false;          // gameplay input is live (false: a menu is open, dead, loading, a bench)
    uint64_t held = 0;             // Act bits held now
    uint64_t pressed = 0;          // went down since the last step (or this frame, for PlayerInput_Frame)
    uint64_t released = 0;         // went up since then
    Vector2 look = { 0, 0 };       // mouse motion (pixels) over the same span
    float wheel = 0.0f;            // wheel clicks over the same span
    int slot = -1;                 // hotbar key (0-9) pressed over the same span, -1 none

    bool down(Act a) const { return enabled && (held >> (int)a & 1u); }
    bool hit(Act a) const { return enabled && (pressed >> (int)a & 1u); }
    bool up(Act a) const { return enabled && (released >> (int)a & 1u); }
    static const PlayerInput& None(); // nothing held or pressed, not enabled
    // Scripted input (benches, later a remote player): hold or press an action.
    void hold(Act a) { enabled = true; held |= 1ull << (int)a; }
    void press(Act a) { enabled = true; pressed |= 1ull << (int)a; }
};

// The local player's devices. Once a frame, before the simulation steps: PlayerInput_BeginFrame reads raylib.
// PlayerInput_Frame(): this frame's input (per-frame code: the tools aim and act with this frame's view).
// PlayerInput_Step(): what the next simulation step sees: held state now, edges, motion and wheel since the last step
// (a frame may run no step, or two: world.h SimClock). PlayerInput_EndStep after each step clears the edges.
void PlayerInput_BeginFrame();
void PlayerInput_EndStep();
const PlayerInput& PlayerInput_Frame();
inline Vector2 gFakeLook = { 0.0f, 0.0f }; // added to the mouse motion every frame (--play-new turning, #387)
inline uint64_t gFakeHeld = 0; // actions held as if their keys were down (bit = 1 << Act; --play-load walking, #509)
const PlayerInput& PlayerInput_Step();
// `in` when enabled, else PlayerInput::None() (the old inputEnabled flags).
inline const PlayerInput& PlayerInput_If(bool enabled, const PlayerInput& in) { return enabled ? in : PlayerInput::None(); }

#endif
