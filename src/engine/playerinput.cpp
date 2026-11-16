#include "playerinput.h"

namespace {
// The bindings: a key, or a mouse button (kMouse + button).
constexpr int kMouse = 1000;
constexpr int kBinding[(int)Act::Count] = {
    KEY_W, KEY_S, KEY_A, KEY_D, KEY_LEFT_SHIFT, KEY_SPACE, KEY_LEFT_CONTROL,
    kMouse + MOUSE_BUTTON_LEFT, kMouse + MOUSE_BUTTON_RIGHT,
    KEY_F, KEY_B, KEY_H, KEY_L,
    KEY_E, KEY_I, KEY_O, KEY_R, KEY_V, KEY_G,
    KEY_X, KEY_M, KEY_P, KEY_PERIOD, KEY_COMMA,
    KEY_K, KEY_T, KEY_U, KEY_Z, KEY_C,
    KEY_LEFT_ALT, KEY_F5, KEY_F9,
};
constexpr int kSlotKeys[10] = { KEY_ONE, KEY_TWO, KEY_THREE, KEY_FOUR, KEY_FIVE, KEY_SIX, KEY_SEVEN, KEY_EIGHT, KEY_NINE, KEY_ZERO };
static_assert((int)Act::Count <= 64, "Act bits fit a uint64_t");

PlayerInput frameIn, stepIn;

bool Down(int b) { return b >= kMouse ? IsMouseButtonDown(b - kMouse) : IsKeyDown(b); }
bool Pressed(int b) { return b >= kMouse ? IsMouseButtonPressed(b - kMouse) : IsKeyPressed(b); }
bool Released(int b) { return b >= kMouse ? IsMouseButtonReleased(b - kMouse) : IsKeyReleased(b); }
}

const PlayerInput& PlayerInput::None() { static const PlayerInput none; return none; }

void PlayerInput_BeginFrame() {
    PlayerInput f;
    f.enabled = true;
    for (int a = 0; a < (int)Act::Count; ++a) {
        const uint64_t bit = 1ull << a;
        if (Down(kBinding[a]) || (gFakeHeld & bit)) f.held |= bit;
        if (Pressed(kBinding[a])) f.pressed |= bit;
        if (Released(kBinding[a])) f.released |= bit;
    }
    f.look = GetMouseDelta();
    f.look.x += gFakeLook.x; f.look.y += gFakeLook.y; // --play-load turning (#509)
    f.wheel = GetMouseWheelMove();
    for (int i = 0; i < 10; ++i) if (IsKeyPressed(kSlotKeys[i])) { f.slot = i; break; }
    frameIn = f;
    // The step's input: held state is current; edges, motion and wheel add up until a step takes them.
    stepIn.enabled = true;
    stepIn.held = f.held;
    stepIn.pressed |= f.pressed;
    stepIn.released |= f.released;
    stepIn.look.x += f.look.x; stepIn.look.y += f.look.y;
    stepIn.wheel += f.wheel;
    if (f.slot >= 0) stepIn.slot = f.slot;
}

void PlayerInput_EndStep() {
    stepIn.pressed = stepIn.released = 0;
    stepIn.look = { 0, 0 };
    stepIn.wheel = 0.0f;
    stepIn.slot = -1;
}

const PlayerInput& PlayerInput_Frame() { return frameIn; }
const PlayerInput& PlayerInput_Step() { return stepIn; }
