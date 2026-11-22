#include "playerinput.h"

#include <cstdio>
#include <cstring>

namespace {
// The default bindings: a key, or a mouse button (kMouse + button). Settings can change each (#821).
constexpr int kMouse = kMouseBinding;
constexpr int kDefault[(int)Act::Count] = {
    KEY_W, KEY_S, KEY_A, KEY_D, KEY_LEFT_SHIFT, KEY_SPACE, KEY_LEFT_CONTROL,
    kMouse + MOUSE_BUTTON_LEFT, kMouse + MOUSE_BUTTON_RIGHT,
    KEY_F, KEY_B, KEY_H, KEY_L,
    KEY_E, KEY_I, KEY_O, KEY_R, KEY_V, KEY_G,
    KEY_X, KEY_M, KEY_P, KEY_PERIOD, KEY_COMMA,
    KEY_K, KEY_T, KEY_U, KEY_Z, KEY_C,
    KEY_LEFT_ALT, KEY_F5, KEY_F9,
    KEY_Z,
    KEY_TAB, KEY_C, KEY_K,
    KEY_T, KEY_J, KEY_Q, KEY_N, KEY_B, KEY_MINUS, KEY_EQUAL,
    KEY_Q,
    KEY_ENTER,
    KEY_Q, KEY_F,
    kMouse + MOUSE_BUTTON_MIDDLE,
    0, // (#1218) the dodge key: unbound
    KEY_Y,
};
int binding[(int)Act::Count] = {}; // the live bindings (filled from kDefault on first use)
bool bindingsReady = false;
void EnsureBindings() { if (!bindingsReady) { memcpy(binding, kDefault, sizeof(binding)); bindingsReady = true; } }
struct ActInfo { const char* id; const char* label; };
constexpr ActInfo kActs[(int)Act::Count] = {
    { "forward", "Move forward" }, { "back", "Move back" }, { "left", "Move left" }, { "right", "Move right" },
    { "sprint", "Sprint" }, { "jump", "Jump / jetpack up / swim up" }, { "down", "Crouch / down" },
    { "primary", "Dig, break, hit" }, { "secondary", "Place, use, eat" },
    { "fly", "Fly (where allowed)" }, { "dampeners", "Rocket dampeners" }, { "hover", "Hover" }, { "headlamp", "Headlamp" },
    { "use", "Use: board or leave a rocket; on foot, the inventory" }, { "insert", "Insert into a machine" }, { "extract", "Take from a machine" },
    { "rotate", "Rotate a part" }, { "flip", "Flip a part" }, { "build", "Build mode" },
    { "brake", "Rocket brake" }, { "map", "Map" }, { "autopilot", "Rocket autopilot" }, { "warp_up", "Time warp faster" }, { "warp_down", "Time warp slower" },
    { "dock", "Rocket: dock" }, { "transfer", "Rocket: transfer" }, { "undock", "Rocket: undock" }, { "roll_left", "Rocket: roll left" }, { "roll_right", "Rocket: roll right" },
    { "free_look", "Rocket: free look" }, { "save", "Quick save" }, { "load", "Quick load" },
    { "prone", "Lie down / get up" },
    { "inventory", "Open the inventory" }, { "crafting", "Open crafting" }, { "skills", "Open skills" },
    { "telescope", "Telescope" }, { "player_list", "Show the players (multiplayer)" }, { "deposit", "Deposit into nearby chests" },
    { "blueprint_cycle", "Next blueprint" }, { "blueprint_save", "Save a blueprint" }, { "zoom_out", "Minimap: zoom out" }, { "zoom_in", "Minimap: zoom in" },
    { "offhand_use", "Use the off hand (food, oxygen)" },
    { "chat", "Command line / chat" },
    { "drop_item", "Inventory: drop from the slot (Ctrl: all)" }, { "swap_offhand", "Inventory: swap the slot with the off hand" },
    { "lock_on", "Lock on to a creature (toggle)" },
    { "dodge", "Dodge (also: double-tap a direction)" },
    { "magic_sight", "Magic sight: tap to show or hide, hold to study" },
};
constexpr int kSlotKeys[10] = { KEY_ONE, KEY_TWO, KEY_THREE, KEY_FOUR, KEY_FIVE, KEY_SIX, KEY_SEVEN, KEY_EIGHT, KEY_NINE, KEY_ZERO };
static_assert((int)Act::Count <= 64, "Act bits fit a uint64_t");

PlayerInput frameIn, stepIn;
bool textFocus = false; // (#1081) a text field has the keyboard

bool Down(int b) { return b <= 0 ? false : b >= kMouse ? IsMouseButtonDown(b - kMouse) : !PlayerInput_KeyBlocked(b) && IsKeyDown(b); }
bool Pressed(int b) { return b <= 0 ? false : b >= kMouse ? IsMouseButtonPressed(b - kMouse) : !PlayerInput_KeyBlocked(b) && IsKeyPressed(b); }
bool Released(int b) { return b <= 0 ? false : b >= kMouse ? IsMouseButtonReleased(b - kMouse) : IsKeyReleased(b); }
}

const PlayerInput& PlayerInput::None() { static const PlayerInput none; return none; }

void PlayerInput_SetTextFocus(bool on) { textFocus = on; }
bool PlayerInput_TextFocus() { return textFocus; }
bool PlayerInput_KeyBlocked(int code) { return textFocus && code > 0 && code < kMouse; }

void PlayerInput_BeginFrame() {
    EnsureBindings();
    PlayerInput f;
    f.enabled = true;
    static uint64_t fakeBefore = 0; // (#1216) scripted presses: a fake hold's first frame is a press, its end a release
    for (int a = 0; a < (int)Act::Count; ++a) {
        const uint64_t bit = 1ull << a;
        if (Down(binding[a]) || (gFakeHeld & bit)) f.held |= bit;
        if (Pressed(binding[a])) f.pressed |= bit;
        if (Released(binding[a])) f.released |= bit;
        if ((gFakeHeld & bit) && !(fakeBefore & bit)) f.pressed |= bit;
        if (!(gFakeHeld & bit) && (fakeBefore & bit)) f.released |= bit;
    }
    fakeBefore = gFakeHeld;
    f.look = GetMouseDelta();
    f.look.x += gFakeLook.x; f.look.y += gFakeLook.y; // --play-load turning (#509)
    f.wheel = GetMouseWheelMove();
    for (int i = 0; i < 10; ++i) if (!PlayerInput_KeyBlocked(kSlotKeys[i]) && IsKeyPressed(kSlotKeys[i])) { f.slot = i; break; }
    if (gFakeSlot >= 0 && gFakeSlot < 10) { f.slot = gFakeSlot; gFakeSlot = -1; } // (#1219) a scripted hotbar key
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

// ---- Key bindings (#821) ----
int PlayerInput_Binding(Act a) { EnsureBindings(); return (int)a >= 0 && a < Act::Count ? binding[(int)a] : -1; }
int PlayerInput_DefaultBinding(Act a) { return (int)a >= 0 && a < Act::Count ? kDefault[(int)a] : -1; }
void PlayerInput_SetBinding(Act a, int code) { EnsureBindings(); if ((int)a >= 0 && a < Act::Count && code > 0) binding[(int)a] = code; }
void PlayerInput_ResetBindings() { memcpy(binding, kDefault, sizeof(binding)); bindingsReady = true; }
const char* PlayerInput_ActId(Act a) { return (int)a >= 0 && a < Act::Count ? kActs[(int)a].id : ""; }
const char* PlayerInput_ActLabel(Act a) { return (int)a >= 0 && a < Act::Count ? kActs[(int)a].label : ""; }
bool PlayerInput_ActFromId(const char* id, Act& out) {
    for (int a = 0; a < (int)Act::Count; ++a) if (strcmp(kActs[a].id, id) == 0) { out = (Act)a; return true; }
    return false;
}
bool PlayerInput_BindingPressed(Act a) { const int b = PlayerInput_Binding(a); return b > 0 && Pressed(b); }

const char* PlayerInput_KeyName(int code) {
    static char buf[24];
    if (code <= 0) return "Unbound"; // (#1218) the dodge key starts unbound
    if (code >= kMouse) {
        const int m = code - kMouse;
        return m == MOUSE_BUTTON_LEFT ? "Left mouse" : m == MOUSE_BUTTON_RIGHT ? "Right mouse" : m == MOUSE_BUTTON_MIDDLE ? "Middle mouse" : (snprintf(buf, sizeof(buf), "Mouse %d", m + 1), buf);
    }
    if ((code >= KEY_A && code <= KEY_Z) || (code >= KEY_ZERO && code <= KEY_NINE)) { buf[0] = (char)code; buf[1] = 0; return buf; }
    if (code >= KEY_F1 && code <= KEY_F12) { snprintf(buf, sizeof(buf), "F%d", code - KEY_F1 + 1); return buf; }
    if (code >= KEY_KP_0 && code <= KEY_KP_9) { snprintf(buf, sizeof(buf), "Keypad %d", code - KEY_KP_0); return buf; }
    switch (code) {
        case KEY_SPACE: return "Space"; case KEY_TAB: return "Tab"; case KEY_ENTER: return "Enter"; case KEY_BACKSPACE: return "Backspace";
        case KEY_LEFT_SHIFT: return "Left Shift"; case KEY_RIGHT_SHIFT: return "Right Shift";
        case KEY_LEFT_CONTROL: return "Left Ctrl"; case KEY_RIGHT_CONTROL: return "Right Ctrl";
        case KEY_LEFT_ALT: return "Left Alt"; case KEY_RIGHT_ALT: return "Right Alt";
        case KEY_PERIOD: return "."; case KEY_COMMA: return ","; case KEY_MINUS: return "-"; case KEY_EQUAL: return "=";
        case KEY_SLASH: return "/"; case KEY_SEMICOLON: return ";"; case KEY_APOSTROPHE: return "'"; case KEY_GRAVE: return "`";
        case KEY_LEFT_BRACKET: return "["; case KEY_RIGHT_BRACKET: return "]"; case KEY_BACKSLASH: return "\\";
        case KEY_UP: return "Up"; case KEY_DOWN: return "Down"; case KEY_LEFT: return "Left"; case KEY_RIGHT: return "Right";
        case KEY_INSERT: return "Insert"; case KEY_DELETE: return "Delete"; case KEY_HOME: return "Home"; case KEY_END: return "End";
        case KEY_PAGE_UP: return "Page Up"; case KEY_PAGE_DOWN: return "Page Down"; case KEY_CAPS_LOCK: return "Caps Lock";
        default: snprintf(buf, sizeof(buf), "Key %d", code); return buf;
    }
}

int PlayerInput_PressedCode() {
    // Esc cancels a rebind (never bound); the first other key or mouse button down this frame.
    for (int k = KEY_SPACE; k <= KEY_KP_EQUAL; ++k) if (k != KEY_ESCAPE && IsKeyPressed(k)) return k;
    for (int m = MOUSE_BUTTON_LEFT; m <= MOUSE_BUTTON_BACK; ++m) if (IsMouseButtonPressed(m)) return kMouse + m;
    return -1;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#895): a held button that broke something strikes nothing more until pressed again ----
#include "unit_helpers.h"

TEST(PlayerInputTextFocus, TypingKeepsKeyboardBindingsFromTheGame) { // (#1081: typing "/gm" in chat opened the map)
    PlayerInput_SetTextFocus(true);
    EXPECT_TRUE(PlayerInput_KeyBlocked(PlayerInput_Binding(Act::Map))) << "M while typing is a letter, not the map";
    EXPECT_TRUE(PlayerInput_KeyBlocked(KEY_ONE)) << "and digits aren't hotbar slots";
    EXPECT_FALSE(PlayerInput_KeyBlocked(PlayerInput_Binding(Act::Primary))) << "mouse buttons stay live";
    PlayerInput_SetTextFocus(false);
    EXPECT_FALSE(PlayerInput_KeyBlocked(PlayerInput_Binding(Act::Map))) << "no text field: the keys are the game's";
}

TEST(HoldLatch, HeldDiggingIsNeverInterrupted) {
    HoldLatch l;
    for (int f = 0; f < 120; ++f) EXPECT_TRUE(l.allows(true)) << "frame " << f;
}

TEST(HoldLatch, ABreakUsesUpTheHeldButtonUntilItIsPressedAgain) {
    HoldLatch l;
    EXPECT_TRUE(l.allows(true));
    l.consume(); // a block broke while the button was held
    for (int f = 0; f < 60; ++f) EXPECT_FALSE(l.allows(true)) << "still held: nothing behind it is struck (frame " << f << ")";
    EXPECT_FALSE(l.allows(false)) << "released";
    EXPECT_TRUE(l.allows(true)) << "pressed again: strikes again";
}

TEST(HoldLatch, AReleaseInTheSameFrameAsTheBreakClearsIt) {
    HoldLatch l;
    l.consume();
    EXPECT_FALSE(l.allows(false));
    EXPECT_FALSE(l.latched);
}

TEST(HoldLatch, EatingKeepsLmbDeadUntilAReleaseAfterIt) { // (#896) eating consumes it every frame
    HoldLatch l;
    for (int f = 0; f < 30; ++f) { l.consume(); EXPECT_FALSE(l.allows(f >= 10)) << "pressed while eating (frame " << f << ")"; }
    EXPECT_FALSE(l.allows(true)) << "eating stopped, LMB still held from before: still nothing";
    EXPECT_FALSE(l.allows(false));
    EXPECT_TRUE(l.allows(true)) << "pressed again after eating";
}
#endif
