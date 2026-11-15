#include "siminput.h"

namespace {
constexpr int kKeys = 400, kButtons = 8;
bool keyPressed[kKeys], keyReleased[kKeys], buttonPressed[kButtons], buttonReleased[kButtons];
Vector2 mouseDelta{ 0, 0 };
float wheel = 0.0f;
}

void SimInput_BeginFrame() {
    for (int k = 1; k < kKeys; ++k) {
        if (IsKeyPressed(k)) keyPressed[k] = true;
        if (IsKeyReleased(k)) keyReleased[k] = true;
    }
    for (int b = 0; b < kButtons; ++b) {
        if (IsMouseButtonPressed(b)) buttonPressed[b] = true;
        if (IsMouseButtonReleased(b)) buttonReleased[b] = true;
    }
    const Vector2 d = GetMouseDelta();
    mouseDelta.x += d.x; mouseDelta.y += d.y;
    wheel += GetMouseWheelMove();
}

void SimInput_EndStep() {
    for (int k = 0; k < kKeys; ++k) keyPressed[k] = keyReleased[k] = false;
    for (int b = 0; b < kButtons; ++b) buttonPressed[b] = buttonReleased[b] = false;
    mouseDelta = { 0, 0 };
    wheel = 0.0f;
}

bool SimKeyPressed(int key) { return key > 0 && key < kKeys && keyPressed[key]; }
bool SimKeyReleased(int key) { return key > 0 && key < kKeys && keyReleased[key]; }
bool SimMousePressed(int button) { return button >= 0 && button < kButtons && buttonPressed[button]; }
bool SimMouseReleased(int button) { return button >= 0 && button < kButtons && buttonReleased[button]; }
Vector2 SimMouseDelta() { return mouseDelta; }
float SimMouseWheel() { return wheel; }
