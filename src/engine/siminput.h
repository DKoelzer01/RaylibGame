#ifndef SIMINPUT_H
#define SIMINPUT_H

#include <raylib.h>

// Input for the fixed simulation steps (#315). The gameplay block runs 0..N times a frame (world.h SimClock), so raylib's
// per-frame edges (IsKeyPressed, IsMouseButtonPressed) would be lost on a frame with no step and seen twice on a frame
// with two. SimInput_BeginFrame latches them once a frame, the first step after that sees them, and SimInput_EndStep
// clears them; mouse motion and the wheel add up the same way. Held state (IsKeyDown, IsMouseButtonDown) needs no
// latch. Code that runs per frame (menus, UI, the free-fly camera) keeps reading raylib directly.
void SimInput_BeginFrame();
void SimInput_EndStep();
bool SimKeyPressed(int key);
bool SimKeyReleased(int key);
bool SimMousePressed(int button);
bool SimMouseReleased(int button);
Vector2 SimMouseDelta();
float SimMouseWheel();

#endif
