#ifndef camera_h
#define camera_h

#include <raylib.h>
#include "raymath.h"
#include "rcamera.h"

extern float cameraMoveSpeed; // units per second (Shift/Alt multiply it)

// One frame of free-fly input: movement along the camera's own axes (-1..1 each) and rotation in
// radians about them. Gathered from the keyboard/mouse, or built by tests.
struct FlyInput {
    float forward = 0.0f, right = 0.0f, up = 0.0f;
    float yaw = 0.0f, pitch = 0.0f, roll = 0.0f; // roll > 0 = roll right
    float fast = 1.0f;                           // speed multiplier
};
FlyInput ReadFlyInput();                                           // keyboard + mouse
void ApplyFlyInput(Camera3D* camera, const FlyInput& in, float dt); // camera-relative 6-DOF
void customUpdateCamera(Camera3D* camera);                         // ReadFlyInput + ApplyFlyInput

#endif
