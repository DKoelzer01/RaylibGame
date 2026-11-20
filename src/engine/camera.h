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
// gravityUp null: camera-relative 6-DOF (space, near-zero g). Given (#232, near a body): yaw about that up, pitch clamped,
// no roll, the camera's up eased onto it; WASD moves level (the view only picks the heading), Space / Ctrl along up.
void ApplyFlyInput(Camera3D* camera, const FlyInput& in, float dt, const Vector3* gravityUp = nullptr);
void customUpdateCamera(Camera3D* camera);                         // ReadFlyInput + ApplyFlyInput (with gFlyGravityUp)
// The local up at the fly camera where gravity is felt (#232: set by the player in fly mode each frame), or zero for the
// free 6-DOF camera.
extern Vector3 gFlyGravityUp;

#endif
