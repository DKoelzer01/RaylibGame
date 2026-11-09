#include "camera.h"

// Free-fly camera. Speeds are per second and scaled by frame time, so movement no longer
// depends on frame rate.
float cameraMoveSpeed = 30.0f;

void customUpdateCamera(Camera3D* camera)
{
    const float dt = GetFrameTime();
    float moveSpeed = cameraMoveSpeed; // units / second
    const float keyTurnSpeed = 1.5f;   // radians / second (arrow keys)
    const float mouseSensitivity = 0.003f; // radians / pixel (already frame-rate independent)

    // Free-fly: rotate around the camera itself, allow looking straight up/down.
    const bool lockView = false;
    const bool rotateAroundTarget = false;
    const bool rotateUp = false;
    const bool moveInWorldPlane = true;    // WASD stays level, Space/Ctrl move vertically (unchanged behaviour)

    if (IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT)) moveSpeed *= 2.0f;
    if (IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT)) moveSpeed *= 2.0f;
    const float step = moveSpeed * dt;
    const float turn = keyTurnSpeed * dt;

    if (IsKeyDown(KEY_DOWN))  CameraPitch(camera, -turn, lockView, rotateAroundTarget, rotateUp);
    if (IsKeyDown(KEY_UP))    CameraPitch(camera,  turn, lockView, rotateAroundTarget, rotateUp);
    if (IsKeyDown(KEY_RIGHT)) CameraYaw(camera, -turn, rotateAroundTarget);
    if (IsKeyDown(KEY_LEFT))  CameraYaw(camera,  turn, rotateAroundTarget);

    const Vector2 mouse = GetMouseDelta();
    CameraYaw(camera, -mouse.x * mouseSensitivity, rotateAroundTarget);
    CameraPitch(camera, -mouse.y * mouseSensitivity, lockView, rotateAroundTarget, rotateUp);

    if (IsKeyDown(KEY_W)) CameraMoveForward(camera,  step, moveInWorldPlane);
    if (IsKeyDown(KEY_S)) CameraMoveForward(camera, -step, moveInWorldPlane);
    if (IsKeyDown(KEY_A)) CameraMoveRight(camera, -step, moveInWorldPlane);
    if (IsKeyDown(KEY_D)) CameraMoveRight(camera,  step, moveInWorldPlane);
    if (IsKeyDown(KEY_SPACE))        CameraMoveUp(camera,  step);
    if (IsKeyDown(KEY_LEFT_CONTROL)) CameraMoveUp(camera, -step);
}
