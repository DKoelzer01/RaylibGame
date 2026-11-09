#include "camera.h"

float cameraMoveSpeed = 30.0f;

// Free-fly (dev / spectator) camera: fully camera-relative 6-DOF.
// - Mouse yaw turns about the camera's own up axis, pitch about its own right axis, and pitch also
//   rotates the up vector, so there is no gimbal flip when looking straight up or down.
// - WASD moves along the view direction / camera right; Space / Ctrl along the camera's up.
// - Z / C roll. Speeds are per second and scaled by frame time.
void customUpdateCamera(Camera3D* camera)
{
    const float dt = GetFrameTime();
    float moveSpeed = cameraMoveSpeed; // units / second
    const float keyTurnSpeed = 1.5f;   // radians / second (arrow keys, roll)
    const float mouseSensitivity = 0.003f; // radians / pixel (already frame-rate independent)

    const bool lockView = false;           // no pitch clamp: we carry the up vector instead
    const bool rotateAroundTarget = false; // rotate about the eye, not the target
    const bool rotateUp = true;            // pitch rotates camera.up too -> no flip at the poles
    const bool moveInWorldPlane = false;   // move along where you look, not a world plane

    if (IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT)) moveSpeed *= 2.0f;
    if (IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT)) moveSpeed *= 2.0f;
    const float step = moveSpeed * dt;
    const float turn = keyTurnSpeed * dt;

    if (IsKeyDown(KEY_DOWN))  CameraPitch(camera, -turn, lockView, rotateAroundTarget, rotateUp);
    if (IsKeyDown(KEY_UP))    CameraPitch(camera,  turn, lockView, rotateAroundTarget, rotateUp);
    if (IsKeyDown(KEY_RIGHT)) CameraYaw(camera, -turn, rotateAroundTarget);
    if (IsKeyDown(KEY_LEFT))  CameraYaw(camera,  turn, rotateAroundTarget);
    if (IsKeyDown(KEY_Z))     CameraRoll(camera, -turn);
    if (IsKeyDown(KEY_C))     CameraRoll(camera,  turn);

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
