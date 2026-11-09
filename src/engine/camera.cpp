#include "camera.h"

float cameraMoveSpeed = 30.0f;

FlyInput ReadFlyInput() {
    FlyInput in;
    if (IsKeyDown(KEY_W)) in.forward += 1.0f;
    if (IsKeyDown(KEY_S)) in.forward -= 1.0f;
    if (IsKeyDown(KEY_D)) in.right += 1.0f;
    if (IsKeyDown(KEY_A)) in.right -= 1.0f;
    if (IsKeyDown(KEY_SPACE)) in.up += 1.0f;
    if (IsKeyDown(KEY_LEFT_CONTROL)) in.up -= 1.0f;
    const float keyTurn = 1.5f;        // radians / second (arrow keys, roll)
    const float mouse = 0.003f;        // radians / pixel (already frame-rate independent)
    const float dt = GetFrameTime();
    const Vector2 m = GetMouseDelta();
    in.yaw = -m.x * mouse;
    in.pitch = -m.y * mouse;
    if (IsKeyDown(KEY_LEFT)) in.yaw += keyTurn * dt;
    if (IsKeyDown(KEY_RIGHT)) in.yaw -= keyTurn * dt;
    if (IsKeyDown(KEY_UP)) in.pitch += keyTurn * dt;
    if (IsKeyDown(KEY_DOWN)) in.pitch -= keyTurn * dt;
    if (IsKeyDown(KEY_Q) || IsKeyDown(KEY_Z)) in.roll -= keyTurn * dt; // roll left
    if (IsKeyDown(KEY_E) || IsKeyDown(KEY_C)) in.roll += keyTurn * dt; // roll right
    in.fast = (IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT) ? 2.0f : 1.0f) *
              (IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT) ? 2.0f : 1.0f);
    return in;
}

void ApplyFlyInput(Camera3D* camera, const FlyInput& in, float dt) {
    // The camera's own frame, re-orthonormalised every frame so it never drifts.
    Vector3 f = Vector3Normalize(Vector3Subtract(camera->target, camera->position));
    Vector3 u = Vector3Normalize(Vector3Subtract(camera->up, Vector3Scale(f, Vector3DotProduct(camera->up, f))));
    if (Vector3LengthSqr(u) < 0.5f) u = Vector3Normalize(Vector3Perpendicular(f));
    // Yaw about the camera's up, pitch about its right, roll about its forward: all relative to the
    // current view, so there is no gimbal flip at any orientation.
    const Quaternion q = QuaternionMultiply(QuaternionMultiply(
        QuaternionFromAxisAngle(f, in.roll),
        QuaternionFromAxisAngle(Vector3CrossProduct(f, u), in.pitch)),
        QuaternionFromAxisAngle(u, in.yaw));
    f = Vector3Normalize(Vector3RotateByQuaternion(f, q));
    u = Vector3Normalize(Vector3RotateByQuaternion(u, q));
    u = Vector3Normalize(Vector3Subtract(u, Vector3Scale(f, Vector3DotProduct(u, f))));
    const Vector3 r = Vector3CrossProduct(f, u);

    // Move along the camera's own axes: W/S forward, A/D right, Space/Ctrl up.
    const float step = cameraMoveSpeed * in.fast * dt;
    const Vector3 move = Vector3Add(Vector3Add(Vector3Scale(f, in.forward), Vector3Scale(r, in.right)), Vector3Scale(u, in.up));
    camera->position = Vector3Add(camera->position, Vector3Scale(move, step));
    camera->target = Vector3Add(camera->position, f);
    camera->up = u;
}

// Free-fly (dev / spectator) camera: fully camera-relative 6-DOF. Mouse / arrow keys yaw and
// pitch, Q / E (or Z / C) roll, WASD moves along the view and camera right, Space / Ctrl along the
// camera's up; Shift / Alt double the speed.
void customUpdateCamera(Camera3D* camera) {
    ApplyFlyInput(camera, ReadFlyInput(), GetFrameTime());
}
