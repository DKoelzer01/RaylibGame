#include "camera.h"

void customUpdateCamera(Camera3D* camera)
{
    // Camera movement speed
    float cameraMoveSpeed = 0.5f;
    float cameraRotationSpeed = 0.1f;
    float cameraPanSpeed = 0.05f;
    float cameraMouseMoveSensitivity = 0.003f;

    // Camera mode
    const int mode = CAMERA_CUSTOM; // CAMERA_FIRST_PERSON, CAMERA_THIRD_PERSON, CAMERA_FREE, etc.

    // Camera flags
    const bool lockView = (mode == CAMERA_FIRST_PERSON || mode == CAMERA_THIRD_PERSON);
    const bool rotateAroundTarget = (mode == CAMERA_THIRD_PERSON);
    const bool rotateUp = (mode == CAMERA_FREE || mode == CAMERA_THIRD_PERSON);

    // Mouse position delta
    Vector2 mousePositionDelta = GetMouseDelta();
    bool moveInWorldPlane = (mode != CAMERA_FIRST_PERSON);

    // Camera reference
    Camera3D& cameraRef = *camera;

    UpdateCamera(&cameraRef, CAMERA_CUSTOM);

    if (IsKeyDown(KEY_LEFT_SHIFT) || IsKeyDown(KEY_RIGHT_SHIFT))
    {
        cameraMoveSpeed *= 2.0f; // Increase speed when shift is held
    }
    if (IsKeyDown(KEY_LEFT_ALT) || IsKeyDown(KEY_RIGHT_ALT))
    {
        cameraMoveSpeed *= 2.0f; // Increase speed speed when alt is held
    }

    // Camera rotation
    if (IsKeyDown(KEY_DOWN)) CameraPitch(camera, -cameraRotationSpeed, lockView, rotateAroundTarget, rotateUp);
    if (IsKeyDown(KEY_UP)) CameraPitch(camera, cameraRotationSpeed, lockView, rotateAroundTarget, rotateUp);
    if (IsKeyDown(KEY_RIGHT)) CameraYaw(camera, -cameraRotationSpeed, rotateAroundTarget);
    if (IsKeyDown(KEY_LEFT)) CameraYaw(camera, cameraRotationSpeed, rotateAroundTarget);
    // if (IsKeyDown(KEY_Q)) CameraRoll(camera, -cameraRotationSpeed);
    // if (IsKeyDown(KEY_E)) CameraRoll(camera, cameraRotationSpeed);

    // Camera movement
    // Camera pan (for CAMERA_FREE)
    if ((mode == CAMERA_FREE) && (IsMouseButtonDown(MOUSE_BUTTON_MIDDLE)))
    {
        const Vector2 mouseDelta = GetMouseDelta();
        if (mouseDelta.x > 0.0f) CameraMoveRight(camera, cameraPanSpeed, moveInWorldPlane);
        if (mouseDelta.x < 0.0f) CameraMoveRight(camera, -cameraPanSpeed, moveInWorldPlane);
        if (mouseDelta.y > 0.0f) CameraMoveUp(camera, -cameraPanSpeed);
        if (mouseDelta.y < 0.0f) CameraMoveUp(camera, cameraPanSpeed);
    }
    else
    {
        // Mouse support
        CameraYaw(camera, -mousePositionDelta.x*cameraMouseMoveSensitivity, rotateAroundTarget);
        CameraPitch(camera, -mousePositionDelta.y*cameraMouseMoveSensitivity, lockView, rotateAroundTarget, rotateUp);
    }

    // Keyboard support
    if (IsKeyDown(KEY_W)) CameraMoveForward(camera, cameraMoveSpeed, moveInWorldPlane);
    if (IsKeyDown(KEY_A)) CameraMoveRight(camera, -cameraMoveSpeed, moveInWorldPlane);
    if (IsKeyDown(KEY_S)) CameraMoveForward(camera, -cameraMoveSpeed, moveInWorldPlane);
    if (IsKeyDown(KEY_D)) CameraMoveRight(camera, cameraMoveSpeed, moveInWorldPlane);
    if (IsKeyDown(KEY_SPACE)) CameraMoveUp(camera, cameraMoveSpeed);
    if (IsKeyDown(KEY_LEFT_CONTROL)) CameraMoveUp(camera, -cameraMoveSpeed);
}