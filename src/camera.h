#ifndef camera_h
#define camera_h

#include <raylib.h>
#include "raymath.h"
#include "rcamera.h"

extern float cameraMoveSpeed; // units per second (Shift/Alt multiply it)
void customUpdateCamera(Camera3D* camera);

#endif