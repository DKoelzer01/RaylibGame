#ifndef SCENETARGET_H
#define SCENETARGET_H

// Scene target with reverse-Z (#135): the 3D world draws into an offscreen framebuffer with a 32-bit float
// depth buffer (+ stencil for the water), the same sample count as the window, using glClipControl
// (GL 4.5 / ARB_clip_control) with a [0, 1] depth range, an infinite far plane with depth 1 at the near
// plane and 0 at infinity, and GL_GREATER. Float depth is close to logarithmic that way: centimetres up
// close and metres at 100 km+ in one range, so the far / mid / near passes with their own depth ranges
// and the depth clears between them are not needed (Scene::drawScene draws one pass when this is on).
// Then the colour is resolved to a single-sample texture and drawn to the window by a fullscreen pass (works
// whatever the window's format and sample count); the HUD and the hands draw on the window as before.
//
// Off (old passes) without clip control or the framebuffer, or with RAYLIBGAME_REVZ=0.

#include <raylib.h>

bool SceneTarget_Available();
// Binds the target (created / resized to the window), sets reverse-Z depth state and clears it.
bool SceneTarget_Begin(Color clear);   // false: not available (or the target failed): draw the old passes
// Replaces the current projection (after BeginMode3D) with the infinite reversed one of the same field
// of view, aspect and near plane.
void SceneTarget_ApplyProjection();
// Restores the standard depth state and copies the colour to the window.
void SceneTarget_End();
void SceneTarget_Unload();

#endif
