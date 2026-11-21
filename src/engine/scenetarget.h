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
// Contact shadows (#137): the sun direction in view space and the strength (0: none) for this frame's composite.
void SceneTarget_SetContact(Vector3 sunView, float strength);
// Hi-Z occlusion culling (#134): this frame's view-projection (render space, after SceneTarget_ApplyProjection); at
// SceneTarget_End the resolved depth becomes a min-depth pyramid. SceneTarget_HiZ gives the latest pyramid and the
// view-projection it was drawn with (false: none yet, or RAYLIBGAME_HIZ=0).
void SceneTarget_SetHiZView(Matrix viewProj);
bool SceneTarget_HiZ(unsigned& tex, int& w, int& h, int& levels, Matrix& viewProj);
long long SceneTarget_Frame();    // frames ended so far (the one being drawn has this number)
long long SceneTarget_HiZFrame(); // the frame the current pyramid was built in (-1: none)
// Restores the standard depth state and copies the colour to the window.
void SceneTarget_End();
// (#1164) The opaque scene so far (colour, depth) resolved for screen-space reflections; drawing continues into the target.
bool SceneTarget_ResolveOpaque(unsigned& colorTex, unsigned& depthTex, int& w, int& h);
void SceneTarget_Unload();
// #131 the tonemap curve (0 knee - the default, today's look -, 1 ACES, 2 clip): Settings sets it; RAYLIBGAME_TONEMAP wins.
void SceneTarget_SetTonemap(int curve);
int SceneTarget_Tonemap();

#endif
