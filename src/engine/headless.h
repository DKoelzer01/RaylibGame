#ifndef HEADLESS_H
#define HEADLESS_H

// Headless mode (#311, first slice): the simulation runs without drawing. `--headless` sets it; main.cpp then opens a
// hidden 320x180 window (raylib, the loaders and the GPU terrain still need a GL context: removing it is the next slice)
// and skips every draw: the 3D scene, the viewmodel, the HUD and ImGui. Gameplay steps (SimClock), streaming,
// animals, machines and the fluid simulation keep running, so a dedicated server can host from this mode. Scene::drawScene
// returns before any GPU work; it counts the frames it was asked to draw (gHeadlessSkippedDraws) for --test-headless.
//
// Still to split (see #311): engine_core / engine_render libraries, GL uploads inside the world update (chunk meshes,
// tree VBOs, plot sections, ring VAOs), the window-less clock (GetTime), and HUD code mixed into src/space logic.

inline bool gHeadless = false;
inline long long gHeadlessSkippedDraws = 0;
inline bool Headless() { return gHeadless; }

// No-render (#2130): a scenario that checks no picture skips the drawing and nothing else. Unlike --headless (a dedicated
// server: no sound, no local player, the frame rate capped), the game runs as a player's would: the camera and the
// player's view still update and the pre-draw uploads still run, so streaming and WorldIdle() behave the same. Only the
// shadow cascades, the main pass (with the Hi-Z pyramid and onDraw3D), the viewmodel and the HUD are left out. Set by
// RAYLIBGAME_NO_RENDER=1 (the scenario harness sets it for scenarios that take no shot and keep no recording).
inline bool gNoRender = false;
inline long long gNoRenderSkippedDraws = 0;
inline bool NoRender() { return gNoRender; }

#endif
