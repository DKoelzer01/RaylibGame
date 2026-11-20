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

#endif
