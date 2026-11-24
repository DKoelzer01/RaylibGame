#ifndef LIGHTLEVEL_H
#define LIGHTLEVEL_H

// (#1739 PR 2) The game's static light-level grid (engine lightgrid.h) for the body the player is on: placed torches and
// lamp blocks are its sources (no count limit, unlike the point lights), terrain and solid building blocks its walls.
// Rebuilt from the world (nothing saved); terrain edits and plot changes relight only what they touch, a budget of cells
// a frame. The GPU window round the camera is uploaded through LocalLights_UploadGrid. RAYLIBGAME_LIGHT_GRID=0: off
// (torches and lamps go back to point lights).

#include "lightgrid.h"
#include <raylib.h>

class Scene;
class Planetoid;

bool LightLevel_On();
// Each frame on the main thread, after the camera is placed and before the scene is drawn.
void LightLevel_Update(Scene& world, bool inGame);
// The grid's level at a body-local point (0 where it has none, or on another body). Creatures and crops read this later.
LightCell LightLevel_At(const Planetoid& body, Vector3 local);
// Stats for tests and the profiler: sources, bricks, queued work, cells relit last frame, ms spent last frame.
struct LightLevelStats { int sources = 0; int bricks = 0; bool busy = false; int uploaded = 0; float ms = 0.0f; };
LightLevelStats LightLevel_Stats();

#endif
