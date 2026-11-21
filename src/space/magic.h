#pragma once

// Magic (#1001, the magic doc; steps #1248-#1259). Everything magic is behind RAYLIBGAME_MAGIC (off by default) until the
// first playable slice. Step 1 (#1248): the mana field of the start planet and a debug overlay.

#include "manafield.h"

class Scene;

bool Magic_Enabled(); // RAYLIBGAME_MAGIC=1

// The start planet's mana field (built on first use from its terrain seed; nullptr without a planet with air, or with
// magic off).
const ManaField* Magic_Field(Scene& world);
// What the place adds at a body-local direction on the start planet now: height, climate, biome, wind and storms (the
// weather is a pure function of place and world time, so clients agree), depth = m under the surface.
ManaEnv Magic_EnvAt(Scene& world, Vector3 dir, float depth = 0.0f);

// The debug overlay (RAYLIBGAME_MAGIC=1 and RAYLIBGAME_MANA_DEBUG=1): ley lines within ~300 m of the camera coloured by kind,
// nodes as beacons, and a HUD line with the densities at the camera and the nearest node.
void Magic_DrawDebug3D(Scene& world);
void Magic_DrawDebugHUD(Scene& world);

// Magic sight's look (#1249): ribbons within `range` m of the camera, nodes within 1.5 x range.
void Magic_DrawSight3D(Scene& world, float range, bool throughRock);
// The node whose glow is under the crosshair within `range` m (of those magic sight last drew): false if none.
bool Magic_SightNodeAhead(Scene& world, float range, ManaNode& out, float& dist);
