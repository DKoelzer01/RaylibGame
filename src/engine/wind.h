#ifndef WIND_H
#define WIND_H

// The wind field (#650): one field that grass, bushes and trees bend in and gameplay (wind turbines, later sails, fire,
// snow drift) reads.
//
//  - Mean wind: the weather's regional wind (set each frame by the game, Wind_SetFrame) times a daily cycle (calm at
//    dawn, building through the afternoon, easing at night) and the weather (storms stronger and gustier).
//  - Gusts: a tileable noise texture (512 m period) in the plane at an anchor near the player, sampled at the position
//    minus the distance the air has travelled, so gust fronts roll downwind across the grass. Two parts: rolling waves
//    (~60 m fronts with ~15 m puffs) and rare burst cells that briefly flatten the grass (more often in storms), plus
//    a few degrees of direction swing.
//  - The CPU (Wind_At) and the shaders (Wind_Bind: grass.vs, trees.vs) read the same texture data and formula.
//  - Gameplay extras on the CPU only: terrain exposure (ridges more, hollows and lee less) and a log height profile.
//  - RAYLIBGAME_WIND=0: off (the old regular sway).

#include <raylib.h>
#include "dvec3.h"

class Planetoid;

// Each frame, from the game (weather.cpp): the body the player is on, its regional mean wind (body-local, m/s along
// the surface), storminess (0 calm .. 1 storm) and the daily factor; dt advances the drift (sim time: stops when paused).
void Wind_SetFrame(const Planetoid* body, Vector3 playerLocal, Vector3 meanLocal, float storm, float daily, float dt);
// The daily factor from the sun (body-local directions): ~0.45 at dawn, ~1.3 mid-afternoon, ~0.7 at night.
float Wind_DailyFactor(Vector3 upLocal, Vector3 sunLocal, Vector3 spinAxisLocal);

// The wind at a body-local point (m/s, body-local, along the surface). heightAboveGround (m) scales it by a log profile
// (10 m = 1); exposure from the terrain applies when `withTerrain` (a few height-map lookups).
Vector3 Wind_At(const Planetoid& body, Vector3 local, float heightAboveGround = 10.0f, bool withTerrain = true);
float Wind_Gust01(const Planetoid& body, Vector3 local); // 0 calm .. 1 a burst here (tests, sound)
bool Wind_Enabled();

// Bind the wind texture and uniforms on a foliage shader (grass, trees): windTex on texture slot `slot`.
void Wind_Bind(Shader s, DVec3 renderOrigin, int slot);
// (#888 A/B) RAYLIBGAME_WIND_888=0: the pre-#888 foliage wind (full weather wind at the ground, faster flutter, harder gusts)
inline bool Wind_Old888() { static const bool old = getenv("RAYLIBGAME_WIND_888") && getenv("RAYLIBGAME_WIND_888")[0] == '0'; return old; }
void Wind_Unload();

#endif
