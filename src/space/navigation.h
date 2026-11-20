#ifndef NAVIGATION_H
#define NAVIGATION_H

// Navigation (plan 2, milestone 9): a compass on the HUD, and with a Map in the inventory a
// minimap and a full-screen map (M) of the ground you have explored.
//
// Exploration: each planet or moon keeps a set of explored 256 m tiles (a cube-sphere grid, so
// the poles behave); tiles within 150 m of the player are marked as they go. Map images are drawn
// on the worker threads from the height field (the same biome colours as the far sphere,
// hillshaded, sea and rivers in blue); unexplored ground stays dark. Waypoints are placed on the
// full map and show on the compass and the minimap. Exploration, waypoints and the start are
// saved with the game.
//
// North is the body's spin axis (+Y, body-local); near a pole the compass shows "pole".

#include <raylib.h>
#include <cstdio>
#include <string>
#include <cstdint>
#include <unordered_set>

class Scene;

void Nav_Update(Scene& world, bool inputEnabled, bool gameInput); // once a frame: exploration, the map's input (inputEnabled), images; gameInput: the minimap zoom keys (#278)
// (#338) the map tiles a player at body-local `up` (unit) reveals on a body of radius R, added to `set` (the same tiles and radius
// the local map uses); how many were new. A server counts its players' exploring for their teams with it.
int Nav_RevealTiles(Vector3 up, float R, std::unordered_set<uint64_t>& set);
void Nav_DrawHUD(Scene& world);                   // compass (always), minimap (with a Map)
void Nav_DrawMap(Scene& world);                   // the full-screen map, when open (draw last)
bool Nav_MapOpen();
void Nav_ToggleMap();                             // M (needs a Map in the inventory)
void Nav_CloseMap();
void Nav_Reset();                                 // a new world: nothing explored
void Nav_SetStart(const std::string& body, Vector3 localDir); // the start marker
void Nav_Save(FILE* f);
bool Nav_LoadLine(const char* tok, FILE* f);      // true if it consumed a nav line

// Tests and HUD helpers.
int Nav_ExploredTiles(const std::string& body);
float Nav_Heading(Scene& world);                  // degrees clockwise from north; -1 near a pole
float Nav_Bearing(Scene& world, Vector3 localDir);// to a body-local direction, degrees (-1 none)
void Nav_AddWaypoint(const std::string& body, Vector3 localDir, const std::string& name, bool byPlayer = true); // byPlayer: counts for the quest
void Nav_RemoveWaypoint(const std::string& body, const std::string& name);
int Nav_WaypointCount();
bool Nav_MinimapReady();
bool Nav_MinimapShown(); // (#932) the minimap is on screen (it drew within the last 0.25 s)
bool Nav_MapImageReady();

#endif
