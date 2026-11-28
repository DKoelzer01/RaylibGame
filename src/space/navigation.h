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
#include <functional>
#include <map>
#include <unordered_set>
#include <vector>

class Scene;

// (#1886) explored map tiles by body name (a team's on a server, a player's own), and their text form for the team
// progression store: "tiles <body> <n>" then n hex keys, 16 to a line.
using NavTiles = std::map<std::string, std::unordered_set<uint64_t>>;
std::string Nav_TilesText(const NavTiles& tiles);
void Nav_TilesFromText(const std::string& text, NavTiles& into); // adds to `into` (a union)
// A client: the team's explored map added to the local one (the player keeps what they explored themselves).
void Nav_MergeTilesText(const std::string& text);
// (#1886) A team's shared map pins (the host keeps them per team, TeamQuests_SetPin; a client shows them beside its own
// waypoints). Text: "pin <body> <x> <y> <z> <name>" a line, the direction body-local; names are one word (Nav_PinName).
struct NavPin { std::string body, name; Vector3 dir; };
std::string Nav_PinsText(const std::vector<NavPin>& pins);
std::vector<NavPin> Nav_PinsFromText(const std::string& text);
std::string Nav_PinName(const std::string& raw); // printable, spaces as '_', at most 48 bytes ("" if nothing is left)
// A client: the team's pins as the host has them, replacing the team pins it showed ("" clears them, e.g. on leaving).
void Nav_SetTeamPins(const std::string& text);
int Nav_TeamPinCount();
// Set by the game while joined to a server: a pin the player places (a click on the map) or removes is sent to the team
// (true: sent; it shows as a team pin at once, and the host's store replaces it). Unset (single player, a host): local.
extern std::function<bool(bool add, const NavPin& pin)> gNavSharePin;

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
Rectangle Nav_MinimapRect(); // (#2072) the minimap's box on screen, rim included; zero size when it isn't shown
bool Nav_MapImageReady();

#endif
