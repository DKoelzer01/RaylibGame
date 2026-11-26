#ifndef BLUEPRINT_H
#define BLUEPRINT_H

// Blueprints (design doc section 7 follow-up): a plot's parts saved to a text file and placed again
// anywhere. Presets ship in resources/blueprints/, the player's own go to blueprints/.
// In build mode: B saves the plot you look at; N cycles blueprints (a ghost outline at the aim, its
// parts list on screen; N past the last turns it off); RMB places every part you can afford, bottom
// up, and reports what is still missing (place again later to fill in the rest).
// (#1886) Each team also has a blueprint library, kept in the world (the team progression store's "blueprints"
// section); the blueprints/ folder stays the player's personal library. B saves to both; N cycles the team's, then the
// player's own (one already in the team's is shown once), then the presets.

#include <cstdio>
#include <string>
#include <vector>
#include "plots.h"

class Planetoid;
class Scene;

struct Blueprint {
    std::string name, path;    // path: "" for a team library blueprint
    std::vector<Part> parts;   // origins relative to the blueprint's lowest corner
    int team = 0;              // (#1886) from this team's library (0: a file)
    std::string text;          // the .bp text it was read from
};

const std::vector<Blueprint>& Blueprint_List();     // the team's, the player's and the presets, loaded on first use
void Blueprint_Reload();
bool Blueprint_SaveFile(const Plot& plot, const std::string& path, const std::string& name);
bool Blueprint_LoadFile(const std::string& path, Blueprint& out);
std::string Blueprint_Text(const Plot& plot, const std::string& name); // the .bp text of a plot's parts ("" when it has none)
bool Blueprint_ParseText(const std::string& text, Blueprint& out);    // keeps out.name when the text has no name line

// ---- Team libraries (#1886) ----
void Blueprint_AddToTeam(int team, const std::string& text);  // a .bp text (the same text twice is kept once)
int Blueprint_TeamCount(int team);
std::string Blueprint_TeamText(int team);                      // the library as store text: "blueprint" .bp lines "blueprint_end"
void Blueprint_SetTeamText(int team, const std::string& text); // replaces the library ("" empties it)
void Blueprint_ApplyTeamText(int team, const std::string& text); // a client: its team's library from the host (N shows it)
void Blueprint_ResetTeams();
// The team whose library N shows and B adds to: a client's from the host's store, else the local player's team.
int Blueprint_LocalTeam();

const Blueprint* Blueprint_Active();                // the one being placed, or null
void Blueprint_Cycle();                             // N
void Blueprint_SetActive(int index);                // -1: none
// Place the active blueprint with its corner at `origin` (plot cell): parts that fit and can be paid
// for. Returns the number placed; `missing` lists what ran out. (#1018) With `ghosts`, every part that
// could not be paid for stays as a ghost there, filled in later from the inventory (Blueprint_UpdateGhosts).
int Blueprint_Place(Planetoid& body, Plot& plot, const Blueprint& bp, Int3 origin, bool charge, std::string& missing, bool ghosts = false);

// ---- Ghosts (#1018, Factorio-style): parts of a placed blueprint not built yet ----
// Drawn see-through wherever they are; within kGhostReach of the player they build themselves from the inventory, a few
// a second (later: drones). A ghost whose cell is taken by anything else is dropped. Saved with the world.
constexpr float kGhostReach = 8.0f;   // m from the player's feet
int Blueprint_GhostCount();
int Blueprint_GhostsOn(const Plot& plot, const Blueprint& bp, Int3 origin);     // the blueprint's ghosts at that corner
int Blueprint_RemoveGhosts(const Plot& plot, const Blueprint& bp, Int3 origin); // cancels them; how many went
// One update a sim step: builds at most one ghost within reach (paced to 4 a second); the number built.
int Blueprint_UpdateGhosts(Planetoid* body, Vector3 feetLocal, double now);
struct DVec3;
void Blueprint_DrawGhosts(Planetoid& body, const DVec3& renderOrigin, Vector3 cameraLocal); // inside BeginMode3D
void Blueprint_ResetGhosts();
void Blueprint_SaveGhosts(FILE* f);
bool Blueprint_LoadGhostLine(const char* tok, FILE* f);
// Build-mode hooks: B / N keys (true when a key was used), the ghost and the parts list.
bool Blueprint_Keys(Planetoid* body, Plot* lookingAt);
void Blueprint_DrawGhost(Planetoid& body, const Plot& frame, Int3 origin, bool ok, const DVec3& renderOrigin);
void Blueprint_DrawHUD();

#endif
