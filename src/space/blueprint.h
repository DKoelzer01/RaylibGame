#ifndef BLUEPRINT_H
#define BLUEPRINT_H

// Blueprints (design doc section 7 follow-up): a plot's parts saved to a text file and placed again
// anywhere. Presets ship in resources/blueprints/, the player's own go to blueprints/.
// In build mode: B saves the plot you look at; N cycles blueprints (a ghost outline at the aim, its
// parts list on screen; N past the last turns it off); RMB places every part you can afford, bottom
// up, and reports what is still missing (place again later to fill in the rest).

#include <cstdio>
#include <string>
#include <vector>
#include "plots.h"

class Planetoid;
class Scene;

struct Blueprint {
    std::string name, path;
    std::vector<Part> parts;   // origins relative to the blueprint's lowest corner
};

const std::vector<Blueprint>& Blueprint_List();     // presets and the player's, loaded on first use
void Blueprint_Reload();
bool Blueprint_SaveFile(const Plot& plot, const std::string& path, const std::string& name);
bool Blueprint_LoadFile(const std::string& path, Blueprint& out);

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
