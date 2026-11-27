#ifndef BELTTOOL_H
#define BELTTOOL_H

// (#2054, belts #1735) The belt tool: build mode with the conveyor item lays thin belts (Plot::belts) instead of cell
// conveyors (Ctrl/Alt+wheel switches between the two). Hold RMB on a plot's floor (the top of a part, or the ground by a
// plot) and drag: the run goes from where RMB went down to where it comes up, items flowing that way; one corner (R puts
// it on the other side), V reverses the flow. An end on a machine or chest snaps to that cell's side-face item port. LMB
// on a belt takes it up (refunded). Placing checks Conduit_CheckPlace: green can be placed, red says why.
//
// Units: belt points are conduit slots (1/8 cell, relative to the plot origin), x and z even, y a floor layer (cell * 8).

#include "conduit.h"
#include <functional>
#include <string>
#include <vector>

class Scene;
class Planetoid;
class Plot;
struct PlayerInput;

bool BeltTool_On();      // build mode, the conveyor selected, belts chosen
// Once a frame from BuildTool_Update. True while the crosshair is on a belt (or a run is being dragged): LMB is the belt
// tool's then, so the part break is skipped.
bool BeltTool_Update(Scene& world, const PlayerInput& in);
void BeltTool_Draw3D(Scene& world);
void BeltTool_DrawHUD();
void BeltTool_ForgetWorld(); // drop the aim and the drag (pointers into a world being unloaded)
void BeltTool_Cancel();      // drop a drag in progress (the build tool is unusable: a menu, death), without laying it

// ---- the parts below are pure (unit tests, scenarios) ----
// The run from a to b (one floor layer): straight when they line up on x or z, else one corner (xFirst: along x first).
// {} when a == b.
// firstAxis / lastAxis (0 x, 2 z): the axes of the first and last segment as routed before any snapping.
std::vector<CInt3> BeltTool_Route(CInt3 a, CInt3 b, bool xFirst, int* firstAxis = nullptr, int* lastAxis = nullptr);
// An end inside a machine (or chest) cell moves to the slot of the side face it crosses: `axis` (0 x, 2 z) and `sign`
// (+1: the belt travels toward +axis at this end) give the face; the cell is the machine's outermost along the travel.
// Any other point comes back as it is.
CInt3 BeltTool_SnapEnd(const Plot& plot, CInt3 p, int axis, int sign, bool start);
// The run planned from the drag's two points: route, snap both ends, route again (a straight run's free end follows the
// snapped one across).
std::vector<CInt3> BeltTool_Plan(const Plot& plot, CInt3 a, CInt3 b, bool xFirst);
// Every outer side-face centre of every machine cell, as the factory wires them (factory.cpp BuildBelts, GUESS in #2052).
std::vector<ConduitPortAt> BeltTool_Ports(const Plot& plot);
// The placement check against the plot's belts, parts and the terrain (`ground(cell)`: solid terrain at the cell's centre).
ConduitPlaceResult BeltTool_Check(const Plot& plot, const ConduitRun& run, const std::function<bool(CInt3 cell)>& ground);
std::function<bool(CInt3 cell)> BeltTool_Ground(const Planetoid& body, const Plot& plot);
int BeltTool_RunAt(const Plot& plot, CInt3 slot);   // the belt whose slots hold `slot` (or the slot above it), 0 none
int BeltTool_Cost(const ConduitRun& run);           // conveyor items (conduits.txt's basic belt)
// Places it (checked, charged unless creative or !charge; research gates the conveyor part); false and `why` when not.
bool BeltTool_Place(Planetoid& body, Plot& plot, ConduitRun run, bool charge, std::string* why);
bool BeltTool_Remove(Planetoid& body, Plot& plot, int id, bool refund);
// (scenarios, captures) A drag held from a to b on `plot` of `body`, planned now: the ghost and the HUD show it as if
// dragged; nothing places. A null plot ends it.
void BeltTool_BenchDrag(Planetoid* body, Plot* plot, CInt3 a, CInt3 b, bool xFirst = true);
const ConduitRun& BeltTool_Planned();               // the drag's run (empty points without a drag)
const ConduitPlaceResult& BeltTool_PlannedCheck();

#endif
