#ifndef BUILDTOOL_H
#define BUILDTOOL_H

#include <string>

class Scene;
class Planetoid;
class Plot;
struct Part;
struct Int3;
struct Vector3;

// Placing parts on build plots (milestone 3). G toggles build mode; the wheel picks the part,
// R turns it about the plot's up axis and V about its east axis; RMB places (on a part face, or on
// the ground: snapping to a nearby plot's grid or starting a new plot), LMB removes a part and
// refunds it. Parts cost concrete or iron plate. A part placed against the ground clears the
// terrain inside its cells and fills the ground under them, so it neither floats nor sinks.
struct BuildToolState {
    bool active = false;
    int part = 0;           // PartType
    int orient = 0;         // 0..23
    bool orientManual = false; // #196: R/V pressed: the orientation sticks; otherwise parts face the player automatically
    int autoPart = -1;          // the part the automatic facing was last reset for (another part: automatic again)
    int material = 0;       // BuildMaterial index for shapes
    float reach = 8.0f;
};
extern BuildToolState& buildTool; // the local PlayerSession's (#313)

struct PlayerInput;
void BuildTool_Update(Scene& world, const PlayerInput& in, bool usable);
void BuildTool_Draw3D(Scene& world);
void BuildTool_DrawHUD();
float BuildTool_BreakProgress();   // 0..1 while LMB is held on a part (the HUD bar under the crosshair)
bool BuildTool_AimingPart();       // the crosshair is on a built part in reach (LMB breaks it)
int BuildTool_AimOrient();        // the aimed placement's orientation (0..23), -1 without an aim (#196 test)
// The aimed placement for tests (#550/#551): the plot it goes into (null: a new plot at newOrigin, body-local) and its
// origin cell.
const Plot* BuildTool_AimPlacement(Vector3* newOrigin, Int3* origin);
// Place a part on a plot (charging its cost unless creative or charge is false) and lay its
// foundation against the ground. Used by the tool and by scripted tests.
bool BuildTool_Place(Planetoid& body, Plot& plot, const Part& part, bool charge);
// Clear terrain inside a part's cells and fill small gaps under it (placing; a rocket landing).
void BuildTool_LayFoundation(Planetoid& body, Plot& plot, const Part& part);
// Machines: looking at one shows its state (factory snapshot); I inserts items, O takes outputs.
void Machine_Update(Scene& world, const PlayerInput& in);
void Machine_DrawHUD();
bool IsMachineItem(const std::string& id); // an item that builds a machine part (chests don't take these)

#endif
