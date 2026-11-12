#ifndef BUILDTOOL_H
#define BUILDTOOL_H

class Scene;
class Planetoid;
class Plot;
struct Part;

// Placing parts on build plots (milestone 3). G toggles build mode; the wheel picks the part,
// R turns it about the plot's up axis and V about its east axis; RMB places (on a part face, or on
// the ground: snapping to a nearby plot's grid or starting a new plot), LMB removes a part and
// refunds it. Parts cost concrete or iron plate. A part placed against the ground clears the
// terrain inside its cells and fills the ground under them, so it neither floats nor sinks.
struct BuildToolState {
    bool active = false;
    int part = 0;           // PartType
    int orient = 0;         // 0..23
    int material = 0;       // BuildMaterial index for shapes
    float reach = 8.0f;
};
extern BuildToolState buildTool;

void BuildTool_Update(Scene& world, bool inputEnabled, bool usable);
void BuildTool_Draw3D(Scene& world);
void BuildTool_DrawHUD();
float BuildTool_BreakProgress();   // 0..1 while LMB is held on a part (the HUD bar under the crosshair)
bool BuildTool_AimingPart();       // the crosshair is on a built part in reach (LMB breaks it)
// Place a part on a plot (charging its cost unless creative or charge is false) and lay its
// foundation against the ground. Used by the tool and by scripted tests.
bool BuildTool_Place(Planetoid& body, Plot& plot, const Part& part, bool charge);
// Clear terrain inside a part's cells and fill small gaps under it (placing; a rocket landing).
void BuildTool_LayFoundation(Planetoid& body, Plot& plot, const Part& part);
// Machines: looking at one shows its state (factory snapshot); I inserts items, O takes outputs.
void Machine_Update(Scene& world, bool inputEnabled);
void Machine_DrawHUD();

#endif
