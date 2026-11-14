#ifndef MACHINEMODELS_H
#define MACHINEMODELS_H

#include <raylib.h>
#include <cstddef>

// Machine models (Plan 5, after the animals): details built from boxes and cylinders over each machine's
// patterned casing (plots.cpp draws the casing itself): chimneys, hoppers, drums, tanks, pipes, doors,
// coils, drills, dishes - one set per machine type - and a working animation (spinning drums and fans,
// glowing furnace doors, a bobbing hammer, a humming fusion torus) while the factory reports it working.
// Drawn for machines within 32 m of the camera (the patterned casings carry them beyond).

class Scene;
void Machines_Draw3D(Scene& world);
int Machines_DrawnCount(); // tests: machines given models last frame
// Where a working machine's smoke comes out (#220): its chimney tops, or the middle of its top (campfire: over the
// flames), body-local. Returns how many (at most max).
class Plot;
int Machine_SmokePoints(const Plot& plot, size_t part, Vector3* outBody, int max);
// The machine details outside the casing's cells (chimneys) along a ray in plot metres (#285): the nearest part hit
// within maxDist.
bool Machine_RaycastDetails(const Plot& plot, Vector3 origin, Vector3 dir, float maxDist, int& part, float& dist);

#endif
