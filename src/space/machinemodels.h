#ifndef MACHINEMODELS_H
#define MACHINEMODELS_H

// Machine models (Plan 5, after the animals): details built from boxes and cylinders over each machine's
// patterned casing (plots.cpp draws the casing itself): chimneys, hoppers, drums, tanks, pipes, doors,
// coils, drills, dishes - one set per machine type - and a working animation (spinning drums and fans,
// glowing furnace doors, a bobbing hammer, a humming fusion torus) while the factory reports it working.
// Drawn for machines within 32 m of the camera (the patterned casings carry them beyond).

class Scene;
void Machines_Draw3D(Scene& world);
int Machines_DrawnCount(); // tests: machines given models last frame

#endif
