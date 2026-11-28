#ifndef BELTMODELS_H
#define BELTMODELS_H

// (#2053, belts #1735) Thin belts as the player sees them: every belt run of the plots near the camera (Plot::belts) as a
// low slab between two rails, 1/8 cell tall and half a cell wide, with chevrons that scroll at the belt's speed, and the
// items riding it (Factory::BeltItems, the factory's last tick, carried forward between ticks) as small cubes in their
// icon's colour with the icon on top. The boxes go through the instance batch (lit and shadowed like the machines' details);
// all of it under the plot's double-precision model (#271).

class Scene;

void Belts_Draw3D(Scene& world);      // before InstanceBatch_Flush
void Belts_DrawIcons(Scene& world);   // after it: the icon decals on the items drawn this frame

struct BeltDrawStats { int runs = 0, segments = 0, chevrons = 0, items = 0, icons = 0; };
const BeltDrawStats& Belts_LastStats(); // the last frame's (tests, captures)

#endif
