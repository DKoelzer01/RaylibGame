#ifndef MULTIBLOCKMODELS_H
#define MULTIBLOCKMODELS_H

// Multiblock animation (#1429 PR 6): over the plot mesh's casings, a status lamp on each controller's front (dim red:
// unformed; green: formed, pulsing while it works; blinking amber: stalled), a lamp in the hatch kind's colour on each
// formed hatch's outer face (lit while the structure works), and a seam glow along the formed structure's outer edges.
// The state is the plot's (Multiblock_Looks) and the factory's working and stalled views, so a joined client sees it too.
// Structures within 32 m of the camera, as the machine details.

class Scene;
void MultiblockModels_Draw3D(Scene& world);
int MultiblockModels_DrawnCount(); // tests: structures (controllers) drawn last frame

#endif
