#ifndef TELESCOPE_H
#define TELESCOPE_H

class Scene;

// Hand telescope: T raises/lowers it (needs one in the inventory, on foot), the mouse wheel zooms
// 4x..64x. Narrows the camera's field of view, slows mouse look by the zoom, and tells the terrain
// LOD to refine what is in view as if it were that much closer (SetLodZoom).
struct TelescopeState {
    bool active = false;
    float zoom = 8.0f;
    float minZoom = 4.0f, maxZoom = 64.0f;
    float baseFov = 45.0f;   // the camera's field of view while lowered
};
extern TelescopeState telescope;

void Telescope_Update(Scene& world, bool inputEnabled, bool usable); // after the camera is placed
void Telescope_SetActive(Scene& world, bool on);
void Telescope_DrawHUD(Scene& world);

#endif
