#ifndef WATERFALLS_H
#define WATERFALLS_H

// Waterfalls (#540): where a river steps down a cliff (RiverNetwork::waterfalls), a falling sheet of water from the lip
// down a throw arc to the pool, white water at its foot and mist and spray around it (scaled by the fall's power, #614). Near
// the camera the sheets are animated; out to 2.5 km static, fading in (#612). The falls within that range
// are drawn; each is built in the planet's frame (doubles) every frame, so it stays put on the spinning ground.

class Scene;

void Waterfalls_Update(Scene& world, float dt); // picks the falls near the camera, emits mist and spray
void Waterfalls_Draw3D(Scene& world);           // inside BeginMode3D, after the water (onDraw3D)
void Waterfalls_Reset();                        // a new or unloaded world
int Waterfalls_Visible();                       // tests: the falls drawn last frame

#endif
