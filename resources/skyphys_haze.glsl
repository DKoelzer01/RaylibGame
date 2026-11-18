// #142: aerial perspective with the physical sky (on by default; RAYLIBGAME_SKY_PHYS=0: off). The haze keeps today's amount (fogDensity over the
// path inside fogShell: ~35 km visibility); only its colour changes, to the physical sky toward the fragment, just above
// the horizon (warm toward a low sun, blue away from it), pulled toward today's fogColor by skyPhysTint.
uniform float skyHaze;     // 0: off (switch off, under water, no air) .. 1; scene.cpp scales it by the air and the weather
uniform vec3 skyHazeUp;    // the local up at the camera
uniform vec3 skyHazeSun;   // toward the sun
#include "skyphys.glsl"

vec3 SkyPhysFog(vec3 fragPos, vec3 eye, vec3 fogC)
{
    if (skyHaze <= 0.0) return fogC;
    vec3 rd = normalize(fragPos - eye);
    vec3 phys = SkyPhysSky(rd, skyHazeUp, skyHazeSun, true);
    return mix(fogC, mix(phys, fogC, skyPhysTint), skyHaze);
}
