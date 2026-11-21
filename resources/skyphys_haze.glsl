// #142: aerial perspective with the physical sky (on by default; RAYLIBGAME_SKY_PHYS=0: off). The haze keeps today's amount (fogDensity over the
// path inside fogShell: ~35 km visibility); only its colour changes, to the physical sky toward the fragment, just above
// the horizon (warm toward a low sun, blue away from it), pulled toward today's fogColor by skyPhysTint.
// #819: that colour comes from a small LUT built with the sky view (skyphys.cpp kFsHaze), one fetch a fragment: indexed
// by the cosine of the azimuth from the sun (u) and the view's height, packed toward the horizon (v = 0.5 + 0.5 sqrt|h|).
uniform float skyHaze;     // 0: off (switch off, under water, no air) .. 1; scene.cpp scales it by the air and the weather
uniform vec3 skyHazeUp;    // the local up at the camera
uniform vec3 skyHazeSun;   // toward the sun
uniform sampler2D skyHazeLut;
uniform float skyPhysTint; // 0 = the physical colour as it is .. 1 = today's fogColor

vec3 SkyPhysFog(vec3 fragPos, vec3 eye, vec3 fogC)
{
    if (skyHaze <= 0.0) return fogC;
    vec3 rd = normalize(fragPos - eye);
    float h = dot(rd, skyHazeUp), s = dot(skyHazeSun, skyHazeUp);
    vec3 vt = rd - skyHazeUp * h, st = skyHazeSun - skyHazeUp * s;
    float d = dot(vt, vt) * dot(st, st);
    float lvc = d > 1e-10 ? dot(vt, st) * inversesqrt(d) : 1.0;
    float v = 0.5 + 0.5 * sign(h) * sqrt(abs(h));
    vec3 phys = texture(skyHazeLut, vec2(lvc * 0.5 + 0.5, v)).rgb;
    return mix(fogC, mix(phys, fogC, skyPhysTint), skyHaze);
}

// (#1164) The same physical sky colour toward a direction rd (a reflection off water): the LUT covers every view height;
// rays under the horizon read just above it (water reflects no ground but the sky's rim). Off (no air, switch off): fallbackC.
vec3 SkyPhysToward(vec3 rd, vec3 fallbackC)
{
    if (skyHaze <= 0.0) return fallbackC;
    float h = max(dot(rd, skyHazeUp), 0.02), s = dot(skyHazeSun, skyHazeUp);
    vec3 vt = rd - skyHazeUp * dot(rd, skyHazeUp), st = skyHazeSun - skyHazeUp * s;
    float d = dot(vt, vt) * dot(st, st);
    float lvc = d > 1e-10 ? dot(vt, st) * inversesqrt(d) : 1.0;
    float v = 0.5 + 0.5 * sqrt(h);
    vec3 phys = texture(skyHazeLut, vec2(lvc * 0.5 + 0.5, v)).rgb;
    return mix(fallbackC, mix(phys, fallbackC, skyPhysTint), skyHaze);
}
