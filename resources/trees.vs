#version 330

// Trees and bushes (trees.h): one prism-built variant mesh, instanced. Each instance stands the
// variant on the body's local up at its position (relative to the batch anchor), turns it by its
// hash and scales it; leaves sway a little in the wind, more toward the top.
layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexNormal;
layout(location = 2) in vec4 vertexColor;      // r: leaf (1) / bark (0), g: shade variation
layout(location = 10) in vec4 instPosScale;    // anchor-relative body-local position, scale
layout(location = 11) in vec4 instUpHash;      // body-local up, hash in [0, 1)

uniform mat4 matModel;     // anchor -> render space
uniform mat4 matViewProj;
uniform float time;
uniform vec3 viewPos;     // the main camera (also in the depth passes, so shadows match)
uniform float farForest;  // #745: 1 = far trees grow and darken a little so distant forests read (RAYLIBGAME_FAR_FOREST)
uniform float treeVariationOff; // #1979: 1 = every tree of a variant alike (RAYLIBGAME_TREE_VARIATION=0, A/B)
#include "wind.glsl"

out vec3 fragPosition;
out vec3 fragNormal;
out float fragLeaf;
out float fragShade;
out float fragStone;
out float fragHash;
out float fragLatitude;
out vec3 fragClimate;      // Plan 5 B: temperature (C), moisture, 1 if known (packed in the up vector's length)
out float fragHeight;
out float fragFar;         // #745: 0 near .. 1 far (0 with the switch off)
out float fragAo;          // (#933 layer 2) baked ambient occlusion (trees.cpp BakeFoliageAo): 1 open

void main()
{
    vec3 up = normalize(instUpHash.xyz);
    float h = instUpHash.w;
    vec3 ref = abs(up.y) < 0.9 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 t0 = normalize(cross(up, ref));
    vec3 b0 = cross(up, t0);
    float a = h * 6.2831853;
    vec3 t = t0 * cos(a) + b0 * sin(a);
    vec3 b = cross(t, up);   // right-handed (t, up, b): keeps the mesh winding (front faces CCW)
    float s = instPosScale.w;
    // #1979: no two trees alike. Sub-hashes of the instance hash give each tree its height stretch (0.85-1.2x), crown
    // width (0.92-1.08x) and a lean of up to 4 degrees, a shear that keeps the base on the ground (GUESS ranges, tuned in
    // playtests). The cull bound pads for the stretch (trees.cpp kVariationPad).
    float sy = 1.0, sxz = 1.0;
    vec3 leanDir = vec3(0.0);
    if (treeVariationOff < 0.5) {
        sy = 0.85 + 0.35 * fract(h * 131.7 + 0.13);
        sxz = 0.92 + 0.16 * fract(h * 251.3 + 0.71);
        float la = 6.2831853 * fract(h * 57.9 + 0.29);
        leanDir = (t0 * cos(la) + b0 * sin(la)) * 0.07 * fract(h * 389.1 + 0.53); // tan(4 deg) = 0.07
    }
    // #745: far trees grow up to +70% from 300 m to 1.4 km, so far forests keep their mass as they thin out.
    float far = farForest * smoothstep(300.0, 1400.0, length(vec3(matModel * vec4(instPosScale.xyz, 1.0)) - viewPos));
    s *= 1.0 + 0.7 * far; // (stronger, user 2026-11-19: +70% by 1.4 km)
    vec3 v = vertexPosition * vec3(sxz, sy, sxz);
    // Wind on leaves: a sway growing with height. With the wind field (#650) it grows with the local wind, and strong
    // wind (storms, bursts) leans the whole tree downwind, more toward the top, swaying with the gusts. No damage.
    float leaf = vertexColor.r;
    vec3 sway;
    if (windOn > 0.5) {
        vec4 w = windAt(vec3(matModel * vec4(instPosScale.xyz, 1.0)));
        vec3 wl = transpose(mat3(matModel)) * w.xyz; // render -> body-local directions
        wl -= up * dot(wl, up);
        float speed = length(wl);
        vec3 wdir = speed > 1e-3 ? wl / speed : t0;
        // (#1083) Fixed sway rates (a rate growing with the speed sped the cycle up in wind and jerked the phase at every
        // gust); the wind sets only the size, capped so storms stay subtle.
        float amp = min(0.02 + 0.004 * speed, 0.06);
        sway = (t0 * sin(time * 0.9 + h * 30.0) + b0 * sin(time * 0.7 + h * 17.0)) * amp * max(v.y, 0.0) * leaf;
        float lean = 0.05 * clamp((speed - 8.0) / 12.0, 0.0, 1.2) * (1.0 + 0.35 * sin(time * 1.2 + h * 9.0));
        sway += wdir * lean * max(v.y, 0.0) * max(v.y, 0.0);
    } else {
        sway = (t0 * sin(time * 0.9 + h * 30.0) + b0 * sin(time * 0.7 + h * 17.0)) * 0.04 * max(v.y, 0.0) * leaf;
    }
    vec3 local = instPosScale.xyz + (t * v.x + up * v.y + b * v.z + leanDir * max(v.y, 0.0)) * s + sway * s;
    vec3 vn = vertexNormal / vec3(sxz, sy, sxz); // the scale's inverse-transpose (the small lean shear is left out)
    vec3 n = t * vn.x + up * vn.y + b * vn.z;
    fragPosition = vec3(matModel * vec4(local, 1.0));
    fragNormal = normalize(mat3(matModel) * n);
    fragLeaf = leaf;
    fragShade = vertexColor.g;
    fragAo = vertexColor.a;
    fragStone = vertexColor.b;  // loose stones: grey rock
    fragHash = h;
    fragLatitude = abs(up.y);
    float cl = length(instUpHash.xyz);
    fragClimate = cl > 1.5 ? vec3((floor(cl) - 2.0) / 62.0 * 70.0 - 30.0, fract(cl), 1.0) : vec3(0.0); // codes from 2 (#254)
    fragHeight = v.y;
    fragFar = far;
    gl_Position = matViewProj * vec4(fragPosition, 1.0);
}
