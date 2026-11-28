#version 330

// Grass tufts (foliage.h): three shared tuft/clump meshes, instanced per chunk. Each instance places
// the tuft on the terrain (chunk-local position), sizes it, stands it along the body's local up,
// turns it by its hash and bends the blade tips in the wind.
layout(location = 0) in vec3 vertexPosition;   // tuft space: y = 0 root .. 1 tip
layout(location = 10) in vec4 instPosScale;    // chunk-local position, height (m) (half floats, #1603)
layout(location = 11) in vec4 instUpHash;      // body-local up (unit), hash x 2 - 1 (signed bytes, #1603)
layout(location = 12) in vec4 instClimate;     // (#1603) x: (temperature step + 1) / 255, 0 none; y: moisture; z: biome variant id / 255 (#2032) (bytes)

uniform mat4 matModel;     // chunk -> render space (Planetoid::chunkModel)
uniform mat4 matViewProj;
uniform float time;
uniform vec3 viewPos;      // render space (shared with the fragment shader)
uniform float fadeEnd;     // grass draw distance (m), shared with the fragment shader
#include "wind.glsl"

out vec3 fragPosition;     // render space
out vec3 fragUp;           // render space
out float fragHeight;      // 0 root .. 1 tip
out float fragHash;
out float fragLatitude;    // |body-local up.y|: 0 equator, 1 pole
out vec3 fragClimate;      // Plan 5 B: temperature (C), moisture, 1 if known (packed in the up vector's length)
flat out float fragBiome; // (#2032) the biome variant id (biometint.glsl), 0 if none
out float fragBlade;       // per-blade shade variation in [0, 1)

void main()
{
    vec3 root = vec3(matModel * vec4(instPosScale.xyz, 1.0));
    vec3 up = normalize(mat3(matModel) * instUpHash.xyz);
    float h = clamp(instUpHash.w * 0.5 + 0.5, 0.0, 0.999); // (#1603) stored as hash x 2 - 1
    // A tangent frame around up, turned by the hash.
    vec3 ref = abs(up.y) < 0.9 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 t0 = normalize(cross(up, ref));
    vec3 b0 = cross(up, t0);
    float a = h * 6.2831853;
    vec3 t = t0 * cos(a) + b0 * sin(a);
    vec3 b = cross(up, t);
    float s = mod(instPosScale.w, 10.0);  // the tens carry the variant (foliage.cpp)
    // Shrink toward the ground over the fade band too (#241), so the far edge sinks away rather than
    // ending in a line of full-height blades.
    s *= mix(1.0, 0.3, smoothstep(0.6 * fadeEnd, fadeEnd, length(root - viewPos)));
    float y = vertexPosition.y;
    vec3 sway;
    if (windOn > 0.5) {
        // The wind field (#650): bend by speed^2 toward the wind (bursts flatten the grass for a moment), a flutter whose
        // rate grows with the speed, and a per-blade phase so no two clumps move in step.
        vec4 w = windAt(root);
        vec3 wf = w.xyz - up * dot(w.xyz, up);
        float speed = length(wf);
        vec3 wdir = speed > 1e-3 ? wf / speed : t0;
        // (#1083, user: stronger wind should bend more, not sway faster; keep it minimal) The flutter runs at a fixed rate:
        // a rate that grew with the speed also made every gust jerk the phase (time x speed), faster the longer the game ran.
        // The wind sets only how far: the lean (speed^2, capped lower) and the flutter's size (capped), so storms stay
        // subtle; gusts change it slowly as the wind texture drifts over.
        float bend = windOld > 0.5 ? min(speed * speed / 144.0, 1.7) : min(speed * speed / 144.0, 1.0);   // 12 m/s = 1
        float flutter = windOld > 0.5 ? sin(time * (2.5 + 0.35 * speed) + h * 40.0) * (0.02 + 0.012 * speed) : sin(time * 2.0 + h * 40.0) * min(0.015 + 0.006 * speed, 0.06);
        vec3 side = cross(up, wdir);
        sway = (wdir * (0.06 + 0.32 * bend + flutter) + side * flutter * 0.5) * y * y;
        sway -= up * (0.22 * bend * bend) * y * y;                    // a bent blade is shorter, not just leaning
    } else {
        // No wind field (RAYLIBGAME_WIND=0): the old regular sway.
        float gust = sin(time * 1.3 + dot(root, vec3(0.11, 0.07, 0.13))) * 0.5 + 0.5;
        float flutter = sin(time * 4.0 + h * 40.0);
        vec3 windDir = normalize(t0 + b0 * 0.4);
        sway = windDir * (0.12 + 0.18 * gust + 0.04 * flutter) * y * y;
    }
    vec3 p = root + (t * vertexPosition.x + up * y + b * vertexPosition.z + sway) * s;
    fragPosition = p;
    fragUp = up;
    fragHeight = y;
    fragHash = h;
    fragLatitude = abs(normalize(instUpHash.xyz).y);
    float step1 = floor(instClimate.x * 255.0 + 0.5); // (#1603) temperature step + 1, 0: no climate (it rode in the up vector's length)
    fragClimate = step1 > 0.5 ? vec3((step1 - 1.0) / 62.0 * 70.0 - 30.0, instClimate.y, 1.0) : vec3(0.0); // (#254)
    fragBlade = fract(float(gl_VertexID / 3) * 0.618034 + h * 7.0);
    fragBiome = step1 > 0.5 ? floor(instClimate.z * 255.0 + 0.5) : 0.0; // (#2032) the variant id, its own byte
    gl_Position = matViewProj * vec4(p, 1.0);
    // Thin clumps with distance: past 30% of the draw distance drop a growing share of the blades (down to a third by
    // 85%), the same blades for all three corners so a dropped blade is clipped away whole.
    float keep = mix(1.0, 0.33, smoothstep(0.3 * fadeEnd, 0.85 * fadeEnd, length(root - viewPos)));
    if (fract(float(gl_VertexID / 3) * 0.618034) > keep) gl_Position = vec4(2.0, 2.0, 2.0, 1.0);
}
