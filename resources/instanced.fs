#version 330

// InstanceBatch: sun light with shadows (shadow.glsl), the headlamp, sky ambient (more from above), a soft
// specular, the terrain's haze; glowing instances (furnace doors, coils) unlit.
in vec3 fragPosition;
in vec3 fragNormal;
in vec4 fragColor;
in float fragGlow;
in float fragSky;

uniform vec3 viewPos;
uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;
uniform vec4 sunOccluder; // the nearest planet blocks the sun (night side), as for parts (#284)

#include "shadow.glsl"
#include "hemisphere.glsl"

out vec4 finalColor;

void main()
{
    vec3 albedo = pow(fragColor.rgb, vec3(2.2));
    vec3 color;
    if (fragGlow > 0.5) {
        color = albedo * 1.6;
    } else {
        vec3 n = normalize(fragNormal);
        if (!gl_FrontFacing) n = -n;
        // The planet's night side and how buried the instance is (#284): as parts.fs, so machine details match
        // the casings and parts around them.
        float sunVis = 1.0;
        if (sunOccluder.w > 0.0) {
            vec3 oc = fragPosition - sunOccluder.xyz;
            float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
            if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
        }
        float ndl = max(dot(n, sunDir), 0.0);
        float sun = ndl > 0.0 ? sunVis * ndl * sunShadowNormal(fragPosition, n, sunDir, 1) * fragSky : 0.0;
        vec3 v = normalize(viewPos - fragPosition);
        float spec = pow(max(dot(n, normalize(sunDir + v)), 0.0), 32.0) * 0.1 * sun;
        color = albedo * ((hemiOn() ? hemiAmbient(ambient, n) * 1.4 : ambient * (1.4 + 0.3 * n.y)) * fragSky + vec3(1.0, 0.96, 0.88) * sun + lampLight(fragPosition, n)) + vec3(spec);
    }
    float dist = length(viewPos - fragPosition), fogDist = dist;
    if (fogShell.w > 0.0) {
        vec3 rd = (fragPosition - viewPos) / max(dist, 1e-3);
        vec3 oc = viewPos - fogShell.xyz;
        float b = dot(oc, rd), cc = dot(oc, oc) - fogShell.w * fogShell.w, disc = b * b - cc;
        float t0 = disc > 0.0 ? -b - sqrt(disc) : 0.0, t1 = disc > 0.0 ? -b + sqrt(disc) : 0.0;
        fogDist = max(0.0, min(t1, fogDist) - max(t0, 0.0));
    }
    color = mix(color, fogColor, clamp(1.0 - exp(-fogDensity * fogDist), 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), fragColor.a);
}
