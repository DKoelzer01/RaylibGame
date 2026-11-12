#version 330

// Parts: sRGB vertex colour, sun light with shadows (shadow.glsl) and the planet's night side, the
// headlamp, sky ambient (more from above), a soft specular, and the terrain's haze.
in vec3 fragPosition;
in vec3 fragNormal;
in vec4 fragColor;

uniform vec3 viewPos;
uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec4 sunOccluder;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;

#include "shadow.glsl"

out vec4 finalColor;

void main()
{
    vec3 albedo = pow(fragColor.rgb, vec3(2.2));
    vec3 n = normalize(fragNormal);
    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        // Soft-edged: a penumbra 4% of the body's radius, widening behind it (a hard edge looked cut out).
        if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
    }
    float ndl = max(dot(n, sunDir), 0.0);
    float sun = sunVis * ndl * (ndl > 0.0 ? sunShadowNormal(fragPosition, n, sunDir, 4) : 0.0);
    vec3 v = normalize(viewPos - fragPosition);
    float spec = pow(max(dot(n, normalize(sunDir + v)), 0.0), 32.0) * 0.12 * sun;
    vec3 color = albedo * (ambient * (1.4 + 0.3 * n.y) + vec3(1.0, 0.96, 0.88) * sun + lampLight(fragPosition, n)) + vec3(spec);
    float dist = length(viewPos - fragPosition), fogDist = dist;
    if (fogShell.w > 0.0) {
        vec3 rd = (fragPosition - viewPos) / max(dist, 1e-3);
        vec3 oc = viewPos - fogShell.xyz;
        float b = dot(oc, rd), cc = dot(oc, oc) - fogShell.w * fogShell.w, disc = b * b - cc;
        float t0 = disc > 0.0 ? -b - sqrt(disc) : 0.0, t1 = disc > 0.0 ? -b + sqrt(disc) : 0.0;
        fogDist = max(0.0, min(t1, fogDist) - max(t0, 0.0));
    }
    color = mix(color, fogColor, clamp(1.0 - exp(-fogDensity * fogDist), 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0);
}
