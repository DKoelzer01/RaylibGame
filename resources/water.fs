#version 330

// Sea surface: translucent deep blue, reflecting the sky at grazing angles (Fresnel), small ripples,
// a sun glint, the planet's own shadow at night, and the same distance fog as the terrain.
in vec3 fragPosition;
in vec3 fragNormal;
in vec3 fragTexPos;
in vec3 fragChunkPos;
in float fragDepth;            // depth to the sea floor (m)
flat in int fragKind;          // fluid kind (fluids.txt; 0 water)
in vec4 fragFlow;              // the current (#361): xyz direction (body-aligned), w speed m/s (0: still water)

uniform vec3 viewPos;
uniform vec3 sunDir;           // towards the sun
uniform vec3 skyColor;         // linear horizon sky colour (dark at night)
uniform vec3 chunkClipMax;     // chunk-local overlap clip (as for terrain)
uniform vec4 sunOccluder;      // planet: centre + radius blocking the sun (w = 0: none)
uniform vec3 fogColor;         // linear
uniform float fogDensity;
uniform float time;
uniform vec4 fluidColors[8];   // per kind (Fluid::BindKindUniforms): sRGB colour, opacity (a > 0: not water)
uniform float fluidGlows[8];   // per kind: emissive strength (lava)
uniform int fluidDebug;         // #695: 1 = water coloured by its speed, stripes running along the current

#include "shadow.glsl"

out vec4 finalColor;

void main()
{
    if (any(greaterThan(fragChunkPos, chunkClipMax))) discard;
    vec3 n = normalize(fragNormal);
    // Ripples: a few sine waves on the wrapped body-local position. Their wavelengths divide the
    // 20 m wrap period, so neighbouring chunks line up.
    const float k = 6.2831853 / 20.0;
    vec3 p = fragTexPos;
    float dist = length(viewPos - fragPosition);
    float rippleFade = 1.0 / (1.0 + dist / 60.0); // ripples finer than a pixel far away only alias
    // Rivers carry the ripples downstream (#361): two copies of the pattern advected along the current, each
    // restarting every 4 s half a period apart and cross-faded, so a speed that varies across the channel never
    // shears the pattern. Plus crests across the current running downstream, stronger where it is fast.
    vec3 fdir = fragFlow.xyz;
    float speed = fragFlow.w;
    const float T = 4.0;
    float ph0 = fract(time / T), ph1 = fract(time / T + 0.5);
    vec3 p0 = p - fdir * speed * ph0 * T, p1 = p - fdir * speed * ph1 * T;
    vec3 b0 = vec3(sin(p0.x * k * 7.0 + time * 1.3) + sin(p0.z * k * 5.0 - time * 1.1), 0.0,
                   sin(p0.z * k * 9.0 + time * 0.9) + sin(p0.y * k * 6.0 + time * 1.7));
    vec3 b1 = vec3(sin(p1.x * k * 7.0 + time * 1.3) + sin(p1.z * k * 5.0 - time * 1.1), 0.0,
                   sin(p1.z * k * 9.0 + time * 0.9) + sin(p1.y * k * 6.0 + time * 1.7));
    float w1 = abs(ph0 * 2.0 - 1.0);           // 1 when b0 restarts (ph0 wraps): b1 alone then, mid-cycle
    vec3 bump = mix(b0, b1, w1) * 0.025;
    if (speed > 0.01) {
        float s = dot(p, fdir) * k * 13.0 - time * speed * k * 13.0; // crests moving downstream at the current's speed
        bump += fdir * cos(s) * 0.03 * min(speed, 1.5);
    }
    bump *= rippleFade;
    n = normalize(n + bump - n * dot(bump, n));

    vec3 v = normalize(viewPos - fragPosition);
    float nv = abs(dot(n, v));                                    // also seen from below
    float fresnel = 0.02 + 0.98 * pow(1.0 - nv, 5.0);

    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        // Soft-edged: a penumbra 4% of the body's radius, widening behind it (a hard edge looked cut out).
        if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
    }
    sunVis *= sunShadowSimple(fragPosition, 1); // terrain, trees and builds shade the surface
    float light = sunVis * max(dot(n, sunDir), 0.0);
    vec3 deep = vec3(0.01, 0.05, 0.13) * (0.15 + 0.85 * light); // dark but clearly blue, also from orbit
    // Clear turquoise in the shallows, deep blue offshore.
    float clarity = exp(-fragDepth / 5.0);
    vec3 shallow = vec3(0.03, 0.22, 0.24) * (0.15 + 0.85 * light);
    vec3 color = mix(mix(deep, shallow, clarity), skyColor, fresnel);
    vec3 h = normalize(sunDir + v);
    float glint = pow(max(dot(n, h), 0.0), 800.0) * 1.2 * sunVis;
    color += vec3(1.0, 0.95, 0.85) * glint;
    // Translucent: see into the shallows looking down, mostly sky reflection at grazing angles
    // (Fresnel). Seen from below, the surface is a bright, mostly opaque ceiling. The scene's stencil
    // blends one water layer per pixel, so overlapping LOD copies don't double up.
    bool below = dot(normalize(fragNormal), v) < 0.0;
    float alpha = below ? 0.8 : clamp(max(mix(1.0, 0.3, clarity), fresnel) + glint, 0.0, 1.0);
    vec4 fluidColor = fluidColors[clamp(fragKind, 0, 7)];
    float fluidGlow = fluidGlows[clamp(fragKind, 0, 7)];
    if (fluidColor.a > 0.0) {
        // Another fluid kind (fluid.cpp, fluids.txt): its own colour and opacity, a sky reflection at grazing
        // angles, a faint rainbow sheen on dark ones (oil), a hard sun glint, and its glow (lava).
        vec3 base = pow(fluidColor.rgb, vec3(2.2));
        float dark = 1.0 - clamp(dot(base, vec3(0.333)) * 5.0, 0.0, 1.0);
        vec3 sheen = 0.5 + 0.5 * cos(6.2831853 * (vec3(0.0, 0.33, 0.67) + dot(fragTexPos, vec3(0.13, 0.07, 0.11)) + nv * 2.0));
        color = base * (0.3 + 0.7 * light) + skyColor * fresnel * 0.3 + sheen * 0.012 * light * dark
              + vec3(1.0, 0.95, 0.85) * glint * 1.5 + base * fluidGlow * 3.0;
        alpha = max(fluidColor.a, fresnel);
    }

    if (fluidDebug != 0) {
        // The debug view (#695): opaque, coloured by speed on a log scale (still deep blue, 0.5 m/s green, 2 m/s yellow,
        // >= 5 m/s red) so slow flow shows too; bands across the current, 2 a metre, drifting downstream at its speed.
        float v = fragFlow.w;
        float t = clamp(log2(1.0 + v / 0.05) / log2(101.0), 0.0, 1.0);
        vec3 ramp = t < 0.5 ? mix(vec3(0.02, 0.05, 0.45), vec3(0.05, 0.75, 0.15), t * 2.0)
                            : mix(vec3(0.95, 0.85, 0.05), vec3(0.95, 0.05, 0.02), t * 2.0 - 1.0);
        if (t >= 0.5 && t < 0.75) ramp = mix(vec3(0.05, 0.75, 0.15), vec3(0.95, 0.85, 0.05), (t - 0.5) * 4.0);
        if (v > 0.01) {
            vec3 fd = normalize(fragFlow.xyz);
            float band = fract(dot(fragChunkPos, fd) * 2.0 - time * v * 2.0);
            ramp *= band < 0.5 ? 1.0 : 0.6;
        }
        float lit = 0.55 + 0.45 * max(dot(normalize(fragNormal), sunDir), 0.0);
        finalColor = vec4(pow(ramp * lit, vec3(1.0 / 2.2)), 1.0);
        return;
    }
    float fog = 1.0 - exp(-fogDensity * dist);
    color = mix(color, fogColor, clamp(fog, 0.0, 1.0));
    finalColor = vec4(pow(max(color, vec3(0.0)), vec3(1.0 / 2.2)), alpha);
}
