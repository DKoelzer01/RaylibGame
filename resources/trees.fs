#version 330

// Tree shading: bark brown, leaves coloured by latitude (dry olive near the equator, lush green in
// the temperate band, dark blue-green toward the tundra) with per-block and per-tree variation,
// sunlight wrapped a little through the leaves, the planet's night side, sky ambient, haze and a
// dithered fade at the draw distance.
in vec3 fragPosition;
in vec3 fragNormal;
in float fragLeaf;
in float fragShade;
in float fragStone;
in float fragHash;
in float fragLatitude;
in vec3 fragClimate;   // temperature (C), moisture, 1 if known (Plan 5 B biomes)
in float fragHeight;

uniform vec3 viewPos;
uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec4 sunOccluder;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;
uniform float fadeEnd;
uniform vec3 barkTint;   // per species (trees.h WoodSpeciesInfo); 1 for bushes and stones
uniform vec3 leafTint;

// Sun shadows (two cascades, hardware PCF) and the headlamp.
#include "shadow.glsl"
out vec4 finalColor;
// Foliage colour by climate (Plan 5 B): temperate green, straw-olive where warm and dry, deep in the
// jungle, dark blue-green in the cold, olive in swamps. Without a climate: by latitude as before.
vec3 climateFoliage(vec3 temperate, vec3 dry, vec3 cold, vec3 lush, vec3 swampC, float lat)
{
    if (fragClimate.z < 0.5)
        return lat < 0.35 ? mix(dry, temperate, smoothstep(0.15, 0.35, lat)) : mix(temperate, cold, smoothstep(0.45, 0.7, lat));
    float t = fragClimate.x, m = fragClimate.y;
    // Temperate land has its own gentle gradient, as the ground (#254): cooler blue-green, warmer yellower, wetter deeper.
    vec3 c = temperate * mix(vec3(0.9, 1.0, 1.02), vec3(1.1, 1.05, 0.86), smoothstep(8.0, 17.0, t));
    c *= mix(vec3(1.08, 1.04, 0.9), vec3(0.88, 1.0, 0.9), smoothstep(0.3, 0.68, m));
    c = mix(c, dry, smoothstep(0.66, 0.5, m) * smoothstep(16.0, 20.0, t));
    c = mix(c, lush, smoothstep(18.0, 23.0, t) * smoothstep(0.56, 0.7, m));
    c = mix(c, swampC, smoothstep(0.68, 0.8, m) * smoothstep(10.0, 18.0, t));
    c = mix(c, cold, smoothstep(9.0, 3.0, t));
    return c;
}

void main()
{
    float dist = length(viewPos - fragPosition);
    float fade = clamp((fadeEnd - dist) / (fadeEnd * 0.25), 0.0, 1.0);
    vec2 cell = floor(gl_FragCoord.xy);
    if (fract(sin(dot(cell, vec2(12.9898, 78.233))) * 43758.5453) > fade) discard;

    vec3 n = normalize(fragNormal);
    float lat = fragLatitude;
    vec3 dry = vec3(0.22, 0.27, 0.07), temperate = vec3(0.09, 0.27, 0.05), cold = vec3(0.05, 0.17, 0.08);
    vec3 leafC = climateFoliage(temperate, dry, cold, vec3(0.05, 0.24, 0.04), vec3(0.12, 0.18, 0.05), lat);
    leafC *= 0.75 + 0.5 * fragShade * (0.8 + 0.4 * fragHash);
    leafC = min(leafC * leafTint, vec3(0.9));
    vec3 bark = min(vec3(0.20, 0.13, 0.08) * barkTint, vec3(0.85)) * (0.8 + 0.4 * fragShade);
    vec3 albedo = mix(bark, leafC, fragLeaf);
    albedo = mix(albedo, vec3(0.36, 0.35, 0.33) * (0.7 + 0.6 * fragShade), fragStone); // loose stones

    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        // Soft-edged: a penumbra 4% of the body's radius, widening behind it (a hard edge looked cut out).
        if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
    }
    float ndl = dot(n, sunDir);
    float wrap = mix(0.0, 0.3, fragLeaf);
    // Normal-offset + slope-biased lookup (#188: the fixed-bias one striped the steep canopy faces at low sun). Leaves
    // are lit from both sides, so their normal faces the sun for the lookup.
    vec3 ns = (fragLeaf > 0.5 && ndl < 0.0) ? -n : n;
    float direct = sunVis * sunShadowNormal(fragPosition, ns, sunDir, 4) * clamp((ndl + wrap) / (1.0 + wrap), 0.0, 1.0);
    float sky = 0.7 + 0.3 * n.y; // a little more sky from above
    vec3 color = albedo * (ambient * 1.6 * sky + vec3(1.0, 0.96, 0.88) * direct + lampLight(fragPosition, n));

    float fogDist = dist;
    if (fogShell.w > 0.0) {
        vec3 rd = (fragPosition - viewPos) / max(dist, 1e-3);
        vec3 oc = viewPos - fogShell.xyz;
        float b = dot(oc, rd), c = dot(oc, oc) - fogShell.w * fogShell.w, disc = b * b - c;
        float t0 = disc > 0.0 ? -b - sqrt(disc) : 0.0, t1 = disc > 0.0 ? -b + sqrt(disc) : 0.0;
        fogDist = max(0.0, min(t1, fogDist) - max(t0, 0.0));
    }
    color = mix(color, fogColor, clamp(1.0 - exp(-fogDensity * fogDist), 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0);
}
