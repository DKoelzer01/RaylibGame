#version 330

// Grass shading: colour by latitude (lush temperate green, drier and yellower toward the equator,
// olive toward the tundra) with per-tuft variation, darker roots, wrapped sun light, the planet's
// night-side shadow, sky ambient, the same haze as the terrain and a dithered fade with distance.
in vec3 fragPosition;
in vec3 fragUp;
in float fragHeight;
in float fragHash;
in float fragLatitude;
in vec3 fragClimate;   // temperature (C), moisture, 1 if known (Plan 5 B biomes)
in float fragBlade;

uniform vec3 viewPos;
uniform vec3 sunDir;        // towards the sun
uniform vec3 ambient;       // linear
uniform vec4 sunOccluder;   // planet: centre + radius blocking the sun
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;      // atmosphere: centre + outer radius (w = 0: fog the whole distance)
uniform float fadeEnd;      // fully gone here (m)
uniform float alphaToCoverage; // 1: the fade is written as alpha (GL_SAMPLE_ALPHA_TO_COVERAGE on)

// Sun shadows (two cascades, hardware PCF) and the headlamp.
#include "shadow.glsl"
#include "hemisphere.glsl"
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
    // Fade over the last 40% of the draw distance, without blending (so no sorting): under MSAA the
    // alpha is turned into sample coverage (alpha-to-coverage, #241); single-sampled, a dithered discard.
    float fade = clamp((fadeEnd - dist) / (fadeEnd * 0.4), 0.0, 1.0);
    fade = fade * fade * (3.0 - 2.0 * fade);
    if (alphaToCoverage > 0.5) {
        if (fade <= 0.0) discard;
    } else {
        vec2 cell = floor(gl_FragCoord.xy);
        float dither = fract(sin(dot(cell, vec2(12.9898, 78.233))) * 43758.5453);
        if (dither > fade) discard;
    }

    vec3 temperate = vec3(0.10, 0.30, 0.05);
    vec3 dry = vec3(0.40, 0.35, 0.11);   // straw
    vec3 cold = vec3(0.20, 0.24, 0.10);
    float lat = fragLatitude;
    vec3 albedo = climateFoliage(temperate, dry, cold, vec3(0.06, 0.26, 0.04), vec3(0.16, 0.22, 0.06), lat);
    albedo *= 0.8 + 0.4 * fragHash;
    albedo *= 0.85 + 0.3 * fragBlade;      // blades in a clump differ a little
    // Self-shadowed lower half (the clump and its neighbours shade it), sunlit tips.
    float lit = smoothstep(0.0, 0.65, fragHeight);
    albedo *= mix(0.28, 1.15, lit * lit * (3.0 - 2.0 * lit) * 0.5 + lit * 0.5);

    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        // Soft-edged: a penumbra 4% of the body's radius, widening behind it (a hard edge looked cut out).
        if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
    }
    // Thin blades: light from the ground's angle to the sun, wrapped so backlit grass still glows.
    float ndl = dot(normalize(fragUp), sunDir);
    float direct = sunVis * sunShadowSimple(fragPosition, 1) * clamp(ndl * 0.6 + 0.4, 0.0, 1.0) * smoothstep(-0.05, 0.1, ndl);
    vec3 color = albedo * (hemiAmbient(ambient, normalize(fragUp)) * 1.5 + vec3(1.0, 0.96, 0.88) * direct + lampLight(fragPosition, normalize(fragUp)));

    float fogDist = dist;
    if (fogShell.w > 0.0) {
        vec3 rd = (fragPosition - viewPos) / max(dist, 1e-3);
        vec3 oc = viewPos - fogShell.xyz;
        float b = dot(oc, rd), c = dot(oc, oc) - fogShell.w * fogShell.w, disc = b * b - c;
        float t0 = disc > 0.0 ? -b - sqrt(disc) : 0.0, t1 = disc > 0.0 ? -b + sqrt(disc) : 0.0;
        fogDist = max(0.0, min(t1, fogDist) - max(t0, 0.0));
    }
    float fog = 1.0 - exp(-fogDensity * fogDist);
    color = mix(color, fogColor, clamp(fog, 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), alphaToCoverage > 0.5 ? fade : 1.0);
}
