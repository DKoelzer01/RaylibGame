#version 330

// Grass shading: colour by latitude (lush temperate green, drier and yellower toward the equator,
// olive toward the tundra) with per-tuft variation, darker roots, wrapped sun light, the planet's
// night-side shadow, sky ambient, the same haze as the terrain and a dithered fade with distance.
in vec3 fragPosition;
in vec3 fragUp;
in float fragHeight;
in float fragHash;
in float fragLatitude;

uniform vec3 viewPos;
uniform vec3 sunDir;        // towards the sun
uniform vec3 ambient;       // linear
uniform vec4 sunOccluder;   // planet: centre + radius blocking the sun
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;      // atmosphere: centre + outer radius (w = 0: fog the whole distance)
uniform float fadeEnd;      // fully gone here (m)

// Shadows from the sun's shadow map (terrain, trees and builds), and the headlamp (a spot light).
uniform sampler2D shadowMap;
uniform mat4 lightSpaceMatrix;
uniform int shadowMapResolution;
uniform vec4 lampPos;     // render space; w = intensity (0: off)
uniform vec4 lampTarget;  // the beam's end; w unused
float sunShadow(vec3 p)
{
    if (shadowMapResolution <= 0) return 1.0;
    vec4 ls = lightSpaceMatrix * vec4(p, 1.0);
    vec3 q = ls.xyz / ls.w * 0.5 + 0.5;
    if (q.x < 0.0 || q.x > 1.0 || q.y < 0.0 || q.y > 1.0 || q.z > 1.0) return 1.0;
    vec2 texel = vec2(1.0 / float(shadowMapResolution));
    float lit = 0.0;
    for (int x = -1; x <= 1; x++)
        for (int y = -1; y <= 1; y++)
            lit += (q.z - 0.003 > texture(shadowMap, q.xy + texel * vec2(x, y)).r) ? 0.0 : 1.0;
    return lit / 9.0;
}
vec3 lampLight(vec3 p, vec3 n)
{
    if (lampPos.w <= 0.0) return vec3(0.0);
    vec3 beam = lampTarget.xyz - lampPos.xyz;
    float range = max(length(beam), 1.0);
    vec3 l = normalize(lampPos.xyz - p);
    float cone = smoothstep(0.82, 0.95, dot(-l, beam / range));
    float dl = length(lampPos.xyz - p) / range;
    return vec3(1.0, 0.95, 0.85) * lampPos.w * 3.0 * cone / (1.0 + dl * dl * 4.0) * max(dot(n, l), 0.35);
}
out vec4 finalColor;

void main()
{
    float dist = length(viewPos - fragPosition);
    // Dithered fade over the last 40% of the draw distance (no blending, so no sorting).
    float fade = clamp((fadeEnd - dist) / (fadeEnd * 0.4), 0.0, 1.0);
    vec2 cell = floor(gl_FragCoord.xy);
    float dither = fract(sin(dot(cell, vec2(12.9898, 78.233))) * 43758.5453);
    if (dither > fade) discard;

    vec3 temperate = vec3(0.10, 0.30, 0.05);
    vec3 dry = vec3(0.32, 0.33, 0.08);
    vec3 cold = vec3(0.20, 0.24, 0.10);
    float lat = fragLatitude;
    vec3 albedo = lat < 0.35 ? mix(dry, temperate, smoothstep(0.15, 0.35, lat))
                             : mix(temperate, cold, smoothstep(0.5, 0.7, lat));
    albedo *= 0.8 + 0.4 * fragHash;
    albedo *= mix(0.45, 1.15, fragHeight); // self-shadowed roots, sunlit tips

    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        if (c < 0.0 || (b < 0.0 && b * b > c)) sunVis = 0.0;
    }
    // Thin blades: light from the ground's angle to the sun, wrapped so backlit grass still glows.
    float ndl = dot(normalize(fragUp), sunDir);
    float direct = sunVis * sunShadow(fragPosition) * clamp(ndl * 0.6 + 0.4, 0.0, 1.0) * smoothstep(-0.05, 0.1, ndl);
    vec3 color = albedo * (ambient * 1.5 + vec3(1.0, 0.96, 0.88) * direct + lampLight(fragPosition, normalize(fragUp)));

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
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0);
}
