#version 330

// Tree shading: bark brown, leaves coloured by latitude (dry olive near the equator, lush green in
// the temperate band, dark blue-green toward the tundra) with per-block and per-tree variation,
// sunlight wrapped a little through the leaves, the planet's night side, sky ambient, haze and a
// dithered fade at the draw distance.
in vec3 fragPosition;
in vec3 fragNormal;
in float fragLeaf;
in float fragShade;
in float fragHash;
in float fragLatitude;
in float fragHeight;

uniform vec3 viewPos;
uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec4 sunOccluder;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;
uniform float fadeEnd;

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
    float fade = clamp((fadeEnd - dist) / (fadeEnd * 0.25), 0.0, 1.0);
    vec2 cell = floor(gl_FragCoord.xy);
    if (fract(sin(dot(cell, vec2(12.9898, 78.233))) * 43758.5453) > fade) discard;

    vec3 n = normalize(fragNormal);
    float lat = fragLatitude;
    vec3 dry = vec3(0.22, 0.27, 0.07), temperate = vec3(0.09, 0.27, 0.05), cold = vec3(0.05, 0.17, 0.08);
    vec3 leafC = lat < 0.35 ? mix(dry, temperate, smoothstep(0.15, 0.35, lat)) : mix(temperate, cold, smoothstep(0.45, 0.7, lat));
    leafC *= 0.75 + 0.5 * fragShade * (0.8 + 0.4 * fragHash);
    vec3 bark = vec3(0.20, 0.13, 0.08) * (0.8 + 0.4 * fragShade);
    vec3 albedo = mix(bark, leafC, fragLeaf);

    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        // Soft-edged: a penumbra 4% of the body's radius, widening behind it (a hard edge looked cut out).
        if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
    }
    float ndl = dot(n, sunDir);
    float wrap = mix(0.0, 0.3, fragLeaf);
    float direct = sunVis * sunShadow(fragPosition) * clamp((ndl + wrap) / (1.0 + wrap), 0.0, 1.0);
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
