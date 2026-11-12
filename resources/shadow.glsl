// Sun shadows, shared by lighting.fs, trees.fs, grass.fs and parts.fs (pulled in by `#include` in
// HotShader). Two cascades: a near box (~3 cm texels) and a wide one (~20 cm).
//
// The maps are depth textures with hardware comparison (GL_COMPARE_REF_TO_TEXTURE, linear filter):
// one texture() call returns the bilinear-filtered compare of the four texels around the point, so a
// 2x2 of taps half a texel apart gives the same soft ~2-texel edge the old shader built from 16 manual
// compares, in 4 fetches. The near cascade fades into the far one over the last 8% of its box (no
// visible seam where the sharp shadows end), and the far one fades to lit at its own edge.
uniform sampler2DShadow shadowMap;      // far cascade
uniform sampler2DShadow shadowMapNear;  // near cascade
uniform mat4 lightSpaceMatrix;
uniform mat4 lightSpaceNear;
uniform int shadowMapResolution;        // 0: shadows off
uniform vec2 shadowTexel;               // metres per texel: near, far
uniform float shadowDepthRange;         // metres the light's depth spans

float shadowPcf(sampler2DShadow map, vec3 q, int taps)
{
    if (taps == 1) return texture(map, q);
    float o = 0.5 / float(shadowMapResolution);
    return 0.25 * (texture(map, vec3(q.xy + vec2(-o, -o), q.z)) + texture(map, vec3(q.xy + vec2(o, -o), q.z))
                 + texture(map, vec3(q.xy + vec2(-o, o), q.z)) + texture(map, vec3(q.xy + vec2(o, o), q.z)));
}

// One cascade at p (already offset along the normal by the caller), depth bias in metres.
// edge: 1 well inside the box, falling to 0 at its border (0 also beyond its depth range).
float shadowCascade(sampler2DShadow map, mat4 m, vec3 p, float biasM, int taps, out float edge)
{
    vec4 ls = m * vec4(p, 1.0);
    vec3 q = ls.xyz / ls.w * 0.5 + 0.5;
    vec2 d = min(q.xy, 1.0 - q.xy);
    edge = clamp((min(d.x, d.y) - 0.01) / 0.07, 0.0, 1.0);
    if (edge <= 0.0 || q.z > 1.0) { edge = 0.0; return 1.0; }
    q.z -= biasM / shadowDepthRange;
    return shadowPcf(map, q, taps);
}

// Both cascades: positions and biases per cascade (the normal offset and bias scale with each texel).
float sunShadowCascades(vec3 pNear, vec3 pFar, float biasNear, float biasFar, int taps)
{
    if (shadowMapResolution <= 0) return 1.0;
    float eNear, eFar;
    float sNear = shadowCascade(shadowMapNear, lightSpaceNear, pNear, biasNear, taps, eNear);
    if (eNear >= 1.0) return sNear;
    float sFar = mix(1.0, shadowCascade(shadowMap, lightSpaceMatrix, pFar, biasFar, taps, eFar), eFar);
    return mix(sFar, sNear, eNear);
}

// Foliage and other surfaces without a reliable normal: a depth bias of two texels.
float sunShadowSimple(vec3 p, int taps)
{
    return sunShadowCascades(p, p, 2.0 * shadowTexel.x, 2.0 * shadowTexel.y, taps);
}

// Surfaces with a normal n (lit from l): the lookup moved about a texel along the normal (no
// self-shadowing on flat-shaded triangles) plus a slope-scaled bias of a texel or two.
float sunShadowNormal(vec3 p, vec3 n, vec3 l, int taps)
{
    float ndl = clamp(dot(n, l), 0.0, 1.0);
    float slope = min(sqrt(max(1.0 - ndl * ndl, 0.0)) / max(ndl, 0.2), 4.0);
    float k = 1.5 - 0.5 * ndl;
    return sunShadowCascades(p + n * shadowTexel.x * k, p + n * shadowTexel.y * k,
                             shadowTexel.x * (0.5 + slope), shadowTexel.y * (0.5 + slope), taps);
}

// The headlamp (a spot light): render-space position (w = intensity, 0 off) and the beam's end.
uniform vec4 lampPos;
uniform vec4 lampTarget;
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
