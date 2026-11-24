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
// Low-sun tint (#276): 1 - the sun light's colour (warm and dim near the horizon); unset = 0 = white sun.
uniform vec3 sunFade;
vec3 sunTint() { return vec3(1.0) - clamp(sunFade, 0.0, 1.0); }

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

// Local point lights (#139, locallights.h): lamps, campfires, glowing machines. Per 16 px screen tile a count and up to
// 15 light indices; per light two texels (render-space position + range, linear colour + intensity). The light falls to
// zero at its range (a smooth square falloff), with a little wrap so faces turned slightly away still catch some.
uniform sampler2D localLightTex;
uniform isampler2D localTileTex;
uniform ivec4 localTileInfo;   // tiles across, tile size (px), lights this frame (0: none), ints per tile
vec3 localLights(vec3 p, vec3 n)
{
    if (localTileInfo.z <= 0) return vec3(0.0);
    ivec2 t = ivec2(gl_FragCoord.xy) / localTileInfo.y;
    int base = t.x * localTileInfo.w;
    int count = texelFetch(localTileTex, ivec2(base, t.y), 0).r;
    vec3 sum = vec3(0.0);
    for (int k = 0; k < count; ++k) {
        int i = texelFetch(localTileTex, ivec2(base + 1 + k, t.y), 0).r;
        vec4 a = texelFetch(localLightTex, ivec2(2 * i, 0), 0);
        vec3 d = a.xyz - p;
        float dist = length(d);
        if (dist >= a.w) continue;
        vec4 c = texelFetch(localLightTex, ivec2(2 * i + 1, 0), 0);
        float f = 1.0 - dist / a.w;
        float wrap = clamp((dot(n, d / max(dist, 1e-3)) + 0.25) / 1.25, 0.0, 1.0);
        sum += c.rgb * c.a * f * f * wrap;
    }
    return sum;
}

// (#1739) Static lights (torches, lamps, fires) through the light-level grid (lightgrid.h): a 128 m window of 1 m cells
// round the camera, body-local, in a wrap-around 3D texture. Sampled half a metre out along the normal (in the open cell in
// front of the surface, not the solid one behind it). Texels: the light's colour (rgb, its largest channel 1) and its
// level (a: 8 per metre of reach left, / 255), so colour x (level / 128)^2 falls off as the point lights'
// (1 - d / range)^2 did; info.w is the gain (0: off). Fades out at the window's edge.
uniform sampler3D lightGridTex;
uniform mat4 lightGridFromRender;   // render space -> texture coordinates
uniform vec4 lightGridInfo;         // the window's centre (texture coordinates), gain
vec3 gridLight(vec3 p, vec3 n)
{
    if (lightGridInfo.w <= 0.0) return vec3(0.0);
    vec3 q = (lightGridFromRender * vec4(p + n * 0.5, 1.0)).xyz;
    vec3 off = abs(q - lightGridInfo.xyz);
    float edge = max(off.x, max(off.y, off.z));
    if (edge > 0.49) return vec3(0.0);
    vec4 t = texture(lightGridTex, q);
    float level = t.a * (255.0 / 128.0);
    return t.rgb * level * level * lightGridInfo.w * (1.0 - smoothstep(0.42, 0.49, edge));
}

// A held torch (#218), a point light: render-space position and range (w, 0 off); colour and intensity (a).
// Terrain gets it through lights[3] (lighting.fs), with the same falloff.
uniform vec4 torchPos;
uniform vec4 torchColor;
vec3 torchLight(vec3 p, vec3 n)
{
    if (torchPos.w <= 0.0) return vec3(0.0);
    vec3 l = normalize(torchPos.xyz - p);
    float dl = length(torchPos.xyz - p) / torchPos.w;
    return torchColor.rgb * torchColor.a * 3.0 / (1.0 + dl * dl * 4.0) * max(dot(n, l), 0.35);
}

// The headlamp (a spot light): render-space position (w = intensity, 0 off) and the beam's end.
uniform vec4 lampPos;
uniform vec4 lampTarget;
// lampLight() also adds the held torch and the local lights: every shader that calls it gets them all.
vec3 lampLight(vec3 p, vec3 n)
{
    vec3 other = torchLight(p, n) + localLights(p, n) + gridLight(p, n); // with the headlamp off too (it returned early before: no torch on grass, #218)
    if (lampPos.w <= 0.0) return other;
    vec3 beam = lampTarget.xyz - lampPos.xyz;
    float range = max(length(beam), 1.0);
    vec3 l = normalize(lampPos.xyz - p);
    float cone = smoothstep(0.82, 0.95, dot(-l, beam / range));
    float dl = length(lampPos.xyz - p) / range;
    return vec3(1.0, 0.95, 0.85) * lampPos.w * 3.0 * cone / (1.0 + dl * dl * 4.0) * max(dot(n, l), 0.35) + other;
}
