#version 330

in vec3 fragPosition;   // world space
in vec2 fragTexCoord;    // chunks: material weights 4-5 (sand, gravel)
in vec2 fragTexCoord2;   // chunks: material weights 6-7 (grass, snow)
in vec4 fragColor;
in vec3 fragNormal;     // world space
in vec3 fragTexPos;     // body-aligned, wrapped to the texture repeat (texturing)
in vec3 fragLocalNormal;
in vec3 fragChunkPos;   // chunk-local (overlap clipping)
in float fragSky;       // 1 = open to the sky, 0 = buried (caves, deep pits)
in vec2 fragCanopy;     // planets: forest canopy cover and latitude (painted in beyond the real trees)
in float fragOreKind;   // planets: which ore, (1 + kind) / 8 (0 none)
in float fragMacro;     // planets: slow brightness / warmth variation (70-250 m), per vertex (lighting.vs)

// Material maps: bound by DrawMesh to the units of the matching raylib material slots.
uniform sampler2D colorTexture;    // MATERIAL_MAP_ALBEDO
uniform sampler2D normalTexture;   // MATERIAL_MAP_NORMAL (OpenGL convention, tangent space)
uniform sampler2D roughAoTexture;  // MATERIAL_MAP_ROUGHNESS: R = roughness, G = ambient occlusion
uniform sampler2DArray terrainLayers; // planet palette: dirt, sand, gravel, grass, snow, ice, stone (terrainlayers.h)
uniform sampler2DArray terrainDetail; // its detail maps: normal XY, roughness, AO (unit 12)

flat in vec3 fragClipMax;     // chunk-local overlap clip (lighting.vs: a uniform or the arena record)
uniform mat4 matModel;          // body-local -> world (DrawMesh, or the arena flush: all chunks of a body share the rotation)
uniform float texScale;            // world units -> texture repeats

#define MAX_LIGHTS          4
#define LIGHT_DIRECTIONAL   0
#define LIGHT_POINT         1
#define LIGHT_SPOT          2   // position -> target is the beam; its length is the range

struct Light {
    int enabled;
    int type;
    vec3 position;
    vec3 target;
    vec4 color;
};

uniform Light lights[MAX_LIGHTS];
uniform vec4 ambient;
uniform vec3 viewPos;
uniform int palette;            // 0 = space materials, 1 = fantasy, 2 = planet (texture layers)
uniform vec4 sunOccluder;       // planet: render-space centre + radius of a sphere that blocks the sun (w = 0: none)
uniform vec3 fogColor;          // linear
uniform float fogDensity;       // 0 = off
uniform vec4 fogShell;          // atmosphere: render-space centre + outer radius (w = 0: fog the whole distance)
uniform float treeDistance;     // real trees are drawn to here (m); the canopy tint takes over near it

out vec4 finalColor;

// Triplanar projection: one UV set per axis, blended by the surface normal.
struct Triplanar { vec2 uvX; vec2 uvY; vec2 uvZ; vec3 w; };

Triplanar makeTriplanar(vec3 p, vec3 n)
{
    Triplanar t;
    t.uvX = p.zy * texScale;
    t.uvY = p.xz * texScale + vec2(0.33, 0.67);  // offsets decorrelate the three projections
    t.uvZ = p.xy * texScale + vec2(0.71, 0.13);
    vec3 w = pow(abs(n), vec3(4.0));
    t.w = w / (w.x + w.y + w.z);
    return t;
}

vec4 sampleTriplanar(sampler2D tex, Triplanar t)
{
    return texture(tex, t.uvX) * t.w.x + texture(tex, t.uvY) * t.w.y + texture(tex, t.uvZ) * t.w.z;
}

// One layer of the terrain texture array, triplanar, sRGB -> linear.
vec3 sampleLayer(Triplanar t, float layer)
{
    // Triplanar weights are sharpened (pow 4), so most pixels take nearly everything from one axis:
    // skip axes under 2% and renormalise.
    vec3 c = vec3(0.0);
    float sum = 0.0;
    if (t.w.x > 0.02) { c += texture(terrainLayers, vec3(t.uvX, layer)).rgb * t.w.x; sum += t.w.x; }
    if (t.w.y > 0.02) { c += texture(terrainLayers, vec3(t.uvY, layer)).rgb * t.w.y; sum += t.w.y; }
    if (t.w.z > 0.02) { c += texture(terrainLayers, vec3(t.uvZ, layer)).rgb * t.w.z; sum += t.w.z; }
    return pow(c / max(sum, 1e-4), vec3(2.2));
}

// Terrain detail maps (terrainlayers.h, Plan 5 D): RG tangent normal XY, B roughness, A ambient occlusion.
// Adds one layer's slopes, k-weighted, to the per-axis sums and returns its (roughness, ao).
vec2 sampleDetail(Triplanar t, float layer, float k, inout vec2 dX, inout vec2 dY, inout vec2 dZ)
{
    vec2 ra = vec2(0.0);
    float sum = 0.0;
    if (t.w.x > 0.02) { vec4 d = texture(terrainDetail, vec3(t.uvX, layer)); dX += (d.xy * 2.0 - 1.0) * k; ra += d.zw * t.w.x; sum += t.w.x; }
    if (t.w.y > 0.02) { vec4 d = texture(terrainDetail, vec3(t.uvY, layer)); dY += (d.xy * 2.0 - 1.0) * k; ra += d.zw * t.w.y; sum += t.w.y; }
    if (t.w.z > 0.02) { vec4 d = texture(terrainDetail, vec3(t.uvZ, layer)); dZ += (d.xy * 2.0 - 1.0) * k; ra += d.zw * t.w.z; sum += t.w.z; }
    return ra / max(sum, 1e-4);
}

// The same, from the dominant axis only (one fetch): the far-scale blend, where a second full triplanar
// lookup cost more than it showed (#183).
vec3 sampleLayerDominant(Triplanar t, float layer)
{
    vec2 uv = (t.w.x >= t.w.y && t.w.x >= t.w.z) ? t.uvX : (t.w.y >= t.w.z ? t.uvY : t.uvZ);
    return pow(texture(terrainLayers, vec3(uv, layer)).rgb, vec3(2.2));
}

// Triplanar normal mapping, "whiteout" blend (Ben Golus). Uses the sign of the geometric normal
// so faces pointing down an axis get correctly oriented normals, and perturbs the geometric
// normal instead of replacing it.
vec3 triplanarNormal(sampler2D tex, Triplanar t, vec3 n)
{
    vec3 tnX = texture(tex, t.uvX).xyz * 2.0 - 1.0;
    vec3 tnY = texture(tex, t.uvY).xyz * 2.0 - 1.0;
    vec3 tnZ = texture(tex, t.uvZ).xyz * 2.0 - 1.0;
    tnX = vec3(tnX.xy + n.zy, abs(tnX.z) * n.x);
    tnY = vec3(tnY.xy + n.xz, abs(tnY.z) * n.y);
    tnZ = vec3(tnZ.xy + n.xy, abs(tnZ.z) * n.z);
    return normalize(tnX.zyx * t.w.x + tnY.xzy * t.w.y + tnZ.xyz * t.w.z);
}

// Sun shadows: two cascades, hardware PCF (shadow.glsl). The far map is MATERIAL_MAP_BRDF, the near one
// MATERIAL_MAP_HEIGHT.
#include "shadow.glsl"
float shadowFactor(vec3 worldPos, vec3 n, vec3 l) { return sunShadowNormal(worldPos, n, l, 4); }

// Biome tints (Plan 5 B) from the climate packed in fragCanopy.y by the mesher: 1 + temperature step
// (0-62 over -30..40 C) + moisture as the fraction; 0 = no climate (asteroids, caves, old meshes).
bool decodeClimate(out float temp, out float moist)
{
    float y = fragCanopy.y;
    temp = (floor(y) - 1.0) / 62.0 * 70.0 - 30.0;
    moist = fract(y);
    return y >= 1.0;
}
// Multiplier for one material layer's colour: grass dries to straw on the savanna and in deserts,
// deepens in the jungle, dulls in the cold and turns olive in swamps; dirt darkens to mud where wet
// and reddens in the badlands, as does sand.
vec3 biomeTint(int material, float temp, float moist)
{
    float red = smoothstep(16.0, 21.0, temp) * smoothstep(0.2, 0.27, moist) * (1.0 - smoothstep(0.36, 0.44, moist));
    if (material == 6) {   // grass
        float dry = smoothstep(0.66, 0.5, moist) * smoothstep(16.0, 20.0, temp);
        float lush = smoothstep(18.0, 23.0, temp) * smoothstep(0.56, 0.7, moist);
        float cold = smoothstep(9.0, 3.0, temp);
        float swamp = smoothstep(0.68, 0.8, moist);
        vec3 t = vec3(1.0);
        t = mix(t, vec3(2.0, 1.3, 0.45), dry);
        t = mix(t, vec3(0.72, 1.02, 0.68), lush);
        t = mix(t, vec3(0.92, 0.92, 0.95), cold);
        t = mix(t, vec3(0.82, 0.86, 0.55), swamp);
        return t;
    }
    if (material == 1) return mix(mix(vec3(1.0), vec3(1.3, 0.82, 0.62), red), vec3(0.55, 0.5, 0.45), smoothstep(0.72, 0.88, moist)); // dirt -> red / mud
    if (material == 4) return mix(vec3(1.0), vec3(1.22, 0.72, 0.5), red);                                                                // sand -> red sand
    if (material == 0) return mix(vec3(1.0), vec3(1.25, 0.85, 0.7), red * 0.8);                                                         // banded red rock
    return vec3(1.0);
}

// Smooth value noise in [0, 1] (macro variation of the terrain colour).
float hash3(vec3 p) { p = mod(p, 289.0); return fract(sin(dot(p, vec3(127.1, 311.7, 74.7))) * 43758.5453); } // mod: sin() of huge arguments bands
float vnoise(vec3 p)
{
    vec3 i = floor(p), f = fract(p);
    f = f * f * (3.0 - 2.0 * f);
    float x00 = mix(hash3(i), hash3(i + vec3(1, 0, 0)), f.x), x10 = mix(hash3(i + vec3(0, 1, 0)), hash3(i + vec3(1, 1, 0)), f.x);
    float x01 = mix(hash3(i + vec3(0, 0, 1)), hash3(i + vec3(1, 0, 1)), f.x), x11 = mix(hash3(i + vec3(0, 1, 1)), hash3(i + vec3(1, 1, 1)), f.x);
    return mix(mix(x00, x10, f.y), mix(x01, x11, f.y), f.z);
}

void main()
{
    if (any(greaterThan(fragChunkPos, fragClipMax))) discard;
    vec3 geomNormal = normalize(fragNormal);
    // Texture in body-local space so the pattern is fixed to the rock, then bring the bumped
    // normal into world space for lighting.
    vec3 localNormal = normalize(fragLocalNormal);
    Triplanar tp = makeTriplanar(fragTexPos, localNormal);
    float farMix = smoothstep(60.0, 280.0, length(viewPos - fragPosition));
    Triplanar tpFar = makeTriplanar(fragTexPos * 0.21 + vec3(17.3, 5.1, 9.7), localNormal);

    // The base rock / roughness maps feed the asteroid and fantasy palettes; the planet palette takes
    // its colours from the texture layers instead, so it skips these six fetches.
    vec3 albedo = vec3(0.5);
    vec2 roughAo = vec2(0.8, 1.0);
    // Planets bump from their own layers below (the crater normal map is not used: three fetches saved).
    vec3 normal = geomNormal;
    if (palette != 2) {
        albedo = pow(sampleTriplanar(colorTexture, tp).rgb, vec3(2.2)); // sRGB -> linear
        roughAo = sampleTriplanar(roughAoTexture, tp).rg;
        normal = normalize(mat3(matModel) * triplanarNormal(normalTexture, tp, localNormal));
    }

    float rough = roughAo.r;
    float ao = roughAo.g;

    // Voxel material blend weights from the vertex colour: R rock, G regolith, B ice, A ore.
    // (Meshes without vertex colours get white = equal weights, which normalises to a neutral mix.)
    vec4 w = fragColor;
    w /= max(w.r + w.g + w.b + w.a, 0.001);
    float lum = dot(albedo, vec3(0.299, 0.587, 0.114));
    vec3 rock, regolith, ice, ore;
    vec3 emissive = vec3(0.0);
    if (palette == 1) {
        // Fantasy: stone, grass/soil, snow, glowing crystal.
        rock     = albedo * vec3(0.95, 0.95, 1.0);
        regolith = mix(albedo * vec3(0.45, 0.85, 0.3), albedo * vec3(0.9, 0.7, 0.45), 0.25); // grass, a touch of soil
        ice      = mix(vec3(0.85, 0.88, 0.95), vec3(1.0), lum * 1.5);                   // snow
        ore      = vec3(0.35, 0.45, 1.0) * (0.6 + lum);                                 // crystal
        emissive = vec3(0.15, 0.25, 0.9) * w.a * 0.8;                                   // crystals glow
    } else {
        rock     = albedo;
        regolith = albedo * vec3(1.45, 1.25, 1.0);                      // dusty, warmer, lighter
        ice      = mix(vec3(0.55, 0.72, 0.9), vec3(0.85, 0.95, 1.0), lum * 2.0);
        ore      = albedo * 0.5 + vec3(0.55, 0.3, 0.08) * (0.4 + lum);  // rusty metal
    }
    float shinyMaterials = w.b + w.a; // ice + ore: extra specular
    if (palette == 2) {
        // Planet: eight materials from the texture layers (weights 0-3 in the vertex colour, 4-5 in
        // texcoord, 6-7 in texcoord2). Only layers with weight are sampled (usually one or two).
        float w8[8] = float[8](fragColor.r, fragColor.g, fragColor.b, fragColor.a,
                               fragTexCoord.x, fragTexCoord.y, fragTexCoord2.x, fragTexCoord2.y);
        float total = 0.0;
        float ctemp, cmoist;
        bool hasClimate = decodeClimate(ctemp, cmoist);
        for (int i = 0; i < 8; i++) total += w8[i];
        // material -> texture layer: stone, dirt, ice, ore (stone, rusted), sand, gravel, grass, snow
        const float layerOf[8] = float[8](6.0, 0.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0);
        const float roughOf[8] = float[8](0.75, 0.95, 0.15, 0.35, 0.9, 0.85, 0.95, 0.7);
        // Strength of the (asteroid rock) normal map per material: rocky ones only; on grass, dirt and
        // sand its craters read as the asteroid texture showing through.
        const float bumpOf[8] = float[8](1.0, 0.8, 0.6, 1.0, 0.7, 1.0, 0.7, 0.5); // detail normal strength
        vec3 a = vec3(0.0);
        float r = 0.0, aoSum = 0.0;
        // Detail maps fade out by ~120 m (mips flatten them before that; #183: sampling them out to 280 m cost 0.2 ms)
        float near = 1.0 - smoothstep(40.0, 120.0, length(viewPos - fragPosition));
        vec2 dX = vec2(0.0), dY = vec2(0.0), dZ = vec2(0.0);
        for (int i = 0; i < 8; i++) {
            float wi = w8[i] / max(total, 0.001);
            if (wi < 0.01) continue;
            vec3 c = sampleLayer(tp, layerOf[i]);
            // Far away the repeat of one scale shows: blend in the same layer ~5x larger.
            if (farMix > 0.01) c = mix(c, sampleLayerDominant(tpFar, layerOf[i]), 0.55 * farMix); // one fetch: the dominant axis (#183)
            if (hasClimate) c *= biomeTint(i, ctemp, cmoist);
            if (i == 3) {
                // Ore: flecks of the mineral's colour in the stone, by kind (iron rust, copper green,
                // tin pale, coal black, bauxite red, ilmenite dark violet).
                int kind = int(fragOreKind * 8.0 + 0.5) - 1;
                vec3 mineral = kind == 1 ? vec3(0.10, 0.42, 0.30) : kind == 2 ? vec3(0.62, 0.62, 0.64) : kind == 3 ? vec3(0.03, 0.03, 0.035)
                             : kind == 4 ? vec3(0.62, 0.22, 0.16) : kind == 5 ? vec3(0.13, 0.08, 0.16) : vec3(0.58, 0.27, 0.08);
                // Veins: thin warped bands of the mineral through the rock, plus crystal glints between them.
                float fleck = smoothstep(0.30, 0.55, dot(c, vec3(0.333)) + 0.15 * sin(dot(fragTexPos, vec3(7.1, 5.3, 6.7))));
                float vp = dot(fragTexPos, vec3(0.37, 0.81, 0.45)) * 1.4 + vnoise(fragTexPos * 0.7) * 1.8;
                float vein = smoothstep(0.16, 0.04, abs(fract(vp) - 0.5));
                float glint = step(0.94, hash3(floor(fragTexPos * 7.0)));
                c = mix(c * 0.6, mineral, 0.2 + 0.35 * fleck + 0.45 * vein);
                c += mineral * glint * 0.6;
            }
            a += c * wi;
            vec2 ra = vec2(roughOf[i], 1.0);
            if (near > 0.01) ra = mix(ra, sampleDetail(tp, layerOf[i], wi * bumpOf[i] * near, dX, dY, dZ), near);
            r += ra.x * wi;
            aoSum += ra.y * wi;
        }
        albedo = a;
        // Slow brightness and warmth changes over ~70-250 m break up what tiling is left.
        float macro = fragMacro; // per vertex (70-250 m features; #183)
        albedo *= mix(vec3(0.84, 0.86, 0.9), vec3(1.14, 1.1, 1.02), macro);
        if (fragCanopy.x > 0.01) {
            // Distant forests: where the real trees (drawn to treeDistance) fade out, the ground takes the
            // colour of their crowns, by latitude as the trees do, broken up a little by the texture.
            float c = fragCanopy.x * smoothstep(0.65 * treeDistance, 0.95 * treeDistance, length(viewPos - fragPosition));
            float ct, cm;
            bool cl = decodeClimate(ct, cm);
            vec3 dry = vec3(0.22, 0.27, 0.07), temperate = vec3(0.09, 0.27, 0.05), cold = vec3(0.05, 0.17, 0.08);
            // By the climate as the trees are coloured: dry olive where warm and dry, lush in the temperate
            // and wet warm lands, dark blue-green in the cold.
            vec3 crowns = cl ? mix(mix(temperate, dry, smoothstep(0.5, 0.25, cm) * smoothstep(12.0, 22.0, ct)), cold, smoothstep(8.0, 0.0, ct)) : temperate;
            vec2 crown = floor(fragTexPos.xz / 5.0); // ~5 m dapple, like the prism crowns
            crowns *= 0.75 + 0.5 * fract(sin(dot(crown, vec2(12.9898, 78.233))) * 43758.5453);
            albedo = mix(albedo, crowns * 0.5, c); // crowns in their own shade read darker than open grass
            r = mix(r, 0.95, c);
        }
        rough = r;
        ao = aoSum;
        albedo *= mix(1.0, ao, 0.35);            // cavities read in direct light too (stylized)
        // Normals from the layers' detail maps, triplanar "whiteout" blend (as triplanarNormal) of the
        // weighted per-axis slopes, in body-local space, then to world.
        if (near > 0.01) {
            vec3 n = localNormal;
            vec3 tX = vec3(dX + n.zy, abs(n.x)), tY = vec3(dY + n.xz, abs(n.y)), tZ = vec3(dZ + n.xy, abs(n.z));
            tX.z *= sign(n.x); tY.z *= sign(n.y); tZ.z *= sign(n.z);
            vec3 bumped = normalize(tX.zyx * tp.w.x + tY.xzy * tp.w.y + tZ.xyz * tp.w.z);
            normal = normalize(mat3(matModel) * bumped);
        } else {
            normal = normalize(mat3(matModel) * localNormal);
        }
        shinyMaterials = (w8[2] + w8[3] + 0.3 * w8[7]) / max(total, 0.001);
    } else {
        albedo = rock * w.r + regolith * w.g + ice * w.b + ore * w.a;
        rough = rough * (w.r + w.g) + 0.15 * w.b + 0.3 * w.a;
        rough = max(rough, w.g);                                             // regolith is matte
    }

    float shininess = mix(64.0, 4.0, rough);
    float specStrength = (1.0 - rough) * (0.5 + 0.8 * shinyMaterials);

    vec3 viewD = normalize(viewPos - fragPosition);
    vec3 direct = vec3(0.0);
    for (int i = 0; i < MAX_LIGHTS; i++)
    {
        if (lights[i].enabled == 0) continue;
        vec3 l = (lights[i].type == LIGHT_DIRECTIONAL)
            ? normalize(lights[i].position - lights[i].target)
            : normalize(lights[i].position - fragPosition);
        float NdotL = max(dot(normal, l), 0.0);
        vec3 h = normalize(l + viewD);
        float spec = (NdotL > 0.0) ? pow(max(dot(normal, h), 0.0), shininess) * specStrength : 0.0;
        float shadow = 1.0;
        if (i == 0) {
            // The sun: the local shadow map, how buried this point is (caves get no sun even when
            // the occluding rock lies beyond the shadow map), and the planet itself (the night side,
            // at any distance): no sun if the ray toward it enters the occluder sphere, or starts inside.
            shadow = shadowFactor(fragPosition, geomNormal, l) * fragSky;
            if (sunOccluder.w > 0.0) {
                vec3 oc = fragPosition - sunOccluder.xyz;
                float b = dot(oc, l), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
                if (b < 0.0) shadow *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0))); // soft-edged night side
            }
        }
        // Point lights fall off with distance; their range is passed in target.x (target is unused
        // for point lights). Colour alpha scales intensity.
        float atten = 1.0;
        if (lights[i].type == LIGHT_POINT) {
            float dl = length(lights[i].position - fragPosition) / max(lights[i].target.x, 1.0);
            atten = lights[i].color.a * 3.0 / (1.0 + dl * dl * 4.0);
        } else if (lights[i].type == LIGHT_SPOT) {
            vec3 beam = lights[i].target - lights[i].position;
            float range = max(length(beam), 1.0);
            float cone = smoothstep(0.82, 0.95, dot(-l, beam / range)); // ~35 deg outer, ~18 deg inner
            float dl = length(lights[i].position - fragPosition) / range;
            atten = lights[i].color.a * 3.0 * cone / (1.0 + dl * dl * 4.0);
        }
        direct += (albedo * NdotL + spec) * lights[i].color.rgb * shadow * atten;
    }

    vec3 color = ambient.rgb * albedo * ao * fragSky + direct + emissive;
    // Fog over the part of the view ray inside the air (fogShell: render-space centre + outer radius;
    // w = 0: the whole distance), so a moon or ring far beyond the atmosphere stays visible.
    float fogDist = length(viewPos - fragPosition);
    if (fogShell.w > 0.0) {
        vec3 rd = (fragPosition - viewPos) / max(fogDist, 1e-3);
        vec3 oc = viewPos - fogShell.xyz;
        float b = dot(oc, rd), c = dot(oc, oc) - fogShell.w * fogShell.w, disc = b * b - c;
        float t0 = disc > 0.0 ? -b - sqrt(disc) : 0.0, t1 = disc > 0.0 ? -b + sqrt(disc) : 0.0;
        fogDist = max(0.0, min(t1, fogDist) - max(t0, 0.0));
    }
    float fog = 1.0 - exp(-fogDensity * fogDist);
    color = mix(color, fogColor, clamp(fog, 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0); // linear -> sRGB
}
