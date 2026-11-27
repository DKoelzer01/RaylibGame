#version 330

in vec3 fragPosition;   // world space
in vec2 fragTexCoord;    // chunks: material weights 4-5 (sand, gravel)
in vec2 fragTexCoord2;   // chunks: material weights 6-7 (grass, snow)
in vec4 fragColor;
in vec3 fragNormal;     // world space
in vec3 fragTexPos;     // body-aligned, wrapped to the texture repeat (texturing)
in vec3 fragLocalNormal;
in vec3 fragChunkPos;   // chunk-local (overlap clipping)
centroid in float fragSky;       // 1 = open to the sky, 0 = buried (caves, deep pits); centroid: #402
centroid in float fragOcclusion; // short-range ambient occlusion (#132): 0 open .. 1 enclosed
in vec2 fragCanopy;     // planets: forest canopy cover and latitude (painted in beyond the real trees)
flat in float fragClimate; // the packed climate, per triangle (#254: interpolated, it striped the biome tints)
in float fragOreKind;   // planets: which ore, (1 + kind) / 8 (0 none)
in float fragMacro;     // planets: slow brightness / warmth variation (70-250 m), per vertex (lighting.vs)

// Material maps: bound by DrawMesh to the units of the matching raylib material slots.
uniform sampler2D colorTexture;    // MATERIAL_MAP_ALBEDO
uniform sampler2D normalTexture;   // MATERIAL_MAP_NORMAL (OpenGL convention, tangent space)
uniform sampler2D roughAoTexture;  // MATERIAL_MAP_ROUGHNESS: R = roughness, G = ambient occlusion
uniform sampler2DArray terrainLayers; // planet palette: dirt, sand, gravel, grass, snow, ice, stone, clay, granite, basalt, sandstone, limestone, marble (terrainlayers.h)
uniform sampler2DArray terrainDetail; // its detail maps: normal XY, roughness, AO (unit 12)

flat in vec3 fragClipMax;
flat in float fragLodFade;   // #243: 0 none; (0, 1) fading in; (-1, 0) fading out (the share already given to its replacement)     // chunk-local overlap clip (lighting.vs: a uniform or the arena record)
flat in uint fragMatIds;       // #480 (arena): the triangle's three material ids, 8 bits each (fragColor.xyz: their barycentrics)
uniform mat4 matModel;          // body-local -> world (DrawMesh, or the arena flush: all chunks of a body share the rotation)
uniform float texScale;            // world units -> texture repeats
uniform int arenaDraw;             // 1: a terrain arena multi-draw (lighting.vs; #480: the triangle's palette slots)

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
uniform int debugSky;        // RAYLIBGAME_SHOW_SKY=1: draw fragSky as grey (red = the #898 cave-mouth floor range, 0.55-0.65)
uniform int bounceOn;         // (#898) the sunlit ground's bounce onto walls (RAYLIBGAME_BOUNCE=0: off)
uniform vec4 toolSphere;     // (#1089) the dig / fill sphere: render-space centre, radius (w = 0: none)
uniform vec4 toolSphereColor;
uniform vec3 viewPos;
uniform int palette;            // 0 = space materials, 1 = fantasy, 2 = planet (texture layers)
uniform vec4 sunOccluder;       // planet: render-space centre + radius of a sphere that blocks the sun (w = 0: none)
uniform vec3 fogColor;          // linear
uniform float fogDensity;       // 0 = off
uniform vec4 fogShell;          // atmosphere: render-space centre + outer radius (w = 0: fog the whole distance)
uniform float treeDistance;     // real trees are drawn to here (m); the canopy tint takes over near it
uniform float forestFloorAo;   // (#933 layer 2) strength: the forest floor under the real trees takes less sky light (0: off; RAYLIBGAME_FOLIAGE_AO=0, _STRENGTH)
uniform float farForest;        // #745: 1 = deeper canopy tint, starting earlier and denser (RAYLIBGAME_FAR_FOREST)

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
#include "skyphys_haze.glsl" // #142: the physical haze colour (off by default)
#include "hemisphere.glsl"
#include "biometint.glsl" // #2032
// #138: past the far cascade (+-200 m) the wide terrain-only cascade (+-2 km, ~2 m texels, MATERIAL_MAP_EMISSION), so
// mountains shadow valleys at dawn and dusk. It fills in exactly where the far cascade fades to lit at its edge.
uniform sampler2DShadow shadowMapWide;
uniform mat4 lightSpaceWide;
uniform vec2 shadowWide;     // texel (m), depth range (m)
uniform int wideShadowOn;
float shadowFactor(vec3 worldPos, vec3 n, vec3 l)
{
    float s = sunShadowNormal(worldPos, n, l, 4);
    if (wideShadowOn == 0 || shadowMapResolution <= 0) return s;
    // The far cascade's edge weight here (as shadowCascade works it out): 1 inside, 0 beyond.
    vec4 ls = lightSpaceMatrix * vec4(worldPos + n * shadowTexel.y, 1.0);
    vec3 q = ls.xyz / ls.w * 0.5 + 0.5;
    vec2 d = min(q.xy, 1.0 - q.xy);
    float eFar = (q.z > 1.0) ? 0.0 : clamp((min(d.x, d.y) - 0.01) / 0.07, 0.0, 1.0);
    if (eFar >= 1.0) return s;
    float ndl = clamp(dot(n, l), 0.0, 1.0);
    float slope = min(sqrt(max(1.0 - ndl * ndl, 0.0)) / max(ndl, 0.2), 4.0);
    vec3 pw = worldPos + n * shadowWide.x * (1.5 - 0.5 * ndl);
    vec4 lw = lightSpaceWide * vec4(pw, 1.0);
    vec3 qw = lw.xyz / lw.w * 0.5 + 0.5;
    vec2 dw = min(qw.xy, 1.0 - qw.xy);
    float eW = (qw.z > 1.0) ? 0.0 : clamp((min(dw.x, dw.y) - 0.01) / 0.07, 0.0, 1.0);
    if (eW <= 0.0) return s;
    qw.z -= shadowWide.x * (0.5 + slope) / shadowWide.y;
    float sw = mix(1.0, shadowPcf(shadowMapWide, qw, 4), eW);
    // Where the far cascade fades out (s there = mix(1, sFar, eFar)), the wide one takes over: mix(sw, sFar, eFar).
    return s - (1.0 - eFar) * (1.0 - sw);
}

// Biome tints (Plan 5 B) from the climate packed by the mesher: 2 + temperature step (0-62 over
// -30..40 C) + moisture as the fraction; below 1.5 = no climate (asteroids, caves, old meshes). Read flat:
// interpolated across a triangle the packing swept through whole steps and striped the tints (#254).
bool decodeClimate(out float temp, out float moist)
{
    float y = fragClimate;
    temp = biomeTempOf(y); // (#2032: the variant id rides on the step)
    moist = fract(y);
    return y >= 1.5;
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
        // Temperate land has its own gentle gradient (#254: it all read as the default green): cooler a little
        // blue-green, warmer a little yellow; wetter deeper, drier paler.
        vec3 t = mix(vec3(0.9, 1.0, 1.02), vec3(1.1, 1.05, 0.86), smoothstep(8.0, 17.0, temp));
        t *= mix(vec3(1.08, 1.04, 0.9), vec3(0.88, 1.0, 0.9), smoothstep(0.3, 0.68, moist));
        t = mix(t, vec3(2.0, 1.3, 0.45), dry);
        t = mix(t, vec3(0.72, 1.02, 0.68), lush);
        t = mix(t, vec3(0.92, 0.92, 0.95), cold);
        t = mix(t, vec3(0.82, 0.86, 0.55), swamp);
        return t;
    }
    if (material == 1) return mix(mix(vec3(1.0), vec3(1.3, 0.82, 0.62), red), vec3(0.55, 0.5, 0.45), smoothstep(0.72, 0.88, moist)); // dirt -> red / mud
    if (material == 4) return mix(vec3(1.0), vec3(1.22, 0.72, 0.5), red);                                                                // sand -> red sand
    if (material == 8) return mix(vec3(1.0), vec3(0.72, 0.84, 0.92), 0.8 * smoothstep(0.72, 0.9, moist));                              // clay -> grey-blue (gley) in swamps
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

// Value noise in [0, 1] at `scale` features per metre of fragTexPos, rounded so it repeats with the terrain texture:
// fragTexPos jumps by whole texture periods between chunks, so this noise shows no chunk seams (vnoise above would).
float periodicNoise(vec3 p, float scale)
{
    float per = max(1.0, floor(scale / max(texScale, 1e-4) + 0.5));
    vec3 q = p * (per * max(texScale, 1e-4)), i = floor(q), f = fract(q);
    f = f * f * (3.0 - 2.0 * f);
    float x00 = mix(hash3(mod(i, per)), hash3(mod(i + vec3(1, 0, 0), per)), f.x), x10 = mix(hash3(mod(i + vec3(0, 1, 0), per)), hash3(mod(i + vec3(1, 1, 0), per)), f.x);
    float x01 = mix(hash3(mod(i + vec3(0, 0, 1), per)), hash3(mod(i + vec3(1, 0, 1), per)), f.x), x11 = mix(hash3(mod(i + vec3(0, 1, 1), per)), hash3(mod(i + vec3(1, 1, 1), per)), f.x);
    return mix(mix(x00, x10, f.y), mix(x01, x11, f.y), f.z);
}

void main()
{
    if (any(greaterThan(fragChunkPos, fragClipMax))) discard;
    // LOD cross-fade (#243): a node and the one replacing it are both drawn for a moment and split the pixels by the same
    // screen noise (interleaved gradient noise), so each pixel shows exactly one of them and the swap dissolves.
    if (fragLodFade != 0.0) {
        float n = fract(52.9829189 * fract(dot(gl_FragCoord.xy, vec2(0.06711056, 0.00583715))));
        float t = abs(fragLodFade);
        if (fragLodFade > 0.0 ? n >= t : n < t) discard;
    }
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
        // Planet: materials from the texture layers, by id. Terrain arena (#480): the triangle's three materials (ids
        // from the chunk's palette), weighed by the barycentrics in fragColor (corners of one material merged).
        // Per-chunk draws: the 8 legacy weights (0-3 in the vertex colour, 4-5 in texcoord, 6-7 in texcoord2), ids 0-7.
        // Only materials with weight are sampled (usually one or two).
        float ctemp, cmoist;
        bool hasClimate = decodeClimate(ctemp, cmoist);
        // material id -> texture layer: stone, dirt, ice, ore (stone, rusted), sand, gravel, grass, snow, clay (#206),
        // then the rock strata (#587): granite, basalt, sandstone, limestone, marble, then pebbles, silt, mud (#595), then
        // farmland (#1085: tilled dirt, the dirt layer with furrows below)
        const int kMaterials = 18;
        const float layerOf[kMaterials] = float[kMaterials](6.0, 0.0, 5.0, 6.0, 1.0, 2.0, 3.0, 4.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 0.0);
        const float roughOf[kMaterials] = float[kMaterials](0.75, 0.95, 0.15, 0.35, 0.9, 0.85, 0.95, 0.7, 0.6, 0.65, 0.75, 0.85, 0.75, 0.45, 0.7, 0.9, 0.45, 0.95);
        // Strength of the (asteroid rock) normal map per material: rocky ones only; on grass, dirt and
        // sand its craters read as the asteroid texture showing through.
        const float bumpOf[kMaterials] = float[kMaterials](1.0, 0.8, 0.6, 1.0, 0.7, 1.0, 0.7, 0.5, 0.6, 1.0, 1.0, 0.8, 0.8, 0.5, 1.0, 0.5, 0.4, 0.5); // detail normal strength
        const float shinyOf[kMaterials] = float[kMaterials](0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.3, 0.25, 0.15, 0.1, 0.0, 0.05, 0.4, 0.2, 0.0, 0.45, 0.0); // ice, ore, a little snow, wet clay, polished-ish marble, wet stones and mud
        vec3 a = vec3(0.0);
        float r = 0.0, aoSum = 0.0, shiny = 0.0;
        // Detail maps fade out by ~120 m (mips flatten them before that; #183: sampling them out to 280 m cost 0.2 ms)
        float near = 1.0 - smoothstep(40.0, 120.0, length(viewPos - fragPosition));
        vec2 dX = vec2(0.0), dY = vec2(0.0), dZ = vec2(0.0);
        vec4 aw = vec4(0.0);
        ivec4 aid = ivec4(0);
        int n = 0;
        float total = 1.0;
        if (arenaDraw == 1) {
            int m0 = int(fragMatIds & 255u), m1 = int((fragMatIds >> 8) & 255u), m2 = int((fragMatIds >> 16) & 255u);
            vec3 b = fragColor.xyz;
            if (m1 == m0) { b.x += b.y; b.y = 0.0; }
            if (m2 == m0) { b.x += b.z; b.z = 0.0; } else if (m2 == m1) { b.y += b.z; b.z = 0.0; }
            aw = vec4(b, 0.0);
            aid = ivec4(m0, m1, m2, 0);
            n = 3;
        } else {
            vec4 w0 = fragColor, w1 = vec4(fragTexCoord, fragTexCoord2);
            total = dot(w0, vec4(1.0)) + dot(w1, vec4(1.0));
            float cut = 0.01 * max(total, 0.001);
            // Constant indices only (a loop variable indexing the inputs spilled them to local memory).
#define TAKE(w, c, id) { if (w >= cut && n < 4) { aw[n] = w; aid[n] = id; n++; } }
            TAKE(w0.x, 0, 0) TAKE(w0.y, 1, 1) TAKE(w0.z, 2, 2) TAKE(w0.w, 3, 3) TAKE(w1.x, 4, 4) TAKE(w1.y, 5, 5) TAKE(w1.z, 6, 6) TAKE(w1.w, 7, 7)
#undef TAKE
        }
        // Height blend (#609): where materials meet, each one's texture "height" (its brightness: stones over the sand between
        // them, turf over soil) plus a fine noise decides which shows, so an edge is a ragged, soft line through the triangle
        // instead of the straight cross-fade of the barycentrics (which drew hard, road-like bands along the river banks).
        vec3 cs[4];
        float hw[4];
        float hmax = -1.0;
        // Two octaves of value noise (~1.8 m and ~0.6 m), phased per material id (consistent across triangles), strong
        // enough to move an edge across a whole triangle: the edge then wanders instead of following the mesh's stair steps.
        float blendNoise = 0.6 * periodicNoise(fragTexPos, 0.55) + 0.4 * periodicNoise(fragTexPos + vec3(7.3), 1.6);
        for (int k = 0; k < n; k++) {
            cs[k] = vec3(0.0); hw[k] = -1.0;
            float wi = aw[k] / max(total, 0.001);
            if (wi < 0.01) continue;
            int i = min(aid[k], kMaterials - 1);
            vec3 c = sampleLayer(tp, layerOf[i]);
            // Far away the repeat of one scale shows: blend in the same layer ~5x larger.
            if (farMix > 0.01) c = mix(c, sampleLayerDominant(tpFar, layerOf[i]), 0.55 * farMix); // one fetch: the dominant axis (#183)
            cs[k] = c;
            hw[k] = wi + 0.25 * (dot(sqrt(c), vec3(0.333)) - 0.5) + 0.5 * sin(6.2832 * (2.0 * blendNoise + float(i) * 0.37));
            hmax = max(hmax, hw[k]);
        }
        float hsum = 0.0;
        for (int k = 0; k < n; k++) { hw[k] = hw[k] < -0.9 ? 0.0 : max(hw[k] - hmax + 0.45, 0.0); hsum += hw[k]; }
        for (int k = 0; k < n; k++) {
            float wi = hw[k] / max(hsum, 1e-4);
            if (wi < 0.01) continue;
            int i = min(aid[k], kMaterials - 1);
            shiny += shinyOf[i] * wi;
            vec3 c = cs[k];
            if (hasClimate) c *= i == 6 ? mix(vec3(1.0), biomeTint(i, ctemp, cmoist), biomeClimateWeight) * biomeGrass(biomeIdOf(fragClimate)) : biomeTint(i, ctemp, cmoist); // #2032: grass by biome
            if (i == 3) {
                // Ore: flecks of the mineral's colour in the stone, by kind (iron rust, copper green,
                // tin pale, coal black, bauxite red, ilmenite dark violet).
                int kind = int(fragOreKind * 8.0 + 0.5) - 1;
                vec3 mineral = kind == 1 ? vec3(0.10, 0.42, 0.30) : kind == 2 ? vec3(0.62, 0.62, 0.64) : kind == 3 ? vec3(0.03, 0.03, 0.035)
                             : kind == 4 ? vec3(0.62, 0.22, 0.16) : kind == 5 ? vec3(0.13, 0.08, 0.16) : vec3(0.58, 0.27, 0.08);
                // Each kind also has its own pattern (#1433): iron and bauxite, coal and ilmenite are close in colour for
                // everyone (and copper and bauxite for red-green colour blindness), so the shape tells them apart:
                // iron thin veins, copper veins with green crust patches, tin a fine speckle, coal broad dull seams,
                // bauxite round nodules, ilmenite dense metallic glints.
                float fleck = smoothstep(0.30, 0.55, dot(c, vec3(0.333)) + 0.15 * sin(dot(fragTexPos, vec3(7.1, 5.3, 6.7))));
                float vp = dot(fragTexPos, vec3(0.37, 0.81, 0.45)) * 1.4 + vnoise(fragTexPos * 0.7) * 1.8;
                float vein = smoothstep(0.16, 0.04, abs(fract(vp) - 0.5));
                float glint = step(0.94, hash3(floor(fragTexPos * 7.0)));
                float amount = 0.2 + 0.35 * fleck; // the mineral's share
                if (kind == 0) amount += 0.45 * vein;                                                   // iron: veins
                else if (kind == 1) amount += 0.45 * vein + 0.5 * smoothstep(0.55, 0.68, vnoise(fragTexPos * 1.3)); // copper: + crust
                else if (kind == 2) { amount = 0.2 + 0.45 * step(0.8, hash3(floor(fragTexPos * 12.0))); glint = 0.0; } // tin: speckle
                else if (kind == 3) { // coal: broad layered seams, dull (no glints)
                    float seam = smoothstep(0.30, 0.18, abs(fract(fragTexPos.y * 0.9 + vnoise(fragTexPos * 0.3) * 0.5) - 0.5));
                    amount = 0.15 + 0.8 * seam; glint = 0.0;
                } else if (kind == 4) { // bauxite: round nodules (pisoliths) a few cm to a hand across
                    vec3 q = fragTexPos * 1.8;
                    vec3 cell = floor(q);
                    vec3 o = vec3(hash3(cell), hash3(cell + 17.0), hash3(cell + 31.0)) * 0.4 - 0.2;
                    float nodule = smoothstep(0.32, 0.24, length(fract(q) - 0.5 - o));
                    amount = 0.18 + 0.75 * nodule;
                } else if (kind == 5) { // ilmenite: dense small metallic crystals
                    glint = 0.5 * step(0.86, hash3(floor(fragTexPos * 9.0)));
                    amount += 0.3 + 0.2 * vein; // the dark mineral as the ground, the crystals on it
                }
                c = mix(c * 0.6, mineral, clamp(amount, 0.0, 1.0));
                c += (kind == 5 ? vec3(0.55, 0.55, 0.6) : mineral * 0.6) * glint; // ilmenite's crystals shine steel-grey
            }
            if (i == 17) {
                // Farmland (#1085): darker turned soil in furrows every 25 cm (the old farmland part's pattern, parts.fs).
                // Tilling is a free brush, so the furrows follow a body axis across the ground: the coordinate of the
                // triplanar projection's widest-facing plane that is not the up-most axis (fragTexPos repeats in whole
                // texture periods, a multiple of 25 cm, so the stripes have no chunk seams).
                vec3 an = abs(localNormal);
                float along = an.y >= an.x && an.y >= an.z ? fragTexPos.x : an.x >= an.z ? fragTexPos.z : fragTexPos.x;
                float f = abs(fract(along / 0.25) - 0.5) * 2.0;
                c *= vec3(0.78, 0.72, 0.68) * (0.65 + 0.45 * f);
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
        // (#1022) Strata in buried rock: layers by height along the planet's up (distance from its centre, the fog shell's),
        // so cave walls show beds cutting across a passage instead of the texture reading as rings around the bore. Uneven
        // bed thickness and tone, thin dark bedding seams, fine laminations; the boundaries wander a little.
        if (fogShell.w > 0.0 && fragSky < 0.95) {
            float rockW = 0.0;
            for (int i = 0; i < n; ++i) { int m = aid[i]; if (m == 0 || m == 3 || (m >= 9 && m <= 13)) rockW += aw[i]; }
            rockW = clamp(rockW / max(total, 1e-4), 0.0, 1.0) * (1.0 - smoothstep(0.55, 0.95, fragSky));
            if (rockW > 0.01) {
                float h = length(fragPosition - fogShell.xyz) + 0.7 * (periodicNoise(fragTexPos, 0.12) - 0.5) + 0.25 * (periodicNoise(fragTexPos, 0.5) - 0.5);
                float s = h / 1.9, bed = floor(s), f = fract(s);
                float t1 = hash3(vec3(bed, 7.0, 3.0)), t2 = hash3(vec3(bed, 1.0, 11.0));
                vec3 tone = mix(vec3(0.6, 0.59, 0.58), vec3(1.22, 1.17, 1.08), t1);
                tone *= mix(vec3(1.0), vec3(1.12, 0.96, 0.82), step(0.68, t2));          // some beds warmer (iron-stained)
                float seam = 1.0 - 0.45 * smoothstep(0.07, 0.0, min(f, 1.0 - f) * (0.7 + 0.6 * t2));
                float lam = 0.94 + 0.12 * hash3(vec3(floor(h / 0.35), 5.0, 9.0));
                albedo *= mix(vec3(1.0), tone * seam * lam, rockW);
            }
        }
        if (fragCanopy.x > 0.01) {
            // Distant forests: where the real trees (drawn to treeDistance) fade out, the ground takes the
            // colour of their crowns, by latitude as the trees do, broken up a little by the texture.
            float cdist = length(viewPos - fragPosition);
            float c = fragCanopy.x * smoothstep(0.65 * treeDistance, 0.95 * treeDistance, cdist);
            // #745: thin cover still reads as forest (cover^0.6), the tint starts earlier and its shade is deeper.
            if (farForest > 0.5) c = pow(fragCanopy.x, 0.5) * smoothstep(0.35 * treeDistance, 0.75 * treeDistance, cdist);
            float ct, cm;
            bool cl = decodeClimate(ct, cm);
            vec3 dry = vec3(0.22, 0.27, 0.07), temperate = vec3(0.09, 0.27, 0.05), cold = vec3(0.05, 0.17, 0.08);
            // By the climate as the trees are coloured: dry olive where warm and dry, lush in the temperate
            // and wet warm lands, dark blue-green in the cold.
            vec3 crowns = cl ? mix(mix(temperate, dry, smoothstep(0.5, 0.25, cm) * smoothstep(12.0, 22.0, ct)), cold, smoothstep(8.0, 0.0, ct)) : temperate;
            vec2 crown = floor(fragTexPos.xz / 5.0); // ~5 m dapple, like the prism crowns
            if (cl) crowns = biomeClimate(crowns, temperate) * biomeLeaf(biomeIdOf(fragClimate)); // #2032: as the trees
            crowns *= 0.75 + 0.5 * fract(sin(dot(crown, vec2(12.9898, 78.233))) * 43758.5453);
            albedo = mix(albedo, crowns * (farForest > 0.5 ? 0.32 : 0.5), c); // crowns in their own shade read darker than open grass
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
        shinyMaterials = shiny;
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
    // (#898) Bounce: the sunlit ground reflects light onto what faces sideways or down - a pit's or shaft's walls, a cliff's
    // shaded face (user: "light that behaves like light"; a hand-dug shaft's walls read near black at noon: they faced
    // no sun and took only the sky's small ambient). The ground's irradiance x ~0.25 albedo x how far the surface faces
    // down (0 up, 0.5 sideways, 1 down) x the sun's height; only where the sky term says open ground, a pit or a shaft
    // (fragSky >= ~0.85 there: #1067's open column), so caves and cave mouths (<= ~0.6) stay dark. bounceOn 0: off (A/B).
    if (bounceOn != 0 && lights[0].enabled != 0 && hemiOn()) {
        vec3 l0 = normalize(lights[0].position - lights[0].target);
        float sunUp = max(dot(ambientUp, l0), 0.0);
        float facingDown = 0.5 - 0.5 * dot(normal, ambientUp);
        float open = smoothstep(0.55, 0.85, fragSky);
        direct += albedo * lights[0].color.rgb * (0.25 * sunUp * facingDown * open);
    }
    direct += albedo * (localLights(fragPosition, normal) + gridLight(fragPosition, normal)); // lamps, campfires, glowing machines (#139); placed lights through the grid (#1739)

    // Ambient occlusion (#132) darkens the ambient term only: dug holes and crevices read under flat sky light.
    float occlusionAO = 1.0 - 0.75 * clamp(fragOcclusion, 0.0, 1.0);
    // (#933 layer 2) Under the real trees (to where the painted canopy takes over), the forest floor's sky is mostly crowns:
    // its ambient drops with the canopy cover (direct sun there comes from the shadow maps already).
    float floorAo = 1.0;
    if (forestFloorAo > 0.0 && fragCanopy.x > 0.01) floorAo = 1.0 - min(0.45 * forestFloorAo, 0.9) * fragCanopy.x * (1.0 - smoothstep(0.65 * treeDistance, 0.95 * treeDistance, length(viewPos - fragPosition)));
    vec3 color = hemiAmbient(ambient.rgb, normal) * albedo * ao * fragSky * occlusionAO * floorAo + direct + emissive;
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
    color = mix(color, SkyPhysFog(fragPosition, viewPos, fogColor), clamp(fog, 0.0, 1.0));
    finalColor = vec4(pow(max(color, vec3(0.0)), vec3(1.0 / 2.2)), 1.0); // linear -> sRGB
    // (#1089) The dig / fill preview: where the sphere cuts the terrain, a line on the real surface (over bumps, ledges and
    // cave walls, depth-tested by the terrain itself) and a faint tint inside. Only the planet terrain draws it.
    if (toolSphere.w > 0.0 && palette == 2) {
        float sd = length(fragPosition - toolSphere.xyz) - toolSphere.w;
        float fw = max(fwidth(sd), 1e-4);
        float line = 1.0 - smoothstep(0.75 * fw, 2.0 * fw, abs(sd));
        float inside = sd < 0.0 ? 0.18 : 0.0;
        finalColor.rgb = mix(finalColor.rgb, toolSphereColor.rgb, max(line * toolSphereColor.a, inside * toolSphereColor.a));
    }
    if (debugSky == 1) finalColor = vec4(fragSky, abs(fragSky - 0.6) < 0.05 ? 0.0 : fragSky, abs(fragSky - 0.6) < 0.05 ? 0.0 : fragSky, 1.0);
}
