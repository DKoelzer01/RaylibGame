#version 330

in vec3 fragPosition;   // world space
in vec2 fragTexCoord;
in vec4 fragColor;
in vec3 fragNormal;     // world space
in vec3 fragLocalPos;   // body-local (texturing)
in vec3 fragLocalNormal;

// Material maps: bound by DrawMesh to the units of the matching raylib material slots.
uniform sampler2D colorTexture;    // MATERIAL_MAP_ALBEDO
uniform sampler2D normalTexture;   // MATERIAL_MAP_NORMAL (OpenGL convention, tangent space)
uniform sampler2D roughAoTexture;  // MATERIAL_MAP_ROUGHNESS: R = roughness, G = ambient occlusion
uniform sampler2D shadowMap;       // MATERIAL_MAP_BRDF: light depth texture

uniform mat4 lightSpaceMatrix;
uniform mat4 matModel;          // body-local -> world (set by DrawMesh); see lighting.vs about matNormal
uniform int shadowMapResolution;
uniform float texScale;            // world units -> texture repeats

#define MAX_LIGHTS          4
#define LIGHT_DIRECTIONAL   0
#define LIGHT_POINT         1

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
uniform int palette;            // 0 = space materials, 1 = fantasy materials
uniform vec3 chunkClipMax;      // body-local: discard beyond this (see Chunk::clipMax)
uniform vec3 fogColor;          // linear
uniform float fogDensity;       // 0 = off

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

// 1.0 = fully lit, 0.0 = fully shadowed. 3x3 PCF.
float shadowFactor(vec3 worldPos, vec3 n, vec3 l)
{
    vec4 ls = lightSpaceMatrix * vec4(worldPos, 1.0);
    vec3 p = ls.xyz / ls.w * 0.5 + 0.5;
    if (p.x < 0.0 || p.x > 1.0 || p.y < 0.0 || p.y > 1.0 || p.z > 1.0) return 1.0; // outside the shadow box
    float bias = max(0.002 * (1.0 - dot(n, l)), 0.0005);
    vec2 texel = vec2(1.0 / float(shadowMapResolution));
    float lit = 0.0;
    for (int x = -1; x <= 1; x++)
        for (int y = -1; y <= 1; y++)
            lit += (p.z - bias > texture(shadowMap, p.xy + texel * vec2(x, y)).r) ? 0.0 : 1.0;
    return lit / 9.0;
}

void main()
{
    if (any(greaterThan(fragLocalPos, chunkClipMax))) discard;
    vec3 geomNormal = normalize(fragNormal);
    // Texture in body-local space so the pattern is fixed to the rock, then bring the bumped
    // normal into world space for lighting.
    vec3 localNormal = normalize(fragLocalNormal);
    Triplanar tp = makeTriplanar(fragLocalPos, localNormal);

    vec3 albedo = pow(sampleTriplanar(colorTexture, tp).rgb, vec3(2.2)); // sRGB -> linear
    vec2 roughAo = sampleTriplanar(roughAoTexture, tp).rg;
    vec3 normal = normalize(mat3(matModel) * triplanarNormal(normalTexture, tp, localNormal));

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
    albedo = rock * w.r + regolith * w.g + ice * w.b + ore * w.a;
    rough = rough * (w.r + w.g) + 0.15 * w.b + 0.3 * w.a;
    rough = max(rough, w.g);                                             // regolith is matte

    float shininess = mix(64.0, 4.0, rough);
    float specStrength = (1.0 - rough) * (0.5 + 0.8 * (w.b + w.a));

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
        float shadow = (i == 0) ? shadowFactor(fragPosition, geomNormal, l) : 1.0;
        // Point lights fall off with distance; their range is passed in target.x (target is unused
        // for point lights). Colour alpha scales intensity.
        float atten = 1.0;
        if (lights[i].type == LIGHT_POINT) {
            float dl = length(lights[i].position - fragPosition) / max(lights[i].target.x, 1.0);
            atten = lights[i].color.a * 3.0 / (1.0 + dl * dl * 4.0);
        }
        direct += (albedo * NdotL + spec) * lights[i].color.rgb * shadow * atten;
    }

    vec3 color = ambient.rgb * albedo * ao + direct + emissive;
    float fog = 1.0 - exp(-fogDensity * length(viewPos - fragPosition));
    color = mix(color, fogColor, clamp(fog, 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0); // linear -> sRGB
}
