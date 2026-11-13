#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;   // chunks: material weights 4-5 (sand, gravel)
in vec2 vertexTexCoord2;  // chunks: material weights 6-7 (grass, snow)
in vec3 vertexNormal;
in vec4 vertexColor;
in vec4 vertexTangent;    // chunks: w = sky exposure (default attribute (0,0,0,1): fully exposed)
// Terrain arena (terrainarena.h): this draw's record, instanced attributes picked by the command's baseInstance.
layout(location = 9) in vec4 arenaRow0;   // model rows 0-2 (render space)
layout(location = 10) in vec4 arenaRow1;
layout(location = 11) in vec4 arenaRow2;
layout(location = 12) in vec4 arenaTexOrigin;
layout(location = 13) in vec4 arenaClip;

// Set by DrawMesh() from the material's shader locations
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
// Set per chunk draw: the chunk's body-local origin wrapped to the texture repeat (so texture
// coordinates stay small and precise even 60 km from a body's centre). Zero for other meshes.
uniform vec3 chunkTexOrigin;
uniform vec3 chunkClipMax;      // chunk-local: the fragment shader discards beyond this (Chunk::clipMax)
// Terrain arena (Plan 6 step 1): one multi-draw for all chunks of a body (per-draw data in the arena attributes above).
uniform int arenaDraw;
uniform mat4 matViewProj;

out vec3 fragPosition;     // world (render) space: lighting, shadows, fog
out vec2 fragTexCoord;
out vec2 fragTexCoord2;
out vec4 fragColor;
out vec3 fragNormal;       // world space
out vec3 fragTexPos;       // body-aligned texturing position (wrapped): textures stay attached to a
out vec3 fragLocalNormal;  // moving/spinning body and don't swim when the floating origin shifts
out vec3 fragChunkPos;     // chunk-local position (overlap clipping)
out float fragSky;         // 1 = open to the sky, 0 = buried (scales sunlight)
out vec2 fragCanopy;       // planets: distant forest canopy cover (x) and latitude (y)
flat out float fragClimate; // planets: the packed climate (2 + temperature step + moisture), per triangle (#254)
out float fragOreKind;     // planets: which ore (1 + kind) / 8, 0 none
out float fragMacro;       // slow (70-250 m) brightness variation for the planet palette (per vertex: #183)
flat out vec3 fragClipMax;

// Smooth value noise in [0, 1] (the planet palette's macro variation; was per pixel in lighting.fs).
float hash3(vec3 p) { p = mod(p, 289.0); return fract(sin(dot(p, vec3(127.1, 311.7, 74.7))) * 43758.5453); }
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
    mat4 model = matModel;
    vec3 texOrigin = chunkTexOrigin;
    fragClipMax = chunkClipMax;
    if (arenaDraw == 1) {
        model = transpose(mat4(arenaRow0, arenaRow1, arenaRow2, vec4(0.0, 0.0, 0.0, 1.0)));
        texOrigin = arenaTexOrigin.xyz;
        fragClipMax = arenaClip.xyz;
    }
    vec4 world = model * vec4(vertexPosition, 1.0);
    fragPosition = world.xyz;
    fragTexCoord = vertexTexCoord;
    fragTexCoord2 = vertexTexCoord2;
    fragColor = vertexColor;
    // w = 0: normals are directions and must not pick up the model translation.
    // Normal to world space via the model matrix. raylib also provides matNormal, but in testing its value
    // did not follow per-draw model rotations (normals stayed in body space on spinning bodies), so we
    // don't rely on it. Valid for rotation + uniform scale, which is all the engine uses.
    fragNormal = normalize(mat3(model) * vertexNormal);
    fragTexPos = vertexPosition + texOrigin;
    fragChunkPos = vertexPosition;
    fragSky = vertexTangent.w;
    fragCanopy = vertexTangent.xy;
    fragClimate = vertexTangent.y;
    fragOreKind = vertexTangent.z;
    fragMacro = vnoise(fragTexPos / 70.0) * 0.6 + vnoise(fragTexPos / 250.0) * 0.4;
    fragLocalNormal = vertexNormal;
    gl_Position = arenaDraw == 1 ? matViewProj * world : mvp * vec4(vertexPosition, 1.0);
}
