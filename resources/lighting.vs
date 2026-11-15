#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec4 vertexTexCoord;   // chunks: material weights 4-5 (sand, gravel); arena: palette slot weights 0-3 (#480)
in vec4 vertexTexCoord2;  // chunks: material weights 6-7 (grass, snow); arena: slot weights 8-11
in vec3 vertexNormal;
in vec4 vertexColor;      // chunks: material weights 0-3; arena: slot weights 4-7
in vec4 vertexTangent;    // chunks: w = sky exposure (default attribute (0,0,0,1): fully exposed)
layout(location = 14) in float vertexOcclusion; // terrain arena: ambient occlusion 0 open .. 1 enclosed (#132); other meshes: 0
// Terrain arena (terrainarena.h): this draw's record, instanced attributes picked by the command's baseInstance.
layout(location = 9) in vec4 arenaRow0;   // model rows 0-2 (render space)
layout(location = 10) in vec4 arenaRow1;
layout(location = 11) in vec4 arenaRow2;
layout(location = 12) in vec4 arenaTexOrigin;
layout(location = 13) in vec4 arenaClip;
layout(location = 6) in vec4 arenaMatW3;      // arena: palette slot weights 12-15 (#480)
layout(location = 15) in uvec4 arenaPalette;  // arena: the chunk's 16 material ids, a byte each

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
uniform int palette;            // lighting.fs's: 2 = planet (reads the palette slots), else the asteroid/fantasy base materials
uniform mat4 matViewProj;

out vec3 fragPosition;     // world (render) space: lighting, shadows, fog
out vec2 fragTexCoord;
out vec2 fragTexCoord2;
out vec4 fragColor;
out vec3 fragNormal;       // world space
out vec3 fragTexPos;       // body-aligned texturing position (wrapped): textures stay attached to a
out vec3 fragLocalNormal;  // moving/spinning body and don't swim when the floating origin shifts
out vec3 fragChunkPos;     // chunk-local position (overlap clipping)
centroid out float fragSky;  // 1 = open to the sky, 0 = buried (scales sunlight). Centroid (#402): with MSAA, samples past a triangle edge
                            // extrapolated it below 0 (negative light, NaN through the gamma pow)
centroid out float fragOcclusion; // #132 (centroid: as fragSky)
out vec2 fragCanopy;       // planets: distant forest canopy cover (x) and latitude (y)
flat out float fragClimate; // planets: the packed climate (2 + temperature step + moisture), per triangle (#254)
out float fragOreKind;     // planets: which ore (1 + kind) / 8, 0 none
out float fragMacro;       // slow (70-250 m) brightness variation for the planet palette (per vertex: #183)
flat out vec3 fragClipMax;
// #480: material weights per palette slot (16) and the slots' material ids (bytes; lighting.fs reads materials by id).
// Per-chunk draws (no arena) carry the 8 legacy weights with the identity palette 0-7.
out vec4 fragMatW[4];
flat out uvec4 fragPalette;
flat out int fragSlots;    // the palette slots in use (arena record clip.w): lighting.fs skips the unused groups

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
    fragTexCoord = vertexTexCoord.xy;
    fragTexCoord2 = vertexTexCoord2.xy;
    fragColor = vertexColor;
    if (arenaDraw == 1) {
        fragMatW[0] = vertexTexCoord; fragMatW[1] = vertexColor; fragMatW[2] = vertexTexCoord2; fragMatW[3] = arenaMatW3;
        fragPalette = arenaPalette;
        fragSlots = int(arenaClip.w + 0.5);
        // The 4 base materials' weights for the asteroid palettes (lighting.fs: rock, regolith, ice, ore in fragColor).
        vec4 base = vec4(0.0);
        if (palette != 2) // the planet palette reads the slots; skipping this loop there saves ~0.65 ms of near pass
        for (int s = 0; s < 16; ++s) {
            uint id = (arenaPalette[s >> 2] >> (8u * uint(s & 3))) & 255u;
            if (id < 4u) base[id] += fragMatW[s >> 2][s & 3];
        }
        fragColor = base;
        fragTexCoord = vec2(0.0);
        fragTexCoord2 = vec2(0.0);
    } else {
        fragMatW[0] = vertexColor; fragMatW[1] = vec4(vertexTexCoord.xy, vertexTexCoord2.xy);
        fragMatW[2] = vec4(0.0); fragMatW[3] = vec4(0.0);
        fragPalette = uvec4(0x03020100u, 0x07060504u, 0u, 0u);
        fragSlots = 8;
    }
    // w = 0: normals are directions and must not pick up the model translation.
    // Normal to world space via the model matrix. raylib also provides matNormal, but in testing its value
    // did not follow per-draw model rotations (normals stayed in body space on spinning bodies), so we
    // don't rely on it. Valid for rotation + uniform scale, which is all the engine uses.
    // The arena's packed vertex (#302): an octahedral normal, the climate scaled down by 66 (terrainarena.cpp Attr).
    vec3 nrm = vertexNormal;
    vec4 tng = vertexTangent;
    if (arenaDraw == 1) {
        vec2 e = vertexNormal.xy;
        nrm = vec3(e.x, 1.0 - abs(e.x) - abs(e.y), e.y);
        if (nrm.y < 0.0) nrm.xz = (1.0 - abs(nrm.zx)) * vec2(nrm.x >= 0.0 ? 1.0 : -1.0, nrm.z >= 0.0 ? 1.0 : -1.0);
        nrm = normalize(nrm);
        { float raw = floor(tng.y * 65535.0 + 0.5); tng.y = raw < 0.5 ? 0.0 : 1.0 + floor(raw / 1024.0) + mod(raw, 1024.0) / 1023.0; } // the climate's bit fields: 2 + step + moisture
    }
    fragNormal = normalize(mat3(model) * nrm);
    fragTexPos = vertexPosition + texOrigin;
    fragChunkPos = vertexPosition;
    fragSky = tng.w;
    fragOcclusion = vertexOcclusion;
    fragCanopy = tng.xy;
    fragClimate = tng.y;
    fragOreKind = tng.z;
    fragMacro = vnoise(fragTexPos / 70.0) * 0.6 + vnoise(fragTexPos / 250.0) * 0.4;
    fragLocalNormal = nrm;
    gl_Position = arenaDraw == 1 ? matViewProj * world : mvp * vec4(vertexPosition, 1.0);
}
