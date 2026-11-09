#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;   // chunks: material weights 4-5 (sand, gravel)
in vec2 vertexTexCoord2;  // chunks: material weights 6-7 (grass, snow)
in vec3 vertexNormal;
in vec4 vertexColor;

// Set by DrawMesh() from the material's shader locations
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
// Set per chunk draw: the chunk's body-local origin wrapped to the texture repeat (so texture
// coordinates stay small and precise even 60 km from a body's centre). Zero for other meshes.
uniform vec3 chunkTexOrigin;

out vec3 fragPosition;     // world (render) space: lighting, shadows, fog
out vec2 fragTexCoord;
out vec2 fragTexCoord2;
out vec4 fragColor;
out vec3 fragNormal;       // world space
out vec3 fragTexPos;       // body-aligned texturing position (wrapped): textures stay attached to a
out vec3 fragLocalNormal;  // moving/spinning body and don't swim when the floating origin shifts
out vec3 fragChunkPos;     // chunk-local position (overlap clipping)

void main()
{
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    fragTexCoord = vertexTexCoord;
    fragTexCoord2 = vertexTexCoord2;
    fragColor = vertexColor;
    // w = 0: normals are directions and must not pick up the model translation.
    // Normal to world space via matModel. raylib also provides matNormal, but in testing its value
    // did not follow per-draw model rotations (normals stayed in body space on spinning bodies), so we
    // don't rely on it. Valid for rotation + uniform scale, which is all the engine uses.
    fragNormal = normalize(mat3(matModel) * vertexNormal);
    fragTexPos = vertexPosition + chunkTexOrigin;
    fragChunkPos = vertexPosition;
    fragLocalNormal = vertexNormal;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
