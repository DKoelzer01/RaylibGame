#version 330

// Input vertex attributes
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec4 vertexColor;

// Set by DrawMesh() from the material's shader locations
uniform mat4 mvp;
uniform mat4 matModel;
uniform mat4 matNormal;
// Set per chunk draw: the chunk's offset inside its body (voxels). Zero for other meshes.
uniform vec3 chunkOrigin;

out vec3 fragPosition;     // world (render) space: lighting, shadows, fog
out vec2 fragTexCoord;
out vec4 fragColor;
out vec3 fragNormal;       // world space
out vec3 fragLocalPos;     // body-local space: texturing, so textures stay attached to a
out vec3 fragLocalNormal;  // moving/spinning body and don't swim when the floating origin shifts

void main()
{
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    fragTexCoord = vertexTexCoord;
    fragColor = vertexColor;
    // w = 0: normals are directions and must not pick up the model translation.
    // Normal to world space via matModel. raylib also provides matNormal, but in testing its value
    // did not follow per-draw model rotations (normals stayed in body space on spinning bodies), so we
    // don't rely on it. Valid for rotation + uniform scale, which is all the engine uses.
    fragNormal = normalize(mat3(matModel) * vertexNormal);
    fragLocalPos = vertexPosition + chunkOrigin;
    fragLocalNormal = vertexNormal;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
