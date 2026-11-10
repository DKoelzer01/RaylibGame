#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;

uniform mat4 mvp;
uniform mat4 matModel;
uniform vec3 chunkTexOrigin; // chunk origin wrapped to the texture repeat (as in lighting.vs)

out vec3 fragPosition;  // world (render) space
out vec3 fragNormal;    // world space
out vec3 fragTexPos;    // body-aligned, wrapped (ripples stay put and seamless)
out vec3 fragChunkPos;  // chunk-local (overlap clipping)

void main()
{
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    fragNormal = normalize(mat3(matModel) * vertexNormal);
    fragTexPos = vertexPosition + chunkTexOrigin;
    fragChunkPos = vertexPosition;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
