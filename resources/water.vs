#version 330

in vec3 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexTexCoord;     // x: depth to the floor (m), y: fluid kind (fluids.txt; 0 water)
in vec4 vertexColor;        // the current (#361): rgb = body-local direction * 0.5 + 0.5, a = speed / 4 m/s

uniform mat4 mvp;
uniform mat4 matModel;
uniform vec3 chunkTexOrigin; // chunk origin wrapped to the texture repeat (as in lighting.vs)

out vec3 fragPosition;  // world (render) space
out vec3 fragNormal;    // world space
out vec3 fragTexPos;    // body-aligned, wrapped (ripples stay put and seamless)
out vec3 fragChunkPos;  // chunk-local (overlap clipping)
out float fragDepth;     // depth to the sea floor (m)
flat out int fragKind;   // fluid kind
out vec4 fragFlow;       // the current: xyz direction (body-aligned, as fragTexPos), w speed (m/s)
out vec3 fragFlowW;      // the current's direction in world (render) space, as fragNormal: for tilting the normal
out float fragGen;       // (#1164) 1: undisturbed generated water (sea, lake): its shores foam

void main()
{
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    fragNormal = normalize(mat3(matModel) * vertexNormal);
    fragTexPos = vertexPosition + chunkTexOrigin;
    fragChunkPos = vertexPosition;
    fragDepth = vertexTexCoord.x;
    fragKind = int(vertexTexCoord.y + 0.5);
    fragGen = fract(vertexTexCoord.y) > 0.1 ? 1.0 : 0.0; // the mesher adds 0.25 for generated water
    fragFlow = vec4(vertexColor.rgb * 2.0 - 1.0, vertexColor.a * 4.0);
    fragFlowW = mat3(matModel) * fragFlow.xyz;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
