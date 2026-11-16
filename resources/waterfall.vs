#version 330
// Waterfall sheets and their foam (#540, src/space/waterfalls.cpp): vertices come in render space.
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec2 vertexTexCoord;   // sheet: (across 0..1, metres fallen); foam disc: (0..1, 0..1) round its centre
in vec4 vertexColor;      // r: the fraction fallen (sheet); a: 1 sheet, 0 foam disc

uniform mat4 mvp;

out vec3 fragPosition;
out vec3 fragNormal;
out vec2 fragUV;
out vec4 fragColor;

void main() {
    fragPosition = vertexPosition;
    fragNormal = vertexNormal;
    fragUV = vertexTexCoord;
    fragColor = vertexColor;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
