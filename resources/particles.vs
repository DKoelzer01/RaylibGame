#version 330

// Particles (particles.cpp, #2087): smoke, dust, debris and splash billboards in render space (rlgl's batch). The normal
// carries the particle's day light, not a direction: x = its sun-angle light times the sky it sees (0 in caves).
in vec3 vertexPosition;
in vec2 vertexTexCoord;
in vec3 vertexNormal;
in vec4 vertexColor;

uniform mat4 mvp;

out vec3 fragPosition;
out vec2 fragTexCoord;
out vec4 fragColor;
out float fragDay;

void main()
{
    fragPosition = vertexPosition;
    fragTexCoord = vertexTexCoord;
    fragColor = vertexColor;
    fragDay = clamp(vertexNormal.x, 0.0, 1.0);
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
