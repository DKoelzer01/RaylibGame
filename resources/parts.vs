#version 330

// Build plot sections (plots.h): merged part meshes with flat normals and per-face vertex colours
// (alpha = the surface pattern id, see parts.fs). Positions are section-local plot space, so the
// patterns stay fixed to the 0.5 m build grid.
in vec3 vertexPosition;
in vec3 vertexNormal;
in vec4 vertexColor;
in vec2 vertexTexCoord;   // x: how buried the cell is (1 - sky exposure, plots.cpp; #216); meshes without it: 0 = open

uniform mat4 mvp;
uniform mat4 matModel;

out vec3 fragPosition;
out vec3 fragNormal;
out vec4 fragColor;
out vec3 fragLocal;
out vec3 fragLocalNormal;
out float fragSky;        // 1 open to the sky .. 0 buried ~10 m+ (scales ambient and sun, as for the terrain)

void main()
{
    fragSky = 1.0 - clamp(vertexTexCoord.x, 0.0, 1.0);
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    fragNormal = normalize(mat3(matModel) * vertexNormal);
    fragColor = vertexColor;
    fragLocal = vertexPosition;
    fragLocalNormal = vertexNormal;
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
