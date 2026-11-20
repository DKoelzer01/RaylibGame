#version 330
in vec3 vertexPosition;
// Terrain arena (terrainarena.h): this draw's record, instanced attributes picked by the command's baseInstance.
layout(location = 9) in vec4 arenaRow0;   // model rows 0-2 (render space)
layout(location = 10) in vec4 arenaRow1;
layout(location = 11) in vec4 arenaRow2;
layout(location = 12) in vec4 arenaTexOrigin;
layout(location = 13) in vec4 arenaClip;
uniform mat4 mvp;   // set by DrawMesh from the light camera's matrices
// Terrain arena (Plan 6 step 1): the model comes from the arena attributes.
uniform int arenaDraw;
uniform mat4 matViewProj;
void main() {
    if (arenaDraw == 1) {
        mat4 model = transpose(mat4(arenaRow0, arenaRow1, arenaRow2, vec4(0.0, 0.0, 0.0, 1.0)));
        gl_Position = matViewProj * (model * vec4(vertexPosition, 1.0));
    } else {
        gl_Position = mvp * vec4(vertexPosition, 1.0);
    }
}
