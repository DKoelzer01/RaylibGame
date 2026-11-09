#version 330
in vec3 vertexPosition;
uniform mat4 mvp;   // set by DrawMesh from the light camera's matrices
void main() {
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
