#version 330

// Distance sphere (farsphere.h): a unit sphere scaled to the body; the fragment shader looks the
// textures up by body-local direction.
in vec3 vertexPosition;

uniform mat4 mvp;
uniform mat4 matModel;

out vec3 fragLocal;     // body-local direction (unnormalised)
out vec3 fragPosition;  // render space

void main()
{
    fragLocal = vertexPosition;
    fragPosition = vec3(matModel * vec4(vertexPosition, 1.0));
    gl_Position = mvp * vec4(vertexPosition, 1.0);
}
