#version 330

// Input vertex attributes
in vec3 vertexPosition;

// Input uniform values
uniform mat4 matProjection;
uniform mat4 matView;

// Output vertex attributes (to fragment shader)
out vec3 fragPosition;

void main()
{
    // Calculate fragment position based on model transformations
    fragPosition = vertexPosition;

    // Remove translation from the view matrix
    mat4 rotView = mat4(mat3(matView));
    // Scale the unit cube well past the near plane (0.5); translation is already removed.
    vec4 clipPos = matProjection*rotView*vec4(vertexPosition*10.0, 1.0);

    // Calculate final vertex position
    gl_Position = clipPos;
}
