#version 330

// InstanceBatch (instancebatch.h): unit cubes and tapered cylinders, one instanced draw per shape a frame.
// Per instance: a model matrix (render space), an sRGB colour, and params (x: top radius scale for
// cylinders, y: glow - unlit and bright).
layout(location = 0) in vec3 vertexPosition;   // cube: [-0.5, 0.5]^3; cylinder: radius 1 at y = 0, y in [0, 1]
layout(location = 2) in vec3 vertexNormal;
layout(location = 10) in vec4 instRow0;         // rows of the model matrix (row-major, 3x4)
layout(location = 11) in vec4 instRow1;
layout(location = 12) in vec4 instRow2;
layout(location = 13) in vec4 instColor;
layout(location = 14) in vec4 instParams;

uniform mat4 matViewProj;
uniform int shape;                               // 0 cube, 1 cylinder

out vec3 fragPosition;
out vec3 fragNormal;
out vec4 fragColor;
out float fragGlow;
out float fragSky;     // 1 open .. 0 buried (IB_SetBuried, #284)

void main()
{
    vec3 v = vertexPosition;
    vec3 n = vertexNormal;
    if (shape == 1) {
        float r = mix(1.0, instParams.x, v.y);         // taper: bottom radius 1, top instParams.x
        v.xz *= r;
        n = vec3(n.x, (1.0 - instParams.x) * length(n.xz), n.z) + vec3(0.0, n.y, 0.0); // the cone's slope
    }
    mat4 model = transpose(mat4(instRow0, instRow1, instRow2, vec4(0.0, 0.0, 0.0, 1.0)));
    vec4 world = model * vec4(v, 1.0);
    fragPosition = world.xyz;
    fragNormal = normalize(transpose(inverse(mat3(model))) * n); // non-uniform scales
    fragColor = instColor;
    fragGlow = instParams.y;
    fragSky = 1.0 - clamp(instParams.z, 0.0, 1.0);
    gl_Position = matViewProj * world;
}
