#version 330

// Ring rocks (ring.h): one of a few lumpy rock meshes, instanced. Each instance is placed in render
// space, scaled to its mean radius, turned by its spin (a quaternion) and tinted.
layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexNormal;
layout(location = 10) in vec4 instPosScale; // render-space centre, mean radius
layout(location = 11) in vec4 instRot;      // quaternion
layout(location = 12) in vec4 instTint;

uniform mat4 matViewProj;

out vec3 fragPosition;
out vec3 fragNormal;
out vec3 fragTint;

vec3 rotate(vec4 q, vec3 v) { return v + 2.0 * cross(q.xyz, cross(q.xyz, v) + q.w * v); }

void main()
{
    fragPosition = instPosScale.xyz + rotate(instRot, vertexPosition * instPosScale.w);
    fragNormal = rotate(instRot, vertexNormal);
    fragTint = instTint.rgb;
    gl_Position = matViewProj * vec4(fragPosition, 1.0);
}
