#version 330

// The ring as a band (ring.h): a unit annulus (x = cos, z = sin, y = 0 inner / 1 outer edge) scaled
// to the ring's radii around its planet, in the planet's equatorial plane.
layout(location = 0) in vec3 vertexPosition;

uniform mat4 matViewProj;
uniform vec4 ringCentre;   // render-space centre, inner radius
uniform float outerRadius;

out vec3 fragPosition;
out float fragRadial;      // 0 inner edge .. 1 outer edge

void main()
{
    float r = mix(ringCentre.w, outerRadius, vertexPosition.y);
    fragPosition = ringCentre.xyz + vec3(vertexPosition.x * r, 0.0, vertexPosition.z * r);
    fragRadial = vertexPosition.y;
    gl_Position = matViewProj * vec4(fragPosition, 1.0);
}
