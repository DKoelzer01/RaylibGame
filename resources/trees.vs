#version 330

// Trees and bushes (trees.h): one prism-built variant mesh, instanced. Each instance stands the
// variant on the body's local up at its position (relative to the batch anchor), turns it by its
// hash and scales it; leaves sway a little in the wind, more toward the top.
layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexNormal;
layout(location = 2) in vec4 vertexColor;      // r: leaf (1) / bark (0), g: shade variation
layout(location = 10) in vec4 instPosScale;    // anchor-relative body-local position, scale
layout(location = 11) in vec4 instUpHash;      // body-local up, hash in [0, 1)

uniform mat4 matModel;     // anchor -> render space
uniform mat4 matViewProj;
uniform float time;

out vec3 fragPosition;
out vec3 fragNormal;
out float fragLeaf;
out float fragShade;
out float fragStone;
out float fragHash;
out float fragLatitude;
out float fragHeight;

void main()
{
    vec3 up = normalize(instUpHash.xyz);
    float h = instUpHash.w;
    vec3 ref = abs(up.y) < 0.9 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 t0 = normalize(cross(up, ref));
    vec3 b0 = cross(up, t0);
    float a = h * 6.2831853;
    vec3 t = t0 * cos(a) + b0 * sin(a);
    vec3 b = cross(t, up);   // right-handed (t, up, b): keeps the mesh winding (front faces CCW)
    float s = instPosScale.w;
    vec3 v = vertexPosition;
    // Wind on leaves: a slow sway growing with height.
    float leaf = vertexColor.r;
    vec3 sway = (t0 * sin(time * 0.9 + h * 30.0) + b0 * sin(time * 0.7 + h * 17.0)) * 0.04 * max(v.y, 0.0) * leaf;
    vec3 local = instPosScale.xyz + (t * v.x + up * v.y + b * v.z) * s + sway * s;
    vec3 n = t * vertexNormal.x + up * vertexNormal.y + b * vertexNormal.z;
    fragPosition = vec3(matModel * vec4(local, 1.0));
    fragNormal = normalize(mat3(matModel) * n);
    fragLeaf = leaf;
    fragShade = vertexColor.g;
    fragStone = vertexColor.b;  // loose stones: grey rock
    fragHash = h;
    fragLatitude = abs(up.y);
    fragHeight = v.y;
    gl_Position = matViewProj * vec4(fragPosition, 1.0);
}
