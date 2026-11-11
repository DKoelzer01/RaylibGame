#version 330

// Grass tufts (foliage.h): one shared five-blade tuft, instanced per chunk. Each instance places
// the tuft on the terrain (chunk-local position), sizes it, stands it along the body's local up,
// turns it by its hash and bends the blade tips in the wind.
layout(location = 0) in vec3 vertexPosition;   // tuft space: y = 0 root .. 1 tip
layout(location = 10) in vec4 instPosScale;    // chunk-local position, height (m)
layout(location = 11) in vec4 instUpHash;      // body-local up, hash in [0, 1)

uniform mat4 matModel;     // chunk -> render space (Planetoid::chunkModel)
uniform mat4 matViewProj;
uniform float time;

out vec3 fragPosition;     // render space
out vec3 fragUp;           // render space
out float fragHeight;      // 0 root .. 1 tip
out float fragHash;
out float fragLatitude;    // |body-local up.y|: 0 equator, 1 pole

void main()
{
    vec3 root = vec3(matModel * vec4(instPosScale.xyz, 1.0));
    vec3 up = normalize(mat3(matModel) * instUpHash.xyz);
    float h = instUpHash.w;
    // A tangent frame around up, turned by the hash.
    vec3 ref = abs(up.y) < 0.9 ? vec3(0.0, 1.0, 0.0) : vec3(1.0, 0.0, 0.0);
    vec3 t0 = normalize(cross(up, ref));
    vec3 b0 = cross(up, t0);
    float a = h * 6.2831853;
    vec3 t = t0 * cos(a) + b0 * sin(a);
    vec3 b = cross(up, t);
    float s = instPosScale.w;
    float y = vertexPosition.y;
    // Wind: a slow gust field over position plus a quick flutter, bending tips more than roots.
    float gust = sin(time * 1.3 + dot(root, vec3(0.11, 0.07, 0.13))) * 0.5 + 0.5;
    float flutter = sin(time * 4.0 + h * 40.0);
    vec3 windDir = normalize(t0 + b0 * 0.4);
    vec3 sway = windDir * (0.12 + 0.18 * gust + 0.04 * flutter) * y * y;
    vec3 p = root + (t * vertexPosition.x + up * y + b * vertexPosition.z + sway) * s;
    fragPosition = p;
    fragUp = up;
    fragHeight = y;
    fragHash = h;
    fragLatitude = abs(instUpHash.y);
    gl_Position = matViewProj * vec4(p, 1.0);
}
