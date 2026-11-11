#version 330

// Ring rocks: grey-brown rock lit by the sun (none in the planet's shadow) and a little sky.
in vec3 fragPosition;
in vec3 fragNormal;
in vec3 fragTint;

uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec4 sunOccluder;

out vec4 finalColor;

void main()
{
    vec3 n = normalize(fragNormal);
    vec3 albedo = pow(fragTint, vec3(2.2)) * 0.9;
    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        if (c < 0.0 || (b < 0.0 && b * b > c)) sunVis = 0.0;
    }
    float direct = sunVis * max(dot(n, sunDir), 0.0);
    vec3 color = albedo * (ambient * 1.5 + vec3(0.015) + vec3(1.0, 0.97, 0.9) * direct);
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0);
}
