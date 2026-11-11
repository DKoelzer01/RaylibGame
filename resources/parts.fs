#version 330

// Parts: sRGB vertex colour, sun light with the planet's night side, sky ambient, the terrain's haze.
in vec3 fragPosition;
in vec3 fragNormal;
in vec4 fragColor;

uniform vec3 viewPos;
uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec4 sunOccluder;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;

out vec4 finalColor;

void main()
{
    vec3 albedo = pow(fragColor.rgb, vec3(2.2));
    vec3 n = normalize(fragNormal);
    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        if (c < 0.0 || (b < 0.0 && b * b > c)) sunVis = 0.0;
    }
    float direct = sunVis * max(dot(n, sunDir), 0.0);
    vec3 color = albedo * (ambient * (1.4 + 0.3 * n.y) + vec3(1.0, 0.96, 0.88) * direct);
    float dist = length(viewPos - fragPosition), fogDist = dist;
    if (fogShell.w > 0.0) {
        vec3 rd = (fragPosition - viewPos) / max(dist, 1e-3);
        vec3 oc = viewPos - fogShell.xyz;
        float b = dot(oc, rd), cc = dot(oc, oc) - fogShell.w * fogShell.w, disc = b * b - cc;
        float t0 = disc > 0.0 ? -b - sqrt(disc) : 0.0, t1 = disc > 0.0 ? -b + sqrt(disc) : 0.0;
        fogDist = max(0.0, min(t1, fogDist) - max(t0, 0.0));
    }
    color = mix(color, fogColor, clamp(1.0 - exp(-fogDensity * fogDist), 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), 1.0);
}
