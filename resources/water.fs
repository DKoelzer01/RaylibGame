#version 330

// Sea surface: translucent deep blue, reflecting the sky at grazing angles (Fresnel), small ripples,
// a sun glint, the planet's own shadow at night, and the same distance fog as the terrain.
in vec3 fragPosition;
in vec3 fragNormal;
in vec3 fragTexPos;
in vec3 fragChunkPos;
in float fragDepth;            // depth to the sea floor (m)

uniform vec3 viewPos;
uniform vec3 sunDir;           // towards the sun
uniform vec3 skyColor;         // linear horizon sky colour (dark at night)
uniform vec3 chunkClipMax;     // chunk-local overlap clip (as for terrain)
uniform vec4 sunOccluder;      // planet: centre + radius blocking the sun (w = 0: none)
uniform vec3 fogColor;         // linear
uniform float fogDensity;
uniform float time;

out vec4 finalColor;

void main()
{
    if (any(greaterThan(fragChunkPos, chunkClipMax))) discard;
    vec3 n = normalize(fragNormal);
    // Ripples: a few sine waves on the wrapped body-local position. Their wavelengths divide the
    // 20 m wrap period, so neighbouring chunks line up.
    const float k = 6.2831853 / 20.0;
    vec3 p = fragTexPos;
    float dist = length(viewPos - fragPosition);
    float rippleFade = 1.0 / (1.0 + dist / 60.0); // ripples finer than a pixel far away only alias
    vec3 bump = vec3(sin(p.x * k * 7.0 + time * 1.3) + sin(p.z * k * 5.0 - time * 1.1),
                     0.0,
                     sin(p.z * k * 9.0 + time * 0.9) + sin(p.y * k * 6.0 + time * 1.7)) * 0.025 * rippleFade;
    n = normalize(n + bump - n * dot(bump, n));

    vec3 v = normalize(viewPos - fragPosition);
    float nv = abs(dot(n, v));                                    // also seen from below
    float fresnel = 0.02 + 0.98 * pow(1.0 - nv, 5.0);

    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        if (c < 0.0 || (b < 0.0 && b * b > c)) sunVis = 0.0;
    }
    float light = sunVis * max(dot(n, sunDir), 0.0);
    vec3 deep = vec3(0.01, 0.05, 0.13) * (0.15 + 0.85 * light); // dark but clearly blue, also from orbit
    // Clear turquoise in the shallows, deep blue offshore.
    float clarity = exp(-fragDepth / 5.0);
    vec3 shallow = vec3(0.03, 0.22, 0.24) * (0.15 + 0.85 * light);
    vec3 color = mix(mix(deep, shallow, clarity), skyColor, fresnel);
    vec3 h = normalize(sunDir + v);
    float glint = pow(max(dot(n, h), 0.0), 800.0) * 1.2 * sunVis;
    color += vec3(1.0, 0.95, 0.85) * glint;
    // Translucent: see into the shallows looking down, mostly sky reflection at grazing angles
    // (Fresnel). Seen from below, the surface is a bright, mostly opaque ceiling. The scene's stencil
    // blends one water layer per pixel, so overlapping LOD copies don't double up.
    bool below = dot(normalize(fragNormal), v) < 0.0;
    float alpha = below ? 0.8 : clamp(max(mix(1.0, 0.3, clarity), fresnel) + glint, 0.0, 1.0);

    float fog = 1.0 - exp(-fogDensity * dist);
    color = mix(color, fogColor, clamp(fog, 0.0, 1.0));
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), alpha);
}
