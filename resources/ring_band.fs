#version 330

// The ring band: the rock band (banded density, soft edges, gaps) inside a faint dust halo, lit by
// the sun from either side, dark in the planet's shadow, fading out as you get close to it (the
// rocks take over).
in vec3 fragPosition;
in float fragRadial;

uniform vec3 viewPos;
uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec4 sunOccluder;
uniform vec4 ringCentre;   // render-space centre, inner radius of the drawn band
uniform float outerRadius;
uniform vec2 coreRadii;    // the rock band
uniform float seed;

out vec4 finalColor;

float hash(float x) { return fract(sin(x * 127.1 + seed * 311.7) * 43758.5453); }
float vnoise(float x) { return mix(hash(floor(x)), hash(floor(x) + 1.0), smoothstep(0.0, 1.0, fract(x))); }

void main()
{
    float r = length(fragPosition.xz - ringCentre.xz);
    float core = (r - coreRadii.x) / (coreRadii.y - coreRadii.x);      // 0..1 across the rocks
    float inCore = smoothstep(-0.05, 0.05, core) * (1.0 - smoothstep(0.95, 1.05, core));
    float off = max(0.0, abs(r - 0.5 * (coreRadii.x + coreRadii.y)) - 0.5 * (coreRadii.y - coreRadii.x)); // m outside the rocks
    float halo = exp(-off / 18000.0);
    float bands = 0.5 * vnoise(r / 140.0) + 0.5 * vnoise(r / 47.0 + 50.0);
    float density = inCore * (0.35 + 0.65 * bands) + 0.35 * halo * (0.55 + 0.45 * vnoise(r / 2500.0 + 9.0));
    float edges = smoothstep(0.0, 0.05, fragRadial) * (1.0 - smoothstep(0.95, 1.0, fragRadial));
    density *= edges;
    float dist = length(fragPosition - viewPos);
    float fade = smoothstep(20000.0, 60000.0, dist);
    // Seen edge-on the band is thin; seen face-on it is a sheet: more opaque at grazing angles.
    vec3 v = (fragPosition - viewPos) / max(dist, 1.0);
    float grazing = 1.0 / max(abs(v.y), 0.08);
    float alpha = clamp(density * 0.22 * grazing, 0.0, 0.9) * fade;
    if (alpha < 0.003) discard;
    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        // Soft-edged: a penumbra 4% of the body's radius, widening behind it (a hard edge looked cut out).
        if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
    }
    // Lit from either face; the side away from the sun still glows through a little.
    float lit = sunVis * (0.35 + 0.65 * abs(sunDir.y));
    vec3 albedo = vec3(0.30, 0.27, 0.24);
    vec3 color = albedo * (ambient * 2.0 + vec3(1.0, 0.97, 0.9) * lit);
    finalColor = vec4(pow(color, vec3(1.0 / 2.2)), alpha);
}
