#version 330
// Waterfall sheets and their foam (#540). Like the other world shaders it writes display-encoded colour.
in vec3 fragPosition;
in vec3 fragNormal;
in vec2 fragUV;
in vec4 fragColor;

uniform float time;
uniform vec3 sunDir;      // towards the sun
uniform vec3 ambient;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec3 viewPos;

out vec4 finalColor;

float hash(vec2 p) { return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453); }
float noise(vec2 p) {
    vec2 i = floor(p), f = fract(p);
    vec2 u = f * f * (3.0 - 2.0 * f);
    return mix(mix(hash(i), hash(i + vec2(1, 0)), u.x), mix(hash(i + vec2(0, 1)), hash(i + vec2(1, 1)), u.x), u.y);
}

void main() {
    vec3 n = normalize(fragNormal);
    vec3 v = normalize(viewPos - fragPosition);
    float light = 0.45 + 0.55 * max(dot(n, sunDir), 0.0) + 0.25 * max(dot(-n, sunDir), 0.0); // thin water: lit through
    vec3 water = vec3(0.55, 0.70, 0.78), foam = vec3(0.93, 0.96, 0.98);
    vec3 col;
    float alpha;
    if (fragColor.a > 0.5) {
        // The sheet: streaks running down, faster the further the water has fallen (v grows ~ t^2).
        // The pattern rides on the time water takes to fall this far (thrown at v0, pulled by g), less the clock: a
        // feature then always moves down, faster as it falls. (Scrolling by fallen - time * speed(fallen) ran upwards
        // once the clock grew: the speed term's slope along the sheet outgrew the distance term's.)
        float fallen = max(fragUV.y, 0.0), s = fragColor.r;
        const float v0 = 2.2, g = 9.8;
        float tFall = (sqrt(v0 * v0 + 2.0 * g * fallen) - v0) / g; // s
        float phase = tFall - time;
        float streak = noise(vec2(fragUV.x * 22.0, phase * 3.0));
        float fine = noise(vec2(fragUV.x * 61.0, phase * 9.0));
        float white = clamp(0.35 + 0.55 * streak + 0.25 * fine + 0.6 * smoothstep(0.75, 1.0, s) + 0.3 * (1.0 - smoothstep(0.0, 0.12, s)), 0.0, 1.0);
        col = mix(water, foam, white);
        float edge = smoothstep(0.0, 0.12, fragUV.x) * smoothstep(1.0, 0.88, fragUV.x);      // ragged-ish banks
        edge *= 0.75 + 0.25 * noise(vec2(fragUV.x * 9.0, phase * 2.0));
        alpha = clamp(0.55 + 0.4 * white, 0.0, 0.95) * edge * smoothstep(0.0, 0.03, s + 0.03);
    } else {
        // The foam where it lands: churning rings fading out from the middle.
        vec2 c = fragUV - 0.5;
        float r = length(c) * 2.0;
        float churn = noise(c * 9.0 + vec2(time * 0.7, -time * 0.5)) * 0.6 + noise(c * 23.0 - time * 1.3) * 0.4;
        col = foam;
        alpha = (1.0 - smoothstep(0.35, 1.0, r)) * (0.45 + 0.5 * churn) * 0.85;
    }
    col *= ambient + vec3(light);
    col = clamp(col, 0.0, 1.2);
    float dist = length(viewPos - fragPosition);
    float fogF = 1.0 - exp(-dist * fogDensity);
    col = mix(col, fogColor, clamp(fogF, 0.0, 1.0));
    finalColor = vec4(col, alpha);
}
