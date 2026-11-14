#version 330

// Parts: sRGB vertex colour with a procedural surface pattern (the id in the vertex alpha, set by
// plots.cpp PartPatternOf), drawn in plot-local metres so it stays fixed to the 0.5 m build grid:
// cut stone, brick, planks, metal panels with rivets, glass, tiles, polished stone and marble, logs,
// gems, machine casings with vents, conveyor belts, cables, farmland, rocket hull plating. Then sun
// light with shadows (shadow.glsl) and the planet's night side, the headlamp, sky ambient, a
// per-pattern specular, and the terrain's haze. Pattern detail fades out once a pixel covers more
// than a few cm (no shimmer far away).
in vec3 fragPosition;
in vec3 fragNormal;
in vec4 fragColor;
in vec3 fragLocal;
in vec3 fragLocalNormal;
in float fragSky;         // 1 open .. 0 buried (parts.vs; #216: parts underground took full sky light)

uniform vec3 viewPos;
uniform vec3 sunDir;
uniform vec3 ambient;
uniform vec4 sunOccluder;
uniform vec3 fogColor;
uniform float fogDensity;
uniform vec4 fogShell;

#include "shadow.glsl"
#include "hemisphere.glsl"

out vec4 finalColor;

float h21(vec2 p) { p = mod(p, 289.0); return fract(sin(dot(p, vec2(127.1, 311.7))) * 43758.5453); }
float vn(vec2 p)
{
    vec2 i = floor(p), f = fract(p);
    f = f * f * (3.0 - 2.0 * f);
    return mix(mix(h21(i), h21(i + vec2(1.0, 0.0)), f.x), mix(h21(i + vec2(0.0, 1.0)), h21(i + vec2(1.0, 1.0)), f.x), f.y);
}
float fbm(vec2 p) { return 0.5 * vn(p) + 0.3 * vn(p * 2.13 + 7.1) + 0.2 * vn(p * 4.37 + 3.3); }
// 0 inside a groove of half-width w (m) at distance d from its centre line, 1 outside; aa = m per pixel.
float groove(float d, float w, float aa) { return smoothstep(w, w + aa, d); }
// Distance to the nearest edge of a running-bond block pattern (rows of height rh, blocks of length bl,
// odd rows shifted half a block); cell = the block's id.
float bond(vec2 uv, float bl, float rh, out vec2 cell)
{
    float row = floor(uv.y / rh);
    float u = uv.x + mod(row, 2.0) * bl * 0.5;
    cell = vec2(floor(u / bl), row);
    float fx = fract(u / bl) * bl, fy = fract(uv.y / rh) * rh;
    return min(min(fx, bl - fx), min(fy, rh - fy));
}
float gridEdge(vec2 uv, float period) { vec2 g = abs(fract(uv / period + 0.5) - 0.5) * period; return min(g.x, g.y); }
float rivets(vec2 uv, float period, float inset, float r)
{
    vec2 q = abs(fract(uv / period) - 0.5) * period;   // from the panel centre
    vec2 c = vec2(period * 0.5 - inset) - q;           // to the rivet at the nearest corner
    return 1.0 - smoothstep(r * 0.6, r, length(c));
}

// Pattern -> albedo (linear), specular strength, emission.
vec3 surface(int pat, vec3 base, vec2 uv, vec3 ln, float aa, out float spec, out float emit)
{
    spec = 0.08; emit = 0.0;
    vec2 cell;
    if (pat == 1) {           // cut stone: 0.5 x 0.25 m blocks
        float e = bond(uv, 0.5, 0.25, cell);
        vec3 a = base * (0.84 + 0.3 * h21(cell)) * (0.88 + 0.24 * fbm(uv * 9.0));
        return mix(base * 0.5, a * (e < 0.03 ? 1.08 : 1.0), groove(e, 0.008, aa));
    }
    if (pat == 2) {           // metal plate: 0.5 m panels, seams, rivets, brushed
        float e = gridEdge(uv, 0.5);
        vec3 a = base * (0.9 + 0.12 * fbm(vec2(uv.x * 3.0, uv.y * 90.0)));
        a = mix(a * 0.45, a, groove(e, 0.004, aa));
        a = mix(a, base * 1.25, rivets(uv, 0.5, 0.035, 0.012));
        spec = 0.55;
        return a;
    }
    if (pat == 3) {           // glass: a thin frame per cell, streaks of reflection
        float e = gridEdge(uv, 0.5);
        float s = fract((uv.x + uv.y) * 0.9);
        float streak = smoothstep(0.05, 0.0, abs(s - 0.35)) * 0.3 + smoothstep(0.02, 0.0, abs(s - 0.45)) * 0.2;
        spec = 0.9;
        return mix(base * 0.45, base * 1.05 + vec3(streak), groove(e, 0.018, aa));
    }
    if (pat == 4 || pat == 12) {   // packed dirt, clay
        float k = pat == 4 ? 0.5 : 0.18;
        return base * (1.0 - k * 0.5 + k * fbm(uv * 14.0)) * (1.0 + (pat == 4 ? 0.2 : 0.05) * (h21(floor(uv * 60.0)) - 0.5));
    }
    if (pat == 5) {           // heat shield tiles: 0.25 m, dark gaps
        float e = gridEdge(uv, 0.25);
        cell = floor(uv / 0.25);
        spec = 0.25;
        return mix(vec3(0.01), base * (0.8 + 0.35 * h21(cell)), groove(e, 0.006, aa));
    }
    if (pat == 6) {           // planks: 12.5 cm boards along u, ends staggered
        float row = floor(uv.y / 0.125);
        float u = uv.x + h21(vec2(row, 3.0));
        float fx = fract(u), fy = fract(uv.y / 0.125) * 0.125;
        float e = min(min(fx, 1.0 - fx), min(fy, 0.125 - fy));
        float grain = fbm(vec2(u * 2.5, uv.y * 70.0 + row * 13.0));
        vec3 a = base * (0.86 + 0.24 * h21(vec2(row, floor(u)))) * (0.8 + 0.4 * grain);
        float nail = 1.0 - smoothstep(0.004, 0.007, length(vec2(min(fx, 1.0 - fx) - 0.04, fy - 0.0625)));
        a = mix(a, vec3(0.08), nail * 0.8);
        return mix(base * 0.35, a, groove(e, 0.003, aa));
    }
    if (pat == 7) {           // granite, basalt, sandstone: speckled slabs, 1 m seams, sandstone layers
        float e = gridEdge(uv, 1.0);
        vec3 a = base * (0.78 + 0.35 * h21(floor(uv * 90.0))) * (0.9 + 0.2 * fbm(uv * 4.0));
        a *= 1.0 + 0.08 * sin(uv.y * 28.0 + fbm(uv * 2.0) * 5.0);
        spec = 0.3;
        return mix(base * 0.55, a, groove(e, 0.003, aa));
    }
    if (pat == 16) {          // marble, limestone: soft veins, 1 m slabs
        float e = gridEdge(uv, 1.0);
        float v = abs(sin((uv.x * 0.8 + uv.y) * 3.5 + fbm(uv * 2.5) * 7.0));
        vec3 a = base * (0.95 + 0.08 * fbm(uv * 6.0));
        a = mix(a * 0.62, a, smoothstep(0.0, 0.12, v));
        spec = 0.5;
        return mix(base * 0.6, a, groove(e, 0.002, aa));
    }
    if (pat == 8) {           // brick: 25 x 8.3 cm, light mortar
        float e = bond(uv, 0.25, 0.5 / 6.0, cell);
        vec3 a = base * (0.8 + 0.35 * h21(cell)) * (0.9 + 0.2 * fbm(uv * 25.0));
        return mix(vec3(0.42, 0.4, 0.37), a, groove(e, 0.006, aa));
    }
    if (pat == 9) return base * (0.95 + 0.08 * fbm(uv * 5.0));   // plaster
    if (pat == 10) {          // polished gem blocks: faceted 0.25 m tiles, sparkles
        float e = gridEdge(uv, 0.25);
        vec2 f = fract(uv / 0.25);
        float facet = (f.x + f.y < 1.0 ? 1.12 : 0.88) * (f.x > f.y ? 1.04 : 0.96);
        vec3 a = base * facet;
        a += vec3(step(0.985, h21(floor(uv * 80.0)))) * 0.6;
        spec = 1.0;
        return mix(base * 0.5, a, groove(e, 0.004, aa));
    }
    if (pat == 11) {          // logs: bark on the sides, growth rings on the ends
        if (abs(ln.y) < 0.5) {
            float b = fbm(vec2(uv.x * 16.0, uv.y * 1.6));
            return base * (0.55 + 0.75 * b) * (1.0 - 0.35 * smoothstep(0.32, 0.25, b));
        }
        float r = length(fract(uv / 0.5) - 0.5) * 0.5;
        vec3 wood = mix(base * 1.8, vec3(0.62, 0.45, 0.27), 0.6);
        return r > 0.235 ? base : wood * (0.85 + 0.15 * sin(r * 160.0 + fbm(uv * 8.0) * 3.0));
    }
    if (pat == 13) {          // ice: pale streaks
        spec = 0.9;
        return base * 1.05 + vec3(0.12) * smoothstep(0.08, 0.0, abs(fract(uv.x * 1.3 + fbm(uv * 3.0)) - 0.5));
    }
    if (pat == 14) {          // coal: dark, glinting
        spec = 0.4;
        return base * (0.7 + 0.6 * fbm(uv * 10.0)) + vec3(step(0.97, h21(floor(uv * 70.0)))) * 0.2;
    }
    if (pat == 15) { emit = 1.0; return base; }   // lamp
    if (pat == 20 || pat == 24) {   // machine casing / rocket hull: panels, bevels, rivets, vents
        float period = pat == 24 ? 1.0 : 0.5;
        float e = gridEdge(uv, period);
        vec3 a = base * (0.92 + 0.1 * fbm(uv * 6.0));
        a *= e < 0.02 ? 1.1 : 1.0;                             // bevelled panel edges catch the light
        a = mix(a * 0.4, a, groove(e, 0.004, aa));
        a = mix(a, base * 1.3, rivets(uv, period, 0.03, 0.01));
        if (pat == 20 && abs(ln.y) < 0.5) {
            // A louvred vent across the middle of each 1 m of height, on the sides.
            float fy = fract(fragLocal.y), fx = fract(uv.x);
            if (fy > 0.58 && fy < 0.8 && fx > 0.15 && fx < 0.85)
                a = mix(a * 0.3, a, groove(abs(fract(uv.y / 0.035) - 0.5) * 0.035, 0.007, aa));
        }
        spec = pat == 24 ? 0.45 : 0.3;
        return a;
    }
    if (pat == 21) {          // conveyor belt: ridges across the run
        return mix(base * 0.5, base * (0.95 + 0.1 * fbm(uv * 20.0)), groove(abs(fract(uv.x / 0.06) - 0.5) * 0.06, 0.008, aa));
    }
    if (pat == 22) { spec = 0.5; return base * (0.85 + 0.2 * sin((uv.x + uv.y) * 90.0)); }   // twisted copper
    if (pat == 23) {          // farmland: furrows every 25 cm
        float f = abs(fract(uv.y / 0.25) - 0.5) * 2.0;
        return base * (0.65 + 0.45 * f) * (0.85 + 0.3 * fbm(uv * 20.0));
    }
    // Concrete (0, and anything unknown): mottled, faint 1 m formwork lines, a few pits.
    float e = gridEdge(uv, 1.0);
    vec3 a = base * (0.9 + 0.18 * fbm(uv * 6.0));
    a *= 1.0 - 0.18 * step(0.975, h21(floor(uv * 40.0)));
    return mix(a * 0.8, a, groove(e, 0.003, aa));
}

void main()
{
    vec3 base = pow(fragColor.rgb, vec3(2.2));
    int pat = int(fragColor.a * 255.0 + 0.5);
    vec3 ln = normalize(fragLocalNormal);
    vec3 an = abs(ln);
    vec2 uv = (an.x > an.y && an.x > an.z) ? fragLocal.zy : (an.y > an.z ? fragLocal.xz : fragLocal.xy);
    float aa = max(max(fwidth(uv.x), fwidth(uv.y)), 1e-5);
    float spec, emit;
    vec3 albedo = surface(pat, base, uv, ln, aa, spec, emit);
    albedo = mix(albedo, base, clamp((aa - 0.01) / 0.03, 0.0, 1.0)); // detail fades past ~1 cm per pixel

    vec3 n = normalize(fragNormal);
    float sunVis = 1.0;
    if (sunOccluder.w > 0.0) {
        vec3 oc = fragPosition - sunOccluder.xyz;
        float b = dot(oc, sunDir), c = dot(oc, oc) - sunOccluder.w * sunOccluder.w;
        // Soft-edged: a penumbra 4% of the body's radius, widening behind it (a hard edge looked cut out).
        if (b < 0.0) sunVis *= smoothstep(sunOccluder.w - max(sunOccluder.w * 0.04, -b * 0.03), sunOccluder.w + max(sunOccluder.w * 0.04, -b * 0.03), sqrt(max(c + sunOccluder.w * sunOccluder.w - b * b, 0.0)));
    }
    float ndl = max(dot(n, sunDir), 0.0);
    float sun = sunVis * ndl * (ndl > 0.0 ? sunShadowNormal(fragPosition, n, sunDir, 4) : 0.0) * fragSky;
    vec3 v = normalize(viewPos - fragPosition);
    float specular = pow(max(dot(n, normalize(sunDir + v)), 0.0), mix(16.0, 96.0, spec)) * spec * 0.6 * sun;
    vec3 color = albedo * ((hemiOn() ? hemiAmbient(ambient, n) * 1.4 : ambient * (1.4 + 0.3 * n.y)) * fragSky + vec3(1.0, 0.96, 0.88) * sun + lampLight(fragPosition, n)) + vec3(specular)
               + albedo * emit * 1.2;
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
