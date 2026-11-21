#version 330

// Sky of the planet the camera is near, blended over the star skybox: blue overhead, paler at the
// horizon, dusk colours around a low sun, a sun disc, fading to the stars at night and in space.
in vec3 fragPosition;          // view direction (cube vertex, from skybox.vs)

uniform vec3 sunDir;           // towards the sun
uniform vec3 upDir;            // away from the planet centre, at the camera
uniform float density;
uniform float overcast;        // weather: 0 clear .. 1 a thick grey sky (no sun disc)
// Clouds: a shell over the planet, sampled in the planet's own frame (they stay over the ground
// and drift with the wind). Cover never drops to nothing: a few clouds are always out.
uniform vec3 camLocal;         // camera, planet-local (m)
uniform vec3 toLocalX, toLocalY, toLocalZ; // the planet-local images of the world axes
uniform float cloudShell;      // the cloud layer's radius (m); 0 = no clouds
uniform vec3 cloudDrift;       // wind drift so far (m, planet-local)
uniform float cloudCover;      // 0..1 (weather), with a floor applied here
uniform vec4 lightning;        // (#1162) xyz: toward the strike (world), w: the flash 0..1 - lights the clouds and the sky there
// #142: the physical sky (on by default; RAYLIBGAME_SKY_PHYS=0: off): the sky-view LUT (skyphys.cpp), pulled toward the palette below by skyPhysTint.
uniform float skyPhys;         // 1: on (the LUT is bound)
#include "skyphys.glsl"
// #142 later slices: the transmittance LUT for the limb seen from orbit and for sunlight on the clouds; the top of the
// planet's air (m; with skyRealHorizon.x its radius); and which of them are on (x: limb, y: clouds lit from the LUTs).
uniform sampler2D skyTransLut;
uniform float skyAirTop;
uniform vec2 skyExtras;

// Height in the planet's air (m) -> in the LUTs' Earth-sized shell (km), as the camera's (SkyPhys_LutHeight).
float ShellHeight(float hm) { return clamp(hm / max(skyAirTop - skyRealHorizon.x, 1.0), 0.0, 1.0) * (skyRt - skyR); }
// Sunlight reaching a planet-local point (m) at shell height hk (km): the LUT's transmittance toward the sun, and none
// in the planet's shadow (terra's real geometry: its terminator is where it is).
vec3 SunlightAt(vec3 pm, float hk, vec3 sunL)
{
    if (SkyRaySphere(pm / 1000.0, sunL, skyRealHorizon.x / 1000.0) > 0.0) return vec3(0.0);
    return texture(skyTransLut, SkyTransmittanceUv(skyR + hk, dot(normalize(pm), sunL))).rgb;
}
// The glow of the air around the planet's rim seen from outside it (linear, exposed): the view ray's closest pass to
// the centre (its tangent point) sets the air column along the chord (exponential layers: sqrt(2 pi r H) of each at
// that height); single scattering of the sunlight reaching the tangent point, by both phase functions.
vec3 LimbGlow(vec3 cam, vec3 vl, vec3 sunL)
{
    float b = dot(cam, vl);
    if (b >= 0.0) return vec3(0.0);                       // looking away from the planet
    float p = sqrt(max(dot(cam, cam) - b * b, 0.0));      // the closest pass (m)
    float R = skyRealHorizon.x;
    if (p < R || p > skyAirTop) return vec3(0.0);          // the planet's disc (it draws over), or clear of the air
    float hk = ShellHeight(p - R), rs = skyR + hk;
    vec3 tauR = kRayleighScat * exp(-hk / kRayleighH) * sqrt(2.0 * kPi * rs * kRayleighH);
    float colM = exp(-hk / kMieH) * sqrt(2.0 * kPi * rs * kMieH);
    vec3 tau = tauR + vec3(kMieExt * colM);
    float c = dot(vl, sunL);
    vec3 scat = (tauR * SkyPhaseRayleigh(c) + vec3(kMieScat * colM * SkyPhaseMie(c))) / max(tau, vec3(1e-4));
    return SunlightAt(cam + vl * (-b), hk, sunL) * scat * (1.0 - exp(-tau)) * skyPhysExposure; // (unit sun illuminance, as the sky-view LUT)
}

out vec4 finalColor;

float hash3(vec3 p) {
    p = fract(p * 0.3183099 + 0.1);
    p *= 17.0;
    return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
}
float vnoise(vec3 x) {
    vec3 i = floor(x), f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix(hash3(i + vec3(0, 0, 0)), hash3(i + vec3(1, 0, 0)), f.x), mix(hash3(i + vec3(0, 1, 0)), hash3(i + vec3(1, 1, 0)), f.x), f.y),
               mix(mix(hash3(i + vec3(0, 0, 1)), hash3(i + vec3(1, 0, 1)), f.x), mix(hash3(i + vec3(0, 1, 1)), hash3(i + vec3(1, 1, 1)), f.x), f.y), f.z);
}
float fbm(vec3 p) {
    float a = 0.5, s = 0.0;
    for (int i = 0; i < 5; i++) { s += a * vnoise(p); p = p * 2.03 + vec3(1.7, 9.2, 3.1); a *= 0.5; }
    return s;
}

void main()
{
    vec3 v = normalize(fragPosition);
    float s = dot(upDir, sunDir);                       // sun height above the horizon (sine)
    // Sunrise and sunset (#276). The day sky is full from just before sunrise to just after sunset (user: switch much
    // earlier), with a short deep-blue twilight (the blue hour) under it; SkyDayFactor in scene.cpp matches.
    float day = smoothstep(-0.10, 0.02, s);             // ~-5.7 deg .. +1.1 deg
    float twilight = smoothstep(-0.22, -0.08, s) * (1.0 - day);
    float h = dot(v, upDir);                            // view height above the horizon
    float horizon = pow(1.0 - clamp(h, 0.0, 1.0), 4.0);
    float lowSun = 1.0 - smoothstep(0.04, 0.35, s);     // the sun low: a deeper zenith, warm light
    vec3 zenith = mix(vec3(0.18, 0.42, 0.85), vec3(0.09, 0.22, 0.62), lowSun);
    vec3 sky = mix(zenith, mix(vec3(0.62, 0.75, 0.92), vec3(0.72, 0.70, 0.78), lowSun * 0.6), horizon);
    float mu = dot(v, sunDir);
    sky *= day;
    sky += vec3(0.05, 0.09, 0.24) * twilight * (0.6 + 0.4 * horizon);    // the blue hour
    // The twilight glow: strongest with the sun at the horizon (not dimmed by `day`: it lasts past sunset), toward the
    // sun and low in the sky; orange-red at the horizon, gold above it.
    float glow = exp(-(s * s) / (2.0 * 0.07 * 0.07));  // ~4 deg either side of the horizon
    float toward = pow(max(mu, 0.0), 2.5) * (0.25 + 0.75 * horizon);
    vec3 glowCol = mix(vec3(1.0, 0.32, 0.10), vec3(1.0, 0.62, 0.30), smoothstep(0.0, 0.25, h));
    sky += glowCol * glow * toward * 0.9;
    // The pink anti-solar band (the "belt of Venus") a little above the opposite horizon.
    float band = smoothstep(0.0, 0.08, h) * (1.0 - smoothstep(0.10, 0.30, h));
    sky += vec3(0.55, 0.32, 0.45) * glow * pow(max(-mu, 0.0), 1.5) * band * 0.35;
    float sunset = glow;                                // (the clouds' dusk tint below)
    if (skyPhys > 0.5) {
        vec3 phys = pow(max(SkyPhysSky(v, upDir, sunDir, false), vec3(0.0)), vec3(1.0 / 2.2)); // display colours, as the palette
        sky = mix(phys, sky, skyPhysTint);
    }
    sky = mix(sky, vec3(0.52, 0.55, 0.60) * max(day, 0.15 * twilight), overcast * 0.85); // cloud cover greys the sky
    // Glare around the sun, warmer as it sinks.
    sky += mix(vec3(1.0, 0.9, 0.7), vec3(1.0, 0.55, 0.25), lowSun) * pow(max(mu, 0.0), 300.0) * max(day, glow) * density * (1.0 - overcast);
    float alpha = density * (0.06 + 0.94 * max(day, max(twilight * 0.85, glow * toward)));  // night: the stars show through

    float cloudA = 0.0;
    // Clouds: where the view ray meets the cloud shell (from below; from above, its near side).
    if (cloudShell > 0.0 && density > 0.01) {
        vec3 vl = normalize(toLocalX * v.x + toLocalY * v.y + toLocalZ * v.z);
        float b = dot(camLocal, vl), c = dot(camLocal, camLocal) - cloudShell * cloudShell, disc2 = b * b - c;
        float t = disc2 < 0.0 ? -1.0 : (c < 0.0 ? -b + sqrt(disc2) : -b - sqrt(disc2));
        if (t > 0.0) {
            vec3 p = camLocal + vl * t + cloudDrift;
            float cover = clamp(max(cloudCover, 0.22), 0.0, 1.0);
            float thr = 0.66 - 0.32 * cover;                              // a few clouds at the floor, a full deck at 1
            float n = fbm(p / 2600.0) * 0.8 + vnoise(p / 450.0) * 0.2;
            float a = smoothstep(thr, thr + 0.12, n) * exp(-t / 45000.0) * smoothstep(-0.02, 0.10, dot(vl, normalize(camLocal)));
            float lit = mix(0.72, 1.0, smoothstep(thr, thr + 0.3, n));    // denser cores read brighter on top, greyer below
            vec3 cloud = vec3(lit) * mix(vec3(1.0), vec3(0.62, 0.64, 0.68), overcast);
            cloud = mix(cloud, vec3(1.0, 0.62, 0.42), sunset * 0.45);    // dusk tints them
            cloud *= max(day, max(0.04, 0.7 * glow));                    // lit by the glow after sunset
            if (skyPhys > 0.5 && skyExtras.y > 0.5) {
                // #142: lit by the sunlight that reaches the cloud's height (the transmittance LUT: white at noon, golden low,
                // red at the end, none once the planet shadows it) and by the sky above it (the sky-view LUT's zenith).
                vec3 sunL = normalize(toLocalX * sunDir.x + toLocalY * sunDir.y + toLocalZ * sunDir.z);
                vec3 pc = camLocal + vl * t;
                vec3 sunlit = SunlightAt(pc, ShellHeight(length(pc) - skyRealHorizon.x), sunL) * mix(0.55, 1.0, smoothstep(thr, thr + 0.3, n));
                sunlit *= 1.15 + 2.2 * pow(max(mu, 0.0), 6.0); // the droplets scatter forward: clouds toward the sun glow (their silver and gold edges)
                vec3 skylight = SkyPhysSky(upDir, upDir, sunDir, false) * 0.35;
                vec3 phys = pow(max(sunlit + skylight, vec3(0.0)), vec3(1.0 / 2.2)) * mix(vec3(1.0), vec3(0.62, 0.64, 0.68), overcast);
                cloud = mix(phys, cloud, skyPhysTint);
            }
            if (lightning.w > 0.0) {   // (#1162) lit from within: the denser cores most, strongest toward the strike
                float lobe = pow(max(dot(v, lightning.xyz), 0.0), 3.0);
                cloud += vec3(0.82, 0.86, 1.0) * lightning.w * (0.12 + 0.88 * lobe) * lit;
            }
            cloudA = a;
            sky = mix(sky, cloud, a * 0.92);
            alpha = max(alpha, a * 0.9 * density);                       // they hide the stars at night too
        }
    }

    // #142: the limb seen from outside the air (from orbit): the glow of the air around the planet's rim, over the stars.
    if (skyPhys > 0.5 && skyExtras.x > 0.5) {
        float out_ = smoothstep(0.97, 1.03, length(camLocal) / max(skyAirTop, 1.0)); // (inside the air the sky is there)
        if (out_ > 0.0) {
            vec3 vl = normalize(toLocalX * v.x + toLocalY * v.y + toLocalZ * v.z);
            vec3 sunL = normalize(toLocalX * sunDir.x + toLocalY * sunDir.y + toLocalZ * sunDir.z);
            vec3 limb = pow(max(LimbGlow(camLocal, vl, sunL), vec3(0.0)), vec3(1.0 / 2.2)) * out_;
            float la = clamp(max(limb.r, max(limb.g, limb.b)), 0.0, 1.0);
            if (la > 0.001) { // an emissive veil: as bright as it is, it hides the stars behind
                float na = alpha + la * (1.0 - alpha);
                sky = (sky * alpha + limb * (1.0 - alpha)) / max(na, 1e-4);
                alpha = na;
            }
        }
    }
    // (#1162) Lightning lights the air too (the clear sky between clouds, the whole dome a little): an emissive veil like the
    // limb's, so it shows over the night's stars, brightest toward the strike.
    if (lightning.w > 0.0 && density > 0.01) {
        float lobe = pow(max(dot(v, lightning.xyz), 0.0), 3.0);
        float la = clamp(lightning.w * density * (0.15 + 0.55 * lobe), 0.0, 1.0);
        vec3 glowC = vec3(0.72, 0.76, 0.92);
        float na = alpha + la * (1.0 - alpha);
        sky = (sky * alpha + glowC * la) / max(na, 1e-4); // light added over what was there (blended: col x na = sky x alpha + glow x la)
        alpha = na;
    }
    // The sun disc (also in space, where the atmosphere is gone).
    float disc = smoothstep(0.99985, 0.99995, mu) * (1.0 - smoothstep(0.4, 0.8, overcast)) * (1.0 - cloudA);
    vec3 col = mix(sky, mix(vec3(1.0, 0.97, 0.9), vec3(1.0, 0.45, 0.18), 1.0 - smoothstep(0.0, 0.15, s)), disc); // reddened near the horizon
    alpha = max(alpha, disc);
    finalColor = vec4(col, alpha); // already display (sRGB-ish) colours, like the skybox
}
