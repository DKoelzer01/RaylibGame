// #142: physical sky, Hillaire 2020 ("A Scalable and Production Ready Sky and Atmosphere Rendering Technique") look-up
// tables. Shared by the LUT passes (skyphys.cpp prepends "#version 430" and appends one pass's main) and atmosphere.fs
// (the sky-view lookup). Planet-local frame with the centre at the origin; lengths in km (float precision at Earth size).
// Earth's air (Rayleigh, Mie, ozone) over an Earth-sized shell (skyphys.cpp): terra's 60 km radius and 25 km of air give a
// horizon path ~15x shorter than Earth's, too short for a red sunset; the camera height and view angles are mapped onto it.

uniform float skyR;    // the LUT shell's sea level (km)
uniform float skyRt;   // its top of the air (km)

const float kPi = 3.14159265;
const vec3 kRayleighScat = vec3(5.802, 13.558, 33.1) * 1e-3;
const float kRayleighH = 8.0;
const float kMieScat = 3.996e-3, kMieExt = 4.440e-3, kMieH = 1.2, kMieG = 0.8;
const vec3 kOzoneAbs = vec3(0.650, 1.881, 0.085) * 1e-3;

struct SkyMedium { vec3 rayleigh; float mie; vec3 scat; vec3 ext; };

SkyMedium SkyMediumAt(float h) // h: height above sea level
{
    h = max(h, 0.0);
    SkyMedium m;
    m.rayleigh = kRayleighScat * exp(-h / kRayleighH);
    float dm = exp(-h / kMieH);
    m.mie = kMieScat * dm;
    float oz = max(0.0, 1.0 - abs(h - 25.0) / 15.0); // ozone: Earth's tent (25 km +- 15 km), the blue zenith at dusk
    m.scat = m.rayleigh + vec3(m.mie);
    m.ext = m.rayleigh + vec3(kMieExt * dm) + kOzoneAbs * oz;
    return m;
}

// The nearest hit at t >= 0 of the ray with the sphere of radius r about the origin, else -1.
float SkyRaySphere(vec3 ro, vec3 rd, float r)
{
    float b = dot(ro, rd), c = dot(ro, ro) - r * r, d = b * b - c;
    if (d < 0.0) return -1.0;
    float s = sqrt(d), t0 = -b - s, t1 = -b + s;
    if (t0 >= 0.0) return t0;
    if (t1 >= 0.0) return t1;
    return -1.0;
}

float SkyPhaseRayleigh(float c) { return 3.0 / (16.0 * kPi) * (1.0 + c * c); }
float SkyPhaseMie(float c) // Cornette-Shanks
{
    float g = kMieG, k = 3.0 / (8.0 * kPi) * (1.0 - g * g) / (2.0 + g * g);
    return k * (1.0 + c * c) / pow(max(1.0 + g * g - 2.0 * g * c, 1e-4), 1.5);
}

// Transmittance LUT parameterisation (Bruneton 2017): height and the cosine of the view zenith <-> uv.
vec2 SkyTransmittanceUv(float vh, float mu)
{
    float H = sqrt(max(0.0, skyRt * skyRt - skyR * skyR));
    float rho = sqrt(max(0.0, vh * vh - skyR * skyR));
    float disc = vh * vh * (mu * mu - 1.0) + skyRt * skyRt;
    float d = max(0.0, -vh * mu + sqrt(max(disc, 0.0)));
    float dMin = skyRt - vh, dMax = rho + H;
    return vec2((d - dMin) / max(dMax - dMin, 1e-3), rho / H);
}
void SkyTransmittanceParams(vec2 uv, out float vh, out float mu)
{
    float H = sqrt(max(0.0, skyRt * skyRt - skyR * skyR));
    float rho = H * uv.y;
    vh = sqrt(rho * rho + skyR * skyR);
    float dMin = skyRt - vh, dMax = rho + H, d = dMin + uv.x * (dMax - dMin);
    mu = d == 0.0 ? 1.0 : (H * H - rho * rho - d * d) / (2.0 * vh * d);
    mu = clamp(mu, -1.0, 1.0);
}

// Sky-view LUT parameterisation (Hillaire 2020, 5.3): u = the azimuth from the sun (squared toward it), v = the view
// zenith, packed around the horizon (half above, half below).
vec2 SkyViewUv(float vh, float viewCosZenith, float lightViewCos)
{
    float vHorizon = sqrt(max(vh * vh - skyR * skyR, 0.0));
    float beta = acos(clamp(vHorizon / vh, -1.0, 1.0));
    float zenithHorizon = kPi - beta;
    float a = acos(clamp(viewCosZenith, -1.0, 1.0));
    float v;
    if (a < zenithHorizon) { float c = a / zenithHorizon; c = 1.0 - sqrt(max(1.0 - c, 0.0)); v = c * 0.5; }
    else { float c = (a - zenithHorizon) / max(beta, 1e-4); v = sqrt(clamp(c, 0.0, 1.0)) * 0.5 + 0.5; }
    float u = sqrt(clamp(-lightViewCos * 0.5 + 0.5, 0.0, 1.0));
    return vec2(u, v);
}

// ---- The sky-view lookup (atmosphere.fs for the sky, skyphys_haze.glsl for the haze on lit surfaces) ----
uniform sampler2D skyViewLut;
uniform float skyViewHeight;   // the camera's height in the shell (km, SkyPhys_LutHeight)
uniform vec2 skyRealHorizon;   // terra's radius and the camera's distance from its centre (m): the real horizon's dip
uniform float skyPhysExposure; // sun illuminance x exposure into display
uniform float skyPhysTint;     // 0 = the physical sky as it is .. 1 = today's palette

// The in-scattered light toward v (unit; up = local up at the camera, sun = toward the sun), linear, exposed. The shell
// is Earth-sized, so its horizon dips less than terra's: the view zenith angle is remapped piecewise from the real horizon
// to the shell's (zenith and nadir stay put). aboveHorizon: hold v just over the real horizon (haze on the ground).
vec3 SkyPhysSky(vec3 v, vec3 up, vec3 sun, bool aboveHorizon)
{
    float h = dot(v, up), s = dot(sun, up);
    vec3 vt = v - up * h, st = sun - up * s;
    float lt = length(vt), ls = length(st);
    float lightViewCos = (lt > 1e-4 && ls > 1e-4) ? dot(vt, st) / (lt * ls) : 1.0;
    float vh = clamp(skyViewHeight, skyR + 0.001, skyRt - 0.001);
    float zr = kPi - asin(clamp(skyRealHorizon.x / skyRealHorizon.y, 0.0, 1.0));
    float zv = kPi - asin(clamp(skyR / vh, 0.0, 1.0));
    float za = acos(clamp(h, -1.0, 1.0));
    if (aboveHorizon) za = min(za, zr - 0.02);
    za = za < zr ? za * zv / zr : zv + (za - zr) * (kPi - zv) / max(kPi - zr, 1e-4);
    // Brightness held near today's as the eye adapts (captures, #142): 0.7x from 17 deg up, a 5 deg sun 1x, dusk up to 3.5x.
    float adapt = s < 0.087 ? mix(3.5, 1.0, smoothstep(-0.05, 0.087, s)) : mix(1.0, 0.7, smoothstep(0.087, 0.3, s));
    return texture(skyViewLut, SkyViewUv(vh, cos(za), lightViewCos)).rgb * skyPhysExposure * adapt;
}
