// Planet terrain on the GPU (Plan 6 step 3, #146): a port of SimplexNoise (3D noise, fractal,
// fractalLimited) and of PlanetField::land / climateAt / biomeShape / elevationBase (voxel.cpp), with
// the same arithmetic in the same order, so the results match the CPU within float rounding
// (--test-terrain-gpu). `precise` keeps the compiler from fusing multiply-adds the CPU does separately.
//
// Uniforms the includer provides: noiseParams[N_COUNT] (frequency, seed x, y, z per field, from
// PlanetField::gpuNoiseParams), planetRadius, planetMoon (0/1).

#define N_CONTINENTS 0
#define N_RIDGES 1
#define N_HILLS 2
#define N_VALLEYS 3
#define N_DETAIL 4
#define N_CAVES_A 5
#define N_CAVES_B 6
#define N_VEINS 7
#define N_CLIMATE 8
#define N_CAVERNS 9
#define N_RAVINES 10
#define N_RAVINE_MASK 11
#define N_STREAMS 12
#define N_ROCKS 13
#define N_ROCK_MASK 14
#define N_RUGGED 15
#define N_CLIFFS 16
#define N_WARP3 17
#define N_TEMP 18
#define N_MOIST 19
#define N_BIOME_JITTER 20
#define N_CAVE_MOUTHS 21
#define N_COUNT 22

uniform vec4 noiseParams[N_COUNT];
uniform float planetRadius;
uniform int planetMoon;

const int kPerm[256] = int[256](
    151, 160, 137, 91, 90, 15,
    131, 13, 201, 95, 96, 53, 194, 233, 7, 225, 140, 36, 103, 30, 69, 142, 8, 99, 37, 240, 21, 10, 23,
    190, 6, 148, 247, 120, 234, 75, 0, 26, 197, 62, 94, 252, 219, 203, 117, 35, 11, 32, 57, 177, 33,
    88, 237, 149, 56, 87, 174, 20, 125, 136, 171, 168, 68, 175, 74, 165, 71, 134, 139, 48, 27, 166,
    77, 146, 158, 231, 83, 111, 229, 122, 60, 211, 133, 230, 220, 105, 92, 41, 55, 46, 245, 40, 244,
    102, 143, 54, 65, 25, 63, 161, 1, 216, 80, 73, 209, 76, 132, 187, 208, 89, 18, 169, 200, 196,
    135, 130, 116, 188, 159, 86, 164, 100, 109, 198, 173, 186, 3, 64, 52, 217, 226, 250, 124, 123,
    5, 202, 38, 147, 118, 126, 255, 82, 85, 212, 207, 206, 59, 227, 47, 16, 58, 17, 182, 189, 28, 42,
    223, 183, 170, 213, 119, 248, 152, 2, 44, 154, 163, 70, 221, 153, 101, 155, 167, 43, 172, 9,
    129, 22, 39, 253, 19, 98, 108, 110, 79, 113, 224, 232, 178, 185, 112, 104, 218, 246, 97, 228,
    251, 34, 242, 193, 238, 210, 144, 12, 191, 179, 162, 241, 81, 51, 145, 235, 249, 14, 239, 107,
    49, 192, 214, 31, 181, 199, 106, 157, 184, 84, 204, 176, 115, 121, 50, 45, 127, 4, 150, 254,
    138, 236, 205, 93, 222, 114, 67, 29, 24, 72, 243, 141, 128, 195, 78, 66, 215, 61, 156, 180);

// Division and square root as IEEE (correctly rounded) like the CPU: GPUs use faster approximations (a few
// ulp), so the quotient / root gets one FMA refinement step.
float xdiv(float a, float b) {
    precise float q = a / b;
    precise float r = fma(-q, b, a);
    precise float res = fma(r, 1.0 / b, q);
    return res;
}
float xsqrt(float x) {
    if (x <= 0.0) return 0.0;
    precise float s = sqrt(x);
    precise float r = fma(-s, s, x);
    precise float res = fma(r, 0.5 / s, s);
    return res;
}

int nHash(int i) { return kPerm[i & 255]; }
int nFloor(float f) { int i = int(f); return (f < float(i)) ? i - 1 : i; }
float nGrad(int h, float x, float y, float z) {
    h &= 15;
    float u = h < 8 ? x : y;
    float v = h < 4 ? y : (h == 12 || h == 14) ? x : z;
    return (((h & 1) != 0) ? -u : u) + (((h & 2) != 0) ? -v : v);
}

float snoise(float x, float y, float z) {
    const float F3 = 1.0 / 3.0;
    const float G3 = 1.0 / 6.0;
    precise float s = (x + y + z) * F3;
    int i = nFloor(x + s), j = nFloor(y + s), k = nFloor(z + s);
    precise float t = float(i + j + k) * G3;
    precise float X0 = float(i) - t, Y0 = float(j) - t, Z0 = float(k) - t;
    precise float x0 = x - X0, y0 = y - Y0, z0 = z - Z0;
    int i1, j1, k1, i2, j2, k2;
    if (x0 >= y0) {
        if (y0 >= z0)      { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 1; k2 = 0; }
        else if (x0 >= z0) { i1 = 1; j1 = 0; k1 = 0; i2 = 1; j2 = 0; k2 = 1; }
        else               { i1 = 0; j1 = 0; k1 = 1; i2 = 1; j2 = 0; k2 = 1; }
    } else {
        if (y0 < z0)       { i1 = 0; j1 = 0; k1 = 1; i2 = 0; j2 = 1; k2 = 1; }
        else if (x0 < z0)  { i1 = 0; j1 = 1; k1 = 0; i2 = 0; j2 = 1; k2 = 1; }
        else               { i1 = 0; j1 = 1; k1 = 0; i2 = 1; j2 = 1; k2 = 0; }
    }
    precise float x1 = x0 - float(i1) + G3, y1 = y0 - float(j1) + G3, z1 = z0 - float(k1) + G3;
    precise float x2 = x0 - float(i2) + 2.0 * G3, y2 = y0 - float(j2) + 2.0 * G3, z2 = z0 - float(k2) + 2.0 * G3;
    precise float x3 = x0 - 1.0 + 3.0 * G3, y3 = y0 - 1.0 + 3.0 * G3, z3 = z0 - 1.0 + 3.0 * G3;
    int gi0 = nHash(i + nHash(j + nHash(k)));
    int gi1 = nHash(i + i1 + nHash(j + j1 + nHash(k + k1)));
    int gi2 = nHash(i + i2 + nHash(j + j2 + nHash(k + k2)));
    int gi3 = nHash(i + 1 + nHash(j + 1 + nHash(k + 1)));
    precise float n0 = 0.0, n1 = 0.0, n2 = 0.0, n3 = 0.0;
    precise float t0 = 0.6 - x0 * x0 - y0 * y0 - z0 * z0;
    if (t0 >= 0.0) { t0 *= t0; n0 = t0 * t0 * nGrad(gi0, x0, y0, z0); }
    precise float t1 = 0.6 - x1 * x1 - y1 * y1 - z1 * z1;
    if (t1 >= 0.0) { t1 *= t1; n1 = t1 * t1 * nGrad(gi1, x1, y1, z1); }
    precise float t2 = 0.6 - x2 * x2 - y2 * y2 - z2 * z2;
    if (t2 >= 0.0) { t2 *= t2; n2 = t2 * t2 * nGrad(gi2, x2, y2, z2); }
    precise float t3 = 0.6 - x3 * x3 - y3 * y3 - z3 * z3;
    if (t3 >= 0.0) { t3 *= t3; n3 = t3 * t3 * nGrad(gi3, x3, y3, z3); }
    precise float r = 32.0 * (n0 + n1 + n2 + n3);
    return r;
}

float fractal(int f, int octaves, float x, float y, float z) {
    vec4 p = noiseParams[f];
    precise float output_ = 0.0, denom = 0.0;
    precise float frequency = p.x, amplitude = 1.0;
    for (int o = 0; o < octaves; ++o) {
        precise float xs = (x + p.y) * frequency, ys = (y + p.z) * frequency, zs = (z + p.w) * frequency;
        output_ += amplitude * snoise(xs, ys, zs);
        denom += amplitude;
        frequency *= 2.0;
        amplitude *= 0.5;
    }
    precise float r = xdiv(output_, denom);
    return r;
}

float fractalLimited(int f, int octaves, float footprint, float x, float y, float z) {
    vec4 p = noiseParams[f];
    precise float output_ = 0.0, denom = 0.0;
    precise float frequency = p.x, amplitude = 1.0;
    for (int o = 0; o < octaves; ++o) {
        denom += amplitude;
        precise float w = footprint > 0.0 ? xdiv(1.0, frequency * 2.0 * footprint) - 1.0 : 1.0;
        w = w >= 1.0 ? 1.0 : (w <= 0.0 ? 0.0 : w);
        if (w > 0.0) {
            precise float xs = (x + p.y) * frequency, ys = (y + p.z) * frequency, zs = (z + p.w) * frequency;
            output_ += w * amplitude * snoise(xs, ys, zs);
        }
        frequency *= 2.0;
        amplitude *= 0.5;
    }
    precise float r = xdiv(output_, denom);
    return r;
}

float tSmooth(float a, float b, float x) {
    precise float t = clamp(xdiv(x - a, b - a), 0.0, 1.0);
    precise float r = t * t * (3.0 - 2.0 * t);
    return r;
}

// PlanetField::Land: h, c, mountain, valley, valleyWeight, detail, stream, ridge.
struct Land { float h, c, mountain, valley, valleyWeight, detail, stream, ridge; };

Land landAt(vec3 dir, float fp) {
    precise float sx = dir.x * planetRadius, sy = dir.y * planetRadius, sz = dir.z * planetRadius;
    Land l = Land(0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0);
    if (planetMoon != 0) {
        precise float base = fractalLimited(N_CONTINENTS, 4, fp, sx, sy, sz) * 60.0;
        precise float r = 1.0 - abs(fractalLimited(N_RIDGES, 4, fp, sx, sy, sz));
        precise float hm = base + r * r * 250.0 + fractalLimited(N_HILLS, 3, fp, sx, sy, sz) * 25.0 + fractalLimited(N_DETAIL, 2, fp, sx, sy, sz) * 2.5;
        l.h = hm;
        return l;
    }
    precise float c = fractalLimited(N_CONTINENTS, 4, fp, sx, sy, sz);
    precise float base = c < 0.0 ? c * 650.0 : c * 200.0;
    precise float landMask = tSmooth(-0.02, 0.06, c);
    precise float mountainMask = tSmooth(0.08, 0.35, c);
    precise float h = base;
    if (mountainMask > 0.0) {
        precise float r = 1.0 - abs(fractalLimited(N_RIDGES, 5, fp, sx, sy, sz));
        l.ridge = r * r * r * 1400.0 * mountainMask;
        h += l.ridge;
    }
    h += fractalLimited(N_HILLS, 3, fp, sx, sy, sz) * (25.0 + 60.0 * mountainMask);
    precise float valleyWeight = landMask * clamp(2.0 - xdiv(fp, 64.0), 0.0, 1.0);
    l.valley = 1.0;
    precise float ravineWeight = landMask * clamp(2.0 - xdiv(fp, 24.0), 0.0, 1.0);
    if (ravineWeight > 0.0) {
        precise float region = tSmooth(0.25, 0.4, fractal(N_RAVINE_MASK, 2, sx, sy, sz));
        if (region > 0.0) {
            precise float r = abs(fractal(N_RAVINES, 2, sx, sy, sz));
            if (r < 0.006) h -= xsqrt(1.0 - xdiv(r, 0.006)) * 90.0 * region * ravineWeight;
        }
    }
    l.detail = fractalLimited(N_DETAIL, 2, fp, sx, sy, sz) * 2.5;
    precise float iceCap = tSmooth(0.86, 0.92, abs(dir.y) + 0.03 * fractalLimited(N_CLIMATE, 2, fp, sx, sy, sz));
    precise float hd = h + l.detail;
    if (iceCap > 0.0 && hd < 2.0) { hd += (2.0 - hd) * iceCap; h = hd - l.detail; }
    l.h = h;
    l.c = c;
    l.mountain = mountainMask;
    l.valleyWeight = valleyWeight * (1.0 - iceCap);
    return l;
}

// PlanetField::climateAt: (temperature C, moisture 0..1).
vec2 climateAt(vec3 dir, float elev) {
    if (planetMoon != 0) return vec2(-20.0, 0.0);
    precise float lat = abs(dir.y);
    precise float x = dir.x * planetRadius, y = dir.y * planetRadius, z = dir.z * planetRadius;
    precise float temp = 30.0 - 44.0 * (lat * xsqrt(lat)) - xdiv(6.5 * max(elev, 0.0), 1000.0) + 7.0 * fractal(N_TEMP, 2, x, y, z);
    precise float moist = 0.5 + 0.85 * fractal(N_MOIST, 3, x, y, z) + 0.08 * (1.0 - tSmooth(10.0, 150.0, elev));
    return vec2(temp, clamp(moist, 0.0, 1.0));
}

// PlanetField::biomeShape: dunes, mesas, swamp floors.
float biomeShape(vec3 dir, float h, float fp) {
    if (planetMoon != 0 || h < 2.0) return 0.0;
    vec2 cl = climateAt(dir, h);
    precise float hot = tSmooth(16.0, 20.0, cl.x);
    precise float desert = hot * tSmooth(0.27, 0.2, cl.y);
    precise float badlands = hot * tSmooth(0.2, 0.27, cl.y) * (1.0 - tSmooth(0.36, 0.42, cl.y));
    precise float swamp = tSmooth(0.7, 0.82, cl.y) * tSmooth(9.0, 15.0, cl.x) * (1.0 - tSmooth(60.0, 110.0, h));
    if (desert + badlands + swamp < 1e-3) return 0.0;
    precise float x = dir.x * planetRadius, y = dir.y * planetRadius, z = dir.z * planetRadius;
    precise float dh = 0.0;
    if (desert > 0.0 && fp < 48.0) {
        precise float warp = 2.2 * fractal(N_BIOME_JITTER, 2, x * 0.25, y * 0.25, z * 0.25);
        precise float u = xdiv(x * 0.8 + z * 0.6, 55.0) + warp;
        precise float fu = u - floor(u), tri = 1.0 - abs(2.0 * fu - 1.0), sm = tri * tri * (3.0 - 2.0 * tri);
        precise float crest = sm * sm;
        precise float swell = 0.55 + 0.45 * fractal(N_MOIST, 1, x * 6.0, y * 6.0, z * 6.0);
        dh += desert * crest * 9.0 * swell * tSmooth(4.0, 20.0, h) * clamp(2.0 - xdiv(fp, 24.0), 0.0, 1.0);
    }
    if (badlands > 0.0) {
        const float step_ = 14.0;
        precise float t = xdiv(h, step_), f = t - floor(t);
        precise float terraced = (floor(t) + tSmooth(0.72, 0.92, f)) * step_;
        dh += (terraced - h) * badlands;
    }
    if (swamp > 0.0) {
        precise float floorH = 3.5 + 1.5 * fractal(N_BIOME_JITTER, 1, x, y, z);
        if (h > floorH) dh += (floorH - h) * 0.85 * swamp;
    }
    return dh;
}

// PlanetField::elevationBase: the height map without rivers, lakes and ponds.
float elevationBase(vec3 dir, float fp) {
    Land l = landAt(dir, fp);
    if (planetMoon != 0) return l.h;
    precise float h = l.h + l.detail;
    precise float r = h + biomeShape(dir, h, fp);
    return r;
}

// ---------------------------------------------------------------------------------------------
// Density near the surface (Plan 6 step 4, #147): PlanetField::densityNear and its parts (terrain3D,
// rockFeature, caves, floatingIsland) from voxel.cpp. More uniforms: mantleRadius, islandCount; and the
// includer defines islandData(i, k): vec4 k of island i: (cell x, y, z, ok), (dir, waist), (r, depth, -, -).

#define ISO_LEVEL 0.5
#define DENSITY_GRADIENT 0.25
#define K_TUNNEL_LEVEL 0.096
#define K_CAVERN_LEVEL 0.68
#define K_BAND_BELOW 16.0
#define K_BAND_ABOVE 50.0
#define K_BOULDER_CELL 24.0
#define K_BOULDER_MAX 5.0
#define K_ISLAND_CELL 4000.0

uniform float mantleRadius;
int islandCount = 0;   // set per block by the includer

float waveSin(float x) {
    precise float u = x * 0.15915494, f = u - floor(u);
    precise float r = f < 0.5 ? 16.0 * f * (0.5 - f) : -16.0 * (f - 0.5) * (1.0 - f);
    return r;
}
float vlen(vec3 v) { precise float s = v.x * v.x + v.y * v.y + v.z * v.z; return xsqrt(s); }
float vdot(vec3 a, vec3 b) { precise float s = a.x * b.x + a.y * b.y + a.z * b.z; return s; }

// PlanetField::terrain3D with the column given (mask, f).
float terrain3D(vec3 p, float d, float surface, float fp, float mask, float f) {
    if (planetMoon != 0 || d < 1e-3) return 0.0;
    precise float h = d - surface;
    if (h < -K_BAND_BELOW || h > K_BAND_ABOVE) return 0.0;
    vec3 dir = vec3(xdiv(p.x, d), xdiv(p.y, d), xdiv(p.z, d));
    if (mask <= 0.0) return 0.0;
    precise float amp = 16.0 * (0.5 + 1.5 * f), top = 20.0 + 26.0 * f;
    precise float hs = h * 0.6;
    precise vec3 q = vec3(p.x - dir.x * hs, p.y - dir.y * hs, p.z - dir.z * hs);
    precise float w = 9.0 * fractalLimited(N_WARP3, 2, fp, q.x, q.y, q.z);
    q = vec3(q.x + w, q.y - 0.7 * w, q.z + 0.4 * w);
    precise float n = fractalLimited(N_CLIFFS, 3, fp, q.x, q.y, q.z);
    precise float shape = n + 0.3 * (1.0 - 2.0 * abs(n)) + 0.15;
    precise float fade = (1.0 - tSmooth(3.0, top, h)) * tSmooth(-16.0, -8.0, h);
    precise float r = max(amp * mask * shape * fade, -3.0);
    return r;
}

uint pondHash(int x, int y, int z, uint seed) {
    uint h = uint(x) * 73856093u ^ uint(y) * 19349663u ^ uint(z) * 83492791u ^ seed * 2654435761u;
    h ^= h >> 16; h *= 0x7feb352du; h ^= h >> 15; h *= 0x846ca68bu; h ^= h >> 16;
    return h;
}
float pondUnit(inout uint h) { h = h * 1664525u + 1013904223u; precise float r = float(h >> 8) * (1.0 / 16777216.0); return r; }

// PlanetField::rockFeature: metres inside the nearest tor or boulder (< 0 outside).
float rockFeature(vec3 p, float d, float surface, float fp) {
    if (planetMoon != 0 || d < 1e-3) return -1e9;
    precise float w = clamp(2.0 - xdiv(fp, 4.0), 0.0, 1.0);
    precise float h = d - surface;
    if (w <= 0.0 || h > 14.0 || h < -3.0) return -1e9; // tops reach 12 m (#242)
    precise float elev = surface - planetRadius;
    if (elev < 4.0) return -1e9;
    precise float best = -1e9;
    precise float mask = tSmooth(0.5, 0.75, fractal(N_ROCK_MASK, 2, p.x, p.y, p.z) + 0.2 * tSmooth(150.0, 600.0, elev));
    vec3 dir = vec3(xdiv(p.x, d), xdiv(p.y, d), xdiv(p.z, d));
    precise vec3 s = vec3(dir.x * planetRadius, dir.y * planetRadius, dir.z * planetRadius);
    if (mask > 0.0) {
        precise float n = fractal(N_ROCKS, 2, s.x, s.y, s.z);
        precise float top = (3.0 + 9.0 * mask) * (0.6 + 0.4 * snoise(s.x * 0.01, s.y * 0.01, s.z * 0.01));
        precise float ff = clamp(xdiv(h, max(top, 0.5)), 0.0, 1.0);
        precise float thr = 0.3 - 0.06 * ff * ff + 0.04 * waveSin(h * 2.1 + n * 9.0) * (1.0 - tSmooth(0.35, 0.7, ff));
        precise float side = (n - thr) * 14.0;
        best = max(best, min(side, top - h) * w);
    }
    precise float reach = K_BOULDER_MAX * 1.3;
    int x0 = int(floor(xdiv(s.x - reach, K_BOULDER_CELL))), x1 = int(floor(xdiv(s.x + reach, K_BOULDER_CELL)));
    int y0 = int(floor(xdiv(s.y - reach, K_BOULDER_CELL))), y1 = int(floor(xdiv(s.y + reach, K_BOULDER_CELL)));
    int z0 = int(floor(xdiv(s.z - reach, K_BOULDER_CELL))), z1 = int(floor(xdiv(s.z + reach, K_BOULDER_CELL)));
    for (int cz = z0; cz <= z1; ++cz)
        for (int cy = y0; cy <= y1; ++cy)
            for (int cx = x0; cx <= x1; ++cx) {
                uint hh = pondHash(cx, cy, cz, uint(planetRadius * 3.0) ^ 0xb01du);
                if (!(pondUnit(hh) < 0.16)) continue;
                precise float ux = pondUnit(hh);
                precise float uy = pondUnit(hh);
                precise float uz = pondUnit(hh);
                precise vec3 c = vec3((float(cx) + ux) * K_BOULDER_CELL, (float(cy) + uy) * K_BOULDER_CELL, (float(cz) + uz) * K_BOULDER_CELL);
                precise float len = vlen(c);
                if (!(len > 1.0)) continue;
                precise float il = xdiv(1.0, len);
                precise vec3 bdir = vec3(c.x * il, c.y * il, c.z * il);
                precise vec3 q = vec3(bdir.x * planetRadius, bdir.y * planetRadius, bdir.z * planetRadius);
                if (!(int(floor(xdiv(q.x, K_BOULDER_CELL))) == cx && int(floor(xdiv(q.y, K_BOULDER_CELL))) == cy && int(floor(xdiv(q.z, K_BOULDER_CELL))) == cz)) continue;
                precise float u = pondUnit(hh);
                precise float br = 1.2 + (K_BOULDER_MAX - 1.2) * u * u;
                precise float lift = br * (0.1 + 0.3 * pondUnit(hh));
                precise vec3 dd = vec3(bdir.x - dir.x, bdir.y - dir.y, bdir.z - dir.z);
                precise float horiz = vlen(dd) * planetRadius;
                if (horiz > br * 1.3) continue;
                precise float vert = xdiv(h - lift, 0.75);
                precise float lump = 0.25 * br * fractal(N_ROCKS, 2, p.x * 1.7, p.y * 1.7, p.z * 1.7);
                precise float rr = horiz * horiz + vert * vert;
                best = max(best, (br + lump - xsqrt(rr)) * w);
            }
    return best;
}

float islandTop(vec3 p, float q) {
    precise float hills = fractal(N_DETAIL, 2, p.x * 1.4, p.y * 1.4, p.z * 1.4);
    precise float t = 4.0 + 10.0 * max(hills, -0.2) * (1.0 - q * q);
    const float step_ = 5.0;
    precise float terr = xdiv(t, step_), fr = terr - floor(terr);
    precise float terraced = (floor(terr) + tSmooth(0.65, 0.9, fr)) * step_;
    precise float cliffMask = tSmooth(0.0, 0.3, fractal(N_RUGGED, 1, p.x * 8.0, p.y * 8.0, p.z * 8.0));
    t += (terraced - t) * cliffMask;
    precise float r = t + 2.0 * (1.0 - q);
    return r;
}

vec4 islandData(int i, int k);

// PlanetField::floatingIsland: metres inside an island (< 0 outside).
float floatingIsland(vec3 p, float d) {
    precise float sc = xdiv(planetRadius, d);
    precise vec3 s = vec3(p.x * sc, p.y * sc, p.z * sc);
    ivec3 cell = ivec3(int(floor(xdiv(s.x, K_ISLAND_CELL))), int(floor(xdiv(s.y, K_ISLAND_CELL))), int(floor(xdiv(s.z, K_ISLAND_CELL))));
    int found = -1;
    for (int i = 0; i < islandCount; ++i) {
        vec4 c = islandData(i, 0);
        if (ivec3(c.xyz) == cell) { found = i; break; }
    }
    if (found < 0) return -1e9;
    vec4 dw = islandData(found, 1), rd = islandData(found, 2);
    vec3 idir = dw.xyz;
    float waist = dw.w, ir = rd.x, idepth = rd.y;
    precise float pd = vdot(p, idir);
    precise float along = pd - waist;
    if (along > 30.0 || along < -idepth - 6.0) return -1e9;
    precise vec3 lat3 = vec3(p.x - idir.x * pd, p.y - idir.y * pd, p.z - idir.z * pd);
    precise float lateral = vlen(lat3);
    if (lateral > ir * 1.5) return -1e9;
    precise float wobble = 0.35 * fractal(N_DETAIL, 2, p.x * 0.8, p.y * 0.8, p.z * 0.8) + 0.1 * fractal(N_DETAIL, 2, p.x * 2.5, p.y * 2.5, p.z * 2.5);
    precise float lump = fractal(N_DETAIL, 2, p.x * 2.0 + 31.0, p.y * 2.0, p.z * 2.0);
    precise float rr = ir * (1.0 + wobble), q = clamp(xdiv(lateral, rr), 0.0, 1.0);
    precise float top = islandTop(p, q);
    precise float oq = 1.0 - q;
    precise float bottom = -idepth * (0.75 + 0.25 * lump) * (oq * xsqrt(xsqrt(oq))) - 2.0;
    return min(min(top - along, along - bottom), rr - lateral);
}

// PlanetField::densityNear with the column given.
float densityNear(vec3 p, float d, float surface, float fp, float cmask, float cf) {
    precise float v = ISO_LEVEL + DENSITY_GRADIENT * (surface - d);
    precise float depth = surface - d;
    if (depth < K_BAND_BELOW && depth > -K_BAND_ABOVE) v += DENSITY_GRADIENT * terrain3D(p, d, surface, fp, cmask, cf);
    if (abs(depth) < 6.0) v += DENSITY_GRADIENT * 1.5 * fractalLimited(N_DETAIL, 2, fp, p.x, p.y * 1.3, p.z);
    if (depth < 3.0 && depth > -14.0) {
        precise float m = rockFeature(p, d, surface, fp);
        if (m > -4.0) v = max(v, ISO_LEVEL + DENSITY_GRADIENT * m);
    }
    precise float caveWeight = clamp(2.0 - xdiv(fp, 8.0), 0.0, 1.0);
    precise float caveTop = 5.0;
    if (depth < 14.0 && depth > -3.0 && caveWeight > 0.0 && planetMoon == 0) caveTop = 5.0 - 7.5 * tSmooth(0.3, 0.45, fractal(N_CAVE_MOUTHS, 1, p.x, p.y, p.z));
    if (caveWeight > 0.0 && depth > caveTop && d > mantleRadius + 20.0) {
        precise float near_ = 2.0 * max(fp, 1.0) + 2.0;
        precise float m = 1e9;
        precise float a = abs(fractal(N_CAVES_A, 2, p.x, p.y, p.z));
        if ((a - K_TUNNEL_LEVEL) * 100.0 <= near_) {
            precise float fade = min(1.0, xdiv(depth - caveTop, 8.0)) * min(1.0, xdiv(d - mantleRadius - 20.0, 10.0));
            precise float t = max(a, abs(fractal(N_CAVES_B, 2, p.x, p.y, p.z)));
            m = (t - K_TUNNEL_LEVEL) * 100.0 + (1.0 - fade * caveWeight) * 20.0;
        }
        if (depth > 30.0 && d > mantleRadius + 40.0) {
            precise float fade = min(1.0, xdiv(depth - 30.0, 20.0)) * min(1.0, xdiv(d - mantleRadius - 40.0, 20.0));
            precise float c = fractal(N_CAVERNS, 2, p.x, p.y * 1.4, p.z);
            m = min(m, (K_CAVERN_LEVEL - c) * 120.0 + (1.0 - fade * caveWeight) * 30.0);
        }
        if (m < near_) v = min(v, ISO_LEVEL + DENSITY_GRADIENT * m);
    }
    if (depth < -15.0 && planetMoon == 0) {
        precise float m = floatingIsland(p, d);
        if (m > -4.0) v = max(v, ISO_LEVEL + DENSITY_GRADIENT * m);
    }
    return v;
}
