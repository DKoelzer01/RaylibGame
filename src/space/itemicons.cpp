#include "itemicons.h"
#include "inventory.h"
#include "buildtool.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <functional>
#include <thread>
#include <unordered_map>
#include <vector>

// Icons are painted into a float canvas from signed distance functions: each shape is anti-aliased
// (coverage from the distance), shaded (light from the top left), and edged with a darker outline of
// its own colour. Painters are chosen by the item's id (suffixes like _ore, _log, _pickaxe, or the
// exact id), coloured by items.txt.

namespace {

constexpr int kIcon = 64;     // pixels per icon
constexpr int kCols = 16;     // icons per atlas row

struct C3 { float r, g, b; };
C3 operator*(C3 a, float k) { return { a.r * k, a.g * k, a.b * k }; }
C3 operator+(C3 a, C3 b) { return { a.r + b.r, a.g + b.g, a.b + b.b }; }
C3 Mix(C3 a, C3 b, float t) { return { a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t }; }
C3 FromColor(Color c) { return { c.r / 255.0f, c.g / 255.0f, c.b / 255.0f }; }
const C3 kWhite{ 1, 1, 1 }, kWood{ 0.55f, 0.36f, 0.2f }, kSteel{ 0.62f, 0.64f, 0.68f }, kDark{ 0.16f, 0.17f, 0.19f };
const C3 kGlow{ 1.0f, 0.62f, 0.15f }, kWater{ 0.25f, 0.5f, 0.9f }, kLeaf{ 0.3f, 0.58f, 0.2f }, kGold{ 0.9f, 0.72f, 0.25f };

float Clamp01(float v) { return v < 0 ? 0 : v > 1 ? 1 : v; }
float Hash(int x, int y, int seed) {
    unsigned h = (unsigned)x * 374761393u + (unsigned)y * 668265263u + (unsigned)seed * 2246822519u;
    h = (h ^ (h >> 13)) * 1274126177u;
    return ((h ^ (h >> 16)) & 0xFFFF) / 65535.0f;
}

using Sdf = std::function<float(float, float)>;
using Shade = std::function<C3(float, float)>;
struct V2 { float x, y; };

Sdf Circle(float cx, float cy, float r) { return [=](float x, float y) { return hypotf(x - cx, y - cy) - r; }; }
Sdf Ellipse(float cx, float cy, float rx, float ry) {
    return [=](float x, float y) { return (hypotf((x - cx) / rx, (y - cy) / ry) - 1.0f) * std::min(rx, ry); };
}
Sdf Box(float cx, float cy, float hx, float hy, float r = 0, float ang = 0) {
    const float cs = cosf(ang), sn = sinf(ang);
    return [=](float x, float y) {
        const float dx = x - cx, dy = y - cy, lx = cs * dx + sn * dy, ly = -sn * dx + cs * dy;
        const float qx = fabsf(lx) - hx + r, qy = fabsf(ly) - hy + r;
        return hypotf(std::max(qx, 0.0f), std::max(qy, 0.0f)) + std::min(std::max(qx, qy), 0.0f) - r;
    };
}
Sdf Capsule(float ax, float ay, float bx, float by, float r) {
    return [=](float x, float y) {
        const float px = x - ax, py = y - ay, ex = bx - ax, ey = by - ay;
        const float t = Clamp01((px * ex + py * ey) / std::max(ex * ex + ey * ey, 1e-6f));
        return hypotf(px - ex * t, py - ey * t) - r;
    };
}
Sdf Poly(std::vector<V2> v) {
    return [v](float x, float y) {
        float d = (x - v[0].x) * (x - v[0].x) + (y - v[0].y) * (y - v[0].y), s = 1.0f;
        for (size_t i = 0, j = v.size() - 1; i < v.size(); j = i++) {
            const float ex = v[j].x - v[i].x, ey = v[j].y - v[i].y, wx = x - v[i].x, wy = y - v[i].y;
            const float t = Clamp01((wx * ex + wy * ey) / std::max(ex * ex + ey * ey, 1e-6f));
            const float bx = wx - ex * t, by = wy - ey * t;
            d = std::min(d, bx * bx + by * by);
            const bool c1 = y >= v[i].y, c2 = y < v[j].y, c3 = ex * wy > ey * wx;
            if ((c1 && c2 && c3) || (!c1 && !c2 && !c3)) s = -s;
        }
        return s * sqrtf(d);
    };
}
Sdf Ring(float cx, float cy, float r, float w) { return [=](float x, float y) { return fabsf(hypotf(x - cx, y - cy) - r) - w; }; }
Sdf Union(Sdf a, Sdf b) { return [=](float x, float y) { return std::min(a(x, y), b(x, y)); }; }
Sdf Cut(Sdf a, Sdf b) { return [=](float x, float y) { return std::max(a(x, y), -b(x, y)); }; }
Sdf Inter(Sdf a, Sdf b) { return [=](float x, float y) { return std::max(a(x, y), b(x, y)); }; }
// The shape turned by ang (radians) about (cx, cy).
Sdf Rot(Sdf a, float cx, float cy, float ang) {
    const float cs = cosf(-ang), sn = sinf(-ang);
    return [=](float x, float y) { const float dx = x - cx, dy = y - cy; return a(cx + cs * dx - sn * dy, cy + sn * dx + cs * dy); };
}

// Light from the top left, over the icon.
Shade Lin(C3 c, float k = 0.32f) { return [=](float x, float y) { return c * (1.16f - k * (x + y) / kIcon); }; }
// A rounded surface centred at (cx, cy) with radius r: brighter toward the top left, a soft rim.
Shade Round(C3 c, float cx, float cy, float r) {
    return [=](float x, float y) {
        const float nx = (x - cx) / r, ny = (y - cy) / r;
        return c * (0.82f + 0.38f * Clamp01(-0.6f * nx - 0.8f * ny + 0.3f) - 0.2f * Clamp01(nx * nx + ny * ny - 0.6f));
    };
}
Shade Flat(C3 c) { return [=](float, float) { return c; }; }
Shade Speckle(Shade s, float amt, int seed, int cell = 2) {
    return [=](float x, float y) { return s(x, y) * (1.0f + amt * (Hash((int)x / cell, (int)y / cell, seed) * 2.0f - 1.0f)); };
}

struct Canvas {
    std::vector<float> px = std::vector<float>(kIcon * kIcon * 4, 0.0f); // premultiplied RGBA
    int seed = 1;
    // Fill a shape: anti-aliased edge, shaded body, an outline of the shade at 35% brightness.
    void paint(const Sdf& sdf, const Shade& shade, float outline = 1.4f, float alpha = 1.0f) {
        for (int y = 0; y < kIcon; ++y)
            for (int x = 0; x < kIcon; ++x) {
                const float fx = x + 0.5f, fy = y + 0.5f, d = sdf(fx, fy);
                if (d >= 0.5f) continue;
                const float a = Clamp01(0.5f - d) * alpha;
                C3 c = shade(fx, fy);
                if (outline > 0) c = Mix(c, c * 0.32f, Clamp01(d + outline + 0.5f) * 0.9f);
                float* p = &px[(y * kIcon + x) * 4];
                const float ia = 1.0f - a;
                p[0] = Clamp01(c.r) * a + p[0] * ia; p[1] = Clamp01(c.g) * a + p[1] * ia; p[2] = Clamp01(c.b) * a + p[2] * ia;
                p[3] = a + p[3] * ia;
            }
    }
    void line(float ax, float ay, float bx, float by, float w, C3 c, float alpha = 1.0f) { paint(Capsule(ax, ay, bx, by, w), Flat(c), 0, alpha); }
    void shine(float cx, float cy, float rx, float ry, float alpha = 0.55f) { paint(Ellipse(cx, cy, rx, ry), Flat(kWhite), 0, alpha); }
};

// ---------------------------------------------------------------------------------------------
// Painters
// ---------------------------------------------------------------------------------------------

std::vector<V2> LumpPoints(int seed, float cx, float cy, float r, int n = 9) {
    std::vector<V2> v;
    for (int i = 0; i < n; ++i) {
        const float a = 6.2831853f * i / n + 0.3f * Hash(i, 3, seed);
        const float k = r * (0.78f + 0.3f * Hash(i, 7, seed));
        v.push_back({ cx + cosf(a) * k, cy + sinf(a) * k * 0.85f });
    }
    return v;
}
void Lump(Canvas& cv, C3 c, float speckle = 0.12f) {
    cv.paint(Poly(LumpPoints(cv.seed, 32, 35, 24)), Speckle(Round(c, 30, 30, 26), speckle, cv.seed));
    cv.shine(24, 24, 5, 3, 0.25f);
}
void Ore(Canvas& cv, C3 c) {
    cv.paint(Poly(LumpPoints(cv.seed, 32, 35, 24)), Speckle(Round({ 0.45f, 0.44f, 0.43f }, 30, 30, 26), 0.12f, cv.seed));
    const V2 spots[5] = { { 22, 28 }, { 38, 24 }, { 34, 42 }, { 44, 36 }, { 24, 42 } };
    for (int i = 0; i < 5; ++i) {
        const V2 s = spots[i];
        const float r = 3.5f + 2.5f * Hash(i, 1, cv.seed);
        cv.paint(Poly({ { s.x, s.y - r * 1.3f }, { s.x + r, s.y }, { s.x, s.y + r * 1.3f }, { s.x - r, s.y } }), Lin(c * 1.15f), 0.9f);
        cv.shine(s.x - r * 0.3f, s.y - r * 0.4f, r * 0.3f, r * 0.25f, 0.6f);
    }
}
void Ingot(Canvas& cv, C3 c) {
    // A cast bar seen from above-left: sloped top, front and side faces.
    cv.paint(Poly({ { 8, 36 }, { 34, 50 }, { 34, 56 }, { 8, 42 } }), Flat(c * 0.75f), 1.2f);
    cv.paint(Poly({ { 34, 50 }, { 56, 34 }, { 56, 40 }, { 34, 56 } }), Flat(c * 0.55f), 1.2f);
    cv.paint(Poly({ { 8, 36 }, { 30, 22 }, { 56, 34 }, { 34, 50 } }), Flat(c * 0.9f), 1.2f);
    cv.paint(Poly({ { 16, 35 }, { 31, 26 }, { 47, 34 }, { 33, 44 } }), Lin(c * 1.2f), 0.8f);
    cv.line(20, 34, 30, 28, 1.1f, kWhite, 0.45f);
}
void Shards(Canvas& cv, C3 c) {
    const float s[3][4] = { { 22, 38, 8, -0.3f }, { 36, 30, 11, 0.15f }, { 44, 44, 7, 0.5f } };
    for (auto& q : s) {
        Sdf d = Rot(Poly({ { q[0], q[1] - q[2] * 1.8f }, { q[0] + q[2] * 0.6f, q[1] }, { q[0], q[1] + q[2] }, { q[0] - q[2] * 0.6f, q[1] } }), q[0], q[1], q[3]);
        cv.paint(d, Lin(Mix(c, kWhite, 0.15f)), 1.2f);
        cv.paint(Rot(Poly({ { q[0], q[1] - q[2] * 1.8f }, { q[0] + q[2] * 0.6f, q[1] }, { q[0], q[1] } }), q[0], q[1], q[3]), Flat(c * 0.7f), 0);
    }
}
void Pile(Canvas& cv, C3 c) {
    cv.paint(Union(Ellipse(32, 46, 26, 10), Ellipse(32, 36, 17, 16)), Speckle(Round(c, 30, 34, 26), 0.18f, cv.seed, 1));
}
void Pebbles(Canvas& cv, C3 c) {
    const float p[6][3] = { { 20, 42, 9 }, { 40, 44, 10 }, { 30, 30, 9 }, { 46, 30, 7 }, { 16, 28, 6 }, { 30, 50, 6 } };
    for (int i = 0; i < 6; ++i) {
        const C3 k = c * (0.8f + 0.4f * Hash(i, 2, cv.seed));
        cv.paint(Ellipse(p[i][0], p[i][1], p[i][2], p[i][2] * 0.8f), Round(k, p[i][0], p[i][1], p[i][2]));
    }
}
void Gem(Canvas& cv, C3 c) {
    cv.paint(Poly({ { 12, 24 }, { 22, 12 }, { 42, 12 }, { 52, 24 }, { 32, 54 } }), Flat(c), 1.6f);
    cv.paint(Poly({ { 22, 12 }, { 42, 12 }, { 38, 24 }, { 26, 24 } }), Flat(c * 1.35f), 0);
    cv.paint(Poly({ { 12, 24 }, { 26, 24 }, { 32, 54 } }), Flat(c * 0.75f), 0);
    cv.paint(Poly({ { 38, 24 }, { 52, 24 }, { 32, 54 } }), Flat(c * 0.6f), 0);
    cv.paint(Poly({ { 26, 24 }, { 38, 24 }, { 32, 54 } }), Flat(c * 1.0f), 0);
    cv.paint(Poly({ { 12, 24 }, { 22, 12 }, { 26, 24 } }), Flat(c * 1.1f), 0);
    cv.paint(Poly({ { 42, 12 }, { 52, 24 }, { 38, 24 } }), Flat(c * 0.85f), 0);
    cv.shine(27, 16, 4, 2, 0.7f);
    cv.paint(Poly({ { 12, 24 }, { 22, 12 }, { 42, 12 }, { 52, 24 }, { 32, 54 } }), Flat(c), 0, 0.0f);
}

// An isometric block; pattern(u, v, face) returns a brightness factor for face-local u, v in [0, 1].
enum class Pat { Plain, Stone, Brick, Planks, Tiles, Glass, Soil, Gem };
float Pattern(Pat p, float u, float v, int face, int seed) {
    switch (p) {
    case Pat::Stone: return 1.0f + 0.16f * (Hash((int)(u * 12), (int)(v * 12), seed + face) - 0.5f);
    case Pat::Brick: {
        const int row = (int)(v * 4); const float off = (row & 1) ? 0.25f : 0.0f;
        const float fu = fmodf(u * 2.0f + off, 1.0f), fv = fmodf(v * 4.0f, 1.0f);
        return (fu < 0.07f || fv < 0.12f) ? 0.55f : 1.0f + 0.1f * (Hash((int)(u * 2 + off), row, seed) - 0.5f);
    }
    case Pat::Planks: { const float fv = fmodf(v * 3.0f, 1.0f); return fv < 0.08f ? 0.6f : 1.0f + 0.06f * sinf(u * 40.0f + v * 9.0f); }
    case Pat::Tiles: { const float fu = fmodf(u * 2.0f, 1.0f), fv = fmodf(v * 2.0f, 1.0f); return (fu < 0.06f || fv < 0.06f) ? 0.6f : 1.0f + 0.25f * (0.5f - fu) * (0.5f - fv); }
    case Pat::Soil: { const float fu = fmodf(u * 4.0f, 1.0f); return (face == 0 ? (fu < 0.3f ? 0.75f : 1.0f) : 1.0f) * (1.0f + 0.2f * (Hash((int)(u * 16), (int)(v * 16), seed) - 0.5f)); }
    case Pat::Gem: { const float fu = fmodf(u * 2.0f, 1.0f), fv = fmodf(v * 2.0f, 1.0f); return (fu < 0.06f || fv < 0.06f) ? 0.7f : 1.0f + 0.35f * Clamp01(0.6f - fu - fv); }
    default: return 1.0f + 0.05f * (Hash((int)(u * 8), (int)(v * 8), seed + face) - 0.5f);
    }
}
void Cube(Canvas& cv, C3 c, Pat pat, C3 top = { -1, 0, 0 }) {
    const float alpha = pat == Pat::Glass ? 0.6f : 1.0f;
    const C3 topC = top.r >= 0 ? top : c;
    const int seed = cv.seed;
    cv.paint(Poly({ { 32, 6 }, { 56, 18 }, { 32, 30 }, { 8, 18 } }), [=](float x, float y) {
        const float u = ((x - 32) / 24 + (y - 6) / 12) * 0.5f, v = ((y - 6) / 12 - (x - 32) / 24) * 0.5f;
        return topC * (1.18f * Pattern(pat, u, v, 0, seed));
    }, 1.2f, alpha);
    cv.paint(Poly({ { 8, 18 }, { 32, 30 }, { 32, 58 }, { 8, 46 } }), [=](float x, float y) {
        const float u = (x - 8) / 24, v = (y - (18 + (x - 8) * 0.5f)) / 28;
        return c * (0.86f * Pattern(pat, u, v, 1, seed));
    }, 1.2f, alpha);
    cv.paint(Poly({ { 32, 30 }, { 56, 18 }, { 56, 46 }, { 32, 58 } }), [=](float x, float y) {
        const float u = (x - 32) / 24, v = (y - (30 - (x - 32) * 0.5f)) / 28;
        return c * (0.62f * Pattern(pat, u, v, 2, seed));
    }, 1.2f, alpha);
    if (pat == Pat::Glass) { cv.line(14, 30, 22, 26, 1.2f, kWhite, 0.7f); cv.line(38, 44, 50, 30, 1.2f, kWhite, 0.5f); }
}
void Log(Canvas& cv, C3 c) {
    const int seed = cv.seed;
    cv.paint(Box(28, 34, 20, 13, 3, -0.35f), [=](float x, float y) { return c * ((0.85f + 0.3f * Clamp01((52 - y) / 30)) * (1.0f + 0.18f * (Hash((int)(x * 0.6f + y * 0.25f), 0, seed) - 0.5f))); });
    const C3 wood = Mix(c, { 0.85f, 0.68f, 0.45f }, 0.65f);
    cv.paint(Rot(Ellipse(47, 27, 6.5f, 13), 47, 27, -0.35f), [=](float x, float y) {
        const float r = hypotf((x - 47) / 6.5f, (y - 27) / 13.0f);
        return wood * (0.9f + 0.15f * sinf(r * 18.0f));
    });
}
void Planks(Canvas& cv, C3 c) {
    for (int i = 0; i < 3; ++i) {
        const float y = 22 + i * 11;
        const C3 k = c * (1.0f + 0.08f * (1 - i));
        cv.paint(Box(32, y, 25, 4.5f, 1, -0.3f), [=](float x, float yy) { return k * (1.05f + 0.08f * sinf(x * 1.3f + yy * 0.4f) - 0.15f * (x / kIcon)); });
    }
}
void Plate(Canvas& cv, C3 c) {
    cv.paint(Poly({ { 8, 32 }, { 32, 44 }, { 32, 50 }, { 8, 38 } }), Flat(c * 0.8f), 1.2f);
    cv.paint(Poly({ { 32, 44 }, { 56, 30 }, { 56, 36 }, { 32, 50 } }), Flat(c * 0.6f), 1.2f);
    cv.paint(Poly({ { 8, 32 }, { 32, 18 }, { 56, 30 }, { 32, 44 } }), Lin(c * 1.12f), 1.2f);
    const V2 r[4] = { { 18, 31 }, { 32, 24 }, { 46, 31 }, { 32, 38 } };
    for (const V2& p : r) cv.paint(Ellipse(p.x, p.y, 2, 1.4f), Flat(c * 0.75f), 0);
    cv.line(20, 30, 30, 25, 1.0f, kWhite, 0.35f);
}
void Bottle(Canvas& cv, C3 liquid, C3 glass = { 0.85f, 0.92f, 0.95f }) {
    cv.paint(Box(32, 40, 15, 17, 7), Flat(glass), 1.2f, 0.45f);
    cv.paint(Inter(Box(32, 40, 13.5f, 15.5f, 6), [](float, float y) { return 30.0f - y; }), Lin(liquid), 0.6f);
    cv.paint(Box(32, 19, 6, 6, 1.5f), Flat(glass), 1.0f, 0.5f);
    cv.paint(Box(32, 11, 7, 3.5f, 1.5f), Lin({ 0.55f, 0.38f, 0.22f }), 1.0f);
    cv.line(23, 30, 23, 48, 1.6f, kWhite, 0.55f);
}
void Canister(Canvas& cv, C3 c) {
    cv.paint(Box(32, 37, 12, 21, 7), Round(c, 26, 30, 22));
    cv.paint(Box(32, 29, 12.5f, 2, 0), Flat(c * 0.55f), 0);
    cv.paint(Box(32, 13, 6, 4, 1.5f), Lin(kSteel));
    cv.line(25, 22, 25, 50, 1.5f, kWhite, 0.4f);
}
void Cell(Canvas& cv, C3 c) {
    cv.paint(Box(32, 13, 4, 3, 1), Lin(kSteel));
    cv.paint(Box(32, 35, 12, 20, 3), Lin(c));
    cv.paint(Box(32, 20, 12, 4, 2), Lin(kDark * 1.6f), 0.8f);
    cv.line(32, 32, 32, 42, 1.4f, kWhite, 0.9f); cv.line(27, 37, 37, 37, 1.4f, kWhite, 0.9f);
}
void Berries(Canvas& cv, C3 c) {
    cv.paint(Rot(Ellipse(40, 14, 9, 4), 40, 14, -0.5f), Lin(kLeaf));
    const float b[6][2] = { { 24, 30 }, { 38, 28 }, { 31, 40 }, { 18, 42 }, { 44, 41 }, { 30, 52 } };
    for (auto& p : b) { cv.paint(Circle(p[0], p[1], 7.5f), Round(c, p[0], p[1], 7.5f)); cv.shine(p[0] - 2.5f, p[1] - 2.5f, 1.8f, 1.4f, 0.6f); }
}
void Bread(Canvas& cv, C3 c) {
    cv.paint(Ellipse(32, 38, 25, 14), Round(c, 28, 32, 26));
    for (int i = 0; i < 3; ++i) cv.line(20.0f + i * 10, 32, 26.0f + i * 10, 40, 1.5f, Mix(c, kWhite, 0.4f), 0.8f);
}
void Meat(Canvas& cv, C3 c) {
    cv.paint(Capsule(46, 18, 54, 10, 3), Flat({ 0.92f, 0.9f, 0.84f }), 1.0f);
    cv.paint(Circle(55, 9, 4), Flat({ 0.92f, 0.9f, 0.84f }), 1.0f);
    cv.paint(Poly({ { 10, 34 }, { 18, 18 }, { 36, 14 }, { 50, 22 }, { 52, 40 }, { 40, 54 }, { 20, 52 } }), Round(c, 30, 32, 24), 1.6f);
    cv.paint(Poly({ { 14, 34 }, { 20, 22 }, { 34, 19 }, { 46, 25 }, { 47, 38 }, { 38, 49 }, { 22, 47 } }), Flat(Mix(c, kWhite, 0.15f)), 0, 0.35f);
}
void Egg(Canvas& cv, C3 c) { cv.paint(Ellipse(32, 35, 15, 20), Round(c, 28, 28, 20)); cv.shine(26, 26, 3, 5, 0.5f); }
void Wheat(Canvas& cv, C3 c) {
    cv.line(32, 58, 32, 14, 1.4f, c * 0.8f);
    for (int i = 0; i < 5; ++i) {
        const float y = 14.0f + i * 7;
        cv.paint(Rot(Ellipse(27, y, 5, 2.6f), 27, y, 0.6f), Lin(c));
        cv.paint(Rot(Ellipse(37, y, 5, 2.6f), 37, y, -0.6f), Lin(c));
    }
    cv.paint(Ellipse(32, 9, 2.6f, 5), Lin(c));
}
void Seeds(Canvas& cv, C3 c) {
    const float s[7][3] = { { 22, 26, 0.4f }, { 36, 22, -0.6f }, { 44, 34, 0.9f }, { 28, 38, -0.2f }, { 18, 46, 0.7f }, { 38, 48, 0.2f }, { 48, 50, -0.8f } };
    for (auto& p : s) cv.paint(Rot(Ellipse(p[0], p[1], 5.5f, 3.2f), p[0], p[1], p[2]), Round(c, p[0], p[1], 5));
}
void Root(Canvas& cv, C3 c) {
    cv.paint(Rot(Ellipse(18, 14, 7, 2.5f), 18, 14, -0.8f), Lin(kLeaf));
    cv.paint(Rot(Ellipse(26, 10, 7, 2.5f), 26, 10, -1.3f), Lin(kLeaf));
    cv.paint(Poly({ { 16, 20 }, { 30, 14 }, { 36, 22 }, { 52, 54 } }), Round(c, 26, 22, 24));
    for (int i = 0; i < 3; ++i) cv.line(26.0f + i * 6, 24.0f + i * 8, 31.0f + i * 6, 22.0f + i * 8, 0.8f, c * 0.6f);
}
void Sack(Canvas& cv, C3 c) {
    cv.paint(Poly({ { 22, 20 }, { 42, 20 }, { 52, 40 }, { 48, 56 }, { 16, 56 }, { 12, 40 } }), Round(c, 28, 34, 26));
    cv.paint(Poly({ { 24, 20 }, { 40, 20 }, { 44, 10 }, { 20, 10 } }), Lin(c * 0.95f));
    cv.paint(Box(32, 20, 10, 2, 1), Flat({ 0.55f, 0.38f, 0.22f }), 0.8f);
}
void Fibre(Canvas& cv, C3 c) {
    for (int i = 0; i < 5; ++i) {
        const float x = 18.0f + i * 7;
        cv.line(x, 56, x + 6 - i * 2.0f, 10 + (i & 1) * 4, 1.6f, c * (0.85f + 0.08f * i));
    }
    cv.paint(Box(32, 38, 16, 2.5f, 1), Flat(c * 0.6f), 0.6f);
}
void Coil(Canvas& cv, C3 c, bool copper = false) {
    cv.paint(Ring(30, 34, 17, 3.5f), Lin(c));
    cv.paint(Ring(30, 34, 9.5f, 3.5f), Lin(c * 0.9f));
    cv.line(46, 36, 56, 52, 3.5f, c * 0.95f);
    if (copper) cv.paint(Circle(56, 52, 3), Flat({ 0.85f, 0.5f, 0.25f }), 0.8f);
}
void Wool(Canvas& cv, C3 c) {
    const float p[7][3] = { { 22, 32, 10 }, { 36, 26, 11 }, { 44, 38, 10 }, { 30, 42, 11 }, { 18, 44, 8 }, { 40, 48, 8 }, { 28, 22, 8 } };
    for (auto& q : p) cv.paint(Circle(q[0], q[1], q[2]), Round(c, q[0] - 2, q[1] - 2, q[2]), 1.0f);
}
void Hide(Canvas& cv, C3 c) {
    cv.paint(Poly({ { 14, 12 }, { 26, 18 }, { 38, 18 }, { 50, 12 }, { 46, 30 }, { 52, 50 }, { 38, 46 }, { 32, 56 }, { 26, 46 }, { 12, 50 }, { 18, 30 } }), Speckle(Lin(c), 0.1f, cv.seed, 3));
}
void Feather(Canvas& cv, C3 c) {
    cv.paint(Rot(Ellipse(32, 30, 9, 22), 32, 30, 0.5f), Lin(c));
    cv.line(20, 54, 44, 8, 1.2f, Mix(c, kWhite, 0.5f));
    for (int i = 0; i < 4; ++i) cv.line(26.0f + i * 5, 42.0f - i * 9, 22.0f + i * 5, 36.0f - i * 9, 0.6f, c * 0.7f);
}
void Cloth(Canvas& cv, C3 c) {
    cv.paint(Box(32, 34, 20, 18, 2, 0.1f), Lin(c));
    cv.paint(Poly({ { 30, 52 }, { 52, 30 }, { 52, 52 } }), Lin(c * 0.8f), 1.0f);
    for (int i = 0; i < 6; ++i) cv.line(16.0f + i * 5, 20.5f + i * 0.5f, 18.0f + i * 5, 20.7f + i * 0.5f, 0.6f, kWhite, 0.7f);
}
void Tube(Canvas& cv, C3 c) {
    cv.paint(Box(32, 50, 10, 6, 1.5f), Lin(kDark * 2.0f));
    for (int i = -1; i <= 1; ++i) cv.line(32.0f + i * 5, 56, 32.0f + i * 5, 61, 1.0f, kSteel);
    cv.paint(Ellipse(32, 28, 12, 19), Flat(Mix(c, kWhite, 0.4f)), 1.2f, 0.45f);
    cv.line(28, 38, 30, 22, 0.9f, kGlow); cv.line(36, 38, 34, 22, 0.9f, kGlow); cv.line(30, 22, 34, 22, 0.9f, kGlow);
    cv.shine(26, 20, 2, 6, 0.5f);
}
void Circuit(Canvas& cv, C3 c, int chips) {
    cv.paint(Box(32, 32, 25, 19, 2), Lin(c));
    for (int i = 0; i < 4; ++i) cv.line(10.0f, 18.0f + i * 9, 54.0f, 18.0f + i * 9, 0.7f, kGold, 0.8f);
    for (int i = 0; i < chips; ++i) {
        const float x = 18.0f + (i % 3) * 14, y = 26.0f + (i / 3) * 12;
        cv.paint(Box(x, y, 5, 4, 0.5f), Lin(kDark), 0.6f);
    }
}
void Hull(Canvas& cv, C3 c) {
    cv.paint(Box(32, 32, 23, 23, 3), Lin(c));
    cv.paint(Box(32, 32, 15, 15, 2), Lin(c * 0.7f), 1.0f);
    for (int i = 0; i < 4; ++i) cv.paint(Circle(i & 1 ? 50.0f : 14.0f, i & 2 ? 50.0f : 14.0f, 2.2f), Flat(c * 1.2f), 0.6f);
}

enum class Sym { None, Flame, Bolt, Gear, Drop, Flask, Teeth, Ring, Arrow, Hammer, Plus, Drill, Snow, Screen, Ingots };
void Symbol(Canvas& cv, Sym s, float cx, float cy, float r, C3 col) {
    switch (s) {
    case Sym::Flame:
        cv.paint(Poly({ { cx, cy - r }, { cx + r * 0.7f, cy + r * 0.2f }, { cx + r * 0.5f, cy + r * 0.8f }, { cx - r * 0.5f, cy + r * 0.8f }, { cx - r * 0.7f, cy + r * 0.2f } }), Flat(kGlow), 0.6f);
        cv.paint(Poly({ { cx, cy - r * 0.2f }, { cx + r * 0.35f, cy + r * 0.4f }, { cx, cy + r * 0.75f }, { cx - r * 0.35f, cy + r * 0.4f } }), Flat({ 1.0f, 0.92f, 0.5f }), 0);
        break;
    case Sym::Bolt: cv.paint(Poly({ { cx + r * 0.2f, cy - r }, { cx - r * 0.5f, cy + r * 0.1f }, { cx, cy + r * 0.1f }, { cx - r * 0.2f, cy + r }, { cx + r * 0.5f, cy - r * 0.15f }, { cx, cy - r * 0.15f } }), Flat({ 1.0f, 0.85f, 0.2f }), 0.6f); break;
    case Sym::Gear: {
        Sdf g = Circle(cx, cy, r * 0.65f);
        for (int i = 0; i < 8; ++i) { const float a = i * 0.785f; g = Union(g, Box(cx + cosf(a) * r * 0.75f, cy + sinf(a) * r * 0.75f, r * 0.16f, r * 0.2f, 0, a + 1.5708f)); }
        cv.paint(Cut(g, Circle(cx, cy, r * 0.28f)), Flat(col), 0.6f);
        break;
    }
    case Sym::Drop: cv.paint(Union(Circle(cx, cy + r * 0.3f, r * 0.55f), Poly({ { cx, cy - r }, { cx + r * 0.5f, cy + r * 0.1f }, { cx - r * 0.5f, cy + r * 0.1f } })), Flat(kWater * 1.2f), 0.6f); break;
    case Sym::Flask:
        cv.paint(Poly({ { cx - r * 0.2f, cy - r }, { cx + r * 0.2f, cy - r }, { cx + r * 0.2f, cy - r * 0.3f }, { cx + r * 0.8f, cy + r * 0.8f }, { cx - r * 0.8f, cy + r * 0.8f }, { cx - r * 0.2f, cy - r * 0.3f } }), Flat(col), 0.6f);
        cv.paint(Poly({ { cx - r * 0.45f, cy + r * 0.2f }, { cx + r * 0.45f, cy + r * 0.2f }, { cx + r * 0.75f, cy + r * 0.75f }, { cx - r * 0.75f, cy + r * 0.75f } }), Flat({ 0.4f, 0.9f, 0.5f }), 0);
        break;
    case Sym::Teeth: for (int i = 0; i < 4; ++i) { const float x = cx - r * 0.75f + i * r * 0.5f; cv.paint(Poly({ { x, cy - r * 0.6f }, { x + r * 0.5f, cy - r * 0.6f }, { x + r * 0.25f, cy + r * 0.3f } }), Flat(col), 0.5f);
                                                   cv.paint(Poly({ { x, cy + r * 0.7f }, { x + r * 0.5f, cy + r * 0.7f }, { x + r * 0.25f, cy - r * 0.1f } }), Flat(col * 0.8f), 0.5f); } break;
    case Sym::Ring: cv.paint(Ring(cx, cy, r * 0.65f, r * 0.18f), Flat(col), 0.5f); cv.paint(Circle(cx, cy, r * 0.22f), Flat(kGlow), 0.5f); break;
    case Sym::Arrow: cv.paint(Union(Poly({ { cx, cy - r }, { cx + r * 0.7f, cy - r * 0.2f }, { cx - r * 0.7f, cy - r * 0.2f } }), Box(cx, cy + r * 0.3f, r * 0.25f, r * 0.55f)), Flat(col), 0.5f); break;
    case Sym::Hammer: cv.paint(Union(Box(cx, cy - r * 0.5f, r * 0.75f, r * 0.3f, 1), Box(cx, cy + r * 0.3f, r * 0.15f, r * 0.65f)), Flat(col), 0.5f); break;
    case Sym::Plus: cv.paint(Union(Box(cx, cy, r * 0.7f, r * 0.2f), Box(cx, cy, r * 0.2f, r * 0.7f)), Flat(col), 0.5f); break;
    case Sym::Drill: cv.paint(Poly({ { cx - r * 0.6f, cy - r * 0.8f }, { cx + r * 0.6f, cy - r * 0.8f }, { cx, cy + r } }), Flat(col), 0.5f);
                     for (int i = 0; i < 3; ++i) cv.line(cx - r * 0.45f + i * r * 0.15f, cy - r * 0.4f + i * r * 0.45f, cx + r * 0.45f - i * r * 0.15f, cy - r * 0.6f + i * r * 0.45f, 0.8f, col * 0.6f);
                     break;
    case Sym::Snow: for (int i = 0; i < 3; ++i) { const float a = i * 1.047f; cv.line(cx - cosf(a) * r * 0.8f, cy - sinf(a) * r * 0.8f, cx + cosf(a) * r * 0.8f, cy + sinf(a) * r * 0.8f, 1.2f, { 0.7f, 0.9f, 1.0f }); } break;
    case Sym::Screen: cv.paint(Box(cx, cy, r * 0.8f, r * 0.6f, 1), Flat({ 0.15f, 0.55f, 0.35f }), 0.6f); cv.line(cx - r * 0.5f, cy, cx + r * 0.5f, cy - r * 0.2f, 0.8f, { 0.6f, 1.0f, 0.7f }); break;
    case Sym::Ingots: cv.paint(Box(cx - r * 0.35f, cy + r * 0.2f, r * 0.4f, r * 0.3f, 1), Flat({ 0.85f, 0.55f, 0.3f }), 0.6f);
                      cv.paint(Box(cx + r * 0.35f, cy + r * 0.2f, r * 0.4f, r * 0.3f, 1), Flat({ 0.8f, 0.82f, 0.86f }), 0.6f);
                      cv.paint(Box(cx, cy - r * 0.35f, r * 0.4f, r * 0.3f, 1), Flat(kGold), 0.6f); break;
    default: break;
    }
}
void Machine(Canvas& cv, C3 c, Sym s) {
    cv.paint(Box(32, 34, 25, 25, 3), Lin(c));
    cv.paint(Box(32, 12, 25, 3, 1), Flat(c * 0.7f), 0.6f);
    cv.paint(Box(32, 37, 16, 15, 2), Lin(kDark * 1.4f), 1.0f);
    Symbol(cv, s, 32, 37, 12, { 0.85f, 0.87f, 0.9f });
    for (int i = 0; i < 2; ++i) cv.paint(Circle(i ? 51.0f : 13.0f, 56, 1.8f), Flat(c * 1.3f), 0.4f);
}
void Campfire(Canvas& cv) {
    for (int i = 0; i < 5; ++i) { const float a = 3.1416f * (0.15f + 0.175f * i); cv.paint(Circle(32 + cosf(a) * 22, 50 + sinf(a) * 4 - 2, 4.5f), Round({ 0.5f, 0.5f, 0.52f }, 32 + cosf(a) * 22, 48, 5)); }
    cv.line(12, 54, 50, 42, 4, kWood); cv.line(14, 42, 52, 54, 4, kWood * 0.85f);
    Symbol(cv, Sym::Flame, 32, 30, 18, kGlow);
}
void Chest(Canvas& cv, C3 c) {
    cv.paint(Box(32, 40, 25, 16, 2), Lin(c));
    cv.paint(Box(32, 20, 25, 7, 4), Lin(c * 1.1f));
    for (int i = 0; i < 2; ++i) cv.paint(Box(i ? 48.0f : 16.0f, 34, 3, 22), Flat(kSteel * 0.8f), 0.6f);
    cv.paint(Box(32, 28, 4, 5, 1), Flat(kGold), 0.8f);
}
void RocketEngine(Canvas& cv, C3 c) {
    cv.paint(Box(32, 12, 10, 6, 1), Lin(kSteel));
    cv.paint(Poly({ { 22, 18 }, { 42, 18 }, { 54, 56 }, { 10, 56 } }), Lin(c));
    cv.paint(Ellipse(32, 56, 22, 4), Flat(kDark), 0.6f);
    for (int i = 0; i < 3; ++i) cv.line(24.0f - i * 4, 28.0f + i * 10, 40.0f + i * 4, 28.0f + i * 10, 0.8f, c * 0.6f);
}
void Tank(Canvas& cv, C3 c) {
    cv.paint(Box(32, 32, 14, 26, 12), Round(c, 26, 26, 24));
    cv.paint(Box(32, 20, 14.5f, 1.8f), Flat(c * 0.6f), 0); cv.paint(Box(32, 44, 14.5f, 1.8f), Flat(c * 0.6f), 0);
    cv.line(24, 14, 24, 50, 1.5f, kWhite, 0.35f);
}
void Cockpit(Canvas& cv, C3 c) {
    cv.paint(Poly({ { 32, 4 }, { 52, 48 }, { 52, 56 }, { 12, 56 }, { 12, 48 } }), Round(c, 26, 30, 30));
    cv.paint(Poly({ { 32, 18 }, { 42, 38 }, { 22, 38 } }), Lin({ 0.3f, 0.6f, 0.9f }), 1.0f);
    cv.shine(29, 28, 2, 4, 0.6f);
}
void Capsule3(Canvas& cv, C3 c) {
    cv.paint(Poly({ { 22, 12 }, { 42, 12 }, { 54, 50 }, { 10, 50 } }), Round(c, 26, 26, 30));
    cv.paint(Box(32, 52, 22, 4, 2), Lin(kDark * 1.5f));
    cv.paint(Circle(26, 32, 4), Lin({ 0.3f, 0.6f, 0.9f }), 0.8f); cv.paint(Circle(38, 32, 4), Lin({ 0.3f, 0.6f, 0.9f }), 0.8f);
}
void LandingLeg(Canvas& cv, C3 c) {
    cv.line(20, 10, 40, 48, 3.5f, c); cv.line(44, 14, 36, 40, 2.5f, c * 0.8f);
    cv.paint(Ellipse(42, 52, 14, 4.5f), Lin(c * 0.9f));
}
void Computer(Canvas& cv, C3 c) {
    cv.paint(Box(32, 34, 24, 20, 3), Lin(c));
    Symbol(cv, Sym::Screen, 32, 31, 14, c);
    for (int i = 0; i < 4; ++i) cv.paint(Circle(18.0f + i * 9, 48, 2), Flat(i == 0 ? C3{ 1, 0.3f, 0.3f } : kGold), 0.4f);
}

enum class ToolKind { Pickaxe, Shovel, Axe, Hoe };
void Tool(Canvas& cv, ToolKind t, C3 head, C3 handle) {
    const float ang = 0.72f; // the tool leans to the right
    auto R = [&](Sdf s) { return Rot(s, 32, 32, ang); };
    cv.paint(R(Capsule(32, 58, 32, 12, 3)), Lin(handle));
    Sdf h;
    switch (t) {
    case ToolKind::Pickaxe: h = Poly({ { 6, 20 }, { 18, 11 }, { 32, 8 }, { 46, 11 }, { 58, 20 }, { 55, 22 }, { 44, 16 }, { 32, 14 }, { 20, 16 }, { 9, 22 } }); break;
    case ToolKind::Shovel: h = Union(Poly({ { 22, 24 }, { 22, 10 }, { 32, 2 }, { 42, 10 }, { 42, 24 } }), Box(32, 26, 4, 3)); break;
    case ToolKind::Axe: h = Poly({ { 30, 8 }, { 50, 3 }, { 56, 15 }, { 52, 28 }, { 30, 22 } }); break;
    case ToolKind::Hoe: h = Poly({ { 30, 8 }, { 52, 8 }, { 54, 22 }, { 48, 22 }, { 46, 14 }, { 30, 14 } }); break;
    }
    cv.paint(R(h), Lin(head, 0.4f), 1.6f);
}
void Shears(Canvas& cv) {
    cv.paint(Rot(Poly({ { 30, 32 }, { 34, 32 }, { 50, 6 }, { 46, 6 } }), 32, 32, 0.25f), Lin(kSteel));
    cv.paint(Rot(Poly({ { 30, 32 }, { 34, 32 }, { 18, 6 }, { 14, 6 } }), 32, 32, -0.25f), Lin(kSteel * 0.9f));
    cv.paint(Ring(22, 46, 7, 2.2f), Flat({ 0.75f, 0.2f, 0.2f })); cv.paint(Ring(42, 46, 7, 2.2f), Flat({ 0.75f, 0.2f, 0.2f }));
}
void Bucket(Canvas& cv, C3 c, bool water) {
    cv.paint(Cut(Ring(32, 24, 16, 1.5f), Box(32, 34, 20, 10)), Flat(kSteel * 0.8f), 0.6f);
    cv.paint(Poly({ { 15, 24 }, { 49, 24 }, { 45, 56 }, { 19, 56 } }), Lin(c));
    cv.paint(Ellipse(32, 24, 17, 4.5f), Flat(water ? kWater : c * 0.45f), 1.0f);
    if (water) cv.shine(27, 23, 5, 1.2f, 0.6f);
}
void Telescope(Canvas& cv, C3 c) {
    cv.line(32, 36, 20, 58, 1.5f, kWood); cv.line(32, 36, 44, 58, 1.5f, kWood); cv.line(32, 36, 32, 58, 1.5f, kWood * 0.8f);
    cv.paint(Capsule(12, 42, 50, 18, 5.5f), Lin(c));
    cv.paint(Capsule(44, 22, 54, 15, 6.5f), Lin(c * 0.85f));
    cv.paint(Ellipse(55, 15, 2.5f, 5), Flat({ 0.5f, 0.75f, 0.95f }), 0.6f);
}
void MapIcon(Canvas& cv) {
    cv.paint(Poly({ { 8, 14 }, { 24, 10 }, { 40, 14 }, { 56, 10 }, { 56, 52 }, { 40, 56 }, { 24, 52 }, { 8, 56 } }), Lin({ 0.9f, 0.84f, 0.65f }));
    cv.paint(Ellipse(22, 30, 9, 7), Flat({ 0.45f, 0.7f, 0.4f }), 0, 0.8f);
    cv.paint(Ellipse(42, 40, 8, 6), Flat({ 0.4f, 0.6f, 0.85f }), 0, 0.8f);
    cv.line(40, 22, 46, 28, 1.3f, { 0.8f, 0.15f, 0.15f }); cv.line(46, 22, 40, 28, 1.3f, { 0.8f, 0.15f, 0.15f });
    cv.line(24, 10, 24, 52, 0.6f, { 0.5f, 0.45f, 0.3f }, 0.6f); cv.line(40, 14, 40, 56, 0.6f, { 0.5f, 0.45f, 0.3f }, 0.6f);
}
void Suit(Canvas& cv, C3 c) {
    cv.line(16, 30, 12, 50, 4.5f, c * 0.9f); cv.line(48, 30, 52, 50, 4.5f, c * 0.9f);
    cv.paint(Box(32, 42, 15, 16, 5), Lin(c));
    cv.paint(Circle(32, 17, 11), Round(c, 30, 14, 11));
    cv.paint(Ellipse(32, 18, 7, 5), Lin({ 0.15f, 0.25f, 0.4f }), 0.8f);
    cv.shine(29, 16, 2.5f, 1.4f, 0.6f);
    cv.paint(Box(32, 40, 6, 4, 1), Flat(c * 0.6f), 0.5f);
}
void Module(Canvas& cv, C3 c, Sym s) {
    cv.paint(Box(32, 32, 21, 21, 6), Lin(c));
    cv.paint(Box(32, 32, 14, 14, 4), Lin(kDark * 1.5f), 0.8f);
    Symbol(cv, s, 32, 32, 11, { 0.9f, 0.9f, 0.92f });
}
void Crucible(Canvas& cv, C3 c) {
    cv.paint(Poly({ { 14, 22 }, { 50, 22 }, { 44, 52 }, { 20, 52 } }), Round(c, 28, 30, 26));
    cv.paint(Ellipse(32, 22, 18, 5), Flat(c * 0.4f), 1.0f);
}
void Fence(Canvas& cv, C3 c, bool gate) {
    cv.paint(Box(14, 34, 4, 22, 1), Lin(c)); cv.paint(Box(50, 34, 4, 22, 1), Lin(c));
    cv.paint(Box(32, 24, 16, 3), Lin(c * 0.9f)); cv.paint(Box(32, 42, 16, 3), Lin(c * 0.9f));
    if (gate) cv.line(18, 42, 46, 24, 2.5f, c * 0.8f);
}
void Lantern(Canvas& cv, C3 c) {
    cv.paint(Ring(32, 10, 6, 1.3f), Flat(kDark * 2.0f), 0.5f);
    cv.paint(Box(32, 18, 12, 3, 1), Lin(kDark * 2.0f));
    cv.paint(Box(32, 36, 10, 15, 3), Flat(Mix(c, kGlow, 0.5f) * 1.2f), 1.0f, 0.95f);
    cv.paint(Ellipse(32, 36, 5, 8), Flat({ 1.0f, 0.95f, 0.7f }), 0);
    cv.paint(Box(32, 53, 13, 3, 1), Lin(kDark * 2.0f));
}
void Dye(Canvas& cv, C3 c) {
    cv.paint(Box(32, 40, 15, 16, 6), Round({ 0.62f, 0.42f, 0.3f }, 28, 34, 18));
    cv.paint(Ellipse(32, 25, 14, 4), Flat(c), 1.0f);
    cv.paint(Union(Capsule(40, 26, 40, 38, 2.2f), Circle(40, 39, 2.8f)), Flat(c), 0.6f);
}
void Pane(Canvas& cv, C3 c) {
    cv.paint(Box(32, 32, 20, 25, 1.5f), Flat(Mix(c, kWhite, 0.15f)), 1.4f, 0.55f);
    cv.line(18, 30, 30, 14, 1.3f, kWhite, 0.75f); cv.line(26, 50, 44, 26, 1.0f, kWhite, 0.55f);
}
void Conveyor(Canvas& cv, C3 c) {
    cv.paint(Box(32, 36, 26, 9, 8), Lin(c));
    for (int i = 0; i < 4; ++i) cv.paint(Circle(13.0f + i * 12.7f, 36, 3.5f), Flat(kSteel), 0.6f);
    Symbol(cv, Sym::Arrow, 32, 18, 8, { 0.95f, 0.8f, 0.2f });
}
void Pipe(Canvas& cv, C3 c, bool outlet) {
    cv.paint(Union(Capsule(8, 22, 40, 22, 7), Capsule(40, 22, 40, 50, 7)), Round(c, 30, 20, 30));
    cv.paint(Box(14, 22, 2.5f, 9), Lin(c * 0.75f)); cv.paint(Box(40, 46, 9, 2.5f), Lin(c * 0.75f));
    if (outlet) Symbol(cv, Sym::Drop, 40, 58, 6, kWater);
}
void Repair(Canvas& cv) {
    cv.paint(Cut(Ring(32, 20, 9, 2), Box(32, 28, 12, 6)), Flat(kDark * 2.0f), 0.5f);
    cv.paint(Box(32, 38, 24, 16, 3), Lin({ 0.8f, 0.2f, 0.18f }));
    Symbol(cv, Sym::Plus, 32, 38, 11, kWhite);
}

// One item -> its painter.
void PaintItem(Canvas& cv, const ItemDef& d) {
    const std::string& id = d.id;
    const C3 c = FromColor(d.color);
    auto ends = [&](const char* s) { const size_t n = strlen(s); return id.size() >= n && id.compare(id.size() - n, n, s) == 0; };
    auto starts = [&](const char* s) { return id.rfind(s, 0) == 0; };
    auto is = [&](std::initializer_list<const char*> ids) { for (const char* s : ids) if (id == s) return true; return false; };

    // Tools: <material>_<tool>, the head in the item's colour; stone tools have a wooden handle.
    if (ends("_pickaxe") || ends("_shovel") || ends("_axe") || ends("_hoe")) {
        const ToolKind t = ends("_pickaxe") ? ToolKind::Pickaxe : ends("_shovel") ? ToolKind::Shovel : ends("_axe") ? ToolKind::Axe : ToolKind::Hoe;
        Tool(cv, t, c, starts("titanium") || starts("steel") ? C3{ 0.3f, 0.3f, 0.34f } : kWood);
        return;
    }
    if (id == "campfire") { Campfire(cv); return; }
    if (id == "chest") { Chest(cv, c); return; }
    if (IsMachineItem(id)) {
        Sym s = Sym::Gear;
        if (is({ "kiln", "coke_oven", "blast_furnace", "electric_furnace", "arc_furnace" })) s = Sym::Flame;
        else if (id == "alloy_smelter") s = Sym::Ingots;
        else if (ends("generator")) s = Sym::Bolt;
        else if (id == "battery") s = Sym::Plus;
        else if (ends("macerator")) s = Sym::Teeth;
        else if (is({ "ore_washer", "water_pump", "boiler" })) s = Sym::Drop;
        else if (is({ "electrolyzer", "chemical_reactor", "refinery" })) s = Sym::Flask;
        else if (id == "centrifuge" || id == "fusion_reactor") s = Sym::Ring;
        else if (id == "miner") s = Sym::Drill;
        else if (id == "forge_hammer") s = Sym::Hammer;
        else if (id == "mass_driver") s = Sym::Arrow;
        Machine(cv, c, s);
        return;
    }
    if (ends("_ore")) { Ore(cv, c); return; }
    // Material forms (materials.txt): crushed, purified, dust, ingot, plate, rod, wire, gear.
    if (ends("_crushed")) { Pebbles(cv, Mix(c, { 0.5f, 0.5f, 0.5f }, 0.35f)); return; }
    if (ends("_purified")) { Shards(cv, c); return; }
    if (ends("_dust")) { Pile(cv, c); return; }
    if (ends("_ingot")) { Ingot(cv, c); return; }
    if (ends("_plate")) { Plate(cv, c); return; }
    if (ends("_rod")) { cv.paint(Capsule(14, 50, 50, 14, 4.5f), Round(c, 28, 28, 30)); cv.line(18, 44, 44, 18, 1.0f, kWhite, 0.4f); return; }
    if (ends("_wire")) { Coil(cv, c, true); return; }
    if (ends("_gear")) { Symbol(cv, Sym::Gear, 32, 32, 26, c); cv.shine(24, 20, 5, 2.5f, 0.35f); return; }
    if (ends("_log") || id == "any_log") { Log(cv, c); return; }
    if (ends("_planks")) { Planks(cv, c); return; }
    if (ends("_dye")) { Dye(cv, c); return; }
    if (ends("_plaster")) { Cube(cv, c, Pat::Plain); return; }
    if (ends("_glass") || id == "glass") { Pane(cv, c); return; }
    if (ends("_block")) { Cube(cv, c, Pat::Gem); return; }
    if (ends("_hull")) { Hull(cv, c); return; }
    if (ends("_circuit")) { Circuit(cv, c, starts("advanced") ? 6 : starts("good") ? 4 : 2); return; }
    if (ends("_suit")) { Suit(cv, c); return; }
    if (ends("_seeds")) { Seeds(cv, c); return; }
    if (is({ "quartz", "amethyst", "ruby", "sapphire", "emerald" })) { Gem(cv, c); return; }
    if (is({ "granite", "basalt", "sandstone", "limestone", "marble" })) { Cube(cv, c, Pat::Stone); return; }
    if (id == "brick") { Cube(cv, c, Pat::Brick); return; }
    if (is({ "concrete" })) { Cube(cv, c, Pat::Stone); return; }
    if (id == "heat_shield_tile") { Cube(cv, c, Pat::Tiles); return; }
    if (id == "farmland") { Cube(cv, { 0.4f, 0.27f, 0.16f }, Pat::Soil); return; }
    if (is({ "sand", "regolith", "dirt", "helium3_dust" })) { Pile(cv, c); return; }
    if (id == "gravel") { Pebbles(cv, c); return; }
    if (id == "ice") { cv.paint(Poly(LumpPoints(cv.seed, 32, 35, 24, 6)), Lin(c), 1.4f, 0.85f); cv.shine(24, 24, 6, 3, 0.6f); return; }
    if (is({ "coal", "charcoal", "coke" })) { Lump(cv, c, 0.35f); return; }
    if (is({ "stone", "clay", "anorthite" })) { Lump(cv, c, id == "clay" ? 0.04f : 0.12f); return; }
    if (is({ "iron_plate", "polyethylene" })) { Plate(cv, c); return; }
    if (id == "superconductor") { Coil(cv, c); return; }
    if (is({ "water", "milk", "rocket_fuel", "crude_oil" })) { Bottle(cv, c); return; }
    if (id == "peat") { Pile(cv, c); return; }
    if (id == "cactus") {
        cv.paint(Union(Box(32, 36, 7, 22, 6), Union(Union(Box(20, 32, 3.5f, 8, 3), Box(23, 38, 5, 3, 2)), Union(Box(44, 26, 3.5f, 8, 3), Box(41, 32, 5, 3, 2)))), Lin(c));
        for (int i = 0; i < 5; ++i) cv.line(30.0f, 18.0f + i * 7, 30.0f, 21.0f + i * 7, 0.6f, Mix(c, kWhite, 0.5f), 0.7f);
        return;
    }
    if (is({ "oxygen_canister", "hydrogen", "oxygen", "helium3", "oxygen_tank_module" })) { Canister(cv, c); return; }
    if (is({ "power_cell", "battery_pack" })) { Cell(cv, c); return; }
    if (id == "berries") { Berries(cv, c); return; }
    if (id == "bread") { Bread(cv, c); return; }
    if (ends("_meat")) { Meat(cv, c); return; }
    if (id == "egg") { Egg(cv, c); return; }
    if (id == "grain") { Wheat(cv, c); return; }
    if (ends("roots")) { Root(cv, c); return; }
    if (is({ "flour", "food_pouch" })) { Sack(cv, c); return; }
    if (id == "fibre") { Fibre(cv, c); return; }
    if (id == "rope") { Coil(cv, c); return; }
    if (id == "cable") { Coil(cv, c, true); return; }
    if (id == "wool") { Wool(cv, c); return; }
    if (id == "hide") { Hide(cv, c); return; }
    if (id == "feather") { Feather(cv, c); return; }
    if (id == "cloth") { Cloth(cv, c); return; }
    if (id == "vacuum_tube") { Tube(cv, c); return; }
    if (id == "conveyor") { Conveyor(cv, c); return; }
    if (id == "pipe" || id == "water_outlet") { Pipe(cv, c, id == "water_outlet"); return; }
    if (id == "rocket_engine") { RocketEngine(cv, c); return; }
    if (id == "fuel_tank") { Tank(cv, c); return; }
    if (id == "cockpit") { Cockpit(cv, c); return; }
    if (id == "guidance_computer" || id == "ore_scanner") { Computer(cv, c); return; }
    if (id == "landing_leg" || id == "servo_legs") { LandingLeg(cv, c); return; }
    if (is({ "cabin", "station_core", "cargo_capsule", "launched_capsule" })) { Capsule3(cv, c); return; }
    if (id == "shears") { Shears(cv); return; }
    if (id == "bucket" || id == "water_bucket") { Bucket(cv, kSteel, id == "water_bucket"); return; }
    if (id == "telescope") { Telescope(cv, c); return; }
    if (id == "map") { MapIcon(cv); return; }
    if (id == "repair_kit") { Repair(cv); return; }
    if (id == "crucible_mold" || id == "fired_crucible") { Crucible(cv, c); return; }
    if (id == "fence" || id == "gate") { Fence(cv, c, id == "gate"); return; }
    if (id == "lamp" || id == "lamp_upgrade") { Lantern(cv, c); return; }
    if (id == "any_planks") { Planks(cv, c); return; }
    if (id == "insulation_module") { Module(cv, c, Sym::Plus); return; }
    if (id == "rebreather") { Module(cv, c, Sym::Drop); return; }
    if (id == "thruster_pack") { Module(cv, c, Sym::Arrow); return; }
    if (id == "thermal_regulator") { Module(cv, c, Sym::Snow); return; }
    if (id == "grip_boots") { cv.paint(Poly({ { 20, 10 }, { 36, 10 }, { 36, 40 }, { 54, 44 }, { 54, 54 }, { 18, 54 } }), Lin(c)); cv.paint(Box(36, 55, 19, 2.5f), Flat(kDark), 0.5f); return; }
    // By kind.
    switch (d.kind) {
    case ItemKind::Block: Cube(cv, c, Pat::Plain); break;
    case ItemKind::Tool: Module(cv, c, Sym::Gear); break;
    default: Lump(cv, c); break;
    }
}

Texture2D atlas{};
Image atlasImage{};
std::unordered_map<std::string, int> cellOf;

} // namespace

void ItemIcons_Init() {
    ItemIcons_Unload();
    const auto& items = itemDB.items();
    const int n = (int)items.size();
    if (n == 0) return;
    const int rows = (n + kCols - 1) / kCols;
    atlasImage = GenImageColor(kCols * kIcon, rows * kIcon, BLANK);
    Color* out = (Color*)atlasImage.data;
    // Paint in parallel: each icon writes only its own cell.
    const int threads = std::max(1, std::min(8, (int)std::thread::hardware_concurrency()));
    std::vector<std::thread> pool;
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t] {
            for (int i = t; i < n; i += threads) {
                Canvas cv;
                cv.seed = (int)std::hash<std::string>{}(items[i].id) & 0x7FFF;
                PaintItem(cv, items[i]);
                const int ox = (i % kCols) * kIcon, oy = (i / kCols) * kIcon;
                for (int y = 0; y < kIcon; ++y)
                    for (int x = 0; x < kIcon; ++x) {
                        const float* p = &cv.px[(y * kIcon + x) * 4];
                        const float a = p[3], ia = a > 1e-4f ? 1.0f / a : 0.0f;
                        out[(oy + y) * atlasImage.width + ox + x] = { (unsigned char)(Clamp01(p[0] * ia) * 255.0f + 0.5f), (unsigned char)(Clamp01(p[1] * ia) * 255.0f + 0.5f),
                                                                      (unsigned char)(Clamp01(p[2] * ia) * 255.0f + 0.5f), (unsigned char)(Clamp01(a) * 255.0f + 0.5f) };
                    }
            }
        });
    for (auto& th : pool) th.join();
    for (int i = 0; i < n; ++i) cellOf[items[i].id] = i;
    atlas = LoadTextureFromImage(atlasImage);
    GenTextureMipmaps(&atlas);
    SetTextureFilter(atlas, TEXTURE_FILTER_TRILINEAR);
}

void ItemIcons_Unload() {
    if (atlas.id > 0) UnloadTexture(atlas);
    if (atlasImage.data) UnloadImage(atlasImage);
    atlas = Texture2D{};
    atlasImage = Image{};
    cellOf.clear();
}

bool ItemIcons_Has(const std::string& id) { return atlas.id > 0 && cellOf.count(id) > 0; }
Texture2D ItemIcons_Texture() { return atlas; }

Rectangle ItemIcons_Source(const std::string& id) {
    auto it = cellOf.find(id);
    if (it == cellOf.end()) return { 0, 0, 0, 0 };
    return { (float)((it->second % kCols) * kIcon), (float)((it->second / kCols) * kIcon), (float)kIcon, (float)kIcon };
}

void ItemIcons_Draw(const std::string& id, Rectangle dst, Color tint) {
    if (!ItemIcons_Has(id)) return;
    DrawTexturePro(atlas, ItemIcons_Source(id), dst, { 0, 0 }, 0.0f, tint);
}

void ItemIcons_UV(const std::string& id, float uv0[2], float uv1[2]) {
    const Rectangle r = ItemIcons_Source(id);
    const float w = (float)std::max(atlas.width, 1), h = (float)std::max(atlas.height, 1);
    uv0[0] = r.x / w; uv0[1] = r.y / h;
    uv1[0] = (r.x + r.width) / w; uv1[1] = (r.y + r.height) / h;
}

bool ItemIcons_SaveSheet(const char* path) {
    if (!atlasImage.data) return false;
    // Over a mid-grey so outlines and translucent glass read in the sheet.
    Image sheet = GenImageColor(atlasImage.width, atlasImage.height, Color{ 70, 74, 80, 255 });
    ImageDraw(&sheet, atlasImage, { 0, 0, (float)atlasImage.width, (float)atlasImage.height }, { 0, 0, (float)atlasImage.width, (float)atlasImage.height }, WHITE);
    const bool ok = ExportImage(sheet, path);
    UnloadImage(sheet);
    return ok;
}
