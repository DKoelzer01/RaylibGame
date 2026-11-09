#include "sealrecognizer.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

namespace {

// ---------------------------------------------------------------------------------------------
// $P point-cloud recogniser
// ---------------------------------------------------------------------------------------------
constexpr int kPoints = 32;

struct PPoint { float x, y; int id; };
using Cloud = std::vector<PPoint>;

float Dist(const PPoint& a, const PPoint& b) { return sqrtf((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y)); }

float PathLength(const Cloud& pts) {
    float d = 0.0f;
    for (size_t i = 1; i < pts.size(); ++i) if (pts[i].id == pts[i - 1].id) d += Dist(pts[i - 1], pts[i]);
    return d;
}

Cloud Resample(Cloud pts, int n) {
    const float I = PathLength(pts) / (n - 1);
    Cloud out;
    if (pts.empty()) return out;
    out.push_back(pts[0]);
    if (I <= 1e-6f) { while ((int)out.size() < n) out.push_back(pts[0]); return out; }
    float D = 0.0f;
    for (size_t i = 1; i < pts.size(); ++i) {
        if (pts[i].id != pts[i - 1].id) continue;
        const float d = Dist(pts[i - 1], pts[i]);
        if (D + d >= I) {
            const float t = (I - D) / d;
            PPoint q{ pts[i - 1].x + t * (pts[i].x - pts[i - 1].x), pts[i - 1].y + t * (pts[i].y - pts[i - 1].y), pts[i].id };
            out.push_back(q);
            pts.insert(pts.begin() + i, q); // q becomes the next start point
            D = 0.0f;
        } else {
            D += d;
        }
    }
    while ((int)out.size() < n) out.push_back(pts.back()); // rounding
    out.resize(n);
    return out;
}

Cloud Normalize(const std::vector<Stroke>& strokes, bool flipY) {
    Cloud pts;
    for (size_t s = 0; s < strokes.size(); ++s)
        for (const Vector2& p : strokes[s]) pts.push_back({ p.x, flipY ? -p.y : p.y, (int)s });
    pts = Resample(pts, kPoints);
    float minX = 1e9f, minY = 1e9f, maxX = -1e9f, maxY = -1e9f;
    for (auto& p : pts) { minX = std::min(minX, p.x); maxX = std::max(maxX, p.x); minY = std::min(minY, p.y); maxY = std::max(maxY, p.y); }
    const float size = std::max(std::max(maxX - minX, maxY - minY), 1e-6f);
    float cx = 0, cy = 0;
    for (auto& p : pts) { p.x = (p.x - minX) / size; p.y = (p.y - minY) / size; cx += p.x; cy += p.y; }
    cx /= pts.size(); cy /= pts.size();
    for (auto& p : pts) { p.x -= cx; p.y -= cy; }
    return pts;
}

float CloudDistance(const Cloud& a, const Cloud& b, int start) {
    const int n = (int)a.size();
    std::vector<bool> matched(n, false);
    float sum = 0.0f;
    int i = start;
    do {
        int index = -1;
        float best = 1e30f;
        for (int j = 0; j < n; ++j) {
            if (matched[j]) continue;
            const float d = Dist(a[i], b[j]);
            if (d < best) { best = d; index = j; }
        }
        matched[index] = true;
        const float weight = 1.0f - (float)((i - start + n) % n) / n;
        sum += weight * best;
        i = (i + 1) % n;
    } while (i != start);
    return sum;
}

float GreedyCloudMatch(const Cloud& a, const Cloud& b) {
    const int n = (int)a.size();
    const int step = (int)floorf(powf((float)n, 0.5f));
    float best = 1e30f;
    for (int i = 0; i < n; i += step) best = std::min({ best, CloudDistance(a, b, i), CloudDistance(b, a, i) });
    return best;
}

// ---------------------------------------------------------------------------------------------
// Templates (math coordinates: +y up, roughly within [-1, 1])
// ---------------------------------------------------------------------------------------------
struct Template { std::string name; std::vector<Stroke> strokes; Cloud cloud; };

Stroke Polyline(std::initializer_list<Vector2> pts, int perSegment = 12) {
    Stroke s;
    std::vector<Vector2> v(pts);
    for (size_t i = 0; i + 1 < v.size(); ++i)
        for (int k = 0; k < perSegment; ++k) s.push_back(Vector2Lerp(v[i], v[i + 1], (float)k / perSegment));
    s.push_back(v.back());
    return s;
}

Stroke Circle(float r, float start = 0.0f, float sweep = 2 * PI, int n = 48) {
    Stroke s;
    for (int i = 0; i <= n; ++i) { const float a = start + sweep * i / n; s.push_back({ r * cosf(a), r * sinf(a) }); }
    return s;
}

Stroke Rotate(Stroke s, float angle) {
    for (auto& p : s) p = Vector2Rotate(p, angle);
    return s;
}

std::vector<Template> BuildTemplates() {
    std::vector<Template> t;
    auto add = [&](const std::string& name, std::vector<Stroke> strokes) {
        Template tp{ name, strokes, {} };
        tp.cloud = Normalize(strokes, false);
        t.push_back(tp);
    };
    // Core signs
    add("flame", { Polyline({ { 0, 1 }, { 0.87f, -0.6f }, { -0.87f, -0.6f }, { 0, 1 } }) });               // triangle, apex up
    add("earth", { Polyline({ { -0.8f, 0.8f }, { 0.8f, 0.8f }, { 0.8f, -0.8f }, { -0.8f, -0.8f }, { -0.8f, 0.8f } }) }); // square
    {
        Stroke wave;
        for (int i = 0; i <= 48; ++i) { const float x = -1.0f + 2.0f * i / 48; wave.push_back({ x, 0.35f * sinf(x * 1.5f * PI) }); } // 1.5 periods
        add("water", { wave });                                                                            // wave
    }
    {
        Stroke spiral;
        for (int i = 0; i <= 64; ++i) { const float u = (float)i / 64, a = u * 4.0f * PI, r = 0.1f + 0.9f * u; spiral.push_back({ r * cosf(a), r * sinf(a) }); }
        add("wind", { spiral });                                                                           // spiral
    }
    add("light", { Polyline({ { 0, 1 }, { 0, -1 } }), Polyline({ { -0.87f, 0.5f }, { 0.87f, -0.5f } }),
                   Polyline({ { -0.87f, -0.5f }, { 0.87f, 0.5f } }) });                                    // six-point star
    // Modifier signs
    const Stroke chevron = Polyline({ { -0.5f, 0.8f }, { 0.5f, 0.0f }, { -0.5f, -0.8f } });
    for (int k = 0; k < 4; ++k) add("arrow", { Rotate(chevron, k * PI / 2) });                            // any direction
    add("orb", { Circle(1.0f) });                                                                          // small circle
    for (int k = 0; k < 4; ++k) add("line", { Rotate(Polyline({ { -1, 0 }, { 1, 0 } }), k * PI / 4) });   // stroke, any angle
    return t;
}

const std::vector<Template>& Templates() {
    static const std::vector<Template> t = BuildTemplates();
    return t;
}

const Template* FindTemplate(const std::string& name) {
    for (const Template& t : Templates()) if (t.name == name) return &t;
    return nullptr;
}

// Kasa algebraic circle fit: minimise sum (x^2 + y^2 + D x + E y + F)^2.
bool FitCircle(const Stroke& s, Vector2& c, float& r) {
    double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sz = 0, sxz = 0, syz = 0;
    const double n = (double)s.size();
    for (const Vector2& p : s) {
        const double x = p.x, y = p.y, z = x * x + y * y;
        sx += x; sy += y; sxx += x * x; syy += y * y; sxy += x * y; sz += z; sxz += x * z; syz += y * z;
    }
    // Normal equations for [D E F]: [sxx sxy sx; sxy syy sy; sx sy n] * [D E F] = -[sxz syz sz]
    const double A[3][3] = { { sxx, sxy, sx }, { sxy, syy, sy }, { sx, sy, n } };
    const double b[3] = { -sxz, -syz, -sz };
    const double det = A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
    if (fabs(det) < 1e-9) return false;
    auto solve = [&](int col) {
        double M[3][3];
        for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) M[i][j] = (j == col) ? b[i] : A[i][j];
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])) / det;
    };
    const double D = solve(0), E = solve(1), F = solve(2);
    c = { (float)(-D / 2), (float)(-E / 2) };
    const double rr = c.x * c.x + c.y * c.y - F;
    if (rr <= 0) return false;
    r = (float)sqrt(rr);
    return true;
}

RingFit EvaluateRing(const Stroke& s, int index) {
    RingFit f;
    if (s.size() < 16 || !FitCircle(s, f.center, f.radius) || f.radius < 20.0f) return f;
    double err = 0.0;
    std::vector<float> angles;
    for (const Vector2& p : s) {
        const float d = Vector2Distance(p, f.center) - f.radius;
        err += d * d;
        angles.push_back(atan2f(p.y - f.center.y, p.x - f.center.x));
    }
    const float rms = (float)sqrt(err / s.size());
    std::sort(angles.begin(), angles.end());
    float largestGap = angles.front() + 2 * PI - angles.back();
    for (size_t i = 1; i < angles.size(); ++i) largestGap = std::max(largestGap, angles[i] - angles[i - 1]);
    f.coverage = 1.0f - largestGap / (2 * PI);
    f.roundness = std::clamp(1.0f - (rms / f.radius) / 0.12f, 0.0f, 1.0f);
    const float gap = Vector2Distance(s.front(), s.back());
    f.closure = std::clamp(1.0f - gap / (0.35f * f.radius), 0.0f, 1.0f);
    f.quality = std::clamp(0.6f * f.roundness + 0.4f * f.closure, 0.0f, 1.0f) * std::clamp((f.coverage - 0.8f) / 0.15f, 0.0f, 1.0f);
    f.found = f.coverage > 0.85f && rms / f.radius < 0.2f;
    f.stroke = index;
    return f;
}

} // namespace

const std::vector<std::string>& CoreSigns() { static const std::vector<std::string> v = { "flame", "water", "wind", "light", "earth" }; return v; }
const std::vector<std::string>& ModifierSigns() { static const std::vector<std::string> v = { "arrow", "orb", "line" }; return v; }

GlyphMatch RecognizeGlyph(const std::vector<Stroke>& strokes) {
    GlyphMatch g;
    const Cloud c = Normalize(strokes, true); // screen space is y-down; templates are y-up
    float best = 1e30f;
    for (const Template& t : Templates()) {
        const float d = GreedyCloudMatch(c, t.cloud);
        if (d < best) { best = d; g.sign = t.name; }
    }
    g.score = std::max((best - 2.0f) / -2.0f, 0.0f);
    return g;
}

SealAnalysis AnalyzeSeal(const std::vector<Stroke>& strokes) {
    SealAnalysis a;
    // 1. The ring: the largest stroke that fits a circle well and sweeps most of the way round.
    for (size_t i = 0; i < strokes.size(); ++i) {
        const RingFit f = EvaluateRing(strokes[i], (int)i);
        if (f.found && (!a.ring.found || f.radius > a.ring.radius)) a.ring = f;
    }
    if (!a.ring.found) { a.problem = "No seal ring: draw a closed circle around your signs"; return a; }
    const float R = a.ring.radius;

    // 2. Group the strokes inside the ring into glyphs (bounding boxes closer than 8% of R).
    struct Box { float x0, y0, x1, y1; };
    std::vector<int> inside;
    std::vector<Box> boxes(strokes.size());
    for (size_t i = 0; i < strokes.size(); ++i) {
        if ((int)i == a.ring.stroke || strokes[i].empty()) continue;
        Box b{ 1e9f, 1e9f, -1e9f, -1e9f };
        Vector2 c{};
        for (const Vector2& p : strokes[i]) { b.x0 = std::min(b.x0, p.x); b.y0 = std::min(b.y0, p.y); b.x1 = std::max(b.x1, p.x); b.y1 = std::max(b.y1, p.y); c = Vector2Add(c, p); }
        c = Vector2Scale(c, 1.0f / strokes[i].size());
        boxes[i] = b;
        if (Vector2Distance(c, a.ring.center) < R * 1.1f) inside.push_back((int)i);
    }
    std::vector<int> group(strokes.size());
    std::iota(group.begin(), group.end(), 0);
    std::function<int(int)> root = [&](int x) { return group[x] == x ? x : group[x] = root(group[x]); };
    const float gapLimit = 0.08f * R;
    for (size_t i = 0; i < inside.size(); ++i)
        for (size_t j = i + 1; j < inside.size(); ++j) {
            const Box& p = boxes[inside[i]];
            const Box& q = boxes[inside[j]];
            const float dx = std::max({ 0.0f, p.x0 - q.x1, q.x0 - p.x1 });
            const float dy = std::max({ 0.0f, p.y0 - q.y1, q.y0 - p.y1 });
            if (sqrtf(dx * dx + dy * dy) <= gapLimit) group[root(inside[i])] = root(inside[j]);
        }
    std::vector<std::vector<int>> groups;
    std::vector<int> groupOf(strokes.size(), -1);
    for (int s : inside) {
        const int r = root(s);
        if (groupOf[r] < 0) { groupOf[r] = (int)groups.size(); groups.push_back({}); }
        groups[groupOf[r]].push_back(s);
    }

    // 3. Classify each group: tiny ones are power dots, the rest are signs.
    for (const auto& g : groups) {
        Box b{ 1e9f, 1e9f, -1e9f, -1e9f };
        std::vector<Stroke> gs;
        for (int s : g) { gs.push_back(strokes[s]); b.x0 = std::min(b.x0, boxes[s].x0); b.y0 = std::min(b.y0, boxes[s].y0); b.x1 = std::max(b.x1, boxes[s].x1); b.y1 = std::max(b.y1, boxes[s].y1); }
        const float size = sqrtf((b.x1 - b.x0) * (b.x1 - b.x0) + (b.y1 - b.y0) * (b.y1 - b.y0));
        if (size < 0.07f * R) { a.dots++; continue; }
        GlyphMatch m = RecognizeGlyph(gs);
        m.center = { (b.x0 + b.x1) * 0.5f, (b.y0 + b.y1) * 0.5f };
        m.size = size;
        m.polarR = Vector2Distance(m.center, a.ring.center) / R;
        m.angleDeg = atan2f(-(m.center.y - a.ring.center.y), m.center.x - a.ring.center.x) * RAD2DEG;
        m.strokes = g;
        a.glyphs.push_back(m);
    }

    // 4. Grammar: the most central core sign is the core; signs out toward the ring modify it.
    float bestR = 1e9f;
    for (size_t i = 0; i < a.glyphs.size(); ++i) {
        const GlyphMatch& g = a.glyphs[i];
        const bool isCore = std::find(CoreSigns().begin(), CoreSigns().end(), g.sign) != CoreSigns().end();
        if (isCore && g.polarR < 0.45f && g.polarR < bestR) { bestR = g.polarR; a.core = (int)i; }
    }
    for (size_t i = 0; i < a.glyphs.size(); ++i) if ((int)i != a.core && a.glyphs[i].polarR >= 0.3f) a.modifiers.push_back((int)i);
    if (a.core < 0) { a.problem = "The seal has no core sign at its centre"; return a; }
    a.quality = a.ring.quality * (0.5f + 0.5f * a.glyphs[a.core].score);
    return a;
}

// ---------------------------------------------------------------------------------------------
// Synthetic seals (tests)
// ---------------------------------------------------------------------------------------------

std::vector<Stroke> SynthesizeSeal(const std::string& core, const std::vector<std::string>& modifiers, int dots,
                                   float sloppiness, Vector2 center, float radius, std::mt19937& rng) {
    std::normal_distribution<float> n(0.0f, 1.0f);
    std::uniform_real_distribution<float> u(0.0f, 1.0f);
    std::vector<Stroke> out;
    auto place = [&](const Stroke& s, Vector2 at, float scale, float wobble) {
        // Template math coords (y up) -> screen (y down), plus hand jitter.
        Stroke r;
        for (const Vector2& p : s) r.push_back({ at.x + (p.x + n(rng) * wobble) * scale, at.y - (p.y + n(rng) * wobble) * scale });
        return r;
    };
    // Ring: a wobbly circle that may not quite close.
    {
        Stroke ring;
        const float start = u(rng) * 2 * PI;
        const float sweep = 2 * PI * (1.0f - 0.04f * sloppiness * u(rng));
        const float squash = 1.0f + 0.12f * sloppiness * (u(rng) - 0.5f);
        const int count = 96;
        float drift = 0.0f;
        for (int i = 0; i <= count; ++i) {
            const float a = start + sweep * i / count;
            drift += n(rng) * 0.004f * sloppiness;
            const float rr = radius * (1.0f + drift + n(rng) * 0.006f * sloppiness);
            ring.push_back({ center.x + rr * cosf(a) * squash, center.y + rr * sinf(a) / squash });
        }
        out.push_back(ring);
    }
    const float wobble = 0.03f + 0.05f * sloppiness;
    if (const Template* t = FindTemplate(core)) for (const Stroke& s : t->strokes) out.push_back(place(s, center, radius * 0.24f, wobble));
    for (size_t i = 0; i < modifiers.size(); ++i) {
        const float a = PI / 2 + (float)i * 2 * PI / std::max<size_t>(1, modifiers.size()) + 0.3f;
        const Vector2 at = { center.x + cosf(a) * radius * 0.66f, center.y - sinf(a) * radius * 0.66f };
        if (const Template* t = FindTemplate(modifiers[i])) for (const Stroke& s : t->strokes) out.push_back(place(s, at, radius * 0.15f, wobble));
    }
    for (int d = 0; d < dots; ++d) {
        const float a = -PI / 2 + (d - (dots - 1) * 0.5f) * 0.25f;
        const Vector2 at = { center.x + cosf(a) * radius * 0.82f, center.y - sinf(a) * radius * 0.82f };
        out.push_back(place(Circle(1.0f, 0, 2 * PI, 8), at, radius * 0.02f, 0.1f));
    }
    return out;
}
