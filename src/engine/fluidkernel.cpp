#include "fluidkernel.h"
#include "jobsystem.h"
#include "world.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <thread>

// A line-by-line port of resources/compute/fluid_step.comp (#654). Comments on the rules live there.
namespace FluidKernel {

void Grid::resize(int slots) {
    a.assign((size_t)slots * 4096, 0u); b.assign((size_t)slots * 4096, 0u);
    faceP.assign((size_t)slots * 4096, Vec4{}); faceN.assign((size_t)slots * 4096 * 2, Vec4{});
    nbr.assign((size_t)slots * 6, -1); processed.assign(slots, 0); flags.assign(slots, 0u);
}

namespace {

constexpr uint32_t SOLID = 1u, SOURCE = 2u, DUG = 8u, SINK = 16u;
constexpr float G = 9.81f;
inline uint32_t F(uint32_t c) { return c & 255u; }
inline uint32_t FL(uint32_t c) { return (c >> 8) & 255u; }
inline uint32_t TG(uint32_t c) { return (c >> 16) & 255u; }
inline uint32_t KIND(uint32_t c) { return (FL(c) >> 5) & 7u; }
inline uint32_t SH(uint32_t c) { return c >> 24; }
inline uint32_t CAP(uint32_t c) { return 255u - SH(c); }
inline uint32_t Pack(uint32_t f, uint32_t fl, uint32_t c) { return (f & 255u) | ((fl & 255u) << 8) | (c & 0xFFFF0000u); }
inline float Sign(float v) { return v > 0.0f ? 1.0f : v < 0.0f ? -1.0f : 0.0f; }
// A float rounded through IEEE half precision (nearest even), as the shader stores its face flows (#675: packHalf2x16).
inline float H(float f) {
    if (std::fabs(f) < 6.103515625e-05f) return std::nearbyint(f * 16777216.0f) / 16777216.0f; // half subnormals: 2^-24 steps
    if (std::fabs(f) > 65504.0f) return std::copysign(INFINITY, f);
    uint32_t x; std::memcpy(&x, &f, 4);
    x += 0x0FFFu + ((x >> 13) & 1u); x &= ~0x1FFFu; // 23 -> 10 mantissa bits
    std::memcpy(&f, &x, 4);
    return f;
}
inline Vec4 H(Vec4 v) { return { H(v.x), H(v.y), H(v.z), H(v.w) }; }

struct Pass {
    Grid& g;
    const Params& p;
    int pass;
    const uint32_t* inC;
    uint32_t* outC;

    int neighbour(int slot, int x, int y, int z, int dx, int dy, int dz, bool& processed, int& face) const {
        x += dx; y += dy; z += dz;
        face = -1;
        if (x < 0) { face = 0; x += 16; } else if (x >= 16) { face = 1; x -= 16; }
        else if (y < 0) { face = 2; y += 16; } else if (y >= 16) { face = 3; y -= 16; }
        else if (z < 0) { face = 4; z += 16; } else if (z >= 16) { face = 5; z -= 16; }
        int s = slot;
        if (face >= 0) { s = g.nbr[(size_t)slot * 6 + face]; if (s < 0) { processed = false; return -1; } }
        processed = g.processed[s] != 0;
        return s * 4096 + x + 16 * y + 256 * z;
    }
    static void cellPos(int gi, int& slot, int& x, int& y, int& z) {
        slot = gi / 4096; const int l = gi % 4096;
        x = l % 16; y = (l / 16) % 16; z = l / 256;
    }
    uint32_t readCell(int gi, bool processed) const { return (pass >= 1 && !processed) ? outC[gi] : inC[gi]; }

    static uint32_t downFlux(uint32_t u, uint32_t l, bool lProc, bool lExists) {
        if (!lExists || !lProc) return 0u;
        if ((FL(u) & SOLID) || (FL(l) & SOLID)) return 0u;
        if ((FL(l) & (SOURCE | SINK)) == SOURCE) return 0u;
        if ((FL(u) & SOURCE) && !(FL(l) & DUG)) return 0u;
        if (F(l) > 0u && F(u) > 0u && KIND(l) != KIND(u)) return 0u;
        return std::min(F(u), CAP(l) > F(l) ? CAP(l) - F(l) : 0u);
    }
    bool rests(int slot, int x, int y, int z, uint32_t c) const {
        bool proc; int face;
        const int gi = neighbour(slot, x, y, z, 0, 0, -1, proc, face);
        if (gi < 0) return false;
        const uint32_t b = readCell(gi, proc);
        if (FL(b) & SOLID) return true;
        return F(b) >= CAP(b) || (F(b) > 0u && KIND(b) != KIND(c));
    }
    bool onGround(int slot, int x, int y, int z) const {
        bool proc; int face;
        const int gi = neighbour(slot, x, y, z, 0, 0, -1, proc, face);
        return gi >= 0 && (FL(readCell(gi, proc)) & SOLID);
    }
    uint32_t kindFor(int slot, int x, int y, int z) const {
        uint32_t best = 0u, kind = 0u;
        for (int k = 0; k < 4; ++k) {
            const int dx = k == 0 ? -1 : k == 1 ? 1 : 0, dy = k == 2 ? -1 : k == 3 ? 1 : 0;
            bool proc; int face;
            const int gi = neighbour(slot, x, y, z, dx, dy, 0, proc, face);
            if (gi < 0) continue;
            const uint32_t c = readCell(gi, proc);
            if ((FL(c) & SOLID) || F(c) <= best) continue;
            best = F(c); kind = KIND(c);
        }
        return kind;
    }
    static bool sameWater(uint32_t c, uint32_t kind) { return !(FL(c) & SOLID) && F(c) > 0u && KIND(c) == kind; }
    int headOf(int gi, uint32_t c) const {
        int h = (int)(SH(c) + F(c));
        if (F(c) < CAP(c) || F(c) == 0u) return h;
        const uint32_t kind = KIND(c);
        int s, x, y, z;
        cellPos(gi, s, x, y, z);
        for (int k = 0; k < p.columnCells; ++k) {
            bool pr; int face;
            const int ga = neighbour(s, x, y, z, 0, 0, 1, pr, face);
            if (ga < 0) break;
            const uint32_t a = readCell(ga, pr);
            if (!sameWater(a, kind)) break;
            h += (int)(SH(a) + F(a));
            if (F(a) < CAP(a)) break;
            cellPos(ga, s, x, y, z);
        }
        return h;
    }
    bool submerged(int gi, uint32_t c) const {
        int s, x, y, z;
        cellPos(gi, s, x, y, z);
        bool pr; int face;
        const int ga = neighbour(s, x, y, z, 0, 0, 1, pr, face);
        return ga >= 0 && sameWater(readCell(ga, pr), KIND(c));
    }
    float columnDepth(int gi, uint32_t c) const {
        float d = (float)F(c) / p.hu;
        const uint32_t kind = KIND(c);
        int s, x, y, z;
        cellPos(gi, s, x, y, z);
        for (int k = 0; k < p.columnCells; ++k) {
            bool pr; int face;
            const int gb = neighbour(s, x, y, z, 0, 0, -1, pr, face);
            if (gb < 0) break;
            const uint32_t b = readCell(gb, pr);
            if (!sameWater(b, kind) || F(b) < CAP(b)) break;
            d += (float)F(b) / p.hu;
            cellPos(gb, s, x, y, z);
        }
        return d;
    }
    void faceFlow(int ga, uint32_t a, int gb, uint32_t b, float qp, float& qTot, float& qMom) const {
        qTot = 0.0f; qMom = 0.0f;
        if ((FL(a) & SOLID) || (FL(b) & SOLID)) return;
        if (F(a) >= CAP(a) && F(b) >= CAP(b)) return; // both full: the taker has no room either way (before the head walks)
        const int ha = headOf(ga, a), hb = headOf(gb, b);
        if (ha == hb) return;
        const bool aGives = ha > hb;
        const uint32_t d = aGives ? a : b, r = aGives ? b : a;
        const int gd = aGives ? ga : gb, gr = aGives ? gb : ga;
        if ((FL(r) & (SOURCE | SINK)) == SOURCE) return;
        if ((FL(d) & SOURCE) && !(FL(r) & DUG)) return;
        const uint32_t kd = KIND(d);
        const bool under = submerged(gd, d);
        const float colDepth = under ? (float)F(d) / p.hu : columnDepth(gd, d);
        if (colDepth * p.hu < (float)p.minSpread[kd]) return;
        int ds, dx, dy, dz, rs, rx, ry, rz;
        cellPos(gd, ds, dx, dy, dz); cellPos(gr, rs, rx, ry, rz);
        if (F(r) > 0u) { if (KIND(r) != kd) return; }
        else if (kindFor(rs, rx, ry, rz) != kd) return;
        if (!rests(ds, dx, dy, dz, d)) return;
        if (F(r) >= CAP(r)) return;
        const float sgn = aGives ? 1.0f : -1.0f;
        const int localDiff = (int)(SH(d) + F(d)) - (int)(SH(r) + F(r));
        float share = 0.0f;
        if (F(r) == 0u && !rests(rs, rx, ry, rz, r)) share = (float)F(d);
        else if (localDiff >= p.visc[kd]) share = (float)(localDiff / p.visc[kd]);
        const float depth = colDepth;
        const float dh = (float)std::abs(ha - hb) / p.hu;
        float q = std::max(qp * sgn, 0.0f);
        const float speed = q / (p.cellSize * std::max(depth, 0.02f)); // m/s through the face (a cell wide)
        const float accel = std::abs(ha - hb) >= p.visc[kd] ? p.dt * G * depth * dh : 0.0f;
        q = (q + accel) / (1.0f + p.dt * (p.frictionLin[kd] + p.frictionQuad[kd] * speed / std::max(depth, 0.02f)));
        if (accel == 0.0f && q * p.dt * p.units < 0.05f) q = 0.0f;
        qMom = q * sgn;
        qTot = std::max(q, share / (p.dt * p.units)) * sgn;
    }
    void outward(int gi, float amt[4]) const {
        int s, x, y, z;
        cellPos(gi, s, x, y, z);
        const Vec4& own = g.faceN[(size_t)2 * gi];
        amt[1] = own.x * p.dt * p.units + own.z;
        amt[3] = own.y * p.dt * p.units + own.w;
        bool pr; int face;
        const int gx = neighbour(s, x, y, z, -1, 0, 0, pr, face);
        amt[0] = (gx >= 0 && pr) ? -(g.faceN[(size_t)2 * gx].x * p.dt * p.units + g.faceN[(size_t)2 * gx].z) : 0.0f;
        const int gy = neighbour(s, x, y, z, 0, -1, 0, pr, face);
        amt[2] = (gy >= 0 && pr) ? -(g.faceN[(size_t)2 * gy].y * p.dt * p.units + g.faceN[(size_t)2 * gy].w) : 0.0f;
    }
    float giveScale(int gi, uint32_t c) const {
        float a[4]; outward(gi, a);
        const float total = std::max(a[0], 0.0f) + std::max(a[1], 0.0f) + std::max(a[2], 0.0f) + std::max(a[3], 0.0f);
        return total > (float)F(c) ? (float)F(c) / total : 1.0f;
    }
    uint32_t roomShare(int gi, uint32_t c) const {
        float a[4]; outward(gi, a);
        const uint32_t n = (uint32_t)(a[0] < 0.0f) + (uint32_t)(a[1] < 0.0f) + (uint32_t)(a[2] < 0.0f) + (uint32_t)(a[3] < 0.0f);
        const uint32_t room = CAP(c) > F(c) ? CAP(c) - F(c) : 0u;
        return n > 0u ? room / n : room;
    }

    void cell(int slot, int local) const {
        const int x = local % 16, y = (local / 16) % 16, z = local / 256;
        const int gi = slot * 4096 + local;
        const uint32_t c = inC[gi];
        const uint32_t fl = FL(c);
        if (pass == 1) {
            const Vec4 P = g.faceP[gi];
            Vec4 N{}, M{};
            for (int k = 0; k < 2; ++k) {
                bool pN; int face;
                const int gn = neighbour(slot, x, y, z, k == 0 ? 1 : 0, k == 1 ? 1 : 0, 0, pN, face);
                if (gn < 0 || !pN) continue;
                if (F(c) == 0u && F(inC[gn]) == 0u) continue; // two dry cells: nothing flows (#669)
                float q, qm;
                faceFlow(gi, c, gn, inC[gn], k == 0 ? P.x : P.y, q, qm);
                if (q == 0.0f) continue;
                if (k == 0) { N.x = q; N.z = Sign(P.z) == Sign(q) ? P.z : 0.0f; M.x = qm; }
                else { N.y = q; N.w = Sign(P.w) == Sign(q) ? P.w : 0.0f; M.y = qm; }
            }
            g.faceN[(size_t)2 * gi] = H(N);
            g.faceN[(size_t)2 * gi + 1] = H(M);
            return;
        }
        if (fl & SOLID) { outC[gi] = c; if (pass == 2) g.faceP[gi] = Vec4{}; return; }
        const int f = (int)F(c);
        uint32_t kind = KIND(c);
        if (pass == 0) {
            bool pB, pA; int fB, fA;
            const int gb = neighbour(slot, x, y, z, 0, 0, -1, pB, fB);
            const int ga = neighbour(slot, x, y, z, 0, 0, 1, pA, fA);
            const uint32_t below = gb >= 0 ? inC[gb] : 0u, above = ga >= 0 ? inC[ga] : 0u;
            const uint32_t outF = f > 0 ? downFlux(c, below, pB, gb >= 0) : 0u;
            const uint32_t inF = ga >= 0 ? downFlux(above, c, true, true) * (uint32_t)pA : 0u;
            if (gb < 0 && f > 0 && !(fl & SOURCE)) g.flags[slot] |= 1u << (1 + fB);
            const int nf = f - (int)outF + (int)inF;
            if (f == 0 && inF > 0u) kind = KIND(above);
            outC[gi] = Pack((uint32_t)std::clamp(nf, 0, 255), (fl & 31u) | (kind << 5), c);
            return;
        }
        float mine[4];
        outward(gi, mine);
        float myScale = -1.0f; // computed when a face first needs them (#669)
        int myRoom = -1;
        int net = 0;
        uint32_t moving = 0u, inKind = kind;
        int best = 0;
        Vec4 applied{};
        for (int k = 0; k < 4; ++k) {
            const int dx = k == 0 ? -1 : k == 1 ? 1 : 0, dy = k == 2 ? -1 : k == 3 ? 1 : 0;
            bool pN; int face;
            const int gn = neighbour(slot, x, y, z, dx, dy, 0, pN, face);
            if (gn < 0) {
                if (f >= p.minSpread[kind] && !(fl & SOURCE) && rests(slot, x, y, z, c)) g.flags[slot] |= 1u << (1 + face);
                continue;
            }
            const float a = mine[k];
            if (!pN || a == 0.0f) continue;
            const uint32_t n = readCell(gn, pN);
            if (a > 0.0f && myScale < 0.0f) myScale = giveScale(gi, c);
            if (a < 0.0f && myRoom < 0) myRoom = (int)roomShare(gi, c);
            const float sc = a > 0.0f ? myScale : giveScale(gn, n);
            const uint32_t share = a > 0.0f ? roomShare(gn, n) : (uint32_t)myRoom;
            const float mag = std::fabs(a) * sc;
            const int units = (int)std::min(std::floor(mag), (float)share);
            const bool clipped = (float)units < std::floor(mag);
            const float carry = clipped ? 0.0f : mag - std::floor(mag);
            const int signedUnits = a > 0.0f ? -units : units;
            net += signedUnits;
            moving += (uint32_t)units;
            if (signedUnits > best) { best = signedUnits; inKind = KIND(n); }
            const float dirAxis = (a > 0.0f ? 1.0f : -1.0f) * ((k == 1 || k == 3) ? 1.0f : -1.0f);
            if (k == 1 || k == 3) {
                const Vec4& M = g.faceN[(size_t)2 * gi + 1];
                const float qMom = k == 1 ? M.x : M.y;
                float flowAxis = qMom * sc;
                if (clipped) flowAxis = Sign(flowAxis) * std::min(std::fabs(flowAxis), (float)units / (p.dt * p.units));
                if (k == 1) { applied.x = flowAxis; applied.z = dirAxis * carry; }
                else { applied.y = flowAxis; applied.w = dirAxis * carry; }
            }
        }
        g.faceP[gi] = H(applied);
        int nf = std::clamp(f + net, 0, 255);
        if (f == 0 && nf > 0) kind = inKind;
        if (fl & SOURCE) nf = (int)TG(c);
        else if (p.evaporate && nf > 0 && nf <= 2 && moving == 0u && onGround(slot, x, y, z)) nf = 0;
        const uint32_t old = outC[gi];
        const uint32_t res = Pack((uint32_t)nf, (fl & 31u) | (kind << 5), c);
        if (F(old) != (uint32_t)nf) g.flags[slot] |= 1u;
        outC[gi] = res;
    }
};

// One pass over every listed slot, the slots split across the job system's workers (this thread takes the first share).
void RunPass(Grid& g, const Params& p, const std::vector<int>& list, int pass) {
    const uint32_t* in = pass == 0 ? g.a.data() : g.b.data();
    uint32_t* out = pass == 0 ? g.b.data() : g.a.data();
    const Pass ps{ g, p, pass, in, out };
    const int n = (int)list.size();
    if (n == 0) return;
    const int parts = std::max(1, std::min(n, (int)std::thread::hardware_concurrency()));
    auto run = [&ps, &list, n, parts](int part) {
        for (int i = part; i < n; i += parts)
            for (int local = 0; local < 4096; ++local) ps.cell(list[i], local);
    };
    std::atomic<int> left{ parts - 1 };
    for (int part = 1; part < parts; ++part) {
        JobSystem::tagNext(JobTag::Fluid);
        WorldJobs().submit(-3.0f, [&run, &left, part] { run(part); left.fetch_sub(1); });
    }
    run(0);
    while (left.load() > 0) std::this_thread::yield();
}

} // namespace

void Substep(Grid& g, const Params& p, const std::vector<int>& list) {
    RunPass(g, p, list, 0);
    RunPass(g, p, list, 1);
    RunPass(g, p, list, 2);
}

} // namespace FluidKernel
