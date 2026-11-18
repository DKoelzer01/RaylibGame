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
std::atomic<long long> dbgWhy[12]; // #680 diagnosis: why column pairs of deep runs did not move (DBWHY)
std::atomic<long long> dbgPair[3]; long long dbgPairPlan = 0, dbgPairUnits = 0; long long dbgPairFric = 0, dbgPairDH = 0, dbgPairH = 0, dbgPairN = 0; // TEMP #680 (DBWHY): deep pairs limited by plan / avail / room
#define WHY(i) do { if (A.surf - A.bed > 128) ++dbgWhy[i]; } while (0)

void Grid::resize(int slots) {
    a.assign((size_t)slots * 4096, 0u); b.assign((size_t)slots * 4096, 0u);
    faceP.assign((size_t)slots * 4096, Vec4{}); faceN.assign((size_t)slots * 4096 * 2, Vec4{});
    nbr.assign((size_t)slots * 6, -1); processed.assign(slots, 0); flags.assign(slots, 0u); delta.assign((size_t)slots * 4096, 0); sinkOut.assign(slots, 0u);
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
        const int gate = colDepth < 0.05f ? p.visc[kd] : 2; // films hold a slight slope; deeper water levels on 2 units (#680)
        const float accel = std::abs(ha - hb) >= gate ? p.dt * G * depth * dh : 0.0f;
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

    // ---- Column-pair flow (#680, FLUID_COLUMNS) ----------------------------------------------------------------------
    // Water moves sideways between columns, not cells: a run of water cells resting on support has a bed and a surface;
    // each pair of neighbouring runs carries one discharge Q (m^3/s, kept on the owning run's bottom cell) driven by the
    // difference of their surfaces over the whole wetted depth. Its units leave the giver's top cell and arrive at the
    // taker's head (atomically, into g.delta); the apply pass spills a head's overflow into the cell above and covers a
    // top's shortfall from the cell below. Per cell faces passed water only through the one layer at the downstream
    // surface (#680): a quarter of the depth at 0.25 m cells.
    bool waterOf(uint32_t c, uint32_t kind) const { return !(FL(c) & SOLID) && F(c) > 0u && KIND(c) == kind; }
    int stepZ(int gi, int dz, bool& proc) const {
        int s, x, y, z; cellPos(gi, s, x, y, z);
        int face;
        return neighbour(s, x, y, z, 0, 0, dz, proc, face);
    }
    // A run from its bottom cell up: its top cell, its surface (units over the reference row's floor), the cell under the
    // top, the head (where arriving water goes: the top if it has room, else the open cell over it) and the cell over the
    // head; ok when the cells it changes are in processed slots (the apply pass runs only there).
    struct ColRun { int bottom = -1, top = -1, below = -1, head = -1, above = -1; int bed = 0, surf = 0; bool ok = false; };
    ColRun runUp(int bottom, int rowOffset) const {
        ColRun r;
        const uint32_t kind = KIND(inC[bottom]);
        r.bottom = bottom; r.bed = rowOffset * 255 + (int)SH(inC[bottom]);
        int cur = bottom, dz = 0;
        bool proc = true;
        for (int k = 0; k < 64; ++k) {
            bool pr;
            const int up = stepZ(cur, 1, pr);
            if (F(inC[cur]) >= CAP(inC[cur]) && up >= 0 && waterOf(inC[up], kind)) { cur = up; ++dz; proc = proc && pr; continue; }
            break;
        }
        r.top = cur;
        r.surf = (rowOffset + dz) * 255 + (int)(SH(inC[cur]) + F(inC[cur]));
        bool pb = true, ph = true, pa = true;
        r.below = cur == bottom ? -1 : stepZ(cur, -1, pb);
        if (r.below >= 0 && !pb) r.below = -1; // its apply must run to give the shortfall: only in a processed slot
        r.head = cur;
        if (F(inC[cur]) >= CAP(inC[cur])) {
            const int up = stepZ(cur, 1, ph);
            r.head = (up >= 0 && ph && !(FL(inC[up]) & SOLID) && (F(inC[up]) == 0u || KIND(inC[up]) == kind)) ? up : -1;
        }
        r.above = r.head >= 0 ? stepZ(r.head, 1, pa) : -1;
        if (r.above >= 0 && (!pa || (FL(inC[r.above]) & SOLID) || (F(inC[r.above]) > 0u && KIND(inC[r.above]) != kind))) r.above = -1;
        int s, x, y, z;
        cellPos(cur, s, x, y, z);
        r.ok = proc && g.processed[s];
        return r;
    }
    // What a run can give this substep (its top cell and the full cell under it) and take (its head and the cell above).
    int runAvail(const ColRun& r) const { return (int)F(inC[r.top]) + (r.below >= 0 ? (int)F(inC[r.below]) : 0); }
    int runRoom(const ColRun& r) const {
        if (r.head < 0) return 0;
        return (int)(CAP(inC[r.head]) - F(inC[r.head])) + (r.above >= 0 ? (int)(CAP(inC[r.above]) - F(inC[r.above])) : 0);
    }
    // Plans the pair between this run (bottom cell gi) and the neighbour column on side (dx, dy): moves its units into
    // g.delta and returns the new Q (+: from this run toward the neighbour). A wet pair belongs to the -x / -y side's run.
    float colPair(int gi, int dx, int dy, float qPrev, bool owner) const {
        const uint32_t c = inC[gi];
        const uint32_t kind = KIND(c);
        const ColRun A = runUp(gi, 0);
        if (!A.ok) return 0.0f;
        int s, x, y, z;
        cellPos(gi, s, x, y, z);
        bool pn; int face;
        const int nb = neighbour(s, x, y, z, dx, dy, 0, pn, face);
        if (nb < 0) { WHY(0); if (face >= 0 && F(c) >= (uint32_t)p.minSpread[kind]) g.flags[s] |= 1u << (1 + face); return 0.0f; }
        // A sleeping neighbour brick: ask for it to be woken (the face rules wake theirs; without this a hole cut into a wall
        // or a channel beyond it stayed asleep and the outflow started late, #680).
        if (!pn) { WHY(1); if (face >= 0 && F(c) >= (uint32_t)p.minSpread[kind]) g.flags[s] |= 1u << (1 + face); return 0.0f; }
        const int topRows = (A.surf + 254) / 255; // rows A's water reaches over its bottom row
        ColRun B;
        bool wet = false, spill = false;
        int takerSpill = -1;
        const uint32_t n0 = inC[nb];
        if (FL(n0) & SOLID) {
            int cur = nb, k = 0;
            bool pr = true;
            while (k < topRows) {
                const int up = stepZ(cur, 1, pr);
                if (up < 0 || !pr) return 0.0f;
                cur = up; ++k;
                if (!(FL(inC[cur]) & SOLID)) break;
            }
            if (FL(inC[cur]) & SOLID) { WHY(2); return 0.0f; } // a wall over A's water
            if (waterOf(inC[cur], kind)) { B = runUp(cur, k); wet = true; }
            else if (F(inC[cur]) > 0u) return 0.0f;          // another fluid
            else {                                           // a dry floor higher than A's bed
                B.bottom = B.top = B.head = cur;
                B.bed = B.surf = k * 255 + (int)SH(inC[cur]);
                bool pa;
                B.above = stepZ(cur, 1, pa);
                if (B.above >= 0 && (!pa || (FL(inC[B.above]) & SOLID) || F(inC[B.above]) > 0u)) B.above = -1;
                B.ok = true;
            }
        } else if (waterOf(n0, kind)) {
            int cur = nb, k = 0;
            bool pr = true;
            while (k > -64) {                                // down to the run's bottom
                const int dn = stepZ(cur, -1, pr);
                if (dn < 0 || !pr || !waterOf(inC[dn], kind) || F(inC[dn]) < CAP(inC[dn])) break;
                cur = dn; --k;
            }
            B = runUp(cur, k); wet = true;
        } else if (F(n0) > 0u) {
            return 0.0f;                                     // another fluid
        } else {
            spill = true; takerSpill = nb;                   // air beside A's bottom: over the edge (it falls from there)
        }
        if (wet && !owner) { WHY(3); return 0.0f; }
        // A pair that just turned wet was the other run's while this side was dry (a spill toward -x / -y, kept in its z / w):
        // its momentum carries over (that side stores nothing now), or every front moving toward -x restarted from rest.
        if (wet && owner) qPrev -= dx == 1 ? g.faceP[B.bottom].z : g.faceP[B.bottom].w;
        if (!spill && !B.ok) { WHY(4); return 0.0f; }
        // Surfaces and the flow (+: from A toward the neighbour).
        const int bedMax = spill ? A.bed : std::max(A.bed, B.bed);
        const int surfB = spill ? A.bed : B.surf;
        const bool aUp = A.surf >= surfB;
        const int surfU = aUp ? A.surf : surfB, surfD = aUp ? surfB : A.surf;
        const int h = surfU - bedMax, dH = surfU - std::max(surfD, bedMax);
        const float hm = (float)h / p.hu, dHm = (float)dH / p.hu;
        float q = qPrev;
        if (h >= p.minSpread[kind] && dH >= 2) q += p.dt * G * hm * dHm * (aUp ? 1.0f : -1.0f);
        if (!wet && q < 0.0f) q = 0.0f;                      // a dry or lower neighbour never gives back
        if (h < p.minSpread[kind]) q = 0.0f;                 // films rest
        // Momentum is carried with the water (upwind advection, -u dQ/dx): a pair's flow eases toward the flow of the pair
        // upstream of it at u / cs a second. Without it water arriving in a column brought no momentum into the next pair,
        // which had to be pushed from rest by a rising surface: the front piled into a mound and crawled (#680).
        if (p.advect && h >= p.minSpread[kind]) { // also from rest: a front is pushed by the water behind it
            const int k = dx != 0 ? 0 : 1, dir = dx + dy;  // axis; +1 toward +x / +y
            float qUp = 0.0f; bool have = false;
            const float way = q != 0.0f ? q : (aUp ? 1.0f : -1.0f); // the flow's direction (+: A gives)
            if (way > 0.0f) {                                // A gives: the pair behind A (on its far side)
                if (dir > 0) {
                    bool pl; int fl2;
                    int lb = neighbour(s, x, y, z, -dx, -dy, 0, pl, fl2);
                    if (lb >= 0 && waterOf(inC[lb], kind)) {
                        for (int i = 0; i < 64; ++i) { bool pd; const int dn = stepZ(lb, -1, pd); if (dn < 0 || !waterOf(inC[dn], kind) || F(inC[dn]) < CAP(inC[dn])) break; lb = dn; }
                        qUp = k == 0 ? g.faceP[lb].x : g.faceP[lb].y; have = true;
                    } else { qUp = -(k == 0 ? g.faceP[gi].z : g.faceP[gi].w); have = qUp != 0.0f; }
                } else { qUp = -(k == 0 ? g.faceP[gi].x : g.faceP[gi].y); have = true; } // toward -x: A's own +x pair, reversed
            } else if (wet && dir > 0) {                      // the neighbour gives: the pair beyond it
                qUp = k == 0 ? g.faceP[B.bottom].x : g.faceP[B.bottom].y; have = true;
            }
            if (have && qUp * way > 0.0f) {
                // Uncapped: with a full reservoir the 0.25 m front runs at Ritter's speed (1.03x to 10 m); capping u by
                // sqrt(g h) at the pair's or the giver's depth stalled it (#680; the earlier 1.78x came from a half-filled pool).
                const float u = 0.5f * (std::fabs(q) + std::fabs(qUp)) / (p.cellSize * std::max(hm, 0.02f));
                q += std::min(1.0f, p.dt * u / p.cellSize) * (qUp - q);
            }
        }
        const float speed = std::fabs(q) / (p.cellSize * std::max(hm, 0.02f));
        const float fricDiv = 1.0f + p.dt * (p.frictionLin[kind] + p.frictionQuad[kind] * speed / std::max(hm, 0.02f));
        q /= fricDiv;
        if (A.surf - A.bed > 128) { // TEMP #680 (DBWHY): per deep pair the friction divisor, the head drop and the depth (x1000)
            __atomic_fetch_add(&dbgPairFric, (long long)(fricDiv * 1000.0f), __ATOMIC_RELAXED);
            __atomic_fetch_add(&dbgPairDH, (long long)(dHm * 1000.0f), __ATOMIC_RELAXED);
            __atomic_fetch_add(&dbgPairH, (long long)(hm * 1000.0f), __ATOMIC_RELAXED);
            __atomic_fetch_add(&dbgPairN, 1LL, __ATOMIC_RELAXED);
        }
        const float plan = std::fabs(q) * p.dt * p.units;
        if (plan < 1.0f) { WHY(5); return q; }
        const bool aGives = q > 0.0f;
        const int giverTop = aGives ? A.top : B.top;
        const int avail = aGives ? runAvail(A) : runAvail(B);
        int taker = -1, room = 0;
        if (aGives) {
            if (spill) { taker = takerSpill; room = (int)(CAP(inC[taker]) - F(inC[taker])); }
            else { taker = B.head; room = runRoom(B); }
        } else { taker = A.head; room = runRoom(A); }
        if (taker < 0) { WHY(6); return 0.0f; }
        // Sources (generated water, held at their level each step) as in the per-face rules: one never takes, and gives
        // only into dug cells (infinite water reaches a dig, not the open ground); a sink takes.
        if ((FL(inC[taker]) & (SOURCE | SINK)) == SOURCE) return 0.0f;
        if ((FL(inC[giverTop]) & SOURCE) && !(FL(inC[taker]) & DUG)) return 0.0f;
        const int units = (int)std::floor(std::min({ plan, (float)avail / 4.0f, (float)room / 4.0f }));
        if (A.surf - A.bed > 128) { // TEMP #680 (DBPAIR): which limit set a deep pair's units this substep
            const float lim = std::min({ plan, (float)avail / 4.0f, (float)room / 4.0f });
            ++dbgPair[lim == plan ? 0 : lim == (float)avail / 4.0f ? 1 : 2];
            __atomic_fetch_add(&dbgPairPlan, (long long)plan, __ATOMIC_RELAXED); __atomic_fetch_add(&dbgPairUnits, (long long)units, __ATOMIC_RELAXED);
        }
        if (units <= 0) { WHY(7); return q * 0.5f; }
        WHY(spill ? 8 : wet ? 9 : 10);
        if ((float)units < std::floor(plan)) q = (aGives ? 1.0f : -1.0f) * (float)units / (p.dt * p.units);
        __atomic_fetch_sub(&g.delta[giverTop], units, __ATOMIC_RELAXED);
        __atomic_fetch_add(&g.delta[taker], units, __ATOMIC_RELAXED);
        return q;
    }
    // A cell's new fill from the column flows: its own delta (clamped), a head's overflow from the cell below, a top's
    // shortfall taken from it when it lies under one.
    int colApply(int gi, uint32_t c, uint32_t& kindOut) const {
        const int f = (int)F(c), cap = (int)CAP(c);
        int nf = std::clamp(f + g.delta[gi], 0, cap);
        bool pb, pa;
        const int bl = stepZ(gi, -1, pb), ab = stepZ(gi, 1, pa);
        if (bl >= 0) {
            const uint32_t b = inC[bl];
            const int over = (int)F(b) + g.delta[bl] - (int)CAP(b);
            if (over > 0 && !(FL(b) & SOLID)) { if (nf == 0) kindOut = KIND(b); nf += over; }
        }
        if (ab >= 0) {
            const uint32_t a = inC[ab];
            const int under = (int)F(a) + g.delta[ab];
            if (under < 0 && !(FL(a) & SOLID)) nf += under;
        }
        return std::clamp(nf, 0, cap);
    }
    void cell(int slot, int local) const {
        const int x = local % 16, y = (local / 16) % 16, z = local / 256;
        const int gi = slot * 4096 + local;
        const uint32_t c = inC[gi];
        const uint32_t fl = FL(c);
        if (pass == 1 && p.columns) {
            // Column-pair flow (#680): a run's bottom cell plans its pairs (+x, +y always; -x, -y only toward a dry or lower
            // column: a wet pair belongs to the -side run); Q per side in the face record (x +x, y +y, z -x, w -y).
            Vec4 Q{};
            bool pb; int fb;
            const int gb = neighbour(slot, x, y, z, 0, 0, -1, pb, fb);
            const bool bottom = !(fl & SOLID) && F(c) > 0u && !(fl & SOURCE) && gb >= 0 &&
                                ((FL(inC[gb]) & SOLID) || (F(inC[gb]) > 0u && KIND(inC[gb]) != KIND(c)));
            if (bottom) {
                const Vec4 P = g.faceP[gi];
                Q.x = colPair(gi, 1, 0, P.x, true);
                Q.y = colPair(gi, 0, 1, P.y, true);
                Q.z = colPair(gi, -1, 0, P.z, false);
                Q.w = colPair(gi, 0, -1, P.w, false);
                if (Q.x != 0.0f || Q.y != 0.0f || Q.z != 0.0f || Q.w != 0.0f) g.flags[slot] |= 1u;
            }
            g.faceN[(size_t)2 * gi] = Q;
            return;
        }
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
            g.delta[gi] = 0; // the column flows planned next (#680)
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
        if (p.columns) {
            // Column-pair flow (#680): the planned units, a head's overflow and a top's shortfall.
            g.faceP[gi] = H(g.faceN[(size_t)2 * gi]);
            uint32_t inKind = kind;
            int nf = colApply(gi, c, inKind);
            if (f == 0 && nf > 0) kind = inKind;
            if (fl & SOURCE) nf = (int)TG(c);
            else if (p.evaporate && nf > 0 && nf <= 2 && g.delta[gi] == 0 && onGround(slot, x, y, z)) nf = 0;
            if (f != nf) g.flags[slot] |= 1u;
            outC[gi] = Pack((uint32_t)nf, (fl & 31u) | (kind << 5), c);
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
        if ((fl & (SOURCE | SINK)) == (SOURCE | SINK) && nf > (int)TG(c)) g.sinkOut[slot] += (uint32_t)(nf - (int)TG(c)); // a released reach's sink takes in (#656)
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
    static const bool why = EnvFlag("DBWHY");
    static int nWhy = 0;
    if (why && ++nWhy % (int)(60.0f / p.cellSize) == 0) {
        printf("dbwhy:");
        static const char* names[11] = { "noBrick", "unproc", "wall", "notOwner", "Bnotok", "plan<1", "noTaker", "units0", "SPILL", "WET", "DRY" };
        for (int i = 0; i < 11; ++i) { printf(" %s %lld", names[i], dbgWhy[i].load()); dbgWhy[i] = 0; }
        printf("\ndbpair: limited by plan %lld, avail/4 %lld, room/4 %lld; planned %lld units, moved %lld\n", dbgPair[0].load(), dbgPair[1].load(), dbgPair[2].load(),
               __atomic_load_n(&dbgPairPlan, __ATOMIC_RELAXED), __atomic_load_n(&dbgPairUnits, __ATOMIC_RELAXED));
        for (auto& d : dbgPair) d = 0;
        { const long long n = std::max(1LL, __atomic_load_n(&dbgPairN, __ATOMIC_RELAXED)); printf("dbpair: mean friction divisor %.4f, head drop %.4f m, depth %.3f m over %lld pair-substeps; dt %.5f s\n", __atomic_load_n(&dbgPairFric, __ATOMIC_RELAXED) / 1000.0 / n, __atomic_load_n(&dbgPairDH, __ATOMIC_RELAXED) / 1000.0 / n, __atomic_load_n(&dbgPairH, __ATOMIC_RELAXED) / 1000.0 / n, n, p.dt); __atomic_store_n(&dbgPairFric, 0, __ATOMIC_RELAXED); __atomic_store_n(&dbgPairDH, 0, __ATOMIC_RELAXED); __atomic_store_n(&dbgPairH, 0, __ATOMIC_RELAXED); __atomic_store_n(&dbgPairN, 0, __ATOMIC_RELAXED); }
        __atomic_store_n(&dbgPairPlan, 0, __ATOMIC_RELAXED); __atomic_store_n(&dbgPairUnits, 0, __ATOMIC_RELAXED);
    }
}

} // namespace FluidKernel
