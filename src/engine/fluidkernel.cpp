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

constexpr float kMatchPlan = 3.0f; // column pairs planning fewer units a substep carry only the momentum they moved (#855)

void Grid::resize(int slots) {
    a.assign((size_t)slots * 4096, 0u); b.assign((size_t)slots * 4096, 0u);
    faceP.assign((size_t)slots * 4096, Vec4{}); faceN.assign((size_t)slots * 4096 * 2, Vec4{});
    nbr.assign((size_t)slots * 6, -1); processed.assign(slots, 0); flags.assign(slots, 0u); delta.assign((size_t)slots * 4096, 0); sinkOut.assign(slots, 0u);
    dem.assign((size_t)slots * 4096, 0u); cap.assign((size_t)slots * 4096, 0u);
    pairCode.assign((size_t)slots * 4096 * 4, (uint16_t)0);
    gmask.assign((size_t)slots * 6, 0u); gon.assign(slots, ~0ull); groupsFull = true;
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
        // Inside the brick (#747: the column passes step up and down constantly): the cell 256 on, the slot's own flag.
        const int zl = (gi & 4095) >> 8;
        if ((dz == 1 && zl < 15) || (dz == -1 && zl > 0)) { proc = g.processed[gi >> 12] != 0; return gi + dz * 256; }
        int s, x, y, z; cellPos(gi, s, x, y, z);
        int face;
        return neighbour(s, x, y, z, 0, 0, dz, proc, face);
    }
    // A run from its bottom cell up: its top cell, its surface (units over the reference row's floor), the cell under the
    // top and the head (where arriving water goes: the top if it has room, else the open cell over it; headRoom counts the
    // room over it too), with their rows; ok when the cells it changes are in processed slots (the apply pass runs only there).
    // topRow / headRow: their rows over the reference row (the move pass walks back to them, #680 GPU layout).
    struct ColRun { int bottom = -1, top = -1, below = -1, head = -1; int bed = 0, surf = 0; int topRow = 0, headRow = 0; bool ok = false; };
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
        r.topRow = r.headRow = rowOffset + dz;
        r.surf = (rowOffset + dz) * 255 + (int)(SH(inC[cur]) + F(inC[cur]));
        bool pb = true, ph = true;
        r.below = cur == bottom ? -1 : stepZ(cur, -1, pb);
        if (r.below >= 0 && !pb) r.below = -1; // its apply must run to give the shortfall: only in a processed slot
        r.head = cur;
        if (F(inC[cur]) >= CAP(inC[cur])) {
            const int up = stepZ(cur, 1, ph); // a missing or sleeping brick over it: no head (askAbove, when water comes)
            r.head = (up >= 0 && ph && !(FL(inC[up]) & SOLID) && (F(inC[up]) == 0u || KIND(inC[up]) == kind)) ? up : -1;
            r.headRow = r.topRow + 1;
        }
        int s, x, y, z;
        cellPos(cur, s, x, y, z);
        r.ok = proc && g.processed[s];
        return r;
    }
    // What a run can give this substep (its top cell and the full cell under it) and take (its head and the cell above).
    int runAvail(const ColRun& r) const { return (int)F(inC[r.top]) + (r.below >= 0 ? (int)F(inC[r.below]) : 0); }
    // What a taker cell can take this substep: its room and the room of the open cell over it (the apply lifts its
    // overflow there). One value a cell, whoever asks - a run's head, a dry floor, a spill (#680 GPU layout: the room is
    // OR-ed into one word a cell; a spill once counted only its own cell's room).
    int headRoom(int cell, uint32_t kind) const {
        if (cell < 0) return 0;
        int room = (int)(CAP(inC[cell]) - F(inC[cell]));
        bool pa;
        const int up = stepZ(cell, 1, pa);
        if (up >= 0 && pa && !(FL(inC[up]) & SOLID) && (F(inC[up]) == 0u || KIND(inC[up]) == kind)) room += (int)(CAP(inC[up]) - F(inC[up]));
        return room;
    }
    int runRoom(const ColRun& r, uint32_t kind) const { return r.head < 0 ? 0 : headRoom(r.head, kind); }
    // A cell `row` rows over (under, if negative) a column cell, or -1 (the move pass, from a pair's code).
    int cellAtRow(int cell, int row) const {
        bool pr;
        for (; row > 0 && cell >= 0; --row) cell = stepZ(cell, 1, pr);
        for (; row < 0 && cell >= 0; ++row) cell = stepZ(cell, -1, pr);
        return cell;
    }
    // Whether water would move between run A and a sleeping neighbour's column B (#855): a spill over A's edge (A as deep as
    // a film that spreads), or surfaces apart by wakeUnits or more over the higher bed.
    bool wakesNeighbour(const ColRun& A, const ColRun& B, bool spill, uint32_t kind) const {
        if (spill) return A.surf - A.bed >= p.minSpread[kind];
        if (B.bottom < 0) return false;
        const int bedMax = std::max(A.bed, B.bed);
        const int surfU = std::max(A.surf, B.surf), surfD = std::min(A.surf, B.surf);
        return surfU - bedMax >= p.minSpread[kind] && surfU - std::max(surfD, bedMax) >= p.wakeUnits;
    }
    // A run whose top is full and has no head (the brick over it is missing or asleep) asks for that brick when a pair wants
    // to put water in (#680: a run could not rise past its brick's top - the dam-break's channel held 0.3 m past x = 5 m while
    // a 1.7 m bore stood behind it). Only then (#855: asking from every full top kept the brick over a still pond awake).
    // Atomic: the run may lie in a slot another thread is planning.
    void askAbove(const ColRun& r) const {
        if (r.top < 0 || F(inC[r.top]) < CAP(inC[r.top])) return;
        int s2, x2, y2, z2; cellPos(r.top, s2, x2, y2, z2);
        bool pu; int fu;
        neighbour(s2, x2, y2, z2, 0, 0, 1, pu, fu);
        if (fu >= 0 && !pu) __atomic_fetch_or(&g.flags[s2], 1u << (1 + fu), __ATOMIC_RELAXED);
    }
    // Plans the pair between this run (bottom cell gi) and the neighbour column on side (dx, dy): moves its units into
    // g.delta and returns the new Q (+: from this run toward the neighbour). A wet pair belongs to the -x / -y side's run.
    // k: the pair's slot in the run's record (0 +x, 1 +y, 2 -x, 3 -y); planned: the units it plans to move (0: none).
    // A: this run (runUp(gi, 0)), walked once for its four pairs (#747: each pair walked it again).
    float colPair(int gi, const ColRun& A, int dx, int dy, float qPrev, bool owner, int k, float& planned) const {
        planned = 0.0f;
        const uint32_t c = inC[gi];
        const uint32_t kind = KIND(c);
        if (!A.ok) return 0.0f;
        int s, x, y, z;
        cellPos(gi, s, x, y, z);
        bool pn; int face;
        const int nb = neighbour(s, x, y, z, dx, dy, 0, pn, face);
        if (nb < 0) { if (face >= 0 && F(c) >= (uint32_t)p.minSpread[kind]) g.flags[s] |= 1u << (1 + face); return 0.0f; }
        // A sleeping neighbour brick (its cells are packed, read only): woken when water would move between the two (without
        // that a hole cut into a wall or a channel beyond it stayed asleep and the outflow started late, #680), not by resting
        // water beside it - asking every substep woke every neighbour of an awake brick and no pond ever slept (#855).
        const bool sleepNb = !pn;
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
                if (up < 0 || (!pr && !sleepNb)) return 0.0f;
                cur = up; ++k;
                if (!(FL(inC[cur]) & SOLID)) break;
            }
            if (FL(inC[cur]) & SOLID) { return 0.0f; } // a wall over A's water
            if (waterOf(inC[cur], kind)) { B = runUp(cur, k); wet = true; }
            else if (F(inC[cur]) > 0u) return 0.0f;          // another fluid
            else {                                           // a dry floor higher than A's bed
                B.bottom = B.top = B.head = cur;
                B.bed = B.surf = k * 255 + (int)SH(inC[cur]);
                B.topRow = B.headRow = k;
                B.ok = true;
            }
        } else if (waterOf(n0, kind)) {
            int cur = nb, k = 0;
            bool pr = true;
            while (k > -64) {                                // down to the run's bottom
                const int dn = stepZ(cur, -1, pr);
                if (dn < 0 || (!pr && !sleepNb) || !waterOf(inC[dn], kind) || F(inC[dn]) < CAP(inC[dn])) break;
                cur = dn; --k;
            }
            B = runUp(cur, k); wet = true;
        } else if (F(n0) > 0u) {
            return 0.0f;                                     // another fluid
        } else if (SH(n0) > 0u) {
            // An empty cell with ground in it (a wall's or a floor's partly solid cell): a dry floor at its ground's height, not a
            // drop - spilling into it ignored the ground, and water chained through a wall's partial cells (#680: the dam-break's
            // wall leaked 1.1 m^3 and the leak test's 2.3 m^3 during the fill).
            B.bottom = B.top = B.head = nb;
            B.bed = B.surf = (int)SH(n0);
            B.topRow = B.headRow = 0;
            B.ok = true;
        } else {
            spill = true; takerSpill = nb;                   // air beside A's bottom: over the edge (it falls from there)
        }
        if (sleepNb) {
            if (face >= 0 && wakesNeighbour(A, B, spill, kind)) g.flags[s] |= 1u << (1 + face);
            return 0.0f;
        }
        if (wet && !owner) { return 0.0f; }
        // A pair that just turned wet was the other run's while this side was dry (a spill toward -x / -y, kept in its z / w):
        // its momentum carries over (that side stores nothing now), or every front moving toward -x restarted from rest.
        if (wet && owner) qPrev -= dx == 1 ? g.faceP[B.bottom].z : g.faceP[B.bottom].w;
        if (!spill && !B.ok) { return 0.0f; }
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
                bool pl; int fl2;
                int lb = neighbour(s, x, y, z, -dx, -dy, 0, pl, fl2);
                if (dir > 0) {
                    if (lb >= 0 && waterOf(inC[lb], kind)) {
                        for (int i = 0; i < 64; ++i) { bool pd; const int dn = stepZ(lb, -1, pd); if (dn < 0 || !waterOf(inC[dn], kind) || F(inC[dn]) < CAP(inC[dn])) break; lb = dn; }
                        qUp = k == 0 ? g.faceP[lb].x : g.faceP[lb].y; have = true;
                    } else { qUp = -(k == 0 ? g.faceP[gi].z : g.faceP[gi].w); have = qUp != 0.0f; }
                } else { qUp = -(k == 0 ? g.faceP[gi].x : g.faceP[gi].y); have = true; } // toward -x: A's own +x pair, reversed
            } else if (wet && dir > 0) {                      // the neighbour gives: the pair beyond it
                qUp = k == 0 ? g.faceP[B.bottom].x : g.faceP[B.bottom].y; have = true;
            }
            if (have && qUp * way > 0.0f) {
                // The water carries its velocity, not its discharge (#680): blending a thin pair toward the deep upstream
                // pair's discharge gave the front u = Q / (cs h) -> huge as h -> 0 (a dam break's sheet ran 16-20 m/s, ~2x
                // Ritter's 2 sqrt(g h0), the frictionless maximum). The target is the upstream velocity at this pair's depth,
                // never more than its discharge, the giver run's depth standing in for the upstream pair's (its true depth,
                // the deeper of its two runs, left a thin tip no target: the front crawled). With one substep scale the
                // dam-break's front runs 0.86x Ritter. FLUID_ADVECT_Q=1: the discharge blend (A/B).
                static const bool advectQ = EnvFlag("FLUID_ADVECT_Q");
                float qT = qUp;
                if (!advectQ) {
                    const int gDepth = way > 0.0f ? A.surf - A.bed : B.surf - B.bed; // the giver run's depth (units)
                    qT *= std::min(1.0f, hm / std::max((float)gDepth / p.hu, 0.02f));
                }
                const float u = 0.5f * (std::fabs(q) + std::fabs(qT)) / (p.cellSize * std::max(hm, 0.02f));
                q += std::min(1.0f, p.dt * u / p.cellSize) * (qT - q);
            }
        }
        const float speed = std::fabs(q) / (p.cellSize * std::max(hm, 0.02f));
        const float fricDiv = 1.0f + p.dt * (p.frictionLin[kind] + p.frictionQuad[kind] * speed / std::max(hm, 0.02f));
        q /= fricDiv;
        const float plan = std::fabs(q) * p.dt * p.units;
        if (plan < 1.0f) { return q; }
        const bool aGives = q > 0.0f;
        const int giverTop = aGives ? A.top : B.top;
        const int avail = aGives ? runAvail(A) : runAvail(B);
        int taker = -1, room = 0, giverRow = aGives ? A.topRow : B.topRow, takerRow = 0;
        if (aGives) {
            if (spill) { taker = takerSpill; room = headRoom(taker, kind); takerRow = 0; }
            else { taker = B.head; room = runRoom(B, kind); takerRow = B.headRow; }
        } else { taker = A.head; room = runRoom(A, kind); takerRow = A.headRow; }
        if (taker < 0) { askAbove(aGives ? B : A); return 0.0f; }
        // Sources (generated water, held at their level each step) as in the per-face rules: one never takes, and gives
        // only into dug cells (infinite water reaches a dig, not the open ground); a sink takes.
        if ((FL(inC[taker]) & (SOURCE | SINK)) == SOURCE) return 0.0f;
        if ((FL(inC[giverTop]) & SOURCE) && !(FL(inC[taker]) & DUG)) return 0.0f;
        // How much this pair may move. Each of a run's up to four pairs once took at most a quarter of what the run could
        // give or take, so none could overdraw it - but a thin run's quarter (one partial cell: ~33 units a substep) was a
        // third of a dam break's discharge, and a bore stood behind the choke and backed upstream (#680). Now a pair plans
        // up to the whole avail / room, adds it to its giver's and taker's demand, and the move pass scales every pair by
        // min(1, avail / out, room / in): exact conservation, any order (integer atomics), no choke.
        // FLUID_COL_QUARTER=1: the old quarter limit, moved here (A/B).
        static const bool quarter = EnvFlag("FLUID_COL_QUARTER");
        const int units = (int)std::floor(std::min({ plan, (float)avail / (quarter ? 4.0f : 1.0f), (float)room / (quarter ? 4.0f : 1.0f) }));
        if (units <= 0) { return q * 0.5f; }
        // A slow, nearly level pair's momentum is what it moved (#855): below kMatchPlan units a substep and 2 x wakeUnits of
        // head (3 cm), keeping the fraction floor() left behind carried momentum the water never had - it stuck, then slipped
        // a unit too far - and a basin at rest rocked by a few cm for good (no pond slept). Other pairs keep it: a thin front
        // moves a few units a substep (dropping up to one every substep stalled the dam-break's front at 5 m), and a stream on
        // a slope plans as few (dropping it there, the flow meter read a 1 m^3/s stream as 1.4).
        if ((float)units < std::floor(plan) || ((float)units < plan && plan < kMatchPlan && dH < 2 * p.wakeUnits)) q = (aGives ? 1.0f : -1.0f) * (float)units / (p.dt * p.units);
        if (quarter) {
            __atomic_fetch_sub(&g.delta[giverTop], units, __ATOMIC_RELAXED);
            __atomic_fetch_add(&g.delta[taker], units, __ATOMIC_RELAXED);
            return q;
        }
        const int rowDiff = takerRow - giverRow;
        if (giverRow < -64 || giverRow > 191 || rowDiff < -110 || rowDiff > 17) return q; // out of the 16-bit code's reach: no move
        planned = (float)units;
        g.pairCode[(size_t)4 * gi + k] = (uint16_t)((aGives ? 1u : 0u) | ((uint32_t)(giverRow + 64) << 1) | ((uint32_t)(rowDiff + 110) << 9));
        __atomic_fetch_add(&g.dem[giverTop], (uint32_t)units << 16, __ATOMIC_RELAXED); // out
        __atomic_fetch_add(&g.dem[taker], (uint32_t)units, __ATOMIC_RELAXED);         // in
        __atomic_fetch_or(&g.cap[giverTop], (uint32_t)avail << 16, __ATOMIC_RELAXED);  // one avail a cell (its run's)
        __atomic_fetch_or(&g.cap[taker], (uint32_t)room, __ATOMIC_RELAXED);            // one room a cell (headRoom)
        return q;
    }
    // A cell's new fill from the column flows: its own delta (clamped), a head's overflow from the cell below, a top's
    // shortfall taken from it when it lies under one.
    int colApply(int gi, uint32_t c, uint32_t& kindOut) const {
        // A cell holding more than its room before this substep keeps it (the GPU twin can find one: a CPU change placed by
        // a mirror a tick old; clamping it deleted the excess). The overflow passed up counts only what this substep's flows
        // push past that bound. On the CPU no cell is over its room, so this changes nothing here.
        const int f = (int)F(c), cap = std::max((int)CAP(c), (int)F(c));
        int nf = std::clamp(f + g.delta[gi], 0, cap);
        bool pb, pa;
        const int bl = stepZ(gi, -1, pb), ab = stepZ(gi, 1, pa);
        // Only processed cells move water (an unprocessed brick's scratch holds whatever it last had: on the GPU slots persist).
        if (bl >= 0 && pb) {
            const uint32_t b = inC[bl];
            const int over = (int)F(b) + g.delta[bl] - std::max((int)CAP(b), (int)F(b));
            if (over > 0 && !(FL(b) & SOLID)) { if (nf == 0) kindOut = KIND(b); nf += over; }
        }
        if (ab >= 0 && pa) {
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
            Vec4 Q{}, U{}; // U: the units each pair plans (the move pass scales and moves them)
            bool pb; int fb;
            const int gb = neighbour(slot, x, y, z, 0, 0, -1, pb, fb);
            // A source's run plans its pairs too (#982: a column of generated water never gave sideways, so a cave dug out to a
            // river stayed dry); the pair's source rules keep it to giving into dug cells only, and it holds its level.
            const bool bottom = !(fl & SOLID) && F(c) > 0u && gb >= 0 &&
                                ((FL(inC[gb]) & SOLID) || (F(inC[gb]) > 0u && KIND(inC[gb]) != KIND(c)));
            if (bottom) {
                const Vec4 P = g.faceP[gi];
                const ColRun A = runUp(gi, 0);
                Q.x = colPair(gi, A, 1, 0, P.x, true, 0, U.x);
                Q.y = colPair(gi, A, 0, 1, P.y, true, 1, U.y);
                Q.z = colPair(gi, A, -1, 0, P.z, false, 2, U.z);
                Q.w = colPair(gi, A, 0, -1, P.w, false, 3, U.w);
                if (Q.x != 0.0f || Q.y != 0.0f || Q.z != 0.0f || Q.w != 0.0f) g.flags[slot] |= 1u;
            }
            g.faceN[(size_t)2 * gi] = H(Q); // in half precision, as the GPU stores it (#680 GPU layout)
            g.faceN[(size_t)2 * gi + 1] = U;
            return;
        }
        if (pass == 3) {
            // The column pairs' move (#680): each planned pair moves its units scaled by its giver's and taker's demand, and
            // its flow follows what moved. Floors keep every cell's sum within its avail / room.
            const Vec4 U = g.faceN[(size_t)2 * gi + 1];
            if (U.x == 0.0f && U.y == 0.0f && U.z == 0.0f && U.w == 0.0f) return;
            Vec4& Q = g.faceN[(size_t)2 * gi];
            const float u4[4] = { U.x, U.y, U.z, U.w };
            float* q4[4] = { &Q.x, &Q.y, &Q.z, &Q.w };
            for (int k = 0; k < 4; ++k) {
                const int planned = (int)u4[k];
                if (planned <= 0) continue;
                // Back to the cells from the pair's code: rows from this run's bottom, in A's column or the neighbour's.
                const uint32_t code = g.pairCode[(size_t)4 * gi + k];
                const bool aGives = (code & 1u) != 0u;
                const int giverRow = (int)((code >> 1) & 255u) - 64, takerRow = giverRow + (int)((code >> 9) & 127u) - 110;
                bool pn; int face;
                const int nb = neighbour(slot, x, y, z, k == 0 ? 1 : k == 2 ? -1 : 0, k == 1 ? 1 : k == 3 ? -1 : 0, 0, pn, face);
                const int giver = cellAtRow(aGives ? gi : nb, giverRow), taker = cellAtRow(aGives ? nb : gi, takerRow);
                if (giver < 0 || taker < 0) continue; // can't happen: the plan walked these cells
                const uint32_t demG = g.dem[giver], demT = g.dem[taker], capG = g.cap[giver], capT = g.cap[taker];
                const int out = (int)(demG >> 16), in = (int)(demT & 0xFFFFu), av = (int)(capG >> 16), rm = (int)(capT & 0xFFFFu);
                int units = planned;
                if (out > av) units = std::min(units, planned * av / out);
                if (in > rm) units = std::min(units, planned * rm / in);
                if (units < planned) *q4[k] = H(*q4[k] * ((float)units / (float)planned));
                if (units <= 0) continue;
                __atomic_fetch_sub(&g.delta[giver], units, __ATOMIC_RELAXED);
                __atomic_fetch_add(&g.delta[taker], units, __ATOMIC_RELAXED);
            }
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
            g.delta[gi] = 0; // the column flows planned next (#680), and their limiter's demand
            g.dem[gi] = 0u; g.cap[gi] = 0u;
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
            if ((fl & (SOURCE | SINK)) == (SOURCE | SINK)) g.sinkOut[slot] += (uint32_t)(nf - (int)TG(c)); // a sink's net intake (#656, #875: the columns counted none)
            if (fl & SOURCE) nf = (int)TG(c);
            else if (p.evaporate && nf > 0 && nf <= 2 && g.delta[gi] == 0 && onGround(slot, x, y, z)) nf = 0;
            // Changed against the cell before this substep (A, as the face rules compare), not after the vertical pass: water
            // that only fell inside a brick left it unflagged, so the unpack skipped it and the fall was lost (#680: the leak
            // test's sim volume fell 103 -> 55 m^3 in 60 s with nothing reaching the channel).
            if (F(outC[gi]) != (uint32_t)nf) g.flags[slot] |= 1u;
            outC[gi] = Pack((uint32_t)nf, (fl & 31u) | (kind << 5), c);
            if (nf > 0) markWet(slot, local);
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
        // What a released reach's sink takes in (#656), net (#875): a sink is a source held at the natural level, so once the
        // reach drains below it, it gives into the reach beside it and takes the same water back - counting only the intake read
        // that loop as outflow, and a cut upstream never reached the springs downstream. Signed (two's complement): what it gives
        // counts against.
        if ((fl & (SOURCE | SINK)) == (SOURCE | SINK)) g.sinkOut[slot] += (uint32_t)(nf - (int)TG(c));
        if (fl & SOURCE) nf = (int)TG(c);
        else if (p.evaporate && nf > 0 && nf <= 2 && moving == 0u && onGround(slot, x, y, z)) nf = 0;
        const uint32_t old = outC[gi];
        const uint32_t res = Pack((uint32_t)nf, (fl & 31u) | (kind << 5), c);
        if (F(old) != (uint32_t)nf) g.flags[slot] |= 1u;
        outC[gi] = res;
        if (nf > 0) markWet(slot, local);
    }
    // The dry-cell skip (#747): this substep's groups that hold water. A slot's cells run on one thread: no atomics.
    void markWet(int slot, int local) const { g.gmask[(size_t)slot * 6 + g.gcur * 2 + (local >> 11)] |= 1u << ((local >> 6) & 31); }
};

// One pass over every listed slot, the slots split across the job system's workers (this thread takes the first share).
void RunPass(Grid& g, const Params& p, const std::vector<int>& list, int pass) {
    const uint32_t* in = pass == 0 ? g.a.data() : g.b.data();
    uint32_t* out = pass == 0 ? g.b.data() : g.a.data();
    const Pass ps{ g, p, pass, in, out };
    const int n = (int)list.size();
    if (n == 0) return;
    const int parts = std::max(1, std::min(n, (int)std::thread::hardware_concurrency()));
    auto run = [&ps, &list, &g, n, parts](int part) {
        for (int i = part; i < n; i += parts) {
            const int slot = list[i];
            const uint64_t on = g.gon[slot];
            if (on == ~0ull) { for (int local = 0; local < 4096; ++local) ps.cell(slot, local); continue; }
            for (int grp = 0; grp < 64; ++grp) // the dry-cell skip (#747): only the groups near water
                if ((on >> grp) & 1ull)
                    for (int local = grp * 64; local < grp * 64 + 64; ++local) ps.cell(slot, local);
        }
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

void BeginTick(Grid& g, const std::vector<int>& list, int slots) {
    // Every slot's deltas (a ring slot's could still hold last tick's: colApply reads the cells over and under a head
    // without asking whether they are processed), then per processed slot its face records and its wet groups.
    std::fill(g.delta.begin(), g.delta.begin() + (size_t)slots * 4096, 0);
    const int n = (int)list.size();
    if (n == 0) { g.groupsFull = false; return; }
    const int parts = std::max(1, std::min(n, (int)std::thread::hardware_concurrency()));
    auto run = [&g, &list, n, parts](int part) {
        for (int i = part; i < n; i += parts) {
            const int s = list[i];
            std::fill(g.faceN.begin() + (size_t)s * 4096 * 2, g.faceN.begin() + (size_t)(s + 1) * 4096 * 2, Vec4{});
            uint32_t w[2] = { 0u, 0u };
            const uint32_t* c = &g.a[(size_t)s * 4096];
            for (int grp = 0; grp < 64; ++grp) {
                uint32_t any = 0u;
                for (int k = grp * 64; k < grp * 64 + 64; ++k) any |= c[k] & 255u;
                if (any) w[grp >> 5] |= 1u << (grp & 31);
            }
            uint32_t* m = &g.gmask[(size_t)s * 6]; // the host put the brick's history in 0 (last substep) and 2 (the one before)
            m[0] |= w[0]; m[1] |= w[1];            // and the groups wet now: the first substep lists them and their borders
        }
    };
    std::atomic<int> left{ parts - 1 };
    for (int part = 1; part < parts; ++part) {
        JobSystem::tagNext(JobTag::Fluid);
        WorldJobs().submit(-3.0f, [&run, &left, part] { run(part); left.fetch_sub(1); });
    }
    run(0);
    while (left.load() > 0) std::this_thread::yield();
    g.gcur = 0;
    g.groupsFull = false;
}

void Substep(Grid& g, const Params& p, const std::vector<int>& list) {
    // The dry-cell skip (#747, the CPU twin of fluid_groups.comp, #675): run only the 64-cell groups that held water in
    // either of the last two substeps or border one (in the brick, or across a face to a processed brick). Water moves at
    // most a cell a pass, so a skipped group can't change: the result is the same as running every cell. BeginTick seeds a
    // tick's first substep (slots are new each tick); FLUID_CPU_FULLFIRST=1 runs that one over every cell instead.
    // RAYLIBGAME_FLUID_GROUPS=0: off (A/B), as on the GPU.
    static const bool groupsOn = !(getenv("RAYLIBGAME_FLUID_GROUPS") && atoi(getenv("RAYLIBGAME_FLUID_GROUPS")) == 0);
    if (!groupsOn || g.groupsFull) {
        for (int s : list) { g.gon[s] = ~0ull; for (int w = 0; w < 6; ++w) g.gmask[(size_t)s * 6 + w] = 0u; }
        g.gcur = 0;
        g.groupsFull = false;
    } else {
        const int p1 = g.gcur, p2 = (g.gcur + 2) % 3; // the last two substeps' masks
        g.gcur = (g.gcur + 1) % 3;
        auto wet = [&](int s) -> uint64_t {
            const uint32_t* m = &g.gmask[(size_t)s * 6];
            return ((uint64_t)(m[p1 * 2 + 1] | m[p2 * 2 + 1]) << 32) | (uint64_t)(m[p1 * 2] | m[p2 * 2]);
        };
        auto wetNb = [&](int s, int face) -> uint64_t { const int o = g.nbr[(size_t)s * 6 + face]; return (o >= 0 && g.processed[o]) ? wet(o) : 0ull; };
        constexpr uint64_t q0 = 0x1111111111111111ull, q3 = q0 << 3;          // groups at y / 4 = 0 / 3 in every layer
        constexpr uint64_t z0 = 0xFull, z15 = 0xFull << 60;                    // the bottom / top layer's groups
        for (int s : list) {
            const uint64_t w = wet(s);
            uint64_t on = w | ((w << 1) & ~q0) | ((w >> 1) & ~q3) | (w << 4) | (w >> 4); // in the brick: y +-1 group, z +-1 layer
            on |= wetNb(s, 0) | wetNb(s, 1);                                    // across x: the same groups
            const uint64_t ym = wetNb(s, 2), yp = wetNb(s, 3);                  // across y: a layer's last / first group
            on |= ((ym & q3) >> 3) | ((yp & q0) << 3);
            const uint64_t zm = wetNb(s, 4), zp = wetNb(s, 5);                  // across z: the top / bottom layer
            on |= ((zm & z15) >> 60) | ((zp & z0) << 60);
            g.gon[s] = on;
            g.gmask[(size_t)s * 6 + g.gcur * 2] = 0u; g.gmask[(size_t)s * 6 + g.gcur * 2 + 1] = 0u; // the apply pass marks it
        }
    }
    static const bool quarter = EnvFlag("FLUID_COL_QUARTER");
    RunPass(g, p, list, 0);
    RunPass(g, p, list, 1);
    if (p.columns && !quarter) RunPass(g, p, list, 3); // the column pairs' scaled move (#680)
    RunPass(g, p, list, 2);
}

} // namespace FluidKernel

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): the CPU fluid kernel on a few bricks: exact conservation, determinism, the 255 clamp ----
#include "unit_helpers.h"

namespace {
using namespace FluidKernel;
constexpr uint32_t kSolidFlag = 1u; // flags bit 0 (fluid_step.comp's SOLID)
inline uint32_t Cell(uint32_t fill, uint32_t flags = 0u) { return (fill & 255u) | (flags << 8); }
inline int Idx(int x, int y, int z) { return x + 16 * y + 256 * z; }

// Params as the game sets them for water (fluid.cpp: kindVisc 5, kindMinSpread 10, KindFriction) on 0.25 m cells.
Params WaterParams() {
    Params p;
    for (int k = 0; k < 8; ++k) { p.minSpread[k] = 10; p.visc[k] = 5; p.frictionLin[k] = 0.05f; p.frictionQuad[k] = 0.005f; }
    p.dt = 1.0f / 120.0f;
    p.setCell(0.25f, 32);
    p.evaporate = false;
    return p;
}

// `n` bricks in a row along x (each one's +x face to the next), every one processed; a solid floor at z = 0 and a dam
// of solid cells at x = 8 of brick 0 (rows y 0..15) holding a block of water behind it, which `removeDam` takes away.
Grid DamBreak(int n, bool removeDam) {
    Grid g;
    g.resize(n);
    for (int s = 0; s < n; ++s) {
        g.processed[s] = 1;
        if (s > 0) g.nbr[(size_t)s * 6 + 0] = s - 1;
        if (s + 1 < n) g.nbr[(size_t)s * 6 + 1] = s + 1;
        for (int y = 0; y < 16; ++y) for (int x = 0; x < 16; ++x) g.a[(size_t)s * 4096 + Idx(x, y, 0)] = Cell(0, kSolidFlag);
    }
    for (int z = 1; z < 12; ++z) for (int y = 0; y < 16; ++y) for (int x = 0; x < 8; ++x) g.a[Idx(x, y, z)] = Cell(255);
    if (!removeDam) for (int z = 1; z < 16; ++z) for (int y = 0; y < 16; ++y) g.a[Idx(8, y, z)] = Cell(0, kSolidFlag);
    g.b = g.a;
    return g;
}
uint64_t TotalFill(const Grid& g) { uint64_t t = 0; for (uint32_t c : g.a) t += c & 255u; return t; }
std::vector<int> AllSlots(int n) { std::vector<int> l(n); for (int i = 0; i < n; ++i) l[i] = i; return l; }
void RunSubsteps(Grid& g, const Params& p, int substeps) {
    const std::vector<int> list = AllSlots((int)g.processed.size());
    BeginTick(g, list, (int)g.processed.size());
    for (int k = 0; k < substeps; ++k) Substep(g, p, list);
}
}

// Water moves but is neither made nor lost: the total fill is the same to the unit after every substep (#840).
TEST(FluidKernel, DamBreakConservesWaterExactly) {
    const Params p = WaterParams();
    Grid g = DamBreak(3, true);
    const uint64_t before = TotalFill(g);
    ASSERT_EQ(before, 11ull * 16 * 8 * 255);
    const std::vector<int> list = AllSlots(3);
    BeginTick(g, list, 3);
    int moved = 0;
    for (int k = 0; k < 600; ++k) { // 5 s
        Substep(g, p, list);
        ASSERT_EQ(TotalFill(g), before) << "substep " << k;
    }
    for (int x = 9; x < 16; ++x) moved += (int)(g.a[Idx(x, 8, 1)] & 255u);
    uint64_t inLater = 0;
    for (size_t i = 4096; i < g.a.size(); ++i) inLater += g.a[i] & 255u;
    EXPECT_GT(moved, 0) << "the water left the dam";
    EXPECT_GT(inLater, 0u) << "and crossed into the next brick";
}

// The same state stepped twice gives the same bits (the kernel splits slots across threads; the order must not matter).
TEST(FluidKernel, IsDeterministic) {
    const Params p = WaterParams();
    Grid g1 = DamBreak(3, true), g2 = DamBreak(3, true);
    RunSubsteps(g1, p, 240);
    RunSubsteps(g2, p, 240);
    EXPECT_TRUE(g1.a == g2.a) << "cells differ";
    bool facesSame = g1.faceP.size() == g2.faceP.size();
    for (size_t i = 0; facesSame && i < g1.faceP.size(); ++i)
        facesSame = memcmp(&g1.faceP[i], &g2.faceP[i], sizeof(Vec4)) == 0;
    EXPECT_TRUE(facesSame) << "face momenta differ";
}

// Water held by a dam stays put: still water at rest doesn't creep through solid cells or drift.
TEST(FluidKernel, DammedWaterStaysBehindTheDam) {
    const Params p = WaterParams();
    Grid g = DamBreak(2, false);
    const uint64_t before = TotalFill(g);
    RunSubsteps(g, p, 360);
    EXPECT_EQ(TotalFill(g), before);
    uint64_t past = 0;
    for (int s = 0; s < 2; ++s) for (int z = 0; z < 16; ++z) for (int y = 0; y < 16; ++y) for (int x = 0; x < 16; ++x)
        if (s == 1 || x >= 8) past += g.a[(size_t)s * 4096 + Idx(x, y, z)] & 255u;
    EXPECT_EQ(past, 0u);
}

// No cell ever holds more than its capacity (255 less its ground share), even under a tall column's pressure.
TEST(FluidKernel, CellsNeverOverfill) {
    const Params p = WaterParams();
    Grid g = DamBreak(1, true);
    for (int z = 1; z < 16; ++z) for (int y = 6; y < 10; ++y) for (int x = 10; x < 14; ++x) g.a[Idx(x, y, z)] = Cell(255);
    g.b = g.a;
    const uint64_t before = TotalFill(g);
    const std::vector<int> list = AllSlots(1);
    BeginTick(g, list, 1);
    for (int k = 0; k < 300; ++k) {
        Substep(g, p, list);
        for (uint32_t c : g.a) {
            const uint32_t fill = c & 255u, share = c >> 24;
            ASSERT_LE(fill, 255u - share);
            if ((c >> 8) & kSolidFlag) ASSERT_EQ(fill, 0u) << "water inside a solid cell";
        }
    }
    EXPECT_EQ(TotalFill(g), before);
}
#endif
