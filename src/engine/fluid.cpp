#include "fluid.h"
#include "cubemarch.h"
#include "gpucompute.h"
#include "jobsystem.h"
#include "object.h"
#include "terraindeltas.h"
#include "voxel.h"
#include "waterfield.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <deque>
#include <map>
#include <memory>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <thread>

namespace {

constexpr int B = 16, B3 = B * B * B;
using WaterField::kSolid; using WaterField::kSource; using WaterField::kStatic; // kStatic: inside infinite generated water
using WaterField::kDug;                // open only because of an edit: generated infinite water may flow into it
using WaterField::kKindShift;          // flag bits 5-7: the cell's fluid kind (0 water .. 7; Fluid::Kind, fluids.txt)
using WaterField::kKindMask;
using WaterField::CellKind;
// The kind table, read by the step job (written only at startup): sideways divisor and film limit per kind.
std::vector<Fluid::Kind> kinds(1);
int kindVisc[8] = { 5, 5, 5, 5, 5, 5, 5, 5 }, kindMinSpread[8] = { 10, 10, 10, 10, 10, 10, 10, 10 };
constexpr int kSiteBricks = 16;        // a site reaches +-16 bricks (256 m) from its origin
// Tests only: RAYLIBGAME_FLUID_SITE_BRICKS=N shrinks a site so a pond straddles its edge (#384).
inline int SiteBricks() { static const int n = getenv("RAYLIBGAME_FLUID_SITE_BRICKS") ? std::max(1, atoi(getenv("RAYLIBGAME_FLUID_SITE_BRICKS"))) : kSiteBricks; return n; }
constexpr float kSiteJoin = 180.0f;    // a disturbance within this of a site's origin joins it
constexpr int kEvaporate = 3;          // resting specks this thin (~1 cm) dry up
constexpr int kSubsteps = 3;           // CPU simulation steps per 50 ms (water moves a cell a step)
constexpr double kTickSeconds = 0.025; // a tick (a job + meshes) every 25 ms: water moved in 1-3 cell jumps every 50-100 ms, hard to follow (#366)
constexpr int kSleepSteps = 120;       // 2 s without change: a brick sleeps
constexpr int kSettleSteps = 180;      // 3 s window: a brick whose cells all moved less than kSettleDiff in it sleeps (#366)
constexpr int kSettleDiff = 4;         // fill units (1.6 cm)
constexpr int kLevelEvery = 15;        // levelling pass every 0.25 s
constexpr int kInitPerStep = 24;       // bricks classified per step, in parallel on the worker pool (each a few ms)
constexpr int kMaxBricks = 4096;       // per site

inline int Idx(int x, int y, int z) { return x + y * B + z * B * B; }
inline int64_t Key(int i, int j, int k) { return ((int64_t)(i + 1024) << 22) | ((int64_t)(j + 1024) << 11) | (int64_t)(k + 1024); }
inline int FloorDiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }

struct Brick {
    int bi = 0, bj = 0, bk = 0;
    std::array<uint8_t, B3> fill{}, flags{}, target{};
    Brick* nbr[6] = {};             // -x +x -y +y -z +z
    bool ready = false;             // cells classified (solid / source)
    bool keepFill = false;          // loaded from a save: the classification keeps the stored water
    bool basin = false;             // holds finite generated water (a pond): saved even once drained (#360)
    bool active = true;
    int idle = 0;
    bool changed = true;            // since the last publish
    std::array<uint8_t, B3> check{}; long long checkStep = -1; // the cells a settle window ago (#366: resting water sleeps)
    bool settled = false;           // put to sleep resting (#366): small changes (the GPU step at its edge) do not wake it
    std::array<uint8_t, B3> drawn{}, drawnFlags{}; // the cells as the chunks last redrew them (visible changes only)
    bool drawnSet = false;
    int slot = -1;                  // GPU mode: its place in the GPU brick pool (-1: none yet)
    bool upload = false;            // GPU mode: the whole brick (cells, flags, targets) goes up next tick
};

struct Site {
    uint64_t body = 0;
    int index = 0;
    Vector3 origin{}, e{}, n{}, u{};
    std::unordered_map<int64_t, std::unique_ptr<Brick>> bricks;
    std::deque<Brick*> pending;     // waiting for classification
    // Sources with a rate and a head (springs, seeps): site cell, kind, fill units added per step at most,
    // and the level (site cells, fractional) the fluid surface over it is held to.
    struct Spring { int ci, cj, ck, kind; float perStep, head; float acc = 0.0f; }; // acc: fill units owed (a slow source adds < 1 a step)
    std::vector<Spring> springs;

    Vector3 cellCentre(int ci, int cj, int ck) const {
        return Vector3Add(origin, Vector3Add(Vector3Add(Vector3Scale(e, ci + 0.5f), Vector3Scale(n, cj + 0.5f)), Vector3Scale(u, ck + 0.5f)));
    }
    void cellOf(Vector3 local, int& ci, int& cj, int& ck) const {
        const Vector3 d = Vector3Subtract(local, origin);
        ci = (int)floorf(Vector3DotProduct(d, e)); cj = (int)floorf(Vector3DotProduct(d, n)); ck = (int)floorf(Vector3DotProduct(d, u));
    }
    bool inRange(int bi, int bj, int bk) const { const int r = SiteBricks(); return abs(bi) <= r && abs(bj) <= r && abs(bk) <= r; }
    Brick* find(int bi, int bj, int bk) const { auto it = bricks.find(Key(bi, bj, bk)); return it == bricks.end() ? nullptr : it->second.get(); }
    // Allocates (unclassified, queued) if missing and in range; nullptr if out of range or full.
    Brick* ensure(int bi, int bj, int bk) {
        if (Brick* b = find(bi, bj, bk)) return b;
        if (!inRange(bi, bj, bk) || (int)bricks.size() >= kMaxBricks) return nullptr;
        auto nb = std::make_unique<Brick>();
        nb->bi = bi; nb->bj = bj; nb->bk = bk;
        Brick* b = nb.get();
        bricks[Key(bi, bj, bk)] = std::move(nb);
        const int d[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
        for (int f = 0; f < 6; ++f)
            if (Brick* o = find(bi + d[f][0], bj + d[f][1], bk + d[f][2])) { b->nbr[f] = o; o->nbr[f ^ 1] = b; }
        pending.push_back(b);
        return b;
    }
    // A cell by site coordinates: its brick (ready) and index; false if missing or unclassified.
    bool cell(int ci, int cj, int ck, Brick*& b, int& idx) const {
        b = find(FloorDiv(ci, B), FloorDiv(cj, B), FloorDiv(ck, B));
        if (!b || !b->ready) return false;
        idx = Idx(ci - b->bi * B, cj - b->bj * B, ck - b->bk * B);
        return true;
    }
};

struct BodyRef {
    std::shared_ptr<const DensityField> field;
    std::shared_ptr<const TerrainDeltas> edits;
    const PlanetField* pf = nullptr;
};

struct Command { int kind; uint64_t body; Vector3 c; float r, m3; int fluid = 0; float head = 0.0f; }; // 0 edit, 1 take, 2 pour, 3 add a source (r = rate m^3/s); fluid = the kind

// Simulation state: owned by the step job while it runs (simMx).
std::mutex simMx;
std::vector<std::unique_ptr<Site>> sites;
std::unordered_map<uint64_t, BodyRef> refs;     // snapshot handed to the job
long long stepCount = 0;
Fluid::Stats stats;
long long settledTotal = 0;     // bricks put to sleep resting (#366), since the start
long long dbgDown = 0, dbgSide = 0, dbgSrc = 0, dbgLevel = 0; // tests: changes by kind (since the last publish)

// Main thread <-> job.
std::mutex cmdMx;
std::vector<Command> commands;
std::atomic<bool> running{ false };
std::atomic<int> speed{ 1 };
double accum = 0.0;

// Published results: the cells go out as WaterField snapshots (waterfield.h: the chunks draw them, queries read them).
std::mutex pubMx;
Fluid::Stats pubStats;

// ---------------------------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------------------------

// Cells of a brick (#360: the one water field): solid from the terrain density (edits included); open cells under
// generated water hold it. Infinite basins (the sea, rivers, large lakes) give sources held at their level, which feed
// only cells an edit opened; finite ones (ponds) give ordinary water, once: a brick classified again (an edit) or
// loaded keeps what it holds, so a drained pond stays drained. There is no take-over: the chunks draw the published
// cells and the generated water as one surface. `faces` (bits -x +x -y +y -z +z): faces holding finite generated
// water; RunJob pulls those neighbours in, so a pond joins the simulation whole and drains as one body.
// Runs in parallel over bricks (pool jobs): it writes only its own brick.
void Classify(const Site& s, Brick& b, const BodyRef& ref, uint8_t& faces) {
    const PlanetField* pf = ref.pf;
    const float sea = pf ? pf->seaLevelRadius() : 0.0f;
    b.basin = false;
    for (int y = 0; y < B; ++y)
        for (int x = 0; x < B; ++x) {
            // The column's ground and generated water (from its middle: a column tilts < 0.3 deg in a site).
            const int ci = b.bi * B + x, cj = b.bj * B + y;
            const Vector3 mid = s.cellCentre(ci, cj, b.bk * B + B / 2);
            const float ml = Vector3Length(mid);
            if (ml < 1.0f) continue;
            const Vector3 dir = Vector3Scale(mid, 1.0f / ml);
            float bed = 0.0f, level = -1e9f;
            int kind = 0;                               // 1 river, 2 lake, 4 sea: infinite; 5 a pond: finite
            if (pf) {
                bed = pf->radius + pf->elevation(dir);
                if (bed < sea) { level = sea; kind = 4; } // the sea
                const PlanetField::WaterSource ws = pf->waterSourceAt(dir);
                if (ws.kind == 3) {
                    if (bed < sea + ws.level && sea + ws.level > level) { level = sea + ws.level; kind = 5; } // a pond (never an infinite source, even its rim)
                } else if (kind != 4 && ws.kind && sea + ws.level > level) { // at sea the sea's level wins (#190: a river mouth's higher level fed the sea columns: overfilled channels that never slept)
                    level = sea + ws.level; kind = ws.kind; }
            }
            const float top = pf ? bed + 14.0f : 1e9f;  // tors, boulders and detail stay below this
            for (int z = 0; z < B; ++z) {
                const int i = Idx(x, y, z);
                const Vector3 c = s.cellCentre(ci, cj, b.bk * B + z);
                const float r = Vector3Length(c);
                const bool solid = r < top && ref.field->density(c, ref.edits.get()) >= ISO_LEVEL;
                uint8_t fl = solid ? kSolid : 0;
                // Dug: open only because of an edit. Infinite water feeds only these (the real ground has
                // crags the height map doesn't: natural dips by a river are not to be flooded).
                if (!solid && r < bed + 4.0f && r < top && ref.field->density(c, nullptr) >= ISO_LEVEL) fl |= kDug;
                uint8_t tg = 0;
                if (!solid && kind && r - 0.5f < level && r > bed - 0.5f) {
                    tg = (uint8_t)std::clamp((int)lroundf((level - (r - 0.5f)) * 255.0f), 0, 255);
                    if (tg > 0 && kind != 5) fl |= kSource | kStatic;
                    if (tg > 0 && kind == 5) {
                        b.basin = true;
                        faces |= (x == 0 ? 1 : 0) | (x == B - 1 ? 2 : 0) | (y == 0 ? 4 : 0) | (y == B - 1 ? 8 : 0) | (z == 0 ? 16 : 0) | (z == B - 1 && tg == 255 ? 32 : 0);
                    }
                }
                if (b.keepFill && !solid && !(fl & kSource) && b.fill[i] > 0) fl |= (b.flags[i] & kKindMask); // kept fluid stays its kind
                b.flags[i] = fl;
                b.target[i] = (fl & kSource) ? tg : 0;
                if (b.keepFill) { if (solid) b.fill[i] = 0; else if (fl & kSource) b.fill[i] = tg; }
                else b.fill[i] = solid ? 0 : tg;
            }
        }
    // Kept (loaded or reclassified) bricks pull their pond neighbours in too (#382 review): ensure() is a no-op for those
    // that exist; one never classified (saved mid pull-in, an old save's claimed pond) comes in fresh with its generated
    // water, so a drained brick never borders the generated pond (a wall of water at the brick face).
    b.ready = true;
    b.upload = true; // GPU mode: the classified brick goes up whole
    b.keepFill = false;
    b.active = true; b.idle = 0; b.changed = true;
}

// A pond's water at face f whose neighbour brick cannot join (outside the site, or the site is full, #384): the
// face's water cells hold the pond's generated level as sources, so the pond keeps that level like a lake and feeds
// a dug channel, instead of draining inside the site against the generated pond's wall at the brick face. Runs right
// after Classify, before the brick goes up to the GPU (b.upload), on its generated fill.
long long pinnedFaces = 0;
void PinPondFace(Brick& b, int f) {
    ++pinnedFaces;
    for (int u = 0; u < B; ++u)
        for (int v = 0; v < B; ++v) {
            const int m = (f & 1) ? B - 1 : 0;
            const int x = f < 2 ? m : u, y = f < 2 ? u : (f < 4 ? m : v), z = f < 4 ? v : m;
            const int i = Idx(x, y, z);
            if ((b.flags[i] & kSolid) || b.fill[i] == 0) continue;
            b.flags[i] |= kSource;
            b.target[i] = b.fill[i];
        }
}

// ---------------------------------------------------------------------------------------------
// Step
// ---------------------------------------------------------------------------------------------

struct Neighbour { Brick* b; int i; };
// A cell's neighbour along an axis (dx, dy, dz one of +-1): crosses into the next brick.
inline bool Next(Brick* b, int x, int y, int z, int dx, int dy, int dz, Neighbour& out) {
    x += dx; y += dy; z += dz;
    int f = -1;
    if (x < 0) { f = 0; x += B; } else if (x >= B) { f = 1; x -= B; }
    else if (y < 0) { f = 2; y += B; } else if (y >= B) { f = 3; y -= B; }
    else if (z < 0) { f = 4; z += B; } else if (z >= B) { f = 5; z -= B; }
    if (f >= 0) { b = b->nbr[f]; if (!b || !b->ready) return false; }
    out = { b, Idx(x, y, z) };
    return true;
}

long long wakeGpu = 0, wakeTouch = 0; // settled bricks woken (stats)
void Touch(Brick* b) { if (b->settled) ++wakeTouch; b->active = true; b->idle = 0; b->changed = true; b->settled = false; }

// GPU mode (fluid on compute shaders): CPU-side changes to the mirror (buckets, pumps, springs, levelling) go
// to the GPU as deltas next tick; bricks waiting for a full upload take theirs from the mirror instead.
bool gpuMode = false;
struct GpuDeltaRec { Brick* b; int i; int d; uint8_t kind; };
std::vector<GpuDeltaRec> gpuDeltas;
inline void GpuDelta(Brick* b, int i, int d) {
    if (!gpuMode || d == 0 || b->slot < 0 || b->upload) return;
    gpuDeltas.push_back({ b, i, d, (uint8_t)CellKind(b->flags[i]) });
}

// Kinds don't mix: fluid moves only into an empty cell or one of its own kind; an empty cell
// takes the kind of what first flows in.
inline bool SameKind(const Brick* a, int ia, const Brick* b, int ib) { return b->fill[ib] == 0 || ((a->flags[ia] ^ b->flags[ib]) & kKindMask) == 0; }
inline void TakeKind(Brick* to, int it, uint8_t kindBits) { if (to->fill[it] == 0) to->flags[it] = (uint8_t)((to->flags[it] & ~kKindMask) | kindBits); }

// Water moving out of a cell into a missing brick: queue that brick (it takes water from the next
// step on, once classified).
void RequestBeyond(Site& s, Brick* b, int x, int y, int z, int dx, int dy, int dz) {
    const int ci = b->bi * B + x + dx, cj = b->bj * B + y + dy, ck = b->bk * B + z + dz;
    s.ensure(FloorDiv(ci, B), FloorDiv(cj, B), FloorDiv(ck, B));
}

// A source with a rate and a head: it adds to the surface cell of its column while that surface is
// below the head. A basin over it fills to the head and the source rests; a lower rim lets the overflow
// run off (the step carries it away), so the source keeps adding and a stream forms.
void StepSpring(Site& s, Site::Spring& sp) {
    sp.acc = std::min(sp.acc + sp.perStep, std::max(16.0f, sp.perStep * 30.0f)); // owed at most ~half a second's worth
    if (sp.acc < 8.0f) return; // in drops of 3 cm+: thinner specks would soak away (kEvaporate) before they spread
    for (int k = 0; k < 64; ++k) {
        const int ck = sp.ck + k;
        if ((float)ck >= sp.head) return;                       // the surface is at the head
        Brick* b; int i;
        if (!s.cell(sp.ci, sp.cj, ck, b, i)) { s.ensure(FloorDiv(sp.ci, B), FloorDiv(sp.cj, B), FloorDiv(ck, B)); return; }
        if (b->flags[i] & kSolid) { if (k < 3) continue; return; }  // a buried source wells up through the ground just above it
        if (b->fill[i] > 0 && CellKind(b->flags[i]) != sp.kind) return; // another fluid lies on it
        if (b->fill[i] >= 255) continue;                          // full: the surface is higher up
        const float surface = ck + b->fill[i] / 255.0f;
        const int add = (int)std::min({ sp.acc, (sp.head - surface) * 255.0f, 255.0f - b->fill[i] });
        if (add <= 0) return;
        TakeKind(b, i, (uint8_t)(sp.kind << kKindShift));
        b->fill[i] = (uint8_t)(b->fill[i] + add);
        GpuDelta(b, i, add);
        sp.acc -= (float)add;
        Touch(b); ++dbgSrc;
        return;
    }
}

void StepSite(Site& s, int parity) {
    for (Site::Spring& sp : s.springs) StepSpring(s, sp); // springs and seeps
    std::vector<Brick*> act;
    for (auto& kv : s.bricks) if (kv.second->ready && kv.second->active) act.push_back(kv.second.get());
    std::sort(act.begin(), act.end(), [](const Brick* a, const Brick* b) { return a->bk != b->bk ? a->bk < b->bk : a->bj != b->bj ? a->bj < b->bj : a->bi < b->bi; });
    static const int side[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
    for (Brick* b : act) {
        bool moved = false;
        for (int z = 0; z < B; ++z)
            for (int y = 0; y < B; ++y)
                for (int x = 0; x < B; ++x) {
                    const int i = Idx(x, y, z);
                    int f = b->fill[i];
                    if (f == 0 || (b->flags[i] & kSolid)) continue;
                    const bool src = b->flags[i] & kSource;
                    // Down.
                    Neighbour d;
                    bool belowOpen = false, belowFull = true;
                    if (Next(b, x, y, z, 0, 0, -1, d)) {
                        if (!(d.b->flags[d.i] & kSolid)) {
                            belowOpen = true;
                            const int room = SameKind(b, i, d.b, d.i) ? 255 - d.b->fill[d.i] : 0; // oil rests on water (and water on oil)
                            if (room > 0 && !(d.b->flags[d.i] & kSource) && (!src || (d.b->flags[d.i] & kDug))) { // static water takes nothing in
                                const int m = std::min(f, room);
                                TakeKind(d.b, d.i, b->flags[i] & kKindMask);
                                d.b->fill[d.i] = (uint8_t)(d.b->fill[d.i] + m);
                                f -= m;
                                if (d.b != b) Touch(d.b);
                                moved = true; ++dbgDown;
                            }
                            belowFull = d.b->fill[d.i] == 255 || !SameKind(b, i, d.b, d.i);
                        }
                    } else if (z == 0 && !b->nbr[4]) {
                        if (!src) RequestBeyond(s, b, x, y, z, 0, 0, -1); // static water never asks for more bricks
                        belowOpen = true; belowFull = false; // wait for it
                    }
                    // Sideways, once the cell rests on something.
                    const int kd = CellKind(b->flags[i]), visc = kindVisc[kd], minSpread = kindMinSpread[kd]; // per fluid kind
                    if (f >= minSpread && (!belowOpen || belowFull)) {
                        const int f0 = f;
                        for (int k = 0; k < 4; ++k) {
                            const int kk = (k + parity) & 3;
                            Neighbour o;
                            if (!Next(b, x, y, z, side[kk][0], side[kk][1], 0, o)) {
                                // Out of the allocated region: ask for the brick.
                                const int nx = x + side[kk][0], ny = y + side[kk][1];
                                if (!src && (nx < 0 || nx >= B || ny < 0 || ny >= B)) RequestBeyond(s, b, x, y, z, side[kk][0], side[kk][1], 0);
                                continue;
                            }
                            if (o.b->flags[o.i] & kSolid) continue;
                            if ((o.b->flags[o.i] & kSource) || (src && !(o.b->flags[o.i] & kDug))) continue; // sources feed dug cells only, and take nothing in
                            if (!SameKind(b, i, o.b, o.i)) continue;
                            const int diff = f0 - o.b->fill[o.i];
                            if (diff < visc) continue;
                            const int m = std::min(diff / visc, f); // thicker fluids spread slower
                            if (o.b->fill[o.i] == 0 && m < minSpread) continue; // no trickles into dry cells: a film's edge holds
                            if (m <= 0) continue;
                            TakeKind(o.b, o.i, b->flags[i] & kKindMask);
                            o.b->fill[o.i] = (uint8_t)(o.b->fill[o.i] + m);
                            f -= m;
                            if (o.b != b) Touch(o.b);
                            moved = true; ++dbgSide;
                        }
                    }
                    if (f > 0 && f <= kEvaporate && !src && (!belowOpen || belowFull)) { f = 0; moved = true; } // a resting speck soaks away
                    b->fill[i] = (uint8_t)f;
                }
        // Sources hold their level.
        for (int i = 0; i < B3; ++i)
            if ((b->flags[i] & kSource) && b->fill[i] != b->target[i]) { b->fill[i] = b->target[i]; moved = true; ++dbgSrc; }
        if (moved) Touch(b);
        else if (++b->idle > kSleepSteps) b->active = false;
    }
}

// Levelling: bodies of water joined through full cells communicate. Each body's open surfaces (the
// partial cells and the empty cells over full ones) ease toward one level: its sources' level if it
// has any, else the level its volume fills to.
void LevelSite(Site& s) {
    std::unordered_set<int64_t> seen;
    struct C { Brick* b; int i; int ck; bool leaf; };
    for (auto& kv : s.bricks) {
        Brick* b0 = kv.second.get();
        if (!b0->ready || !b0->active) continue;
        for (int i0 = 0; i0 < B3; ++i0) {
            if (b0->fill[i0] != 255 || (b0->flags[i0] & kSolid)) continue;
            const int x0 = i0 % B, y0 = (i0 / B) % B, z0 = i0 / (B * B);
            const int64_t k0 = Key(b0->bi * B + x0, b0->bj * B + y0, b0->bk * B + z0);
            if (seen.count(k0)) continue;
            const uint8_t kind0 = b0->flags[i0] & kKindMask; // a body is all one kind
            // BFS through full cells; partial and empty-above neighbours are leaves.
            std::vector<C> comp;
            std::deque<std::array<int, 3>> q;
            q.push_back({ b0->bi * B + x0, b0->bj * B + y0, b0->bk * B + z0 });
            seen.insert(k0);
            bool hasSource = false;
            float srcLevel = -1e9f, srcSurface = 1e9f; // highest source top; lowest partial (surface) source
            std::unordered_set<int64_t> leafSeen;
            while (!q.empty() && comp.size() < 150000) {
                const auto c = q.front(); q.pop_front();
                Brick* b; int i;
                if (!s.cell(c[0], c[1], c[2], b, i)) continue;
                comp.push_back({ b, i, c[2], false });
                if (b->flags[i] & kSource) { hasSource = true; srcLevel = std::max(srcLevel, c[2] + b->target[i] / 255.0f); if (b->target[i] < 255) srcSurface = std::min(srcSurface, c[2] + b->target[i] / 255.0f); }
                const int d[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
                for (int k = 0; k < 6; ++k) {
                    const int nc[3] = { c[0] + d[k][0], c[1] + d[k][1], c[2] + d[k][2] };
                    Brick* nb; int ni;
                    if (!s.cell(nc[0], nc[1], nc[2], nb, ni) || (nb->flags[ni] & kSolid)) continue;
                    if (nb->fill[ni] > 0 && (nb->flags[ni] & kKindMask) != kind0) continue; // the other fluid
                    const int64_t nk = Key(nc[0], nc[1], nc[2]);
                    if (nb->fill[ni] == 255) {
                        if (!seen.count(nk)) { seen.insert(nk); q.push_back({ nc[0], nc[1], nc[2] }); }
                    } else if ((nb->fill[ni] > 0 || k == 4) && !leafSeen.count(nk)) {
                        leafSeen.insert(nk);
                        comp.push_back({ nb, ni, nc[2], true });
                        if (nb->flags[ni] & kSource) { hasSource = true; srcLevel = std::max(srcLevel, nc[2] + nb->target[ni] / 255.0f); if (nb->target[ni] < 255) srcSurface = std::min(srcSurface, nc[2] + nb->target[ni] / 255.0f); }
                    }
                }
            }
            if (comp.size() >= 150000) continue;
            float lo = 1e9f, hi = -1e9f, vol = 0.0f;
            int kmin = INT32_MAX, kmax = INT32_MIN;
            for (const C& c : comp) {
                vol += c.b->fill[c.i] / 255.0f;
                kmin = std::min(kmin, c.ck); kmax = std::max(kmax, c.ck);
                if (c.leaf) { const float h = c.ck + c.b->fill[c.i] / 255.0f; lo = std::min(lo, h); hi = std::max(hi, h); }
            }
            if (hi < lo) continue; // no open surface (a sealed cave full to the roof)
            float level;
            if (hasSource) {
                level = srcSurface < 1e8f ? srcSurface : srcLevel; // a sloping river: its lowest surface, never overfilling
                if (lo > level - 0.25f && hi < level + 0.25f) continue;
            } else {
                if (hi - lo < 0.3f) continue;
                float a = (float)kmin, z = (float)kmax + 1.0f;
                for (int it = 0; it < 24; ++it) {
                    const float m = 0.5f * (a + z);
                    float cap = 0.0f;
                    for (const C& c : comp) cap += std::clamp(m - c.ck, 0.0f, 1.0f);
                    if (cap < vol) a = m; else z = m;
                }
                level = 0.5f * (a + z);
            }
            // Ease toward the level (30% a pass): a flooded cave rises, a U-bend's arms meet. A body within
            // a few 255ths everywhere is left alone (so it sleeps); otherwise the rounding is carried from
            // cell to cell so a finite body keeps its volume.
            float worst = 0.0f;
            for (const C& c : comp) {
                if (c.b->flags[c.i] & kSource) continue;
                const float want = std::clamp(level - c.ck, 0.0f, 1.0f) * 255.0f, now = c.b->fill[c.i];
                if (hasSource && want <= now) continue; // fed by static water: only raise (a sloping river feeds at several levels)
                worst = std::max(worst, fabsf(want - now));
            }
            if (worst < 6.0f) continue;
            float carry = 0.0f;
            for (const C& c : comp) {
                if (c.b->flags[c.i] & kSource) continue;
                const float want = std::clamp(level - c.ck, 0.0f, 1.0f) * 255.0f, now = c.b->fill[c.i];
                if (hasSource && want <= now) continue;
                const float d = (want - now) * 0.3f + (hasSource ? 0.0f : carry);
                const int di = (int)lroundf(d);
                const int nf = std::clamp((int)now + di, 0, 255);
                carry = d - (float)(nf - (int)now);
                if (nf != c.b->fill[c.i]) { if (nf > 0) TakeKind(c.b, c.i, kind0); GpuDelta(c.b, c.i, nf - (int)c.b->fill[c.i]); c.b->fill[c.i] = (uint8_t)nf; if (c.b->settled && abs(nf - (int)c.b->check[c.i]) <= 2 * kSettleDiff) c.b->changed = true; else Touch(c.b); ++dbgLevel; } // a resting brick: small levelling nudges do not wake it (#366)
            }
        }
    }
}
// ---------------------------------------------------------------------------------------------
// The job
// ---------------------------------------------------------------------------------------------

Site* SiteFor(uint64_t body, Vector3 c, bool create) {
    Site* best = nullptr;
    float bd = kSiteJoin;
    for (auto& s : sites) {
        if (s->body != body) continue;
        const float d = Vector3Distance(s->origin, c);
        if (d < bd) { bd = d; best = s.get(); }
    }
    if (best || !create) return best;
    auto s = std::make_unique<Site>();
    s->body = body;
    s->index = (int)sites.size();
    const Vector3 up = Vector3Normalize(c);
    s->u = up;
    s->e = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    s->n = Vector3CrossProduct(up, s->e);
    s->origin = { floorf(c.x), floorf(c.y), floorf(c.z) };
    sites.push_back(std::move(s));
    return sites.back().get();
}

// Static water (any kind) within reach of a point: a site is worth starting.
bool StaticWaterNear(const BodyRef& ref, Vector3 c, float r) {
    if (!ref.pf) return false;
    const float sea = ref.pf->seaLevelRadius();
    const Vector3 up = Vector3Normalize(c);
    const Vector3 t1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
    const Vector3 t2 = Vector3CrossProduct(up, t1);
    const Vector3 offs[9] = { { 0, 0, 0 }, t1, Vector3Negate(t1), t2, Vector3Negate(t2),
                              Vector3Normalize(Vector3Add(t1, t2)), Vector3Normalize(Vector3Subtract(t1, t2)), Vector3Normalize(Vector3Add(Vector3Negate(t1), t2)), Vector3Normalize(Vector3Negate(Vector3Add(t1, t2))) };
    for (const Vector3& o : offs) {
        const Vector3 p = Vector3Add(c, Vector3Scale(o, r + 1.5f));
        const float pl = Vector3Length(p);
        const Vector3 dir = Vector3Scale(p, 1.0f / pl);
        float level = sea;
        const PlanetField::WaterSource ws = ref.pf->waterSourceAt(dir);
        if (ws.kind) level = std::max(level, sea + ws.level);
        if (pl - r - 1.0f < level && ref.pf->radius + ref.pf->elevation(dir) < level + 0.5f) return true;
    }
    return false;
}

bool SiteHasWaterNear(Site& s, Vector3 c, float r) {
    int ci, cj, ck;
    s.cellOf(c, ci, cj, ck);
    const int R = (int)ceilf(r) + 2;
    for (int dk = -R; dk <= R; ++dk)
        for (int dj = -R; dj <= R; ++dj)
            for (int di = -R; di <= R; ++di) {
                Brick* b; int i;
                if (s.cell(ci + di, cj + dj, ck + dk, b, i) && b->fill[i] > 0) return true;
            }
    return false;
}

void EnsureAround(Site& s, Vector3 c, float r) {
    int ci, cj, ck;
    s.cellOf(c, ci, cj, ck);
    const int R = (int)ceilf(r) + 2;
    for (int bk = FloorDiv(ck - R, B); bk <= FloorDiv(ck + R, B); ++bk)
        for (int bj = FloorDiv(cj - R, B); bj <= FloorDiv(cj + R, B); ++bj)
            for (int bi = FloorDiv(ci - R, B); bi <= FloorDiv(ci + R, B); ++bi)
                if (Brick* b = s.ensure(bi, bj, bk)) {
                    if (b->ready) { // reclassify (solidity changed); keep the water
                        b->ready = false; b->keepFill = true;
                        if (std::find(s.pending.begin(), s.pending.end(), b) == s.pending.end()) s.pending.push_front(b);
                    }
                }
}

void ApplyCommand(const Command& c) {
    auto rit = refs.find(c.body);
    if (rit == refs.end()) return;
    const BodyRef& ref = rit->second;
    if (c.kind == 0) {
        // A terrain edit: refresh the bricks it reaches; start a site if it opens space by water.
        Site* s = SiteFor(c.body, c.c, false);
        const bool near = (s && SiteHasWaterNear(*s, c.c, c.r)) || StaticWaterNear(ref, c.c, c.r);
        if (!s && !near) return;
        if (!s) s = SiteFor(c.body, c.c, true);
        if (near) EnsureAround(*s, c.c, c.r);
        else {
            // Only bricks that exist: their solidity changed.
            int ci, cj, ck; s->cellOf(c.c, ci, cj, ck);
            const int R = (int)ceilf(c.r) + 2;
            for (int bk = FloorDiv(ck - R, B); bk <= FloorDiv(ck + R, B); ++bk)
                for (int bj = FloorDiv(cj - R, B); bj <= FloorDiv(cj + R, B); ++bj)
                    for (int bi = FloorDiv(ci - R, B); bi <= FloorDiv(ci + R, B); ++bi)
                        if (Brick* b = s->find(bi, bj, bk); b && b->ready) { b->ready = false; b->keepFill = true; s->pending.push_front(b); }
        }
        return;
    }
    if (c.kind == 3) {
        // A source (spring, seep): rate c.r m^3/s, head c.head m above the point. One per cell.
        Site* s = SiteFor(c.body, c.c, true);
        int ci, cj, ck;
        s->cellOf(c.c, ci, cj, ck);
        const float fz = Vector3DotProduct(Vector3Subtract(c.c, s->origin), s->u);
        Site::Spring sp{ ci, cj, ck, std::clamp(c.fluid, 0, 7), c.r * 255.0f / (20.0f * kSubsteps), fz + c.head };
        bool replaced = false;
        for (Site::Spring& o : s->springs) if (!replaced && abs(o.ci - ci) <= 2 && abs(o.cj - cj) <= 2 && abs(o.ck - ck) <= 2) { o = sp; replaced = true; } // the same source again (a reload)
        if (!replaced) s->springs.push_back(sp);
        EnsureAround(*s, c.c, 2.0f);
        return;
    }
    Site* s = SiteFor(c.body, c.c, true);
    int ci, cj, ck;
    s->cellOf(c.c, ci, cj, ck);
    if (c.kind == 2) {
        // Pour: into the first open cell at or above the point, spilling upward.
        float left = c.m3 * 255.0f;
        Brick* b; int i;
        if (!s->cell(ci, cj, ck, b, i)) { EnsureAround(*s, c.c, 1.0f); std::lock_guard<std::mutex> lock(cmdMx); commands.push_back(c); return; } // after classification
        for (int k = 0; k < 6 && left > 0.0f; ++k) {
            if (!s->cell(ci, cj, ck + k, b, i)) break;
            if (b->flags[i] & kSolid) continue;
            if (b->fill[i] > 0 && CellKind(b->flags[i]) != c.fluid) continue; // another fluid: spill past it
            TakeKind(b, i, (uint8_t)(c.fluid << kKindShift));
            const int add = (int)std::min(left, 255.0f - b->fill[i]);
            b->fill[i] = (uint8_t)(b->fill[i] + add);
            GpuDelta(b, i, add);
            left -= add;
            Touch(b);
        }
        return;
    }
    // Take: from the wettest cells within the radius (sources are free: they refill).
    float want = c.m3 * 255.0f;
    const int R = (int)ceilf(c.r);
    EnsureAround(*s, c.c, c.r);
    for (int dk = R; dk >= -R && want > 0.0f; --dk)
        for (int dj = -R; dj <= R && want > 0.0f; ++dj)
            for (int di = -R; di <= R && want > 0.0f; ++di) {
                Brick* b; int i;
                if (!s->cell(ci + di, cj + dj, ck + dk, b, i) || b->fill[i] == 0) continue;
                if (CellKind(b->flags[i]) != c.fluid) continue; // only the asked-for fluid
                const int t = (int)std::min(want, (float)b->fill[i]);
                if (!(b->flags[i] & kSource)) { b->fill[i] = (uint8_t)(b->fill[i] - t); GpuDelta(b, i, -t); }
                want -= t;
                Touch(b);
            }
}

void RunJob() {
    const auto t0 = std::chrono::steady_clock::now();
    std::lock_guard<std::mutex> lock(simMx);
    std::vector<Command> cmds;
    { std::lock_guard<std::mutex> c(cmdMx); cmds.swap(commands); }
    for (const Command& c : cmds) ApplyCommand(c);
    // Classify queued bricks, in parallel on the worker pool (each costs ~4k density samples, a few ms):
    // this job waits for them, so nothing else touches the site meanwhile.
    struct Work { Site* s; Brick* b; const BodyRef* ref; uint8_t faces = 0; };
    std::vector<Work> work;
    for (auto& s : sites) {
        auto rit = refs.find(s->body);
        if (rit == refs.end()) continue;
        while (!s->pending.empty() && (int)work.size() < kInitPerStep) {
            Brick* b = s->pending.front(); s->pending.pop_front();
            if (b->ready) continue;
            bool dup = false;
            for (const Work& w : work) dup = dup || w.b == b;
            if (!dup) work.push_back({ s.get(), b, &rit->second, 0 });
        }
    }
    if (!work.empty()) {
        std::atomic<int> left{ (int)work.size() - 1 };
        for (size_t k = 1; k < work.size(); ++k) {
            Work* w = &work[k];
            WorldJobs().submit(-3.0f, [w, &left] { Classify(*w->s, *w->b, *w->ref, w->faces); left.fetch_sub(1); });
        }
        Classify(*work[0].s, *work[0].b, *work[0].ref, work[0].faces); // this thread takes one too
        while (left.load() > 0) std::this_thread::yield();
        static const int d[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
        for (Work& w : work) {
            for (Brick* o : w.b->nbr) if (o) Touch(o); // neighbours wake: water may now flow in or out
            // A pond's water at a face: the brick next to it holds more of the pond (#360: it joins whole).
            for (int f = 0; f < 6; ++f)
                if ((w.faces & (1 << f)) && !w.s->ensure(w.b->bi + d[f][0], w.b->bj + d[f][1], w.b->bk + d[f][2])) PinPondFace(*w.b, f);
        }
    }
    // CPU steps: kSubsteps per 50 ms on average (spring rates assume it), whatever the tick (#366: ticks every 25 ms).
    static double stepAcc = 0.0;
    stepAcc += kSubsteps * (kTickSeconds / 0.05);
    int n = (int)stepAcc; stepAcc -= n;
    n *= std::max(1, speed.load());
    for (int k = 0; k < n; ++k) {
        ++stepCount;
        if (gpuMode) { for (auto& s : sites) for (Site::Spring& sp : s->springs) StepSpring(*s, sp); } // the GPU runs the cells (GpuTick)
        else for (auto& s : sites) StepSite(*s, (int)(stepCount & 3));
        if (stepCount % kLevelEvery == 0) for (auto& s : sites) LevelSite(*s);
    }
    const auto t1 = std::chrono::steady_clock::now();
    // Publish: a new snapshot of the water field (waterfield.h) sharing every unchanged brick; the changed bricks'
    // bounds go out with it, so the chunks there redraw their water (the chunk mesher draws all fluid, #360).
    Fluid::Stats st;
    const WaterField::SnapshotPtr prev = WaterField::Current();
    std::shared_ptr<WaterField::Snapshot> snap;
    std::vector<WaterField::Dirty> dirty;
    for (auto& s : sites) {
        ++st.sites;
        for (auto& kv : s->bricks) {
            Brick& b = *kv.second;
            ++st.bricks;
            if (!b.ready) continue;
            if (b.active) ++st.activeBricks;
            if (b.basin) ++st.basinBricks;
            uint8_t present = 0; // kinds in the brick (bit k)
            for (int i = 0; i < B3; ++i)
                if (b.fill[i] && !(b.flags[i] & kSolid)) {
                    const int k = CellKind(b.flags[i]);
                    present |= (uint8_t)(1 << k);
                    if (b.flags[i] & kSource) continue;
                    st.wetCells += 1.0;
                    if (k) st.otherVolume += b.fill[i] / 255.0;
                    else st.volume += b.fill[i] / 255.0;
                }
            // Resting water sleeps (#366): an awake brick whose cells all stayed within kSettleDiff of where they were
            // kSettleSteps ago is put to sleep, though films still trickle back and forth in it (a drained pond kept
            // 30-50 bricks awake for minutes, and their changes kept the chunks redrawing, #381). Water moving into it
            // wakes it again (Touch); a real flow changes some cell by more than that within the window.
            if (b.active) {
                if (b.checkStep < 0) { b.check = b.fill; b.checkStep = stepCount; }
                else if (stepCount - b.checkStep >= kSettleSteps) {
                    int maxd = 0;
                    for (int i = 0; i < B3 && maxd <= kSettleDiff; ++i) maxd = std::max(maxd, abs((int)b.fill[i] - (int)b.check[i]));
                    if (maxd <= kSettleDiff) { b.active = false; b.idle = kSleepSteps + 1; b.settled = true; ++settledTotal; }
                    b.check = b.fill; b.checkStep = stepCount;
                }
            } else b.checkStep = -1;
            if (!b.changed) continue;
            b.changed = false;
            if (!snap) {
                snap = std::make_shared<WaterField::Snapshot>();
                snap->seq = prev->seq + 1;
                snap->sites.resize(sites.size());
                for (size_t k = 0; k < sites.size() && k < prev->sites.size(); ++k) snap->sites[k] = prev->sites[k]; // shares the bricks
                for (auto& o : sites) {
                    WaterField::Site& ps = snap->sites[o->index];
                    ps.body = o->body; ps.origin = o->origin; ps.e = o->e; ps.n = o->n; ps.u = o->u;
                }
            }
            WaterField::Site& ps = snap->sites[s->index];
            auto pb = std::make_shared<WaterField::Brick>();
            pb->fill = b.fill; pb->flags = b.flags; pb->kinds = present;
            ps.bricks[kv.first] = std::move(pb);
            Vector3 lo, hi;
            ps.brickBounds(b.bi, b.bj, b.bk, lo, hi);
            ps.lo = Vector3Min(ps.lo, lo); ps.hi = Vector3Max(ps.hi, hi);
            // A cell's surface reaches its neighbours' samples (solid cells carry it, empty ones look below). Only a
            // visible change redraws the chunks, and only around the cells that changed: some cell 2.5 cm+ off what
            // they last drew, or its solidity or fluid kind (a body easing by a 255th at a time kept ~2000 chunk
            // redraws a second going, #360). Those cells are drawn anew; the rest keep what was drawn.
            // Per 8^3 octant of the brick (#381: one box over scattered changes marked a whole brick's chunks).
            int c0[8][3], c1[8][3], mag[8] = {}; // mag: the largest change in each (coarse chunks wait for a bigger one, world.cpp)
            for (int o = 0; o < 8; ++o) { c0[o][0] = c0[o][1] = c0[o][2] = B; c1[o][0] = c1[o][1] = c1[o][2] = -1; }
            for (int i = 0; i < B3; ++i) {
                const bool v = !b.drawnSet || abs((int)b.fill[i] - (int)b.drawn[i]) >= 6 || ((b.flags[i] ^ b.drawnFlags[i]) & kSolid) ||
                               (b.fill[i] && ((b.flags[i] ^ b.drawnFlags[i]) & kKindMask));
                if (!v) continue;
                const int x = i % B, y = (i / B) % B, z = i / (B * B);
                const int o = (x >> 3) | ((y >> 3) << 1) | ((z >> 3) << 2);
                mag[o] = std::max(mag[o], (!b.drawnSet || ((b.flags[i] ^ b.drawnFlags[i]) & (kSolid | kKindMask))) ? 255 : abs((int)b.fill[i] - (int)b.drawn[i]));
                c0[o][0] = std::min(c0[o][0], x); c0[o][1] = std::min(c0[o][1], y); c0[o][2] = std::min(c0[o][2], z);
                c1[o][0] = std::max(c1[o][0], x); c1[o][1] = std::max(c1[o][1], y); c1[o][2] = std::max(c1[o][2], z);
                b.drawn[i] = b.fill[i]; b.drawnFlags[i] = b.flags[i];
            }
            b.drawnSet = true;
            for (int o = 0; o < 8; ++o) {
                if (c1[o][0] < 0) continue;
                Vector3 dlo, dhi;
                ps.cellBounds(b.bi * B + c0[o][0] - 1, b.bj * B + c0[o][1] - 1, b.bk * B + c0[o][2] - 1, b.bi * B + c1[o][0] + 2, b.bj * B + c1[o][1] + 2, b.bk * B + c1[o][2] + 2, dlo, dhi);
                dirty.push_back({ s->body, dlo, dhi, mag[o] });
            }
        }
    }
    if (snap) WaterField::Publish(std::move(snap), std::move(dirty));
    const auto t2 = std::chrono::steady_clock::now();
    {
        std::lock_guard<std::mutex> p(pubMx);
        const double stepMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
        st.stepMs = stepMs;
        st.maxStepMs = std::max(pubStats.maxStepMs, stepMs);
        st.meshMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
        st.steps = stepCount;
        static long long tot[4] = {}; tot[0] += dbgDown; tot[1] += dbgSide; tot[2] += dbgSrc; tot[3] += dbgLevel;
        for (int k = 0; k < 4; ++k) st.dbgTotal[k] = tot[k];
        st.wakes[0] = wakeGpu; st.wakes[1] = wakeTouch;
        st.pinnedFaces = pinnedFaces;
        st.dbg[0] = dbgDown; st.dbg[1] = dbgSide; st.dbg[2] = dbgSrc; st.dbg[3] = dbgLevel; dbgDown = dbgSide = dbgSrc = dbgLevel = 0;
        st.settled = (int)settledTotal;
        pubStats = st;
    }
}

// ---------------------------------------------------------------------------------------------
// GPU mode: the cells live in a GPU brick pool; once a tick (main thread, no job running) the last
// tick's results come back for the bricks that changed, new and reclassified bricks and the CPU's deltas
// go up, and the next substeps run as compute dispatches (resources/compute/fluid_step.comp).
// ---------------------------------------------------------------------------------------------

struct GpuState {
    bool init = false, ok = false;
    ComputeProgram step, delta;
    ComputeBuffer a, b, nbrs, list, inS, flags, deltas;
    int capacity = 0;                      // bricks
    std::vector<Brick*> owner;             // slot -> brick
    std::vector<Site*> ownerSite;
    std::vector<int> lastList;             // slots processed by the last dispatch
    int lastN = 0;                         // substeps in it
    bool nbrDirty = false;
} G;

constexpr unsigned kCellsPerBrick = (unsigned)B3;

bool GpuInit() {
    if (G.init) return G.ok;
    G.init = true;
    G.ok = Compute_Available() && G.step.load("resources/compute/fluid_step.comp") && G.delta.load("resources/compute/fluid_delta.comp");
    if (G.ok) {
        G.list.create(64 * 4); G.deltas.create(64 * 16); G.inS.create(64 * 4); G.flags.create(64 * 4); G.nbrs.create(64 * 24);
        G.a.create(64 * kCellsPerBrick * 4); G.b.create(64 * kCellsPerBrick * 4);
        G.capacity = 64;
        G.owner.assign(64, nullptr); G.ownerSite.assign(64, nullptr);
        const std::vector<unsigned> zeros(64, 0u);
        G.flags.upload(zeros.data(), 64 * 4);
    }
    TraceLog(LOG_WARNING, "Fluid: %s", G.ok ? "GPU compute" : "CPU (no compute)");
    return G.ok;
}

void GpuGrow(int need) {
    if (need <= G.capacity) return;
    int cap = G.capacity;
    while (cap < need) cap *= 2;
    auto grow = [](ComputeBuffer& buf, unsigned bytesOld, unsigned bytesNew) {
        ComputeBuffer nb; nb.create(bytesNew);
        rlCopyShaderBuffer(nb.id, buf.id, 0, 0, bytesOld);
        buf.destroy(); buf = nb;
    };
    grow(G.a, G.capacity * kCellsPerBrick * 4, cap * kCellsPerBrick * 4);
    grow(G.b, G.capacity * kCellsPerBrick * 4, cap * kCellsPerBrick * 4);
    grow(G.flags, G.capacity * 4, cap * 4);
    G.inS.resize(cap * 4);
    G.nbrs.resize(cap * 24);
    const std::vector<unsigned> zeros(cap - G.capacity, 0u);
    G.flags.upload(zeros.data(), (unsigned)zeros.size() * 4, G.capacity * 4);
    G.owner.resize(cap, nullptr); G.ownerSite.resize(cap, nullptr);
    G.capacity = cap;
    G.nbrDirty = true;
}

void GpuUploadBrick(const Brick& b) {
    static std::vector<unsigned> packed(kCellsPerBrick);
    for (int i = 0; i < B3; ++i) packed[i] = b.fill[i] | ((unsigned)b.flags[i] << 8) | ((unsigned)b.target[i] << 16);
    G.a.upload(packed.data(), kCellsPerBrick * 4, b.slot * kCellsPerBrick * 4);
    G.b.upload(packed.data(), kCellsPerBrick * 4, b.slot * kCellsPerBrick * 4);
}

void GpuTick(int n) {
    // 1. Results of the last dispatch: which bricks changed (and want neighbours), then their cells.
    if (!G.lastList.empty()) {
        std::vector<unsigned> flags(G.capacity);
        G.flags.download(flags.data(), G.capacity * 4);
        static std::vector<unsigned> cells(kCellsPerBrick);
        static const int d[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
        for (int slot : G.lastList) {
            Brick* br = G.owner[slot];
            Site* site = G.ownerSite[slot];
            if (!br || !site) continue;
            const unsigned fl = flags[slot];
            // A brick the CPU reclassified since that dispatch (upload pending) or is reclassifying now (not ready:
            // a worker reads its fill) keeps the CPU state: the GPU's older fill would put water back into cells
            // that became solid (PR #129 review).
            if ((fl & 1u) && br->ready && !br->upload) {
                G.a.download(cells.data(), kCellsPerBrick * 4, slot * kCellsPerBrick * 4);
                int drift = 0; // how far a settled brick's cells moved from where it settled
                for (int i = 0; i < B3; ++i) {
                    br->fill[i] = (uint8_t)(cells[i] & 255u);
                    br->flags[i] = (uint8_t)((br->flags[i] & ~kKindMask) | (((cells[i] >> 8) & 255u) & kKindMask));
                    if (br->settled) drift = std::max(drift, abs((int)br->fill[i] - (int)br->check[i]));
                }
                if (br->settled && drift <= kSettleDiff) { br->changed = true; continue; } // a resting brick at the edge of a dispatch: stays asleep (#366)
                if (br->settled) { ++wakeGpu; br->settled = false; }
                Touch(br);
                for (Brick* o : br->nbr) if (o) { o->active = true; o->idle = 0; }
            } else if (br->active && (br->idle += G.lastN) > kSleepSteps) {
                br->active = false;
            }
            for (int face = 0; face < 6; ++face)
                if (fl & (2u << face)) site->ensure(br->bi + d[face][0], br->bj + d[face][1], br->bk + d[face][2]);
        }
        std::fill(flags.begin(), flags.end(), 0u);
        G.flags.upload(flags.data(), G.capacity * 4);
        G.lastList.clear();
    }
    // 2. Slots for classified bricks; whole bricks up.
    int used = 0;
    for (Brick* o : G.owner) used += o ? 1 : 0;
    for (auto& s : sites)
        for (auto& kv : s->bricks) {
            Brick& br = *kv.second;
            if (!br.ready) continue;
            if (br.slot < 0) {
                GpuGrow(used + 1);
                int slot = 0;
                while (G.owner[slot]) ++slot;
                br.slot = slot; G.owner[slot] = &br; G.ownerSite[slot] = s.get();
                ++used;
                G.nbrDirty = true;
                br.upload = true;
            }
            if (br.upload) { GpuUploadBrick(br); br.upload = false; br.active = true; br.idle = 0; G.nbrDirty = true; }
        }
    if (G.nbrDirty) {
        std::vector<int> nb(G.capacity * 6, -1);
        for (int slot = 0; slot < G.capacity; ++slot)
            if (Brick* br = G.owner[slot])
                for (int f = 0; f < 6; ++f) nb[slot * 6 + f] = (br->nbr[f] && br->nbr[f]->ready && br->nbr[f]->slot >= 0) ? br->nbr[f]->slot : -1;
        G.nbrs.upload(nb.data(), G.capacity * 24);
        G.nbrDirty = false;
    }
    // 3. The CPU's deltas (summed per cell).
    if (!gpuDeltas.empty()) {
        struct D { unsigned cell; int delta; unsigned kind; unsigned pad; };
        std::unordered_map<unsigned, D> merged;
        for (const GpuDeltaRec& r : gpuDeltas) {
            if (r.b->slot < 0) continue;
            const unsigned cell = (unsigned)r.b->slot * kCellsPerBrick + (unsigned)r.i;
            D& dd = merged[cell];
            dd.cell = cell; dd.delta += r.d; dd.kind = r.kind;
        }
        gpuDeltas.clear();
        std::vector<D> list;
        for (auto& kv : merged) if (kv.second.delta) list.push_back(kv.second);
        if (!list.empty()) {
            if (G.deltas.bytes < list.size() * sizeof(D)) G.deltas.resize((unsigned)(list.size() * 2 * sizeof(D)));
            G.deltas.upload(list.data(), (unsigned)(list.size() * sizeof(D)));
            G.delta.use();
            G.delta.setUint("count", (unsigned)list.size());
            G.a.bind(0); G.deltas.bind(1); G.flags.bind(2);
            G.delta.dispatch(((unsigned)list.size() + 63) / 64);
            for (const D& dd : list) if (Brick* br = G.owner[dd.cell / kCellsPerBrick]) { if (!br->settled) { br->active = true; br->idle = 0; } } // a settled brick took a CPU nudge (levelling) asleep
        }
    }
    // 4. The bricks to run: the awake ones and their neighbours (flux only between bricks both in the set).
    std::vector<unsigned> inS(G.capacity, 0u);
    std::vector<int> list;
    for (int slot = 0; slot < G.capacity; ++slot) {
        Brick* br = G.owner[slot];
        if (!br || !br->active) continue;
        if (!inS[slot]) { inS[slot] = 1; list.push_back(slot); }
        for (Brick* o : br->nbr) if (o && o->slot >= 0 && !inS[o->slot]) { inS[o->slot] = 1; list.push_back(o->slot); }
    }
    if (list.empty()) return;
    G.inS.upload(inS.data(), G.capacity * 4);
    if (G.list.bytes < list.size() * 4) G.list.resize((unsigned)(list.size() * 8));
    G.list.upload(list.data(), (unsigned)(list.size() * 4));
    // 5. n substeps: vertical A -> B, sideways B -> A.
    G.step.use();
    G.step.setUint("count", (unsigned)list.size());
    int visc[8], spread[8];
    for (int k = 0; k < 8; ++k) { visc[k] = kindVisc[k]; spread[k] = kindMinSpread[k]; }
    rlSetUniform(G.step.uniform("visc"), visc, RL_SHADER_UNIFORM_INT, 8);
    rlSetUniform(G.step.uniform("minSpread"), spread, RL_SHADER_UNIFORM_INT, 8);
    G.nbrs.bind(2); G.list.bind(3); G.inS.bind(4); G.flags.bind(5);
    const unsigned bricks = (unsigned)list.size(); // 2D dispatch: 64 groups of 64 cells per brick, one row per brick
    for (int k = 0; k < n; ++k) {
        G.step.setInt("pass", 0); G.a.bind(0); G.b.bind(1); G.step.dispatch(kCellsPerBrick / 64u, bricks);
        G.step.setInt("pass", 1); G.b.bind(0); G.a.bind(1); G.step.dispatch(kCellsPerBrick / 64u, bricks);
    }
    G.lastList = list;
    G.lastN = n;
}

void GpuReset() {
    if (G.capacity > 0) { std::vector<unsigned> zeros(G.capacity, 0u); G.flags.upload(zeros.data(), G.capacity * 4); } // no stale face requests on reused slots
    for (Brick*& o : G.owner) o = nullptr;
    for (Site*& o : G.ownerSite) o = nullptr;
    G.lastList.clear();
    gpuDeltas.clear();
    G.nbrDirty = true;
}

} // namespace

// ---------------------------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------------------------

void Fluid::Reset() {
    WaitIdle();
    std::lock_guard<std::mutex> lock(simMx);
    GpuReset(); // the GPU pool forgets its bricks (before the sites holding them go)
    sites.clear();
    refs.clear();
    stepCount = 0;
    { std::lock_guard<std::mutex> c(cmdMx); commands.clear(); }
    {
        std::lock_guard<std::mutex> p(pubMx);
        pubStats = Stats{};
    }
    WaterField::Clear(); // only generated water until the simulation publishes again
    accum = 0.0;
}

void Fluid::OnTerrainEdit(Planetoid& body, Vector3 localCentre, float radius) {
    if (!dynamic_cast<const PlanetField*>(body.field.get())) return;
    std::lock_guard<std::mutex> c(cmdMx);
    commands.push_back({ 0, body.id, localCentre, radius, 0.0f });
}

float Fluid::Take(Planetoid& body, Vector3 local, float m3, float radius, int kind) {
    const float have = kind ? FluidAt(body, local, radius, kind) : WaterAt(body, local, radius);
    const float t = std::min(m3, have);
    if (t <= 0.0f) return 0.0f;
    std::lock_guard<std::mutex> c(cmdMx);
    commands.push_back({ 1, body.id, local, radius, t, kind });
    return t;
}

float Fluid::Pour(Planetoid& body, Vector3 local, float m3, int kind) {
    if (!dynamic_cast<const PlanetField*>(body.field.get())) return 0.0f;
    std::lock_guard<std::mutex> c(cmdMx);
    commands.push_back({ 2, body.id, local, 1.0f, m3, kind });
    return m3;
}

void Fluid::AddSource(Planetoid& body, Vector3 local, int kind, float rate, float head) {
    if (!dynamic_cast<const PlanetField*>(body.field.get()) || rate <= 0.0f) return;
    std::lock_guard<std::mutex> c(cmdMx);
    Command cm{ 3, body.id, local, rate, 0.0f, kind };
    cm.head = head;
    commands.push_back(cm);
}

int Fluid::SourceCount() {
    WaitIdle();
    std::lock_guard<std::mutex> lock(simMx);
    int n = 0;
    for (auto& s : sites) n += (int)s->springs.size();
    return n;
}

// fluids.txt: id | name | colour r,g,b | opacity | glow | viscosity | min spread | bucket item count
void Fluid::LoadKinds(const char* path) {
    kinds.assign(1, Kind{});
    FILE* f = fopen(path, "r");
    if (f) {
        char line[512];
        while (fgets(line, sizeof(line), f)) {
            if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
            char id[64] = {}, name[64] = {}, item[64] = {};
            int r = 0, g = 0, b = 0, visc = 5, minSpread = 10, count = 0;
            float opacity = 0.9f, glow = 0.0f;
            if (sscanf(line, " %63[^| ] | %63[^|] | %d , %d , %d | %f | %f | %d | %d | %63s %d", id, name, &r, &g, &b, &opacity, &glow, &visc, &minSpread, item, &count) < 9) continue;
            Kind k;
            k.id = id; k.name = name;
            while (!k.name.empty() && k.name.back() == ' ') k.name.pop_back();
            k.color = { (unsigned char)r, (unsigned char)g, (unsigned char)b, 255 };
            k.opacity = opacity; k.glow = glow; k.viscosity = std::max(1, visc); k.minSpread = std::max(1, minSpread);
            k.bucketItem = strcmp(item, "-") == 0 ? "" : item; k.bucketCount = count;
            if (k.id == "water") kinds[0] = k;
            else if (kinds.size() < 8) kinds.push_back(k);
        }
        fclose(f);
    }
    for (int i = 0; i < 8; ++i) {
        const Kind& k = kinds[std::min(i, (int)kinds.size() - 1)];
        // The GPU step conserves water only with viscosity >= 4 (four neighbours may each take diff / visc).
        if (k.viscosity < 4 && i < (int)kinds.size()) logger.logf("[Fluid] %s: viscosity %d below 4, using 4\n", k.id.c_str(), k.viscosity);
        kindVisc[i] = std::max(k.viscosity, 4); kindMinSpread[i] = k.minSpread;
    }
}
int Fluid::KindCount() { return (int)kinds.size(); }
const Fluid::Kind& Fluid::KindInfo(int kind) { return kinds[std::clamp(kind, 0, (int)kinds.size() - 1)]; }
int Fluid::KindByName(const std::string& id) {
    for (size_t i = 0; i < kinds.size(); ++i) if (kinds[i].id == id) return (int)i;
    return -1;
}

void Fluid::Update(Scene& world, float dt) {
    // A step every 50 ms when the last one has finished (the job owns the state while it runs).
    accum += dt;
    if (accum < kTickSeconds || running.load()) return;
    accum = std::min(accum - kTickSeconds, 0.1);
    bool work;
    { std::lock_guard<std::mutex> c(cmdMx); work = !commands.empty(); }
    {
        std::lock_guard<std::mutex> lock(simMx); // free: no job is running
        static bool decided = false;
        if (!decided) { decided = true; gpuMode = !getenv("FLUID_CPU") && GpuInit(); }
        if (gpuMode && !G.lastList.empty()) work = true; // the last dispatch's results still to read back
        for (auto& s : sites) {
            if (!s->pending.empty()) work = true;
            for (auto& kv : s->bricks) if (kv.second->active || kv.second->changed) { work = true; break; }
            if (work) break;
        }
        if (!work) return;
        // Snapshot the bodies the job may touch (fields are immutable; edits are copy-on-write).
        std::unordered_set<uint64_t> ids;
        for (auto& s : sites) ids.insert(s->body);
        { std::lock_guard<std::mutex> c(cmdMx); for (auto& cm : commands) ids.insert(cm.body); }
        refs.clear();
        for (uint64_t id : ids)
            if (Planetoid* p = Planetoid::find(id)) refs[id] = { p->field, p->edits, dynamic_cast<const PlanetField*>(p->field.get()) };
    }
    // GPU substeps: twice the CPU's - the order-independent GPU step moves water one cell a substep, where the
    // CPU's in-place sweep carried it several cells along the sweep.
    if (gpuMode) { std::lock_guard<std::mutex> lock(simMx); GpuTick((int)lround(2 * kSubsteps * (kTickSeconds / 0.05)) * std::max(1, speed.load())); } // 120 a second, as before
    (void)world;
    running = true;
    WorldJobs().submit(-0.5f, [] { RunJob(); running = false; });
}

void Fluid::BindKindUniforms(const Shader& shader) {
    // water.fs: per fluid kind (a chunk vertex's second texcoord) its colour (sRGB) + opacity, and glow (lava).
    float col[8 * 4] = {}, glow[8] = {};
    for (int k = 0; k < (int)kinds.size() && k < 8; ++k) {
        const Fluid::Kind& kd = kinds[k];
        col[k * 4 + 0] = kd.color.r / 255.0f; col[k * 4 + 1] = kd.color.g / 255.0f; col[k * 4 + 2] = kd.color.b / 255.0f;
        col[k * 4 + 3] = k == 0 ? 0.0f : std::max(kd.opacity, 0.01f); // water keeps its own look (a = 0)
        glow[k] = kd.glow;
    }
    // Locations looked up once per shader program (a hot reload makes a new id).
    static unsigned cachedId = 0; static int colLoc = -1, glowLoc = -1;
    if (shader.id != cachedId) { cachedId = shader.id; colLoc = GetShaderLocation(shader, "fluidColors"); glowLoc = GetShaderLocation(shader, "fluidGlows"); }
    if (colLoc >= 0) rlSetUniform(colLoc, col, RL_SHADER_UNIFORM_VEC4, 8);
    if (glowLoc >= 0) rlSetUniform(glowLoc, glow, RL_SHADER_UNIFORM_FLOAT, 8);
}

namespace {
    // A published cell at a body-local point: its fill (0-255) and flags; false if none.
    bool PubCell(const WaterField::Site& s, Vector3 local, int dk, uint8_t& fill, uint8_t& flags, float* frac = nullptr) {
        float fx, fy, fz;
        s.frame(local, fx, fy, fz);
        if (!s.cell((int)floorf(fx), (int)floorf(fy), (int)floorf(fz) + dk, fill, flags)) return false;
        if (frac) *frac = fz - floorf(fz);
        return true;
    }
}

bool Fluid::Wet(const Planetoid& body, Vector3 local) {
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        uint8_t f, fl; float frac;
        if (!PubCell(s, local, 0, f, fl, &frac) || (fl & kSolid)) continue;
        if (f / 255.0f > frac) return true;
        // A full cell above counts too (we are in the gap at a surface's edge).
    }
    return false;
}

float Fluid::SurfaceAbove(const Planetoid& body, Vector3 local, float reach) {
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        uint8_t f, fl; float frac;
        if (!PubCell(s, local, 0, f, fl, &frac) || f == 0) continue;
        float h = f / 255.0f - frac;
        for (int k = 1; k <= (int)reach && f == 255; ++k) {
            if (!PubCell(s, local, k, f, fl) || (fl & kSolid)) break;
            h += f / 255.0f;
        }
        return h;
    }
    return -1.0f;
}

float Fluid::WaterAt(const Planetoid& body, Vector3 local, float radius) {
    float total = 0.0f;
    {
        const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
            if (s.body != body.id) continue;
            const int R = (int)ceilf(radius);
            for (int dk = -R; dk <= R; ++dk)
                for (int dj = -R; dj <= R; ++dj)
                    for (int di = -R; di <= R; ++di) {
                        const Vector3 q = Vector3Add(local, Vector3Add(Vector3Add(Vector3Scale(s.e, (float)di), Vector3Scale(s.n, (float)dj)), Vector3Scale(s.u, (float)dk)));
                        uint8_t f, fl;
                        if (PubCell(s, q, 0, f, fl) && CellKind(fl) == 0) total += (fl & kSource) ? 100.0f : f / 255.0f;
                    }
        }
    }
    // Generated water where no simulated cell decides (#360: a drained pond's cells are simulated and dry).
    if (total <= 0.0f && WaterField::CellFill(*WaterField::Current(), body.id, local, 0) < 0) {
        const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
        if (pf) {
            const float d = Vector3Length(local);
            const Vector3 dir = Vector3Scale(local, 1.0f / std::max(d, 1.0f));
            float level = pf->seaLevelRadius();
            const PlanetField::WaterSource ws = pf->waterSourceAt(dir);
            if (ws.kind) level = std::max(level, pf->seaLevelRadius() + ws.level);
            if (d - radius < level && pf->radius + pf->elevation(dir) < level) total = ws.kind == 3 ? 10.0f : 100.0f;
        }
    }
    return total;
}

namespace {
    // Moving fluid within radius, m^3 per kind.
    void FluidVolumes(const Planetoid& body, Vector3 local, float radius, float out[8]) {
        for (int k = 0; k < 8; ++k) out[k] = 0.0f;
        const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
            if (s.body != body.id) continue;
            const int R = (int)ceilf(radius);
            for (int dk = -R; dk <= R; ++dk)
                for (int dj = -R; dj <= R; ++dj)
                    for (int di = -R; di <= R; ++di) {
                        const Vector3 q = Vector3Add(local, Vector3Add(Vector3Add(Vector3Scale(s.e, (float)di), Vector3Scale(s.n, (float)dj)), Vector3Scale(s.u, (float)dk)));
                        uint8_t f, fl;
                        if (PubCell(s, q, 0, f, fl) && !(fl & (kSolid | kSource))) out[CellKind(fl)] += f / 255.0f;
                    }
        }
    }
}

float Fluid::FluidAt(const Planetoid& body, Vector3 local, float radius, int kind) {
    float v[8];
    FluidVolumes(body, local, radius, v);
    return v[std::clamp(kind, 0, 7)];
}

int Fluid::KindAt(const Planetoid& body, Vector3 local, float radius) {
    float v[8];
    FluidVolumes(body, local, radius, v);
    int best = -1;
    for (int k = 1; k < 8; ++k) if (v[k] > 0.05f && (best < 0 || v[k] > v[best])) best = k;
    return best;
}


bool Fluid::StandingWaterNear(const Planetoid& body, Vector3 local, float range) {
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        const int R = (int)ceilf(range);
        for (int dk = -3; dk <= 1; ++dk)
            for (int dj = -R; dj <= R; dj += 1)
                for (int di = -R; di <= R; di += 1) {
                    const Vector3 q = Vector3Add(local, Vector3Add(Vector3Add(Vector3Scale(s.e, (float)di), Vector3Scale(s.n, (float)dj)), Vector3Scale(s.u, (float)dk)));
                    uint8_t f, fl;
                    if (PubCell(s, q, 0, f, fl) && f >= 64 && !(fl & kSolid)) return true;
                }
    }
    return false;
}

void Fluid::Save(FILE* f, uint64_t bodyId) {
    WaitIdle();
    std::lock_guard<std::mutex> lock(simMx);
    std::vector<Site*> mine;
    for (auto& s : sites) if (s->body == bodyId) mine.push_back(s.get());
    fprintf(f, "water %zu\n", mine.size());
    for (Site* s : mine) {
        std::vector<Brick*> wet;
        for (auto& kv : s->bricks) {
            Brick& b = *kv.second;
            bool any = false;
            for (int i = 0; i < B3 && !any; ++i) any = b.fill[i] && !(b.flags[i] & kSource);
            if (any || b.basin) wet.push_back(&b); // a drained pond's bricks too: loaded fresh they would hold its water again (#360)
        }
        fprintf(f, "site %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f %.6f %zu %zu", s->origin.x, s->origin.y, s->origin.z, s->e.x, s->e.y, s->e.z, s->u.x, s->u.y, s->u.z, (size_t)0, wet.size()); // 0: claimed ponds are gone (#360)
        fprintf(f, "\n");
        fprintf(f, "springs %zu\n", s->springs.size());
        for (const Site::Spring& sp : s->springs) fprintf(f, "%d %d %d %d %.4f %.3f\n", sp.ci, sp.cj, sp.ck, sp.kind, sp.perStep, sp.head);
        for (Brick* b : wet) {
            // Run-length fill (sources and solid cells as 0: they are re-derived on load).
            fprintf(f, "%d %d %d", b->bi, b->bj, b->bk);
            int i = 0;
            while (i < B3) {
                // A cell's fill (0-255) + 256 x its fluid kind.
            auto code = [&](int c) { return (b->flags[c] & (kSource | kSolid)) ? 0 : b->fill[c] + (b->fill[c] ? 256 * CellKind(b->flags[c]) : 0); };
            const int v = code(i);
                int run = 1;
                while (i + run < B3 && code(i + run) == v) ++run;
                fprintf(f, " %d*%d", v, run);
                i += run;
            }
            fprintf(f, "\n");
        }
    }
}

void Fluid::Load(FILE* f, uint64_t bodyId) {
    int n = 0;
    if (fscanf(f, " water %d", &n) != 1) return;
    WaitIdle();
    std::lock_guard<std::mutex> lock(simMx);
    for (int k = 0; k < n; ++k) {
        auto s = std::make_unique<Site>();
        size_t np = 0, nb = 0;
        if (fscanf(f, " site %f %f %f %f %f %f %f %f %f %zu %zu", &s->origin.x, &s->origin.y, &s->origin.z, &s->e.x, &s->e.y, &s->e.z, &s->u.x, &s->u.y, &s->u.z, &np, &nb) != 11) return;
        s->n = Vector3CrossProduct(s->u, s->e);
        s->body = bodyId;
        s->index = (int)sites.size();
        for (size_t i = 0; i < np; ++i) {
            unsigned long long id = 0;
            if (fscanf(f, " %llu", &id) != 1) break; // old saves: ponds taken over (#360: no more claims; their drained bricks were not stored: the pond refills, pulled in whole)
        }
        int ns = 0; // sources (springs, seeps): "springs N" then one per line
        if (fscanf(f, " springs %d", &ns) == 1)
            for (int j = 0; j < ns; ++j) {
                Site::Spring sp{};
                if (fscanf(f, " %d %d %d %d %f %f", &sp.ci, &sp.cj, &sp.ck, &sp.kind, &sp.perStep, &sp.head) == 6) s->springs.push_back(sp);
            }
        for (size_t j = 0; j < nb; ++j) {
            int bi, bj, bk;
            if (fscanf(f, " %d %d %d", &bi, &bj, &bk) != 3) return;
            Brick* b = s->ensure(bi, bj, bk);
            int i = 0;
            while (i < B3) {
                int v, run;
                if (fscanf(f, " %d*%d", &v, &run) != 2) return;
                for (int r = 0; r < run && i < B3; ++r) {
                    if (b) { b->fill[i] = (uint8_t)(v & 255); b->flags[i] = (uint8_t)((b->flags[i] & ~kKindMask) | ((v >> 8) << kKindShift)); }
                    ++i;
                }
            }
            if (b) b->keepFill = true;
        }
        sites.push_back(std::move(s));
    }
}

Fluid::Stats Fluid::GetStats() { std::lock_guard<std::mutex> p(pubMx); return pubStats; }

void Fluid::WaitIdle() {
    while (running.load()) std::this_thread::yield();
}

void Fluid::SetSpeed(int steps) { speed = std::max(1, steps); }

bool Fluid::Asleep() {
    if (running.load()) return false;
    { std::lock_guard<std::mutex> c(cmdMx); if (!commands.empty()) return false; }
    std::lock_guard<std::mutex> lock(simMx);
    for (auto& s : sites) {
        if (!s->pending.empty()) return false;
        for (auto& kv : s->bricks) if (kv.second->active) return false;
    }
    return true;
}

std::string Fluid::DebugCell(const Planetoid& body, Vector3 local) {
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        uint8_t f, fl;
        if (!PubCell(s, local, 0, f, fl)) return "no brick";
        return TextFormat("fill %d flags %s%s%s%s", f, (fl & kSolid) ? "S" : "", (fl & kSource) ? "R" : "", (fl & kStatic) ? "T" : "", (fl & kDug) ? "D" : "");
    }
    return "no site";
}

bool Fluid::Simulated(const Planetoid& body, Vector3 local) {
    return WaterField::CellFill(*WaterField::Current(), body.id, local, 0) >= 0;
}

bool Fluid::GeneratedWaterValid(uint64_t body, Vector3 local) {
    const int f = WaterField::CellFill(*WaterField::Current(), body, local, 0);
    return f < 0 || f > 0;
}
