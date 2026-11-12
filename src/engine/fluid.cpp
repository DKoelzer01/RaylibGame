#include "fluid.h"
#include "cubemarch.h"
#include "jobsystem.h"
#include "object.h"
#include "terraindeltas.h"
#include "voxel.h"
#include "world.h"
#include "worldtypes.h"

#include <raymath.h>
#include <rlgl.h>
#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
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
constexpr uint8_t kSolid = 1, kSource = 2, kStatic = 4; // kStatic: inside static water (the chunks draw its surface)
constexpr uint8_t kDug = 8;            // open only because of an edit: static water may flow into it
constexpr int kSiteBricks = 16;        // a site reaches +-16 bricks (256 m) from its origin
constexpr float kSiteJoin = 180.0f;    // a disturbance within this of a site's origin joins it
constexpr int kMinSpread = 10;         // films thinner than this (~4 cm) don't spread sideways
constexpr int kEvaporate = 3;          // resting specks this thin (~1 cm) dry up
constexpr int kSubsteps = 3;           // simulation steps per 20 Hz tick (water moves a cell a step)
constexpr int kSleepSteps = 120;       // 2 s without change: a brick sleeps
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
    bool active = true;
    int idle = 0;
    bool meshDirty = true;
    int lastMesh = -100;            // step of the last mesh build
    bool published = false;
    bool changed = true;            // since the last publish
};

struct Site {
    uint64_t body = 0;
    int index = 0;
    Vector3 origin{}, e{}, n{}, u{};
    std::unordered_map<int64_t, std::unique_ptr<Brick>> bricks;
    std::deque<Brick*> pending;     // waiting for classification
    std::vector<uint64_t> ponds;    // ponds this site took over
    struct PondInfo { Vector3 dir; float r, level; };
    std::vector<PondInfo> pondInfo; // their extent (this session: a load keeps the stored water instead)

    Vector3 cellCentre(int ci, int cj, int ck) const {
        return Vector3Add(origin, Vector3Add(Vector3Add(Vector3Scale(e, ci + 0.5f), Vector3Scale(n, cj + 0.5f)), Vector3Scale(u, ck + 0.5f)));
    }
    void cellOf(Vector3 local, int& ci, int& cj, int& ck) const {
        const Vector3 d = Vector3Subtract(local, origin);
        ci = (int)floorf(Vector3DotProduct(d, e)); cj = (int)floorf(Vector3DotProduct(d, n)); ck = (int)floorf(Vector3DotProduct(d, u));
    }
    bool inRange(int bi, int bj, int bk) const { return abs(bi) <= kSiteBricks && abs(bj) <= kSiteBricks && abs(bk) <= kSiteBricks; }
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

struct Command { int kind; uint64_t body; Vector3 c; float r, m3; }; // 0 edit, 1 take, 2 pour

struct MeshOut {
    int site; int64_t key; uint64_t body; DVec3 anchor;
    std::vector<float> pos, nrm, uv;
    std::vector<unsigned short> idx;
};

struct PubBrick { std::array<uint8_t, B3> fill{}, flags{}; };
struct PubSite {
    uint64_t body = 0; Vector3 origin{}, e{}, n{}, u{};
    std::unordered_map<int64_t, std::shared_ptr<PubBrick>> bricks;
};

struct GpuMesh { uint64_t body = 0; DVec3 anchor{}; Mesh mesh{}; };

// Simulation state: owned by the step job while it runs (simMx).
std::mutex simMx;
std::vector<std::unique_ptr<Site>> sites;
std::unordered_map<uint64_t, BodyRef> refs;     // snapshot handed to the job
long long stepCount = 0;
Fluid::Stats stats;
long long dbgDown = 0, dbgSide = 0, dbgSrc = 0, dbgLevel = 0; // tests: changes by kind (since the last publish)

// Main thread <-> job.
std::mutex cmdMx;
std::vector<Command> commands;
std::atomic<bool> running{ false };
std::atomic<int> speed{ 1 };
double accum = 0.0;

std::mutex pubMx;                                // published results
std::vector<PubSite> published;
std::vector<MeshOut> meshQueue;
std::vector<std::pair<int, int64_t>> meshDrops;  // bricks whose mesh became empty
std::vector<std::pair<uint64_t, Vector4>> remeshQueue; // claimed ponds: chunks to rebuild (body, local centre + radius)
Fluid::Stats pubStats;

std::map<std::pair<int, int64_t>, GpuMesh> gpu;  // main thread

// ---------------------------------------------------------------------------------------------
// Classification
// ---------------------------------------------------------------------------------------------

// Cells of a brick: solid from the terrain density (edits included); open cells inside static water
// are sources (infinite) or, for ponds, take the pond over.
struct PondClaim { uint64_t id; Vector3 dir; float r, level; };
// Runs in parallel over bricks (pool jobs): it writes only its own brick; ponds it finds are reported
// in `claims` and taken over after the batch (RunJob), never here.
void Classify(const Site& s, Brick& b, const BodyRef& ref, std::vector<PondClaim>& claims) {
    std::vector<Site::PondInfo> ponds = s.pondInfo; // taken over before, and found in this brick
    const PlanetField* pf = ref.pf;
    const float sea = pf ? pf->seaLevelRadius() : 0.0f;
    for (int y = 0; y < B; ++y)
        for (int x = 0; x < B; ++x) {
            // The column's ground and static water (from its middle: a column tilts < 0.3 deg in a site).
            const int ci = b.bi * B + x, cj = b.bj * B + y;
            const Vector3 mid = s.cellCentre(ci, cj, b.bk * B + B / 2);
            const float ml = Vector3Length(mid);
            if (ml < 1.0f) continue;
            const Vector3 dir = Vector3Scale(mid, 1.0f / ml);
            float bed = 0.0f, level = -1e9f;
            int kind = 0;                               // 1 river, 2 lake, 4 sea: infinite; 5 a pond taken over: finite
            if (pf) {
                bed = pf->radius + pf->elevation(dir);
                if (bed < sea) { level = sea; kind = 4; } // the sea
                const PlanetField::WaterSource ws = pf->waterSourceAt(dir);
                if (ws.kind == 3 && bed < sea + ws.level) {
                    // A pond: the simulation takes it over (finite water from here on); reported to RunJob,
                    // which hides its static surface and queues every brick of its water.
                    bool known = false;
                    for (const PondClaim& c : claims) known = known || c.id == ws.pond;
                    if (!known) { claims.push_back({ ws.pond, ws.pondDir, ws.pondR, ws.level }); ponds.push_back({ ws.pondDir, ws.pondR, ws.level }); }
                } else if (ws.kind && sea + ws.level > level) { level = sea + ws.level; kind = ws.kind; }
                for (const Site::PondInfo& p : ponds)
                    if (Vector3Distance(dir, p.dir) * pf->radius < p.r * 1.05f && sea + p.level > level) { level = sea + p.level; kind = 5; }
            }
            const float top = pf ? bed + 14.0f : 1e9f;  // tors, boulders and detail stay below this
            for (int z = 0; z < B; ++z) {
                const int i = Idx(x, y, z);
                const Vector3 c = s.cellCentre(ci, cj, b.bk * B + z);
                const float r = Vector3Length(c);
                const bool solid = r < top && ref.field->density(c, ref.edits.get()) >= ISO_LEVEL;
                uint8_t fl = solid ? kSolid : 0;
                // Dug: open only because of an edit. Static water feeds only these (the real ground has
                // crags the height map doesn't: natural dips by a river are not to be flooded).
                if (!solid && r < bed + 4.0f && r < top && ref.field->density(c, nullptr) >= ISO_LEVEL) fl |= kDug;
                uint8_t tg = 0;
                if (!solid && kind && r - 0.5f < level && r > bed - 0.5f) {
                    tg = (uint8_t)std::clamp((int)lroundf((level - (r - 0.5f)) * 255.0f), 0, 255);
                    if (tg > 0 && kind != 5) fl |= kSource | kStatic;
                }
                b.flags[i] = fl;
                b.target[i] = (fl & kSource) ? tg : 0;
                if (b.keepFill) { if (solid) b.fill[i] = 0; else if (fl & kSource) b.fill[i] = tg; }
                else b.fill[i] = solid ? 0 : tg;
            }
        }
    b.ready = true;
    b.keepFill = false;
    b.active = true; b.idle = 0; b.meshDirty = true; b.changed = true;
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

void Touch(Brick* b) { b->active = true; b->idle = 0; b->meshDirty = true; b->changed = true; }

// Water moving out of a cell into a missing brick: queue that brick (it takes water from the next
// step on, once classified).
void RequestBeyond(Site& s, Brick* b, int x, int y, int z, int dx, int dy, int dz) {
    const int ci = b->bi * B + x + dx, cj = b->bj * B + y + dy, ck = b->bk * B + z + dz;
    s.ensure(FloorDiv(ci, B), FloorDiv(cj, B), FloorDiv(ck, B));
}

void StepSite(Site& s, int parity) {
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
                            const int room = 255 - d.b->fill[d.i];
                            if (room > 0 && !(d.b->flags[d.i] & kSource) && (!src || (d.b->flags[d.i] & kDug))) { // static water takes nothing in
                                const int m = std::min(f, room);
                                d.b->fill[d.i] = (uint8_t)(d.b->fill[d.i] + m);
                                f -= m;
                                if (d.b != b) Touch(d.b);
                                moved = true; ++dbgDown;
                            }
                            belowFull = d.b->fill[d.i] == 255;
                        }
                    } else if (z == 0 && !b->nbr[4]) {
                        if (!src) RequestBeyond(s, b, x, y, z, 0, 0, -1); // static water never asks for more bricks
                        belowOpen = true; belowFull = false; // wait for it
                    }
                    // Sideways, once the cell rests on something.
                    if (f >= kMinSpread && (!belowOpen || belowFull)) {
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
                            const int diff = f0 - o.b->fill[o.i];
                            if (diff < 5) continue;
                            const int m = std::min(diff / 5, f);
                            if (o.b->fill[o.i] == 0 && m < kMinSpread) continue; // no trickles into dry cells: a film's edge holds
                            if (m <= 0) continue;
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
                if (nf != c.b->fill[c.i]) { c.b->fill[c.i] = (uint8_t)nf; Touch(c.b); ++dbgLevel; }
            }
        }
    }
}

// ---------------------------------------------------------------------------------------------
// Meshing
// ---------------------------------------------------------------------------------------------

void MeshBrick(Site& s, Brick& b, MeshOut& out) {
    constexpr int S = B + 1;
    std::vector<float> D((size_t)S * S * S, -0.5f);
    std::vector<uint8_t> stat((size_t)S * S * S, 0);
    auto at = [&](int x, int y, int z, float& d, uint8_t& fl) {
        Brick* bb; int i;
        if (!s.cell(b.bi * B + x, b.bj * B + y, b.bk * B + z, bb, i)) { d = -0.5f; fl = 0; return false; }
        fl = bb->flags[i];
        d = bb->fill[i] / 255.0f - 0.5f;
        return true;
    };
    bool any = false;
    for (int z = 0; z < S; ++z)
        for (int y = 0; y < S; ++y)
            for (int x = 0; x < S; ++x) {
                float d; uint8_t fl;
                at(x, y, z, d, fl);
                if (fl & kSolid) {
                    // Carry the water surface into the rock so it meets the terrain without a gap.
                    float best = -0.5f;
                    const int side[4][2] = { { 1, 0 }, { -1, 0 }, { 0, 1 }, { 0, -1 } };
                    for (auto& sd : side) { float od; uint8_t ofl; if (at(x + sd[0], y + sd[1], z, od, ofl) && !(ofl & kSolid)) best = std::max(best, od); }
                    d = best;
                }
                const size_t k = ((size_t)z * S + y) * S + x;
                D[k] = d;
                stat[k] = (fl & kStatic) ? 1 : 0;
                if (d > 0.0f && !(fl & kStatic) && !(fl & kSolid)) any = true;
            }
    out.site = s.index; out.key = Key(b.bi, b.bj, b.bk); out.body = s.body;
    const Vector3 o = s.cellCentre(b.bi * B, b.bj * B, b.bk * B);
    out.anchor = DVec3(o);
    if (!any) return;
    std::vector<Vector3> verts;
    std::vector<int> idx;
    std::vector<int> cache((size_t)S * S * S * 3, -1);
    for (int z = 0; z < B; ++z)
        for (int y = 0; y < B; ++y)
            for (int x = 0; x < B; ++x) marchCube(x, y, z, D, S, verts, idx, cache, 0.0f);
    if (idx.empty() || verts.size() > 65000) return;
    auto sample = [&](Vector3 g) { // trilinear D at grid coordinates
        const int x0 = std::clamp((int)floorf(g.x), 0, S - 2), y0 = std::clamp((int)floorf(g.y), 0, S - 2), z0 = std::clamp((int)floorf(g.z), 0, S - 2);
        const float fx = g.x - x0, fy = g.y - y0, fz = g.z - z0;
        auto v = [&](int a, int bq, int c) { return D[((size_t)(z0 + c) * S + (y0 + bq)) * S + (x0 + a)]; };
        const float c00 = v(0, 0, 0) * (1 - fx) + v(1, 0, 0) * fx, c10 = v(0, 1, 0) * (1 - fx) + v(1, 1, 0) * fx;
        const float c01 = v(0, 0, 1) * (1 - fx) + v(1, 0, 1) * fx, c11 = v(0, 1, 1) * (1 - fx) + v(1, 1, 1) * fx;
        return (c00 * (1 - fy) + c10 * fy) * (1 - fz) + (c01 * (1 - fy) + c11 * fy) * fz;
    };
    std::vector<int> remap(verts.size(), -1);
    std::vector<Vector3> nacc;
    for (size_t t = 0; t + 2 < idx.size(); t += 3) {
        Vector3 a = verts[idx[t]], bq = verts[idx[t + 1]], c = verts[idx[t + 2]];
        const Vector3 cen = Vector3Scale(Vector3Add(Vector3Add(a, bq), c), 1.0f / 3.0f);
        // Static water (the chunks draw its surface) and faces deep in rock are left out.
        const int sx = std::clamp((int)lroundf(cen.x), 0, S - 1), sy = std::clamp((int)lroundf(cen.y), 0, S - 1), sz = std::clamp((int)lroundf(cen.z), 0, S - 1);
        if (stat[((size_t)sz * S + sy) * S + sx]) continue;
        // Outward = down the fill gradient; fix the winding to match.
        const float h = 0.25f;
        const Vector3 grad = { sample({ cen.x + h, cen.y, cen.z }) - sample({ cen.x - h, cen.y, cen.z }),
                               sample({ cen.x, cen.y + h, cen.z }) - sample({ cen.x, cen.y - h, cen.z }),
                               sample({ cen.x, cen.y, cen.z + h }) - sample({ cen.x, cen.y, cen.z - h }) };
        Vector3 fn = Vector3CrossProduct(Vector3Subtract(bq, a), Vector3Subtract(c, a));
        int i0 = idx[t], i1 = idx[t + 1], i2 = idx[t + 2];
        if (Vector3DotProduct(fn, grad) > 0.0f) { std::swap(i1, i2); fn = Vector3Negate(fn); }
        const int tri[3] = { i0, i1, i2 };
        for (int e = 0; e < 3; ++e) {
            int& r = remap[tri[e]];
            if (r < 0) {
                r = (int)nacc.size();
                const Vector3 g = verts[tri[e]];
                const Vector3 lp = Vector3Add(Vector3Add(Vector3Scale(s.e, g.x), Vector3Scale(s.n, g.y)), Vector3Scale(s.u, g.z));
                out.pos.insert(out.pos.end(), { lp.x, lp.y, lp.z });
                // Depth for the shader's tint: water cells below this point in the column.
                int depth = 1;
                for (int z = (int)floorf(g.z); z >= 0 && depth < 8; --z) {
                    const int gx = std::clamp((int)lroundf(g.x), 0, S - 1), gy = std::clamp((int)lroundf(g.y), 0, S - 1);
                    if (D[((size_t)z * S + gy) * S + gx] > -0.4f) ++depth; else break;
                }
                out.uv.insert(out.uv.end(), { (float)depth, 0.0f });
                nacc.push_back({ 0, 0, 0 });
            }
            nacc[r] = Vector3Add(nacc[r], fn);
            out.idx.push_back((unsigned short)r);
        }
    }
    for (const Vector3& nv : nacc) {
        // Body-local normal from grid space.
        Vector3 g = Vector3Normalize(nv);
        const Vector3 l = Vector3Add(Vector3Add(Vector3Scale(s.e, g.x), Vector3Scale(s.n, g.y)), Vector3Scale(s.u, g.z));
        out.nrm.insert(out.nrm.end(), { l.x, l.y, l.z });
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
            const int add = (int)std::min(left, 255.0f - b->fill[i]);
            b->fill[i] = (uint8_t)(b->fill[i] + add);
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
                const int t = (int)std::min(want, (float)b->fill[i]);
                if (!(b->flags[i] & kSource)) b->fill[i] = (uint8_t)(b->fill[i] - t);
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
    // this job waits for them, so nothing else touches the site meanwhile. Ponds found are taken over
    // after the batch.
    struct Work { Site* s; Brick* b; const BodyRef* ref; std::vector<PondClaim> claims; };
    std::vector<Work> work;
    for (auto& s : sites) {
        auto rit = refs.find(s->body);
        if (rit == refs.end()) continue;
        while (!s->pending.empty() && (int)work.size() < kInitPerStep) {
            Brick* b = s->pending.front(); s->pending.pop_front();
            if (b->ready) continue;
            bool dup = false;
            for (const Work& w : work) dup = dup || w.b == b;
            if (!dup) work.push_back({ s.get(), b, &rit->second, {} });
        }
    }
    if (!work.empty()) {
        std::atomic<int> left{ (int)work.size() - 1 };
        for (size_t k = 1; k < work.size(); ++k) {
            Work* w = &work[k];
            WorldJobs().submit(-3.0f, [w, &left] { Classify(*w->s, *w->b, *w->ref, w->claims); left.fetch_sub(1); });
        }
        Classify(*work[0].s, *work[0].b, *work[0].ref, work[0].claims); // this thread takes one too
        while (left.load() > 0) std::this_thread::yield();
        for (Work& w : work) {
            for (Brick* o : w.b->nbr) if (o) Touch(o); // neighbours wake: water may now flow in or out
            const PlanetField* pf = w.ref->pf;
            for (const PondClaim& c : w.claims) {
                if (!pf || std::find(w.s->ponds.begin(), w.s->ponds.end(), c.id) != w.s->ponds.end()) continue;
                // Take the pond over: no static surface (the chunks around it rebuild); queue its bricks.
                pf->claimPond(c.id);
                w.s->ponds.push_back(c.id);
                w.s->pondInfo.push_back({ c.dir, c.r, c.level });
                const Vector3 pc = Vector3Scale(c.dir, pf->radius + c.level);
                { std::lock_guard<std::mutex> lock(pubMx); remeshQueue.push_back({ w.s->body, { pc.x, pc.y, pc.z, c.r * 1.3f + 8.0f } }); }
                int pci, pcj, pck;
                w.s->cellOf(pc, pci, pcj, pck);
                const int rc = (int)ceilf(c.r * 1.1f) + 1;
                for (int bk = FloorDiv(pck - 8, B); bk <= FloorDiv(pck, B); ++bk)
                    for (int bj = FloorDiv(pcj - rc, B); bj <= FloorDiv(pcj + rc, B); ++bj)
                        for (int bi = FloorDiv(pci - rc, B); bi <= FloorDiv(pci + rc, B); ++bi) w.s->ensure(bi, bj, bk);
            }
        }
    }
    const int n = kSubsteps * std::max(1, speed.load());
    for (int k = 0; k < n; ++k) {
        ++stepCount;
        for (auto& s : sites) StepSite(*s, (int)(stepCount & 3));
        if (stepCount % kLevelEvery == 0) for (auto& s : sites) LevelSite(*s);
    }
    const auto t1 = std::chrono::steady_clock::now();
    // Publish changed bricks; mesh the dirty ones (each at most every other step).
    std::vector<MeshOut> meshes;
    std::vector<std::pair<int, int64_t>> drops;
    Fluid::Stats st;
    for (auto& s : sites) {
        ++st.sites;
        for (auto& kv : s->bricks) {
            Brick& b = *kv.second;
            ++st.bricks;
            if (!b.ready) continue;
            if (b.active) ++st.activeBricks;
            for (int i = 0; i < B3; ++i) if (b.fill[i] && !(b.flags[i] & (kSource | kSolid))) { st.wetCells += 1.0; st.volume += b.fill[i] / 255.0; }
            if (b.meshDirty && stepCount - b.lastMesh >= 2) {
                MeshOut m;
                MeshBrick(*s, b, m);
                b.meshDirty = false; b.lastMesh = (int)stepCount;
                if (m.idx.empty()) drops.push_back({ s->index, kv.first });
                else meshes.push_back(std::move(m));
            }
        }
        st.claimedPonds += (int)s->ponds.size();
    }
    const auto t2 = std::chrono::steady_clock::now();
    {
        std::lock_guard<std::mutex> p(pubMx);
        published.resize(sites.size());
        for (auto& s : sites) {
            PubSite& ps = published[s->index];
            ps.body = s->body; ps.origin = s->origin; ps.e = s->e; ps.n = s->n; ps.u = s->u;
            for (auto& kv : s->bricks) {
                Brick& b = *kv.second;
                if (!b.ready || !b.changed) continue;
                b.changed = false;
                auto pb = std::make_shared<PubBrick>();
                pb->fill = b.fill; pb->flags = b.flags;
                ps.bricks[kv.first] = pb;
            }
        }
        for (auto& m : meshes) meshQueue.push_back(std::move(m));
        for (auto& d : drops) meshDrops.push_back(d);
        const double stepMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
        st.stepMs = stepMs;
        st.maxStepMs = std::max(pubStats.maxStepMs, stepMs);
        st.meshMs = std::chrono::duration<double, std::milli>(t2 - t1).count();
        st.steps = stepCount;
        st.dbg[0] = dbgDown; st.dbg[1] = dbgSide; st.dbg[2] = dbgSrc; st.dbg[3] = dbgLevel; dbgDown = dbgSide = dbgSrc = dbgLevel = 0;
        st.meshes = (int)gpu.size();
        pubStats = st;
    }
}

} // namespace

// ---------------------------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------------------------

void Fluid::Reset() {
    WaitIdle();
    std::lock_guard<std::mutex> lock(simMx);
    // Ponds taken over come back as static water (a load re-claims the ones it stored).
    for (auto& s : sites)
        if (Planetoid* p = Planetoid::find(s->body))
            if (const auto* pf = dynamic_cast<const PlanetField*>(p->field.get())) pf->clearClaimedPonds();
    sites.clear();
    refs.clear();
    stepCount = 0;
    { std::lock_guard<std::mutex> c(cmdMx); commands.clear(); }
    {
        std::lock_guard<std::mutex> p(pubMx);
        published.clear(); meshQueue.clear(); meshDrops.clear(); remeshQueue.clear();
        pubStats = Stats{};
    }
    for (auto& kv : gpu) if (kv.second.mesh.vaoId) UnloadMesh(kv.second.mesh);
    gpu.clear();
    accum = 0.0;
}

void Fluid::OnTerrainEdit(Planetoid& body, Vector3 localCentre, float radius) {
    if (!dynamic_cast<const PlanetField*>(body.field.get())) return;
    std::lock_guard<std::mutex> c(cmdMx);
    commands.push_back({ 0, body.id, localCentre, radius, 0.0f });
}

float Fluid::Take(Planetoid& body, Vector3 local, float m3, float radius) {
    const float have = WaterAt(body, local, radius);
    const float t = std::min(m3, have);
    if (t <= 0.0f) return 0.0f;
    std::lock_guard<std::mutex> c(cmdMx);
    commands.push_back({ 1, body.id, local, radius, t });
    return t;
}

float Fluid::Pour(Planetoid& body, Vector3 local, float m3) {
    if (!dynamic_cast<const PlanetField*>(body.field.get())) return 0.0f;
    std::lock_guard<std::mutex> c(cmdMx);
    commands.push_back({ 2, body.id, local, 1.0f, m3 });
    return m3;
}

void Fluid::Update(Scene& world, float dt) {
    // Uploads (main thread: GL).
    std::vector<MeshOut> meshes;
    std::vector<std::pair<int, int64_t>> drops;
    std::vector<std::pair<uint64_t, Vector4>> remesh;
    {
        std::lock_guard<std::mutex> p(pubMx);
        meshes.swap(meshQueue); drops.swap(meshDrops); remesh.swap(remeshQueue);
    }
    for (auto& d : drops) { auto it = gpu.find(d); if (it != gpu.end()) { if (it->second.mesh.vaoId) UnloadMesh(it->second.mesh); gpu.erase(it); } }
    for (auto& m : meshes) {
        GpuMesh& g = gpu[{ m.site, m.key }];
        if (g.mesh.vaoId) UnloadMesh(g.mesh);
        g = GpuMesh{};
        g.body = m.body; g.anchor = m.anchor;
        Mesh& mesh = g.mesh;
        mesh.vertexCount = (int)(m.pos.size() / 3);
        mesh.triangleCount = (int)(m.idx.size() / 3);
        mesh.vertices = (float*)MemAlloc((unsigned)(m.pos.size() * sizeof(float))); memcpy(mesh.vertices, m.pos.data(), m.pos.size() * sizeof(float));
        mesh.normals = (float*)MemAlloc((unsigned)(m.nrm.size() * sizeof(float))); memcpy(mesh.normals, m.nrm.data(), m.nrm.size() * sizeof(float));
        mesh.texcoords = (float*)MemAlloc((unsigned)(m.uv.size() * sizeof(float))); memcpy(mesh.texcoords, m.uv.data(), m.uv.size() * sizeof(float));
        mesh.indices = (unsigned short*)MemAlloc((unsigned)(m.idx.size() * sizeof(unsigned short))); memcpy(mesh.indices, m.idx.data(), m.idx.size() * sizeof(unsigned short));
        UploadMesh(&mesh, false);
    }
    for (auto& r : remesh)
        if (Planetoid* p = Planetoid::find(r.first)) RequestRemesh(*p, { r.second.x, r.second.y, r.second.z }, r.second.w);

    // A step every 50 ms when the last one has finished (the job owns the state while it runs).
    accum += dt;
    if (accum < 0.05 || running.load()) return;
    accum = std::min(accum - 0.05, 0.1);
    bool work;
    { std::lock_guard<std::mutex> c(cmdMx); work = !commands.empty(); }
    {
        std::lock_guard<std::mutex> lock(simMx); // free: no job is running
        for (auto& s : sites) {
            if (!s->pending.empty()) work = true;
            for (auto& kv : s->bricks) if (kv.second->active || kv.second->meshDirty) { work = true; break; }
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
    (void)world;
    running = true;
    WorldJobs().submit(-0.5f, [] { RunJob(); running = false; });
}

void Fluid::Draw(const Planetoid& body, const RenderContext& ctx, ChunkBatch& batch) {
    if (gpu.empty()) return;
    const double period = 1.0 / std::max(ctx.texScale, 1e-6f);
    auto wrap = [period](double v) { return (float)(fmod(fmod(v, period) + period, period)); };
    const Vector3 clip = { 1e9f, 1e9f, 1e9f };
    batch.setUniform(ctx.waterClipLoc, &clip, SHADER_UNIFORM_VEC3);
    for (auto& kv : gpu) {
        const GpuMesh& g = kv.second;
        if (g.body != body.id || !g.mesh.vaoId) continue;
        const Matrix m = body.localModel(g.anchor, ctx.renderOrigin);
        const Vector3 c = Vector3Transform({ 8, 8, 8 }, m);
        if (Vector3DistanceSqr(c, ctx.cameraPos) > 1500.0f * 1500.0f) continue;
        const Vector3 texOrigin = { wrap(g.anchor.x), wrap(g.anchor.y), wrap(g.anchor.z) };
        batch.setUniform(ctx.waterTexOriginLoc, &texOrigin, SHADER_UNIFORM_VEC3);
        batch.draw(g.mesh, m);
        if (ctx.drawCalls) ++*ctx.drawCalls;
    }
}

namespace {
    // A published cell at a body-local point: its fill (0-255) and flags; false if none.
    bool PubCell(const PubSite& s, Vector3 local, int dk, uint8_t& fill, uint8_t& flags, float* frac = nullptr) {
        const Vector3 d = Vector3Subtract(local, s.origin);
        const float fx = Vector3DotProduct(d, s.e), fy = Vector3DotProduct(d, s.n), fz = Vector3DotProduct(d, s.u);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz) + dk;
        auto it = s.bricks.find(Key(FloorDiv(ci, B), FloorDiv(cj, B), FloorDiv(ck, B)));
        if (it == s.bricks.end()) return false;
        const int i = Idx(ci - FloorDiv(ci, B) * B, cj - FloorDiv(cj, B) * B, ck - FloorDiv(ck, B) * B);
        fill = it->second->fill[i]; flags = it->second->flags[i];
        if (frac) *frac = fz - floorf(fz);
        return true;
    }
}

bool Fluid::Wet(const Planetoid& body, Vector3 local) {
    std::lock_guard<std::mutex> p(pubMx);
    for (const PubSite& s : published) {
        if (s.body != body.id) continue;
        uint8_t f, fl; float frac;
        if (!PubCell(s, local, 0, f, fl, &frac) || (fl & kSolid)) continue;
        if (f / 255.0f > frac) return true;
        // A full cell above counts too (we are in the gap at a surface's edge).
    }
    return false;
}

float Fluid::SurfaceAbove(const Planetoid& body, Vector3 local, float reach) {
    std::lock_guard<std::mutex> p(pubMx);
    for (const PubSite& s : published) {
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
        std::lock_guard<std::mutex> p(pubMx);
        for (const PubSite& s : published) {
            if (s.body != body.id) continue;
            const int R = (int)ceilf(radius);
            for (int dk = -R; dk <= R; ++dk)
                for (int dj = -R; dj <= R; ++dj)
                    for (int di = -R; di <= R; ++di) {
                        const Vector3 q = Vector3Add(local, Vector3Add(Vector3Add(Vector3Scale(s.e, (float)di), Vector3Scale(s.n, (float)dj)), Vector3Scale(s.u, (float)dk)));
                        uint8_t f, fl;
                        if (PubCell(s, q, 0, f, fl)) total += (fl & kSource) ? 100.0f : f / 255.0f;
                    }
        }
    }
    // Static water not yet in a site.
    if (total <= 0.0f) {
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

bool Fluid::StandingWaterNear(const Planetoid& body, Vector3 local, float range) {
    std::lock_guard<std::mutex> p(pubMx);
    for (const PubSite& s : published) {
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
            if (any) wet.push_back(&b);
        }
        fprintf(f, "site %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f %.6f %zu %zu", s->origin.x, s->origin.y, s->origin.z, s->e.x, s->e.y, s->e.z, s->u.x, s->u.y, s->u.z, s->ponds.size(), wet.size());
        for (uint64_t id : s->ponds) fprintf(f, " %llu", (unsigned long long)id);
        fprintf(f, "\n");
        for (Brick* b : wet) {
            // Run-length fill (sources and solid cells as 0: they are re-derived on load).
            fprintf(f, "%d %d %d", b->bi, b->bj, b->bk);
            int i = 0;
            while (i < B3) {
                const int v = (b->flags[i] & (kSource | kSolid)) ? 0 : b->fill[i];
                int run = 1;
                while (i + run < B3 && ((b->flags[i + run] & (kSource | kSolid)) ? 0 : b->fill[i + run]) == v) ++run;
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
    Planetoid* body = Planetoid::find(bodyId);
    const auto* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
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
            if (fscanf(f, " %llu", &id) == 1) { s->ponds.push_back(id); if (pf) pf->claimPond(id); }
        }
        for (size_t j = 0; j < nb; ++j) {
            int bi, bj, bk;
            if (fscanf(f, " %d %d %d", &bi, &bj, &bk) != 3) return;
            Brick* b = s->ensure(bi, bj, bk);
            int i = 0;
            while (i < B3) {
                int v, run;
                if (fscanf(f, " %d*%d", &v, &run) != 2) return;
                for (int r = 0; r < run && i < B3; ++r) if (b) b->fill[i++] = (uint8_t)v; else ++i;
            }
            if (b) b->keepFill = true;
        }
        sites.push_back(std::move(s));
    }
}

Fluid::Stats Fluid::GetStats() { std::lock_guard<std::mutex> p(pubMx); Stats s = pubStats; s.meshes = (int)gpu.size(); return s; }

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
    std::lock_guard<std::mutex> p(pubMx);
    for (const PubSite& s : published) {
        if (s.body != body.id) continue;
        uint8_t f, fl;
        if (!PubCell(s, local, 0, f, fl)) return "no brick";
        return TextFormat("fill %d flags %s%s%s%s", f, (fl & kSolid) ? "S" : "", (fl & kSource) ? "R" : "", (fl & kStatic) ? "T" : "", (fl & kDug) ? "D" : "");
    }
    return "no site";
}
