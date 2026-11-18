#include "fluid.h"
#include "envflag.h"
#include "gpumem.h"
#include "cubemarch.h"
#include "gpucompute.h"
#include "jobsystem.h"
#include "object.h"
#include "plots.h"
#include "rivers.h"
#include "terraindeltas.h"
#include "voxel.h"
#include "waterfield.h"
#include "watersurface.h"
#include "fluidkernel.h"
#include "glcaps.h"
#include "gputimer.h"
#include "external/glad.h"
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
#include <optional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <thread>

namespace {

constexpr int B = 16, B3 = B * B * B;
using WaterField::kSolid; using WaterField::kSource; using WaterField::kStatic; // kStatic: inside infinite generated water
using WaterField::kDug; using WaterField::kSink; // kDug: open only because of an edit: generated infinite water may flow into it
using WaterField::kKindShift;          // flag bits 5-7: the cell's fluid kind (0 water .. 7; Fluid::Kind, fluids.txt)
using WaterField::kKindMask;
using WaterField::CellKind;
// The kind table, read by the step job (written only at startup): sideways divisor and film limit per kind.
std::vector<Fluid::Kind> kinds(1);
int kindVisc[8] = { 5, 5, 5, 5, 5, 5, 5, 5 }, kindMinSpread[8] = { 10, 10, 10, 10, 10, 10, 10, 10 };
constexpr float kSiteReach = 256.0f;   // m: a site reaches this far from its origin (16 bricks of 1 m cells)
// Tests only: RAYLIBGAME_FLUID_SITE_BRICKS=N shrinks a site so a pond straddles its edge (#384).
inline int SiteBricks() {
    static const int n = getenv("RAYLIBGAME_FLUID_SITE_BRICKS") ? std::max(1, atoi(getenv("RAYLIBGAME_FLUID_SITE_BRICKS")))
                                                                 : (int)lroundf(kSiteReach / (16.0f * WaterField::CellSize()));
    return n;
}
constexpr float kSiteJoin = 180.0f;    // a disturbance within this of a site's origin joins it
constexpr int kEvaporate = 2;          // resting specks this thin (~8 mm) dry up (user 2026-11-16: puddles should last a while)
constexpr unsigned kEvapEvery = 1200;   // specks soak away on one step in this many (~10 s on the GPU): only the ones that stay put (#454)
constexpr int kSubsteps = 3;           // CPU simulation steps per 50 ms (water moves a cell a step)
constexpr double kTickSeconds = 0.025; // a tick (a job + meshes) every 25 ms: water moved in 1-3 cell jumps every 50-100 ms, hard to follow (#366)
constexpr float kGpuSubstepSeconds = 0.05f / (2 * kSubsteps); // the GPU runs twice the CPU's substeps: 120 a second (GpuTick)
// Smaller cells take shorter substeps (W4c #642): a cell gives at most what it holds a substep, so a quarter-metre surface
// cell capped a deep channel at ~4 m/s with 120 a second. 1 at 1 m, 4 at 0.25 m.
inline int SubstepScale() {
    // The column model (#680) moves a run's top two cells a pair a substep, at any speed the limiter allows: it needs no
    // extra substeps at small cells (the gravity-wave limit, cs / sqrt(g h), is ~56 ms for 2 m of water at 0.25 m). At
    // 0.25 m one scale ran a quarter of the work and a dam break's front came closer to Ritter: 0.86x to 10 m (4x the
    // substeps: 1.50x), conservation exact (#747). RAYLIBGAME_FLUID_SUBSTEP_SCALE=N forces N (A/B).
    static const int forced = getenv("RAYLIBGAME_FLUID_SUBSTEP_SCALE") ? atoi(getenv("RAYLIBGAME_FLUID_SUBSTEP_SCALE")) : 0;
    static const bool columns = EnvFlag("FLUID_COLUMNS");
    if (forced > 0) return forced;
    return columns ? 1 : std::max(1, (int)lroundf(1.0f / WaterField::CellSize()));
}
// How far down a surface cell's column depth looks (the momentum flow's wetted height): ~8 m, at most 32 cells (W4c #642).
inline int ColumnCells() { return std::clamp((int)lroundf(8.0f / WaterField::CellSize()), 1, 32); }
// Face-flow friction per kind (#642), from its viscosity (water 5): linear (1/s: a pool's last ripples die out) and
// quadratic (drag on the speed: a river's speed levels out with its slope; water ~1.4 m/s on a 1:1000 slope 1 m deep).
inline void KindFriction(int kind, float& lin, float& quad) {
    const float v = (float)kindVisc[kind] / 5.0f;
    lin = 0.05f * v;
    quad = 0.005f * v * v;
}
constexpr int kSleepSteps = 120;       // 2 s without change: a brick sleeps
constexpr int kSettleSteps = 180;      // 3 s window: a brick whose cells all moved less than kSettleDiff in it sleeps (#366)
// Fill units for 1.6 cm (4 at 1 m, 16 at 0.25 m: in units a quarter-metre cell's resting water never settled, W4c #674).
const int kSettleDiff = std::max(4, (int)lroundf(0.016f * 255.0f / WaterField::CellSize()));
// A brick with its levels held is steady flow (kept running) only with a net current of at least this (cm/s, #680); below it,
// at rest. RAYLIBGAME_STEADY_NET=0: any GPU change counts as flow (the old rule), for A/B.
const float kSteadyNet = [] { const char* e = getenv("RAYLIBGAME_STEADY_NET"); return e ? (float)atof(e) : 5.0f; }();
constexpr long long kSteadySampleTicks = 40; // a steady brick is read back once in this many ticks (1 s; staggered by slot, #674)
constexpr long long kFilmHoldSteps = 300; // 5 s of steps: a brick's films stay drawn this long after it stops moving (#660)
inline long long FilmHoldSteps() { static const long long h = getenv("RAYLIBGAME_FILM_HOLD") && atoi(getenv("RAYLIBGAME_FILM_HOLD")) == 0 ? 1 : kFilmHoldSteps; return h; } // A/B: =0 no hold
// Levelling pass interval (steps), by how much water is awake: 4x a second for a few hundred bricks (a dug channel fills
// at W4a's speed), down to once a second for thousands (a big pool's nudges cost ~10 ms a tick, #669).
int awakeBricksLast = 0; // as of the last publish
inline int LevelEvery() { return awakeBricksLast < 300 ? 15 : awakeBricksLast < 1500 ? 30 : 60; }
constexpr int kInitPerStep = 24;       // bricks classified per step, in parallel on the worker pool (each a few ms)
// Bricks per site: 4096 of 1 m cells; smaller cells need more for the same water surface (1/cell^2), at most 32768 (~1.2 GB
// of CPU bricks: W4c #642).
inline int MaxBricks() { static const int n = std::min(32768, (int)lroundf(4096.0f / (WaterField::CellSize() * WaterField::CellSize()))); return n; }

inline int Idx(int x, int y, int z) { return x + y * B + z * B * B; }
inline int64_t Key(int i, int j, int k) { return ((int64_t)(i + 1024) << 22) | ((int64_t)(j + 1024) << 11) | (int64_t)(k + 1024); }
inline int FloorDiv(int a, int b) { return a >= 0 ? a / b : -((-a + b - 1) / b); }
// Partly solid cells (#454): RAYLIBGAME_WATER_SHUT=0 treats every open cell as wholly open (the old model, for A/B).
inline bool ShutOn() { static const bool on = !(getenv("RAYLIBGAME_WATER_SHUT") && getenv("RAYLIBGAME_WATER_SHUT")[0] == '0'); return on; }
constexpr int kMaxShut = 240;          // a cell keeps at least ~6% room

struct Brick {
    int bi = 0, bj = 0, bk = 0;
    std::array<uint8_t, B3> fill{}, flags{}, target{};
    // The share of an open cell that is ground (0-240 of 255, #454), taken as lying at the cell's bottom: the cell holds
    // 255 - shut, and its water stands at shut + fill. Measured at classification; not saved (reclassified on load).
    std::array<uint8_t, B3> shut{};
    Brick* nbr[6] = {};             // -x +x -y +y -z +z
    bool ready = false;             // cells classified (solid / source)
    bool keepFill = false;          // loaded from a save: the classification keeps the stored water
    bool basin = false;             // holds finite generated water (a pond): saved even once drained (#360)
    bool reach = false;             // inside a released river reach (#361): saved even once drained
    bool active = true;
    int idle = 0;
    bool changed = true;            // since the last publish
    std::array<uint8_t, B3> check{}; long long checkStep = -1; // the cells a settle window ago (#366: resting water sleeps)
    bool settled = false;           // put to sleep resting (#366): small changes (the GPU step at its edge) do not wake it
    // CPU mode's dry-cell skip (#747): the wet groups of its last two substeps, kept across ticks (the kernel's masks belong
    // to slots, handed out afresh each tick). Without them a group that dried in a tick's last substep was skipped from the
    // next tick on and kept a stale momentum in its face records.
    uint32_t cpuMask[4] = { ~0u, ~0u, ~0u, ~0u }; // last substep (2 words), the one before; all set until first run
    std::array<uint8_t, B3> drawn{}, drawnFlags{}; // the cells as the chunks last redrew them (visible changes only)
    bool drawnSet = false;
    bool drawnMoving = false;       // published as moving (#649: its thin films are drawn while it is awake)
    std::array<uint8_t, B3> disp{}; // the cells as published: easing toward fill (#366)
    bool dispSet = false, easing = false;
    long long movingStep = -1000000; // the last step it was awake or easing (#660: its films stay drawn a while after)
    // Its counts for the stats (#669: recounted only when it changes - every brick's 4096 cells every tick cost ~10 ms at 3000 bricks).
    bool statSet = false;
    uint8_t statKinds = 0;
    float statWet = 0.0f, statVolume = 0.0f, statOther = 0.0f;
    int statSolid = 0; // solid cells (all solid: nothing to run, #669)
    int slot = -1;                  // GPU mode: its place in the GPU brick pool (-1: none yet)
    std::vector<FluidKernel::Vec4> faces; // CPU mode: the face flows as applied (#654: momentum between ticks; empty: at rest)
    std::vector<int16_t> colFlow;         // the drawn current (#677): per column (x + 16 y) its flow along east, north (cm/s); empty: none yet
    std::vector<int16_t> drawnFlow; int flowStep = -100000; // the current the chunks last drew, and when (#677: a steady stream's levels never change, so it redraws on its flow)
    int cpuSlot = -1;               // CPU mode: its place in this tick's grid
    bool upload = false;            // GPU mode: the whole brick (cells, flags, targets) goes up next tick
    int sinkRelease = -1;           // its cells hold a released reach's downstream sink: that release's index (river flow, #656)
    int blockRelease = -1, blocked = 0; // cells of ground added into that released reach's river (the dam: #656 follow-up)
    long long slotNeeded = 0;       // GPU mode: the last tick it was awake or near awake water (#675: else its slot is given back)
    bool hadSlot = false;           // GPU mode: it held a slot before (back from release: uploaded without waking)
    // Steady flow (#674): levels stable over a settle window while water keeps moving through (the GPU reported changes in at
    // least half its ticks). It keeps running - the water still reaches what is downstream - but is read back and published
    // once a second, not every tick, and never falls asleep (a steady stream's bricks slept and woke every second, #660).
    bool steady = false;
    int winTicks = 0, winFlow = 0;  // this settle window: ticks published, ticks the GPU changed it
    float checkVol = 0.0f;          // its water (m^3) at the window's start: a slow drain or fill changes it, jitter does not (#680)
};
inline int Cap(const Brick* b, int i) { return 255 - b->shut[i]; }          // what the cell holds (#454)
inline int Height(const Brick* b, int i) { return b->shut[i] + b->fill[i]; } // where its water stands, in 255ths of the cell

static Vector3 dbgAxO{}, dbgAxX{}, dbgAxY{}, dbgAxU{}; static bool dbgAxSet = false; // TEMP #680 (DBTRACE): the test's channel frame

struct Site {
    uint64_t body = 0;
    int index = 0;
    Vector3 origin{}, e{}, n{}, u{};
    std::unordered_map<int64_t, std::unique_ptr<Brick>> bricks;
    std::deque<Brick*> pending;     // waiting for classification
    // Sources with a rate and a head (springs, seeps): site cell, kind, fill units added per step at most,
    // and the level (site cells, fractional) the fluid surface over it is held to.
    struct Spring { int ci, cj, ck, kind; float perStep, head; float acc = 0.0f; int rel = -1; float scale = 1.0f; }; // acc: fill units owed (a slow source adds < 1 a step); rel: the released reach it feeds, scale: the river's flow arriving there (#656)
    std::vector<Spring> springs;
    // Released river reaches (#361): ground added into a river (a dam) turns its static source cells within `r` cells
    // (horizontally) of the edit into ordinary water. The upstream edge is fed by springs at the reach's discharge
    // (in `springs`), the downstream edge is a sink held at the generated level; beyond `r` the river stays generated.
    struct Release { int ci, cj; float r; Vector3 centre, down; float natQ = 0.0f; int piece = -1; double sinkUnits = 0.0; float share = 1.0f; int upPiece = -1; int peakBlocked = 0; }; // natQ: the springs' m^3/s; piece: the river link below it; sinkUnits: taken in since the last measure; share: outflow / natQ (#656)
    std::vector<Release> releases;

    Vector3 cellCentre(int ci, int cj, int ck) const { return cellPoint(ci + 0.5f, cj + 0.5f, ck + 0.5f); }
    // A point at site coordinates (cells, continuous), body-local (W4c #642: the axes are unit, cells CellSize() m).
    Vector3 cellPoint(float fi, float fj, float fk) const {
        const float cs = WaterField::CellSize();
        return Vector3Add(origin, Vector3Add(Vector3Add(Vector3Scale(e, fi * cs), Vector3Scale(n, fj * cs)), Vector3Scale(u, fk * cs)));
    }
    // A body-local point's height over the origin, in cells (continuous).
    float heightCells(Vector3 local) const { return Vector3DotProduct(Vector3Subtract(local, origin), u) / WaterField::CellSize(); }
    void cellOf(Vector3 local, int& ci, int& cj, int& ck) const {
        const Vector3 d = Vector3Scale(Vector3Subtract(local, origin), 1.0f / WaterField::CellSize());
        ci = (int)floorf(Vector3DotProduct(d, e)); cj = (int)floorf(Vector3DotProduct(d, n)); ck = (int)floorf(Vector3DotProduct(d, u));
    }
    bool inRange(int bi, int bj, int bk) const { const int r = SiteBricks(); return abs(bi) <= r && abs(bj) <= r && abs(bk) <= r; }
    Brick* find(int bi, int bj, int bk) const { auto it = bricks.find(Key(bi, bj, bk)); return it == bricks.end() ? nullptr : it->second.get(); }
    // Allocates (unclassified, queued) if missing and in range; nullptr if out of range or full.
    Brick* ensure(int bi, int bj, int bk) {
        if (Brick* b = find(bi, bj, bk)) return b;
        if (!inRange(bi, bj, bk) || (int)bricks.size() >= MaxBricks()) return nullptr;
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

// Build parts as solid for water (#391): an immutable copy of each plot's solid cells, made on the main thread when a
// plot's version changes and read by the classification jobs. A fluid cell is solid when any of seven points in it
// (its centre and +-0.35 m along the plot's axes) lies in a solid part cell: plot cells are 0.5 m and turned against
// the site's grid, so one block's wall at any angle still holds water.
struct PlotSnap {
    uint32_t id = 0, version = 0;
    Vector3 origin{}, east{}, up{}, north{}, centre{};
    float radius = 0.0f;                 // a sphere around its solid cells (body-local), for a quick reject
    std::unordered_set<Int3> solid;      // plot cells whose shape holds the cell's centre
    bool cellSolid(Vector3 q) const {
        const Vector3 d = Vector3Subtract(q, origin);
        const Int3 c{ (int)floorf(Vector3DotProduct(d, east) / kCellSize), (int)floorf(Vector3DotProduct(d, up) / kCellSize), (int)floorf(Vector3DotProduct(d, north) / kCellSize) };
        return solid.count(c) > 0;
    }
};
struct PartSolid {
    std::vector<PlotSnap> plots;
    bool solidAt(Vector3 c) const {
        for (const PlotSnap& p : plots) {
            if (p.solid.empty() || Vector3DistanceSqr(c, p.centre) > (p.radius + 1.0f) * (p.radius + 1.0f)) continue;
            if (p.cellSolid(c)) return true;
            // Probes inside the fluid cell, 0.35 of its edge from the centre along the plot's axes (#694: a fixed 0.35 m reached past
            // 0.25 m cells - water stood off every block face and a 1 m hole in a block dam passed ~9% of its flow).
            // Cells smaller than a plot cell (0.5 m): solid when a part covers the centre or most of the cell (4 of 6 probes) -
            // exact enough at that size, and a one-block hole opens. Larger cells: any probe, so a half-block panel or pane stays
            // watertight; the cost is that a one-block hole half a block off the fluid grid stays shut at 1 m cells.
            const float cs = WaterField::CellSize(), off = 0.35f * cs;
            const Vector3 ax[3] = { p.east, p.up, p.north };
            int hits = 0;
            for (const Vector3& a : ax)
                hits += (int)p.cellSolid(Vector3Add(c, Vector3Scale(a, off))) + (int)p.cellSolid(Vector3Subtract(c, Vector3Scale(a, off)));
            if (hits >= (cs < 0.5f ? 4 : 1)) return true;
        }
        return false;
    }
};

struct BodyRef {
    std::shared_ptr<const DensityField> field;
    std::shared_ptr<const TerrainDeltas> edits;
    const PlanetField* pf = nullptr;
    std::shared_ptr<const PartSolid> parts; // build parts (#391), or null
};
std::unordered_map<uint64_t, std::shared_ptr<const PartSolid>> partSnaps; // per body (main thread): its plots' solid cells (#391)

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
    b.basin = false; b.reach = false; b.sinkRelease = -1; b.blockRelease = -1; b.blocked = 0;
    const float cs = WaterField::CellSize(), hc = 0.5f * cs; // a cell, half a cell (m: W4c #642)
    std::vector<float> corner((B + 1) * (B + 1) * (B + 1), NAN); // the density at the cell corners, sampled when first needed
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
                // Inside a released reach the natural river (#656): the reach is the simulated river there, and a cut
                // further up must not classify it dry (no sink: its outflow could never be measured again).
                bool inReach = false;
                for (const Site::Release& rel : s.releases)
                    if ((float)((ci - rel.ci) * (ci - rel.ci) + (cj - rel.cj) * (cj - rel.cj)) < rel.r * rel.r) { inReach = true; break; }
                std::optional<RiverNetwork::NaturalScope> natural;
                if (inReach) natural.emplace();
                const PlanetField::WaterSource ws = pf->waterSourceAt(dir);
                if (ws.kind == 3) {
                    if (bed < sea + ws.level && sea + ws.level > level) { level = sea + ws.level; kind = 5; } // a pond (never an infinite source, even its rim)
                } else if (kind != 4 && ws.kind && sea + ws.level > level) { // at sea the sea's level wins (#190: a river mouth's higher level fed the sea columns: overfilled channels that never slept)
                    level = sea + ws.level; kind = ws.kind; }
            }
            const float top = pf ? bed + 14.0f : 1e9f;  // tors, boulders and detail stay below this
            // A released river reach (#361): 1 inside (ordinary water), 2 its downstream edge (a sink at the level).
            int released = 0, relIdx = -1;
            if (kind == 1)
                for (size_t ri = 0; ri < s.releases.size(); ++ri) {
                    const Site::Release& rel = s.releases[ri];
                    const float d = sqrtf((float)((ci - rel.ci) * (ci - rel.ci) + (cj - rel.cj) * (cj - rel.cj)));
                    if (d >= rel.r) continue;
                    released = (d >= rel.r - 2.0f && Vector3DotProduct(Vector3Subtract(mid, rel.centre), rel.down) > 0.0f) ? 2 : 1;
                    if (released == 2) b.sinkRelease = (int)ri; // its outflow is measured here (river flow, #656)
                    relIdx = (int)ri;
                    b.reach = true;
                    break;
                }
            // The column's ground first: generated water reaches a cell only through open cells from the water's surface
            // down (#559: a cave dug under a river, roofed by rock, filled with the river's water when a torch placed in
            // it reclassified its brick). Dug through, the brick is reclassified and the river pours in.
            float dcCol[B];
            bool solidCol[B];
            for (int z = 0; z < B; ++z) {
                const Vector3 c = s.cellCentre(ci, cj, b.bk * B + z);
                const float r = Vector3Length(c);
                dcCol[z] = r < top ? ref.field->density(c, ref.edits.get()) : 0.0f;
                solidCol[z] = (r < top && dcCol[z] >= ISO_LEVEL) || (ref.parts && ref.parts->solidAt(c)); // build parts hold water too (#391)
            }
            int coveredTo = -1; // cells at or below this z have ground between them and the water's surface
            if (kind) {
                bool roof = false; // ground above the brick, under the level
                for (int k = B; !roof && k < 2 * B; ++k) {
                    const Vector3 q = s.cellCentre(ci, cj, b.bk * B + k);
                    const float rq = Vector3Length(q);
                    if (rq - hc >= level || rq >= top) break;
                    roof = ref.field->density(q, ref.edits.get()) >= ISO_LEVEL || (ref.parts && ref.parts->solidAt(q));
                }
                if (roof) coveredTo = B - 1;
                else
                    for (int z = B - 1; z >= 0; --z)
                        if (solidCol[z] && Vector3Length(s.cellCentre(ci, cj, b.bk * B + z)) - hc < level) { coveredTo = z - 1; break; }
            }
            for (int z = 0; z < B; ++z) {
                const int i = Idx(x, y, z);
                const Vector3 c = s.cellCentre(ci, cj, b.bk * B + z);
                const float r = Vector3Length(c);
                const float dc = dcCol[z];
                const bool solid = solidCol[z];
                uint8_t fl = solid ? kSolid : 0;
                // Ground added into the reach's river (#656 follow-up): a solid cell under the natural level that the
                // natural ground leaves open - the dam. None left: the river runs as generated again (UpdateRiverFlow).
                if (solid && relIdx >= 0 && r - hc < level && r > bed - hc && ref.field->density(c, nullptr) < ISO_LEVEL) {
                    b.blockRelease = relIdx; ++b.blocked;
                }
                // How much of an open cell is ground (#454), in cells near the surface: 16 points (2 x 2 across, 4 heights)
                // interpolated from the density at the cell's 8 corners - as marching cubes draws the ground. The corners
                // are sampled once and shared with the neighbouring cells (16 field samples a cell cost ~40 ms of a step).
                int sh = 0;
                if (!solid && ShutOn() && r < top && dc > ISO_LEVEL - 0.75f * DENSITY_GRADIENT / kVoxelSize) {
                    float d8[8];
                    for (int k = 0; k < 8; ++k) {
                        const int ax = x + (k & 1), ay = y + ((k >> 1) & 1), az = z + (k >> 2);
                        float& d = corner[(az * (B + 1) + ay) * (B + 1) + ax];
                        if (std::isnan(d)) {
                            const Vector3 q = s.cellPoint((float)(ci - x + ax), (float)(cj - y + ay), (float)(b.bk * B + az));
                            d = ref.field->density(q, ref.edits.get());
                        }
                        d8[k] = d;
                    }
                    int n = 0;
                    for (int a = 0; a < 2; ++a)
                        for (int bb = 0; bb < 2; ++bb)
                            for (int h = 0; h < 4; ++h) {
                                const float tx = a ? 0.75f : 0.25f, ty = bb ? 0.75f : 0.25f, tz = 0.125f + 0.25f * h;
                                const float d00 = d8[0] + (d8[1] - d8[0]) * tx, d10 = d8[2] + (d8[3] - d8[2]) * tx;
                                const float d01 = d8[4] + (d8[5] - d8[4]) * tx, d11 = d8[6] + (d8[7] - d8[6]) * tx;
                                const float lo = d00 + (d10 - d00) * ty, hi = d01 + (d11 - d01) * ty;
                                n += lo + (hi - lo) * tz >= ISO_LEVEL ? 1 : 0;
                            }
                    sh = std::min(kMaxShut, n * 255 / 16);
                }
                b.shut[i] = (uint8_t)sh;
                // Dug: open only because of an edit. Infinite water feeds only these (the real ground has
                // crags the height map doesn't: natural dips by a river are not to be flooded).
                if (!solid && r < bed + 4.0f && r < top && ref.field->density(c, nullptr) >= ISO_LEVEL) fl |= kDug;
                // An open cell in the water's band under a roof (#559) gets no generated water, but the water may reach it
                // from the side - a tunnel into a bank, a hollow in built rock: fed like a dug cell, and a sealed cave stays dry.
                if (!solid && kind && kind != 5 && z <= coveredTo && r - hc < level && r > bed - hc) fl |= kDug;
                uint8_t tg = 0;
                if (!solid && kind && z > coveredTo && r - hc < level && r > bed - hc) {
                    tg = (uint8_t)std::clamp((int)lroundf((level - (r - hc)) / cs * 255.0f), 0, 255);
                    if ((int)tg - sh > 0 && kind != 5) fl |= released == 1 ? kDug : released == 2 ? (kSource | kStatic | kSink) : (kSource | kStatic);
                    if (tg > 0 && kind == 5) {
                        b.basin = true;
                        faces |= (x == 0 ? 1 : 0) | (x == B - 1 ? 2 : 0) | (y == 0 ? 4 : 0) | (y == B - 1 ? 8 : 0) | (z == 0 ? 16 : 0) | (z == B - 1 && tg == 255 ? 32 : 0);
                    }
                }
                if (b.keepFill && !solid && !(fl & kSource) && b.fill[i] > 0) fl |= (b.flags[i] & kKindMask); // kept fluid stays its kind
                b.flags[i] = fl;
                const int tf = std::max(0, (int)tg - sh); // the generated level is a height: the ground under it holds no water (#454)
                b.target[i] = (fl & kSource) ? (uint8_t)tf : 0;
                if (b.keepFill) { if (solid) b.fill[i] = 0; else if (fl & kSource) b.fill[i] = (uint8_t)tf; else b.fill[i] = (uint8_t)std::min((int)b.fill[i], 255 - sh); }
                else b.fill[i] = solid ? 0 : (uint8_t)tf;
            }
        }
    // Kept (loaded or reclassified) bricks pull their pond neighbours in too (#382 review): ensure() is a no-op for those
    // that exist; one never classified (saved mid pull-in, an old save's claimed pond) comes in fresh with its generated
    // water, so a drained brick never borders the generated pond (a wall of water at the brick face).
    b.ready = true;
    b.upload = true; // GPU mode: the classified brick goes up whole
    b.faces.clear(); // CPU mode: a (re)classified brick starts at rest (#654)
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

long long wakeGpu = 0, wakeTouch = 0; // settled bricks woken (stats)
std::atomic<long long> pourDroppedUnits{ 0 }; // fill units pours found no room for, even widened (#766)
long long filmFlips = 0;              // bricks whose films switched shown/hidden (stats, #660)
long long slotReleases = 0;           // GPU slots given back by sleeping bricks (stats, #675)
int slotsUsed = 0;                    // GPU slots held now (stats, #675)
int G_Capacity();                     // the GPU brick pool (below)
constexpr long long kSlotKeepTicks = 80; // 2 s of ticks: a brick away from awake water keeps its slot this long (#675)
void Touch(Brick* b) { if (b->settled) ++wakeTouch; b->active = true; b->idle = 0; b->changed = true; b->settled = false; b->steady = false; }
void CpuTick(int n); // CPU mode: the GPU's step on the CPU (#654), below
// Whether a neighbour (on face fc of a brick whose water moved) needs waking and running (#669): not a brick of solid rock
// (water never enters it) nor a dry brick above (water rises into it only through sources, pours or levelling, which wake
// it). Under and over a big pool they were most of every dispatch.
inline bool RunsBeside(const Brick& o, int fc) { return !(o.statSet && (o.statSolid == B3 || (fc == 5 && o.statKinds == 0))); }

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
inline void TakeKind(Brick* to, int it, uint8_t kindBits) { if (to->fill[it] == 0) to->flags[it] = (uint8_t)((to->flags[it] & ~kKindMask) | kindBits); }

// A source with a rate and a head: it adds to the surface cell of its column while that surface is
// below the head. A basin over it fills to the head and the source rests; a lower rim lets the overflow
// run off (the step carries it away), so the source keeps adding and a stream forms.
// The pressurised outlet (W4c #642, user: a small source pushes out more than its own cell holds): the spring's
// accumulator is the pressure. It wells up its column (through up to ~3 m of ground over a buried source) to the first
// open cell; once that is full, the head pushes the rest on through the full cells of its fluid into the open cells beside
// and above them, out to SpringReach() cells sideways - 0 for 1 m cells (the column alone, as before), 1 for quarter-metre
// cells (a 3 x 3 outlet: one cell takes 1/64 m^3 a step, so a source of a few m^3/s could not get its water out).
// open(b, i, ck) takes from the walk an open (not full) cell under the head and returns true to stop; missing() is told of
// a brick still to come (and stops the walk).
inline int SpringReach() { return (int)(0.25f / WaterField::CellSize() + 1e-3f); }
template <class Open, class Missing>
void WalkOutlet(const Site& s, const Site::Spring& sp, Open&& open, Missing&& missing) {
    const int R = SpringReach();
    // The start: the column's first cell that is not ground (a buried source wells up), then up through its full cells.
    int k0 = -1;
    const int buried = std::max(3, (int)lroundf(3.0f / WaterField::CellSize())); // ~3 m of ground over a buried source
    for (int k = 0; k < buried && k0 < 0; ++k) {
        Brick* b; int i;
        if (!s.cell(sp.ci, sp.cj, sp.ck + k, b, i)) { missing(sp.ci, sp.cj, sp.ck + k); return; }
        if (!(b->flags[i] & kSolid)) k0 = sp.ck + k;
    }
    if (k0 < 0) return;
    struct C { int dx, dy, ck; };
    std::vector<C> q{ { 0, 0, k0 } };
    std::vector<int64_t> seen{ Key(0, 0, 0) }; // relative to the start (cell keys hold +-1024)
    for (size_t h = 0; h < q.size() && q.size() < 4096; ++h) {
        const C c = q[h];
        if ((float)c.ck >= sp.head || c.ck - k0 >= 64) continue;  // the surface is at the head
        Brick* b; int i;
        if (!s.cell(sp.ci + c.dx, sp.cj + c.dy, c.ck, b, i)) { missing(sp.ci + c.dx, sp.cj + c.dy, c.ck); return; }
        if (b->flags[i] & kSolid) continue;
        if (b->fill[i] > 0 && CellKind(b->flags[i]) != sp.kind) continue; // another fluid lies there
        if (b->fill[i] < Cap(b, i)) {
            if (open(b, i, c.ck)) return;
            if (R == 0 || b->fill[i] < Cap(b, i)) { if (R == 0) return; continue; } // still open: its surface holds the pressure
        }
        // Full: the pressure passes on - up, and sideways within the outlet.
        static const int d[5][3] = { { 0, 0, 1 }, { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 } };
        for (const auto& dd : d) {
            const C n{ c.dx + dd[0], c.dy + dd[1], c.ck + dd[2] };
            if (abs(n.dx) > R || abs(n.dy) > R) continue;
            const int64_t key = Key(n.dx, n.dy, n.ck - k0);
            if (std::find(seen.begin(), seen.end(), key) != seen.end()) continue;
            seen.push_back(key); q.push_back(n);
        }
    }
}

void StepSpring(Site& s, Site::Spring& sp) {
    sp.acc = std::min(sp.acc + sp.perStep * sp.scale, std::max(16.0f, sp.perStep * 30.0f)); // owed at most ~half a second's worth
    if (sp.acc < 8.0f) return; // in drops of 3 cm+: thinner specks would soak away (kEvaporate) before they spread
    WalkOutlet(s, sp,
        [&](Brick* b, int i, int ck) {
            const float surface = ck + Height(b, i) / 255.0f;
            const int add = (int)std::min({ sp.acc, (sp.head - surface) * 255.0f, (float)(Cap(b, i) - b->fill[i]) });
            if (add <= 0) return SpringReach() == 0; // at the head here (another outlet cell may be lower)
            TakeKind(b, i, (uint8_t)(sp.kind << kKindShift));
            b->fill[i] = (uint8_t)(b->fill[i] + add);
            GpuDelta(b, i, add);
            sp.acc -= (float)add;
            Touch(b); ++dbgSrc;
            return sp.acc < 1.0f;
        },
        [&](int ci, int cj, int ck) { s.ensure(FloorDiv(ci, B), FloorDiv(cj, B), FloorDiv(ck, B)); });
}

// Whether a spring would still add (an outlet cell's surface is under its head): read-only, the same walk as StepSpring.
// Fluid::Update steps while one does, though every brick sleeps (#645: a dammed river's inflow stopped once the pool
// behind the dam settled, so it never rose to the crest).
bool SpringWants(const Site& s, const Site::Spring& sp) {
    if (sp.scale <= 0.0f) return false; // its river is cut upstream (#656)
    bool wants = false;
    WalkOutlet(s, sp,
        [&](Brick* b, int i, int ck) { wants = sp.head - (ck + Height(b, i) / 255.0f) > 1.0f / 255.0f; return wants || SpringReach() == 0; },
        [&](int, int, int) { wants = true; }); // its brick is still to come
    return wants;
}

// Levelling: bodies of water joined through full cells communicate. Each body's open surfaces (the
// partial cells and the empty cells over full ones) ease toward one level: its sources' level if it
// has any, else the level its volume fills to.
void LevelSite(Site& s) {
    // In metres (W4c: in cells, a quarter-metre channel's sloping surface - more than 7.5 cm end to end - was levelled, moving its
    // water at once instead of the face flows carrying it: a flow meter read a third of the stream). Cells: / CellSize().
    const float kLevelSpread = 0.3f / WaterField::CellSize(), kLevelNear = 0.25f / WaterField::CellSize();
    std::unordered_set<int64_t> seen;
    struct C { Brick* b; int i; int ck; bool leaf; };
    for (auto& kv : s.bricks) {
        Brick* b0 = kv.second.get();
        if (!b0->ready || !b0->active) continue;
        for (int i0 = 0; i0 < B3; ++i0) {
            if (b0->fill[i0] < Cap(b0, i0) || (b0->flags[i0] & kSolid)) continue; // full: at its capacity (#454)
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
                if (b->flags[i] & kSource) { const int th = b->shut[i] + b->target[i]; hasSource = true; srcLevel = std::max(srcLevel, c[2] + th / 255.0f); if (th < 255) srcSurface = std::min(srcSurface, c[2] + th / 255.0f); }
                const int d[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
                for (int k = 0; k < 6; ++k) {
                    const int nc[3] = { c[0] + d[k][0], c[1] + d[k][1], c[2] + d[k][2] };
                    Brick* nb; int ni;
                    if (!s.cell(nc[0], nc[1], nc[2], nb, ni) || (nb->flags[ni] & kSolid)) continue;
                    if (nb->fill[ni] > 0 && (nb->flags[ni] & kKindMask) != kind0) continue; // the other fluid
                    const int64_t nk = Key(nc[0], nc[1], nc[2]);
                    if (nb->fill[ni] >= Cap(nb, ni)) {
                        if (!seen.count(nk)) { seen.insert(nk); q.push_back({ nc[0], nc[1], nc[2] }); }
                    } else if ((nb->fill[ni] > 0 || k == 4) && !leafSeen.count(nk)) {
                        leafSeen.insert(nk);
                        comp.push_back({ nb, ni, nc[2], true });
                        if (nb->flags[ni] & kSource) { const int th = nb->shut[ni] + nb->target[ni]; hasSource = true; srcLevel = std::max(srcLevel, nc[2] + th / 255.0f); if (th < 255) srcSurface = std::min(srcSurface, nc[2] + th / 255.0f); }
                    }
                }
            }
            if (comp.size() >= 150000) continue;
            // Only bodies whose open surfaces are apart - a U-bend's arms, a flooded chamber under its shaft - are levelled. One
            // surface joined side to side (a step up or down for a slope) is the face flows' to carry: levelling it moved a
            // channel's water at once (W4c: a quarter-metre trench's meter read a third of its stream, the rest teleported).
            {
                std::unordered_map<int64_t, int> leafAt;
                std::vector<int> parent;
                std::vector<std::array<int, 3>> leafPos;
                for (const C& c : comp)
                    if (c.leaf) {
                        const int li = c.i, lx = c.b->bi * B + li % B, ly = c.b->bj * B + (li / B) % B;
                        leafAt.emplace(Key(lx, ly, c.ck), (int)parent.size());
                        leafPos.push_back({ lx, ly, c.ck });
                        parent.push_back((int)parent.size());
                    }
                auto root = [&](int a) { while (parent[a] != a) a = parent[a] = parent[parent[a]]; return a; };
                int groups = (int)parent.size();
                for (size_t a = 0; a < leafPos.size(); ++a)
                    for (int dz = -1; dz <= 1; ++dz)
                        for (int dy = -1; dy <= 1; ++dy)
                            for (int dx = -1; dx <= 1; ++dx) {
                                if ((dx != 0) + (dy != 0) + (dz != 0) == 0 || (dx != 0 && dy != 0)) continue; // side, up/down, side +- a layer
                                auto it = leafAt.find(Key(leafPos[a][0] + dx, leafPos[a][1] + dy, leafPos[a][2] + dz));
                                if (it == leafAt.end()) continue;
                                const int ra = root((int)a), rb = root(it->second);
                                if (ra != rb) { parent[ra] = rb; --groups; }
                            }
                if (groups <= 1) continue;
                // A surface counts from 1 m^2 (one cell at 1 m, 16 at 0.25 m): a stray wet cell under a lip or a drip is not an
                // arm of a U-bend, and with it the whole trench was levelled.
                const int minLeaves = std::max(1, (int)lroundf(1.0f / (WaterField::CellSize() * WaterField::CellSize())));
                std::unordered_map<int, int> size;
                for (size_t a = 0; a < leafPos.size(); ++a) ++size[root((int)a)];
                int real = 0;
                for (const auto& kv2 : size) real += kv2.second >= minLeaves ? 1 : 0;
                if (real <= 1) continue;
            }
            float lo = 1e9f, hi = -1e9f, vol = 0.0f;
            int kmin = INT32_MAX, kmax = INT32_MIN;
            for (const C& c : comp) {
                vol += c.b->fill[c.i] / 255.0f;
                kmin = std::min(kmin, c.ck); kmax = std::max(kmax, c.ck);
                if (c.leaf) { const float h = c.ck + Height(c.b, c.i) / 255.0f; lo = std::min(lo, h); hi = std::max(hi, h); }
            }
            if (hi < lo) continue; // no open surface (a sealed cave full to the roof)
            float level;
            if (hasSource) {
                level = srcSurface < 1e8f ? srcSurface : srcLevel; // a sloping river: its lowest surface, never overfilling
                if (lo > level - kLevelNear && hi < level + kLevelNear) continue;
            } else {
                if (hi - lo < kLevelSpread) continue;
                float a = (float)kmin, z = (float)kmax + 1.0f;
                for (int it = 0; it < 24; ++it) {
                    const float m = 0.5f * (a + z);
                    float cap = 0.0f;
                    for (const C& c : comp) cap += std::clamp((m - c.ck) * 255.0f - c.b->shut[c.i], 0.0f, (float)Cap(c.b, c.i)) / 255.0f; // the water under level m
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
                const float want = std::clamp((level - c.ck) * 255.0f - c.b->shut[c.i], 0.0f, (float)Cap(c.b, c.i)), now = c.b->fill[c.i];
                if (hasSource && want <= now) continue; // fed by static water: only raise (a sloping river feeds at several levels)
                worst = std::max(worst, fabsf(want - now));
            }
            if (worst < 6.0f) continue;
            float carry = 0.0f;
            for (const C& c : comp) {
                if (c.b->flags[c.i] & kSource) continue;
                const float want = std::clamp((level - c.ck) * 255.0f - c.b->shut[c.i], 0.0f, (float)Cap(c.b, c.i)), now = c.b->fill[c.i];
                if (hasSource && want <= now) continue;
                const float d = (want - now) * 0.3f + (hasSource ? 0.0f : carry);
                const int di = (int)lroundf(d);
                const int nf = std::clamp((int)now + di, 0, Cap(c.b, c.i));
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
    const int R = (int)ceilf(r / WaterField::CellSize()) + 2;
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
    const int R = (int)ceilf(r / WaterField::CellSize()) + 2;
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

// Ground added into a river's water (a dam, #361): release the reach around it. Its bricks are (re)classified (the river
// cells inside become ordinary water, the downstream edge a sink), and the upstream edge gets springs that add the
// inflow the simulated channel can carry: ~0.25 m^3/s per metre of width (#432: the cellular step has no momentum, so a
// gentle surface barely flows; a river's real discharge, 17 m^3/s on a 10 m river, mounded up and flooded the banks), with a
// head 3 m over the level so a backed-up pool keeps rising until it spills over a low dam. Beyond the reach the river
// stays generated: a pool that rises past the upstream edge shows a step there.
constexpr float kReleaseMetres = 40.0f;
void MaybeRelease(Site& s, const BodyRef& ref, Vector3 c, float r) {
    const PlanetField* pf = ref.pf;
    if (!pf || !pf->rivers()) return;
    const float cl = Vector3Length(c);
    if (cl < 1.0f) return;
    const Vector3 dir = Vector3Scale(c, 1.0f / cl);
    const PlanetField::WaterSource ws = pf->waterSourceAt(dir);
    if (ws.kind != 1) return;                                   // rivers only (lakes and the sea stay as they are)
    const float sea = pf->seaLevelRadius(), level = sea + ws.level;
    if (cl - r > level) return;                                 // the edit is above the water
    int ci, cj, ck;
    s.cellOf(c, ci, cj, ck);
    for (const Site::Release& o : s.releases)
        if ((ci - o.ci) * (ci - o.ci) + (cj - o.cj) * (cj - o.cj) < (int)(o.r * o.r * 0.36f)) return; // already released
    const RiverNetwork::Hit h = pf->rivers()->nearest(dir, 1e4f);
    if (!h.found) return;
    Site::Release rel{ ci, cj, kReleaseMetres / WaterField::CellSize(), c, h.down }; // r in cells
    // The river link just below the reach (#656): its flow follows what the reach lets out.
    rel.piece = pf->rivers()->pieceBeyond(h.piece, c, kReleaseMetres + 5.0f);
    rel.upPiece = h.piece; // the link through the reach: the flow arriving from upstream (a cut further up feeds less, #656)
    s.releases.push_back(rel);
    // The reach's bricks: per brick column, the river's depth under its own level down, 8 m over it up (#645: a reach
    // that climbs upstream - 16 m over its 40 m on seed 42's test river - had no bricks at its upstream edge).
    auto levelZ = [&](Vector3 d, float lv) { return s.heightCells(Vector3Scale(d, sea + lv)); }; // site height (cells)
    const int R = (int)rel.r;
    for (int bj = FloorDiv(cj - R, B); bj <= FloorDiv(cj + R, B); ++bj)
        for (int bi = FloorDiv(ci - R, B); bi <= FloorDiv(ci + R, B); ++bi) {
            float lo = 1e30f, hi = -1e30f; // the river levels over the brick column's corners and middle (else the edit's)
            for (int k = 0; k < 5; ++k) {
                const int x = bi * B + (k == 4 ? B / 2 : (k & 1) * (B - 1)), y = bj * B + (k == 4 ? B / 2 : (k >> 1) * (B - 1));
                const Vector3 md = Vector3Normalize(s.cellCentre(x, y, 0));
                const PlanetField::WaterSource wq = pf->waterSourceAt(md);
                const float z = levelZ(md, wq.kind == 1 ? wq.level : ws.level);
                lo = std::min(lo, z); hi = std::max(hi, z);
            }
            for (int bk = FloorDiv((int)floorf(lo - (h.depth + 3.0f) / WaterField::CellSize()), B); bk <= FloorDiv((int)ceilf(hi + 8.0f / WaterField::CellSize()), B); ++bk)
                if (Brick* b = s.ensure(bi, bj, bk); b && b->ready) {
                    b->ready = false; b->keepFill = true;
                    if (std::find(s.pending.begin(), s.pending.end(), b) == s.pending.end()) s.pending.push_front(b);
                }
        }
    int lci, lcj, lck;
    s.cellOf(Vector3Scale(dir, level), lci, lcj, lck);
    // Upstream edge: one ring of river columns behind the edit along the flow, each with a head 3 m over its own level
    // (#645: the edit's level, 16 m under the upstream bed on a climbing reach, held every spring at its head from the start).
    struct Up { int ci, cj, ck; float head; };
    std::vector<Up> up;
    for (int dy = -R; dy <= R; ++dy)
        for (int dx = -R; dx <= R; ++dx) {
            const float d = sqrtf((float)(dx * dx + dy * dy));
            if (d < rel.r - 2.0f || d >= rel.r - 1.0f) continue;
            const Vector3 m = s.cellCentre(ci + dx, cj + dy, lck);
            if (Vector3DotProduct(Vector3Subtract(m, c), rel.down) >= 0.0f) continue;
            const Vector3 md = Vector3Normalize(m);
            const PlanetField::WaterSource wq = pf->waterSourceAt(md);
            if (wq.kind != 1) continue;
            int bci, bcj, bck;
            s.cellOf(Vector3Scale(md, pf->radius + pf->elevation(md)), bci, bcj, bck);
            up.push_back({ ci + dx, cj + dy, bck, levelZ(md, wq.level) + 3.0f / WaterField::CellSize() });
        }
    if (up.empty()) return;
    const float q = std::max(1.0f, std::min({ 0.25f * 2.0f * h.halfWidth, 40.0f, 2.0f * h.halfWidth * h.depth })); // m^3/s (#432)
    s.releases.back().natQ = q; // what flows in when nothing blocks the river (#656: its outflow is measured against it)
    for (const Up& u : up)
        { s.springs.push_back({ u.ci, u.cj, u.ck, 0, q / (float)up.size() * 255.0f / WaterField::CellVolume() / (20.0f * kSubsteps), u.head }); s.springs.back().rel = (int)s.releases.size() - 1; }
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
        if (near) {
            EnsureAround(*s, c.c, c.r);
            if (c.m3 > 0.0f) MaybeRelease(*s, ref, c.c, c.r); // ground added into a river: release the reach (#361)
        } else {
            // Only bricks that exist: their solidity changed.
            int ci, cj, ck; s->cellOf(c.c, ci, cj, ck);
            const int R = (int)ceilf(c.r / WaterField::CellSize()) + 2;
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
        const float fz = s->heightCells(c.c);
        Site::Spring sp{ ci, cj, ck, std::clamp(c.fluid, 0, 7), c.r * 255.0f / WaterField::CellVolume() / (20.0f * kSubsteps), fz + c.head / WaterField::CellSize() };
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
        // Pour: into the open cells at or above the point, layer by layer up to 6 m, over a ~1 m footprint (W4c #642: one
        // column of quarter-metre cells held a tenth of a bucket).
        float left = c.m3 * 255.0f / WaterField::CellVolume();
        Brick* b; int i;
        if (!s->cell(ci, cj, ck, b, i)) { EnsureAround(*s, c.c, 1.0f); std::lock_guard<std::mutex> lock(cmdMx); commands.push_back(c); return; } // after classification
        // A full footprint widens (to 4x) before any water is given up (#766: the rest was dropped silently, and Pour had
        // already reported it poured); whatever still finds no room is counted (PourDropped).
        const int hw0 = (int)(0.5f / WaterField::CellSize()), layers = (int)ceilf(6.0f / WaterField::CellSize()); // half-width (cells), layers
        for (int hw = hw0, ring = 0; ring < 4 && left > 0.0f; ++ring, hw = std::max(hw + 1, hw * 2))
            for (int k = 0; k < layers && left > 0.0f; ++k)
                for (int dj = -hw; dj <= hw && left > 0.0f; ++dj)
                    for (int di = -hw; di <= hw && left > 0.0f; ++di) {
                        if (!s->cell(ci + di, cj + dj, ck + k, b, i)) continue;
                        if (b->flags[i] & kSolid) continue;
                        if (b->fill[i] > 0 && CellKind(b->flags[i]) != c.fluid) continue; // another fluid: spill past it
                        TakeKind(b, i, (uint8_t)(c.fluid << kKindShift));
                        const int add = (int)std::min(left, (float)std::max(0, Cap(b, i) - b->fill[i]));
                        if (add <= 0) continue;
                        b->fill[i] = (uint8_t)(b->fill[i] + add);
                        GpuDelta(b, i, add);
                        left -= add;
                        Touch(b);
                        if (dbgAxSet && EnvFlag("DBTRACE")) { // TEMP #680: pour water placed at or past the wall
                            const Vector3 d = Vector3Subtract(s->cellCentre(ci + di, cj + dj, ck + k), dbgAxO);
                            const float along = Vector3DotProduct(d, dbgAxX);
                            if (along > 0.0f) printf("dbtrace: pour ring %d (hw %d) put %d units at along %.2f m, h %.2f m (pour at along %.2f)\n", ring, hw, add, along,
                                                       Vector3DotProduct(d, dbgAxU), Vector3DotProduct(Vector3Subtract(c.c, dbgAxO), dbgAxX));
                        }
                    }
        if (left > 0.0f) pourDroppedUnits += (long long)left;
        return;
    }
    // Take: from the wettest cells within the radius (sources are free: they refill).
    float want = c.m3 * 255.0f / WaterField::CellVolume();
    const int R = (int)ceilf(c.r / WaterField::CellSize());
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
            JobSystem::tagNext(JobTag::Fluid); WorldJobs().submit(-3.0f, [w, &left] { Classify(*w->s, *w->b, *w->ref, w->faces); left.fetch_sub(1); });
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
        for (auto& s : sites) for (Site::Spring& sp : s->springs) StepSpring(*s, sp); // the cells: the GPU (GpuTick) or CpuTick below
        static const bool levelOff = getenv("RAYLIBGAME_FLUID_LEVEL") && atoi(getenv("RAYLIBGAME_FLUID_LEVEL")) == 0; // A/B: no levelling pass
        if (!levelOff && stepCount % LevelEvery() == 0) for (auto& s : sites) LevelSite(*s);
    }
    // CPU mode: the GPU's step on the CPU (#654), at its substep count (120 a second).
    if (!gpuMode) CpuTick((int)lround(2 * kSubsteps * (kTickSeconds / 0.05)) * std::max(1, speed.load()));
    const auto t1 = std::chrono::steady_clock::now();
    // Publish: a new snapshot of the water field (waterfield.h) sharing every unchanged brick; the changed bricks'
    // bounds go out with it, so the chunks there redraw their water (the chunk mesher draws all fluid, #360).
    Fluid::Stats st;
    const WaterField::SnapshotPtr prev = WaterField::Current();
    std::shared_ptr<WaterField::Snapshot> snap;
    std::vector<WaterField::Dirty> dirty;
    // Per brick on the workers (#676: one thread took ~26 ms a tick over 3500 quarter-metre bricks): its counts, settling,
    // easing, film hold and - when it changed - its published copy and the cells to redraw. Each touches only its own brick;
    // the sums, the snapshot and the redraw boxes are merged in order below.
    struct PubOut {
        Site* s; int64_t key; Brick* b;
        std::shared_ptr<WaterField::Brick> pb;
        int box[8][7]; int boxes = 0; // per 8^3 octant: cell range (x0 y0 z0 x1 y1 z1, brick-local) and the largest change
        bool settledNow = false, flip = false;
    };
    static std::vector<PubOut> outs;
    outs.clear();
    for (auto& s : sites) {
        ++st.sites;
        for (auto& kv : s->bricks) {
            ++st.bricks;
            if (kv.second->ready) outs.push_back({ s.get(), kv.first, kv.second.get() });
        }
    }
    auto publishBrick = [&](PubOut& out) {
            Brick& b = *out.b;
            if (b.changed || !b.statSet) { // kinds in the brick (bit k), wet cells, volumes: recounted when it changed
                uint8_t kindsHere = 0;
                int wet = 0, vol = 0, other = 0, solid = 0;
                for (int i = 0; i < B3; ++i)
                    if (b.fill[i] && !(b.flags[i] & kSolid)) {
                        const int k = CellKind(b.flags[i]);
                        kindsHere |= (uint8_t)(1 << k);
                        if (b.flags[i] & kSource) continue;
                        ++wet;
                        (k ? other : vol) += b.fill[i];
                    }
                for (int i = 0; i < B3; ++i) solid += (b.flags[i] & kSolid) ? 1 : 0;
                b.statSolid = solid; b.statKinds = kindsHere; b.statWet = (float)wet; b.statVolume = vol / 255.0f * WaterField::CellVolume(); b.statOther = other / 255.0f * WaterField::CellVolume(); b.statSet = true;
            }
            const uint8_t present = b.statKinds;
            // Resting water sleeps (#366): an awake brick whose cells all stayed within kSettleDiff of where they were
            // kSettleSteps ago is put to sleep, though films still trickle back and forth in it (a drained pond kept
            // 30-50 bricks awake for minutes, and their changes kept the chunks redrawing, #381). Water moving into it
            // wakes it again (Touch); a real flow changes some cell by more than that within the window.
            if (b.active) {
                ++b.winTicks;
                if (b.checkStep < 0) { b.check = b.fill; b.checkStep = stepCount; b.winTicks = b.winFlow = 0; b.checkVol = b.statVolume + b.statOther; }
                else if (stepCount - b.checkStep >= kSettleSteps) {
                    int maxd = 0;
                    for (int i = 0; i < B3 && maxd <= 2 * kSettleDiff; ++i) maxd = std::max(maxd, abs((int)b.fill[i] - (int)b.check[i]));
                    // Levels held: at rest it sleeps; with water still passing through (the GPU changed it in half its ticks
                    // or more) it is steady (#674: keeps running, read back once a second). Levels moved: neither.
                    // Passing through means a net current too (#680): at 0.25 m cells resting water jitters a fill unit nearly
                    // every tick, so "the GPU changed it" held ~140 bricks of a still pond awake as steady flow. The net current
                    // (the mean of its columns' signed flows, #677's colFlow, cm/s) averages jitter out; a real stream moves
                    // its columns the same way. No current published yet: the old rule.
                    bool flowing = b.winTicks > 0 && b.winFlow * 2 >= b.winTicks;
                    if (flowing && kSteadyNet > 0.0f && !b.colFlow.empty()) {
                        double sx = 0.0, sy = 0.0; int nCol = 0;
                        for (size_t k = 0; k + 1 < b.colFlow.size(); k += 2)
                            if (b.colFlow[k] || b.colFlow[k + 1]) { sx += b.colFlow[k]; sy += b.colFlow[k + 1]; ++nCol; }
                        flowing = nCol > 0 && sqrt(sx * sx + sy * sy) / nCol >= kSteadyNet;
                    }
                    // Its volume held too (#680): a big pool draining slowly through a hole (a few mm a second) moved every cell
                    // less than the tolerance and fell asleep over the water it still owed; jitter moves units, not volume.
                    // At most 2 mm over the brick's area.
                    const float side = (float)B * WaterField::CellSize();
                    const bool volumeHeld = fabsf(b.statVolume + b.statOther - b.checkVol) <= 0.002f * side * side;
                    if (maxd <= 2 * kSettleDiff && flowing && volumeHeld && gpuMode) b.steady = true; // a moving surface wobbles more than a still one
                    else if (maxd <= kSettleDiff && volumeHeld) { b.active = false; b.idle = kSleepSteps + 1; b.settled = true; b.steady = false; out.settledNow = true; }
                    else b.steady = false;
                    b.check = b.fill; b.checkStep = stepCount; b.winTicks = b.winFlow = 0; b.checkVol = b.statVolume + b.statOther;
                }
            } else { b.checkStep = -1; b.steady = false; }
            // The drawn water eases toward the simulated: 40% of the way a tick (25 ms), so a cell the water reaches fills in
            // over ~0.15 s instead of appearing whole (#366: water read as packets jumping a metre at a time). Solid cells follow
            // at once. A brick keeps publishing while it eases. Only the mesher reads the eased copy (WaterField::Brick::drawn);
            // queries, buckets and pumps read the simulated fill, so the ease never creates or loses water.
            if (b.changed || b.easing) {
                bool still = false;
                if (!b.dispSet) { b.disp = b.fill; b.dispSet = true; }
                else
                    for (int i = 0; i < B3; ++i) {
                        const int f = b.fill[i], d = b.disp[i];
                        if (f == d) continue;
                        const int diff = f - d;
                        // Draining eases slower (1/8 a tick, a unit a tick at the end): a thin film's cell that empties
                        // for a tick or two as the water moves on in packets showed a hole (W4c: quarter-metre streams
                        // flickered near their spring); a film that really dries fades within ~0.25 s.
                        int step = diff > 0 ? diff * 2 / 5 : diff / 8;
                        if ((b.flags[i] & kSolid) || (diff > 0 && diff <= 3)) step = diff;
                        else if (step == 0) step = diff > 0 ? 1 : -1;
                        b.disp[i] = (uint8_t)(d + step);
                        still = still || b.disp[i] != f;
                    }
                if (b.easing || still) b.changed = true;
                b.easing = still;
            }
            // Its thin films are drawn while it moves and for kFilmHoldSteps after (#660: bricks in a steady stream fall asleep
            // and wake again every second or so, and their films blinked out and back each time).
            if (b.active || b.easing) b.movingStep = stepCount;
            const bool movingNow = stepCount - b.movingStep < FilmHoldSteps();
            const bool movedFlip = b.drawnSet && b.drawnMoving != movingNow; // its thin films show or hide (#649)
            out.flip = movedFlip;
            // The current changed (#677): a steady stream's levels hold, so it redraws on its flow - some column 20 cm/s off
            // what was drawn, at most every 2 s (80 ticks).
            bool flowRedraw = false;
            if (!b.colFlow.empty() && stepCount - b.flowStep >= 80) {
                if (b.drawnFlow.size() != b.colFlow.size()) flowRedraw = true;
                else for (size_t k = 0; k < b.colFlow.size() && !flowRedraw; ++k) flowRedraw = abs((int)b.colFlow[k] - (int)b.drawnFlow[k]) > 20;
            }
            if (!b.changed && !movedFlip && !flowRedraw) return;
            b.changed = false;
            if (flowRedraw) { b.drawnFlow = b.colFlow; b.flowStep = stepCount; }
            out.pb = std::make_shared<WaterField::Brick>();
            WaterField::Brick* pb = out.pb.get();
            pb->fill = b.fill; pb->drawn = b.disp; pb->flags = b.flags; pb->shut = b.shut; pb->kinds = present; pb->moving = movingNow; // queries read the simulated water, the mesher the eased
            // The drawn current (#677, #695): per column the solver's flow at its surface, cm/s along the site's east and
            // north - the top water cell's face flows over the column's wetted depth (a column's top cell carries the
            // column's flow). GPU: fluid_flowpack.comp, read back with the cells; CPU kernel: from its faces here.
            if (!b.colFlow.empty()) pb->flow = b.colFlow;
            else if (b.faces.size() == (size_t)B3 && !(getenv("RAYLIBGAME_FLUID_REALFLOW") && atoi(getenv("RAYLIBGAME_FLUID_REALFLOW")) == 0)) {
                pb->flow.assign(512, 0);
                const float cs = WaterField::CellSize();
                for (int y = 0; y < B; ++y)
                    for (int x = 0; x < B; ++x) {
                        float depth = 0.0f; // the water from the column's bottom up to this cell (m), within the brick
                        for (int z = 0; z < B; ++z) {
                            const int i = x + B * y + B * B * z;
                            if (!b.fill[i] || (b.flags[i] & kSolid)) { depth = 0.0f; continue; }
                            depth += b.fill[i] / 255.0f * cs;
                            const int above = z + 1 < B ? i + B * B : -1;
                            if (above >= 0 && b.fill[above] && !(b.flags[above] & kSolid)) continue; // not the surface
                            const float area = cs * std::max(0.02f, depth);
                            pb->flow[2 * (x + B * y)] = (int16_t)std::clamp(b.faces[i].x / area * 100.0f, -32000.0f, 32000.0f);
                            pb->flow[2 * (x + B * y) + 1] = (int16_t)std::clamp(b.faces[i].y / area * 100.0f, -32000.0f, 32000.0f);
                        }
                    }
            }
            // A cell's surface reaches its neighbours' samples (solid cells carry it, empty ones look below). Only a
            // visible change redraws the chunks, and only around the cells that changed: some cell 2.5 cm+ off what
            // they last drew, or its solidity or fluid kind (a body easing by a 255th at a time kept ~2000 chunk
            // redraws a second going, #360). Those cells are drawn anew; the rest keep what was drawn.
            // Per 8^3 octant of the brick (#381: one box over scattered changes marked a whole brick's chunks).
            int c0[8][3], c1[8][3], mag[8] = {}; // mag: the largest change in each (coarse chunks wait for a bigger one, world.cpp)
            for (int o = 0; o < 8; ++o) { c0[o][0] = c0[o][1] = c0[o][2] = B; c1[o][0] = c1[o][1] = c1[o][2] = -1; }
            for (int i = 0; i < B3; ++i) {
                const bool v = !b.drawnSet || abs((int)b.disp[i] - (int)b.drawn[i]) >= 6 || ((b.flags[i] ^ b.drawnFlags[i]) & kSolid) ||
                               (b.disp[i] && ((b.flags[i] ^ b.drawnFlags[i]) & kKindMask)) || (movedFlip && b.disp[i] > 0 && b.disp[i] < 64) ||
                               (flowRedraw && b.disp[i] > 0 && b.disp[i] < 255); // the surface cells, for their current
                if (!v) continue;
                const int x = i % B, y = (i / B) % B, z = i / (B * B);
                const int o = (x >> 3) | ((y >> 3) << 1) | ((z >> 3) << 2);
                mag[o] = std::max(mag[o], (!b.drawnSet || ((b.flags[i] ^ b.drawnFlags[i]) & (kSolid | kKindMask)) || flowRedraw) ? 255 : abs((int)b.disp[i] - (int)b.drawn[i]));
                c0[o][0] = std::min(c0[o][0], x); c0[o][1] = std::min(c0[o][1], y); c0[o][2] = std::min(c0[o][2], z);
                c1[o][0] = std::max(c1[o][0], x); c1[o][1] = std::max(c1[o][1], y); c1[o][2] = std::max(c1[o][2], z);
                b.drawn[i] = b.disp[i]; b.drawnFlags[i] = b.flags[i];
            }
            b.drawnSet = true;
            b.drawnMoving = movingNow;
            for (int o = 0; o < 8; ++o) {
                if (c1[o][0] < 0) continue;
                int* bx = out.box[out.boxes++];
                bx[0] = c0[o][0]; bx[1] = c0[o][1]; bx[2] = c0[o][2]; bx[3] = c1[o][0]; bx[4] = c1[o][1]; bx[5] = c1[o][2]; bx[6] = mag[o];
            }
    };
    {
        const int n = (int)outs.size();
        const int parts = std::max(1, std::min(n / 64, (int)std::thread::hardware_concurrency()));
        std::atomic<int> left{ parts - 1 };
        auto run = [&](int part) { for (int k = part; k < n; k += parts) publishBrick(outs[k]); };
        for (int part = 1; part < parts; ++part) {
            JobSystem::tagNext(JobTag::Fluid);
            WorldJobs().submit(-3.0f, [&run, &left, part] { run(part); left.fetch_sub(1); });
        }
        run(0);
        while (left.load() > 0) std::this_thread::yield();
    }
    for (PubOut& out : outs) {
        const Brick& b = *out.b;
        if (b.active) ++st.activeBricks;
        if (b.steady) ++st.steadyBricks;
        if (b.basin) ++st.basinBricks;
        st.wetCells += b.statWet; st.volume += b.statVolume; st.otherVolume += b.statOther;
        settledTotal += out.settledNow ? 1 : 0;
        filmFlips += out.flip ? 1 : 0;
        if (!out.pb) continue;
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
        WaterField::Site& ps = snap->sites[out.s->index];
        ps.bricks[out.key] = std::move(out.pb);
        Vector3 lo, hi;
        ps.brickBounds(b.bi, b.bj, b.bk, lo, hi);
        ps.lo = Vector3Min(ps.lo, lo); ps.hi = Vector3Max(ps.hi, hi);
        for (int k = 0; k < out.boxes; ++k) {
            const int* bx = out.box[k];
            Vector3 dlo, dhi;
            ps.cellBounds(b.bi * B + bx[0] - 1, b.bj * B + bx[1] - 1, b.bk * B + bx[2] - 1, b.bi * B + bx[3] + 2, b.bj * B + bx[4] + 2, b.bk * B + bx[5] + 2, dlo, dhi);
            dirty.push_back({ out.s->body, dlo, dhi, bx[6] });
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
        st.filmFlips = filmFlips;
        st.gpuSlots = slotsUsed; st.gpuCapacity = G_Capacity(); st.slotReleases = slotReleases;
        st.dbg[0] = dbgDown; st.dbg[1] = dbgSide; st.dbg[2] = dbgSrc; st.dbg[3] = dbgLevel; dbgDown = dbgSide = dbgSrc = dbgLevel = 0;
        st.settled = (int)settledTotal;
        pubStats = st;
        awakeBricksLast = st.activeBricks; // the levelling interval (LevelEvery)
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
    ComputeBuffer faceP, faceN;            // face flows (#642, fluid_step.comp): applied (a vec4 a cell) / new (two vec4 a cell)
    // Async readback (#669): the dispatch's flags and its bricks' cells packed (fluid_pack.comp) are copied into a persistent
    // mapped buffer behind a fence and read the next tick; without buffer storage, the old per-brick downloads.
    ComputeProgram pack, flowPack; // flowPack: per column the solver's current, for drawing (#677)
    ComputeBuffer flowPacked;
    bool realFlow = false;      // the drawn current from the solver (RAYLIBGAME_FLUID_REALFLOW=0: the surface-slope estimate)
    ComputeBuffer packed;
    unsigned rbBuf = 0; unsigned char* rbPtr = nullptr; size_t rbCap = 0; void* rbFence = nullptr; size_t rbFlagsBytes = 0, rbSinkAt = 0, rbSinkBytes = 0; bool sinkOn = false; // sinkOn: some release measures its outflow (#656); rbSinkAt: its counters in the mapped buffer
    size_t rbFlowAt = 0, rbFlowBytes = 0; // the drawn current's columns in the mapped buffer (#677), after the sink counters
    bool asyncRb = false;
    // The dry-cell skip (#675, fluid_groups.comp): per slot 3 substeps x 64 group bits (wet), the listed groups + dispatch args.
    ComputeProgram groups;
    ComputeBuffer mask, glist;
    ComputeBuffer sinkOut; // per slot: fill units its sinks took in (river flow, #656; binding 12)
    bool groupsOn = false;
    int maskCur = 0;
    int capacity = 0;                      // bricks
    std::vector<Brick*> owner;             // slot -> brick
    std::vector<Site*> ownerSite;
    std::vector<int> lastList;             // slots processed by the last dispatch
    int lastN = 0;                         // substeps in it
    bool nbrDirty = false;
    // Flow readout (#392): watched bricks' flux blocks (fluid_step.comp bindings 6-7).
    ComputeBuffer watch, flux;
    std::vector<int> watchOf;              // slot -> flux block, -1 (uploaded to `watch`)
    int fluxBlocks = 0;
    struct WatchCell { int block, cell; Vector3 e, n, u; float wx, wy, wz; }; // w: 1 / the ball's cells on its line along each axis
    std::vector<std::vector<WatchCell>> watchCells; // per flow watch: its cells as resolved for the last dispatch
} G;
int G_Capacity() { return G.capacity; }

// Flow watches (#392): Fluid::FlowAt asks for the flow around a point; GpuTick resolves each to cells, reads their flux
// every tick and publishes a smoothed value. A watch nobody asked about for 3 s is dropped.
struct FlowWatch { uint64_t body = 0; Vector3 local{}; float r = 1.0f; double lastAsk = 0.0; Vector3 value{}; };
std::vector<FlowWatch> flowWatches;
std::mutex flowMx;

constexpr unsigned kCellsPerBrick = (unsigned)B3;
constexpr unsigned kFacePBytes = 8, kFaceNBytes = 16; // half-precision face flows a cell (#675; were 16 and 32)

bool GpuInit() {
    if (G.init) return G.ok;
    G.init = true;
    G.ok = Compute_Available() && G.step.load("resources/compute/fluid_step.comp") && G.delta.load("resources/compute/fluid_delta.comp");
    if (G.ok) {
        G.list.create(64 * 4); G.deltas.create(64 * 16); G.inS.create(64 * 4); G.flags.create(64 * 4); G.nbrs.create(64 * 24);
        G.a.create(64 * kCellsPerBrick * 4); G.b.create(64 * kCellsPerBrick * 4);
        const std::vector<float> noFlow((size_t)64 * kCellsPerBrick * 8, 0.0f);
        G.faceP.create(64 * kCellsPerBrick * kFacePBytes, noFlow.data()); G.faceN.create(64 * kCellsPerBrick * kFaceNBytes, noFlow.data());
        G.capacity = 64;
        G.owner.assign(64, nullptr); G.ownerSite.assign(64, nullptr);
        const std::vector<unsigned> zeros(64, 0u);
        G.flags.upload(zeros.data(), 64 * 4);
        G.sinkOut.create(64 * 4, zeros.data());
        // Async readback (#669): buffer storage for the persistent mapping, the pack pass; else the old downloads.
        G.asyncRb = !EnvFlag("RAYLIBGAME_FLUID_SYNC_READBACK") && GL_Caps().bufferStorage && glad_glFenceSync && G.pack.load("resources/compute/fluid_pack.comp");
        if (G.asyncRb) G.packed.create(64 * kCellsPerBrick * 2);
        G.realFlow = G.asyncRb && !(getenv("RAYLIBGAME_FLUID_REALFLOW") && atoi(getenv("RAYLIBGAME_FLUID_REALFLOW")) == 0) && G.flowPack.load("resources/compute/fluid_flowpack.comp");
        if (G.realFlow) G.flowPacked.create(64 * 256 * 4);
        G.groupsOn = !(getenv("RAYLIBGAME_FLUID_GROUPS") && atoi(getenv("RAYLIBGAME_FLUID_GROUPS")) == 0) && glad_glDispatchComputeIndirect && G.groups.load("resources/compute/fluid_groups.comp");
        if (G.groupsOn) {
            const std::vector<unsigned> zero((size_t)(8 + 64 * 64), 0u);
            G.mask.create(64 * 6 * 4, zero.data()); G.glist.create((8 + 64 * 64) * 4, zero.data());
        }
    }
    TraceLog(LOG_WARNING, "Fluid: %s%s", G.ok ? "GPU compute" : "CPU (no compute)", G.asyncRb ? ", async readback" : "");
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
    grow(G.faceP, G.capacity * kCellsPerBrick * kFacePBytes, cap * kCellsPerBrick * kFacePBytes);
    grow(G.faceN, G.capacity * kCellsPerBrick * kFaceNBytes, cap * kCellsPerBrick * kFaceNBytes);
    grow(G.flags, G.capacity * 4, cap * 4);
    grow(G.sinkOut, G.capacity * 4, cap * 4);
    if (G.groupsOn) {
        grow(G.mask, G.capacity * 6 * 4, cap * 6 * 4);
        const std::vector<unsigned> zero((size_t)(8 + cap * 64), 0u);
        G.glist.destroy(); G.glist.create((unsigned)zero.size() * 4, zero.data());
    }
    G.inS.resize(cap * 4);
    G.nbrs.resize(cap * 24);
    const std::vector<unsigned> zeros(cap - G.capacity, 0u);
    G.flags.upload(zeros.data(), (unsigned)zeros.size() * 4, G.capacity * 4);
    G.sinkOut.upload(zeros.data(), (unsigned)zeros.size() * 4, G.capacity * 4);
    G.owner.resize(cap, nullptr); G.ownerSite.resize(cap, nullptr);
    G.capacity = cap;
    G.nbrDirty = true;
}

void GpuUploadBrick(const Brick& b) {
    static std::vector<unsigned> packed(kCellsPerBrick);
    for (int i = 0; i < B3; ++i) packed[i] = b.fill[i] | ((unsigned)b.flags[i] << 8) | ((unsigned)b.target[i] << 16) | ((unsigned)b.shut[i] << 24); // shut in bits 24-31 (#454)
    G.a.upload(packed.data(), kCellsPerBrick * 4, b.slot * kCellsPerBrick * 4);
    G.b.upload(packed.data(), kCellsPerBrick * 4, b.slot * kCellsPerBrick * 4);
    static const std::vector<float> noFlow((size_t)kCellsPerBrick * 4, 0.0f); // a (re)classified brick starts at rest (#642)
    G.faceP.upload(noFlow.data(), kCellsPerBrick * kFacePBytes, b.slot * kCellsPerBrick * kFacePBytes);
    if (G.groupsOn) { static const unsigned all[6] = { ~0u, ~0u, ~0u, ~0u, ~0u, ~0u }; G.mask.upload(all, 24, b.slot * 24); } // every group listed once: the next substeps narrow it
}

// Flow readout (#392). Read: the flux the last dispatch added in the watched cells, summed per watch into body-local
// m^3/s (two faces per axis: halved), each cell weighted 1 / the cells of the ball on its line along that axis: the mean flux
// through the ball's parallel slices, so a watch spanning a channel reads its discharge at any cell size (W4c: the old
// division by the diameter under-read 3x with quarter-metre cells), smoothed over ~0.5 s; the buffer is cleared for the next dispatch.
void GpuFlowRead() {
    if (G.fluxBlocks <= 0 || G.watchCells.empty()) return;
    std::vector<int> fx((size_t)G.fluxBlocks * kCellsPerBrick * 4);
    G.flux.download(fx.data(), (unsigned)(fx.size() * 4));
    const float perSecond = WaterField::CellVolume() / (255.0f * (float)kTickSeconds * (float)std::max(1, speed.load())); // fill units a tick -> m^3/s
    std::lock_guard<std::mutex> lk(flowMx);
    for (size_t w = 0; w < G.watchCells.size() && w < flowWatches.size(); ++w) {
        Vector3 sum{};
        for (const auto& c : G.watchCells[w]) {
            const int* f = &fx[((size_t)c.block * kCellsPerBrick + c.cell) * 4];
            sum = Vector3Add(sum, Vector3Add(Vector3Add(Vector3Scale(c.e, (float)f[0] * c.wx), Vector3Scale(c.n, (float)f[1] * c.wy)), Vector3Scale(c.u, (float)f[2] * c.wz)));
        }
        const Vector3 q = Vector3Scale(sum, 0.5f * perSecond); // each line through the ball weighted 1 / its cells: the flux through its disc (W4c)
        flowWatches[w].value = Vector3Lerp(flowWatches[w].value, q, 0.1f);
    }
    std::fill(fx.begin(), fx.end(), 0);
    G.flux.upload(fx.data(), (unsigned)(fx.size() * 4));
}

// Resolve: each live watch to the cells within its radius in bricks on the GPU; one flux block per watched brick.
void GpuFlowResolve() {
    {   // Nothing watched now or last time, buffers sized: nothing to rebuild (the common case: no meters or wheels).
        std::lock_guard<std::mutex> lk(flowMx);
        if (flowWatches.empty() && G.fluxBlocks == 0 && (int)G.watchOf.size() == G.capacity) return;
    }
    std::vector<int> watchOf(G.capacity, -1);
    std::vector<std::vector<GpuState::WatchCell>> cells;
    int blocks = 0;
    {
        std::lock_guard<std::mutex> lk(flowMx);
        const double now = GetTime();
        flowWatches.erase(std::remove_if(flowWatches.begin(), flowWatches.end(), [&](const FlowWatch& w) { return now - w.lastAsk > 3.0; }), flowWatches.end());
        cells.resize(flowWatches.size());
        for (size_t w = 0; w < flowWatches.size(); ++w) {
            const FlowWatch& fw = flowWatches[w];
            Site* s = SiteFor(fw.body, fw.local, false);
            if (!s) continue;
            int ci, cj, ck;
            s->cellOf(fw.local, ci, cj, ck);
            const float rc = fw.r / WaterField::CellSize(); // the watch radius in cells
            const int R = (int)ceilf(rc);
            for (int dk = -R; dk <= R; ++dk)
                for (int dj = -R; dj <= R; ++dj)
                    for (int di = -R; di <= R; ++di) {
                        if (di * di + dj * dj + dk * dk > rc * rc) continue;
                        Brick* b; int i;
                        if (!s->cell(ci + di, cj + dj, ck + dk, b, i) || b->slot < 0) continue;
                        if (watchOf[b->slot] < 0) watchOf[b->slot] = blocks++;
                        auto line = [&](int a, int b2) { return 1.0f / (2.0f * floorf(sqrtf(std::max(0.0f, rc * rc - (float)(a * a + b2 * b2)))) + 1.0f); };
                        cells[w].push_back({ watchOf[b->slot], i, s->e, s->n, s->u, line(dj, dk), line(di, dk), line(di, dj) });
                    }
        }
    }
    if ((int)G.watchOf.size() != G.capacity || watchOf != G.watchOf) {
        if (G.watch.bytes < (unsigned)G.capacity * 4) { if (G.watch.id) G.watch.destroy(); G.watch.create((unsigned)G.capacity * 4); }
        G.watch.upload(watchOf.data(), (unsigned)G.capacity * 4);
        G.watchOf = watchOf;
    }
    const unsigned need = (unsigned)std::max(1, blocks) * kCellsPerBrick * 16;
    if (G.flux.bytes < need) { // grown (cleared)
        if (G.flux.id) { G.flux.destroy(); } G.flux.create(need);
        const std::vector<int> zeros(need / 4, 0);
        G.flux.upload(zeros.data(), need);
    }
    G.fluxBlocks = blocks;
    G.watchCells = std::move(cells);
}

void GpuTick(int n) {
    // 0. Released reaches' outflow (#656): what their sinks took in, credited to each sink brick's release. It rides the
    // async readback (copied with the flags each dispatch, read a tick later); without it, a stalling download once a
    // second (with 1800 bricks at 0.25 m cells that stall was a 35 ms frame every second).
    auto creditSinks = [](const unsigned* units, int slots) {
        for (int slot = 0; slot < slots && slot < (int)G.owner.size(); ++slot) {
            const Brick* br = G.owner[slot];
            Site* st = G.ownerSite[slot];
            if (units[slot] && br && st && br->sinkRelease >= 0 && br->sinkRelease < (int)st->releases.size()) st->releases[br->sinkRelease].sinkUnits += units[slot];
        }
    };
    G.sinkOn = false;
    for (auto& s : sites) for (const Site::Release& r : s->releases) G.sinkOn = G.sinkOn || r.piece >= 0;
    static int sinkTick = 0;
    if (!G.asyncRb && G.sinkOn && ++sinkTick % 20 == 0 && G.capacity > 0) {
        static std::vector<unsigned> units, zeros;
        units.resize(G.capacity); zeros.assign(G.capacity, 0u);
        G.sinkOut.download(units.data(), G.capacity * 4);
        G.sinkOut.upload(zeros.data(), G.capacity * 4);
        creditSinks(units.data(), G.capacity);
    }
    // 1. Results of the last dispatch: which bricks changed (and want neighbours), then their cells.
    if (!G.lastList.empty()) {
        // The flags and packed cells: from the mapped buffer the last dispatch copied into (#669: no stall; the copy was
        // queued a tick ago), else downloaded.
        std::vector<unsigned> flagsV;
        const unsigned* flags = nullptr;
        const uint32_t* packed = nullptr;
        if (G.asyncRb && G.rbFence) {
            glClientWaitSync((GLsync)G.rbFence, GL_SYNC_FLUSH_COMMANDS_BIT, 100000000ull); // at most 100 ms
            glDeleteSync((GLsync)G.rbFence); G.rbFence = nullptr;
            flags = (const unsigned*)G.rbPtr;
            packed = (const uint32_t*)(G.rbPtr + G.rbFlagsBytes);
            if (G.rbSinkBytes) creditSinks((const unsigned*)(G.rbPtr + G.rbSinkAt), (int)(G.rbSinkBytes / 4)); // #656
            // The drawn current (#677): per column the solver's flow, kept on the brick until it next publishes.
            if (G.rbFlowBytes) {
                const uint32_t* fp = (const uint32_t*)(G.rbPtr + G.rbFlowAt);
                const size_t n = std::min(G.lastList.size(), G.rbFlowBytes / (256 * 4));
                for (size_t li = 0; li < n; ++li) {
                    Brick* br = G.owner[G.lastList[li]];
                    if (!br) continue;
                    br->colFlow.resize(512);
                    for (int c = 0; c < 256; ++c) { const uint32_t v = fp[li * 256 + c]; br->colFlow[2 * c] = (int16_t)(v & 0xFFFFu); br->colFlow[2 * c + 1] = (int16_t)(v >> 16); }
                }
            }
        } else {
            flagsV.resize(G.capacity);
            G.flags.download(flagsV.data(), G.capacity * 4);
            flags = flagsV.data();
        }
        const int flagSlots = packed ? (int)(G.rbFlagsBytes / 4) : G.capacity;
        static std::vector<unsigned> cells(kCellsPerBrick);
        static const int d[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
        // A brick the CPU reclassified since that dispatch (upload pending) or is reclassifying now (not ready: a worker reads
        // its fill) keeps the CPU state: the GPU's older fill would put water back into cells that became solid (PR #129 review).
        // A steady brick (#674) is read back on one tick in kSteadySampleTicks (staggered by slot); it keeps running between.
        static long long rbTick = 0;
        ++rbTick;
        auto sampled = [&](int slot, const Brick* br) { return !br->steady || (rbTick + slot) % kSteadySampleTicks == 0; };
        auto takes = [&](size_t li) {
            const int slot = G.lastList[li];
            const Brick* br = G.owner[slot];
            return br && G.ownerSite[slot] && slot < flagSlots && (flags[slot] & 1u) && br->ready && !br->upload && sampled(slot, br);
        };
        // The packed cells into the bricks, on the workers (#669: ~3000 bricks x 4096 cells on the main thread took ~10 ms a tick);
        // drift: how far a settled brick's cells moved from where it settled.
        std::vector<int> drift(G.lastList.size(), 0);
        if (packed) {
            const int n = (int)G.lastList.size();
            const int parts = std::max(1, std::min(n, (int)std::thread::hardware_concurrency()));
            auto unpack = [&](int part) {
                for (int li = part; li < n; li += parts) {
                    if (!takes((size_t)li)) continue;
                    Brick* br = G.owner[G.lastList[li]];
                    const uint32_t* pc = packed + (size_t)li * (kCellsPerBrick / 2);
                    int dr = 0;
                    for (int i = 0; i < B3; ++i) {
                        const uint32_t v = (pc[i >> 1] >> ((i & 1) * 16)) & 0xFFFFu; // fill | kind << 8 (fluid_pack.comp)
                        br->fill[i] = (uint8_t)(v & 255u);
                        br->flags[i] = (uint8_t)((br->flags[i] & ~kKindMask) | (((v >> 8) & 7u) << kKindShift));
                        if (br->settled) dr = std::max(dr, abs((int)br->fill[i] - (int)br->check[i]));
                    }
                    drift[li] = dr;
                }
            };
            std::atomic<int> left{ parts - 1 };
            for (int part = 1; part < parts; ++part) {
                JobSystem::tagNext(JobTag::Fluid);
                WorldJobs().submit(-3.0f, [&unpack, &left, part] { unpack(part); left.fetch_sub(1); });
            }
            unpack(0);
            while (left.load() > 0) std::this_thread::yield();
        }
        for (size_t li = 0; li < G.lastList.size(); ++li) {
            const int slot = G.lastList[li];
            Brick* br = G.owner[slot];
            Site* site = G.ownerSite[slot];
            if (!br || !site || slot >= flagSlots) continue;
            const unsigned fl = flags[slot];
            if (fl & 1u) ++br->winFlow; // the GPU changed it this tick (steady flow, #674)
            if (br->steady && !takes(li)) {
                // Steady, not read this tick: it keeps running (no idle count: the water still passes through).
            } else if (br->steady) {
                br->changed = true; // read back: published (its levels are checked at the end of the settle window)
            } else if (takes(li)) {
                if (!packed) {
                    G.a.download(cells.data(), kCellsPerBrick * 4, slot * kCellsPerBrick * 4);
                    for (int i = 0; i < B3; ++i) {
                        br->fill[i] = (uint8_t)(cells[i] & 255u);
                        br->flags[i] = (uint8_t)((br->flags[i] & ~kKindMask) | (((cells[i] >> 8) & 255u) & kKindMask));
                        if (br->settled) drift[li] = std::max(drift[li], abs((int)br->fill[i] - (int)br->check[i]));
                    }
                }
                if (br->settled && drift[li] <= kSettleDiff) { br->changed = true; continue; } // a resting brick at the edge of a dispatch: stays asleep (#366)
                if (br->settled) { ++wakeGpu; br->settled = false; }
                Touch(br);
                for (int fc = 0; fc < 6; ++fc) if (Brick* o = br->nbr[fc]; o && RunsBeside(*o, fc)) { o->active = true; o->idle = 0; }
            } else if (br->active && (br->idle += G.lastN) > kSleepSteps) {
                br->active = false;
            }
            for (int face = 0; face < 6; ++face)
                if (fl & (2u << face)) site->ensure(br->bi + d[face][0], br->bj + d[face][1], br->bk + d[face][2]);
        }
        static std::vector<unsigned> zeros;
        zeros.assign(G.capacity, 0u);
        G.flags.upload(zeros.data(), G.capacity * 4);
        G.lastList.clear();
    }
    GpuTimer_Mark("fluid-read"); // benches: the readback above (#669)
    GpuFlowRead(); // 1b. the watched cells' flux since the last read (#392)
    // 2. Slots for classified bricks; whole bricks up.
    int used = 0;
    for (Brick* o : G.owner) used += o ? 1 : 0;
    // Only bricks near awake water hold a slot (#675): the awake ones, and two rings of neighbours (the step runs a ring and
    // reads across its faces; the column depth looks down further). One not needed for kSlotKeepTicks gives its slot back:
    // its CPU copy is current (results come back every tick it runs), and it goes up whole when it is needed again, without
    // waking (a read-only ring would wake the next ring: a whole lake). A settled lake then holds no GPU memory.
    static long long slotTick = 0;
    ++slotTick;
    static const bool keepAll = EnvFlag("RAYLIBGAME_FLUID_KEEP_SLOTS"); // A/B: every brick keeps its slot (before #675)
    for (auto& s : sites)
        for (auto& kv : s->bricks) {
            Brick& br = *kv.second;
            if (!br.ready || !(br.active || br.upload || keepAll)) continue;
            br.slotNeeded = slotTick;
            for (Brick* o : br.nbr) {
                if (!o) continue;
                o->slotNeeded = slotTick;
                for (Brick* o2 : o->nbr) if (o2) o2->slotNeeded = slotTick;
            }
        }
    for (auto& s : sites)
        for (auto& kv : s->bricks) {
            Brick& br = *kv.second;
            if (!br.ready) continue;
            const bool needed = slotTick - br.slotNeeded <= kSlotKeepTicks;
            if (br.slot >= 0 && !needed && !br.active && !br.upload) { // released
                G.owner[br.slot] = nullptr; G.ownerSite[br.slot] = nullptr;
                br.slot = -1; --used; G.nbrDirty = true; ++slotReleases;
                continue;
            }
            if (br.slot < 0) {
                if (!needed) continue;
                GpuGrow(used + 1);
                int slot = 0;
                while (G.owner[slot]) ++slot;
                br.slot = slot; G.owner[slot] = &br; G.ownerSite[slot] = s.get();
                ++used;
                G.nbrDirty = true;
                if (br.hadSlot && !br.upload) { GpuUploadBrick(br); continue; } // back from release: up as it is, asleep
                br.upload = true;
            }
            if (br.upload) { GpuUploadBrick(br); br.upload = false; br.active = true; br.idle = 0; G.nbrDirty = true; }
            br.hadSlot = true;
        }
    slotsUsed = used;
    if (G.nbrDirty) {
        std::vector<int> nb(G.capacity * 6, -1);
        for (int slot = 0; slot < G.capacity; ++slot)
            if (Brick* br = G.owner[slot])
                for (int f = 0; f < 6; ++f) nb[slot * 6 + f] = (br->nbr[f] && br->nbr[f]->ready && br->nbr[f]->slot >= 0) ? br->nbr[f]->slot : -1;
        G.nbrs.upload(nb.data(), G.capacity * 24);
        G.nbrDirty = false;
    }
    GpuTimer_Mark("fluid-upload"); // benches: slots and whole-brick uploads, neighbours (#669)
    // 3. The CPU's deltas (summed per cell).
    if (!gpuDeltas.empty()) {
        struct D { unsigned cell; int delta; unsigned kind; unsigned pad; };
        // Summed per cell by a sort (#669: a hash map of a big pool's levelling nudges took ~10 ms a tick).
        static std::vector<D> raw;
        raw.clear();
        for (const GpuDeltaRec& r : gpuDeltas)
            if (r.b->slot >= 0) raw.push_back({ (unsigned)r.b->slot * kCellsPerBrick + (unsigned)r.i, r.d, r.kind, 0u });
        gpuDeltas.clear();
        std::stable_sort(raw.begin(), raw.end(), [](const D& a, const D& b) { return a.cell < b.cell; }); // later deltas' kinds win
        std::vector<D> list;
        for (size_t i = 0; i < raw.size();) {
            D sum = raw[i];
            size_t j = i + 1;
            for (; j < raw.size() && raw[j].cell == sum.cell; ++j) { sum.delta += raw[j].delta; sum.kind = raw[j].kind; }
            if (sum.delta) list.push_back(sum);
            i = j;
        }
        if (!list.empty()) {
            if (G.deltas.bytes < list.size() * sizeof(D)) G.deltas.resize((unsigned)(list.size() * 2 * sizeof(D)));
            G.deltas.upload(list.data(), (unsigned)(list.size() * sizeof(D)));
            G.delta.use();
            G.delta.setUint("count", (unsigned)list.size());
            G.delta.setInt("useGroups", G.groupsOn ? 1 : 0);
            G.a.bind(0); G.deltas.bind(1); G.flags.bind(2); if (G.groupsOn) G.mask.bind(10); // a change wets its group (#675)
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
        for (int fc = 0; fc < 6; ++fc) {
            Brick* o = br->nbr[fc];
            if (!o || o->slot < 0 || inS[o->slot]) continue;
            // Not run (#669): a brick of solid rock (water never enters it) or a dry brick above (water only rises into it
            // through sources, pours or levelling, which wake it): under and over a big pool they were most of the dispatch.
            if (!o->active && !RunsBeside(*o, fc)) continue;
            inS[o->slot] = 1; list.push_back(o->slot);
        }
    }
    if (list.empty()) return;
    G.inS.upload(inS.data(), G.capacity * 4);
    if (G.list.bytes < list.size() * 4) G.list.resize((unsigned)(list.size() * 8));
    G.list.upload(list.data(), (unsigned)(list.size() * 4));
    GpuTimer_Mark("fluid-list"); // benches: deltas and the brick list (#669)
    // 5. n substeps: vertical A -> B, face flows B -> N, apply B, N -> A, P (#642).
    G.step.use();
    G.step.setUint("count", (unsigned)list.size());
    int spread[8], visc[8];
    float fricLin[8], fricQuad[8];
    for (int k = 0; k < 8; ++k) { spread[k] = kindMinSpread[k]; visc[k] = kindVisc[k]; KindFriction(k, fricLin[k], fricQuad[k]); }
    rlSetUniform(G.step.uniform("minSpread"), spread, RL_SHADER_UNIFORM_INT, 8);
    rlSetUniform(G.step.uniform("visc"), visc, RL_SHADER_UNIFORM_INT, 8);
    rlSetUniform(G.step.uniform("frictionLin"), fricLin, RL_SHADER_UNIFORM_FLOAT, 8);
    rlSetUniform(G.step.uniform("frictionQuad"), fricQuad, RL_SHADER_UNIFORM_FLOAT, 8);
    n *= SubstepScale();
    const float substepDt = kGpuSubstepSeconds / SubstepScale();
    rlSetUniform(G.step.uniform("dt"), &substepDt, RL_SHADER_UNIFORM_FLOAT, 1);
    const float cellSize = WaterField::CellSize();
    rlSetUniform(G.step.uniform("cellSize"), &cellSize, RL_SHADER_UNIFORM_FLOAT, 1);
    G.step.setInt("kColumnCells", ColumnCells());
    static const bool noEvap = EnvFlag("RAYLIBGAME_FLUID_NOEVAP"); // tests: specks kept (water accounting, #454)
    static unsigned evapTick = 0;
    GpuFlowResolve(); // which bricks the flow watches need this dispatch (#392)
    G.nbrs.bind(2); G.list.bind(3); G.inS.bind(4); G.flags.bind(5); G.watch.bind(6); G.flux.bind(7); G.faceP.bind(8); G.faceN.bind(9); G.sinkOut.bind(12);
    const unsigned bricks = (unsigned)list.size(); // 2D dispatch: 64 groups of 64 cells per brick, one row per brick
    // The dry-cell skip (#675, fluid_groups.comp): each substep runs only the 64-cell groups near water, by indirect dispatch.
    G.step.setInt("useGroups", G.groupsOn ? 1 : 0);
    if (G.groupsOn) {
        G.groups.use();
        G.groups.setUint("count", bricks);
        G.mask.bind(10); G.glist.bind(11);
    }
    auto pass = [&](int p) {
        G.step.setInt("pass", p);
        if (G.groupsOn) G.step.dispatchIndirect(G.glist, 0);
        else G.step.dispatch(kCellsPerBrick / 64u, bricks);
    };
    for (int k = 0; k < n; ++k) {
        if (G.groupsOn) {
            G.maskCur = (G.maskCur + 1) % 3;
            G.groups.use();
            G.groups.setInt("cur", G.maskCur);
            G.groups.setInt("mode", 0); G.groups.dispatch(bricks);
            G.groups.setInt("mode", 1); G.groups.dispatch(1);
            G.step.use();
            G.step.setInt("maskCur", G.maskCur);
        }
        G.a.bind(0); G.b.bind(1); pass(0);
        G.b.bind(0); G.a.bind(1); pass(1);
        // Specks soak away on one substep in kEvapEvery: one must rest across it. Every substep, a spreading film lost the
        // edge it had just laid down, again and again (#454: 7-12% of a draining pond).
        G.step.setInt("evaporate", !noEvap && ++evapTick % (kEvapEvery * (unsigned)SubstepScale()) == 0 ? 1 : 0);
        pass(2);
    }
    G.lastList = list;
    G.lastN = n;
    // 6. Async readback (#669): pack the processed bricks' cells, copy them and the flags into the mapped buffer, fence.
    if (G.asyncRb) {
        const size_t packedBytes = (size_t)list.size() * kCellsPerBrick * 2, flagsBytes = (size_t)G.capacity * 4;
        if (G.packed.bytes < packedBytes) G.packed.resize((unsigned)(packedBytes * 2));
        G.pack.use();
        G.pack.setUint("count", (unsigned)list.size());
        G.a.bind(0); G.list.bind(3); G.packed.bind(10);
        G.pack.dispatch(kCellsPerBrick / 2 / 64u, bricks);
        const size_t sinkBytes = G.sinkOn ? (size_t)G.capacity * 4 : 0; // released reaches' sink counters (#656)
        // The drawn current (#677): per column the solver's flow at its surface, 1 KB a brick (fluid_flowpack.comp).
        const size_t flowBytes = G.realFlow ? (size_t)list.size() * 256 * 4 : 0;
        if (flowBytes) {
            if (G.flowPacked.bytes < flowBytes) G.flowPacked.resize((unsigned)(flowBytes * 2));
            G.flowPack.use();
            G.flowPack.setUint("count", (unsigned)list.size());
            const float cs = WaterField::CellSize();
            rlSetUniform(G.flowPack.uniform("cellSize"), &cs, RL_SHADER_UNIFORM_FLOAT, 1);
            G.a.bind(0); G.list.bind(3); G.faceN.bind(9); G.flowPacked.bind(13);
            G.flowPack.dispatch(1, bricks);
        }
        const size_t total = flagsBytes + packedBytes + sinkBytes + flowBytes; // mapped: flags, cells, sink counters, current
        if (G.rbCap < total) {
            if (G.rbBuf) { glDeleteBuffers(1, (GLuint*)&G.rbBuf); G.rbBuf = 0; G.rbPtr = nullptr; }
            size_t cap = 1u << 20;
            while (cap < total) cap <<= 1;
            GLuint buf = 0;
            glGenBuffers(1, &buf);
            glBindBuffer(GL_COPY_WRITE_BUFFER, buf);
            const GLbitfield mapFlags = GL_MAP_READ_BIT | GL_MAP_PERSISTENT_BIT | GL_MAP_COHERENT_BIT;
            GpuMemTag gm_("fluid readback"); // VRAM accounting (#302)
            glBufferStorage(GL_COPY_WRITE_BUFFER, (GLsizeiptr)cap, nullptr, mapFlags | GL_CLIENT_STORAGE_BIT);
            G.rbPtr = (unsigned char*)glMapBufferRange(GL_COPY_WRITE_BUFFER, 0, (GLsizeiptr)cap, mapFlags);
            glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
            if (!G.rbPtr) { glDeleteBuffers(1, &buf); G.asyncRb = false; TraceLog(LOG_WARNING, "Fluid: persistent map failed: synchronous readback"); return; }
            G.rbBuf = buf; G.rbCap = cap;
        }
        glBindBuffer(GL_COPY_WRITE_BUFFER, G.rbBuf);
        glBindBuffer(GL_COPY_READ_BUFFER, G.flags.id);
        glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, 0, (GLsizeiptr)flagsBytes);
        glBindBuffer(GL_COPY_READ_BUFFER, G.packed.id);
        glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, (GLintptr)flagsBytes, (GLsizeiptr)packedBytes);
        if (sinkBytes) {
            glBindBuffer(GL_COPY_READ_BUFFER, G.sinkOut.id);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, (GLintptr)(flagsBytes + packedBytes), (GLsizeiptr)sinkBytes);
        }
        if (flowBytes) {
            glBindBuffer(GL_COPY_READ_BUFFER, G.flowPacked.id);
            glCopyBufferSubData(GL_COPY_READ_BUFFER, GL_COPY_WRITE_BUFFER, 0, (GLintptr)(flagsBytes + packedBytes + sinkBytes), (GLsizeiptr)flowBytes);
        }
        glBindBuffer(GL_COPY_READ_BUFFER, 0);
        glBindBuffer(GL_COPY_WRITE_BUFFER, 0);
        if (sinkBytes) { static std::vector<unsigned> zeros; zeros.assign(G.capacity, 0u); G.sinkOut.upload(zeros.data(), (unsigned)sinkBytes); } // counted from here on
        G.rbFlagsBytes = flagsBytes;
        G.rbSinkAt = flagsBytes + packedBytes; G.rbSinkBytes = sinkBytes;
        G.rbFlowAt = flagsBytes + packedBytes + sinkBytes; G.rbFlowBytes = flowBytes;
        G.rbFence = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
        glFlush();
    }
}

// CPU mode (#654: no compute shaders - a headless server, FLUID_CPU=1): the same step as the GPU's (fluidkernel.h) on a grid
// packed each tick from the awake bricks, their neighbours (processed too, as on the GPU) and the ring past them (read only,
// for the columns and the cells beside), then unpacked with the GPU readback's rules for waking and sleeping.
FluidKernel::Grid cpuGrid;
void CpuTick(int n) {
    // TEMP #747 (DBSTEPPROF): where a CPU tick's time goes - pack, the kernel's phases, unpack - every 40 ticks.
    static const bool stepProf = EnvFlag("DBSTEPPROF");
    static double profPack = 0.0, profKernel = 0.0, profUnpack = 0.0; static int profTicks = 0; static long long profSlots = 0;
    using ProfClock = std::chrono::steady_clock;
    const auto profT0 = ProfClock::now();
    std::vector<Brick*> slots;
    std::vector<Site*> slotSite;
    auto add = [&](Brick* b, Site* s) {
        if (!b || !b->ready || b->cpuSlot >= 0) return;
        b->cpuSlot = (int)slots.size();
        slots.push_back(b); slotSite.push_back(s);
    };
    std::vector<int> list;
    for (auto& s : sites)
        for (auto& kv : s->bricks) {
            Brick* b = kv.second.get();
            if (!b->ready || !b->active) continue;
            add(b, s.get());
            for (Brick* o : b->nbr) add(o, s.get());
        }
    const int processedCount = (int)slots.size();
    for (int i = 0; i < processedCount; ++i) list.push_back(i);
    for (int i = 0; i < processedCount; ++i) for (Brick* o : slots[i]->nbr) add(o, slotSite[i]); // read-only ring
    if (list.empty()) { for (Brick* b : slots) b->cpuSlot = -1; return; }
    FluidKernel::Grid& g = cpuGrid;
    const int S = (int)slots.size();
    if ((int)g.processed.size() < S) g.resize(std::max(S, (int)g.processed.size() * 2));
    std::fill(g.flags.begin(), g.flags.begin() + S, 0u);
    std::fill(g.sinkOut.begin(), g.sinkOut.begin() + S, 0u);
    for (int i = 0; i < S; ++i) {
        Brick* b = slots[i];
        uint32_t* c = &g.a[(size_t)i * B3];
        for (int k = 0; k < B3; ++k) c[k] = b->fill[k] | ((uint32_t)b->flags[k] << 8) | ((uint32_t)b->target[k] << 16) | ((uint32_t)b->shut[k] << 24);
        std::copy(c, c + B3, &g.b[(size_t)i * B3]);
        if (b->faces.size() == (size_t)B3) std::copy(b->faces.begin(), b->faces.end(), &g.faceP[(size_t)i * B3]);
        else std::fill(&g.faceP[(size_t)i * B3], &g.faceP[(size_t)i * B3] + B3, FluidKernel::Vec4{});
        for (int f = 0; f < 6; ++f) g.nbr[(size_t)i * 6 + f] = b->nbr[f] && b->nbr[f]->cpuSlot >= 0 ? b->nbr[f]->cpuSlot : -1;
        g.processed[i] = i < processedCount ? 1 : 0;
        // The skip's history (#747), as the masks of substeps "0" (the last) and "2" (the one before): BeginTick adds the
        // groups wet now, the first substep (gcur 0 -> 1) reads 0 and 2.
        uint32_t* m = &g.gmask[(size_t)i * 6];
        m[0] = b->cpuMask[0]; m[1] = b->cpuMask[1]; m[2] = 0u; m[3] = 0u; m[4] = b->cpuMask[2]; m[5] = b->cpuMask[3];
    }
    FluidKernel::Params p;
    for (int k = 0; k < 8; ++k) { p.minSpread[k] = kindMinSpread[k]; p.visc[k] = kindVisc[k]; KindFriction(k, p.frictionLin[k], p.frictionQuad[k]); }
    n *= SubstepScale();
    p.dt = kGpuSubstepSeconds / SubstepScale();
    p.setCell(WaterField::CellSize(), ColumnCells());
    static const bool noEvap = EnvFlag("RAYLIBGAME_FLUID_NOEVAP");
    static unsigned evapTick = 0;
    // #680 diagnosis (DBVOL): water in the processed slots' plain cells before and after (sources, sinks and evaporation aside).
    static const bool dbVol = EnvFlag("DBVOL");
    auto plainSum = [&]() { long long v = 0; for (int i = 0; i < processedCount; ++i) for (int c = 0; c < B3; ++c) { const uint32_t u = g.a[(size_t)i * B3 + c]; if (!((u >> 8) & (kSource | kSink | kSolid))) v += u & 255u; } return v; };
    const long long volBefore = dbVol ? plainSum() : 0;
    // The dry-cell skip (#747): slots are new this tick. FLUID_CPU_FULLFIRST=1: its first substep runs every cell (the old
    // way, A/B); else the kernel clears the slots' scratch and seeds the wet groups from the cells, and every substep skips.
    static const bool fullFirst = EnvFlag("FLUID_CPU_FULLFIRST");
    if (fullFirst) g.groupsFull = true;
    else FluidKernel::BeginTick(g, list, S);
    const auto profT1 = ProfClock::now();
    static const bool dbTrace = EnvFlag("DBTRACE");
    // TEMP #680 (DBTRACE): transfers into partial / solid-capped cells near the test's wall, in the channel frame (0.25 m
    // buckets): (taker along, taker height, giver along, branch * 2 + capped) -> units.
    static std::map<std::array<int, 5>, long long> traceSum;
    static std::map<std::array<int, 3>, std::array<long long, 2>> pairSum;
    auto chan = [&](int cell, float& along, float& across, float& h) {
        const Brick* b = slots[cell / B3];
        const int l = cell % B3;
        const Vector3 d = Vector3Subtract(slotSite[cell / B3]->cellCentre(b->bi * B + l % B, b->bj * B + (l / B) % B, b->bk * B + l / (B * B)), dbgAxO);
        along = Vector3DotProduct(d, dbgAxX); across = Vector3DotProduct(d, dbgAxY); h = Vector3DotProduct(d, dbgAxU);
    };
    if (dbTrace && dbgAxSet) { // the slots whose brick lies within ~a brick of the wall (x 0-2 m)
        g.dbgMark.assign(S, 0);
        for (int i = 0; i < S; ++i) {
            const Vector3 d = Vector3Subtract(slotSite[i]->cellCentre(slots[i]->bi * B + B / 2, slots[i]->bj * B + B / 2, slots[i]->bk * B + B / 2), dbgAxO);
            const float a = Vector3DotProduct(d, dbgAxX), c = Vector3DotProduct(d, dbgAxY);
            g.dbgMark[i] = a > -3.0f && a < 12.0f && std::fabs(c) < 6.0f;
        }
    }
    for (int k = 0; k < n; ++k) {
        p.evaporate = !noEvap && ++evapTick % (kEvapEvery * (unsigned)SubstepScale()) == 0;
        FluidKernel::Substep(g, p, list);
        if (dbTrace && dbgAxSet) {
            const int m = std::min(FluidKernel::dbgTraceN.exchange(0), FluidKernel::kDbgTrace);
            for (int i = 0; i < m; ++i) {
                const FluidKernel::DbgTransfer& t = FluidKernel::dbgTrace[i];
                float ta, tc, th, ga, gc, gh;
                chan(t.taker, ta, tc, th); chan(t.giver, ga, gc, gh);
                if (ta < -0.5f || ta > 0.25f || std::fabs(tc) > 2.5f || th < -1.0f) continue; // at the wall's face, over the pool's surface
                traceSum[{ (int)std::floor(ta * 4.0f), (int)std::floor(th * 4.0f), (int)std::floor(ga * 4.0f), (int)std::floor(gh * 4.0f), t.branch * 2 + t.capped }] += t.units;
            }
            // The pairs runs at x 3-8 m planned: (along 0.5 m, direction 0 down / 1 up / 2 across, exit code) -> count, depth.
            const int mp = std::min(FluidKernel::dbgPairsN.exchange(0), FluidKernel::kDbgTrace);
            for (int i = 0; i < mp; ++i) {
                const FluidKernel::DbgPairRec& r = FluidKernel::dbgPairs[i];
                float aa, ac, ah;
                chan(r.a, aa, ac, ah);
                if (aa < 0.0f || aa > 4.0f || std::fabs(ac) > 2.5f) continue;
                const Site* st = slotSite[r.a / B3];
                const float da = Vector3DotProduct(Vector3Add(Vector3Scale(st->e, (float)r.dx), Vector3Scale(st->n, (float)r.dy)), dbgAxX);
                const int dir = da > 0.3f ? 0 : da < -0.3f ? 1 : 2;
                auto& e = pairSum[{ (int)std::floor(aa * 2.0f), dir, r.why }];
                ++e[0]; e[1] += r.depth;
            }
        } else if (dbTrace) { FluidKernel::dbgTraceN = 0; FluidKernel::dbgPairsN = 0; }
    }
    if (dbTrace) {
        static int pticks = 0;
        if (++pticks % 10 == 0 && !pairSum.empty()) {
            printf("dbpairs: runs at x 3-8 m (along, d/u/x, exit: count, mean depth m):");
            for (auto& kv : pairSum)
                printf(" %.1f%c%d:%lld/%.2f", kv.first[0] / 2.0f, "dux"[kv.first[1]], kv.first[2], kv.second[0], kv.second[1] / (double)kv.second[0] / (255.0 / WaterField::CellSize()));
            printf("\n");
            pairSum.clear();
        }
    }
    if (dbTrace) {
        static int ticks = 0;
        if (++ticks % 40 == 0 && !traceSum.empty()) {
            printf("dbtrace: into partial/capped cells by the wall (taker along m, h m <- giver along m: branch s/w/f, c = under solid, units):\n");
            for (auto& kv : traceSum)
                printf("  %.2f h%.2f <- %.2f h%.2f %c%s %lld\n", kv.first[0] / 4.0f, kv.first[1] / 4.0f, kv.first[2] / 4.0f, kv.first[3] / 4.0f, "swf"[kv.first[4] / 2], kv.first[4] % 2 ? "c" : "", kv.second);
            traceSum.clear();
        }
    }
    if (dbVol) {
        static long long drift = 0; static int ticks = 0;
        drift += plainSum() - volBefore;
        if (++ticks % 40 == 0) { printf("dbvol: plain-cell water changed by %+.2f m^3 over the last second (%d slots)\n", drift * WaterField::CellVolume() / 255.0, processedCount); drift = 0; }
    }
    const auto profT2 = ProfClock::now();
    // Unpack, as the GPU readback (GpuTick 1.).
    static const int d[6][3] = { { -1, 0, 0 }, { 1, 0, 0 }, { 0, -1, 0 }, { 0, 1, 0 }, { 0, 0, -1 }, { 0, 0, 1 } };
    for (int i = 0; i < processedCount; ++i) {
        Brick* br = slots[i];
        const uint32_t fl = g.flags[i];
        br->faces.assign(&g.faceP[(size_t)i * B3], &g.faceP[(size_t)i * B3] + B3);
        { // the skip's history for the next tick (#747): this tick's last two substeps' wet groups
            const uint32_t* m = &g.gmask[(size_t)i * 6];
            const int c0 = g.gcur, c1 = (g.gcur + 2) % 3;
            br->cpuMask[0] = m[c0 * 2]; br->cpuMask[1] = m[c0 * 2 + 1]; br->cpuMask[2] = m[c1 * 2]; br->cpuMask[3] = m[c1 * 2 + 1];
        }
        if (g.sinkOut[i] && br->sinkRelease >= 0 && br->sinkRelease < (int)slotSite[i]->releases.size()) slotSite[i]->releases[br->sinkRelease].sinkUnits += g.sinkOut[i]; // #656
        if (fl & 1u) {
            const uint32_t* c = &g.a[(size_t)i * B3];
            int drift = 0;
            for (int k = 0; k < B3; ++k) {
                br->fill[k] = (uint8_t)(c[k] & 255u);
                br->flags[k] = (uint8_t)((br->flags[k] & ~kKindMask) | (((c[k] >> 8) & 255u) & kKindMask));
                if (br->settled) drift = std::max(drift, abs((int)br->fill[k] - (int)br->check[k]));
            }
            if (br->settled && drift <= kSettleDiff) br->changed = true; // a resting brick at the edge: stays asleep (#366)
            else {
                if (br->settled) { ++wakeGpu; br->settled = false; }
                Touch(br);
                for (int fc = 0; fc < 6; ++fc) if (Brick* o = br->nbr[fc]; o && RunsBeside(*o, fc)) { o->active = true; o->idle = 0; }
            }
        } else if (br->active && (br->idle += n) > kSleepSteps) {
            br->active = false;
        }
        // A face request: the neighbour brick is made if missing, woken if asleep (#680: the column pairs ask for sleeping ones).
        for (int face = 0; face < 6; ++face)
            if (fl & (2u << face)) {
                if (Brick* o = br->nbr[face]) { o->active = true; o->idle = 0; }
                else slotSite[i]->ensure(br->bi + d[face][0], br->bj + d[face][1], br->bk + d[face][2]);
            }
    }
    for (Brick* b : slots) b->cpuSlot = -1;
    if (stepProf) {
        const auto profT3 = ProfClock::now();
        auto ms = [](ProfClock::time_point a, ProfClock::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); };
        profPack += ms(profT0, profT1); profKernel += ms(profT1, profT2); profUnpack += ms(profT2, profT3);
        profSlots += (long long)slots.size();
        if (++profTicks == 40) {
            const double k = 1.0 / profTicks;
            namespace FK = FluidKernel;
            printf("dbstep: per tick (40): pack %.1f ms, kernel %.1f ms [select %.1f, vertical %.1f, plan/faces %.1f, move %.1f, apply %.1f], unpack %.1f ms; "
                   "%.0f slots, %.1f substeps (%.2f full), groups run %.1f%%\n",
                   profPack * k, profKernel * k, FK::dbgPhaseMs[0] * k, FK::dbgPhaseMs[1] * k, FK::dbgPhaseMs[2] * k, FK::dbgPhaseMs[3] * k, FK::dbgPhaseMs[4] * k,
                   profUnpack * k, profSlots * k, FK::dbgSubsteps * k, FK::dbgFullSubsteps * k, 100.0 * FK::dbgGroupsRun / std::max(1LL, FK::dbgGroupsAll));
            profPack = profKernel = profUnpack = 0.0; profTicks = 0; profSlots = 0;
            for (double& v : FK::dbgPhaseMs) v = 0.0;
            FK::dbgGroupsRun = FK::dbgGroupsAll = FK::dbgSubsteps = FK::dbgFullSubsteps = 0;
        }
    }
}

void GpuReset() {
    if (G.capacity > 0) { std::vector<unsigned> zeros(G.capacity, 0u); G.flags.upload(zeros.data(), G.capacity * 4); } // no stale face requests on reused slots
    for (Brick*& o : G.owner) o = nullptr;
    for (Site*& o : G.ownerSite) o = nullptr;
    G.lastList.clear();
    if (G.rbFence) { glDeleteSync((GLsync)G.rbFence); G.rbFence = nullptr; } // the last dispatch's readback is dropped (#669)
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
    partSnaps.clear(); // build parts are read again from the new world (#391)
    sites.clear();
    refs.clear();
    stepCount = 0;
    { std::lock_guard<std::mutex> c(cmdMx); commands.clear(); }
    {
        std::lock_guard<std::mutex> p(pubMx);
        pubStats = Stats{};
    }
    WaterField::Clear(); // only generated water until the simulation publishes again
    WaterSurface::Clear(); // the old world's water surface map (#639): a new field may reuse its address
    accum = 0.0;
}

void Fluid::OnTerrainEdit(Planetoid& body, Vector3 localCentre, float radius, bool adds) {
    if (!dynamic_cast<const PlanetField*>(body.field.get())) return;
    std::lock_guard<std::mutex> c(cmdMx);
    commands.push_back({ 0, body.id, localCentre, radius, adds ? 1.0f : 0.0f }); // m3 > 0: the edit adds ground
}

Vector3 Fluid::FlowAt(const Planetoid& body, Vector3 local, float radius) {
    std::lock_guard<std::mutex> lk(flowMx);
    const double now = GetTime();
    radius = std::max(0.5f, radius); // as stored: a smaller radius matched no watch and added one every call (#392 review)
    for (FlowWatch& w : flowWatches)
        if (w.body == body.id && fabsf(w.r - radius) < 0.01f && Vector3DistanceSqr(w.local, local) < 0.25f) { w.lastAsk = now; return w.value; }
    FlowWatch w; w.body = body.id; w.local = local; w.r = radius; w.lastAsk = now;
    flowWatches.push_back(w);
    return { 0, 0, 0 };
}

int Fluid::ReleasedReaches() {
    std::lock_guard<std::mutex> lock(simMx);
    int n = 0;
    for (auto& s : sites) n += (int)s->releases.size();
    return n;
}

float Fluid::ReleaseInflowScale() {
    std::lock_guard<std::mutex> lock(simMx);
    float lo = -1.0f;
    for (auto& s : sites)
        for (const Site::Spring& sp : s->springs)
            if (sp.rel >= 0) lo = lo < 0.0f ? sp.scale : std::min(lo, sp.scale);
    return lo;
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

int Fluid::DbgCellShut(const Planetoid& body, Vector3 local) {
    Site* s = SiteFor(body.id, local, false);
    if (!s) return -1;
    int ci, cj, ck; s->cellOf(local, ci, cj, ck);
    Brick* b; int i;
    if (!s->cell(ci, cj, ck, b, i) || !b->ready) return -1;
    return (b->flags[i] & kSolid) ? 255 : (int)b->shut[i];
}
void Fluid::DbgTraceAxes(Vector3 o, Vector3 ax, Vector3 ay, Vector3 up) { dbgAxO = o; dbgAxX = ax; dbgAxY = ay; dbgAxU = up; dbgAxSet = true; }
float Fluid::PourDropped() { return (float)(pourDroppedUnits.load() * (double)WaterField::CellVolume() / 255.0); }

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
        if (k.viscosity < 4 && i < (int)kinds.size()) LOG_WARN("Fluid", "%s: viscosity %d below 4, using 4\n", k.id.c_str(), k.viscosity);
        kindVisc[i] = std::max(k.viscosity, 4); kindMinSpread[i] = k.minSpread;
    }
}
int Fluid::KindCount() { return (int)kinds.size(); }
const Fluid::Kind& Fluid::KindInfo(int kind) { return kinds[std::clamp(kind, 0, (int)kinds.size() - 1)]; }
int Fluid::KindByName(const std::string& id) {
    for (size_t i = 0; i < kinds.size(); ++i) if (kinds[i].id == id) return (int)i;
    return -1;
}

// Build parts (#391), main thread: each body's plots are copied into an immutable PartSolid when a plot's version
// changes; the cells that became solid or open are sent as terrain-edit commands (2 m buckets), so the bricks there
// are reclassified, a site starts where parts go into water, and solid added into a river releases the reach (#361).
static void UpdatePartSnaps(Scene& world) {
    for (auto& o : world.objects) {
        Planetoid* p = dynamic_cast<Planetoid*>(o.get());
        if (!p || !p->plots || !dynamic_cast<const PlanetField*>(p->field.get())) continue;
        const auto old = partSnaps.find(p->id);
        bool same = old != partSnaps.end() && old->second->plots.size() == p->plots->plots.size();
        if (same)
            for (size_t i = 0; i < p->plots->plots.size() && same; ++i)
                same = old->second->plots[i].id == p->plots->plots[i]->id && old->second->plots[i].version == p->plots->plots[i]->version;
        if (same) continue;
        auto snap = std::make_shared<PartSolid>();
        for (const auto& pl : p->plots->plots) {
            PlotSnap ps;
            ps.id = pl->id; ps.version = pl->version;
            ps.origin = pl->origin.toFloat(); ps.east = pl->east; ps.up = pl->up; ps.north = pl->north;
            Vector3 lo{ 1e30f, 1e30f, 1e30f }, hi{ -1e30f, -1e30f, -1e30f };
            for (const auto& [cell, ref] : pl->cells) {
                bool in = ref.shape.full;
                if (!in) { // a sloped or cut shape: solid where it holds the cell's centre
                    in = true;
                    for (int k = 0; k < ref.shape.planes; ++k) in = in && Vector3DotProduct(ref.shape.n[k], { 0.5f, 0.5f, 0.5f }) <= ref.shape.d[k];
                }
                if (!in) continue;
                ps.solid.insert(cell);
                const Vector3 c = pl->toBody({ (cell.x + 0.5f) * kCellSize, (cell.y + 0.5f) * kCellSize, (cell.z + 0.5f) * kCellSize });
                lo = Vector3Min(lo, c); hi = Vector3Max(hi, c);
            }
            // Every plot keeps an entry, an empty one too: the version check above compares them one to one (an
            // emptied plot left out made the counts differ for good, and the whole snapshot was rebuilt every tick).
            if (ps.solid.empty()) { ps.centre = ps.origin; ps.radius = 0.0f; } // solidAt skips it
            else {
                ps.centre = Vector3Scale(Vector3Add(lo, hi), 0.5f);
                ps.radius = Vector3Distance(lo, hi) * 0.5f + kCellSize;
            }
            snap->plots.push_back(std::move(ps));
        }
        // What changed: solid cells added or removed, per plot, as edit commands bucketed to 2 m.
        std::unordered_map<int64_t, std::pair<Vector3, bool>> buckets;
        auto note = [&](const PlotSnap& ps, const Int3& cell, bool added) {
            const Vector3 c = Vector3Add(ps.origin, Vector3Add(Vector3Add(Vector3Scale(ps.east, (cell.x + 0.5f) * kCellSize), Vector3Scale(ps.up, (cell.y + 0.5f) * kCellSize)), Vector3Scale(ps.north, (cell.z + 0.5f) * kCellSize)));
            const int64_t k = Key((int)floorf(c.x / 2.0f), (int)floorf(c.y / 2.0f), (int)floorf(c.z / 2.0f));
            auto it = buckets.find(k);
            if (it == buckets.end()) buckets[k] = { c, added }; else it->second.second = it->second.second || added;
        };
        const PartSolid* prev = old != partSnaps.end() ? old->second.get() : nullptr;
        for (const PlotSnap& ps : snap->plots) {
            const PlotSnap* was = nullptr;
            if (prev) for (const PlotSnap& q : prev->plots) if (q.id == ps.id) was = &q;
            for (const Int3& c : ps.solid) if (!was || !was->solid.count(c)) note(ps, c, true);
            if (was) for (const Int3& c : was->solid) if (!ps.solid.count(c)) note(ps, c, false);
        }
        if (prev)
            for (const PlotSnap& q : prev->plots) {
                bool kept = false;
                for (const PlotSnap& ps : snap->plots) kept = kept || ps.id == q.id;
                if (!kept) for (const Int3& c : q.solid) note(q, c, false); // a whole plot gone
            }
        partSnaps[p->id] = snap;
        if (old == partSnaps.end() && buckets.size() > 4096) continue; // a load: nothing moved yet, sites classify with it
        std::lock_guard<std::mutex> c(cmdMx);
        for (const auto& [k, v] : buckets) commands.push_back({ 0, p->id, v.first, 1.6f, v.second ? 1.0f : 0.0f });
    }
}

// River flow (#656): once a second, each released reach's outflow (what its downstream sink took in) against the inflow its
// springs give, as a share of the river's flow below it; then the network carries the changes down and the water along the
// links whose flow changed is rebuilt (surface map tiles and chunk meshes). Under simMx with no job running.
const RiverNetwork* RiversOf(uint64_t bodyId) {
    Planetoid* body = Planetoid::find(bodyId);
    const PlanetField* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
    const RiverNetwork* net = pf ? pf->rivers() : nullptr;
    return net && net->bakeNative() ? net : nullptr;
}

void UpdateRiverFlow() {
    // The links whose water changed, refreshed one a tick (40 a second): a change running down a river remeshed every link
    // it reached in the same frame, a 30-40 ms hitch once a second while it ran.
    static std::deque<std::pair<uint64_t, int>> refreshQueue;
    static std::set<std::pair<uint64_t, int>> refreshQueued;
    if (!refreshQueue.empty()) {
        const auto [bodyId, piece] = refreshQueue.front();
        refreshQueue.pop_front(); refreshQueued.erase({ bodyId, piece });
        Planetoid* body = Planetoid::find(bodyId);
        const PlanetField* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
        if (const RiverNetwork* net = pf ? pf->rivers() : nullptr; net && net->pieceBounds(piece).radius > 0.0f) {
            const RiverNetwork::PieceBounds b = net->pieceBounds(piece);
            WaterSurface::Invalidate(*pf, Vector3Normalize(b.mid), b.radius);
            RequestRemesh(*body, b.mid, b.radius);
        }
    }
    static double last = -1.0;
    static long long lastStep = 0;
    const double now = GetTime();
    if (last < 0.0) { last = now; lastStep = stepCount; return; }
    const float dt = (float)(now - last);
    if (dt < 1.0f) return;
    last = now;
    // The outflow per second of SIMULATED time (#656 follow-up): the springs give natQ per simulated second (kSubsteps steps a
    // 50 ms), but a busy sim drops ticks (a refilling river kept ~1800 bricks awake), so per wall second the same flow read low.
    const long long steps = stepCount - lastStep;
    lastStep = stepCount;
    static const bool wallTime = EnvFlag("RAYLIBGAME_RIVERFLOW_WALL"); // A/B: the old per-wall-second measure
    const float simDt = wallTime ? dt : (float)steps * 0.05f / (float)kSubsteps;
    static const bool flowLog = EnvFlag("RAYLIBGAME_RIVERFLOW_LOG"); // diagnostics: each reach's measure once a second
    static std::vector<int> changed;
    for (auto& s : sites) {
        Planetoid* body = Planetoid::find(s->body);
        const PlanetField* pf = body ? dynamic_cast<const PlanetField*>(body->field.get()) : nullptr;
        const RiverNetwork* net = pf ? pf->rivers() : nullptr;
        if (!net || !net->bakeNative()) continue;
        bool anyPending = !s->pending.empty();
        for (auto& kv : s->bricks) anyPending = anyPending || !kv.second->ready;
        for (size_t ri = 0; ri < s->releases.size(); ++ri) {
            Site::Release& r = s->releases[ri];
            if (r.piece < 0 || r.natQ <= 0.0f) continue;
            const float q = simDt > 0.05f ? (float)(r.sinkUnits * WaterField::CellVolume() / 255.0 / simDt) : 0.0f; // m^3/s (no steps: nothing moved)
            const double units = r.sinkUnits;
            r.sinkUnits = 0.0;
            // Bricks waiting for (re)classification (a load, a dig) hold no ground counts and their sinks measure nothing yet:
            // keep the last state until the site is classified (a load then read every reach open for a second). The same
            // after a load reset the step count (no simulated time to measure against).
            if (anyPending || steps < 0) continue;
            // The dam gone (#656 follow-up): (almost) no ground added into the reach's river any more - the river runs as
            // generated again at once. The measured outflow lags the opening by minutes: the pool must refill the reach's
            // drained channel first, and the cellular step carries little on a gentle slope (#680), so the river below stayed
            // low long after both dams were dug out. Open: at most 15% of the most ground the dam put in the river (seed 42's
            // dug-out dam left 3 of 31 cells) or a couple of cells; a partial dig-out leaves more: measured, as before.
            int blocked = 0;
            for (auto& kv : s->bricks) if (kv.second->blockRelease == (int)ri) blocked += kv.second->blocked;
            r.peakBlocked = std::max(r.peakBlocked, blocked);
            const bool open = (float)blocked * WaterField::CellVolume() <= std::max(2.0f * WaterField::CellVolume(), 0.5f) || blocked <= r.peakBlocked * 15 / 100;
            if (open) r.share = 1.0f;
            else r.share += (std::clamp(q / r.natQ, 0.0f, 1.5f) - r.share) * 0.3f; // smoothed: a few seconds of trend
            // Most of the inflow leaving: the river as generated (the springs' spread and the reach's side edges lose some).
            if (r.share > 0.6f) net->clearReachFlow(r.piece);
            else net->setReachFlow(r.piece, r.share / 0.6f);
            if (flowLog) {
                int nb = 0, act = 0, std_ = 0, slept = 0, siteAct = 0, siteSlot = 0, notReady = 0;
                float vol = 0.0f;
                for (auto& kv : s->bricks) {
                    const Brick& b = *kv.second;
                    siteAct += b.active ? 1 : 0; siteSlot += b.slot >= 0 ? 1 : 0; notReady += b.ready ? 0 : 1; vol += b.statVolume;
                    if (b.sinkRelease != (int)ri) continue;
                    ++nb; act += b.active ? 1 : 0; std_ += b.steady ? 1 : 0; slept += b.settled ? 1 : 0;
                }
                int sp = 0, want = 0; float sc = 1.0f;
                for (const Site::Spring& g : s->springs)
                    if (g.rel == (int)ri) { ++sp; want += SpringWants(*s, g) ? 1 : 0; sc = g.scale; }
                printf("riverflow: t=%.0f rel %zu: q %.2f (wall %.2f, sim/wall %.2f) natQ %.2f share %.2f blocked %d flow(piece %d) %.2f flow(up %d) %.2f; sink bricks %d (active %d steady %d asleep %d); springs %d want %d scale %.2f; site %d: bricks %zu active %d slots %d unready %d water %.0f m3\n",
                       now, ri, q, (float)(units * WaterField::CellVolume() / 255.0 / dt), simDt / dt, r.natQ, r.share, blocked, r.piece, net->flowOf(r.piece), r.upPiece, net->flowOf(r.upPiece), nb, act, std_, slept, sp, want, sc, s->index, s->bricks.size(), siteAct, siteSlot, notReady, vol);
            }
        }
        // A reach's springs bring what arrives from upstream: a cut further up the river feeds it less (#656).
        for (Site::Spring& sp : s->springs)
            if (sp.rel >= 0 && sp.rel < (int)s->releases.size()) {
                const int up = s->releases[sp.rel].upPiece;
                sp.scale = up >= 0 ? std::clamp(net->flowOf(up), 0.0f, 1.5f) : 1.0f;
            }
        changed.clear();
        if (!net->stepFlow(std::min(dt, 2.0f), changed)) continue; // a long gap (a load, a stall) must not jump the river: it eases at its speed
        for (int piece : changed)
            if (refreshQueued.insert({ s->body, piece }).second) refreshQueue.push_back({ s->body, piece });
    }
}

void Fluid::Update(Scene& world, float dt) {
    // A step every 50 ms when the last one has finished (the job owns the state while it runs).
    accum += dt;
    // RAYLIBGAME_FLUID_LOCKSTEP=1 (tests): wait for a running step instead of skipping the tick, so the water keeps the game's
    // clock however slow the step is (#747: the CPU kernel at 0.25 m ran at ~25% of real time and fronts looked stalled).
    static const bool lockstep = EnvFlag("RAYLIBGAME_FLUID_LOCKSTEP");
    if (accum < kTickSeconds) return;
    if (running.load()) {
        if (!lockstep) return;
        const double w0 = GetTime();
        WaitIdle();
        // A progress line every 40 waits (~1 s of water at the CPU's pace): how much of the run went to waiting for it.
        static double waited = 0.0; static long long lsTicks = 0;
        waited += GetTime() - w0;
        if (++lsTicks % 40 == 0) { printf("fluid lockstep: %lld waits, %.1f s waited in all\n", lsTicks, waited); fflush(stdout); }
    }
    accum = std::min(accum - kTickSeconds, 0.1);
    UpdatePartSnaps(world); // build parts placed or removed: refresh the water there (#391)
    bool work;
    { std::lock_guard<std::mutex> c(cmdMx); work = !commands.empty(); }
    {
        std::lock_guard<std::mutex> lock(simMx); // free: no job is running
        static bool decided = false;
        if (!decided) { decided = true; gpuMode = !EnvFlag("FLUID_CPU") && GpuInit(); }
        if (gpuMode && !G.lastList.empty()) work = true; // the last dispatch's results still to read back
        for (auto& s : sites) {
            if (!s->pending.empty()) work = true;
            for (auto& kv : s->bricks) if (kv.second->active || kv.second->changed) { work = true; break; }
            for (const Site::Spring& sp : s->springs) if (!work && SpringWants(*s, sp)) work = true; // #645: inflow continues while everything sleeps
            if (work) break;
        }
        UpdateRiverFlow(); // #656
        if (!work) return;
        // Snapshot the bodies the job may touch (fields are immutable; edits are copy-on-write).
        std::unordered_set<uint64_t> ids;
        for (auto& s : sites) ids.insert(s->body);
        { std::lock_guard<std::mutex> c(cmdMx); for (auto& cm : commands) ids.insert(cm.body); }
        refs.clear();
        for (uint64_t id : ids)
            if (Planetoid* p = Planetoid::find(id)) {
                auto ps = partSnaps.find(id);
                refs[id] = { p->field, p->edits, dynamic_cast<const PlanetField*>(p->field.get()), ps != partSnaps.end() ? ps->second : nullptr };
            }
    }
    // GPU substeps: twice the CPU's - the order-independent GPU step moves water one cell a substep, where the
    // CPU's in-place sweep carried it several cells along the sweep.
    if (gpuMode) { std::lock_guard<std::mutex> lock(simMx); GpuTick((int)lround(2 * kSubsteps * (kTickSeconds / 0.05)) * std::max(1, speed.load())); } // 120 a second, as before
    (void)world;
    running = true;
    JobSystem::tagNext(JobTag::Fluid); WorldJobs().submit(-0.5f, [] { RunJob(); running = false; });
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
    // The debug view (#695): RAYLIBGAME_FLUID_DEBUG=1, or toggled in game (Fluid::SetDebugView).
    static int dbgLoc = -1; static unsigned dbgId = 0;
    if (shader.id != dbgId) { dbgId = shader.id; dbgLoc = GetShaderLocation(shader, "fluidDebug"); }
    const int dbg = DebugView() ? 1 : 0;
    if (dbgLoc >= 0) rlSetUniform(dbgLoc, &dbg, RL_SHADER_UNIFORM_INT, 1);
}

static bool debugView = getenv("RAYLIBGAME_FLUID_DEBUG") && atoi(getenv("RAYLIBGAME_FLUID_DEBUG")) != 0;
bool Fluid::DebugView() { return debugView; }
void Fluid::SetDebugView(bool on) { debugView = on; }

namespace {
    // A published cell at a body-local point: its fill (0-255) and flags (and its ground share, #454); false if none.
    bool PubCell(const WaterField::Site& s, Vector3 local, int dk, uint8_t& fill, uint8_t& flags, float* frac = nullptr, uint8_t* shut = nullptr) {
        float fx, fy, fz;
        s.frame(local, fx, fy, fz);
        uint8_t sh = 0;
        if (!s.cell((int)floorf(fx), (int)floorf(fy), (int)floorf(fz) + dk, fill, flags, sh)) return false;
        if (frac) *frac = fz - floorf(fz);
        if (shut) *shut = sh;
        return true;
    }
}

float Fluid::PressureAt(const Planetoid& body, Vector3 local) {
    constexpr float kRhoG = 1000.0f * 9.81f; // Pa per metre of water
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        float fx, fy, fz;
        s.frame(local, fx, fy, fz);
        const int ci = (int)floorf(fx), cj = (int)floorf(fy), ck = (int)floorf(fz);
        uint8_t f, fl, sh;
        if (!s.cell(ci, cj, ck, f, fl, sh)) continue;      // another site, or the generated water decides (below)
        if ((fl & kSolid) || f == 0) return 0.0f;
        if (f + sh < 255) return std::max(0.0f, ck + (sh + f) / 255.0f - fz) * WaterField::CellSize() * kRhoG; // the surface is in this cell (over its ground, #454)
        // A full cell: the body's free surface, through connected full cells of this kind (as the levelling pass joins
        // them): the highest partial cell or empty cell resting on a full one; a sealed body: its top.
        const uint8_t kind = fl & kKindMask;
        std::unordered_set<int64_t> seen;
        std::deque<std::array<int, 3>> q;
        q.push_back({ ci, cj, ck }); seen.insert(Key(ci, cj, ck));
        float level = -1e9f, top = (float)ck + 1.0f;
        static const int d[6][3] = { { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 }, { 0, 0, 1 }, { 0, 0, -1 } };
        while (!q.empty() && seen.size() < 20000) {
            const auto c = q.front(); q.pop_front();
            top = std::max(top, (float)c[2] + 1.0f);
            for (const auto& dd : d) {
                const int nx = c[0] + dd[0], ny = c[1] + dd[1], nz = c[2] + dd[2];
                uint8_t nf, nfl, nsh;
                if (!s.cell(nx, ny, nz, nf, nfl, nsh) || (nfl & kSolid)) continue;
                if (nf > 0 && (nfl & kKindMask) != kind) continue;
                if (nf + nsh >= 255) { if (seen.insert(Key(nx, ny, nz)).second) q.push_back({ nx, ny, nz }); }
                else if (nf > 0) level = std::max(level, nz + (nsh + nf) / 255.0f); // a partial cell: the surface
                else if (dd[2] == 1) level = std::max(level, (float)nz);      // open air right over a full cell
            }
        }
        if (level < -1e8f) level = top; // sealed (or past the search): its top
        return std::max(0.0f, level - fz) * WaterField::CellSize() * kRhoG;
    }
    // Generated water: the sea, or a river, lake or pond over it (as IsUnderwater).
    const auto* pf = dynamic_cast<const PlanetField*>(body.field.get());
    const float dl = Vector3Length(local);
    if (!pf || dl < 1.0f || pf->seaLevelRadius() <= 0.0f) return 0.0f;
    const Vector3 dir = Vector3Scale(local, 1.0f / dl);
    const float sea = pf->seaLevelRadius();
    float surface = sea;
    if (dl >= sea - 50.0f && pf->inlandWaterNear(dir, 0.0f)) {
        const float w = pf->waterElevation(dir);
        if (w != PlanetField::kNoWater) surface = std::max(surface, sea + w);
    }
    if (dl >= surface || dl < pf->radius + pf->elevation(dir)) return 0.0f; // over the water, or in the bed
    return (surface - dl) * kRhoG;
}

bool Fluid::Wet(const Planetoid& body, Vector3 local) {
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        uint8_t f, fl, sh; float frac;
        if (!PubCell(s, local, 0, f, fl, &frac, &sh) || (fl & kSolid)) continue;
        if (f > 0 && (sh + f) / 255.0f > frac) return true; // the water stands over the cell's ground (#454)
        // A full cell above counts too (we are in the gap at a surface's edge).
    }
    return false;
}

float Fluid::SurfaceAbove(const Planetoid& body, Vector3 local, float reach) {
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        uint8_t f, fl, sh; float frac;
        if (!PubCell(s, local, 0, f, fl, &frac, &sh) || f == 0) continue;
        float h = (sh + f) / 255.0f - frac;
        for (int k = 1; k <= (int)(reach / WaterField::CellSize()) && sh + f >= 255; ++k) {
            if (!PubCell(s, local, k, f, fl, nullptr, &sh) || (fl & kSolid)) break;
            h += f > 0 ? (sh + f) / 255.0f : 0.0f;
        }
        return h * WaterField::CellSize(); // m
    }
    return -1.0f;
}

float Fluid::VolumeAbove(const Planetoid& body, Vector3 local, float radius, float above) {
    const WaterField::SnapshotPtr snapNow = WaterField::Current();
    double vol = 0.0;
    for (const WaterField::Site& s : snapNow->sites) {
        if (s.body != body.id) continue;
        float fx, fy, fz;
        s.frame(local, fx, fy, fz);
        const int R = (int)ceilf(radius / WaterField::CellSize());
        const float aboveC = above / WaterField::CellSize(), rC = radius / WaterField::CellSize(); // cells
        for (int dk = (int)floorf(aboveC); dk <= R; ++dk)
            for (int dj = -R; dj <= R; ++dj)
                for (int di = -R; di <= R; ++di) {
                    const int ci = (int)floorf(fx) + di, cj = (int)floorf(fy) + dj, ck = (int)floorf(fz) + dk;
                    if ((float)ck < fz + aboveC) continue;                                    // the cell's bottom under the line
                    const float ex = ci + 0.5f - fx, ey = cj + 0.5f - fy, ez = ck + 0.5f - fz;
                    if (ex * ex + ey * ey + ez * ez > rC * rC) continue;
                    uint8_t f, fl;
                    uint8_t sh = 0;
                    if (s.cell(ci, cj, ck, f, fl, sh) && !(fl & (kSolid | kSource)) && CellKind(fl) == 0) { vol += f / 255.0 * WaterField::CellVolume(); if (f && EnvFlag("VOLDBG")) { uint8_t fb = 0, flb = 0, shb = 0; const bool hb = s.cell(ci, cj, ck - 1, fb, flb, shb); printf("voldbg: cell %d %d %d (+%.1f over) fill %d shut %d; below %s fill %d shut %d\n", di, dj, dk, ck - fz, f, sh, !hb ? "none" : (flb & kSolid) ? "SOLID" : "open", fb, shb); } }
                }
    }
    return (float)vol;
}

float Fluid::WaterAt(const Planetoid& body, Vector3 local, float radius) {
    float total = 0.0f;
    {
        const WaterField::SnapshotPtr snapNow = WaterField::Current();
    for (const WaterField::Site& s : snapNow->sites) {
            if (s.body != body.id) continue;
            const int R = (int)ceilf(radius / WaterField::CellSize());
            for (int dk = -R; dk <= R; ++dk)
                for (int dj = -R; dj <= R; ++dj)
                    for (int di = -R; di <= R; ++di) {
                        const Vector3 q = Vector3Add(local, Vector3Add(Vector3Add(Vector3Scale(s.e, di * WaterField::CellSize()), Vector3Scale(s.n, dj * WaterField::CellSize())), Vector3Scale(s.u, dk * WaterField::CellSize())));
                        uint8_t f, fl;
                        if (PubCell(s, q, 0, f, fl) && CellKind(fl) == 0) total += (fl & kSource) ? 100.0f : f / 255.0f * WaterField::CellVolume();
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
            const int R = (int)ceilf(radius / WaterField::CellSize());
            for (int dk = -R; dk <= R; ++dk)
                for (int dj = -R; dj <= R; ++dj)
                    for (int di = -R; di <= R; ++di) {
                        const Vector3 q = Vector3Add(local, Vector3Add(Vector3Add(Vector3Scale(s.e, di * WaterField::CellSize()), Vector3Scale(s.n, dj * WaterField::CellSize())), Vector3Scale(s.u, dk * WaterField::CellSize())));
                        uint8_t f, fl;
                        if (PubCell(s, q, 0, f, fl) && !(fl & (kSolid | kSource))) out[CellKind(fl)] += f / 255.0f * WaterField::CellVolume();
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
        const int R = (int)ceilf(range / WaterField::CellSize());
        // From 3 m under the point to 1 m over it, in metres (#730: it was 3 cells down and 1 up - at the 0.25 m default
        // (#714) a field more than 0.75 m over its water no longer counted as watered).
        const int down = (int)ceilf(3.0f / WaterField::CellSize()), upC = (int)ceilf(1.0f / WaterField::CellSize());
        for (int dk = -down; dk <= upC; ++dk)
            for (int dj = -R; dj <= R; dj += 1)
                for (int di = -R; di <= R; di += 1) {
                    const Vector3 q = Vector3Add(local, Vector3Add(Vector3Add(Vector3Scale(s.e, di * WaterField::CellSize()), Vector3Scale(s.n, dj * WaterField::CellSize())), Vector3Scale(s.u, dk * WaterField::CellSize())));
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
    fprintf(f, "water %zu cell %.4f\n", mine.size(), WaterField::CellSize()); // cell edge (m): a site saved at another size is dropped on load (W4c)
    for (Site* s : mine) {
        std::vector<Brick*> wet;
        for (auto& kv : s->bricks) {
            Brick& b = *kv.second;
            bool any = false;
            for (int i = 0; i < B3 && !any; ++i) any = b.fill[i] && !(b.flags[i] & kSource);
            if (any || b.basin || b.reach) wet.push_back(&b); // a drained pond's bricks too: loaded fresh they would hold its water again (#360)
        }
        fprintf(f, "site %.4f %.4f %.4f %.6f %.6f %.6f %.6f %.6f %.6f %zu %zu", s->origin.x, s->origin.y, s->origin.z, s->e.x, s->e.y, s->e.z, s->u.x, s->u.y, s->u.z, (size_t)0, wet.size()); // 0: claimed ponds are gone (#360)
        fprintf(f, "\n");
        fprintf(f, "springs %zu\n", s->springs.size());
        for (const Site::Spring& sp : s->springs) fprintf(f, "%d %d %d %d %.4f %.3f r %d\n", sp.ci, sp.cj, sp.ck, sp.kind, sp.perStep, sp.head, sp.rel);
        fprintf(f, "releases %zu\n", s->releases.size()); // released river reaches (#361)
        for (const Site::Release& r : s->releases)
            fprintf(f, "%d %d %.2f %.4f %.4f %.4f %.6f %.6f %.6f q %.3f %d %.4f u %d b %d\n", r.ci, r.cj, r.r, r.centre.x, r.centre.y, r.centre.z, r.down.x, r.down.y, r.down.z, r.natQ, r.piece, r.share, r.upPiece, r.peakBlocked);
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
    // River flow (#656): the links off their natural flow and the reaches setting one.
    std::vector<std::pair<int, float>> cur, set;
    if (const RiverNetwork* net = RiversOf(bodyId)) net->flowState(cur, set);
    fprintf(f, "riverflow %zu %zu\n", cur.size(), set.size());
    for (const auto& [i, v] : cur) fprintf(f, "%d %.4f\n", i, v);
    for (const auto& [i, v] : set) fprintf(f, "%d %.4f\n", i, v);
}

void Fluid::Load(FILE* f, uint64_t bodyId) {
    int n = 0;
    if (fscanf(f, " water %d", &n) != 1) return;
    float cell = 1.0f; // saves before W4c: 1 m cells
    if (fscanf(f, " cell %f", &cell) != 1) cell = 1.0f;
    const bool keep = fabsf(cell - WaterField::CellSize()) < 1e-3f; // else read and dropped: the sites start over from the world
    const auto waitT0 = std::chrono::steady_clock::now();
    WaitIdle();
    if (getenv("RAYLIBGAME_LOAD_PROFILE") && getenv("RAYLIBGAME_LOAD_PROFILE")[0] == '1')
        printf("load: fluid waited %.0f ms for the running step (%d sites in the save)\n", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - waitT0).count(), n);
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
                if (fscanf(f, " %d %d %d %d %f %f", &sp.ci, &sp.cj, &sp.ck, &sp.kind, &sp.perStep, &sp.head) == 6) {
                    const long at = ftell(f); // its released reach (#656; older saves: none - and "releases" must not lose its r)
                    if (fscanf(f, " r %d", &sp.rel) != 1) { sp.rel = -1; fseek(f, at, SEEK_SET); }
                    s->springs.push_back(sp);
                }
            }
        int nr = 0; // released river reaches (#361; older saves have none)
        if (fscanf(f, " releases %d", &nr) == 1)
            for (int j = 0; j < nr; ++j) {
                Site::Release r{};
                if (fscanf(f, " %d %d %f %f %f %f %f %f %f", &r.ci, &r.cj, &r.r, &r.centre.x, &r.centre.y, &r.centre.z, &r.down.x, &r.down.y, &r.down.z) == 9) {
                    if (fscanf(f, " q %f %d %f", &r.natQ, &r.piece, &r.share) != 3) { r.natQ = 0.0f; r.piece = -1; r.share = 1.0f; } // #656 (older saves: none)
                    else { const long at = ftell(f); if (fscanf(f, " u %d", &r.upPiece) != 1) { r.upPiece = -1; fseek(f, at, SEEK_SET); } else { const long at2 = ftell(f); if (fscanf(f, " b %d", &r.peakBlocked) != 1) { r.peakBlocked = 0; fseek(f, at2, SEEK_SET); } } } // b: the most ground the dam put in the river (#656 follow-up)
                    s->releases.push_back(r);
                }
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
        if (keep) sites.push_back(std::move(s));
    }
    // River flow (#656; older saves: none - back to where the line would start).
    const long at = ftell(f);
    int nc = 0, ns = 0;
    if (fscanf(f, " riverflow %d %d", &nc, &ns) != 2) { fseek(f, at, SEEK_SET); return; }
    std::vector<std::pair<int, float>> cur, set;
    for (int i = 0; i < nc + ns; ++i) {
        int k = 0; float v = 0.0f;
        if (fscanf(f, " %d %f", &k, &v) != 2) break;
        (i < nc ? cur : set).push_back({ k, v });
    }
    // Only with flow state to restore: RiversOf builds the whole river network on first use (~470 ms), which an empty
    // "riverflow 0 0" line made every load pay up front (#737; it is built later anyway, when water first needs it).
    if (!cur.empty() || !set.empty())
        if (const RiverNetwork* net = RiversOf(bodyId)) net->setFlowState(cur, set);
}

Fluid::Stats Fluid::GetStats() { std::lock_guard<std::mutex> p(pubMx); return pubStats; }
bool Fluid::GpuMode() { return gpuMode; }

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
        if (!PubCell(s, local, 0, f, fl)) continue;
        return TextFormat("fill %d flags %s%s%s%s%s", f, (fl & kSolid) ? "S" : "", (fl & kSource) ? "R" : "", (fl & kStatic) ? "T" : "", (fl & kDug) ? "D" : "", (fl & WaterField::kSink) ? "K" : "");
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
