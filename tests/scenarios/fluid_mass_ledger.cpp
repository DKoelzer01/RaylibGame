// fluid/mass-ledger (#1885): water is neither made nor lost by the simulation. In the sandbox, a closed 2 x 2 m pen of stone-brick
// cubes (plot cells, 0.5 m; walls 1 m high) on a floating floor; 1.9 m^3 is poured into one corner over 2 s (it sloshes across),
// 0.5 m^3 is taken from the middle at 10 s (Fluid::Take), and it settles until 40 s. Evaporation is off
// (RAYLIBGAME_FLUID_NOEVAP=1). Twice a second the simulation's volume (Fluid::GetStats().volume) is compared with the CPU
// ledger (Fluid::GetLedger: poured - taken - dried): their difference may drift by at most 0.5% of what was poured. Samples
// within 1 s of a pour or a take are skipped (the ledger counts a command when it is applied, the stats when its bricks are
// next published). The -cpu twin runs the CPU kernel; "volume at the end" in both is the CPU/GPU parity scalar.
// Not in the ledger: the kernels' own source refills and sink intake (no sources or sinks here).
#include "scenario.h"
#include "fluid.h"
#include "plots.h"
#include "player.h"
#include "buildtool.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <memory>

namespace {
constexpr int kFloor = 8; // floor cells across (4 m)
constexpr int kPen = 5;   // the pen's walls at cells 0 and kPen (inside: 1..kPen-1, 2 m)
struct State {
    Plot* plot = nullptr;
    Int3 base{};
    double nextPour = 0.0, lastOp = 0.0;
    float poured = 0.0f, taken = 0.0f;
    bool tookOnce = false;
    double v0 = 0.0, worst = 0.0;
    Fluid::Ledger l0{};
    int samples = 0;
    Vector3 at(float x, float y, float z) const { return plot->toBody({ x, y, z }); } // plot metres -> body-local
    float floorTop() const { return (base.y + 1) * kCellSize; }
    float cell(int c, int b) const { return (b + c + 0.5f) * kCellSize; } // a pen cell's centre (plot metres)
} S;
double LedgerNet(const Fluid::Ledger& l) { return l.poured - l.taken - l.dried; }

void Build(ScenarioDef& s) {
    s.tags = { "fluid", "water" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 180.0;
    s.step("build the closed pen", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        S = State{};
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        S.plot = c.body->plots->plotFor(feet, 0.0f, true);
        S.base = S.plot->cellAt(S.plot->toPlot(feet));
        S.base.x += 4; S.base.y += 2; S.base.z -= kFloor / 2 + 3; // ahead of the player, 1 m up (clear of the ground)
        const uint8_t mat = (uint8_t)std::max(0, BuildMaterialForItem("stone_brick"));
        int placed = 0, want = 0;
        auto put = [&](int x, int y, int z) { ++want; placed += S.plot->place(Part{ PartType::Cube, { S.base.x + x, S.base.y + y, S.base.z + z }, 0, mat }) >= 0; };
        for (int x = 0; x < kFloor; ++x)
            for (int z = 0; z < kFloor; ++z) put(x, 0, z);
        for (int y = 1; y <= 2; ++y)
            for (int i = 0; i <= kPen; ++i) {
                put(i, y, 0); put(i, y, kPen);                    // two walls, corners included
                if (i > 0 && i < kPen) { put(0, y, i); put(kPen, y, i); } // the other two
            }
        c.eq("cubes placed", placed, want);
        // The camera over the pen, looking down into it (no player physics here, so c.camera's fly mode does no harm).
        const float mid = S.cell(kPen / 2, S.base.x) + 0.25f, midz = S.cell(kPen / 2, S.base.z) + 0.25f;
        c.camera(S.at(mid + 2.5f, S.floorTop() + 3.0f, midz + 2.5f), S.at(mid, S.floorTop(), midz), Vector3Normalize(Vector3Subtract(S.at(0, 1, 0), S.at(0, 0, 0))));
        return true;
    });
    s.step("pour 1.9 m^3 into a corner, take 0.5 m^3 at 10 s, settle to 40 s", 60.0, [](ScenarioCtx& c) {
        const float cx = S.cell(1, S.base.x), cz = S.cell(1, S.base.z);
        if (c.t >= S.nextPour && c.t < 2.0) {
            S.nextPour = c.t + 0.1;
            S.poured += Fluid::Pour(*c.body, S.at(cx, S.floorTop() + 0.6f, cz), 0.1f);
            S.lastOp = c.t;
        }
        if (!S.tookOnce && c.t >= 10.0) {
            S.tookOnce = true;
            const float mid = (S.base.x + kPen * 0.5f) * kCellSize, midz = (S.base.z + kPen * 0.5f) * kCellSize;
            S.taken = Fluid::Take(*c.body, S.at(mid, S.floorTop() + 0.1f, midz), 0.5f, 1.5f);
            S.lastOp = c.t;
        }
        if (c.t > 40.0) return true;
        if ((int)(c.t * 60.0) % 30 != 0 || c.t - S.lastOp < 1.0) return false; // twice a second, 1 s clear of a pour or take
        const double v = Fluid::GetStats().volume;
        const Fluid::Ledger l = Fluid::GetLedger();
        if (S.samples++ == 0) { S.v0 = v; S.l0 = l; } // (the sandbox's own water, if any, is in both)
        const double err = (v - S.v0) - (LedgerNet(l) - LedgerNet(S.l0));
        S.worst = std::max(S.worst, std::fabs(err));
        return false;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        const Fluid::Ledger l = Fluid::GetLedger();
        const double v = Fluid::GetStats().volume;
        c.note(TextFormat("asked: %.3f m^3 poured, %.3f taken; ledger: poured %.4f, taken %.4f (free %.4f), spring %.4f, dried %.4f, "
                          "dropped %.4f; volume %.4f m^3; %d samples, largest drift from the ledger %.5f m^3",
                          S.poured, S.taken, l.poured, l.taken, l.takenFree, l.spring, l.dried, Fluid::PourDropped(), v, S.samples, S.worst));
        c.metric("volume at the end", v, "m^3", false);
        c.metric("largest drift from the ledger", S.worst, "m^3");
        c.gt("water poured", l.poured, 1.5);
        c.gt("water taken", l.taken, 0.2);
        c.gt("samples", S.samples, 40);
        c.lt("drift from the ledger under 0.5% of the poured", S.worst, 0.005 * l.poured);
        return true;
    });
}
}

SCENARIO("fluid/mass-ledger") { s.env = { { "RAYLIBGAME_FLUID_NOEVAP", "1" } }; Build(s); }
SCENARIO("fluid/mass-ledger-cpu") { s.env = { { "RAYLIBGAME_FLUID_NOEVAP", "1" }, { "FLUID_CPU", "1" } }; Build(s); }
// At Luna's fluid gravity (#1942: 0.98 / 5.886): the scaled kernels keep the ledger too, and (--scenario-parity) the GPU
// and the CPU agree bit for bit with k != 1.
SCENARIO("fluid/mass-ledger-luna") { s.env = { { "RAYLIBGAME_FLUID_NOEVAP", "1" }, { "RAYLIBGAME_FLUID_GSCALE", "0.16650" } }; Build(s); }
