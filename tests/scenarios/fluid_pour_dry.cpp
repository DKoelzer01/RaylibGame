// fluid/pour-dry (#1944, user 2026-11-27: "slow water film drying down to minutes"): 0.5 m^3 poured on the sandbox's flat
// ground stands as a puddle a few cm deep with films at its edge; both (shallower than 5 cm, one cell deep) dry a unit (~1 mm)
// every 20 s, so it stays a puddle for minutes. Reported: the water held at 10 s, 60 s, 5 and 8 min of fluid time (ticks x
// 25 ms; wall speed varies). Checked: most of it at 60 s, drying at 5 min, mostly dry by 8 min, and the held water against
// the mass ledger (#1885: what dried is counted). Before #1944 a film lost a unit every 2 s and puddles never dried.
#include "scenario.h"
#include "fluid.h"
#include "player.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
struct State {
    Vector3 at{};
    long long tick0 = -1;
    float queued = 0.0f;
    double v0 = 0.0, held[4] = { -1.0, -1.0, -1.0, -1.0 };
    Fluid::Ledger l0{};
} S;
constexpr double kMarks[4] = { 10.0, 60.0, 300.0, 480.0 }; // s of fluid time
double FluidT() { return (double)(Fluid::GetStats().ticks - S.tick0) * 0.025; }
}

SCENARIO("fluid/pour-dry") {
    s.tags = { "fluid", "water" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 780.0;
    s.step("pour 0.5 m^3 on flat ground", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        S = State{};
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        const Vector3 up = { 0, 1, 0 }, ahead = { 0, 0, 1 }; // the sandbox is flat: body-local +y is up
        S.at = Vector3Add(feet, Vector3Add(Vector3Scale(ahead, 6.0f), Vector3Scale(up, 0.5f)));
        S.v0 = Fluid::GetStats().volume;
        S.l0 = Fluid::GetLedger();
        S.tick0 = Fluid::GetStats().ticks;
        S.queued = Fluid::Pour(*c.body, S.at, 0.5f);
        // Looking down at the puddle from 6 m (no player physics here: c.camera's fly mode does no harm).
        c.camera(Vector3Add(S.at, Vector3Add(Vector3Scale(up, 5.0f), Vector3Scale(ahead, -4.0f))), S.at, up);
        return true;
    });
    s.step("wait 8 min of fluid time", 720.0, [](ScenarioCtx& c) {
        const double t = FluidT();
        for (int k = 0; k < 4; ++k)
            if (S.held[k] < 0.0 && t >= kMarks[k]) {
                S.held[k] = Fluid::GetStats().volume - S.v0;
                c.note(TextFormat("at %.0f s: %.3f m^3 held (%.0f%% of the pour)", kMarks[k], S.held[k], 100.0 * S.held[k] / std::max(0.01f, S.queued)));
            }
        c.waitingFor(TextFormat("fluid time %.0f of 480 s (%.3f m^3 held)", t, Fluid::GetStats().volume - S.v0));
        return S.held[3] >= 0.0;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        const Fluid::Ledger l = Fluid::GetLedger(), &l0 = S.l0;
        const double net = (l.poured - l0.poured) - (l.taken - l0.taken) - (l.dried - l0.dried);
        const double poured = l.poured - l0.poured;
        c.note(TextFormat("0.5 m^3 poured (%.3f queued): held %.3f at 10 s, %.3f at 60 s, %.3f at 5 min, %.3f at 8 min; ledger poured %.3f, dried %.3f",
                          S.queued, S.held[0], S.held[1], S.held[2], S.held[3], poured, l.dried - l0.dried));
        c.metric("held at 10 s", S.held[0], "m^3", false);
        c.metric("held at 60 s", S.held[1], "m^3", false);
        c.metric("held at 5 min", S.held[2], "m^3", false);
        c.metric("held at 8 min", S.held[3], "m^3", true);
        c.gt("water poured (m^3)", poured, 0.45);
        c.gt("held at 10 s (share of the pour)", S.held[0] / poured, 0.9);
        c.gt("held at 60 s: most of it (share)", S.held[1] / poured, 0.5);
        c.lt("held at 5 min: drying (share)", S.held[2] / poured, 0.7);
        c.lt("held at 8 min: mostly dry (share)", S.held[3] / poured, 0.2);
        c.near("water held matches the ledger (m^3)", Fluid::GetStats().volume - S.v0, net, 0.01);
        return true;
    });
}
