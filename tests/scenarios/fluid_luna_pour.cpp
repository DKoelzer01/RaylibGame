// fluid/luna-pour (#1942, Plan 11 #667): water flows by its body's gravity. The kernels' g is 9.81 x the body's surface
// gravity over Terra's (Fluid::GravityScale: exactly 1 on Terra, 0.98 / 5.886 on Luna). The same 1 m^3 is poured 1 m over
// the ground beside the player on Terra, then on Luna; after 1 s of fluid time (ticks x 25 ms, wall speed varies) the share
// still within 1.5 m of the pour is noted for both (Luna's water spreads slower: it should keep more). Checked: the scales,
// and the water within 30 m of each pour against the mass ledger (#1885) after 10 s. Not a parity scenario (the pours wait on landing and
// streaming, wall-clock): fluid/mass-ledger-luna compares the GPU and the CPU at Luna's scale.
// No c.camera(): it turns fly mode on, and the player's own placement is used.
#include "scenario.h"
#include "fluid.h"
#include "player.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <cmath>

namespace {
struct Pour {
    Planetoid* body = nullptr;
    Vector3 at{};
    long long tick0 = -1;
    float queued = 0.0f, near1s = -1.0f;
    double vol0 = 0.0;
    Fluid::Ledger l0{};
};
struct State { Pour terra, luna; } S;

// The pour beside the player (3 m to the side, 1 m up), body-local.
Vector3 BesidePlayer(Planetoid& b) {
    const Vector3 feet = b.toLocalD(player.position).toFloat(), up = Vector3Normalize(feet);
    const Vector3 side = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.z) < 0.9f ? Vector3{ 0, 0, 1 } : Vector3{ 1, 0, 0 }));
    return Vector3Add(feet, Vector3Add(Vector3Scale(side, 3.0f), up));
}

// The simulated water within 30 m of a pour (m^3). Each body's water is judged by its own box: the world's total also
// moved with the other body's water (Terra's pour, then left behind while the player stood on Luna).
double NearPour(const Pour& p) {
    const Vector3 up = Vector3Normalize(p.at);
    const Vector3 ax = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.z) < 0.9f ? Vector3{ 0, 0, 1 } : Vector3{ 1, 0, 0 }));
    const Vector3 ay = Vector3CrossProduct(up, ax);
    return Fluid::VolumeIn(*p.body, p.at, ax, ay, up, { 30.0f, 30.0f, 15.0f });
}

void AddSteps(ScenarioDef& s, Pour State::*which, const char* name) {
    s.step(TextFormat("stand on %s", name), 90.0, [which](ScenarioCtx& c) {
        Pour& p = S.*which;
        c.waitingFor("standing in an idle world");
        player.flyMode = false;
        return c.t > 2.0 && player.grounded && WorldIdle() && player.nearest == p.body;
    });
    s.step(TextFormat("pour 1 m^3 on %s", name), 240.0, // (fluid time runs at ~0.1 x sim on the headless CPU kernel) [which](ScenarioCtx& c) {
        Pour& p = S.*which;
        if (p.tick0 < 0) {
            p.at = BesidePlayer(*p.body);
            p.vol0 = NearPour(p);
            p.l0 = Fluid::GetLedger();
            p.tick0 = 0;
            p.queued = Fluid::Pour(*p.body, p.at, 1.0f);
            return false;
        }
        // The clock starts once the pour is in the ledger: read before it, the tick count can be stale (no water was
        // active on Terra yet) and jumped 10 s at once, judging a pour that never ran (#1942).
        if (p.tick0 == 0) {
            c.waitingFor("the pour in the simulation");
            if (Fluid::GetLedger().poured - p.l0.poured < 0.5) return false;
            p.tick0 = Fluid::GetStats().ticks;
        }
        const double t = (double)(Fluid::GetStats().ticks - p.tick0) * 0.025;
        if (p.near1s < 0.0f && t >= 1.0) p.near1s = Fluid::WaterAt(*p.body, p.at, 1.5f) + Fluid::WaterAt(*p.body, Vector3Subtract(p.at, Vector3Scale(Vector3Normalize(p.at), 1.0f)), 1.5f);
        c.waitingFor(TextFormat("10 s of fluid time (%.1f s)", t));
        return t >= 10.0;
    });
    s.step(TextFormat("judge %s", name), 1.0, [which, name](ScenarioCtx& c) {
        Pour& p = S.*which;
        const Fluid::Ledger l = Fluid::GetLedger(), &l0 = p.l0;
        const double held = NearPour(p) - p.vol0;
        const double net = (l.poured - l0.poured) - (l.taken - l0.taken) - (l.dried - l0.dried);
        c.note(TextFormat("%s: gravity %.3f m/s^2, fluid scale %.4f; %.2f m^3 queued, near the pour after 1 s %.3f m^3; after 10 s "
                          "held %.3f, ledger poured %.3f, dried %.3f, net %.3f", name, p.body->surfaceGravity, Fluid::GravityScale(*p.body),
                          p.queued, p.near1s, held, l.poured - l0.poured, l.dried - l0.dried, net));
        c.metric(TextFormat("%s water near the pour after 1 s", name), p.near1s, "m^3", true);
        c.gt(TextFormat("water poured on %s (m^3)", name), l.poured - l0.poured, 0.8);
        c.near(TextFormat("water on %s matches the ledger (m^3)", name), held, net, 0.01);
        return true;
    });
}
}

SCENARIO("fluid/luna-pour") {
    s.tags = { "fluid", "space", "frames" };
    // The ledger leaves out the kernels' evaporation (fluid.cpp, #454): Luna's pour spreads into films, and 0.03-0.16 m^3 of it
    // dried uncounted in 10 s. Without evaporation every unit is accounted for, as in fluid/mass-ledger.
    s.env = { { "RAYLIBGAME_FLUID_NOEVAP", "1" } };
    s.world = ScenarioWorld::Planet;
    s.timeout = 840.0;
    s.step("find Terra and Luna", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        S = State{};
        for (Planetoid* b : Planetoid::All()) {
            if (b->kind == BodyKind::Moon && b->name == "luna") S.luna.body = b;
            if (b->kind == BodyKind::Planet && !b->orbitParent && !S.terra.body) S.terra.body = b;
        }
        if (!S.luna.body || !S.terra.body) { c.fail("no Terra or no body named luna"); return true; }
        // Terra's water is the tuned reference: its scale is exactly 1 (bit-identical to before #1942).
        c.eq("Terra's fluid gravity scale", (double)Fluid::GravityScale(*S.terra.body), 1.0);
        c.near("Luna's fluid gravity scale", (double)Fluid::GravityScale(*S.luna.body),
               (double)(S.luna.body->surfaceGravity / S.terra.body->surfaceGravity), 1e-6);
        player.spawnOn(c.world, *S.terra.body, { 0, 1, 0 });
        player.flyMode = false;
        return true;
    });
    AddSteps(s, &State::terra, "Terra");
    s.step("go to Luna", 5.0, [](ScenarioCtx& c) {
        player.spawnOn(c.world, *S.luna.body, { 0, 1, 0 });
        player.flyMode = false;
        return true;
    });
    AddSteps(s, &State::luna, "Luna");
    s.step("compare", 1.0, [](ScenarioCtx& c) {
        c.note(TextFormat("near the pour after 1 s: Terra %.3f m^3, Luna %.3f m^3 (Luna's slower water should keep more)",
                          S.terra.near1s, S.luna.near1s));
        return true;
    });
}
