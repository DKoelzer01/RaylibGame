// space/luna-walk (#1942, Plan 11 #667): the player stands, jumps and pours water on Luna under Luna's own gravity, not
// Terra's. Luna (worlds/luna.txt) is 10 km with 0.98 m/s^2; Terra 5.9. The player is put on Luna's surface (spawnOn), waits
// to land, then is jumps (autoJump: a Space press, 5 m/s): the apex over the take-off radius gives the gravity felt,
// g = v^2 / (2 h) (Luna: ~12.8 m, Terra's pull would give ~2.1 m), which must be within 5% of Luna's surfaceGravity. Then
// 0.5 m^3 of water is poured beside the player: the simulation must hold what the mass ledger says (#1885).
// No c.camera(): it turns fly mode on, and the player's physics is what is tested.
#include "scenario.h"
#include "fluid.h"
#include "player.h"
#include "survival.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <string>

namespace {
struct State {
    Planetoid* luna = nullptr;
    double r0 = 0.0, apex = 0.0, pourAt = -1.0;
    double spawnR = 0.0;                   // the feet's radius where spawnOn put them
    double settleR = -1.0, settleAt = 0.0; // 'land': the feet's radius and since when it held still
    int spawns0 = 0, deaths0 = 0;          // after the move to Luna: a change later is a respawn (#2000)
    // (#2000) The landing as it went, for a death: the lowest feet radius, the fastest speed against the ground, and the
    // time and speed at the first touchdown.
    double minR = 1e30, maxRel = 0.0, touchAt = -1.0, touchRel = 0.0;
    bool wasGrounded = false;
    float vol0 = 0.0f, poured = 0.0f;
    Fluid::Ledger l0{}; // the mass ledger (#1885) at the pour
    bool jumped = false, airborne = false;
    std::string trace; // the jump's height over time (s:m)
} S;
double FeetRadius() { return S.luna->toLocalD(player.position).length(); } // m from Luna's centre
Vector3 UpLocal() { return Vector3Normalize(S.luna->toLocalD(player.position).toFloat()); }
// (#2000) A respawn puts the player back on Terra, ~600 km from Luna's centre, and the jump then never "lands". Say so.
bool Respawned(ScenarioCtx& c) {
    if (player.spawns == S.spawns0 && !survival.dead) return false;
    c.fail(TextFormat("the player died or was respawned (deaths %d -> %d, last harm '%s' at %.1f m/s, health %.0f, oxygen %.0f, vacuum %.0f s); "
                      "spawned at r %.2f m, lowest %.2f, now %.2f; fastest against the ground %.2f m/s; first touchdown at %.2f s, %.2f m/s",
                      S.deaths0, survival.deaths, survival.lastHarm, survival.lastFallSpeed, survival.health, survival.oxygen, survival.vacuumTime,
                      S.spawnR, S.minR, FeetRadius(), S.maxRel, S.touchAt, S.touchRel));
    return true;
}
// (#2000) Once a frame while landing: the radius and the speed against Luna's turning, orbiting ground.
void TrackLanding(ScenarioCtx& c) {
    // (F8b-3) The player's own reading: on the frame path its velocity is against the ground already (subtracting the
    // ground's motion again reported Luna's ~1 km/s orbit as the landing speed).
    const double rel = Vector3Length(player.groundVelocity(c.world));
    S.minR = std::min(S.minR, FeetRadius());
    S.maxRel = std::max(S.maxRel, rel);
    if (player.grounded && !S.wasGrounded && S.touchAt < 0.0) { S.touchAt = c.t; S.touchRel = rel; }
    S.wasGrounded = player.grounded;
}
}

SCENARIO("space/luna-walk") {
    s.tags = { "space", "physics", "frames" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;
    s.step("go to Luna", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        S = State{};
        for (Planetoid* b : Planetoid::All()) if (b->kind == BodyKind::Moon && b->name == "luna") S.luna = b;
        if (!S.luna) { c.fail("no body named luna"); return true; }
        c.note(TextFormat("Luna: radius %.0f m, surface gravity %.3f m/s^2", S.luna->radius * S.luna->scale, S.luna->surfaceGravity));
        player.spawnOn(c.world, *S.luna, { 0, 1, 0 });
        player.flyMode = false;
        S.spawns0 = player.spawns;
        S.deaths0 = survival.deaths;
        S.spawnR = FeetRadius();
        return true;
    });
    // spawnOn used to raycast the terrain built so far: on the runner a coarse mesh put the feet 2 m inside Luna's ground,
    // and the push out when the fine chunks came was a fatal "fall" (#2000; spawnOn now marches the field). Wait until the
    // feet hold still on the ground for 1 s.
    s.step("land", 90.0, [](ScenarioCtx& c) {
        player.flyMode = false;
        TrackLanding(c);
        if (Respawned(c)) return true;
        const double r = FeetRadius();
        if (!(player.grounded && WorldIdle() && player.nearest == S.luna) || std::fabs(r - S.settleR) > 0.05) { S.settleR = r; S.settleAt = c.t; }
        c.waitingFor(TextFormat("standing still on Luna in an idle world for 1 s (grounded %d, idle %d, still for %.1f s)", (int)player.grounded,
                                (int)WorldIdle(), c.t - S.settleAt));
        if (c.t <= 2.0 || c.t - S.settleAt <= 1.0) return false;
        c.note(TextFormat("spawnOn put the feet at %.1f m from Luna's centre; they stand at %.1f m", S.spawnR, r));
        return true;
    });
    s.step("jump (5 m/s up)", 30.0, [](ScenarioCtx& c) {
        if (Respawned(c)) return true;
        if (!S.jumped) {
            S.jumped = true;
            S.r0 = FeetRadius();
            player.autoJump = true; // a Space press: the jump key's own path (jumpSpeed 5 m/s)
            return false;
        }
        const double h = FeetRadius() - S.r0;
        S.apex = std::max(S.apex, h);
        if ((int)(c.t * 60.0) % 15 == 0) S.trace += TextFormat("%.2f:%.2f ", c.t, h); // the height every 0.25 s
        if (!player.grounded) S.airborne = true;
        c.waitingFor(TextFormat("back on the ground (apex %.2f m, now %.2f m, airborne %d, grounded %d, fly %d, nearest luna %d)", S.apex, h, (int)S.airborne, (int)player.grounded, (int)player.flyMode, (int)(player.nearest == S.luna)));
        return S.airborne && player.grounded && c.t > 1.0;
    });
    s.step("pour 0.5 m^3 beside the player", 20.0, [](ScenarioCtx& c) {
        if (S.pourAt < 0.0) {
            S.pourAt = c.t;
            S.vol0 = (float)Fluid::GetStats().volume;
            S.l0 = Fluid::GetLedger();
            const Vector3 up = UpLocal(), side = Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 }));
            const Vector3 at = Vector3Add(S.luna->toLocalD(player.position).toFloat(), Vector3Add(Vector3Scale(side, 3.0f), Vector3Scale(up, 1.0f)));
            S.poured = Fluid::Pour(*S.luna, at, 0.5f);
            return false;
        }
        if (Respawned(c)) return true;
        c.waitingFor("10 s of water on the moon");
        return c.t - S.pourAt > 10.0;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        const double g = S.apex > 0.01 ? 25.0 / (2.0 * S.apex) : 0.0;
        const double held = Fluid::GetStats().volume - S.vol0;
        const Fluid::Ledger l = Fluid::GetLedger(), &l0 = S.l0;
        const double net = (l.poured - l0.poured) - (l.taken - l0.taken) - (l.dried - l0.dried); // what the ledger says is there
        c.note(TextFormat("jump apex %.2f m (Luna's gravity predicts %.2f m, Terra's 5.9 would give %.2f m); felt %.3f m/s^2; "
                          "water: %.2f m^3 queued, %.3f m^3 held after 10 s; ledger: poured %.3f, dried %.3f, net %.3f",
                          S.apex, 25.0 / (2.0 * S.luna->surfaceGravity), 25.0 / (2.0 * 5.9), g, S.poured, held,
                          l.poured - l0.poured, l.dried - l0.dried, net));
        c.note("height over the jump (s:m): " + S.trace);
        c.metric("jump apex", S.apex, "m", false);
        c.near("gravity felt on Luna (m/s^2)", g, S.luna->surfaceGravity, 0.05 * S.luna->surfaceGravity);
        c.gt("water poured on Luna (m^3)", l.poured - l0.poured, 0.4);
        // Kept by the simulation as its ledger says (#1885: thin films dry, but nothing else may vanish).
        c.near("water on Luna matches the ledger (m^3)", held, net, 0.01);
        return true;
    });
}
