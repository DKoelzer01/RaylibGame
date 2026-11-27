// space/luna-walk (#1942, Plan 11 #667): the player stands, jumps and pours water on Luna under Luna's own gravity, not
// Terra's. Luna (worlds/luna.txt) is 10 km with 0.98 m/s^2; Terra 5.9. The player is put on Luna's surface (spawnOn), waits
// to land, then is pushed up at 5 m/s as a jump does (pendingPush): the apex over the take-off radius gives the gravity felt,
// g = v^2 / (2 h) (Luna: ~12.8 m, Terra's pull would give ~2.1 m), which must be within 5% of Luna's surfaceGravity. Then
// 0.5 m^3 of water is poured beside the player: the simulation must hold it (water works on the moon).
// No c.camera(): it turns fly mode on, and the player's physics is what is tested.
#include "scenario.h"
#include "fluid.h"
#include "player.h"
#include "world.h"
#include "worldtypes.h"
#include <raymath.h>
#include <algorithm>
#include <cmath>

namespace {
struct State {
    Planetoid* luna = nullptr;
    double r0 = 0.0, apex = 0.0, pourAt = -1.0;
    float vol0 = 0.0f, poured = 0.0f;
    bool jumped = false, airborne = false;
} S;
double FeetRadius() { return S.luna->toLocalD(player.position).length(); } // m from Luna's centre
Vector3 UpLocal() { return Vector3Normalize(S.luna->toLocalD(player.position).toFloat()); }
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
        return true;
    });
    s.step("land", 90.0, [](ScenarioCtx& c) {
        c.waitingFor("standing on Luna in an idle world");
        player.flyMode = false;
        return c.t > 2.0 && player.grounded && WorldIdle() && player.nearest == S.luna;
    });
    s.step("jump (5 m/s up)", 30.0, [](ScenarioCtx& c) {
        if (!S.jumped) {
            S.jumped = true;
            S.r0 = FeetRadius();
            player.pendingPush = S.luna->dirToWorld(Vector3Scale(UpLocal(), 5.0f));
            return false;
        }
        const double h = FeetRadius() - S.r0;
        S.apex = std::max(S.apex, h);
        if (!player.grounded) S.airborne = true;
        c.waitingFor("back on the ground");
        return S.airborne && player.grounded && c.t > 1.0;
    });
    s.step("pour 0.5 m^3 beside the player", 20.0, [](ScenarioCtx& c) {
        if (S.pourAt < 0.0) {
            S.pourAt = c.t;
            S.vol0 = (float)Fluid::GetStats().volume;
            const Vector3 up = UpLocal(), side = Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 }));
            const Vector3 at = Vector3Add(S.luna->toLocalD(player.position).toFloat(), Vector3Add(Vector3Scale(side, 3.0f), Vector3Scale(up, 1.0f)));
            S.poured = Fluid::Pour(*S.luna, at, 0.5f);
            return false;
        }
        c.waitingFor("10 s of water on the moon");
        return c.t - S.pourAt > 10.0;
    });
    s.step("report", 1.0, [](ScenarioCtx& c) {
        const double g = S.apex > 0.01 ? 25.0 / (2.0 * S.apex) : 0.0;
        const double held = Fluid::GetStats().volume - S.vol0;
        c.note(TextFormat("jump apex %.2f m (Luna's gravity predicts %.2f m, Terra's 5.9 would give %.2f m); felt %.3f m/s^2; "
                          "water: %.2f m^3 queued, %.2f m^3 held after 10 s",
                          S.apex, 25.0 / (2.0 * S.luna->surfaceGravity), 25.0 / (2.0 * 5.9), g, S.poured, held));
        c.metric("jump apex", S.apex, "m", false);
        c.near("gravity felt on Luna (m/s^2)", g, S.luna->surfaceGravity, 0.05 * S.luna->surfaceGravity);
        c.gt("water poured on Luna is held (m^3)", held, 0.4);
        return true;
    });
}
