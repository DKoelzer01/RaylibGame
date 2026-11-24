// world/spawn-by-water (#1787): a spawn (placeAt, as a teleport or a bench does) while water is moving releases the player
// promptly, and the player is under gravity after it. --test-pond's swim spawn hit a loading screen that never saw an idle
// world (20 s, #1793) and then judged its swim on the first frame after it, with gravity still 0 (the player isn't updated
// while the screen is up). On the sandbox, 300 m^3 are poured on flat ground (it spreads for a long while), then the player
// is put down 40 m off, on dry ground, 1 m up.
#include "scenario.h"
#include "fluid.h"
#include "player.h"
#include "playersession.h"
#include "world.h"
#include <raymath.h>
#include <cmath>

namespace {
struct State {
    Vector3 home{}, side{}, spot{};
    double placedWall = 0.0; // the wall clock when the player was put down
    double h0 = 0.0;         // height over the centre then (m)
} S;
double Height(ScenarioCtx& c) { return Vector3Length(c.body->toLocalOf(LocalSession().player.position, c.world.renderOrigin)); }
}

SCENARIO("world/spawn-by-water") {
    s.tags = { "world", "fluid" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("settle", 30.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("pour 300 m^3 on flat ground, then put the player down 40 m off", 2.0, [](ScenarioCtx& c) {
        S = State{};
        S.home = c.field->findSpawnDir(c.body->seed);
        S.side = Vector3Normalize(Vector3CrossProduct(S.home, { 0, 1, 0 }));
        const Vector3 dir = Vector3Normalize(Vector3Add(S.home, Vector3Scale(S.side, 8.0f / c.field->radius)));
        S.spot = Vector3Scale(dir, c.field->radius + c.field->elevation(dir) + 0.5f);
        Fluid::Pour(*c.body, S.spot, 300.0f);
        const Vector3 at = Vector3Normalize(Vector3Add(dir, Vector3Scale(S.side, -40.0f / c.field->radius)));
        const Vector3 feet = Vector3Scale(at, c.field->radius + c.field->elevation(at) + 1.0f); // 1 m up: it falls onto the ground
        Player& p = LocalSession().player;
        p.flyMode = false;
        p.placeAt(c.world, *c.body, c.body->toWorld(feet)); // a spawn: the loading screen comes up
        S.placedWall = GetTime();
        S.h0 = Vector3Length(feet);
        return true; // (no c.camera: a scenario camera flies the player with it)
    });
    // Steps don't run while a loading screen is up: the first frame here comes when it has closed.
    s.step("the loading screen closes once the ground is built", 2.0, [](ScenarioCtx& c) {
        const double held = GetTime() - S.placedWall;
        c.check("the water is moving (the world isn't idle)", !Fluid::Asleep());
        c.lt("the spawn's loading screen beside moving water (s; it timed out at 20 s)", held, 10.0);
        return true;
    });
    s.step("the player falls onto the ground (under gravity, not held)", 3.0, [](ScenarioCtx& c) {
        if (c.t < 1.5) return false;
        const Player& p = LocalSession().player;
        c.gt("gravity on the player (m/s^2)", p.gravity, 1.0);
        c.check("the player has landed", p.grounded);
        c.gt("the player came down from 1 m up (m)", S.h0 - Height(c), 0.3);
        return true;
    });
}
