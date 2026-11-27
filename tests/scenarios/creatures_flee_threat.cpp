// creatures/flee-threat (#1850, bug hunt 2026-11-26): an animal fleeing a predator ran from the player instead - every step
// the threat it fled was replaced by the player's position - so with the player on the far side it ran into the predator.
// In the sandbox a woolback stands 30 m from the player, a threat point 10 m beyond it (away from the player); it is made
// to flee that point. It must open the gap to the threat (it runs toward the player's side), not close it.
#include "scenario.h"
#include "animals.h"
#include "player.h"
#include "world.h"
#include <raymath.h>

namespace {
uint64_t net = 0;
Vector3 threat{};
float gap0 = 0.0f;

Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocalD(player.position).toFloat(); }
Vector3 Ground(ScenarioCtx& c, Vector3 near, Vector3 dirFlat, float d) { // d metres from `near` along the ground
    const Vector3 dir = Vector3Normalize(Vector3Add(near, Vector3Scale(dirFlat, d)));
    return Vector3Scale(dir, c.field->radius + c.field->elevation(dir));
}
}

SCENARIO("creatures/flee-threat") {
    s.tags = { "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("a woolback 30 m off, a threat 10 m beyond it", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false);
        player.flyMode = false;
        const Vector3 feet = FeetLocal(c), up = Vector3Normalize(feet);
        const Vector3 flat = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.z) < 0.9f ? Vector3{ 0, 0, 1 } : Vector3{ 1, 0, 0 }));
        const Vector3 at = Ground(c, feet, flat, 30.0f);
        threat = Ground(c, feet, flat, 40.0f);
        Animals_Spawn(c.world, *c.body, "woolback", at, 1, false, 0.01f);
        const int i = Animals_Nearest(at, 3.0f);
        if (i < 0) { c.fail("no woolback spawned"); return true; }
        net = Animals_NetId(i);
        Animals_TestFlee(i, threat, 3.0f);
        gap0 = Vector3Distance(Animals_Position(i), threat);
        return true;
    });
    s.step("it runs from the threat, not from the player", 10.0, [](ScenarioCtx& c) {
        const int i = Animals_FindNet(net);
        if (i < 0) { c.fail("the woolback is gone"); return true; }
        c.waitingFor("2 s of fleeing");
        if (c.t < 2.0) return false;
        const float gap = Vector3Distance(Animals_Position(i), threat);
        c.note(TextFormat("gap to the threat %.1f -> %.1f m", gap0, gap));
        c.gt("metres gained from the threat in 2 s", gap - gap0, 3.0);
        return true;
    });
}
