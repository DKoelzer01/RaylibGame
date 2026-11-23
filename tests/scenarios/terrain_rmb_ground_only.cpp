// terrain/rmb-ground-only (#1577, rule #1098): RMB places ground only with a ground item in hand (or a block stack in the off
// hand, #836); bare hands and a tool never place it, in creative too (the playtest: creative placed ground from an empty hand).
// Placing strokes at the ground 2 m ahead (fresh spots), in creative: bare hands, a pickaxe, then dirt in hand (free).
#include "scenario.h"
#include "inventory.h"
#include "playerinput.h"
#include "playersession.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
bool AimAhead(ScenarioCtx& c, float side) { // as terrain/place-dig-back: the camera at the eyes, the ground 2 m ahead
    Player& p = LocalSession().player;
    const Vector3 feet = (p.position - c.world.renderOrigin).toFloat();
    const Vector3 up = Vector3Normalize(p.up);
    Vector3 fwd = Vector3Subtract(p.forward, Vector3Scale(up, Vector3DotProduct(p.forward, up)));
    fwd = Vector3Length(fwd) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(fwd);
    c.world.camera.position = Vector3Add(feet, Vector3Scale(up, 1.6f));
    c.world.camera.target = Vector3Add(Vector3Add(feet, Vector3Scale(fwd, 2.0f)), Vector3Scale(Vector3Normalize(Vector3CrossProduct(fwd, up)), side));
    c.world.camera.up = up;
    TerrainTool_Update(c.world, PlayerInput::None(), true);
    Planetoid* body = nullptr; Vector3 local{}, n{};
    return TerrainTool_LastHit(body, local, n);
}
}

SCENARIO("terrain/rmb-ground-only") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    s.step("creative: bare hands and a tool place no ground; dirt in hand does, for free", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("the player standing on the sandbox");
        if (c.T < 1.0) return false;
        terrainTool.creative = true;
        terrainTool.terraform = false;
        c.check("nothing in the off hand (its block stack would place, #836)", inventory.offHand().id.empty(), inventory.offHand().id);
        struct Case { const char* held; bool places; } cases[] = { { "", false }, { "stone_pickaxe", false }, { "dirt", true } };
        float side = -3.0f;
        for (const Case& k : cases) {
            terrainTool.heldTool = k.held;
            if (k.places) inventory.add("dirt", 5);
            const double dirt0 = inventory.amount("dirt");
            bool aimed = AimAhead(c, side);
            const bool placed = aimed && TerrainTool_Stroke(false);
            const char* what = k.held[0] ? k.held : "bare hands";
            c.check(TextFormat("aimed at the ground (%s)", what), aimed);
            c.check(TextFormat("%s: %s", what, k.places ? "places ground" : "places no ground"), placed == k.places);
            if (k.places) c.near("creative: placing cost nothing", inventory.amount("dirt"), dirt0, 1e-6);
            side += 3.0f;
        }
        terrainTool.heldTool = "";
        return true;
    });
}
