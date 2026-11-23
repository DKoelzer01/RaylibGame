// terrain/mine-fragment-yield (#1718): mining a floating fragment yields what the strokes take out, once. On the sandbox, a
// fragment of placed gravel floats 3 m ahead of the player (a sphere of 1.2 m radius, nothing under it). Bare hands dig it
// away stroke by stroke. The gravel gained must match the fragment's volume (the placed-ground rate, 1.8 an m^3, within 25 %:
// the yield grid samples the strokes' spheres), a further 20 swings and a yield probe over the old spot must give nothing
// (a swing that removes nothing pays nothing), and the aim must no longer hit anything there.
#include "scenario.h"
#include "playerinput.h"
#include "playersession.h"
#include "terraintool.h"
#include "inventory.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
// The camera at the player's eyes looking at the fragment's centre; then the tool's aim from it.
bool AimAtFragment(ScenarioCtx& c, Vector3 fragmentRender) {
    Player& p = LocalSession().player;
    const Vector3 feet = (p.position - c.world.renderOrigin).toFloat();
    const Vector3 up = Vector3Normalize(p.up);
    c.world.camera.position = Vector3Add(feet, Vector3Scale(up, 1.6f));
    c.world.camera.target = fragmentRender;
    c.world.camera.up = up;
    TerrainTool_Update(c.world, PlayerInput::None(), true);
    Planetoid* body = nullptr; Vector3 local{}, n{};
    return TerrainTool_LastHit(body, local, n);
}
}

SCENARIO("terrain/mine-fragment-yield") {
    s.tags = { "world", "mining" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 90.0;
    s.step("mine a floating gravel fragment away: yield = its volume, then nothing", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("the player standing on the sandbox");
        if (c.T < 1.0) return false;
        terrainTool.creative = false;
        terrainTool.heldTool = ""; // bare hands: soft ground only
        Player& pl = LocalSession().player;
        const Vector3 feet = (pl.position - c.world.renderOrigin).toFloat();
        const Vector3 up = Vector3Normalize(pl.up);
        Vector3 fwd = Vector3Subtract(pl.forward, Vector3Scale(up, Vector3DotProduct(pl.forward, up)));
        fwd = Vector3Length(fwd) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(fwd);
        const Vector3 centreRender = Vector3Add(Vector3Add(feet, Vector3Scale(fwd, 3.0f)), Vector3Scale(up, 2.2f)); // 2.2 m up: floating
        const float radius = 1.2f;
        TerrainEdit fill{};
        fill.fill = true; fill.center = c.body->toLocal(centreRender); fill.radius = radius; fill.paint = (int)VoxelMaterial::Gravel;
        ApplyTerrainEdit(*c.body, fill);
        const double expected = 1.8 * (4.0 / 3.0) * PI * radius * radius * radius; // the placed-ground rate a m^3
        c.check("the fragment is aimable", AimAtFragment(c, centreRender));

        const double g0 = inventory.amount("gravel");
        int swings = 0;
        for (int n = 0; n < 60 && AimAtFragment(c, centreRender); ++n) if (TerrainTool_Stroke(true)) ++swings;
        const double total = inventory.amount("gravel") - g0;
        c.note(TextFormat("%d swings took the fragment: %.2f gravel (fragment volume %.2f m^3 -> %.2f expected)", swings, total, expected / 1.8, expected));
        c.near("total yield matches the fragment's volume", total, expected, 0.25 * expected);

        // Nothing left to take: 20 more swings (wherever the aim lands) and a yield probe over the old spot give nothing.
        const double g1 = inventory.amount("gravel");
        for (int n = 0; n < 20; ++n) if (AimAtFragment(c, centreRender)) TerrainTool_Stroke(true);
        c.eq("20 more swings yield nothing", inventory.amount("gravel") - g1, 0.0, 1e-6);
        const DigYield probe = TerrainTool_DigYield(*c.body, fill.center, radius, fill.center, true, 0.0f);
        float left = 0.0f;
        for (const auto& kv : probe.items) left += kv.second;
        for (const auto& kv : probe.placed) left += kv.second;
        c.eq("a stroke over the old spot would yield nothing", left, 0.0, 1e-6);

        // The aim no longer finds anything there (the ray runs on into the sandbox's open air).
        c.check("nothing is aimable where the fragment was", !AimAtFragment(c, centreRender));
        return true;
    });
}
