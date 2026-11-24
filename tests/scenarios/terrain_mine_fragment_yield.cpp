// terrain/mine-fragment-yield (#1718): mining a floating fragment yields what the strokes take out, once. On the sandbox, a
// fragment of placed gravel floats 3 m ahead of the player (a sphere of 1.2 m radius, nothing under it). A stone shovel digs it
// away: each pass sweeps the aim over the fragment's face on a grid and strokes wherever the aim hits it, until a whole pass
// finds nothing to aim at. The gravel gained must be the fragment's volume (the placed-ground rate, 1.8 an m^3; 60-105 %:
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
// The camera at the player's eyes looking at a point near the fragment; then the tool's aim from it. True: the aim hit
// something within `near` metres of the fragment (render space) - the sandbox's ground is 2 m below and never in the way.
bool AimAt(ScenarioCtx& c, Vector3 targetRender, Vector3 fragmentRender, float near) {
    Player& p = LocalSession().player;
    const Vector3 feet = (p.position - c.world.renderOrigin).toFloat();
    const Vector3 up = Vector3Normalize(p.up);
    c.world.camera.position = Vector3Add(feet, Vector3Scale(up, 1.6f));
    c.world.camera.target = targetRender;
    c.world.camera.up = up;
    TerrainTool_Update(c.world, PlayerInput::None(), true);
    Planetoid* body = nullptr; Vector3 local{}, n{};
    if (!TerrainTool_LastHit(body, local, n)) return false;
    return Vector3Distance(body->toWorld(local), fragmentRender) < near;
}
}

SCENARIO("terrain/mine-fragment-yield") {
    s.tags = { "world", "mining" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    s.step("mine a floating gravel fragment away: yield = its volume, then nothing", 100.0, [](ScenarioCtx& c) {
        c.waitingFor("the player standing on the sandbox, the world built");
        // and an idle world: on a slow start (cold shader compiles, #1769) the ground isn't built after 1 s of play
        if (c.T < 1.0 || !WorldIdle()) return false;
        terrainTool.creative = false;
        terrainTool.heldTool = "stone_shovel";
        Player& pl = LocalSession().player;
        const Vector3 feet = (pl.position - c.world.renderOrigin).toFloat();
        const Vector3 up = Vector3Normalize(pl.up);
        Vector3 fwd = Vector3Subtract(pl.forward, Vector3Scale(up, Vector3DotProduct(pl.forward, up)));
        fwd = Vector3Length(fwd) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(fwd);
        const Vector3 right = Vector3Normalize(Vector3CrossProduct(fwd, up));
        const Vector3 centreRender = Vector3Add(Vector3Add(feet, Vector3Scale(fwd, 3.0f)), Vector3Scale(up, 2.2f)); // 2.2 m up: floating
        const float radius = 1.2f;
        TerrainEdit fill{};
        fill.fill = true; fill.center = c.body->toLocal(centreRender); fill.radius = radius; fill.paint = (int)VoxelMaterial::Gravel;
        ApplyTerrainEdit(*c.body, fill);
        const double expected = 1.8 * (4.0 / 3.0) * PI * radius * radius * radius; // the placed-ground rate a m^3
        const float near = radius + 0.6f;
        // The aim grid over the fragment's face: -1.0 .. 1.0 m across and up, every 0.5 m.
        auto sweep = [&](bool stroke) {
            int hits = 0;
            for (int iy = -2; iy <= 2; ++iy)
                for (int ix = -2; ix <= 2; ++ix) {
                    const Vector3 target = Vector3Add(centreRender, Vector3Add(Vector3Scale(right, 0.5f * ix), Vector3Scale(up, 0.5f * iy)));
                    if (!AimAt(c, target, centreRender, near)) continue;
                    ++hits;
                    if (stroke) TerrainTool_Stroke(true);
                }
            return hits;
        };
        c.check("the fragment is aimable", AimAt(c, centreRender, centreRender, near));

        const double g0 = inventory.amount("gravel");
        int passes = 0, swings = 0;
        for (; passes < 40; ++passes) {
            const int hits = sweep(true);
            swings += hits;
            if (hits == 0) break;
        }
        const double total = inventory.amount("gravel") - g0;
        c.note(TextFormat("%d passes, %d swings: %.2f gravel (fragment volume %.2f m^3 -> %.2f expected)", passes, swings, total, expected / 1.8, expected));
        c.lt("the fragment was mined away in the pass limit", passes, 40);
        c.between("total yield is the fragment's volume (60-105 %: the rim's ground that is not marked as placed digs as natural ground)", total, 0.6 * expected, 1.05 * expected);

        // Nothing left to take: 20 more sweeps' strokes and a yield probe over the old spot give nothing.
        const double g1 = inventory.amount("gravel");
        for (int n = 0; n < 20; ++n) if (AimAt(c, centreRender, centreRender, near)) TerrainTool_Stroke(true);
        c.eq("20 more swings yield nothing", inventory.amount("gravel") - g1, 0.0, 1e-6);
        const DigYield probe = TerrainTool_DigYield(*c.body, fill.center, radius, fill.center, true, 0.0f);
        float left = 0.0f;
        for (const auto& kv : probe.items) left += kv.second;
        for (const auto& kv : probe.placed) left += kv.second;
        c.lt("what a stroke over the old spot could still take is only unaimable slivers", left, 1.0);
        c.lt("yield and what is left together never exceed the fragment", total + left, 1.05 * expected);

        // The aim no longer finds anything there (the ray runs on into the sandbox's open air).
        c.eq("nothing is aimable where the fragment was", sweep(false), 0, 0);
        return true;
    });
}
