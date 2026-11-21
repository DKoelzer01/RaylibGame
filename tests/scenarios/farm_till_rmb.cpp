// farm/till-rmb (#1085 playtest fail): a hoe's or shovel's RMB tills soil into farmland. #1141's rule "RMB places ground
// only with a ground item in hand" kept the stroke from running with a tool at all, so nothing tilled (--test-farm calls
// Farm_Till directly and missed it). Through the real input path (TerrainTool_Update, RMB held): a stone hoe and a stone
// shovel each till the sandbox's natural surface and a patch of placed dirt; a pickaxe (the control) tills and places
// nothing.
#include "scenario.h"
#include "farming.h"
#include "playerinput.h"
#include "playersession.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
// The camera at the player's eyes, looking at the ground 2 m ahead (side: metres to the right); the tool's aim from it.
bool Aim(ScenarioCtx& c, float side, Vector3& local) {
    Player& p = LocalSession().player;
    const Vector3 feet = (p.position - c.world.renderOrigin).toFloat();
    const Vector3 up = Vector3Normalize(p.up);
    Vector3 fwd = Vector3Subtract(p.forward, Vector3Scale(up, Vector3DotProduct(p.forward, up)));
    fwd = Vector3Length(fwd) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(fwd);
    c.world.camera.position = Vector3Add(feet, Vector3Scale(up, 1.6f));
    c.world.camera.target = Vector3Add(Vector3Add(feet, Vector3Scale(fwd, 2.0f)), Vector3Scale(Vector3Normalize(Vector3CrossProduct(fwd, up)), side));
    c.world.camera.up = up;
    TerrainTool_Update(c.world, PlayerInput::None(), true);
    Planetoid* body = nullptr; Vector3 n{};
    return TerrainTool_LastHit(body, local, n);
}
// The ground under the aimed point is farmland (the tilled state, #1085): its material a hand's breadth down. (Farm_Aims
// also says yes near a crop field: the wrong probe here.)
bool Farmland(ScenarioCtx& c, Vector3 spot) {
    return c.body->materialAt(Vector3Subtract(spot, Vector3Scale(Vector3Normalize(spot), 0.12f))) == VoxelMaterial::Farmland;
}
// RMB held for one update with `held` in hand (the caller waited out the cooldown): what the player's RMB does there.
void Rmb(ScenarioCtx& c, const char* held) {
    if (inventory.count(held) <= 0) inventory.add(held, 1); // (the tool in hand must be carried: HeldTool reads the inventory)
    terrainTool.heldTool = held;
    terrainTool.holdingItem = false; // (as the hotbar's take() does when a tool or ground item is selected: the start's hotbar
                                     // slot held food, which keeps RMB from striking at all)
    c.note(TextFormat("RMB with %s: holdingItem %d rmbBlocks %d cooldown %.2f", held, terrainTool.holdingItem ? 1 : 0, terrainTool.rmbBlocks ? 1 : 0, TerrainTool_Cooldown()));
    PlayerInput in; in.hold(Act::Secondary);
    TerrainTool_Update(c.world, in, true);
}
}

SCENARIO("farm/till-rmb") {
    s.tags = { "world", "farm" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    static int stage = 0;
    static float side = -12.0f; static double rmbAt = 0.0;
    static int strokes0 = 0;
    static double dirt0 = 0.0;
    static Vector3 spot{};
    struct Case { const char* tool; bool placedDirt; bool tills; };
    static const Case cases[] = { { "stone_hoe", false, true }, { "stone_shovel", false, true }, { "stone_hoe", true, true },
                                  { "stone_shovel", true, true }, { "stone_pickaxe", false, false } };
    s.step("RMB with a hoe or shovel tills; with a pickaxe nothing", 40.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = false;
        const int i = stage / 3, phase = stage % 3;
        if (i >= (int)(sizeof(cases) / sizeof(cases[0]))) return true;
        const Case& k = cases[i];
        if (TerrainTool_Cooldown() > 0.0f) return false; // one stroke per cooldown: wait it out in sim time
        if (phase == 0) { // a fresh spot (3 m apart); with placedDirt, a patch of dirt laid there first
            if (!Aim(c, side, spot)) { c.fail("no ground ahead"); return true; }
            c.check(TextFormat("case %d: not farmland before", i), !Farmland(c, spot));
            if (k.placedDirt) { inventory.add("dirt", 20); Rmb(c, "dirt"); }
            ++stage;
            return false;
        }
        if (phase == 1) {
            Aim(c, side, spot);
            strokes0 = terrainTool.swings;
            dirt0 = inventory.amount("dirt");
            Rmb(c, k.tool); rmbAt = c.T;
            ++stage;
            return false;
        }
        if (c.T - rmbAt < 0.5) return false; // (the till's ground edit lands over the next frames)
        Aim(c, side, spot);
        const std::string what = std::string(k.tool) + (k.placedDirt ? " on placed dirt" : " on the natural ground");
        if (k.tills) {
            c.check(("RMB tills: " + what).c_str(), Farmland(c, spot), TextFormat("strokes %d -> %d", strokes0, terrainTool.swings));
        } else {
            c.check("RMB with a pickaxe tills nothing (the control)", !Farmland(c, spot));
            c.eq("RMB with a pickaxe strikes nothing", terrainTool.swings - strokes0, 0);
        }
        c.eq(("no ground placed or paid: " + what).c_str(), inventory.amount("dirt"), dirt0);
        side += 6.0f; // (6 m apart: a till brush and the farm aim reach past 3 m)
        ++stage;
        return false;
    });
}
