// farm/till-rmb (#1085 playtest fail): a hoe's RMB tills soil into farmland. #1141's rule "RMB places ground only with a
// ground item in hand" kept the stroke from running with a tool at all, so nothing tilled (--test-farm calls Farm_Till
// directly and missed it). Through the real input path (TerrainTool_Update, RMB held): a stone hoe tills the sandbox's
// natural surface and a patch of placed dirt, its preview the size it tills; a shovel (#1085, the user: it digs, it doesn't
// till) and a pickaxe till and place nothing.
#include "scenario.h"
#include "farming.h"
#include "playerinput.h"
#include "playersession.h"
#include "scene.h"
#include "terraintool.h"
#include "tools.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
// The camera at the player's eyes, looking at the ground 2 m ahead (side: metres to the right); the tool's aim from it.
// (#1176) From a pose taken once, body-local: the player turns to the camera's view, so building each aim from the
// player's facing turned it further with every sideways look - the stroke and the check 3.8-8.3 m apart.
bool Aim(ScenarioCtx& c, float side, Vector3& local) {
    static bool taken = false;
    static Vector3 feetL{}, upL{}, fwdL{};
    Player& p = LocalSession().player;
    Planetoid* b = c.body;
    if (!taken) {
        taken = true;
        feetL = b->toLocal((p.position - c.world.renderOrigin).toFloat());
        upL = Vector3Normalize(feetL);
        const Vector3 f = b->dirToLocal(p.forward);
        fwdL = Vector3Subtract(f, Vector3Scale(upL, Vector3DotProduct(f, upL)));
        fwdL = Vector3Length(fwdL) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(upL, { 0, 0, 1 })) : Vector3Normalize(fwdL);
    }
    const Vector3 rightL = Vector3Normalize(Vector3CrossProduct(fwdL, upL));
    c.world.camera.position = b->toWorld(Vector3Add(feetL, Vector3Scale(upL, 1.6f)));
    c.world.camera.target = b->toWorld(Vector3Add(Vector3Add(feetL, Vector3Scale(fwdL, 2.0f)), Vector3Scale(rightL, side)));
    c.world.camera.up = b->dirToWorld(upL);
    TerrainTool_Update(c.world, PlayerInput::None(), true);
    Planetoid* body = nullptr; Vector3 n{};
    return TerrainTool_LastHit(body, local, n);
}
// RMB held for one update with `held` in hand (the caller waited out the cooldown): what the player's RMB does there.
void Rmb(ScenarioCtx& c, const char* held) {
    // In hand as a player has it (#1176): carried (the terrain tool uses a tool only while it is in the inventory) and
    // nothing else held (holdingItem, left on by the hotbar's start selection, blocks every RMB). Setting heldTool alone
    // worked only while the sandbox's starting kit happened to carry these tools.
    if (inventory.count(held) <= 0) inventory.add(held, 1);
    terrainTool.holdingItem = false;
    terrainTool.heldTool = held;
    PlayerInput in; in.hold(Act::Secondary);
    TerrainTool_Update(c.world, in, true);
}
}

SCENARIO("farm/till-rmb") {
    s.tags = { "world", "farm" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    static int stage = 0;
    static float side = -6.0f;
    static int strokes0 = 0;
    static double dirt0 = 0.0;
    static Vector3 spot{};
    struct Case { const char* tool; bool placedDirt; bool tills; };
    static const Case cases[] = { { "stone_hoe", false, true }, { "stone_shovel", false, false }, { "stone_hoe", true, true },
                                  { "stone_shovel", true, false }, { "stone_pickaxe", false, false } };
    s.step("RMB with a hoe tills (its preview the size it tills); with a shovel or a pickaxe nothing", 40.0, [](ScenarioCtx& c) {
        if (c.T < 1.0) return false;
        terrainTool.creative = false;
        const int i = stage / 3, phase = stage % 3;
        if (i >= (int)(sizeof(cases) / sizeof(cases[0]))) return true;
        const Case& k = cases[i];
        if (TerrainTool_Cooldown() > 0.0f) return false; // one stroke per cooldown: wait it out in sim time
        if (phase == 0) { // a fresh spot (3 m apart); with placedDirt, a patch of dirt laid there first
            if (!Aim(c, side, spot)) { c.fail("no ground ahead"); return true; }
            c.check(TextFormat("case %d: not farmland before", i), !Farm_Aims(*c.body, spot));
            if (k.placedDirt) { inventory.add("dirt", 20); Rmb(c, "dirt"); }
            ++stage;
            return false;
        }
        if (phase == 1) {
            if (inventory.count(k.tool) <= 0) inventory.add(k.tool, 1); // (in hand for the aim: the preview is the held tool's)
            terrainTool.holdingItem = false;
            terrainTool.heldTool = k.tool;
            Aim(c, side, spot);
            if (k.tills) // (#1085, the user: the hoe's preview showed the dig's size, not the patch it tills)
                c.near(TextFormat("case %d: the preview is the till's size (m)", i), gToolSphere.radius, TillRadius(ToolOf(k.tool)), 1e-4);
            strokes0 = terrainTool.swings;
            dirt0 = inventory.amount("dirt");
            Rmb(c, k.tool);
            ++stage;
            return false;
        }
        Aim(c, side, spot);
        const std::string what = std::string(k.tool) + (k.placedDirt ? " on placed dirt" : " on the natural ground");
        if (k.tills) {
            c.check(("RMB tills: " + what).c_str(), Farm_Aims(*c.body, spot), TextFormat("strokes %d -> %d", strokes0, terrainTool.swings));
        } else {
            c.check(("RMB tills nothing: " + what).c_str(), !Farm_Aims(*c.body, spot));
            c.eq(("RMB strikes nothing: " + what).c_str(), terrainTool.swings - strokes0, 0);
        }
        c.eq(("no ground placed or paid: " + what).c_str(), inventory.amount("dirt"), dirt0);
        side += 3.0f;
        ++stage;
        return false;
    });
}
