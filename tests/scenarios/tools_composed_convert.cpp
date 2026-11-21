// tools/composed-convert (#1449): old per-material tools in a save become composed tools on load, and wear in their own record.
// On the sandbox: the player holds a stone pickaxe, axe and hoe, the pickaxe 20 strokes worn (the old per-id wear); the world
// is saved and loaded back. Then: the inventory holds tool_pickaxe, tool_axe and tool_hoe with stone heads, the pickaxe's
// record carries w=20 and the old per-id wear is gone; six dig strokes with it in hand raise its own wear to 26, and its wear
// bar reads below full. (Chopping and tilling read the same ToolOf hook; trees aren't on the sandbox.)
#include "scenario.h"
#include "composed.h"
#include "inventory.h"
#include "itemdata.h"
#include "playerinput.h"
#include "playersession.h"
#include "terraintool.h"
#include "tools.h"
#include "toolwear.h"
#include "hotbar.h"
#include "skills.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct State { int pick = -1; } S;

bool AimDown(ScenarioCtx& c) { // the eye over the feet, looking at the ground 2 m ahead; then the tool's aim
    Player& p = LocalSession().player;
    const Vector3 feet = (p.position - c.world.renderOrigin).toFloat();
    const Vector3 up = Vector3Normalize(p.up);
    Vector3 fwd = Vector3Subtract(p.forward, Vector3Scale(up, Vector3DotProduct(p.forward, up)));
    fwd = Vector3Length(fwd) < 1e-3f ? Vector3Normalize(Vector3CrossProduct(up, { 0, 0, 1 })) : Vector3Normalize(fwd);
    c.world.camera.position = Vector3Add(feet, Vector3Scale(up, 1.6f));
    c.world.camera.target = Vector3Add(feet, Vector3Scale(fwd, 2.0f));
    c.world.camera.up = up;
    TerrainTool_Update(c.world, PlayerInput::None(), true);
    Planetoid* body = nullptr; Vector3 local{}, n{};
    return TerrainTool_LastHit(body, local, n);
}
int CellOf(const char* id) { for (int i = 0; i < Inventory::kCells; ++i) if (inventory.grid()[i].id == id) return i; return -1; }
}

SCENARIO("tools/composed-convert") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    s.step("old tools, saved and loaded back", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("the player standing on the sandbox");
        if (c.T < 1.0) return false;
        terrainTool.creative = false;
        inventory.add("stone_pickaxe", 1); inventory.add("stone_axe", 1); inventory.add("stone_hoe", 1);
        LocalSession().toolWear.used["stone_pickaxe"] = 20;
        const bool saved = SaveWorld(c.world);
        const bool loaded = saved && LoadWorldSave(c.world);
        c.check("saved and loaded", loaded);
        return true;
    });
    s.step("the old tools are composed ones now", 10.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        c.eq("old stone pickaxes left", inventory.count("stone_pickaxe"), 0);
        c.eq("composed pickaxes", inventory.countInstances("tool_pickaxe"), 1);
        c.eq("composed axes", inventory.countInstances("tool_axe"), 1);
        c.eq("composed hoes", inventory.countInstances("tool_hoe"), 1);
        S.pick = CellOf("tool_pickaxe");
        if (!c.check("a composed pickaxe held", S.pick >= 0)) { c.fail("no pickaxe"); return true; }
        const std::string d = inventory.grid()[S.pick].data;
        c.check("its head is stone", ItemData::Get(d, "h") == "stone", d);
        c.eq("its record carries the old wear (w)", ItemData::GetInt(d, "w", 0), 20);
        c.eq("the old per-id wear went", (int)LocalSession().toolWear.used.count("stone_pickaxe"), 0);
        return true;
    });
    s.step("dig with it: its own wear rises", 20.0, [](ScenarioCtx& c) {
        if (S.pick < 0) return true;
        inventory.preferred = S.pick;
        terrainTool.heldTool = "tool_pickaxe";
        c.eq("ToolOf the held composed pickaxe: tier", ToolOf("tool_pickaxe").tier, 1);
        int strokes = 0;
        for (int n = 0; n < 6; ++n) if (AimDown(c) && TerrainTool_Stroke(true)) ++strokes; // (each stroke counts a swing: the main loop wears the held tool by them)
        c.eq("dig strokes made", strokes, 6);
        return true;
    });
    s.step("the main loop wore the held copy", 5.0, [](ScenarioCtx& c) {
        if (S.pick < 0) return true;
        if (c.t < 0.2) return false; // the frame's wear pass (main.cpp) has run
        c.note(TextFormat("held '%s', hotbar selected %d, preferred %d, swings %d, pick_wear stat %.2f, creative %d", terrainTool.heldTool.c_str(),
                          Hotbar_Selected(), inventory.preferred, terrainTool.swings, Skills_Stat(LocalSession().skills, "pick_wear"), terrainTool.creative ? 1 : 0));
        const int cell = CellOf("tool_pickaxe");
        c.check("the pickaxe is still held", cell >= 0);
        if (cell >= 0) c.eq("its wear after six strokes", ItemData::GetInt(inventory.grid()[cell].data, "w", 0), 26);
        c.lt("its wear bar (left)", ToolWear_Left("tool_pickaxe"), 1.0);
        c.keepCapture("composed tools after a load (#1449): the held pickaxe digging");
        return true;
    });
}
