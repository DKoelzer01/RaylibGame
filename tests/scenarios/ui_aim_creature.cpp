// ui/aim-creature (#1564): with the crosshair on a creature the terrain tool showed its dig/place circle and text for the
// ground behind it. A hornhide is pinned 4 m in front of the player, who stays on foot (the creature aim is off in flight)
// and looks through the game's own camera: aimed at it, the creature is the target and the terrain tool has none (no
// preview); aimed at the bare ground beside it, the terrain tool has its target again.
#include "scenario.h"
#include "animals.h"
#include "player.h"
#include "terraintool.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>

namespace {
int idx = -1;
Vector3 groundAt(ScenarioCtx& c, Vector3 dir) { return Vector3Scale(dir, c.field->radius + c.field->elevation(dir)); }
struct View { Vector3 eye, creature, beside, up; };
View Frame(ScenarioCtx& c) { // from where the player stands, body-local
    const Vector3 feet = c.body->toLocal((player.position - c.world.renderOrigin).toFloat()), up = Vector3Normalize(feet);
    const Vector3 east = Vector3Normalize(Vector3CrossProduct(up, Vector3{ 0, 0, 1 }));
    const Vector3 north = Vector3CrossProduct(east, up);
    const float R = c.field->radius;
    View v;
    v.up = up;
    v.eye = Vector3Add(feet, Vector3Scale(up, 1.6f));
    v.creature = groundAt(c, Vector3Normalize(Vector3Add(Vector3Scale(up, R), Vector3Scale(east, 4.0f))));
    v.beside = groundAt(c, Vector3Normalize(Vector3Add(Vector3Scale(up, R), Vector3Add(Vector3Scale(east, 4.0f), Vector3Scale(north, 3.0f)))));
    return v;
}
void Look(ScenarioCtx& c, const View& v, Vector3 at) { // the player's own view (the game puts the camera at its eye)
    player.flyMode = false;
    player.forward = c.body->dirToWorld(Vector3Normalize(Vector3Subtract(at, v.eye)));
}
}

SCENARIO("ui/aim-creature") {
    s.tags = { "creatures" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 60.0;
    s.step("settle on foot, a hornhide pinned 4 m ahead", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play, an idle world, on the ground");
        player.flyMode = false;
        if (c.T < 1.0 || !WorldIdle() || !player.grounded) return false;
        const View v = Frame(c);
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "hornhide", v.creature, 1, false, 0.5f);
        idx = Animals_Count() > before ? Animals_Count() - 1 : -1;
        if (!c.check("a hornhide spawned", idx >= 0)) return true;
        Animals_TestPin(idx, v.creature, v.eye);
        return true;
    });
    s.step("aimed at the creature: it is the target, no terrain preview", 5.0, [](ScenarioCtx& c) {
        if (idx < 0) return true;
        const View v = Frame(c);
        Animals_TestPin(idx, v.creature, v.eye);
        Look(c, v, Vector3Add(v.creature, Vector3Scale(v.up, 0.6f)));
        if (c.t < 0.5) return false;
        c.check("the creature is aimed", Animals_Aiming());
        c.check("the terrain tool aims at nothing (no preview, no text)", !TerrainTool_AimedGround(nullptr, nullptr));
        c.grab("on-creature");
        return true;
    });
    s.step("aimed at the ground beside it: the terrain tool's target again", 5.0, [](ScenarioCtx& c) {
        if (idx < 0) return true;
        const View v = Frame(c);
        Animals_TestPin(idx, v.creature, v.eye);
        Look(c, v, v.beside);
        if (c.t < 0.5) return false;
        c.check("no creature aimed", !Animals_Aiming());
        c.check("the terrain tool aims at the ground", TerrainTool_AimedGround(nullptr, nullptr));
        return true;
    });
}
