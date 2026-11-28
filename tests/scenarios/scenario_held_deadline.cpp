// harness/held-deadline (#2119): a step that waits while the world builds doesn't spend its deadline on the build. The
// headless sim runs many times real time while chunk rebuilds take wall time, so on a loaded runner a sim deadline ran out
// before the world was idle again (11 of the 20 flaky scenarios of #2107). Here: dig six strokes on the sandbox, then wait
// for the rebuild with a deadline of 0.05 s sim, a few frames: before #2119 that step failed on any machine.
#include "scenario.h"
#include "playerinput.h"
#include "playersession.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>

namespace {
// The camera at the player's eyes, looking at the ground 2 m ahead and `side` metres to the right; true when it aims at ground.
bool AimAt(ScenarioCtx& c, float side) {
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
int gStrokes = 0;
}

SCENARIO("harness/held-deadline") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 60.0;
    s.step("settle: the start built", 40.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        return c.T >= 1.0 && WorldIdle();
    });
    s.step("dig six strokes", 5.0, [](ScenarioCtx& c) {
        terrainTool.creative = true;
        terrainTool.heldTool = "";
        gStrokes = 0;
        for (float side = -5.0f; side <= 5.0f; side += 2.0f)
            if (AimAt(c, side) && TerrainTool_Stroke(true)) ++gStrokes;
        c.eq("strokes dug", gStrokes, 6);
        return true;
    });
    // The rebuild takes frames of wall time; only the frames after the world is idle again count toward 0.05 s.
    s.step("wait for the rebuild on a 0.05 s deadline", 0.05, [](ScenarioCtx& c) {
        c.waitingFor("the dug chunks rebuilt (an idle world)");
        if (!WorldIdle()) return false;
        c.note(TextFormat("rebuilt after %.2f s sim, %.2f s wall; the deadline counted %.2f s of it", c.t, c.wall, c.clock));
        c.check("the world was busy after the digs (the step had something to wait for)", c.t > c.clock + 0.01,
                TextFormat("%.3f s sim, %.3f s counted", c.t, c.clock));
        c.lt("the deadline counted only the idle frames (s)", c.clock, 0.05);
        return true;
    });
}
