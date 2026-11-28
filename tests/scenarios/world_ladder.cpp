// world/ladder (#1105, the user: "ladders for traversing vertically"): a 10 m stone-brick wall with a column of ladders on
// its face, in the sandbox. The player, facing it: climbs 1.5 m and pushes off (Space); climbs with W (~2 m/s) and holds
// with no input; climbs down with S to the ground, unhurt; climbs to the top and steps off onto the wall's top.
#include "scenario.h"
#include "player.h"
#include "plots.h"
#include "survival.h"
#include "playerinput.h"
#include "world.h"
#include <raymath.h>

namespace {
constexpr int kWallCells = 20;                 // 10 m of wall (0.5 m cells)
constexpr uint64_t kW = 1ull << (int)Act::Forward, kS = 1ull << (int)Act::Back, kJump = 1ull << (int)Act::Jump;
Plot* plot = nullptr;
Vector3 upL{}, faceL{}, wallAt{};               // body-local: up, the way the player faces (the wall), the wall's face
double startH = 0.0, startHealth = 0.0, t0 = 0.0, h0 = 0.0, pushFrom = 0.0;
double Height(ScenarioCtx& c) { return Vector3DotProduct(c.body->toLocalD(player.position).toFloat(), upL) - startH; }
double Off(ScenarioCtx& c) { // how far from the wall the player stands (horizontal, body-local)
    const Vector3 p = c.body->toLocalD(player.position).toFloat(), w = wallAt;
    Vector3 d = Vector3Subtract(p, w);
    d = Vector3Subtract(d, Vector3Scale(upL, Vector3DotProduct(d, upL)));
    return Vector3Length(d);
}
void FollowCam(ScenarioCtx& c) { // at eye height 2.5 m behind the player, on the ladder at their height (as they see it)
    const Vector3 feet = c.body->toLocalD(player.position).toFloat();
    const Vector3 eye = Vector3Add(Vector3Subtract(feet, Vector3Scale(faceL, 2.5f)), Vector3Scale(upL, 1.6f));
    c.world.camera.position = c.body->toWorld(eye);
    c.world.camera.target = c.body->toWorld(Vector3Add(Vector3Add(feet, Vector3Scale(faceL, 0.6f)), Vector3Scale(upL, 1.2f)));
    c.world.camera.up = c.body->dirToWorld(upL);
}
void Trace(ScenarioCtx& c) { // once a second of each step: where the climb stands (a failure's story in results.json)
    static int lastStep = -1, lastSecond = -1;
    if (c.stepIndex == lastStep && (int)c.t == lastSecond) return;
    lastStep = c.stepIndex; lastSecond = (int)c.t;
    const Vector3 surface = c.body->surfaceVelocity(c.body->toWorld(c.body->toLocalD(player.position).toFloat()));
    c.note(TextFormat("step %d t=%d s: height %.2f m, on ladder %d, grounded %d, climb input %.0f, off the wall %.2f m, vertical speed %.2f m/s",
                      c.stepIndex, lastSecond, Height(c), (int)player.onLadder, (int)player.grounded, player.climbInput, Off(c),
                      Vector3DotProduct(Vector3Subtract(player.velocity, surface), player.up)));
}
void FaceWall(ScenarioCtx& c) { // the player looks at the wall (horizontally) whatever turned them
    Vector3 d = Vector3Subtract(wallAt, c.body->toLocalD(player.position).toFloat());
    d = Vector3Subtract(d, Vector3Scale(upL, Vector3DotProduct(d, upL)));
    if (Vector3Length(d) > 0.05f) player.forward = c.body->dirToWorld(Vector3Normalize(d));
}
void Hold(uint64_t keys) { gFakeHeld = (gFakeHeld & ~(kW | kS | kJump)) | keys; }
}

SCENARIO("world/ladder") {
    s.tags = { "world", "build" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    s.args = { "--day" }; // (lit, for its recording)
    s.step("build a 10 m wall with ladders, stand facing it", 20.0, [](ScenarioCtx& c) {
        player.flyMode = false; // (a scenario starts the player flying over the start: down to the ground first)
        if (c.T < 1.0 || !player.grounded) return false;
        // Around where the player stands (moving them re-anchored them at the start): the plot under their feet; the
        // ladders one cell along the plot's -x from theirs (their frame on the cell's low-x face: against the wall), the
        // wall 2 cells thick beyond; the player faces -x.
        const Vector3 feet = c.body->toLocalD(player.position).toFloat();
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        const Int3 me = plot->cellAt(plot->toPlot(feet));
        const int mat = std::max(0, BuildMaterialForItem("stone_brick"));
        int built = 0, ladders = 0;
        for (int y = 0; y < kWallCells; ++y) {
            for (int x = me.x - 3; x <= me.x - 2; ++x) built += plot->place(Part{ PartType::Cube, { x, me.y + y, me.z }, 0, (uint8_t)mat }) >= 0;
            ladders += plot->place(Part{ PartType::Ladder, { me.x - 1, me.y + y, me.z }, 0, 0 }) >= 0;
        }
        c.eq("wall cubes placed", built, 2 * kWallCells);
        c.eq("ladders placed", ladders, kWallCells);
        upL = Vector3Normalize(Vector3Subtract(plot->toBody({ 0, 1, 0 }), plot->toBody({ 0, 0, 0 })));
        faceL = Vector3Normalize(Vector3Subtract(plot->toBody({ 0, 0, 0 }), plot->toBody({ 1, 0, 0 }))); // -x: the wall
        wallAt = plot->toBody({ (me.x - 1) * kCellSize, me.y * kCellSize, (me.z + 0.5f) * kCellSize }); // the wall's face
        player.forward = c.body->dirToWorld(faceL);
        startH = Vector3DotProduct(feet, upL);
        startHealth = survival.health;
        gScenarioInput = true;
        gScenarioCamera = true; // (the view from behind the player, for the recording)
        return true;
    });
    s.step("settle on the ground", 5.0, [](ScenarioCtx& c) { FollowCam(c); return c.t > 1.0 && player.grounded; });
    s.step("climb 1.5 m and push off (Space)", 8.0, [](ScenarioCtx& c) {
        FollowCam(c); Trace(c); FaceWall(c);
        if (Height(c) < 1.5) { Hold(kW); return false; }
        if (pushFrom == 0.0) { pushFrom = Off(c); Hold(kJump); t0 = c.t; return false; }
        Hold(0);
        if (c.t - t0 < 0.6) return false;
        c.check("pushed off: no longer on the ladder", !player.onLadder);
        c.gt("pushed off: away from the wall (m)", Off(c) - pushFrom, 0.3);
        return true;
    });
    s.step("land, then walk back to the ladder", 8.0, [](ScenarioCtx& c) {
        FollowCam(c); Trace(c); FaceWall(c);
        static bool landed = false;
        landed = landed || player.grounded;
        if (!landed) return false; // (after the push-off: down first; then W walks back and up the ladder)
        Hold(kW);
        return player.onLadder;
    });
    s.step("climb to 4 m with W (~2 m/s)", 8.0, [](ScenarioCtx& c) {
        FollowCam(c); Trace(c); FaceWall(c);
        static bool started = false;
        if (!started) { started = true; t0 = c.t; h0 = Height(c); }
        Hold(kW);
        if (Height(c) > 2.5) c.grab("climbing"); // (#1071) a slideshow picture, when asked for
        if (Height(c) < 4.0) return false;
        c.between("climb speed (m/s)", (Height(c) - h0) / std::max(0.1, c.t - t0), 1.5, 2.5);
        Hold(0); t0 = c.t; h0 = Height(c);
        return true;
    });
    s.step("no input holds on the ladder", 3.0, [](ScenarioCtx& c) {
        FollowCam(c); Trace(c); FaceWall(c);
        Hold(0);
        // (#2086) This step's own start: c.t restarts each step, so the previous step's t0 is on another clock.
        static bool started = false;
        if (!started) { started = true; t0 = c.t; h0 = Height(c); }
        if (c.t - t0 < 1.0) return false;
        c.lt("height change while holding (m)", fabs(Height(c) - h0), 0.1);
        c.check("still on the ladder", player.onLadder);
        return true;
    });
    s.step("climb down with S to the ground, unhurt", 10.0, [](ScenarioCtx& c) {
        FollowCam(c); Trace(c); FaceWall(c);
        static bool started = false; static double offMax = 0.0;
        if (!started) { started = true; t0 = c.t; h0 = Height(c); }
        Hold(kS);
        offMax = std::max(offMax, Off(c));
        if (!(player.grounded && Height(c) < 0.3)) return false;
        Hold(0);
        c.between("climb-down speed (m/s): on the ladder all the way, not dropped", (h0 - Height(c)) / std::max(0.1, c.t - t0), 1.5, 2.5);
        c.lt("stayed at the ladder going down (m from the wall at most)", offMax, 0.8);
        c.eq("health after climbing down", survival.health, startHealth, 0.01);
        return true;
    });
    s.step("climb to the top and step off onto the wall", 14.0, [](ScenarioCtx& c) {
        FollowCam(c); Trace(c); FaceWall(c);
        Hold(kW);
        if (!(player.grounded && !player.onLadder && Height(c) > kWallCells * kCellSize - 0.2)) return false;
        Hold(0);
        c.gt("standing on the wall's top (m above the ground)", Height(c), kWallCells * kCellSize - 0.2);
        c.eq("health after the whole climb", survival.health, startHealth, 0.01);
        gScenarioInput = false; gScenarioCamera = false;
        c.keepCapture("the ladder and the climb (a new part the user sees)");
        return true;
    });
}
