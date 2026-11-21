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
Vector3 upL{}, faceL{};                         // body-local: up, and the way the player faces (the wall)
double startH = 0.0, startHealth = 0.0, t0 = 0.0, h0 = 0.0, pushFrom = 0.0;
double Height(ScenarioCtx& c) { return Vector3DotProduct(c.body->toLocalD(player.position).toFloat(), upL) - startH; }
double Off(ScenarioCtx& c) { // how far from the wall the player stands (horizontal, body-local)
    const Vector3 p = c.body->toLocalD(player.position).toFloat(), w = plot->toBody({ 1.0f * kCellSize, 0, 0.5f * kCellSize });
    Vector3 d = Vector3Subtract(p, w);
    d = Vector3Subtract(d, Vector3Scale(upL, Vector3DotProduct(d, upL)));
    return Vector3Length(d);
}
void FollowCam(ScenarioCtx& c) { // from 7 m out and 1 m up, at the player
    const Vector3 feet = c.body->toLocalD(player.position).toFloat();
    const Vector3 eye = Vector3Add(Vector3Subtract(feet, Vector3Scale(faceL, 7.0f)), Vector3Scale(upL, 1.0f));
    c.world.camera.position = c.body->toWorld(eye);
    c.world.camera.target = c.body->toWorld(Vector3Add(feet, Vector3Scale(upL, 0.9f)));
    c.world.camera.up = c.body->dirToWorld(upL);
}
void Hold(uint64_t keys) { gFakeHeld = (gFakeHeld & ~(kW | kS | kJump)) | keys; }
}

SCENARIO("world/ladder") {
    s.tags = { "world", "build" };
    s.world = ScenarioWorld::Sandbox;
    s.timeout = 120.0;
    s.step("build a 10 m wall with ladders, stand facing it", 10.0, [](ScenarioCtx& c) {
        if (c.T < 0.5) return false;
        // The plot 3 m from the start along the ground; the wall 2 cells thick at x 0-1, the ladders at x 2 (their frame on
        // the cell's low-x face: against the wall), the player at x 3.5 facing -x.
        const Vector3 dir = c.spawnDirLocal();
        const Vector3 east = Vector3Normalize(Vector3CrossProduct({ 0, 1, 0 }, dir));
        const Vector3 at = Vector3Normalize(Vector3Add(Vector3Scale(dir, c.field->radius), Vector3Scale(east, 3.0f)));
        const Vector3 ground = Vector3Scale(at, c.field->radius + c.field->elevation(at));
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(ground, 0.0f, true);
        const int mat = std::max(0, BuildMaterialForItem("stone_brick"));
        int built = 0, ladders = 0;
        for (int y = 0; y < kWallCells; ++y) {
            for (int x = 0; x < 2; ++x) built += plot->place(Part{ PartType::Cube, { x, y, 0 }, 0, (uint8_t)mat }) >= 0;
            ladders += plot->place(Part{ PartType::Ladder, { 2, y, 0 }, 0, 0 }) >= 0;
        }
        c.eq("wall cubes placed", built, 2 * kWallCells);
        c.eq("ladders placed", ladders, kWallCells);
        upL = Vector3Normalize(Vector3Subtract(plot->toBody({ 0, 1, 0 }), plot->toBody({ 0, 0, 0 })));
        faceL = Vector3Normalize(Vector3Subtract(plot->toBody({ 0, 0, 0 }), plot->toBody({ 1, 0, 0 }))); // -x: the wall
        const Vector3 stand = plot->toBody({ 3.5f * kCellSize, 0.05f, 0.5f * kCellSize });
        player.flyMode = false;
        player.placeAt(c.world, *c.body, c.body->toWorld(stand));
        player.forward = c.body->dirToWorld(faceL);
        startH = Vector3DotProduct(plot->toBody({ 0, 0, 0 }), upL);
        startHealth = survival.health;
        gScenarioInput = true;
        return true;
    });
    s.step("settle on the ground", 5.0, [](ScenarioCtx& c) { FollowCam(c); return c.t > 1.0 && player.grounded; });
    s.step("climb 1.5 m and push off (Space)", 8.0, [](ScenarioCtx& c) {
        FollowCam(c);
        if (Height(c) < 1.5) { Hold(kW); return false; }
        if (pushFrom == 0.0) { pushFrom = Off(c); Hold(kJump); t0 = c.t; return false; }
        Hold(0);
        if (c.t - t0 < 0.6) return false;
        c.check("pushed off: no longer on the ladder", !player.onLadder);
        c.gt("pushed off: away from the wall (m)", Off(c) - pushFrom, 0.3);
        return true;
    });
    s.step("land, then walk back to the ladder", 8.0, [](ScenarioCtx& c) {
        FollowCam(c);
        if (!player.grounded) return false;
        Hold(kW);
        return player.onLadder;
    });
    s.step("climb to 4 m with W (~2 m/s)", 8.0, [](ScenarioCtx& c) {
        FollowCam(c);
        static bool started = false;
        if (!started) { started = true; t0 = c.t; h0 = Height(c); }
        Hold(kW);
        if (Height(c) < 4.0) return false;
        c.between("climb speed (m/s)", (Height(c) - h0) / std::max(0.1, c.t - t0), 1.5, 2.5);
        Hold(0); t0 = c.t; h0 = Height(c);
        return true;
    });
    s.step("no input holds on the ladder", 3.0, [](ScenarioCtx& c) {
        FollowCam(c);
        Hold(0);
        if (c.t - t0 < 1.0) return false;
        c.lt("height change while holding (m)", fabs(Height(c) - h0), 0.1);
        c.check("still on the ladder", player.onLadder);
        return true;
    });
    s.step("climb down with S to the ground, unhurt", 10.0, [](ScenarioCtx& c) {
        FollowCam(c);
        Hold(kS);
        if (!(player.grounded && Height(c) < 0.3)) return false;
        Hold(0);
        c.eq("health after climbing down", survival.health, startHealth, 0.01);
        return true;
    });
    s.step("climb to the top and step off onto the wall", 14.0, [](ScenarioCtx& c) {
        FollowCam(c);
        Hold(kW);
        if (!(player.grounded && !player.onLadder && Height(c) > kWallCells * kCellSize - 0.2)) return false;
        Hold(0);
        c.gt("standing on the wall's top (m above the ground)", Height(c), kWallCells * kCellSize - 0.2);
        c.eq("health after the whole climb", survival.health, startHealth, 0.01);
        gScenarioInput = false;
        return true;
    });
}
