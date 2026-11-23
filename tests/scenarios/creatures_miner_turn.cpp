// creatures/miner-turn (#1717, playtest: hunting Lost Miners "turning slowly one way, flipping 180, turning back, flipping").
// A Lost Miner 4.5 m from the player with a brick wall (5 wide, 3 high) between them, provoked so it hunts; three cases: the
// player behind it, to its left, to its right. Close in (the turn-rate rule's range) a blocked hunter used to detour or turn
// away by up to ~180 degrees at once, then turn back at 180 deg/s into the same wall and flip again. Checked per case over 3 s:
// no turn of more than 45 degrees in one update, and it ends facing the player (within 35 degrees).
#include "scenario.h"
#include "animals.h"
#include "player.h"
#include "plots.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>

namespace {
struct St { uint64_t net = 0; int caseNo = 0; float lastHeading = 0.0f; float worstJump = 0.0f; bool have = false; } S;
Plot* plot = nullptr;
Int3 base{};
constexpr const char* kCase[3] = { "behind it", "to its left", "to its right" };
constexpr float kOffset[3] = { PI, 0.5f * PI, -0.5f * PI }; // its heading minus the bearing to the player

Vector3 Feet(ScenarioCtx& c) { return c.body->toLocalD(player.position).toFloat(); }
int Idx() { return S.net ? Animals_FindNet(S.net) : -1; }
float Wrap(float a) { return atan2f(sinf(a), cosf(a)); }
// The bearing from the miner to the player in its local tangent frame (the heading's convention: from east toward north).
float BearingToPlayer(ScenarioCtx& c) {
    const Vector3 at = Animals_Position(Idx()), up = Vector3Normalize(at);
    const Vector3 cr = Vector3CrossProduct({ 0, 1, 0 }, up); // animals.cpp TangentFrame
    const Vector3 east = Vector3Length(cr) < 1e-3f ? Vector3{ 1, 0, 0 } : Vector3Normalize(cr), north = Vector3Normalize(Vector3CrossProduct(up, east));
    const Vector3 to = Vector3Subtract(Feet(c), at);
    return atan2f(Vector3DotProduct(to, north), Vector3DotProduct(to, east));
}
} // namespace

SCENARIO("creatures/miner-turn") {
    s.tags = { "creatures", "combat" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;

    s.step("settle, a wall and a miner behind it", 20.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        Animals_SetNaturalSpawns(false);
        Animals_SetThreats(true);
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        const Vector3 feet = Feet(c);
        if (!c.body->plots) c.body->plots = std::make_unique<PlotSet>();
        plot = c.body->plots->plotFor(feet, 0.0f, true);
        base = plot->cellAt(plot->toPlot(feet));
        const int mat = std::max(0, BuildMaterialForItem("stone_brick"));
        int built = 0;
        for (int z = -2; z <= 2; ++z)
            for (int y = 0; y < 3; ++y) built += plot->place(Part{ PartType::Cube, { base.x + 2, base.y + y, base.z + z }, 0, (uint8_t)mat }) >= 0;
        c.eq("wall cells placed", built, 15);
        const Vector3 minerAt = plot->toBody(Vector3{ (float)base.x + 4.5f, (float)base.y + 0.5f, (float)base.z + 0.5f });
        const int before = Animals_Count();
        Animals_Spawn(c.world, *c.body, "lost_miner", minerAt, 1, false, 0.01f);
        if (!c.check("a lost miner", Animals_Count() > before)) { c.fail("no miner"); return true; }
        S.net = Animals_NetId(Animals_Count() - 1);
        return true;
    });
    for (int k = 0; k < 3; ++k)
        s.step(std::string("the player ") + kCase[k] + ": it turns to face them, no flip", 6.0, [](ScenarioCtx& c) {
            survival.health = 100.0f; survival.dead = false;
            const int i = Idx();
            if (i < 0) { c.check("the miner is there", false); c.fail("gone"); return true; }
            if (S.caseNo == 0 && c.t < 1.0) { Animals_TestProvokeAll(); return false; } // the first case: awake and hunting first
            static double startAt = -1.0;
            if (startAt < 0.0 || c.t < startAt) { // the case's start: its facing set so the player is behind / left / right of it
                startAt = c.t;
                S.have = false; S.worstJump = 0.0f;
                Animals_TestSetHeading(i, Wrap(BearingToPlayer(c) + kOffset[S.caseNo]));
                Animals_TestProvokeAll();
                return false;
            }
            Animals_TestProvokeAll(); // kept on it (hurt: it hunts the player)
            const float h = Animals_Heading(i);
            if (S.have) S.worstJump = std::max(S.worstJump, fabsf(Wrap(h - S.lastHeading)));
            S.lastHeading = h; S.have = true;
            if (c.t - startAt < 3.0) return false;
            const float err = fabsf(Wrap(h - BearingToPlayer(c)));
            c.lt(std::string("largest turn in one update, deg (") + kCase[S.caseNo] + ")", S.worstJump * RAD2DEG, 45.0);
            c.lt(std::string("facing error at the end, deg (") + kCase[S.caseNo] + ")", err * RAD2DEG, 35.0);
            ++S.caseNo;
            return true;
        });
}
