// Cave gas pockets (#1472, the underground plan #1465).
//   hazards/cave-gas (sandbox, pockets placed by hand ahead of the player): walking into a CO2 pocket, the suit warns and the
//     hiss starts before the air turns bad; inside, the air is not breathable and the oxygen drains 3x. Walking into methane
//     with a lit torch in hand sets it off once: a blast that hurts, then the pocket is burnt off (no gas, no warning).
#include "scenario.h"
#include "cavegas.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "suit.h"
#include "survival.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
struct St {
    Vector3 start{}, up{}, ahead{};
    double warnAt = -1.0, hissAt = -1.0, gasAt = -1.0;
    float warnDist = -1.0f, oxygenIn = 0.0f, health0 = 0.0f;
    int blasts0 = 0;
    bool harmedBeforeGas = false;
} S;

Vector3 FeetLocal(ScenarioCtx& c) { return c.body->toLocal((player.position - c.world.renderOrigin).toFloat()); }
// Stand `d` m along the walk from the start (on the sandbox's flat ground). Not placeAt: each one counts as a spawn, and
// the loading wait after a spawn pauses the survival step and the gas.
void StandAt(ScenarioCtx& c, float d) {
    const Vector3 feet = c.body->toWorld(Vector3Add(S.start, Vector3Scale(S.ahead, d)));
    player.position = c.world.renderOrigin + DVec3(feet);
    player.velocity = c.body->surfaceVelocity(feet);
}
// A pocket of this kind `d` m along the walk, its centre `lift` m over the start's feet.
void AddPocket(CaveGasKind kind, float d, float lift, float radius) {
    CaveGasPocket p;
    p.kind = kind; p.radius = radius; p.up = S.up;
    p.centre = Vector3Add(Vector3Add(S.start, Vector3Scale(S.ahead, d)), Vector3Scale(S.up, lift));
    CaveGas_TestAdd(p);
}
}

SCENARIO("hazards/cave-gas") {
    s.tags = { "caves", "hazards", "survival" }; // (caves: the runner's headless part 2)
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" }; // (the haze shows in the capture)
    s.timeout = 90.0;

    s.step("settle", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    // (#2109) The settle step's spawn raises the loading screen, which holds the gas update until the terrain round the
    // spawn is built: on the hosted null-GL loop that takes longer than the walk's deadline in sim time. So the walk
    // starts only once the update sees the pocket (12 m away, inside the 40 m near range).
    s.step("a CO2 pocket 12 m ahead", 30.0, [](ScenarioCtx& c) {
        if (c.t < 0.5) return false;
        auto place = [&] {
            CaveGas_Reset();
            S.start = FeetLocal(c); S.up = Vector3Normalize(c.body->dirToLocal(player.up)); // (the sandbox's ground isn't radial)
            const Vector3 f = c.body->dirToLocal(player.forward);
            S.ahead = Vector3Normalize(Vector3Subtract(f, Vector3Scale(S.up, Vector3DotProduct(f, S.up))));
            // r 6 m, its centre 1 m over the ground: the gas reaches 2.8 m up, over the head (1.6 m).
            AddPocket(CaveGasKind::CO2, 12.0f, 1.0f, 6.0f);
        };
        static bool added = false;
        if (!added) { added = true; place(); }
        c.waitingFor("the gas update to see the pocket (the spawn's loading screen holds it)");
        if (CaveGas_Stats().pockets == 0) return false;
        place(); // again from where the player stands now the world runs (the spawn may have settled them meanwhile)
        survival.oxygen = 100.0f; survival.health = 100.0f;
        S.health0 = survival.health;
        c.check("clean air at the start", !survival.inGas && survival.breathable);
        return true;
    });
    s.step("walk in: the warning and the hiss come first", 12.0, [](ScenarioCtx& c) {
        StandAt(c, std::min(7.0f, (float)c.t * 0.8f)); // 0.8 m/s, up to 1 m inside the edge (12 - 6 m)
        const CaveGasStats& g = CaveGas_Stats();
        if (g.warnKind == (int)CaveGasKind::CO2 && S.warnAt < 0.0) { S.warnAt = c.t; S.warnDist = g.warnDist; }
        if (g.hissing && S.hissAt < 0.0) S.hissAt = c.t;
        if (survival.health < S.health0 - 0.01f && S.gasAt < 0.0) S.harmedBeforeGas = true;
        if (survival.caveGas == (uint8_t)CaveGasKind::CO2 && S.gasAt < 0.0) { S.gasAt = c.t; S.oxygenIn = survival.oxygen; }
        if (S.gasAt < 0.0 || c.t < S.gasAt + 0.1) {
            c.waitingFor("the head in the CO2: " + std::to_string(g.pockets) + " pockets near, the head " +
                         std::to_string(g.edgeDist) + " m from the gas, warning " + std::to_string(g.warnKind));
            return false;
        }
        c.check("warned", S.warnAt >= 0.0);
        c.check("the hiss heard", S.hissAt >= 0.0);
        c.lt("the warning before the gas (s)", S.warnAt, S.gasAt - 0.5);
        c.lt("the hiss before the gas (s)", S.hissAt, S.gasAt - 0.5);
        c.lt("warned within 3 m of the edge (m)", S.warnDist, 3.01);
        c.check("no harm before the gas", !S.harmedBeforeGas);
        c.check("not breathable in it", survival.inGas && !survival.breathable);
        return true;
    });
    s.step("inside: oxygen drains 3x", 4.0, [](ScenarioCtx& c) {
        StandAt(c, 9.0f); // 3 m into the pocket
        static double t0 = -1.0; static float o0 = 0.0f;
        if (c.t < 0.2) return false;
        if (t0 < 0.0) { t0 = c.t; o0 = survival.oxygen; }
        if (c.t < t0 + 2.0) return false;
        const float rate = (o0 - survival.oxygen) / (float)(c.t - t0);
        c.near("oxygen drain (per s): 3 x 0.35 / the suit's capacity", rate, 3.0f * survival.oxygenDrainPerSec / Suit_Stats().oxygen, 0.15f * 3.0f * survival.oxygenDrainPerSec);
        c.check("the warning says we're in it", CaveGas_Stats().warnKind == (int)CaveGasKind::CO2 && CaveGas_Stats().warnDist == 0.0f);
        return true;
    });
    s.step("back out; methane ahead and a torch in hand", 4.0, [](ScenarioCtx& c) {
        StandAt(c, 0.0f);
        if (c.t < 0.5) return false;
        CaveGas_Reset();
        survival.oxygen = 100.0f; survival.health = 100.0f;
        // Methane rises: its centre 2 m up, r 5 m (down to 0.5 m over the ground, the head in it).
        AddPocket(CaveGasKind::Methane, 10.0f, 2.0f, 5.0f);
        inventory.add("torch", 3);
        Hotbar_Bind(4, "torch");
        Hotbar_Press(4, false);
        S.blasts0 = CaveGas_Stats().blasts; S.health0 = survival.health;
        c.check("torch in hand", Hotbar_Slot(Hotbar_Selected()) == "torch" && inventory.count("torch") > 0);
        return true;
    });
    s.step("walk in with the torch: one blast", 12.0, [](ScenarioCtx& c) {
        StandAt(c, std::min(7.0f, (float)c.t * 0.8f));
        static double blastAt = -1.0;
        if (CaveGas_Stats().blasts == S.blasts0) { c.waitingFor("the blast"); return false; }
        if (blastAt < 0.0) blastAt = c.t;
        if (c.t < blastAt + 0.3) return false; // (the next survival steps take it)
        c.near("one blast", CaveGas_Stats().blasts - S.blasts0, 1.0, 0.0);
        c.gt("it hurt (hp)", S.health0 - survival.health, 10.0);
        return true;
    });
    s.step("burnt off: no gas, no second blast", 4.0, [](ScenarioCtx& c) {
        StandAt(c, 10.0f); // the pocket's middle
        if (c.t < 2.0) return false;
        c.near("still one blast", CaveGas_Stats().blasts - S.blasts0, 1.0, 0.0);
        c.check("no methane left in it", survival.caveGas == 0);
        c.check("no warning", CaveGas_Stats().warnKind == 0);
        survival.health = 100.0f;
        CaveGas_Reset();
        return true;
    });
}
