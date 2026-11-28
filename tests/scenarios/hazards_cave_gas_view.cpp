// hazards/cave-gas-view (#1472): the cave gas pockets the seed puts in a real cave (the planet world, seed 42), seen from
// the player's eyes in the dark: a CO2 pocket ahead in a passage (the suit warning, the grey haze sinking against the rock),
// a walk into it (the warning turns to "in CO2"), then a methane pocket walked into with a lit torch (one blast, then the
// pocket burnt off). The checks repeat hazards/cave-gas's on seed pockets (not placed by hand), so the real path is tested
// too; the capture is kept for a person to judge. The light: the suit's headlamp (and the torch), nothing else. The walk
// stands 2 m outside the gas: inside the basic suit's 3 m warning.
#include "scenario.h"
#include "cavegas.h"
#include "cliptelemetry.h"
#include "hotbar.h"
#include "inventory.h"
#include "player.h"
#include "suit.h"
#include "survival.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <cmath>
#include <string>
#include <vector>

namespace {
constexpr float kEye = 1.6f;   // m, feet to eye
// A walk to a pocket: from `from` (eye height, outside the gas) to `to` (eye height, in it), on the passage floor.
struct Walk { bool ok = false; CaveGasPocket p; Vector3 from{}, to{}; };
struct St { Walk co2, ch4; int blasts0 = 0; float health0 = 0.0f; double warnAt = -1.0, gasAt = -1.0, blastAt = -1.0; } S;

Vector3 Up(Vector3 local) { return Vector3Normalize(local); }
bool Air(ScenarioCtx& c, Vector3 q) { return c.body->densityAt(q) < ISO_LEVEL; }
// The floor under a point in the air (body-local feet), within `depth` m; false: none.
bool FloorUnder(ScenarioCtx& c, Vector3 q, float depth, Vector3& feet) {
    const Vector3 down = Vector3Negate(Up(q));
    if (!Air(c, q)) return false;
    for (float d = 0.1f; d <= depth; d += 0.1f) {
        const Vector3 s = Vector3Add(q, Vector3Scale(down, d));
        if (!Air(c, s)) { feet = Vector3Add(q, Vector3Scale(down, d - 0.1f)); return true; }
    }
    return false;
}
// A pocket of this kind the player can walk into: its centre in the air over a floor, the eye there in the gas, and a
// level stretch of passage from 3 m outside the gas to the middle.
bool FindWalk(ScenarioCtx& c, const std::vector<CaveGasPocket>& ps, CaveGasKind kind, Walk& w) {
    for (const CaveGasPocket& p : ps) {
        if (p.kind != kind) continue;
        Vector3 fc;
        if (!FloorUnder(c, p.centre, p.radius, fc)) continue;
        const Vector3 up = Up(fc), eyeIn = Vector3Add(fc, Vector3Scale(up, kEye));
        if (!CaveGas_InPocket(p, eyeIn) || !Air(c, eyeIn)) continue;
        const Vector3 e1 = Vector3Normalize(Vector3CrossProduct(up, fabsf(up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        const Vector3 e2 = Vector3CrossProduct(up, e1);
        for (int k = 0; k < 24; ++k) {
            const float a = 6.2831853f * k / 24.0f;
            const Vector3 d = Vector3Add(Vector3Scale(e1, cosf(a)), Vector3Scale(e2, sinf(a)));
            const Vector3 side = Vector3CrossProduct(up, d);
            float edge = 0.0f; // where this line leaves the gas (m from the middle)
            while (edge < 2.0f * p.radius && CaveGas_InPocket(p, Vector3Add(eyeIn, Vector3Scale(d, edge)))) edge += 0.1f;
            const float out = edge + 2.0f;
            bool clear = true; Vector3 lastFeet = fc, feetOut = fc;
            for (float s = 0.0f; s <= out && clear; s += 0.5f) {
                Vector3 f;
                clear = FloorUnder(c, Vector3Add(eyeIn, Vector3Scale(d, s)), kEye + 1.0f, f) &&
                        fabsf(Vector3DotProduct(Vector3Subtract(f, lastFeet), up)) < 0.4f && // no steps over 0.4 m
                        fabsf(Vector3DotProduct(Vector3Subtract(f, fc), up)) < 0.8f;         // and level: within 0.8 m
                const Vector3 eye = Vector3Add(f, Vector3Scale(up, kEye));
                for (const CaveGasPocket& q : ps) clear = clear && (q.id == p.id || !CaveGas_InPocket(q, eye)); // no other gas on the way
                feetOut = f;
                // Headroom: the eye clear of the walls by more than the near plane's reach (the clip telemetry's 0.15 m).
                clear = clear && Air(c, Vector3Add(f, Vector3Scale(up, kEye + 0.5f))) &&
                        Air(c, Vector3Add(Vector3Add(f, Vector3Scale(up, kEye)), Vector3Scale(side, 0.35f))) &&
                        Air(c, Vector3Add(Vector3Add(f, Vector3Scale(up, kEye)), Vector3Scale(side, -0.35f)));
                lastFeet = f;
            }
            const Vector3 eyeOut = Vector3Add(feetOut, Vector3Scale(up, kEye)); // the eye over the floor where the walk starts
            if (!clear || CaveGas_InPocket(p, eyeOut) || !Air(c, eyeOut)) continue;
            w.ok = true; w.p = p; w.from = eyeOut; w.to = eyeIn;
            return true;
        }
    }
    return false;
}
// Face along the walk (level), toward the pocket's middle. The player walks it by physics (autoWalk: W held), not placed
// each frame: the capsule's own collision keeps it off the passage walls.
void Face(ScenarioCtx& c, const Walk& w) {
    const Vector3 up = c.body->dirToLocal(player.up);
    Vector3 look = Vector3Subtract(w.to, c.body->toLocalD(player.position).toFloat());
    look = Vector3Subtract(look, Vector3Scale(up, Vector3DotProduct(look, up)));
    if (Vector3Length(look) > 0.3f) player.forward = c.body->dirToWorld(Vector3Normalize(look));
    // The suit's headlamp the only light (the basic suit has no battery: power topped up by hand).
    survival.headlamp = true; survival.power = 100.0f;
}
void PlaceAt(ScenarioCtx& c, Vector3 eye) {
    player.flyMode = false;
    player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Subtract(eye, Vector3Scale(Up(eye), kEye))));
}
// Go to a pocket found from the start, let the world there load, then choose the walk again: the terrain far from the start
// changes once it loads (the erosion bake is applied where it is loaded, and the caves follow the surface).
void AddGoTo(ScenarioDef& s, Walk* w, CaveGasKind kind, const std::string& name) {
    s.step("go to the " + name, 60.0, [w](ScenarioCtx& c) {
        if (c.t < 0.05) PlaceAt(c, w->from);
        c.waitingFor("the world loaded there");
        return c.t > 2.0 && WorldIdle();
    });
    s.step("the walk into the " + name + ", on the loaded cave", 60.0, [w, kind, name](ScenarioCtx& c) {
        if (c.t < 0.05) {
            std::vector<CaveGasPocket> ps;
            CaveGas_PocketsNear(*c.field, w->p.centre, 80.0f, ps);
            Walk again;
            if (!c.check("a " + name + " pocket to walk into, there", FindWalk(c, ps, kind, again))) { c.fail("no pocket"); return true; }
            *w = again;
            PlaceAt(c, w->from);
        }
        Face(c, *w);
        survival.oxygen = 100.0f; survival.health = 100.0f;
        c.waitingFor("the cave built round the player (an idle world)");
        return c.t > 2.0 && WorldIdle();
    });
}
}

SCENARIO("hazards/cave-gas-view") {
    s.tags = { "hazards", "survival", "visual" }; // (visual: local only; hazards/cave-gas runs on the runner)
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;

    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("find a CO2 and a methane pocket", 1.0, [](ScenarioCtx& c) {
        CaveGas_Reset();
        inventory.add("battery_pack", 1); // suit power: the basic suit's headlamp needs one (#217)
        c.check("a battery pack fitted", Suit_Fit("battery_pack"));
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        std::vector<CaveGasPocket> ps;
        for (float depth : { 90.0f, 200.0f, 400.0f }) {
            ps.clear();
            CaveGas_PocketsNear(*c.field, Vector3Scale(home, c.field->radius - depth), 1500.0f, ps);
            if (!S.co2.ok) FindWalk(c, ps, CaveGasKind::CO2, S.co2);
            if (!S.ch4.ok) FindWalk(c, ps, CaveGasKind::Methane, S.ch4);
        }
        if (!c.check("a CO2 pocket to walk into", S.co2.ok) | !c.check("a methane pocket to walk into", S.ch4.ok)) { c.fail("no pocket"); return true; }
        c.metric("CO2 pocket radius (m)", S.co2.p.radius, "m");
        c.metric("CO2 pocket depth (m)", c.field->radius - Vector3Length(S.co2.p.centre), "m");
        c.metric("methane pocket radius (m)", S.ch4.p.radius, "m");
        return true;
    });
    AddGoTo(s, &S.co2, CaveGasKind::CO2, "CO2");
    s.step("look: the haze ahead", 3.0, [](ScenarioCtx& c) {
        Face(c, S.co2);
        const CaveGasStats& g = CaveGas_Stats();
        if (g.warnKind == (int)CaveGasKind::CO2 && g.warnDist > 0.0f && S.warnAt < 0.0) S.warnAt = c.T;
        if (c.t < 2.5) return false;
        c.check("the camera out of the rock", !ClipTelemetry_CameraInside());
        c.check("breathable here", !survival.inGas);
        c.gt("pockets near the player", CaveGas_Stats().pockets, 0.0);
        c.gt("haze puffs so far", CaveGas_Stats().puffs, 20.0);
        return true;
    });
    s.step("walk in: the warning first", 20.0, [](ScenarioCtx& c) {
        Face(c, S.co2);
        player.autoWalk = true;
        const CaveGasStats& g = CaveGas_Stats();
        if (g.warnKind == (int)CaveGasKind::CO2 && g.warnDist > 0.0f && S.warnAt < 0.0) S.warnAt = c.T;
        if (survival.caveGas == (uint8_t)CaveGasKind::CO2 && S.gasAt < 0.0) S.gasAt = c.T;
        c.waitingFor("the eye in the CO2 and 2 s in it (" + std::to_string(g.edgeDist) + " m from it)");
        if (S.gasAt < 0.0 || c.T < S.gasAt + 2.0) return false;
        player.autoWalk = false;
        c.check("warned before the gas", S.warnAt >= 0.0);
        c.lt("the warning ahead of the gas (s)", S.warnAt, S.gasAt - 0.3);
        c.check("not breathable in it", survival.inGas && !survival.breathable);
        c.check("the warning says we're in it", g.warnDist == 0.0f);
        return true;
    });
    AddGoTo(s, &S.ch4, CaveGasKind::Methane, "methane");
    s.step("by the methane, a torch in hand", 3.0, [](ScenarioCtx& c) {
        if (c.t < 0.05) { inventory.add("torch", 3); Hotbar_Bind(4, "torch"); Hotbar_Press(4, false); }
        Face(c, S.ch4);
        S.blasts0 = CaveGas_Stats().blasts; S.health0 = survival.health; S.warnAt = -1.0;
        if (c.t <= 2.0) return false;
        c.check("methane: the camera out of the rock", !ClipTelemetry_CameraInside());
        c.check("torch in hand", Hotbar_Slot(Hotbar_Selected()) == "torch");
        return true;
    });
    s.step("walk in with the torch: one blast", 20.0, [](ScenarioCtx& c) {
        Face(c, S.ch4);
        player.autoWalk = S.blastAt < 0.0;
        const CaveGasStats& g = CaveGas_Stats();
        if (g.warnKind == (int)CaveGasKind::Methane && g.warnDist > 0.0f) S.warnAt = c.t;
        if (g.blasts == S.blasts0) { c.waitingFor("the blast (" + std::to_string(g.edgeDist) + " m from the gas)"); return false; }
        if (S.blastAt < 0.0) { S.blastAt = c.t; c.check("warned of methane before the blast", S.warnAt >= 0.0); }
        if (c.t < S.blastAt + 3.0) return false; // the fire and the smoke
        c.near("one blast", CaveGas_Stats().blasts - S.blasts0, 1.0, 0.0);
        c.gt("it hurt (hp)", S.health0 - survival.health, 5.0);
        c.check("no methane left in it", survival.caveGas == 0);
        return true;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) {
        CaveGas_Reset();
        survival.health = 100.0f;
        c.keepCapture("cave gas in a real cave (#1472): for a person to judge");
        return true;
    });
}
