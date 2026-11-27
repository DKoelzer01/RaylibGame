// world/ponds (#1965, user 2026-11-27: "ponds don't seem to generate at all anymore, just actual large lakes"): the bake's
// small ponds (#360) through the terrain's own rules (PlanetField::pondCensus -> bakePond, as pondOfCell): how many become
// ponds, and why the rest do not (two in one 260 m cell, by a lake or river, level too low, the bowl's depth), and how many
// show water by any rule (most stand in the river network as small lakes). Checked: some ponds exist, and one within 5 km
// of the spawn. Then the camera looks down on the nearest hollow showing water (the capture is kept for a person to judge).
// Run on seeds 42 and 1025335508 (--seed).
#include "scenario.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>

namespace {
struct State {
    bool wet = false;
    Vector3 at{}, up{}, side{};
} S;
}

SCENARIO("world/ponds") {
    s.tags = { "world", "water", "terrain" };
    s.world = ScenarioWorld::Planet;
    s.args = { "--day" };
    s.timeout = 180.0;
    s.step("count the ponds", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const PlanetField* pf = c.field;
        if (!pf) { c.fail("no planet field"); return true; }
        const Vector3 spawn = Vector3Normalize(c.body->toLocalD(player.position).toFloat());
        const PlanetField::PondCensus n = pf->pondCensus(spawn, 5000.0f);
        c.note(TextFormat("bake ponds %d: kept %d, shadowed (same 260 m cell) %d, by a lake or river %d, level too low %d, bowl depth %d; "
                          "kept within 5 km of the spawn %d, nearest %.0f m, smallest radius %.1f m",
                          n.bake, n.kept, n.shadowed, n.nearWater, n.low, n.depth, n.within, n.nearest, n.smallest));
        c.note(TextFormat("by water: in a lake's water %d (at the pond's own level %d, the largest of those lakes %d cells), on a lake's bank %d, by a river %d; "
                          "bake hollows with water at their middle (any rule) %d of %d, the nearest %.0f m from the spawn",
                          n.inLake, n.lakeIsPond, n.lakeCellsMax, n.lakeBank, n.river, n.wet, n.bake, n.nearestWet));
        c.metric("ponds kept", n.kept, "", true);
        c.metric("ponds within 5 km of the spawn", n.within, "", true);
        c.metric("nearest pond", n.nearest, "m", false);
        c.metric("bake hollows showing water", n.wet, "", true);
        c.gt("bake ponds", n.bake, 0);
        c.gt("ponds kept", n.kept, 0);
        c.gt("ponds within 5 km of the spawn", n.within, 0);
        S.wet = n.wet > 0;
        if (S.wet) {
            S.up = n.nearestWetDir;
            S.at = Vector3Scale(S.up, pf->radius + pf->elevation(S.up, 0.0f));
            S.side = Vector3Normalize(Vector3CrossProduct(S.up, fabsf(S.up.y) < 0.9f ? Vector3{ 0, 1, 0 } : Vector3{ 1, 0, 0 }));
        }
        return true;
    });
    s.step("look down on the nearest wet hollow", 90.0, [](ScenarioCtx& c) {
        if (!S.wet) return true;
        c.camera(Vector3Add(S.at, Vector3Add(Vector3Scale(S.up, 70.0f), Vector3Scale(S.side, 60.0f))), S.at, S.up);
        c.waitingFor("the terrain there built (an idle world), then 2 s");
        static double idleAt = -1.0;
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) { idleAt = c.t; c.cue(); }
        if (idleAt >= 0.0 && c.t > idleAt + 1.5) c.grab("above");
        return idleAt >= 0.0 && c.t > idleAt + 2.0;
    });
    s.step("done", 1.0, [](ScenarioCtx& c) {
        c.keepCapture("the nearest bake hollow showing water, from above (#1965): for a person to judge");
        return true;
    });
}
