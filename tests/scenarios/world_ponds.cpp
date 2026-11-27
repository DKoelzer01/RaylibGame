// world/ponds (#1965, user 2026-11-27: "ponds don't seem to generate at all anymore, just actual large lakes"): the bake's
// small ponds (#360) through the terrain's own rules (PlanetField::pondCensus -> bakePond, as pondOfCell): how many become
// ponds, and why the rest do not (two in one 260 m cell, by a lake or river, level too low, the bowl's depth). Checked: some
// ponds exist, and one within 5 km of the spawn. Run on seeds 42 and 1025335508 (--seed). Reported: the census, the nearest
// pond and the smallest.
#include "scenario.h"
#include "player.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>

SCENARIO("world/ponds") {
    s.tags = { "world", "water", "terrain" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 120.0;
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
        c.metric("ponds kept", n.kept, "", true);
        c.metric("ponds within 5 km of the spawn", n.within, "", true);
        c.metric("nearest pond", n.nearest, "m", false);
        c.gt("bake ponds", n.bake, 0);
        c.gt("ponds kept", n.kept, 0);
        c.gt("ponds within 5 km of the spawn", n.within, 0);
        return true;
    });
}
