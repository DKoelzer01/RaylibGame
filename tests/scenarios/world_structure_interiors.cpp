// world/structure-interiors (#1637, the user: generated structures kept terrain inside their rooms): every set that stands on
// the ground (ruins, huts, camps, towers...; not the wonders, which shape their own ground) is built forced on the side of a
// hill raised beside it on the sandbox, so the hill reaches into its rooms. Every cell inside its walls with no part in it
// must be air (rooms, corridors, doorways); the outside may meet the hill.
#include "scenario.h"
#include "player.h"
#include "structures.h"
#include "world.h"
#include <raymath.h>

SCENARIO("world/structure-interiors") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--day" };
    s.timeout = 90.0;
    s.step("every ground set built on a hillside: its interior is air", 60.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 here = c.body->toLocal((player.position - c.world.renderOrigin).toFloat());
        const std::vector<StructureInterior> rs = Structures_TestInteriors(c.world, *c.body, here);
        c.gt("ground sets tested", (double)rs.size(), 0.0);
        for (const StructureInterior& r : rs) {
            if (r.cells < 0) { c.check(TextFormat("%s: built", r.set.c_str()), false); continue; }
            c.note(TextFormat("%s: %d interior cells, %d solid", r.set.c_str(), r.cells, r.solid));
            c.eq(TextFormat("%s: interior cells that are solid ground", r.set.c_str()), r.solid, 0);
        }
        c.grab("structures");
        return true;
    });
}
