// world/structure-interiors (#1637, the user: generated structures kept terrain inside their rooms): every set that stands on
// the ground (ruins, huts, camps, towers...; not the wonders, which shape their own ground) is built forced on the side of a
// hill raised beside it (seed 42, round the start), so the hill reaches into its rooms. Every cell inside its walls with no part in it
// must be air (rooms, corridors, doorways); the outside may meet the hill.
#include "scenario.h"
#include "player.h"
#include "structures.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>

SCENARIO("world/structure-interiors") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Planet; // (the sandbox's plane is at sea level: no site is built there)
    s.args = { "--day" };
    s.timeout = 300.0;
    s.step("every ground set built on a hillside: its interior is air", 200.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        const Vector3 home = c.spawnDirLocal(); // (round the start: the player may not be placed yet on a slow machine, #1555)
        const Vector3 here = Vector3Scale(home, c.field->radius + c.field->elevation(home));
        const std::vector<StructureInterior> rs = Structures_TestInteriors(c.world, *c.body, here);
        int built = 0;
        for (const StructureInterior& r : rs) {
            if (r.cells < 0) { c.note(TextFormat("%s: not built here (%s)", r.set.c_str(), r.skip.c_str())); continue; } // (water at its spot)
            ++built;
            c.note(TextFormat("%s: %d interior cells, %d solid", r.set.c_str(), r.cells, r.solid));
            c.eq(TextFormat("%s: interior cells that are solid ground", r.set.c_str()), r.solid, 0);
        }
        c.gt("sets built and tested (of the ground sets)", built, (double)rs.size() * 0.6);
        c.grab("structures");
        return true;
    });
}
