// creatures/cave-spawn-region (#1740, the underground plan #1465): cave spawns read the underground region map
// (underground.h). Two test rows make the hornhide a cave dweller of the region at seed 42's nearest chamber and the
// woolback one of a different region; the player stands in the chamber for a few spawn passes. Every new animal is checked
// the frame it appears (where it was placed): each hornhide stands in its region, each woolback in its own, and some
// hornhides spawned (the check is not empty).
#include "scenario.h"
#include "animals.h"
#include "creatures.h"
#include "player.h"
#include "survival.h"
#include "underground.h"
#include "voxel.h"
#include "world.h"
#include <raymath.h>
#include <string>
#include <unordered_set>

namespace {
struct St {
    Vector3 chamber{}, out{};
    UndergroundRegion in = UndergroundRegion::None, other = UndergroundRegion::None;
    std::unordered_set<uint64_t> checked; // net ids already checked (the ones there before the setup too)
    int hornIn = 0, hornOut = 0, woolIn = 0, woolOut = 0;
} S;
}

SCENARIO("creatures/cave-spawn-region") {
    s.tags = { "creatures", "caves" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 200.0;
    s.record = false;

    s.step("settle", 120.0, [](ScenarioCtx& c) { c.waitingFor("2 s of play and an idle world"); return c.T >= 2.0 && WorldIdle(); });
    s.step("two region-bound cave dwellers, in the nearest chamber", 5.0, [](ScenarioCtx& c) {
        const Vector3 home = c.field->findSpawnDir(c.body->seed);
        if (!c.check("a chamber within 3 km", c.field->nearestCaveMouth(home, 3000.0f, S.chamber, S.out, 3))) { c.fail("no chamber"); return true; }
        const Vector3 up = Vector3Normalize(S.chamber);
        Vector3 feet = S.chamber;
        for (float d = 0.0f; d < 40.0f; d += 0.25f) { const Vector3 q = Vector3Subtract(S.chamber, Vector3Scale(up, d)); if (c.body->densityAt(q) >= ISO_LEVEL) { feet = Vector3Add(q, Vector3Scale(up, 0.3f)); break; } }
        S.in = Underground_RegionAt(*c.field, feet);
        S.other = (UndergroundRegion)(((int)S.in + 1) % (int)UndergroundRegion::Count);
        c.note(TextFormat("chamber region: %s; the woolback's: %s", UndergroundRegion_Name(S.in), UndergroundRegion_Name(S.other)));
        std::string err;
        const std::string horn = std::string("hornhide | livestock | aggro=passive habitat=cave depth=6-500 dark=0.1 cap=2 respawn=30 region=") + UndergroundRegion_Name(S.in);
        const std::string wool = std::string("woolback | livestock | aggro=passive habitat=cave depth=6-500 dark=0.1 cap=2 respawn=30 region=") + UndergroundRegion_Name(S.other);
        if (!c.check("the test rows parse", Creatures_TestAdd(horn, &err) && Creatures_TestAdd(wool, &err), err)) { c.fail(err); return true; }
        Animals_SetNaturalSpawns(false);
        Animals_SetThreats(true);
        worldDifficulty = 1;
        Animals_ResetCaveSpawns();
        player.flyMode = false;
        player.placeAt(c.world, *c.body, c.body->toWorld(feet));
        for (int i = 0; i < Animals_Count(); ++i) S.checked.insert(Animals_NetId(i)); // the ones already there are not cave spawns
        return true;
    });
    s.step("spawn passes: each in its region", 60.0, [](ScenarioCtx& c) {
        static double idleAt = -1.0;
        // Check each hornhide and woolback the first frame it is there (by its net id: it has not walked yet).
        for (int i = 0; i < Animals_Count(); ++i) {
            if (Animals_Body(i) != c.body->name) continue;
            const std::string sp = Animals_SpeciesOf(i).first;
            if ((sp != "hornhide" && sp != "woolback") || !S.checked.insert(Animals_NetId(i)).second) continue;
            const UndergroundRegion r = Underground_RegionAt(*c.field, Animals_Position(i));
            if (sp == "hornhide") (r == S.in ? S.hornIn : S.hornOut)++;
            else (r == S.other ? S.woolIn : S.woolOut)++;
        }
        c.waitingFor("the chamber built, then 6 s (three passes)");
        if (idleAt < 0.0 && c.t > 1.0 && WorldIdle()) { idleAt = c.t; Animals_CaveSpawnNow(); }
        if (idleAt < 0.0 || c.t < idleAt + 6.0) return false;
        c.metric("spawn sites in the last pass", Animals_CaveSites(), "", false);
        c.metric("hornhides spawned in their region", S.hornIn, "", false);
        c.metric("woolbacks spawned in their region", S.woolIn, "", false);
        c.gt("hornhides spawned (in the chamber's region)", S.hornIn, 0.5);
        c.eq("hornhides outside their region", S.hornOut, 0);
        c.eq("woolbacks outside their region", S.woolOut, 0);
        return true;
    });
}
