// caves/islands (#1465, the user's rule): no cave kind may add floating rock. Every chamber within 800 m of seed 42's start
// is sampled on a 1 m grid (CaveIslands_Measure, as --cave-islands 800): the solid pieces not joined to the surrounding rock
// float. Today (ladder c7305e02, 2026-11-22): 20 chambers, 426 pieces, 722 m^3 (#1116 brought it down from 3349). The
// budget is that plus ~5%: a new kind or wall rule that adds floating pieces fails here. ~60 s on the runner (CPU density).
#include "scenario.h"
#include "caveislands.h"
#include "voxel.h"
#include "world.h"
#include <string>

namespace {
constexpr float kRange = 800.0f;          // m from the start
constexpr long long kBudget = 760;        // m^3 of floating rock in all (today 722; GUESS margin, the user's rule is "none added")
constexpr int kMinChambers = 10;          // fewer: the check would say little (today 20)
}

SCENARIO("caves/islands") {
    s.tags = { "caves" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 120.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("floating rock in the chambers", 10.0, [](ScenarioCtx& c) {
        if (!c.field || !c.body) { c.fail("no planet"); return true; }
        const CaveIslands r = CaveIslands_Measure(*c.field, c.body->seed, kRange, false);
        c.metric("chambers within 800 m", r.chambers, "");
        c.metric("floating pieces", (double)r.pieces, "");
        c.metric("floating rock (m^3)", (double)r.volume, "m^3");
        c.metric("measure time (s)", r.seconds, "s");
        std::string h;
        for (int b = 0; b < 8; ++b) h += std::to_string(r.byHeight[b]) + (b < 7 ? " " : "");
        c.note("floating m^3 by height from the chamber centre, -28..+28 m in 7 m bins: " + h);
        c.check("enough chambers to judge (>= 10)", r.chambers >= kMinChambers, std::to_string(r.chambers));
        c.check("floating rock within today's budget (<= 760 m^3)", r.volume <= kBudget, std::to_string(r.volume) + " m^3");
        return true;
    });
}
