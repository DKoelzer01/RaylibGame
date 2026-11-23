// caves/islands (#1465, the user's rule): no cave kind may add floating rock. Every chamber within 800 m of seed 42's start
// is sampled on a 1 m grid (CaveIslands_Measure, as --cave-islands 800): the solid pieces not joined to the surrounding rock
// float. Today (ladder c7305e02, 2026-11-22): 20 chambers, 426 pieces, 722 m^3 (#1116 brought it down from 3349). The
// budget is that plus ~5%: a new kind or wall rule that adds floating pieces fails here. Large caverns (kind 4) within 1200 m
// are measured too, on a 2 m grid, against their own budget. ~75 s on the runner (CPU density).
#include "scenario.h"
#include "caveislands.h"
#include "voxel.h"
#include "world.h"
#include <string>

namespace {
constexpr float kRange = 800.0f;          // m from the start
constexpr long long kBudget = 760;        // m^3 of floating rock in all (today 722; GUESS margin, the user's rule is "none added")
constexpr int kMinChambers = 10;          // fewer: the check would say little (today 20)
constexpr float kCavernRange = 1200.0f;   // (#1465) large caverns within this (seed 42: 2, r 86 and 91 m, 2 m grid)
constexpr long long kCavernBudget = 830;  // m^3 floating in them (784 when caverns landed; GUESS margin, the same rule)
constexpr long long kSheetBudget = 460;   // (#1465) m^3 floating in the sheet caves within the same range (432 when they landed: 36 sheets; GUESS)
}

SCENARIO("caves/islands") {
    s.tags = { "caves" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 120.0;
    s.step("settle", 60.0, [](ScenarioCtx& c) { c.waitingFor("an idle world"); return c.T >= 1.0 && WorldIdle(); });
    s.step("floating rock in the chambers", 10.0, [](ScenarioCtx& c) {
        if (!c.field || !c.body) { c.fail("no planet"); return true; }
        const CaveIslands r = CaveIslands_Measure(*c.field, c.body->seed, kRange, false, kCavernRange);
        c.metric("chambers within 800 m", r.chambers, "");
        c.metric("floating pieces", (double)r.pieces, "");
        c.metric("floating rock (m^3)", (double)r.volume, "m^3");
        c.metric("measure time (s)", r.seconds, "s");
        std::string h;
        for (int b = 0; b < 8; ++b) h += std::to_string(r.byHeight[b]) + (b < 7 ? " " : "");
        c.note("floating m^3 by height from the chamber centre, -28..+28 m in 7 m bins: " + h);
        c.check("enough chambers to judge (>= 10)", r.chambers >= kMinChambers, std::to_string(r.chambers));
        c.check("floating rock within today's budget (<= 760 m^3)", r.volume <= kBudget, std::to_string(r.volume) + " m^3");
        c.metric("caverns within 1200 m", r.caverns, "");
        c.metric("floating rock in caverns (m^3)", (double)r.cavernVolume, "m^3");
        c.check("a cavern to judge (>= 1)", r.caverns >= 1, std::to_string(r.caverns));
        c.check("floating rock in caverns within budget (<= 830 m^3)", r.cavernVolume <= kCavernBudget, std::to_string(r.cavernVolume) + " m^3");
        c.metric("sheet caves within 1200 m", r.sheets, "");
        c.metric("floating rock in sheet caves (m^3)", (double)r.sheetVolume, "m^3");
        c.check("floating rock in sheet caves within budget (<= 460 m^3)", r.sheetVolume <= kSheetBudget, std::to_string(r.sheetVolume) + " m^3");
        return true;
    });
}
