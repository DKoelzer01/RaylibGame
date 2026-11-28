// caves/islands (#1465, the user's rule; #1545 the gate): no cave kind may add floating rock. Every chamber within 800 m of
// seed 42's start is sampled on a 1 m grid (CaveIslands_Measure, as --cave-islands 800): the solid pieces not joined to
// the surrounding rock float. Today (ladder c7305e02, 2026-11-22): 20 chambers, 426 pieces, 722 m^3 (#1116 brought it down
// from 3349). Large caverns (kind 4) within 1200 m are measured too, on a 2 m grid, against their own budget. ~75 s on the
// runner (CPU density).
//
// The acceptance rule (#2025): the chamber and sheet budgets are per room, not totals. A terrain change moves the sample
// (the erosion bake v27, #1884, brought two more chambers within 800 m: 22 chambers, 825 m^3, yet 37.5 m^3 a chamber
// against 36.1 before), so a total failed on rooms the caves didn't change. The budget is today's rock per room plus ~5%:
// a new kind or wall rule that adds floating pieces to the rooms fails here, more rooms in the sample don't.
#include "scenario.h"
#include "caveislands.h"
#include "voxel.h"
#include "world.h"
#include <string>

namespace {
constexpr float kRange = 800.0f;          // m from the start
constexpr double kPerChamber = 38.0;      // m^3 of floating rock a chamber (722 / 20 = 36.1 today; GUESS margin, the user's rule is "none added")
constexpr int kMinChambers = 10;          // fewer: the check would say little (today 20-22)
constexpr float kCavernRange = 1200.0f;   // (#1465) large caverns within this (seed 42: 2, r 86 and 91 m, 2 m grid)
constexpr long long kCavernBudget = 830;  // m^3 floating in them (784 when caverns landed; GUESS margin, the same rule). A total: only two
// (#1465) m^3 floating a sheet cave within the same range: 432 / 36 = 12.0 when they landed; 472 / 36 = 13.1 since the
// bake v27 window (the same 36 sheets; the chambers' jump is the bake, measured: 20 chambers and 722 m^3 just before it;
// the sheets' share is inferred from the same window). Budget 13.1 + ~5% (GUESS).
constexpr double kPerSheet = 13.8;
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
        const double perChamber = r.chambers > 0 ? (double)r.volume / r.chambers : 0.0;
        c.metric("floating rock a chamber (m^3)", perChamber, "m^3");
        c.check("floating rock a chamber within today's budget (<= 38 m^3)", perChamber <= kPerChamber,
                std::to_string(perChamber) + " m^3 (" + std::to_string(r.volume) + " in " + std::to_string(r.chambers) + ")");
        c.metric("caverns within 1200 m", r.caverns, "");
        c.metric("floating rock in caverns (m^3)", (double)r.cavernVolume, "m^3");
        c.check("a cavern to judge (>= 1)", r.caverns >= 1, std::to_string(r.caverns));
        c.check("floating rock in caverns within budget (<= 830 m^3)", r.cavernVolume <= kCavernBudget, std::to_string(r.cavernVolume) + " m^3");
        c.metric("sheet caves within 1200 m", r.sheets, "");
        c.metric("floating rock in sheet caves (m^3)", (double)r.sheetVolume, "m^3");
        const double perSheet = r.sheets > 0 ? (double)r.sheetVolume / r.sheets : 0.0;
        c.metric("floating rock a sheet cave (m^3)", perSheet, "m^3");
        c.check("floating rock a sheet cave within budget (<= 13.8 m^3)", perSheet <= kPerSheet,
                std::to_string(perSheet) + " m^3 (" + std::to_string(r.sheetVolume) + " in " + std::to_string(r.sheets) + ")");
        return true;
    });
}
