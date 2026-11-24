// world/structure-site-stall (#1791, the user: structure generation on a worker, no main-thread stall walking into range):
// on the planet (seed 42), the player is set down 200 m from a site not built yet (inside the 250 m build range), so the
// world placement builds it as in play - its ground cuts (pad, interior, doorways, a dungeon's rooms) and its parts'
// foundations folded on a worker, its parts a few a frame. Checks: the site is built and its job finishes, its strokes are
// in the body's edits, and no frame spends more than a few ms in Structures_Update (it was 150-420 ms for the fold alone).
// Then at a second site, the world is saved while its ground is still folding: the save puts it in at once (nothing lost).
#include "scenario.h"
#include "perf_budget.h"
#include "player.h"
#include "structures.h"
#include "terraindeltas.h"
#include "world.h"
#include <raymath.h>
#include <algorithm>
#include <string>
#include <vector>

namespace {
int sites0 = 0;
size_t strokes0 = 0;
double worstMs = 0.0, worstAt = 0.0;
bool placed = false; // (the last step: the player set down at the second site)
BudgetFrames frames;

// The player 200 m from the which-th site (nearest first) 700-4000 m from here that is not a wonder (wonders shape their
// own ground), along the great circle toward here: inside the build range, so its first BuildAround has the player there
// (as walking in does, a step later). False: no such site.
bool SetDownNear(ScenarioCtx& c, int which) {
    const Vector3 here = c.body->toLocalD(player.position).toFloat();
    std::vector<SiteInfo> near;
    Structures_SitesNear(*c.body, here, 4000.0f, near);
    const SiteInfo* pick = nullptr;
    for (const SiteInfo& si : near)
        if (si.distance > 700.0f && si.kind != "wonder" && which-- == 0) { pick = &si; break; }
    if (!pick) return false;
    c.note(TextFormat("site: %s (%s%s), %.0f m away", pick->set.c_str(), pick->kind.c_str(), pick->dungeon ? ", a dungeon under it" : "", pick->distance));
    const Vector3 sd = Vector3Normalize(pick->at), hd = Vector3Normalize(here);
    const Vector3 perp = Vector3Normalize(Vector3Subtract(hd, Vector3Scale(sd, Vector3DotProduct(hd, sd))));
    const float a = 200.0f / c.field->radius;
    const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(sd, cosf(a)), Vector3Scale(perp, sinf(a))));
    sites0 = Structures_SitesBuilt();
    strokes0 = c.body->edits->strokes;
    player.flyMode = false;
    player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Scale(d, c.field->radius + c.field->elevation(d) + 1.0f)));
    return true;
}
}

SCENARIO("world/structure-site-stall") {
    s.tags = { "world", "perf" };
    s.world = ScenarioWorld::Planet; // (the sandbox's plane is at sea level: no site is built there)
    s.args = { "--day" };
    s.timeout = 300.0;
    s.record = false; // (a perf check: the recorder costs frame time)
    s.step("set the player down 200 m from a site not built yet", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle() || !player.nearest) return false;
        if (!SetDownNear(c, 0)) c.fail("no site 700-4000 m from the start");
        worstMs = 0.0;
        return true;
    });
    s.step("the site is built with no main-thread stall", 90.0, [](ScenarioCtx& c) {
        frames.frame();
        const double ms = Structures_LastUpdateMs();
        if (ms > worstMs) { worstMs = ms; worstAt = c.t; }
        c.waitingFor(TextFormat("the site built (%d -> %d) and its job done (%d pending)", sites0, Structures_SitesBuilt(), Structures_JobsPending()));
        if (c.t < 3.0 || Structures_SitesBuilt() <= sites0 || Structures_JobsPending() > 0 || frames.ms.size() <= 60) return false;
        c.note(TextFormat("built in %.1f s; worst Structures_Update %.2f ms at %.1f s", c.t, worstMs, worstAt));
        c.gt("sites built by walking in", Structures_SitesBuilt() - sites0, 0.0);
        c.gt("site strokes in the body's edits", (double)(c.body->edits->strokes - strokes0), 0.0);
        c.metric("Structures_Update max", worstMs, "ms");
        // (GUESS 10: 0.9-2.2 ms measured on the sandbox, 6.1 once with another job on it; the fold on the main thread was 150-420)
        c.timingLt("Structures_Update max ms", worstMs, 10.0);
        frames.report(c);
        return true;
    });
    s.step("a save while the next site's ground folds puts it in", 90.0, [](ScenarioCtx& c) {
        if (!placed) {
            placed = true;
            if (!SetDownNear(c, 1)) { c.fail("no second site 700-4000 m away"); return true; }
            return false;
        }
        c.waitingFor(TextFormat("the second site's job (%d pending, %d built)", Structures_JobsPending(), Structures_SitesBuilt() - sites0));
        if (Structures_SitesBuilt() <= sites0 || Structures_JobsPending() == 0) return false; // (a job with nothing to fold ends at once: wait for one)
        const bool folding = c.body->edits->strokes == strokes0; // its cuts not in yet
        c.note(folding ? "saving with the site's ground still folding" : "saving with the site's parts still going down");
        c.check("the world saved", SaveWorld(c.world));
        c.eq("build jobs left after the save", Structures_JobsPending(), 0.0);
        c.gt("the site's strokes in the body's edits after the save", (double)(c.body->edits->strokes - strokes0), 0.0);
        return true;
    });
}
