// world/structure-site-stall (#1791, the user: structure generation on a worker, no main-thread stall walking into range):
// on the planet (seed 42), the player is set down 200 m from a site not built yet (inside the 250 m build range), so the
// world placement builds it as in play - its ground cuts (pad, interior, doorways, a dungeon's rooms) folded on a worker,
// then its parts a few a frame. Checks: the site is built and its job finishes (ground and parts in), its strokes are in
// the body's edits, and no frame spends more than a few ms in Structures_Update (it was 150-420 ms for the fold alone).
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
std::string siteName;
BudgetFrames frames;
}

SCENARIO("world/structure-site-stall") {
    s.tags = { "world", "perf" };
    s.world = ScenarioWorld::Planet; // (the sandbox's plane is at sea level: no site is built there)
    s.args = { "--day" };
    s.timeout = 240.0;
    s.record = false; // (a perf check: the recorder costs frame time)
    s.step("set the player down 200 m from a site not built yet", 120.0, [](ScenarioCtx& c) {
        c.waitingFor("1 s of play and an idle world");
        if (c.T < 1.0 || !WorldIdle() || !player.nearest) return false;
        const Vector3 here = c.body->toLocalD(player.position).toFloat();
        std::vector<SiteInfo> near;
        Structures_SitesNear(*c.body, here, 4000.0f, near);
        const SiteInfo* pick = nullptr;
        for (const SiteInfo& si : near) // the nearest ground set past the start's build range (a wonder shapes its own ground)
            if (si.distance > 700.0f && si.kind != "wonder") { pick = &si; break; }
        if (!pick) { c.fail("no site 700-4000 m from the start"); return true; }
        siteName = pick->set;
        c.note(TextFormat("site: %s (%s%s), %.0f m from the start", pick->set.c_str(), pick->kind.c_str(), pick->dungeon ? ", a dungeon under it" : "", pick->distance));
        // 200 m from it along the great circle toward the start: inside the build range, and the site's first BuildAround
        // happens with the player there (as walking in does, a step later).
        const Vector3 sd = Vector3Normalize(pick->at), hd = Vector3Normalize(here);
        const Vector3 perp = Vector3Normalize(Vector3Subtract(hd, Vector3Scale(sd, Vector3DotProduct(hd, sd))));
        const float a = 200.0f / c.field->radius;
        const Vector3 d = Vector3Normalize(Vector3Add(Vector3Scale(sd, cosf(a)), Vector3Scale(perp, sinf(a))));
        sites0 = Structures_SitesBuilt();
        strokes0 = c.body->edits->strokes;
        worstMs = 0.0;
        player.flyMode = false;
        player.placeAt(c.world, *c.body, c.body->toWorld(Vector3Scale(d, c.field->radius + c.field->elevation(d) + 1.0f)));
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
        c.timingLt("Structures_Update max ms (GUESS 4)", worstMs, 4.0);
        frames.report(c);
        return true;
    });
}
