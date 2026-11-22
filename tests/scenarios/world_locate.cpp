// world/locate (#1171): /locate and /tp on the real planet (seed 42), from the start. The playtest skipped a volcano, lava
// tubes, ruins and dungeons, cave mouths and karst because nobody could find them. Every feature is asked for; the
// common ones (a cave mouth, a ruin, a spring, a lake, a biome) must be found; the rare ones are reported. Then /tp to the
// nearest cave mouth and to the nearest ruin: the player must arrive there (the camera follows) and stand on the ground.
#include "scenario.h"
#include "commandline.h"
#include "locate.h"
#include "playersession.h"
#include "terraintool.h"
#include "world.h"
#include <raymath.h>
#include <string>

namespace {
Vector3 target{};
double tpAt = -1.0;
std::string tpWhat;
bool Arrived(ScenarioCtx& c, const char* what) {
    Player& p = LocalSession().player;
    if (!p.nearest) return false;
    const Vector3 here = p.nearest->toLocalD(p.position).toFloat();
    const float d = Vector3Distance(here, target);
    // The subject in view: the camera a few metres behind and above the player, looking at the target.
    const Vector3 up = Vector3Normalize(target);
    c.camera(Vector3Add(here, Vector3Add(Vector3Scale(up, 4.0f), Vector3Scale(Vector3Normalize(Vector3CrossProduct(up, { 0, 1, 0 })), 6.0f))), target, up);
    c.waitingFor(TextFormat("the player at the %s (%.0f m away, idle %d)", what, d, WorldIdle() ? 1 : 0));
    if (!WorldIdle() || c.t < 3.0) return false;
    c.lt(TextFormat("/tp %s: the player arrived (m from the target)", what), d, 30.0);
    c.note(TextFormat("/tp %s: %s", what, tpWhat.c_str()));
    return true;
}
}

SCENARIO("world/locate") {
    s.tags = { "world" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 240.0;
    s.step("locate every feature from the start", 180.0, [](ScenarioCtx& c) { // (most of it waiting for the spawn screen on a slow machine)
        Player& p = LocalSession().player;
        if (c.T < 1.0) return false;
        terrainTool.creative = true; // (the commands are creative's)
        // From the player where it stands; before its first update (no body yet) from the spawn point on the planet
        Planetoid* const body = p.nearest ? p.nearest : c.body;
        const Vector3 from = p.nearest ? p.nearest->toLocalD(p.position).toFloat() : Vector3Scale(c.field->findSpawnDir(c.body->seed), c.field->radius + 2.0f);
        if (!p.nearest) c.note("the player has no body yet: searching from the spawn point");
        const char* common[] = { "cave", "ruin", "spring", "lake", "biome:desert", "ore", "ore:iron" };
        const char* rare[] = { "hillside", "sinkhole", "resurgence", "karst", "tube", "chamber", "volcano", "dungeon", "camp", "wreck", "wonder", "oilseep", "river" };
        for (const char* f : common) {
            std::vector<LocateHit> hits; std::string why;
            const bool ok = Locate_Find(*body, from, f, 3, hits, why);
            c.check(TextFormat("/locate %s finds one", f), ok, ok ? TextFormat("%s at %.0f m", hits[0].what.c_str(), hits[0].distance) : why);
        }
        for (const char* f : rare) {
            std::vector<LocateHit> hits; std::string why;
            if (Locate_Find(*body, from, f, 1, hits, why)) c.note(TextFormat("/locate %s: %s, %.1f km %s", f, hits[0].what.c_str(), hits[0].distance / 1000.0f, Locate_CompassPoint(hits[0].bearing)));
            else c.note(TextFormat("/locate %s: %s", f, why.c_str()));
        }
        // The command itself (a waypoint, its reply) and an unknown feature's answer.
        c.note(("/locate cave 2 -> " + CommandLine_Run("/locate cave 2")).c_str());
        c.check("an unknown feature lists the features", CommandLine_Run("/locate unicorn").find("Features:") != std::string::npos);
        terrainTool.creative = false;
        c.check("survival: /tp is refused", CommandLine_Run("/tp cave").find("creative") != std::string::npos);
        terrainTool.creative = true;
        return true;
    });
    s.step("/tp cave: the player goes to the nearest cave mouth", 90.0, [](ScenarioCtx& c) {
        if (tpAt < 0.0) {
            Player& p = LocalSession().player;
            if (!p.nearest) return false;
            std::vector<LocateHit> hits; std::string why;
            if (!Locate_Find(*p.nearest, p.nearest->toLocalD(p.position).toFloat(), "cave", 1, hits, why)) { c.fail("no cave to go to: " + why); return true; }
            target = hits[0].at;
            tpWhat = CommandLine_Run("/tp cave");
            tpAt = c.T;
            return false;
        }
        if (!Arrived(c, "cave mouth")) return false;
        tpAt = -1.0;
        return true;
    });
    s.step("/tp ruin: the player goes to the nearest ruin", 90.0, [](ScenarioCtx& c) {
        if (tpAt < 0.0) {
            Player& p = LocalSession().player;
            if (!p.nearest) return false;
            std::vector<LocateHit> hits; std::string why;
            if (!Locate_Find(*p.nearest, p.nearest->toLocalD(p.position).toFloat(), "ruin", 1, hits, why)) { c.fail("no ruin to go to: " + why); return true; }
            target = hits[0].at;
            tpWhat = CommandLine_Run("/tp ruin");
            tpAt = c.T;
            return false;
        }
        if (!Arrived(c, "ruin")) return false;
        return true;
    });
}
