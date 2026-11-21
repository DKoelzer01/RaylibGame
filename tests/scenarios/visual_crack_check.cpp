// visual/crack-check (#870 visual tier; moved from --test-stand 3 --look-down --crack-check, #725): with the skybox off and
// a magenta clear colour (--crack-check: no HUD either), holes between terrain chunks show magenta. The player stands at
// the start looking ~50 deg down (--look-down) at terrain filling the view; after the world settles, the magenta pixels
// enclosed by drawn ones (CrackHoles: not the open sky) must stay under 50 (the share is a metric). The flag printed the
// count only; 50 is the level its fly check already treats as a hole worth logging.
#include "scenario.h"
#include "crackcheck.h"
#include "player.h"
#include "world.h"

SCENARIO("visual/crack-check") {
    s.tags = { "visual", "terrain", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = 120.0;
    s.args = { "--crack-check", "--look-down" };
    s.legacyFlag = "--test-stand 3 --look-down --crack-check";
    s.step("stand at the start", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle()) return false;
        player.flyMode = false;
        player.spawnOn(c.world, *c.body, c.spawnDirWorld());
        return true;
    });
    s.step("settle (3 s standing, the world idle)", 30.0, [](ScenarioCtx& c) { c.waitingFor("3 s and an idle world"); return c.t > 3.0 && WorldIdle(); });
    s.step("count the holes", 5.0, [](ScenarioCtx& c) {
        Image shot = LoadImageFromScreen();
        Color* px = LoadImageColors(shot);
        long sky = 0;
        const long holes = CrackHoles(px, shot.width, shot.height, 1, &sky);
        const double total = (double)shot.width * shot.height;
        UnloadImageColors(px);
        UnloadImage(shot);
        c.metric("hole share", 100.0 * holes / total, "%");
        c.metric("open sky pixels (not counted)", (double)sky, "");
        c.lt("hole pixels", (double)holes, 50.0);
        return true;
    });
}
