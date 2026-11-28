// saves/server-player (#2157): a dedicated server (--host without --host-draw: headless) keeps no record of its own player,
// so a load left that player at its default, body-local zero. Streaming, cave gas (#2153: the --host hang) and the render
// origin all centre on it. Now a server's load puts it at the world's start, as a first run does. The scenario moves the
// player near the core (where the bug left it), saves, loads the save back and checks the player is at the start again.
// No GameNetworkingSockets needed: the host loads and saves whether or not it can listen.
#include "scenario.h"
#include "headless.h"
#include "player.h"
#include "world.h"
#include "worldtypes.h"

namespace {
struct State { DVec3 start{}; } S;
}

SCENARIO("saves/server-player") {
    s.tags = { "saves", "net" }; // net: the runner's part2 step (a server load)
    s.world = ScenarioWorld::Sandbox;
    s.args = { "--host", "saves/scenario-saves-server-player.sav", "--port", "27957" };
    s.record = false; // headless: nothing to record
    s.step("the server's player is at the start", 30.0, [](ScenarioCtx& c) {
        c.waitingFor("an idle world");
        if (c.T < 1.0 || !WorldIdle() || !c.body) return false;
        c.check("a dedicated server (headless)", Headless());
        S.start = c.body->toLocalD(player.position);
        c.gt("the player is away from the core (m)", S.start.length(), 1000.0);
        return true;
    });
    s.step("saved with the player at the core and loaded back: at the start again", 30.0, [](ScenarioCtx& c) {
        player.position = c.body->toWorldD(DVec3(0.0, 40.0, 0.0)); // where the bug left it (#2157: local 0, 40, 0)
        const bool saved = SaveWorld(c.world), loaded = saved && LoadWorldSave(c.world);
        c.check("saved and loaded", saved && loaded, LastLoadError());
        const DVec3 now = c.body->toLocalD(player.position);
        c.gt("the player is away from the core after the load (m)", now.length(), 1000.0);
        c.lt("the player is at the start (m from the first spot)", (now - S.start).length(), 50.0);
        return true;
    });
}
