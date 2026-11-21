// visual/shot-* (#870 visual tier): every --shot-* screenshot bench as a visual scenario. Each runs its legacy flag
// unchanged (legacyDrives: the bench drives the world, the scenario never holds it) on the planet at seed 42, then every
// bench*.png the bench wrote is compared with its baseline on the share (Scenario_BaselineDir()/visual_shot-NAME/), as
// ScenarioCtx::shot does: a failure keeps the shot, a difference image and the recording. The first run records the
// baselines; RAYLIBGAME_BASELINE_UPDATE=1 re-records them after an intended look change.
//
// They move the benches into the registry without rewriting them; a bench that later becomes scenario steps (as
// visual/lamps and visual/torch already are) replaces its wrapper here.
#include "scenario.h"

namespace {
void Shot(ScenarioDef& s, const char* flag, double timeout = 240.0, std::vector<std::string> extra = {}) {
    s.tags = { "visual", "shot", "planet" };
    s.world = ScenarioWorld::Planet;
    s.timeout = timeout;
    s.legacyFlag = flag;
    s.legacyDrives = true;
    s.args = { flag };
    for (auto& a : extra) s.args.push_back(a);
}
}

// The world and its look
SCENARIO("visual/shot-terrain") { Shot(s, "--shot-terrain"); }
SCENARIO("visual/shot-mountain") { Shot(s, "--shot-mountain"); }
SCENARIO("visual/shot-cave-mouth") { Shot(s, "--shot-cave-mouth"); }
SCENARIO("visual/shot-hillshade") { Shot(s, "--shot-hillshade"); }
SCENARIO("visual/shot-biomes") { Shot(s, "--shot-biomes", 420.0); }
SCENARIO("visual/shot-forest") { Shot(s, "--shot-forest"); }
SCENARIO("visual/shot-far-forest") { Shot(s, "--shot-far-forest"); }
SCENARIO("visual/shot-giant") { Shot(s, "--shot-giant"); }
SCENARIO("visual/shot-shadow") { Shot(s, "--shot-shadow"); }
SCENARIO("visual/shot-particles") { Shot(s, "--shot-particles"); }
SCENARIO("visual/shot-clouds") { Shot(s, "--shot-clouds"); }
SCENARIO("visual/shot-limb") { Shot(s, "--shot-limb"); }
SCENARIO("visual/shot-dig") { Shot(s, "--shot-dig"); }
// Things built and placed
SCENARIO("visual/shot-campfire") { Shot(s, "--shot-campfire"); }
SCENARIO("visual/shot-shelter") { Shot(s, "--shot-shelter"); }
SCENARIO("visual/shot-room") { Shot(s, "--shot-room"); }
SCENARIO("visual/shot-machines") { Shot(s, "--shot-machines"); }
SCENARIO("visual/shot-glow") { Shot(s, "--shot-glow"); }
SCENARIO("visual/shot-chest") { Shot(s, "--shot-chest"); }
SCENARIO("visual/shot-preview") { Shot(s, "--shot-preview"); }
SCENARIO("visual/shot-graves") { Shot(s, "--shot-graves"); }
SCENARIO("visual/shot-claims") { Shot(s, "--shot-claims"); }
// Creatures and players
SCENARIO("visual/shot-animals") { Shot(s, "--shot-animals"); }
SCENARIO("visual/shot-avatar") { Shot(s, "--shot-avatar"); }
// The interface
SCENARIO("visual/shot-hud") { Shot(s, "--shot-hud"); }
SCENARIO("visual/shot-inventory") { Shot(s, "--shot-inventory"); }
SCENARIO("visual/shot-game-window") { Shot(s, "--shot-game-window"); }
SCENARIO("visual/shot-tabs") { Shot(s, "--shot-tabs"); }
SCENARIO("visual/shot-team") { Shot(s, "--shot-team"); }
SCENARIO("visual/shot-icons") { Shot(s, "--shot-icons"); }
SCENARIO("visual/shot-settings") { Shot(s, "--shot-settings"); }
SCENARIO("visual/shot-mp-menu") { Shot(s, "--shot-mp-menu"); }
