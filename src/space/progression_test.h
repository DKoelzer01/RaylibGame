#ifndef PROGRESSION_TEST_H
#define PROGRESSION_TEST_H

class Scene;
// Milestone 5 gate (--test-progression): one step a frame; 1 done (pass), -1 failed, 0 running.
int Progression_Update(Scene& world, bool fullLadder, bool evLadder = false); // fullLadder: on through MV and HV to titanium alloy (milestone 8); evLadder: tiers 9-10 from a moon/ring haul (M24)

#endif
