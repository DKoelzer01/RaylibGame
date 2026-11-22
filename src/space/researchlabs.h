#pragma once

// Research labs at work (#1193, #1030): every 0.5 s (single player and a host's own team for now; joined clients read the
// server's later), each lab is set to run the first node of its team's queue it can (its tier at most the lab's: LV 4,
// MV 5, HV 7), and the research units it made are taken out and added to that node; a node done leaves the queue with a
// toast. Labs on plots the factory doesn't simulate (far ones) wait for the summary sim's step (doc section 8).

class Scene;

int ResearchLabs_Update(Scene& world, double now, bool force = false); // units counted this call
int ResearchLabs_UnitsCounted();                                         // (tests) since the reset
void ResearchLabs_Reset();
