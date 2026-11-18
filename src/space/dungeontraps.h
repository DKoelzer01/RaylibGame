#ifndef DUNGEONTRAPS_H
#define DUNGEONTRAPS_H

// Dungeon traps (#719 slice 2): pressure plates and gas vents, plot parts the dungeon planner places (structures.cpp), so
// they save with the dungeon and need no state of their own.
// - A plate (PressurePlate, material 0) drops a stone block from the ceiling: a click and a puff of dust warn, the block lands
//   0.5 s later on whoever is still under it and stays as rubble (a stone cube); the plate is spent.
// - A dart plate (material 1) shoots a dart across the corridor from its nearer wall; it re-arms after 2.5 s.
// - A gas vent (GasVent) fills the room round it with poison: the air there is not breathable and the suit's oxygen drains
//   three times as fast (the filters choke), then suffocation hurts as anywhere.
// They spring for the local player in single player and on a host; a joined client's traps wait for the server (#337).
#include <raylib.h>

class Scene;
class Planetoid;

void DungeonTraps_Update(Scene& world);            // per frame, after the player moved (timed on the sim clock)
void DungeonTraps_Draw(Scene& world);              // inside BeginMode3D (render space): falling blocks and darts
// Whether a render-space point is in a vent's gas (the survival step's breathable air).
bool DungeonTraps_InGas(const Scene& world, Vector3 renderPos);
void DungeonTraps_Reset();                         // a new world

// Tests (#719): counters since the reset.
struct DungeonTrapStats { int plates = 0, vents = 0, blocksDropped = 0, blockHits = 0, dartsShot = 0, dartHits = 0; };
const DungeonTrapStats& DungeonTraps_Stats();
void DungeonTraps_Rescan();                        // look for traps now (they are rescanned every 0.5 s)
int DungeonTraps_Test(Scene& world, Planetoid& body, int frame); // --test-dungeon-traps: 0 running, 1 pass, -1 fail

#endif
