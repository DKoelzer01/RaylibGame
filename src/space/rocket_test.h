#ifndef ROCKET_TEST_H
#define ROCKET_TEST_H

class Scene;
// Milestone 7 gate (--test-rocket): 1 pass, -1 fail, 0 running.
int RocketTest_Update(Scene& world);
// Plan 3 M20 gate (--test-station): deploy a station in orbit, dock, move fuel, undock, save and load.
int StationTest_Update(Scene& world);

#endif
