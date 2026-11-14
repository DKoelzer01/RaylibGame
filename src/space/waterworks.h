#ifndef WATERWORKS_H
#define WATERWORKS_H

// Waterworks (Plan 4, F6): water pumps, pipes and outlets on build plots.
//
// A water pump is a powered LV machine (factory.h): while it has power its recipe fills its outputs
// with water units. Twice a second each pump's units go down its pipe network - pipe parts joined face
// to face, starting next to the pump - to every water outlet the network touches: the water is taken
// from under the pump (a river, a pond, the sea, moving water: fluid.h) and poured out of the outlets,
// split evenly. A pump with no outlet on its network keeps its water (its outputs fill and it stops);
// one with nothing under it pumps air (the units are lost).

class Scene;

constexpr float kWaterPerUnit = 0.5f; // m^3 per water unit (a pump makes one a second)

void Waterworks_Update(Scene& world, float dt);
// Hand pump (#280): the player works it (hold RMB on it): dt seconds of pumping, 50 litres a second, sent down its
// pipes at the next pass like a water pump's units (no power; no water under it: nothing).
constexpr float kHandPumpRate = 0.05f; // m^3 a second
#include <cstdint>
void Waterworks_HandPump(uint64_t body, uint32_t plot, int part, float dt);

struct WaterworksStats { int pumps = 0, networks = 0, outlets = 0; double pumped = 0.0; int dry = 0; };
WaterworksStats Waterworks_Stats(); // totals since the start (pumped m^3), and the last pass (pumps, networks, outlets, dry)

#endif
