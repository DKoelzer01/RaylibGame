#ifndef ROTATIONSOURCES_H
#define ROTATIONSOURCES_H

// Water wheels and windmills (#159, #1270; the Factory and Progression doc, sections 5.1 and 2.9): the main thread samples
// each one's river or wind twice a second and hands the factory its speed (Factory::SetSourceRpm); the factory's rotation
// step turns it into speed and stress capacity like any source.
//
// - Water wheel (1.5 m across): the water's speed across its rim at the bottom. Generated rivers: RiverCurrent (from the
//   bake, valid anywhere). Simulated water (a dug channel, a dammed or diverted river, a pond's outflow): Fluid::FlowAt,
//   which wins wherever a simulated cell decides (so a dammed river's wheel slows). Still water: nothing.
// - Windmill: the wind's speed at its hub (Wind_At: the mean wind x the daily cycle x the terrain's exposure x height;
//   gusts only near the player), from the cut-in speed up to full speed at kWindFullMps. Calm: it stops.
// - Far plots (decided, 2.9): the sampled speed is what the full simulation would see at most: a wheel on simulated water
//   uses the lowest speed seen over the last minute (the fluid sim does not run there; the lowest sample keeps the factory
//   at or under the full run); windmills far from the player have no gusts anyway.

#include <raylib.h>
#include <deque>

class Scene;

namespace RotationSources {
// Tuning (the doc's guesses, as named constants).
constexpr float kWheelRadius = 0.75f;   // m: 3 cells across
constexpr float kWheelSlip = 0.8f;      // the rim turns at 80% of the water's speed
constexpr float kWheelDipArea = 0.75f;  // m^2 dipped (1.5 m wide x 0.5 m deep): simulated flow (m^3/s) to speed (m/s)
constexpr float kWindCutInMps = 1.5f;   // under this a windmill stands
constexpr float kWindFullMps = 8.0f;    // at and over this it turns at full speed
constexpr double kSampleSeconds = 0.5;
constexpr int kFarWindow = 120;         // samples (60 s) whose lowest a far wheel on simulated water uses
constexpr float kNearMetres = 1000.0f;  // as the factory's "near"

// A water wheel's speed (signed RPM about +axis) from the water's velocity at the bottom of its rim; all vectors in one frame.
int WaterWheelRpm(Vector3 waterVelocity, Vector3 axis, Vector3 up);
// A windmill's speed (RPM, >= 0) in a wind of this speed (m/s).
int WindmillRpm(float windSpeed);
// The lowest-magnitude speed in the window (0 if they turn both ways): what a far wheel on simulated water gets.
int LowestRpm(const std::deque<int>& window);
}

void RotationSources_Update(Scene& world);   // main thread, every frame (samples every kSampleSeconds of sim time)
void RotationSources_Reset();                // a new world / a load
// Tests: every windmill sees this wind speed (m/s) instead of Wind_At (on = false: back to the real wind).
void RotationSources_SetWindOverride(bool on, float speed = 0.0f);

#endif
