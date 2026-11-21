#ifndef MACHINESOUNDS_H
#define MACHINESOUNDS_H

// Machine sounds (#1243, audio plan #1003 step 7): a hum per machine type and tier (fire, steam, LV, MV, HV, EV) while it
// works, quieter while an electric machine idles; state cues as one-shots (stalled: a falling tone then a click every few
// seconds; a batch done: a soft chime); and the factory bed: the nearest kIndividual machines play on their own, the rest
// of the working ones within range merge into one bed whose loudness grows with their count only slowly, so a big
// factory stays tame (user decision, with the Machines slider). Everything on the Machines bus. Audio v2 only.

#include "plots.h"

#include <vector>

class Scene;

constexpr int kMachineIndividual = 8;    // machines heard one by one (the nearest)

// The hum a machine type makes ("hum_lv", ...), or nullptr (chests, parts that make no sound of their own).
const char* MachineSound_Hum(PartType t);
// (#1269) Mechanical power: shafts, cogs, gearboxes and cranks (a plot's sound as one clatter at their middle), and the
// loops' pitch from the network's speed (made at 64 RPM; 0.5x..2x).
bool MachineSound_IsTransmission(PartType t);
float MachineSound_RpmPitch(float rpm);
// Whether it hums softly when idle (powered electric machines do; fires and steam engines are silent).
bool MachineSound_IdleHums(PartType t);

// 0 idle, 1 working, 2 stalled (no power, output full, no fuel, no water).
enum class MachineCue : uint8_t { None, Stall, Done };
MachineCue MachineSound_Cue(int was, int now);   // a change of state -> the cue to play
// The factory bed's volume for `count` working machines beyond the individual ones: grows like a log, at most 1.
float MachineSound_BedVolume(int count);

// Pure selection: which machines (by index into `distances`, sorted nearest first) get their own voice: the nearest `max`
// of those that make a sound (`audible`).
void MachineSound_Pick(const std::vector<float>& distances, const std::vector<bool>& audible, int max, std::vector<int>& out);

void MachineSounds_Update(Scene& world, float dt);   // once a frame
void MachineSounds_Reset();                          // a new or unloaded world

#endif
