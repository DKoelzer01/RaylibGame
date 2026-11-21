#ifndef FOOTSTEPS_H
#define FOOTSTEPS_H

// Footsteps by the surface underfoot and landings by fall speed (#1240, audio plan #1003 step 4). The local player and every
// remote player (their avatars, from what multiplayer already replicates): one step per stride on the ground, the
// sound of what is underfoot (terrain material, or the part stood on), louder running, quiet crouching or prone.
// The same loudness is a gameplay number: every step is also a noise event that creature hearing (#1223) can read.
// Steps are audio v2 sounds (procedural "step_<surface>"; resources/audio/step_<surface>.wav etc. replace them); with the
// old engine nothing plays, the noise events still happen.

#include "dvec3.h"
#include "player.h"
#include "voxel.h"
#include "plots.h"

#include <cstdint>
#include <vector>

class Scene;
struct AvatarState;

enum class StepSurface : uint8_t { Stone, Dirt, Grass, Sand, Gravel, Snow, Wood, Metal, Water, Count };
const char* StepSurface_Sound(StepSurface s);         // "step_stone", ...
StepSurface StepSurfaceOfVoxel(VoxelMaterial m);
StepSurface StepSurfaceOfPart(const Part& p);
// How loud a step is (1 = walking on stone): posture and pace, times the surface (grass and sand soft, metal and gravel loud).
float StepLoudness(StepSurface s, Player::Posture posture, bool running);

// One walker's strides (pure): distance along the ground since the last step, and touchdowns.
struct StepTracker {
    static constexpr float kStride = 1.8f;           // one step per 1.8 m standing (the hand bob's side swing, viewmodel.cpp)
    static constexpr float kLandingSpeed = 2.5f;     // m/s down at touchdown: a landing, not a step
    float distance = kStride * 0.6f;                 // the first step comes soon after starting off
    bool wasGrounded = true;
    float lastDownSpeed = 0.0f;
    enum Event { None, Step, Landing };
    // moved: metres along the ground this frame; downSpeed: m/s into the ground (for remote walkers, whose touchdown
    // speed is not sent: the speed the frame before they landed). Returns what to play; `landing` = the touchdown speed.
    Event update(bool grounded, float moved, float downSpeed, Player::Posture posture, float& landing);
};

struct StepNoise { DVec3 pos; float loudness; double time; uint32_t who; }; // who: 0 the local player, else a session id
// Steps and landings of the last 2 s (creature hearing, #1223).
const std::vector<StepNoise>& Footsteps_RecentNoises();

// Once a frame, after the player's update (and the avatars' interpolation).
void Footsteps_Update(Scene& world, const Player& self, const std::vector<AvatarState>& others, float dt);

#endif
