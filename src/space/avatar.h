#ifndef AVATAR_H
#define AVATAR_H

// Remote player avatars (#328, multiplayer): a procedural humanoid in a suit, code-built like the animal models (Plan 5 A),
// lit and batched through the instance batch (sun, shadows, headlamp). It walks, runs, crouches, lies prone, swims, falls
// and flies; holds its item or tool, swings it when it digs or chops; its helmet lamp glows when the headlamp is on; a
// nameplate floats over it. Single player never draws one (except RAYLIBGAME_AVATAR_SELF=1, a debug view of yourself).
// The network side (src/space/netclient.cpp, #327) fills an AvatarState per remote player each frame, already
// interpolated; this file only draws from it. A glTF body can replace the procedural one later: RAYLIBGAME_AVATAR_GLTF.

#include <raylib.h>
#include <cstdint>
#include <string>
#include <vector>
#include "dvec3.h"
#include "player.h"

class Scene;

struct AvatarState {
    uint32_t session = 0;                 // the server's session id (stable per connection): keys the animation state
    std::string name;                     // nameplate ("" = none)
    DVec3 position{};                     // feet, world (interpolated)
    Vector3 up{ 0, 1, 0 };                // body up (world, unit)
    Vector3 facing{ 0, 0, 1 };            // facing along the ground (world, unit)
    Vector3 velocity{};                   // world m/s relative to the ground (drives the gait)
    float lookPitch = 0.0f;               // head and arms tilt, radians (+ up)
    Player::Posture posture = Player::Posture::Stand;
    bool grounded = true, flying = false, swimming = false, jetpack = false, headlamp = false;
    std::string held;                     // item or tool id in the hand ("" = empty)
    int action = 0;                       // 0 none, 1 dig / mine swing, 2 chop, 3 place / use
    float actionPhase = 0.0f;             // 0..1 within the current swing
    int suitTier = 0;                     // suit upgrades show as accent colours
    Color suitColor{ 0, 0, 0, 0 };        // the colour the player picked when joining (#341); alpha 0: the default suit
};

// Draw inside the 3D pass, before InstanceBatch_Flush (the boxes go through the instance batch).
void Avatar_Draw3D(Scene& world, const std::vector<AvatarState>& avatars);
// After the 3D pass (2D): names over the heads, fading with distance.
void Avatar_DrawNameplates(Scene& world, const std::vector<AvatarState>& avatars);
// A Player's own state as an avatar (the local player for RAYLIBGAME_AVATAR_SELF, bench players); the extras (held item,
// action, headlamp, name) are filled by the caller.
AvatarState Avatar_FromPlayer(const Player& p, uint32_t session);

#endif
