#ifndef PLAYER_H
#define PLAYER_H

#include <raylib.h>
#include "dvec3.h"
#include <cstdint>
#include <algorithm>

class Scene;
class Planetoid;

// First-person character for small bodies:
// - gravity is summed over every body (GravityAt) and "up" follows it, so you can walk all
//   the way around an asteroid;
// - collision is resolved directly against the density field (two spheres: feet and head), so
//   dug tunnels and built terrain are walkable with no separate collision mesh;
// - standing on a spinning body carries you with it;
// - jump, and a jetpack while airborne (Space up, Ctrl down).
struct Player {
    DVec3 position{};                 // feet, world space (double: survives the floating origin)
    Vector3 velocity{};               // world units / s
    Vector3 forward = { 0, 0, 1 };    // view direction
    Vector3 up = { 0, 1, 0 };
    bool grounded = false;
    bool flyMode = false;             // F: free-fly dev camera (no gravity/collision)
    float jetpackFuel = 1.0f;         // 0..1, recharges on the ground
    bool jetpackFiring = false;       // fuel was burnt this update (survival / audio / HUD)
    Planetoid* nearest = nullptr;     // body we collide with (gravity comes from every body)
    float gravity = 0.0f;             // current acceleration magnitude (HUD/debug)
    bool autoWalk = false;            // tests: walk forward without input
    float gravityScale = 1.0f;        // spells/effects (e.g. featherfall)
    bool dampeners = true;            // B: inertial dampeners (deep space only)
    bool autoJump = false;            // tests: jump once on the next update
    float lastContactDot = 0.0f;      // debug: dot(contact normal, up) of the last collision
    float impactSpeed = 0.0f;         // hardest collision speed this frame (fall damage)

    // Tunables
    float eyeHeight = 1.6f;
    float radius = 0.45f;
    float walkSpeed = 6.0f;
    float runMultiplier = 2.0f;
    float jumpSpeed = 5.0f;
    float jetpackAccel = 9.0f;
    float jetpackFuelPerSec = 0.1f;   // full tank = 10 s of thrust
    float airControlAccel = 12.0f;    // m/s^2 toward the wished walking velocity while airborne in gravity
    float mouseSensitivity = 0.003f;

    // Stand on the surface of `body` in direction `dir` from its centre (default: its +Y pole).
    void spawnOn(Scene& world, Planetoid& body, Vector3 dir = { 0, 1, 0 });
    void update(Scene& world, bool inputEnabled, float dt);
    void applyToCamera(Scene& world) const;  // camera at the eye, looking along forward

private:
    DVec3 groundLocal{};              // feet in the ground body's local space, double (follows its spin exactly)
    Planetoid* groundBody = nullptr;
    uint64_t groundId = 0;            // validates groundBody (bodies can be destroyed, e.g. world regen)
    Vector3 groundNormal = { 0, 1, 0 };
    Quaternion groundRotPrev = { 0, 0, 0, 1 }; // ground body's rotation at the end of last update
    DVec3 lastEndPosition{};                    // where update() left us (detects outside moves)
    void stepPhysics(Scene& world, float dt, Vector3 wishMove, bool jetUp, bool jetDown);
    float jumpBuffer = 0.0f; // seconds a jump press stays pending
    float coyote = 0.0f;     // seconds a jump is still allowed after leaving the ground
    float jumpLock = 0.0f;   // seconds after a jump during which the player cannot re-ground
    bool resolveCollisions(Scene& world, Vector3& feet);
};

extern Player player;

#endif
