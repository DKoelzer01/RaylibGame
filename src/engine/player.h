#ifndef PLAYER_H
#define PLAYER_H

#include <raylib.h>
#include "dvec3.h"

class Scene;
class Planetoid;

// First-person character for small bodies:
// - gravity pulls toward the nearest planetoid's centre and "up" follows it, so you can walk all
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
    Planetoid* nearest = nullptr;     // body whose gravity applies
    float gravity = 0.0f;             // current acceleration magnitude (HUD/debug)
    bool autoWalk = false;            // tests: walk forward without input
    float gravityScale = 1.0f;        // spells/effects (e.g. featherfall)
    float lastContactDot = 0.0f;      // debug: dot(contact normal, up) of the last collision
    float impactSpeed = 0.0f;         // hardest collision speed this frame (fall damage)

    // Tunables
    float eyeHeight = 1.6f;
    float radius = 0.45f;
    float walkSpeed = 6.0f;
    float runMultiplier = 2.0f;
    float jumpSpeed = 5.0f;
    float jetpackAccel = 9.0f;
    float mouseSensitivity = 0.003f;

    void spawnOn(Scene& world, Planetoid& body);
    void update(Scene& world, bool inputEnabled, float dt);
    void applyToCamera(Scene& world) const;  // camera at the eye, looking along forward

private:
    Vector3 groundLocal{};            // feet in the ground body's local space (to follow its spin)
    Planetoid* groundBody = nullptr;
    Vector3 groundNormal = { 0, 1, 0 };
    void stepPhysics(Scene& world, float dt, Vector3 wishMove, bool jump, bool jetUp, bool jetDown);
    bool resolveCollisions(Scene& world, Vector3& feet);
};

extern Player player;

#endif
