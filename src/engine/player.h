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
    int spawns = 0;                   // counts spawnOn/placeAt: the game raises the loading screen on a new one
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
    Vector3 groundNormal = { 0, 1, 0 }; // contact normal while grounded
    float eyeOffset = 0.0f;           // camera height smoothing (metres along up, see update)
    float eyeVz = 0.0f;               // low-passed vertical speed of the feet in the ground frame
    bool eyeSmoothing = true;         // tests: --no-eye-smooth turns it off (A/B)
    float eyeStep = 0.0f;             // tests: the eye's height change in the last update (ground frame, m)
    bool eyeStepKnown = false;

    // Tunables
    float eyeHeight = 1.6f;
    float radius = 0.45f;
    float walkSpeed = 6.0f;
    float runMultiplier = 2.0f;
    float jumpSpeed = 5.0f;
    float jetpackAccel = 9.0f;
    float jetpackFuelPerSec = 0.1f;   // full tank = 10 s of thrust
    float airControlAccel = 12.0f;    // m/s^2 toward the wished walking velocity while airborne in gravity
    float walkableCos = 0.55f;        // ground steeper than this (cos of its slope, ~57 deg) is a wall; grip boots lower it
    float mouseSensitivity = 0.003f;
    bool smoothMotion = false;        // Settings "Reduce motion": eased walking and a smoothed camera up
    mutable Vector3 cameraUp = { 0, 1, 0 };
    float lookScale = 1.0f;        // mouse look multiplier (the telescope slows it by its zoom)

    // Stand on the surface of `body` in direction `dir` from its centre (default: its +Y pole).
    void spawnOn(Scene& world, Planetoid& body, Vector3 dir = { 0, 1, 0 });
    // Put the feet at a render-space point, moving with the body there (tests, respawns).
    void placeAt(Scene& world, Planetoid& body, Vector3 renderFeet);
    // Ride along with body from the next update on (else one slow frame leaves you behind its spin).
    void anchorTo(Planetoid& body);
    void update(Scene& world, bool inputEnabled, float dt);
    void applyToCamera(Scene& world) const;  // camera at the eye, looking along forward
    // Fly mode follows the reference body (the nearest whose sphere of influence holds the camera):
    // its spin and motion within 2 radii, its motion only further out. recordFly() runs after the
    // camera's own input each frame (Scene::drawScene); update() then carries it with the body.
    bool flyFollow = true;
    void recordFly(Scene& world);
    Planetoid* flyReference(Scene& world, bool& rotate) const;

private:
    DVec3 groundLocal{};              // feet in the ground body's local space, double (follows its spin exactly)
    Planetoid* groundBody = nullptr;
    uint64_t groundId = 0;            // validates groundBody (bodies can be destroyed, e.g. world regen)

    Quaternion groundRotPrev = { 0, 0, 0, 1 }; // ground body's rotation at the end of last update
    DVec3 lastEndPosition{};                    // where update() left us (detects outside moves)
    // Airborne near a body: carried with its motion each frame and moved by velocity relative to its
    // surface, so substep collisions line up with the (already moved) terrain. Same inertial path to
    // first order; without it you ran up to a frame of spin (86 m/s at 60 km) ahead of the ground.
    Planetoid* airBody = nullptr;
    uint64_t airId = 0;
    DVec3 airLocal{};
    // Airborne near an orbiting body we move in its free-falling frame: its velocity at the end of
    // last frame (its change is added to ours), and whether this frame uses felt gravity.
    Vector3 airFrameVel{};
    bool airFelt = false;
    uint64_t flyBodyId = 0;           // body the fly camera was recorded against (0: none)
    bool flyRotate = false;
    DVec3 flyLocal{}, flyBodyPos{};
    Vector3 flyLocalFwd{}, flyLocalUp{};
    DVec3 flyCarried{};               // the carried camera in double (float render space rounds it)
    bool flyCarriedValid = false;
    void stepPhysics(Scene& world, float dt, Vector3 wishMove, bool jetUp, bool jetDown);
    float jumpBuffer = 0.0f; // seconds a jump press stays pending
    float coyote = 0.0f;     // seconds a jump is still allowed after leaving the ground
    float jumpLock = 0.0f;   // seconds after a jump during which the player cannot re-ground
    bool resolveCollisions(Scene& world, Vector3& feet);
    void collideTrunks(Vector3& feet); // tree trunks: vertical capsules on the nearest body
    // Velocity of the terrain being collided with, in the frame `velocity` is in: zero while grounded
    // (velocity is then relative to the surface), the surface velocity while airborne (world frame).
    Vector3 collideFrameVel = { 0, 0, 0 };
};

extern Player player;

#endif
