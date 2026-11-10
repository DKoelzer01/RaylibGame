#ifndef SHIP_H
#define SHIP_H

#include <raylib.h>
#include <raymath.h>
#include "dvec3.h"
#include <cstdint>

class Scene;
class Planetoid;

// A small pilotable ship. Newtonian 6-DOF flight with assisted damping, chase camera, terrain
// collision against the density field, an energy store recharged with power cells, and a
// life-support field that refills the player's suit.
struct Ship {
    DVec3 position{};                      // world (double)
    Quaternion orientation = QuaternionIdentity();
    Vector3 velocity{};
    float energy = 100.0f;                 // 0..100
    bool piloted = false;
    bool spawned = false;
    bool autopilot = false;                // fly to autopilotTarget (tests / assist)
    DVec3 autopilotTarget{};
    float lastImpact = 0.0f;               // hardest terrain hit this frame (reset by caller)
    bool hovering = false;                 // flight assist is holding the ship against gravity (HUD/tests)
    // Parked on a body: carried along with its spin (see Ship::update).
    Planetoid* anchorBody = nullptr;
    uint64_t anchorId = 0;
    DVec3 anchorLocal{};                   // body-local, double (rides large bodies exactly)
    Quaternion anchorRotPrev = { 0, 0, 0, 1 };
    DVec3 anchorEndPos{};

    // Tunables
    float thrust = 18.0f;                  // m/s^2
    float boostMultiplier = 6.0f;
    float turnRate = 1.6f;                 // rad/s (keys); mouse uses sensitivity
    float mouseSensitivity = 0.0025f;
    float dampingPerSec = 0.6f;            // flight assist
    float energyPerThrustSec = 1.2f;
    float lifeSupportRadius = 8.0f;
    float hullRadius = 2.2f;

    Vector3 forward() const { return Vector3RotateByQuaternion({ 0, 0, 1 }, orientation); }
    Vector3 up() const { return Vector3RotateByQuaternion({ 0, 1, 0 }, orientation); }
    Vector3 right() const { return Vector3RotateByQuaternion({ 1, 0, 0 }, orientation); }
    Vector3 renderPos(const Scene& world) const;

    void parkNear(Scene& world, Vector3 renderPoint, Vector3 surfaceUp); // place on the ground
    void update(Scene& world, bool inputEnabled, float dt);
    void applyChaseCamera(Scene& world) const;
    void draw(const Scene& world) const;   // inside BeginMode3D
    bool playerInRange(const Scene& world, Vector3 playerRender, float range) const;
};

extern Ship ship;

#endif
