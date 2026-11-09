#include "player.h"
#include "world.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>

Player player;


static Vector3 RotateAround(Vector3 v, Vector3 axis, float angle) {
    return Vector3RotateByQuaternion(v, QuaternionFromAxisAngle(axis, angle));
}

void Player::spawnOn(Scene& world, Planetoid& body, Vector3 dir) {
    const float R = body.radius * body.scale;
    const Vector3 c = body.position;
    if (body.flat) dir = { 0, 1, 0 };
    dir = Vector3Normalize(dir);
    const Vector3 above = body.flat ? Vector3Add(c, { 0, body.maxHeight + 20.0f, 0 }) : Vector3Add(c, Vector3Scale(dir, R + 60.0f));
    TerrainHit hit;
    Vector3 feet = Vector3Add(c, Vector3Scale(dir, R + 5.0f));
    if (RaycastTerrain(world, Ray{ above, Vector3Negate(dir) }, body.flat ? body.maxHeight - body.minHeight + 40.0f : R + 80.0f, hit)) feet = Vector3Add(hit.point, Vector3Scale(hit.normal, 0.2f));
    position = world.renderOrigin + feet;
    up = body.flat ? Vector3{ 0, 1, 0 } : Vector3Normalize(Vector3Subtract(feet, c));
    forward = Vector3Normalize(Vector3CrossProduct(up, { 1, 0, 0 }));
    if (Vector3Length(forward) < 0.1f) forward = { 0, 0, 1 };
    velocity = { 0, 0, 0 };
    grounded = false;
    groundBody = nullptr;
    flyMode = false;
}

void Player::update(Scene& world, bool inputEnabled, float dt) {
    PROFILE_FUNCTION();
    if (inputEnabled && IsKeyPressed(KEY_F)) {
        flyMode = !flyMode;
        velocity = { 0, 0, 0 };
        groundBody = nullptr;
    }
    impactSpeed = 0.0f; // per-frame; cleared before the fly-mode early return too
    if (flyMode) {
        // The free-fly camera (Scene / customUpdateCamera) drives; keep the body under the eye.
        position = world.renderOrigin + Vector3Subtract(world.camera.position, Vector3Scale(up, eyeHeight));
        forward = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        up = Vector3Normalize(world.camera.up); // keep a consistent basis for leaving fly mode
        return;
    }
    dt = std::min(dt, 0.1f); // don't explode after a hitch

    // Drop a ground anchor whose body was destroyed (world regen, a load) whether or not we are
    // still flagged grounded; everything below may dereference groundBody.
    if (groundBody && Planetoid::find(groundId) != groundBody) groundBody = nullptr;

    // Ride along with a spinning body we are standing on: position AND orientation turn with it,
    // so standing still on a spinning asteroid looks still (the view no longer drifts).
    if (grounded && groundBody) {
        if ((position - lastEndPosition).length() > 1e-6 || Planetoid::find(groundId) != groundBody) {
            // Something else moved us since last frame (respawn, leaving the ship, a load), or the
            // body was destroyed (world regen): the old ground anchor no longer applies.
            grounded = false;
            groundBody = nullptr;
        } else {
            position = world.renderOrigin + groundBody->toWorld(groundLocal);
            // Turn the view by however much the body turned since last frame (a delta, so view
            // changes made elsewhere between frames are kept).
            const Quaternion now = groundBody->rotationQuat();
            const Quaternion delta = QuaternionMultiply(now, QuaternionInvert(groundRotPrev));
            forward = Vector3Normalize(Vector3RotateByQuaternion(forward, delta));
            up = Vector3Normalize(Vector3RotateByQuaternion(up, delta));
        }
    }
    Vector3 feet = (position - world.renderOrigin).toFloat();

    nearest = NearestPlanetoid(world, feet);
    if (nearest) {
        const Vector3 newUp = nearest->flat ? Vector3{ 0, 1, 0 } : Vector3Normalize(Vector3Subtract(feet, nearest->position));
        // Carry the view along as "up" changes so walking around the body feels continuous. A zero
        // "up" (a save without one) has nothing to carry from: take the new up directly.
        if (Vector3LengthSqr(up) < 0.25f) up = newUp;
        forward = Vector3Normalize(Vector3RotateByQuaternion(forward, QuaternionFromVector3ToVector3(up, newUp)));
        up = newUp;
    } else if (Vector3LengthSqr(up) < 0.25f) {
        up = { 0, 1, 0 };
    }

    Vector3 wish = { 0, 0, 0 };
    bool jump = false, jetUp = false, jetDown = false;
    if (autoWalk) wish = Vector3Normalize(Vector3Subtract(forward, Vector3Scale(up, Vector3DotProduct(forward, up))));
    if (inputEnabled) {
        const Vector2 m = GetMouseDelta();
        forward = RotateAround(forward, up, -m.x * mouseSensitivity);
        const Vector3 right = Vector3Normalize(Vector3CrossProduct(forward, up));
        const Vector3 pitched = RotateAround(forward, right, -m.y * mouseSensitivity);
        if (fabsf(Vector3DotProduct(pitched, up)) < 0.985f) forward = pitched; // clamp ~10 deg from straight up/down

        const Vector3 flat = Vector3Normalize(Vector3Subtract(forward, Vector3Scale(up, Vector3DotProduct(forward, up))));
        const Vector3 side = Vector3CrossProduct(flat, up);
        // On the ground WASD walks level; in the air it fires the suit thrusters along the view.
        const Vector3 f = grounded ? flat : forward;
        if (IsKeyDown(KEY_W)) wish = Vector3Add(wish, f);
        if (IsKeyDown(KEY_S)) wish = Vector3Subtract(wish, f);
        if (IsKeyDown(KEY_D)) wish = Vector3Add(wish, side);
        if (IsKeyDown(KEY_A)) wish = Vector3Subtract(wish, side);
        if (Vector3Length(wish) > 0.0f) wish = Vector3Normalize(wish);
        if (IsKeyDown(KEY_LEFT_SHIFT)) wish = Vector3Scale(wish, runMultiplier);
        jump = IsKeyPressed(KEY_SPACE);
        jetUp = IsKeyDown(KEY_SPACE);
        jetDown = IsKeyDown(KEY_LEFT_CONTROL);
        if (IsKeyPressed(KEY_B)) dampeners = !dampeners;
    }
    if (autoJump) { jump = true; autoJump = false; }

    // Fixed substeps keep collision stable at any frame rate, and are shortened further at speed
    // (<= 0.2 m of travel each, under the body radius) so a fast player can't pass through thin walls.
    const float h = 1.0f / 120.0f;
    const int steps = std::clamp(std::max((int)ceilf(dt / h), (int)ceilf(Vector3Length(velocity) * dt / 0.2f)), 1, 1024);
    const float sdt = dt / steps;
    for (int i = 0; i < steps; ++i) { stepPhysics(world, sdt, wish, jump && i == 0, jetUp, jetDown); }

    feet = (position - world.renderOrigin).toFloat();
    if (grounded && nearest) {
        groundBody = nearest;
        groundId = nearest->id;
        groundLocal = nearest->toLocal(feet);
        groundRotPrev = nearest->rotationQuat();
    } else {
        groundBody = nullptr;
    }
    lastEndPosition = position;
}

void Player::stepPhysics(Scene& world, float dt, Vector3 wish, bool jump, bool jetUp, bool jetDown) {
    Vector3 feet = (position - world.renderOrigin).toFloat();

    const bool wasGrounded = grounded;

    // Gravity: inverse-square toward the nearest body, all the way out (no cutoff), so orbits
    // work. (Summing every body's pull is milestone 1 in the design doc.)
    Vector3 gdir = { 0, 0, 0 };
    gravity = 0.0f;
    if (nearest) {
        const Vector3 toC = Vector3Subtract(nearest->position, feet);
        const float dist = Vector3Length(toC);
        const float R = nearest->radius * nearest->scale;
        if (nearest->flat) {
            gdir = { 0, -1, 0 };
            gravity = 9.81f; // open world: ordinary gravity
        } else if (dist > 1e-3f) {
            gdir = Vector3Scale(toC, 1.0f / dist);
            const float g0 = SurfaceGravity(R);
            gravity = g0 * std::min(1.0f, (R * R) / (dist * dist));
        }
    }

    const float vUp = Vector3DotProduct(velocity, up);
    if (grounded) {
        // Direct ground control along the contact plane, pressed lightly into it so crests and dips
        // keep contact. Collision-induced upward speed is discarded (only jumps leave the ground).
        const Vector3 n = groundNormal;
        Vector3 along = Vector3Subtract(wish, Vector3Scale(n, Vector3DotProduct(wish, n)));
        const float wl = Vector3Length(wish), al = Vector3Length(along);
        if (al > 1e-4f) along = Vector3Scale(along, wl / al);
        velocity = Vector3Subtract(Vector3Scale(along, walkSpeed), Vector3Scale(n, 1.5f));
        (void)vUp;
        if (jump) { velocity = Vector3Add(velocity, Vector3Scale(up, jumpSpeed)); grounded = false; }
        jetpackFuel = std::min(1.0f, jetpackFuel + dt * 0.5f);
    } else {
        // Airborne: Newtonian. Suit thrusters push along the view (WASD) and up/down (Space/Ctrl)
        // while there is fuel; nothing clamps speed, so you can build up orbital velocity.
        Vector3 thrust = Vector3Length(wish) > 0.0f ? Vector3Normalize(wish) : wish;
        if (jetUp) thrust = Vector3Add(thrust, up);
        if (jetDown) thrust = Vector3Subtract(thrust, up);
        const float tl = Vector3Length(thrust);
        if (tl > 0.0f && jetpackFuel > 0.0f) {
            velocity = Vector3Add(velocity, Vector3Scale(thrust, jetpackAccel * dt / std::max(tl, 1.0f)));
            jetpackFuel = std::max(0.0f, jetpackFuel - dt * jetpackFuelPerSec);
        } else if (tl == 0.0f && dampeners && gravity < 0.5f) {
            // Inertial dampeners (B toggles) only in deep space, so they never fight a fall or an orbit.
            velocity = Vector3Scale(velocity, 1.0f - std::min(1.0f, 0.6f * dt));
        }
    }
    velocity = Vector3Add(velocity, Vector3Scale(gdir, gravity * gravityScale * dt));

    feet = Vector3Add(feet, Vector3Scale(velocity, dt));
    grounded = resolveCollisions(world, feet);
    // Leaving a spinning body: on the ground velocity is relative to the surface, so add the
    // surface's own motion; you keep the body's momentum, as you would really.
    if (wasGrounded && !grounded && groundBody && Planetoid::find(groundId) == groundBody) velocity = Vector3Add(velocity, groundBody->surfaceVelocity(feet));
    position = world.renderOrigin + feet;
}

// Push the feet/head spheres out of solid density, then probe just below the feet for ground.
bool Player::resolveCollisions(Scene&, Vector3& feet) {
    if (!nearest) return false;
    Planetoid& p = *nearest;
    bool onGround = false;
    const float offsets[2] = { radius, eyeHeight - radius };
    for (int iter = 0; iter < 4; ++iter) {
        bool moved = false;
        for (float off : offsets) {
            float dist; Vector3 n;
            if (!TerrainSurfaceDistance(p, Vector3Add(feet, Vector3Scale(up, off)), dist, n) || dist >= radius) continue;
            const float push = std::min(radius - dist, 1.0f); // limit per iteration if deeply embedded
            feet = Vector3Add(feet, Vector3Scale(n, push));
            const float vUpBefore = Vector3DotProduct(velocity, up);
            const float vn = Vector3DotProduct(velocity, n);
            if (vn < 0.0f) { impactSpeed = std::max(impactSpeed, -vn); velocity = Vector3Subtract(velocity, Vector3Scale(n, vn)); }
            lastContactDot = Vector3DotProduct(n, up);
            // Steep walls must not convert a push into upward motion (no wall climbing).
            if (lastContactDot <= 0.55f) {
                const float vUpAfter = Vector3DotProduct(velocity, up);
                if (vUpAfter > vUpBefore) velocity = Vector3Subtract(velocity, Vector3Scale(up, vUpAfter - vUpBefore));
            }
            if (lastContactDot > 0.55f) { onGround = true; groundNormal = n; } // walkable slope (< ~57 deg)
            moved = true;
        }
        if (!moved) break;
    }
    // Ground probe: the push-out is approximate and can leave a small gap; treat anything within a
    // few centimetres below the feet sphere as standing (unless moving away, e.g. jumping).
    if (!onGround && Vector3DotProduct(velocity, up) <= 0.5f) {
        float dist; Vector3 n;
        const float skin = 0.12f;
        if (TerrainSurfaceDistance(p, Vector3Add(feet, Vector3Scale(up, radius)), dist, n) && dist < radius + skin && Vector3DotProduct(n, up) > 0.55f) {
            feet = Vector3Subtract(feet, Vector3Scale(n, std::max(0.0f, dist - radius))); // settle onto it
            onGround = true;
            groundNormal = n;
        }
    }
    return onGround;
}

void Player::applyToCamera(Scene& world) const {
    if (flyMode) return;
    const Vector3 feet = (position - world.renderOrigin).toFloat();
    world.camera.position = Vector3Add(feet, Vector3Scale(up, eyeHeight));
    world.camera.target = Vector3Add(world.camera.position, forward);
    world.camera.up = up;
}
