#include "ship.h"
#include "world.h"

#include <rlgl.h>
#include <algorithm>
#include <cmath>

Ship ship;

Vector3 Ship::renderPos(const Scene& world) const { return (position - world.renderOrigin).toFloat(); }

void Ship::parkNear(Scene& world, Vector3 point, Vector3 surfaceUp) {
    // Drop onto the terrain below `point` so the hull rests on the surface.
    Vector3 rest = Vector3Add(point, Vector3Scale(surfaceUp, hullRadius));
    TerrainHit hit;
    if (RaycastTerrain(world, Ray{ Vector3Add(point, Vector3Scale(surfaceUp, 30.0f)), Vector3Negate(surfaceUp) }, 60.0f, hit))
        rest = Vector3Add(hit.point, Vector3Scale(hit.normal, hullRadius));
    position = world.renderOrigin + rest;
    orientation = QuaternionFromVector3ToVector3({ 0, 1, 0 }, surfaceUp);
    velocity = { 0, 0, 0 };
    spawned = true;
}

bool Ship::playerInRange(const Scene& world, Vector3 playerRender, float range) const {
    return spawned && Vector3Distance(playerRender, renderPos(world)) < range;
}

void Ship::update(Scene& world, bool inputEnabled, float dt) {
    if (!spawned) return;
    PROFILE_FUNCTION();
    dt = std::min(dt, 0.1f);

    // Parked on a body: ride its spin (position and orientation), skip flight physics.
    if (anchorBody) {
        if ((position - anchorEndPos).length() > 1e-6 || Planetoid::find(anchorId) != anchorBody) {
            anchorBody = nullptr; // moved by something else (parkNear, a load) or the body is gone
        } else if (!piloted || !(autopilot || (inputEnabled && (IsKeyDown(KEY_W) || IsKeyDown(KEY_S) || IsKeyDown(KEY_A) || IsKeyDown(KEY_D) || IsKeyDown(KEY_SPACE) || IsKeyDown(KEY_LEFT_CONTROL))))) {
            position = world.renderOrigin + anchorBody->toWorld(anchorLocal);
            const Quaternion now = anchorBody->rotationQuat();
            orientation = QuaternionNormalize(QuaternionMultiply(QuaternionMultiply(now, QuaternionInvert(anchorRotPrev)), orientation));
            anchorRotPrev = now;
            anchorEndPos = position;
            velocity = { 0, 0, 0 };
            return;
        } else {
            // Taking off from a spinning body: start with the surface's own motion.
            velocity = anchorBody->surfaceVelocity(renderPos(world));
            anchorBody = nullptr;
        }
    }
    Vector3 accel = { 0, 0, 0 };
    bool braking = false;

    if (piloted && inputEnabled) {
        const Vector2 m = GetMouseDelta();
        Quaternion q = orientation;
        q = QuaternionMultiply(q, QuaternionFromAxisAngle({ 0, 1, 0 }, -m.x * mouseSensitivity));
        q = QuaternionMultiply(q, QuaternionFromAxisAngle({ 1, 0, 0 }, m.y * mouseSensitivity));
        float roll = 0.0f;
        if (IsKeyDown(KEY_Z)) roll -= turnRate * dt;
        if (IsKeyDown(KEY_C)) roll += turnRate * dt;
        q = QuaternionMultiply(q, QuaternionFromAxisAngle({ 0, 0, 1 }, roll));
        orientation = QuaternionNormalize(q);

        Vector3 dir = { 0, 0, 0 };
        if (IsKeyDown(KEY_W)) dir = Vector3Add(dir, forward());
        if (IsKeyDown(KEY_S)) dir = Vector3Subtract(dir, forward());
        if (IsKeyDown(KEY_D)) dir = Vector3Subtract(dir, right());
        if (IsKeyDown(KEY_A)) dir = Vector3Add(dir, right());
        if (IsKeyDown(KEY_SPACE)) dir = Vector3Add(dir, up());
        if (IsKeyDown(KEY_LEFT_CONTROL)) dir = Vector3Subtract(dir, up());
        braking = IsKeyDown(KEY_X);
        if (Vector3Length(dir) > 0.0f && energy > 0.0f) {
            const float boost = IsKeyDown(KEY_LEFT_SHIFT) ? boostMultiplier : 1.0f;
            accel = Vector3Scale(Vector3Normalize(dir), thrust * boost);
            energy = std::max(0.0f, energy - energyPerThrustSec * boost * dt * 0.5f);
        }
    }
    if (autopilot && piloted) {
        // Test/assist: turn toward the target and fly there, braking on approach.
        const Vector3 to = Vector3Subtract((autopilotTarget - world.renderOrigin).toFloat(), renderPos(world));
        const float dist = Vector3Length(to);
        const Vector3 want = Vector3Scale(to, 1.0f / std::max(dist, 1e-3f));
        orientation = QuaternionSlerp(orientation, QuaternionMultiply(QuaternionFromVector3ToVector3(forward(), want), orientation), std::min(1.0f, 2.0f * dt));
        const float speed = Vector3Length(velocity);
        const float stopDist = speed * speed / (2.0f * thrust) + 25.0f;
        // Same energy rule as manual flight: no energy, no thrust.
        if (dist > stopDist && Vector3DotProduct(forward(), want) > 0.95f) {
            if (energy > 0.0f) {
                accel = Vector3Scale(forward(), thrust * boostMultiplier);
                energy = std::max(0.0f, energy - energyPerThrustSec * boostMultiplier * dt * 0.5f);
            }
        } else braking = true;
    }

    velocity = Vector3Add(velocity, Vector3Scale(accel, dt));
    // Flight assist: gentle damping always, strong when braking or parked.
    const float damp = braking ? 2.5f : piloted ? dampingPerSec * (Vector3Length(accel) > 0 ? 0.2f : 1.0f) : 4.0f;
    velocity = Vector3Scale(velocity, std::max(0.0f, 1.0f - damp * dt));

    Vector3 p = renderPos(world);

    // Move in substeps of at most ~1 m (less than the 1.2 m hull clearance) so a fast ship can't
    // tunnel through a small asteroid between frames; each substep pushes the hull samples out of
    // the density field.
    const int moveSteps = std::clamp((int)ceilf(Vector3Length(velocity) * dt / 1.0f), 1, 256);
    const float h = dt / (float)moveSteps;
    for (int step = 0; step < moveSteps; ++step) {
        p = Vector3Add(p, Vector3Scale(velocity, h));
        Planetoid* body = NearestPlanetoid(world, p);
        if (!body) continue;
        const Vector3 samples[3] = { p, Vector3Add(p, Vector3Scale(forward(), 2.2f)), Vector3Subtract(p, Vector3Scale(forward(), 2.2f)) };
        for (int iter = 0; iter < 3; ++iter) {
            for (const Vector3& s0 : samples) {
                const Vector3 s = Vector3Add(s0, Vector3Subtract(p, samples[0]));
                float dist; Vector3 n;
                if (!TerrainSurfaceDistance(*body, s, dist, n) || dist >= 1.2f) continue;
                p = Vector3Add(p, Vector3Scale(n, std::min(1.2f - dist, 2.0f)));
                const float vn = Vector3DotProduct(velocity, n);
                if (vn < 0.0f) {
                    lastImpact = std::max(lastImpact, -vn);
                    velocity = Vector3Subtract(velocity, Vector3Scale(n, 1.2f * vn)); // small bounce
                    velocity = Vector3Scale(velocity, 0.7f);                           // scrape friction
                }
            }
        }
    }
    position = world.renderOrigin + p;

    // Come to rest on the ground (parked, or landed with no thrust) -> anchor to that body so it
    // carries us with its spin; thrusting releases it.
    if (Vector3Length(accel) == 0.0f && Vector3Length(velocity) < 0.5f) {
        if (Planetoid* body = NearestPlanetoid(world, p)) {
            float dist; Vector3 n;
            if (TerrainSurfaceDistance(*body, p, dist, n) && dist < hullRadius + 1.0f) {
                anchorBody = body;
                anchorId = body->id;
                anchorLocal = body->toLocal(p);
                anchorRotPrev = body->rotationQuat();
                anchorEndPos = position;
                velocity = { 0, 0, 0 };
            }
        }
    }
}

void Ship::applyChaseCamera(Scene& world) const {
    const Vector3 p = renderPos(world);
    world.camera.position = Vector3Add(Vector3Subtract(p, Vector3Scale(forward(), 14.0f)), Vector3Scale(up(), 4.5f));
    world.camera.target = Vector3Add(p, Vector3Scale(forward(), 8.0f));
    world.camera.up = up();
}

void Ship::draw(const Scene& world) const {
    if (!spawned) return;
    const Vector3 p = renderPos(world);
    const Matrix rot = QuaternionToMatrix(orientation);
    rlPushMatrix();
    rlTranslatef(p.x, p.y, p.z);
    rlMultMatrixf(MatrixToFloat(rot));
    // Simple kit-bashed hull (local +Z is forward, +Y up).
    DrawCube({ 0, 0, 0 }, 2.4f, 1.2f, 5.0f, { 150, 155, 165, 255 });
    DrawCubeWires({ 0, 0, 0 }, 2.4f, 1.2f, 5.0f, { 60, 60, 70, 255 });
    DrawCube({ 0, 0.75f, 1.0f }, 1.4f, 0.7f, 1.8f, { 90, 170, 230, 255 });       // canopy
    DrawCube({ 0, -0.1f, -0.6f }, 6.5f, 0.18f, 1.8f, { 110, 115, 125, 255 });    // wings
    DrawCubeWires({ 0, -0.1f, -0.6f }, 6.5f, 0.18f, 1.8f, { 50, 50, 60, 255 });
    const bool firing = piloted && (IsKeyDown(KEY_W) || autopilot);
    const Color glow = firing ? Color{ 255, 170, 60, 255 } : Color{ 120, 80, 50, 255 };
    DrawCube({ -0.8f, 0, -2.7f }, 0.7f, 0.7f, 0.5f, glow);                      // engines
    DrawCube({ 0.8f, 0, -2.7f }, 0.7f, 0.7f, 0.5f, glow);
    DrawCube({ -0.9f, -0.9f, 1.4f }, 0.15f, 0.8f, 0.15f, DARKGRAY);              // landing legs
    DrawCube({ 0.9f, -0.9f, 1.4f }, 0.15f, 0.8f, 0.15f, DARKGRAY);
    DrawCube({ 0, -0.9f, -1.8f }, 0.15f, 0.8f, 0.15f, DARKGRAY);
    rlPopMatrix();
    // Life-support field hint when the player is outside.
    if (!piloted) DrawSphereWires(p, lifeSupportRadius, 6, 12, Fade(SKYBLUE, 0.15f));
}
