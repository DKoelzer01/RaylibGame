#include "player.h"
#include "world.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdio>

Player player;


static Vector3 RotateAround(Vector3 v, Vector3 axis, float angle) {
    return Vector3RotateByQuaternion(v, QuaternionFromAxisAngle(axis, angle));
}

void Player::spawnOn(Scene& world, Planetoid& body, Vector3 dir) {
    const float R = body.radius * body.scale;
    const Vector3 c = body.position;
    if (body.flat) dir = { 0, 1, 0 };
    dir = Vector3Normalize(dir);
    const float outer = body.field->maxRadius() * body.scale; // above the highest terrain (planets: mountains)
    const Vector3 above = body.flat ? Vector3Add(c, { 0, body.maxHeight + 20.0f, 0 }) : Vector3Add(c, Vector3Scale(dir, outer + 60.0f));
    TerrainHit hit;
    Vector3 feet = Vector3Add(c, Vector3Scale(dir, R + 5.0f));
    if (RaycastTerrain(world, Ray{ above, Vector3Negate(dir) }, body.flat ? body.maxHeight - body.minHeight + 40.0f : outer + 80.0f, hit)) feet = Vector3Add(hit.point, Vector3Scale(hit.normal, 0.2f));
    position = world.renderOrigin + feet;
    up = body.flat ? Vector3{ 0, 1, 0 } : Vector3Normalize(Vector3Subtract(feet, c));
    forward = Vector3Normalize(Vector3CrossProduct(up, { 1, 0, 0 }));
    if (Vector3Length(forward) < 0.1f) forward = { 0, 0, 1 };
    velocity = body.surfaceVelocity(feet); // move with the ground (a spinning planet's surface is fast)
    grounded = false;
    groundBody = nullptr;
    anchorTo(body);
    flyMode = false;
}

void Player::update(Scene& world, bool inputEnabled, float dt) {
    PROFILE_FUNCTION();
    if (inputEnabled && IsKeyPressed(KEY_F)) {
        flyMode = !flyMode;
        velocity = { 0, 0, 0 };
        if (!flyMode) {
            // Drop in moving with the ground below (a spinning planet's surface is fast).
            const Vector3 feet = (position - world.renderOrigin).toFloat();
            if (Planetoid* p = NearestPlanetoid(world, feet)) velocity = p->surfaceVelocity(feet);
        }
        groundBody = nullptr;
    }
    impactSpeed = 0.0f; // per-frame; cleared before the fly-mode early return too
    jetpackFiring = false;
    if (flyMode) {
        // Carry the camera with its reference body's motion since recordFly() (after last frame's input).
        bool rotate = false;
        Planetoid* ref = flyFollow ? flyReference(world, rotate) : nullptr;
        if (ref && ref->id == flyBodyId && rotate == flyRotate) {
            if (rotate) {
                flyCarried = ref->toWorldD(flyLocal);
                flyCarriedValid = true;
                world.camera.position = (flyCarried - world.renderOrigin).toFloat();
                world.camera.target = Vector3Add(world.camera.position, ref->dirToWorld(flyLocalFwd));
                world.camera.up = ref->dirToWorld(flyLocalUp);
            } else {
                const Vector3 move = (ref->worldPos - flyBodyPos).toFloat();
                world.camera.position = Vector3Add(world.camera.position, move);
                world.camera.target = Vector3Add(world.camera.target, move);
            }
        }
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
            position = groundBody->toWorldD(groundLocal);
            // Turn the view by however much the body turned since last frame (a delta, so view
            // changes made elsewhere between frames are kept).
            const Quaternion now = groundBody->rotationQuat();
            const Quaternion delta = QuaternionMultiply(now, QuaternionInvert(groundRotPrev));
            forward = Vector3Normalize(Vector3RotateByQuaternion(forward, delta));
            up = Vector3Normalize(Vector3RotateByQuaternion(up, delta));
        }
    }
    // Airborne: carry along with the body we were near (its spin and orbit since last frame); the
    // motion of our own (world-frame) velocity is then integrated relative to its surface below.
    if (airBody && Planetoid::find(airId) != airBody) airBody = nullptr;
    if (!grounded && airBody) {
        if ((position - lastEndPosition).length() > 1e-6) airBody = nullptr; // moved from outside (respawn, load, ship)
        else position = airBody->toWorldD(airLocal);
    }
    Vector3 feet = (position - world.renderOrigin).toFloat();

    nearest = NearestPlanetoid(world, feet);
    // "Up" is against the summed gravity. Where gravity is too weak to define it (deep space, the
    // balance point between bodies) it is left alone so the view doesn't spin; in the air it turns
    // at most 2 rad/s. On the ground it follows exactly, so walking around a body stays continuous.
    {
        // Felt gravity: near a body on an orbit, it and everything on it fall freely around its parent,
        // so subtract the body's own acceleration (the outside gravity at its centre). Otherwise the
        // parent's pull would tilt "up" (~4 deg on the moon) and jumps would go out sideways.
        Vector3 g = GravityAt(world, feet);
        if (nearest && nearest->orbitParent) g = Vector3Subtract(g, GravityAt(world, nearest->position));
        const float gl = Vector3Length(g);
        Vector3 newUp = up;
        if (nearest && nearest->flat) newUp = { 0, 1, 0 };
        else if (gl > 0.05f) newUp = Vector3Scale(g, -1.0f / gl);
        if (Vector3LengthSqr(up) < 0.25f) {
            // A zero "up" (a save without one) has nothing to carry from: take the new up directly.
            up = Vector3LengthSqr(newUp) > 0.25f ? newUp : Vector3{ 0, 1, 0 };
            newUp = up;
        }
        if (!grounded) {
            const float angle = acosf(std::clamp(Vector3DotProduct(up, newUp), -1.0f, 1.0f));
            const float maxTurn = 2.0f * dt;
            if (angle > maxTurn)
                newUp = Vector3Normalize(Vector3RotateByQuaternion(up, QuaternionSlerp(QuaternionIdentity(), QuaternionFromVector3ToVector3(up, newUp), maxTurn / angle)));
        }
        // Carry the view along as "up" changes so walking around the body feels continuous.
        forward = Vector3Normalize(Vector3RotateByQuaternion(forward, QuaternionFromVector3ToVector3(up, newUp)));
        up = newUp;
    }

    Vector3 wish = { 0, 0, 0 };
    bool jump = false, jetUp = false, jetDown = false;
    if (autoWalk) wish = Vector3Normalize(Vector3Subtract(forward, Vector3Scale(up, Vector3DotProduct(forward, up))));
    if (inputEnabled) {
        const Vector2 m = GetMouseDelta();
        forward = RotateAround(forward, up, -m.x * mouseSensitivity * lookScale);
        const Vector3 right = Vector3Normalize(Vector3CrossProduct(forward, up));
        const Vector3 pitched = RotateAround(forward, right, -m.y * mouseSensitivity * lookScale);
        if (fabsf(Vector3DotProduct(pitched, up)) < 0.985f) forward = pitched; // clamp ~10 deg from straight up/down

        const Vector3 flat = Vector3Normalize(Vector3Subtract(forward, Vector3Scale(up, Vector3DotProduct(forward, up))));
        const Vector3 side = Vector3CrossProduct(flat, up);
        // On the ground and in the air under gravity WASD is level (walking / air control); in deep
        // space it fires the suit thrusters along the view.
        const Vector3 f = (grounded || gravity >= 0.5f) ? flat : forward;
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
    // Jump buffer: a press is remembered briefly, so it is not lost on a frame where "grounded"
    // flickers off over bumpy ground (see coyote time in stepPhysics).
    if (jump) jumpBuffer = 0.15f;

    // Fixed substeps keep collision stable at any frame rate, and are shortened further at speed
    // (<= 0.2 m of travel each, under the body radius) so a fast player can't pass through thin walls.
    const float h = 1.0f / 120.0f;
    // Speed relative to the ground (airborne velocity is world-frame and carries the spin).
    const Vector3 relVel = (!grounded && nearest) ? Vector3Subtract(velocity, nearest->surfaceVelocity((position - world.renderOrigin).toFloat())) : velocity;
    const int steps = std::clamp(std::max((int)ceilf(dt / h), (int)ceilf(Vector3Length(relVel) * dt / 0.2f)), 1, 1024);
    const float sdt = dt / steps;
    for (int i = 0; i < steps; ++i) { stepPhysics(world, sdt, wish, jetUp, jetDown); }

    feet = (position - world.renderOrigin).toFloat();
    if (grounded && nearest) {
        groundBody = nearest;
        groundId = nearest->id;
        groundLocal = nearest->toLocalD(position);
        groundRotPrev = nearest->rotationQuat();
    } else {
        groundBody = nullptr;
    }
    if (!grounded && nearest) {
        airBody = nearest;
        airId = nearest->id;
        airLocal = nearest->toLocalD(position);
    } else {
        airBody = nullptr;
    }
    lastEndPosition = position;
}

void Player::stepPhysics(Scene& world, float dt, Vector3 wish, bool jetUp, bool jetDown) {
    Vector3 feet = (position - world.renderOrigin).toFloat();

    const bool wasGrounded = grounded;

    // Gravity summed over every body (inverse-square, no cutoff, so orbits work).
    const Vector3 gvec = GravityAt(world, feet);
    gravity = Vector3Length(gvec);
    const Vector3 gdir = gravity > 1e-6f ? Vector3Scale(gvec, 1.0f / gravity) : Vector3{ 0, 0, 0 };

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
        // A jump starts from the walking velocity only: keeping the press into a sloped ground would
        // push you downhill on takeoff.
        if (jumpBuffer > 0.0f) {
            velocity = Vector3Add(Vector3Scale(along, walkSpeed), Vector3Scale(up, jumpSpeed));
            grounded = false;
            jumpBuffer = coyote = 0.0f;
            jumpLock = 0.15f;
        } else {
            coyote = 0.12f;
        }
        jetpackFuel = std::min(1.0f, jetpackFuel + dt * 0.5f);
    } else {
        // Airborne: Newtonian; nothing clamps speed, so you can build up orbital velocity.
        const bool inGravity = gravity >= 0.5f;
        if (jumpBuffer > 0.0f && coyote > 0.0f) {
            // Coyote time: just left the ground (or "grounded" flickered off for a frame) - still jump,
            // setting the upward speed relative to the surface below.
            const Vector3 surface = nearest ? nearest->surfaceVelocity(feet) : Vector3{ 0, 0, 0 };
            velocity = Vector3Add(velocity, Vector3Scale(up, jumpSpeed - Vector3DotProduct(Vector3Subtract(velocity, surface), up)));
            jumpBuffer = coyote = 0.0f;
            jumpLock = 0.15f;
        }
        Vector3 thrust = { 0, 0, 0 };
        if (inGravity) {
            // Air control (free, no fuel): steer toward the wished walking velocity relative to the
            // ground below - enough to reach a ledge mid-jump. With no input, momentum is untouched.
            if (Vector3Length(wish) > 0.0f) {
                const Vector3 surface = nearest ? nearest->surfaceVelocity(feet) : Vector3{ 0, 0, 0 };
                const Vector3 rel = Vector3Subtract(velocity, surface);
                const Vector3 relLevel = Vector3Subtract(rel, Vector3Scale(up, Vector3DotProduct(rel, up)));
                const Vector3 dv = Vector3Subtract(Vector3Scale(wish, walkSpeed), relLevel);
                const float dvl = Vector3Length(dv), maxDv = airControlAccel * dt;
                if (dvl > 1e-5f) velocity = Vector3Add(velocity, Vector3Scale(dv, std::min(1.0f, maxDv / dvl)));
            }
        } else if (Vector3Length(wish) > 0.0f) {
            thrust = Vector3Normalize(wish); // deep space: suit thrusters along the view
        }
        // Jetpack (Space up / Ctrl down, and WASD in deep space) burns fuel; recharging it on the ground
        // is what costs suit power (see the survival code).
        if (jetUp) thrust = Vector3Add(thrust, up);
        if (jetDown) thrust = Vector3Subtract(thrust, up);
        const float tl = Vector3Length(thrust);
        if (tl > 0.0f && jetpackFuel > 0.0f) {
            velocity = Vector3Add(velocity, Vector3Scale(thrust, jetpackAccel * dt / std::max(tl, 1.0f)));
            jetpackFuel = std::max(0.0f, jetpackFuel - dt * jetpackFuelPerSec);
            jetpackFiring = true;
        } else if (tl == 0.0f && dampeners && !inGravity) {
            // Inertial dampeners (B toggles) only in deep space, so they never fight a fall or an orbit.
            velocity = Vector3Scale(velocity, 1.0f - std::min(1.0f, 0.6f * dt));
        }
    }
    Vector3 gAccel = Vector3Scale(gdir, gravity * gravityScale * dt);
    // Static friction: standing on walkable ground, only gravity's push into the ground applies.
    // (The downhill part made a standing player creep ~1.6 cm/s down slopes at fixed substeps.)
    if (grounded) gAccel = Vector3Scale(groundNormal, std::min(0.0f, Vector3DotProduct(gAccel, groundNormal)));
    velocity = Vector3Add(velocity, gAccel);

    // Airborne, velocity is world-frame but the terrain already stands where it will be at the end of
    // this frame (and update() carried us with it): move relative to its surface.
    const Vector3 frameVel = (!wasGrounded && nearest) ? nearest->surfaceVelocity(feet) : Vector3{ 0, 0, 0 };
    feet = Vector3Add(feet, Vector3Scale(Vector3Subtract(velocity, frameVel), dt));
    collideFrameVel = (wasGrounded || !nearest) ? Vector3{ 0, 0, 0 } : nearest->surfaceVelocity(feet);
    grounded = resolveCollisions(world, feet);
    collideTrunks(feet);
    // Just jumped: the probe below the feet would still find the slope (they have risen only a few cm)
    // and re-ground the player, cancelling the jump. Collisions still push out; only grounding waits.
    if (jumpLock > 0.0f) { grounded = false; jumpLock = std::max(0.0f, jumpLock - dt); }
    // Just landed: on the ground velocity is relative to the surface (as the walking code expects), so
    // take the ground's own motion out now rather than leave a world-frame velocity for a frame.
    if (!wasGrounded && grounded) velocity = Vector3Subtract(velocity, collideFrameVel);
    // Leaving a spinning body: on the ground velocity is relative to the surface, so add the
    // surface's own motion; you keep the body's momentum, as you would really.
    if (wasGrounded && !grounded && groundBody && Planetoid::find(groundId) == groundBody) velocity = Vector3Add(velocity, groundBody->surfaceVelocity(feet));
    position = world.renderOrigin + feet;
    jumpBuffer = std::max(0.0f, jumpBuffer - dt);
    if (!grounded) coyote = std::max(0.0f, coyote - dt);
}

// Push the feet/head spheres out of solid density, then probe just below the feet for ground.
bool Player::resolveCollisions(Scene&, Vector3& feet) {
    if (!nearest) return false;
    Planetoid& p = *nearest;
    bool onGround = false;
    const float offsets[2] = { radius, eyeHeight - radius };
    // The linear estimate (iso - density) / |gradient| is only a distance where the density behaves
    // like one; dig and build edits (big density steps) break that, and it reported contacts in open
    // air that pushed up to 1 m per iteration, sometimes into the ground (the "wall launch"). So each
    // contact is confirmed by sampling along the normal, the push is the measured penetration, and the
    // total push in a frame is capped.
    auto solid = [&](Vector3 c) { return p.densityAt(p.toLocal(c)) >= ISO_LEVEL; };
    float pushBudget = 0.6f;
    for (int iter = 0; iter < 4 && pushBudget > 0.0f; ++iter) {
        bool moved = false;
        for (float off : offsets) {
            float dist; Vector3 n;
            const Vector3 c = Vector3Add(feet, Vector3Scale(up, off));
            if (!TerrainSurfaceDistance(p, c, dist, n) || dist >= radius) continue;
            float push = 0.0f;
            if (!solid(c)) {
                // Centre in the open: is there solid ground within the radius along -n, and how far?
                if (!solid(Vector3Subtract(c, Vector3Scale(n, radius)))) continue; // a false contact
                float lo = 0.0f, hi = radius;                                       // lo open, hi solid
                for (int k = 0; k < 6; ++k) {
                    const float mid = 0.5f * (lo + hi);
                    if (solid(Vector3Subtract(c, Vector3Scale(n, mid)))) hi = mid; else lo = mid;
                }
                push = radius - hi;
            } else {
                // Embedded: step out along the normal to open space.
                float s = 0.05f;
                while (s < 1.0f && solid(Vector3Add(c, Vector3Scale(n, s)))) s += 0.05f;
                push = s + radius;
            }
            push = std::min(push, pushBudget);
            if (push <= 1e-4f) continue;
            pushBudget -= push;
            feet = Vector3Add(feet, Vector3Scale(n, push));
            const float vUpBefore = Vector3DotProduct(velocity, up);
            // Impact and push-out relative to the (moving) terrain: against the world-frame velocity a
            // spinning planet's hole walls hit like 86 m/s impacts and stopped you dead while they moved on.
            const float vn = Vector3DotProduct(Vector3Subtract(velocity, collideFrameVel), n);
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
    // Build plot parts: the same two spheres pushed out of part shapes (exact convex cells).
    if (p.plots) {
        for (int iter = 0; iter < 2; ++iter)
            for (float off : offsets)
                for (auto& pl : p.plots->plots) {
                    const Vector3 before = p.toLocal(Vector3Add(feet, Vector3Scale(up, off)));
                    Vector3 c = before, nl;
                    if (!pl->collideSphere(c, radius / p.scale, nl)) continue;
                    const Vector3 moved = Vector3Subtract(c, before);
                    const float len = Vector3Length(moved) * p.scale;
                    if (len > 1e-6f) feet = Vector3Add(feet, Vector3Scale(p.dirToWorld(Vector3Normalize(moved)), len));
                    const Vector3 n = p.dirToWorld(nl);
                    const float vn = Vector3DotProduct(Vector3Subtract(velocity, collideFrameVel), n);
                    if (vn < 0.0f) { impactSpeed = std::max(impactSpeed, -vn); velocity = Vector3Subtract(velocity, Vector3Scale(n, vn)); }
                    if (Vector3DotProduct(n, up) > 0.55f) { onGround = true; groundNormal = n; }
                }
        // Standing on a part: a probe just below the feet sphere.
        if (!onGround && Vector3DotProduct(velocity, up) <= 0.5f)
            for (auto& pl : p.plots->plots) {
                Vector3 c = p.toLocal(Vector3Add(feet, Vector3Scale(up, radius - 0.12f))), nl;
                if (pl->collideSphere(c, radius / p.scale, nl) && Vector3DotProduct(p.dirToWorld(nl), up) > 0.55f) { onGround = true; groundNormal = p.dirToWorld(nl); break; }
            }
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

// Tree trunks are vertical capsules (bushes are soft): push the body cylinder out sideways and
// drop the velocity into the trunk, relative to the (moving) ground the tree stands on.
void Player::collideTrunks(Vector3& feet) {
    if (!nearest || !nearest->trees) return;
    static thread_local std::vector<Trunk> near;
    near.clear();
    const Vector3 local = nearest->toLocal(feet);
    nearest->trees->trunksNear(local, radius + 2.0f, near);
    for (const Trunk& t : near) {
        const Vector3 d = Vector3Subtract(local, t.base);
        const float along = Vector3DotProduct(d, t.up);
        if (along < -eyeHeight || along > t.height) continue;      // below the roots or over the top
        const Vector3 side = Vector3Subtract(d, Vector3Scale(t.up, along));
        const float dist = Vector3Length(side);
        const float minDist = t.radius + radius;
        if (dist >= minDist || dist < 1e-4f) continue;
        const Vector3 n = nearest->dirToWorld(Vector3Scale(side, 1.0f / dist));
        feet = Vector3Add(feet, Vector3Scale(n, minDist - dist));
        const float vn = Vector3DotProduct(Vector3Subtract(velocity, collideFrameVel), n);
        if (vn < 0.0f) velocity = Vector3Subtract(velocity, Vector3Scale(n, vn));
    }
}

Planetoid* Player::flyReference(Scene& world, bool& rotate) const {
    // Sphere of influence: a = orbit radius, r = a (m / M)^0.4 with masses ~ g R^2. Bodies with no
    // parent reach everywhere. Several hold the camera: the nearest surface wins.
    const DVec3 cam = world.renderOrigin + world.camera.position;
    Planetoid* best = nullptr;
    double bestDist = 1e300;
    for (auto& obj : world.objects) {
        Planetoid* b = dynamic_cast<Planetoid*>(obj.get());
        if (!b || b->flat) continue;
        const double d = (cam - b->worldPos).length();
        const double R = (double)b->radius * b->scale;
        double soi = 1e300;
        if (b->orbitParent && Planetoid::find(b->orbitParentId) == b->orbitParent) {
            const Planetoid* m = b->orbitParent;
            const double Rm = (double)m->radius * m->scale;
            const double ratio = (b->surfaceGravity * R * R) / std::max(1e-9, (double)m->surfaceGravity * Rm * Rm);
            soi = std::max(2.0 * R, b->orbitRadius * pow(std::max(ratio, 0.0), 0.4));
        } else if (b->kind == BodyKind::Asteroid) {
            soi = 10.0 * R; // loose asteroids (no orbit): their neighbourhood
        }
        if (d > soi) continue;
        if (d - R < bestDist) { bestDist = d - R; best = b; }
    }
    if (best) rotate = (cam - best->worldPos).length() < 2.0 * best->radius * best->scale;
    return best;
}

void Player::recordFly(Scene& world) {
    flyBodyId = 0;
    if (!flyMode || !flyFollow) return;
    bool rotate = false;
    Planetoid* ref = flyReference(world, rotate);
    if (!ref) return;
    flyBodyId = ref->id;
    flyRotate = rotate;
    flyBodyPos = ref->worldPos;
    // Unmoved by input since the carry: the double pose, not its float rounding (which drifted ~0.2 m in 3 s).
    const DVec3 cam = world.renderOrigin + world.camera.position;
    flyLocal = ref->toLocalD(flyCarriedValid && (cam - flyCarried).length() < 0.01 ? flyCarried : cam);
    flyCarriedValid = false;
    flyLocalFwd = ref->dirToLocal(Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position)));
    flyLocalUp = ref->dirToLocal(world.camera.up);
}

void Player::anchorTo(Planetoid& body) {
    airBody = &body;
    airId = body.id;
    airLocal = body.toLocalD(position);
    lastEndPosition = position;
}

void Player::placeAt(Scene& world, Planetoid& body, Vector3 renderFeet) {
    position = world.renderOrigin + renderFeet;
    velocity = body.surfaceVelocity(renderFeet);
    grounded = false;
    groundBody = nullptr;
    flyMode = false;
    anchorTo(body);
}
