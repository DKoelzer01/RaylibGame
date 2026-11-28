#include "camera.h" // gFlyGravityUp (#232)
#include "player.h"
#include <chrono>
#include "world.h"
#include "playerinput.h"

#include <raymath.h>
#include <algorithm>
#include <cmath>
#include <cstdio>

// (player: the local PlayerSession's, src/space/playersession.cpp, #313)


static Vector3 RotateAround(Vector3 v, Vector3 axis, float angle) {
    return Vector3RotateByQuaternion(v, QuaternionFromAxisAngle(axis, angle));
}

void Player::spawnOn(Scene& world, Planetoid& body, Vector3 dir) {
    ++spawns;
    const float R = body.radius * body.scale;
    const Vector3 c = body.position;
    if (body.flat) dir = { 0, 1, 0 };
    dir = Vector3Normalize(dir);
    const float outer = body.field->maxRadius() * body.scale; // above the highest terrain (planets: mountains)
    const Vector3 above = body.flat ? Vector3Add(c, { 0, body.maxHeight + 20.0f, 0 }) : Vector3Add(c, Vector3Scale(dir, outer + 60.0f));
    TerrainHit hit;
    Vector3 feet = Vector3Add(c, Vector3Scale(dir, R + 5.0f));
    bool found = false;
    if (!body.flat) {
        // (#2000) The surface from the field itself, edits included: march down from maxRadius in 2 m steps, then halve to
        // ~1 cm. A mesh raycast only sees what is built so far: with no mesh there (a far body) R + 5 m could be a deadly
        // fall or underground, and a coarse LOD mesh sits up to ~2 m under the real ground, so the player stood inside it
        // when the fine chunks came and was thrown out ("fall" deaths on Luna, runner).
        const DensityField& f = *body.field;
        const TerrainDeltas* ed = body.edits.get();
        const Vector3 local = body.dirToLocal(dir); // the field is in the body's own (turning) axes
        const float top = f.maxRadius(), bottom = std::max(0.0f, f.minRadius());
        for (float r = top; r > bottom; r -= 2.0f) {
            if (f.density(Vector3Scale(local, r - 2.0f), ed) < ISO_LEVEL) continue;
            float lo = r - 2.0f, hi = r; // lo solid, hi air
            for (int i = 0; i < 8; ++i) { const float m = 0.5f * (lo + hi); (f.density(Vector3Scale(local, m), ed) >= ISO_LEVEL ? lo : hi) = m; }
            feet = Vector3Add(c, Vector3Scale(dir, hi * body.scale + 0.2f));
            found = true;
            break;
        }
    }
    if (!found && RaycastTerrain(world, Ray{ above, Vector3Negate(dir) }, body.flat ? body.maxHeight - body.minHeight + 40.0f : outer + 80.0f, hit))
        feet = Vector3Add(hit.point, Vector3Scale(hit.normal, 0.2f));
    position = world.renderOrigin + feet;
    up = body.flat ? Vector3{ 0, 1, 0 } : Vector3Normalize(Vector3Subtract(feet, c));
    forward = Vector3Normalize(Vector3CrossProduct(up, { 1, 0, 0 }));
    if (Vector3Length(forward) < 0.1f) forward = { 0, 0, 1 };
    velocity = body.surfaceVelocity(feet); // move with the ground (a spinning planet's surface is fast)
    grounded = false;
    groundBody = nullptr;
    // (#2000) The body put on is the nearest now, before the next update finds it: the spawn's loading screen holds the
    // player in player.nearest's frame, and the old one (Terra, for a move to Luna) held them still against Terra while Luna
    // moved on at ~1 km/s: released, its ground hit them at hundreds of m/s, a fatal "fall".
    nearest = &body;
    anchorTo(body);
    flyMode = false;
}


// Felt gravity: near a body on an orbit, it and everything on it fall freely around its parent, so
// its parent's pull cancels in its frame - subtract the body's own acceleration (the outside
// gravity at its centre). A small body only counts once its own pull here is felt (0.1 m/s^2):
// flying past ring rocks the planet stays the reference instead of "up" snapping between rocks.
// (F8b-4b) These take world points in double: the player's float space is its own (Player::spaceOrigin).
static bool FeltIsRelative(const DVec3& feet, const Planetoid* nearest) {
    return nearest && nearest->orbitParent && BodyGravityAtD(*nearest, feet).length() >= 0.1;
}
// (#2149) Another body's surface motion at `at`, in the terms of the player's velocity. That velocity is measured against
// the nearest body's ground less `nearestFrameVel` (collideFrameVel: the nearest's surface motion when it is world-frame,
// airborne off the frame path; 0 when it is ground-relative, grounded or on the frame path). Passing the other body's
// world-frame motion against a ground-relative velocity read a ring rock's orbit (km/s) as an impact.
static Vector3 ContactFrameVelocity(const Planetoid& other, const Planetoid& nearest, const DVec3& at, Vector3 nearestFrameVel) {
    return Vector3Add((Frames_SurfaceVelocity(other, at) - Frames_SurfaceVelocity(nearest, at)).toFloat(), nearestFrameVel);
}
static Vector3 FeltGravity(Scene& world, const DVec3& feet, Planetoid* nearest) {
    const DVec3 g = GravityAtD(world, feet);
    if (!FeltIsRelative(feet, nearest)) return g.toFloat();
    return (g - GravityAtD(world, nearest->worldPos)).toFloat();
}
void Player::update(Scene& world, const PlayerInput& in, float dt) {
    PROFILE_FUNCTION();
    // #768: the collision cache lives for one update only. Carried across updates it made the physics depend on history:
    // a client replaying from a corrected state skipped differently from the server, and predictions drifted.
    for (bool& v : collideCache.valid) v = false;
    // (F2b, #1146) Moved from outside since the last update (a spawn, placeAt, a respawn, a load): the velocity was set in
    // world terms then, so the frame path converts it again instead of restoring the last flight's stored frame velocity.
    // A spawn or placeAt anchors the player (anchorTo sets lastEndPosition), so the spawn count tells those apart. In the game
    // the loading hold then moves the player anyway; with no hold (a server's remote players) the Fixed path would rebuild the
    // last flight's velocity after a placement mid-flight (F8b-2, #2113).
    movedOutside = (position - lastEndPosition).length() > 1e-6 || spawns != updateSpawns;
    if (movedOutside) { fxActive = false; fcActive = false; fcBody = 0; airVelBody = 0; viewBody = 0; }
    // (#1562) In creative, a double-tap of Space toggles flying as F does: the first tap still jumps (or, in flight, rises); a
    // second within kDoubleTapFly toggles; a third starts a new pair.
    bool doubleTap = false;
    sinceJumpTap += dt;
    if (in.hit(Act::Jump)) {
        doubleTap = doubleTapFly && flyAllowed && sinceJumpTap < kDoubleTapFly;
        sinceJumpTap = doubleTap ? 1e9f : 0.0f;
    }
    if (((in.hit(Act::Fly) || doubleTap) && (flyAllowed || flyMode)) || (flyMode && !flyAllowed)) { // out of creative: F only lands (#531)
        flyMode = !flyMode;
        velocity = { 0, 0, 0 };
        if (!flyMode) {
            // Drop in moving with the ground below (a spinning planet's surface is fast).
            if (Planetoid* p = NearestPlanetoidD(world, position)) velocity = Frames_SurfaceVelocity(*p, position).toFloat();
        }
        groundBody = nullptr;
    }
    impactSpeed = landingSpeed = 0.0f; // per-frame; cleared before the fly-mode early return too
    jetpackFiring = false;
    if (flyMode) {
        // (#1555) The body under the player, flying too: this branch returns before the walking path sets it, so a player who
        // started flying (every scenario does) had none - the build tool, /locate and /tp found "no world" and /tp never landed.
        nearest = NearestPlanetoidD(world, position);
        // Carry the camera with its reference body's motion since recordFly() (after last frame's input).
        bool rotate = false;
        Planetoid* ref = flyFollow ? flyReference(world, rotate) : nullptr;
        // (F4, #1412) Carried in the frame recordFly() saw, even when the reference has changed since: crossing a body's
        // altitude line (2 radii: spinning frame <-> centred) or into another body's SOI skipped the carry that frame and
        // left the camera a frame of the body's motion behind (frames/fly-cross: 3.2 m at the moon's line). The next
        // recordFly() takes the new frame.
        if (ref && ref->id != flyBodyId && flyBodyId != 0) {
            ref = nullptr;
            for (auto& obj : world.objects)
                if (auto* q = dynamic_cast<Planetoid*>(obj.get())) if (q->id == flyBodyId) { ref = q; break; }
        }
        if (ref && ref->id == flyBodyId) {
            if (flyRotate) {
                flyCarried = ref->toWorldD(flyLocal);
                flyCarriedValid = true;
                world.camera.position = (flyCarried - world.renderOrigin).toFloat();
                world.camera.target = Vector3Add(world.camera.position, ref->dirToWorld(flyLocalFwd));
                world.camera.up = ref->dirToWorld(flyLocalUp);
            } else { // (F4, #1412) centred on the body: its centre + the double offset, the view kept
                const Vector3 view = Vector3Subtract(world.camera.target, world.camera.position);
                world.camera.position = (ref->worldPos + flyOffset - world.renderOrigin).toFloat();
                world.camera.target = Vector3Add(world.camera.position, view);
            }
        }
        // Near a body the fly camera keeps the horizon level (#232): its local up while gravity there is felt (> 0.5 m/s^2),
        // else the free 6-DOF camera (space, near-zero g).
        {
            bool rot2 = false;
            Planetoid* near = flyReference(world, rot2);
            gFlyGravityUp = { 0, 0, 0 };
            if (near) {
                const DVec3 rel = (world.renderOrigin + world.camera.position) - near->worldPos;
                const double d = rel.length(), R = (double)near->radius * near->scale;
                const double g = near->surfaceGravity * (R * R) / std::max(d * d, 1.0);
                if (g > 0.5 && d > 1.0) gFlyGravityUp = { (float)(rel.x / d), (float)(rel.y / d), (float)(rel.z / d) };
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
            // changes made elsewhere between frames are kept). (On the frame path F2c turns it below, grounded or not.)
            if (!fxActive) {
                const Quaternion now = groundBody->rotationQuat();
                const Quaternion delta = QuaternionMultiply(now, QuaternionInvert(groundRotPrev));
                forward = Vector3Normalize(Vector3RotateByQuaternion(forward, delta));
                up = Vector3Normalize(Vector3RotateByQuaternion(up, delta));
            }
        }
    }
    // Airborne: carry along with the body we were near (its spin and orbit since last frame); the
    // motion of our own (world-frame) velocity is then integrated relative to its surface below.
    if (airBody && Planetoid::find(airId) != airBody) airBody = nullptr;
    if (!grounded && airBody) {
        if ((position - lastEndPosition).length() > 1e-6) airBody = nullptr; // moved from outside (respawn, load, ship)
        else position = airBody->toWorldD(airLocal);
    }
    spaceOrigin = position; // (F8b-4b) the float space of this update's physics: offsets from here
    Vector3 feet = toSpace(position);
    // Eye smoothing (below) compares the feet in the ground's frame before and after this update.
    Planetoid* const eyeBody = (grounded && groundBody) ? groundBody : (airBody ? airBody : nullptr);
    const DVec3 eyeFrom = (grounded && groundBody) ? groundLocal : airLocal;

    nearest = NearestPlanetoidD(world, position);
    // Plan 11 F2b (#667, #1146), RAYLIBGAME_FRAMES=1: in a Fixed frame the velocity is always the speed against the turning
    // ground - airborne too, not only on the ground - so nothing converts at landing or lift-off and the move is velocity *
    // dt. An airborne frame velocity is kept in the body's axes between updates (it turns with the body, as the position
    // does); entering or leaving the frame path converts once. (Centred frames keep today's path until they are wired.)
    {
        const bool fx = nearest && frame.frame.kind == FrameKind::Fixed && frame.frame.body == nearest->netId;
        // (F8b-2) Above the altitude line the frame is Centred: the velocity is against the body's centre (world axes).
        const bool fc = !grounded && nearest && frame.frame.kind == FrameKind::Centred && frame.frame.body == nearest->netId;
        if (fcActive && !(fc && fcBody == nearest->id)) { // leaving it: back to world-frame first (the Fixed entry below expects that)
            if (const Planetoid* was = Planetoid::find(fcBody)) velocity = Vector3Add(velocity, Frames_Pose(*was).vel.toFloat());
        }
        if (fx && !grounded) {
            if (fxActive && airVelBody == nearest->id) velocity = Vector3Scale(nearest->dirToWorld(airVelDir), airVelLen); // turned with the body
            else velocity = Vector3Subtract(velocity, surfaceVel(*nearest, feet));                                          // into the frame, once
        } else if (!fx && fxActive && !grounded && nearest) {
            velocity = Vector3Add(velocity, surfaceVel(*nearest, feet));                                                    // back to world-frame
        }
        // F2c (#1146, #1143): the view turns with the frame's body, grounded or airborne, by its rotation since the last
        // update - R_now R_prev^T in double from the same pose as every transform (the float quaternion delta above was
        // grounded only). A delta, so mouse look between updates is kept.
        if (fxActive && nearest && viewBody == nearest->id) {
            double R[3][3]; DVec3 w;
            BodyRotationAndSpin(*nearest, R, w);
            double M[3][3]; // R_now * R_prev^T
            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j) M[i][j] = R[i][0] * viewR[j][0] + R[i][1] * viewR[j][1] + R[i][2] * viewR[j][2];
            auto turn = [&M](Vector3 v) {
                return Vector3Normalize({ (float)(M[0][0] * v.x + M[0][1] * v.y + M[0][2] * v.z), (float)(M[1][0] * v.x + M[1][1] * v.y + M[1][2] * v.z),
                                          (float)(M[2][0] * v.x + M[2][1] * v.y + M[2][2] * v.z) });
            };
            forward = turn(forward);
            up = turn(up);
        }
        fxActive = fx;
        if (fc) {
            const BodyPose p = Frames_Pose(*nearest);
            // Carried with the centre only, not turned with the spin (the airBody carry above did both): the frame state is
            // the offset from the centre at the end of the last update (updateFrame).
            if (!movedOutside) { position = p.centre + frame.pos; feet = toSpace(position); }
            if (!(fcActive && fcBody == nearest->id)) velocity = Vector3Subtract(velocity, p.vel.toFloat()); // into the frame, once
        }
        fcActive = fc;
        fcBody = fc ? nearest->id : 0;
    }
    // In water when the body's middle is (static or flowing water).
    const DVec3 middle = fromSpace(Vector3Add(feet, Vector3Scale(up, 0.9f)));
    inWater = !flyMode && IsUnderwaterD(world, middle);
    riverCurrent = inWater ? RiverCurrentD(nearest, middle) : Vector3{ 0, 0, 0 }; // a river pushes (#361)
    // "Up" is against the summed gravity. Where gravity is too weak to define it (deep space, the
    // balance point between bodies) it is left alone so the view doesn't spin; in the air it turns
    // at most 2 rad/s. On the ground it follows exactly, so walking around a body stays continuous.
    {
        // Felt gravity: near a body on an orbit, it and everything on it fall freely around its parent,
        // so subtract the body's own acceleration (the outside gravity at its centre). Otherwise the
        // parent's pull would tilt "up" (~4 deg on the moon) and jumps would go out sideways.
        const Vector3 g = FeltGravity(world, fromSpace(feet), nearest);
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
    if (in.enabled) {
        // (mouse look: look(), once per frame)
        const Vector3 flat = Vector3Normalize(Vector3Subtract(forward, Vector3Scale(up, Vector3DotProduct(forward, up))));
        const Vector3 side = Vector3CrossProduct(flat, up);
        // On the ground and in the air under gravity WASD is level (walking / air control); in deep
        // space it fires the suit thrusters along the view.
        const Vector3 f = (grounded || gravity >= 0.5f) ? flat : forward;
        if (in.down(Act::Forward)) wish = Vector3Add(wish, f);
        if (in.down(Act::Back) && !onLadder) wish = Vector3Subtract(wish, f); // (#1105: on a ladder S only climbs down; it walked you off it)
        if (in.down(Act::Right)) wish = Vector3Add(wish, side);
        if (in.down(Act::Left)) wish = Vector3Subtract(wish, side);
        if (Vector3Length(wish) > 0.0f) wish = Vector3Normalize(wish);
        if (in.down(Act::Sprint)) wish = Vector3Scale(wish, runMultiplier);
        jump = in.hit(Act::Jump);
        jetUp = in.down(Act::Jump);
        jetDown = in.down(Act::Down);
        if (in.hit(Act::Dampeners)) dampeners = !dampeners;
        if (in.hit(Act::Hover) && canHover) hovering = !hovering; // flight pack (#163)
        if (in.hit(Act::Prone)) proneWanted = !proneWanted; // Z (#578)
    }
    // Posture (#578): prone while wanted, crouched while Ctrl is held on the ground (Ctrl still swims down and descends
    // in the air), else upright; Space from prone gets up. Lying down or crouching always works; getting taller needs
    // the room above, else the lower posture stays.
    {
        const bool onFoot = grounded || coyote > 0.0f; // (grounded flickers over bumps: the coyote time bridges it)
        if (inWater || hovering || flyMode) proneWanted = false;
        if (proneWanted && jump) { proneWanted = false; jump = false; }
        const bool crouchWanted = in.enabled && in.down(Act::Down) && onFoot && !inWater && !hovering;
        Posture want = proneWanted ? Posture::Prone : crouchWanted ? Posture::Crouch : Posture::Stand;
        if (want < posture && !roomAbove(world, postureHeight(want))) { // taller: room?
            want = want == Posture::Stand && posture == Posture::Prone && roomAbove(world, kCrouchHeight) ? Posture::Crouch : posture;
            if (posture == Posture::Prone && want != Posture::Prone) proneWanted = false;
        }
        posture = want;
        height = postureHeight(posture);
        eyeHeight += (postureEye(posture) - eyeHeight) * std::min(1.0f, dt * 12.0f); // eased (~0.1 s)
        if (posture != Posture::Stand) {
            if (Vector3LengthSqr(wish) > 0.0f) wish = Vector3Scale(Vector3Normalize(wish), posture == Posture::Prone ? kProneSpeed : kCrouchSpeed); // no sprint
            if (posture == Posture::Prone) { jump = false; jetUp = false; }
        }
    }
    if (!canHover || !hoverPowered) hovering = false;
    // Fall safety (#163): a charged flight pack catches a fall faster than 12 m/s against the ground below.
    if (canHover && hoverPowered && !hovering && !grounded && !inWater && gravity >= 0.5f && nearest) {
        const float vUp = Vector3DotProduct(Vector3Subtract(velocity, groundMotion(world, toSpace(position))), up);
        if (vUp < -12.0f) hovering = true;
    }
    // (#1105) A ladder's input for the physics steps: W climbs, S descends; Space pushes off (not also a ground jump)
    climbInput = in.enabled && !flyMode ? (in.down(Act::Forward) ? 1.0f : in.down(Act::Back) ? -1.0f : 0.0f) : 0.0f;
    climbJump = onLadder && jump;
    if (climbJump) jump = false;
    // Jump buffer: a press is remembered briefly, so it is not lost on a frame where "grounded"
    // flickers off over bumpy ground (see coyote time in stepPhysics).
    if (jump) jumpBuffer = 0.15f;

    // Fixed substeps keep collision stable at any frame rate, and are shortened further at speed
    // (<= 0.2 m of travel each, under the body radius) so a fast player can't pass through thin walls.
    const float h = 1.0f / 120.0f;
    // Speed relative to the ground (airborne velocity is world-frame and carries the spin).
    const Vector3 relVel = !grounded ? Vector3Subtract(velocity, groundMotion(world, toSpace(position))) : velocity;
    const int steps = std::clamp(std::max((int)ceilf(dt / h), (int)ceilf(Vector3Length(relVel) * dt / 0.2f)), 1, 1024);
    const float sdt = dt / steps;
    // Airborne by an orbiting body: move in its free-falling frame - take its change of velocity
    // since last frame, then fall with the felt gravity. A rail orbit ignores other bodies' pull (a
    // big rock next to a small one), so world gravity drifted the player off a 5 m rock as fast as
    // it pulled them in: they hovered a metre up and never landed.
    airFelt = false;
    if (!fxActive && !fcActive && !grounded && airBody && airBody == nearest && FeltIsRelative(position, nearest)) {
        velocity = Vector3Add(velocity, Vector3Subtract(nearest->linearVelocity, airFrameVel));
        airFelt = true;
    }
    // In water Space and Ctrl swim instead of firing the jetpack.
    swimUp = inWater && jetUp; swimDown = inWater && jetDown;
    if (inWater) jetUp = jetDown = false;
    for (int i = 0; i < steps; ++i) { stepPhysics(world, sdt, wish, jetUp, jetDown); }

    if (grounded && nearest) {
        groundBody = nearest;
        groundId = nearest->id;
        groundLocal = nearest->toLocalD(position);
        groundRotPrev = nearest->rotationQuat();
    } else {
        groundBody = nullptr;
    }
    // Eye smoothing: walking, the collision push-out moves the feet up and down by millimetres to
    // centimetres from one frame to the next (marching-cubes facets, the contact normal changing),
    // and the camera showed every step as shake - worst climbing steep slopes. The eye follows a
    // low-passed vertical *speed* of the feet in the ground's frame, so a steady climb has no lag and
    // only the frame-to-frame noise is filtered; a slow pull keeps it within 25 cm of the body.
    {
        float dz = 0.0f;
        const float offBefore = eyeOffset;
        const bool known = eyeBody && nearest == eyeBody && dt > 0.0f;
        if (known) {
            const DVec3 d = eyeBody->toLocalD(position) - eyeFrom;
            const Vector3 u = eyeBody->dirToLocal(up);
            dz = (float)(d.x * u.x + d.y * u.y + d.z * u.z) * eyeBody->scale;
        }
        if (eyeSmoothing && known && grounded && fabsf(dz) < 1.0f) {
            const float vz = dz / dt;
            const float sc = std::max(0.1f, eyeSmoothScale); // #168: the setting scales the filter and the leash
            eyeVz += (vz - eyeVz) * (1.0f - expf(-dt / (0.10f * sc)));
            eyeOffset += eyeVz * dt - dz;
            eyeOffset *= expf(-dt * 5.0f / sc);
            eyeOffset = std::clamp(eyeOffset, -0.25f * sc, 0.25f * sc);
        } else {
            // In the air (jumps, falls, a frame of lost contact) the eye follows the body and eases
            // out what is left; the speed is kept so landing again on a slope doesn't lag.
            eyeVz = known ? dz / dt : 0.0f;
            eyeOffset *= expf(-dt * 15.0f);
        }
        eyeStep = known ? dz + (eyeOffset - offBefore) : 0.0f; // tests: the eye's height change this update, ground frame
        eyeStepKnown = known;
        // Sideways (#255), in the air too (wedged in a crevice the player is often not grounded): the same low-passed
        // velocity, a short time constant and a 3 cm leash so walking off or stopping never trails visibly.
        Vector3 dh = { 0, 0, 0 };
        if (known) {
            const DVec3 d = eyeBody->toLocalD(position) - eyeFrom;
            const Vector3 u = eyeBody->dirToLocal(up);
            dh = Vector3Subtract(Vector3Scale(Vector3{ (float)d.x, (float)d.y, (float)d.z }, eyeBody->scale), Vector3Scale(u, dz));
        }
        if (eyeSmoothing && known && eyeOffsetBody == eyeBody->id && Vector3Length(dh) < 0.5f) {
            eyeVh = Vector3Lerp(eyeVh, Vector3Scale(dh, 1.0f / dt), 1.0f - expf(-dt / 0.05f));
            eyeOffsetH = Vector3Scale(Vector3Add(eyeOffsetH, Vector3Subtract(Vector3Scale(eyeVh, dt), dh)), expf(-dt * 10.0f));
            const float len = Vector3Length(eyeOffsetH);
            if (len > 0.03f) eyeOffsetH = Vector3Scale(eyeOffsetH, 0.03f / len);
        } else {
            eyeVh = known ? Vector3Scale(dh, 1.0f / dt) : Vector3{ 0, 0, 0 };
            eyeOffsetH = { 0, 0, 0 };
        }
        eyeOffsetBody = known ? eyeBody->id : 0;
    }
    if (!grounded && nearest) {
        airBody = nearest;
        airId = nearest->id;
        airLocal = nearest->toLocalD(position);
        airFrameVel = nearest->linearVelocity;
        if (fxActive) { // (F2b) the frame velocity in the body's axes, turned back with it next update
            airVelLen = Vector3Length(velocity);
            airVelDir = airVelLen > 1e-9f ? nearest->dirToLocal(Vector3Scale(velocity, 1.0f / airVelLen)) : Vector3{ 0, 0, 0 };
            airVelBody = nearest->id;
        }
    } else {
        airBody = nullptr;
        airVelBody = 0; // (F2b) on the ground: the next flight starts from the grounded velocity, not the last one stored
    }
    lastEndPosition = position;
    updateSpawns = spawns;
    updateFrame(world);
    if (fxActive && nearest) { DVec3 w; BodyRotationAndSpin(*nearest, viewR, w); viewBody = nearest->id; } // (F2c) the pose the view was last turned to
    else viewBody = 0;
    CacheNetState(world);
    // Render interpolation (#315): the feet after this step and the one before, in the nearest body's frame. A new
    // body, a teleport (spawns) or a jump of over 5 m (placed from outside) starts again from here.
    if (nearest) {
        const DVec3 l = nearest->toLocalD(position);
        // (F2d) the update knows when the player was moved from outside, so the interpolation restarts on that, not on a
        // guess from a 5 m step (which also cut it for anything faster than 300 m/s).
        const bool cont = nearest->id == interpBody && interpSpawns == spawns && !movedOutside;
        interpPrev = cont ? interpCur : l;
        interpCur = l;
        interpBody = nearest->id;
        interpSpawns = spawns;
    } else interpBody = 0;
}

void Player::look(Vector2 m) {
    forward = RotateAround(forward, up, -m.x * mouseSensitivity * lookScale);
    const Vector3 right = Vector3Normalize(Vector3CrossProduct(forward, up));
    const Vector3 pitched = RotateAround(forward, right, -m.y * mouseSensitivity * lookScale);
    if (fabsf(Vector3DotProduct(pitched, up)) < 0.985f) forward = pitched; // clamp ~10 deg from straight up/down
}

Vector3 Player::groundVelocity(Scene& world) const {
    if (fxActive || grounded) return velocity;
    if (fcActive) return Vector3Subtract(velocity, groundMotion(world, toSpace(position)));
    Planetoid* body = nearest ? nearest : NearestPlanetoidD(world, position);
    return body ? Vector3Subtract(velocity, Frames_SurfaceVelocity(*body, position).toFloat()) : velocity;
}

DVec3 Player::renderPosition(float alpha) const {
    if (flyMode || !interpBody || interpSpawns != spawns) return position;
    const Planetoid* b = Planetoid::find(interpBody);
    if (!b) return position;
    // Only while nothing moved the player since the last step (benches, respawns, the death hold): then the body-local
    // feet still match `position`.
    if ((b->toWorldD(interpCur) - position).length() > 1e-3) return position;
    const double a = std::clamp((double)alpha, 0.0, 1.0);
    return b->toWorldD(interpPrev + (interpCur - interpPrev) * a);
}

void Player::stepPhysics(Scene& world, float dt, Vector3 wish, bool jetUp, bool jetDown) {
    Vector3 feet = toSpace(position);

    const bool wasGrounded = grounded;

    // Gravity summed over every body (inverse-square, no cutoff, so orbits work).
    const Vector3 gvec = GravityAtD(world, position).toFloat();   // world frame: airborne velocity is world-frame
    const Vector3 gFelt = FeltGravity(world, position, nearest); // the ground's frame: grounded velocity is relative to it
    gravity = Vector3Length(gFelt);

    // (#1105) A ladder: the capsule's lower or upper sphere over a ladder's cell. Climbing up from the ground starts as a
    // jump does (the ground branch takes it off, with its frame); in the air the ladder holds the climb speed (below).
    ladderCooldown = std::max(0.0f, ladderCooldown - dt);
    onLadder = false;
    if (!flyMode && ladderCooldown <= 0.0f && nearest && nearest->plots) {
        const float r = 0.35f / nearest->scale;
        const Vector3 lo = localOf(*nearest, Vector3Add(feet, Vector3Scale(up, 0.45f))), hi = localOf(*nearest, Vector3Add(feet, Vector3Scale(up, 1.25f)));
        for (const auto& pl : nearest->plots->plots) if (pl->touchesPart(lo, r, PartType::Ladder) || pl->touchesPart(hi, r, PartType::Ladder)) { onLadder = true; break; }
    }
    if (onLadder && grounded && climbInput > 0.0f) jumpBuffer = 0.15f;

    const float vUp = Vector3DotProduct(velocity, up);
    bool resting = false; // grounded with no walking input
    if (grounded) {
        // Direct ground control along the contact plane, pressed lightly into it so crests and dips
        // keep contact. Collision-induced upward speed is discarded (only jumps leave the ground).
        const Vector3 n = groundNormal;
        Vector3 along = Vector3Subtract(wish, Vector3Scale(n, Vector3DotProduct(wish, n)));
        const float wl = Vector3Length(wish), al = Vector3Length(along);
        if (al > 1e-4f) along = Vector3Scale(along, wl / al);
        Vector3 walk = Vector3Scale(along, walkSpeed * statusSpeed * (inWater ? 0.55f : 1.0f)); // wading is slow; (#1215) chill
        walk = Vector3Add(walk, Vector3Scale(Vector3Subtract(riverCurrent, Vector3Scale(n, Vector3DotProduct(riverCurrent, n))), 0.4f)); // wading in a river: pushed downstream (#361)
        if (dodgeLeft > 0.0f) { // (#1218) a dodge: a fixed step along the ground, whatever the keys
            const Vector3 dd = Vector3Subtract(dodgeDir, Vector3Scale(n, Vector3DotProduct(dodgeDir, n)));
            const float dl = Vector3Length(dd);
            if (dl > 1e-4f) walk = Vector3Scale(dd, dodgeSpeed / dl);
        }
        if (smoothMotion && dodgeLeft <= 0.0f) {
            // Reduce motion: ease into and out of walking (instant starts and stops are hard to watch).
            const Vector3 cur = Vector3Subtract(velocity, Vector3Scale(n, Vector3DotProduct(velocity, n)));
            const Vector3 dv = Vector3Subtract(walk, cur);
            const float dl = Vector3Length(dv), maxDv = 20.0f * dt;
            walk = dl > maxDv ? Vector3Add(cur, Vector3Scale(dv, maxDv / dl)) : walk;
        }
        // At rest, no press into the ground: the push back out is found to a few mm, so pressing in
        // every substep made a standing camera (and the shadow map with it) vibrate by 5-9 mm a frame.
        // The ground probe keeps a resting player grounded.
        resting = Vector3LengthSqr(walk) < 1e-6f;
        velocity = resting ? Vector3{ 0, 0, 0 } : Vector3Subtract(walk, Vector3Scale(n, 1.5f));
        (void)vUp;
        // A jump starts from the walking velocity only: keeping the press into a sloped ground would
        // push you downhill on takeoff.
        if (Vector3LengthSqr(pendingPush) > 1e-6f) {
            // Knocked back (#165): off the ground like a jump, at the push (surface-relative, as grounded velocity is).
            velocity = pendingPush;
            grounded = false;
            jumpBuffer = coyote = 0.0f;
            jumpLock = 0.15f;
        } else if (jumpBuffer > 0.0f) {
            velocity = Vector3Add(Vector3Scale(along, walkSpeed * statusSpeed), Vector3Scale(up, jumpSpeed * (statusSpeed > 0.0f ? 1.0f : 0.0f)));
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
        velocity = Vector3Add(velocity, pendingPush); // knocked back in the air (#165): the push adds to the flight
        if (jumpBuffer > 0.0f && coyote > 0.0f) {
            // Coyote time: just left the ground (or "grounded" flickered off for a frame) - still jump,
            // setting the upward speed relative to the surface below.
            const Vector3 surface = groundMotion(world, feet);
            velocity = Vector3Add(velocity, Vector3Scale(up, jumpSpeed - Vector3DotProduct(Vector3Subtract(velocity, surface), up)));
            jumpBuffer = coyote = 0.0f;
            jumpLock = 0.15f;
        }
        Vector3 thrust = { 0, 0, 0 };
        if (hovering && inGravity && !inWater) {
            // Hover (#163): the velocity eases fast to the wished one against the ground below (no drift when the keys
            // are let go); Space/Ctrl climb and descend. Gravity is cancelled below; no jetpack fuel.
            const Vector3 surface = groundMotion(world, feet);
            Vector3 target = Vector3Add(surface, Vector3Scale(wish, hoverSpeed));
            if (jetUp) target = Vector3Add(target, Vector3Scale(up, hoverClimb));
            if (jetDown) target = Vector3Subtract(target, Vector3Scale(up, hoverClimb));
            velocity = Vector3Add(target, Vector3Scale(Vector3Subtract(velocity, target), expf(-10.0f * dt)));
        } else if (inGravity) {
            // Air control (free, no fuel): steer toward the wished walking velocity relative to the
            // ground below - enough to reach a ledge mid-jump. With no input, momentum is untouched.
            if (Vector3Length(wish) > 0.0f) {
                const Vector3 surface = groundMotion(world, feet);
                const Vector3 rel = Vector3Subtract(velocity, surface);
                const Vector3 relLevel = Vector3Subtract(rel, Vector3Scale(up, Vector3DotProduct(rel, up)));
                const Vector3 dv = Vector3Subtract(Vector3Scale(wish, walkSpeed * statusSpeed), relLevel);
                const float dvl = Vector3Length(dv), maxDv = airControlAccel * dt;
                if (dvl > 1e-5f) velocity = Vector3Add(velocity, Vector3Scale(dv, std::min(1.0f, maxDv / dvl)));
            }
        } else if (Vector3Length(wish) > 0.0f) {
            thrust = Vector3Normalize(wish); // deep space: suit thrusters along the view
        }
        // Jetpack (Space up / Ctrl down, and WASD in deep space) burns fuel; recharging it on the ground
        // is what costs suit power (see the survival code).
        if (jetUp && !hovering) thrust = Vector3Add(thrust, up);
        if (jetDown && !hovering) thrust = Vector3Subtract(thrust, up);
        const float tl = Vector3Length(thrust);
        if (tl > 0.0f && jetpackFuel > 0.0f) {
            velocity = Vector3Add(velocity, Vector3Scale(thrust, jetpackAccel * dt / std::max(tl, 1.0f)));
            jetpackFuel = std::max(0.0f, jetpackFuel - dt * jetpackFuelPerSec);
            jetpackFiring = true;
        } else if (tl == 0.0f && dampeners && !inGravity) {
            // Inertial dampeners (B toggles) only in deep space, so they never fight a fall or an orbit.
            // Relative to the nearest body: damping the world-frame velocity (a ring rock's ~3 km/s
            // orbit) dragged the player off small rocks at hundreds of m/s^2.
            const Vector3 frame = groundMotion(world, feet);
            velocity = Vector3Add(frame, Vector3Scale(Vector3Subtract(velocity, frame), 1.0f - std::min(1.0f, 0.6f * dt)));
        }
    }
    pendingPush = { 0, 0, 0 };
    if (dodgeLeft > 0.0f) dodgeLeft = std::max(0.0f, dodgeLeft - dt); // (#1218)
    // On the ground (surface-relative velocity) the felt gravity; in the air the world's (a rock's
    // orbit is then followed by the frame carry): the world's planet pull on the ground of a ring rock
    // shoved the player off it.
    Vector3 gAccel = Vector3Scale(grounded ? gFelt : (fxActive || fcActive) ? frameGravity(world, feet) : airFelt ? gFelt : gvec, gravityScale * dt);
    if (hovering && !grounded && !inWater && gravity >= 0.5f) gAccel = { 0, 0, 0 }; // the flight pack holds the height (#163)
    if (inWater && !grounded) {
        // Swimming: buoyancy all but cancels gravity, the water drags (relative to the ground's frame),
        // Space swims up and Ctrl down (no jetpack under water).
        gAccel = Vector3Scale(gAccel, 0.12f);
        Vector3 frame = groundMotion(world, feet);
        frame = Vector3Add(frame, riverCurrent); // a river carries a swimmer downstream (#361)
        velocity = Vector3Add(frame, Vector3Scale(Vector3Subtract(velocity, frame), expf(-1.8f * dt)));
        if (swimUp) velocity = Vector3Add(velocity, Vector3Scale(up, 9.0f * dt));
        if (swimDown) velocity = Vector3Subtract(velocity, Vector3Scale(up, 6.0f * dt));
    }
    if (onLadder && !grounded && !wasGrounded) { // (not the step that leaves the ground: its velocity is still the ground's)
        // (#1105) On a ladder: no gravity; up and down at the climb speed against the surface (W, S; nothing holds); sideways
        // only as the walking input asks (W into the ladder, A/D along it), never a leftover drift (it carried the player off
        // the ladder going down); Space pushes off, back from where the player faces (the ladder), a little up.
        const Vector3 surface = groundMotion(world, feet);
        Vector3 side = Vector3Subtract(wish, Vector3Scale(up, Vector3DotProduct(wish, up)));
        if (const float sl = Vector3Length(side); sl > 1e-4f) side = Vector3Scale(side, kClimbSpeed / std::max(1.0f, sl));
        velocity = Vector3Add(surface, Vector3Add(side, Vector3Scale(up, climbInput * kClimbSpeed)));
        gAccel = { 0, 0, 0 };
        if (climbJump) {
            Vector3 back = Vector3Subtract(forward, Vector3Scale(up, Vector3DotProduct(forward, up)));
            back = Vector3LengthSqr(back) > 1e-6f ? Vector3Scale(Vector3Normalize(back), -kLadderPushOff) : Vector3{ 0, 0, 0 };
            velocity = Vector3Add(velocity, Vector3Add(back, Vector3Scale(up, 2.0f)));
            climbJump = false; ladderCooldown = 0.4f; onLadder = false;
        }
    }
    // Static friction: standing on walkable ground, only gravity's push into the ground applies.
    // (The downhill part made a standing player creep ~1.6 cm/s down slopes at fixed substeps.)
    if (grounded) gAccel = resting ? Vector3{ 0, 0, 0 } : Vector3Scale(groundNormal, std::min(0.0f, Vector3DotProduct(gAccel, groundNormal)));
    velocity = Vector3Add(velocity, gAccel);

    // Airborne, velocity is world-frame but the terrain already stands where it will be at the end of
    // this frame (and update() carried us with it): move relative to its surface.
    // (F8b-2) In a Centred frame the position is carried with the centre: the move is the velocity too.
    const Vector3 frameVel = (!wasGrounded && nearest && !fxActive && !fcActive) ? surfaceVel(*nearest, feet) : Vector3{ 0, 0, 0 };
    const Vector3 move = Vector3Scale(Vector3Subtract(velocity, frameVel), dt);
    const Vector3 feetBefore = feet;
    feet = Vector3Add(feet, move);
    collideFrameVel = wasGrounded ? Vector3{ 0, 0, 0 } : groundMotion(world, feet); // impacts against the terrain's motion
    const Vector3 vPre = velocity; // before the push-out and the ground probe touch it
    { const auto t0 = std::chrono::steady_clock::now(); grounded = resolveCollisions(world, feet); gPlayerCollideSeconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count(); } // (#768)
    collideTrunks(feet);
    // Just jumped: the probe below the feet would still find the slope (they have risen only a few cm)
    // and re-ground the player, cancelling the jump. Collisions still push out; only grounding waits.
    if (jumpLock > 0.0f) { grounded = false; jumpLock = std::max(0.0f, jumpLock - dt); }
    // Touchdown (#263): the speed into the floor, from the velocity before collisions, however the landing was
    // found (a push-out, or the ground probe's skin, which at fall speeds is most landings and recorded nothing).
    // Only walkable floors land you (groundNormal); walls and ceilings never count as landings (#255).
    if (!wasGrounded && grounded) landingSpeed = std::max(landingSpeed, -Vector3DotProduct(Vector3Subtract(vPre, collideFrameVel), groundNormal));
    // Wedged (#255): airborne, meant to fall but held up by steep walls (a V, an overhang corner), the fall speed
    // must not build up while the walls stop it (it reached lethal speeds in a second or two).
    if (!grounded) {
        const float meantDown = -Vector3DotProduct(move, up);
        const float wentDown = -Vector3DotProduct(Vector3Subtract(feet, feetBefore), up);
        if (meantDown > 0.005f && wentDown < 0.25f * meantDown) {
            const float vDown = Vector3DotProduct(Vector3Subtract(velocity, frameVel), up);
            if (vDown < -1.0f) velocity = Vector3Add(velocity, Vector3Scale(up, -1.0f - vDown));
        }
    }
    // Just landed: on the ground velocity is relative to the surface (as the walking code expects), so
    // take the ground's own motion out now rather than leave a world-frame velocity for a frame.
    // (F8b-2) A Centred flight that lands is ground-relative from here, as any grounded velocity.
    if (!wasGrounded && grounded) { velocity = Vector3Subtract(velocity, collideFrameVel); fcActive = false; fcBody = 0; }
    // Leaving a spinning body: on the ground velocity is relative to the surface, so add the
    // surface's own motion; you keep the body's momentum, as you would really.
    // Landed and lifted off within one frame (contact flicker in a tight corner) leaves groundBody
    // unset (it is recorded at the end of the frame): the landing took out the nearest body's motion,
    // so give that back. Skipping it lost the planet's spin (~87 m/s) and the walls threw the player.
    if (wasGrounded && !grounded && nearest) {
        // Walked off an edge (#240): the walking velocity followed a contact normal tilting out over
        // the lip and carried the 1.5 m/s press into the ground, so leaving it threw the player down
        // (more when sprinting). If there is no ground a walkable slope's drop ahead and below, keep
        // only the level part and let gravity start from rest. Contact flickers on bumpy ground and
        // walking down a slope (ground still there) keep their velocity.
        const float vu = Vector3DotProduct(velocity, up);
        if (vu < 0.0f && jumpLock <= 0.0f) {
            Planetoid& p = *nearest;
            const Vector3 level = Vector3Subtract(velocity, Vector3Scale(up, vu));
            const float ll = Vector3Length(level);
            const float ahead = std::min(ll * 0.1f, 1.5f); // where the feet will be in 0.1 s
            const Vector3 probe = ll > 1e-4f ? Vector3Add(feet, Vector3Scale(level, ahead / ll)) : feet;
            const float depth = 0.4f + 1.6f * ahead; // ~58 degrees: steeper than walkable
            bool ground = false;
            for (float k = 0.1f; k <= depth && !ground; k += 0.15f) {
                const Vector3 at = Vector3Subtract(probe, Vector3Scale(up, k));
                ground = p.densityAt(localOf(p, at)) >= ISO_LEVEL;
                if (!ground && p.plots) // built parts count as ground too (stairs, platform steps)
                    for (auto& pl : p.plots->plots) {
                        Vector3 c = localOf(p, at), nl;
                        if (pl->collideSphere(c, 0.08f / p.scale, nl)) { ground = true; break; }
                    }
            }
            if (!ground) velocity = level;
        }
    }
    if (wasGrounded && !grounded && !fxActive) { // (F2b: lifting off keeps the frame velocity as it is)
        Planetoid* from = (groundBody && Planetoid::find(groundId) == groundBody) ? groundBody : nearest;
        if (from) velocity = Vector3Add(velocity, surfaceVel(*from, feet));
    }
    position = fromSpace(feet);
    jumpBuffer = std::max(0.0f, jumpBuffer - dt);
    if (!grounded) coyote = std::max(0.0f, coyote - dt);
}

// Room to get taller (#578): the capsule's spheres from its top now up to toHeight, each a little smaller than the
// capsule's (a shallow touch is fine: the collision pass pushes out of it), clear of solid ground and built parts.
bool Player::roomAbove(Scene& world, float toHeight) const {
    if (!nearest || toHeight <= height) return true;
    const Planetoid& p = *nearest;
    const Vector3 feet = toSpace(position);
    const float r = radius * 0.85f;
    for (float o = height - radius; ; o = std::min(o + 0.2f, toHeight - radius)) {
        const Vector3 c = Vector3Add(feet, Vector3Scale(up, o));
        if (p.densityAt(localOf(p, Vector3Add(c, Vector3Scale(up, r)))) >= ISO_LEVEL) return false; // the sphere's top
        if (p.plots)
            for (auto& pl : p.plots->plots) { Vector3 lc = localOf(p, c), nl; if (pl->collideSphere(lc, r / p.scale, nl)) return false; }
        if (o >= toHeight - radius) return true;
    }
}

// Push the feet/head spheres out of solid density, then probe just below the feet for ground.
// One body's terrain and build parts against the capsule spheres at `offsets` (#255): deepest contact first, one push a pass,
// each push measured and the frame's total capped (pushBudget, shared across bodies). `frameVel` is the body's surface motion
// there (impacts are against it). `cache`: the nearest body's sphere cache (#768), or null. True if standing on it.
bool Player::collideBody(Planetoid& p, const float* offsets, int kSamples, Vector3& feet, Vector3 frameVel, float& pushBudget, CollideCache* cache) {
    bool onGround = false;
    auto solid = [&](Vector3 c) { return p.densityAt(localOf(p, c)) >= ISO_LEVEL; };
    constexpr float kClearMargin = 0.05f; // m
    for (int iter = 0; iter < 6 && pushBudget > 0.0f; ++iter) {
        bool moved = false;
        // Contacts deepest first; one push per pass (the next pass re-measures).
        struct Cand { float d, off; Vector3 n; } cand[8]; int nc = 0;
        for (int si = 0; si < kSamples; ++si) { const float o2 = offsets[si];
            const Vector3 c2 = Vector3Add(feet, Vector3Scale(up, o2));
            const Vector3 l2 = localOf(p, c2);
            if (cache && cache->valid[si] && cache->dist[si] - Vector3Distance(l2, cache->centre[si]) * p.scale > radius + kClearMargin) continue; // clear: can't touch yet
            float d2; Vector3 n2;
            const bool ok = TerrainSurfaceDistanceD(p, fromSpace(c2), d2, n2);
            if (cache) { cache->valid[si] = ok; cache->dist[si] = d2; cache->centre[si] = l2; }
            if (ok && d2 < radius) cand[nc++] = { d2, o2, n2 };
        }
        std::sort(cand, cand + nc, [](const Cand& a, const Cand& b) { return a.d < b.d; });
        for (int ci = 0; ci < nc && !moved; ++ci) {
            const float off = cand[ci].off;
            const Vector3 n = cand[ci].n; // (measured just above, at this same centre)
            const Vector3 c = Vector3Add(feet, Vector3Scale(up, off));
            float push = 0.0f;
            if (!solid(c)) {
                // Centre in the open: is there solid ground within the radius along -n, and how far?
                if (!solid(Vector3Subtract(c, Vector3Scale(n, radius)))) continue; // a false contact
                float lo = 0.0f, hi = radius;                                       // lo open, hi solid
                for (int k = 0; k < 10; ++k) { // to ~0.4 mm
                    const float m = 0.5f * (lo + hi);
                    if (solid(Vector3Subtract(c, Vector3Scale(n, m)))) hi = m; else lo = m;
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
            const float vn = Vector3DotProduct(Vector3Subtract(velocity, frameVel), n);
            if (vn < 0.0f) { impactSpeed = std::max(impactSpeed, -vn); velocity = Vector3Subtract(velocity, Vector3Scale(n, vn)); }
            lastContactDot = Vector3DotProduct(n, up);
            // Steep walls must not convert a push into upward motion (no wall climbing).
            if (lastContactDot <= 0.55f) {
                const float vUpAfter = Vector3DotProduct(velocity, up);
                if (vUpAfter > vUpBefore) velocity = Vector3Subtract(velocity, Vector3Scale(up, vUpAfter - vUpBefore));
            }
            if (lastContactDot > walkableCos) { onGround = true; groundNormal = n; } // walkable slope (< ~57 deg)
            moved = true;
        }
        if (!moved) break;
    }
    // Build plot parts: the same capsule spheres pushed out of part shapes (exact convex cells). Only the first
    // kSamples offsets are set (the array is sized for the tallest override).
    if (p.plots)
        for (int iter = 0; iter < 2; ++iter)
            for (int si = 0; si < kSamples; ++si)
                for (auto& pl : p.plots->plots) {
                    // From a double body-local centre (#1732): a float one 60 km out steps ~4 mm, and its round trip jittered
                    // the pushes by that much a frame (kicks of 0.1-0.4 m/s against a wall).
                    const DVec3 centre = p.toLocalD(fromSpace(Vector3Add(feet, Vector3Scale(up, offsets[si]))));
                    Vector3 moved, nl;
                    Plot::Contacts hits;
                    // The spheres above the feet meet a part's top edge as a wall (#1732); the feet sphere keeps its slope
                    // over a lip, which steps up onto a slab.
                    if (!pl->collideSphereD(centre, radius / p.scale, moved, nl, &hits, si > 0 ? walkableCos : 0.0f)) continue;
                    const float len = Vector3Length(moved) * p.scale;
                    if (len > 1e-6f) feet = Vector3Add(feet, Vector3Scale(p.dirToWorld(Vector3Normalize(moved)), len));
                    // Each surface touched on its own (#1732): the sum of a floor's and a wall's normals is a 45-degree
                    // "surface" that turned the walk into the wall, and the press into the floor, into an upward and
                    // backward kick, and made it ground (dot 0.71), tilting the walking plane up the wall.
                    const float vUpBefore = Vector3DotProduct(velocity, up);
                    bool steep = false;
                    for (int k = 0; k < hits.count; ++k) {
                        const Vector3 n = p.dirToWorld(hits.n[k]);
                        const float vn = Vector3DotProduct(Vector3Subtract(velocity, frameVel), n);
                        if (vn < 0.0f) { impactSpeed = std::max(impactSpeed, -vn); velocity = Vector3Subtract(velocity, Vector3Scale(n, vn)); }
                        const float nu = Vector3DotProduct(n, up);
                        if (nu > walkableCos) { if (!onGround || nu > Vector3DotProduct(groundNormal, up)) groundNormal = n; onGround = true; }
                        else steep = true;
                    }
                    // A wall doesn't turn a push into upward motion (no wall climbing; as the terrain pass above).
                    if (steep) { const float vUpAfter = Vector3DotProduct(velocity, up); if (vUpAfter > vUpBefore && vUpAfter > 0.0f) velocity = Vector3Subtract(velocity, Vector3Scale(up, vUpAfter - std::max(vUpBefore, 0.0f))); }
                }
    return onGround;
}

bool Player::resolveCollisions(Scene& world, Vector3& feet) {
    if (!nearest) return false;
    Planetoid& p = *nearest;
    bool onGround = false;
    // A capsule (#255): sphere centres from `radius` up to `height - radius`, at most 0.23 m apart (no waist: the two
    // spheres at 0.45 and 1.15 m left a 0.28 m neck that jagged rock poked into). The player is `height` (1.7 m) tall, so
    // 2 m doorways and tunnels pass (user, 2026-11-14: the first capsule ran up to the eye and made the player 2.05 m;
    // RAYLIBGAME_PLAYER_HEIGHT=2.05 is that one, for A/B); the eye (1.6 m) sits inside the top sphere, and the eye
    // keep-out below lowers the camera where a ceiling comes closer. Each pass pushes out of the deepest contact only:
    // two spheres pushing along different normals on jagged rock fought each other.
    static const float heightOverride = getenv("RAYLIBGAME_PLAYER_HEIGHT") ? (float)atof(getenv("RAYLIBGAME_PLAYER_HEIGHT")) : 0.0f;
    if (heightOverride > 2.0f * radius) { standHeight = heightOverride; if (posture == Posture::Stand) height = standHeight; } // (crouching and prone keep theirs, #578)
    constexpr int kMaxSamples = 8;
    const int kSamples = std::clamp((int)ceilf((height - 2.0f * radius) / 0.23f - 1e-3f) + 1, 2, kMaxSamples); // 5 at 1.7 m, 6 at 2.05 m
    float offsets[kMaxSamples];
    for (int i = 0; i < kSamples; ++i) offsets[i] = radius + (height - 2.0f * radius) * (float)i / (float)(kSamples - 1);
    // The linear estimate (iso - density) / |gradient| is only a distance where the density behaves
    // like one; dig and build edits (big density steps) break that, and it reported contacts in open
    // air that pushed up to 1 m per iteration, sometimes into the ground (the "wall launch"). So each
    // contact is confirmed by sampling along the normal, the push is the measured penetration, and the
    // total push in a frame is capped.
    auto solid = [&](Vector3 c) { return p.densityAt(localOf(p, c)) >= ISO_LEVEL; };
    float pushBudget = 0.6f;
    // #768: skip spheres that were measured clear of the terrain and can't have reached it since (see CollideCache).
    static const bool useCache = !(getenv("RAYLIBGAME_COLLIDE_CACHE") && atoi(getenv("RAYLIBGAME_COLLIDE_CACHE")) == 0);
    CollideCache& cc = collideCache;
    if (cc.body != &p || cc.edits != p.edits.get()) { for (bool& v : cc.valid) v = false; cc.body = &p; cc.edits = p.edits.get(); }
    onGround = collideBody(p, offsets, kSamples, feet, collideFrameVel, pushBudget, useCache ? &cc : nullptr) || onGround;
    // (Plan 11 F6, #1698) Every other body the capsule reaches, not only the nearest: two promoted ring rocks side by side, a
    // rock beside a station - the player went through the one that wasn't nearest. Each with its own surface motion for
    // the impact; the sphere cache stays the nearest body's (#768). A body is in reach when the capsule touches its bounding
    // sphere (radius x scale: the voxel field never reaches past it).
    const Vector3 mid = Vector3Add(feet, Vector3Scale(up, 0.5f * height));
    for (auto& obj : world.objects) {
        auto* other = dynamic_cast<Planetoid*>(obj.get());
        if (!other || other == &p || other->flat || (other->vehicle && !other->freeFall)) continue; // not a flying rocket (as frames.cpp)
        if ((fromSpace(mid) - other->worldPos).length() > other->radius * other->scale + 0.5f * height + radius) continue;
        onGround = collideBody(*other, offsets, kSamples, feet, ContactFrameVelocity(*other, p, fromSpace(mid), collideFrameVel), pushBudget, nullptr) || onGround;
    }
    if (p.plots) {
        // Standing on a part: a probe just below the feet sphere.
        if (!onGround && Vector3DotProduct(Vector3Subtract(velocity, collideFrameVel), up) <= 0.5f)
            for (auto& pl : p.plots->plots) {
                Vector3 c = localOf(p, Vector3Add(feet, Vector3Scale(up, radius - 0.12f))), nl;
                Plot::Contacts hits;
                if (!pl->collideSphere(c, radius / p.scale, nl, &hits)) continue;
                for (int k = 0; k < hits.count; ++k) { // the floor's own normal, not its sum with a wall beside the feet (#1732)
                    const Vector3 n = p.dirToWorld(hits.n[k]);
                    if (Vector3DotProduct(n, up) > walkableCos && (!onGround || Vector3DotProduct(n, up) > Vector3DotProduct(groundNormal, up))) { onGround = true; groundNormal = n; }
                }
                if (onGround) break;
            }
    }
    // Ground probe: the push-out is approximate and can leave a small gap; treat anything within a
    // few centimetres below the feet sphere as standing (unless moving away, e.g. jumping).
    if (!onGround && Vector3DotProduct(Vector3Subtract(velocity, collideFrameVel), up) <= 0.5f) {
        float dist; Vector3 n;
        const float skin = 0.12f;
        if (TerrainSurfaceDistanceD(p, fromSpace(Vector3Add(feet, Vector3Scale(up, radius))), dist, n) && dist < radius + skin && Vector3DotProduct(n, up) > walkableCos) { // (#1855 M5: grip boots lower walkableCos; 0.55 here kept them off steep ground)
            // Settle onto it, past a 2 cm dead band: the distance here is a linear estimate and the
            // push-out a bisection; settling every frame see-sawed between the two by 5-9 mm.
            const float gap = dist - radius;
            if (gap > 0.02f) feet = Vector3Subtract(feet, Vector3Scale(n, gap - 0.01f));
            onGround = true;
            groundNormal = n;
        } else if (solid(Vector3Subtract(feet, Vector3Scale(up, skin)))) {
            // The estimate overshoots on small, sharply curved bodies: solid right under the feet is
            // standing too (a resting player no longer presses into the ground to make contact).
            onGround = true;
            groundNormal = up;
        }
    }
    return onGround;
}

void Player::aimRay(const Scene& world, const Planetoid& body, DVec3& originLocal, Vector3& dirLocal, Vector3& dirWorld) const {
    if (!flyMode && lastEyeD.length() > 0.0) {
        dirWorld = Vector3Normalize(forward);
        originLocal = body.toLocalD(lastEyeD);
    } else {
        dirWorld = Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position));
        originLocal = body.toLocalRenderD(world.camera.position, world.renderOrigin);
    }
    dirLocal = body.dirToLocal(dirWorld);
}

void Player::applyToCamera(Scene& world) const {
    if (flyMode) return;
    const Vector3 feet = (renderPosition(SimAlpha()) - world.renderOrigin).toFloat(); // between the last two steps (#315)
    // The eye's offset from the feet, kept on its own (small numbers): (F4, #1412) the double eye below adds it to the
    // double feet; eye - feet in render space carried both points' rounding (~61 um 1 km from the render origin).
    Vector3 eyeOff = Vector3Scale(up, eyeHeight + eyeOffset);
    if (nearest && eyeOffsetBody == nearest->id) { // sideways smoothing (#255); dirToWorld normalises, so scale back
        const float len = Vector3Length(eyeOffsetH);
        if (len > 1e-6f) eyeOff = Vector3Add(eyeOff, Vector3Scale(nearest->dirToWorld(eyeOffsetH), len));
    }
    Vector3 eye = Vector3Add(feet, eyeOff);
    lastEyePush = 0.0f;
    // The eye sits at the head sphere's top (and the smoothing can lift it 0.25 m more), outside the collision
    // shape: keep it 0.25 m clear of the terrain (the near plane's reach) so it never looks through a ceiling or an
    // overhang (#255). It was pushed out along the surface normal each frame: against jagged walls and overhangs
    // that normal and the push size swung between facets frame to frame and the camera shook (and inside rock the
    // normal could flip: metre jumps). Now the eye only moves down the body's axis toward the head sphere's centre,
    // which the collision keeps clear: the highest point of that segment with the clearance (a short search). It
    // drops there at once (never into the rock) and comes back up eased (~0.15 s): one direction, one number.
    if (nearest) {
        const float clearance = 0.25f;
        auto clear = [&](Vector3 q) {
            float dist; Vector3 n;
            if (!TerrainSurfaceDistance(*nearest, q, dist, n)) return true;
            return dist >= clearance && nearest->densityAt(nearest->toLocal(q)) < ISO_LEVEL;
        };
        const Vector3 headCentre = Vector3Add(feet, Vector3Scale(up, height - radius)); // the capsule's top sphere
        const float span = Vector3Distance(headCentre, eye); // ~0.45 m (+ the eye offset)
        float want = 0.0f;                                   // how far below the nominal eye
        if (!clear(eye)) {
            float lo = 0.0f, hi = span;                      // lo: blocked, hi: clear (or the centre)
            for (int i = 0; i < 7; ++i) {
                const float mid = 0.5f * (lo + hi);
                if (clear(Vector3Subtract(eye, Vector3Scale(up, mid)))) hi = mid; else lo = mid;
            }
            want = hi;
        }
        const float dt = GameFrameTime();
        if (dt <= 0.0f || dt > 0.25f) eyeDrop = want;
        else if (want > eyeDrop) {
            // Down eased too (~40 ms) while the eye where it is now keeps 12 cm; closer than that, at once.
            float dist; Vector3 n;
            const Vector3 now = Vector3Subtract(eye, Vector3Scale(up, eyeDrop));
            const bool roomy = !TerrainSurfaceDistance(*nearest, now, dist, n) || (dist >= 0.12f && nearest->densityAt(nearest->toLocal(now)) < ISO_LEVEL);
            eyeDrop = roomy ? eyeDrop + (want - eyeDrop) * (1.0f - expf(-dt / 0.04f)) : want;
        } else eyeDrop += (want - eyeDrop) * (1.0f - expf(-dt / 0.15f)); // back up eased
        eye = Vector3Subtract(eye, Vector3Scale(up, eyeDrop));
        eyeOff = Vector3Subtract(eyeOff, Vector3Scale(up, eyeDrop));
        lastEyePush = eyeDrop;
    }
    world.camera.position = eye;
    world.camera.target = Vector3Add(world.camera.position, forward);
    // (F4, #1412) the eye in double: the interpolated feet (double) plus the eye's small offset from them. The render-space
    // camera rounds it (~61 um steps 1 km from the render origin), and aims taken from it wobbled while the player stood still.
    lastEyeD = renderPosition(SimAlpha()) + DVec3(eyeOff);
    // Reduce motion: the camera's up follows the ground's gently (no lurch when it changes).
    if (smoothMotion && Vector3Length(cameraUp) > 0.5f) cameraUp = Vector3Normalize(Vector3Lerp(cameraUp, up, std::min(1.0f, 4.0f * GameFrameTime())));
    else cameraUp = up;
    world.camera.up = cameraUp;
}

// Tree trunks are vertical capsules (bushes are soft): push the body cylinder out sideways and
// drop the velocity into the trunk, relative to the (moving) ground the tree stands on.
void Player::collideTrunks(Vector3& feet) {
    if (!nearest || !nearest->trees) return;
    static thread_local std::vector<Trunk>* nearP = nullptr; if (!nearP) nearP = new std::vector<Trunk>(); auto& near = *nearP; // (#1133, #485: MinGW frees thread_local storage before running its destructor at thread exit: never destroyed)
    near.clear();
    const Vector3 local = localOf(*nearest, feet);
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
    // Unmoved by input since the carry: the double pose, not its float rounding (which drifted ~0.2 m in 3 s).
    const DVec3 cam = world.renderOrigin + world.camera.position;
    flyLocal = ref->toLocalD(flyCarriedValid && (cam - flyCarried).length() < 0.01 ? flyCarried : cam);
    // (F4, #1412) unmoved by input since the centred carry: keep the double offset (the float camera rounds it each frame)
    const DVec3 offsetNow = cam - ref->worldPos;
    if (!(!rotate && ref->id == flyBodyIdPrev && (offsetNow - flyOffset).length() < 0.01)) flyOffset = offsetNow;
    flyBodyIdPrev = ref->id;
    flyCarriedValid = false;
    flyLocalFwd = ref->dirToLocal(Vector3Normalize(Vector3Subtract(world.camera.target, world.camera.position)));
    flyLocalUp = ref->dirToLocal(world.camera.up);
}

void Player::anchorTo(Planetoid& body) {
    airBody = &body;
    airId = body.id;
    airLocal = body.toLocalD(position);
    airFrameVel = body.linearVelocity;
    lastEndPosition = position;
}

// The state as the last step left it, taken there (a server that steps a player only by its inputs reads it later, when
// the body has moved on: the world position and velocity are then stale, the body-local ones are not).
void Player::CacheNetState(const Scene& world) {
    netCached = false;
    Planetoid* b = (grounded && groundBody) ? groundBody : nearest ? nearest : airBody;
    if (!b) return;
    const Quaternion toLocal = QuaternionInvert(b->rotationQuat()); // (dirToLocal normalizes: a speed needs the rotation)
    const Vector3 rel = (grounded || fxActive) ? velocity : Vector3Subtract(velocity, (fcActive && b == nearest) ? groundMotion(world, toSpace(position)) : Frames_SurfaceVelocity(*b, position).toFloat()); // (#1855 M3) the frame path flies against the ground already
    netCache.body = b->netId;
    netCache.local = b->toLocalD(position);
    netCache.vel = Vector3RotateByQuaternion(rel, toLocal);
    netCacheSpawns = spawns;
    netCached = true;
}

bool Player::netState(const Scene& world, NetState& s) const {
    if (!netCached || netCacheSpawns != spawns) { // never stepped since a spawn: as it stands now
        Planetoid* b = (grounded && groundBody) ? groundBody : nearest ? nearest : airBody;
        if (!b) return false;
        const Vector3 rel = (grounded || fxActive) ? velocity : Vector3Subtract(velocity, (fcActive && b == nearest) ? groundMotion(world, toSpace(position)) : Frames_SurfaceVelocity(*b, position).toFloat()); // (M3)
        s.body = b->netId;
        // (#856) The anchor's body-local feet, not the world position: a player spawned and not stepped yet (a server's
        // joiner still loading its world) keeps its spawn-time world position while the planet turns under it (~105 m/s at
        // the surface), so that position read later lay kilometres off; the client started there and was snapped back.
        if (grounded && groundBody == b) s.local = groundLocal;
        else if (!grounded && airBody == b) s.local = airLocal;
        else s.local = b->toLocalD(position);
        s.vel = Vector3RotateByQuaternion(rel, QuaternionInvert(b->rotationQuat()));
    } else {
        s.body = netCache.body; s.local = netCache.local; s.vel = netCache.vel;
    }
    s.grounded = grounded; s.hovering = hovering;
    s.posture = (uint8_t)posture;
    s.jumpBuffer = jumpBuffer; s.coyote = coyote; s.jumpLock = jumpLock; s.jetpackFuel = jetpackFuel;
    return true;
}

bool Player::applyNetState(Scene& world, const NetState& s) {
    Planetoid* b = Planetoid::findNet(s.body);
    if (!b) return false;
    position = b->toWorldD(s.local);
    const Vector3 rel = Vector3RotateByQuaternion(s.vel, b->rotationQuat());
    velocity = s.grounded ? rel : Vector3Add(rel, Frames_SurfaceVelocity(*b, position).toFloat());
    // (#1855 M4) The state is the server's now: drop the frame path's stored flight (its velocity, kept turned with the body, and
    // the view carry), or the next step rebuilt the old velocity from it and the correction was lost. The frame path starts over
    // from this world-frame velocity, as after a teleport.
    fxActive = false; airVelBody = 0; viewBody = 0;
    fcActive = false; fcBody = 0;
    frame = FramePos{}; // (F8b-2) not "moved from outside" (anchored): the next update must not restore the Centred offset from before
    grounded = s.grounded; hovering = s.hovering;
    posture = (Posture)std::min<uint8_t>(s.posture, (uint8_t)Posture::Prone);
    height = postureHeight(posture);
    jumpBuffer = s.jumpBuffer; coyote = s.coyote; jumpLock = s.jumpLock; jetpackFuel = s.jetpackFuel;
    nearest = b;
    if (grounded) {
        groundBody = b; groundId = b->id;
        groundLocal = s.local;
        groundRotPrev = b->rotationQuat();
        airBody = nullptr;
        lastEndPosition = position;
    } else {
        groundBody = nullptr;
        anchorTo(*b); // (sets lastEndPosition)
    }
    netCache = s; netCacheSpawns = spawns; netCached = true;
    return true;
}

// Plan 11 F2a (#667, #1146): the player's state in its frame, from today's state. The world-frame velocity: on the ground
// the velocity is relative to the surface (add the ground's motion), airborne it is world-frame already. The frame: the
// smallest SOI holding the feet, Fixed below the body's line (contact forces it), with hysteresis against the last one.
void Player::updateFrame(Scene& world) {
    // The velocity is against the ground: on the ground, and airborne on the frame path (F2b) - add the ground's motion.
    const Planetoid* ground = grounded ? ((groundBody && Planetoid::find(groundId) == groundBody) ? groundBody : nearest) : fxActive ? nearest : nullptr;
    Vector3 vw = ground ? Vector3Add(velocity, Frames_SurfaceVelocity(*ground, position).toFloat()) : velocity;
    if (!grounded && fcActive && nearest) vw = Vector3Add(velocity, Frames_Pose(*nearest).vel.toFloat()); // (F8b-2) against the centre
    const FrameRef to = Frames_Choose(Planetoid::All(), position, frame.frame, grounded);
    if (to.body != 0) Frames_FromGlobal(position, DVec3(vw.x, vw.y, vw.z), to, frame);
}

bool Player::restoreFrame(const FramePos& f) {
    Planetoid* b = Planetoid::findNet(f.frame.body);
    DVec3 pg, vg;
    if (!b || !Frames_ToGlobal(f, pg, vg)) return false;
    position = pg;
    frame = f;
    nearest = b;
    grounded = false; groundBody = nullptr;
    // Carried with the body from here as any airborne player (update's airBody path), and not "moved from outside": the
    // velocity below is already in the terms the next update expects.
    airBody = b; airId = b->id; airLocal = b->toLocalD(position);
    airFrameVel = b->linearVelocity;
    lastEndPosition = position;
    viewBody = 0;
    const double vl = f.vel.length();
    if (f.frame.kind == FrameKind::Fixed) { // the frame velocity (against the ground), kept in the body's axes
        const DVec3 d = vl > 1e-12 ? f.vel * (1.0 / vl) : DVec3{};
        airVelDir = { (float)d.x, (float)d.y, (float)d.z };
        airVelLen = (float)vl;
        airVelBody = b->id;
        velocity = Vector3Scale(b->dirToWorld(airVelDir), airVelLen);
        fxActive = true;
    } else { // (F8b-2) Centred: the velocity against the centre, world axes; the next update carries the offset in `frame`
        velocity = f.vel.toFloat();
        airVelBody = 0;
        fxActive = false;
        fcActive = true;
        fcBody = b->id;
    }
    updateSpawns = spawns;
    return true;
}

// Plan 11 F2b (#1146): the acceleration an airborne player feels in the nearest body's Fixed frame (felt gravity, Coriolis
// and centrifugal: Frames_Acceleration), in world axes. The velocity is the frame velocity (against the ground).
Vector3 Player::frameGravity(Scene& world, Vector3 feet) const {
    if (!nearest) return { 0, 0, 0 };
    if (fcActive) { // (F8b-2) Centred: the felt gravity (the centre's own fall taken out) in world axes, no spin forces
        FramePos fc;
        fc.frame = FrameRef{ nearest->netId, FrameKind::Centred };
        fc.pos = fromSpace(feet) - Frames_Pose(*nearest).centre;
        fc.vel = DVec3(velocity.x, velocity.y, velocity.z);
        return Frames_Acceleration(Planetoid::All(), fc).toFloat();
    }
    FramePos fp;
    fp.frame = FrameRef{ nearest->netId, FrameKind::Fixed };
    fp.pos = nearest->toLocalD(fromSpace(feet));
    const float vl = Vector3Length(velocity);
    if (vl > 1e-9f) { const Vector3 d = Vector3Scale(nearest->dirToLocal(Vector3Scale(velocity, 1.0f / vl)), vl); fp.vel = DVec3(d.x, d.y, d.z); }
    const DVec3 a = Frames_Acceleration(Planetoid::All(), fp);
    const double al = a.length();
    if (al < 1e-12) return { 0, 0, 0 };
    return Vector3Scale(nearest->dirToWorld(Vector3{ (float)(a.x / al), (float)(a.y / al), (float)(a.z / al) }), (float)al);
}

Vector3 Player::groundMotion(const Scene& world, Vector3 feet) const {
    if (!nearest || fxActive) return { 0, 0, 0 };
    if (fcActive) {
        const BodyPose p = Frames_Pose(*nearest);
        const DVec3 r = fromSpace(feet) - p.centre;
        return DVec3(p.omega.y * r.z - p.omega.z * r.y, p.omega.z * r.x - p.omega.x * r.z, p.omega.x * r.y - p.omega.y * r.x).toFloat(); // w x r
    }
    return surfaceVel(*nearest, feet);
}

Vector3 Player::surfaceVel(const Planetoid& p, Vector3 s) const { return Frames_SurfaceVelocity(p, fromSpace(s)).toFloat(); }
Vector3 Player::localOf(const Planetoid& p, Vector3 s) const { return p.toLocalD(fromSpace(s)).toFloat(); }

void Player::placeAt(Scene& world, Planetoid& body, Vector3 renderFeet) {
    ++spawns;
    position = world.renderOrigin + renderFeet;
    velocity = Frames_SurfaceVelocity(body, position).toFloat();
    grounded = false;
    groundBody = nullptr;
    nearest = &body; // (#2000) the loading screen's hold uses it (spawnOn)
    flyMode = false;
    anchorTo(body);
}

void Player::poisonCollideCache() {
    // At the player's own feet (body-local), so the "moved since" term is small and every sphere would be skipped.
    const Vector3 here = nearest ? nearest->toLocalD(position).toFloat() : Vector3{};
    for (int si = 0; si < 8; ++si) { collideCache.valid[si] = true; collideCache.dist[si] = 1000.0f; collideCache.centre[si] = here; }
    collideCache.body = nearest; collideCache.edits = nearest ? nearest->edits.get() : nullptr;
}

#ifdef RAYLIBGAME_UNIT_TESTS
// ---- Unit tests (#869): the multiplayer state (netState / applyNetState) in the body's frame ----
#include "unit_helpers.h"

namespace {
struct PlayerNetTest : ::testing::Test {
    Scene& world = UnitTest::FieldOnlyScene(); // reads renderOrigin only (zero)
    std::unique_ptr<Planetoid> body;
    void SetUp() override {
        body = std::make_unique<Planetoid>("unit", DVec3{ 0.0, 0.0, 0.0 }, Vector3{ 0, 40, 0 }, 1.0f, 42u, 60000.0f);
        body->rotation = { 0, 40, 0 };
        body->spin = { 0, 0.1f, 0 }; // ~105 m/s at the surface, as terra
        body->netId = 0x1000000000000001ull;
    }
    Player::NetState State(bool grounded) const {
        Player::NetState s;
        s.body = body->netId; s.local = { 1200.25, 59950.5, -640.125 }; s.vel = { 1.5f, grounded ? 0.0f : -3.0f, 0.25f };
        s.grounded = grounded; s.posture = 1; s.jumpBuffer = 0.1f; s.coyote = 0.05f; s.jumpLock = 0.02f; s.jetpackFuel = 0.6f;
        return s;
    }
    static void ExpectSame(const Player::NetState& a, const Player::NetState& b) {
        EXPECT_EQ(a.body, b.body);
        EXPECT_LT((a.local - b.local).length(), 1e-6) << "local position";
        EXPECT_LT(Vector3Distance(a.vel, b.vel), 2e-3f) << "local velocity " << b.vel.x << " " << b.vel.y << " " << b.vel.z;
        EXPECT_EQ(a.grounded, b.grounded);
        EXPECT_EQ(a.posture, b.posture);
        EXPECT_FLOAT_EQ(a.jumpBuffer, b.jumpBuffer); EXPECT_FLOAT_EQ(a.coyote, b.coyote);
        EXPECT_FLOAT_EQ(a.jumpLock, b.jumpLock); EXPECT_FLOAT_EQ(a.jetpackFuel, b.jetpackFuel);
    }
};
}

TEST_F(PlayerNetTest, AppliedStateReadsBack) {
    for (bool grounded : { true, false }) {
        Player p;
        const Player::NetState s = State(grounded);
        ASSERT_TRUE(p.applyNetState(world, s));
        Player::NetState back;
        ASSERT_TRUE(p.netState(world, back));
        ExpectSame(s, back);
    }
}

// The same through the frame conversions: a player not stepped since a spawn reports what it stands on now.
TEST_F(PlayerNetTest, UnsteppedStateGoesThroughTheBodyFrame) {
    for (bool grounded : { true, false }) {
        Player p;
        const Player::NetState s = State(grounded);
        ASSERT_TRUE(p.applyNetState(world, s));
        ++p.spawns; // the cached state is a spawn old: netState converts from the anchors and the world velocity
        Player::NetState back;
        ASSERT_TRUE(p.netState(world, back));
        ExpectSame(s, back);
    }
}

// #856: a joiner the server never steps (still loading) keeps its spawn-time world position while the planet turns; its
// reported state must stay where it was put on the planet, not slide by spin x time (12 km in #856's run).
TEST_F(PlayerNetTest, UnsteppedPlayerStaysPutOnATurningPlanet) {
    for (bool grounded : { true, false }) {
        Player p;
        const Player::NetState s = State(grounded);
        ASSERT_TRUE(p.applyNetState(world, s));
        ++p.spawns;
        body->rotation.y += 60.0f * body->spin.y; // a minute of spin, the player not stepped
        Player::NetState later;
        ASSERT_TRUE(p.netState(world, later));
        EXPECT_LT((later.local - s.local).length(), 1e-6) << (grounded ? "grounded" : "airborne") << " player slid on the body";
    }
}

TEST_F(PlayerNetTest, PlacedPlayerStaysPutOnATurningPlanet) {
    Player p;
    const DVec3 local{ -800.0, 59990.0, 300.0 };
    p.placeAt(world, *body, body->toWorldD(local).toFloat()); // render origin = world origin
    body->rotation.y += 90.0f * body->spin.y;
    Player::NetState s;
    ASSERT_TRUE(p.netState(world, s));
    EXPECT_LT((s.local - local).length(), 0.02) << "float render feet: within 2 cm at 60 km";
}

// #2149: a second body's motion at a contact, in the player's velocity terms. Two ring rocks on near orbits, 3 km/s each:
// ground-relative (grounded, or the frame path) the contact sees only their difference, never the orbit itself.
TEST(PlayerContact, OtherBodyMotionIsInThePlayersTerms) {
    Planetoid a("rockA", DVec3{ 0.0, 0.0, 0.0 }, Vector3{ 0, 0, 0 }, 1.0f, 7u, 20.0f);
    Planetoid b("rockB", DVec3{ 45.0, 0.0, 0.0 }, Vector3{ 0, 0, 0 }, 1.0f, 8u, 20.0f);
    a.linearVelocity = { 0, 0, 3000.0f };
    b.linearVelocity = { 0, 0.5f, 3001.0f };
    a.spin = b.spin = { 0, 0, 0 };
    const Vector3 at = { 22.5f, 0, 0 }; // (render origin 0: the world point too)
    const Vector3 rel = ContactFrameVelocity(b, a, DVec3(at), Vector3{ 0, 0, 0 }); // ground-relative velocity
    EXPECT_LT(Vector3Distance(rel, Vector3{ 0, 0.5f, 1.0f }), 1e-3f) << rel.x << " " << rel.y << " " << rel.z;
    const Vector3 world = ContactFrameVelocity(b, a, DVec3(at), a.surfaceVelocity(at)); // world-frame (airborne, frames off)
    EXPECT_LT(Vector3Distance(world, b.surfaceVelocity(at)), 1e-3f);
}

TEST_F(PlayerNetTest, UnknownBodyIsRefused) {
    Player p;
    Player::NetState s = State(true);
    s.body = 0x1000000000000999ull;
    EXPECT_FALSE(p.applyNetState(world, s));
    Player::NetState out;
    EXPECT_FALSE(Player{}.netState(world, out)) << "no body yet";
}
#endif
