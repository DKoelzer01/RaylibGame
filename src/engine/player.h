#ifndef PLAYER_H
#define PLAYER_H

#include <raylib.h>
#include "dvec3.h"
#include "frames.h"
#include <cstdint>
#include <algorithm>

class Scene;
class Planetoid;
struct PlayerInput;

// First-person character for small bodies:
// - gravity is summed over every body (GravityAt) and "up" follows it, so you can walk all
//   the way around an asteroid;
// - collision is resolved directly against the density field (two spheres: feet and head), so
//   dug tunnels and built terrain are walkable with no separate collision mesh;
// - standing on a spinning body carries you with it;
// - jump, and a jetpack while airborne (Space up, Ctrl down).
inline double gPlayerCollideSeconds = 0.0; // (#768) time in resolveCollisions since the start (the player-cost bench)

struct Player {
    DVec3 position{};                 // feet, world space (double: survives the floating origin)
    Vector3 velocity{};               // world units / s
    Vector3 forward = { 0, 0, 1 };    // view direction
    Vector3 up = { 0, 1, 0 };
    bool grounded = false;
    bool flyMode = false;             // F: free-fly dev camera (no gravity/collision)
    bool flyAllowed = true;           // F works (#531: creative mode only; the game sets it each frame); off: leaves fly mode
    bool doubleTapFly = false;        // (#1562) a double-tap of Space toggles flying (creative; the game sets it; benches don't)
    float sinceJumpTap = 1e9f;        // (#1562) s since the last Space press: a second within kDoubleTapFly toggles flying
    static constexpr float kDoubleTapFly = 0.3f;
    int spawns = 0;                   // counts spawnOn/placeAt: the game raises the loading screen on a new one
    float jetpackFuel = 1.0f;         // 0..1, recharges on the ground
    // (#1105) Ladders: touching one (the capsule overlaps a ladder's cell) turns gravity off; W climbs, S descends, nothing
    // holds; Space pushes off. climbInput / climbJump are this frame's input, set by update() for the physics steps.
    bool onLadder = false;
    float climbInput = 0.0f;          // +1 W, -1 S, 0 holds
    bool climbJump = false;
    float ladderCooldown = 0.0f;      // s after pushing off before a ladder catches again
    static constexpr float kClimbSpeed = 2.0f;     // m/s up or down a ladder
    static constexpr float kLadderPushOff = 3.0f;  // m/s away from it (and 2 up) on Space
    // Flight pack (#163): H toggles hover - height held, no drift; Space/Ctrl climb/descend, Shift boosts. It runs
    // on suit power (survival drains it, hoverPowered off at zero) and switches itself on in a fast fall.
    bool canHover = false, hovering = false, hoverPowered = true;
    float hoverSpeed = 6.0f, hoverClimb = 5.0f;
    bool jetpackFiring = false;       // fuel was burnt this update (survival / audio / HUD)
    Planetoid* nearest = nullptr;     // body we collide with (gravity comes from every body)
    float gravity = 0.0f;             // current acceleration magnitude (HUD/debug)
    bool autoWalk = false;            // tests: walk forward (main.cpp turns it into PlayerInput: W held, #312)
    float gravityScale = 1.0f;        // spells/effects (e.g. featherfall)
    bool dampeners = true;            // B: inertial dampeners (deep space only)
    bool autoJump = false;            // tests: jump once on the next update (a Space press in PlayerInput)
    bool autoCrouch = false;          // tests: Ctrl held (crouch, #578)
    bool autoSprint = false;          // tests: Shift held with autoWalk (#338: the explore test walks 600 m)
    float lastContactDot = 0.0f;      // debug: dot(contact normal, up) of the last collision
    float impactSpeed = 0.0f;         // hardest collision speed this frame (any contact: tests, launches)
    float landingSpeed = 0.0f;        // speed into the floor at touchdown this frame (fall damage, #263)
    // Plan 11 F2a (#667, #1146): with RAYLIBGAME_FRAMES=1, the player's state in its body frame (body NetId, Fixed below
    // the line, Centred above; double position and the frame velocity: in Fixed the speed against the ground, airborne
    // too), refreshed at the end of every update. Observed only for now (tests, net); F2b makes it the state physics runs on.
    FramePos frame{};
    // (F3, #1277) The velocity against the ground of the body the player is near, for readers (HUD, perf log, benches):
    // on the frame path that is the velocity itself, airborne too; else airborne velocities are world-frame and take the
    // surface's motion out (grounded ones are surface-relative already). Fly mode has no nearest: the nearest body is found.
    Vector3 groundVelocity(Scene& world) const;
    Vector3 groundNormal = { 0, 1, 0 }; // contact normal while grounded
    bool inWater = false;             // the body's middle is in water (swimming: buoyancy, drag, Space/Ctrl swim)
    Vector3 riverCurrent{ 0, 0, 0 };  // in a river: its current (world m/s, #361): carries a swimmer, pushes a wader
    bool swimUp = false, swimDown = false;
    float eyeOffset = 0.0f;           // camera height smoothing (metres along up, see update)
    float eyeVz = 0.0f;               // low-passed vertical speed of the feet in the ground frame
    // The same sideways (#255): walking into jagged walls, the collision push-out moved the feet back and forth by
    // centimetres each frame. Body-local, at most 3 cm, so starting and stopping never visibly trail.
    Vector3 eyeOffsetH = { 0, 0, 0 }, eyeVh = { 0, 0, 0 };
    uint64_t eyeOffsetBody = 0;
    bool eyeSmoothing = true;         // tests: --no-eye-smooth turns it off (A/B); the Camera smoothing setting Off too
    float eyeSmoothScale = 1.0f;      // Camera smoothing (#168): Low 0.5, Default 1, High 2.4 (time constant and leash scale)
    bool speedVignette = false;       // Speed vignette setting (#168): darken the screen edges when moving or turning fast
    float eyeStep = 0.0f;             // tests: the eye's height change in the last update (ground frame, m)
    bool eyeStepKnown = false;

    // Tunables
    float eyeHeight = 1.6f;
    float height = 1.7f;     // collision capsule top (#255): 2 m doorways pass
    float radius = 0.45f;
    // Posture (#578): crouching (hold Ctrl on the ground) and prone (Z) shorten the capsule and slow the walk; getting up
    // needs the room above. `height` and the eye ease follow the posture; standHeight / standEye are the upright ones.
    enum class Posture { Stand, Crouch, Prone };
    Posture posture = Posture::Stand;
    bool proneWanted = false;
    float standHeight = 1.7f, standEye = 1.6f;
    static constexpr float kCrouchHeight = 1.25f, kCrouchEye = 1.15f, kProneHeight = 0.9f, kProneEye = 0.6f;
    static constexpr float kCrouchSpeed = 0.5f, kProneSpeed = 0.25f; // of the walk
    float postureHeight(Posture p) const { return p == Posture::Prone ? kProneHeight : p == Posture::Crouch ? kCrouchHeight : standHeight; }
    float postureEye(Posture p) const { return p == Posture::Prone ? kProneEye : p == Posture::Crouch ? kCrouchEye : standEye; }
    bool roomAbove(Scene& world, float toHeight) const; // nothing solid between the capsule's top now and toHeight
    float walkSpeed = 6.0f;
    float statusSpeed = 1.0f;         // (#1215) status effects on the walk: chilled 0.7, frozen 0 (Combat_Step sets it)
    // (#1218) A dodge under way (Combat_Step starts it): a fixed step along dodgeDir (world, level) for dodgeLeft s.
    Vector3 dodgeDir = { 0, 0, 0 };
    float dodgeLeft = 0.0f, dodgeSpeed = 0.0f;
    DVec3 dodgeOrbit{};             // (#1218) a sideways dodge while locked on: circles this point (world)
    float dodgeOrbitSign = 0.0f;     //   +1 / -1 the way round (n x (here - target)), 0 a straight dodge
    bool dodgeOrbitDone = false;     //   this dodge was looked at (animals.cpp)
    float runMultiplier = 2.0f;
    // Knockback (#165): a push relative to the ground, applied at the next physics step. On the ground it lifts
    // you off it like a jump does (grounded velocity is overwritten by walking each step); in the air it adds.
    Vector3 pendingPush = { 0, 0, 0 };
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
    // Multiplayer (#327): everything a step depends on, in the nearest body's frame (time-independent: the server and a
    // predicting client step the same input at different times). vel is relative to the surface under the feet, on body
    // axes. applyNetState sets the player to it with the ground or air anchor update() expects, so the next step goes on
    // from there as if it had got there itself (setting position alone dropped the ground contact for a step).
    struct NetState {
        uint64_t body = 0;              // NetId
        DVec3 local{};                  // feet
        Vector3 vel{};
        bool grounded = false, hovering = false;
        uint8_t posture = 0;            // Posture
        float jumpBuffer = 0.0f, coyote = 0.0f, jumpLock = 0.0f, jetpackFuel = 1.0f;
    };
    bool netState(const Scene& world, NetState& s) const; // false: no body yet
    bool applyNetState(Scene& world, const NetState& s);  // false: the body isn't loaded
    void update(Scene& world, const PlayerInput& in, float dt); // one simulation step (SimDt, #315) with that step's input (#312)
    // Mouse look, once per frame (not per simulation step, so the view turns smoothly at any frame rate; #315).
    void look(Vector2 mouseDelta);
    // The feet to draw from (#315): between the last two simulation steps by alpha (SimAlpha), in the frame of the body
    // they were on, so walking looks smooth while the simulation steps at 60 Hz. `position` when nothing to blend.
    DVec3 renderPosition(float alpha) const;
    void applyToCamera(Scene& world) const;  // camera at the eye (renderPosition), looking along forward
    // (F4, #1412) the eye applyToCamera last placed the camera at, in absolute world doubles (the camera's float rounds it).
    // Aim rays start here with RAYLIBGAME_FRAMES=1 (AimOrigin): a still player's aim then holds still.
    mutable DVec3 lastEyeD{};
    // The aim ray for a body: its body-local double origin and body-local direction. With the switch on, the eye above in
    // double and `forward`; off (or flying), the render camera as before. dirWorld: the same direction in world axes.
    void aimRay(const Scene& world, const Planetoid& body, DVec3& originLocal, Vector3& dirLocal, Vector3& dirWorld) const;
    mutable float lastEyePush = 0.0f;         // how far applyToCamera pushed the eye out of the terrain last frame (m; #255 tests)
    mutable float eyeDrop = 0.0f;              // how far the eye is lowered toward the head sphere's centre to keep it out of the terrain (#255)
    // Fly mode follows the reference body (the nearest whose sphere of influence holds the camera):
    // its spin and motion within 2 radii, its motion only further out. recordFly() runs after the
    // camera's own input each frame (Scene::drawScene); update() then carries it with the body.
    bool flyFollow = true;
    void recordFly(Scene& world);
    Planetoid* flyReference(Scene& world, bool& rotate) const;

    // Tests (#768): make every collision-cache entry claim "clear by a kilometre" - a stale cache's worst case. A player
    // whose physics depends on history would then sink through the ground; the determinism check steps it beside a copy.
    void poisonCollideCache();
private:
    // Render interpolation (#315): the feet after the last two steps, local to interpBody; interpSpawns detects teleports.
    uint64_t interpBody = 0;
    DVec3 interpPrev{}, interpCur{};
    int interpSpawns = -1;
    DVec3 groundLocal{};              // feet in the ground body's local space, double (follows its spin exactly)
    Planetoid* groundBody = nullptr;
    uint64_t groundId = 0;            // validates groundBody (bodies can be destroyed, e.g. world regen)

    void updateFrame(Scene& world);
    Quaternion groundRotPrev = { 0, 0, 0, 1 }; // ground body's rotation at the end of last update
    DVec3 lastEndPosition{};                    // where update() left us (detects outside moves)
    // Airborne near a body: carried with its motion each frame and moved by velocity relative to its
    // surface, so substep collisions line up with the (already moved) terrain. Same inertial path to
    // first order; without it you ran up to a frame of spin (86 m/s at 60 km) ahead of the ground.
    Planetoid* airBody = nullptr;
    uint64_t airId = 0;
    DVec3 airLocal{};
    // Plan 11 F2b (#1146): the frame path is in use (RAYLIBGAME_FRAMES=1, a Fixed frame of the nearest body), and an
    // airborne frame velocity kept in that body's axes between updates (direction, length, the body's runtime id).
    bool fxActive = false;
    bool movedOutside = false; // (F2d) moved from outside since the last update (spawn, placeAt, respawn, load)
    double viewR[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } }; // F2c: the frame body's rotation at the end of the last update
    uint64_t viewBody = 0;                                            // its runtime id (0: no view carry yet)
    Vector3 frameGravity(Scene& world, Vector3 feet) const; // the frame's acceleration (felt gravity + spin forces), world axes
    Vector3 airVelDir{};
    float airVelLen = 0.0f;
    uint64_t airVelBody = 0;
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
    // (Plan 11 F4, #1412, RAYLIBGAME_FRAMES=1) beyond 2 radii the camera's frame is centred on its body (it moves with it,
    // without its spin): the camera as a double offset from the body's centre, not a float delta added to a float camera
    // every frame (which accumulated rounding as drift).
    DVec3 flyOffset{};
    uint64_t flyBodyIdPrev = 0;       // the body flyOffset was kept against
    void stepPhysics(Scene& world, float dt, Vector3 wishMove, bool jetUp, bool jetDown);
    float jumpBuffer = 0.0f; // seconds a jump press stays pending
    float coyote = 0.0f;     // seconds a jump is still allowed after leaving the ground
    float jumpLock = 0.0f;   // seconds after a jump during which the player cannot re-ground
    bool resolveCollisions(Scene& world, Vector3& feet);
    // #768: per capsule sphere, the last measured terrain surface distance (m) and where it was measured (body-local).
    // Within one update (cleared at its start: physics must not depend on history - client replays and the server must
    // agree), a sphere whose distance, less how far it has moved since, is still beyond its radius + 5 cm can't touch: it is not
    // measured again (each measurement is 7 density samples). Any edit (a new snapshot) or another body clears it.
    struct CollideCache { Vector3 centre[8]{}; float dist[8]{}; bool valid[8]{}; const void* body = nullptr; const void* edits = nullptr; } collideCache;
    void collideTrunks(Vector3& feet); // tree trunks: vertical capsules on the nearest body
    // Velocity of the terrain being collided with, in the frame `velocity` is in: zero while grounded
    // (velocity is then relative to the surface), the surface velocity while airborne (world frame).
    Vector3 collideFrameVel = { 0, 0, 0 };
    // Multiplayer (#327, #750): the body-local state at the end of the last step (netState reads it later).
    NetState netCache;
    bool netCached = false;
    int netCacheSpawns = -1;
    void CacheNetState(const Scene& world);
};

extern Player& player;         // the local PlayerSession's (src/space/playersession.h, #313)

#endif
