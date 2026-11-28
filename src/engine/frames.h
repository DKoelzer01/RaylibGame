#ifndef FRAMES_H
#define FRAMES_H
// Plan 11 (#667): body-relative physics frames. Every simulated thing near a body lives in that body's frame, in double:
//   Fixed(body)    - the body's centre and its own axes, turning with it (the voxel space): terrain, parts and water are
//                    static there, so collision sees no spin or orbit. Used on and near a body, below its altitude line.
//   Centred(body)  - the body's centre, axes parallel to the global ones (no spin), moving with the body: orbits and
//                    flight above the line.
// The global frame (the system's origin, non-rotating) is only the hub between frames and the start of drawing; the
// floating render origin never enters here. Bodies are referred to by NetId (#317), never a pointer or a runtime id.
// Frames switch at a body's altitude line (Fixed <-> Centred) and at sphere-of-influence boundaries, with hysteresis.
#include "dvec3.h"
#include <cstdint>
#include <vector>

class Planetoid;

// The switch (RAYLIBGAME_FRAMES=1): off by default until F2's gate passes; tests may set it.
bool Frames_On();
void Frames_SetOn(bool on);

enum class FrameKind : uint8_t { Fixed = 0, Centred = 1 };
struct FrameRef {
    uint64_t body = 0;                // the body's NetId
    FrameKind kind = FrameKind::Fixed;
    bool operator==(const FrameRef& o) const { return body == o.body && kind == o.kind; }
    bool operator!=(const FrameRef& o) const { return !(*this == o); }
};
struct FramePos {
    FrameRef frame;
    DVec3 pos{};                      // m, in the frame
    DVec3 vel{};                      // m/s, in the frame (in Fixed: the speed against the ground)
};

// A body's state at the current sim step, all double: centre, its velocity and acceleration (orbits from the clock),
// orientation (rows of R map body-local into the global axes) and angular velocity (global axes, rad/s).
struct BodyPose {
    DVec3 centre{}, vel{}, acc{}, omega{};
    double R[3][3] = { { 1, 0, 0 }, { 0, 1, 0 }, { 0, 0, 1 } };
};
BodyPose Frames_Pose(const Planetoid& body);
// (Plan 11 F6, #1698) A frame body's pose by NetId: a loaded body's, or a ring rock's while it is an instance (demoted): its
// orbit and spin from the clock, as when promoted, so a frame on a rock never dangles. False: no such body.
bool Frames_PoseOf(uint64_t netId, BodyPose& out);
bool Frames_RingRockPose(uint64_t netId, BodyPose& out); // the analytic pose alone (tests compare it with the live rock's)

// Conversions through the global frame: position and velocity (Fixed: x_g = c + R x, v_g = v_c + w x (R x) + R v;
// Centred: x_g = c + x, v_g = v_c + v). False if the frame's body is not loaded.
bool Frames_ToGlobal(const FramePos& f, DVec3& posGlobal, DVec3& velGlobal);
bool Frames_FromGlobal(const DVec3& posGlobal, const DVec3& velGlobal, FrameRef to, FramePos& out);
bool Frames_Convert(const FramePos& f, FrameRef to, FramePos& out);

// Gravity summed over every body (inverse-square outside, linear inside, as BodyGravityAt) at a global point, in double.
DVec3 Frames_GravityGlobal(const std::vector<Planetoid*>& bodies, const DVec3& posGlobal);
// The acceleration an unsupported thing feels in its frame: gravity (into the frame's axes) minus the frame body's own
// acceleration, plus in Fixed the fictitious terms: Coriolis -2 w x v and centrifugal -w x (w x x) (the user's decision,
// #667: spin forces for anything airborne). spinForces = false leaves them out (A/B).
DVec3 Frames_Acceleration(const std::vector<Planetoid*>& bodies, const FramePos& f, bool spinForces = true);

// The lines: a body's altitude line above its surface radius (m; Terra-like planets: the top of the air, 35 km; moons:
// max(2 km, 0.2 R); asteroids and ring rocks: their whole SOI - no Centred zone), and its sphere of influence (m from
// its centre; today's Laplace rule r = a (m / M)^0.4, at least 2 R; loose asteroids 10 R; root bodies: everywhere).
double Frames_AltitudeLine(const Planetoid& body);
double Frames_Soi(const Planetoid& body);
// The frame a global point belongs in: the smallest SOI that holds it (a child's before its parent's), Fixed below that
// body's line, Centred above. With hysteresis against `current` (a body is kept out to 1.03 x its SOI and entered inside
// 0.97 x; the line is crossed upward at 1.05 x and downward at 0.95 x); `contact` forces Fixed.
FrameRef Frames_Choose(const std::vector<Planetoid*>& bodies, const DVec3& posGlobal, FrameRef current, bool contact = false);

#endif
