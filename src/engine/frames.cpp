#include "frames.h"
#include "world.h"
#include "worldtypes.h"
#include "bodies.h"
#include <algorithm>
#include <cmath>

// Plan 11 (#667). See frames.h.

// Read on first use, not in a static initialiser: a scenario's s.env is set inside main() (Scenario_PrepareArgs), after
// static initialisers ran - read there, every "-pose" scenario ran with the switch off.
static int gFramesOn = -1; // -1: not read yet
bool Frames_On() {
    if (gFramesOn < 0) { const char* e = getenv("RAYLIBGAME_FRAMES"); gFramesOn = (e && e[0] == '1') ? 1 : 0; }
    return gFramesOn == 1;
}
void Frames_SetOn(bool on) { gFramesOn = on ? 1 : 0; }

namespace {
constexpr double kInf = 1e300;
DVec3 Cross(const DVec3& a, const DVec3& b) { return { a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x }; }
double Dot(const DVec3& a, const DVec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
DVec3 Mul(const double R[3][3], const DVec3& v) { // R v: body-local -> global axes
    return { R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z, R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z, R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z };
}
DVec3 MulT(const double R[3][3], const DVec3& v) { // R^T v: global axes -> body-local (R is orthonormal)
    return { R[0][0] * v.x + R[1][0] * v.y + R[2][0] * v.z, R[0][1] * v.x + R[1][1] * v.y + R[2][1] * v.z, R[0][2] * v.x + R[1][2] * v.y + R[2][2] * v.z };
}
double SurfaceRadius(const Planetoid& b) { return (double)b.radius * b.scale; }
bool OnOrbit(const Planetoid& b) { return b.orbitParent && Planetoid::find(b.orbitParentId) == b.orbitParent; }
}

BodyPose Frames_Pose(const Planetoid& b) {
    BodyPose p;
    p.centre = b.worldPos;
    BodyRotationAndSpin(b, p.R, p.omega);
    if (OnOrbit(b)) {
        // A circular orbit from the clock in the parent's XZ plane (PlaceBodyOnOrbit): its velocity and its acceleration
        // (the felt-gravity term of anything in its frame) are exact functions of the offset from the parent.
        const BodyPose q = Frames_Pose(*b.orbitParent);
        const DVec3 rel = b.worldPos - q.centre;
        const DVec3 inPlane{ rel.x, 0.0, rel.z };
        const double w = b.orbitRate;
        p.vel = q.vel + DVec3{ -inPlane.z * w, 0.0, inPlane.x * w };
        p.acc = q.acc - inPlane * (w * w);
    } else {
        // Root bodies, and bodies moved by other code (stations: Kepler elements): their velocity as kept. A free-falling
        // one (freeFall: stations) accelerates as its code moves it (freeFallAcc), so anything in its frame feels only the
        // difference, not the fall a second time (Plan 11 F5, #1457: stations as Fixed(station) frames).
        p.vel = DVec3{ b.linearVelocity.x, b.linearVelocity.y, b.linearVelocity.z };
        // Its own code's acceleration, not this file's gravity model: a station flies a two-body Kepler orbit, and the
        // all-bodies model differs from it (Luna's pull, the inside model): 1.4 m of drift in 30 s beside one
        // (rocket/station-hold), since a framed rocket is drawn at the station's real place plus its offset.
        if (b.freeFall) p.acc = b.freeFallAcc;
    }
    return p;
}

bool Frames_ToGlobal(const FramePos& f, DVec3& pg, DVec3& vg) {
    const Planetoid* b = Planetoid::findNet(f.frame.body);
    if (!b) return false;
    const BodyPose p = Frames_Pose(*b);
    if (f.frame.kind == FrameKind::Centred) { pg = p.centre + f.pos; vg = p.vel + f.vel; return true; }
    const DVec3 rx = Mul(p.R, f.pos);
    pg = p.centre + rx;
    vg = p.vel + Cross(p.omega, rx) + Mul(p.R, f.vel);
    return true;
}

bool Frames_FromGlobal(const DVec3& pg, const DVec3& vg, FrameRef to, FramePos& out) {
    const Planetoid* b = Planetoid::findNet(to.body);
    if (!b) return false;
    const BodyPose p = Frames_Pose(*b);
    out.frame = to;
    const DVec3 d = pg - p.centre;
    if (to.kind == FrameKind::Centred) { out.pos = d; out.vel = vg - p.vel; return true; }
    out.pos = MulT(p.R, d);
    out.vel = MulT(p.R, vg - p.vel - Cross(p.omega, d));
    return true;
}

bool Frames_Convert(const FramePos& f, FrameRef to, FramePos& out) {
    if (f.frame == to) { out = f; return true; }
    DVec3 pg, vg;
    return Frames_ToGlobal(f, pg, vg) && Frames_FromGlobal(pg, vg, to, out);
}

DVec3 Frames_GravityGlobal(const std::vector<Planetoid*>& bodies, const DVec3& pg) {
    DVec3 g{ 0, 0, 0 };
    for (const Planetoid* b : bodies) {
        if (!b) continue;
        if (b->flat) { g.y -= b->surfaceGravity; continue; } // open-world terrain: uniform, straight down
        const DVec3 toC = b->worldPos - pg;
        const double d = toC.length();
        if (d < 1e-4) continue;
        const double R = b->gravityRadius > 0.0f ? (double)b->gravityRadius : SurfaceRadius(*b);
        const double a = d >= R ? b->surfaceGravity * (R * R) / (d * d) : b->surfaceGravity * d / R; // as BodyGravityAt
        g += toC * (a / d);
    }
    return g;
}

DVec3 Frames_Acceleration(const std::vector<Planetoid*>& bodies, const FramePos& f, bool spinForces) {
    const Planetoid* b = Planetoid::findNet(f.frame.body);
    if (!b) return DVec3{ 0, 0, 0 };
    const BodyPose p = Frames_Pose(*b);
    DVec3 pg, vg;
    Frames_ToGlobal(f, pg, vg);
    const DVec3 felt = Frames_GravityGlobal(bodies, pg) - p.acc; // the frame's own fall taken out (felt gravity, exact)
    if (f.frame.kind == FrameKind::Centred) return felt;
    DVec3 a = MulT(p.R, felt);
    if (spinForces) {
        const DVec3 w = MulT(p.R, p.omega);
        a -= Cross(w, f.vel) * 2.0;           // Coriolis
        a -= Cross(w, Cross(w, f.pos));       // centrifugal (a constant spin: no Euler term)
    }
    return a;
}

double Frames_AltitudeLine(const Planetoid& b) {
    switch (b.kind) {
        case BodyKind::Planet: return 35000.0; // the top of the air (rocket.cpp kAirTop): orbits (periapsis above it) are Centred
        case BodyKind::Moon: return std::max(2000.0, 0.2 * SurfaceRadius(b));
        default: return kInf;                   // asteroids, ring rocks: Fixed over their whole SOI
    }
}

double Frames_Soi(const Planetoid& b) {
    const double R = SurfaceRadius(b);
    if (OnOrbit(b)) { // (Player::flyReference's rule)
        const Planetoid* m = b.orbitParent;
        const double Rm = SurfaceRadius(*m);
        const double ratio = (b.surfaceGravity * R * R) / std::max(1e-9, (double)m->surfaceGravity * Rm * Rm);
        return std::max(2.0 * R, b.orbitRadius * std::pow(std::max(ratio, 0.0), 0.4));
    }
    if (b.kind == BodyKind::Asteroid) return 10.0 * R; // loose asteroids (no orbit): their neighbourhood
    return kInf;                                        // root bodies reach everywhere
}

FrameRef Frames_Choose(const std::vector<Planetoid*>& bodies, const DVec3& pg, FrameRef current, bool contact) {
    const Planetoid* best = nullptr;
    double bestSoi = kInf * 10.0, bestGap = kInf;
    for (const Planetoid* b : bodies) {
        if (!b || (b->vehicle && !b->freeFall) || b->netId == 0) continue; // a flying rocket isn't a frame; a station (also a vehicle body, EmptyField) is (F5)
        const double d = (pg - b->worldPos).length();
        const double soi = Frames_Soi(*b);
        const double reach = soi >= kInf ? kInf : soi * (b->netId == current.body ? 1.03 : 0.97);
        if (d > reach) continue;
        const double gap = d - SurfaceRadius(*b);
        // the smallest SOI that holds the point (a child before its parent); among equals the nearest surface
        if (soi < bestSoi || (soi == bestSoi && gap < bestGap)) { best = b; bestSoi = soi; bestGap = gap; }
    }
    if (!best) return current;
    FrameRef out{ best->netId, FrameKind::Fixed };
    if (contact) return out;
    const double line = Frames_AltitudeLine(*best);
    if (line >= kInf) return out;
    const double alt = (pg - best->worldPos).length() - SurfaceRadius(*best);
    if (best->netId == current.body)
        out.kind = current.kind == FrameKind::Fixed ? (alt > 1.05 * line ? FrameKind::Centred : FrameKind::Fixed)
                                                    : (alt < 0.95 * line ? FrameKind::Fixed : FrameKind::Centred);
    else out.kind = alt < line ? FrameKind::Fixed : FrameKind::Centred;
    return out;
}

#ifdef RAYLIBGAME_UNIT_TESTS
#include <gtest/gtest.h>
#include <memory>

namespace {
// A body at a global position, turned and spinning: a planet (Terra-like) or a child on a circular orbit.
std::unique_ptr<Planetoid> Body(const char* name, uint64_t netId, DVec3 pos, float radius, Vector3 rotDeg, Vector3 spinDeg, BodyKind kind) {
    auto b = std::make_unique<Planetoid>(name, pos, rotDeg, 1.0f, 42u, radius);
    b->rotation = rotDeg; b->spin = spinDeg; b->netId = netId; b->kind = kind;
    b->surfaceGravity = 5.886f; b->gravityRadius = radius; // (the constructor takes it from its terrain field)
    return b;
}
double Dist(const DVec3& a, const DVec3& b) { return (a - b).length(); }
}

// Round trips through every frame keep position and velocity to sub-micrometre / nano-m/s at planet scale.
TEST(FramesCore, ConversionsRoundTrip) {
    auto terra = Body("terra", 0x9001, { 1000.0, -2000.0, 500.0 }, 60000.0f, { 10.0f, 237.25f, -4.0f }, { 0.0f, 0.1f, 0.0f }, BodyKind::Planet);
    const FrameRef fixed{ 0x9001, FrameKind::Fixed }, centred{ 0x9001, FrameKind::Centred };
    FramePos f{ fixed, { 1234.5, 59000.25, -812.75 }, { 3.0, -9.5, 1.25 } };
    FramePos c, back;
    ASSERT_TRUE(Frames_Convert(f, centred, c));
    ASSERT_TRUE(Frames_Convert(c, fixed, back));
    EXPECT_LT(Dist(back.pos, f.pos), 1e-6);
    EXPECT_LT(Dist(back.vel, f.vel), 1e-9);
    // Centred = global minus the centre (same axes); a point at rest on the ground moves at w x r there.
    FramePos rest{ fixed, f.pos, { 0, 0, 0 } };
    FramePos rc;
    ASSERT_TRUE(Frames_Convert(rest, centred, rc));
    const DVec3 wr = Cross(Frames_Pose(*terra).omega, rc.pos); // w x r (the spin's axis is tilted by the body's X angle)
    EXPECT_LT(Dist(rc.vel, wr), 1e-9 * wr.length() + 1e-12);
    EXPECT_NEAR(Frames_Pose(*terra).omega.length(), 2.0 * PI / 3600.0, 1e-9); // Terra's day: an hour
    EXPECT_NEAR(Dot(rc.vel, rc.pos), 0.0, 1e-6);
}

// The fictitious terms: at rest on the equator the centrifugal acceleration is w^2 r outward (Terra: ~0.183 m/s^2);
// moving east at 10 m/s adds Coriolis 2 w v, outward too.
TEST(FramesCore, SpinForcesInTheFixedFrame) {
    auto terra = Body("terra", 0x9002, { 0.0, 0.0, 0.0 }, 60000.0f, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.1f, 0.0f }, BodyKind::Planet);
    const std::vector<Planetoid*> bodies{ terra.get() };
    const double w = 2.0 * PI / 3600.0, r = 60000.0;
    FramePos rest{ { 0x9002, FrameKind::Fixed }, { r, 0.0, 0.0 }, { 0, 0, 0 } };
    const DVec3 a0 = Frames_Acceleration(bodies, rest);
    const DVec3 g0 = Frames_Acceleration(bodies, rest, false);
    EXPECT_NEAR(g0.x, -5.886, 1e-3);                 // gravity alone: straight down
    EXPECT_NEAR(a0.x - g0.x, w * w * r, 1e-6);       // + centrifugal outward
    EXPECT_NEAR(w * w * r, 0.1827, 1e-3);
    // Moving along +z (east for a spin about +Y? the sign is the cross product's: check Coriolis = -2 w x v)
    FramePos moving = rest;
    moving.vel = { 0.0, 0.0, 10.0 };
    const DVec3 a1 = Frames_Acceleration(bodies, moving);
    const DVec3 wv = Cross(DVec3{ 0, w, 0 }, moving.vel) * -2.0;
    EXPECT_NEAR(a1.x - a0.x, wv.x, 1e-9);
    EXPECT_NEAR(std::fabs(a1.x - a0.x), 2.0 * w * 10.0, 1e-9);
}

// A moon on a circular orbit: its frame's acceleration is its orbit's (rate^2 r toward the parent), so a thing falling
// with it feels only the moon (and the parent's tide), and its pose velocity is the orbit's.
TEST(FramesCore, OrbitingFramePose) {
    auto terra = Body("terra", 0x9003, { 0.0, 0.0, 0.0 }, 60000.0f, { 0, 0, 0 }, { 0.0f, 0.1f, 0.0f }, BodyKind::Planet);
    auto luna = Body("luna", 0x9004, { 560000.0, 0.0, 0.0 }, 10000.0f, { 0, 0, 0 }, { 0, 0, 0 }, BodyKind::Moon);
    luna->orbitParent = terra.get(); luna->orbitParentId = terra->id;
    luna->orbitRadius = 560000.0; luna->orbitRate = std::sqrt(5.886 * 60000.0 * 60000.0 / std::pow(560000.0, 3));
    const BodyPose p = Frames_Pose(*luna);
    EXPECT_NEAR(p.vel.z, luna->orbitRate * 560000.0, 1e-9);
    EXPECT_NEAR(p.acc.x, -luna->orbitRate * luna->orbitRate * 560000.0, 1e-12);
    // In Centred(luna) at its centre's far side 20 km out, the acceleration is the moon's pull plus Terra's tide only.
    const std::vector<Planetoid*> bodies{ terra.get(), luna.get() };
    FramePos f{ { 0x9004, FrameKind::Centred }, { 20000.0, 0.0, 0.0 }, { 0, 0, 0 } };
    const DVec3 a = Frames_Acceleration(bodies, f);
    const double moon = 5.886 * 10000.0 * 10000.0 / (20000.0 * 20000.0);
    const double GM = 5.886 * 60000.0 * 60000.0; // Terra's pull at the point less at the moon's centre (the frame's fall): the tide, outward
    const double tide = GM / (560000.0 * 560000.0) - GM / (580000.0 * 580000.0);
    EXPECT_NEAR(a.x, -moon + tide, 1e-6); // (surfaceGravity is a float: 5.886f)
}

// Frame choice: the moon's SOI holds a point near it (Centred above its 2 km line, Fixed below); Terra's frame elsewhere,
// Fixed under 35 km, Centred above; hysteresis keeps a frame across its boundary's band.
TEST(FramesCore, ChooseFramesWithHysteresis) {
    auto terra = Body("terra", 0x9005, { 0.0, 0.0, 0.0 }, 60000.0f, { 0, 0, 0 }, { 0.0f, 0.1f, 0.0f }, BodyKind::Planet);
    auto luna = Body("luna", 0x9006, { 560000.0, 0.0, 0.0 }, 10000.0f, { 0, 0, 0 }, { 0, 0, 0 }, BodyKind::Moon);
    luna->orbitParent = terra.get(); luna->orbitParentId = terra->id; luna->orbitRadius = 560000.0; luna->orbitRate = 1e-4;
    const std::vector<Planetoid*> bodies{ terra.get(), luna.get() };
    const double soi = Frames_Soi(*luna);
    EXPECT_GT(soi, 2.0 * 10000.0);
    const FrameRef none{};
    EXPECT_EQ(Frames_Choose(bodies, { 60100.0, 0.0, 0.0 }, none), (FrameRef{ 0x9005, FrameKind::Fixed }));
    EXPECT_EQ(Frames_Choose(bodies, { 100000.0, 0.0, 0.0 }, none), (FrameRef{ 0x9005, FrameKind::Centred }));
    EXPECT_EQ(Frames_Choose(bodies, { 561000.0 + 10000.0, 0.0, 0.0 }, none), (FrameRef{ 0x9006, FrameKind::Fixed }));
    EXPECT_EQ(Frames_Choose(bodies, { 560000.0 + 30000.0, 0.0, 0.0 }, none), (FrameRef{ 0x9006, FrameKind::Centred }));
    // Hysteresis at Terra's 35 km line: 36 km keeps Fixed (switches above 36.75 km), 34 km keeps Centred (below 33.25 km).
    const FrameRef tf{ 0x9005, FrameKind::Fixed }, tc{ 0x9005, FrameKind::Centred };
    EXPECT_EQ(Frames_Choose(bodies, { 60000.0 + 36000.0, 0.0, 0.0 }, tf), tf);
    EXPECT_EQ(Frames_Choose(bodies, { 60000.0 + 37000.0, 0.0, 0.0 }, tf), tc);
    EXPECT_EQ(Frames_Choose(bodies, { 60000.0 + 34000.0, 0.0, 0.0 }, tc), tc);
    EXPECT_EQ(Frames_Choose(bodies, { 60000.0 + 33000.0, 0.0, 0.0 }, tc), tf);
    // Contact forces Fixed; the moon's SOI keeps its frame out to 1.03 x, is entered only inside 0.97 x.
    EXPECT_EQ(Frames_Choose(bodies, { 60000.0 + 50000.0, 0.0, 0.0 }, tc, true), tf);
    const DVec3 edge{ 560000.0 + soi * 1.01, 0.0, 0.0 };
    EXPECT_EQ(Frames_Choose(bodies, edge, FrameRef{ 0x9006, FrameKind::Centred }).body, 0x9006u);
    EXPECT_EQ(Frames_Choose(bodies, edge, tc).body, 0x9005u);
}

// (Plan 11 F5, #1457) A station: an Asteroid-kind body of 4 m with no gravity of its own, 45 km over a Terra-like root.
// A point 20 m from it, coming from Centred(Terra), belongs in Fixed(station) (its SOI is 10 x 4 m); and as a free
// faller its frame accelerates with the gravity at its centre.
TEST(FramesCore, StationFrameAndFreeFall) {
    auto terra = Body("terra", 0x9101, { 0.0, 0.0, 0.0 }, 60000.0f, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.1f, 0.0f }, BodyKind::Planet);
    auto st = Body("station", 0x9102, { 0.0, 105000.0, 0.0 }, 4.0f, { 0.0f, 0.0f, 0.0f }, { 0.0f, 0.0f, 0.0f }, BodyKind::Asteroid);
    st->surfaceGravity = 0.0f;
    st->freeFall = true;
    st->vehicle = true; // as in the game: stations are EmptyField bodies, which the constructor marks as vehicles
    st->freeFallAcc = Frames_GravityGlobal({ terra.get() }, st->worldPos); // (set by station.cpp from its orbit in the game)
    const std::vector<Planetoid*> bodies{ terra.get(), st.get() };
    const DVec3 p{ 0.0, 105000.0, 20.0 };
    const FrameRef pick = Frames_Choose(bodies, p, FrameRef{ 0x9101, FrameKind::Centred });
    EXPECT_EQ(pick.body, 0x9102u);
    EXPECT_EQ((int)pick.kind, (int)FrameKind::Fixed);
    const DVec3 g = st->freeFallAcc;
    EXPECT_GT(g.length(), 1.0); // ~1.9 m/s^2 at 105 km
    EXPECT_LT(Dist(Frames_Pose(*st).acc, g), 1e-9);
}
#endif
